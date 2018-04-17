#include "wfcp.h"
#include <time.h>
   
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp);
void print_error(const char *err);
double second();       
void debug(const char *err);       
int time_limit_expired(instance *inst);   

void mip_timelimit(CPXENVptr env, double timelimit, instance *inst);
int mip_update_incumbent(CPXENVptr env, CPXLPptr lp, instance *inst);
int compute_nocross_cut(instance *inst, int i, int j, int k, int *index, double *value);
int nocross_separation(CPXENVptr env, CPXLPptr lp, instance *inst);
int installplotfunction( CPXENVptr env, CPXLPptr lp, instance *inst);

FILE *gp;
char color[20][20];
int plt[20];
int xpos(int i, int j, int k, instance *inst) 
{ 
	return inst->xstart + i * (inst->nturbines * inst->ncables) + j * inst->ncables + k; 
}                                              
int ypos(int i, int j, instance *inst) 
{ 
	return inst->ystart + i * inst->nturbines + j; 
}                                              
int fpos(int i, int j, instance *inst) 
{ 
	return inst->fstart + i * inst->nturbines + j; 
}               
int setColor()
{
	sprintf(color[0],"#FF0000");
	sprintf(color[1],"#00FF00");
	sprintf(color[2],"#0000FF");
	sprintf(color[3],"#FFFF00");
	sprintf(color[4],"#00FFFF");
	sprintf(color[5],"#FF00FF");
	sprintf(color[6],"#C0C0C0");
	sprintf(color[7],"#800000");
	sprintf(color[8],"#808000");
	sprintf(color[9],"#008000");
	sprintf(color[10],"#800080");
	sprintf(color[11],"#008080");
	sprintf(color[12],"#000080");
	sprintf(color[13],"#808080");
	sprintf(color[14],"#000000");
	return 0;
}
int makeScript(instance *inst)
{
	FILE *s;	
	setColor();
	int flag = 0;
	s = fopen("plot/script_plot.p","w");
	fprintf(s,"set autoscale\n");
	for(int k = 0; k < inst->ncables; k++)
	{
		if(plt[k] > 0)
		{
			if(k != 0 && flag != 0)
				fprintf(s, "\nre" );	
			fprintf(s, "plot \\\n\t\'plot/plot%d.dat\' using 1:2 with lines lc rgb \"%s\" lw 2 title \"Cable %d\",\\\n",k,color[k],k+1);
			fprintf(s, "\t\'plot/plot%d.dat\' using 1:2:(0.6) with circles fill solid lc rgb \"black\" notitle,\\\n",k);
			fprintf(s, "\t\'plot/plot%d.dat\' using 1:2:3     with labels tc rgb \"black\" offset (0,0) font \'Arial Bold\' notitle",k );
			flag++;
		}
	}
	fclose(s);
	return 0;
}
int plotGraph(CPXENVptr env, CPXLPptr lp, instance *inst)
{
	double *x;
	x = (double *) calloc(CPXgetnumcols(env, lp), sizeof(double)); 
	CPXgetx(env, lp, x, 0, CPXgetnumcols(env, lp) -1);
	FILE *f;
	char filename[30];
	printf("PLOT GRAPH\n");

	for(int k = 0; k < inst->ncables; k++)
	{
		plt[k] = 0;
		sprintf(filename,"plot/plot%d.dat",k);
		f = fopen(filename,"w");
		for(int i = 0; i < inst->nturbines; i++)
		{
			for(int j = 0; j < inst->nturbines; j++)
			{
				if(x[xpos(i,j,k,inst)] > 0.5)
				{
					fprintf(f, "%lf %lf %d\n", inst->xcoord[i], inst->ycoord[i], i );
					fprintf(f, "%lf %lf %d\n\n", inst->xcoord[j], inst->ycoord[j], j );
					plt[k]++;
				}
			}
		}
		fclose(f);
	}
	makeScript(inst);
	fprintf(gp, "load 'plot/script_plot.p'\n");
	return 0;
}

int noCross(int p1, int p2, int p3, int p4, instance *inst) // p1 < - > p2 cross with p3 < - > p4 ?
{
	double lambda = 0;
	double mu = 0;
	double x1 = inst->xcoord[p1];
	double y1 = inst->ycoord[p1];

	double x2 = inst->xcoord[p2];
	double y2 = inst->ycoord[p2];

	double x3 = inst->xcoord[p3];
	double y3 = inst->ycoord[p3];

	double x4 = inst->xcoord[p4];
	double y4 = inst->ycoord[p4];
	double det = (((x2-x1)*(y4-y3)) - ((x4-x3)*(y2-y1)));
	if(fabs(det)>= EPSILON)
	{
		double dmu = (x1-x3)*(y1-y2) - (y1-y3) * (x1-x2);
		double dlambda = (x4-x3)*(y1-y3) - (y4-y3)*(x1-x3);
		mu=dmu/det;
		lambda=dlambda/det;
		if(mu>=XSMALL && mu<=1-XSMALL && lambda >= XSMALL && lambda<= 1-XSMALL)
			return 1;
	}
	return 0;
	//printf("Non C'è Crossing\n");
	
}

double dist(int i, int j, instance *inst)
{
	double dx = inst->xcoord[i] - inst->xcoord[j];
	double dy = inst->ycoord[i] - inst->ycoord[j]; 
	return sqrt(dx*dx+dy*dy);
	
}        

int is_fractional(double x) 						// it works for x in [0,1] only
{
	return ( (x > XSMALL) && (x < 1-XSMALL) );
}    

int is_all_integer(int n, const double *x) 			// it works for x_j in [0,1] only
{
	for ( int j = 0; j < n; j++ ) 
	{
		if ( is_fractional(x[j]) ) 
			return 0; 
	}
	return 1;
}                                                                                                                               
                         

int time_limit_expired(instance *inst)	 
{
	double tspan = second() - inst->tstart;
	if (  tspan > inst->timelimit ) 
	{
		if ( VERBOSE >= 100 ) 
			printf("\n\n$$$ time limit of %10.1lf sec.s expired after %10.1lf sec.s $$$\n\n", inst->timelimit, tspan);
		//exit(0); 
		return 1;
	}  
	return 0;
}
void mip_timelimit(CPXENVptr env, double timelimit, instance *inst)
{
	double residual_time = inst->tstart + inst->timelimit - second();
	if ( residual_time < 0.0 ) residual_time = 0.0;
	CPXsetintparam(env, CPX_PARAM_CLOCKTYPE, 2);
	CPXsetdblparam(env, CPX_PARAM_TILIM, residual_time); 							// real time
	CPXsetdblparam(env, CPX_PARAM_DETTILIM, TICKS_PER_SECOND*timelimit);			// ticks
}

void double_vector_copy(int n, const double *from, double *to) // vector copy
{
	for ( int j = 0; j < n; j++ ) 
		to[j] = from[j];
}

double mip_value(CPXENVptr env, CPXLPptr lp)
{
	double zz;
	if ( CPXgetobjval(env, lp, &zz) ) zz = CPX_INFBOUND;
	return zz;
}

int mip_update_incumbent(CPXENVptr env, CPXLPptr lp, instance *inst)
{
	int ncols = CPXgetnumcols(env, lp);

	int newsol = 0;

	inst->tbest = second() - inst->tstart;
	inst->zbest = mip_value(env, lp);
	CPXgetx(env, lp, inst->best_sol, 0, ncols-1);
	plotGraph(env, lp, inst);
	if ( VERBOSE >= 40 ) printf("\n >>>>>>>>>> incumbent update of value %lf at time %7.2lf <<<<<<<<\n", inst->zbest, inst->tbest);
	newsol = 1;
    
	
	// save the solution in a file (does not work if the callbacks changed it...)
	if ( newsol && (VERBOSE >= 10) )
	{
		if ( VERBOSE >= 100 ) CPXwritemipstarts(env, lp, "model.mst", 0, 0);
		printf("... New incumbent of value %20.5lf found after %7.2lf sec.s \n", inst->zbest,  inst->tbest);
		fflush(NULL);
	}

	return newsol;
}

/**************************************************************************************************************************/
int CableOpt(instance *inst)
/**************************************************************************************************************************/
{  
/* 1. initialization ------------------------------------------------- */

	//time_t tt_1 = time(NULL); 		// seconds since the epoch    
	inst->tstart = second();   
	inst->best_lb = -CPX_INFBOUND;   
	int done = 0;
	int solved = 2;
	// open cplex model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "WFCP"); 
	if(inst->rins != -1)
		if(CPXsetintparam(env,CPXPARAM_MIP_Strategy_RINSHeur,inst->rins))
			print_error("Error set of rins");
	if(inst->polishing_time > 0.0 )
		if(CPXsetdblparam(env,CPXPARAM_MIP_PolishAfter_Time,inst->polishing_time))
			print_error("Error set of polishing time");
	if ( inst->randomseed != 0 ) 
		CPXsetintparam(env, CPX_PARAM_RANDOMSEED, fabs(inst->randomseed));
	if( inst->gap >= 0.0 && inst->gap < 1.0)
		if(CPXsetdblparam(env,CPX_PARAM_EPGAP,inst->gap))
			print_error("Error set of minimum gap");
/* 2. build initial model  ------------------------------------------------- */
	gp = popen("gnuplot -p","w");
	build_model(inst, env, lp);
	
	
	int ncols = CPXgetnumcols(env, lp);
	inst->best_sol = (double *) calloc(ncols, sizeof(double)); 	// all entries to zero  
	inst->zbest = CPX_INFBOUND;  
/* 3. final MIP run ------------------------------------------------------ */


	mip_timelimit(env, CPX_INFBOUND, inst);
	if ( VERBOSE >= 50 ) 
		CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); 
	if ( time_limit_expired(inst) ) 
	{
		CPXgetbestobjval(env, lp, &inst->best_lb); 
		mip_update_incumbent(env, lp, inst);                 

		// free pools and close cplex model
		CPXfreeprob(env, &lp);
		CPXcloseCPLEX(&env); 	
		
		return 0;
	}         
	if ( VERBOSE >= 100 ) 
		CPXwriteprob(env, lp, "model/final.lp", NULL);  
	     
	//installplotfunction(env, lp, inst);
	while(!done && solved)
	{
		CPXmipopt(env,lp);     
		CPXgetbestobjval(env, lp, &inst->best_lb); 
		//plotGraph(env, lp, inst);
		mip_update_incumbent(env, lp, inst);
		done = 1;
		if(nocross_separation(env, lp, inst))
		{
			/*
			inst->timelimit = inst->timelimit - (second() - inst->tstart);
			mip_timelimit(env, CPX_INFBOUND, inst);
			inst->tstart = second();
			*/
			done = 0;
		}
		else
		{
			/*
			inst->timelimit = inst->timelimit - (second() - inst->tstart);
			mip_timelimit(env, CPX_INFBOUND, inst);
			inst->tstart = second();
			*/
			printf("NO CROSS CONSTRAINTS TO ADD\n");
			mip_update_incumbent(env, lp, inst);
			if(CPXsetdblparam(env,CPX_PARAM_EPGAP,0.0))
				print_error("Error set of minimum gap");
			 done = 0;
			 solved --;
		}
	}	
	
	printf("\n     Time,Solution,Lower B,GAP%%,Seed,Rins,Relax,Polishing Time\n");
	printf("STAT,%.0lf,%.0lf,%.0lf,%.2lf%%,%d,%d,%d,%.0lf\n\n",inst->tbest,inst->zbest,inst->best_lb,((inst->zbest-inst->best_lb)*100/inst->zbest),inst->randomseed,inst->rins,inst->relax,inst->polishing_time);
	
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env); 	
	fclose(gp);
	return 0;
}  
 
/**************************************************************************************************************************/
void build_model0(instance *inst, CPXENVptr env, CPXLPptr lp) 
/**************************************************************************************************************************/
{
  
	double zero = 0.0; // one = 1.0; 	
	char binary = 'B'; 
	char continuous = 'C';
	//char integer = 'I';

	char **cname = (char **) calloc(1, sizeof(char *));		// (char **) required by cplex...
	cname[0] = (char *) calloc(100, sizeof(char));

	//////////////////////////////////////////////////////////////////////////////////////////////// add binary var.s y(i,j) 
	if ( inst->ystart != -1 ) 
		print_error(" ... error in build_model(): var. y cannot be redefined!");
	inst->ystart = CPXgetnumcols(env,lp); 		// position of the first q(,) variable   
	for ( int i = 0; i < inst->nturbines; i++ )
	{
		for ( int j = 0; j < inst->nturbines; j++ )
		{
			sprintf(cname[0], "y(%d,%d)", i+1,j+1);
			double obj = 0.0;
			double ub = ( i == j ) ?  0.0 : 1.0;
			if ( CPXnewcols(env, lp, 1, &obj, &zero, &ub, &binary, cname) ) 
				print_error(" wrong CPXnewcols on y var.s");
			//if ( CPXnewcols(env, lp, 1, &obj, &zero, &ub, &integer, cname) ) print_error(" wrong CPXnewcols on q var.s");
    		if ( CPXgetnumcols(env,lp)-1 != ypos(i,j, inst) ) 
    			print_error(" wrong position for y var.s");
		}
	} 

	//////////////////////////////////////////////////////////////////////////////////////////// add binary var.s x(i,j)   
	if ( inst->xstart != -1 ) 
		print_error(" ... error in build_model(): var. x cannot be redefined!");
	
	inst->xstart = CPXgetnumcols(env,lp); 		// position of the first x(,) variable   
	for ( int i = 0; i < inst->nturbines; i++ )
	{
		for ( int j = 0; j < inst->nturbines; j++ )
		{
			for(int k = 0 ; k < inst->ncables; k++)
			{
				sprintf(cname[0], "x(%d,%d,%d)", i+1,j+1,k+1);
				double obj = dist(i,j,inst)*(inst->cablecost[k]);   
				double ub =( i == j) ? 0.0 : 1.0;
				if ( CPXnewcols(env, lp, 1, &obj, &zero, &ub, &binary, cname) ) 
					print_error(" wrong CPXnewcols on x var.s");
	    		if ( CPXgetnumcols(env,lp)-1 != xpos(i,j,k, inst) ) 
	    		{
	    			printf("I: %d J:%d K:%d ---- xpos: %d getCols: %d\n", i, j, k, xpos(i,j,k,inst), CPXgetnumcols(env,lp)-1 );
	    			print_error(" wrong position for x var.s");
	    		}

	    	}
		}
	} 

	//////////////////////////////////////////////////////////////////////////////////////////////// add continuous var.s f(i,j) = flux throws (i,j) edge
	if ( inst->fstart != -1 ) 
		print_error(" ... error in build_model(): var. f cannot be redefined!");
	inst->fstart = CPXgetnumcols(env,lp); 		// position of the first f() variable   
	for ( int i = 0; i < inst->nturbines; i++ )
	{
		for ( int j = 0; j < inst->nturbines; j++ )
		{
			sprintf(cname[0], "f(%d,%d)", i+1,j+1);
			double obj = 0.0;
			double ub =  CPX_INFBOUND;
			if ( CPXnewcols(env, lp, 1, &obj, &zero, &ub, &continuous, cname) ) 
				print_error(" wrong CPXnewcols on f var.s");
	    	if ( CPXgetnumcols(env,lp)-1 != fpos(i,j, inst) ) 
	    		print_error(" wrong position for f var.s");
	    }
	} 
	
	///////////////////////////////////////////////////////////////////////////////////////////////// add s var
	if(inst->relax == 1)
	{
		sprintf(cname[0], "s");
		double obj = 1e9;
		double ub = CPX_INFBOUND;
		if ( CPXnewcols(env, lp, 1, &obj, &zero, &ub, &continuous, cname) ) 
				print_error(" wrong CPXnewcols on s var.s");
		inst->sstart = CPXgetnumcols(env,lp)-1;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////// 
	for ( int h = 0; h < inst->nturbines; h++ )  // out edges constraints
	{
		int lastrow = CPXgetnumrows(env,lp);
		double rhs = ( inst->power[h] >= 0 ) ? 1.0 : 0.0; 
		char sense = 'E'; 
		sprintf(cname[0], "OutEdges(%d)", h+1);   
		if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) 
			print_error(" wrong CPXnewrows [x1]");
		for ( int j = 0; j < inst->nturbines; j++ )
		{
			if ( CPXchgcoef(env, lp, lastrow, ypos(h,j, inst), 1.0) ) 
				print_error(" wrong CPXchgcoef [x1]");
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	for ( int h = 0; h < inst->nturbines; h++ )  // in root edges constraints
	{
		if(inst->power[h] <= 0.0)
		{
			int lastrow = CPXgetnumrows(env,lp);
			double rhs = inst->C; 
			char sense = 'L'; 
			sprintf(cname[0], "inRoot[%d]",h+1);   
			if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) 
				print_error(" wrong CPXnewrows [x1]");
			for ( int i = 0; i < inst->nturbines; i++ )
			{
				if ( CPXchgcoef(env, lp, lastrow, ypos(i,h, inst), 1.0) ) 
					print_error(" wrong CPXchgcoef [x1]");
			}
			if(inst->relax == 1)
			{
				if ( CPXchgcoef(env, lp, lastrow, inst->sstart, -1.0) ) 
					print_error(" wrong CPXchgcoef [x1]");
			}
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	for ( int h = 0; h < inst->nturbines; h++ )  // flux constraints
	{
		if(inst->power[h] >= 0.0)
		{
			int lastrow = CPXgetnumrows(env,lp);
			double rhs = inst->power[h]; 
			char sense = 'E'; 
			sprintf(cname[0], "fluxThrows(%d)", h+1);   
			if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) 
				print_error(" wrong CPXnewrows [x1]");
			for ( int j = 0; j < inst->nturbines; j++ )
			{
				if ( CPXchgcoef(env, lp, lastrow, fpos(h,j, inst), 1.0) ) 
					print_error(" wrong CPXchgcoef [x1]");
				if ( CPXchgcoef(env, lp, lastrow, fpos(j,h, inst), -1.0) ) 
					print_error(" wrong CPXchgcoef [x1]");
			}
			
		}
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	for ( int i = 0; i < inst->nturbines; i++ )  // choice of cable constraints
	{
		for ( int j = 0; j < inst->nturbines; j++ )  
		{
			int lastrow = CPXgetnumrows(env,lp);
			double rhs = 0.0; 
			char sense = 'E'; 
			sprintf(cname[0], "EdgesCable(%d, %d)", i+1, j+1);   
			if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) 
				print_error(" wrong CPXnewrows [x1]");
			for ( int k = 0; k < inst->ncables; k++ )
			{
				if ( CPXchgcoef(env, lp, lastrow, xpos(i, j, k, inst), 1.0) ) 
					print_error(" wrong CPXchgcoef [x1]");
			}
			if ( CPXchgcoef(env, lp, lastrow, ypos(i, j, inst), -1.0) ) 
				print_error(" wrong CPXchgcoef [x1]");
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	for ( int i = 0; i < inst->nturbines; i++ )  // total capacity of edges constraints
	{
		for ( int j = 0; j < inst->nturbines; j++ )  
		{
			int lastrow = CPXgetnumrows(env,lp);
			double rhs = 0.0; 
			char sense = 'G'; 
			sprintf(cname[0], "EdgesCapacity(%d, %d)", i+1, j+1);   
			if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) 
				print_error(" wrong CPXnewrows [x1]");
			for ( int k = 0; k < inst->ncables; k++ )
			{
				if ( CPXchgcoef(env, lp, lastrow, xpos(i, j, k, inst), inst->cablecapacity[k]) ) 
					print_error(" wrong CPXchgcoef [x1]");
			}
			if ( CPXchgcoef(env, lp, lastrow, fpos(i, j, inst), -1.0) ) 
				print_error(" wrong CPXchgcoef [x1]");
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	if(inst->noCross == 0)
	{
		for ( int i = 0; i < inst->nturbines; i++ )  // no Cross condition
		{
			for ( int j = i+1; j < inst->nturbines; j++ )  
			{
				for ( int k = 0; k < inst->nturbines; k++ )
				{
					int lastrow = CPXgetnumrows(env,lp);
					double rhs = 1.0; 
					char sense = 'L'; 
					sprintf(cname[0], "CrossON(%d, %d) by OutEdgeFrom(%d)", i+1, j+1,k+1);   
					if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) 
						print_error(" wrong CPXaddlazyconstraints [x1]");
					if ( CPXchgcoef(env, lp, lastrow, ypos(i, j, inst), 1.0) ) 
						print_error(" wrong CPXchgcoef [x1]");
					if ( CPXchgcoef(env, lp, lastrow, ypos(j, i, inst), 1.0) ) 
						print_error(" wrong CPXchgcoef [x1]");
					for ( int l = k+1; l < inst->nturbines; l++ )
					{
						if(!noCross(i,j,k,l, inst))
						{
							if ( CPXchgcoef(env, lp, lastrow, ypos(k, l, inst), 1.0) ) 
								print_error(" wrong CPXchgcoef [x1]");

						}
					}
				}
			}
		}
	}
	else if(inst->noCross == 1)
	{
		for ( int i = 0; i < inst->nturbines; i++ )  // no Cross condition add to pool of lazy constraints
		{
			for ( int j = i+1; j < inst->nturbines; j++ )  
			{
				for ( int k = 0; k < inst->nturbines; k++ )
				{
					if(k == i || k == j)
						continue;

					int vectInd[inst->nturbines];			// specifies the column index of the corresponding coefficient
					double vectVal[inst->nturbines];		// corrisponding coefficient
					int vectBag[inst->nturbines];			// The nonzero elements of every lazy constraint must be stored in sequential locations in this array from position rmatbeg[i] to rmatbeg[i+1]-1
					int nzcnt = 0;							// An integer that specifies the number of nonzero constraint coefficients to be added to the constraint matrix. This specifies the length of the arrays rmatind and rmatval.
									
					vectBag[0] = 0;

					vectInd[0] = ypos(i,j,inst);
					vectVal[0] = 1.0;
					nzcnt ++;
					
					vectInd[1] = ypos(j,i,inst);
					vectVal[1] = 1.0;
					nzcnt ++;

					for ( int l = k+1; l < inst->nturbines; l++ )
					{
						if(!noCross(i,j,k,l, inst))
						{
							vectInd[nzcnt] = ypos(k,l,inst);
							vectVal[nzcnt] = 1.0;
							nzcnt ++;
						}
					}
					
					vectBag[1] = nzcnt;
					if(nzcnt > 2)
					{
						double rhs = 1.0; 
						char sense = 'L'; 
						sprintf(cname[0], "CrossON(%d, %d) by OutEdgeFrom(%d)", i+1, j+1,k+1);  
						if ( CPXaddlazyconstraints (env, lp, 1, nzcnt, &rhs, &sense, vectBag, vectInd, vectVal, cname) ) 
							print_error(" wrong CPXaddlazyconstraints [x1]");
					}
				}
			}
		}
	} 
	free(cname[0]);
	free(cname);
}

/***************************************************************************************************************************/
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp)
/**************************************************************************************************************************/
{    
	

	switch (inst->model_type) 	
	{
		case 0 : 								
			build_model0(inst, env,lp);
			break;
		default: 
			print_error(" model type unknown!!"); 
			break;
		
	}  

	
	if ( VERBOSE >= -100 ) CPXwriteprob(env, lp, "model/model.lp", NULL); 
	
}  

int CPXPUBLIC plotcallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, double *objval_p,  double *x, int *checkfeas_p,  int *useraction_p)
{

	instance *inst = (instance *)cbhandle;
	inst->zbest = objval_p[0];
	inst->tbest = second() - inst->tstart;
	FILE *f;
	char filename[30];


	for(int k = 0; k < inst->ncables; k++)
	{
		plt[k] = 0;
		sprintf(filename,"plot/plot%d.dat",k);
		f = fopen(filename,"w");
		for(int i = 0; i < inst->nturbines; i++)
		{
			for(int j = 0; j < inst->nturbines; j++)
			{
				if(x[xpos(i,j,k,inst)] > 0.5)
				{
					fprintf(f, "%lf %lf %d\n", inst->xcoord[i], inst->ycoord[i], i );
					fprintf(f, "%lf %lf %d\n\n", inst->xcoord[j], inst->ycoord[j], j );
					plt[k]++;
				}
			}
		}
		fclose(f);
	}
	makeScript(inst);
	fprintf(gp, "load 'plot/script_plot.p'\n"); 
	return 0;
}

int installplotfunction( CPXENVptr env, CPXLPptr lp, instance *inst)
{
	if(CPXsetheuristiccallbackfunc(env, plotcallback, (void *)inst))
		print_error("Error lazy callback function");
	return 0;
}

int compute_nocross_cut(instance *inst, double *x,int i, int j, int k, int *index, double *value)
{
	int count = 0;
	index[count] =  ypos(j,i,inst);
	value[count] = 1.0;
	count++;
	index[count] = ypos(i,j,inst);
	value[count] = 1.0;
	count++;
	for ( int l = k+1; l < inst->nturbines; l++ )
	{
		if(l == i || l == j)continue;
		//if(x[ypos(k,l,inst)] < 0.5)continue;
		if(!noCross(i,j,k,l, inst))
		{
			index[count] = ypos(k,l,inst);
			value[count] = 1.0;
			count++;
		}
	}
	if(count < 3)
		count = 0;
	return count;
}
int nocross_separation(CPXENVptr env, CPXLPptr lp, instance *inst)
{
	int count = 0;
	double *x;
	x = (double *) calloc(CPXgetnumcols(env, lp), sizeof(double)); 
	CPXgetx(env, lp, x, 0, CPXgetnumcols(env, lp) -1);
	for ( int i = 0; i < inst->nturbines; i++ )  // no Cross condition
	{
		for ( int j = i+1; j < inst->nturbines; j++ )  
		{

			if(x[ypos(i,j,inst)] < 0.5 && x[ypos(j,i,inst)] < 0.5)
				continue;
			for ( int k = i+1; k < inst->nturbines; k++ )
			{
				if(i == k || j == k)continue;
				int *index;
				double *values;
				index = (int *) calloc(inst->nturbines - k + 3, sizeof(int ));
				values = (double *) calloc(inst->nturbines- k + 3, sizeof(double ));
				int nnz = compute_nocross_cut(inst, x, i, j, k, index, values);
				if(nnz == 0)
					break;
				double rhs = 1.0; 
				char sense = 'L';
				int vectBag[1]; 
				vectBag[0] = 0;
				char **cname = (char **) calloc(1, sizeof(char *));		// (char **) required by cplex...
				cname[0] = (char *) calloc(100, sizeof(char));
				sprintf(cname[0], "CrossON(%d, %d) by OutEdgeFrom(%d)", i+1, j+1,k+1); 
				if ( CPXaddrows(env, lp,0, 1, nnz, &rhs, &sense, vectBag, index, values, NULL,cname) ) 
					print_error(" wrong addRows [x1]");
				
				free(index);
				free(values);
				count++;
			}
		}
	}
	return count;
}
/*
	callback per mettere il separatore dentro il branch and bound


*/