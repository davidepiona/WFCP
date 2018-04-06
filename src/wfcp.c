#include "wfcp.h"
#include <time.h>
   
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp);
void print_error(const char *err);
double second();       
void debug(const char *err);       
int time_limit_expired(instance *inst);   

void mip_timelimit(CPXENVptr env, double timelimit, instance *inst);
int mip_update_incumbent(CPXENVptr env, CPXLPptr lp, instance *inst);

FILE *gp;
char color[20][20];
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
int makeScript(CPXENVptr env, CPXLPptr lp, instance *inst)
{
	FILE *s;	
	setColor();

	s = fopen("plot/script_plot.p","w");
	fprintf(s,"set autoscale\n");
	for(int k = 0; k < inst->ncables; k++)
	{
		if(k != 0)
			fprintf(s, "\nre" );	
		
		fprintf(s, "plot \\\n\t\'plot/plot%d.dat\' using 1:2 with lines lc rgb \"%s\" lw 2 title \"Cable %d\",\\\n",k,color[k],k+1);
		fprintf(s, "\t\'plot/plot%d.dat\' using 1:2:(0.6) with circles fill solid lc rgb \"black\" notitle,\\\n",k);
		fprintf(s, "\t\'plot/plot%d.dat\' using 1:2:3     with labels tc rgb \"black\" offset (0,0) font \'Arial Bold\' notitle",k );
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

	makeScript(env, lp, inst);

	for(int k = 0; k < inst->ncables; k++)
	{
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

				}
			}
		}
		fclose(f);
	}
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

	double a = x2 - x1;
	double b = x3 - x4;
	double c = x3 - x1;
	double d = y2 - y1;
	double e = y3 - y4;
	double f = y3 - y1;

	if( abs(a) < EPSILON || abs(e) < EPSILON)
		return 1;

	lambda = (c - (f * b / e)) / (a - (d * b / e));
	mu = (f - (c * d / a)) / (e - (b * d / a));

	if(((lambda > XSMALL) && (lambda < 1 - XSMALL)) || ((mu > XSMALL) && (mu < 1 - XSMALL)))
		return 0;

	return 1;
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

	if ( mip_value(env,lp) < inst->zbest - EPSILON )
	{
		inst->tbest = second() - inst->tstart;
		inst->zbest = mip_value(env, lp);
		CPXgetx(env, lp, inst->best_sol, 0, ncols-1);
		plotGraph(env, lp, inst);
		if ( VERBOSE >= 40 ) printf("\n >>>>>>>>>> incumbent update of value %lf at time %7.2lf <<<<<<<<\n", inst->zbest, inst->tbest);
		newsol = 1;
	}     
	
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
	     
	CPXmipopt(env,lp);     


/* 99. final statistics ------------------------------------------------- */
	
	
	CPXgetbestobjval(env, lp, &inst->best_lb); 
	mip_update_incumbent(env, lp, inst);                 

	printf("\n     Time,Solution,Lower B,GAP%%,Seed,Rins,Relax,Polishing Time\n");
	printf("STAT,%.0lf,%.0lf,%.0lf,%.2lf%%,%d,%d,%d,%.0lf\n\n",inst->tbest,inst->zbest,inst->best_lb,((inst->zbest-inst->best_lb)*100/inst->zbest),inst->randomseed,inst->rins,inst->relax,inst->polishing_time);
	
	// free pools and close cplex model
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
	free(cname[0]);
	free(cname);
	   
}

/***************************************************************************************************************************/
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp)
/**************************************************************************************************************************/
{    
	

	switch (inst->model_type) 	
	{
		case 0 : 								// basic model with asymmetric x and q
			build_model0(inst, env,lp);
			break;
		default: 
			print_error(" model type unknown!!"); 
			break;
		
	}  

	
	if ( VERBOSE >= -100 ) CPXwriteprob(env, lp, "model/model.lp", NULL); 
	
}  
