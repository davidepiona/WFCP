#include "wfcp.h"
#include <time.h>


void build_model(instance *inst, CPXENVptr env, CPXLPptr lp);
void print_error(const char *err);
double second();       
void debug(const char *err);       
int time_limit_expired(instance *inst, double time);   

void mip_timelimit(CPXENVptr env, double timelimit, instance *inst);
int mip_update_incumbent(CPXENVptr env, CPXLPptr lp, instance *inst);
int compute_nocross_cut(instance *inst, double *x, int i, int j, int k, int *index, double *value);
int nocross_separation(CPXENVptr env, CPXLPptr lp, instance *inst);
int installLazyCallback( CPXENVptr env, CPXLPptr lp, instance *inst);
int compute_nocross_cut_all(instance *inst, double *x,int i, int j, int k, int *index, double *value);
int resetHardFix(CPXENVptr env, CPXLPptr lp, instance *inst, int *index, char *lu, double *bd, int count );
int setHardFix(CPXENVptr env, CPXLPptr lp, instance *inst, int *index,char *lu, double *bd);

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
void fileNames(instance* inst, char* p)
{
	if(inst->names)	
	{
		char* path = inst->cables_file;
		memcpy(p, path+strlen(path)-11, 7); 
	}
	else
	{	
		memcpy(p, "datastd", 7);
	}
}
int makeScript(instance *inst)
{
	FILE *s;	
	setColor();
	int flag = 0;
	s = fopen("plot/script_plot.p","w");
	fprintf(s,"set autoscale\n");
	char p[8];
	fileNames(inst, p);
	//printf("\n\n----%s----\n\n", p);
	for(int k = 0; k < inst->ncables; k++)
	{
		//printf("Script cavo < %d >\n",k);
		if(plt[k] > 0)
		{
			//fprintf(s,"\nset term png\n");
			//fprintf(s,"set output '%s.png'\n", p);
			//fprintf(s,"set datafile separator \" \"\n");
			if(k != 0 && flag != 0)
				fprintf(s, "\nre" );	
			fprintf(s, "plot \\\n\t\'plot/plot_%s_%d.dat\' using 1:2 with lines lc rgb \"%s\" lw 2 title \"Cable %d\",\\\n", p,k,color[k],k+1);
			fprintf(s, "\t\'plot/plot_%s_%d.dat\' using 1:2:(0.6) with circles fill solid lc rgb \"black\" notitle,\\\n",p ,k);
			fprintf(s, "\t\'plot/plot_%s_%d.dat\' using 1:2:3     with labels tc rgb \"black\" offset (0,0) font \'Arial Bold\' notitle",p, k );
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
	char p[8];
	fileNames(inst, p);
	for(int k = 0; k < inst->ncables; k++)
	{
		//printf("File data < %d >\n",k);
		plt[k] = 0;
		sprintf(filename,"plot/plot_%s_%d.dat",p ,k);
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
	fflush(NULL);
	return 0;
}
int plotGraph(instance *inst, double *x)
{
	FILE *f;
	char filename[30];
	char p[8];
	fileNames(inst, p);

	for(int k = 0; k < inst->ncables; k++)
	{
		//printf("File data < %d >\n",k);
		plt[k] = 0;
		sprintf(filename,"plot/plot_%s_%d.dat",p,k);
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
	fflush(NULL);
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
                         
int num_cores(instance *inst)
{                          
	int ncores = inst->num_threads;
	if ( ncores == 0 ) 
	{
		int error;
		CPXENVptr env = CPXopenCPLEX(&error);
		CPXgetnumcores(env, &ncores);  
		CPXcloseCPLEX(&env);  
	}
	return ncores;
}

int time_limit_expired(instance *inst, double time)	 
{
	double tspan = second() - inst->tstart;
	if (  tspan > time ) 
	{
		if ( VERBOSE >= 100 ) 
			printf("\n\n$$$ time limit of %10.1lf sec.s expired after %10.1lf sec.s $$$\n\n", time, tspan);
		//exit(0); 
		return 1;
	}  
	return 0;
}

void mip_timelimit(CPXENVptr env, double timelimit, instance *inst, double time)
{
	double residual_time = inst->tstart + time - second();
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
	int newsol = 0;
	plotGraph(env,lp,inst);
	if ( mip_value(env,lp) < inst->zbest - EPSILON )
	{
		inst->tbest = second() - inst->tstart;
		inst->zbest = mip_value(env, lp);
		CPXgetx(env, lp, inst->best_sol, 0, inst->ncols-1);
		if ( VERBOSE >= 40 ) 
			printf("\n >>>>>>>>>> incumbent update of value %lf at time %7.2lf <<<<<<<<\n", inst->zbest, inst->tbest);
		newsol = 1;
	}     
	
	// save the solution in a file (does not work if the callbacks changed it...)
	if ( newsol && (VERBOSE >= 10) )
	{
		if ( VERBOSE >= 100 ) CPXwritemipstarts(env, lp, "model.mst", 0, 0);
		printf("... New incumbent of value %20.5lf found after %7.2lf sec.s \n", inst->zbest, inst->tbest);
		fflush(NULL);
	}

	return newsol;
}
int mip_update_incumbent(instance *inst, double *x, double z)
{
	int newsol = 0;
	plotGraph(inst, x);
	if ( z < inst->zbest - EPSILON )
	{
		inst->tbest = second() - inst->tstart;
		inst->zbest = z;
		inst->best_sol = x;
		if ( VERBOSE >= 40 ) 
			printf("\n >>>>>>>>>> incumbent update of value %lf at time %7.2lf <<<<<<<<\n", inst->zbest, inst->tbest);
		newsol = 1;
	}     
	
	// save the solution in a file (does not work if the callbacks changed it...)
	if ( newsol && (VERBOSE >= 10) )
	{
		if ( VERBOSE >= 100 )
		printf("... New incumbent of value %20.5lf found after %7.2lf sec.s \n", inst->zbest, inst->tbest);
		fflush(NULL);
	}

	return newsol;
}
/**********************************************************************************************************/
double mip_getobj(CPXENVptr env, CPXLPptr lp, const double *x)
/**********************************************************************************************************/
{
	int ncols = CPXgetnumcols(env, lp);
	double *obj = (double *) calloc(ncols, sizeof(double));
	CPXgetobj(env, lp, obj, 0, ncols-1);
	double value = 0.0;
	for ( int j = 0; j < ncols; j++ ) value += obj[j] * x[j];
	free(obj);
	return value;
}

void mip_reload_solution(CPXENVptr env, CPXLPptr lp, int ncols, double *xstar)
{
	int nmipstart = CPXgetnummipstarts(env, lp);
	if ( nmipstart > 0 && CPXdelmipstarts (env, lp, 0, nmipstart-1) ) 
		print_error("CPXdelmipstarts error");	

	char *ctype = (char *) calloc(ncols, sizeof(char));
	CPXgetctype(env, lp, ctype, 0, ncols-1);
	int j;
	for ( j = 0; j < ncols; j++ )
	{
		if ( ctype[j] == 'B' ) xstar[j] =  (xstar[j] < 0.5) ? 0.0 : 1.0;
	}
	free(ctype);

	int *indices = (int *) malloc(ncols*sizeof(int));
	for ( j = 0; j < ncols; j++ ) indices[j] = j;

	int beg = 0;
	int effortlevel = CPX_MIPSTART_SOLVEFIXED; 										// it was CPX_MIPSTART_SOLVEMIP;
	CPXaddmipstarts(env, lp, 1, ncols, &beg, indices, xstar, &effortlevel, NULL); 	// it was CPXcopymipstart(env, lp, ncols, indices, xstar);
	free(indices);

	double myval = mip_getobj(env, lp, xstar);
	if ( VERBOSE >= 10 ) 
		printf("Solution of value %lf (%lf) reloaded!\n\n", mip_value(env,lp),myval);
}
/**************************************************************************************************************************/
int CableOpt(instance *inst)
/**************************************************************************************************************************/
{  
/* 1. initialization ------------------------------------------------- */
  
	inst->tstart = second();   
	inst->best_lb = -CPX_INFBOUND; 
	inst->num_threads = num_cores(inst);  
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
		if(inst->noCross != 4)
		{
			if(CPXsetdblparam(env,CPXPARAM_MIP_PolishAfter_Time,inst->polishing_time))
				print_error("Error set of polishing time");
			else
				if(CPXsetdblparam(env,CPXPARAM_MIP_PolishAfter_Time,inst->timelimit - inst->polishing_time))
					print_error("Error set of polishing time");
		}
	if ( inst->randomseed != 0 ) 
		CPXsetintparam(env, CPX_PARAM_RANDOMSEED, fabs(inst->randomseed));
	if( inst->gap >= 0.0 && inst->gap < 1.0)
		if(CPXsetdblparam(env,CPX_PARAM_EPGAP,inst->gap))
			print_error("Error set of minimum gap");
	if(inst->num_threads > 1)
	{
		if ( VERBOSE >= 2 ) printf(" ... Cplex in opportunistic mode with %d thread(s)\n", inst->num_threads);
		CPXsetintparam(env, CPX_PARAM_PARALLELMODE, -1); 	
	} 
/* 2. build initial model  ------------------------------------------------- */
	gp = popen("gnuplot -p","w");
	build_model(inst, env, lp);
		
	
	inst->ncols = CPXgetnumcols(env, lp);
	inst->best_sol = (double *) calloc(inst->ncols, sizeof(double)); 	// all entries to zero  
	inst->zbest = CPX_INFBOUND;  
/* 3. final MIP run ------------------------------------------------------ */


	if ( VERBOSE >= 50 ) 
		CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); 
	if ( time_limit_expired(inst, inst->timelimit) ) 
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
	 
	if( inst->noCross == 4) 		// HardFixing
	{
		int done = 0;
		if(inst->timeStartSol > 0)
			mip_timelimit(env, CPX_INFBOUND, inst, inst->timeStartSol);
		else if(inst->timeloop > 0)
			mip_timelimit(env, CPX_INFBOUND, inst, inst->timeloop);
		else
			printf("Set time to find start solution\n");

		CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);	
		installLazyCallback(env,lp,inst);
		CPXmipopt(env,lp); 
		CPXgetbestobjval(env, lp, &inst->best_lb);    
		printf("Hard Fixing start\n");
		int times = 1;
		while(!done)
		{
			
			if(inst->timeloop > 0 && (inst->timelimit - second() + inst->tstart) > inst->timeloop)
			{
				mip_timelimit(env, CPX_INFBOUND, inst, inst->timeloop*times);
				if(inst->polishing_time > 0.0 )
					if(CPXsetdblparam(env,CPXPARAM_MIP_PolishAfter_Time,(inst->timeloop*times - inst->polishing_time)))
						print_error("Error set of polishing time");
			}
			else
			{
				mip_timelimit(env, CPX_INFBOUND, inst, inst->timelimit - (second() - inst->tstart));
				if(inst->polishing_time > 0.0 )
				if(CPXsetdblparam(env,CPXPARAM_MIP_PolishAfter_Time,(inst->timelimit - (second() - inst->tstart) - inst->polishing_time)))
					print_error("Error set of polishing time");
			}
			
			times++;
			int *index;
			char *lu;
			double *bd;

			index = (int *) calloc(inst->nturbines * inst->nturbines * inst->nturbines * inst->nturbines , sizeof(int ));
			lu = (char *) calloc(inst->nturbines * inst->nturbines * inst->nturbines * inst->nturbines , sizeof(char ));
			bd = (double *) calloc(inst->nturbines * inst->nturbines * inst->nturbines * inst->nturbines , sizeof(double ));
			
			int nzcnt = setHardFix( env,  lp, inst, index, lu, bd);

			mip_reload_solution(env, lp, inst->ncols, inst->best_sol);

			installLazyCallback(env,lp,inst);
			CPXmipopt(env,lp); 
			
			
			resetHardFix( env,  lp, inst, index, lu, bd, nzcnt );
			free(index);
			free(lu);
			free(bd);
			
			
			if(time_limit_expired(inst, inst->timelimit) || times > 1000)
				done = 1;
		}
		plotGraph( inst, inst->best_sol);
	}   
	else if(inst->noCross == 3)		// CallBack
	{
		mip_timelimit(env, CPX_INFBOUND, inst, inst->timelimit);
		CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);	
		installLazyCallback(env,lp,inst);
		CPXmipopt(env,lp);     
		CPXgetbestobjval(env, lp, &inst->best_lb); 
		mip_update_incumbent(env, lp, inst);
	}
	else if(inst->noCross == 2)	// Loop Method
	{
		int times = 1;
		while(!done && solved)
		{
			if(inst->timeloop > 0)
				mip_timelimit(env, CPX_INFBOUND, inst, inst->timeloop*times);
			CPXmipopt(env,lp);     
			CPXgetbestobjval(env, lp, &inst->best_lb); 
			plotGraph(env, lp, inst);
			done = 1;
			if(nocross_separation(env, lp, inst))
			{
				done = 0;
			}
			else
			{
				printf("NO CROSS CONSTRAINTS TO ADD\n");
				mip_update_incumbent(env, lp, inst);
				double t = inst->timelimit - second() + inst->tstart;
				printf("\n\n time %lf \n\n",t);
				mip_timelimit(env, CPX_INFBOUND, inst, t);
				if(CPXsetdblparam(env,CPX_PARAM_EPGAP,0.0))
					print_error("Error set of minimum gap");
				 done = 0;
				 solved --;
			}
			if(time_limit_expired(inst, inst->timelimit))
				solved = 0;
		}
	}	
	else 					// normal execution if noCross = 10 add it as constraints
	{
		CPXmipopt(env,lp);     
		CPXgetbestobjval(env, lp, &inst->best_lb); 
		mip_update_incumbent(env, lp, inst);
	}
	printf("\n     Time,Solution,Lower B,GAP%%,Seed,Rins,Relax,Polishing Time\n");
	printf("STAT,%.0lf,%.0lf,%.0lf,%.2lf%%,%d,%d,%d,%.0lf\n\n",inst->tbest,inst->zbest,inst->best_lb,((inst->zbest-inst->best_lb)*100/inst->zbest),inst->randomseed,inst->rins,inst->relax,inst->polishing_time);
	fprintf(gp, "exit\n");
	fclose(gp);
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env); 	
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
	if(inst->noCross == 10)
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

int myseparation(instance *inst, double *x, CPXCENVptr env, void *cbdata, int wherefrom)
{
	
	int count = 0;
	//printf("\n\nCONSTRAINT ADD :\n");
	for ( int i = 0; i < inst->nturbines; i++ )  // no Cross condition
	{
		for ( int j = i+1; j < inst->nturbines; j++ )  
		{
			if(i == j)continue;

			if(x[ypos(i,j,inst)] < 0.1 && x[ypos(j,i,inst)] < 0.1)
				continue;
			//printf("i< %d > j< %d >\n",i,j );
			for ( int k = i+1; k < inst->nturbines; k++ )
			{
				if(i == k || j == k)continue;
				int *index;
				double *values;
				index = (int *) calloc(inst->nturbines, sizeof(int ));
				values = (double *) calloc(inst->nturbines , sizeof(double ));
				int nnz = compute_nocross_cut(inst, x, i, j, k, index, values);
				if(nnz == 0)
				{
					free(index);
					free(values);
					continue;
				}
				double rhs = 1.0; 
				char sense = 'L';
				if(CPXcutcallbackadd(env,cbdata,wherefrom,nnz,rhs,sense,index,values,0))
					print_error("USER_SEPARATION_CUT_ERROR");
				
				free(index);
				free(values);
				count++;
			}
		}
	}
	//printf("\n\n");
	return count;
}
static int CPXPUBLIC lazyCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p)
{
	*useraction_p = CPX_CALLBACK_DEFAULT;
	instance *inst = (instance*)cbhandle;
	double z = CPX_INFBOUND;
	CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &z);
	double *xstar = (double*) malloc(inst->ncols * sizeof(double));
	if ( CPXgetcallbacknodex(env, cbdata, wherefrom, xstar, 0, inst->ncols-1) )
		return 1;
	
	int ncuts = myseparation(inst, xstar, env, cbdata, wherefrom);
	printf("Cuts added : < %d >\n",ncuts );
	if(ncuts >= 1)
	{
		plotGraph(inst, xstar);
		*useraction_p = CPX_CALLBACK_SET;
	}
	else
	{
		mip_update_incumbent(inst, xstar, z);
	}
	return 0;

}
int installLazyCallback( CPXENVptr env, CPXLPptr lp, instance *inst)
{
	if(CPXsetlazyconstraintcallbackfunc(env, lazyCallback, (void *)inst))
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
	
	for ( int l = 0; l < inst->nturbines; l++ )
	{
		if(l == i || l == j || l == k)continue;
		if(x[ypos(k,l,inst)] < 0.1)continue;
		if(noCross(i,j,k,l, inst))
		{
			index[count] = ypos(k,l,inst);
			value[count] = 1.0;
			count++;
			//printf("ADD constraints %d %d %d %d\n",i,j,k,l);
		}
	}
	
	if(count < 3)
		count = 0;
	return count;
}
int compute_nocross_cut_all(instance *inst, double *x,int i, int j, int k, int *index, double *value)
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
		if(l == i || l == j || l == k)continue;
		//if(x[ypos(k,l,inst)] < 0.1)continue;
		if(noCross(i,j,k,l, inst))
		{
			index[count] = ypos(k,l,inst);
			value[count] = 1.0;
			count++;
			//printf("ADD constraints %d %d %d %d\n",i,j,k,l);
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
	//printf("\n\nCONSTRAINT ADD :\n");
	for ( int i = 0; i < inst->nturbines; i++ )  // no Cross condition
	{
		for ( int j = i+1; j < inst->nturbines; j++ )  
		{
			if(i == j)continue;

			if(x[ypos(i,j,inst)] < 0.1 && x[ypos(j,i,inst)] < 0.1)
				continue;
			//printf("i< %d > j< %d >\n",i,j );
			for ( int k = i+1; k < inst->nturbines; k++ )
			{
				if(i == k || j == k)continue;
				int *index;
				double *values;
				index = (int *) calloc(inst->nturbines + 3, sizeof(int ));
				values = (double *) calloc(inst->nturbines + 3, sizeof(double ));
				int nnz = compute_nocross_cut(inst, x, i, j, k, index, values);
				if(nnz == 0)
				{
					free(index);
					free(values);
					continue;
				}
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

int setHardFix(CPXENVptr env, CPXLPptr lp, instance *inst, int *index,char *lu, double *bd)
{
	int count = 0;
	for(int i = 0; i < inst->nturbines; i++)
		for(int j = i + 1; j < inst->nturbines; j++)
		{
			if(i == j) continue;
			
			if((inst->best_sol[ypos(i,j,inst)] > 0.1 ) && (rand() % 100) > 50)
			{
			
				index[count] = ypos(i,j,inst);
				lu[count] = 'L';
				bd[count] = 1;
			}
		}
	
	CPXchgbds (env, lp, count, index, lu, bd);		
	return count;
}
int resetHardFix(CPXENVptr env, CPXLPptr lp, instance *inst, int *index, char *lu, double *bd, int count )
{
	for(int i = 0; i < count ; i++)
		if(lu[i] == 'L')
			bd[i] = 0;
		else if(lu[i] == 'U')
			bd[i] = 1;
	CPXchgbds (env, lp, count, index, lu, bd);	
	return 0;
}
/*
	Algoritmi euristici
	Cercano una soluzione ottima senza dimostrarne l'ottimalità.
			________________________
  input	--->|	     Algoritmo  	| --->	Soluzione 
	?	--->|         esatto        | --->   ottima
		--->|_______________________| --->

	Hard fixing
	Supponiamo di avere una soluzione del nostro problema, avremmo un vettore di soluzioni (x, f, y),  ci concentriamo sulle y!
	Posso fissare alcune variabili. Scelgo un certo numero di archi (50%) e decido di fissarli, ovvero aggiungere al mio modello le variabili da fissare
	con lower bound e upper bound fissati al valore di riferimento. (mi conviene fissare le variabili a 1 ).

	Codice:
	build_model(); install callback; set CPX params;
	
	--> exact_solve()
			o
	--> hard_fix() 
	
	CPX_getx()
	grafico()
	stampa()


	hardFix()
	{
		while()
		{
			fisso un timelimit ragionevole
			prendo una soluzione di riferimento

			scelgo gli archi da fissare
			CPXchgbound()
			
			CPXmipopt() -> con warmStart(bestSol)
			aggiorno la inst->best_sol
			unfixing delle variabili ---> CPXchgbound()

		}
	}
*/