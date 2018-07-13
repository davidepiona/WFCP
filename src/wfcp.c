#include "wfcp.h"
#include <time.h>
#include "float.h" 

void execute(instance *inst, CPXENVptr env, CPXLPptr lp);
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp);
void print_error(const char *err);
double second();       
void debug(const char *err);       

/*****************************************************************************************************************/
// 													utilities.c function
/*****************************************************************************************************************/
int xpos(int i, int j, int k, instance *inst);
int ypos(int i, int j, instance *inst);
int fpos(int i, int j, instance *inst);
int spos(int s, instance *inst);
void wait(int time);
int time_limit_expired(instance *inst, double time);   
void mip_timelimit(CPXENVptr env, double timelimit, instance *inst);
int mip_update_incumbent(CPXENVptr env, CPXLPptr lp, instance *inst);
int compute_nocross_cut(instance *inst, double *x, int i, int j, int k, int *index, double *value);
void mip_timelimit(CPXENVptr env, double timelimit, instance *inst, double time);
double mip_value(CPXENVptr env, CPXLPptr lp);
int mip_update_incumbent(CPXENVptr env, CPXLPptr lp, instance *inst);
int mip_update_incumbent(instance *inst, double *x, double z);
double mip_getobj(CPXENVptr env, CPXLPptr lp, const double *x);
void mip_reload_solution(CPXENVptr env, CPXLPptr lp, int ncols, instance *inst, double *xstar);
int num_cores(instance *inst);
int time_limit_expired(instance *inst, double time);
int cableregularize(instance *inst, double *x, long double z );
int myseparation(instance *inst, double *x, CPXCENVptr env, void *cbdata, int wherefrom);
int noCross(int p1, int p2, int p3, int p4, instance *inst);
double dist(int i, int j, instance *inst);
int is_all_integer(int n, const double *x);
int compute_nocross_cut_all(instance *inst, double *x,int i, int j, int k, int *index, double *value);
int compute_nocross_cut(instance *inst, double *x,int i, int j, int k, int *index, double *value);
int simmetricHamming(double *yr, int s, int l, int *index);
int asimmetricHamming(double *yr, int s, int l, int *index);
int RinsAHamming(double *yr1, double *yr2, int s, int l, int *index);
int RinsSHamming(double *yr1, double *yr2, int s, int l, int *index);
int myseparation(instance *inst, double *x, CPXCENVptr env, void *cbdata, int wherefrom);
int cableregularize(instance *inst, double *x, long double z );
int nocross_separation(CPXENVptr env, CPXLPptr lp, instance *inst);
int PrimDijkstra(double *mat, int nodes, int *pred, int r);
int fluxCalculator(int *suc, double *flux, int nodes);
int cableregularize(instance *inst, double *x, double *flux );
int PrimDijkstraGrasp(double*mat, int nodes, int *pred, int r);
int arraycontainsINT(int *array, int length, int val);
int tabucontains(int *arrayout, int* arrayin, int n, int valout, int valin);
int antFindPath(double*mat, int nodes, int *pred, double *pheromones);
int antFindPathKruskal(double*mat, int nodes, int *pred, double *pheromones, instance* inst);
int updatepheromones(double *pheromones,double* x, double* flux,double *cost, int nodes, double *mat, double zsol);
int removeloop(double*mat, int nodes, int *suc);
/*****************************************************************************************************************/

int installLazyCallback( CPXENVptr env, CPXLPptr lp, instance *inst);
int installheuristicCallback( CPXENVptr env, CPXLPptr lp, instance *inst);
int resetHardFix(CPXENVptr env, CPXLPptr lp, instance *inst, int *index, char *lu, double *bd, int count );
int setHardFix(CPXENVptr env, CPXLPptr lp, instance *inst, int *index,char *lu, double *bd);
int setHardRins(CPXENVptr env, CPXLPptr lp, instance *inst, int *index,char *lu, double *bd);
int resetSoftFix(CPXENVptr env, CPXLPptr lp);
int setSoftFix(CPXENVptr env, CPXLPptr lp, instance *inst, int K);
int setSoftRinsA(CPXENVptr env, CPXLPptr lp, instance *inst, int K);
int setSoftFixA(CPXENVptr env, CPXLPptr lp, instance *inst, int K);
int setSoftFixS(CPXENVptr env, CPXLPptr lp, instance *inst, int K);
int PrimDijkstraMat(instance *inst);
int unoOpt(instance *inst, int* suc, int *best);
int tabuOpt(instance *inst, int* suc, int* best, int* tabulistout, int* tabulistin, int &etabu, int tabuthresh);


FILE *gp;
char color[20][20];
int plt[20];

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
	
	fprintf(s,"\nset term wxt title '%Le'\n",(long double)inst->zbest);
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
	int ncols = CPXgetnumcols(env, lp);
	double x[ncols]; 

	CPXgetx(env, lp, x, 0, ncols -1);
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
int plotGraph(instance *inst, double *x, double *flux, double cost)
{
	FILE *f;
	char filename[60];
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
				if(x[i+j*inst->nturbines] == k && flux[i+j*inst->nturbines] > 0.5)
				{
					//printf("collego il cavo ( %d ) - ( %d ) con il cavo %d\n",i,j,k );
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

double objectiveFunction(instance *inst, double*x, double*flux)
{
	long double cost = 0;
	double Mcrossing = BIG_M_CROSS;
	double Minstation = BIG_M_SUBSTATION;
	for(int i = 0; i < inst->nturbines; i++){
		for(int j = 0; j < inst->nturbines; j++)
		{
			if(x[i+j*inst->nturbines] < -0.5)continue;
			int cable = x[i+j*inst->nturbines];
			int d = dist(i,j,inst);
			if(cable == (inst->ncables-1))
			{
				d = d*flux[i+j*inst->nturbines];
			}
		
			cost = cost + inst->cablecost[cable]*d;
			
		}
	}
	//printf("\n\n\nLa soluzione trovata costa di lunghezza dei cavi: %Le\n", cost);
	int n = 0;
	for(int i = 0; i < inst->nturbines; i++){
		for(int j = 0; j < inst->nturbines; j++){
			if(x[i+j*inst->nturbines] < 0 || i == j) continue;
			for(int k = 0; k < inst->nturbines; k++){
				for(int l = 0; l < inst->nturbines; l++){
					if(x[k+l*inst->nturbines] < 0 || l == k) continue;
					//printf("provo < %d - %d > -- < %d - %d >\n",i,j,k,l);
					if(!noCross(i,j,k,l,inst))
						n++;
				}
			}
		}
	}
	cost = cost + n * Mcrossing;
	//printf("Ci sono < %d > cross quindi ora costa : %Le\n",n,cost );
	n = 0;
	for(int i = 0; i < inst->nturbines; i++){
		if( x[i] != -1 ) n++;
	}
	//printf("Cavi entranti nella substation: %d\n",n );
	if(n > inst->C ) 
	{
		cost = cost + (n - inst->C) * Minstation;	
	}
	//printf("La soluzione trovata costa : %Le\n", cost);
	return cost;
}


/**************************************************************************************************************************/
int CableOpt(instance *inst)
/**************************************************************************************************************************/
{  
/* 1. initialization ------------------------------------------------- */
  
	inst->tstart = second();   
	inst->best_lb = -CPX_INFBOUND; 
	inst->num_threads = num_cores(inst);  
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
	inst->second_best_sol = (double *) calloc(inst->ncols, sizeof(double)); 	// all entries to zero  
	inst->zbest = CPX_INFBOUND;  
	inst->second_zbest = CPX_INFBOUND;  
/* 3. final MIP run ------------------------------------------------------ */


	if ( VERBOSE >= 50 ) 
		CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); 
	if ( time_limit_expired(inst, inst->timelimit) ) 
	{
		CPXgetbestobjval(env, lp, &inst->best_lb); 
		mip_update_incumbent(env, lp, inst);                 
		CPXfreeprob(env, &lp);
		CPXcloseCPLEX(&env); 	
		return 0;
	}         

	if ( VERBOSE >= 100 ) 
		CPXwriteprob(env, lp, "model/final.lp", NULL);  
	 
	execute(inst, env, lp );

	printf("\n     Time,Solution,Lower B,GAP%%,Seed,Rins,Relax,Polishing Time, modello, cc,softF,hardF\n");
	printf("STAT,%.0lf,%.0lf,%.0lf,%.2lf%%,%d,%d,%d,%.0lf,%d,%d,%d,%d\n\n",inst->tbest,inst->zbest,inst->best_lb,((inst->zbest-inst->best_lb)*100/inst->zbest),inst->randomseed,inst->rins,inst->relax,inst->polishing_time,inst->model_type,inst->noCross,inst->softF, inst->hardF);
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
		double obj = BIG_M_CABLE;
		double ub = CPX_INFBOUND;
		if ( CPXnewcols(env, lp, 1, &obj, &zero, &ub, &continuous, cname) ) 
				print_error(" wrong CPXnewcols on s var.s");
		inst->sstart = CPXgetnumcols(env,lp)-1;
	}
	else if(inst->relax == 2 || inst->relax == 3)
	{
		if ( inst->sstart != -1 ) 
			print_error(" ... error in build_model(): var. f cannot be redefined!");
		inst->sstart = CPXgetnumcols(env,lp);
		for(int i = 0; i < inst->nturbines ; i++)
		{
			sprintf(cname[0], "s ( %d)", i);
			double obj = BIG_M_CABLE;
			double ub = CPX_INFBOUND;
			if ( CPXnewcols(env, lp, 1, &obj, &zero, &ub, &continuous, cname) ) 
				print_error(" wrong CPXnewcols on s var.s");
			if ( CPXgetnumcols(env,lp)-1 != spos(i, inst) ) 
	    		print_error(" wrong position for f var.s");
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////// 
	for ( int h = 0; h < inst->nturbines; h++ )  // out edges constraints
	{
		int lastrow = CPXgetnumrows(env,lp);
		double rhs = ( inst->power[h] >= 0 ) ? 1.0 : 0.0; 
		char sense = 'E';
		if(inst->relax == 3) 
			if(inst->power[h] >= 0)
				sense = 'L';
		
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
			char sense =  'E';
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
			if(inst->relax == 2 || inst->relax == 3)
				if ( CPXchgcoef(env, lp, lastrow, spos(h, inst), 1.0) ) 
					print_error(" wrong CPXchgcoef [x1]");	
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
/**************************************************************************************************************************/
void build_model1(instance *inst, CPXENVptr env, CPXLPptr lp) 
/**************************************************************************************************************************/
{
	inst->mat = (double *) calloc(inst->nturbines * inst->nturbines, sizeof(double ));	
	inst->cablecapacity[inst->ncables] = DBL_MAX;
	inst->cablecost[inst->ncables] = BIG_M_CABLE;
	inst->ncables++;

	printf("Creo la matrice dei costi ( distanze )\n");
	PrimDijkstraMat(inst);

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
		case 1 : 								
			build_model1(inst, env,lp);
			break;
		default: 
			print_error(" model type unknown!!"); 
			break;
		
	}  

	
	if ( VERBOSE >= -100 ) CPXwriteprob(env, lp, "model/model.lp", NULL); 	
} 
 
/***************************************************************************************************************************/
void execute0(instance *inst, CPXENVptr env, CPXLPptr lp)
/***************************************************************************************************************************/
{
	mip_timelimit(env, CPX_INFBOUND, inst, inst->timelimit);
	CPXmipopt(env,lp);     
	CPXgetbestobjval(env, lp, &inst->best_lb); 
	mip_update_incumbent(env, lp, inst);
}
/***************************************************************************************************************************/
void execute1(instance *inst, CPXENVptr env, CPXLPptr lp)
/***************************************************************************************************************************/
{
	mip_timelimit(env, CPX_INFBOUND, inst, inst->timelimit);
	CPXmipopt(env,lp);     
	CPXgetbestobjval(env, lp, &inst->best_lb); 
	mip_update_incumbent(env, lp, inst);
}
/***************************************************************************************************************************/
void execute2(instance *inst, CPXENVptr env, CPXLPptr lp)
/***************************************************************************************************************************/
{
	int done = 0;
	int solved = 2;
	int times = 1;
	while(!done && solved)
	{
		if(inst->cableReg == 1) {CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);	installheuristicCallback(env,lp,inst);}
		if(inst->timeloop > 0 && inst->timeloop*times < inst->timelimit)
			mip_timelimit(env, CPX_INFBOUND, inst, inst->timeloop*times );
		else{
			printf("time limit\n");
			mip_timelimit(env, CPX_INFBOUND, inst, inst->timelimit - inst->timeloop*(times-1));
		}
		
		CPXmipopt(env,lp);     
		CPXgetbestobjval(env, lp, &inst->best_lb); 
		mip_update_incumbent(env, lp, inst);
		plotGraph(env, lp, inst);
		
		if(!nocross_separation(env, lp, inst))
		{
			printf("NO CROSS CONSTRAINTS TO ADD\n");
		}
		if(time_limit_expired(inst, inst->timelimit))
			solved = 0;
		times ++;
	}
}
/***************************************************************************************************************************/
void execute3(instance *inst, CPXENVptr env, CPXLPptr lp)
/***************************************************************************************************************************/
{
	mip_timelimit(env, CPX_INFBOUND, inst, inst->timelimit);
	CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);		
	installLazyCallback(env,lp,inst);
	if(inst->cableReg == 1) installheuristicCallback(env,lp,inst);
	CPXmipopt(env,lp);     
	CPXgetbestobjval(env, lp, &inst->best_lb); 
	mip_update_incumbent(env, lp, inst);
}
/***************************************************************************************************************************/
void execute4(instance *inst, CPXENVptr env, CPXLPptr lp)
/***************************************************************************************************************************/
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
	if(inst->cableReg == 1) installheuristicCallback(env,lp,inst);
	CPXmipopt(env,lp); 
	CPXgetbestobjval(env, lp, &inst->best_lb);  
	mip_update_incumbent(env, lp, inst);  
	printf("Hard Fixing start\n");
	int times = 1;
	
	while(!done)
	{
		if(inst->timeloop > 0 && (inst->timelimit - (second() - inst->tstart )- inst->polishing_time) > inst->timeloop)
		{
			//printf("---------------------------------------------------------------- time loop %lf\n", inst->timeloop);
			mip_timelimit(env, CPX_INFBOUND, inst, inst->timeloop + second() - inst->tstart );
			if(inst->polishing_time > 0.0 )
				if(CPXsetdblparam(env,CPXPARAM_MIP_PolishAfter_Time,(inst->timeloop + second() - inst->tstart - inst->polishing_time)))
					print_error("Error set of polishing time");
		}
		else
		{
			mip_timelimit(env, CPX_INFBOUND, inst, inst->timelimit - (second() - inst->tstart));
			if(inst->polishing_time > 0.0 )
				if(CPXsetdblparam(env,CPXPARAM_MIP_PolishAfter_Time,(inst->timelimit - (second() - inst->tstart) - inst->polishing_time)))
					print_error("Error set of polishing time");
			done = 1;
		}
		
		int index[inst->nturbines * inst->nturbines * inst->nturbines * inst->nturbines ];
		char lu[inst->nturbines * inst->nturbines * inst->nturbines * inst->nturbines ];
		double bd[inst->nturbines * inst->nturbines * inst->nturbines * inst->nturbines ];
		for (int i = 0; i < inst->nturbines * inst->nturbines * inst->nturbines * inst->nturbines ; i++)
		{
			index[i] = 0;
			lu[i] = 0;
			bd[i] = 0;
		}

		//index = (int *) calloc(inst->nturbines * inst->nturbines * inst->nturbines * inst->nturbines , sizeof(int ));
		//lu = (char *) calloc(inst->nturbines * inst->nturbines * inst->nturbines * inst->nturbines , sizeof(char ));
		//bd = (double *) calloc(inst->nturbines * inst->nturbines * inst->nturbines * inst->nturbines , sizeof(double ));
		
		mip_reload_solution(env, lp, inst->ncols,inst, inst->best_sol);
		int nzcnt = setHardFix( env,  lp, inst, index, lu, bd);


		installLazyCallback(env,lp,inst);
		if(inst->cableReg == 1) installheuristicCallback(env,lp,inst);
		//printf("-----------------------------------START MIP OPT----------------------------------\n");
		CPXmipopt(env,lp);
		mip_update_incumbent(env, lp, inst);
		//printf("-----------------------------------END MIP OPT------------------------------------\n");
		
		
		resetHardFix( env,  lp, inst, index, lu, bd, nzcnt );
		
		
		
		if(time_limit_expired(inst, inst->timelimit) || times > 1000)
		{
			printf("----------------------------------------------------------------- ULTIMA ESECUZIONE PRIMA DI USCIRE\n");
			mip_timelimit(env, CPX_INFBOUND, inst, inst->timelimit - (second() - inst->tstart));
			if(inst->polishing_time > 0.0 )
				if(CPXsetdblparam(env,CPXPARAM_MIP_PolishAfter_Time,(inst->timelimit - (second() - inst->tstart) - inst->polishing_time)))
					print_error("Error set of polishing time");
			mip_reload_solution(env, lp, inst->ncols, inst, inst->best_sol);
			installLazyCallback(env,lp,inst);
			CPXmipopt(env,lp);
			done = 1;
		}
		times++;
	}
	plotGraph( inst, inst->best_sol);
	printf("Zbest : %lf\n",inst->zbest);
}
/***************************************************************************************************************************/
void execute5(instance *inst, CPXENVptr env, CPXLPptr lp)
/***************************************************************************************************************************/
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
	if(inst->cableReg == 1) installheuristicCallback(env,lp,inst);
	CPXmipopt(env,lp); 
	CPXgetbestobjval(env, lp, &inst->best_lb);
	mip_update_incumbent(env, lp, inst);
	//printf("Hard Fixing start\n");
	int times = 1;
	int K = 3;
	double z = inst->zbest;
	while(!done)
	{
		if(inst->timeloop > 0 && (inst->timelimit - (second() - inst->tstart )- inst->polishing_time) > inst->timeloop)
		{
			//printf("---------------------------------------------------------------- time loop %lf\n", inst->timeloop);
			mip_timelimit(env, CPX_INFBOUND, inst, inst->timeloop + second() - inst->tstart );
			if(inst->polishing_time > 0.0 )
				if(CPXsetdblparam(env,CPXPARAM_MIP_PolishAfter_Time,(inst->timeloop + second() - inst->tstart - inst->polishing_time)))
					print_error("Error set of polishing time");
		}
		else
		{
			mip_timelimit(env, CPX_INFBOUND, inst, inst->timelimit - (second() - inst->tstart));
			if(inst->polishing_time > 0.0 )
				if(CPXsetdblparam(env,CPXPARAM_MIP_PolishAfter_Time,(inst->timelimit - (second() - inst->tstart) - inst->polishing_time)))
					print_error("Error set of polishing time");
			done = 1;
		}
		
		mip_reload_solution(env, lp, inst->ncols,inst, inst->best_sol);
		
		setSoftFix( env,  lp, inst, K);


		installLazyCallback(env,lp,inst);
		if(inst->cableReg == 1) installheuristicCallback(env,lp,inst);
		//printf("-----------------------------------START MIP OPT----------------------------------\n");
		CPXmipopt(env,lp);
		mip_update_incumbent(env, lp, inst); 
		//printf("-----------------------------------END MIP OPT------------------------------------\n");
		
		resetSoftFix( env,  lp);
		
		
		if(time_limit_expired(inst, inst->timelimit) || times > 1000)
		{
			printf("----------------------------------------------------------------- ULTIMA ESECUZIONE PRIMA DI USCIRE\n");
			mip_timelimit(env, CPX_INFBOUND, inst, inst->timelimit - (second() - inst->tstart));
			if(inst->polishing_time > 0.0 )
				if(CPXsetdblparam(env,CPXPARAM_MIP_PolishAfter_Time,(inst->timelimit - (second() - inst->tstart) - inst->polishing_time)))
					print_error("Error set of polishing time");
			mip_reload_solution(env, lp, inst->ncols, inst, inst->best_sol);
			installLazyCallback(env,lp,inst);
			CPXmipopt(env,lp);
			done = 1;
		}
		if(K < 20 && z == inst->zbest)
		{
			K = K + 2;
		}
		else 
		{
			if(K > 20)
				K = 3;
			z = inst->zbest;
		}
		times++;
	}
	plotGraph( inst, inst->best_sol);
	printf("Zbest : %lf\n",inst->zbest);
}
/***************************************************************************************************************************/
void execute6(instance *inst, CPXENVptr env, CPXLPptr lp)
/***************************************************************************************************************************/
{

	int *suc;
	suc = (int *) calloc(inst->nturbines, sizeof(int));
	double *flux;
	flux = (double *) calloc((inst->nturbines*inst->nturbines), sizeof(double ));
	double *x;
	x = (double *) calloc((inst->nturbines*inst->nturbines), sizeof(double ));
	double cost = 0;


	printf("Eseguo PrimDijkstra\n");
	PrimDijkstra(inst->mat, inst->nturbines, suc, 0);

	printf("Calcolo il flusso\n");
	fluxCalculator(suc, flux, inst->nturbines);
	
	printf("Metto i cavi\n");
	cableregularize(inst, x, flux );

	printf("Calcolo quanto costa la solutione\n");
	cost = objectiveFunction(inst, x, flux);
	inst->zbest = cost;
	printf("La soluzione trovata costa : %Le\n",(long double)cost );

	inst->cablecost[inst->ncables-1] = 0;
	inst->best_lb = objectiveFunction(inst, x, flux);
	printf("La soluzione, senza considerare i cavi inventati costa : %Le\n",(long double)inst->best_lb );
	inst->cablecost[inst->ncables-1] = 10e10;


	printf("Plotto\n");
	plotGraph(inst, x, flux, cost);

	fprintf(gp, "exit\n");
	fclose(gp);
		
	gp = popen("gnuplot -p","w");

	free(x);
	free(flux);
	free(suc);
}
/***************************************************************************************************************************/
void execute7(instance *inst, CPXENVptr env, CPXLPptr lp)
/***************************************************************************************************************************/
{

	for(int i = 1; i <= inst->times; i++)
	{

		int *suc;
		suc = (int *) calloc(inst->nturbines, sizeof(int));
		double *flux;
		flux = (double *) calloc(inst->nturbines*inst->nturbines, sizeof(double ));
		double *x;
		x = (double *) calloc(inst->nturbines*inst->nturbines, sizeof(double ));

		printf("Eseguo PrimDijkstra\n");
		inst->randomseed = time(NULL)*i;
		PrimDijkstraGrasp(inst->mat, inst->nturbines, suc, inst->randomseed);

		printf("Calcolo il flusso\n");
		fluxCalculator(suc, flux, inst->nturbines);
		
		printf("Metto i cavi\n");
		cableregularize(inst, x, flux );

		printf("Calcolo quanto costa la solutione\n");
		double cost = objectiveFunction(inst, x, flux);
		inst->zbest = cost;
		printf("La soluzione trovata costa : %Le\n",(long double)cost );

		inst->cablecost[inst->ncables-1] = 0;
		inst->best_lb = objectiveFunction(inst, x, flux);
		printf("La soluzione, senza considerare i cavi inventati costa : %Le\n",(long double)inst->best_lb );
		inst->cablecost[inst->ncables-1] = 10e10;

		printf("Plotto\n");
		plotGraph(inst, x, flux, cost);
		fprintf(gp, "exit\n");
		fclose(gp);
		
		gp = popen("gnuplot -p","w");
		//wait(5);
		
		printf("------------------------------- [ %d ] --------------------------------\n", i);
	}

}
/***************************************************************************************************************************/
void execute8(instance *inst, CPXENVptr env, CPXLPptr lp)
/***************************************************************************************************************************/
{
	double cost;

	int *suc;
	suc = (int *) calloc(inst->nturbines, sizeof(int));
	double *flux;
	flux = (double *) calloc((inst->nturbines*inst->nturbines), sizeof(double ));
	double *x;
	x = (double *) calloc((inst->nturbines*inst->nturbines), sizeof(double ));


	double newcost = DBL_MAX;

	int *best;
	best = (int *) calloc((inst->nturbines+1), sizeof(int));
	double *flux_best;
	flux_best = (double *) calloc((inst->nturbines*inst->nturbines), sizeof(double ));
	double *x_best;
	x_best = (double *) calloc((inst->nturbines*inst->nturbines), sizeof(double ));
	

	printf("Eseguo PrimDijkstra\n");
	PrimDijkstra(inst->mat, inst->nturbines, suc, 0);

	printf("Calcolo il flusso\n");
	fluxCalculator(suc, flux, inst->nturbines);
	
	printf("Metto i cavi\n");
	cableregularize(inst, x, flux );

	printf("Calcolo quanto costa la solutione\n");
	cost = objectiveFunction(inst, x, flux);
	inst->zbest = cost;
	printf("La soluzione trovata costa : %Le\n",(long double)cost );

	/*inst->cablecost[inst->ncables-1] = 0;
	inst->best_lb = objectiveFunction(inst, x);
	printf("La soluzione, senza considerare i cavi inventati costa : %Le\n",(long double)inst->best_lb );
	inst->cablecost[inst->ncables-1] = 10e10;*/


	printf("Plotto\n");
	plotGraph(inst, x, flux, cost);

	fprintf(gp, "exit\n");
	fclose(gp);
		
	gp = popen("gnuplot -p","w");
	free(x);
	free(flux);


	//-------------------------------------------------------------- uno opt --------------------
	for( int count = 1; cost != newcost; count++){
		for(int i = 0; i < inst->nturbines && count != 1; i++){suc[i] = best[i];best[i] = 0;}
		for(int i = 0; i < inst->nturbines*inst->nturbines; i++){x_best[i] = 0;flux_best[i] = 0;}
		
		printf("\n\n----------------------- RUN #%d ----------------------------------------\n", count );
		cost = (count == 1) ? cost : newcost;
		
		unoOpt(inst, suc, best);

		printf("Calcolo il flusso\n");
		fluxCalculator(best, flux_best, inst->nturbines);

		printf("Metto i cavi\n");
		cableregularize(inst, x_best, flux_best );

		printf("Calcolo quanto costa la solutione\n");
		newcost = objectiveFunction(inst, x_best, flux_best);
		inst->zbest = newcost;
		printf("La soluzione trovata costa : %Le\n",(long double)newcost );

		printf("Plotto\n");
		//plotGraph(inst, x_best, flux_best, cost);

		//fprintf(gp, "exit\n");
		//fclose(gp);
		//gp = popen("gnuplot -p","w");
		//wait(1);
	}

	plotGraph(inst, x_best, flux_best, cost);
	printf("Ottimo Locale trovato, costo : %Le\n",(long double )newcost );
	//free(suc);
	//free(best);
	//free(x_best);
	//free(flux_best);
}
/***************************************************************************************************************************/
void execute9(instance *inst, CPXENVptr env, CPXLPptr lp)
/***************************************************************************************************************************/
{
	double cost;

	int *suc;
	suc = (int *) calloc(inst->nturbines, sizeof(int));
	double *flux;
	flux = (double *) calloc((inst->nturbines*inst->nturbines), sizeof(double));
	double *x;
	x = (double *) calloc((inst->nturbines*inst->nturbines), sizeof(double));


	double newcost = DBL_MAX;
	double best_cost = DBL_MAX;

	int *best;
	best = (int *) calloc((inst->nturbines+1), sizeof(int));

	int *temp;
	temp = (int *) calloc((inst->nturbines+1), sizeof(int));

	double *flux_best;
	flux_best = (double *) calloc((inst->nturbines*inst->nturbines), sizeof(double));
	double *x_best;
	x_best = (double *) calloc((inst->nturbines*inst->nturbines), sizeof(double));

	int *tabulistin;
	tabulistin = (int *) calloc((inst->nturbines*inst->nturbines), sizeof(int ));
	int *tabulistout;
	tabulistout = (int *) calloc((inst->nturbines*inst->nturbines), sizeof(int ));
	int etabu = 0;
	int thresh = 10;
	//int changethresh = 40;

	printf("Eseguo PrimDijkstra\n");
	PrimDijkstra(inst->mat, inst->nturbines, suc, 0);

	printf("Calcolo il flusso\n");
	fluxCalculator(suc, flux, inst->nturbines);
	
	printf("Metto i cavi\n");
	cableregularize(inst, x, flux );

	printf("Calcolo quanto costa la solutione\n");
	cost = objectiveFunction(inst, x,flux);
	inst->zbest = cost;
	printf("La soluzione trovata costa : %Le\n",(long double)cost );

	/*inst->cablecost[inst->ncables-1] = 0;
	inst->best_lb = objectiveFunction(inst, x);
	printf("La soluzione, senza considerare i cavi inventati costa : %Le\n",(long double)inst->best_lb );
	inst->cablecost[inst->ncables-1] = 10e10;*/


	/*printf("Plotto\n");
	plotGraph(inst, x, flux, cost);

	fprintf(gp, "exit\n");
	fclose(gp);
		
	gp = popen("gnuplot -p","w");*/
	free(x);
	free(flux);

	best_cost = cost;
	for(int h = 0; h < inst->nturbines; h++){best[h]=suc[h];}

	//-------------------------------------------------------------- tabu opt --------------------
	for( int count = 1;(inst->tstart + inst->timelimit > second()); count++){
		for(int i = 0; i < inst->nturbines*inst->nturbines; i++){x_best[i] = 0;flux_best[i] = 0;}
		if(count%70 && thresh >= 70)
		{
			thresh = 10;
		}
		printf("\n\n----------------------- RUN #%d ----------------------------------------\n", count );
		if(cost == newcost)
		{
			if(count%20)thresh = thresh + 10;
			printf("FERMO\n");
		}
		else if(cost < newcost)
		{
			printf("RISALGO\n");

		}
		else if(newcost < cost)
		{
			printf("SCENDO\n");
		}
		else {
			printf("-----------------------------------------------------------------------------------------------------------------------cost : %f  new : %f\n", cost,newcost);
		}
		//printf("------------------------------------------------ suc[%d] = %d \n", 0,suc[0]);
		cost = (count == 1) ? cost : newcost;
		
		//printf("Tabu opt\n");
		tabuOpt(inst, suc, temp, tabulistout, tabulistin, etabu, thresh);

		//printf("Calcolo il flusso\n");
		fluxCalculator(temp, flux_best, inst->nturbines);

		//printf("Metto i cavi\n");
		cableregularize(inst, x_best, flux_best );

		//printf("Calcolo quanto costa la solutione\n");
		newcost = objectiveFunction(inst, x_best,flux_best);
		//printf("La soluzione trovata costa : %Le\n",(long double)newcost );

		//printf("Plotto\n");
		//plotGraph(inst, x_best, flux_best, newcost);

		//fprintf(gp, "exit\n");
		//fclose(gp);
		//gp = popen("gnuplot -p","w");

		
		if(newcost < best_cost)
		{
			best_cost = newcost;
			for(int h = 0; h < inst->nturbines; h++){best[h]=temp[h];}
			//printf("best[%d] = %d \n", 0,best[0]);
			printf("new cost %f\n",best_cost );
			plotGraph(inst, x_best, flux_best, newcost);

			/*fprintf(gp, "exit\n");
			fclose(gp);
			gp = popen("gnuplot -p","w");*/

		}
		//for(int h = 0; h < etabu; h++){printf("Tabu beetween ( %d ) - ( %d )\n",tabulistout[h], tabulistin[h]);}
		//printf("------------------------------------- TOT tabu list : %d \n",ntabu );
		//wait(1);
		for(int h = 0; h < inst->nturbines; h++){suc[h]=temp[h];}
		if(newcost < 1)return;
	}

	printf("Calcolo il flusso\n");
	fluxCalculator(best, flux_best, inst->nturbines);

	printf("Metto i cavi\n");
	cableregularize(inst, x_best, flux_best );

	printf("Calcolo quanto costa la solutione\n");
	best_cost = objectiveFunction(inst, x_best, flux_best);
	inst->zbest = best_cost;
	printf("La soluzione trovata costa : %Le\n",(long double)best_cost );

	printf("Plotto\n");

	plotGraph(inst, x_best, flux_best, best_cost);
	
	//free(suc);
	//free(best);
	//free(x_best);
	//free(flux_best);
}
/***************************************************************************************************************************/
void execute10(instance *inst, CPXENVptr env, CPXLPptr lp)
/***************************************************************************************************************************/
{

	int best[inst->nturbines];
	double flux_best[inst->nturbines*inst->nturbines];
	double x_best[inst->nturbines*inst->nturbines];
	double best_cost = DBL_MAX;

	for(int i = 1; (inst->tstart + inst->timelimit > second()); i++)
	{

		printf("------------------------------- Cerco la soluzione [ %d ] --------------------------------\n", i);
		int suc[inst->nturbines];
		double flux[inst->nturbines*inst->nturbines];
		double x[inst->nturbines*inst->nturbines];


		double newcost = DBL_MAX;

		int best_temp[inst->nturbines+1];
		double flux_best_temp[inst->nturbines*inst->nturbines];
		double x_best_temp[inst->nturbines*inst->nturbines];
	

		printf("Eseguo PrimDijkstra\n");
		inst->randomseed = time(NULL)*i;
		PrimDijkstraGrasp(inst->mat, inst->nturbines, suc, inst->randomseed);

		printf("Calcolo il flusso\n");
		fluxCalculator(suc, flux, inst->nturbines);
		
		printf("Metto i cavi\n");
		cableregularize(inst, x, flux );

		printf("Calcolo quanto costa la solutione\n");
		double cost = objectiveFunction(inst, x, flux);
		inst->zbest = cost;
		printf("La soluzione trovata costa : %Le\n",(long double)cost );

		for( int count = 1; cost != newcost; count++){
			for(int i = 0; i < inst->nturbines && count != 1; i++){suc[i] = best_temp[i];best_temp[i] = 0;}
			for(int i = 0; i < inst->nturbines*inst->nturbines; i++){x_best_temp[i] = 0;flux_best_temp[i] = 0;}
			
			printf("\n\n----------------------- RUN 1 OPT #%d ----------------------------------------\n", count );
			cost = (count == 1) ? cost : newcost;
			
			unoOpt(inst, suc, best_temp);

			//printf("Calcolo il flusso\n");
			fluxCalculator(best_temp, flux_best_temp, inst->nturbines);

			//printf("Metto i cavi\n");
			cableregularize(inst, x_best_temp, flux_best_temp );

			//printf("Calcolo quanto costa la solutione\n");
			newcost = objectiveFunction(inst, x_best_temp, flux_best_temp);
			//printf("La soluzione trovata costa : %Le\n",(long double)newcost );

			//printf("Plotto\n");
			//plotGraph(inst, x_best_tmp, flux_best_tmp, cost);

			//fprintf(gp, "exit\n");
			//fclose(gp);
			//gp = popen("gnuplot -p","w");
			//wait(1);
		}
		//plotGraph(inst, x_best_temp, flux_best_temp, newcost);

		fprintf(gp, "exit\n");
		fclose(gp);
		gp = popen("gnuplot -p","w");
		if(newcost < best_cost)
		{
			printf("new %lf\n", newcost);
			
			best_cost = newcost;
			inst->zbest = newcost;
			plotGraph(inst, x_best, flux_best, best_cost);
			for(int h = 0; h < inst->nturbines; h++){best[h] = best_temp[h];}
			for(int h = 0; h < inst->nturbines*inst->nturbines; h++){
				flux_best[h] = flux_best_temp[h];
				x_best[h] = x_best_temp[h];
			}

		}
		//wait(5);
		
	}
	//plotGraph(inst, x_best, flux_best, best_cost);

	fprintf(gp, "exit\n");
	fclose(gp);
	gp = popen("gnuplot -p","w");
	//-------------------------------------------------------------- uno opt --------------------
	


	printf("La soluzione migliore trovata ha costo : %Le\n",(long double )best_cost );
	//free(suc);
	//free(best);
	//free(x_best);
	//free(flux_best);
}
/***************************************************************************************************************************/
int execute11(instance *inst, CPXENVptr env, CPXLPptr lp)
/***************************************************************************************************************************/
{
	printf("huntscolonyOpt\n");
	double pheromones[inst->nturbines*inst->nturbines];
	
	int best[inst->nturbines];
	double best_flux[inst->nturbines*inst->nturbines];
	double best_x[inst->nturbines*inst->nturbines];
	double best_cost = DBL_MAX;

	int tmp[inst->nturbines];
	double tmp_flux[inst->nturbines*inst->nturbines];
	double tmp_x[inst->nturbines*inst->nturbines];
	double tmp_cost = 0;

	for(int i = 0; i<inst->nturbines*inst->nturbines;i++)
	{
		pheromones[i]=1;
		best_flux[i]=0;
		best_x[i]=0;
		tmp_flux[i]=0;
		tmp_x[i]=0;
	}
	for(int i = 0; i<inst->nturbines;i++)
	{
		best[i]=0;
		tmp[i]=0;
	}
	if(inst->times == -1)inst->times = INT_MAX;
	for(int i = 0 ;  (inst->tstart + inst->timelimit > second()); i++)
	{
		printf("\n\n----------------------------------------Iterazione %d--------------------------------\n",i );
		printf("----------------------------------------Tempo %.0f--------------------------------\n",(second()-inst->tstart) );
		if(i == INT_MAX-1)i=0;
		antFindPathKruskal(inst->mat, inst->nturbines, tmp, pheromones, inst);
		
		/*for(int h = 0; h < inst->nturbines; h++)
		{
			printf("tmp[ %d ] = %d\n",h,tmp[h] );
		}*/

		//printf("Calcolo il flusso\n");
		fluxCalculator(tmp, tmp_flux, inst->nturbines);

		//printf("Metto i cavi\n");
		cableregularize(inst, tmp_x, tmp_flux );

		//printf("Calcolo quanto costa la solutione: ");
		tmp_cost = objectiveFunction(inst, tmp_x, tmp_flux);


		//plotGraph(inst, tmp_x, tmp_flux, tmp_cost);
		//plotGraph(inst, tmp_x, tmp_flux, tmp_cost);

		//printf(": %f\nAggiorno i feromoni\n",tmp_cost);
		updatepheromones(pheromones, tmp_x, tmp_flux, inst->cablecost, inst->nturbines, inst->mat, tmp_cost);

		if(tmp_cost < best_cost)
		{
			for(int j = 0; j < inst->nturbines; j++){
				best[j]=tmp[j];
			}
			fluxCalculator(best, best_flux, inst->nturbines);
			cableregularize(inst, best_x, best_flux );
			best_cost = objectiveFunction(inst, best_x, best_flux);
			inst->zbest=best_cost;

			printf("NEW BEST SOLUTION : %lf\n",best_cost );
			//plotGraph(inst, best_x, best_flux, best_cost);
			/*
			printf("Plotto migliore soluzione trovata\n");
			fprintf(gp, "exit\n");
			fclose(gp);
			gp = popen("gnuplot -p","w");
			*/
		}
		//wait(1);
		//exit(0);
	}

	//printf("Plotto\n");
	plotGraph(inst, best_x, best_flux, best_cost);
	fprintf(gp, "exit\n");
	fclose(gp);
	gp = popen("gnuplot -p","w");

	return 0;
}
/***************************************************************************************************************************/
int execute12(instance *inst, CPXENVptr env, CPXLPptr lp)
/***************************************************************************************************************************/
{
	printf("huntscolonyOpt\n");
	double pheromones[inst->nturbines*inst->nturbines];
	
	int best[inst->nturbines];
	double best_flux[inst->nturbines*inst->nturbines];
	double best_x[inst->nturbines*inst->nturbines];
	double best_cost = DBL_MAX;

	int tmp[inst->nturbines];
	double tmp_flux[inst->nturbines*inst->nturbines];
	double tmp_x[inst->nturbines*inst->nturbines];
	double tmp_cost = 0;

	for(int i = 0; i<inst->nturbines*inst->nturbines;i++)
	{
		pheromones[i]=1;
		best_flux[i]=0;
		best_x[i]=0;
		tmp_flux[i]=0;
		tmp_x[i]=0;
	}
	for(int i = 0; i<inst->nturbines;i++)
	{
		best[i]=0;
		tmp[i]=0;
	}
	if(inst->times == -1)inst->times = INT_MAX;
	for(int i = 0 ;  (inst->tstart + inst->timelimit > second()); i++)
	{
		printf("\n\n----------------------------------------Iterazione %d--------------------------------\n",i );
		printf("----------------------------------------Tempo %.0f--------------------------------\n",(second()-inst->tstart) );
		if(i == INT_MAX-1)i=0;
		antFindPathKruskal(inst->mat, inst->nturbines, tmp, pheromones, inst);
		
		/*for(int h = 0; h < inst->nturbines; h++)
		{
			printf("tmp[ %d ] = %d\n",h,tmp[h] );
		}*/

		//printf("Calcolo il flusso\n");
		fluxCalculator(tmp, tmp_flux, inst->nturbines);

		//printf("Metto i cavi\n");
		cableregularize(inst, tmp_x, tmp_flux );

		//printf("Calcolo quanto costa la solutione: ");
		tmp_cost = objectiveFunction(inst, tmp_x, tmp_flux);


		int best_temp[inst->nturbines];
		double flux_best_temp[inst->nturbines*inst->nturbines];
		double x_best_temp[inst->nturbines*inst->nturbines];
		double newcost = 0;

		for( int count = 1; tmp_cost != newcost; count++){
			for(int i = 0; i < inst->nturbines && count != 1; i++){tmp[i] = best_temp[i];best_temp[i] = 0;}
			for(int i = 0; i < inst->nturbines*inst->nturbines; i++){x_best_temp[i] = 0;flux_best_temp[i] = 0;}
			
			printf("\n\n----------------------- RUN 1 OPT #%d ----------------------------------------\n", count );
			tmp_cost = (count == 1) ? tmp_cost : newcost;
			
			unoOpt(inst, tmp, best_temp);

			//printf("Calcolo il flusso\n");
			fluxCalculator(best_temp, flux_best_temp, inst->nturbines);

			//printf("Metto i cavi\n");
			cableregularize(inst, x_best_temp, flux_best_temp );

			//printf("Calcolo quanto costa la solutione\n");
			newcost = objectiveFunction(inst, x_best_temp, flux_best_temp);
		}
		
		//plotGraph(inst, tmp_x, tmp_flux, tmp_cost);

		//printf(": %f\nAggiorno i feromoni\n",tmp_cost);
		updatepheromones(pheromones, tmp_x, tmp_flux, inst->cablecost, inst->nturbines, inst->mat, tmp_cost);

		if(tmp_cost < best_cost)
		{
			for(int j = 0; j < inst->nturbines; j++){
				best[j]=tmp[j];
			}
			fluxCalculator(best, best_flux, inst->nturbines);
			cableregularize(inst, best_x, best_flux );
			best_cost = objectiveFunction(inst, best_x, best_flux);
			inst->zbest=best_cost;

			
			printf("Plotto migliore soluzione trovata\n");
			plotGraph(inst, best_x, best_flux, best_cost);
			fprintf(gp, "exit\n");
			fclose(gp);
			gp = popen("gnuplot -p","w");
			
		}
		//wait(1);
		//exit(0);
	}

	//printf("Plotto\n");
	plotGraph(inst, best_x, best_flux, best_cost);
	fprintf(gp, "exit\n");
	fclose(gp);
	gp = popen("gnuplot -p","w");

	return 0;
}

/***************************************************************************************************************************/
void execute(instance *inst, CPXENVptr env, CPXLPptr lp)
/**************************************************************************************************************************/
{    
	
	if(inst->model_type == 0)
	{
		printf("Execution with Cplex\n");
		switch (inst->noCross) 	
		{
			case 0 : 								
				execute0(inst, env,lp); // normal execution
				break;
			case 1 : 								
				execute1(inst, env,lp); // normal execution + cplex lazy constraints
				break;
			case 2 : 								
				execute2(inst, env,lp); // loop method
				break;
			case 3 : 								
				execute3(inst, env,lp); // with lazy callback
				break;
			case 4 : 								
				execute4(inst, env,lp); // hard fixing
				break;
			case 5 : 								
				execute5(inst, env,lp); // soft fixing
				break;
			default: 
				print_error(" execution type unknown!!"); 
				break;
		}
	} 
	else if(inst->model_type == 1)
	{
		switch (inst->noCross) 	
		{
			case 6:
				execute6(inst, env, lp);	// heuristic
				break;
			case 7:
				execute7(inst,env,lp);		// heuristic loop to have more solution
				break;
			case 8:
				execute8(inst, env, lp);	// heuristic 1 Opt
				break;
			case 9:
				execute9(inst, env, lp);	// heuristic tabu search
				break;
			case 10:
				execute10(inst, env, lp);	// heuristic multi start
				break;
			case 11:
				execute11(inst, env, lp);	// heuristic ant colony
				break;
			case 12:
				execute12(inst, env, lp);	// heuristic ant colony
				break;
		}
	}
	
}
static int CPXPUBLIC lazyCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p)
{
	*useraction_p = CPX_CALLBACK_DEFAULT;
	instance *inst = (instance*)cbhandle;
	double z = CPX_INFBOUND;
	CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &z);
	double xstar[inst->ncols];
	if ( CPXgetcallbacknodex(env, cbdata, wherefrom, xstar, 0, inst->ncols-1) )
		return 1;
	
	int ncuts = myseparation(inst, xstar, env, cbdata, wherefrom);
	printf("Cuts added : < %d >\n",ncuts );
	if(ncuts >= 1)
	{
		//plotGraph(inst, xstar);
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
static int CPXPUBLIC heuristicCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, double *objval_p, double *x, int *checkfeas_p, int *useraction_p)
{
	//printf("HEURISTIC CALLBACK START\n");
	instance *inst = (instance*)cbhandle;
	*checkfeas_p = 0;
	*useraction_p = CPX_CALLBACK_DEFAULT;
	long double zbest; 
	CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_INTEGER, &zbest);
	//if(zbest < 1) return 0;
	if(inst->cableRegFA == 0 )
	{	
		int count = cableregularize(inst, x, zbest);
		if(count > 0)
		{
			printf("Cables switched < %d >\n",count);
			*checkfeas_p = 1;
			*useraction_p = CPX_CALLBACK_SET;
		}
		inst->cableRegFA = inst->cableRegF; 
	}
	else
	{
		inst->cableRegFA--; 
	}
	return 0;
}

int installheuristicCallback( CPXENVptr env, CPXLPptr lp, instance *inst)
{
	if(CPXsetheuristiccallbackfunc(env, heuristicCallback, (void *)inst))
	{
		print_error("Error heuristic Callback function");
	}
	
		
	return 0;
}

int nocross_separation(CPXENVptr env, CPXLPptr lp, instance *inst)
{
	int count = 0;
	double x[CPXgetnumcols(env, lp)]; 
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
				values = (double *) calloc(inst->nturbines + 3, sizeof(double));
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
	printf("count %d\n",count );
	return count;
}

int setHardFixRandom(CPXENVptr env, CPXLPptr lp, instance *inst, int *index,char *lu, double *bd)
{
	int count = 0;
	for(int i = 0; i < inst->nturbines; i++)
		for(int j = 0; j < inst->nturbines; j++)
		{
			if(i == j) continue;
			
			if((inst->best_sol[ypos(i,j,inst)] > 0.5 ) && (rand() % 100) > 50)
			{
				printf("Fixing beetween ( %d ) - ( %d )\n", i,j);
				index[count] = ypos(i,j,inst);
				lu[count] = 'L';
				bd[count] = 1;
				inst->best_sol[ypos(i,j,inst)] = 1.0;
				count++;

			}
		}
	
	CPXchgbds (env, lp, count, index, lu, bd);		
	return count;
}
int setHardRins(CPXENVptr env, CPXLPptr lp, instance *inst, int *index,char *lu, double *bd)
{
	printf("Hard Fixing - Rins\n");
	int count = 0;
	for(int i = 0; i < inst->nturbines; i++)
		for(int j = i + 1; j < inst->nturbines; j++)
		{
			if(i == j) continue;
			
			if((inst->best_sol[ypos(i,j,inst)] > 0.5 ) && (inst->second_best_sol[ypos(i,j,inst)] > 0.5 ) && (rand() % 100) > 50)
			{
				//printf("Fixing beetween ( %d ) - ( %d )\n", i,j);
				index[count] = ypos(i,j,inst);
				lu[count] = 'U';
				bd[count] = 1;
				inst->best_sol[ypos(i,j,inst)] = 1.0;
				count++;
			}
			else if((inst->best_sol[ypos(i,j,inst)] < 0.5 ) && (inst->second_best_sol[ypos(i,j,inst)] < 0.5 ) && (rand() % 100) > 50)
			{
				//printf("Fixing beetween ( %d ) - ( %d )\n", i,j);
				index[count] = ypos(i,j,inst);
				lu[count] = 'L';
				bd[count] = 0;
				inst->best_sol[ypos(i,j,inst)] = 0.0;
				count++;				
			}
		}
	
	CPXchgbds (env, lp, count, index, lu, bd);
	
	
	return count;
}
int setHardFix(CPXENVptr env, CPXLPptr lp, instance *inst, int *index,char *lu, double *bd)
{
	int count = 0;
	if(inst->hardF == 1)
		count = setHardFixRandom(env, lp, inst, index, lu, bd);
	else if(inst->hardF == 2)
		count = setHardRins(env, lp, inst, index, lu, bd);
	printf("Setting of  < %d > edges\n", count );
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
int setSoftFixA(CPXENVptr env, CPXLPptr lp, instance *inst, int K)		//Sum (for every (i,j) of yr = 1)( 1 - yij )>= (n - 1) - K
{
	printf("Asimmetric Local Branching with k = < %d >\n",K);
	int l = xpos(0,0,0,inst);
	int *index = (int *) calloc(inst->nturbines * inst->nturbines, sizeof(int ));
	int count = asimmetricHamming(inst->best_sol,ypos(0,0,inst), l, index);

	char **cname = (char **) calloc(1, sizeof(char *));		// (char **) required by cplex...
	cname[0] = (char *) calloc(100, sizeof(char));


	int lastrow = CPXgetnumrows(env,lp);
	double rhs = count - 1 - K ; 
	char sense = 'G'; 
	sprintf(cname[0], "Local Branching");   
	if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) 
		print_error(" wrong CPXnewrows [x1]");
	for ( int i = 0; i < count; i++ )
	{

		if ( CPXchgcoef(env, lp, lastrow, index[i], 1.0) ) 
			print_error(" wrong CPXchgcoef [x1]");
	}	

	free(index);
	free(cname);
	return count;
}
int setSoftFixS(CPXENVptr env, CPXLPptr lp, instance *inst, int K)	//Sum (for every (i,j) of yr = 0)( yij ) + Sum (for every (i,j) of yr = 1)( 1 - yij )<= K 
{
	printf("Simmetric Local Branching with k = < %d >\n",K);
	int l = xpos(0,0,0,inst);
	int *index = (int *) calloc(inst->nturbines * inst->nturbines, sizeof(int ));
	int count = asimmetricHamming(inst->best_sol,ypos(0,0,inst), l, index);

	char **cname = (char **) calloc(1, sizeof(char *));		// (char **) required by cplex...
	cname[0] = (char *) calloc(100, sizeof(char));

	int lastrow = CPXgetnumrows(env,lp);
	double rhs = count - K ; 
	char sense = 'L'; 
	sprintf(cname[0], "Local Branching");   
	if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) 
		print_error(" wrong CPXnewrows [x1]");
	for ( int i = 0; i < count; i++ )
	{
		if ( CPXchgcoef(env, lp, lastrow, index[i], 1.0) ) 
			print_error(" wrong CPXchgcoef [x1]");

	}
	free(index);
	free(cname);
	return count;
}
int setSoftRinsA(CPXENVptr env, CPXLPptr lp, instance *inst, int K)
{
	printf("Rins Local Branching with k = < %d >\n",K);
	int l = xpos(0,0,0,inst);
	int *index = (int *) calloc(inst->nturbines * inst->nturbines, sizeof(int ));
	int count = RinsAHamming(inst->best_sol, inst->second_best_sol,ypos(0,0,inst), l, index);

	char **cname = (char **) calloc(1, sizeof(char *));		// (char **) required by cplex...
	cname[0] = (char *) calloc(100, sizeof(char));


	int lastrow = CPXgetnumrows(env,lp);
	double rhs = count - 1 - K ; 
	char sense = 'G'; 
	sprintf(cname[0], "Local Branching");   
	if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) 
		print_error(" wrong CPXnewrows [x1]");
	for ( int i = 0; i < count; i++ )
	{

		if ( CPXchgcoef(env, lp, lastrow, index[i], 1.0) ) 
			print_error(" wrong CPXchgcoef [x1]");
	}	

	free(index);
	free(cname);
	return count;
}
int setSoftRinsS(CPXENVptr env, CPXLPptr lp, instance *inst, int K)	
{
	printf("Simmetric Local Branching with k = < %d >\n",K);
	int l = xpos(0,0,0,inst);
	int *index = (int *) calloc(inst->nturbines * inst->nturbines, sizeof(int ));
	int count = RinsSHamming(inst->best_sol,inst->second_best_sol,ypos(0,0,inst), l, index);

	char **cname = (char **) calloc(1, sizeof(char *));		// (char **) required by cplex...
	cname[0] = (char *) calloc(100, sizeof(char));

	int lastrow = CPXgetnumrows(env,lp);
	double rhs = count - K; 
	char sense = 'L'; 
	sprintf(cname[0], "Local Branching");   
	if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) 
		print_error(" wrong CPXnewrows [x1]");
	for ( int i = 0; i < count; i++ )
	{
		if ( CPXchgcoef(env, lp, lastrow, index[i], 1.0) ) 
			print_error(" wrong CPXchgcoef [x1]");

	}
	free(index);
	free(cname);
	return count;
}
int setSoftFix(CPXENVptr env, CPXLPptr lp, instance *inst, int K)
{
	int res = 0;
	if(inst->softF == 1)
		res = setSoftFixA( env,  lp, inst, K);
	else if(inst->softF == 2)
		res = setSoftFixS( env,  lp, inst, K);
	else if(inst->softF == 3)
		res = setSoftRinsA(env, lp, inst, K);
	else if(inst->softF == 4)
		res = setSoftRinsS(env, lp, inst, K);
	CPXwriteprob(env, lp, "model/model_softFix.lp", NULL);
	return res;
}
int resetSoftFix(CPXENVptr env, CPXLPptr lp)
{
	int begin = CPXgetnumrows(env,lp)-1;
	int end = begin;	
	CPXdelrows(env, lp, begin,  end);

	return 0;
}
int PrimDijkstraMat(instance *inst)
{
	for(int i = 0; i < inst->nturbines; i++)
	{
		for(int j = 0; j < inst->nturbines; j++)
		{
			if(i == j)
			{
				inst->mat[i+j*inst->nturbines] = DBL_MAX; 
			}
			else
			{
				inst->mat[i+j*inst->nturbines] = dist(i,j,inst);
			}
			//printf("[ %d ] [ %d ]",i,j);
		}
	}
	return 0;
}


int unoOpt(instance *inst, int* suc, int* best)
{
	//printf("unoOpt\n");
	int predec[inst->nturbines];

	int newsuc[inst->nturbines];

	double flux[inst->nturbines*inst->nturbines];
	
	double x[inst->nturbines*inst->nturbines];
	
	double minobj = inst->zbest;
	
	for(int i = 0; i < inst->nturbines; i++){best[i] = suc[i];}
	for(int i = 0; i < inst->nturbines; i++){newsuc[i] = suc[i];}

	for(int i = 0; i < inst->nturbines; i++){
		if(suc[i] != -1){
			//printf("Ramo uscente da ----------------< %d >-----------------------\n",i);
			for(int h = 0; h < inst->nturbines; h++){predec[h] = 0;}
			predec[0] = i;
			int npredec = 0;
			int iter = 0;
			int flag = 1;
			//printf("<- %d -",predec[0] );
			for(int k = 0; flag > 0;k++){
				flag--;
				for(int j = 0; j < inst->nturbines; j++){
					if(suc[j] == predec[iter])
					{
						flag ++;
						npredec++;
						predec[npredec] = j;
						//printf("- %d -",predec[npredec]);
					}
				}
				iter++;
			}
			//printf(">\n");
			flag = 0;
			npredec++;

			for(int j = 0; j < inst->nturbines; j++){
				flag = 0;

				if(i == j || arraycontainsINT(predec, npredec, j) == 1)
				{
					//printf("Non posso collegarlo al ramo : %d\n",j);
					flag = 1;
				}

				if(suc[i] != -1 && flag == 0){
					//printf("\tRamo entrante in ----------------< %d >-----------------------\n",j);
					newsuc[i] = j;
					//printf("Flusso\n");
					//printf("prima di flusso predec[10] = %d\n", predec[10]);
					fluxCalculator(newsuc, flux, inst->nturbines);
					//printf("dopo flusso predec[10] = %d\n", predec[10]);
					//printf("Cavi\n");
					cableregularize(inst, x, flux );
					//printf("Costo\n");
					double cost = objectiveFunction(inst, x, flux);
					if(cost < minobj){
						inst->zbest = cost;
						inst->second_zbest = minobj;
						//printf("New incumbment: %f \nOld incumbment: %f\n\n", cost,minobj );
						minobj = cost;
						for(int h = 0; h < inst->nturbines; h++){best[h] = newsuc[h];}
						//plotGraph(inst, x, flux, cost);
						//wait(5);
					}

					for(int h = 0; h < inst->nturbines*inst->nturbines; h++){x[h] = 0;}
					for(int h = 0; h < inst->nturbines*inst->nturbines; h++){flux[h] = 0;}
					for(int h = 0; h < inst->nturbines; h++){newsuc[h] = suc[h];}
			
				}
			}
			
		}
	}

	return 0;
}

int tabuOpt(instance *inst, int* suc, int* best, int* tabulistout, int* tabulistin, int &etabu, int tabuthresh)
{
	printf("tabuopt\n");
	//printf("------------------------------------------------ suc[%d] = %d \n", 0,suc[0]);
	int *predec;//[inst->nturbines];
	predec = (int *) calloc(inst->nturbines, sizeof(int));

	int *newsuc;//[inst->nturbines];
	newsuc = (int *) calloc(inst->nturbines, sizeof(int));

	double *flux;//[inst->nturbines*inst->nturbines];
	flux = (double *) calloc((inst->nturbines*inst->nturbines), sizeof(double));
	
	double *x;//[inst->nturbines*inst->nturbines]; 
	x = (double *) calloc((inst->nturbines*inst->nturbines), sizeof(double));
	
	double minobj = DBL_MAX;
	
	for(int i = 0; i < inst->nturbines; i++){best[i] = suc[i];}
	for(int i = 0; i < inst->nturbines; i++){newsuc[i] = suc[i];}

	for(int i = 1; i < inst->nturbines; i++){
		//printf("suc[%d] = %d \n", 0,suc[0]);
		if(suc[i] != -1){
			//printf("Ramo uscente da ----------------< %d >-----------------------\n",i);
			for(int h = 0; h < inst->nturbines; h++){predec[h] = 0;}
			predec[0] = i;
			int npredec = 0;
			int iter = 0;
			int flag = 1;
			//printf("<- %d -",predec[0] );
			for(int k = 0; flag > 0;k++){
				flag--;
				for(int j = 0; j < inst->nturbines; j++){
					if(suc[j] == predec[iter])
					{
						flag ++;
						npredec++;
						predec[npredec] = j;
						//printf("- %d -",predec[npredec]);
					}
				}
				iter++;
			}
			//printf(">\n");
			flag = 0;
			npredec++;

			for(int j = 0; j < inst->nturbines; j++){
				flag = 0;
				if( i == j || suc[i] == j )
				{
					//printf("Non posso collegarlo al ramo : %d\n",j);
					flag = 1;
				}else if(arraycontainsINT(predec, npredec, j) == 1 )
				{
					//printf("Non posso collegarlo al ramo : %d perchè creerebbe un ciclo\n",j);
					flag = 1;
				}else if( tabucontains(tabulistout, tabulistin, tabuthresh, i ,j) == 1)
				{
					//printf("Non posso collegarlo al ramo : %d prechè è tabu\n",j);
					flag = 1;
				}
				if(suc[i] != -1 && flag == 0){
					//printf("\tRamo entrante in ----------------< %d >-----------------------\n",j);
					newsuc[i] = j;
					//printf("Flusso\n");
					fluxCalculator(newsuc, flux, inst->nturbines);
					//printf("Cavi\n");
					cableregularize(inst, x, flux );
					//printf("Costo\n");
					double cost = objectiveFunction(inst, x, flux);
					if(cost < minobj){
						inst->zbest = cost;
						inst->second_zbest = minobj;
						//printf("New incumbment: %f \nOld incumbment: %f\n\n", cost,minobj );
						minobj = cost;
						for(int h = 0; h < inst->nturbines; h++){best[h] = newsuc[h];}
						//plotGraph(inst, x, flux, cost);
						//wait(5);
					}

					for(int h = 0; h < inst->nturbines*inst->nturbines; h++){x[h] = 0;}
					for(int h = 0; h < inst->nturbines*inst->nturbines; h++){flux[h] = 0;}
					for(int h = 0; h < inst->nturbines; h++){newsuc[h] = suc[h];}
			
				}
			}
			
		}
	}
	for(int h = 0; h < inst->nturbines; h++)
	{
		if(best[h] != suc[h])
		{	
			tabulistout[etabu] = h;
			tabulistin[etabu] = best[h];
			etabu++;
			if(etabu > tabuthresh)
			{
				etabu = 0;
			}
			continue;
		}
	}
	//printf("------------------------------------------------ suc[%d] = %d \n", 0,suc[0]);
	free(newsuc);
	free(predec);
	free(flux);
	free(x);
	return 0;
}
