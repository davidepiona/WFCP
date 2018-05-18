#include <cplex.h>
#include <time.h>
#include "wfcp.h"
#include "float.h" 
#include <stdlib.h>

double second();  
void print_error(const char *err);

int plotGraph(CPXENVptr env, CPXLPptr lp, instance *inst);
int plotGraph(instance *inst, double *x);

void wait(int time)  
{	
	double t = second();	
	while(second() - t < time);
}

int xpos(int i, int j, int k, instance *inst){ 	return inst->xstart + i * (inst->nturbines * inst->ncables) + j * inst->ncables + k; }                                
int ypos(int i, int j, instance *inst) {	return inst->ystart + i * inst->nturbines + j; }                                              
int fpos(int i, int j, instance *inst) { 	return inst->fstart + i * inst->nturbines + j; }    
int spos(int s, instance *inst) { 	return inst->sstart + s; }  

int noCross(int p1, int p2, int p3, int p4, instance *inst)
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
			return 0;
	}
	return 1;
	//printf("Non C'è Crossing\n");
	
}

double dist(int i, int j, instance *inst)
{
	double dx = inst->xcoord[i] - inst->xcoord[j];
	double dy = inst->ycoord[i] - inst->ycoord[j]; 
	return sqrt(dx*dx+dy*dy);
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
	//printf("----------------------------------------------------time %lf tspan %lf\n",time,tspan );
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
	if ( residual_time < 0.0 ) 
		residual_time = 0.0;
	//CPXsetintparam(env, CPX_PARAM_CLOCKTYPE, 2);
	CPXsetdblparam(env, CPX_PARAM_TILIM, residual_time); 							// real time
	CPXsetdblparam(env, CPX_PARAM_DETTILIM, TICKS_PER_SECOND*timelimit);			// ticks
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
	if ( mip_value(env,lp) < inst->second_zbest - EPSILON )
	{
		printf("------------------------------- ADD SOLUTION\n");
		if( mip_value(env, lp) < inst->zbest)
		{
			inst->second_zbest = inst->zbest;
			inst->second_best_sol = inst->best_sol;

			inst->tbest = second() - inst->tstart;
			inst->zbest = mip_value(env, lp);
			CPXgetx(env, lp, inst->best_sol, 0, inst->ncols-1);
		}
		else
		{
			inst->second_zbest = mip_value(env, lp);
			CPXgetx(env, lp, inst->second_best_sol, 0, inst->ncols-1);	
		}

		if ( VERBOSE >= 40 ) 
			printf("\n >>>>>>>>>> incumbent update of value %lf at time %7.2lf , the second best solution is %lf <<<<<<<<\n", inst->zbest, inst->tbest, inst->second_zbest);
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
	if ( z < inst->second_zbest )
	{
		printf("------------------------------- ADD SOLUTION\n");
		
		if(z < inst->zbest )
		{
			inst->second_zbest = inst->zbest;
			inst->second_best_sol = inst->best_sol; 	
			
			inst->zbest = z;
			inst->best_sol = x;
			inst->tbest = second() - inst->tstart;	
		}
		else
		{
			inst->second_zbest = z;
			inst->second_best_sol = x; 	
		}

		if ( VERBOSE >= 40 ) 
			printf("\n >>>>>>>>>> incumbent update of value %lf at time %7.2lf , the second best solution is %lf <<<<<<<<\n", inst->zbest, inst->tbest, inst->second_zbest);
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

double mip_getobj(CPXENVptr env, CPXLPptr lp, const double *x)
{
	int ncols = CPXgetnumcols(env, lp);
	double *obj = (double *) calloc(ncols, sizeof(double));
	CPXgetobj(env, lp, obj, 0, ncols-1);
	double value = 0.0;
	for ( int j = 0; j < ncols; j++ ) value += obj[j] * x[j];
	free(obj);
	return value;
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
		if(!noCross(i,j,k,l, inst))
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
		if(!noCross(i,j,k,l, inst))
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
int cableregularize(instance *inst, double *x, long double z )
{
	int count = 0;
	double sol = z;
	for(int i = 0 ; i < inst->nturbines ; i++)
	{
		for(int j = 0 ; j < inst->nturbines ; j++)
		{
			if(x[ypos(i,j,inst)] < 0.5) continue;
			double flux = x[fpos(i,j, inst)];
			if(flux < 1) continue;
			int cable = 0;
			for(int k = 0 ; k < inst->ncables ; k++)
			{
				if(x[xpos(i,j,k, inst)] > 0.5)
				{
					cable = k;
				}
			}
			for(int k = 0 ; k < inst->ncables  ; k++)
			{
				if(inst->cablecost[k] < inst->cablecost[cable] && flux < inst->cablecapacity[k])
				{
					x[xpos(i,j,cable, inst)] = 0.0;
					x[xpos(i,j,k, inst)] = 1.0;
					sol = sol - (dist(i,j,inst)*inst->cablecost[cable]) + (dist(i,j,inst)*inst->cablecost[k]);
					//printf("Cable switched < %d > with < %d > in edge ( %d , %d )\n",cable+1, k+1, i ,j);
					z = sol;
					cable = k;
					count ++;
				}
			}
		}
	}
	
	
	z = sol;

	return count;
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
void mip_reload_solution(CPXENVptr env, CPXLPptr lp, int ncols, instance *inst, double *xstar)
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

	cableregularize(inst,xstar,0);

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
int simmetricHamming(double *yr, int s, int l, int *index)
{
	int k = 0;

	for(int i = 0; i < l; i++)
	{
		if( yr[i] > 0.5 )
			k ++;
		
		index[i] = i+s;
		
	}
	return k;
}
int asimmetricHamming(double *yr, int s, int l, int *index)
{
	int k = 0;
	for(int i = 0; i < l; i++)
	{
		if( yr[i] > 0.5 )
		{
			index[k] = i+s;
			k++;
			//printf("add value 1 to cols %d - %d\n",i, l);
		}
		
	}
	return k;
}
int RinsAHamming(double *yr1, double *yr2, int s, int l, int *index)
{
	int k = 0;
	for(int i = 0; i < l; i++)
	{
		if( yr1[i] > 0.5 && yr2[i] > 0.5)
		{
			index[k] = i+s;
			k++;
			//printf("add value 1 to cols %d - %d\n",i, l);
		}
		
	}
	return k;
}
int RinsSHamming(double *yr1, double *yr2, int s, int l, int *index)
{
	int k = 0;

	for(int i = 0; i < l; i++)
	{
		if( yr1[i] > 0.5 && yr2[i] > 0.5 )
			k ++;
		
		index[i] = i+s;
		
	}
	return k;
}

int PrimDijkstra(double*mat, int nodes, int *pred, int r)
{
	int *flag;
	flag = (int *) calloc(nodes, sizeof(int));	
	double *L;
	L = (double *) calloc(nodes, sizeof(double));	
	double *P;
	P = (double *) calloc(nodes, sizeof(double));	
	
	flag[0] = 1;
	pred[0] = -1;
	P[0] = -600;
	if(r != 0)
	{
		srand(r);
	}

	for(int j = 1; j < nodes ; j++)
	{
		flag[j] = 0;
		L[j] = mat[0+j*nodes] + P[0];
		pred[j] = 0;
	}

	for(int k = 0; k < nodes-1; k++)
	{
		double min = DBL_MAX;
		int pool[nodes];
		int npool = 0;
		int h = 0;
		for(int j = 1; j < nodes; j++)
		{
			if(flag[j] == 0 && L[j] + P[j] < min)
			{
				min = L[j] + P[j];
				h = j;
			}
		}
		for(int j = 1; j < nodes; j++)
		{
			if(L[j] + P[j] == min )
			{
				pool[npool] = j;
				npool++;
			}
		}
		if( r != 0 && npool > 0)
		{
			h = pool[rand()%npool];
			P[h] = P[h] + 50;
		}
		flag[h] = 1;
		
		for(int j = 1; j < nodes; j++)
		{
			if(flag[j] == 0 && mat[h+j*nodes] < L[j])
			{
				L[j] = mat[h+j*nodes] + P[h];
				pred[j] = h;
				P[h] = P[h] - 10;
				P[j] = P[j] + 10;
			}
		}
	}
	free(flag);
	free(L);
	free(P);
	return 0;
}
int PrimDijkstraGrasp(double*mat, int nodes, int *pred, int r)
{
	int *flag;
	flag = (int *) calloc(nodes, sizeof(int));	
	double *L;
	L = (double *) calloc(nodes, sizeof(double));	
	int ngrasp = 10;

	flag[0] = 1;
	pred[0] = -1;
	if(r != 0)
	{
		srand(r);
	}

	for(int j = 1; j < nodes ; j++)
	{
		flag[j] = 0;
		L[j] = mat[0+j*nodes];
		pred[j] = 0;
	}

	for(int k = 0; k < nodes-1; k++)
	{
		double min = DBL_MAX;
		int pool[nodes];
		int npool = 0;
		int v = 0;
		int h = 0;
		for(int i = 0; i < ngrasp; i++)
		{
			h = 0;
			min = DBL_MAX;
			for(int j = 1; j < nodes; j++)
			{
				v = 0;
				for(int t = 0; t < npool; t++)
				{
					if(j == pool[t])
					{
						v = 1;
					}
				}
				if(flag[j] == 0 && L[j] < min && v != 1)
				{

					min = L[j];
					h = j;
				}
			}
			pool[npool] = h;
			npool++;
		}
		for(int t = 0; t < npool; t++)
		{
			printf("- %d -",pool[t]);
		}
		printf("\n----------------------------------------------------------------------\n");
		if( r != 0 && npool > 0)
		{
			int index = rand()%(2*npool);
			if(index < npool)
				h = pool[0];
			else
				h = pool[index - npool];

			//h = pool[rand()%npool];
		}
		flag[h] = 1;

		
		for(int j = 1; j < nodes; j++)
		{
			if(flag[j] == 0 && mat[h+j*nodes] < L[j])
			{
				L[j] = mat[h+j*nodes];
				pred[j] = h;
			}
		}
	}
	free(flag);
	free(L);
	return 0;
}

int fluxCalculator(int *suc, double*flux, int nodes)
{
	double *acc;
	acc = (double *) calloc(nodes, sizeof(double ));
	double *start;
	start = (double *) calloc(nodes, sizeof(double ));
	int n = 0;
	
	for( int i = 0; i < nodes; i++)
	{
		acc[suc[i]]++;
	}
	
	for(int j = 0; j < nodes-1; j++)
	{
		n = 0;
		for( int i = 0; i < nodes; i++)
		{
			if(acc[i] == j )
			{
				start[n] = i;
				n++;
			}
		}
		//printf("ci sono %d nodi con %d rami entranti\n",n,j);	
		for(int i = 0; i < n; i++)
		{
			int s = start[i];
			while(s != -1)
			{
				flux[s+suc[s]*nodes]++;
				s = suc[s];
			}
		}
		
	}/*
	for(int i= 0; i<nodes;i++)
		for(int j=0;j<nodes;j++)
			if(flux[j+i*nodes] > 0)printf("(j-i) (%d - %d) flux : %.0f\n",j,i, flux[j+i*nodes] );
	*/
	return 0;
}

int cableregularize(instance *inst, double*x, double*flux )
{
	int count = 0;
	int cable_max = 0;
	for(int k = 0 ; k < inst->ncables ; k++)
	{
		if(inst->cablecost[k] > inst->cablecost[cable_max])
			cable_max = k;
	}
	//printf("cable max %d\n", cable_max);
	for(int i = 0 ; i < inst->nturbines ; i++)
	{
		for(int j = 0 ; j < inst->nturbines ; j++)
		{
			if(flux[i+j*inst->nturbines] < 0.5) continue;

			int cable = cable_max;
			x[i+j*inst->nturbines] = cable;

			for(int k = 0 ; k < inst->ncables ; k++)
			{
				if((inst->cablecost[k] < inst->cablecost[cable]) && (flux[i+j*inst->nturbines] < inst->cablecapacity[k]) && (flux[i+j*inst->nturbines] > 0.5))
				{
					x[i+j*inst->nturbines] = k;
					//sol = sol - (dist(i,j,inst)*inst->cablecost[cable]) + (dist(i,j,inst)*inst->cablecost[k]);
					cable = k;
					count ++;
				}
			}
		}
	}
	
	
	
	return count;
}
