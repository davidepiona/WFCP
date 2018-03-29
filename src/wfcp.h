#ifndef WFCP_H_  

#define WFCP_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <stdio.h>  

#include <cplex.h>  
#include <pthread.h>  

#define VERBOSE				    50		// printing level  (=10 only incumbent, =20 little output, =50-60 good, =70 verbose, >=100 cplex log)

//hard-wired parameters

#define TICKS_PER_SECOND 	  1000.0  	// cplex's ticks on Intel Core i7 quadcore @2.3GHZ
#define EPSILON		  		  1e-9		// 1e-9		// very small numerical tolerance 
#define XSMALL		  		  1e-5 		// 1e-4*	// tolerance used to decide ingerality of 0-1 var.s
//data structures  

typedef struct {   
	
	//input data
	int nturbines; 	
	int ncables; 	
	double *power;   
	double *xcoord;
	double *ycoord;
	double *cablecapacity;
	double *cablecost;

	// parameters 
	int model_type; 
	double timelimit;						// overall time limit, in sec.s
	char cables_file[1000];		  			// input file
	char turbines_file[1000];	  			// input file
	int randomseed;							// random seed
	int rins;								// rins
	int relax;								// relax on nturbines in root
	double polishing_time;					// polishing time

	//global data
	double tstart;								
	double zbest;							// best sol. available  
	double tbest;							// time for the best sol. available  
	double *best_sol;						// best sol. available        
	double	best_lb;						// best lower bound available 
	
	// model;     
	int xstart;
	int ystart;
	int fstart;
	int sstart;
	int C;

} instance;        


#endif   /* VRP_H_ */ 
