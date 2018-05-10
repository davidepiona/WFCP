#include "wfcp.h"           
#include "limits.h"

double second();
void print_error(const char *err);  
void check_dblsort();       
double random01();     
void read_input(instance *inst);
void read_cables(instance *inst);
void read_turbines(instance *inst);
void parse_command_line(int argc, char** argv, instance *inst); 

void debug(const char *err) 
{ 
	printf("\nDEBUG: %s \n", err); 
	fflush(NULL); 
}
void print_error(const char *err) 
{ 
	printf("\n\n ERROR: %s \n\n", err); 
	fflush(NULL); 
	exit(1); 
}   

void free_instance(instance *inst)
{     
	free(inst->power);
	free(inst->xcoord);
	free(inst->ycoord);
	free(inst->cablecapacity);
	free(inst->cablecost);
	free(inst->best_sol);
	free(inst->second_best_sol);
	free(inst->mat);
}

int CableOpt(instance *inst);    

int main(int argc, char **argv) 
{ 
	if( argc < 2 ) 
	{ 
		printf("Usage: %s -help for help\n", argv[0]); 
		exit(1); 
	}       
	if( VERBOSE >= 2 ) 
	{ 
		for(int a = 0; a < argc; a++) 
			printf("%s ", argv[a]); 
		printf("\n"); 
	}

	double t1 = second(); 
	instance inst;

	parse_command_line(argc,argv, &inst);     
	
	//printf(" file %s has %d non-empty lines\n", inst.input_file, number_of_nonempty_lines(inst.input_file)); exit(1);
	  
	read_input(&inst);  
	if( CableOpt(&inst) ) 
		print_error(" error within CableOpt()");
	double t2 = second(); 
    
	if( VERBOSE >= 1 )   
	{
		printf("... WFCP problem solved in %lf sec.s\n", t2-t1);  
	}
	
	free_instance(&inst);
	return 0; 
}         

void read_input(instance *inst) // simplified CVRP parser, not all SECTIONs detected  
{
	read_cables(inst);
	read_turbines(inst);
}
void read_cables(instance *inst)
{
                            
	FILE *fin = fopen(inst->cables_file, "r");
	if( fin == NULL ) print_error(" input file not found!");
	
	inst->ncables = -1;
	

	char line[180];
	int num = 0;

	while( fgets(line, sizeof(line), fin) != NULL ) 
	{
		num++;		    
	}         

	inst->ncables = num;
	inst->cablecapacity = (double *) calloc(inst->ncables, sizeof(double)); 	 
	inst->cablecost = (double *) calloc(inst->ncables, sizeof(double));      
	fclose(fin);  

	fin = fopen(inst->cables_file, "r");
	for(int i = 0; fgets(line, sizeof(line), fin) != NULL ; i++) 
	{
		sscanf(line,"%lf %lf", &inst->cablecapacity[i], &inst->cablecost[i]);		    
		if(VERBOSE >= 500)
		{	
			printf("(%d) %s", i, line);
			printf("Cap: %lf cost: %lf\n", inst->cablecapacity[i], inst->cablecost[i]);
		}
	}
	fclose(fin);
}
void read_turbines(instance *inst)
{
                            
	FILE *fin = fopen(inst->turbines_file, "r");
	if( fin == NULL ) print_error(" input file not found!");
	
	inst->nturbines = -1;
	

	char line[180];
	int num = 0;
	int p = 0;

	while( fgets(line, sizeof(line), fin) != NULL ) 
	{
		num++;		    
	}   
	inst->nturbines = num;
	inst->xcoord = (double *) calloc(inst->nturbines, sizeof(double)); 	 
	inst->ycoord = (double *) calloc(inst->nturbines, sizeof(double));    
	inst->power = (double *) calloc(inst->nturbines, sizeof(double)); 	
	fclose(fin); 
	
	fin = fopen(inst->turbines_file, "r");
	for(int i = 0; fgets(line, sizeof(line), fin) != NULL ; i++) 
	{
		sscanf(line,"%lf %lf %lf", &inst->xcoord[i], &inst->ycoord[i], &inst->power[i]);
		if(VERBOSE >= 500)
		{
			printf("(%d) %s\n", i, line);
			printf("X: %lf  Y:%lf  P:%lf\n\n", inst->xcoord[i], inst->ycoord[i], inst->power[i]);
		}
		if(inst->power[i] <= -0.5)
		{
			p = num;
		}
				    
	}  
	if(p == -1 )
	{
		print_error("Error, root not found");
	}  
	fclose(fin);
}
void parse_command_line(int argc, char** argv, instance *inst) 
{ 
	
	if( VERBOSE >= 100 ) printf(" running %s with %d parameters \n", argv[0], argc-1); 
		
	// default   
	strcpy(inst->cables_file, "NULL");
	strcpy(inst->turbines_file, "NULL");
	inst->timelimit = CPX_INFBOUND;
	inst->timeloop = CPX_INFBOUND;
	inst->timeStartSol = CPX_INFBOUND; 
	inst->model_type = 0;
	inst->xstart = -1;
	inst->ystart = -1;
	inst->fstart = -1;
	inst->sstart = -1;
	inst->C = INT_MAX;
	inst->rins = -1;
	inst->relax = 0;
	inst->polishing_time = -1.0;
	inst->randomseed = 0; 
	inst->num_threads = 0;
	inst->noCross = 1;
	inst->cableReg = 0;
	inst->cableRegF = 10;
	inst->cableRegFA = 0;
	inst->names = 0;
	inst->gap = 0.0;
	inst->softF = 1;
	inst->hardF = 1;

    int help = 0; if( argc < 1 ) help = 1;	
	for( int i = 1; i < argc; i++ ) 
	{ 
		if( strcmp(argv[i],"-cables_file") == 0 ) { strcpy(inst->cables_file,argv[++i]); continue; } 			// input cables file
		if( strcmp(argv[i],"-turbines_file") == 0 ) { strcpy(inst->turbines_file,argv[++i]); continue; } 		// input turbines file
		if( strcmp(argv[i],"-fc") == 0 )  { strcpy(inst->cables_file,argv[++i]); continue; }					// input cables file
		if( strcmp(argv[i],"-ft") == 0 ) { strcpy(inst->turbines_file,argv[++i]); continue; }					// input turbines file
		if( strcmp(argv[i],"-C") == 0 ) { inst->C = atoi(argv[++i]); continue; } 								// Capacity of root
		if( strcmp(argv[i],"-time_loop") == 0 ) { inst->timeloop = atof(argv[++i]); continue; }					// time for loop in loop/heuristic method
		if( strcmp(argv[i],"-time_limit") == 0 ) { inst->timelimit = atof(argv[++i]); continue; }				// total time limit
		if( strcmp(argv[i],"-time_start") == 0 ) { inst->timeStartSol = atof(argv[++i]); continue; }			// time start to heuristic method
		if( strcmp(argv[i],"-model_type") == 0 ) { inst->model_type = atoi(argv[++i]); continue; } 				// model type
		if( strcmp(argv[i],"-model") == 0 ) { inst->model_type = atoi(argv[++i]); continue; } 					// model type
		if( strcmp(argv[i],"-rins") == 0 ) { inst->rins = atoi(argv[++i]); continue; } 							// rins
		if( strcmp(argv[i],"-relax") == 0 ) { inst->relax = atoi(argv[++i]); continue; }						// relax
		if( strcmp(argv[i],"-polishing_time") == 0 ) { inst->polishing_time = atof(argv[++i]); continue; }		// polishing time
		if( strcmp(argv[i],"-gap") == 0 ) { inst->gap = atof(argv[++i]); continue; }							// gap to terminate
		if ( strcmp(argv[i],"-seed") == 0 ) { inst->randomseed = abs(atoi(argv[++i])); continue; } 				// random seed
		if ( strcmp(argv[i],"-threads") == 0 ) { inst->num_threads = atoi(argv[++i]); continue; } 				// n. threads
		if ( strcmp(argv[i],"-Cross_Constraints") == 0 ) { inst->noCross = abs(atoi(argv[++i])); continue;}		// Cross Constraints
		if ( strcmp(argv[i],"-CC") == 0 ) { inst->noCross = abs(atoi(argv[++i])); continue;}					// Cross Constraints
		if ( strcmp(argv[i],"-Cable_Regularization") == 0 ) { inst->cableReg = abs(atoi(argv[++i])); continue;}	// Cable regularization
		if ( strcmp(argv[i],"-CR") == 0 ) { inst->cableReg = abs(atoi(argv[++i])); continue;}					// Cable regularization
		if ( strcmp(argv[i],"-CRF") == 0 ) { inst->cableRegF = abs(atoi(argv[++i])); continue;}					// Cable regularization
		if ( strcmp(argv[i],"-soft_fix") == 0 ) { inst->softF = abs(atoi(argv[++i])); continue;}				// Type of soft fixing
		if ( strcmp(argv[i],"-hard_fix") == 0 ) { inst->hardF = abs(atoi(argv[++i])); continue;}				// Type of hard fixing
		if ( strcmp(argv[i],"-names") == 0 ) { inst->names = abs(atoi(argv[++i])); continue;}					// Names
		if( strcmp(argv[i],"-help") == 0 ) { help = 1; continue; } 												// help
		if( strcmp(argv[i],"--help") == 0 ) { help = 1; continue; } 											// help
		help = 1;
    }      

	if( help || (VERBOSE >= 10) )		// print current parameters
	{
		printf("\n\navailable parameters (vers. 0.1-2018) ---------------------------------------------------\n");
		printf("-cables_file %s\n", inst->cables_file);
		printf("-turbines_file %s\n", inst->turbines_file); 
		printf("-model_type %d\n", inst->model_type); 
		printf("-C %d\n", inst->C);
		printf("-rins %d\n", inst->rins);
		printf("-relax %d\n", inst->relax);
		printf("-polishing_time %lf\n", inst->polishing_time);
		printf("-gap %lf\n", inst->gap);
		printf("-seed %d\n", inst->randomseed);
		printf("-threads %d\n", inst->num_threads); 
		printf("-Cross_Constraints %d\n", inst->noCross);
		printf("-Cable_Regularization %d\n", inst->cableReg); 
		printf("-CRF %d\n", inst->cableRegF); 
		printf("-soft_fix %d\n", inst->softF); 
		printf("-hard_fix %d\n", inst->hardF); 
		printf("-Names %d\n", inst->names); 
		printf("-time_limit %lf\n", inst->timelimit); 
		printf("-time_loop %lf\n", inst->timeloop); 
		printf("-time_start %lf\n", inst->timeStartSol); 
		printf("\nenter -help or --help for help\n");
		printf("----------------------------------------------------------------------------------------------\n\n");
	}        
	
	if( help ) exit(1);

}    





