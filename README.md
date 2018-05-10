# WFCP
## INPUT STRING CONDITION
* **fc**                    : input cables file
* **ft**                    : input turbines file
* **C**                     : Capacity of root
* **time_loop**             : time for loop in loop/heuristic method
* **time_limit**            : total time limit
* **time_start**            : time start to heuristic method
* **model**                 : model type 
	* **0** : Cplex model
	* **1** : Matrix model
* **rins**                  : rins
* **relax**                 : relax 	
	* **1** : relax on station capacity
	* **2** : relax on flux
	* **3** : relax on flux + out edges
	* **else** : no relax
* **polishing_time**        : polishing time
* **gap**                   : gap to terminate
* **seed**                  : random seed
* **threads**               : n threads
* **CC**                    : Cross Constraints	
	* **1** : Lazy constraints to the model
	* **2** : loop Method
	* **3** : Normal execution + lazy callback
	* **4** : Hard Fixing
	* **5** : Soft Fixing
	* **else** : Normal Execution
	* **10** : Normal execution with no cross cable as normal constraints
* **soft_fix**              : Type of soft fixing 	
	* **1** : Asimmetric Local Branching
	* **2** : Simmetric Local Branching
	* **3** : rins asimmetric 
	* **4** : rins simmetric 													  
* **hard_fix**              : Type of hard fixing  
	* **1** : Random hard fixing
	* **2** : rins

## EXAMPLE OF STRING 
```
./wfcp -fc '/home/michele/Scrivania/RO2/WFCP/data/data_01.cbl' -ft '/home/michele/Scrivania/RO2/WFCP/data/data_01.turb' -C 10 -rins 5 -relax 3 -seed 9 -CC 5 -time_start 300 -time loop 30 -time_limit 1200 -CR 0 -soft_fix 2
```
## TESTING
* **CR**					: Cable regularization 
* **CRF**					: Cable regularization
