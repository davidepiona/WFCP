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
*//*

	Nuovo modello - soft fixing - euristico

			________________________
  input	--->|	     Algoritmo  	| --->	Soluzione 
	+	--->|         esatto        | --->   ottima
 vincolo--->|_______________________| --->

	il vincolo dice di fissare almeno il 90% degli archi nella soluzione di riferimento

	Soluzione di riferimento : yr = (...)
	Data una soluzione generica y = (...)
	Calcolo la distanza di hamming d
	il numero di flip consentiti è limitato ( <= K)


	E' come guardare solo l'intorno della soluzione generica y
	>>>>>>>>>>>>>>>>>>>> Local Branching >>>>>>>>>>>>> Vincolo simmetrico: siccome vengono contati sia i flip da 0 a 1 che da 1 a 0
										 >>>>>>>>>>>>> Vincolo assimmetrico: vengono contati solo i flip da 1 a 0

	vincolo simmetrico
	Sum (for every (i,j) of yr = 0)( yij ) + Sum (for every (i,j) of yr = 1)( 1 - yij )<= K 
	
	Vincolo Assimmetrico
	//Sum (for every (i,j) of yr = 1)( 1 - yij )>= (n - 1) - K

	Il valore di K va variato secondo politiche di "buon senso" 
		- partendo da k molto piccoli ( 3 / 5, quando non riesco a migliorare la soluzione aumento, massimo 20)
		//- partendo da k molto grandi\\

	loop iniziale
		prendiamo la prima soluzione 
		cambiamo i bound data quella soluzione di riferimento
		fino ad un time limit
	ripeti
	

	Heuristic 
	---> constractive
	---> refining

	Constractive : Greedy
	---> l'algoritmo di kruskal è un algoritmo greedy
	---> il greedy funziona sempre con strutture dati "matroidi"
	
	TSP - Greedy
	Nel tsp io ho un certo numero di città da visitare e mi trovo nella città 1, il primo passo è andare nella città più vicina, e avanti così.
	Funziona con gli animali!

	TO DO
	- Prim-Dijkstra 
	- Regolarizzare flusso
	- scelta cavi (cable regularize)
	- Calcolo funzione obbiettivo + M1(#edges in substation - C )  + M2(#crossing) : M1 ~ 10e9 , M2 ~ 10e8

	Si può provare a far partire l'algoritmo da punti diversi, MULTI-START ---> Da noi non può funzionare
	Si può provare a randomizzare leggermente i costi ( ovvero le distanze )
	Si può provare con il "GRASP" : Greedy randomize adaptive ... 
	Dove è lecito scegliere non scelgo sempre la scelta migliore, altrimenti avrei il greedy, ma scelgo le 5 scelte migliori ad ogni iterazione
	e poi con il 50% di probabilità sclego la scelta migliore o una a caso delle altre
*/