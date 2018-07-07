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


	Eusristici di raffinamento 
	Prende una soluzione e cerca di migliorarla in qualche modo, applicando una "mossa" (move).
	- 1-opt : prendo la soluzione, elimino un arco non valido e ne aggiungo uno valido, ci muovimo in un intorno della soluzione
			  si dice che le soluzioni trovate sono a distanza 1
			  il problema di questo procedimento sono gli ottimi locali
	- 2-opt : posso valutare due mosse in contemporanea per decidere qual'è la migliore 
	- k-opt : valuto k mosse

	MetaEuristici
	Metaeuristico : modo euristico per descrivere euristici
	Si puo uscire da un ottimo locale:
	- Multistart : provo a partire da soluzioni diverse andando a finire, molto probabilmente, in ottimi locali differenti
	- Tabu-Search : Posso fare anche mosse peggiorative se sono in un ottimo locale e tengo traccia di tutte le soluzioni precedenti in modo da non tornare indietro
					quando trovo la possibilità di scendere senza tornare indietro posso cancellare la tabu-list tenendo comunque traccia della soluzione migliore mai trovata.
					La tabu-list viene proposta come un oggetto dinamico che contiene una "tenuta" (tenure) ovvero dopo I iterazioni in cui una soluzione è presente nella lista
					viene rimossa, valore ragionevole 10/40 < del numero di turbine.
					Posso dichiarare tabu un vertice che ho appena cambiato al posto di tenere tutta la soluzione. Tengo un array dove memorizzo che il vertice i-esimo 
					è stato messo nella lista nell'iterazion k e guardo se l'iterazione corrente - k è > della tenure.
					Provo a far oscillare la tenure. Le fasi in cui si lascia poca libertà all'algoritmo (tenure << ) si chiamano di intensification mentre quelle con alta 
					libertà	(tenure >> ) diversification.
	- Simulated Annealing : Invece di cercare ogni volta la migliore soluzione, valuto un possibile scambio tra archi, questo scambio mi porta ad un delta-Costo e scelgo 
					la soluzione in cui muovermi con una certa probabilità ( vedi letteratura ) che dipende dalla temperatura e dal delta
	- Variable Neighborhood Search (VNS) : L'idea è di fare l'1-opt e quando arrivo ad un ottimo locale faccio una volta il 3-opt
	- Algoritmo formiche : 
							Wikipedia
							https://www.dbai.tuwien.ac.at/staff/musliu/ACO_online_final.pdf
							https://ieeexplore.ieee.org/abstract/document/585892/
							https://i11www.iti.kit.edu/_media/teaching/theses/ba-nedlin-17.pdf
	- Algoritmo genetico : Lavora con una popolazione di individui, ogni individuo è una soluzione. Nella prima generazione tutte le soluzioni vengono generate casualmente
					Noi speriamo che la seconda generazione sia mediamente un po migliore. La fitness è una misura di quanto una soluzione è buona ( per noi è - costo ).
					Passare da una soluzione ad una successiva: prendo due soluzioni e creo un figlio. Il patrimonio genetico dei genitori sono l'ordine di visita dei nodi, 
					definisco un punto di taglio e creo il corredo genetico del figlio. Va bene ammazzare figli in base alla loro fitness ma in modo probabilistico.
					Siccome potrei avere i figli "fatti male", devo fare una procedura di refinement in modo da sistemarli.







	Per la tesina parlare dell'articolo sull'ant algorithm su wind farm cable problem, interessa al professore.
	Dobbiamo parlare anche della regolarizzazione dei cavi.


	pytorch -> libreria deep learning

*/


