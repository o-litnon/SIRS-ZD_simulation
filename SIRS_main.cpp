/*

To compile:		g++ SIRS_main.c -lm -lpthread


S = Susceptible to infection or immunisation.
I = currently Infected.
R = currently Recovered, and immune to change.
Z = currently Immune to infection.
D = currently dead.


S -> I with probability p1 if at least one neighbour of i is I; otherwise site i is unchanged.
S -> Z with probability p4 if at least one neighbour of i is Z; otherwise site i is unchanged.
I -> R with probability p2.
R -> S with probability p3.
I -> Z with probability p4 immediately after recovering from infection.
* -> D with probability p0.
D -> R with probability p2.

bash script:
declare -a p0values=(0.0)
declare -a p1values=(0.01 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99)
declare -a p4values=(0.0)
g++ SIRS_main.cpp -lm -lpthread; for i in ${p0values[@]}; do for ii in ${p1values[@]}; do for iii in ${p4values[@]}; do ./a.out -g -n 10 -m 10 -0 $i -2 $ii -4 $iii; done; done; done;

*/

#include <stdio.h>
#include <time.h>

#include "SIRS_simulation_control.h"

//###################################################################################################
//#####################################    Main program loop    #####################################
//###################################################################################################
int main(int argc, char *argv[]){
	clock_t begin = clock();	/* Used to measure program cpu run time */
	time_t begin_t = time(NULL);
	srand((unsigned) time(NULL));	/* Intializes random number generator */

	sirsControl *sirsContr = new sirsControl;	/* Controler object for the simulation */

	(*sirsContr).userInput(argc, argv, true);
	(*sirsContr).initialise();
	(*sirsContr).gatherData();
	(*sirsContr).outputDataText();
	(*sirsContr).outputDataMap();

	delete sirsContr;

	printf("\tComplete...");
	fflush(stdout);

	/* Print out program run time */
	clock_t end = clock();
	time_t end_t = time(NULL);
	printf("\tCPU:%-7fs\tClock:%-7fs\n\n", (double)(end - begin) / CLOCKS_PER_SEC, difftime(end_t, begin_t));
	return 1;
}
