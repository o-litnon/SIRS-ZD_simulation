Requires gnuplot to be installed, and UNIX path linked as 'gnuplot'

To compile:		g++ SIRS_main.cpp -lm -lpthread
To run:			./a.out <options>
				-g 		print out graphs from each sim
				-n <int> 	give an integer number of threads
				-m <int> 	give an integer number for sim size
				-r <int> 	give an integer number map resolution
				-0 <float> 	give a float for p0 <0.0 : 1.0>
				-2 <float> 	give a float for p2 <0.0 : 1.0>
				-4 <float> 	give a float for p4 <0.0 : 1.0>


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

#!/bin/bash
declare -a p0values=(0.0)
declare -a p1values=(0.01 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99)
declare -a p4values=(0.0)
g++ SIRS_main.cpp -lm -lpthread; for i in ${p0values[@]}; do for ii in ${p1values[@]}; do for iii in ${p4values[@]}; do ./a.out -g -n 10 -m 10 -0 $i -2 $ii -4 $iii; done; done; done;
