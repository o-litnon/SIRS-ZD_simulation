#ifndef SIRS_SIM_CONTROL
#define SIRS_SIM_CONTROL

#include "SIRS_simulation.h"

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <math.h>

//###################################################################################################
//####################################    Class 1 declerations   ####################################
//###################################################################################################
class sirsControl {
private:
	pthread_mutex_t mutex1, mutex2;		/* mutex1 lock for reading parameters & mutex2 lock for printing graphs */
	int *possibleThreads; /* Number of accessible threads */
	double ***m_infected, ***m_recovered, ***m_susceptible, ***m_immune, ***m_dead;	/* Pointers for data collection mattrix */
	pthread_t *tid;	/* Declate a thread number for each thread */
	int *pSteps;	/* Number of steps between p1=0.0 and p1=1.0 */
	sirsParameters *simParams;	/*Declare a sirsParameters object to pass into each thread */
	static void* sirs(void *sirsParam);
	void gnuplot_graph_setup_map(FILE * pointerPipe);
public:
	sirsControl();
	~sirsControl();
	void userInput(int argc, char *argv[], int print);
	void userValues(int graphPrintIn,int possiblethreadsIn,int matrixSizeIn,int pStepsIn,double p0In,double p2In,double p4In);
	void initialise();
	void gatherData();
	void outputDataText();
	void outputDataMap();
};

/* Simulation control constructor */
sirsControl::sirsControl()
{
	pthread_mutex_init(&mutex1, NULL);
	pthread_mutex_init(&mutex2, NULL);
	possibleThreads = new int(1);
	pSteps = new int(30);
	simParams = new sirsParameters;
}
/* Simulation control destructor */
sirsControl::~sirsControl()
{
	for (int i=0 ; i<*pSteps ; ++i){
		for (int ii=0 ; ii<*pSteps ; ++ii){
				delete [] m_infected[i][ii];
				delete [] m_recovered[i][ii];
				delete [] m_susceptible[i][ii];
				delete [] m_immune[i][ii];
				delete [] m_dead[i][ii];
		}
		delete [] m_infected[i];
		delete [] m_recovered[i];
		delete [] m_susceptible[i];
		delete [] m_immune[i];
		delete [] m_dead[i];
	}
	delete [] m_infected;
	delete [] m_recovered;
	delete [] m_susceptible;
	delete [] m_immune;
	delete [] m_dead;

	pthread_mutex_destroy(&mutex1);
	pthread_mutex_destroy(&mutex2);
	delete possibleThreads;
	delete [] tid;
	delete pSteps;
	delete simParams;
	//printf("sirsControl destructer called\n");
}

/* Simulations thread build */
void* sirsControl::sirs(void *sirsParam)
{
	sirsParameters *sirsParams = (sirsParameters*)sirsParam;
	sirsSimulation *sirsSim = new sirsSimulation;

	(*sirsSim).setValues((*sirsParams).graphPrint,(*sirsParams).numSamples,(*sirsParams).matrixSize,(*sirsParams).nSweeps,(*sirsParams).p0,(*sirsParams).p1,(*sirsParams).p2,(*sirsParams).p3,(*sirsParams).p4); /* Imported object from SIRS_simulation.h */
	
	double *m_infected = (*sirsParams).p_infected;
	double *m_recovered = (*sirsParams).p_recovered;
	double *m_susceptible = (*sirsParams).p_susceptible;
	double *m_immune = (*sirsParams).p_immune;
	double *m_dead = (*sirsParams).p_dead;
	pthread_mutex_t *mutex1 = (*sirsParams).mutex1;
	pthread_mutex_t *mutex2 = (*sirsParams).mutex2;
	pthread_mutex_unlock(mutex1);

	(*sirsSim).fillMatrix();
	(*sirsSim).simRun();

	pthread_mutex_lock(mutex2);
	(*sirsSim).graphSimulation();
	pthread_mutex_unlock(mutex2);

	(*sirsSim).caluclateAverages();

	/* Return the average values and errors to the global array */
	*(m_infected + 0) =  (*sirsSim).returnP0();
	*(m_recovered + 0) = (*sirsSim).returnP0();
	*(m_susceptible + 0) = (*sirsSim).returnP0();
	*(m_immune + 0) = (*sirsSim).returnP0();
	*(m_dead + 0) = (*sirsSim).returnP0();

	*(m_infected + 1) =  (*sirsSim).returnP1();
	*(m_recovered + 1) = (*sirsSim).returnP1();
	*(m_susceptible + 1) = (*sirsSim).returnP1();
	*(m_immune + 1) = (*sirsSim).returnP1();
	*(m_dead + 1) = (*sirsSim).returnP1();

	*(m_infected + 2) = (*sirsSim).returnP2();
	*(m_recovered + 2) = (*sirsSim).returnP2();
	*(m_susceptible + 2) = (*sirsSim).returnP2();
	*(m_immune + 2) = (*sirsSim).returnP2();
	*(m_dead + 2) = (*sirsSim).returnP2();

	*(m_infected + 3) = (*sirsSim).returnP3();
	*(m_recovered + 3) = (*sirsSim).returnP3();
	*(m_susceptible + 3) = (*sirsSim).returnP3();
	*(m_immune + 3) = (*sirsSim).returnP3();
	*(m_dead + 3) = (*sirsSim).returnP3();

	*(m_infected + 4) = (*sirsSim).returnP4();
	*(m_recovered + 4) = (*sirsSim).returnP4();
	*(m_susceptible + 4) = (*sirsSim).returnP4();
	*(m_immune + 4) = (*sirsSim).returnP4();
	*(m_dead + 4) = (*sirsSim).returnP4();

	*(m_infected + 5) = (*sirsSim).returnInfected();
	*(m_recovered + 5) = (*sirsSim).returnRecovered();
	*(m_susceptible + 5) = (*sirsSim).returnSusceptible();
	*(m_immune + 5) = (*sirsSim).returnImmune();
	*(m_dead + 5) = (*sirsSim).returnDead();

	*(m_infected + 6) = (*sirsSim).returnInfectedError();
	*(m_recovered + 6) = (*sirsSim).returnRecoveredError();
	*(m_susceptible + 6) = (*sirsSim).returnSusceptibleError();
	*(m_immune + 6) = (*sirsSim).returnImmuneError();
	*(m_dead + 6) = (*sirsSim).returnDeadError();

	delete sirsSim;
	pthread_exit(NULL);
}

/* Initialising GNUPLOT settings */
void sirsControl::gnuplot_graph_setup_map(FILE * pointerPipe)
{
	if (pointerPipe == NULL){
		printf("ERROR:\tOpening pipe to gnuplot failed. Check if it exists.\n ");
		exit(0);
	}
	fprintf(pointerPipe, "set terminal png font ',10.0' size 1280,960 \n");
	fprintf(pointerPipe, "set autoscale yfix \n");
	fprintf(pointerPipe, "set autoscale xfix \n");
	fprintf(pointerPipe, "set cbrange [0.0:1.0]\n");
	fprintf(pointerPipe, "set xlabel 'p1(Infection)' \n");
	fprintf(pointerPipe, "set ylabel 'p3(Susceptible)' \n");
	fprintf(pointerPipe, "set view map \n");
	fprintf(pointerPipe, "set pm3d interpolate 1,1\n");
 	fflush(pointerPipe);	/* Flushing the buffer of fprintf */
 	return;
}

void sirsControl::userInput(int argc, char *argv[], int print)
{
	int c;
	while ((c = getopt(argc, argv, "g n:m:r:0:2:4:")) != -1)
		switch (c){
			case 'g':
				(*simParams).graphPrint = 1;
				break;
			case 'n':
				*possibleThreads = atoi(optarg);
				break;
			case 'm':
				(*simParams).matrixSize = atoi(optarg);
				break;
			case 'r':
				*pSteps = atoi(optarg);
				break;
			case '0':
				(*simParams).p0 = atof(optarg);
				break;
			case '2':
				(*simParams).p2 = atof(optarg);
				break;
			case '4':
				(*simParams).p4 = atof(optarg);
				break;
			case '?':
				printf("illegal option...\n\n");
				printf("%s <options>\n", argv[0]);
				printf("\t-g \t\tprint out graphs from each sim\n");
				printf("\t-n <int> \tgive an integer number of threads\n");
				printf("\t-m <int> \tgive an integer number for sim size\n");
				printf("\t-r <int> \tgive an integer number map resolution\n");
				printf("\t-0 <float> \tgive a float for p0 <0.0 : 1.0>\n");
				printf("\t-2 <float> \tgive a float for p2 <0.0 : 1.0>\n");
				printf("\t-4 <float> \tgive a float for p4 <0.0 : 1.0>\n");
				printf("\n");
				exit(0);
				break;	/* Good practice to include 'break' in a each switch case */
			default:
				break;
    }

    if (print){
	    printf("\nSIRS simulation parameters...\tResolution:%ix%i\tN_Threads:%i\tGaphing:%i\tSize:%i^3\tp0:%g\tp2:%g\tp4:%g",*pSteps,*pSteps, *possibleThreads, (*simParams).graphPrint,(*simParams).matrixSize,(*simParams).p0,(*simParams).p2,(*simParams).p4);
		printf("\n");
		fflush(stdout);
	}
}

void sirsControl::userValues(int graphPrintIn,int possiblethreadsIn,int matrixSizeIn,int pStepsIn,double p0In,double p2In,double p4In)
{
	(*simParams).graphPrint = graphPrintIn;
	*possibleThreads = possiblethreadsIn;
	(*simParams).matrixSize = matrixSizeIn;
	*pSteps = pStepsIn;
	(*simParams).p0 = p0In;
	(*simParams).p2 = p2In;
	(*simParams).p4 = p4In;
}

void sirsControl::initialise()
{
	/* Initialise matrixes for simulation and data collection */
	m_infected = new double**[*pSteps];
	m_recovered = new double**[*pSteps];
	m_susceptible = new double**[*pSteps];
	m_immune = new double**[*pSteps];
	m_dead = new double**[*pSteps];
	for (int i=0 ; i<*pSteps ; ++i){
		m_infected[i] = new double*[*pSteps];
		m_recovered[i] = new double*[*pSteps];
		m_susceptible[i] = new double*[*pSteps];
		m_immune[i] = new double*[*pSteps];
		m_dead[i] = new double*[*pSteps];
		for (int ii=0 ; ii<*pSteps ; ++ii){
			m_infected[i][ii] = new double[7];
			m_recovered[i][ii] = new double[7];
			m_susceptible[i][ii] = new double[7];
			m_immune[i][ii] = new double[7];
			m_dead[i][ii] = new double[7];
			for (int iii=0; iii<7; ++iii){
				m_infected[i][ii][iii] = 0.0;
				m_recovered[i][ii][iii] = 0.0;
				m_susceptible[i][ii][iii] = 0.0;
				m_immune[i][ii][iii] = 0.0;
				m_dead[i][ii][iii] = 0.0;
			}
		}
	}

	tid = new pthread_t[*possibleThreads];	/* Declate a thread number for each thread */
	return;
}

void sirsControl::gatherData()
{
	/* Loop over each data point */
	for (int i=0 ; i<*pSteps ; ++i){
		for (int ii=0 ; ii<*pSteps ; ++ii){

			static int thread_counter = 0;

			/* Create parameter settings for each thread */
			int threadID = thread_counter % *possibleThreads;

			if (thread_counter >= *possibleThreads){
				pthread_join(tid[threadID],NULL); /* Wait for last thread of the same ID to finish before changing the sim parameters*/
			}

			pthread_mutex_lock(&mutex1);	/* Lock the parameters so that they are not lost between thread calls */
			(*simParams).p1 = (double) i/(*pSteps-1);
			(*simParams).p3 = (double) ii/(*pSteps-1);
			(*simParams).p_infected = &m_infected[i][ii][0];
			(*simParams).p_recovered = &m_recovered[i][ii][0];
			(*simParams).p_susceptible = &m_susceptible[i][ii][0];
			(*simParams).p_immune = &m_immune[i][ii][0];
			(*simParams).p_dead = &m_dead[i][ii][0];
			(*simParams).mutex1 = &mutex1;
			(*simParams).mutex2 = &mutex2;

			/* Run the simulation and output values to global matrix */
			pthread_create(&tid[threadID], NULL, sirs, simParams);
			++thread_counter;

			/* Print progress */
			double percent = (((double)i * *pSteps) + (double)ii) / (pow(*pSteps,2)-1);
			printf("\t%07.3f%%\r", 100.0*percent);
			fflush(stdout);
		}
	}
	/* Make sure all threads are finished */
	for (int i=0 ; i < *possibleThreads ; ++i){
		pthread_join(tid[i],NULL);
	}
}

void sirsControl::outputDataText()
{
	/* Create output folder */
	system("mkdir -p Graphs");

	/* Output map data to file */
	FILE *filePointer = fopen("Graphs/_output.txt", "w+");	/* Creating the data output file */
	for (int i=0 ; i<*pSteps ; ++i){
		for (int iii=0 ; iii<*pSteps ; ++iii){
			/* [p1][p2][p3][p4][infected][error][recovered][error][susceptible][error][immune][error][dead][error] */
			fprintf(filePointer,"%6f \t%6f \t%6f \t%6f \t%6f \t%6f \t%6f \t%6f \t%6f \t%6f \t%6f \t%6f \t%6f \t%6f \t%6f\n",m_infected[i][iii][0],m_infected[i][iii][1], m_infected[i][iii][2], m_infected[i][iii][3], m_infected[i][iii][4], m_infected[i][iii][5], m_infected[i][iii][6], m_recovered[i][iii][5], m_recovered[i][iii][6], m_susceptible[i][iii][5], m_susceptible[i][iii][6], m_immune[i][iii][5], m_immune[i][iii][6], m_dead[i][iii][5], m_dead[i][iii][6]);
		}
		fprintf(filePointer,"\n");
	}
	fflush(filePointer);
	fclose(filePointer);
}

void sirsControl::outputDataMap()
{
	FILE *gnuPipe = popen("gnuplot", "w"); /* Opening the gnuplot program */
	gnuplot_graph_setup_map(gnuPipe);

	fprintf(gnuPipe, "set output 'Graphs/SIRS_Map_%.3f_%.3f_%.3f.png' \n", m_infected[0][0][0], m_infected[0][0][2], m_infected[0][0][4]);
	fprintf(gnuPipe, "set multiplot layout 5,2\n");

	fprintf(gnuPipe, "set title '(Infected)%%  [p0=%.3f,p2=%.3f,p4=%.3f]' \n",m_infected[0][0][0],m_infected[0][0][2], m_infected[0][0][4]);
	fprintf(gnuPipe, "splot 'Graphs/_output.txt' using 2:4:6 with pm3d notitle\n");
	fprintf(gnuPipe, "unset cbrange\n");
	fprintf(gnuPipe, "set title '(Infected error)%%  [p0=%.3f,p2=%.3f,p4=%.3f]' \n",m_infected[0][0][0],m_infected[0][0][2], m_infected[0][0][4]);
	fprintf(gnuPipe, "splot 'Graphs/_output.txt' using 2:4:7 with pm3d notitle\n");

	fprintf(gnuPipe, "set cbrange [0.0:1.0]\n");
	fprintf(gnuPipe, "set title '(Recovered)%%  [p0=%.3f,p2=%.3f,p4=%.3f]' \n",m_infected[0][0][0],m_infected[0][0][2], m_infected[0][0][4]);
	fprintf(gnuPipe, "splot 'Graphs/_output.txt' using 2:4:8 with pm3d notitle\n");
	fprintf(gnuPipe, "unset cbrange\n");
	fprintf(gnuPipe, "set title '(Recovered error)%%  [p0=%.3f,p2=%.3f,p4=%.3f]' \n",m_infected[0][0][0],m_infected[0][0][2], m_infected[0][0][4]);
	fprintf(gnuPipe, "splot 'Graphs/_output.txt' using 2:4:9 with pm3d notitle\n");

	fprintf(gnuPipe, "set cbrange [0.0:1.0]\n");
	fprintf(gnuPipe, "set title '(Susceptible)%%  [p0=%.3f,p2=%.3f,p4=%.3f]' \n",m_infected[0][0][0],m_infected[0][0][2], m_infected[0][0][4]);
	fprintf(gnuPipe, "splot 'Graphs/_output.txt' using 2:4:10 with pm3d notitle\n");
	fprintf(gnuPipe, "unset cbrange\n");
	fprintf(gnuPipe, "set title '(Susceptible error)%%  [p0=%.3f,p2=%.3f,p4=%.3f]' \n",m_infected[0][0][0],m_infected[0][0][2], m_infected[0][0][4]);
	fprintf(gnuPipe, "splot 'Graphs/_output.txt' using 2:4:11 with pm3d notitle\n");

	fprintf(gnuPipe, "set cbrange [0.0:1.0]\n");
	fprintf(gnuPipe, "set title '(Immune)%%  [p0=%.3f,p2=%.3f,p4=%.3f]' \n",m_infected[0][0][0],m_infected[0][0][2], m_infected[0][0][4]);
	fprintf(gnuPipe, "splot 'Graphs/_output.txt' using 2:4:12 with pm3d notitle\n");
	fprintf(gnuPipe, "unset cbrange\n");
	fprintf(gnuPipe, "set title '(Immune error)%%  [p0=%.3f,p2=%.3f,p4=%.3f]' \n",m_infected[0][0][0],m_infected[0][0][2], m_infected[0][0][4]);
	fprintf(gnuPipe, "splot 'Graphs/_output.txt' using 2:4:13 with pm3d notitle\n");

	fprintf(gnuPipe, "set cbrange [0.0:1.0]\n");
	fprintf(gnuPipe, "set title '(Dead)%%  [p0=%.3f,p2=%.3f,p4=%.3f]' \n",m_infected[0][0][0],m_infected[0][0][2], m_infected[0][0][4]);
	fprintf(gnuPipe, "splot 'Graphs/_output.txt' using 2:4:14 with pm3d notitle\n");
	fprintf(gnuPipe, "unset cbrange\n");
	fprintf(gnuPipe, "set title '(Dead error)%%  [p0=%.3f,p2=%.3f,p4=%.3f]' \n",m_infected[0][0][0],m_infected[0][0][2], m_infected[0][0][4]);
	fprintf(gnuPipe, "splot 'Graphs/_output.txt' using 2:4:15 with pm3d notitle\n");

	fprintf(gnuPipe, "unset multiplot\n");
	fflush(gnuPipe);
	pclose(gnuPipe);
}

#endif
