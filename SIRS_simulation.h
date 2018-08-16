#ifndef SIRS_SIM
#define SIRS_SIM

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>

//###################################################################################################
//####################################    Class 1 declerations   ####################################
//###################################################################################################

class sirsParameters {
public:
	int graphPrint;	/* [1=yes, 0=no] */
	int numSamples;	/* Number of Monte Carlo sweeps [maximum = 1000] */
	int matrixSize;	/* Edge length of 3D matrix [maximum = 1000] */
	int nSweeps; 	/* Number of simulation sweeps to equilibriate */
	double p0;		/* Probability of becoming dead [0.0 : 1.0] */
	double p1;		/* Probability of becoming infected [0.0 : 1.0] */
	double p2;		/* Probability of becoming recovered [0.0 : 1.0] */
	double p3;		/* Probability of becoming susceptible [0.0 : 1.0] */
	double p4;	/* Probability of becoming immune [0.0 : 1.0] */
	double *p_infected, *p_recovered, *p_susceptible, *p_immune, *p_dead;	/* Pointers for data collection output */
	pthread_mutex_t *mutex1, *mutex2;	/* mutex lock pointers */
	sirsParameters();	/* Parameters constructor */
	~sirsParameters();	/* Parameters destructor */
};

/* Simulation parameters constructor */
sirsParameters::sirsParameters()
{
	graphPrint = 0;	/* [1=yes, 0=no] */
	numSamples = 100;	/* Number of Monte Carlo sweeps [maximum = 1000] */
	matrixSize = 5;	/* Edge length of 3D matrix [maximum = 1000] */
	nSweeps = 25; 	/* Number of simulation sweeps to equilibriate */
	p0 = 0.00000;		/* Probability of becoming dead [0.0 : 1.0] */
	p1 = 0.50000;		/* Probability of becoming infected [0.0 : 1.0] */
	p2 = 0.50000;		/* Probability of becoming recovered [0.0 : 1.0] */
	p3 = 0.50000;		/* Probability of becoming susceptible [0.0 : 1.0] */
	p4 = 0.00000;	/* Probability of becoming immune [0.0 : 1.0] */
}

/* Simulation parameters destructor */
sirsParameters::~sirsParameters()
{
	pthread_mutex_destroy(mutex1);
	pthread_mutex_destroy(mutex2);
	//printf("sirsParameters destructer called\n");
}

//###################################################################################################
//####################################    Class 2 declerations   ####################################
//###################################################################################################
class sirsSimulation {
	private:
		int *graphPrint;	/* [1=yes, 0=no] */
		int *numSamples;	/* Number of Monte Carlo sweeps [maximum = 1000] */
		int *matrixSize;	/* Edge length of 3D matrix [maximum = 1000] */
		int *nSweeps; 	/* Number of simulation sweeps to equilibriate */
		double *p0;		/* Probability of becoming dead [0.0 : 1.0] */
		double *p1;		/* Probability of becoming infected [0.0 : 1.0] */
		double *p2;		/* Probability of becoming recovered [0.0 : 1.0] */
		double *p3;		/* Probability of becoming susceptible [0.0 : 1.0] */
		double *p4;		/* Probability of becoming immune [0.0 : 1.0] */

		double *infected_average, *recovered_average, *susceptible_average, *immune_average, *dead_average;
		double *infected_average_error, *recovered_average_error, *susceptible_average_error, *immune_average_error, *dead_average_error;

		/* Declare arrays for collecting samples */
		double *susceptible, *infected, *healthy, *immune, *dead;
		int *j;	/* Counter for filling the above arrays */

		int *sweep;	/* Size of one Monte Carlo sweep */
		int *nRuns;	/* Total number of runs needed */
		int ***matrixS;	/* Matrix decleration for running the simulation inside */

		double randgen(int min, int max);
		void gnuplot_graph_setup_graph(FILE * pointerPipe);
	public:
		sirsSimulation();	/* Example of function overloading constructor */
		~sirsSimulation();	/* Simulation destructor */
		void setValues(int graphPrintIn,int numSamplesIn,int matrixSizeIn,int nSweepsIn,double p0In,double p1In,double p2In,double p3In,double p4In);
		void fillMatrix();
		int accept_reject_change(int *i,int *j,int *l);
		double percent_of_check_3D(int check);
		void graphSimulation();
		double matrix_sum_av1D(int *Cols, double *matrix);
		double matrix_sum_av_error_1D(int *Cols, double *matrix);
		void caluclateAverages();
		void simRun();
		double returnP0();
		double returnP1();
		double returnP2();
		double returnP3();
		double returnP4();
		double returnInfected();
		double returnInfectedError();
		double returnRecovered();
		double returnRecoveredError();
		double returnSusceptible();
		double returnSusceptibleError();
		double returnImmune();
		double returnImmuneError();
		double returnDead();
		double returnDeadError();
};

/* Simulation constructor */
sirsSimulation::sirsSimulation()
{
	graphPrint = new int(0);
	numSamples = new int(100);
	matrixSize = new int(5);
	nSweeps = new int(25);
	p0 = new double(0.00000);
	p1 = new double(0.50000);
	p2 = new double(0.50000);
	p3 = new double(0.50000);
	p4 = new double(0.00000);

	infected_average = new double(0.00000);
	recovered_average = new double(0.00000);
	susceptible_average = new double(0.00000);
	immune_average = new double(0.00000);
	dead_average = new double(0.00000);
	infected_average_error = new double(0.00000);
	recovered_average_error = new double(0.00000);
	susceptible_average_error = new double(0.00000);
	immune_average_error = new double(0.00000);
	dead_average_error = new double(0.00000);

	j = new int(0);
}

/* Simulation destructor */
sirsSimulation::~sirsSimulation()
{
	for (int i=0 ; i < *matrixSize ; ++i){
		for (int ii=0 ; ii < *matrixSize ; ++ii){
			delete [] matrixS[i][ii];
		}
		delete [] matrixS[i];
	}
	delete [] matrixS;

	delete graphPrint;
	delete numSamples;
	delete matrixSize;
	delete nSweeps;
	delete p0;
	delete p1;
	delete p2;
	delete p3;
	delete p4;

	delete infected_average;
	delete recovered_average;
	delete susceptible_average;
	delete immune_average;
	delete dead_average;
	delete infected_average_error;
	delete recovered_average_error;
	delete susceptible_average_error;
	delete immune_average_error;
	delete dead_average_error;

	delete [] susceptible;
	delete [] infected;
	delete [] healthy;
	delete [] immune;
	delete [] dead;
	delete j;

	delete sweep;
	delete nRuns;
	//printf("sirsSimulation destructer called\n");
}

void sirsSimulation::setValues(int graphPrintIn,int numSamplesIn,int matrixSizeIn,int nSweepsIn,double p0In,double p1In,double p2In,double p3In,double p4In)
{
	*graphPrint = graphPrintIn;
	*numSamples = numSamplesIn;
	*matrixSize = matrixSizeIn;
	*nSweeps = nSweepsIn;
	*p0 = p0In;
	*p1 = p1In;
	*p2 = p2In;
	*p3 = p3In;
	*p4 = p4In;

	susceptible = new double[*numSamples];
	infected = new double[*numSamples];
	healthy = new double[*numSamples];
	immune = new double[*numSamples];
	dead = new double[*numSamples];

	sweep = new int(pow(*matrixSize,3));
	nRuns = new int((*nSweeps * *sweep) + (*numSamples * *sweep));

	/* Initialise matrixes for simulation and data collection */
	matrixS = new int**[*matrixSize];
	for (int i=0 ; i < *matrixSize ; ++i){
		matrixS[i] = new int*[*matrixSize];
		for (int ii=0 ; ii < *matrixSize ; ++ii){
			matrixS[i][ii] = new int[*matrixSize];
			for (int iii=0; iii < *matrixSize; ++iii){
				matrixS[i][ii][iii] = 0;
			}
		}
	}
	return;
}

/* Generate random double between min & max */
double sirsSimulation::randgen(int min, int max)
{
	double i, ii, iii;
	i = rand();
	ii = RAND_MAX;
	iii = (i/ii) * (max - min);
	return (iii + min);
}

void sirsSimulation::gnuplot_graph_setup_graph(FILE * pointerPipe)
{
	if (pointerPipe == NULL){
		printf("ERROR:\tOpening pipe to gnuplot failed. Check if it exists.\n ");
		exit(0);
	}
	fprintf(pointerPipe, "set terminal png font ',10.0' size 640,480 \n");
	fprintf(pointerPipe, "set autoscale yfix \n");
	fprintf(pointerPipe, "set autoscale xfix \n");
	fprintf(pointerPipe, "set yrange [0.0:1.0] \n");
	fprintf(pointerPipe, "set xlabel 'Monte Carlo sweeps' \n");
	fprintf(pointerPipe, "set ylabel '%%' \n");
 	fflush(pointerPipe);	/* Flushing the buffer of fprintf */
 	return;
}

void sirsSimulation::fillMatrix()
{
	for (int i=0 ; i < *matrixSize ; ++i){
		for(int j=0 ; j < *matrixSize ; ++j){
			for (int l=0 ; l < *matrixSize ; ++l){
				double value = randgen(0,100);		/* Randomly choose a value */
				double value1 = randgen(0,50);		/* Randomly choose cut off point bellow halfway */
				double value2 = randgen(50,100);	/* Randomly choose cut off point above halfway */
				if (0.0<=value && value<value1){
					matrixS[i][j][l] = 1;	/* 1 = Susceptible */
					continue;
				} else if (value1<=value && value<value2){
					matrixS[i][j][l] = 2;	/* 2 = Infected */
					continue;
				} else if (value2<=value && value<=100.0){
					matrixS[i][j][l] = 3;	/* 3 = Recovered */
					continue;
				}
			}
		}
	}
	return;
}

/* Parameters for changing the lattice point */
int sirsSimulation::accept_reject_change(int *i,int *j,int *l)
{
	/* Use these to check cells on each side in the x,y,z directions */
	int iup = (*i + 1) % *matrixSize;
	int idown = (*i - 1 + *matrixSize) % *matrixSize;
	int jup = (*j + 1) % *matrixSize;
	int jdown = (*j - 1 + *matrixSize) % *matrixSize;
	int lup = (*l + 1) % *matrixSize;
	int ldown = (*l - 1 + *matrixSize) % *matrixSize;

	double p = randgen(0,1);	/* Generate probability number */
	/* [0=dead, 1=susceptible, 2=infected, 3=recovered, 4=immune] */
	if (matrixS[*i][*j][*l] != 0 && p <= *p0){
		matrixS[*i][*j][*l] = 0;	/* Becomes dead */
		return 1;	/* Matrix is changed */
	} else if (matrixS[*i][*j][*l] == 0 && p <= *p2){
		matrixS[*i][*j][*l] = 3;	/* Becomes healthy (rebirth) */
		return 1;	/* Matrix is changed */
	} else if (matrixS[*i][*j][*l] == 1){
		if (p <= *p1 && (matrixS[iup][*j][*l] == 2 || matrixS[idown][*j][*l] == 2 || matrixS[*i][jup][*l] == 2 || matrixS[*i][jdown][*l] == 2 || matrixS[*i][*j][lup] == 2 || matrixS[*i][*j][ldown] == 2)){
				matrixS[*i][*j][*l] = 2;	/* Becomes infected */
				return 1;	/* Matrix is changed */
			} else if (p <= *p4 && (matrixS[iup][*j][*l] == 4 || matrixS[idown][*j][*l] == 4 || matrixS[*i][jup][*l] == 4 || matrixS[*i][jdown][*l] == 4 || matrixS[*i][*j][lup] == 4 || matrixS[*i][*j][ldown] == 4)){
				matrixS[*i][*j][*l] = 4;	/* Becomes immune */
				return 1;	/* Matrix is changed */
			}
	} else 	if (matrixS[*i][*j][*l] == 2){
		if (p <= (*p4 * *p2)){
				matrixS[*i][*j][*l] = 4;	/* Becomes immune */
				return 1;	/* Matrix is changed */
			} else if (p <= *p2){
				matrixS[*i][*j][*l] = 3;	/* Becomes healthy */
				return 1;	/* Matrix is changed */
			}
	} else if (matrixS[*i][*j][*l] == 3){
		if (p <= *p3){
				matrixS[*i][*j][*l] = 1;	/* Becomes susceptible */
				return 1;	/* Matrix is changed */
			}
	}
	return 0; 	/* Matrix is NOT changed */
}

double sirsSimulation::percent_of_check_3D(int check)
{
	double ii = 0.0;
	for (int i=0 ; i < *matrixSize ; ++i){
		for(int j=0; j < *matrixSize ; ++j){
			for(int l=0; l < *matrixSize ; ++l){
				ii += ( matrixS[i][j][l] == check) ? 1 : 0;	/* (test) ? <return this if true> : <return this if false> */
			}
		}
	}
	return (ii / pow(*matrixSize,3.0));
}

void sirsSimulation::graphSimulation()
{
	if (*graphPrint){
		/* Create output folder */
		char graphsFolderName[50];
		sprintf(graphsFolderName,"mkdir -p Graphs/p2_%.2f",*p2);
		system(graphsFolderName);

		/* Output to file */
		char graphsPath[50];
		sprintf(graphsPath,"Graphs/p2_%.2f/_output.txt",*p2);
		FILE *filePointer = fopen(graphsPath, "w+");	/*Creating the data output file */
		for (int i=0 ; i < *numSamples ; ++i){
			fprintf(filePointer,"%6i \t%6f \t%6f \t%6f \t%6f \t%6f\n",*nSweeps+i,infected[i],healthy[i],susceptible[i],immune[i],dead[i]);
		}
		fflush(filePointer);
		fclose(filePointer);

		/* Run gnuPlot commands */
		FILE *gnuPipe1 = popen("gnuplot", "w");	/* Opening the gnuplot program */
		gnuplot_graph_setup_graph(gnuPipe1);
		fprintf(gnuPipe1, "set output 'Graphs/p2_%.2f/SIRS_Graph_%.3f_%.3f_%.3f_%.3f_%.3f.png' \n",*p2,*p0,*p1,*p2,*p3,*p4);
		fprintf(gnuPipe1, "plot 'Graphs/p2_%.2f/_output.txt' using 1:2 with lines title 'infected',",*p2);
		fprintf(gnuPipe1, " 'Graphs/p2_%.2f/_output.txt' using 1:3 with lines title 'healthy',",*p2);
		fprintf(gnuPipe1, " 'Graphs/p2_%.2f/_output.txt' using 1:4 with lines title 'susceptible',",*p2);
		fprintf(gnuPipe1, " 'Graphs/p2_%.2f/_output.txt' using 1:5 with lines title 'immune',",*p2);
		fprintf(gnuPipe1, " 'Graphs/p2_%.2f/_output.txt' using 1:6 with lines title 'dead' \n",*p2);
		fflush(gnuPipe1);
		pclose(gnuPipe1);
	}
	return;
}

/* Calculate average for 1D matrix */
double sirsSimulation::matrix_sum_av1D(int *Cols, double *matrix)
{
	double M=0.0;
	for (int i=0 ; i < *Cols ; ++i){
		M += *(matrix + i);
	}
	double av_M = M / *Cols;
	return av_M;
}

/* Calculate error on the average for 1D matrix */
double sirsSimulation::matrix_sum_av_error_1D(int *Cols, double *matrix)
{
	double k=0.0;
	for (int i=0 ; i < *Cols ; ++i){
		k += *(matrix + i);
	}
	double mean = k / *Cols;
	double variance=0.0;
	for (int i=0 ; i < *Cols ; ++i){
		variance += pow((*(matrix + i) - mean),2.0)/(*Cols -1);
	}
	double sigma= sqrt(variance);
	double error= sigma / sqrt(*Cols);
	return error;
}

void sirsSimulation::caluclateAverages()
{
	/* Computate average values and their errors */
	*infected_average = matrix_sum_av1D(numSamples, &infected[0]);
	*recovered_average = matrix_sum_av1D(numSamples, &healthy[0]);
	*susceptible_average = matrix_sum_av1D(numSamples, &susceptible[0]);
	*immune_average = matrix_sum_av1D(numSamples, &immune[0]);
	*dead_average = matrix_sum_av1D(numSamples, &dead[0]);

	*infected_average_error = matrix_sum_av_error_1D(numSamples, &infected[0]);
	*recovered_average_error = matrix_sum_av_error_1D(numSamples, &healthy[0]);
	*susceptible_average_error = matrix_sum_av_error_1D(numSamples, &susceptible[0]);
	*immune_average_error = matrix_sum_av_error_1D(numSamples, &immune[0]);
	*dead_average_error = matrix_sum_av_error_1D(numSamples, &dead[0]);
	return;
}

void sirsSimulation::simRun()
{
	int iPoint,jPoint,lPoint;	/* i, j, and l matrix points where check occures */

	for (int i=0 ; i < *nRuns ; ++i){
		/* Find coordinates to check */
		iPoint = (int) randgen(0,*matrixSize -1);
		jPoint = (int) randgen(0,*matrixSize -1);
		lPoint = (int) randgen(0,*matrixSize -1);

		/* Check matrix point and see if it needs to be updated to a new state */
		accept_reject_change(&iPoint,&jPoint,&lPoint);

		if(i>=((*nSweeps * *sweep) - 1) && ((i)% *sweep)==0){
			dead[*j] = percent_of_check_3D(0);
			susceptible[*j] = percent_of_check_3D(1);
			infected[*j] = percent_of_check_3D(2);
			healthy[*j] = percent_of_check_3D(3);
			immune[*j] = percent_of_check_3D(4);
			++(*j);
		}
	}
	return;
}
double sirsSimulation::returnP0()
{
	return *p0;
}
double sirsSimulation::returnP1()
{
	return *p1;
}
double sirsSimulation::returnP2()
{
	return *p2;
}
double sirsSimulation::returnP3()
{
	return *p3;
}
double sirsSimulation::returnP4()
{
	return *p4;
}
double sirsSimulation::returnInfected()
{
	return *infected_average;
}
double sirsSimulation::returnInfectedError()
{
	return *infected_average_error;
}
double sirsSimulation::returnRecovered()
{
	return *recovered_average;
}
double sirsSimulation::returnRecoveredError()
{
	return *recovered_average_error;
}
double sirsSimulation::returnSusceptible()
{
	return *susceptible_average;
}
double sirsSimulation::returnSusceptibleError()
{
	return *susceptible_average_error;
}
double sirsSimulation::returnImmune()
{
	return *immune_average;
}
double sirsSimulation::returnImmuneError()
{
	return *immune_average_error;
}
double sirsSimulation::returnDead()
{
	return *dead_average;
}
double sirsSimulation::returnDeadError()
{
	return *dead_average_error;
}

#endif
