#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "integration.h"
#include "plot.h"
#define pi M_PI

#define a 0.
#define b 1.
//double f(double x) {return (x==0.) ? 1 :(sin(x)/x)*(sin(x)/x);}
//double f(double x) {return (1.0/sqrt(pi))*exp(-x*x/4096);}
double f(double x) {return (4.0)/(1+x*x);}


int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	int size, rank, i, n, N;
	double S = 0.;
	double A[2];
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (argc > 2) {sscanf(argv[1],"%d %d", &n, &N);} else {n = 256; N = 4096;}
	//printf("Hi I'm #%d. There is %d of us in total.\n", rank, size);
	double t0 = .0;
	if (!rank) {t0 = MPI_Wtime();}
	if (size == 1) {
		double c[2];
		c[0] = a;
		c[1] = b;
		printf("#%d:\tresult=%.10lf.\n", rank, integrate(f, c, N));
		printf("#%d:\ttime = %f\n", rank, MPI_Wtime()-t0);
		MPI_Finalize();
		return 0;
	}

	MPI_Status st;
	double H = (b-a)/n;
	double I[2];
	I[1] = a;
	if (!rank) {
		for (i = 1; i < size; i++) {
			I[0] = I[1];
			I[1] = I[1] + H;
			MPI_Send((void *) &I, 2, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
		}
	}

	if (rank) {
		while (I[1] < b) {
			MPI_Recv((void *) &I, 2, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &st);
			if (I[0] == b) {continue;}
			A[0] = integrate(f,I,N/n);
			A[1] = (double) rank;
			MPI_Send((void *) &A, 2, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
//			printf("#%d:\tintegrated on [%f,%f] = %.10lf\n", rank, I[0], I[1], A[0]);
		} 
	}
	if (!rank) {
		while (I[1] != b) {
			MPI_Recv((void *) &A, 2, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &st);
			S += A[0];
			I[0] = I[1];
			I[1] = I[1] + H;
			MPI_Send((void *) &I, 2, MPI_DOUBLE,(int) A[1], 1, MPI_COMM_WORLD);
		}
		I[0] = I[1];
		for (i = 1; i < size; i++) {
			if (i == A[1]) continue;
			MPI_Send((void *) &I, 2, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
		}
	}
	if (!rank) printf("#%d:\tresult=%.10lf\n", rank, S);
	if (!rank) printf("#%d:\ttime = %f\n", rank, MPI_Wtime()-t0);
	MPI_Finalize();
	return 0;
}
