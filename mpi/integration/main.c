#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "integration.h"
#include "plot.h"
#define pi M_PI

#define a 0.0
#define b 1.0
double f(double x) {return (1.0/sqrt(pi))*exp(-x*x);}
double g(double x) {return (1.0/pi)/(1+x*x);}


int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	int size, rank, i, N;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	printf("rank = %d\tsize = %d\n", rank, size);
	double t0 = .0;
	if (!rank) {t0 = MPI_Wtime();}
	

	if (argc > 1) {sscanf(argv[1],"%d", &N);} else {N = 2048;}
	double *F = (double *) malloc(N*sizeof(double));
	double *X = (double *) malloc(N*sizeof(double));

	double h = (b-a)/N;
	for (i = 0; i < N; i++) {
		X[i] = a+h*i;
		F[i] = integrate(g, a, X[i], N);
	}
//	FILE *gpp = gpinit();
//	plot(gpp, X, F, N);
//	usleep(3e6);
//	pclose(gpp);
	printf("rank = %d\tI = %lf\n",rank, F[N-1]);
	free(F);
	free(X);
	if (!rank) printf("rank = %d\ttime = %f\n", rank, MPI_Wtime()-t0);
	MPI_Finalize();
	return 0;
}
//	printf("Pi = %.16f\t(%.16f)\n", 4*I*I, pi);
//	printf("Rectangle:\t I = %.8f\t	error = %.8f\n", rectangle(x,f,N), M2*(b-a)*h*h/6);
//	printf("Trapezoid:\t I = %.8f\t	error = %.8f\n", trapezoid(x,f,N), M2*(b-a)*h*h/24);
//	printf("Simpson:\t I = %.8f\t	error = %.8f\n\n", simpson(x,f,N), M4*(b-a)*h*h*h*h/180);

