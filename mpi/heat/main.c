#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "plot.h"

#define N 512
#define M 4*1024
#define pi M_PI
double f(double x, double y) {return sin(2*pi*y)*sin(2*pi*y)*sin(2*pi*x)*sin(2*pi*x);}

int main(int argc, char *argv[]) {
	int np, rk, rt, i, j, k; double start = .0;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &rk);
	if (!rk) {start = MPI_Wtime();}
	
	double **a = malloc(N*sizeof(double*));
	for (i = 0; i < N; i++)	*(a+i) = malloc(N*sizeof(double));
	double **b = malloc(N*sizeof(double*));
	for (i = 0; i < N; i++)	*(b+i) = malloc(N*sizeof(double));
	double **c;
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++) {
			b[i][j] = 0.;
			a[i][j] = f(1.*i/(N-1),1.*j/(N-1));
		}

	FILE *gpp = gpinit();	if (!rk) plot(gpp, a, N);
	for (k = 0; k < M; k++) {
		for (i = rk*N/np; i < (rk+1)*N/np; i++)
			for (j = 0; j < N; j++)
				b[i][j] = .25*(a[(i-1+N)%N][j]+a[i][(j-1+N)%N]+a[(i+1)%N][j]+a[i][(j+1)%N]);
		for (rt = 0; rt < np; rt++) {
			MPI_Bcast((void *) &b[rt*N/np][0], N, MPI_DOUBLE, rt, MPI_COMM_WORLD);
			MPI_Bcast((void *) &b[(rt+1)*N/np-1][0], N, MPI_DOUBLE, rt, MPI_COMM_WORLD);
		}
		c = b; b = a; a = c;
	}
	for (rt = 0; rt < np; rt++) 
		for (i = rt*N/np; i < (rt+1)*N/np; i++)
			MPI_Bcast((void *) *(a+i), N, MPI_DOUBLE, rt, MPI_COMM_WORLD);
	if (!rk) printf("\ttime = %.0lf sec\n", MPI_Wtime()-start);
	if (!rk) plot(gpp, a, N);

	for (i = 0; i < N; i++) {free(*(a+i)); free(*(b+i));}
	free(a); free(b);
	MPI_Finalize();
	return 0;
}
