#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "plot.h"

#define N 1024
#define M 4*1024
#define rank omp_get_thread_num()
#define pi M_PI
double f(double x, double y) {return sin(2*pi*y)*sin(2*pi*y)*sin(2*pi*x)*sin(2*pi*x);}

int main(int argc, char *argv[]) {
	int np, i, j, k;
	if (argc > 1) {sscanf(argv[1],"%d", &np);} else {np = 1;} 
	omp_set_num_threads(np); printf("%d\t", np);
	
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
	FILE *gpp = gpinit();	plot(gpp, a, N);
	double t0 = omp_get_wtime();
#pragma omp parallel default(shared) private(k,i,j) 
{
	for (k = 0; k < M; k++) {
		#pragma omp for nowait
		for (i = 1; i < N-1; i++)
			for (j = 1; j < N-1; j++)
				b[i][j] = .25*(a[i-1][j]+a[i][j-1]+a[i+1][j]+a[i][j+1]);
		#pragma omp barrier
		if (!rank) {c = b; b = a; a = c;}
		if (!rank) plot(gpp, a, N);
	}
}
	printf("time =\t %.0lf sec\n", omp_get_wtime()-t0);

	for (i = 0; i < N; i++) {free(*(a+i)); free(*(b+i));}
	free(a); free(b);
	return 0;
}
