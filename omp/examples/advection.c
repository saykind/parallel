#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

double  linear(double x1, double y1, double x2, double y2, double x0);
double  U(double t, double x);
double 	residual(double *u, int N);
FILE   *gpinit(void);
void 	plot(FILE *gpp, double *a, int N);

#define A -1.0
#define B 1.0
#define C 1
#define K 0.2
#define MN 10	// MN*K ?= 1
#define X(i,N) (A+1.0*i*(B-A)/N) 
#define O(i,N) (i>0 ? (i): (i+N))
#define ABS(a) ((a)>0 ? (a) : (-a))
#define rank omp_get_thread_num()

int main(int argc, char *argv[]) {
	int nt, i, j;
	int const N = 10000;
	int const M = MN*N;
	if (argc > 1) {sscanf(argv[1],"%d", &nt);} else {nt = 4;}
	omp_set_num_threads(nt);
	printf("%d\t",nt);

	double *u = (double *) malloc(N*sizeof(double));
	for (i = 0; i < N; i++) {u[i] = U(0, X(i,N));}
	double const dx = (B-A)/N;	
	//FILE *gpp = NULL;
	//gpp = gpinit();
	//plot(gpp, u, N);
	double time = omp_get_wtime();
#pragma omp parallel\
	default(shared)\
	shared(u) private(j, i)
{
	for (j = 0; j < M; j++) {
		#pragma omp for nowait
		for (i = 0; i < N; i++) {
			u[i] = linear(-dx, u[O(i,N)-1], 0, u[i], -K*dx);
		}
		#pragma omp barrier	
		//if (!rank) plot(gpp, u, N);
	}		
}
	time = omp_get_wtime() - time;
	printf("%.3lf\n", time);
	//printf("%.2lfe-3\n", 1e3*dx*residual(u, N));
	//plot(gpp, u, N);
	//pclose(gpp);
	return 0;
}

double linear(double x1, double y1, double x2, double y2, double x0) {
	double k = (y2-y1)/(x2-x1);
	double b =  y1-k*x1;
	return (k*x0+b);
}
double U(double t, double x) {
	if (t == -1) {
		if ((x>=-0.5)&&(x<=0.5)){return 1;} 
		else 	{return 0;}
	}
	if (t == 1) {return sin(M_PI*x);}
	if (t == 0) {return pow(sin(M_PI*x),4);}
	return 0.0;
}
double 	residual(double *u, int N) {
	double sum = 0.0, x = 0.0;
	int i;
	for (i = 0; i < N; i++) 
		if((x = U(0,X(i,N)-u[i])) > 0) sum += x;
		else sum -= x;
	return sum;
}

/* GNUPlot functions */
FILE *gpinit(void) {
	FILE *gpp = popen("gnuplot --persist", "w");
	if (gpp == NULL) {return NULL;}
	fprintf(gpp, "set term x11\n");
	fprintf(gpp, "unset key\n");
	fprintf(gpp, "unset border\n");
	fprintf(gpp, "set style fill solid\n");
	fprintf(gpp, "set yrange [-0.01:1.01]\n");
	fprintf(gpp, "set xrange [-1.01:1.01]\n");
	fprintf(gpp, "set multiplot\n");
	return gpp;
}
void plot(FILE *gpp, double *a, int N) {
	int i;
	fprintf(gpp, "plot '-' w lines\n");
	for (i = 0; i < N; i++) {
		fprintf(gpp, "%.6f %.6f\n", X(i,N), a[i]);
	}
	fprintf(gpp,"e\n");
	fflush(gpp);
}
