#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi/mpi.h>

double  linear(double x1, double y1, double x2, double y2, double x0);
double  U(double t, double x);
double 	residual(double *u, int N);
FILE   *gpinit(void);
void 	plot(FILE *gpp, double *a, int N);

#define A -1.0
#define B 1.0
#define C 1
#define K 0.2
#define MN 4	// MN*K ?= 1
#define X(i,N) (A+1.0*i*(B-A)/N) 
#define O(i,N) (i>0 ? (i): (i+N))
#define Q(r,s) (r<s ? (r) : (0))
#define W(r,s) (r>=0 ? (r):(r+size))
#define ABS(a) ((a)>0 ? (a) : (-a))

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	int size, rank;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	printf("%d\t", size);
	int N;
	if (argc > 1) {sscanf(argv[1],"%d", &N);} else {N = 10000;}
	int M = MN*N-1;

	double *u = (double *) malloc(N*sizeof(double));
	double dx = (B-A)/N;	
	int i, j;
	//FILE *gpp = NULL;
	//if (!rank) gpp = gpinit();
	for (i = 0; i < N; i++) {u[i] = U(0, X(i,N));}
	//if (!rank) plot(gpp, u, N);
	double time_start = 0, time_end = 0;
        if (!rank) {time_start = MPI_Wtime();}
	/* main body */
	int n = N/size;
	int from = n*rank;
	int to = n*(rank+1);
	if (rank == (size-1)) {to = N;}
	for (j = 0; j < M; j++) {
		for (i = from; i < to; i++) {
			u[i] = linear(-dx, u[O(i,N)-1], 0, u[i], -K*dx);
			//printf("[%d,%d]\n", j, i);
		}	
		MPI_Status st;
		if (size > 1) MPI_Send((void *) &u[to-1], 1, MPI_DOUBLE, Q(rank+1,size), 1, MPI_COMM_WORLD);
		if (size > 1) MPI_Recv((void *) &u[O(from,N)-1], 1, MPI_DOUBLE, W(rank-1,size), 1, MPI_COMM_WORLD, &st);
	}		
	if (rank && size > 1)
		MPI_Send((void *) (u+from), to-from, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	if (!rank && size > 1) {
		for (i = 1; i < size-1; i++)
			MPI_Recv((void *) (u+i*n), n, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv((void *) (u+(size-1)*n), n + N%size, MPI_DOUBLE, size-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	/* main body */
	if (!rank) {time_end = MPI_Wtime(); printf("%.3lf\n", (time_end-time_start));}
	//if (!rank) plot(gpp, u, N);
	//if (!rank) pclose(gpp);
	//if (!rank) printf("%.2lfe-3\t", 1e3*dx*residual(u, N));
	MPI_Finalize();
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
	FILE *gpp = popen("gnuplot", "w");
	if (gpp == NULL) {return NULL;}
	fprintf(gpp, "set term x11\n");
	fprintf(gpp, "unset key\n");
	fprintf(gpp, "unset border\n");
	fprintf(gpp, "set style fill solid\n");
	fprintf(gpp, "set yrange [-0.01:1.01]\n");
	fprintf(gpp, "set xrange [-1.01:1.01]\n");
	return gpp;
}
void plot(FILE *gpp, double *a, int N) {
	int i;
	fprintf(gpp, "plot '-' w boxes\n");
	for (i = 0; i < N; i++) {
		fprintf(gpp, "%.6f %.6f\n", X(i,N), a[i]);
	}
	fprintf(gpp,"e\n");
	fflush(gpp);
}
