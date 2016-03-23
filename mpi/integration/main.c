#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#define pi M_PI
#define abs(x) ((x) > 0) ? (x) : (-(x))

#define a 0.
#define b 1.
//double f(double x) {return (x==0.) ? 1 :(sin(x)/x)*(sin(x)/x);}
//double f(double x) {return (1.0/sqrt(pi))*exp(-x*x/4096);}
double f(double x) {return (4.0)/(1+x*x);}


int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	int size, rank, i, j, n, N;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (argc > 2) {sscanf(argv[1],"%d %d", &n, &N);} else {n = 64; N = 65536;}
	//printf("Hi I'm #%d. There is %d of us in total.\n", rank, size);
	double t0 = .0;
	if (!rank) {t0 = MPI_Wtime();}
	
	double h = (b-a)/N;
	double S = 0.;
	if (size == 1) {
		S = 0.;
		for (i = 1; i < N; i+=2) 
			S += (h/3)*(f(a+h*(i-1))+4*f(a+h*i)+f(a+h*(i+1)));
		printf("#%d:\tresult =\t%.16lf\n\terror = \t%.16f\n\ttime =\t\t%f sec\n", rank, S, abs(S-pi), MPI_Wtime()-t0);
		MPI_Finalize();
		return 0;
	}

	MPI_Status st;
	st.MPI_TAG = 1;
	int M = N/n;
	double H = (b-a)/n;
	double A = a;
	double s = 0.;
	if (!rank) {
		for (j = 1; j < size; j++) {
			MPI_Send((void *) &A, 1, MPI_DOUBLE, j, 1, MPI_COMM_WORLD);
			A += H;
		}
	}
	if (rank) {
		MPI_Recv((void *) &A, 1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
		while (st.MPI_TAG) {
			s = 0.;
			for (i = 1; i < M; i+=2) 
				s += (h/3)*(f(A+h*(i-1))+4*f(A+h*i)+f(A+h*(i+1)));
			//for (i = 0; i < M; i++) 
			//	s += .5*h*(f(A+h*i)+f(A+h*i+h));
			MPI_Send((void *) &s, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
			MPI_Recv((void *) &A, 1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
		} 
	}
	if (!rank) {
		while (j < n+1) {
			MPI_Recv((void *) &s, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
			S += s;
			MPI_Send((void *) &A, 1, MPI_DOUBLE, st.MPI_SOURCE, 1, MPI_COMM_WORLD);
			A += H;
			j++;
		}
		for (i = 1; i < size; i++) {
			MPI_Recv((void *) &s, 1, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
			S += s;
			MPI_Send((void *) &A, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}
	}
	if (!rank) printf("#%d:\tresult =\t%.16lf\n\terror = \t%.16f\n\ttime =\t\t%f sec\n", rank, S, abs(S-pi), MPI_Wtime()-t0);
	MPI_Finalize();
	return 0;
}
