#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include "plot.h"
#define PLOT 
#define DATA 

#define q 3
#define J -1.0
#define H .15
#define T1 4.0

int main(int argc, char *argv[]) {
	int l;	int const N = 64;	int const NN = N*N;
	int const L = 512;	int const K = 1024*1024;
	double *T = malloc(L*sizeof(double));
	double *E = malloc(L*sizeof(double));
	double *M = malloc(L*sizeof(double));
	for (l = 0; l < L; l++) {T[l] = 0.; E[l] = .0; M[l] = .0;}
	char **s;
	int i, j, i0, j0, k, t;
	double e, h; e = .0; h = .0;
	int c, m; c = 0; m = 0;

#pragma omp parallel default(shared) private(s,i,j,k,t,i0,j0,e,h,c,m) num_threads(333)
{
	s = malloc(N*sizeof(char*));
	for (l = 0; l < N; l++)	*(s+l) = malloc(N*sizeof(char));
	
	srand(time(NULL));
#pragma omp for 
	for (t = 0; t < L; t++)	{							// temperature change
		e = .0;
		m = .0;
		T[t] = T1/L*(t+1);
		for (i = 0; i < N; i++)	{						// grid initialization
			for (j = 0; j < N; j++) {
				s[i][j] = rand() % q;
			}
		}
		for (i = 0; i < N; i++)							// initial instant values of energy e and total moment m
			for (j = 0; j < N; j++) {
				if (s[i][j] == s[(i+1)%N][j])	e += J;
				if (s[i][j] == s[i][(j+1)%N])	e += J;
				if (s[i][j] == s[(i-1+N)%N][j])	e += J;
				if (s[i][j] == s[i][(j-1+N)%N])	e += J;
								m += s[i][j];
			}
		e += -H*m;
		E[t] += e/NN;		
		M[t] += 1.*m/NN;
		for (k = 1; k < K; k++)	{						// Markov chain steps
			i0 = rand() % N;
			j0 = rand() % N;
			c = (s[i0][j0] + 1) % q; 					// changed value of random spin c=s'-s
			h = .0;
			if (s[i0][j0] == s[(i0+1)%N][j0])	h -= J;			// change in energy h=e'-e
			else if (c == s[(i0+1)%N][j0])		h += J;
			if (s[i0][j0] == s[i0][(j0+1)%N])	h -= J;
			else if (c == s[i0][(j0+1)%N])		h += J;
			if (s[i0][j0] == s[(i0-1+N)%N][j0])	h -= J;
			else if (c == s[(i0-1+N)%N][j0])	h += J;
			if (s[i0][j0] == s[i0][(j0-1+N)%N])	h -= J;
			else if (c == s[i0][(j0-1+N)%N])	h += J;
			h += -H*(c-s[i0][j0]);			
			if (rand() < 1.0/(1.0+exp(h/T[t]))*RAND_MAX) {			// whether new state is accepted
				e += h; 
				m += c-s[i0][j0]; 
				s[i0][j0] = c;
			}
			E[t] += e/NN;							// averaging instant values
			M[t] += 1.*m/NN;
		}
		E[t] = E[t]/K;
		M[t] = M[t]/K;
	}
#pragma omp barrier
	for (l = 0; l < N; l++)	free(*(s+l));
	free(s);
}
#ifdef PLOT
	FILE *gpp = gpinit();
	fprintf(gpp, "set title 'Energy'\n set xlabel 'T/J'\n set ylabel 'E/J'\n");
	fprintf(gpp, "set label 1 'H = %.2lf' at 2.3, -2.5\n", H);
	plot(gpp, T, E, L, "energy.eps");
	fprintf(gpp, "set title 'Moment'\n set xlabel 'T/J'\n set ylabel 'M/J'\n");
	fprintf(gpp, "set label 1 'H = %.2lf' at 2.3, 1.55\n", H);
	plot(gpp, T, M, L, "moment.eps");
#endif
#ifdef DATA
	char name[20];
	sprintf(name,"E_data_H=%.2lf.dat", H);
	FILE *efp = fopen(name, "a");
	sprintf(name,"M_data_H=%.2lf.dat", H);
	FILE *mfp = fopen(name, "a");
	for (l = 0; l < N; l++)	fprintf(efp,"%.14lf\t%.14lf\n", T[l], E[l]);
	for (l = 0; l < N; l++)	fprintf(mfp,"%.14lf\t%.14lf\n", T[l], M[l]);
#endif
	free(T);
	free(E);
	free(M);
	return 0;
}
