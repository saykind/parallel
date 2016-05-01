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
#define H 0.
#define T0 0.
#define T1 2.

int main(int argc, char *argv[]) {
	int l;	int const N = 64;	int const NN = N*N;
	int const L = 256;		int const K = 8*1024*1024;
	int const KK = K/8;									// equalibration
	double *T = malloc(L*sizeof(double));
	double *E = malloc(L*sizeof(double));
	double *M = malloc(L*sizeof(double));
	for (l = 0; l < L; l++) {T[l] = 0.; E[l] = .0; M[l] = .0;}
	char s[N][N];
	int i, j, i0, j0, k, t;
	double e, h; e = .0; h = .0;
	int c, m; c = 0; m = 0;

#pragma omp parallel default(shared) private(s,i,j,k,t,i0,j0,e,h,c,m) num_threads(1)
{
	srand(time(NULL));
#pragma omp for 
	for (t = 0; t < L; t++)	{							// temperature change
		T[t] = T0+(T1-T0)/L*(t+1);
		e = .0;
		m = .0;
		for (i = 0; i < N; i++)	{						// grid initialization
			for (j = 0; j < N; j++) {
				s[i][j] = q-1;
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
		for (k = 1; k <= K; k++)	{					// Markov chain steps
			i0 = rand() % N;
			j0 = rand() % N;
			c = (s[i0][j0] + 1 + (rand() % (q-1))) % q; 			// changed value of random spin c=s'-s
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
			if (k > KK) {							// equalibration
				E[t] += e/NN;						// averaging instant values
				M[t] += 1.*m/NN;
			}
		}
		E[t] = E[t]/(K-KK);
		M[t] = M[t]/(K-KK);
	}
#pragma omp barrier
}
#ifdef PLOT
	FILE *gpp = gpinit();
	char output[40];
	fprintf(gpp, "set title 'Energy in Potts model (Z = 4, q = %d)'\n", q);
       	fprintf(gpp, "set xlabel 'T/J'\n set ylabel 'E/ZJ'\n");
	sprintf(output,"energy_q=%d_H=%.2lf.eps", q, H);
	plot(gpp, T, E, L, output);
	fprintf(gpp, "set title 'Magnetic Moment in Potts model (Z = 4, q = %d)'\n", q);
	fprintf(gpp, "set xlabel 'T/J'\n set ylabel 'M_{normalized}'\n");
	sprintf(output,"moment_q=%d_H=%.2lf.eps", q, H);
	plot(gpp, T, M, L, output);
#endif
#ifdef DATA
	char name[30];
	sprintf(name,"E_data_q=%d_H=%.2lf.dat", q, H);	FILE *efp = fopen(name, "a");
	sprintf(name,"M_data_q=%d_H=%.2lf.dat", q, H);	FILE *mfp = fopen(name, "a");
	for (l = 0; l < L; l++)	fprintf(efp,"%.14lf\t%.14lf\n", T[l], E[l]);
	for (l = 0; l < L; l++)	fprintf(mfp,"%.14lf\t%.14lf\n", T[l], M[l]);
	fclose(efp); fclose(mfp);
#endif
	free(T); free(E); free(M);
	return 0;
}
