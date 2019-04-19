#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include "params.h"

int simulate(int N, double complex **h, double complex **S, double **g, double **g2, int ***ac);

int main(int argc, char *argv[]) {
    int nt = 1, N = N0, M=M0, i, j, i0, j0, k1, k2, q1, q2;
    char thrmlzd = 0;
    if (argc > 1) {sscanf(argv[1], "%d", &nt);}
    if (argc > 2) {sscanf(argv[2], "%d", &N);}
    if (argc > 3) {sscanf(argv[3], "%d", &M);}
    omp_set_num_threads(nt);
    int L=2*N+1;    
    printf("N =\t%d\nL =\t%d\nthrds =\t%d\n", N, L, nt);

    /* Memory allocation */
    double complex **h = (double complex **)malloc(L*sizeof(double complex*));
    double complex **S = (double complex **)malloc(L*sizeof(double complex*));
	double **g = (double **)malloc(L*sizeof(double *));
	double **g2 = (double **)malloc(L*sizeof(double *));
    int ***c = malloc(2*sizeof(int**));
    *(c+0) = (int **)malloc(L*sizeof(int*));
    *(c+1) = (int **)malloc(L*sizeof(int*));
	for (i = 0; i < L; i++) {
	    *(h+i) = (double complex *)malloc(L*sizeof(double complex));
	    *(S+i) = (double complex *)malloc(L*sizeof(double complex));
	    *(g+i) = (double *)malloc(L*sizeof(double));
	    *(g2+i) = (double *)malloc(L*sizeof(double));
	    *(*(c+0)+i) = (int *)malloc(L*sizeof(int));
	    *(*(c+1)+i) = (int *)malloc(L*sizeof(int));
	}

    
    /* Data initialization */
    char name[10]; 
    sprintf(name,"N=%d.dat", N);
    FILE *file = fopen(name, "r");
    if (file) {
        double re, im;
        for (i = 0; i < L; i++)
            for (j = 0; j < L; j++) {
                fscanf(file,"%lf\t%lf\n", &re, &im );
                h[i][j] = re + im*I;
                fscanf(file,"%lf\t%lf\n", &re, &im );
                S[i][j] = re + im*I;
                fscanf(file,"%lf\t%lf\n", &g[i][j], &g2[i][j] );
                fscanf(file,"%d\t%d\n", &c[0][i][j], &c[1][i][j] );
                //printf("%lf\t%lf\n", creal(h[i][j]), cimag(h[i][j]) );
                //printf("%d\t%d\n", c[0][i][j], c[1][i][j] );
            }        
        fclose(file);
        thrmlzd = 1;
    } else {
        double a = 2*pi/L, x;
        for (i = 0; i < L; i++)
             for (j = 0; j < L; j++) {
                if ( (i==0) && (j==0) ) {
                    h[i][j] = 1./a/a;
                    continue;
                }
                i0 = (i<(N+1)) ? i : (i-L);
                j0 = (j<(N+1)) ? j : (j-L);
                x = a*a*( i0*i0 + j0*j0 )*( i0*i0 + j0*j0 );
                h[i][j] = 1./x;
                g[i][j] = 0;
                g2[i][j] = 0;
                c[0][i][j] = 0;
                c[1][i][j] = 0;
                S[i][j] = 0;
            }
        for (q1 = -N; q1 < N+1; q1++)
            for (q2 = -N; q2 < N+1; q2++) {
                if ( (q1==0) && (q2==0) )
                    continue;
                for (k1 = -N; k1 < N+1; k1++)
                    for (k2 = -N; k2 < N+1; k2++)
                        S[(q1+L)%L][(q2+L)%L] += a*a*(k1*q2-k2*q1)*(k1*q2-k2*q1)/(q1*q1+q2*q2)*h[(k1+L)%L][(k2+L)%L]*conj(h[(k1+q1+L)%L][(k2+q2+L)%L]);
            }
    }
    
	/* Monte Carlo simulation */
	/* Thermalization */
	double t0 = omp_get_wtime();
	double T0 = 24e-8;
	srand(time(NULL));
	if (!thrmlzd) {
        printf("MTH =\t%d\t(est: %.2f min)\n", MT0*L*L, T0*MT0*L*L*L*L/60 );
        for (i = 0; i < MT0*L*L; i++) {
 		    #ifdef VERBOSE
			    double t = (double) omp_get_wtime()-t0, pc = 1.*(i+1)/L/L/MT0;
			    printf(" \t %.0lf sec  \t%.0lf%%\test: %.0lf sec\r", t, 100.*pc, t/pc); 
			    fflush(stdout);
		    #endif
            simulate(N, h, S, g, g2, NULL);
        }
        printf("time =\t%.1lf sec\n", omp_get_wtime()-t0);
        t0 = omp_get_wtime();
    }
    /* Monte Carlo average */
    printf("M =\t%d\t(est: %.2f min)\n", M*L*L, T0*M*L*L*L*L/60);
    for (i = 0; i < M*L*L; i++)    {
		#ifdef VERBOSE
			double t = (double) omp_get_wtime()-t0, pc = 1.*(i+1)/L/L/M;
			printf(" \t %.0lf sec  \t%.0lf%%\test: %.0lf sec\r", t, 100.*pc, t/pc); 
			fflush(stdout);
		#endif
        simulate(N, h, S, g, g2, c);
    }
    printf("time =\t%.1lf sec\n", omp_get_wtime()-t0);

    /* Data dump */
    file = fopen(name, "w+"); if (!file) {printf("Cannot save data\n"); return -1;}
    for (i = 0; i < L; i++)
        for (j = 0; j < L; j++) {
            fprintf(file,"%.15lf\t%.15lf\n", creal(h[i][j]), cimag(h[i][j]) );
            fprintf(file,"%.15lf\t%.15lf\n", creal(S[i][j]), cimag(S[i][j]) );
            fprintf(file,"%lf\t%lf\n", g[i][j], g2[i][j] );
            fprintf(file,"%d\t%d\n", c[0][i][j], c[1][i][j] );
        }
    fclose(file);

    /* Memory dellocation */ 
   	for (i = 0; i < L; i++) {free(*(h+i)); free(*(g+i)); free(*(g2+i)); free(*(S+i)); free(*(*(c+0)+i)); free(*(*(c+1)+i));}
	free(h); free(g); free(g2); free(S); free(*(c+0)); free(*(c+1));
	free(c);
    return 0;
}
