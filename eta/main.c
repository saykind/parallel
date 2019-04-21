#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include "params.h"

int simulate(double p8, int N, double complex **h, double complex **S, double **g, double **g2, int ***ac);

int main(int argc, char *argv[]) {
    int nt = nt0, N = N0, M=M0, MTH=MT0, L=2*N+1, i, j, i0, j0, k1, k2, q1, q2;
    double p8 = p80;
    char thrmlzd = 0, verbose = 0;
    for(i = 1; i < argc; i++ ) {
        if (!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help")) {
            printf("Usage:\n -h\tprint this message\n -v\tverbose mode [not recommended for N<20]\n -d\tprint default values\n");
            printf(" nt=%%d\tset number of threads\n");
            printf(" N=%%d\tset linear lattice size L=2*N+1 total size L*L\n");
            printf(" MTH=%%d\tset number of thermalization MC steps (total: MTH*L*L) [recommended value MTH=2000]\n");
            printf(" M=%%d\tset number of MC steps (total: M*L*L)\n");
            printf(" p8=%%f\tset interaction force (0 < p8 < 6.28)\n");
            return 0;
        } 
        if (!strcmp(argv[i],"-v")) {
            verbose = 1;
            continue;
        }
        if (!strcmp(argv[i],"-d") || !strcmp(argv[i],"--default") ) {
            printf("Default values:\n");
            printf(" N=%d\t L=%d\t M=%d\t MTH=%d thrds=%d\t p8=%.1lf\n", N, L, M, MTH, nt, p8);
            return 0;
        }  
        if (argv[i][0] == '-') {
            printf("Unknown option '%s'\n", argv[i]);
            return 0;
        }           
        sscanf(argv[i], "nt=%d", &nt);
        sscanf(argv[i], "N=%d", &N);
        sscanf(argv[i], "M=%d", &M);
        sscanf(argv[i], "MTH=%d", &MTH);
        sscanf(argv[i], "p8=%lf", &p8);
    }
    omp_set_num_threads(nt);
    L=2*N+1;    
    printf("\nN =\t%d\nL =\t%d\nthrds =\t%d\np8 =\t%.1lf\n", N, L, nt, p8);

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
    double t0 = omp_get_wtime();
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
#pragma omp parallel default(shared) private(q1, q2, k1, k2, i, j, i0, j0, x)
{
        #pragma omp for nowait
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
        #pragma omp barrier   
        #pragma omp for nowait
        for (q1 = -N; q1 < N+1; q1++)
            for (q2 = -N; q2 < N+1; q2++) {
                if ( (q1==0) && (q2==0) )
                    continue;
                for (k1 = -N; k1 < N+1; k1++)
                    for (k2 = -N; k2 < N+1; k2++)   {
                        double p = a*(k1*q2-k2*q1)/(q1*q1+q2*q2);
                        S[(q1+L)%L][(q2+L)%L] += p*p*h[(k1+L)%L][(k2+L)%L]*conj(h[(k1+q1+L)%L][(k2+q2+L)%L]);
                    }
            }
    }
}
    printf("init =\t%.1lf min\n", (omp_get_wtime()-t0)/60);
    
	/* Monte Carlo simulation */
	time_t rawtime;
	/* Thermalization */
	t0 = omp_get_wtime();
	double T0 = 22e-8;
	srand(time(NULL));
	if (!thrmlzd) {
	    time(&rawtime);
	    printf("\nnow =\t%s", asctime (localtime ( &rawtime )) );
        printf("MTH =\t%d*%d\test: %.0f/%d min\n", MTH, L*L, T0*MT0*L*L*L*L/60, nt );
        rawtime += T0*MT0*L*L*L*L;
        printf("fin =\t%s", asctime (localtime ( &rawtime )) );
        for (i = 0; i < MTH; i++) {
		    for (j = 0; j < L*L; j++)
                simulate(p8, N, h, S, NULL, NULL, NULL);
     		if (verbose) {
		        double t = (double) omp_get_wtime()-t0, pc = 1.*(i+1)/MTH;
		        printf("  %.1lf min  \t%.0lf%%\test: %.1f min\tleft: %.1f min\t \r", t/60, 100.*pc, t/pc/60, (1-pc)*(t/pc)/60); 
		        fflush(stdout);
		    }
        }
        printf("\ntime =\t%.2lf min\n", (omp_get_wtime()-t0)/60);
        t0 = omp_get_wtime();
    }
	time(&rawtime);
	printf("\nnow =\t%s", asctime (localtime ( &rawtime )) );
    /* Monte Carlo average */
    printf("M =\t%d*%d\test: %.0f/%d min\n", M, L*L, T0*M*L*L*L*L/60, nt);
    rawtime += T0*M*L*L*L*L;
    printf("fin =\t%s", asctime (localtime ( &rawtime )) );
    for (i = 0; i < M; i++) {    
    	for (j = 0; j < L*L; j++)
        	simulate(p8, N, h, S, g, g2, c);
	    if (verbose) {
	        double t = (double) omp_get_wtime()-t0, pc = 1.*(i+1)/M;
	        printf("  %.1lf min  \t%.0lf%%\test: %.1f min\tleft: %.1f min\t \r", t/60, 100.*pc, t/pc/60, (1-pc)*(t/pc)/60);
	        fflush(stdout);
	    }
    }
    printf("\ntime =\t%.2lf min\n", (omp_get_wtime()-t0)/60);
    time(&rawtime);
	printf("\nnow =\t%s", asctime (localtime ( &rawtime )) );

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
	printf("\n");
    return 0;
}
