#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "params.h"
#include "plot.h"

int main(int argc, char *argv[]) {
    int N = N0, i, j;
    if (argc > 1) {sscanf(argv[1], "%d", &N);}
    int L =2*N+1;

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
                fscanf(file,"%lf\t%lf\n", &re, &im );
                fscanf(file,"%lf\t%lf\n", &g[i][j], &g2[i][j] );
                fscanf(file,"%d\t%d\n", &c[0][i][j], &c[1][i][j] );
                if (c[1][i][j]) {
                    c[0][i][j] /= c[1][i][j];
                    c[0][i][j] /= c[1][i][j];
                    g[i][j]  /= c[1][i][j];
                    g2[i][j] /= c[1][i][j];
                } else {printf("c[1][%d][%d] == 0\n", i, j);}
            }        
        fclose(file);
    } else {
        printf("Nothing to read :(\n");
        return -1;
    }
    
    /* Plots */
	FILE *gpp = gpinit();
    char output[20]; 
    double Q[N], G[N];
    
    sprintf(output,"Lx=%d.eps", L);
	fprintf(gpp, "set title 'Inverse Green function (L = %d)'\n", L);
    fprintf(gpp, "set xlabel 'q_x'\n set ylabel 'G^{-1}'\n");
    fprintf(gpp, "set logscale xy\n");
    for (i = 0; i < N; i++) {
        Q[i] = i+1;
        G[i] = 1./g[i+1][0];
    }
    fprintf(gpp, "set output '%s'\n", output);
	fprintf(gpp, "plot %f*x**4, '-' w points pt 1 ps .5\n", 2*(2.*pi/L)*(2.*pi/L)*(2.*pi/L)*(2.*pi/L));
	for (i = 0; i < N; i++)
		fprintf(gpp, "%.6f %.6f\n", Q[i], G[i]);
	fprintf(gpp,"e\n");
	fflush(gpp);
	
    sprintf(output,"Ly=%d.eps", L);
	fprintf(gpp, "set title 'Inverse Green function (L = %d)'\n", L);
    fprintf(gpp, "set xlabel 'q_y'\n set ylabel 'G^{-1}'\n");
    fprintf(gpp, "set logscale xy\n");
    for (i = 0; i < N; i++) {
        Q[i] = i+1;
        G[i] = 1./g[0][i+1];
    }
    fprintf(gpp, "set output '%s'\n", output);
	fprintf(gpp, "plot %f*x**4, '-' w points pt 1 ps .5\n", 2*(2.*pi/L)*(2.*pi/L)*(2.*pi/L)*(2.*pi/L));
	for (i = 0; i < N; i++)
		fprintf(gpp, "%.6f %.6f\n", Q[i], G[i]);
	fprintf(gpp,"e\n");
	fflush(gpp);
	
	pclose(gpp);

    /* Memory de-allocation */ 
   	for (i = 0; i < L; i++) {free(*(h+i)); free(*(g+i)); free(*(g2+i)); free(*(S+i)); free(*(*(c+0)+i)); free(*(*(c+1)+i));}
	free(h); free(g); free(g2); free(S); free(*(c+0)); free(*(c+1));
	free(c);
    return 0;
}
