#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>
#include "params.h"

FILE   *gpinit(void);

int main(int argc, char *argv[]) {
    int N = N0, i, j, d3 = 0;
    double p8 = p80;
    for(i = 1; i < argc; i++ ) {
        if (!strcmp(argv[i],"-3d")) {
            d3 = 1;
            continue;
        }
        if (argv[i][0] == '-') {
            printf("Unknown option '%s'\n", argv[i]);
            return 0;
        }
        sscanf(argv[i], "N=%d", &N);
        sscanf(argv[i], "p8=%lf", &p8);
    }    
    int L =2*N+1;

    /* Memory allocation */
	double **g = (double **)malloc(L*sizeof(double *));
	double **g2 = (double **)malloc(L*sizeof(double *));
    int ***c = malloc(2*sizeof(int**));
    *(c+0) = (int **)malloc(L*sizeof(int*));
    *(c+1) = (int **)malloc(L*sizeof(int*));
	for (i = 0; i < L; i++) {
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
                    g[i][j]  /= c[1][i][j];
                    g2[i][j] /= c[1][i][j];
                    g2[i][j] = 4*sqrt(fabs((g2[i][j]-g[i][j]*g[i][j])/c[1][i][j]));
                    printf("[%d,%d]\t%.3f +/- %.3f  \t(%.0f %% of %d times)\n", i, j, g[i][j], g2[i][j], 100.*c[0][i][j]/c[1][i][j], c[1][i][j]);
                    g2[i][j] /= g[i][j];
                } else {printf("c[1][%d][%d] == %d\n", i, j, N);}
            }        
        fclose(file);
    } else {
        printf("Nothing to read :(\n");
        return -1;
    }
    
    /* Plots */
	FILE *gpp = gpinit();	
	fprintf(gpp, "set xrange [ %lf : %lf ]\n", 1.8*pi/L, 1.1*pi);
    char output[20]; 
    double Q[N], Gx[N], Gy[N], Gr[N], A=2;
    for (i = 0; i < N; i++) {
        Q[i] = (i+1)*2*pi/L;
        Gx[i] = 1./g[i+1][0];
        Gy[i] = 1./g[0][i+1];
        Gr[i] = 1./g[i+1][i+1];
    }
    fprintf(gpp, "set title 'Inverse Green function (L = %d, p*=%.0lf)'\n", L, p8);
    fprintf(gpp, "set logscale xy\n");
    // x-axis
    sprintf(output,"Lx=%d.eps", L); fprintf(gpp, "set output '%s'\n", output);
    fprintf(gpp, "set xlabel 'q_x'\n set ylabel 'G^{-1}'\n");
	fprintf(gpp, "plot %f*x**4, %f*x**4*(1+(%f/x)**2)**(.2), '-' w errorbars pt 7 ps .5\n", A, A, p8);
	for (i = 0; i < N; i++)
		fprintf(gpp, "%.6f %.6f %.6f\n", Q[i], Gx[i], Gx[i]*g2[i+1][0]);
	fprintf(gpp,"e\n");
	fflush(gpp);
	// y-axis
    sprintf(output,"Ly=%d.eps", L); fprintf(gpp, "set output '%s'\n", output);
    fprintf(gpp, "set xlabel 'q_y'\n set ylabel 'G^{-1}'\n");
	fprintf(gpp, "plot %f*x**4, %f*x**4*(1+(%f/x)**2)**(.2), '-' w errorbars pt 7 ps .5\n", A, A, p8);
	for (i = 0; i < N; i++)
		fprintf(gpp, "%.6f %.6f %.6f\n", Q[i], Gy[i], Gy[i]*g2[0][i+1]);
	fprintf(gpp,"e\n");
	fflush(gpp);
	// r-axis
    sprintf(output,"Lr=%d.eps", L); fprintf(gpp, "set output '%s'\n", output);
    fprintf(gpp, "set xlabel 'q'\n set ylabel 'G^{-1}'\n");
	fprintf(gpp, "plot 4*%f*x**4, 4*%f*x**4*(1+(%f/x)**2)**(.2), '-' w errorbars pt 7 ps .5\n", A, A, p8);
	for (i = 0; i < N; i++)
		fprintf(gpp, "%.6f %.6f %.6f\n", Q[i], Gr[i], Gr[i]*g2[i+1][i+1]);
	fprintf(gpp,"e\n");
	fflush(gpp);	
	// 3d plot
	if (d3) {
	    fprintf(gpp, "set term x11\n");
	    fprintf(gpp, "unset logscale xy\n");
	    fprintf(gpp, "set xrange [ * : * ]\n");
	    fprintf(gpp, "splot %f*(x**2+y**2)**2, '-' w points pt 7 ps 1\n", A);
	    for (i = -N; i < N+1; i++)
	        for (j = -N; j< N+1; j++) {
	            if (!i && !j) {continue;}
		        fprintf(gpp, "%d %d %.6f\n", i, j, 1./g[(i+L)%L][(j+L)%L]);
		    }
	    fprintf(gpp,"e\n");
	    fflush(gpp);	
	    scanf("%c",output);
    }
    /* Memory de-allocation */ 
   	for (i = 0; i < L; i++) {free(*(g+i)); free(*(g2+i)); free(*(*(c+0)+i)); free(*(*(c+1)+i));}
	free(g); free(g2); free(*(c+0)); free(*(c+1));
	free(c);
	pclose(gpp);
    return 0;
}

FILE *gpinit(void) {
	FILE *gpp = popen("gnuplot", "w");
	if (gpp == NULL) {return NULL;}
	fprintf(gpp, "set term postscript eps enhanced color font 'Helvetica, 20'\n");
	fprintf(gpp, "unset key\n");
	fprintf(gpp, "unset border\n");
	fprintf(gpp, "set grid\n");
	fprintf(gpp, "set style fill solid\n");
	fprintf(gpp, "set xrange [ 1 : * ]\n");
	fprintf(gpp, "set yrange [ * : * ]\n");
	return gpp;
}

