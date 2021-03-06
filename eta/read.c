#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>
#include "params.h"

FILE   *gpinit(void);

int main(int argc, char *argv[]) {
    int N = N0, i, j, d3=0, verbose=0;
    double p8 = p80;
    for(i = 1; i < argc; i++ ) {
        if (!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help")) {
            printf("Data plotting. Usage:\n -h\tprint this message\n -v\tverbose mode [not recommended for N<20]\n -d\tprint default values\n");
            printf(" N=%%d\tset linear lattice size L=2*N+1 total size L*L\n");
            printf(" p8=%%f\tset interaction force for fit (0 < p8 < 6.28)\n");
            return 0;
        } 
        if (!strcmp(argv[i],"-3d")) {
            d3 = 1;
            continue;
        }
        if (!strcmp(argv[i],"-v")) {
            verbose = 1;
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
                    if (verbose)
                        printf("[%d,%d]\t%.3f +/- %.3f  \t(%.0f %% of %d times)\n", i, j, g[i][j], g2[i][j], 100.*c[0][i][j]/c[1][i][j], c[1][i][j]);
                    g2[i][j] /= g[i][j];
                } else if(i || j) {printf("c[1][%d][%d] == 0\n", i, j);}
            }        
        fclose(file);
    } else {
        printf("Nothing to read :(\n");
        return -1;
    }
    if (!verbose)
        printf("\n[%d,%d]\t%.0f%%\n[%d,%d]\t%.0f%%\n[%d,%d]\t%.0f%%\n[%d,%d]\t%.0f%%\n[%d,%d]\t%.0f%%\n[%d,%d]\t%.0f%%\n\n",\
         1,1,100.*c[0][1][1]/c[1][1][1], 2,2,100.*c[0][2][2]/c[1][2][2], 4,4,100.*c[0][4][4]/c[1][4][4],\
         N,N,100.*c[0][N][N]/c[1][N][N], L-4,L-4,100.*c[0][L-4][L-4]/c[1][L-4][L-4], L-2,L-2,100.*c[0][L-2][L-2]/c[1][L-2][L-2]);
    
    /* Plots */
	FILE *gpp = gpinit();	
	fprintf(gpp, "set xrange [ %lf : %lf ]\n", 1.9*pi/50, 1.1*pi);
    char output[20]; 
    double Q[N], Gx[N], Gy[N], Gr[N];
    for (i = 0; i < N; i++) {
        Q[i] = (i+1)*2*pi/L;
        Gx[i] = 1./g[i+1][0];
        Gy[i] = 1./g[0][i+1];
        Gr[i] = 1./g[i+1][i+1];
    }
    fprintf(gpp, "set title 'Inverse Green function (L=%d, p*=%.1lf, M=%d)'\n", L, p8, c[1][N][N]);
    fprintf(gpp, "set logscale xy\n");
    // x-axis
    sprintf(output,"Lx=%d.eps", L); fprintf(gpp, "set output '%s'\n", output);
    fprintf(gpp, "set xlabel 'q_x'\n set ylabel 'G^{-1}'\n");
	fprintf(gpp, "plot x**4, x**4*(1+(%f/x)**2)**(.4), x**4*%f**.8/x**.8, '-' w errorbars pt 7 ps .5\n", p8, p8);
	for (i = 0; i < N; i++)
		fprintf(gpp, "%.6f %.6f %.6f\n", Q[i], Gx[i], Gx[i]*g2[i+1][0]);
	fprintf(gpp,"e\n");
	fflush(gpp);
	// y-axis
    sprintf(output,"Ly=%d.eps", L); fprintf(gpp, "set output '%s'\n", output);
    fprintf(gpp, "set xlabel 'q_y'\n set ylabel 'G^{-1}'\n");
	fprintf(gpp, "plot x**4, x**4*(1+(%f/x)**2)**(.4), x**4*%f**.8/x**.8, '-' w errorbars pt 7 ps .5\n", p8, p8);
	for (i = 0; i < N; i++)
		fprintf(gpp, "%.6f %.6f %.6f\n", Q[i], Gy[i], Gy[i]*g2[0][i+1]);
	fprintf(gpp,"e\n");
	fflush(gpp);
	// r-axis
    sprintf(output,"Lr=%d.eps", L); fprintf(gpp, "set output '%s'\n", output);
    fprintf(gpp, "set xlabel 'q'\n set ylabel 'G^{-1}'\n");
	fprintf(gpp, "plot x**4, x**4*(1+%f**2/x**2)**(.4), x**4*%f**.8/(x**2)**.4, '-' w errorbars pt 7 ps .5\n", p8, p8);
	for (i = 0; i < N; i++)
		fprintf(gpp, "%.6f %.6f %.6f\n", sqrt(2)*(i+1)*2*pi/L, Gr[i], Gr[i]*g2[i+1][i+1]);
	fprintf(gpp,"e\n");
	fflush(gpp);	
	// 3d plot
	if (d3) {
	    fprintf(gpp, "set term x11\n");
	    fprintf(gpp, "unset logscale xy\n");
	    fprintf(gpp, "set xrange [ * : * ]\n");
	    fprintf(gpp, "splot (x**2+y**2)**2, '-' w points pt 7 ps 1\n");
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

