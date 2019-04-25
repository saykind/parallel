#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include "params.h"

int init(int N, double complex **h, double complex **S, double **g, double **g2, int ***c) {
    int L=2*N+1, q1, q2, k1, k2;
    char name[20]; 
    sprintf(name,"data/N=%d.dat", N);
    FILE *file = fopen(name, "r");
    if (file) {
        double re, im;
        for (q1 = 0; q1 < L; q1++)
            for (q2 = 0; q2 < L; q2++) {
                fscanf(file,"%lf\t%lf\n", &re, &im );
                h[q1][q2] = re + im*I;
                fscanf(file,"%lf\t%lf\n", &re, &im );
                S[q1][q2] = re + im*I;
                fscanf(file,"%lf\t%lf\n", &g[q1][q2], &g2[q1][q2] );
                fscanf(file,"%d\t%d\n", &c[0][q1][q2], &c[1][q1][q2] );
            }        
        fclose(file);
        return 1;
    }
    double a = 2*pi/L;
    h[0][0] = 1/a/a;
    #pragma omp parallel for collapse(2) private(q1,q2,k1,k2)
    for (q1 = 0; q1 < L; q1++)
         for (q2 = 0; q2 < L; q2++) {
            g[q1][q2] = 0;
            g2[q1][q2] = 0;
            c[0][q1][q2] = 0;
            c[1][q1][q2] = 0;
            S[q1][q2] = 0;
            if (!q1 && !q2) continue;
            k1 = (q1<(N+1)) ? q1 : (q1-L), k2 = (q2<(N+1)) ? q2 : (q2-L);
            h[q1][q2] = 1/a/a/(k1*k1+k2*k2);
        }  
    #pragma omp parallel for collapse(2) private(q1,q2,k1,k2)
    for (q1 = -N; q1 < N+1; q1++)
        for (q2 = -N; q2 < N+1; q2++) {
            if (!q1 && !q2) continue;
            for (k1 = -N; k1 < N+1; k1++)
                for (k2 = -N; k2 < N+1; k2++)   {
                    double p = a*a*(k1*q2-k2*q1)*(k1*q2-k2*q1)/(q1*q1+q2*q2);
                    S[(q1+L)%L][(q2+L)%L] += p*h[(k1+L)%L][(k2+L)%L]*conj(h[(k1+q1+L)%L][(k2+q2+L)%L]);
                }
        }
    return 0;
}

int dump(int N, double complex **h, double complex **S, double **g, double **g2, int ***c) {    
    int L=2*N+1, q1, q2;
    char name[20];
    sprintf(name,"data/N=%d.dat", N); 
    FILE *file = fopen(name, "w+"); if (!file) {printf("Cannot save data\n"); return -1;}
    for (q1 = 0; q1 < L; q1++)
        for (q2 = 0; q2 < L; q2++) {
            fprintf(file,"%.15lf\t%.15lf\n", creal(h[q1][q2]), cimag(h[q1][q2]) );
            fprintf(file,"%.15lf\t%.15lf\n", creal(S[q1][q2]), cimag(S[q1][q2]) );
            fprintf(file,"%lf\t%lf\n", g[q1][q2], g2[q1][q2] );
            fprintf(file,"%d\t%d\n", c[0][q1][q2], c[1][q1][q2] );
        }
    fclose(file);
    return 0;
}
