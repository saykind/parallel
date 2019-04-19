#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include "params.h"

int simulate(int N, double complex **h, double complex **S, double **g, double **g2, int ***c) {
    int L=2*N+1, i;
    double a = 2.*pi/L, Y = (2*pi)*(2*pi)*(2*pi)/3*p8*p8;
    
    int k1 = rand()%L-N, k2 = rand()%L-N, q1, q2;
    if (!k1 && !k2) {simulate(N, h, S, g, g2, c); return 0;}
    
    double A = a*a*(k1*k1+k2*k2), d = d0/A; A = A*A;
    double complex **dS = (double complex **)malloc(L*sizeof(double complex*));
    for (i = 0; i < L; i++) *(dS+i) = (double complex *)malloc(L*sizeof(double complex));

    double complex z = (1.*rand()/RAND_MAX-.5) + (1.*rand()/RAND_MAX-.5)*I;
    double w = 0; dS[0][0] = 0;
    for (q1 = -N; q1 < N+1; q1++)
        for (q2 = -N; q2 < N+1; q2++) {
            if ( (q1==0) && (q2==0) )
                continue;
            int kq1 = (L-k1-q1)%L < N+1 ? (L-k1-q1)%L : (L-k1-q1)%L-L, kq2 = (L-k2-q2)%L < N+1 ? (L-k2-q2)%L : (L-k2-q2)%L-L;
            int qk1 = (L+k1-q1)%L < N+1 ? (L+k1-q1)%L : (L+k1-q1)%L-L, qk2 = (L+k2-q2)%L < N+1 ? (L+k2-q2)%L : (L+k2-q2)%L-L;
            dS[(q1+L)%L][(q2+L)%L] = a*a*(k1*q2-k2*q1)*(k1*q2-k2*q1)*( conj(h[(k1+q1+L)%L][(k2+q2+L)%L])*z );
            dS[(q1+L)%L][(q2+L)%L]+= a*a*(kq1*q2-kq2*q1)*(k1*q2-k2*q1)*( h[(2*L-k1-q1)%L][(2*L-k2-q2)%L]*z );
            dS[(q1+L)%L][(q2+L)%L]+= a*a*(qk1*q2-qk2*q1)*(k1*q2-k2*q1)*( h[(k1-q1+L)%L][(k2-q2+L)%L]*conj(z) );
            dS[(q1+L)%L][(q2+L)%L]+= a*a*(k1*q2-k2*q1)*(k1*q2-k2*q1)*( conj(h[(q1-k1+L)%L][(q2-k2+L)%L])*conj(z) );
            if (!((q1+2*k1)%L) && !((q2+2*k2)%L))
                dS[(q1+L)%L][(q2+L)%L] += a*a*(k1*q2-k2*q1)*(k1*q2-k2*q1)*z*z*d;
            if (!((q1-2*k1)%L) && !((q2-2*k2)%L))
                dS[(q1+L)%L][(q2+L)%L] += a*a*(k1*q2-k2*q1)*(k1*q2-k2*q1)*conj(z)*conj(z)*d;
            dS[(q1+L)%L][(q2+L)%L] *= d/(q1*q1+q2*q2)/(q1*q1+q2*q2);
            w += creal( (2*S[(q1+L)%L][(q2+L)%L]+dS[(q1+L)%L][(q2+L)%L])*conj(dS[(q1+L)%L][(q2+L)%L]) );
        }
    w *= -Y/L/L;
    w -= A*creal((2*h[(k1+L)%L][(k2+L)%L] + d*z)*conj(z))*d;
    if (w > log(1.*rand()/RAND_MAX)) {
        h[(k1+L)%L][(k2+L)%L] += d*z;
        h[(L-k1)%L][(L-k2)%L] += d*conj(z);
        for (q1 = 0; q1 < L; q1++)
            for (q2 = 0; q2 < L; q2++)
                S[q1][q2] += dS[q1][q2];
        if (c) {
            c[0][(k1+L)%L][(k2+L)%L]++;
            c[0][(L-k1)%L][(L-k2)%L]++;
        }
    }
    if (c) {
        a = creal(h[(k1+L)%L][(k2+L)%L]*conj(h[(k1+L)%L][(k2+L)%L]));
        g[(k1+L)%L][(k2+L)%L] += a;
        g2[(k1+L)%L][(k2+L)%L] += a*a;
        c[1][(k1+L)%L][(k2+L)%L]++;
        c[1][(L-k1)%L][(L-k2)%L]++;
    }  
    for (i = 0; i < L; i++) free(*(dS+i));
    free(dS);
    return 0;
}
