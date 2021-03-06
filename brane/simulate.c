#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include "params.h"

int simulate(double p8, int N, double complex **h, double complex **S, double **g, double **g2, int ***c) {
    int L=2*N+1, i;
    double a = 2.*pi/L, Y = (2*pi)/3*p8*p8;
    
    int k1 = rand()%L-N, k2 = rand()%L-N, q, q1, q2;
    if (!k1 && !k2) {simulate(p8, N, h, S, g, g2, c); return 0;}
    
    double A = a*a*(k1*k1+k2*k2), d = d0/A/pow(1+p8*p8/a/a/(k1*k1+k2*k2),.2); A = A*A;
    double complex **dS = (double complex **)malloc(L*sizeof(double complex*));
    for (i = 0; i < L; i++) *(dS+i) = (double complex *)malloc(L*sizeof(double complex));
    dS[0][0] = 0;
    double complex z = (1.*rand()/RAND_MAX-.5) + (1.*rand()/RAND_MAX-.5)*I;
    double w = 0; 
    #pragma omp parallel for private(q,q1,q2) reduction(+:w)
    for (q = 0; q < L*L; q++) {
        q1 = q/L-N; q2= q%L-N;
        if (!q1 && !q2) continue;
        int kq1 = (L-k1-q1)%L < N+1 ? (L-k1-q1)%L : (L-k1-q1)%L-L;
        int kq2 = (L-k2-q2)%L < N+1 ? (L-k2-q2)%L : (L-k2-q2)%L-L;
        int qk1 = (L+k1-q1)%L < N+1 ? (L+k1-q1)%L : (L+k1-q1)%L-L;
        int qk2 = (L+k2-q2)%L < N+1 ? (L+k2-q2)%L : (L+k2-q2)%L-L;
        double complex s;
        s  =  (k1*q2-k2*q1) * (k1*q2-k2*q1) *( conj(h[(k1+q1+L)%L][(k2+q2+L)%L])*z );
        s += (kq1*q2-kq2*q1)*(kq1*q2-kq2*q1)*( h[(L-k1-q1)%L][(L-k2-q2)%L]*z   );
        s += (qk1*q2-qk2*q1)*(qk1*q2-qk2*q1)*( h[(k1-q1+L)%L][(k2-q2+L)%L]*conj(z) );
        s +=  (k1*q2-k2*q1) * (k1*q2-k2*q1) *( conj(h[(q1-k1+L)%L][(q2-k2+L)%L])*conj(z) );
        if (!((q1+2*k1)%L) && !((q2+2*k2)%L))
            s += (k1*q2-k2*q1)*(k1*q2-k2*q1)*z*z*d;
        if (!((q1-2*k1+2*L)%L) && !((q2-2*k2+2*L)%L))
            s += (k1*q2-k2*q1)*(k1*q2-k2*q1)*conj(z)*conj(z)*d;
        s *= a*a*d/(q1*q1+q2*q2);
        dS[(q1+L)%L][(q2+L)%L] = s;
        w += creal( (2*S[(q1+L)%L][(q2+L)%L]+s)*conj(s) );
    }
    w *= -Y/L/L;
    w -= A*creal((2*h[(k1+L)%L][(k2+L)%L] + d*z)*conj(z))*d;
    if (w > log(1.*rand()/RAND_MAX)) {
        h[(k1+L)%L][(k2+L)%L] += d*z;
        h[(L-k1)%L][(L-k2)%L] += d*conj(z);      
        #pragma omp parallel for private(q)
        for (q = 0; q < L*L; q++)
            S[q/L][q%L] += dS[q/L][q%L];
        if (c) {
            c[0][(k1+L)%L][(k2+L)%L]++;
            c[0][(L-k1)%L][(L-k2)%L]++;
        }
    }
    if (c && g && g2) {
        a = creal(h[(k1+L)%L][(k2+L)%L]*conj(h[(k1+L)%L][(k2+L)%L]));
        g[(k1+L)%L][(k2+L)%L] += a;
        g[(L-k1)%L][(L-k2)%L] += a;
        g2[(k1+L)%L][(k2+L)%L] += a*a;
        g2[(L-k1)%L][(L-k2)%L] += a*a;
        c[1][(k1+L)%L][(k2+L)%L]++;
        c[1][(L-k1)%L][(L-k2)%L]++;
    }  
    for (i = 0; i < L; i++) free(*(dS+i));
    free(dS);
    return 0;
}
