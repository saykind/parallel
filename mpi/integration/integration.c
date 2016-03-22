#include <stdio.h>
#include <stdlib.h>
#include "plot.h"

double rectangle(double *x, double *f, int N){
	int i;	double I = .0;
	for (i = 0; i < N; i++) {I += (x[i+1]-x[i])*(f[i]);}
	return I;
}
double trapezoid(double *x, double *f, int N){
	int i;	double I = .0;
	for (i = 0; i < N; i++) {I += .5*(x[i+1]-x[i])*(f[i]+f[i+1]);}
	return I;
}
double simpson(double *x, double *f, int N){
	int i;	double I = .0;
	for (i = 1; i < N; i+=2) {I += (1.0/6)*(x[i+1]-x[i-1])*(f[i-1]+4*f[i]+f[i+1]);}
	return I;
}
double integrate(double (*F)(double), double *a, int N) {
	int i;
	double h = (a[1]-a[0])/N;
	double *x = malloc((N+1)*sizeof(double));
	for (i = 0; i < N+1; i++) {x[i] = a[0] + h*i;}
	double *f = malloc((N+1)*sizeof(double));
	for (i = 0; i< N+1; i++) {f[i]=F(x[i]);}
	double I = simpson(x, f, N);
	free(x);
	free(f);
	return I;
}

