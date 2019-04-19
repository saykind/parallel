#include <stdio.h>
#include <stdlib.h>
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
void plot(FILE *gpp, double *x, double *f, int N) {
	int i;
	fprintf(gpp, "plot '-' w points pt 1 ps .5\n");
	for (i = 0; i < N; i++)
		fprintf(gpp, "%.6f %.6f\n", x[i], f[i]);
	fprintf(gpp,"e\n");
	fflush(gpp);
	return;
}
