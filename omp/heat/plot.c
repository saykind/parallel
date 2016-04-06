#include <stdio.h>
#include <stdlib.h>
FILE *gpinit(void) {
	FILE *gpp = popen("gnuplot", "w");
	if (gpp == NULL) {return NULL;}
	fprintf(gpp, "set term x11\n");
	fprintf(gpp, "unset key\n");
	fprintf(gpp, "unset border\n");
	fprintf(gpp, "set grid\n");
	fprintf(gpp, "set pm3d at bs interpolate 0,0\n");
	fprintf(gpp, "set xrange [ 0 : 1 ]\n");
	fprintf(gpp, "set yrange [ 0 : 1 ]\n");
	fprintf(gpp, "set zrange [ 0 : 1 ]\n");
	fprintf(gpp, "set cbrange [ 0 : 1 ]\n");
	return gpp;
}
void plot(FILE *gpp, double **a,  int N) {
	int i,j;
	fprintf(gpp, "splot '-' w pm3d \n");
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			fprintf(gpp, "%.6f %.6f %.6f\n", 1.0*i/N, 1.0*j/N, a[i][j]);
		}
		fprintf(gpp, "\n");
	}
	fprintf(gpp,"e\n");
	fflush(gpp);
	//usleep(200);
	return;
}
