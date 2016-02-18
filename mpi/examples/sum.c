#include <stdio.h>
#include <stdlib.h>
#include <mpi/mpi.h>

int main(int argc,char *argv[]) {
	MPI_Init(&argc, &argv);
	
	int size, rank, j;
	unsigned long int i=0, n=0, m=0, N=0, S = 0, last=0;
	double start, end;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//printf("[%d]:\tsize = %d\n", rank, size);

	if (argc < 2) {MPI_Finalize(); return 0;}
	N = atoi(argv[1]); 

 	if (rank == 0) {start = MPI_Wtime();}
	
	n = N/size;
	last = (rank+1)*n;
	if (rank == (size-1)) {last = N+1;}

	for (i = (rank)*n; i < last; i++) {S += i;}

	printf("[%d]:\tsum = %lu\n", rank, S);	

	
	
	if (rank != 0) {
		MPI_Send((void *) &S, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
	}
	if (rank == 0) {
		MPI_Status st;
		unsigned long int SUM = S;
		for (j =1 ; j < size; j++) {
			MPI_Recv((void *) &m, 1, MPI_INT, j, 1, MPI_COMM_WORLD, &st);
			printf("[%d]:\tRecived from %d: %lu\n", rank,j, m);
			SUM += m; 
		}
		printf("\n[%d]:\tSUM = %lu\n", rank, SUM);

		end = MPI_Wtime();
		usleep(1000);
		printf("[%d]:\ttime = %lf millisec\n\n", rank, 1e3*(end-start));
		
	}
 
	

	MPI_Finalize();
	return 0;
}
