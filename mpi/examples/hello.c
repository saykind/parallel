#include <stdio.h>
#include <mpi.h>


void main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);

   printf ("Hello world!\n");

   MPI_Finalize();
}
