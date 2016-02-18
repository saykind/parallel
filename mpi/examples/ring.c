#include	<stdio.h>
#include	<mpi.h>

#define to_right 201


void main (int argc, char *argv[])
{
  int my_rank, size;
  int sent_buf, recv_buf;
  int right, left;
  int sum, i;

  MPI_Status  status;
  MPI_Request request;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  right = (my_rank+1)      % size;
  left  = (my_rank-1+size) % size;
/* ... this SPMD-style neighbor computation with modulo has the same meaning as: */
/* right = my_rank + 1;          */
/* if (right == size) right = 0; */
/* left = my_rank - 1;           */
/* if (left == -1) left = size-1;*/

  sum = 0;
  sent_buf = my_rank;

  for( i = 0; i < size; i++) {
    MPI_Issend(&sent_buf, 1, MPI_INT, right, to_right,
                          MPI_COMM_WORLD, &request);
    MPI_Recv(&recv_buf, 1, MPI_INT, left, to_right,
                        MPI_COMM_WORLD, &status);
    MPI_Wait(&request, &status);
    
    sum += recv_buf;
    sent_buf = recv_buf;
  }
  
  printf ("PE%i:\tSum = %i\n", my_rank, sum);
  MPI_Finalize();
}
