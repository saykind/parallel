#include	<stdio.h>
#include	<mpi.h>

#define msg_tag 111
#define COUNT 2


void main (int argc, char *argv[])
{
  int my_rank, size;
  int right, left;

  int   int_send_buf, int_recv_buf, int_sum, i;
  float float_send_buf, float_recv_buf, float_sum;

  int          array_of_blocklengths[COUNT];
  MPI_Aint     array_of_displacements[COUNT],
               first_var_address,
               second_var_address;
  MPI_Datatype array_of_types[COUNT],
               sendtype, recvtype;

  MPI_Status status;
  MPI_Request request;


  /* Get process and neighbour info. */
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  right = my_rank + 1;
  if (right == size) right = 0;

  left = my_rank - 1;
  if (left == -1) left = size-1;

  /* Set MPI datatypes for sending and receiving partial sums. */

  array_of_blocklengths[0] = 1;
  array_of_blocklengths[1] = 1;

  MPI_Address(&int_send_buf, &first_var_address);
  MPI_Address(&float_send_buf, &second_var_address);

  array_of_displacements[0] = (MPI_Aint) 0;
  array_of_displacements[1] = second_var_address - 
                              first_var_address;

  array_of_types[0] = MPI_INT;
  array_of_types[1] = MPI_FLOAT;

  MPI_Type_struct(COUNT, array_of_blocklengths,
                         array_of_displacements, 
                         array_of_types, 
                         &sendtype);

  MPI_Type_commit(&sendtype);

  MPI_Address(&int_recv_buf, &first_var_address);
  MPI_Address(&float_recv_buf, &second_var_address);

  array_of_displacements[0] = (MPI_Aint) 0;
  array_of_displacements[1] = second_var_address - 
                              first_var_address;

  MPI_Type_struct(COUNT, array_of_blocklengths,
                         array_of_displacements, 
                         array_of_types, 
                         &recvtype);

  MPI_Type_commit(&recvtype);

  /* Compute global sum. */
  int_sum = 0;
  float_sum = 0;
  int_send_buf = my_rank;
  float_send_buf = (float) my_rank;

  for( i = 0; i < size; i++) 
  {
    MPI_Issend(&int_send_buf, 1, sendtype, right, msg_tag,
                          MPI_COMM_WORLD, &request);
    
    MPI_Recv(&int_recv_buf, 1, recvtype, left, msg_tag,
                         MPI_COMM_WORLD, &status);
    
    MPI_Wait(&request, &status);
    
    int_sum += int_recv_buf;
    int_send_buf = int_recv_buf;

    float_sum += float_recv_buf;
    float_send_buf = float_recv_buf;
  }

  printf ("PE%i:\tSum = %i\t%f\n", my_rank, int_sum, float_sum);

  MPI_Finalize();
}
