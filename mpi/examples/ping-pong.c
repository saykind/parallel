#include <stdio.h>
#include <mpi.h>

#define proc_A 0
#define proc_B 1
#define ping 100
#define pong 101
#define number_of_messages 50
#define length_of_message   1

void main(int argc, char *argv[])
{
  int my_rank;
  float buffer[length_of_message];
  int i;   
  double start, finish, time;

  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if (my_rank == proc_A) 
  {
    start = MPI_Wtime();
    for (i = 1; i <= number_of_messages; i++)
    {
      MPI_Ssend(buffer, length_of_message, MPI_FLOAT, proc_B, ping,
                        MPI_COMM_WORLD);
      MPI_Recv(buffer, length_of_message, MPI_FLOAT, proc_B, pong,
                       MPI_COMM_WORLD, &status);
    }
    finish = MPI_Wtime();
    time = finish - start;
    printf("Time for one messsage: %f seconds.\n", 
           (float)(time / (2 * number_of_messages)));
  }
  else if (my_rank == proc_B) 
  {
    for (i = 1; i <= number_of_messages; i++)
    { 
      MPI_Recv(buffer, length_of_message, MPI_FLOAT, proc_A, ping,
                       MPI_COMM_WORLD, &status);
      MPI_Ssend(buffer, length_of_message, MPI_FLOAT, proc_A, pong,
                        MPI_COMM_WORLD);
    }
  } 
  MPI_Finalize();
}


