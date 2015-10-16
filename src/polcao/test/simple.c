#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv)
{
  unsigned int rank;
  unsigned int size;
  int mpierr;

  mpierr = MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  printf("I'm %u out of %u\n", rank, size);

  MPI_Finalize();
  return 0;
}
