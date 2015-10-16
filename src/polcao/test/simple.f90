program gaBlockCyclic

  use MPI

  implicit none

!#include "mafdecls.fh"
!#include "global.fh"

  integer :: mpiRank, mpiSize, mpierr

  mpiRank = 2
  mpiSize = 0
  mpierr = 0

  call mpi_init(mpierr)
  print *, "Init: ", mpierr
  call MPI_COMM_RANK(MPI_COMM_WORLD, mpiRank, mpierr)
  print *, "Rank: ", mpierr
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpiSize, mpierr)
  print *, "Size: ", mpierr


  print *, "\n\nStuff: ", mpiRank, mpiSize

  call mpi_finalize(mpierr)
end program
