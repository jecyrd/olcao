program nodeTest
  use MPI

#include "mafdecls.fh"
#include "global.fh"

  integer :: mpierr, mpiRank, gaRank, gastat

  call MPI_INIT(mpierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mpiRank, mpierr)

  call ga_initialize()

  gastat = MA_init(MT_F_DCPL, 3000000,300000)

  gaRank = ga_nodeid()

  print *, gaRank, mpiRank




end program
