program blacsTest
  use HDF5
  use MPI
  use O_Parallel
  use O_ParallelSetup
  use O_bstAtomPair
  
  implicit none

  integer :: mpirank
  integer :: mpisize
  integer :: mpierr
  integer, dimension(MPI_STATUS_SIZE) :: mpistatus

  type(BlacsInfo) :: blcsinfo
  type(ArrayInfo) :: arrinfo
  type(AtomPair), pointer :: atomPairs
  type(bst_atom_pair_node), pointer :: atomTree

  integer :: i,j
  integer :: valeDim, coreDim, numkp
  integer :: numargs
  integer :: whichArr
  character(len=10) :: arg

  call MPI_INIT(mpierr)

  numargs = command_argument_count()
  call get_command_argument(1,arg)
  read (arg,'(I10)') valedim
  call get_command_argument(2,arg)
  read (arg,'(I10)') coredim
  numkp = 1

  




  call MPI_FINALIZE(MPIERR)
end program blacsTest 

