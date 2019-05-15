program blacsTest
  use HDF5
  use MPI
  use O_ParallelSetup

  use O_pOLCAOhdf5
  
  implicit none

  integer :: mpirank
  integer :: mpisize
  integer :: mpierr
  integer, dimension(MPI_STATUS_SIZE) :: mpistatus

  integer, dimension(2) :: pgrid
  type(AtomPair), pointer :: atomList

  integer :: i
  integer :: hdferr


  call MPI_INIT(mpierr)

  call h5open_f(hdferr)
  if (hdferr <0) stop 'Failed to open the HDF library'

  call MPI_COMM_RANK(MPI_COMM_WORLD, mpirank, mpierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpisize, mpierr)


  do i=0,(mpisize-1)
    if (i==mpirank) then
      call accessAtomPairsHDF5(atomList, pgrid)
      print *, 'Process: ', mpirank, 'Has Grid: ', pgrid
      print *, 'Process: ', mpirank, 'Has AtomPairs:'

      do
        print *, '    ', atomList%I, atomList%J
        if (.not. associated(atomList%next)) exit
        atomList => atomList%next
      enddo

    endif
    call MPI_Barrier(MPI_COMM_WORLD,mpierr)
  enddo

  call h5close_f(hdferr)
  call MPI_FINALIZE(MPIERR)
end program blacsTest 

