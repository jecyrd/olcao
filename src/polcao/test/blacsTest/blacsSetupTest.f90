program blacsTest
  use HDF5
  use MPI
  use O_Parallel
  
  implicit none

  integer :: mpirank
  integer :: mpisize
  integer :: mpierr
  integer, dimension(MPI_STATUS_SIZE) :: mpistatus

  type(BlacsInfo) :: blcsinfo
  type(ArrayInfo) :: arrinfo

  integer :: valeDim, coreDim, numkp

  integer :: i,j

  integer :: numargs
  character(len=10) :: arg

  numargs = command_argument_count()

  call get_command_argument(1,arg)
  read (arg,'(I10)') valedim
  call get_command_argument(2,arg)
  read (arg,'(I10)') coredim

  call MPI_INIT(mpierr)
  
  numkp = 1

  call setupBlacs(blcsinfo)

  call setupArrayDesc(arrinfo, blcsinfo, coredim, valedim, numkp)

  print *, blcsinfo%mpirank, blcsinfo%mpisize,blcsinfo%prows,blcsinfo%pcols, &
         & blcsinfo%myprow,blcsinfo%mypcol,arrinfo%I,arrinfo%J,arrinfo%numKP, &
         & arrinfo%mb,arrinfo%nb,arrinfo%nrblocks,arrinfo%ncblocks, &
         & arrinfo%extraRows,arrinfo%extraCols

  call MPI_FINALIZE(MPIERR)

end program blacsTest 

