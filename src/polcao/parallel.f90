module O_Parallel
  
  ! Import necessasry modules
  use O_Kinds

  implicit none

  contains

  ! Notation choice was chosen as a result of the netlib documentation
  ! and can be viewd at: http://netlib.org/scalapack/slug/node76.html
  ! This structure holds the items relevant to BLACS only.
  type BlacsInfo
    ! The BLACS Context we are using for the process grid.
    integer :: context
    
    ! This is the dimensions of the process grid 
    integer :: prows
    integer :: pcols
   
    ! This is the location in the process grid for a particular process
    integer :: myprow
    integer :: mypcol

    ! MPI related information
    integer :: mpirank
    integer :: mpisize
  end type BlacsInfo

  ! This structure holds information about the Global and Local arrays needed
  ! by a process to complete some distriubted operation.
  ! For every matrix being used in some distributed communication should have
  ! a BlacsInfo structure assosciated with it.
  type ArrayInfo

    ! This is the BLACS array descriptor assosciated with some array
    integer, dimension(9) :: desc

    ! These are the size of the dimensions of the global array
    integer :: I
    integer :: J

    ! Togethere mb and nb define the block size used for the block-cyclic
    ! distribution.
    ! This is the blocking factor used to distribute the rows of the array
    integer :: mb
    ! This is the blocking factor used to distribute the columns of the array
    integer :: nb
    
    ! Holds the extra dimensions for a process in the event that the matrix
    ! is not evenly distributable
    integer :: extraRows
    integer :: extraCols

    ! The local array that will contain the distributed data
#ifndef GAMMA
    complex (kind=double), allocatable, dimension(:,:) :: local
#else
    real (kind=double), allocatable, dimension(:,:) :: local
#endif


  end type Array Info

! Placeholder routine containing some notes on what needs to happen before
! we do blacs/pblas/scalapack stuff
subroutine asdfasdfasdf()
  ! First we need to know how many processes are participating in the 
  ! calculation.

  ! Then we need to create a communicator handle for just the participating
  ! processes, incase they need to barrier, or do any other "world" comms
  ! while there are excluded processes.
  
  ! Then we need to create a condition that excludes processes who are not
  ! participating in the communication.
  if (myinfo%mpirank < x ) then
  endif
  call MPI_Barrier(MPI_COMM_WORLD, mpierr)
end subroutine asdfasdfasdf

! This subroutine calculates the process grid dimensions and sets up the 
! blacs context for a distributed operation
subroutine setupBlacs(blcsinfo)
  implicit none

  ! Define passed parameters
  type(BlacsInfo), intent(out) :: blcsinfo

  call MPI_Comm_rank(MPI_COMM_WORLD, blcsinfo%mpirank, mpierr)
  call MPI_Comm_size(MPI_COMM_WORLD, blcsinfo%mpisize, mpierr)

  ! First we need to calculate the processor grid
  calcProcgrid(blcsinfo)

  call BLACS_GET(-1,0, blcsinfo%ctxt)
  call BLACS_GRIDINIT(blcsinfo%ctxt, 'r', blcsinfo%prows, blcsinfo%pcols)
  call BLACS_GRIDINFO(blcsinfo%ctxt, blcsinfo%prows, blcsinfo%pcols, &
    & blcsinfo%myprow, blcsinfo%mypcol)

end subroutine setubBlacs

subroutine calcProcGrid(blcsinfo)
  implicit none

  ! Define passed parameters
  type(BlacsInfo), intent(inout) :: blcsinfo

  blcsinfo%prows = int(sqrt(dble(blcsinfo%mpisize)))
  blcsinfo%pcols = blcsinfo%mpisize/prow
end subroutine calcProcGrid

subroutine setupArrayDesc(arrinfo, blcsinfo, numGlobalRows, numGlobalCols)
  implicit none

  ! Define Passed Parameters
  type(BlacsInfo), intent(in) :: blcsinfo
  type(ArrayInfo), intent(inout) :: arrinfo

  arrinfo%I = numGlobalRows
  arrinfo%J = numGlobalCols

  ! First we need to calculate the blocking factors, mb, nb, and if there are
  ! any extra dimensions
  call getBlockDims(arrinfo)

  ! Now that we have our block dims, we have sufficient information to allocate
  ! the local array assosciated with each process.
  call allocLocalArray(arrinfo, blcsinfo)

  ! With the local array allocated we can now setup our array descriptor
  call getArrDesc(arrinfo, blcsinfo)

end subroutine setupArrayDesc

! This subroutine determines our block size needed to map our process grid onto
! our global matrix
subroutine getBlockDims(arrinfo, blcsinfo)
  implicit none

  ! Define passed parameters
  type(BlacsInfo), intent(inout) :: blcsinfo
  type(ArrayInfo), intent(inout) :: arrinfo

  ! Local variables
  integer :: extraRow, extraCol

  ! Calculate the blocking factors to be used. From the calcProcGrid
  ! routine we know that procGridRows will be <= procGridCols.
  ! Because of this we'll do the integer divide on procGridCols, ensuring
  ! that our blocks are square.
  arrinfo%mb =  int(arrinfo%I / blcsinfo%pcols)
  arrinfo%nb =  int(arrinfo%J / blcsinfo%pcols)

  ! For the case that global matrix doesn't divide perfectly by our choice in
  ! block size, there will be a couple extra rows and columns that need to be 
  ! asssigned to a process.
  ! The total dim I divided by the block dimension, tells us how many blocks
  ! are going to be along a dimensions. WHen we mod this by the rows in the 
  ! process grid, it tells us the row  coordinate in the process grid the 
  ! extra rows belong to. This of course also works for the columns.
  extraRow = mod( arrinfo%I/arrinfo%mb, blcsinfo%prows)
  extraCol = mod( arrinfo%J/arrinfo%nb, blcsinfo%pcols)

  if ( extraRow == myinfo%prow ) then
    myinfo%extraBlockRows = mod(arrinfo%I, arrinfo%nb)
  else
    myinfo%extraBlockRows = 0
  endif
  if ( extraCol == myinfo%pcol ) then
    myinfo%extraBlockCols = mod(arrinfo%J, arrinfo%nb)
  else
    myinfo$extraBlockCols = 0
  endif
end subroutine getblockDims

! This subroutine allocates the proper amount of space needed for each
! processes local array.  The SCALAPACK function NUMROC (number of rows and
! colums) is meant for this purpose.
subroutine allocLocalArray(arrinfo, blcsinfo)
  use O_Kinds

  implicit none

  ! Define passed parameters
  type(BlacsInfo), intent(inout) :: blcsinfo
  type(ArrayInfo), intent(inout) :: arrinfo

  ! Define local variables
  integer :: nlrows, nlcols

  ! External Functions
  integer, external :: numroc

  nlrows = numroc(arrinfo%I, arrinfo%mb, blcsinfo%myrow, 0, blcsinfo%prows)
  nlcols = numroc(arrinfo%J, arrinfo%nb, blcsinfo%mycol, 0, blcsinfo%pcols)

  allocate(arrinfo%local(nlrows,nlcols))

  ! initialize the local array to zero
#ifndef GAMMA
  arrinfo%local(:,:) = cmplx(0.0_double, 0.0_double)
#else
  arrinfo%local(:,:) = 0.0_double
#endif

end subroutine allocLocalArray

subroutine getArrDesc(arrinfo, blcsinfo)
  implicit none
  
  ! Define passed parameters
  type(BlacsInfo), intent(inout) :: blcsinfo
  type(ArrayInfo), intent(inout) :: arrinfo

  ! Define local variables
  integer :: info

  call descinit(arrinfo%desc, arrinfo%I, arrinfo%J, arrinfo%mb, arrinfo%nb, &
    & 0, 0, blcsinfo%context, size(arrinfo%local,1), info)

  if (info) then
    print *, "BLACS Routine DESCINIT did not exit successfully. Exiting"
    exit
  endif

end subroutine getArrDesc

end module O_Parallel
