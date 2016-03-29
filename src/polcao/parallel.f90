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
    integer :: numKP

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
    complex (kind=double), allocatable, dimension(:,:,:) :: local
#else
    real (kind=double), allocatable, dimension(:,:) :: local
#endif
  end type Array Info

! This subroutine calculates the process grid dimensions and sets up the 
! blacs context for a distributed operation
subroutine setupBlacs(blcsinfo)
  implicit none

  ! Define passed parameters
  type(BlacsInfo), intent(out) :: blcsinfo

  call MPI_Comm_rank(MPI_COMM_WORLD, blcsinfo%mpirank, mpierr)
  call MPI_Comm_size(MPI_COMM_WORLD, blcsinfo%mpisize, mpierr)

  ! First we need to calculate the processor grid
  calcProcGrid(blcsinfo)

  call BLACS_GET(-1,0, blcsinfo%ctxt)
  call BLACS_GRIDINIT(blcsinfo%ctxt, 'r', blcsinfo%prows, blcsinfo%pcols)
  call BLACS_GRIDINFO(blcsinfo%ctxt, blcsinfo%prows, blcsinfo%pcols, &
    & blcsinfo%myprow, blcsinfo%mypcol)

end subroutine setubBlacs

! This subroutine calculates the dimensions of the processor grid given 
! the number of available processes
subroutine calcProcGrid(blcsinfo)
  implicit none

  ! Define passed parameters
  type(BlacsInfo), intent(inout) :: blcsinfo

  blcsinfo%prows = int(sqrt(dble(blcsinfo%mpisize)))
  blcsinfo%pcols = blcsinfo%mpisize/prow
end subroutine calcProcGrid

subroutine setupArrayDesc(arrinfo, blcsinfo, numGlobalRows, numGlobalCols, &
    & numKP)
  implicit none

  ! Define Passed Parameters
  type(BlacsInfo), intent(in) :: blcsinfo
  type(ArrayInfo), intent(inout) :: arrinfo
  integer, intent(in) :: numGlobalRows, numGlobalCols
  integer, intent(in) :: numKP

  arrinfo%I = numGlobalRows
  arrinfo%J = numGlobalCols
  arrinfo%numKP = numKP

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
  if (arrinfo%I == arrinfo%J)
    arrinfo%mb =  int(arrinfo%I / blcsinfo%pcols)
    arrinfo%nb =  int(arrinfo%J / blcsinfo%pcols)
  else
    arrinfo%mb =  int(arrinfo%I / blcsinfo%prows)
    arrinfo%nb =  int(arrinfo%J / blcsinfo%pcols)
  endif

  ! For the case that global matrix doesn't divide perfectly by our choice in
  ! block size, there will be a couple extra rows and columns that need to be 
  ! asssigned to a process.
  ! The total dim I divided by the block dimension, tells us how many blocks
  ! are going to be along a dimensions. When we mod this by the rows in the 
  ! process grid, it tells us the row  coordinate in the process grid the 
  ! extra rows belong to. This of course also works for the columns.
  extraRow = mod( arrinfo%I/arrinfo%mb, blcsinfo%prows)
  extraCol = mod( arrinfo%J/arrinfo%nb, blcsinfo%pcols)

  if ( extraRow == arrinfo%myprow ) then
    myinfo%extraBlockRows = mod(arrinfo%I, arrinfo%nb)
  else
    myinfo%extraBlockRows = 0
  endif
  if ( extraCol == arrinfo%mypcol ) then
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

#ifndef GAMMA
  allocate(arrinfo%local(nlrows,nlcols,arrInfo%numKP))
#else
  allocate(arrinfo%local(nlrows,nlcols))
#endif

  ! initialize the local array to zero
#ifndef GAMMA
  arrinfo%local(:,:) = cmplx(0.0_double, 0.0_double)
#else
  arrinfo%local(:,:) = 0.0_double
#endif
end subroutine allocLocalArray

! This subroutine initializes the BLACS array descriptor and 
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

! This subroutine calculates the global array coordinates for a local block
! and will aid us in both Setup and Main when trying to populate arrays.
!
! Please refer to http://netlib.org/scalapack/slug/node76.html
! For detailed information relating to this subroutine.
!
! Using the documentation above we have the follow equations.
! a=l*MB+x           b=m*NB+y   
! l = (I-1)/(Pr*MB)  m = (J-1)/(Pr*NB)
! I=(l*Pr + pr)+MB+x J=(m*Pc + pc)+MB+y
! Where (a,b)   are the coordinates with the local array.
!       (x,y)   are the coordinates within a local block in the local array.
!       (l,m)   are the coordinates of a block within the local array.
!       (I,J)   are the coordinates in the global array.
!       (Pr,Pc) are the size of the process grid.
!       (pr,pc) are the coordinates of a process in the grid.
!
! What we need to calculate is 2 sets of coordinates, first the global 
! coordinates of the "top left" element of a local block. Second the global
! coordinates of the "bottom right" element of a local block. This'll allow us
! to tell which elements of the global array we require to fill a local block.
!
! Using the above equations we can obtain:
!  I = (a-x)*Pr + 1    J = (b-y)*Pc + 1
! By setting x and y to 0. We can obtain the "top left" global coordinate. Then
! by setting x and y to (MB-1) and (NB-1) respectively, we obtain the
! "bottom right" coordinate of the global array.
!
! For this to work the passed parameters a,b should be the coordinate in
! the local array corresponding to the "top left" element of some local block.
subroutine localToGlobalMap(a, b, glo, ghi, arrinfo, blcsinfo)
  implicit none

  ! Define passed parameters
  type(BlacsInfo), intent(in) :: blcsinfo
  type(ArrayInfo), intent(in) :: arrinfo
  integer, intent(in) :: a,b
  integer, dimension(2), intent(out) :: glo, ghi ! target global values

  ! Calculate the upper left corner
  glo(1) = (a-0)*blcsinfo%prows + 1
  glo(2) = (b-0)*blcsinfo%pcols + 1

  ! Calculate the bottom left corner
  ghi(1) = (a-arrinfo%mb)*blcsinfo%prows + 1
  ghi(2) = (a-arrinfo%nb)*blcsinfo%pcols + 1

end subroutine globalToLocalMap

end module O_Parallel

