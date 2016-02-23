module O_Parallel
  implicit none

  contains

  type BlacsInfo
    ! Holds the coordinates for a process in the Blacs process grid
    integer :: prow
    integer :: pcol

    ! Holds local block size information

  end type BlacsInfo

! This subroutine calculates the Global Array indices from local indices
subroutine getGlobalIndices(    )
  implicit none

end subroutine getGlobalIndices

! This subroutine calculates the Local indices from the global array indices
subroutine getLocalIndices(  )
  implicit none

end subroutine getLocalIndices(  )

! This subroutine determines the process grid to be used for the
! for a matrix of size numDim. There are a lot of considerations that need to be
! taken into account in order to have an optimal process grid.
!
! From netlib documentation the recomendation for routines that use 
! Cholesky, is to have near square blocks and process grids. However, the
! ability to do that depends on the number of processes specified by the
! user. Because there's little we can do about that, we'll try to get the
! "most square" processor grid given the number of processses. But in
! the documentation we'll recommend that the user use perfect square number
! of processes.
!
! Because we're working with Hermitian matrices to solve our eigen problem
! we can typically get very square blocking factors. This should be little
! issue to map onto any processor grid.
!
! We also use this in Setup for the orthogonalization process which may have to
! be reviewed later incase a more optimal type of blocking is used
subroutine getBlockFact(myInfo, numprocs, gridDim, numDim, blkFact)
  implicit none

  ! Define passed parameters
  integer, dimension(3), intent(in) :: myInfo
  integer, intent(in) :: numprocs
  integer, dimension(2), intent(out) :: gridDim
  integer, intent(in) :: numDim
  integer, dimension(2), intent(out) :: blkFact ! blocking factors

  ! Local variables
  integer :: modifiedValeDim

  if (mod(numDim,2) /= 0) then
    modifiedValeDim = numDim - 1
  else
    modifiedValeDim = numDim
  endif

  call calcProcGrid(numProcs, gridDim(1), gridDim(2))

  ! Calculate the blocking factors to be used. From the calcProcGrid
  ! routine we know that gridDim1 (prows) will be <= gridDim2 (pcols).
  ! Because of this we'll do the integer divide on gridDim(2), ensuring
  ! that our blocks are square.
  blkFact(1) =  int(modifiedValeDim / gridDim(2))
  blkFact(2) =  int(modifiedValeDim / gridDim(2))
end subroutine getblkFact

end module O_Parallel
! In order to do the getBlockFact subroutine we need to do integer factorization
! of the number of processes, in order to find the "most square" 
! distribution of processors into the grid. To do this we take the, integer
! square root of numprocs, then divide numprocs by prows. From testing 
! integer values of processes in the range of 2 to 1,000,000 the biggest
! difference in abs(prows-pcols) is 2. Meaning at most only 2 processes
! will be left out by the user's choice of processes.
subroutine calcProcGrid(numprocs, prows, pcols)
  implicit none

  ! Define passed parameters
  integer, intent(in) :: numprocs
  integer, intent(out) :: prows, pcols

  prows = int(sqrt(dble(numprocs)))
  pcols = numprocs/prows
end subroutine calcProcGrid
