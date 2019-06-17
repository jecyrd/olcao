!! This module is focused on the parallelization of OLCAOmain. 
!!
!! Author: James E. Currie
!! Email: jecyrd@mail.umkc.edu
!!
module O_ParallelMain

  ! Although each subroutine should specify implicit none, it is placed
  ! in this module to protect against accidentaly omission
  implicit none
#ifndef GAMMA
   ! This structure holds information about the Global and Local arrays needed
   ! by a process to complete some distributed operation.
   type mArrayInfo

      ! This is the BLACS array descriptor assosciated with some array
      integer, dimension(9) :: desc

      ! These are the size of the dimensions of the global array
      integer :: I
      integer :: J
      integer :: numKP

      ! Together mb and nb define the block size used for the block-cyclic
      ! distribution.
      ! This is the blocking factor used to distribute the rows of the array
      integer :: mb
      ! This is the blocking factor used to distribute the columns of the array
      integer :: nb

      ! This is the number of blocks along the rows that the process has
      integer :: nrblocks
      ! This is the number of blocks along the columns that the process has
      integer :: ncblocks

      ! Holds the extra rows/columns in the respective dimensions for a process
      ! in the event that the matrix is not evenly distributable.
      integer :: extraRows
      integer :: extraCols

      ! The local array that will contain the distributed data
      complex (kind=double), allocatable, dimension(:,:,:) :: local
   end type mArrayInfo
#else
   ! Same as above but only used for the gamma case
   type mArrayInfo

      ! This is the BLACS array descriptor assosciated with some array
      integer, dimension(9) :: desc

      ! These are the size of the dimensions of the global array
      integer :: I
      integer :: J
      integer :: numKP

      ! Together mb and nb define the block size used for the block-cyclic
      ! distribution.
      ! This is the blocking factor used to distribute the rows of the array
      integer :: mb
      ! This is the blocking factor used to distribute the columns of the array
      integer :: nb

      ! This is the number of blocks along the rows that the process has
      integer :: nrblocks
      ! This is the number of blocks along the columns that the process has
      integer :: ncblocks

      ! Holds the extra rows/columns in the respective dimensions for a process
      ! in the event that the matrix is not evenly distributable.
      integer :: extraRows
      integer :: extraCols

      ! The local array that will contain the distributed data
      real (kind=double), allocatable, dimension(:,:) :: local
   end type mArrayInfo
#endif


  contains

! This subroutine takes a given block for a process, and calculates the
! mapping between the global and local proceses.
!
! For some process grid determined by the getBlockFact subroutine. We need
! to extract data from the global array in such away to fit our grid
! using the 2d block cyclic mapping specified in the SCALAPACK
! documentation on netlib. http://netlib.org/scalapack/slug/node76.html
subroutine globalToLocalMap(blockRow, blockCol, nbrow, nbcol, gridDim, &
    & myInfo, blkFact, trailDim, glo, ghi, llo, lhi, ld)

  implicit none

  ! Define passed parameters
  integer, intent(in) :: blockRow, blockCol ! Block row and col index to get
  integer, intent(in) :: nbrow, nbcol ! Total Number of blocks
  integer, dimension(2), intent(in) :: gridDim ! Process Grid Dims 
  integer, dimension(3), intent(in) :: myInfo  ! prow, pcol, mpirank
  integer, dimension(2), intent(in) :: blkFact ! Blocking Factors
  integer, dimension(2), intent(in) :: trailDim ! trailing dims
  integer, dimension(2), intent(out) :: glo, ghi, llo, lhi
  integer, dimension(2), intent(out) :: ld

  ! Define local Variables
  integer :: l,m

  ! Using the equations in the documentation linked above, we can obtain
  ! the following set of equations.
  !  a=l*MB+x  I=(l*Pr + pr)*MB+x
  !  b=m*NB+y  J=(m*Pc + pc)*NB+y
  ! Where the coordinates in the local array are llo(1),llo(2), 
  ! coordinates in global array are (I,J), coordinates of local block are 
  ! (l,m), coordinates inside local block are (x,y), and process 
  ! coordinates (pr,pc). For block size MBxNB, and process grid size PrxPc.
  ! With the equation for I we can solve it twice with x=0, and 
  ! x=MB-1, this will give us the upper left and lower right of the block.
  ! This will be sufficient information to pull the block down from
  ! the global array.
  ! 
  ! The upper left coordinate will be held by lo, and the lower right held
  ! by hi. "+ 0" intentionally left in for clarity.
  !
  ! Lo dimensions need to be incremented by 1, because the math is based
  ! off of a 0 index. Hi should calculate correctly.
 
  ! l and m are a relic of a rewrite, left for clarity.
  l = blockRow-1
  m = blockCol-1

  ! Upper left corner.
  glo(1) = ((l*gridDim(1) + myInfo(1))*blkFact(1) + 0) + 1
  glo(2) = ((m*gridDim(2) + myInfo(2))*blkFact(2) + 0) + 1
  
  ! Starting positions of local lo should be unnefected for every case
  ! Read next comment.
  llo(1) = (l*blkFact(1) + 0) + 1
  llo(2) = (m*blkFact(2) + 0) + 1

  ! If we are an irregular sized block we have to do a special get.
  ! The special get has different ending indices than a normal block
  ! of the exact specified dimensions. 
  ! Should try to think about reordering conditions in terms of
  ! frequency.
  !
  ! The first  condition is: If trailing cols and rows. It'll make
  !                          A block smaller in both dimensions.
  ! The second condition is: If trailing cols.
  ! The third  condition is: If trailing rows.
  ! The last   condition is: If normal full block
  if ( sum(trailDim(:))==2 .and. (blockRow==nbrow) &
    &                                        .and. (blockCol==nbcol) ) then
    ! Lower right corner. 
    ghi(1) = (l*gridDim(1) + myInfo(1))*blkFact(1) + trailDim(1)
    ghi(2) = (m*gridDim(2) + myInfo(2))*blkFact(2) + trailDim(2)
    ! Calculate local ending coordinates
    lhi(1) = llo(1) + trailDim(1) 
    lhi(2) = llo(2) + trailDim(2) 

  else if ( (trailDim(2) > 0) .and. (blockCol==nbcol) ) then ! trailing cols
    ! Lower right corner. 
    ghi(1) = (l*gridDim(1) + myInfo(1))*blkFact(1) + blkFact(1)
    ghi(2) = (m*gridDim(2) + myInfo(2))*blkFact(2) + trailDim(2)
    ! Calculate local ending coordinates
    lhi(1) = llo(1) + blkFact(1)
    lhi(2) = llo(2) + trailDim(2)

  else if ( (trailDim(1) > 0) .and. (blockRow==nbrow) ) then ! trailing rows
    ! Lower right corner. 
    ghi(1) = (l*gridDim(1) + myInfo(1))*blkFact(1) + trailDim(1)
    ghi(2) = (m*gridDim(2) + myInfo(2))*blkFact(2) + blkFact(2)
    ! Calculate local ending coordinates
    lhi(1) = llo(1) + trailDim(1) 
    lhi(2) = llo(2) + blkFact(1)  
                           
  else ! Normal Block
    ! Lower right corner. 
    ghi(1) = (l*gridDim(1) + myInfo(1))*blkFact(1) + blkFact(1)
    ghi(2) = (m*gridDim(2) + myInfo(2))*blkFact(2) + blkFact(2)
    ! Calculate local ending coordinates
    lhi(1) = llo(1) + blkFact(1)
    lhi(2) = llo(2) + blkFact(1)

  endif
  ! Set ld
  ld(1) = (lhi(1)-llo(1)) + 1
  ld(2) = (lhi(2)-llo(2)) + 1
end subroutine globalToLocalMap 

! This subroutine takes the computed wave function and write it to HDF5.
! Because parallel writes on NFS are slow, we'll have process 0 store it's
! sections and deallocate it's local array. Then we'll do a series of 
! blocking MPI calls so that process zero can store it's information.
subroutine writeWaveFunction()
  implicit none

  ! Define passed parameters

  ! Define local variables
  integer :: mpisize, mpierr
  integer :: i ! loop variable

  ! Have process zero store it's sections of the wavefunction.
  ! If not process zero, put up barrier and wait till zero is done.
  if (myInfo(3) == 0) then
    call writeMatrixSection()
    call MPI_BARRIER(mpierr)
  else
    call MPI_BARRIER(mpierr)
  endif

  ! Loop over the world (size-1) and pass information to process zero to
  ! store.
  do i=1, (mpisize-1)
    ! if you're process zero receive local array section then write
    if (myInfo(3) == 0) then
      call recvLocalArray()

      call writeMatrixSection()
      
    ! If you're process mpirank==i, send local array section
    elseif (myInfo(3) == i) then
      call sendLocalArray()

    endif
    ! Else if you're anyone else just hangout and wait for the other two
    ! to catch up.
    call MPI_BARRIER(mpierr)
  enddo

end subroutine writeWaveFunction

! Subroutine to write a matrix section to HDF5. Should only be called by
! process zero.
subroutine writeMatrixSection(localArray, rankToWrite)
  implicit none

  ! Define passed parameters
  complex, (kind=double) :: localArray
  integer :: rankToWrite ! Rank that localArray belongs to. 


end subroutine writeMatrixSection

! Subroutine to handle receiving local arrays
subroutine recvLocalArray()
  use MPI

  implicit none

  call MPI_RECV(buffer, count, datatype, source, tag, comm, status, ierror)
  
end subroutine recvLocalArray

! Subroutine to handle sending local array
subroutine sendLocalArray()
  use MPI

  implicit none

  call MPI_SEND(buffer, count, datatype, dest, tag, comm, ierror)
end subroutine sendLocalArray


! This subroutine is responsible for making sure the eigenvalues are
! redistributed to all processes for use later
subroutine distributeEigenVals()
!  call MPI_BROADCAST

end subroutine

! This subroutine cleans up scalapack descriptors and contexts
subroutine cleanUpSl(slctxt, desca, descb, descz)
  implicit none

end subroutine cleanUpSl

end module O_ParallelMain
