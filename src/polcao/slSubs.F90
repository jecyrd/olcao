!! This module and its subroutines were modeled after various scalapack 
!! interfaces in the Global Arrays, and NWChem packages from 
!! Pacific Northwest National Lab (PNNL). The interfaces in GA and NWChem
!! where developed to maintain F77 compatability. But because OLCAO
!! is nearly entirely in F90 or later, the subroutines may look vastly
!! different from the referenced subroutines. Mostly because the 
!! Memory Allocator (MA) library is a relic and the absense of dynamically
!! allocatable memory in F77.
!!
!! The currently implementation is incomplete due to the intial naive
!! calculation of the BLACS process grid, and the maping to such grid.
!! It is the intention to have a working naive implementation for in place
!! The oga_pdzhegvx subroutine should be mostly immune to changes in the
!! oga_tosl and gridCalc subroutines. The oga_tosl and oga_pdzhegvx
!! subroutines will be edited later to have a more optimal calculation of
!! the BLACS parameters required for the specific calculation.
!!
!! Author: James E. Currie
!! Email: jecyrd@mail.umkc.edu
!!
module O_SlSubs

  ! Although each subroutine should specify implicit none, it is placed
  ! in this module to protect against accidentaly omission
  implicit none

  contains

! This subroutine is responsible for the conversion of a Global Array
! to a local array in preparation for the pzhegvx lapack call in the
! subroutine below. 
! For some process grid determined by the gridCalc subroutine. We need
! to extract data from the global array in such away to fit our grid
! using the 2d block cyclic mapping specified in the SCALAPACK
! documentation on netlib. http://netlib.org/scalapack/slug/node76.html
subroutine oga_tosl(ga_arr, local, valeDim, gridDim, numprocs, myinfo, blkFact, trailDim)
  use MPI
  use O_Kinds

  implicit none

#include "mafdecls.fh"
#include "global.fh"

  ! Define passed parameters
  integer, intent(inout) :: ga_arr ! Global array to get data from
  complex (kind=double), intent(inout), allocatable, dimension(:,:) :: local
  integer, intent(in) :: valeDim
  integer, dimension(2), intent(in) :: gridDim ! Process Grid Dims 
  integer, intent(in) :: numprocs
  integer, dimension(3), intent(in) :: myinfo  ! prow, pcol, mpirank
  integer, dimension(2), intent(in) :: blkFact ! Blocking Factors
  integer, dimension(2), intent(in) :: trailDim

  ! Define local Variables
  integer :: l,m
  integer :: i,j
  integer :: a,b,aend,bend, ld
  integer, dimension(2) :: lo, hi
  integer, dimension(3) :: gaInfo
  integer :: nrows, ncols ! Number of local array rows and columns
  integer :: nbrow, nbcol ! Number of blocks along dimensions

  ! External subroutine
  integer, external ::  numroc

  ! Remove later
  integer, dimension(:), allocatable :: lTrack
  integer, dimension(:), allocatable :: mTrack
  integer, dimension(:), allocatable :: rloTrack 
  integer, dimension(:), allocatable :: cloTrack 
  integer, dimension(:), allocatable :: rhiTrack 
  integer, dimension(:), allocatable :: chiTrack 
  integer :: cnt

  ! Use numroc to determine the number of rows and columns needed by the
  ! distributed array.
  nrows = numroc(valeDim, blkFact(1), myInfo(1), 0, gridDim(1))
  ncols = numroc(valeDim, blkFact(2), myInfo(2), 0, gridDim(2))


  allocate(local(nrows,ncols))
  
  ! Initialize the local array to zero
  local(:,:) = cmplx(0.0_double,0.0_double)

  ! First we need to calculate the number of blocks we'll be getting along
  ! each dimension of the local array.
  nbrow = int(nrows/blkFact(1))
  if ( mod(nrows,blkFact(1)) /= 0 ) then
    nbrow = nbrow + 1
  endif
  nbcol = int(ncols/blkFact(2))
  if ( mod(ncols,blkFact(2)) /= 0 ) then
    nbcol = nbcol + 1
  endif

  allocate(lTrack(nbrow*nbcol))
  allocate(mTrack(nbrow*nbcol))
  allocate(rloTrack(nbrow*nbcol))
  allocate(cloTrack(nbrow*nbcol))
  allocate(rhiTrack(nbrow*nbcol))
  allocate(chiTrack(nbrow*nbcol))
  lTrack(:) = -1
  mTrack(:) = -1
  rloTrack(:) = -1
  cloTrack(:) = -1
  rhiTrack(:) = -1
  chiTrack(:) = -1

  !print *, myInfo(3), "MyInfo: ", myInfo(1), myInfo(2), trailDim
  !print *, myInfo(3), "Rows and Cols: ", nrows, ncols
  !print *, myInfo(3), "!!BRows and BCols: ", nbrow, nbcol


  ! Using the equations in the documentation linked above, we can obtain
  ! the following set of equations.
  !  a=l*MB+x  I=(l*Pr + pr)*MB+x
  !  b=m*NB+y  J=(m*Pc + pc)*NB+y
  ! Where the coordinates in the local array are (a,b), coordinates in
  ! global array are (I,J), coordinates of local block are (l,m),
  ! coordinates inside local block are (x,y), and process 
  ! coordinates (pr,pc). For block size MBxNB, and process grid size PrxPc.
  ! With the equation for I we can solve it twice with x=0, and 
  ! x=MB-1, this will give us the upper left and lower right of the block.
  ! This will be sufficient information to pull the block down from
  ! the global array. As we loop over the local blocks.
  ! 
  ! The upper left coordinate will be held by lo, and the lower right held
  ! by hi. "+ 0" intentionally left in for clarity.
  !
  ! Lo dimensions need to be incremented by 1, because the math is based
  ! off of a 0 index. Hi should calculate correctly.
  !
  ! A special case has to be considered for edge blocks of irregular size.
  call ga_sync()
  do ld = 0, 5
    if (ld==myInfo(3)) then
      print *, myInfo(3), "Nbrow,nbCol: ", nbrow, nbcol
      print *, ""
    endif
    call ga_sync()
  enddo

  cnt = 1
  do i=1,nbrow
    do j=1,nbcol
      l = i-1
      m = j-1
      ! Upper left corner.
      lo(1) = ((l*gridDim(1) + myInfo(1))*blkFact(1) + 0) + 1
      lo(2) = ((m*gridDim(2) + myInfo(2))*blkFact(2) + 0) + 1
      
      ! Starting positions of a and b should be unnefected for every case
      ! Read next comment.
      a = (l*blkFact(1) + 0) + 1
      b = (m*blkFact(2) + 0) + 1

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
      if ( sum(trailDim(:))==2 .and. (i==nbrow) .and. (j==nbcol) ) then
        ! Lower right corner. 
        hi(1) = (l*gridDim(1) + myInfo(1))*blkFact(1) + trailDim(1)
        hi(2) = (m*gridDim(2) + myInfo(2))*blkFact(2) + trailDim(2)
        ! Calculate local starting coordinates
        aend = size(local,1)
        bend = size(local,2)

      else if ( (trailDim(2) > 0) .and. (j==nbcol) ) then ! trailing cols
        ! Lower right corner. 
        hi(1) = (l*gridDim(1) + myInfo(1))*blkFact(1) + blkFact(1)
        hi(2) = (m*gridDim(2) + myInfo(2))*blkFact(2) + trailDim(2)
        ! Calculate local starting coordinates
        aend = a + blkFact(1)
        bend = size(local,2) 

      else if ( (trailDim(1) > 0) .and. (i==nbrow) ) then ! trailing rows
        ! Lower right corner. 
        hi(1) = (l*gridDim(1) + myInfo(1))*blkFact(1) + trailDim(1)
        hi(2) = (m*gridDim(2) + myInfo(2))*blkFact(2) + blkFact(2)
        ! Calculate local starting coordinates
        aend = size(local,1) 
        bend = b + blkFact(1)  
                               
      else ! Normal Block
        ! Lower right corner. 
        hi(1) = (l*gridDim(1) + myInfo(1))*blkFact(1) + blkFact(1)
        hi(2) = (m*gridDim(2) + myInfo(2))*blkFact(2) + blkFact(2)
        ! Calculate local starting coordinates
        aend = a + blkFact(1)
        bend = b + blkFact(1)

      endif


      ! Set ld
      ld = (aend - a) + 1

!      lTrack(cnt) = l
!      mTrack(cnt) = m
!      rloTrack(cnt) = lo(1)
!      cloTrack(cnt) = lo(2)
!      rhiTrack(cnt) = hi(1)
!      chiTrack(cnt) = hi(2)
!      cnt = cnt + 1
      ! Get data from global array, storing it into local array at location
      ! a : (a + blkFact(1)), and b : (b + blkFact(2))
      call ga_get(ga_arr,lo(1),hi(1),lo(2),hi(2), &
        & local( a:aend, b:bend ), &
        & ld)
    enddo
  enddo


!  do ld=0, 5
!    call ga_sync()
!    if (ld==myInfo(3)) then
!      print *, myInfo(3), "pr,pc: ",  myInfo(1), myInfo(2) 
!      print *, "    l: ",  lTrack
!      print *, "    m: ",  mTrack
!      print *, "  rlo: ", rloTrack
!      print *, "  clo: ", cloTrack
!      print *, "  rhi: ", rhiTrack
!      print *, "  chi: ", chiTrack
!      print *, ""
!    endif
!  enddo
!  call ga_sync()
!  stop

  ! We should now have the global data in the proper local configuration
  ! for the grid size defined.

end subroutine oga_tosl

! In order to do the gridCalc subroutine we need to do integer factorization
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

! This subroutine determines the process grid to be used for the
! call to scalapack. There are a lot of considerations that need to be
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
subroutine gridCalc(myInfo, numprocs, gridDim, valeDim, blkFact)
  implicit none

  ! Define passed parameters
  integer, dimension(3), intent(in) :: myInfo
  integer, intent(in) :: numprocs
  integer, dimension(2), intent(out) :: gridDim
  integer, intent(in) :: valeDim
  integer, dimension(2), intent(out) :: blkFact ! blocking factors

  ! Local variables
  integer :: modifiedValeDim

  if (mod(valeDim,2) /= 0) then
    modifiedValeDim = valeDim - 1
  else
    modifiedValeDim = valeDim
  endif

  ! For now we'll assum that the user has used a perfect square number
  ! of processors. And we'll then make the processor grid square.
  call calcProcGrid(numProcs, gridDim(1), gridDim(2))

  ! Calculate the blocking factors to be used. From the calcProcGrid
  ! routine we know that gridDim1 (prows) will be <= gridDim2 (pcols).
  ! Because of this we'll do the integer divide on gridDim(2), ensuring
  ! that our blocks are square.
  blkFact(1) =  int(modifiedValeDim / gridDim(2))
  blkFact(2) =  int(modifiedValeDim / gridDim(2))

end subroutine gridCalc

! Subroutine that calculates the trailing dimensions for a process
subroutine getTrailDim(myInfo, valeDim, gridDim, blkFact, trailDim)
  implicit none

  ! Define passed data
  integer, dimension(3), intent(in) :: myInfo
  integer, intent(in) :: valeDim
  integer, dimension(2), intent(in) :: gridDim
  integer, dimension(2), intent(in) :: blkFact ! blocking factors
  integer, dimension(2), intent(out) :: trailDim ! Trailing Dimension

  ! Define local variables
  integer :: rTrailProc, cTrailProc

  ! remove later
  integer :: i
  
  ! Initialize trailDIm
  trailDim(:) = 0

  ! Because we decided that our blocks are going to be square. It makes it
  ! a touch more difficult to determine if a particular process has
  ! any trailing dimensions. From the process grid dimensions, we need to
  ! know whether or not the process calculating this is going to end up
  ! with trailing dimnesions according to the block-cyclic decomposition.
  !
  ! Luckily we can divide and modulo our asses off in order to get this
  ! info. In order to do this we need to know how many times the a process 
  ! row or process column is going to fit into the dimensions of the global
  ! matrix, given by valdim/(# process rows/cols). We then module that 
  ! number by the (# process rows/cols) with respect to the dimension.
  ! That should give us the process that the trailing dimension belongs
  ! to. We then have each process check whether that process is themself
  ! and if it is, they obtain the trailing dimension.
  ! that number by nprow again

  rTrailProc = mod( valeDim/blkFact(1) , gridDim(1))
  cTrailProc = mod( valeDim/blkFact(2) , gridDim(2))

  if ( rTrailProc == myInfo(1) ) then
    trailDim(1) = mod(valeDim, blkFact(1))
  endif
  if ( cTrailProc == myInfo(2) ) then
    trailDim(2) = mod(valeDim, blkFact(2))
  endif
! do i = 0, 5
!   call ga_sync()
!   if (i==myInfo(3)) then
!     print *, myInfo(3), "Proc Inf: ", myInfo(1), myInfo(2) 
!     print *, myInfo(3), "Proc tra: ", rTrailProc, cTrailProc
!     print *, myInfo(3), "TrailDim: ", trailDim
!     print *, ""
!   endif
! enddo
! call ga_sync()
  

end subroutine getTrailDim

! This subroutine servers as the main wrapper to the zhegvx lapack
! routine.
subroutine oga_pdzhegv(ga_a, ga_b, eigenVals, valeDim)
  use MPI
  use O_Kinds

  implicit none

#include "mafdecls.fh"
#include "global.fh"

  ! Define passed parameters
  integer, intent(inout) :: ga_a ! Global Array holding Matrix A
  integer, intent(inout) :: ga_b ! Global Array holding Matrix B
  ! Place to store eigenValues
  complex (kind=double), intent(out), dimension(:,:) :: eigenVals
  integer, intent(in) :: valeDim

  ! Define local matrices to hold parts from global arrays
  complex (kind=double), allocatable, dimension(:,:) :: localA
  complex (kind=double), allocatable, dimension(:,:) :: localB
  complex (kind=double), allocatable, dimension(:,:) :: localZ

  ! Define local variables
  integer :: mpierr, numprocs
  integer :: Atype, Adim1, Adim2
  integer :: Btype, Bdim1, Bdim2

  ! Variables used to describle BLACS and Scalapack
  integer :: slctxt
  integer, dimension(9) :: desca, descb, descz

  ! Allocate Scalapack parameters to calculate for the zhegvx call
! information on netlib.
  integer :: IBTYPE, N, IA, JA, IB, JB, VL, VU, IL, IU 
  integer :: ABSTOL, M, NZ, W, ORFAC, IZ, JZ
  integer :: LWORK, LRWORK, LIWORK
  integer :: IFAIL, ICLUSTR, GAP, INFO
  integer, allocatable, dimension(:) :: WORK, RWORK, IWORK
  CHARACTER :: JOBZ, RNGE, UPLO

  ! Define information arrays that store important information
  ! Grid dim holds the dimensions of the processor grid
  integer, dimension(2) :: gridDim
  ! blkFact holds the blocking factor used along each dimension
  integer, dimension(2) :: blkFact
  ! trailDim, holds the number of trailng rows or columns along the
  ! respective dimension
  integer, dimension(2) :: trailDim
  ! myInfo holds process specific information, myInfo(1) is the row 
  ! coordinate they belong to in the process grid. myInfo(2) is the 
  ! column coordinate they belong to in the process grid. myInfo(3) is
  ! the MPI rank of the process.
  integer, dimension(3) :: myinfo

  ! Define external functions
  real (kind=double) :: PDLAMCH
  external :: PDLAMCH

  ! Get environment parameters
  myInfo(3) = ga_nodeid()
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,mpierr)

  call ga_check_handle(ga_a, 'In subroutine oga_zhegv, passed ga_a is not &
    a valid Global Arrays handle')
  call ga_check_handle(ga_b, 'In subroutine oga_zhegv, passed ga_b is not &
    a valid Global Arrays handle')

  call ga_inquire(ga_a, Atype, Adim1, Adim2)
  call ga_inquire(ga_b, Btype, Bdim1, Bdim2)
  
  ! Check to make sure arrays are properly sized for the linear algebra
  ! operations
  if (Adim1 .ne. Adim2) then
    call ga_error('oga_pzhegv: matrix A not square ',0)
  endif
  if (Bdim1 .ne. Bdim2) then
    call ga_error('oga_pzhegv: matrix B not square ',0)
  endif
  if (Bdim1 .ne. Adim1) then
    call ga_error('oga_pzhegv: size of matrix A and B differ ',0)
  endif

  ! Sync proccesses
  call ga_sync()

  ! Determine the grid size and blocking factors we will use
  call gridCalc(myInfo, numprocs, gridDim, valeDim, blkFact)

  ! Initialize the SCALAPACK interface and great a process grid
  ! calculated with gridCalc
  call BLACS_GET(-1,0, slctxt)
  call BLACS_GRIDINIT(slctxt, 'r', gridDim(1), gridDim(2))
  call BLACS_GRIDINFO(slctxt, gridDim(1), gridDim(2), myInfo(1), myInfo(2))

  ! Determine the trailing dimensions, which has to be done after
  ! gridCalc and gridInit.
  call getTrailDim(myInfo, valeDim, gridDim, blkFact, trailDim)

  ! If not in participating in the scalapack call, call ga_sync and
  ! wait for all the other processes before returning. 
  !!!!!!!!!!!!! This condition should not evaluate to true as long as the
  !!!!!!!!!!!!! naive processor grid creation is still in place.
  if (myInfo(1) == -1) then
    call ga_sync()
    return
  endif

  ! Pull down data from global arrays to local arrays
  print *, 'OGA 1'
  call flush(6)
  call oga_tosl(ga_a, localA, valeDim, gridDim, numprocs, myinfo, blkFact, trailDim)
  print *, 'OGA 2'
  call flush(6)
  call oga_tosl(ga_b, localB, valeDim, gridDim, numprocs, myinfo, blkFact, trailDim)
  print *, 'Done with pull'
  call flush(6)

  ! Setup array Z to store the eigenValues. For now we'll try to get
  ! by without reorthogonalizing any eigenvectors and should keep
  ! this from being referenced
  allocate(localZ(1,1))

  ! Create Array Descriptors
  call descinit(desca, adim1, adim2, blkFact(1), blkFact(2),  & 
    & 0, 0, slctxt, size(localA,2), info)

  call ga_sync()
  call descinit(descb, bdim1, bdim2, blkFact(1), blkFact(2),  & 
    & 0, 0, slctxt, size(localB,2), info)
  ! Create array descriptor for eigenVector storage
  call descinit(descz, 0,0,1,1,0,0, slctxt, 1, info)

  ! Determine other SCALAPACK info defined in the pzhegvx documentation
  ! on netlib.
!http://www.netlib.org/scalapack/explore-html/d7/dff/pzhegvx_8f_source.html
  IBTYPE = 1
  JOBZ = 'V'
  RNGE = 'A'
  UPLO = 'U'
  N = Adim1

  ! Row and column indices of the locals first column/row in the globals
  IA = myinfo(1) * blkFact(1) 
  JA = myinfo(2) * blkFact(2)
  IB = IB 
  JB = JB

  ! These not referenced if RNGE = 'A'
  VL = 0
  VU = 0
  IL = 0
  IU = 0

  ! This defines the precision the eigenvalues should be calculated to.
  ! Specifying it the way it is below is the most accurate eigenvalue
  ! convergence.
  ABSTOL = 2 * PDLAMCH(slctxt,'S')

  ! Global outputs, just initializing to 0
  M = 0
  NZ = 0
  W = 0

  ! Specifies which eigenvectors should be reorthoganalized. For now
  ! we will specify that none of the eigenvectors should be
  ! reorthogonalized
  ORFAC = 0

  ! Row and column indices of the locals first column to the globals
  !IZ = 0
  !JZ = 0
  IZ = (myinfo(1) * blkFact(1))+1
  JZ = (myinfo(2) * blkFact(2))+1
  print *, 'tingy: ', myinfo(3), IZ, JZ
  call flush(6)

  ! Instead of trying to understand the appropriate size of the WORK
  ! arrays we will call PZHEGVX with LWORK, LRWORK, and LIWORK specified
  ! as -1. This will cause scalapack to do a workspace query and output
  ! the optimal sizes of these parameters in the first index of the
  ! assosciated arrays.
  LWORK =  -1
  allocate(WORK(1)) 
  WORK(:) = 0

  LRWORK = -1
  allocate(RWORK(1))
  RWORK (:) = 0
  
  LIWORK = -1 
  allocate(IWORK(1))
  IWORK(:) = 0

  ! Ouputs, just initalizing to 0
  IFAIL = 0
  ICLUSTR = 0  
  GAP = 0
  INFO = 0
  
  print *, 'lapack 1'
  call flush(6)
  call ga_sync()
  ! As stated above we call the PZHEGVX subroutine as a workspace query
  call PZHEGVX(IBTYPE, JOBZ, RNGE, UPLO, N, localA, IA, JA, DESCA, localB, &
    & IB, JB, DESCB, VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, localZ, &
    & IZ, JZ, DESCZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, IFAIL, &
    & ICLUSTR, GAP, INFO)

  stop

  ! Now we set the proper workspace parameters as needed. And resize
  ! the arrays.
  LWORK = WORK(1)
  LRWORK = RWORK(1)
  LIWORK = IWORK(1)

  deallocate(WORK)
  deallocate(RWORK)
  deallocate(IWORK)
  allocate(WORK(LWORK))
  allocate(RWORK(LRWORK))
  allocate(IWORK(LIWORK))
  WORK(:) = 0
  RWORK(:) = 0
  IWORK(:) = 0

  print *, 'lapack 2'
  call flush(6)
  ! Now we we have sufficient workspace for scalapack and can now 
  ! call PZHEGVX to actually solve our eigen problem.
  call PZHEGVX(IBTYPE, 'V', 'A', 'U', N, localA, IA, JA, DESCA, localB, &
    & IB, JB, DESCB, VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, localZ, &
    & IZ, JZ, DESCZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, IFAIL, &
    & ICLUSTR, GAP, INFO)
  call ga_sync()

  stop
  
  ! Redistribute the eigenvalues to all processes for later use.
!  call MPI_BROADCAST

  ! deallocate the local matrix allocated in oga_tosl
  deallocate(localA)
  deallocate(localB)
  deallocate(localZ)

  ! Close the SCALAPACK interface, and any descriptors
  !


  ! Sync all processes before returning to the calling subroutine
  call ga_sync()


end subroutine oga_pdzhegv


! This subroutine is responsible for reading in the interaction matrices
! from disk, and put them in global arrays.
! Its operation should be similar to the writeValeVale subroutine in the
! parallelSubs module, except in reverse.
subroutine readHDF5_to_GA()
  use MPI
  use O_Kinds

  implicit none

#include "mafdecls.fh"
#include "global.fh"



end subroutine readHDF5_to_GA

end module O_SlSubs
