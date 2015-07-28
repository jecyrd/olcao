!! This module is focused on the parallelization of OLCAOmain. 
!!
!! Author: James E. Currie
!! Email: jecyrd@mail.umkc.edu
!!
module O_SlSubs

  ! Although each subroutine should specify implicit none, it is placed
  ! in this module to protect against accidentaly omission
  implicit none

  contains

! This subroutine is responsible for calculating and allocating the
! local array holding sections of the global array.
subroutine allocLocalArray(local, valeDim, blkFact, myInfo, gridDim)
  use O_Kinds

  implicit none

  complex (kind=double), intent(inout), allocatable, dimension(:,:) :: local
  integer, intent(in) :: valeDim
  integer, intent(in) :: numprocs
  integer, dimension(3), intent(in) :: myinfo  ! prow, pcol, mpirank
  integer, dimension(2), intent(in) :: blkFact ! Blocking Factors
  integer, dimension(2), intent(in) :: gridDim ! Process Grid Dims 

  ! External subroutines
  integer, external :: numroc

  ! Use numroc to determine the number of rows and columns needed by the
  ! distributed array.
  nrows = numroc(valeDim, blkFact(1), myInfo(1), 0, gridDim(1))
  ncols = numroc(valeDim, blkFact(2), myInfo(2), 0, gridDim(2))

  allocate(local(nrows,ncols))
  
  ! Initialize the local array to zero
  local(:,:) = cmplx(0.0_double,0.0_double)

end subroutine allocLocalArray

! This subroutine calculates the number of blocks along each row and 
! column of the local array. It can then be used to loop over the
! globalToLocalMap subroutine to get the appropriate lo,hi parameters
! for taking information from a global array
subroutine numRowColBlocks(nbrow, nbcol, blkFact)
  implicit none
 
  integer, intent(out) :: nbrow, nbcol ! Number of blocks along 
  integer, intent(in) :: nrows, ncols
  integer, dimension(2), intent(in) :: blkFact ! Blocking Factors

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

end subroutine numRowColBlocks

! This subroutine takes a given block for a process, and calculates the
! mapping between the global and local proceses.
!
! For some process grid determined by the getBlockFact subroutine. We need
! to extract data from the global array in such away to fit our grid
! using the 2d block cyclic mapping specified in the SCALAPACK
! documentation on netlib. http://netlib.org/scalapack/slug/node76.html
subroutine globalToLocalMap(blockRow, blockCol, nbrow, nbcol, gridDim, &
    & myinfo, blkFact, trailDim, glo, ghi, llo, lhi, ld)

  implicit none

  ! Define passed parameters
  integer, intent(in) :: blockRow, blockCol ! Block row and col index to get
  integer, intent(in) :: nbrow, nbcol ! Total Number of blocks
  integer, dimension(2), intent(in) :: gridDim ! Process Grid Dims 
  integer, dimension(3), intent(in) :: myinfo  ! prow, pcol, mpirank
  integer, dimension(2), intent(in) :: blkFact ! Blocking Factors
  integer, dimension(2), intent(in) :: trailDim ! trailing dims
  integer, dimension(2), intent(out) :: glo, ghi, llo, lhi
  integer, dimension(2), intent(out) :: ld

  ! Define local Variables
  integer :: l,m

  call numRowColBlocks(nbrow, nbcol, blkfact)

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
  l = blockRow
  m = blockCol

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
  if ( sum(trailDim(:))==2 .and. (i==nbrow) .and. (j==nbcol) ) then
    ! Lower right corner. 
    ghi(1) = (l*gridDim(1) + myInfo(1))*blkFact(1) + trailDim(1)
    ghi(2) = (m*gridDim(2) + myInfo(2))*blkFact(2) + trailDim(2)
    ! Calculate local ending coordinates
    lhi(1) = llo(1) + trailDim(1) 
    lhi(2) = llo(2) + trailDim(2) 

  else if ( (trailDim(2) > 0) .and. (j==nbcol) ) then ! trailing cols
    ! Lower right corner. 
    ghi(1) = (l*gridDim(1) + myInfo(1))*blkFact(1) + blkFact(1)
    ghi(2) = (m*gridDim(2) + myInfo(2))*blkFact(2) + trailDim(2)
    ! Calculate local ending coordinates
    lhi(1) = llo(1) + blkFact(1)
    lhi(2) = llo(2) + trailDim(2)

  else if ( (trailDim(1) > 0) .and. (i==nbrow) ) then ! trailing rows
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
subroutine getBlkFact(myInfo, numprocs, gridDim, valeDim, blkFact)
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
end subroutine getblkFact

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

end subroutine getTrailDim

! This subroutine creates the scalapack array descriptors for the local
! Matrices
subroutine getArrDesc(slctxt, desca, desb, descz)
  implicit none

  ! Define passed parameters
  integer, intent(in) :: slctxt
  integer, intent(inout), dimension(9) :: desca, descb, descz

  call descinit(desca, Adim1, Adim2, blkFact(1), blkFact(2),  & 
    & 0, 0, slctxt, size(localA,1), info)

  call descinit(descb, Bdim1, Bdim2, blkFact(1), blkFact(2),  & 
    & 0, 0, slctxt, size(localB,1), info)

  call descinit(descz, Adim1, Adim2, blkFact(1),blkFact(2),0,0, &
    & slctxt, size(localZ,1), info)
end subroutine getArrDesc


! This subroutine creates a scalapack context and initializes the grid
subroutine initSl(slctxt, grimDim, myInfo)
  implicit none

  ! Define passed parameters
  integer, intent(out) :: slctxt
  integer, intent(inout), dimension(2) :: gridDim
  integer, intent(inout), dimension(3) :: myInfo 
  call BLACS_GET(-1,0, slctxt)
  call BLACS_GRIDINIT(slctxt, 'r', gridDim(1), gridDim(2))
  call BLACS_GRIDINFO(slctxt, gridDim(1), gridDim(2), myInfo(1), myInfo(2))
end subroutine


! This subroutine cleans up scalapack descriptors and contexts
subroutine cleanUpSl(slctxt, desca, descb, descz)
  implicit none

end subroutine cleanUpSl


! This subroutine is responsible for obtaining all the necessary information
! needed by the scalapack routines.
subroutine getAllSlInfo(myInfo, numprocs, gridDim, trailDim, valeDim, &
    & blkFact)
  implicit none
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

  call getBlockFact(myInfo, numprocs, gridDim, valeDim, blkFact)
  call getTrailDim(myInfo, valeDim, gridDim, blkFact, trailDim)
  call getBlkFact(myInfo, numprocs, gridDim, valeDim, blkFact)
end subroutine getAllSlInfo


! This subroutine servers as the main wrapper to the zhegvx lapack
! routine.
subroutine solvePZHEGVX(localA, localB, localZ, desca, descb, descz, &
    & myInfo, numprocs, gridDim, trailDim, valeDim, blkFact, eigenVals)
  use MPI
  use O_Kinds

  implicit none

  ! Define Passed Parameters
  complex (kind=double), intent(inout), dimension(:,:) :: localA
  complex (kind=double), intent(inout), dimension(:,:) :: localB
  complex (kind=double), intent(inout), dimension(:,:) :: localZ
  integer, intent(in), dimension(9) :: desca, descb, descz
  integer, intent(in) :: valeDim
  complex (kind=double), intent(out), dimension(:,:) :: eigenVals

  ! Define local variables
  integer :: mpierr, numprocs
  integer :: Atype, Adim1, Adim2
  integer :: Btype, Bdim1, Bdim2

  ! Variables used to describle BLACS and Scalapack
  integer :: slctxt

  ! Allocate Scalapack parameters to calculate for the pzhegvx call
  ! information on netlib.
  integer :: IBTYPE, N, IA, JA, IB, JB, VL, VU, IL, IU 
  integer :: ABSTOL, M, NZ, ORFAC, IZ, JZ
  real (kind=double), allocatable, dimension(:) :: W
  integer :: LWORK, LRWORK, LIWORK
  integer :: IFAIL, ICLUSTR, GAP, INFO
  complex (kind=double), allocatable, dimension(:) :: WORK
  real (kind=double), allocatable, dimension(:) :: RWORK
  integer, allocatable, dimension(:) ::  IWORK
  CHARACTER :: JOBZ, RNGE, UPLO

  ! Define external functions
  real (kind=double) :: PDLAMCH
  external :: PDLAMCH, PZHEGVX
  integer, external :: numroc

  ! Determine other SCALAPACK info defined in the pzhegvx documentation
  ! on netlib.
!http://www.netlib.org/scalapack/explore-html/d7/dff/pzhegvx_8f_source.html
  IBTYPE = 1
  JOBZ = 'V'
  RNGE = 'A'
  UPLO = 'L'
  N = Adim1

  ! Row and column indices of the locals first column/row in the globals
  IA = 1!firstGlobalA(1)
  JA = 1!firstGlobalA(2)
  IB = 1!firstGlobalB(1) 
  JB = 1!firstGlobalB(2) 

  ! These not referenced if RNGE = 'A'
  VL = 0
  VU = 0
  IL = 0
  IU = 0

  ! This defines the precision the eigenvalues should be calculated to.
  ! Specifying it the way it is below is the most accurate eigenvalue
  ! convergence.
  ABSTOL = 2 * PDLAMCH(slctxt,'S')
  print *, myInfo(3), "abs: ", abstol

  ! Global outputs, just initializing to 0
  M = 0
  NZ = 0

  ! Allocate space for eigenvectors and initialize
  allocate(W(valeDim))
  W(:) = 0.0_double


  ! Specifies which eigenvectors should be reorthoganalized. For now
  ! we will specify that none of the eigenvectors should be
  ! reorthogonalized
  ORFAC = 0

  ! Row and column indices of the locals first column to the globals
  IZ = 1!firstGlobalZ(1)
  JZ = 1!firstGlobalZ(2)

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
  
  ! As stated above we call the PZHEGVX subroutine as a workspace query
  call PZHEGVX(IBTYPE, JOBZ, RNGE, UPLO, N, localA, IA, JA, DESCA, localB, &
    & IB, JB, DESCB, VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, localZ, &
    & IZ, JZ, DESCZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, IFAIL, &
    & ICLUSTR, GAP, INFO)
 
  ! Now we set the proper workspace parameters as needed. And resize
  ! the arrays.
  print *, myInfo(3), "Workspace Query: ", work(1), rwork(1), iwork(1)
  LWORK =   WORK(1)
  LRWORK =  RWORK(1)
  LIWORK =  IWORK(1)

  deallocate(WORK)
  deallocate(RWORK)
  deallocate(IWORK)
  allocate(WORK(LWORK))
  allocate(RWORK(LRWORK))
  allocate(IWORK(LIWORK))
  WORK(:) = 0
  RWORK(:) = 0
  IWORK(:) = 0

  ! Now we we have sufficient workspace for scalapack and can now 
  ! call PZHEGVX to actually solve our eigen problem.
  call PZHEGVX(IBTYPE, JOBZ, RNGE, UPLO, N, localA, IA, JA, DESCA, localB, &
    & IB, JB, DESCB, VL, VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, localZ, &
    & IZ, JZ, DESCZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, IFAIL, &
    & ICLUSTR, GAP, INFO)

  stop
  
  ! Redistribute the eigenvalues to all processes for later use.
  ! call distributeEigenVals()


end subroutine oga_pzhegv

! This subroutine is responsible for making sure the eigenvalues are
! redistributed to all processes for use later
subroutine distributeEigenVals()
!  call MPI_BROADCAST

end subroutine


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

! Luckily after parallelizing setup I made routines for reading in the
! new hdf5 datasets, we can utilize these and then just put the resulting
! matrix into a global array
subroutine readMatrixToGA(dataSetID, matrixDims, ga_vv)
  use MPI
  use O_Kinds

  implicit none

#include "mafdecls.fh"
#include "global.fh"

  ! Define Passed parameters
  integer (hid_t), intent(in) :: dataSetID
  integer (hsize_t), dimension(2) :: matrixDims
  integer, intent(inout) :: ga_vv ! Global array identifier

  ! Define local variables
  integer :: dim1, dim2A
  real (kind=double), dimension(matrixDims(1),matrixDims(2) :: temp

  ! Series of routines taken from secularEquation to be edited later
  call readPackedMatrix(atomNucOverlap_did(i), packedValeVale, atomDims, &
    & valeDim, valeDim)

  call readPackedMatrixAccum(atomKEOverlap_did(i), packedValeVale, &
    & tempPackedValeVale, atomDims, 0.0_double, valeDim, valeDim)

  do j=1, potDim
    call readPackedMatrixAccum(atomPotOverlap_did(i,j), packedaleVale, &
        & tempPackedValeVale, atomDims, potCoeffs(j, spinDirection), &
        & valeDim, valeDim)
  enddo

  call unpackMatrix(valeVale(:,:,1,spinDirection, packedValeVale, valeDim,0)

  call readPackedMatrix(atomOverlap_did(i), packedValeVale, atomDims, &
    & valeDim, valeDim)

  call unpackMatrix(valeValeOL(:,:,1,1), packedValeVale,valeDim,0)


end subroutine readMatrix

end module O_SlSubs
