!! This module is focused on the parallelization of OLCAOmain. 
!!
!! Author: James E. Currie
!! Email: jecyrd@mail.umkc.edu
!!
module O_ParallelMain

  ! Although each subroutine should specify implicit none, it is placed
  ! in this module to protect against accidentaly omission
  implicit none

  contains

! This subroutine servers as the wrapper to neatly tie up all the
! subroutines below.
subroutine solve_PZHEGVX(valeDim, numKPoints, potDim, spinDirection, &
    & hdfMatrixDescriptors)
  use O_Kinds

  implicit none
  integer, dimension(3) :: myInfo
  integer, intent(in) :: valeDim, numKPoints, potDim, spinDirection

  ! Read the hamiltonian and overlap matrices from disk into the proper
  ! local distributions
  call readDistributedWaveFunction()
  
  ! Solve the eigenvalue problem now that we're properly distributed.
  call sl_PZHEGVX(localA, localB, localZ, desca, descb, descz, &
       & myInfo, numprocs, gridDim, trailDim, valeDim, blkFact, eigenVals)

  ! First we'll dump the wavefunction to disk.
  call writeWaveFunction()

  ! Now we'll want to distribute the energy eigenalues for later use in
  ! the computation
  call distributeEigenVals()

  ! Now we'll deallocate any memory we've allocated above

  ! Finally we'll clean up scalapack contexts and descriptors
  call cleanUpSl(slctxt, desca, descb, descz)

end subroutine solve_PZHEGVX

! This Subroutine is responsible for reading only the blocks of data for the
! current processes local array. We need to loop over each block in the 
! distributed matrix, calculate the section of data to be read from HDF5.
! Then, because our data is hermitian, we first read in the real part from
! the upper half. The imaginary part from the lower half (transposed slab).
!
! Because parallel reads shouldn't hamper our performance we'll have each
! process define hyperslabs to read only what it needs.
subroutine readDistributedWaveFunction(localArr, nbrow, nbcol, gridDim, 
  & myInfo, blkFact, trailDim)

  use O_Kinds
  use O_PotTypes only .....
  use HDF5 only ....

  implicit none

  ! Define passed parameters
  complex (kind=double), intent(inout), dimension(:,:) :: localArr
  integer, intent(in) :: nbrow, nbcol ! Total Number of blocks
  integer, dimension(2), intent(in) :: gridDim ! Process Grid Dims 
  integer, dimension(3), intent(in) :: myInfo  ! prow, pcol, mpirank
  integer, dimension(2), intent(in) :: blkFact ! Blocking Factors
  integer, dimension(2), intent(in) :: trailDim ! trailing dims

  ! Define local variables
  integer :: blockRow
  integer :: blockCol
  integer :: mpisize, mpierr
  integer, dimension(2) :: glo, ghi, llo, lhi
  integer, dimension(2) :: ld
  real (kind=double), allocatable, dimension(:,:) :: tempReal
  real (kind=double), allocatable, dimension(:,:) :: tempImag
  real (kind=double), allocatable, dimension(:,:) :: tempRead

  ! Obtain world size
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpisize, mpierr)
  
  ! Create Dataspace and open for access
  call h5screate_simple_f(numdims, dimsizes(,), memspace_dsid, hdferr)

  ! We need to figure out the block indices assosciated with our
  ! particular process. Indirect quote from scalapack documentation: 
  !   For an array length N to be stored on P processes. By convention
  !   the array entries are numbered 1 through N and the processes are
  !   numbered 0 through P-1. The array is divided into contiguous blocks of
  !   size NB. When NB does not divide N evenly, the last block of array
  !   elements will only contain mod(N,NB) entries instead of NB. Blocks
  !   are numbered starting from zero and dealt out to the processes like a
  !   deck of cards. In other words, if we assume that the process 0
  !   receives the first block, the kth block is assigned to the process
  !   of coordinate mod(k,P).
  !
  ! Because the above applies to a one-dimensional array we need only apply
  ! it along both the rows and columns. so we only need to perform the
  ! inner two loops if:
  !        mod(i,mpisize)==myInfo(1) .and. mod(j,mpisize)==myInfo(2)
  do i=1, nbrow ! Loop over blockRows
    do j=1, nbcol ! Loop over blockCols
      ! Logical not of the condition above, because I didn't want to retab
      ! everything below.
      if ( mod(i,mpisize)/=myInfo(1) .or. mod(j,mpisize)/=myInfo(2) ) then
        cycle
      endif

      ! Now we can use the globaToLocalMap subroutine to calculate the lo
      ! and high parameters for our hyperslab for the upcomming hdf5 read
      call globalToLocalMap(i, j, nbrow, nbcol, gridDim, &
        & myInfo, blkfact, trailDim, glo, ghi, llo, lhi, ld)

      allocate(tempRead(ld(1),ld(2)))
      allocate(tempReal(ld(1),ld(2)))
      allocate(tempImag(ld(2),ld(1)))

      do j=1, numKPoints ! Loop over kPoints
        do l=1, potDim ! Loop over potential dimensions
          ! Zero out tempRead
          tempRead(:,:) = 0.0_double

          ! Set dimensions for real part
          hslabstart = 
          hslabcount = 

          ! Define hyperslab to read real part
          call h5sselect_hyperslab_f(valeVale_dsid, H5S_SELECT_SET_F, &
            & hslabStart, hslabCount, hdferr)

          ! Read real part from hdf5 to temporary matrix
          call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tempRead, 
            & hslabdims, hdferr, file_space_id=fspaceid, &
            & xfer_prp=transfer_properties)

          ! Accumulate real part into tempReal
          tempReal = tempReal + tempRead

          ! Zero out tempRead
          tempRead(:,:) = 0.0_double

          ! Set transposed dimensions for imaginary part
          hslabstart = 
          hslabcount = 
          
          ! Redefine hyperslab to read imaginary part 
          call h5sselect_hyperslab_f(valeVale_dsid, H5S_SELECT_SET_F, &
            & hslabStart, hslabCount, hdferr)

          ! Read imaginary part from hdf5 to temporary matrix
          call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tempRead, 
            & hslabdims, hdferr, file_space_id=fspaceid, &
            & xfer_prp=transfer_properties)

          ! Accumulate imaginary part into tempImag
          tempImag = tempImag + tempRead

          ! Combine the two temporary double matrices into the local array 
          ! (complex double matrix).
          call storeIntoLocal(localArr, tempReal, tempImag)
        enddo ! pot dim loop
      enddo ! kpoint loop

      deallocate(tempRead)
      deallocate(tempReal)
      deallocate(tempImag)
    enddo ! blockCols
  enddo ! blockRows

end subroutine

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

subroutine storeIntoLocal(localArr, tempReal, tempImag)
  implicit none

  ! Define passed parameters
  complex, (kind=double), intent(inout), dimension(:,:) :: localArr
  real, (kind=double), intent(inout), dimension(:,:) :: tempReal
  real, (kind=double), intent(inout), dimension(:,:) :: tempImag


end subroutine


! This subroutine serves as the main interface/wrapper to the zhegvx lapack
! routine.
subroutine sl_PZHEGVX(localA, localB, localZ, desca, descb, descz, &
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
  integer :: ABSTOL, M, NZ, IZ, JZ
  real (kind=double) :: ORFAC
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

  ! Global outputs, just initializing to 0
  M = 0
  NZ = 0

  ! Allocate space for eigenvectors and initialize
  allocate(W(valeDim))
  W(:) = 0.0_double

  ! Specifies which eigenvectors should be reorthogonalized. Eigenvectors
  ! that correspond to eigenvalues which are within tol=ORFAC*norm(A) of 
  ! each other are to be reorthogonalized.
  ! For now we'll use the default value of 10^-3
  ORFAC = 1E-3

  ! Row and column index in the global array Z indicating the first
  ! row of sub(Z)
  IZ = 1
  JZ = 1

  ! We will first call PZHEGVX with LWORK, LRWORK, and LIWORK specified
  ! as -1. This will cause scalapack to do a workspace query and output
  ! the optimal sizes of these parameters in the first index of the
  ! assosciated arrays. We'll then reallocate to the correct sizes.
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

end subroutine sl_PZHEGVX

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