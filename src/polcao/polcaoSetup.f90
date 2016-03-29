! Documentation here
!
module O_ParallelSetup
  use MPI
  implicit none

  contains

  type AtomPair
    integer :: I
    integer :: J
  end type AtomPair

! This subroutine is used to balance an loop for use with MPI.  The input
! (toBalance) is the number of things that needs to be split up. The
! output (initialVal, finalVal) are the start and stop of array indices.
! The subroutine currently works in the way that if a loop is not evenly
! distributable then the extra elements are added first to the n process,
! then to the n-1 process, then to the n-2 process, and so on.
subroutine loadBalMPI(toBalance,initialVal,finalVal,myProc,numProcs)
  implicit none

  ! Define passed parameters.
  integer, intent(in) :: toBalance,numProcs,myProc
  integer, intent(out) :: initialVal, finalVal

  ! Define local variables.
  integer :: jobsPer, remainder
  integer :: mpiRank, mpiSize, mpiErr

  mpiSize=numProcs
  mpiRank=myProc

  jobsPer = int(toBalance / mpiSize)
  remainder = mod(toBalance,mpisize)

  initialVal = (jobsPer * mpiRank) + 1
  finalVal = (jobsPer * (mpiRank+1))

  if (mpiRank > (mpiSize - remainder)) then
    initialVal = initialVal + (remainder - (mpiSize - mpiRank))
    finalVal = finalVal + (remainder - (mpiSize - (mpiRank+1)))
  endif
  if (mpiRank == (mpiSize - remainder)) then
    finalVal = finalVal + 1
  endif

!    This is only needed for C type array indices. i.e. first index=0
!    initialVal = initialVal - 1
!    finalVal = finalVal - 1

end subroutine loadBalMPI

subroutine elecBalance(toBalance,initialVal,finalVal,numProcs,numPotSites)
  implicit none
  
  integer, intent(in) :: toBalance,numProcs,numPotSites
  integer :: jobsPer, remainder
  integer, intent(out) :: initialVal, finalVal
  integer :: actMpiRank, mpiSize, mpiErr
  integer :: adjMpiRank

  call MPI_COMM_RANK(MPI_COMM_WORLD, actMpiRank, mpiErr)
!  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpiSize, mpiErr)
  adjMpiRank = (actMpiRank)-numPotSites
  mpiSize=numProcs

  jobsPer = int(toBalance / mpiSize)
  remainder = mod(toBalance,mpisize)

  initialVal = (jobsPer * adjMpiRank) + 1
  finalVal = (jobsPer * (adjMpiRank+1))

  if (adjMpiRank > (mpiSize - remainder)) then
    initialVal = initialVal + (remainder - (mpiSize - adjMpiRank))
    finalVal = finalVal + (remainder - (mpiSize - (adjMpiRank+1)))
  endif
  if (adjMpiRank == (mpiSize - remainder)) then
    finalVal = finalVal + 1
  endif

!    This is only needed for C type array indices. i.e. first index=0
!    initialVal = initialVal - 1
!    finalVal = finalVal - 1

end subroutine elecBalance

! This subroutine find the indices of a 2x2 matrix that correspond to a
! packed matrix (as defined by blas/blacs/lapack). This has only been
! tested for NxN matrices. The subroutine takes the index of a
! packed matrix, and outputs the upper triangle (x,y) coordinates for 
! a unpacked matrix.
subroutine findUnpackedIndices(packedIndex, x, y)
  implicit none

  integer, intent(in) :: packedIndex
  integer, intent(out) :: x, y
  integer :: pIndex 

  pIndex = packedIndex
  pIndex = pIndex - 1
  y = int((-1+ int(sqrt(8.0*real(pIndex)+1.0))) / 2)
  x = pIndex - y*(y+1)/2
  x = x + 1
  y = y + 1

  ! note for future, to go the reverse way
  ! packedIndex = x + (y*(y+1))/2
end subroutine findUnpackedIndices

! This subroutine does the reverse of the one above.
subroutine findPackedIndex(packedIndex, x, y)
  implicit none

  integer, intent(out) :: packedIndex
  integer, intent(in) :: x,y
  
  packedIndex = x + (y*(y+1))/2
  
end subroutine findPackedIndex

! This subroutine is responsible for deciding the atom pairs needed to be
! calculated by a process in order to fill their local block regions. It is to
! be used before the Atom-Atom loop in the integralsSCF subroutines.
!
! We start by looping over the local blocks, then using the subroutine
! localToGlobalMap to find our indices. Global to local map will give us
! 4 indices corresponding to the "top left" and "bottom right" elements we need
! from the global arrays. We can then do searches using getValeAtom, getCoreAtom
! which will define the atom pairs necessary to fill our arrays.
subroutine getAtomPairs(valeArrInfo, coreArrInfo, cvArrInfo, blcsinfo, &
    & atomPairs)
  use O_Parallel, only :: localToGlobalMap

  implicit none

  ! Define passed parameters
  type(ArrayInfo), intent(inout) :: valeArrinfo
  type(ArrayInfo), intent(inout) :: coreArrinfo
  type(ArrayInfo), intent(inout) :: cvArrinfo
  type(BlacsInfo), intent(inout) :: blcsinfo
  type(AtomPair),  intent(inout), allocatable, dimension(:) :: atomPairs

  ! Define local variables
  integer, dimension(2) :: glo, ghi ! Returned indices from globalToLocalMAp
  integer, dimension(2) :: vlo, vhi ! vale returned indices from search 
  integer, dimension(2) :: clo, chi ! core returned indices from search 
  integer, dimension(2) :: cvlo, cvhi ! CV returned indices from search 
  integer :: nrBlocks ! The number of blocks along the local row
  integer :: ncBlocks ! The number of blocks along the local column
  integer :: a,b ! Starting block indices to feed to localToGlobalMap
  integer :: i,j

  ! First we need to calculate the blocks along the local row and columns
  ! so we can loop over them. For now we'll neglect any fringe cases where
  ! we have odd shaped blocks.
  nrBlocks = len(arrinfo%local,1)/arrinfo%mb
  ncBlocks = len(arrinfo%local,1)/arrinfo%nb

  ! Now we loop over these blocks and record the results in atomPairs
  do i=1, nrBlocks
    do j=1, ncBlocks
      ! For localToGlobalMap we need the the starting indices of the block (a,b)
      ! This is just the array index minus 1, times the block size, plus 1
      aVale = (i-1)*valeArrinfo%mb+1
      bVale = (j-1)*valeArrinfo%nb+1
      aCore = (i-1)*coreArrinfo%mb+1
      bCore = (j-1)*coreArrinfo%nb+1
      aVC   = (i-1)*cvArrinfo%mb+1
      bVC   = (j-1)*cvArrinfo%nb+1

      ! Now we get the vale indices having to be searched
      call localToGlobalMap(a,b, glo, ghi, valeArrinfo, blcsinfo)
      ! Now we do the searches for the vale atoms
      call getValeAtom(glo(1), vlo(1))
      call getValeAtom(glo(2), vlo(2))
      call getValeAtom(ghi(1), vhi(1))
      call getValeAtom(ghi(2), vhi(2))
      
      ! Now we get the core indices having to be searched
      call localToGlobalMap(a,b, glo, ghi, coreArrinfo, blcsinfo)
      ! Now we do the searches for the core atoms
      call getCoreAtom(glo(1), clo(1))
      call getCoreAtom(glo(2), clo(2))
      call getCoreAtom(ghi(1), chi(1))
      call getCoreAtom(ghi(2), chi(2))

      ! Now we get the coreVale indices having to be searched
      call localToGlobalMap(a,b, glo, ghi, cvArrinfo, blcsinfo)
      ! Now we do the searches for the valeCore atoms
      call getcoreAtom(glo(1), cvlo(1))
      call getValeAtom(glo(2), cvlo(2))
      call getcoreAtom(ghi(1), cvhi(1))
      call getValeAtom(ghi(2), cvhi(2))

      ! Now we have to test whether or not the returned indices have any overlap
      ! so we can avoid doing any working. This'll be done by the testAtomDupe
      ! subroutine below.
      call testAtomDupe(vlo,vhi,clo,chi,cvlo,cvhi, atomPairs)

    enddo
  enddo
end subroutine getAtomPairs

! This function checks if two rectangles overlap. Our rectangles are given by
! lo (upper left corner) and hi (lower right corner). There are 4 conditions we
! can check to see if two rectangles overlap. We'll call them R1 and R2 to refer
! to all 4 coordinates specified and use "lo" and "hi" to illustrate the 
! conditions, with respective 1 or 2. These conditions are applications of
! deMorgan's law.
!  1.) If the R1 left edge is to the right of R2 right edge. Then no overlap.
!         or if  lo1(1) =< hi2(1)  then there is overlap
!
!  2.) If R1 right edge is to the left of R2 left edge. Then no overlap.
!         or if  hi1(1) >= lo2(1)  then there is overlap
!
!  3.) If R1 top edge is below R2 bottom edge. Then no overlap.
!         or if  lo1(2) >= hi2(2)  then there is overlap
!
!  4.) If R1 bottom edge is above R2 top edge. Then no overlap.
!         or if  hi1(2) =< lo2(2)  then there is overlap
!
! Well this is all good. But this only works for cartesian style coordinate
! systems. For matrix notice we are effectively working in cartesian quadrant 4.
! Also with matrix notation our first element corresponds to rows, which is 
! effectively the cartesian Y coordinate. So for the conditions above, we must 
! transpose our indices such that index 1 is swapped with index 2. And the 
! second 2 conditions must be negative ie: condition 3: -lo1(1)>-hi2(2). This
! transform will effectively put the matrix indices into the correct
! cartesian notation.
function checkRectOverlap(lo1, hi1, lo2, hi2)
  implicit none

  ! Define self and passed parameters
  integer :: checkRangeOverlap
  integer, intent(in), dimension(2) :: lo1, hi1, lo2, hi2

  if ( (lo1(2)=<hi2(2)) .and. (hi1(2)>=lo2(2)) .and. &
     & (-lo1(1)>=-hi2(1)) .and. (-hi1(1)<=-lo2(1)) ) then
    return 1
  endif

  ! else    
  return 0
end function checkRangeOverlap

! This subroutine finds the intersection of two rectangles if they overlap.
! Only use this subroutine if checkRectOverlap function returns true.
! Typically all four statements would be max, but as detailed above the 
! matrix notation is a little different than cartesian coordinate system,
! so we use min for the bottom right corner of the overlap.
subroutine getOverlapRect(lo1, hi1, lo2, hi2, loSect, hiSect)
  implicit none

  ! Define passed parameters
  integer, intent(in), dimension(2) :: lo1, hi1, lo2, hi2
  integer, intent(out), dimension(2) :: loSect, hiSect

  loSect(1) = max(lo1(1), lo2(1))
  loSect(2) = max(lo1(2), lo2(2))
  hiSect(1) = min(hi1(1), hi2(1))
  hiSect(2) = min(hi1(2), hi2(2))

end subroutine getOverlapRect

! This function returns true if an i,j pair is inside some rectangle, where
! loOvlp is the top left corner, and hiOvlp is the bottom right corner.
function inRectangle(i,j, loOvlp, hiOvlp)
  implicit none

  ! Define passed parameters and return value
  integer :: inintersection
  integer, intent(in) :: i,j
  integer, intent(in), dimension(2) :: loOvlp, hiOvlp

  if ( ((i>=loOvlp[0]) .and. (i<=hiOvlp[0])) .and. &
    &  ((j>=loOvlp[1]) .and. (j<=hiOvlp[1])) ) then
    return 1
  endif
  ! else
  return 0
end function inRectangle

! This subroutine appends an atom pair to the list of atomPairs
subroutine addAtomPair(i,j, atomPairs)
  implicit none

  ! Define passed parameters
  integer, intent(in) :: i,j
  type(AtomPair),  intent(inout), allocatable, dimension(:) :: atomPairs

  ! Define local variables
  type(AtomPair), allocatable, dimension(:) :: tmpPairs
  integer :: newsize

  newsize = size(atomPairs,1)+1
  ! Allocate some initial space in tmpPairs

  allocate(tmpPairs(newSize))

  tmpPairs = atomPairs
  tmpPairs(size(tmpPairs,1))%I = i
  tmpPairs(size(tmpPairs,1))%J = j

  deallocate(atomPairs)
  allocate(atomPairs(newSize))
  atomPairs = tmpPairs

  deallocate(tmpPairs)

end subroutine addAtomPairs

! This subroutine is responsible for checking the ranges of atoms that need
! to be calculated in order to fill the local matrices, to make sure that no
! extra work will be done. It'll then copy the atomPairs structure into a 
! temporary structure, so it can reallocate the atomPairs and add the necessary
! atom ranges.
subroutine testAtomDupe(vlo, vhi, clo, chi, cvlo, cvhi, atomPairs)
  implicit none

  ! Define passed parameters
  integer, intent(in), dimension(2) :: vlo, vhi, clo, chi, cvlo, cvhi
  type(AtomPair),  intent(inout), allocatable, dimension(:) :: atomPairs

  ! Define local variables
  integer, dimension(2) :: loOvlp1, hiOvlp1, loOvlp2, hiOvlp2
  integer :: i,j
  
  ! Initialize intersections 
  loOvlp1(:) = 0
  hiOvlp1(:) = 0
  loOvlp2(:) = 0
  hiOvlp2(:) = 0

  ! First add all vale atom pairs
  do i=vlo(1), vlo(2)
    do j=vhi(1), vhi(2)
      addAtomPairs(i,j, atomPairs)
    enddo
  enddo

  ! Check if vale and core overlap
  if ( checkRectOverlap(vlo, vhi, clo, chi) ) then
    call getOverlapRect(vlo, vhi, clo, chi, loOvlp1, hiOvlp1)
    ! Add appropriate ones to atom pairs
  endif 

  ! Now add the core pairs
  do i=clo(1), clo(2)
    do j=chi(1), chi(2)
      if (.not. inRectangle(i,j, loOvlp1, hiOvlp1) then
        addAtomPairs(i,j, atomPairs)
      endif
    enddo
  enddo

  ! Now in order to add the valeCore we need the overlaps of core with coreVale
  ! and the overlap of vale with coreVale.
  loOvlp1(:) = 0
  hiOvlp1(:) = 0
  ! Check if vale and coreVale overlap
  if ( checkRectOverlap(vlo, vhi, cvlo, cvhi) ) then
    call getOverlapRect(vlo, vhi, cvlo, cvhi, loOvlp1, hiOvlp1)
    ! Add appropriate ones to atom pairs
  endif

  ! Check of core and coreVale overlap
  if ( checkRectOverlap(clo, chi, cvlo, cvhi) ) then
    call getOverlapRect(clo, chi, cvlo, cvhi, loSect, hiSect)
    ! Add appropriate ones to atom pairs
  endif
  
  ! Now we add the valeCore pairs
  do i=vclo(1), vclo(2)
    do j=vchi(1), vchi(2)
      if ( (.not. inRectangle(i,j, loOvlp1, hiOvlp1)) .and.
           (.not. inRectangle(i,j, loOvlp2, hiOvlp2)) )then
        addAtomPairs(i,j, atomPairs)
      endif
    enddo
  enddo


end subroutine testAtomDupe

! Given an index gIndx=(1,valeDim), along a dimension gArrayDim=(1,2),
! return the atom indx needed so that the local array can be filled.
!
! Future improvement, write out the conditions necessary so that getValeAtom
! can be combined with getCoreAtom in the same loop. Currently I just want
! things to work. Additionally evaluate the necessity of a binary search over
! a linear one.
subroutine getValeAtom(gIndx, valeIndx)

  use O_AtomicSites, only: atomSites

  implicit none

  ! Define Passed Parameters
  integer, dimension(2), intent(in) :: gIndx ! Global Area of interest
  integer, intent(out) :: valeIndx ! What we're looking for

  if ( gIndx <= atomSites(1)%cumulValeStates ) then
    valeIndx == 1
  else
    do i=1, len(atomSites)-1
      if (gIndx >= atomSites(i)%cumulValeStates .and. &
            gIndx < atomSites(i+1)%cumulValeStates) then
        valeindx = i 
      endif
    enddo
  endif
end subroutine getValeAtom

! Given an index gIndx=(1,coreDim), along a dimension gArrayDim=(1,2),
! return the atom indx needed so that the local array can be filled.
subroutine getCoreAtom(gIndx, coreIndx)

  use O_AtomicSites, only: atomSites

  implicit none

  ! Define Passed Parameters
  integer, dimension(2), intent(in) :: gIndx ! Global Area of interest
  integer, intent(out) :: coreIndx

  if ( gIndx <= atomSites(1)%cumulCoreStates ) then
    coreIndx == 1
  else
    do i=1, len(atomSites)-1
      if (gIndx >= atomSites(i)%cumulCoreStates .and. &
            gIndx < atomSites(i+1)%cumulCoreStates) then
        coreIndx = i 
      endif
    enddo
  endif
end subroutine getCoreAtom

! Purpose: This subroutine undoes some of the load balancing that is done
! previously. It would be extremely easy to combine this with the two other
! load balancing routines earlier in this code. But due to time constraints
! for myself, (James) I'm taking the previouosly mentioned route of 
! "this works and it was easy to implement really quick"
subroutine writeResolve(toBalance, initialVal, finalVal,numProcs,tmpRank,valeDim)
  implicit none
  
  integer, intent(in) :: toBalance,numProcs,tmpRank
  integer, intent(in) :: valeDim
  integer :: jobsPer, remainder
  integer, intent(out) :: initialVal, finalVal
  integer :: mpiRank, mpiSize, mpiErr

!  call MPI_COMM_RANK(MPI_COMM_WORLD, mpiRank, mpiErr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpiSize, mpiErr)
!  mpiRank=tmpRank
!  print *, 'VALEDIM',valeDim
  jobsPer = int(toBalance / numProcs)
  remainder = mod(toBalance,mpisize)

  initialVal = (jobsPer * tmpRank) + 1
  finalVal = (jobsPer * (tmpRank+1))
  if (remainder>0) then
    if (tmpRank < (mpiSize-1)) then
      initialVal = initialVal + tmpRank
      finalVal = finalVal + tmpRank + 1
    endif
    if (tmpRank==(mpiSize-1)) then
      initialVal = initialVal + tmpRank
      finalVal = valeDim
    endif
  endif
!  if (tmpRank == (mpiSize-1) .and. remainder>0) then
!    initialVal= initialVal+1
!    finalVal = valeDim
!  endif

!    This is only needed for C type array indices. i.e. first index=0
!    initialVal = initialVal - 1
!    finalVal = finalVal - 1

end subroutine writeResolve

! This subroutine writes the valeVale matrix to disk, for the case that the
! orthogonaliztion is done. I realize it is kind of confusing and maybe
! not so straight forward so I will do my best to explain it.
! What happens is that for a small block of the matrix, we pull it down
! from a global array, to a local one on the single process execcuting
! the write to HDF5 sections from integralSCF. It then writes it to disk,
! then moves on to another small block and writes that to disk.
! The reason this was done instead of PHDF5 was because of
! NFS (Network File System) writes being inneficcient under the PHDF5
! paradigm as of HDF5 version 1.8.12.
subroutine writeValeVale(ga_valeVale,opCode,numKPoints,potDim, & 
    & currPotTypeNumber,currAlphaNumber,valeDim)
  use HDF5
  use O_Kinds
  use O_Constants
  use O_SetupIntegralsHDF5
  use O_PotTypes
  
  implicit none
#include "mafdecls.fh"
#include "global.fh"

  integer, intent(in),dimension(:) :: ga_valeVale
  integer, intent(in) :: opCode,valeDim
  integer, intent(in) :: currPotTypeNumber,currAlphaNumber
  integer, intent(in) :: numKPoints, potDim
  integer, dimension(2) :: hi, lo, blockDims, numBlocks
  integer :: mpiRank, mpiSize, mpiErr, hdferr
  complex (kind=double), allocatable, dimension(:,:) :: valeValeGA
  real    (kind=double), allocatable, dimension(:,:) :: diskVV
  integer :: i,j,k,m,x,y, xDimCnt, yDimCnt
  integer(hid_t) :: memspace_dsid
  integer(hid_t), dimension(numKPoints,potDim) :: datasetToWrite_did
  integer(hsize_t), dimension(2) :: hslabCount,hslabStart,dims
  integer :: valeBlockDim
  integer :: minDim,maxDim,tmpRank
  ! Define small threshold for eliminating resultant values
  real (kind=double) :: smallThresh10

  smallThresh10 = real(1.0d-10,double)

  call MPI_Comm_rank(MPI_COMM_WORLD, mpiRank, mpiErr)
  call MPI_Comm_size(MPI_COMM_WORLD, mpiSize, mpiErr)

  select case (opCode)
  case (1)
    datasetToWrite_did(:,1) = atomOverlap_did(:)
  case (2)
    datasetToWrite_did(:,1) = atomKEOverlap_did(:)
  case (3)
    datasetToWrite_did(:,1) = atomNucOverlap_did(:)
  case (4)
    datasetToWrite_did(:,:) = atomPotOverlap_did(:,:)
  case default
    print *, "wrong opCode passed to writeValeVale"
    stop
  end select
  
  if (valeDim/mpiSize<1) then
    valeBlockDim=valeDim
  elseif(valeDim/mpiSize>1 .and. mod(valeDim,mpiSize)/=0) then
    valeBlockDim=(valeDim/mpiSize) + 1
  else
    valeBlockDim=valeDim/mpiSize
  endif
  blockDims(1) = valeBlockDim
  blockDims(2) = valeBlockDim
  if (mod(valeDim,blockDims(1)) > 0) then
    numBlocks(1) = int(valeDim/blockDims(1)) + 1
  else
    numBlocks(1) = int(valeDim/blockDims(1))
  endif
  numBlocks(2) = numBlocks(1)
!  print *, "blockDims: ", blockDims
!  print *, "numBlocks: ", numBlocks(1)
!  print *, "Divide:    ", valeDim/blockDims(1)
!  print *, "mod:       ", mod(valeDim,blockDims(1))

  do i=1, numKPoints
    x=0
    y=0
    do m=1,numBlocks(1)**2
      if (x /=  numBlocks(1)-1 .and. y /= numBlocks(2)-1) then
        lo(1) = ((x)*blockDims(1))+1
        lo(2) = ((y)*blockDims(2))+1
        hi(1) = (x+1)*blockDims(1)
        hi(2) = (y+1)*blockDims(2)
      elseif (x == numBlocks(1)-1 .and. y /= numBlocks(2)-1) then
        lo(1) = valeDim-blockDims(1)+1
        lo(2) = ((y)*blockDims(2))+1
        hi(1) = valeDim
        hi(2) = (y+1)*blockDims(2)
      elseif (x /= numBlocks(1)-1 .and. y == numBlocks(2)-1) then
        lo(1) = ((x)*blockDims(1))+1
        lo(2) = valeDim-blockDims(2)+1
        hi(1) = (x+1)*blockDims(1)
        hi(2) = valeDim
      else
        lo(1) = valeDim-blockDims(1)+1
        lo(2) = valeDim-blockDims(2)+1
        hi(1) = valeDim
        hi(2) = valeDim
      endif
      
      allocate(valeValeGA(hi(1)-lo(1)+1, hi(2)-lo(2)+1))
      allocate(diskVV(hi(1)-lo(1)+1, hi(2)-lo(2)+1))
      call nga_get(ga_valeVale(i),lo,hi,valeValeGA, size(valeValeGA,1))
       
      yDimCnt=lo(2)
      do k=1,size(valeValeGA,2)
        xDimCnt=lo(1)
        do j=1,size(valeValeGA,1)
          if (yDimCnt<xDimCnt) then
            diskVV(j,k) = aimag(conjg(valeValeGA(j,k)))
          else
            diskVV(j,k) = real(valeValeGA(j,k))
          endif
          xDimCnt = xDimCnt + 1
        enddo
        yDimCnt = yDimCnt + 1
      enddo
      
      do k=1,size(diskVV,2)
        do j=1,size(diskVV,1)
         if(abs(diskVV(j,k))<smallThresh10) then
            diskVV(j,k) = 0.0_double
         endif
        enddo
      enddo

      hslabStart(1) = lo(1)-1
      hslabStart(2) = lo(2)-1
      
      hslabCount(1)=(hi(1)-lo(1))+1
      hslabCount(2)=(hi(2)-lo(2))+1
      call h5screate_simple_f(2,hslabCount,memspace_dsid,hdferr)
        
      ! define the hyperslab to be written to
      call h5sselect_hyperslab_f(valeVale_dsid,H5S_SELECT_SET_F, &
        & hslabStart,hslabCount,hdferr) 
      select case (opCode)
      case(1:3)
        ! Write slabs to disk
        call h5dwrite_f(datasetToWrite_did(i,1),H5T_NATIVE_DOUBLE, &
          & diskVV(:,:), hslabCount, hdferr, &
            & file_space_id=valeVale_dsid, mem_space_id=memspace_dsid)
      case(4)
        ! Write slabs to disk
        call h5dwrite_f(datasetToWrite_did(i, &
        & potTypes(currPotTypeNumber)%cumulAlphaSum+currAlphaNumber), &
        & H5T_NATIVE_DOUBLE, diskVV(:,:), hslabCount, &
        & hdferr, file_space_id=valeVale_dsid, &
        & mem_space_id=memspace_dsid)
      case default
        print *, "shits fucked up yo"
      end select

      x=x+1
      if (x==numBlocks(1)) then
        y = y +1
        x = 0
      endif

      deallocate(valeValeGA)
      deallocate(diskVV)
      call h5sclose_f(memspace_dsid,hdferr)
    enddo
  enddo

end subroutine writeValeVale

end module O_ParallelSetup
