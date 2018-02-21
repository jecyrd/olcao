! Documentation here
!
module O_ParallelSetup
  use MPI
  implicit none

  ! This is a linked list element
  type AtomPair
    integer :: I
    integer :: J
    type(atomPair), pointer :: next
  end type AtomPair

contains

! This subroutine initializes a new element in the atomPair linked list
subroutine initAtomPair(atomNode)
  implicit none

  ! Define passed parameters
  type(AtomPair), pointer :: atomNode

  allocate(atomNode)

  atomNode%next => null()
  atomNode%I = -1
  atomNode%J = -1
end subroutine initAtomPair

recursive subroutine destroyAtomList(root)
  implicit none

  ! Define passed parameters
  type(AtomPair), pointer :: root

  if (associated(root%next)) then
    call destroyAtomList(root%next)
  endif

  deallocate(root)

end subroutine destroyAtomList


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

  ! Define this function
  logical :: checkRectOverlap

  ! Define self and passed parameters
  integer, intent(in), dimension(2) :: lo1, hi1, lo2, hi2

  !if ( (lo1(2)=<hi2(2)) .and. (hi1(2)>=lo2(2)) .and. &
  if ( (lo1(2)<=hi2(2)) .and. (hi1(2)>=lo2(2)) .and. &
     & (-lo1(1)>=-hi2(1)) .and. (-hi1(1)<=-lo2(1)) ) then
    checkRectOverlap = .true.
    return
  endif

  ! else    
  checkRectOverlap = .false.
  return
end function checkRectOverlap

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

  ! Define this function
  integer :: inRectangle

  ! Define passed parameters and return value
  integer, intent(in) :: i,j
  integer, intent(in), dimension(2) :: loOvlp, hiOvlp

  if ( ((i>=loOvlp(0)) .and. (i<=hiOvlp(0))) .and. &
    &  ((j>=loOvlp(1)) .and. (j<=hiOvlp(1))) ) then
    inRectangle = 1
    return
  endif
  ! else
  inRectangle = 0
  return
end function inRectangle


subroutine getAtomPairs(vvinfo, ccinfo, cvinfo, blcsinfo, atomPairs, atomTree)
  use O_Parallel, only: ArrayInfo, BlacsInfo
  use O_23bst, only: tree_init, tree_destroy

  implicit none
 
  ! Define passed parameters
  type(ArrayInfo), intent(inout) :: vvinfo, ccinfo, cvinfo
  type(BlacsInfo), intent(inout) :: blcsinfo
  type(AtomPair), pointer :: atomPairs
  type(bst_node), pointer, intent(inout) :: atomTree ! 2-3 tree to keep 
                                           ! us from having duplicate atoms, 
                                           ! and to keep track of the local 
                                           ! blocks the atom pairs belong to

  ! Define local variables
  integer, dimension(2) :: vvlo, vvhi
  integer, dimension(2) :: cclo, cchi
  integer, dimension(2) :: cvlo, cvhi


  ! Need to figure out a better way to initialize the atomTree but for now
  ! we're just gonna put a zero in it, representing cantor pairing of 0,0
  ! which we should never have anyway.
  call tree_init(atomTree, 0)

  ! We need to call getArrAtomPairs 3 times one for each array
  call getArrAtomPairs(vvinfo, blcsinfo, atomPairs, atomTree, 0)
  call getArrAtomPairs(ccinfo, blcsinfo, atomPairs, atomTree, 1)
  call getArrAtomPairs(cvinfo, blcsinfo, atomPairs, atomTree, 2)

  ! By this time we're finished with our atomTree, so we can destroy it.
  ! All our atom pair data should be stored in the atomPairs linked list
  call tree_destroy(atomTree)

end subroutine getAtomPairs

! This subroutine is responsible for deciding the atom pairs needed to be
! calculated by a process in order to fill their local block regions. It is to
! be used before the Atom-Atom loop in the integralsSCF subroutines.
!
! We start by looping over the local blocks, then using the subroutine
! localToGlobalMap to find our indices. Global to local map will give us
! 4 indices corresponding to the "top left" and "bottom right" elements we need
! from the global arrays. We can then do searches using getAtoms which will 
! define the atom pairs necessary to fill our arrays.
!
! This subroutine can only be run once for each local array. So it'll need to 
! be called 3 times, once each for valeVale, coreCore, valeCore.
subroutine getArrAtomPairs(arrinfo, blcsinfo, atomPairs, atomTree, whichArr)
  use O_Parallel, only: localToGlobalMap, ArrayInfo, BlacsInfo
  use O_23bst, only: bst_node, tree_init, tree_destroy

  implicit none

  ! Define passed parameters
  type(ArrayInfo), intent(inout) :: arrinfo
  type(BlacsInfo), intent(inout) :: blcsinfo
  type(AtomPair),  pointer :: atomPairs
  type(bst_node), pointer :: atomTree
  integer :: whichArr ! whichArr is a control variable to tell the subroutine
                      ! which array it is working with. Options are:
                      ! = 0  valeVale array
                      !   1  coreCore array
                      !   3  coreVale array
  

  ! Define local variables
  integer, dimension(2) :: glo, ghi ! Returned indices from globalToLocalMAp
  integer, dimension(2) :: alo, ahi
  integer :: a,b ! Starting block indices to feed to localToGlobalMap
  integer :: i,j ! Loop Variables
  integer :: extra ! Variable to denote irregular size in one dimension
                   ! see localToGlobalMap in O_Parallel for description
  integer :: firstDim, secondDim  ! These are control variables that are set
                                  ! based on whichArr. These signal the getAtom
                                  ! subroutine which cumul states it should use

  ! First we need to set which cumul states we are working with.
  ! In the case of valeVale both are set to 0. In the case of coreCore both are
  ! set to 1. In the case of coreVale, firstDim=1, secondDim=0
  if (whichArr == 0) then
    firstDim = 0
    secondDim = 0
  elseif (whichArr == 1) then
    firstDim = 1
    secondDim = 1
  elseif (whichArr == 2) then
    firstDim = 1
    secondDim = 0
  endif

  ! If we have extra rows and columns because the blocking factor didn't divide
  ! the dimensions of the large array equally, then we need to consider 1
  ! more block, which is smaller than usual.
  ! if (arrinfo%extraRows > 0) then
  !   arrinfo%nrblocks = arrinfo%nrblocks + 1
  ! endif
  ! if (arrinfo%extraCols > 0) then
  !   arrinfo%ncblocks = arrinfo%ncblocks + 1
  ! endif

  ! Now we loop over these blocks and record the results in atomPairs
  do i=1, arrinfo%nrblocks
    do j=1, arrinfo%ncblocks
      ! Need to set extra if we have an irregularly sized block that needs to
      ! be handled. If we are at the last block and there are extra rows or 
      ! columns we need to set extra to reflect that. Specifically for J
      ! if we have extra columns, and extra was set to 0, then we don't have
      ! extra rows for this block. However, if extra was already set to 1
      ! then we have both extra rows and columns.
      extra = 0
      if ((i == arrinfo%nrblocks) .and. (arrinfo%extraRows>0)) then
        extra = 1
      endif
      if ((j == arrinfo%ncblocks) .and. (arrinfo%extraCols>0)) then
        if (extra == 0) then
          extra = 2
        elseif (extra == 1) then
          extra = 3
        endif
      endif

      ! For localToGlobalMap we need the the starting indices of the block (a,b)
      ! This is just the array index minus 1, times the block size, plus 1
      a = (i-1)*arrinfo%mb+1
      b = (j-1)*arrinfo%nb+1

      ! Now we get the vale indices having to be searched
      call localToGlobalMap(a, b, glo, ghi, arrinfo, blcsinfo, extra)

      ! Now we do the searches for the atoms
      call getAtoms(glo(1), alo(1), firstDim)
      call getAtoms(glo(2), alo(2), secondDim)
      call getAtoms(ghi(1), ahi(1), firstDim)
      call getAtoms(ghi(2), ahi(2), secondDim)
     
      ! Now we need to enumerate the atomPairs and add them to our tree and list
      call addAtomPairRange(alo, ahi, atomPairs, atomTree, i, j)
    enddo
  enddo
  
end subroutine getArrAtomPairs

! This subroutine enumerates a range of atoms and calls addAtomPair to add
! a range of atomPairs to a list
subroutine addAtomPairRange(alo, ahi, atomPairs, atomTree, nrblock, ncblock)
  use O_23bst, only: bst_node

  implicit none

  ! Define passed parameters
  integer, intent(in), dimension(2) :: alo
  integer, intent(in), dimension(2) :: ahi
  type(AtomPair), pointer :: atomPairs
  type(bst_node), pointer :: atomTree
  integer, intent(in) :: nrblock, ncblock

  ! Define local variables
  integer :: i,j  ! loop vars


  do i=alo(1),ahi(1)
    do j=alo(2),ahi(2)
      ! We only need to fill the upper triangle of the vale vale matrix
      ! so we only don't need atom pairs where i>=j
      if (i>=j) then
        call addAtomPair(i,j, atomPairs, atomTree, nrblock, ncblock)
      endif
    enddo
  enddo
end subroutine addAtomPairRange

! This subroutine appends an atom pair to the list of atomPairs
subroutine addAtomPair(i, j, atomPairs, atomTree, nrblock, ncblock)
  use O_23bst, only: bst_node, tree_search, tree_insert
  use O_Parallel, only: modifiedCantor

  implicit none

  ! Define passed parameters
  integer, intent(in) :: i,j
  type(AtomPair), pointer:: atomPairs
  type(bst_node), pointer :: atomTree
  integer, intent(in) :: nrblock, ncblock

  ! Define local variables
  logical :: exists
  integer :: cantorVal
  type(AtomPair), pointer :: newPair
  type(AtomPair), pointer :: lastPair
  type(bst_node), pointer :: targetPair

  call modifiedCantor(i,j,cantorVal)

  call tree_search(atomTree, cantorVal, exists, targetPair)
  ! If the value is already in the tree, add the lock block information to the
  ! tree, then return.
  if (exists) then
    call addBlockToTreeNode(targetPair, 
    return
  endif

  ! If the value is not in the tree, then we first need to add it to the tree
  call tree_insert(atomTree, cantorVal)

  ! Now we need to check and see if this is our first atomPair in the linked
  ! list.
  if ( (atomPairs%I == -1) .and. (atomPairs%J == -1) ) then
    atomPairs%I = i
    atomPairs%J = j
  else 
    ! If we are not the first pair in the list, we need to initialize a new one
    call initAtomPair(newPair)
    newPair%I = i
    newPair%J = i

    ! Now we need to go to the end of the linked list
    lastPair => atomPairs
    do while ( associated(lastPair%next) )
      lastPair => lastPair%next
    enddo

    ! Point the last pair's pointer to the newPair 
    lastPair%next => newPair
  endif

end subroutine addAtomPair

! Given an index gIndx from 1 to valedim/coredim. Return the atom indx needed
! so that the local array can be filled. whichCumul indicates whether we 
! should search through vale or core states.
subroutine getAtoms(gIndx, atomIndx, whichCumul)

  use O_AtomicSites, only: atomSites, numAtomSites

  implicit none

  ! Define Passed Parameters
  integer, intent(in) :: gIndx ! Global Area of interest
  integer, intent(out) :: atomIndx ! What we're looking for
  integer :: whichCumul ! Which cumul states we should look for 0=vale, 1=core
  integer :: i

  if (whichCumul == 0) then
    if ( gIndx <= atomSites(1)%cumulValeStates ) then
      atomIndx = 1
    else
      do i=1, numAtomSites-1
        if (gIndx >= atomSites(i)%cumulValeStates .and. &
              & gIndx < atomSites(i+1)%cumulValeStates) then
          atomIndx = i 
        endif
      enddo
    endif

  else ! whichCumul = 1
    if ( gIndx <= atomSites(1)%cumulCoreStates ) then
      atomIndx = 1
    else
      do i=1, numAtomSites-1
        if (gIndx >= atomSites(i)%cumulCoreStates .and. &
              & gIndx < atomSites(i+1)%cumulCoreStates) then
          atomIndx = i 
        endif
      enddo
    endif
  endif
end subroutine getAtoms

! This subroutine is responsible for checking the ranges of atoms that need
! to be calculated in order to fill the local matrices, to make sure that no
! extra work will be done. It'll then copy the atomPairs structure into a 
! temporary structure, so it can reallocate the atomPairs and add the necessary
! atom ranges.
!subroutine testAtomDupe(vlo,vhi,clo,chi,cvlo,cvhi, atomPairs)
!  implicit none
!
!  ! Define passed parameters
!  integer, intent(in), dimension(2) :: alo, ahi
!  type(AtomPair),  intent(inout), allocatable, dimension(:) :: atomPairs
!
!  ! Define local variables
!  integer, dimension(2) :: loOvlp1, hiOvlp1, loOvlp2, hiOvlp2
!  integer :: i,j
!  
!  ! Initialize intersections 
!  loOvlp1(:) = 0
!  hiOvlp1(:) = 0
!  loOvlp2(:) = 0
!  hiOvlp2(:) = 0
!
!  ! First add all vale atom pairs
!  do i=vlo(1), vhi(1)
!    do j=vlo(2), vhi(2)
!      call addAtomPairs(i,j, atomPairs)
!    enddo
!  enddo
!
!  ! Check if vale and core overlap
!  if ( checkRectOverlap(vlo, vhi, clo, chi) ) then
!    call getOverlapRect(vlo, vhi, clo, chi, loOvlp1, hiOvlp1)
!  endif 
!
!  ! Now add the core pairs
!  do i=clo(1), clo(2)
!    do j=chi(1), chi(2)
!      if (.not. inRectangle(i,j, loOvlp1, hiOvlp1) then
!        addAtomPairs(i,j, atomPairs)
!      endif
!    enddo
!  enddo
!
!  ! Now in order to add the valeCore we need the overlaps of core with coreVale
!  ! and the overlap of vale with coreVale.
!  loOvlp1(:) = 0
!  hiOvlp1(:) = 0
!  ! Check if vale and coreVale overlap
!  if ( checkRectOverlap(vlo, vhi, cvlo, cvhi) ) then
!    call getOverlapRect(vlo, vhi, cvlo, cvhi, loOvlp1, hiOvlp1)
!  endif
!
!  ! Check of core and coreVale overlap
!  if ( checkRectOverlap(clo, chi, cvlo, cvhi) ) then
!    call getOverlapRect(clo, chi, cvlo, cvhi, loOvlp2, hiOvlp2)
!  endif
!  
!  ! Now we add the valeCore pairs
!  do i=vclo(1), vclo(2)
!    do j=vchi(1), vchi(2)
!      if ( .not. (inRectangle(i,j, loOvlp1, hiOvlp1) .or.
!           inRectangle(i,j, loOvlp2, hiOvlp2)) )then
!        addAtomPairs(i,j, atomPairs)
!      endif
!    enddo
!  enddo
!end subroutine testAtomDupe

end module O_ParallelSetup
