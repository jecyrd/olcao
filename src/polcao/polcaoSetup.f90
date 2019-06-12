! Documentation here
!
module O_ParallelSetup
  use MPI
  implicit none

  ! This is a linked list element
  type AtomPair
    integer :: I
    integer :: J
    type(atomPair), pointer :: next => null()
  end type AtomPair

contains

! This subroutine initializes a new element in the atomPair linked list
subroutine initAtomPair(atomNode,i,j)
  implicit none

  ! Define passed parameters
  type(AtomPair), pointer :: atomNode
  integer, intent(in) :: i,j

  allocate(atomNode)

  atomNode%next => null()
  atomNode%I = i
  atomNode%J = j
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

  if ( ((i>=loOvlp(1)) .and. (i<=hiOvlp(1))) .and. &
    &  ((j>=loOvlp(2)) .and. (j<=hiOvlp(2))) ) then
    inRectangle = 1
    return
  endif
  ! else
  inRectangle = 0
  return
end function inRectangle


subroutine getAtomPairs(vvinfo, ccinfo, cvinfo, blcsinfo, atomPairs, atomTree)
  use O_Parallel, only: ArrayInfo, BlacsInfo
  use O_bstAtomPair

  implicit none
 
  ! Define passed parameters
  type(ArrayInfo), intent(inout) :: vvinfo, ccinfo, cvinfo
  type(BlacsInfo), intent(in) :: blcsinfo
  type(AtomPair), pointer :: atomPairs
  type(bst_atom_pair_node), pointer, intent(inout) :: atomTree 
                                           ! 2-3 tree to keep 
                                           ! us from having duplicate atoms, 
                                           ! and to keep track of the local 
                                           ! blocks the atom pairs belong to

  ! Define local variables
  integer, dimension(2) :: vvlo, vvhi
  integer, dimension(2) :: cclo, cchi
  integer, dimension(2) :: cvlo, cvhi
  type(tree_vals), pointer :: initVal => null()

  integer :: mpierr

  ! Need to figure out a better way to initialize the atomTree but for now
  ! we're just gonna put a zero in it, representing cantor pairing of 0,0
  ! which we should never have anyway.
  allocate(initVal)
  initVal%val = 0
  call tree_init(atomTree, initVal)

  ! We need to call getArrAtomPairs 3 times one for each array
  call getArrAtomPairs(vvinfo, blcsinfo, atomPairs, atomTree, 0)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call getArrAtomPairs(ccinfo, blcsinfo, atomPairs, atomTree, 1)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  call getArrAtomPairs(cvinfo, blcsinfo, atomPairs, atomTree, 2)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

  ! We no longer do this because we need the atomTree for saveCurrentPair
  !  ! By this time we're finished with our atomTree, so we can destroy it.
  !  ! All our atom pair data should be stored in the atomPairs linked list
  !  call tree_destroy(atomTree)

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
  use O_bstAtomPair

  implicit none

  ! Define passed parameters
  type(ArrayInfo), intent(inout) :: arrinfo
  type(BlacsInfo), intent(in) :: blcsinfo
  type(AtomPair),  pointer :: atomPairs
  type(bst_atom_pair_node), pointer :: atomTree
  integer :: whichArr ! whichArr is a control variable to tell the subroutine
                      ! which array it is working with. Options are:
                      ! = 0  valeVale array
                      !   1  coreCore array
                      !   2  coreVale array
  

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
      !a = (i-1)*arrinfo%mb+1
      !b = (j-1)*arrinfo%nb+1
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
      call addAtomPairRange(alo, ahi, atomPairs, atomTree, i, j, whichArr)
    enddo
  enddo
  
end subroutine getArrAtomPairs

! This subroutine enumerates a range of atoms and calls addAtomPair to add
! a range of atomPairs to a list
subroutine addAtomPairRange(alo, ahi, atomPairs, atomTree, nrblock, ncblock, &
                          & whichArr)
  use O_bstAtomPair

  implicit none

  ! Define passed parameters
  integer, intent(in), dimension(2) :: alo
  integer, intent(in), dimension(2) :: ahi
  type(AtomPair), pointer :: atomPairs
  type(bst_atom_pair_node), pointer :: atomTree
  integer, intent(in) :: nrblock, ncblock
  integer, intent(in) :: whichArr ! = 0  valeVale array
                                  !   1  coreCore array
                                  !   2  coreVale array

  ! Define local variables
  integer :: i,j  ! loop vars

  do i=alo(1),ahi(1)
    do j=alo(2),ahi(2)
      ! We only need to fill the upper triangle of the vale vale matrix
      ! so we only don't need atom pairs where i>j
      !if (i<=j) then
        call addAtomPair(i,j, atomPairs, atomTree, nrblock, ncblock, whichArr)
      !endif
    enddo
  enddo
end subroutine addAtomPairRange

! This subroutine appends an atom pair to the list of atomPairs
subroutine addAtomPair(i, j, atomPairs, atomTree, nrblock, ncblock, whichArr)
  use O_bstAtomPair
  use O_Parallel, only: modifiedCantor
  use MPI

  implicit none

  ! Define passed parameters
  integer, intent(in) :: i,j
  type(AtomPair), pointer:: atomPairs
  type(bst_atom_pair_node), pointer :: atomTree
  integer, intent(in) :: nrblock, ncblock
  integer, intent(in) :: whichArr ! = 0  valeVale array
                                  !   1  coreCore array
                                  !   2  coreVale array

  ! Define local variables
  logical :: exists
  integer :: cantorVal
  type(AtomPair), pointer :: newPair => null()
  type(AtomPair), pointer :: lastPair => null()
  type(tree_vals), pointer :: targetPair => null()
  type(tree_vals), pointer :: tempVal => null()

  integer :: mpirank, mpierr

  call modifiedCantor(i,j,cantorVal)
  allocate(tempVal)
  tempVal%val = cantorVal

  call MPI_COMM_RANK(MPI_COMM_WORLD,mpirank,mpierr)

  !call tree_search(atomTree, cantorVal, exists, targetPair)
  call tree_search(atomTree, tempVal, exists, targetPair)

  ! If the value is already in the tree, add the lock block information to the
  ! tree, then return.
  if (exists) then
    ! Deallocate tempVal because we no longer need it
    deallocate(tempVal)
    tempVal => targetPair
    ! Add block to 
    !if (whichArr==0) then
    !  print *, "add existing"
    !  call flush(6)
    !  call tree_addVVBlock(targetpair, nrblock, ncblock)
    !else if (whichArr==1) then
    !  call tree_addCVBlock(targetpair, nrblock, ncblock)
    !else if (whichArr==2) then
    !  call tree_addCCBlock(targetpair, nrblock, ncblock)
    !endif
  else
    ! If the value is not in the tree, then we first need to add it to the tree.
    ! First we need to create a new treevals type.
    call tree_insert(atomTree, tempVal)

    ! Then we initialize lBlockCoords Linked lists for our new val.
    !if (whichArr==0) then
    !  call initBlockCoords(tempVal%vvblocks)
    !else if (whichArr==1) then
    !  call initBlockCoords(tempVal%cvblocks)
    !else if (whichArr==2) then
    !  call initBlockCoords(tempVal%ccblocks)
    !endif

    if (.not. associated(atomPairs)) then
      call initAtomPair(atomPairs,i,j)
    else
      ! If we are not the first pair in the list, 
      ! we need to initialize a new one
      call initAtomPair(newPair,i,j)

      ! Now we need to go to the end of the linked list
      lastPair => atomPairs
      do while ( associated(lastPair%next) )
        lastPair => lastPair%next
      enddo

      ! Point the last pair's pointer to the newPair 
      lastPair%next => newPair
    endif
  endif

  ! Lastly add the block indices to the lists
  if (whichArr==0) then
    call tree_addVVBlock(tempVal, nrblock, ncblock)
  else if (whichArr==1) then
    call tree_addCVBlock(tempVal, nrblock, ncblock)
  else if (whichArr==2) then
    call tree_addCCBlock(tempVal, nrblock, ncblock)
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
!    if ( gIndx <= atomSites(1)%cumulValeStates ) then
!      atomIndx = 1
    if ( gIndx > atomSites(numAtomSites)%cumulValeStates ) then
      atomIndx = numAtomSites
    else
      do i=1, numAtomSites-1
        if (gIndx >= atomSites(i)%cumulValeStates .and. &
              & gIndx <= atomSites(i+1)%cumulValeStates) then
          atomIndx = i 
        endif
      enddo
    endif

  else ! whichCumul = 1
    !if ( gIndx <= atomSites(1)%cumulCoreStates ) then
    !  atomIndx = 1
    if ( gIndx > atomSites(numAtomSites)%cumulCoreStates ) then
      atomIndx = numAtomSites
    else
      do i=1, numAtomSites-1
        if (gIndx >= atomSites(i)%cumulCoreStates .and. &
              & gIndx <= atomSites(i+1)%cumulCoreStates) then
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

! This subroutine handles writing out the vale vale matrix to disk.
! It is only to be done by process 0.
subroutine writeValeVale(arrinfo, blcsinfo, numKPoints, potDim, &
                        & currPotTypeNumber, CurrAlphaNumber, opcode)
  use HDF5

  use O_Kinds
  use O_Parallel, only: ArrayInfo, BlacsInfo, localToGlobalMap

  use O_SetupIntegralsHDF5, only: valeVale_dsid, atomOverlap_did, &
                                & atomKEOverlap_did, atomNucOverlap_did, &
                                & atomPotOverlap_did

  use O_PotTypes, only: potTypes

  implicit none
  
  ! Define passed parameters
  type(ArrayInfo) :: arrinfo
  type(BlacsInfo) :: blcsinfo
  integer, intent(in) :: numKPoints, potDim, currPotTypeNumber, currAlphaNumber
  integer, intent(in) :: opcode

  ! Define local variables
  integer :: i, j, k, l, x, y, a,b, kpl, hdferr
  integer, dimension(2) :: lo,hi

  integer(hid_t) :: memspace_dsid
  integer(hid_t), dimension(numKPoints,potDim) :: datasetToWrite_did
  integer(hsize_t), dimension(2) ::hslabCount, hslabStart

  real(kind=double), allocatable, dimension(:,:,:) :: dataOut
  real(kind=double) :: smallThresh10

  smallThresh10 = real(1.0d-10,double)

  select case(opcode)
  case(1)
    datasetToWrite_did(:,1) = atomOverlap_did(:)
  case(2)
    datasetToWrite_did(:,1) = atomKEOverlap_did(:)
  case(3)
    datasetToWrite_did(:,1) = atomNucOverlap_did(:)
  case(4)
    datasetToWrite_did(:,:) = atomPotOverlap_did(:,:)
  end select

  ! Set the third component to numkp rather than setting in loop every time
  !hslabCount(3) = size(arrinfo%local,3)

  do i=0,arrinfo%nrblocks-1
    do j=0,arrinfo%ncblocks-1
      if ((arrinfo%extraRows>0) .and. (i==arrinfo%nrblocks-1) .and. &
        & (arrinfo%extraCols>0) .and. (j==arrinfo%ncblocks-1)) then
        hslabCount(1) = arrinfo%extraRows
        hslabCount(2) = arrinfo%extraCols
      else if ((arrinfo%extraRows>0) .and. (i==arrinfo%nrblocks-1)) then
        hslabCount(1) = arrinfo%extraRows
        hslabCount(2) = arrinfo%nb
      else if ((arrinfo%extraCols>0) .and. (j==arrinfo%ncblocks-1)) then
        hslabCount(1) = arrinfo%mb
        hslabCount(2) = arrinfo%extraCols
      else
        hslabCount(1) = arrinfo%mb
        hslabCount(2) = arrinfo%nb
      endif

      ! Allocate space for this exact block size
      allocate(dataOut(hslabCount(1),hslabCount(2),size(arrinfo%local,3)))

      a = i*arrinfo%mb
      b = j*arrinfo%nb
      call localToGlobalMap(a,b, lo, hi, arrinfo, blcsinfo, 0)
      hslabStart(1) = lo(1)-1
      hslabStart(2) = lo(2)-1
      !hslabStart(3) = 1

      call h5screate_simple_f(2, hslabCount, memspace_dsid, hdferr)
    
      ! Define hyperslab to be written to
      call h5sselect_hyperslab_f(valeVale_dsid, H5S_SELECT_SET_F, hslabStart, &
        & hslabCount, hdferr)

      ! Need to prepare the data so that complex parts are on saved on the
      ! bottom half of matrix, and real parts on the top half.
      x=1
      do k=hslabStart(1)+1,(hslabStart(1)+hslabCount(1))
        y=1
        do l=hslabStart(2)+1,(hslabStart(2)+hslabCount(2))
          if (l>=k) then ! top half
            dataOut(x,y,:) = real(arrinfo%local(a+x,b+y,:))
          else ! Bottom half
            dataOut(x,y,:) = aimag(arrinfo%local(a+x,b+y,:))
          endif
          y = y + 1
        enddo
        x = x + 1
      enddo

      do kpl=1,numKPoints 
        select case (opCode)
        case(1:3)
          ! write slab to disk
          call h5dwrite_f(datasetToWrite_did(kpl,1), H5T_NATIVE_DOUBLE, &
            & dataOut(:,:,kpl),hslabCount,hdferr, &
            & file_space_id=valeVale_dsid, mem_space_id=memspace_dsid)
        case(4)
          call h5dwrite_f(datasetToWrite_did( &
            & kpl,potTypes(currPotTypeNumber)%cumulAlphaSum+currAlphaNumber), &
            H5T_NATIVE_DOUBLE, dataOut(:,:,kpl), &
            & hslabCount, hdferr, file_space_id=valeVale_dsid, &
            & mem_space_id=memspace_dsid)
        case default
          print *, "Something went very wrong in writeValeVale"
        end select
      enddo
      
      deallocate(dataOut)
    enddo
  enddo


end subroutine writeValeVale

end module O_ParallelSetup
