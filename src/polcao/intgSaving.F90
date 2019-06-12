module O_IntgSaving

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine multWithBasisFn1 (currentBasisFns,pairXBasisFn2,pairXBasisFn12,&
      & currentlmIndex,currentNumTotalStates,maxAlpha1Used)

   ! Import the necessary modules
   use O_Kinds

   ! Import the necessary data modules
   use O_AtomicTypes, only: maxNumAtomAlphas, maxNumStates

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define variables passed to this subroutine
   real (kind=double), dimension (maxNumAtomAlphas,&
         & maxNumStates,2), intent(in) :: currentBasisFns
   real (kind=double), dimension (16,maxNumAtomAlphas,&
         & maxNumStates), intent(in) :: pairXBasisFn2
   real (kind=double), dimension (maxNumStates,maxNumStates), &
         & intent(inout) :: pairXBasisFn12
   integer, dimension (maxNumStates,2), intent(in)  :: currentlmIndex
   integer, dimension(2), intent(in) :: currentNumTotalStates
   integer, intent(in) :: maxAlpha1Used

   ! Define local variables
   integer :: l,m

   do l = 1, currentNumTotalStates(2)
      do m = 1, currentNumTotalStates(1)
         pairXBasisFn12(m,l) = &
               & sum(currentBasisFns(:maxAlpha1Used,m,1) * &
               & pairXBasisFn2(currentlmIndex(m,1),:maxAlpha1Used,l))
      enddo
   enddo

end subroutine multWithBasisFn1


#ifndef GAMMA


subroutine applyPhaseFactors (currentPair,pairXBasisFn12,statesDim1,statesDim2,&
      & k,runCode,currentKPoint)

   ! Import the necessary modules
   use O_Kinds

   ! Import the necessary data modules
   use O_KPoints, only: numKPoints, phaseFactor
   use O_AtomicTypes, only: maxNumStates

   ! Make sure no funny variables are defined accidentally.
   implicit none

   ! Define the variables passed to this subroutine
   integer, intent(in) :: statesDim1
   integer, intent(in) :: statesDim2
   complex (kind=double), dimension (maxNumStates,maxNumStates,&
         & numKPoints), intent(inout) :: currentPair
   real (kind=double), dimension (statesDim1,statesDim2), &
         & intent(in) :: pairXBasisFn12
   integer, intent(in) :: k ! Cell loop index
   integer, intent(in) :: runCode
   integer, intent(in) :: currentKPoint  ! Only used for runCode>0

   ! Define the local variables
   integer :: l,m,n

   if (runCode == 0) then
      do l = 1, numKPoints
         do m = 1, statesDim2
            do n = 1, statesDim1
               currentPair(n,m,l) = currentPair(n,m,l) + &
                     & phaseFactor(l,k) * pairXBasisFn12(n,m)
            enddo
         enddo
      enddo
   elseif (runCode <= 2) then
      do m = 1, statesDim2
         do n = 1, statesDim1
            currentPair(n,m,1) = currentPair(n,m,1) + &
                  & phaseFactor(currentKPoint,k) * pairXBasisFn12(n,m)
         enddo
      enddo
   else  ! Include a -i factor.
      do m = 1, statesDim2
         do n = 1, statesDim1
            currentPair(n,m,1) = currentPair(n,m,1) + cmplx(0.0_double,&
                  & -1.0_double) * phaseFactor(currentKPoint,k) * &
                  & pairXBasisFn12(n,m)
         enddo
      enddo
   endif

end subroutine applyPhaseFactors


subroutine kPointLatticeOriginShift (currentNumTotalStates,currentPair,&
      & latticeVector,kPointCount,kPointIndex)

   ! Import the necessary modules
   use O_Kinds
   use O_Parallel

   ! Import the necessary data modules
   use O_KPoints, only: kPoints, phaseFactor
   use O_AtomicTypes, only: maxNumStates

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define variables passed to this subroutine
   integer, dimension(2), intent(in) :: currentNumTotalStates
   integer, intent(in) :: KPointCount
   complex (kind=double), dimension (maxNumStates,maxNumStates,&
         & kPointCount), intent(inout) :: currentPair
   real (kind=double), dimension (3), intent(in) :: latticeVector
   integer, intent(in) :: kPointIndex

   ! Define local variables for loop control.
   integer :: i,j,k

   ! Lattice origin shift kpoint effect variables
   complex (kind=double) :: dotProduct
   real (kind=double) :: kPointShiftDot

   do i = 1, kPointCount

      ! Find the kpoint, lattice vector dot product   
      kPointShiftDot = sum(kPoints(:,kPointIndex+i) * latticeVector(:))

      ! Check that the vectors are not perpendicular.
      if (kPointShiftDot .ne. 0.0_double) then

         ! Get the cosine and sine of the dot product
         dotProduct = cmplx(cos(kPointShiftDot),sin(kPointShiftDot),double)

         ! Loop over both atom wave functions.
         do j = 1, currentNumTotalStates(2)
            do k = 1, currentNumTotalStates(1)
               currentPair(k,j,i) = dotProduct * currentPair(k,j,i)
            enddo
         enddo
      endif
   enddo

end subroutine kPointLatticeOriginShift


subroutine saveCurrentPair (i,j,kPointCount,currentPair, blcsinfo, vvInfo, &
      & ccInfo, cvInfo)

   ! Import the necessary modules
   use O_Kinds
   
   use O_Parallel

   ! Import the necessary data modules
   use O_KPoints, only: numKPoints
   use O_AtomicTypes, only: maxNumStates, atomTypes
   use O_AtomicSites, only: valeDim, coreDim, atomSites

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the variables passed to this subroutine.
   integer, intent(in) :: i,j  ! Atom1 and Atom2
   integer, intent(in) :: kPointCount
   complex (kind=double), dimension (maxNumStates,maxNumStates,&
         & numKPoints), intent(inout) :: currentPair
   type(BlacsInfo), intent(in) :: blcsinfo
   type(ArrayInfo), intent(inout) :: vvInfo, ccInfo, cvInfo

   ! Define the matrix necessary for quickly accessing the complex conjugate
   !   of the currentPair.
   complex (kind=double), allocatable, dimension (:,:,:) :: currentPairDagger
   integer :: kPointCounter


   ! Index variables for tracking the locations in the large matrices where
   !   the accumulated results should be stored.
   integer, dimension (2) :: valeStateIndex
   integer, dimension (2) :: valeStateNum
   integer, dimension (2) :: coreStateIndex
   integer, dimension (2) :: coreStateNum

   ! Other local variables
   integer, dimension (2) :: currentAtomType

   ! Assign the current atom types
   currentAtomType(1) = atomSites(i)%atomTypeAssn
   currentAtomType(2) = atomSites(j)%atomTypeAssn


   ! Create the complex conjugate of the relevant currentPair section.
   allocate (currentPairDagger(maxNumStates,maxNumStates,kPointCount))
   do kPointCounter = 1, kPointCount
      currentPairDagger(:,:,kPointCounter) = &
            & transpose(conjg(currentPair(:,:,kPointCounter)))
   enddo

   ! Get the indices for where the states for these atoms should begin to
   !   be recorded in the core containing matrices.  Also get the number
   !   of states that these atoms occupy.
   valeStateIndex(1) = atomSites(i)%cumulValeStates
   valeStateIndex(2) = atomSites(j)%cumulValeStates
   valeStateNum(1)   = atomTypes(currentAtomType(1))%numValeStates
   valeStateNum(2)   = atomTypes(currentAtomType(2))%numValeStates

   ! Save the valence valence part.  The upper triangle contains the real
   !   part while the strict lower triangle contains the imaginary part.
   call atomAtomSaving (valeStateIndex, valeStateNum, &
      & currentPair(1:valeStateNum(1),1:valeStateNum(2),:), &
      & vvInfo, blcsinfo) 

   ! If there is no core contribution then we don't have to save any
   !   core parts.
   if (coreDim == 0) return

   ! Get the indices for where the states for these atoms should begin to
   !   be recorded in the core containing matrices.  Also get the number
   !   of states that these atoms occupy.
   coreStateIndex(1) = atomSites(i)%cumulCoreStates
   coreStateIndex(2) = atomSites(j)%cumulCoreStates
   coreStateNum(1)   = atomTypes(currentAtomType(1))%numCoreStates
   coreStateNum(2)   = atomTypes(currentAtomType(2))%numCoreStates

   ! Save the overlap from the core of atom 1 with the valence of atom 2.
   !   This code is basically a replication of the old code.
   call atomAtomSaving ( (/coreStateIndex(1),valeStateIndex(2)/), &
      & (/coreStateNum(1),valeStateNum(2)/), &
      & currentPair(valeStateNum(1)+1:valeStateNum(1)+coreStateNum(1), &
      & 1:valeStateNum(2),:), &
      & cvInfo, blcsinfo) 

   ! Now we do the Core Core part just the same as we did for the Vale
   !   Vale part.  The only difference is the shift of valeStateNum(:)
   !   in determining the start and end points of the matrix to copy.

   if ((coreStateNum(1) == 0) .or. (coreStateNum(2) == 0)) return
   
   call atomAtomSaving (coreStateIndex, coreStateNum, &
      & currentPair(valeStateNum(1)+1:valeStateNum(1)+coreStateNum(1), &
      & valeStateNum(2)+1:valeStateNum(2)+coreStateNum(2),:), &
      & ccInfo, blcsinfo) 


   ! Deallocate the dagger of the current pair since it is no longer needed.
   deallocate (currentPairDagger)

end subroutine saveCurrentPair

! We're going to do this potentially one of the worst scaling ways. But we can
! optimize later if it's not working well for us.
!
! For each element in the currentPair, we will calculate the index the, process
! coordinates, local block and, index of the local array that the element should
! be at. Then we will store it in the local array, and disregard any elements
! that don't belong specifically to this process.
!
! Our first global index is given by stateIndex, and the number of states is
! given by stateNum. This with the data contained in blcsinfo, and localArr
! should be enough information to calculate the location of each element.
subroutine atomAtomSaving(stateIndex, stateNum, cpair, localArr, blcsinfo)

   ! Import necessary modules
   use O_Kinds

   use O_Parallel
   use O_ParallelSetup
   
   implicit none
   !currentPair(1:valeStateNum(1),1:valeStateNum(2),i),&
   
   ! Define passed parameters
   integer, dimension (2), intent(in) :: stateIndex
   integer, dimension (2), intent(in) :: stateNum
   complex (kind=double), dimension(:,:,:), intent(in) :: cpair
   type(ArrayInfo), intent(inout) :: localArr
   type(BlacsInfo), intent(in) :: blcsinfo

   ! Define local variables
   integer :: i,j
   integer :: pr, pc
   integer :: l,m
   integer :: x,y
   integer :: a,b

   do i=1,size(cpair,1)
      do j=1,size(cpair,2)
         ! Calculate the pr and pc of the element, and make sure it belongs
         ! to this process.
         pr = mod((stateIndex(1)+i-1)/localArr%mb, blcsinfo%prows)
         pc = mod((stateIndex(2)+j-1)/localArr%nb, blcsinfo%pcols)
         if ( (pr == blcsinfo%myprow) .and. (pc == blcsinfo%mypcol) ) then
            ! Calculate the local block(l,m), position in local block (x,y), 
            ! and then the position in the local matrix (a,b). Then store the 
            ! element.
            l = (stateIndex(1)+i-1)/(blcsinfo%prows * localArr%mb)
            m = (stateIndex(2)+j-1)/(blcsinfo%pcols * localArr%nb)
            x = mod(stateIndex(1)+i-1,localArr%mb)+1
            y = mod(stateIndex(2)+j-1,localArr%nb)+1
            
            a = l*localArr%mb + x
            b = m*localArr%nb + y

            localArr%local(a,b,:) = cpair(i,j,:)
         endif
      enddo
   enddo

end subroutine atomAtomSaving

#else

subroutine applyPhaseFactorsGamma (currentPairGamma,pairXBasisFn12,statesDim1,&
      & statesDim2,runCode)

   ! Import the necessary modules
   use O_Kinds

   ! Import the necessary data modules
   use O_AtomicTypes, only: maxNumStates

   ! Make sure no funny variables are defined
   implicit none

   ! Define the variables passed to this subroutine
   integer, intent(in) :: statesDim1
   integer, intent(in) :: statesDim2
   real (kind=double), dimension (maxNumStates,maxNumStates), &
         & intent(inout) :: currentPairGamma
   real (kind=double), dimension (statesDim1,statesDim2), &
         & intent(in) :: pairXBasisFn12
   integer, intent(in) :: runCode

   if (runCode <= 2) then
      currentPairGamma(:statesDim1,:statesDim2) = &
            & currentPairGamma(:statesDim1,:statesDim2) + &
            & pairXBasisFn12(:statesDim1,:statesDim2) 
   else  ! Inlude a -i factor and store the imaginary part (real=0 now).
      currentPairGamma(:statesDim1,:statesDim2) = &
            & currentPairGamma(:statesDim1,:statesDim2) - &
            & pairXBasisFn12(:statesDim1,:statesDim2)
   endif


end subroutine applyPhaseFactorsGamma

subroutine saveCurrentPairGamma (i,j,currentPairGamma,descriptVV,descriptCC,&
      & descriptCV,descriptVC,localVV,localCC,localCV,localVC,&
      & currentNumTotalStates)

   ! Import the necessary modules
   use O_Kinds

   ! Import the necessary data modules
   use O_AtomicTypes, only: maxNumStates, atomTypes
   use O_AtomicSites, only: coreDim, valeDim, atomSites

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the variables passed to this subroutine.
   integer, intent(in) :: i,j ! Atom1 and Atom2
   real (kind=double), dimension (maxNumStates,maxNumStates), &
         & intent(in) :: currentPairGamma
   integer, dimension(9), intent(in) :: descriptVV
   integer, dimension(9), intent(in) :: descriptCC
   integer, dimension(9), intent(in) :: descriptCV
   integer, dimension(9), intent(in) :: descriptVC
   real (kind=double), dimension (:,:,:), intent(inout) :: localVV
   real (kind=double), dimension (:,:,:), intent(inout) :: localCC
   real (kind=double), dimension (:,:,:), intent(inout) :: localCV
   real (kind=double), dimension (:,:,:), intent(inout) :: localVC
   integer, dimension(2), intent(in) :: currentNumTotalStates

   ! Define the matrix for quickly accessing the transpose of the currentPair.
   real (kind=double), allocatable, dimension (:,:) :: currentPairGammaTranspose

   ! Index variables for tracking the locations in the large matrices where
   !   the accumulated results should be stored.
   integer, dimension (2) :: valeStateIndex
   integer, dimension (2) :: valeStateNum
   integer, dimension (2) :: coreStateIndex
   integer, dimension (2) :: coreStateNum

   ! Other local variables
   integer, dimension (2) :: currentAtomType

   ! Assign the current atom types
   currentAtomType(1) = atomSites(i)%atomTypeAssn
   currentAtomType(2) = atomSites(j)%atomTypeAssn

   ! Create the transpose of the relevant currentPairGamma section.
   allocate (currentPairGammaTranspose(maxNumStates,maxNumStates))
   currentPairGammaTranspose = transpose(currentPairGamma)

   ! Get the indices for where the states for these atoms should begin to
   !   be recorded in the core containing matrices.  Also get the number
   !   of states that these atoms occupy.
   valeStateIndex(1) = atomSites(i)%cumulValeStates
   valeStateIndex(2) = atomSites(j)%cumulValeStates
   valeStateNum(1)   = atomTypes(currentAtomType(1))%numValeStates
   valeStateNum(2)   = atomTypes(currentAtomType(2))%numValeStates

   ! Save the valence valence part.  The upper triangle contains the real
   !   part while the strict lower triangle contains the imaginary part.
   call valeValeSavingGamma (i,j,valeStateIndex,valeStateNum,&
         & currentPairGamma,currentPairGammaTranspose,descriptVV,localVV)

   ! If there is no core contribution then we don't have to save any
   !   core parts.
   if (coreDim == 0) return

   ! Get the indices for where the states for these atoms should begin to
   !   be recorded in the core containing matrices.  Also get the number
   !   of states that these atoms occupy.
   coreStateIndex(1) = atomSites(i)%cumulCoreStates
   coreStateIndex(2) = atomSites(j)%cumulCoreStates
   coreStateNum(1)   = atomTypes(currentAtomType(1))%numCoreStates
   coreStateNum(2)   = atomTypes(currentAtomType(2))%numCoreStates

   ! Save the overlap from the core of atom 1 with the valence of atom 2.
   !   This code is basically a replication of the old code.
   call coreValeSavingGamma (valeStateIndex,valeStateNum,&
         & coreStateIndex,coreStateNum,currentPairGamma,&
         & currentPairGammaTranspose,descriptCV,descriptVC,localCV,localVC)

   ! Now we do the Core Core part just the same as we did for the Vale
   !   Vale part.  The only difference is the shift of valeStateNum(:)
   !   in determining the start and end points of the matrix to copy.

   if ((coreStateNum(1) == 0) .or. (coreStateNum(2) == 0)) return

   call coreCoreSavingGamma (i,j,valeStateNum,coreStateIndex,&
         & coreStateNum,currentPairGamma,currentPairGammaTranspose,descriptCC,&
         & localCC)

   ! Deallocate the transpose of the current pair since it is no longer needed.
   deallocate (currentPairGammaTranspose)

end subroutine saveCurrentPairGamma



subroutine valeValeSavingGamma (atom1,atom2,valeStateIndex,valeStateNum,&
      & currentPairGamma,currentPairGammaTranspose,descriptVV,localVV)


   ! Import the necessary modules
   use O_Kinds

   ! Import the necessary data modules
   use O_AtomicTypes, only: maxNumStates
   use O_AtomicSites, only: valeDim

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define variables passed to this subroutine
   integer, intent(in) :: atom1
   integer, intent(in) :: atom2
   integer, dimension (2), intent(in) :: valeStateIndex
   integer, dimension (2), intent(in) :: valeStateNum
   real (kind=double), dimension (maxNumStates,&
         & maxNumStates), intent(in) :: currentPairGamma
   real (kind=double), dimension (maxNumStates,&
         & maxNumStates), intent(in) :: currentPairGammaTranspose
   integer, dimension(9), intent(in) :: descriptVV
   real (kind=double), dimension (:,:,:), intent(inout) :: localVV

   ! Define local variables for loop control and tracking.
   integer :: i,j
   integer, dimension(2) :: startIndices
   integer, dimension(2) :: endIndices
   integer, dimension(2) :: transStart
   integer, dimension(2) :: transEnd
   integer, dimension(2) :: ldReal
   integer, dimension(2) :: ldRealTrans

   startIndices(1) = valeStateIndex(1)+1
   startIndices(2) = valeStateIndex(2)+1
   endIndices(1) = valeStateIndex(1)+valeStateNum(1)
   endIndices(2) = valeStateIndex(2)+valeStateNum(2)
   transStart(1) = startIndices(2)
   transStart(2) = startIndices(1)
   transEnd(1) = endIndices(2)
   transEnd(2) = endIndices(1)
   ldReal(1) = valeStateNum(1)
   ldReal(2) = valeStateNum(2)
   ldRealTrans(1) = valeStateNum(2)
   ldRealTrans(2) = valeStateNum(1)

   if (atom1 == atom2) then
      ! Plop the full matrix on the diagonal.
       call ga_put(ga_vv(1), startIndices(1),endIndices(1),&
             & startIndices(2), endIndices(2), &
             & currentPairGamma(1:valeStateNum(1),1:valeStateNum(2)),&
             & ldReal(1))
   else
      ! First save to top half
      call ga_put(ga_vv(1),startIndices(1),endIndices(1), &
            & startIndices(2),endIndices(2), &
            & currentPairGamma(1:valeStateNum(1),1:valeStateNum(2),i),&
            & ldReal(1))

      ! Put the transpose into the lower half
      call ga_put(ga_vv(1),transStart(1),transEnd(1), &
            & transStart(2),transEnd(2), &
            & currentPairGammaTranspose(1:valeStateNum(2),1:valeStatenum(1)),&
            & ldRealTrans(1))
   endif
end subroutine valeValeSavingGamma



subroutine coreValeSavingGamma (valeStateIndex,valeStateNum,&
      & coreStateIndex,coreStateNum,currentPairGamma,&
      & currentPairGammaTranspose,descriptCV,descriptVC,localCV,localVC)

   ! Import the necessary modules
   use O_Kinds

   ! Import the necessary data modules
   use O_AtomicTypes, only: maxNumStates
   use O_AtomicSites, only: coreDim, valeDim

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define variables passed to this subroutine
   integer, dimension (2), intent(in) :: valeStateIndex
   integer, dimension (2), intent(in) :: valeStateNum
   integer, dimension (2), intent(in) :: coreStateIndex
   integer, dimension (2), intent(in) :: coreStateNum
   real (kind=double), dimension (maxNumStates,&
         & maxNumStates), intent(in) :: currentPairGamma
   real (kind=double), dimension (maxNumStates,&
         & maxNumStates), intent(in) :: currentPairGammaTranspose
   integer, dimension(9), intent(in) :: descriptCV
   integer, dimension(9), intent(in) :: descriptVC
   real (kind=double), dimension (:,:,:), intent(inout) :: localCV
   real (kind=double), dimension (:,:,:), intent(inout) :: localVC

   ! Define local variables for loop control.
   integer :: i
   integer, dimension(2) :: startIndices
   integer, dimension(2) :: endIndices
   integer, dimension(2) :: transStart
   integer, dimension(2) :: transEnd
   integer, dimension(2) :: ldReal
   integer, dimension(2) :: ldRealTrans

   startIndices(1) = coreStateIndex(1)+1
   startIndices(2) = valeStateIndex(2)+1

   endIndices(1) = coreStateIndex(1) + coreStateNum(1)
   endIndices(2) = valeStateIndex(2) + valeStateNum(2)

   transStart(1) = coreStateIndex(2)+1
   transStart(2) = valeStateIndex(1)+1

   transEnd(1) = coreStateIndex(2)+coreStateNum(2)
   transEnd(2) = valeStateIndex(1)+valeStateNum(1)

   ldReal(1) = endIndices(1)-startIndices(1)+1
   ldReal(2) = endIndices(2)-startIndices(2)+1

   ldRealTrans(1) = transEnd(1)-transStart(1)+1
   ldRealTrans(2) = transEnd(2)-transStart(2)+1

   if (coreStateNum(1) .ne. 0) then
      call ga_put(ga_cv(1),startIndices(1),endIndices(1),&
            & startIndices(2), endIndices(2),&
            & currentPair(valeStateNum(1)+1:valeStateNum(1)+coreStateNum(1),&
            & 1:valeStateNum(2),i),ldReal(1))
   endif
   if (coreStateNum(2) .ne. 0) then
      call ga_put(ga_cv(1),transStart(1),transEnd(1),&
            & transStart(2), transEnd(2), &
            & currentPairGammaTranspose(valeStateNum(2)+1:valeStateNum(2)+&
            & coreStateNum(2),1:valeStateNum(1)), ldRealTrans(1))
   endif

end subroutine coreValeSavingGamma



subroutine coreCoreSavingGamma (atom1,atom2,valeStateNum,coreStateIndex,&
      & coreStateNum,currentPairGamma,currentPairGammaTranspose,descriptCC,
      & localCC)

   ! Import the necessary modules
   use O_Kinds

   ! Import the necessary data modules
   use O_AtomicTypes, only: maxNumStates
   use O_AtomicSites, only: coreDim

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define variables passed to this subroutine
   integer, intent(in) :: atom1
   integer, intent(in) :: atom2
   integer, dimension (2), intent(in) :: valeStateNum
   integer, dimension (2), intent(in) :: coreStateIndex
   integer, dimension (2), intent(in) :: coreStateNum
   real (kind=double), dimension (maxNumStates,&
         & maxNumStates), intent(in) :: currentPairGamma
   real (kind=double), dimension (maxNumStates,&
         & maxNumStates), intent(in) :: currentPairGammaTranspose
   integer, dimension(9), intent(in) :: descriptCC
   real (kind=double), dimension (:,:,:), intent(inout) :: localCC

   ! Define local variables.
   integer, dimension(2) :: startIndices
   integer, dimension(2) :: endIndices
   integer, dimension(2) :: transStart
   integer, dimension(2) :: transEnd
   integer, dimension(2) :: ldReal
   integer, dimension(2) :: ldRealTrans

   call ga_put(ga_cc(1), startIndices(1),endIndices(1),&
         & startIndices(2), endIndices(2), &
         & currentPair(valeStateNum(1)+1:valeStateNum(1)+coreStateNum(1),&
         & valeStateNum(2)+1:valeStateNum(2)+coreStateNum(2)),ldReal(1))
   if (atom1 .ne. atom2) then
      call ga_put(ga_cc(1),transStart(1),transEnd(1),&
            & transStart(2),transEnd(2), &
            & currentPairGammaTranspose(valeStateNum(2)+1:valeStateNum(2)+&
            & coreStateNum(2),valeStateNum(1)+1:valeStateNum(1)+&
            & coreStateNum(1)),ldRealTrans(1))
   endif

end subroutine coreCoreSavingGamma

#endif

end module O_IntgSaving
