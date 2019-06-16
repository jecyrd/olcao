module O_Parallel

   ! Import necessasry modules
   use O_Kinds

   implicit none


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

   ! Begin list of subroutines and functions
   contains

! Subroutine for transposing a complex distributed matrix
subroutine pctrans(aInfo, cInfo)
   use pztrancInterface

   implicit none

   ! Define passed parameters
   type(ArrayInfo) :: aInfo, cInfo

   ! Define local variables
   integer :: i

   do i=1,aInfo%numKP
   call pztranc(aInfo%J,aInfo%I,(1.0_double,0.0_double),aInfo%local(:,:,i), &
      & 1,1,aInfo%desc,(1.0_double,0.0_double),cInfo%local(:,:,i),1,1, &
      & cInfo%desc)
   enddo
end subroutine pctrans

! Subroutine for transposing a real distributed matrix
subroutine pdtrans(aInfo, cInfo)
   use pdtranInterface

   implicit none

   ! Define passed parameters
   type(gArrayInfo) :: aInfo, cInfo

   ! Define local variables
   integer :: i

   do i=1,aInfo%numKP
   call pdtran(aInfo%J,aInfo%I,1.0_double,aInfo%local, &
      & aInfo%mb,aInfo%nb,aInfo%desc,1.0_double, cInfo%local, &
      & 0,0,cInfo%desc)
   enddo
end subroutine pdtrans

! This subroutine uses a modified Cantor pairing function to encode 2 natural
! numbers into 1 natural number. Normally the Cantor pairing function is 
! unique to the ordering of a pair. i.e. (1,3) produces a different result
! than (3,1). However in our case, we're concerned about calculating unique 
! atom pairs in integralsSCF. In our case (1,3) is the same as conjugate 
! inverse of (3,1). So we wrap our cantor pairing subroutine with one that 
! insures that k1 is always the smallest and k2 is always the largest. This
! will result in (1,3) and (3,1) to not have unique results.
subroutine modifiedCantor(i, j, res)
   implicit none

   ! Define Passed parameters
   integer, intent(in) :: i,j
   integer, intent(out) :: res

   if (i<=j) then
      call cantorPairing(i,j,res)
   else
      call cantorPairing(j,i,res)
   endif
end subroutine modifiedCantor

! This subroutine executes the Cantor Pairing function. This function is 
! used to encode two natural numbers into a single number. It is defined as:
!           pi(k1, k2) = 1/2 * (k1 + k2) * (k1 + k2 + 1) + k2
subroutine cantorPairing(k1, k2, res)
   implicit none

   ! Define passed parameters
   integer, intent(in) :: k1,k2
   integer, intent(out) :: res

   res = 0.5 * (k1 + k2) * (k1 + k2 + 1) + k2
end subroutine cantorPairing

! This subroutine inverts the Cantor pairing function. This is done by the 
! following process:
! k1 = w - k2
! k2 = z - t
! w = floor( (sqrt(8*z + 1) - 1) / 2 )
! t = (w^2 + w) / 2
! z = the result of the Cantor pairing function
subroutine cantorInverse(z, k1, k2)
   implicit none

   ! Define passed parameters
   integer, intent(in) :: z
   integer, intent(out) :: k1, k2

   ! Define local variables
   integer :: w, t

   w = floor( (sqrt(dble(8*z + 1)) - 1.0) / 2 )
   t = (w**2 + w) / 2

   k2 = z - t
   k1 = w - k2
end subroutine cantorInverse

! This subroutine calculates the process grid dimensions and sets up the 
! blacs context for a distributed operation
subroutine setupBlacs(blcsinfo, pgrid)
   use MPI

   implicit none

   ! Define passed parameters
   type(BlacsInfo), intent(inout) :: blcsinfo
   integer, dimension(2), intent(in) :: pgrid

   ! Define local varaibles
   integer :: mpierr
   integer :: val

   external :: BLACS_GET, BLACS_GRIDINIT, BLACS_GRIDINFO

   call MPI_Comm_rank(MPI_COMM_WORLD, blcsinfo%mpirank, mpierr)
   call MPI_Comm_size(MPI_COMM_WORLD, blcsinfo%mpisize, mpierr)

   ! First we need to calculate the processor grid
   !call calcProcGrid(blcsinfo)
   blcsinfo%prows = pgrid(1)
   blcsinfo%pcols = pgrid(2)

   !call BLACS_GET(blcsinfo%context, 0, val)
   call BLACS_GET(-1,0,blcsinfo%context)

   call BLACS_GRIDINIT(blcsinfo%context, 'r', blcsinfo%prows, blcsinfo%pcols)

   call BLACS_GRIDINFO(blcsinfo%context, blcsinfo%prows, blcsinfo%pcols, &
      & blcsinfo%myprow, blcsinfo%mypcol)

end subroutine setupBlacs

! This subroutine calculates the dimensions of the processor grid given 
! the number of available processes. This subroutine creates a process grid
! that is as close to square as possible.
subroutine calcProcGrid(blcsinfo)
   implicit none

   ! Define passed parameters
   type(BlacsInfo), intent(inout) :: blcsinfo

   if (blcsinfo%mpisize > 1) then
      blcsinfo%prows = int(sqrt(dble(blcsinfo%mpisize)))
      blcsinfo%pcols = blcsinfo%mpisize / blcsinfo%prows
   else
      blcsinfo%prows = 1
      blcsinfo%pcols = 1
   endif
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
   call getBlockDims(arrinfo, blcsinfo)

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
   type(BlacsInfo), intent(in) :: blcsinfo
   type(ArrayInfo), intent(inout) :: arrinfo

   ! Local variables
   integer :: extraRow, extraCol

   ! Calculate the blocking factors to be used. From the calcProcGrid
   ! routine we know that procGridRows will be <= procGridCols.
   ! Because of this we'll do the integer divide on procGridCols, ensuring
   ! that our blocks are square.
   arrinfo%mb =  int(arrinfo%I / blcsinfo%prows)
   arrinfo%nb =  int(arrinfo%J / blcsinfo%pcols)

   ! For the case that global matrix doesn't divide perfectly by our choice in
   ! block size, there will be a couple extra rows and columns that need to be 
   ! asssigned to a process.
   ! The total dim I, divided by the block dimension, tells us how many blocks
   ! are going to be along a dimensions. When we mod this by the rows in the 
   ! process grid, it tells us the row  coordinate in the process grid the 
   ! extra rows belong to. This of course also works for the columns.
   ! The extra rows and columns don't always go to the first row or column
   ! of the process grid, because we have to maintain the block cyclic
   ! distribution. This is why the calculation is done in this way. For more info
   ! see http://netlib.org/scalapack/slug/node78.html
   extraRow = mod( arrinfo%I/arrinfo%mb, blcsinfo%prows)
   extraCol = mod( arrinfo%J/arrinfo%nb, blcsinfo%pcols)

   ! Calculate the number of blocks along rows and columns
   arrinfo%nrblocks = int(arrinfo%I / arrinfo%mb / blcsinfo%prows)
   arrinfo%ncblocks = int(arrinfo%J / arrinfo%nb / blcsinfo%pcols)

   if ( (extraRow>0) .and. (extraRow == blcsinfo%myprow) ) then
      arrinfo%extraRows = mod(arrinfo%I, arrinfo%mb)
      arrinfo%nrblocks = arrinfo%nrblocks+1
   else
      arrinfo%extraRows = 0
   endif
   if ( (extraCol>0) .and. (extraCol == blcsinfo%mypcol) ) then
      arrinfo%extraCols = mod(arrinfo%J, arrinfo%nb)
      arrinfo%ncblocks = arrinfo%ncblocks+1
   else
      arrinfo%extraCols = 0
   endif
end subroutine getblockDims

! This subroutine allocates the proper amount of space needed for each
! processes local array.  The SCALAPACK function NUMROC (number of rows and
! colums) is meant for this purpose.
subroutine allocLocalArray(arrinfo, blcsinfo)
   use O_Kinds

   implicit none

   ! Define passed parameters
   type(BlacsInfo), intent(in) :: blcsinfo
   type(ArrayInfo), intent(inout) :: arrinfo

   ! Define local variables
   integer :: nlrows, nlcols

   ! External Functions
   integer, external :: numroc

   nlrows = numroc(arrinfo%I, arrinfo%mb, blcsinfo%myprow, 0, blcsinfo%prows)
   nlcols = numroc(arrinfo%J, arrinfo%nb, blcsinfo%mypcol, 0, blcsinfo%pcols)

   !#ifndef GAMMA
   allocate(arrinfo%local(nlrows,nlcols,arrInfo%numKP))
   !#else
   !  allocate(arrinfo%local(nlrows,nlcols))
   !#endif

   ! initialize the local array to zero
   !#ifndef GAMMA
   arrinfo%local(:,:,:) = cmplx(0.0_double, 0.0_double)
   !#else
   !  arrinfo%local(:,:) = 0.0_double
   !#endif
end subroutine allocLocalArray

! This subroutine cleans up our local array which was allocated in the routine
! above
subroutine deallocLocalArray(arrinfo)
   implicit none

   ! Define passed parameters
   type(ArrayInfo), intent(inout) :: arrinfo

   deallocate(arrinfo%local)
end subroutine deallocLocalArray

! This subroutine initializes the BLACS array descriptor and 
subroutine getArrDesc(arrinfo, blcsinfo)
   implicit none

   ! Define passed parameters
   type(BlacsInfo), intent(in) :: blcsinfo
   type(ArrayInfo), intent(inout) :: arrinfo

   ! Define local variables
   integer :: info

   call descinit(arrinfo%desc, arrinfo%I, arrinfo%J, arrinfo%mb, arrinfo%nb, &
      & 0, 0, blcsinfo%context, size(arrinfo%local,1), info)

   if (info < 0) then
      print *, "BLACS Routine DESCINIT did not exit successfully. Exiting"
      stop
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
!       (Pr,Pc) are the dimensions of the process grid.
!       (pr,pc) are the coordinates of a process in the grid.
!
! What we need to calculate is 2 sets of coordinates, first the global 
! coordinates of the "top left" element of a local block. Second the global
! coordinates of the "bottom right" element of a local block. This'll allow us
! to tell which elements of the global array we require to fill a local block.
!
! Additionally we can calculate I, and J with the following equations:
!         I = (l*Pr + pr)*MB + x,  J = (m*Pc + pc)*NB + y
! Using the above for l and m we can obtain:
!   I = (((a-x)/MB)*Pr + pr)*MB + x,  J = (((b-y)/NB)*Pc + pc)*NB + y
!
! By setting x and y to 0. We can obtain the "top left" global coordinate. Then
! by setting x and y to (MB-1) and (NB-1) respectively, we obtain the
! "bottom right" coordinate of the global array.
!
! For this to work the passed parameters a,b should be the coordinate in
! the local array corresponding to the "top left" element of some local block.
subroutine localToGlobalMap(a, b, glo, ghi, arrinfo, blcsinfo, extra)
   implicit none

   ! Define passed parameters
   type(BlacsInfo), intent(in) :: blcsinfo
   type(ArrayInfo), intent(in) :: arrinfo
   integer, intent(in) :: a,b
   integer, dimension(2), intent(out) :: glo, ghi ! target global values
   integer, intent(in) :: extra ! Control variable that denotes if we are doing
   ! a irregular size block. Options
   ! = 0 normal block size
   !   1 nrows is less than mb
   !   2 ncols is less than nb
   !   3 both nrows and ncols is less than nb and mb

   ! Calculate the upper left corner
   !glo(1) = (((a-1)/arrinfo%mb)*blcsinfo%prows+blcsinfo%myprow)*arrinfo%mb + 1
   !glo(2) = (((b-1)/arrinfo%nb)*blcsinfo%pcols+blcsinfo%mypcol)*arrinfo%nb + 1
   glo(1) = ltog(a,1,arrinfo%mb,blcsinfo%prows,blcsinfo%myprow)
   glo(2) = ltog(b,1,arrinfo%nb,blcsinfo%pcols,blcsinfo%mypcol)


   ! Calculate the bottom right corner
   if (extra == 0) then
      ghi(1) = glo(1) + arrinfo%mb - 1
      ghi(2) = glo(2) + arrinfo%nb - 1
   elseif (extra == 1) then
      ghi(1) = glo(1) + arrinfo%extraRows - 1
      ghi(2) = glo(2) + arrinfo%nb - 1
   elseif (extra == 2) then
      ghi(1) = glo(1) + arrinfo%mb - 1
      ghi(2) = glo(2) + arrinfo%extraCols - 1
   elseif (extra == 3) then
      ghi(1) = glo(1) + arrinfo%extraRows - 1
      ghi(2) = glo(2) + arrinfo%extraCols - 1
   endif

end subroutine localToGlobalMap 

function ltog(a, x, xb, np, myp)
   implicit none

   ! Define function
   integer :: ltog

   ! Define passed parameters
   integer, intent(in) :: a, x, xb, np, myp

   ltog = (((a-x)/xb)*np+myp)*xb + x
end function ltog

! Refer to documentation above localToGlobalMap. This subroutine calculates
! the local array indices for a section of the global
subroutine globalToLocalMap(glo, ghi, llo, lhi, arrinfo, blcsinfo, extra)
   implicit none

   ! Define passed parameters
   type(BlacsInfo), intent(in) :: blcsinfo
   type(ArrayInfo), intent(in) :: arrinfo
   integer, dimension(2), intent(in) :: glo, ghi 
   integer, dimension(2), intent(out) :: llo, lhi ! target local values
   integer, intent(in) :: extra ! Control variable that denotes if we are doing
   ! a irregular size block. Options
   ! = 0 normal block size
   !   1 nrows is less than mb
   !   2 ncols is less than nb
   !   3 both nrows and ncols is less than nb and mb
   !print *, "glob to loc", glo, ghi, extra, blcsInfo%prows, blcsInfo%pcols
   !call flush(6)

   ! Note, because of integer division the mb and nb terms cannot be 
   ! canceled out.
   ! ((I-1)/(Pr*MB))*MB + mod(I,MB) + 1

   llo(1) = gtol(glo(1), blcsinfo%prows, arrinfo%mb)
   llo(2) = gtol(glo(2), blcsinfo%pcols, arrinfo%nb)
   !llo(1) = ((glo(1)-1)/(blcsInfo%prows*arrinfo%mb))*arrinfo%mb + &
   !  & mod((glo(1)-1),arrinfo%mb) + 1
   !llo(2) = ((glo(2)-1)/(blcsInfo%pcols*arrinfo%nb))*arrinfo%nb + &
   !  & mod((glo(2)-1),arrinfo%nb) + 1
   !llo(1) = ((glo(1))/(blcsInfo%prows)) + mod(glo(1),arrInfo%mb) + 1
   !llo(2) = ((glo(2))/(blcsInfo%pcols)) + mod(glo(2),arrInfo%nb) + 1

   if (extra == 0) then
      !lhi(1) = ((ghi(1))/(blcsInfo%prows))+mod(ghi(1),arrInfo%mb)+1
      !lhi(2) = ((ghi(2))/(blcsInfo%pcols))+mod(ghi(2),arrInfo%nb)+1
      lhi(1) = gtol(ghi(1), blcsinfo%prows, arrinfo%mb)
      lhi(2) = gtol(ghi(2), blcsinfo%pcols, arrinfo%nb)
   elseif (extra == 1) then
      !lhi(1) = llo(1) + arrinfo%extraRows
      !lhi(2) = ((ghi(2))/(blcsInfo%pcols))+mod(ghi(2),arrInfo%nb)+1
      lhi(1) = llo(1) + arrinfo%extraRows
      lhi(2) = gtol(ghi(2), blcsinfo%pcols, arrinfo%nb)
   elseif (extra == 2) then
      !lhi(1) = ((ghi(1))/(blcsInfo%prows))+mod(ghi(1),arrInfo%mb)+1
      !lhi(2) = llo(2) + arrinfo%extraCols
      lhi(1) = gtol(ghi(1), blcsinfo%prows, arrinfo%mb)
      lhi(2) = llo(2) + arrinfo%extraCols
   elseif (extra == 3) then
      lhi(1) = llo(1) + arrinfo%extraRows
      lhi(2) = llo(2) + arrinfo%extraCols
   endif

   !print *, "gtol lo", llo,lhi, arrinfo%mb, arrinfo%nb
   !call flush(6)

end subroutine globalToLocalMap

function gtol(I,P,xb)
   implicit none

   ! Define function
   integer :: gtol

   ! Define passed parameters
   integer, intent(in) :: I, P, xb

   ! Note, because of integer division the mb and nb terms cannot be canceled out
   ! ((I-1)/(Pr*MB))*MB + mod(I,MB) + 1
   gtol = ((I-1)/(P*xb)*xb) + mod((I-1),xb)+1
end function gtol


end module O_Parallel

