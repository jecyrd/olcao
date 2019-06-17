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

function ltog(a, x, xb, np, myp)
   implicit none

   ! Define function
   integer :: ltog

   ! Define passed parameters
   integer, intent(in) :: a, x, xb, np, myp

   ltog = (((a-x)/xb)*np+myp)*xb + x
end function ltog

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

