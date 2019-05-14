program maptest
  use MPI
  use O_Parallel

  implicit none

  type(BlacsInfo) :: blcsinfo
  type(ArrayInfo) :: arrinfo

  character(len=10) :: arg
  integer :: numargs

  integer :: valeDim, coreDim, numkp
  integer :: i,j,a,b
  integer :: mpierr
  integer :: extra

  integer, dimension(2) :: glo, ghi
  integer, dimension(2) :: llo, lhi

  numargs = command_argument_count()

  call get_command_argument(1,arg)
  read (arg,'(I10)') valedim
  call get_command_argument(2,arg)
  read (arg,'(I10)') coredim

  call MPI_INIT(mpierr)

  numkp = 1
  extra = 0

  call setupBlacs(blcsinfo)
  call setupArrayDesc(arrinfo, blcsinfo, valedim, valedim, numkp)

  !print *, blcsinfo%mpirank, blcsinfo%prows, blcsinfo%pcols, arrinfo%mb, arrinfo%nb, blcsinfo%myprow, blcsinfo%mypcol
  !print *, blcsinfo%mpirank, "extra", arrinfo%extrarows, arrinfo%extracols
  !print *, blcsinfo%mpirank, arrinfo%nrblocks, arrinfo%ncblocks
  if (blcsinfo%mpirank == 0) then
    open(6,file='proc0.out',status='unknown')
  else
    open(6,file='proc1.out',status='unknown')
  endif

  do i=1, arrinfo%nrblocks
    glo(:) = i
    ghi(:) = i
    call globalToLocalMap(glo, ghi, llo, lhi, arrinfo, blcsinfo, extra)
    call localToGlobalMap(llo(1),llo(2),glo,ghi,arrinfo,blcsinfo,extra)
    print *, i, llo(1), llo(2), glo,ghi
  enddo

  call MPI_FINALIZE(MPIERR)
end program maptest

subroutine ltog2(a, l, np, myp, nb, x)
  implicit none

  integer :: a

  ! Define passed parameters
  integer :: l, np, myp, nb, x


  a = (l*np + myp)*nb + x
end subroutine ltog2

subroutine globalToLocalMap2(I, P, NB)
  implicit none

  ! Define passed parameters
  integer, intent(in) :: I, P, NB

  ! define local variables
  integer :: lp, l, x, a

  lp = mod((I-1)/NB,P)
  l = (I-1)/(P*nb)
  x = mod((I-1),NB)+1

  a = l*NB+x

  print *, P, NB, I, lp, l, x, a
end subroutine globalToLocalMap2
