program gaBlockCyclic

  use MPI
  use O_SlSubs

  implicit none

#include "mafdecls.fh"
#include "global.fh"

  integer :: heap, stack
  integer :: ga_a, ga_s, ga_b
  integer :: mpiRank, mpiSize, mpierr
  integer :: valeDim
  integer :: gaStat

  double complex, dimension(2,2) :: a,s,b,evals

  valeDim = 1200
  
  heap = 300000
  stack = 300000

  call mpi_init(mpierr)

  call ga_initialize()

  gastat = ma_init(MT_F_DCPL, stack, heap)


  ! create array handles
  ga_a = ga_create_handle()
  ga_b = ga_create_handle()
  ga_s = ga_create_handle()

  ! set array name
!  call ga_set_array_name(ga_cc, "coreCore")
!  call ga_set_array_name(ga_vc, "valeCore")

  ! set data properties
  call ga_set_data(ga_a, 2, &
    & (/20,20/), MT_F_DCPL)
  call ga_set_data(ga_s, 2, &
    & (/20,20/), MT_F_DCPL)
  call ga_set_data(ga_b, 2, &
    & (/20,20/), MT_F_DCPL)
 
  ! allocate global array
  gastat = ga_allocate(ga_a)
  gastat = ga_allocate(ga_s)
  gastat = ga_allocate(ga_b)

  call ga_fill(ga_a,cmplx(1.0d0,0.0d0))
  call ga_fill(ga_b,cmplx(0.0d0,0.0d0))

!  call ga_print(ga_a)
!  print *, "filled"
!  call ga_print(ga_vc)
!  call ga_print(ga_cc)
 
!  call ga_zgemm('N','N', valeDim,valeDim,valeDim,1.0d0,ga_vc,ga_cc,1.0d0,ga_vc)
!  print *, "mult"
!  call ga_diag_std(ga_a,ga_b,b)
  call oga_pdzhegv(ga_a, ga_b, evals, 20)

  print *, b
  print *, ''
  print *, ''

  call ga_print(ga_b)

  gastat=ga_destroy(ga_a)
  gastat=ga_destroy(ga_b)
  gastat=ga_destroy(ga_s)


  call ga_terminate()
  call mpi_finalize(mpierr)
end program
