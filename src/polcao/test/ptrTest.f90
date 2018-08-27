module ptrfuncs

contains

subroutine allocptr(p)
  implicit none

  integer, pointer, intent(inout) :: p

  allocate(p)

end subroutine allocptr

subroutine deallocptr(p,tr)
  implicit none

  integer, pointer, intent(inout) :: p
  integer, intent(in) :: tr

  if (associated(p)) then
    print *, "dealloc", tr
    deallocate(p)
  endif

end subroutine deallocptr
end module ptrfuncs

program setupWrapper
  !#use ifport
  use ptrfuncs
  
  implicit none

  integer, pointer :: p1
  integer, pointer :: p2

  integer, target :: i
  
  p1=>null()
  p2=>null()

  call allocptr(p1)

  i=1

  p1=1

  !call allocptr(p2)

  p2 => p1

  nullify(p1)

  print *, "assoc 1", associated(p1)
  print *, "assoc 2", associated(p2)

  call deallocptr(p1,1)
  call deallocptr(p2,2)
end program setupWrapper

