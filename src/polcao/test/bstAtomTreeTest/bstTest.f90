program setupWrapper
  !#use ifport
  use O_bstAtomPair
  
  implicit none

  type(bst_atom_pair_node), pointer :: tree => null()

  integer :: i
  integer(4) :: time_array(3)

  type (tree_vals), pointer :: val => null()
  type (tree_vals), pointer :: found => null()
  logical :: exists

  integer, dimension(9) :: values

  values(1) = 8
  values(2) = 13
  values(3) = 19
  values(4) = 26
  values(5) = 8
  values(6) = 0

  !call itime(time_array)
  !call srand(sum(time_array,1))
  !print *, "init", sum(time_array,1)

  allocate(val)
  val%val = 0
  print *, 'assoc', associated(val%vvblocks)!, val%vvblocks%blockrow
  call flush(6)

  call tree_init(tree,val)

  do i=1,(6)
    nullify(val)
    allocate(val)

    !val%val = mod(irand(),10)
    val%val = values(i)
    print *, ""
    print *, ""
    print *, ""
    print *, ""
    print *, "insert",i, val%val
    call flush(20)
   
    call tree_search(tree, val, exists, found)
    print *, "search",i, val%val, exists, found%val
    call flush(20)
    if (.not. exists) then
      call tree_insert(tree, val)
    else
      deallocate(val)
    endif
    !call tree_printTree(tree)
  enddo
 

  !do i=49,46,-1
  !  print *, ""
  !  print *, ""
  !  print *, ""
  !  print *, ""
  !  print *, "insert",i
  !  call flush(20)
  !  call tree_insert(tree, i)
  !  
  !  call tree_printTree(tree)
  !enddo

  print *, "Printing Tree"
  call tree_printTree(tree)
  call flush(6)

  print *, "Dealloc Tree"
  call flush(6)
  call tree_destroy(tree)
end program setupWrapper