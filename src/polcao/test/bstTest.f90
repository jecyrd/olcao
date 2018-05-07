program setupWrapper
  use ifport
  use O_bstAtomPair
  
  implicit none

  type(bst_atom_pair_node), pointer :: tree

  integer :: i
  integer(4) :: time_array(3)

  type (tree_vals), pointer :: val
  type (tree_vals), pointer :: found
  logical :: exists

  call itime(time_array)
  call srand(sum(time_array,1))
  allocate(val)
  val%val = 0

  print *, "init", sum(time_array,1)
  call flush(20)
  call tree_init(tree,val)

  do i=1,(10)
    nullify(val)
    allocate(val)
    nullify(found)
    allocate(found)

    val%val = mod(irand(),10)
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
end program setupWrapper
