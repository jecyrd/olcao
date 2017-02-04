program setupWrapper

  use O_23bst


  type(bst_node), pointer :: tree

  integer :: i

  print *, "init"
  call flush(20)
  call tree_init(tree, 0)
  
  do i=0,(10000)
!    print *, ""
!    print *, ""
!    print *, ""
!    print *, ""
!    print *, "insert",i
!    call flush(20)
    call tree_insert(tree, i)
    
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
  !call tree_printTree(tree)
end program setupWrapper
