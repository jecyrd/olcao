! 2-3 tree implementation specifically for keeping track of atom pairs, and
! the local blocks that corespond to those pairs.
! This is all to reduce the amount of work and complexity when it comes time
! to save the currentPair in the local arrays.
! Author: James Currie
! email: <jecyrd@mail.umkc.edu>
module O_bstAtomPair
  implicit none

  public 
  
  ! For the purposes of setup we need a special bst_atom_pair_node that will 
  ! work with a set of overloaded routines that exist in the reference
  ! implementation in bst23.f90. However we need to have each value, we need
  ! to have a data structure that holds the cantor value, as well as the 
  ! local blocks that correspond to that atom pair. Rather than rewrite
  ! existing types and routines in bst23.f90 we will declare a new type
  ! to operate on. As well as overloaded operators for the comparison
  ! functions on our new types. This'll keep changes small from the reference 
  ! implementation.
  type bst_atom_pair_node
    type(tree_vals), pointer :: lval
    type(tree_vals), pointer :: mval
    type(tree_vals), pointer :: hval
    type(bst_atom_pair_node), pointer :: lchild, mlchild, mrchild, rchild
    type(bst_atom_Pair), pointer :: parent
  end type bst_atom_pair_node

  ! Data structure to replace the integers in the reference implementation of
  ! the bst
  type tree_vals
    integer :: val
    type(lBlockCoords), pointer :: vvblocks
    type(lBlockCoords), pointer :: cvblocks
    type(lBlockCoords), pointer :: ccblocks
  end type tree_vals

  ! New linked list type for readability, to be used with tree above.
  type lBlockCoords
    integer :: blockrow
    integer :: blockcol
    type(lBlockCoords), pointer :: next
  end type lBlockCoords

  ! Define overload for greater than when using tree_vals type
  interface operator(>)
    function bst_gt(tval1, tval2)
      implicit none

      ! Define passed parameters
      type(tree_vals), pointer :: tval1, tval2

      ! Define function return
      logical :: bst_gt

      if (tval1=>val > tval2=>val) then
        return .true.
      else
        return .false.
      endif
    end function bst_gt
    
    function bst_val_gt(tval1, tval2)
      implicit none

      ! Define passed parameters
      type(tree_vals), pointer :: tval1
      integer, intent(in) :: tval2

      ! Define function return
      logical :: bst_gt

      if (tval1=>val > tval2) then
        return .true.
      else
        return .false.
      endif
    end function bst_val_gt
    
    function val_bst_gt(tval1, tval2)
      implicit none

      ! Define passed parameters
      integer, intent(in) :: tval1
      type(tree_vals), pointer :: tval2

      ! Define function return
      logical :: bst_gt

      if (tval1 > tval2=>val) then
        return .true.
      else
        return .false.
      endif
    end function val_bst_gt
  end interface

  ! Define overload for less than when using tree_vals type
  interface operator(<)
    function bst_lt(tval1, tval2)
      implicit none

      ! Define passed parameters
      type(tree_vals), pointer :: tval1, tval2

      ! Define function return
      logical :: bst_gt

      if (tval1=>val < tval2=>val) then
        return .true.
      else
        return .false.
      endif
    end function bst_lt

    function bst_val_lt(tval1, tval2)
      implicit none

      ! Define passed parameters
      type(tree_vals), pointer :: tval1
      integer, intent(in) :: tval2

      ! Define function return
      logical :: bst_lt

      if (tval1=>val < tval2) then
        return .true.
      else
        return .false.
      endif
    end function bst_val_lt
    
    function val_bst_lt(tval1, tval2)
      implicit none

      ! Define passed parameters
      integer, intent(in) :: tval1
      type(tree_vals), pointer :: tval2

      ! Define function return
      logical :: bst_gt

      if (tval1 < tval2=>val) then
        return .true.
      else
        return .false.
      endif
    end function val_bst_lt
  end interface
  
  ! Define overload for equal to when using tree_vals type
  interface operator(==)
    function bst_eq(tval1, tval2)
      implicit none

      ! Define passed parameters
      type(tree_vals), pointer :: tval1, tval2

      ! Define function return
      logical :: bst_gt

      if (tval1=>val == tval2=>val) then
        return .true.
      else
        return .false.
      endif
    end function bst_eq

    function bst_val_eq(tval1, tval2)
      implicit none

      ! Define passed parameters
      type(tree_vals), pointer :: tval1
      integer, intent(in) :: tval2

      ! Define function return
      logical :: bst_gt

      if (tval1=>val == tval2) then
        return .true.
      else
        return .false.
      endif
    end function bst_val_eq
    
    function val_bst_eq(tval1, tval2)
      implicit none

      ! Define passed parameters
      integer, intent(in) :: tval1
      type(tree_vals), pointer :: tval2

      ! Define function return
      logical :: bst_gt

      if (tval1 == tval2=>val) then
        return .true.
      else
        return .false.
      endif
    end function val_bst_eq
  end interface
contains

! Subroutine to initialize the lBlockCoords data structure inside tree_vals
subroutine initBlockCoords(node)
  implicit none

  ! Define passed parameters
  allocate(lBlockCoords)
  type(lBlockCoords), pointer :: node
  node%next => null()
  node%blockrow = -1
  node%blockcol = -1
end subroutine initBlockCoords

! Subroutine to deallocate our lblockCoords list
recursive subroutine destroyBlockCoordsList(root)
  implicit none

  ! Define passed parameters
  type (lBlockCoords), pointer :: root
  if (assosciated(root%next)) then
    call destroyBlockCoordsList(root%next)
  endif

  deallocate(root)
end subroutine destroyBlockCoordsList

! Subroutine to initialize a new tree node
subroutine tree_init(tree, val)
  implicit none

  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: tree
  integer, intent(in) :: val

  allocate(tree)
  tree%lval => null()
  tree%hval => null()
  tree%mval => null()

  tree%lchild =>  null()
  tree%mlchild => null()
  tree%mrchild => null()
  tree%rchild => null()
  tree%parent => null()
 
  call initBlockCoords(tree%vvblocks)
  call initBlockCoords(tree%cvblocks)
  call initBlockCoords(tree%ccblocks)
end subroutine tree_init

recursive subroutine tree_destroy(tree)
  implicit none

  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: tree

  ! Define local variables
  type(bst_atom_pair_node), pointer :: rchild, mlchild, lchild

  rchild => tree%rchild
  mlchild => tree%mlchild
  lchild => tree%lchild

  if ( associated(lchild) ) then
    call tree_destroy(lchild)
  endif

  if ( associated(mlchild) ) then
    call tree_destroy(mlchild)
  endif

  if ( associated(rchild) ) then
    call tree_destroy(rchild)
  endif

  call destroyBlockCoordsList(tree%vvblocks)
  call destroyBlockCoordsList(tree%cvblocks)
  call destroyBlockCoordsList(tree%ccblocks)
  deallocate(tree)
end subroutine tree_destroy

subroutine tree_destroyNode(node)
  implicit none

  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: node

  deallocate(node)
end subroutine tree_destroyNode

recursive subroutine tree_search(root, val, exists)
  implicit none

  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: root
  integer, intent(in) :: val
  logical, intent(out) :: exists

  ! Define local variables
  integer :: lmr

  exists = .false.
  
  if ( (root%lval == val) .or. (root%hval == val) ) then
    exists = .true.
    return
  endif

  call tree_moveLMR(root, val, lmr)
 
  if (lmr == 0) then
    call tree_search(root%lchild, val, exists)
    return
  elseif (lmr == 1) then
    call tree_search(root%mlchild, val, exists)
    return
  elseif (lmr == 2) then
    call tree_search(root%rchild, val, exists)
    return
  endif
end subroutine tree_search

subroutine tree_nodeType(node, ntype)
  implicit none

  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: node
  integer, intent(out) :: ntype

  if ( (node%hval > -1) .and. (node%mval > -1) ) then
    ntype = 4
  elseif ( (node%hval > -1) .and. (node%mval == -1) ) then
    ntype = 3
  else
    ntype = 2
  endif
end subroutine

subroutine tree_moveLMR(parent, val, lmr)
  implicit none

  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: parent
  integer, intent(in) :: val
  integer, intent(out) :: lmr

  ! Define local variables
  integer :: ntype

  call tree_nodeType(parent, ntype)

  call tree_LMR(parent, val, lmr)

  if ( (lmr == 0) .and. associated(parent%lchild) ) then
    return
  endif

  if (ntype == 2) then
    if ( (lmr==2) .and. associated(parent%rchild) ) then
      return
    endif
  else
    if (  (lmr == 1) .and. associated(parent%mlchild) ) then
      return
    endif

    if ( (lmr==2) .and. associated(parent%rchild) ) then
      return
    endif
  endif
  lmr = -1
end subroutine tree_moveLMR

subroutine tree_LMR(parent, val, lmr)
  implicit none


  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: parent
  integer, intent(in) :: val
  integer, intent(out) :: lmr

  ! Define local variables
  integer :: ntype

  call tree_nodeType(parent, ntype)
  lmr = -1

  if ( val < parent%lval ) then
    lmr = 0
    return
  endif

  if (ntype == 2) then
    if ( val > parent%lval ) then
      lmr = 2
      return
    endif
  else
    if ( (val > parent%lval) .and. (val < parent%hval) ) then
      lmr = 1
      return
    endif

    if ( (val > parent%hval) ) then
      lmr = 2
      return
    endif
  endif
end subroutine tree_LMR

recursive subroutine tree_insert(node, val)
  implicit none

  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: node
  integer, intent(in) :: val

  ! Define local variables
  integer :: lmr
  integer :: ntype, ptype
  logical :: exists
  logical :: isLeaf

  ! First we want to make sure the value is not already in the tree
  ! Took this out because we need to manually search before insert
  !call tree_search(node, val, exists)
  !if (exists) then
  !  return
  !endif

  ! If the value doesn't exist we want to find the place to put it
  call tree_moveLMR(node, val, lmr)
  if (lmr == 0) then
     call tree_insert(node%lchild, val)
  elseif (lmr == 1) then
     call tree_insert(node%mlchild, val)
  elseif (lmr == 2) then
     call tree_insert(node%rchild, val)
  endif

  ! We only want to insert the value if we are at a leaf node. This keeps
  ! recursive calls from accidentally adding the value elsewhere
  call tree_isLeaf(node, isLeaf)
  call tree_nodeType(node, ntype)
  if (isLeaf) then
    if ( ntype == 2 ) then
      call tree_twoToThree(node, val)
    else if ( ntype == 3 ) then
      call tree_threeToFour(node, val)
    endif
  endif


  ! Now that node is the node we should insert into we need to check what type
  ! of node it is and do the proper insertion. If the node is a 2-node, then
  ! we simply make it a 3-node. If node is a 3-Node, then we make it a 4 node,
  ! and recursively split the parents until the tree is balanced.
  call tree_nodeType(node, ntype)
  if (ntype == 4) then
    if ( .not. associated(node%parent) ) then
      call tree_splitRoot(node) 
      !call tree_printTree(node)
    endif
    if ( associated(node%parent) ) then
      call tree_moveMvalUp(node)
      if (.not. isLeaf) then
        call tree_split(node)
        !call tree_printTree(node%parent)
      endif
    endif
  endif
 
end subroutine tree_insert

subroutine tree_split(node)
  implicit none

  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: node

  ! Define local variables
  type(bst_atom_pair_node), pointer :: parent
  type(bst_atom_pair_node), pointer :: newNode
  integer :: ntype
  integer :: lmr
  logical :: isLeaf
  

  !call tree_moveMvalUp(node)
  !print *, "After move val up"
  !call flush(20)
  !    

  call tree_isLeaf(node, isLeaf)
  if ( associated(node%parent) .and. .not. isLeaf ) then
    !call tree_printTree(node%parent)
    !print *, "printTree"
    parent => node%parent
    call tree_nodeType(parent, ntype)
    !print *, "ntype"
    
    ! twoToThree and threeToFour both make the children correctly.
    ! What we need to do is adjust the subchildren. We first need to know
    ! which child node is, so we can adjust the subchildren correctly

!   You're on the right track, but these conditions aren't really present 
!   becuase the value is set to -1 in the children when the twoToThree and
!   threeToFour node routines are called. Consider not setting to -1 in the
!   routines, then these conditions will work. If not, then we have to make the
!   redistribution of children happen in this subroutine.
!   JUST THINK ABOUT IT. THIS ROUTINE SHOULD ALMOST COMPLETELY WORK BECAUSE
!   OF WHERE YOU ARE IN THE RECURSION TREE FROM tree_insert.  Just give it a 
!   little bit.

    if ( associated(parent%rchild%mrchild) ) then
      ! we came from the rightchild
      if (.not. associated(parent%mrchild)) then
        newNode => parent%mlchild
      else
        newNode => parent%mrchild
      endif
      newNode%lchild => parent%rchild%lchild 
      parent%rchild%lchild%parent => newNode

      newNode%rchild => parent%rchild%mlchild
      parent%rchild%mlchild%parent => newNode

      parent%rchild%lchild => parent%rchild%mrchild
      parent%rchild%rchild => parent%rchild%rchild
      parent%rchild%mrchild => null()
      parent%rchild%mlchild => null()

      !parent%rchild%lval = parent%rchild%hval
      !parent%rchild%hval = -1
    else if ( associated(parent%lchild%mrchild) ) then
      ! we came from the leftchild
      parent%mlchild%lchild => parent%lchild%mrchild
      parent%lchild%mrchild%parent => parent%mlchild

      parent%mlchild%rchild => parent%lchild%rchild
      parent%lchild%rchild%parent => parent%mlchild

      parent%lchild%rchild => parent%lchild%mlchild

      parent%lchild%mrchild => null()
      parent%lchild%mlchild => null()

    else 
      ! We came form the middle child
      ! the mrchild needs the mlchild's mrchild and rchild
      ! the mlchild then moves it's mlchild to rchild
      parent%mrchild%lchild => parent%mlchild%mrchild
      parent%mlchild%mrchild => parent%mrchild

      parent%mrchild%rchild => parent%mlchild%rchild
      parent%mlchild%rchild => parent%mrchild

      parent%mlchild%rchild => parent%mlchild%mlchild

      parent%mlchild%mlchild => null()
      parent%mlchild%mrchild => null()
    endif
  endif
end subroutine tree_split

subroutine tree_moveMvalUp(node)
  implicit none

  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: node

  ! Define local variables
  integer :: ntype
  integer :: tempVal

  ! Now we have to combine the mlchild and the mrchild. Also deallocate the
  ! mrnode, and dereference the mrchild pointer
  tempVal = node%mval
  node%mval = -1
  !if (associated(node%mrchild)) then
  !  call tree_twoToThree(node%mlchild, node%mrchild%lval)
  !  call tree_destroy(node%mrchild)
  !  node%mrchild => null()
  !endif

  if ( associated(node%parent) ) then
    call tree_nodeType(node%parent, ntype)
    if ( ntype == 2 ) then
      call tree_twoToThree(node%parent, tempVal)
    elseif (ntype == 3) then
      call tree_threeToFour(node%parent, tempVal)
    endif
  endif

end subroutine tree_moveMvalUp

subroutine tree_splitRoot(node)
  implicit none

  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: node

  ! Define local variables
  type(bst_atom_pair_node), pointer :: parent
  type(bst_atom_pair_node), pointer :: newNode
  integer :: lmr
  logical :: isleaf

  ! First we initialize the new node
  call tree_init(newNode, node%mval)

  ! Now we add the two child nodes
  call tree_addNode(newNode, node%lval, 0)
  call tree_addNode(newNode, node%hval, 3)

  ! Now we need the new nodes to point correctly to the sub-childrenA
  ! We also need to make sure to change the parent pointer for these
  !call tree_isLeaf(node, isLeaf)
  if (associated(node%lchild)) then
    newNode%lchild%lchild => node%lchild
    newNode%lchild%lchild%parent => newNode%lchild
  endif
  if (associated(node%mlchild)) then
    newNode%lchild%rchild => node%mlchild
    newNode%lchild%rchild%parent => newNode%lchild
  endif
  if (associated(node%mrchild)) then
    newNode%rchild%lchild => node%mrchild
    newNode%rchild%lchild%parent => newNode%rchild
  endif
  if (associated(node%rchild)) then
    newNode%rchild%rchild => node%rchild
    newNode%rchild%rchild%parent => newNode%rchild
  endif

  call tree_destroyNode(node)
  node => newNode

end subroutine tree_splitRoot

subroutine tree_threeToFour(node, val)
  implicit none

  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: node
  integer, intent(in) :: val

  ! Define local variables
  type(bst_atom_pair_node), pointer :: tempNode
  logical :: isLeaf
  integer :: tempVal

  call tree_isLeaf(node,isLeaf)

  if (val < node%lval) then
    node%mval = node%lval
    node%lval = val
    ! Now we need to split the child nodes if there are any
    ! If it meets the above condition then the value came from the left child
    ! we need to take the hval from the lchild and make a new node in the 
    ! mlchild position. We also need to move the mlchild to the mr position.
    if (.not. isLeaf) then
      node%mrchild => node%mlchild
      node%mlchild =>  null()
      tempVal = node%lchild%hval
      call tree_addNode(node, tempVal, 1)
      node%lchild%hval = -1
    endif
  elseif ( (node%lval < val) .and. (val < node%hval) ) then
    node%mval = val
    ! Now we need to split the child nodes if there are any
    ! If it meets the above condition then the value came from the middle child
    ! we need to take the hval from the mchild and create a new node in the 
    ! mrchild position. Additionally change the hval to -1 to correctly
    ! identify as a 2-node.
    if (.not. isLeaf) then
      tempVal = node%mlchild%hval
      call tree_addNode(node, tempVal, 2)
      node%mlchild%hval = -1
    endif
  else ! val > node%hval
    node%mval = node%hval
    node%hval = val
    ! Now we need to split the child nodes if there are any
    ! If it meets the above condition then the value came from the right child
    ! we need to take the lval from rchild and make a new node in the mrchild
    ! position. We also need to reassign the lval to hval, and hval to -1
    if (.not. isLeaf) then
      tempVal = node%rchild%lval
      call tree_addNode(node, tempVal, 2)
      node%rchild%lval = node%rchild%hval
      node%rchild%hval = -1
    endif
  endif
end subroutine

subroutine tree_twoToThree(node, val)
  implicit none

  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: node
  integer, intent(in) :: val

  ! Define local variables
  logical :: isLeaf
  integer :: tempVal

  call tree_isLeaf(node, isLeaf)

  if (val < node%lval) then
    node%hval = node%lval
    node%lval = val
    ! Now redistribute children if we aren't a leaf
    if (.not. isLeaf) then
      tempVal = node%lchild%hval
      node%lchild%hval = -1
      call tree_addNode(node, tempVal, 1)
    endif
  else
    node%hval = val
    ! Now redistribute children if we aren't a leaf
    if (.not. isLeaf) then
      tempVal = node%rchild%lval
      node%rchild%lval = node%rchild%hval
      node%rchild%hval = -1
      call tree_addNode(node, tempVal, 1)
    endif
  endif
end subroutine tree_twoToThree

subroutine tree_isLeaf(node, isLeaf)
  implicit none

  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: node
  logical, intent(out) :: isLeaf

  isLeaf = .true.
  if ( associated(node%lchild) .or. associated(node%mlchild) .or. &
     & associated(node%rchild) ) then
    isLeaf = .false.
  endif
end subroutine tree_isLeaf

subroutine tree_addNode(parent, val, lmr)
  implicit none
  
  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: parent
  integer, intent(in) :: val
  integer, intent(in) :: lmr

  ! Define local variables
  type(bst_atom_pair_node), pointer :: newChild

  call tree_init(newChild, val)

  newChild%parent => parent
  if ( (lmr == 0) .and. (.not. associated(parent%lchild)) ) then
    parent%lchild => newChild
  elseif ( (lmr == 1) .and. (.not. associated(parent%mlchild)) ) then
    parent%mlchild => newChild
  elseif ( (lmr == 2) .and. (.not. associated(parent%mrchild)) ) then
    parent%mrchild => newChild
  elseif ( (lmr == 3) .and. (.not. associated(parent%rchild)) ) then
    parent%rchild => newChild
  endif
end subroutine tree_addNode

subroutine tree_printTree(tree)
  implicit none

  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: tree

  ! Define local variables
  logical :: isLeaf

  print *, "root ", tree%lval, tree%mval, tree%hval
  call flush(20)

  call tree_printRecur(tree)
end subroutine tree_printTree

recursive subroutine tree_printRecur(tree)
  implicit none

  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: tree
  
  if (associated(tree%lchild)) then
    print *, "lchild", tree%lchild%lval,tree%lchild%mval, tree%lchild%hval
    call flush(20)
    call tree_printRecur(tree%lchild)
  endif
  if (associated(tree%mlchild)) then
    print *, "mlchild", tree%mlchild%lval,tree%mlchild%mval, tree%mlchild%hval
    call flush(20)
    call tree_printRecur(tree%mlchild)
  endif
  if (associated(tree%mrchild)) then
    print *, "mrchild", tree%mrchild%lval,tree%mrchild%mval, tree%mrchild%hval
    call flush(20)
    call tree_printRecur(tree%mrchild)
  endif
  if (associated(tree%rchild)) then
    print *, "rchild", tree%rchild%lval,tree%rchild%mval, tree%rchild%hval
    call flush(20)
    call tree_printRecur(tree%rchild)
  endif
  print *, "_____________"
  call flush(20)

end subroutine tree_printRecur

end module O_bstAtomPair
