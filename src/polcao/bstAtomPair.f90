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
    type(bst_atom_pair_node), pointer :: parent
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
    module procedure bst_gt
    module procedure bst_val_gt
    module procedure val_bst_gt
  end interface
  
  ! Define overload for less than when using tree_vals type
  interface operator(<)
    module procedure bst_lt
    module procedure bst_val_lt
    module procedure val_bst_lt
  end interface
  
  ! Define overload for equal to when using tree_vals type
  interface operator(==)
    module procedure bst_eq
    module procedure bst_val_eq
    module procedure val_bst_eq
  end interface

contains

! Subroutine to initialize the lBlockCoords data structure inside tree_vals
subroutine initBlockCoords(node)
  implicit none

  ! Define passed parameters
  type(lBlockCoords), pointer :: node

  allocate(node)
  node%next => null()
  node%blockrow = -1
  node%blockcol = -1
end subroutine initBlockCoords

! Subroutine to deallocate our lblockCoords list
recursive subroutine destroyBlockCoordsList(root)
  implicit none

  ! Define passed parameters
  type (lBlockCoords), pointer :: root

  if (associated(root%next)) then
    call destroyBlockCoordsList(root%next)
  endif

  deallocate(root)
end subroutine destroyBlockCoordsList

! Subroutine to dealloc allocated block vals for a tree_vals structure
subroutine destroyVals(val)
  implicit none

  ! Define passed parameters
  type(tree_vals), pointer :: val

  if ( associated(val%vvblocks) ) then
    call destroyBlockCoordsList(val%vvblocks)
  endif

  if ( associated(val%cvblocks) ) then
    call destroyBlockCoordsList(val%vvblocks)
  endif
  
  if ( associated(val%ccblocks) ) then
    call destroyBlockCoordsList(val%vvblocks)
  endif
end subroutine destroyVals

! Subroutine to initialize a new tree node
subroutine tree_init(tree, val)
  implicit none

  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: tree
  type(tree_vals), pointer :: val

  allocate(tree)
  tree%lval=>val
  tree%hval=>null()
  tree%mval=>null()

  tree%lchild=>null()
  tree%mlchild=>null()
  tree%mrchild=>null()
  tree%rchild=>null()
  tree%parent=>null()
 
  !call initBlockCoords(tree%vvblocks)
  !call initBlockCoords(tree%cvblocks)
  !call initBlockCoords(tree%ccblocks)
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

  call destroyVals(tree%lval)
  call destroyVals(tree%mval)
  call destroyVals(tree%hval)
  deallocate(tree)
end subroutine tree_destroy


subroutine tree_destroyNode(node)
  implicit none

  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: node

  deallocate(node)
end subroutine tree_destroyNode

recursive subroutine tree_search(root, val, exists, node)
  implicit none

  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: root
  type(tree_vals), pointer :: val
  logical, intent(out) :: exists
  type(tree_vals), pointer, intent(out) :: node

  ! Define local variables
  integer :: lmr

  exists = .false.
 
  if (associated(root%lval)) then
    if (root%lval == val) then
      node = root%lval
      exists = .true.
      return
    endif
  endif
  if (associated(root%hval)) then 
    if (root%hval == val) then
      node = root%hval
      exists = .true.
      return
    endif
  endif

  call tree_moveLMR(root, val, lmr)
 
  if (lmr == 0) then
    call tree_search(root%lchild, val, exists, node)
    return
  elseif (lmr == 1) then
    call tree_search(root%mlchild, val, exists, node)
    return
  elseif (lmr == 2) then
    call tree_search(root%rchild, val, exists, node)
    return
  endif
end subroutine tree_search

subroutine tree_nodeType(node, ntype)
  implicit none

  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: node
  integer, intent(out) :: ntype
  !print *, "fark"
  !call flush(20)
  !print *, "asdf, " , associated(node%lval)
  !call flush(20)
  !print *, "val: ", (node%lval>-1)

  if ( associated(node%hval) .and. associated(node%mval) ) then
  !if ( (node%hval > -1) .and. (node%mval > -1) ) then
    ntype = 4
  elseif ( associated(node%hval) .and. (.not. associated(node%mval)) ) then
  !elseif ( (node%hval > -1) .and. (node%mval == -1) ) then
    ntype = 3
  else
    ntype = 2
  endif
end subroutine

subroutine tree_moveLMR(parent, val, lmr)
  implicit none

  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: parent
  type(tree_vals), pointer :: val
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
  type(tree_vals), pointer :: val
  integer, intent(out) :: lmr

  ! Define local variables
  integer :: ntype

  call tree_nodeType(parent, ntype)
  lmr = -1

!  if (associated(parent%lval)) then
    if ( val < parent%lval ) then
      lmr = 0
      return
    endif
!  endif

  if (ntype == 2) then
!    if (associated(parent%lval)) then
      if ( val%val > parent%lval ) then
        lmr = 2
        return
      endif
!    endif
  else
!    if ( associated(parent%lval) .and. associated(parent%hval)) then
      if ( (val%val > parent%lval) .and. (val%val < parent%hval) ) then
        lmr = 1
        return
      endif
!    endif

!    if (associated(parent%hval)) then
      if ( (val%val > parent%hval) ) then
        lmr = 2
        return
      endif
!    endif
  endif
end subroutine tree_LMR

recursive subroutine tree_insert(node, val)
  implicit none

  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: node
  type(tree_vals), pointer :: val

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
  type(tree_vals), pointer :: tempVal

  ! Now we have to combine the mlchild and the mrchild. Also deallocate the
  ! mrnode, and dereference the mrchild pointer
  tempVal=>node%mval
  node%mval=>null()
  
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
  type(tree_vals), pointer :: val

  ! Define local variables
  type(bst_atom_pair_node), pointer :: tempNode
  type(tree_vals), pointer :: tempVal
  logical :: isLeaf

  call tree_isLeaf(node,isLeaf)

  if (val < node%lval) then
    node%mval=>node%lval
    node%lval=>val
    ! Now we need to split the child nodes if there are any
    ! If it meets the above condition then the value came from the left child
    ! we need to take the hval from the lchild and make a new node in the 
    ! mlchild position. We also need to move the mlchild to the mr position.
    if (.not. isLeaf) then
      node%mrchild=>node%mlchild
      node%mlchild=>null()
      tempVal=>node%lchild%hval
      call tree_addNode(node, tempVal, 1)
      node%lchild%hval=>null()
    endif
  elseif ( (node%lval < val) .and. (val < node%hval) ) then
    node%mval=>val
    ! Now we need to split the child nodes if there are any
    ! If it meets the above condition then the value came from the middle child
    ! we need to take the hval from the mchild and create a new node in the 
    ! mrchild position. Additionally change the hval to -1 to correctly
    ! identify as a 2-node.
    if (.not. isLeaf) then
      tempVal=>node%mlchild%hval
      call tree_addNode(node, tempVal, 2)
      node%mlchild%hval=>null()
    endif
  else ! val > node%hval
    node%mval=>node%hval
    node%hval=>val
    ! Now we need to split the child nodes if there are any
    ! If it meets the above condition then the value came from the right child
    ! we need to take the lval from rchild and make a new node in the mrchild
    ! position. We also need to reassign the lval to hval, and hval to -1
    if (.not. isLeaf) then
      tempVal=>node%rchild%lval
      call tree_addNode(node, tempVal, 2)
      node%rchild%lval=>node%rchild%hval
      node%rchild%hval=>null()
    endif
  endif
end subroutine tree_threeToFour

subroutine tree_twoToThree(node, val)
  implicit none

  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: node
  type(tree_vals), pointer :: val

  ! Define local variables
  logical :: isLeaf
  type(tree_vals), pointer :: tempVal

  call tree_isLeaf(node, isLeaf)

  if (val < node%lval) then
    node%hval=>node%lval
    node%lval=>val
    ! Now redistribute children if we aren't a leaf
    if (.not. isLeaf) then
      tempVal=> node%lchild%hval
      node%lchild%hval=>null()
      call tree_addNode(node, tempVal, 1)
    endif
  else
    node%hval=>val
    ! Now redistribute children if we aren't a leaf
    if (.not. isLeaf) then
      tempVal=>node%rchild%lval
      node%rchild%lval=>node%rchild%hval
      node%rchild%hval=>null()
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
  type(tree_vals), pointer :: val
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

  print *, "root"
  call printCurr(tree)
  call tree_printRecur(tree)
end subroutine tree_printTree

subroutine printCurr(tree)
  implicit none

  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: tree
  
  if (.not. associated(tree%hval)) then
    print *, tree%lval%val
    call flush(20)
  else if (.not. associated(tree%mval)) then
    print *, tree%lval%val, tree%hval%val  
    call flush(20)
  else
    print *, tree%lval%val, tree%mval%val, tree%hval%val  
    call flush(20)
  end if
end subroutine printCurr

recursive subroutine tree_printRecur(tree)
  implicit none

  ! Define passed parameters
  type(bst_atom_pair_node), pointer :: tree
  
  if (associated(tree%lchild)) then
    print *, "lchild"
    call printCurr(tree%lchild)
    call flush(20)
    call tree_printRecur(tree%lchild)
  endif
  if (associated(tree%mlchild)) then
    print *, "mlchild"
    call printCurr(tree%mlchild)
    call flush(20)
    call tree_printRecur(tree%mlchild)
  endif
  if (associated(tree%mrchild)) then
    print *, "mrchild"
    call printCurr(tree%mrchild)
    call flush(20)
    call tree_printRecur(tree%mrchild)
  endif
  if (associated(tree%rchild)) then
    print *, "rchild"
    call printCurr(tree%rchild)
    call flush(20)
    call tree_printRecur(tree%rchild)
  endif
  print *, "_____________"
  call flush(20)

end subroutine tree_printRecur

! Functions for overloading tree_vals operations with greater than
function bst_gt(tval1, tval2)                    
  implicit none                              
                                                 
  ! Define passed parameters                     
  type(tree_vals), pointer, intent(in) :: tval1, tval2   
                                                 
  ! Define function return                       
  logical :: bst_gt                          
                
  if (tval1%val > tval2%val) then              
    bst_gt = .true.                              
    return                                       
  else                                           
    bst_gt = .false.                             
    return                                       
  endif                                          
end function bst_gt                              
                                                 
function bst_val_gt(tval1, tval2)                
  implicit none                              
                                                 
  ! Define passed parameters                     
  type(tree_vals), pointer, intent(in) :: tval1              
  integer, intent(in) :: tval2               
                                                 
  ! Define function return                       
  logical :: bst_val_gt                          
                                                 
  if (.not. associated(tval1)) then
    bst_val_gt = .false.
    return
  endif

  if (tval1%val > tval2) then                   
    bst_val_gt = .true.                          
    return                                       
  else                                           
    bst_val_gt = .false.                         
    return                                       
  endif                                          
end function bst_val_gt                          
                                                 
function val_bst_gt(tval1, tval2)                
  implicit none                              
                                                 
  ! Define passed parameters                     
  integer, intent(in) :: tval1                   
  type(tree_vals), pointer, intent(in) :: tval2          
                                                 
  ! Define function return                       
  logical :: val_bst_gt                          
  
  if (.not. associated(tval2)) then
    val_bst_gt = .false.
    return
  endif
                                                 
  if (tval1 > tval2%val) then                   
    val_bst_gt = .true.                          
    return                                       
  else                                           
    val_bst_gt = .false.                         
    return                                       
  endif                                          
end function val_bst_gt

! Functions for overloading tree_vals operations with less than
function bst_lt(tval1, tval2)
  implicit none

  ! Define passed parameters
  type(tree_vals), pointer, intent(in) :: tval1, tval2

  ! Define function return
  logical :: bst_lt

  if (tval1%val < tval2%val) then
    bst_lt = .true.
    return
  else
    bst_lt = .false.
    return
  endif
end function bst_lt

function bst_val_lt(tval1, tval2)
  implicit none

  ! Define passed parameters
  type(tree_vals), pointer, intent(in) :: tval1
  integer, intent(in) :: tval2

  ! Define function return
  logical :: bst_val_lt

  if (.not. associated(tval1)) then
    bst_val_lt = .false.
    return
  endif

  if (tval1%val < tval2) then
    bst_val_lt = .true.
    return
  else
    bst_val_lt = .false.
    return
  endif
end function bst_val_lt

function val_bst_lt(tval1, tval2)
  implicit none

  ! Define passed parameters
  integer, intent(in) :: tval1
  type(tree_vals), pointer, intent(in) :: tval2

  ! Define function return
  logical :: val_bst_lt

  if (.not. associated(tval2)) then
    val_bst_lt = .false.
    return
  endif

  if (tval1 < tval2%val) then
    val_bst_lt = .true.
    return 
  else
    val_bst_lt = .false.
    return
  endif
end function val_bst_lt


! Functions for overloading tree_vals operations with less than
function bst_eq(tval1, tval2)
  implicit none

  ! Define passed parameters
  type(tree_vals), pointer, intent(in) :: tval1, tval2

  ! Define function return
  logical :: bst_eq

  if (tval1%val == tval2%val) then
    bst_eq = .true.
    return
  else
    bst_eq = .false.
    return
  endif
end function bst_eq

function bst_val_eq(tval1, tval2)
  implicit none

  ! Define passed parameters
  type(tree_vals), pointer, intent(in) :: tval1
  integer, intent(in) :: tval2

  ! Define function return
  logical :: bst_val_eq

  if (.not. associated(tval1)) then
    bst_val_eq = .false.
    return
  endif

  if (tval1%val == tval2) then
    bst_val_eq = .true.
    return
  else
    bst_val_eq = .false.
    return
  endif
end function bst_val_eq

function val_bst_eq(tval1, tval2)
  implicit none

  ! Define passed parameters
  integer, intent(in) :: tval1
  type(tree_vals), pointer, intent(in) :: tval2

  ! Define function return
  logical :: val_bst_eq

  if (.not. associated(tval2)) then
    val_bst_eq = .false.
    return
  endif

  if (tval1 == tval2%val) then
    val_bst_eq = .true.
    return
  else
    val_bst_eq = .false.
    return
  endif
end function val_bst_eq

! Subroutines for overloading assignment operators when using tree_vals
!subroutine bst_assn_bst(tval1, tval2)
!  implicit none
!
!  ! Define passed parameters
!  type(tree_vals), pointer, intent(out) :: tval1
!  type(tree_Vals), pointer, intent(in) ::  tval2
!
!  tval1%val = tval2%val
!  tval1%vvblocks = tval2%vvblocks
!  tval1%vvblocks = tval2%cvblocks
!  tval1%vvblocks = tval2%ccblocks
!end subroutine bst_assn_bst
!
!subroutine bst_assn_int(tval1, tval2)
!  implicit none
!
!  ! Define passed parameters
!  type(tree_vals), pointer, intent(out) :: tval1
!  integer :: tval2
!
!  tval1%val = tval2
!end subroutine bst_assn_int

end module O_bstAtomPair
