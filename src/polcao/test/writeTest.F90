program writeTest
  use HDF5
  use MPI
  
  implicit none

  integer :: mpirank
  integer :: mpisize
  integer :: mpierr
  integer, dimension(MPI_STATUS_SIZE) :: mpistatus

  integer, dimension(6,6) :: sdata
  integer :: dcount

  integer :: i,j,n

  integer :: to, from, tag
 
  ! HDF5 variables
  integer(HID_T) :: file_id
  integer(HID_T) :: dset_id
  integer(HID_T) :: dspace_id

  complex (kind=kind(0.0d0)) :: ctype
  real (kind=kind(0.0d0)) :: rtype

  call MPI_Init(mpierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpirank,mpierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpisize,mpierr)

  if (mpirank == 0) then
    call setupHDF5(mpisize, file_id, dset_id, dspace_id)
  endif

  print *, "kind: ", kind(0.0d0), kind(ctype), kind(rtype), MPI_COMPLEX
  !from = 1
  to = 0
  tag = 0
  sdata(:,:) = 0
  dcount = 36
  call initData(mpirank, sdata)

  if (mpirank == 0) then
    call writeValeVale(sdata, file_id, dset_id, dspace_id, 0)
  endif
  call MPI_Barrier(MPI_COMM_WORLD, mpierr)
  do i=1,mpisize-1
    if (mpirank == i) then
      print *, "Process ",i,"Sending."
      call MPI_Send(sdata(:,:), dcount, MPI_INTEGER, to, tag, &
                & MPI_COMM_WORLD, mpierr)
    endif
    if (mpirank == 0) then
      print *, "Process 0, receiving from ",i
      from = i
      call MPI_Recv(sdata(:,:), dcount, MPI_INTEGER, from, tag, &
        & MPI_COMM_WORLD, mpistatus, mpierr)
      call writeValeVale(sdata, file_id, dset_id, dspace_id, from)
    endif
    call MPI_Barrier(MPI_COMM_WORLD, mpierr)
  end do

  !if (mpirank == 1) then
  !  print *, "Send Rank: ", mpirank
  !  call MPI_Send(sdata(:,:), dcount, MPI_INTEGER, to, tag, &
  !            & MPI_COMM_WORLD, mpierr)

  !else if (mpirank == 0) then
  !  print *, "Recv Rank: ", mpirank
  !  call MPI_Recv(sdata(:,:), dcount, MPI_INTEGER, from, tag, MPI_COMM_WORLD, &
  !    & mpistatus, mpierr)

  !  do i=1,5
  !    do j=1,5
  !      print *, mpirank, sdata(i,j)
  !    enddo
  !  enddo

  !endif

  if (mpirank == 0) then
    call closeHDF5(file_id, dset_id, dspace_id)
  endif

  call MPI_FINALIZE(mpierr)

end program writeTest

subroutine initData(mpirank, sdata)
  implicit none

  ! Define passed parameters
  integer, dimension(6,6), intent(inout) :: sdata
  integer, intent(in) :: mpirank

  ! Define local variables
  integer :: i,j,n

  !n = (mpirank)*36

  !do i=1,6
  !  do j=1,6
  !    sdata(j,i) = n
  !    n = n + 1
  !  enddo
  !enddo
  if (mpirank == 0) then
    sdata(1,:) = (/  0,   1,   2,   6,   7,   8/)
    sdata(2,:) = (/ 12,  13,  14,  18,  19,  20/)
    sdata(3,:) = (/ 24,  25,  26,  30,  31,  32/)
    sdata(4,:) = (/ 72,  73,  74,  78,  79,  80/)
    sdata(5,:) = (/ 84,  85,  86,  90,  91,  92/)
    sdata(6,:) = (/ 96,  97,  98, 102, 103, 104/)
  else if (mpirank == 1) then
    sdata(1,:) = (/  3,   4,   5,   9,  10,  11/)
    sdata(2,:) = (/ 15,  16,  17,  21,  22,  23/)
    sdata(3,:) = (/ 27,  28,  29,  33,  34,  35/)
    sdata(4,:) = (/ 75,  76,  77,  81,  82,  83/)
    sdata(5,:) = (/ 87,  88,  89,  93,  94,  95/)
    sdata(6,:) = (/ 99, 100, 101, 105, 106, 107/)
  else if (mpirank == 2) then
    sdata(1,:) = (/ 36,  37,  38,  42,  43,  44/)
    sdata(2,:) = (/ 48,  49,  50,  54,  55,  56/)
    sdata(3,:) = (/ 60,  61,  62,  66,  67,  68/)
    sdata(4,:) = (/108, 109, 110, 114, 115, 116/)
    sdata(5,:) = (/120, 121, 122, 126, 127, 128/)
    sdata(6,:) = (/132, 133, 134, 138, 139, 140/)
  else if (mpirank == 3) then
    sdata(1,:) = (/ 39,  40,  41,  45,  46,  47/)
    sdata(2,:) = (/ 51,  52,  53,  57,  58,  59/)
    sdata(3,:) = (/ 63,  64,  65,  69,  70,  71/)
    sdata(4,:) = (/111, 112, 113, 117, 118, 119/)
    sdata(5,:) = (/123, 124, 125, 129, 130, 131/)
    sdata(6,:) = (/135, 136, 137, 141, 142, 143/)
  endif

end subroutine initData

subroutine setupHDF5(mpisize, file_id, dset_id, dspace_id)
  use HDF5

  implicit none

  ! Define passed parameters
  integer, intent(in) :: mpisize
  integer(HID_T) :: file_id
  integer(HID_T) :: dset_id
  integer(HID_T) :: dspace_id

  ! Define local variables
  integer :: hdferr
  integer :: rank
  integer(HSIZE_T), dimension(2) :: dims

  rank = 2
  dims(:) = int(sqrt(dble(mpisize)))*6

  ! initialize hdf5 interface
  call h5open_f(hdferr)

  ! Create a new file using default properties
  call h5fcreate_f("test.hdf5", H5F_ACC_TRUNC_F, file_id, hdferr)

  ! Create the dataspace
  call h5screate_simple_f(rank, dims, dspace_id, hdferr)

  ! Create the dataset with default properties
  call h5dcreate_f(file_id, "sdata", H5T_NATIVE_INTEGER, dspace_id, &
    & dset_id, hdferr)


end subroutine setupHDF5

subroutine closeHDF5(file_id, dset_id, dspace_id)
  use HDF5

  implicit none
 
  ! Define passed parameters
  integer(HID_T) :: file_id
  integer(HID_T) :: dset_id
  integer(HID_T) :: dspace_id

  ! Define local variables
  integer :: hdferr

  ! end access to the dataset and release resources used by it.
  call h5dclose_f(dset_id, hdferr)
  
  ! Terminate access to the data space
  call h5sclose_f(dspace_id, hdferr)

  ! Close the file
  call h5fclose_f(file_id, hdferr)

  ! Close the fortran interface
  call h5close_f(hdferr)

end subroutine closeHDF5

subroutine writeValeVale(sdata, file_id, dset_id, dspace_id, mpirank)
  use HDF5

  implicit none
  
  ! Define passed parameters
  integer, dimension(6,6), intent(in) :: sdata
  integer(HID_T) :: file_id
  integer(HID_T) :: dset_id
  integer(HID_T) :: dspace_id
  integer :: mpirank

  ! Define local variables
  integer(hsize_t), dimension(2) :: hslabCount, hslabStart, dims
  integer(hid_t) :: memspace_dsid
  integer :: i, j, a,b, pr, pc, hdferr
  integer, dimension(2) :: lo,hi


  hslabCount(:) = 3

  if (mpirank == 0) then
    pr = 0
    pc = 0
  else if (mpirank == 1) then
    pr = 0
    pc = 1
  else if (mpirank == 2) then
    pr = 1
    pc = 0
  else if (mpirank == 3) then
    pr = 1
    pc = 1
  endif

  a = 0
  b = 0
  do i=0,1
    do j=0,1
      a = i*3
      b = j*3
      call localToGlobalMap(a,b, lo, hi, 2, 2, 3, 3, pr, pc, 0)
      hslabStart(1) = lo(1)
      hslabStart(2) = lo(2)
      print *, lo(1),lo(2), hslabStart(1),hslabStart(2)
      
      call h5screate_simple_f(2, hslabCount, memspace_dsid, hdferr)
    
      ! Define hyperslab to be written to
      call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hslabStart, &
        & hslabCount, hdferr)

      ! write slab to disk
      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, &
        & sdata(a+1:a+3,b+1:b+3), hslabCount, &
        hdferr, file_space_id=dspace_id, mem_space_id=memspace_dsid)
    enddo
  enddo

  !call localToGlobalMap(0,0, lo, hslabCount, 2, 2, 3, 3, pr, pc, 0)
  !print *, "1,1: ", lo!, hslabCount
  !call localToGlobalMap(0,3, lo, hslabCount, 2, 2, 3, 3, pr, pc, 0)
  !print *, "1,4: ", lo!, hslabCount
  !call localToGlobalMap(3,0, lo, hslabCount, 2, 2, 3, 3, pr, pc, 0)
  !print *, "4,1: ", lo!, hslabCount
  !call localToGlobalMap(3,3, lo, hslabCount, 2, 2, 3, 3, pr, pc, 0)
  !print *, "4,4: ", lo!, hslabCount

end subroutine writeValeVale

subroutine localToGlobalMap(a, b, glo, ghi, prows, pcols, mb, nb, &
    & myprow, mypcol, extra)
  implicit none

  ! Define passed parameters
  integer, intent(in) :: prows, pcols 
  integer, intent(in) :: mb, nb
  integer, intent(in) :: myprow, mypcol 
  integer, intent(in) :: a,b
  integer, dimension(2), intent(out) :: glo, ghi ! target global values
  integer, intent(in) :: extra ! Control variable that denotes if we are doing
                               ! a irregular size block. Options
                               ! = 0 normal block size
                               !   1 nrows is less than mb
                               !   2 ncols is less than nb
                               !   3 both nrows and ncols is less than nb and mb

  ! Calculate the upper left corner
  !glo(1) = (a-0)*prows + 1
  !glo(2) = (b-0)*pcols + 1
  glo(1) = (((a-0)/mb)*prows+myprow)*mb + 0
  glo(2) = (((b-0)/nb)*pcols+mypcol)*nb + 0

  ghi(1) = glo(1) + mb
  ghi(2) = glo(2) + nb
  ! Calculate the bottom right corner
!  if (extra == 0) then
!    ghi(1) = (a-mb)*prows + 1
!    ghi(2) = (b-nb)*pcols + 1
!  elseif (extra == 1) then
!    ghi(1) = glo(1) + extraBlockRows
!    ghi(2) = (b-nb)*pcols + 1
!  elseif (extra == 2) then
!    ghi(1) = (a-mb)*prows + 1
!    ghi(2) = glo(2) + extraBlockCols
!  elseif (extra == 3) then
!    ghi(1) = glo(1) + extraBlockRows
!    ghi(2) = glo(2) + extraBlockCols
!  endif

end subroutine localToGlobalMap 
