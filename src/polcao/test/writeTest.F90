program writeTest
  use HDF5
  use MPI
  use O_Parallel
  
  implicit none

  integer :: mpirank
  integer :: mpisize
  integer :: mpierr
  integer, dimension(MPI_STATUS_SIZE) :: mpistatus

  integer :: dcount

  integer :: i,j,n

  integer :: to, from, tag
 
  ! HDF5 variables
  integer(HID_T) :: file_id
  integer(HID_T) :: dset_id
  integer(HID_T) :: dspace_id

  complex (kind=kind(0.0d0)) :: ctype
  real (kind=kind(0.0d0)) :: rtype

  type(BlacsInfo) :: blcsinfo
  type(ArrayInfo) :: arrinfo

  integer, dimension(2) :: sdata

  call MPI_Init(mpierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpirank,mpierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpisize,mpierr)

  if (mpirank == 0) then
    call setupHDF5(mpisize, file_id, dset_id, dspace_id)
  endif

  call MPI_Barrier(MPI_COMM_WORLD, mpierr)

  call setupBlacs(blcsinfo)
  
  call setupArrayDesc(arrinfo, blcsinfo, 13, 13, 1)

  to = 0
  tag = 0
  dcount = size(arrinfo%local,1)*size(arrinfo%local,2)*size(arrinfo%local,3)
  call initData(mpirank, arrinfo)
  call initData2(mpirank, arrinfo)

  if (mpirank == 0) then
    call writeValeVale(arrinfo, blcsinfo, file_id, dset_id, dspace_id, 0)
  endif
  call MPI_Barrier(MPI_COMM_WORLD, mpierr)

  do i=1,mpisize-1
    if (mpirank == i) then
      sdata = (/blcsinfo%myprow, blcsinfo%mypcol/)
      call MPI_Send(sdata, 2, MPI_INTEGER, to, tag, MPI_COMM_WORLD, mpierr)
      call MPI_Send(arrinfo%local(:,:,:), dcount, MPI_DOUBLE_COMPLEX, to, tag, &
                & MPI_COMM_WORLD, mpierr)
    endif
    if (mpirank == 0) then
      from = i
      call MPI_Recv(sdata, 2, MPI_INTEGER, from, tag, MPI_COMM_WORLD, &
        & mpistatus, mpierr)

      blcsinfo%myprow = sdata(1)
      blcsinfo%mypcol = sdata(2)
      call deallocLocalArray(arrinfo)
      call allocLocalArray(arrinfo,blcsinfo)
      call getBlockDims(arrinfo, blcsinfo)

      dcount = size(arrinfo%local,1)*size(arrinfo%local,2)

      call MPI_Recv(arrinfo%local(:,:,:),dcount,MPI_DOUBLE_COMPLEX, from, tag, &
        & MPI_COMM_WORLD, mpistatus, mpierr)
      call writeValeVale(arrinfo, blcsinfo, file_id, dset_id, dspace_id, from)
    endif
    call MPI_Barrier(MPI_COMM_WORLD, mpierr)
  end do

  call deallocLocalArray(arrinfo)

  call MPI_Barrier(MPI_COMM_WORLD, mpierr)
  if (mpirank == 0) then
    call closeHDF5(file_id, dset_id, dspace_id)
  endif

  call MPI_Barrier(MPI_COMM_WORLD, mpierr)
  call MPI_FINALIZE(mpierr)

end program writeTest

subroutine initData(mpirank, arrinfo)
  use O_Parallel

  implicit none

  ! Define passed parameters
  !integer, dimension(6,6), intent(inout) :: sdata
  type(ArrayInfo) :: arrinfo
  integer, intent(in) :: mpirank

  if (mpirank == 0) then
    arrinfo%local(1,:,1) = (/  0,   1,   2,   6,   7,   8,  -1/)
    arrinfo%local(2,:,1) = (/ 12,  13,  14,  18,  19,  20,  -1/)
    arrinfo%local(3,:,1) = (/ 24,  25,  26,  30,  31,  32,  -1/)
    arrinfo%local(4,:,1) = (/ 72,  73,  74,  78,  79,  80,  -1/)
    arrinfo%local(5,:,1) = (/ 84,  85,  86,  90,  91,  92,  -1/)
    arrinfo%local(6,:,1) = (/ 96,  97,  98, 102, 103, 104,  -1/)
    arrinfo%local(7,:,1) = (/ -1,  -1,  -1,  -1,  -1,  -1,  -1/)
  else if (mpirank == 1) then                             
    arrinfo%local(1,:,1) = (/  3,   4,   5,   9,  10,  11/)
    arrinfo%local(2,:,1) = (/ 15,  16,  17,  21,  22,  23/)
    arrinfo%local(3,:,1) = (/ 27,  28,  29,  33,  34,  35/)
    arrinfo%local(4,:,1) = (/ 75,  76,  77,  81,  82,  83/)
    arrinfo%local(5,:,1) = (/ 87,  88,  89,  93,  94,  95/)
    arrinfo%local(6,:,1) = (/ 99, 100, 101, 105, 106, 107/)
    arrinfo%local(7,:,1) = (/  1,   1,   1,   1,   1,   1/)
  else if (mpirank == 2) then
    arrinfo%local(1,:,1) = (/ 36,  37,  38,  42,  43,  44, 2/)
    arrinfo%local(2,:,1) = (/ 48,  49,  50,  54,  55,  56, 2/)
    arrinfo%local(3,:,1) = (/ 60,  61,  62,  66,  67,  68, 2/)
    arrinfo%local(4,:,1) = (/108, 109, 110, 114, 115, 116, 2/)
    arrinfo%local(5,:,1) = (/120, 121, 122, 126, 127, 128, 2/)
    arrinfo%local(6,:,1) = (/132, 133, 134, 138, 139, 140, 2/)
  else if (mpirank == 3) then
    arrinfo%local(1,:,1) = (/ 39,  40,  41,  45,  46,  47/)
    arrinfo%local(2,:,1) = (/ 51,  52,  53,  57,  58,  59/)
    arrinfo%local(3,:,1) = (/ 63,  64,  65,  69,  70,  71/)
    arrinfo%local(4,:,1) = (/111, 112, 113, 117, 118, 119/)
    arrinfo%local(5,:,1) = (/123, 124, 125, 129, 130, 131/)
    arrinfo%local(6,:,1) = (/135, 136, 137, 141, 142, 143/)
  endif

end subroutine initData

subroutine initData2(mpirank, arrinfo)
  use O_Parallel

  implicit none

  ! Define passed parameters
  !integer, dimension(6,6), intent(inout) :: sdata
  type(ArrayInfo) :: arrinfo
  integer, intent(in) :: mpirank

  if (mpirank == 0) then
    arrinfo%local(1,:,1)=(/(1,-1),(1,-1),(1,-1),(1,-1),(1,-1),(1,-1),(1,-1)/)
    arrinfo%local(2,:,1)=(/(1,-1),(1,-1),(1,-1),(1,-1),(1,-1),(1,-1),(1,-1)/)
    arrinfo%local(3,:,1)=(/(1,-1),(1,-1),(1,-1),(1,-1),(1,-1),(1,-1),(1,-1)/)
    arrinfo%local(4,:,1)=(/(1,-1),(1,-1),(1,-1),(1,-1),(1,-1),(1,-1),(1,-1)/)
    arrinfo%local(5,:,1)=(/(1,-1),(1,-1),(1,-1),(1,-1),(1,-1),(1,-1),(1,-1)/)
    arrinfo%local(6,:,1)=(/(1,-1),(1,-1),(1,-1),(1,-1),(1,-1),(1,-1),(1,-1)/)
    arrinfo%local(7,:,1)=(/(1,-1),(1,-1),(1,-1),(1,-1),(1,-1),(1,-1),(1,-1)/)
  else if (mpirank == 1) then                             
    arrinfo%local(1,:,1)=(/(2,-2),(2,-2),(2,-2),(2,-2),(2,-2),(2,-2)/)
    arrinfo%local(2,:,1)=(/(2,-2),(2,-2),(2,-2),(2,-2),(2,-2),(2,-2)/)
    arrinfo%local(3,:,1)=(/(2,-2),(2,-2),(2,-2),(2,-2),(2,-2),(2,-2)/)
    arrinfo%local(4,:,1)=(/(2,-2),(2,-2),(2,-2),(2,-2),(2,-2),(2,-2)/)
    arrinfo%local(5,:,1)=(/(2,-2),(2,-2),(2,-2),(2,-2),(2,-2),(2,-2)/)
    arrinfo%local(6,:,1)=(/(2,-2),(2,-2),(2,-2),(2,-2),(2,-2),(2,-2)/)
    arrinfo%local(7,:,1)=(/(2,-2),(2,-2),(2,-2),(2,-2),(2,-2),(2,-2)/)
  else if (mpirank == 2) then
    arrinfo%local(1,:,1)=(/(3,-3),(3,-3),(3,-3),(3,-3),(3,-3),(3,-3),(3,-3)/)
    arrinfo%local(2,:,1)=(/(3,-3),(3,-3),(3,-3),(3,-3),(3,-3),(3,-3),(3,-3)/)
    arrinfo%local(3,:,1)=(/(3,-3),(3,-3),(3,-3),(3,-3),(3,-3),(3,-3),(3,-3)/)
    arrinfo%local(4,:,1)=(/(3,-3),(3,-3),(3,-3),(3,-3),(3,-3),(3,-3),(3,-3)/)
    arrinfo%local(5,:,1)=(/(3,-3),(3,-3),(3,-3),(3,-3),(3,-3),(3,-3),(3,-3)/)
    arrinfo%local(6,:,1)=(/(3,-3),(3,-3),(3,-3),(3,-3),(3,-3),(3,-3),(3,-3)/)
  else if (mpirank == 3) then
    arrinfo%local(1,:,1)=(/(4,-4),(4,-4),(4,-4),(4,-4),(4,-4),(4,-4)/)
    arrinfo%local(2,:,1)=(/(4,-4),(4,-4),(4,-4),(4,-4),(4,-4),(4,-4)/)
    arrinfo%local(3,:,1)=(/(4,-4),(4,-4),(4,-4),(4,-4),(4,-4),(4,-4)/)
    arrinfo%local(4,:,1)=(/(4,-4),(4,-4),(4,-4),(4,-4),(4,-4),(4,-4)/)
    arrinfo%local(5,:,1)=(/(4,-4),(4,-4),(4,-4),(4,-4),(4,-4),(4,-4)/)
    arrinfo%local(6,:,1)=(/(4,-4),(4,-4),(4,-4),(4,-4),(4,-4),(4,-4)/)
  endif

end subroutine initData2

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
  !dims(:) = int(sqrt(dble(mpisize)))*7
  dims(:) = 13 

  ! initialize hdf5 interface
  call h5open_f(hdferr)

  ! Create a new file using default properties
  call h5fcreate_f("test.hdf5", H5F_ACC_TRUNC_F, file_id, hdferr)

  ! Create the dataspace
  call h5screate_simple_f(rank, dims, dspace_id, hdferr)

  ! Create the dataset with default properties
  call h5dcreate_f(file_id, "sdata", H5T_NATIVE_DOUBLE, dspace_id, &
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

subroutine writeValeVale(arrinfo, blcsinfo, & 
    & file_id, dset_id, dspace_id, mpirank)
  use HDF5

  use O_Kinds
  use O_Parallel

  implicit none
  
  ! Define passed parameters
  type(ArrayInfo) :: arrinfo
  type(BlacsInfo) :: blcsinfo
  integer(HID_T) :: file_id
  integer(HID_T) :: dset_id
  integer(HID_T) :: dspace_id
  integer :: mpirank

  ! Define local variables
  integer(hsize_t), dimension(3) :: hslabCount, hslabStart
  integer(hid_t) :: memspace_dsid
  integer :: i, j, k, l, x, y, a,b, hdferr
  integer, dimension(2) :: lo,hi

  real(kind=double), allocatable, dimension(:,:,:) :: dataOut

  ! Set the third component to numkp rather than setting in loop every time
  hslabCount(3) = size(arrinfo%local,3)

  ! Allocate space for single blocks
  allocate(dataOut(arrinfo%mb, arrinfo%nb, hslabCount(3)))

  do i=0,arrinfo%nrblocks-1
    do j=0,arrinfo%ncblocks-1
      if ((arrinfo%extraRows>0) .and. (i==arrinfo%nrblocks-1) .and. &
        & (arrinfo%extraCols>0) .and. (j==arrinfo%ncblocks-1)) then
        hslabCount(1) = arrinfo%extraRows
        hslabCount(2) = arrinfo%extraCols
      else if ((arrinfo%extraRows>0) .and. (i==arrinfo%nrblocks-1)) then
        hslabCount(1) = arrinfo%extraRows
        hslabCount(2) = arrinfo%nb
      else if ((arrinfo%extraCols>0) .and. (j==arrinfo%ncblocks-1)) then
        hslabCount(1) = arrinfo%mb
        hslabCount(2) = arrinfo%extraCols
      else
        hslabCount(1) = arrinfo%mb
        hslabCount(2) = arrinfo%nb
      endif

      a = i*arrinfo%mb
      b = j*arrinfo%nb
      call localToGlobalMap(a,b, lo, hi, arrinfo, blcsinfo, 0)
      hslabStart(1) = lo(1)
      hslabStart(2) = lo(2)
      hslabStart(3) = 1

      call h5screate_simple_f(3, hslabCount, memspace_dsid, hdferr)
    
      ! Define hyperslab to be written to
      call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hslabStart, &
        & hslabCount, hdferr)

      ! Need to prepare the data so that complex parts are on saved on the
      ! bottom half of matrix, and real parts on the top half.
      x=1
      do k=hslabStart(1),(hslabStart(1)+hslabCount(1)-1)
        y=1
        do l=hslabStart(2),(hslabStart(2)+hslabCount(2)-1)
          if (l>=k) then ! top half
            dataOut(x,y,:) = real(arrinfo%local(a+x,b+y,:))
          else ! Bottom half
            dataOut(x,y,:) = aimag(arrinfo%local(a+x,b+y,:))
          endif
          y = y + 1
        enddo
        x = x + 1
      enddo

      ! write slab to disk
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
        & dataOut(1:hslabCount(1),1:hslabCount(2), :), hslabCount, hdferr, &
        & file_space_id=dspace_id, mem_space_id=memspace_dsid)
    enddo
  enddo

  deallocate(dataOut)

end subroutine writeValeVale
