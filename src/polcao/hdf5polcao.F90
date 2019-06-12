module O_pOLCAOhdf5
  ! Iport any necessary definition modules.
  
  ! Make sure that no funny variables are defined.
  implicit none
  
  contains

subroutine accessAtomPairsHDF5(atomList, pgrid)
  ! Import necessary modules.
  use hdf5
  use O_ParallelSetup

  ! Make sure that no funny variables are defined
  implicit none

  ! Define passed parameters
  type(AtomPair), pointer :: atomList
  integer, dimension(2), intent(inout) :: pgrid

  ! Declare local variables
  integer :: hdferr
  integer(hid_t) :: atomPair_fid
  integer(hid_t) :: atomPair_plid

  ! Create the property list for the hdf5 file
  !call h5pcreate_f(H5P_FILE_ACCESS_F,atomPair_plid,hdferr)
  !if (hdferr /= 0) stop 'Failed to create setup plid in accessAtomPairsHDF5.'

  ! Open the HDF5 file that contains the processor distribution information
  ! for PBLAS/SCALAPACK.
  call h5fopen_f("polcaoDistribution.hdf5",H5F_ACC_RDONLY_F,atomPair_fid,hdferr)
  if (hdferr /= 0) stop 'Could not open polcaoDistribution.hdf5'

  call readGrid(atomPair_fid, pgrid)
  call readAtomListHDF5(atomPair_fid, atomList)

  ! We're finished with the file we can close it now
  call h5fclose_f(atomPair_fid, hdferr)
  if (hdferr /=0) stop 'Could not close atomPair_fid'

end subroutine accessAtomPairsHDF5

subroutine readGrid(atomPair_fid, pgrid)

  ! Import any necessary modules.
  use hdf5
  use MPI

  ! Define passed parameters
  integer(hid_t), intent(in) :: atomPair_fid
  integer, dimension(2), intent(out) :: pgrid

  ! Define local variables
  integer :: hdferr
  integer :: mpirank, mpisize, mpierr
  integer :: nprocs
  integer(hid_t) :: atomPairNprocs_did
  integer(hid_t) :: atomPairPgrid_did
  integer(hsize_t), dimension(1) :: dims
  integer(hid_t) :: didtype

  ! First open the groups for nprocs and pgrid
  call h5dopen_f(atomPair_fid, "nprocs", atomPairNprocs_did, hdferr)
  if (hdferr /=0) stop 'Failed to open nprocs dataset.'
  call h5dopen_f(atomPair_fid, "pgrid", atomPairPgrid_did, hdferr)
  if (hdferr /=0) stop 'Failed to open nprocs dataset.'

  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

  dims(1) = 1
  call h5dread_f(atomPairNprocs_did, H5T_NATIVE_INTEGER, nprocs, dims, hdferr)
  if (hdferr /=0) stop 'Failed to read nprocs.'
  dims(1) = 2
  call h5dread_f(atomPairPgrid_did, H5T_NATIVE_INTEGER, pgrid, dims, hdferr)
  if (hdferr /=0) stop 'Failed to read pgrid.'

  call MPI_COMM_RANK(MPI_COMM_WORLD, mpirank, mpierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpisize, mpierr)

  if (nprocs /= mpisize) stop 'Distribution not created for this MPI World Size'

  ! Close all datasets and groups opened by this subroutine
  call h5dclose_f(atomPairPgrid_did, hdferr)
  if (hdferr /=0) stop 'Failed to close atomPairPgrid_did'
  call h5dclose_f(atomPairNprocs_did, hdferr)
  if (hdferr /=0) stop 'Failed to close atomPairNprocs_did'

end subroutine readGrid

subroutine readAtomListHDF5(atomPair_fid, atomList)
  
  ! Import any necessary modules.
  use hdf5
  use O_ParallelSetup
  use MPI

  ! Define passed parameters
  integer(hid_t), intent(in) :: atomPair_fid
  type(AtomPair), pointer :: atomList

  ! Define local variables
  integer :: hdferr
  integer :: mpirank, mpisize, mpierr
  integer(hid_t) :: pal_gid
  integer(hid_t) :: atomPairLen_gid, atomPairLen_did
  integer(hid_t) :: atomList_gid, atomList_did
  integer :: atomLen
  character(len=20) :: rankString
  integer(hsize_t), dimension(2) :: atomDims
  integer(hsize_t), dimension(1) :: lenDim 

  integer :: i
  integer, allocatable, dimension(:,:) :: tempAtomList
  type(AtomPair), pointer :: newPair => null()

  call MPI_COMM_RANK(MPI_COMM_WORLD, mpirank, mpierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpisize, mpierr)

  call h5gopen_f(atomPair_fid,"ProcessAtomList", pal_gid, hdferr)
  if (hdferr /=0) stop 'Failed to open ProcessAtomList group.'

  ! Open the group that will tell use how many atom pairs this process has
  call h5gopen_f(pal_gid, "AtomLen", atomPairLen_gid, hdferr)
  if (hdferr /=0) stop 'Failed to open AtomLen group.'
  ! Open the group that contains the atom pair information
  call h5gopen_f(pal_gid, "AtomList", atomList_gid, hdferr)
  if (hdferr /=0) stop 'Failed to open the AtomList group.'

  ! Create a string for the dataset id
  write (rankString, "(I0)") mpirank

  ! Open the corresponding atomLen dataset
  call h5dopen_f(atomPairLen_gid, trim(rankString), atomPairLen_did, hdferr)
  if (hdferr /=0) stop 'Failed to open the atomPairLen dataset.'
  
  ! Read how many atom pairs we need
  lenDim(1)=1
  call h5dread_f(atomPairLen_did, H5T_NATIVE_INTEGER, atomLen, lenDim, hdferr)
  if (hdferr /=0) stop 'Failed to read the atomLen dataset.'

  ! Now that we know how many atoms there are we can allocate the tempAtomList
  ! array, and read it in.
  call h5dopen_f(atomList_gid, trim(rankString), atomList_did, hdferr)
  if (hdferr /=0) stop 'Failed to open the atomList data set.'

  allocate(tempAtomList(2,atomLen))
  tempAtomList(:,:)=0
  atomDims(1)=2
  atomDims(2)=atomLen
  call h5dread_f(atomList_did, H5T_NATIVE_INTEGER, tempAtomList, atomDims, hdferr)
  if (hdferr /=0) stop 'Failed to read the atomList data set.'

  ! Now we can create the linked list that integralSCF will use. We do the first
  ! one manually.
  call initAtomPair(atomList,tempAtomList(1,1),tempAtomList(2,1))
  newPair => atomList
  do i=2,atomLen
    call initAtomPair(newPair%next, tempAtomList(1,i),tempAtomList(2,i))
    newPair=>newPair%next
  end do

  ! Close all groups and datasets opened in this subroutine
  call h5dclose_f(atomList_did,hdferr)
  if (hdferr /=0) stop 'Failed to close atomList_did'
  call h5gclose_f(atomList_gid, hdferr)
  if (hdferr /=0) stop 'Failed to close atomList_gid'
  call h5dclose_f(atomPairLen_did, hdferr)
  if (hdferr /=0) stop 'Failed to close atomPairLen_did'
  call h5gclose_f(atomPairLen_gid, hdferr)
  if (hdferr /=0) stop 'Failed to close atomPairLen_gid'
  call h5gclose_f(pal_gid, hdferr)
  if (hdferr /=0) stop 'Failed to close pal_gid'

  ! Deallocation the temporary array
  deallocate(tempAtomList)


end subroutine readAtomListHDF5


end module O_pOLCAOhdf5
