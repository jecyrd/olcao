module O_pMainHDF5
  ! Iport any necessary definition modules.
  
  ! Make sure that no funny variables are defined.
  implicit none
  
  contains

subroutine accessPGridHDF5(pgrid)
  ! Import necessary modules.
  use hdf5
  use O_ParallelSetup

  ! Make sure that no funny variables are defined
  implicit none

  ! Define passed parameters
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

end module O_pMainHDF5
