module O_Setup

contains

subroutine setupSCF

   ! Import the precision variables
   use O_Kinds

   ! Import the necessary modules.
   use O_SetupHDF5,   only: initSetupHDF5, closeSetupHDF5
   use O_CommandLine, only: parseSetupCommandLine
   use O_Input,       only: parseInput
   use O_Lattice,     only: initializeLattice, initializeFindVec, cleanUpLattice
   use O_KPoints,     only: numKPoints, computePhaseFactors, cleanUpKPoints
   use O_Basis,       only: renormalizeBasis, cleanUpBasis
   use O_ExchangeCorrelation, only: maxNumRayPoints, getECMeshParameters, &
                                  & makeECMeshAndOverlap, cleanUpExchCorr
   use O_ElectroStatics,      only: makeElectrostatics
   use O_GaussianRelations,   only: makeAlphaDist, makeAlphaNucDist, &
                                  & makeAlphaPotDist, cleanUpGaussRelations
   use O_IntegralsSCF,        only: allocateIntegralsSCF, gaussOverlapOL, &
                                  & gaussOverlapKE, gaussOverlapNP, &
                                  & elecPotGaussOverlap
   use O_AtomicSites, only: coreDim, valeDim, cleanUpAtomSites
   use O_AtomicTypes, only: cleanUpRadialFns, cleanUpAtomTypes
   use O_PotSites,    only: cleanUpPotSites
   use O_PotTypes,    only: cleanUpPotTypes
   use O_Potential,   only: cleanUpPotential
   use O_CoreCharge,  only: makeCoreRho
   use O_TimeStamps,  only: initOperationLabels

   ! Import the HDF5 module.
   use HDF5

   ! Import the necessary MPI files
   use MPI
   use O_ParallelSubs

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define global mpi parameters
   integer :: mpiRank
   integer :: mpiSize
   
   ! Define cvOL global array file handle
   integer, allocatable, dimension(:) :: localCVoL

   ! Define tau parameters
!   integer profiler(2) / 0, 0 /
!   save profiler

   ! Define error detectors for hdf, tau, and mpi.
   integer :: hdferr
!   integer :: tauerr
   integer :: mpierr


   ! Initialize the MPI interface
   call MPI_INIT (mpierr)
   call MPI_COMM_RANK (MPI_COMM_WORLD,mpiRank,mpierr)
   call MPI_COMM_SIZE (MPI_COMM_WORLD,mpiSize,mpierr)

   ! Initialize tau timer and start it.

!   call TAU_PROFILE_INIT()
!   call TAU_PROFILE_TIMER(profiler, 'setup')
!   call TAU_PROFILE_START(profiler)
!   call TAU_PROFILE_SET_NODE(0)

   ! Initialize the logging labels.
   call initOperationLabels

   ! Parse the command line parameters
   call parseSetupCommandLine


   ! Read in the input to initialize all the key data structure variables.
   call parseInput


   ! Find specific computational parameters not EXPLICITLY given in the input
   !   file.  These values can, however, be easily determined from the input
   !   file.
   call getImplicitInfo


   ! Create real-space and reciprocal-space super lattices out of the primitive
   !   lattice.  These "supercells" must be big enough so as to include all the
   !   points within a sphere bounded by the negligability limit.  Points
   !   outside the sphere are considered negligable.
   call initializeLattice (1)


   ! Setup the necessary data structures so that we can easily find the lattice
   !   vector that is closest to any other arbitrary vector.
   call initializeFindVec


   ! Compute the kpoint phase factors.
   call computePhaseFactors


   ! Renormalize the basis functions
   call renormalizeBasis


   ! Determine the parameters for the exchange correlation mesh.
   call getECMeshParameters


   ! Now, the dimensions of the system are known.  Therefor we can
   !   initialize the HDF5 file structure format, and datasets.
   if (mpiRank == 0) then
      call initSetupHDF5 (maxNumRayPoints)
   endif

   ! Construct the exchange correlation overlap matrix, and sampling field.
   call makeECMeshAndOverlap()
   call MPI_Barrier(mpierr)

   ! Construct matrix operators and integral vectors that are used to later
   !   determine the electrostatic potential.
   call makeElectrostatics()
   call MPI_Barrier(mpierr)

   ! Create the alpha distance matrices.
   call makeAlphaDist


   ! Allocate space to be used for each of the integrals.  The 1 for the
   !   valeVale cases is a place holder for these matrices because they will
   !   later (main.exe) consider spin.
   ! Setup the CVoL local arrays
   allocate(localCVoL)
   call setupCVoL(localCVoL,coreDim,valeDim,numKPoints)

   ! Calculate the matrix elements of the overlap between all LCAO Bloch
   !   wave functions.
   call gaussOverlapOL(localCVoL)
   call MPI_Barrier(mpierr)

   ! Calculate the matrix elements of the kinetic energy between all LCAO Bloch
   !   wave functions.
   call gaussOverlapKE(localCVoL)
   call MPI_Barrier(mpierr)


   ! Create the alpha distance matrix with nuclear alpha factor
   call makeAlphaNucDist()

   ! Calculate the matrix elements of the overlap between all LCAO Bloch
   !   wave functions and the nuclear potentials.
   call gaussOverlapNP(localCVoL)
   call MPI_Barrier(mpierr)

   ! Create the alpha distance matrix with potential alpha factor
   call makeAlphaPotDist()

   ! Calculate the matrix elements of the overlap between all LCAO Bloch
   !   wave functions and the potential site potential alphas.
   call elecPotGaussOverlap(localCVoL)
   call MPI_Barrier(mpierr)

   ! Now that all the matrices are done being made we can deallocate the
   !   data structures that were used in all the above subroutines but are not
   !   necessary now.
   call cleanUpBasis
   call cleanUpGaussRelations


   ! Construct a vector describing the core charge density since it will not
   !   change throughout the SCF iterations because of core orthogonalization.
   if ((coreDim /= 0) .and. (mpiRank == 0)) then
      call makeCoreRho()
   endif
   call MPI_Barrier(mpierr)


   ! Deallocate the data component of the atomType data structure that holds
   !   the radial functions since they are no longer needed.
   call cleanUpRadialFns


   ! Close all the parts of the setup HDF5 file.
   if (mpiRank == 0) then
      call closeSetupHDF5
   endif

   ! Deallocate the coreValeOL global array
   call MPI_Barrier(mpierr)
   deallocate(localCVoL)

   ! Deallocate all the other as of yet un-deallocated arrays.
   call cleanUpAtomTypes
   call cleanUpAtomSites
   call cleanUpPotTypes
   call cleanUpPotSites
   call cleanUpKPoints
   call cleanUpExchCorr
   call cleanUpLattice
   call cleanUpPotential

   ! End the tau timer
!   call TAU_PROFILE_STOP(profiler)

   ! End the MPI interface
!   call MPI_FINALIZE (mpierr)

   ! Close the output file
   close (20)

   ! Open a file to signal completion of the program.
   if (mpiRank == 0) then
      open (unit=2,file='fort.2',status='unknown')
      close (2)
   endif

   ! Close the GA interface
   call MPI_Barrier(mpierr)

   ! End the MPI interface
   call MPI_FINALIZE (mpierr)

end subroutine setupSCF

subroutine getImplicitInfo

   ! Import necessary modules.
   use O_ExchangeCorrelation
   use O_AtomicTypes
   use O_AtomicSites
   use O_PotSites
   use O_PotTypes
   use O_Lattice
   use O_KPoints
   use O_Potential
   use O_TimeStamps

   implicit none

   call timeStampStart(2)

   ! Subroutines need to be called in this order due to data dependencies.
   call makeSampleVectors

   call getAtomicTypeImplicitInfo
   call getAtomicSiteImplicitInfo
   call getPotSiteImplicitInfo
   call getPotTypeImplicitInfo

   call getRecipCellVectors
   call convertKPointsToXYZ

   call initPotStructures

   call timeStampEnd(2)

end subroutine getImplicitInfo

end module O_Setup
