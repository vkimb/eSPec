! ################################################
! ### Module PDE: Routines to set up and run   ###
! ###             the wavepacket propagation   ###
! ################################################
module pde
  use types
  use utils
  use operators
  use omp_lib
  use io, only: writePlainDataset
  implicit none
  integer :: psize, pdims  ! # Size and dimensionality of the problem
  logical :: PSIZE_SET = .FALSE. ! # Is problem size initialized?
  logical :: SOLVER_WF_ALLOCATED = .FALSE. ! # Are the wavefunctions allocated?
  logical :: FIELDS_ALLOCATED = .FALSE.
  logical, allocatable, dimension(:) :: FIELDS_INITIALIZED
  logical, allocatable, dimension(:) :: SOLVER_WF_INITIALIZED  ! # Are the wavefunctions initialized?
  ! # Solver context stores wavefunctions and time-independent potentials,
  ! # relaxation constants, as well as start and stop time for propagation
  type(solver_context) :: solver
  type(field_type)     :: fields
  ! # Define publically-accessible subroutines 
  public :: startPDEsolver
  public :: initializePde
  public :: psize, pdims
  public :: SetRightHandSide
  public :: PrepareWavefunctions, PrepareFields
  public :: InitializeWavefunction, InitializeField
  public :: InitializeTransitionDipoles
  private ! # Set everything non-public to be private
  

contains


  ! #############################################
  ! # Allocate the wavefunctions and potentials #
  ! #############################################
  subroutine PrepareWavefunctions(number_of_wavefunctions, space_size)
    integer, intent(in) :: number_of_wavefunctions, space_size
    integer :: i,j,k,ierr
    ! # Store the number of wavefunctions
    solver%number_of_wavefunctions=number_of_wavefunctions
    solver%space_size=space_size
    if (.not.(SOLVER_WF_ALLOCATED)) then
       ! # Get memory for wavefunction storage
       allocate(solver%f(number_of_wavefunctions),stat=ierr); call chkmerr(ierr)
       ! # Get memory for df/dt to store right-hand side
       allocate(solver%dfdt(number_of_wavefunctions),stat=ierr); call chkmerr(ierr)
       ! # And get memory for the wavefunctions
       do i=1,number_of_wavefunctions
          ! # i'th wavefunction 
          allocate(solver%f(i)%WF(space_size),stat=ierr); call chkmerr(ierr)
          allocate(solver%dfdt(i)%WF(space_size),stat=ierr); call chkmerr(ierr)
          ! # i'th potential, time-independent real part
          allocate(solver%f(i)%VR(space_size),stat=ierr); call chkmerr(ierr)
          ! # i'th potential, time-independent imaginary part
          allocate(solver%f(i)%VI(space_size),stat=ierr); call chkmerr(ierr)
       end do

       ! # Allocate transition dipole moment storage
       allocate(solver%mu(number_of_wavefunctions, number_of_wavefunctions),stat=ierr); call chkmerr(ierr)
       solver%mu(:,:) = 0.0d00
       SOLVER_WF_ALLOCATED=.TRUE.
       allocate(SOLVER_WF_INITIALIZED(number_of_wavefunctions))
       SOLVER_WF_INITIALIZED(:) = .FALSE.
    end if
  end subroutine PrepareWavefunctions

  ! #################################################
  ! # Initialize values for the idx'th wavefunction #
  ! #################################################
  subroutine InitializeWavefunction(&
       idx, space_size,  WF_Real, WF_Imag, &
       V_Real, V_Imag, Gamma)
    ! # idx: index of wavefunction to initialize
    ! # space_size: size of the space along one dimension
    integer, intent(in) :: idx, space_size
    ! # Real and imaginary parts of the wavefunction
    real(kind=dp), intent(in), dimension(space_size) :: WF_Real, WF_Imag
    ! # Real and imaginary parts of the time-independent potential
    real(kind=dp), intent(in), dimension(space_size) :: V_Real, V_Imag
    ! # Relaxation constant
    real(kind=dp), intent(in) :: Gamma
    ! # Check if the storage was allocated before 
    if (.not.(SOLVER_WF_ALLOCATED)) then
       write(*,*) " Error: wavefunction storage is not allocated "
       stop
    end if
    ! # Now, if the wavefunctions were not initialized, initalize them
    if (.not.(SOLVER_WF_INITIALIZED(idx))) then
       ! # Initialize the wavefunction
       solver%f(idx)%WF(:)%Re=WF_Real(:)
       solver%f(idx)%WF(:)%Im=WF_Imag(:)
       ! # Initialize the time-independent potential
       solver%f(idx)%VR(:)=WF_Real(:)
       solver%f(idx)%VI(:)=WF_Imag(:)
       ! # Initialize the relaxation constant
       solver%f(idx)%Gamma=Gamma
       SOLVER_WF_INITIALIZED(idx)=.TRUE.
    end if
  end subroutine InitializeWavefunction

  ! ###############################
  ! # Allocate the applied fields #
  ! ###############################
  subroutine PrepareFields(number_of_fields)
    integer, intent(in) :: number_of_fields
    integer :: ierr
    ! # Store the number of fields
    solver%number_of_fields = number_of_fields
    if (.not.FIELDS_ALLOCATED) then
       allocate(fields%A(number_of_fields),stat=ierr); call chkmerr(ierr)
       allocate(fields%t0(number_of_fields),stat=ierr); call chkmerr(ierr)
       allocate(fields%omega(number_of_fields),stat=ierr); call chkmerr(ierr)
       allocate(fields%width(number_of_fields),stat=ierr); call chkmerr(ierr)
       allocate(fields%phase_shift(number_of_fields),stat=ierr); call chkmerr(ierr)
       fields%number_of_fields = number_of_fields
       FIELDS_ALLOCATED = .TRUE.
       allocate(FIELDS_INITIALIZED(number_of_fields),stat=ierr); call chkmerr(ierr)
    end if
  end subroutine PrepareFields

  ! ##########################
  ! # Set the applied fields #
  ! ##########################
  ! # Initialize the parameters for a given field (index "idx")
  subroutine InitializeField(idx, t0, omega, width, A, phase_shift)
    real(kind=dp), intent(in) :: t0, omega, width, A, phase_shift
    integer, intent(in)  :: idx
    if (.not.FIELDS_INITIALIZED(idx)) then
       fields%A(idx) = A
       fields%t0(idx) = t0
       fields%width(idx) = width
       fields%phase_shift(idx) = phase_shift
       FIELDS_INITIALIZED(idx) = .TRUE.
    end if
  end subroutine InitializeField

  subroutine InitializeTransitionDipoles(mu_in, i, j)
    real(kind=dp), intent(in)  :: mu_in
    integer, intent(in)        :: i,j
    ! # Symmetric coupling
    solver%mu(i,j) = mu_in
    solver%mu(j,i) = mu_in
  end subroutine InitializeTransitionDipoles
  
  ! #####################
  ! # Compute couplings #
  ! #####################
  ! # Computes coupling between a pair of wavefunctions
  ! # in presence of all applied fields at time t
  pure function Coupling(i,j,t)
    integer, intent(in) :: i, j
    real(kind=dp), intent(in) :: t
    integer :: ifield
    real(kind=dp) :: Coupling
    Coupling = 0.0d00
    Coupling = -solver%mu(i,j)*applyFields(t, fields)
  end function Coupling


  ! ##########################################
  ! # Set up the right-hand side for the PDE #
  ! ##########################################
  subroutine SetRightHandSide(T,Y,YDOT) ! 
    real(kind=dp), dimension(psize), intent(in) :: Y
    real(kind=dp), dimension(psize), intent(out) :: YDOT
    real(kind=dp), intent(in)   :: T
    real(kind=dp), dimension(&
         solver%number_of_wavefunctions,&
         solver%number_of_wavefunctions) :: J ! Couplings stored there
    integer :: iwf, jwf, ispace

    ! # First, pack the wavefunctions from the vector
    call packWavefunctions(Y)

    ! # Then let's compute the couplings
    J(:,:) = 0.0d00
    do iwf=1,solver%number_of_wavefunctions
       do jwf=1,solver%number_of_wavefunctions
          if (iwf.ne.jwf) then
             J(iwf,jwf) = Coupling(iwf,jwf,T)
          end if
       end do
    end do
    
    ! # Now, compute the right-hand side for every wavefunction
    do iwf=1,solver%number_of_wavefunctions
       ! # for every point in space
       !$OMP PARALLEL PRIVATE(ispace)
       !$OMP DO
       do ispace=1,solver%space_size
          ! ############
          ! # Real RHS #
          ! ############
          ! # Set right-hand side without couplings
          solver%dfdt(iwf)%WF(ispace)%Re = &
               + KineticEnergy1D(solver%f(iwf)%WF(:)%Im,ispace)  & ! Kinetic energy operator on imaginary part of the WF at point ispace
               + solver%f(iwf)%VR(ispace)*solver%f(iwf)%WF(ispace)%Im & ! Real part of the time-independent potential
               + solver%f(iwf)%VI(ispace)*solver%f(iwf)%WF(ispace)%Re & ! Imaginary part of the time-independent potential
               - solver%f(iwf)%Gamma*solver%f(iwf)%WF(ispace)%Im        ! Relaxation part
          ! # Iterate over every other wavefunction and compute couplings
          do jwf=1,solver%number_of_wavefunctions ! 
             if (jwf.ne.iwf) then
                solver%dfdt(iwf)%WF(ispace)%Re = solver%dfdt(iwf)%WF(ispace)%Re &
                     + J(iwf,jwf)*solver%f(jwf)%WF(ispace)%Im
             end if
          end do
          
          ! ############
          ! # Imag RHS #
          ! ############
          solver%dfdt(iwf)%WF(ispace)%Im = &
               - KineticEnergy1D(solver%f(iwf)%WF(:)%Re,ispace)  & ! Kinetic energy operator on real part of the WF at point ispace
               - solver%f(iwf)%VR(ispace)*solver%f(iwf)%WF(ispace)%Re & ! Real part of the time-independent potential
               + solver%f(iwf)%VI(ispace)*solver%f(iwf)%WF(ispace)%Im & ! Imaginary part of the time-independent potential
               + solver%f(iwf)%Gamma*solver%f(iwf)%WF(ispace)%Re        ! Relaxation part
          ! # Iterate over every other wavefunction and compute couplings
          do jwf=1,solver%number_of_wavefunctions ! 
             if (jwf.ne.iwf) then
                solver%dfdt(iwf)%WF(ispace)%Im = solver%dfdt(iwf)%WF(ispace)%Im &
                     - J(iwf,jwf)*solver%f(jwf)%WF(ispace)%Re
             end if
          end do
       end do
       !$OMP END DO
       !$OMP END PARALLEL
    end do
    
    ! # Finally, unpack the computed right-hand side from dfdt vectors
    call unpackRightHandSide(YDOT)
  end subroutine SetRightHandSide

  ! # Initalize SUNDIALS and start the solution process
  subroutine startPDEsolver()
    use iso_c_binding
    real(kind=dp), dimension(psize) :: Y0
    real(kind=dp), dimension(psize) :: y_current
    real(kind=dp), dimension(1) :: RPAR=0.0d00
    integer(kind=C_LONG), dimension(1)       :: IPAR=0

    real(kind=dp) :: dt, t, tout, cpu_tstart, cpu_tstop
    integer  :: i, j, k
    integer  :: ierr
    integer(kind=C_LONG), dimension(21) :: IOUT
    real(kind=dp), dimension(6) :: ROUT
    
    ! ##############################################
    ! # Check if everything is allocated and so on #
    ! ##############################################
    ! # Check if module is initalized properly
    call chkset(PSIZE_SET)
    call CheckPDESetup()

    ! # Initialize NVECTOR module (serial)
    call FNVINITS(1,psize,ierr); call chkfcv(ierr)


    ! # Unpack the wavefunctions into a vector
    call unpackWavefunctions(Y0)
    ! # Allocate CVODE internals
    call FCVMALLOC(solver%pde_time_start, Y0, 2, & ! 2 = BDF, 1 = Adams
         2, & ! 1 = functional iteration, 2 = Newton iteration
         1, & ! 1 = scalar absolute tolerance (?), 2 = array absolute tolerance (?), 3 = user-defined
         1.0d-8, &! Relative tolerance
         1.0d-8, &! Absolute tolerance
         IOUT, ROUT, & ! Integer and real optional outputs
         IPAR, RPAR, & ! Integer and real optional parameters
         ierr); call chkfcv(ierr)

    ! # Set diagonal Jacobian approximation
    call FCVDIAG(ierr); call chkfcv(ierr)

    ! # Use SPGMR solver
    call FCVSPGMR(0,& ! 0 = for no preconditioner (problematic!)
         1, & ! Modified GS orthogonalization
         5000, & ! Maximum size of Krylov subspace
         1.0d-8, & ! Tolerance
         ierr); call chkfcv(ierr)


    ! #############################
    ! # Prepare to run the solver #
    ! #############################
    dt = solver%pde_time_step
    t  = solver%pde_time_start
    write (*,*) " **************************** "
    write (*,*) " * Starting the propagation * "
    write (*,*) " **************************** "
    write (*,"(A,F10.5,A,F10.5)") " From tstart = ", solver%pde_time_start, " to tstop = ", solver%pde_time_stop
    write (*,"(A,F10.5,A)") " Solutions will be saved at dt = ", solver%pde_time_step, " intervals"
    ! # Save initial data
    call writePlainDataset(t, y0, psize)
    ! # Solve the system of equations
    do while (t<=solver%pde_time_stop)
       ! # Since solution at the initial time is known,  
       ! # we're stepping forward and ask for the solution
       ! # at time step t + dt
       t = t + dt
       cpu_tstart = omp_get_wtime()
       call FCVODE(t, tout, y_current, &
            1,  & ! 1 = overshoot and interpolate to t,  2 = one-step mode
            ierr);
       cpu_tstop = omp_get_wtime()
       write(*,"(A,F10.4,A,F12.8,A,F12.6,A)") " t = ", t, " ; Field = ",applyFields(t, fields), " ( ",cpu_tstop-cpu_tstart," s)"
       
       call writePlainDataset(t, y_current, psize)
       ! do j=1,psize
       !    write(*,"(F12.6,A)",advance="no") y_current(j), " "
       ! end do
       
       !write(*,"(A)") " ;"
       if (ierr<0) then
          write(*,*) " ERROR HAS OCCURED DURING PROPAGATION: ierr = ", ierr
          write(*,*) " * EXITING * "
          stop
       end if
    end do
    
  end subroutine startPDEsolver

  ! # Initialize the module by setting problem size
  subroutine initializePde(problem_size, problem_dimensions, time_start, time_step, time_stop)
    integer, intent(in)        :: problem_size, problem_dimensions
    real(kind=dp), intent(in)  :: time_start, time_step, time_stop
    psize = problem_size
    pdims = problem_dimensions
    if ((psize/pdims)<10) then
       write(*,*) "Number of spatial points cannot be smaller then 10"
       stop
    end if
    
    solver%pde_time_start = time_start
    solver%pde_time_stop  = time_stop
    solver%pde_time_step  = time_step
    PSIZE_SET = .TRUE.
  end subroutine initializePde

  subroutine CheckPDESetup()
    integer :: i
    logical :: passed=.TRUE.
    ! # First, check if wavefunctions are set up
    do i=1,solver%number_of_wavefunctions
       ! # If one wavefunction is unitialized, we stop
       if (.not.(SOLVER_WF_INITIALIZED(i))) then
          passed = .FALSE.
          write(*,*) " Wavefunction at ", i, " is not initialized; exiting"
          stop
       end if
    end do
    do i=1,solver%number_of_fields
       ! # If one wavefunction is unitialized, we stop
       if (.not.(FIELDS_INITIALIZED(i))) then
          passed = .FALSE.
          write(*,*) " Wavefunction at ", i, " is not initialized; exiting"
          stop
       end if
    end do
  end subroutine CheckPDESetup

  subroutine unpackWavefunctions(vector_out)
    real(kind=dp), dimension(psize), intent(out) :: vector_out
    integer :: iwf, ispace, ipart, idx
    idx=0
    do iwf=1,solver%number_of_wavefunctions
       ! # Unpack real part
       do ispace=1,solver%space_size
          idx=idx+1
          vector_out(idx) = solver%f(iwf)%WF(ispace)%Re
       end do
       ! # Unpack imaginary part
       do ispace=1,solver%space_size
          idx=idx+1
          vector_out(idx) = solver%f(iwf)%WF(ispace)%Im
       end do
    end do
    if (idx.ne.psize) then
       write(*,*) " Error in unpackWavefunctions: balance not met (idx.ne.psize) "
       stop
    end if
  end subroutine unpackWavefunctions
  
  subroutine packWavefunctions(vector_in)
    real(kind=dp), dimension(psize), intent(in) :: vector_in
    integer :: iwf, ispace, ipart, idx
    idx=0
    do iwf=1,solver%number_of_wavefunctions
       ! # Pack real part
       do ispace=1,solver%space_size
          idx=idx+1
          solver%f(iwf)%WF(ispace)%Re = vector_in(idx)
       end do
       ! # Pack imaginary part
       do ispace=1,solver%space_size
          idx=idx+1
          solver%f(iwf)%WF(ispace)%Im = vector_in(idx)
       end do
    end do
    if (idx.ne.psize) then
       write(*,*) " Error in packWavefunctions: balance not met (idx.ne.psize) "
       stop
    end if
  end subroutine packWavefunctions

  
  subroutine unpackRightHandSide(vector_out)
    real(kind=dp), dimension(psize), intent(out) :: vector_out
    integer :: iwf, ispace, ipart, idx
    idx=0
    do iwf=1,solver%number_of_wavefunctions
       do ispace=1,solver%space_size
          idx=idx+1
          vector_out(idx) = solver%dfdt(iwf)%WF(ispace)%Re
       end do
       ! # Unpack imaginary part
       do ispace=1,solver%space_size
          idx=idx+1
          vector_out(idx) = solver%dfdt(iwf)%WF(ispace)%Im
       end do
    end do
    if (idx.ne.psize) then
       write(*,*) " Error in unpackRightHandSide: balance not met (idx.ne.psize) "
       stop
    end if
  end subroutine unpackRightHandSide

end module pde

  
  
