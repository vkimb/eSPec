program wavepacket_propagation
  use types, only: dp
  use utils, only: chkmerr
  use operators
  use pde
  implicit none
  real(kind=dp), allocatable, dimension(:) :: y0, yf
  real(kind=dp), allocatable, dimension(:) :: phi0, phi1
  real(kind=dp), allocatable, dimension(:) :: x
  
  real(kind=dp) :: diff_step, xmax, xmin
  integer :: problem_size, problem_dimensions, number_of_wavefunctions, space_size
  integer :: i, j, k
  integer :: ierr
  

  ! # Set number of dimensions, number of wavefunctions and size of space
  problem_dimensions = 1
  number_of_wavefunctions = 3
  space_size = 500

  ! # Compute problem size
  problem_size = problem_dimensions*2*number_of_wavefunctions*space_size
  ! # Define the points in space to study (1D for now)

  xmax = 200
  xmin = 0
  allocate(x(space_size),stat=ierr); call chkmerr(ierr)
  allocate(phi0(space_size),stat=ierr); call chkmerr(ierr)
  allocate(phi1(space_size),stat=ierr); call chkmerr(ierr)
  x = [(xmin+(i-1)*(xmax-xmin)/(space_size-1), i=1,space_size)]
  diff_step = x(2)-x(1)
  phi0 = 10.0*exp(-((x-(xmax/2.0))/(35))**2)
  phi0 = phi0/sqrt(sum(phi0*phi0*diff_step))
  write(*,*) " Norm = ", sum(phi0*phi0*diff_step)
  phi1(:) = 0.0d00
  
  ! # Set up the problem size, dimensions and times
  call initializePde(&
       problem_size=problem_size,&
       problem_dimensions=problem_dimensions,&
       time_start=0.0d00,&
       time_step=0.55d00,&
       time_stop=400.0d00)

  ! # Pick up spatial differentiation accuracy
  call initializeOperators(4, diff_step) ! differentiain accuracy, step size in space

  ! ##############################
  ! ### Prepare wavefunctions  ###
  ! ##############################
  ! # Here we set the number of wavefunctions
  call PrepareWavefunctions(number_of_wavefunctions=number_of_wavefunctions, &
       space_size=space_size)

  ! # And then initalize them, one by one #

  ! # 1st wf
  call InitializeWavefunction(idx=1, space_size=space_size, &
       WF_Real=phi0, WF_Imag=0.0d00*phi0, &
       V_Real=(x-xmax/2.0d0)**2, V_Imag=0.0d00*phi0, Gamma=0.5d00)

  ! # 2nd wf
  call InitializeWavefunction(idx=2, space_size=space_size, &
       WF_Real=phi1, WF_Imag=0.0d00*phi1, &
       V_Real=30.0d00 + (x-xmax/1.5d0)**2, V_Imag=0.0d00*phi0, Gamma=0.3d00) ! We shift the potential up and displace it slightly
  
  call InitializeWavefunction(idx=3, space_size=space_size, &
       WF_Real=phi1, WF_Imag=0.0d00*phi1, &
       V_Real=40.0d00+(x-xmax/1.85d0)**2, V_Imag=0.0d00*phi0, Gamma=0.5d00) ! 
  ! # Initialize transition dipole moments between every pair of wavefunctions
  ! You can call it only once for each pair
  call InitializeTransitionDipoles(mu_in=10.0d00, i=1, j=3)
  call InitializeTransitionDipoles(mu_in=0.0d00, i=1, j=2)
  call InitializeTransitionDipoles(mu_in=1.0d00, i=2, j=3)
  
  ! ######################
  ! ### Prepare fields ###
  ! ######################

  ! # First, set number of fields
  call PrepareFields(number_of_fields=2)

  ! # Second, initalize them one by one
  call InitializeField(idx=1, t0=150.0d00, &
       omega=40d0, width=10.0d00, A=0.04d00, phase_shift=0.0d00)

  call InitializeField(idx=2, t0=250.0d00, &
       omega=10d0, width=10.0d00, A=0.02d00, phase_shift=0.0d00)

  call startPDEsolver()

end program wavepacket_propagation

subroutine FCVFUN(T, Y, YDOT, IPAR, RPAR, IER)
  use types, only: dp
  use pde
  real(kind=dp), dimension(psize), intent(in) :: Y
  real(kind=dp), dimension(psize), intent(out) :: YDOT
  real(kind=dp), dimension(1), intent(in) :: RPAR
  real(kind=dp), intent(in)         :: t
  integer, dimension(1), intent(in) :: IPAR
  integer :: i, j, k
  integer, intent(out) :: IER

  call SetRightHandSide(T,Y,YDOT)
  IER = 0
end subroutine FCVFUN
