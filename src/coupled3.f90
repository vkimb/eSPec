subroutine coupled3(dim,nd,n,np,xp,xi,sh,shm,U0_in,V0_in,VPOT_in,VMINB,VMINC,gamma,&
     tdipol,omega,e0,tp,td,t0,sni,kl,ti,tf,dt,nshot,mxdct,lmtreort,eigA,eigB,eigC,nstates,&
     prteigvc2,m_fourier,t_ierr)
  !------------------------------------------------------
  !     
  !     subroutine to compute wavepacket propagation 
  !     on 3 potential energy surfaces coupled by strong optical or X-ray fields
  !     within the Rotating Wave Approximation
  !     
  !     Vinicius Vaz da Cruz
  !     Stockholm, 23 of November of 2015
  !     
  !------------------------------------------------------

  use rhside
  use iso_c_binding
  !include ’fftw3.f03’


  implicit none


  integer, parameter :: dp=kind(1.0d00)
  !------external arguments (from eSPec)
  integer, intent(in) :: n,nd,nshot,mxdct,lmtreort,nstates,m_fourier
  integer, dimension(3), intent(in) :: np,kl
  character (len=*) :: dim
  real(kind=dp) VMINB,VMINC,t,dt, ti, tf,EF,xp,tout,start_t,stop_t,t_total
  real(kind=dp), dimension(3) :: omega,e0,tp,td,t0,sni,tdipol,shm,sh,xi,gamma
  real(kind=dp), dimension(mxdct), intent(in) :: U0_in,V0_in
  real(kind=dp), dimension(mxdct), intent(in) :: VPOT_in
  real(kind=dp), dimension(3*n) :: VPOT
  real(kind=dp), dimension(6*n) :: Y, F, HY
  real(kind=dp), intent(in), dimension(mxdct,lmtreort) :: eigA,eigB,eigC
  logical, intent(in) :: prteigvc2
  integer, intent(out) :: t_ierr

  !-----local variables
  integer i,j,k,ntime,ierr,dim_int
  integer(kind=C_LONG), dimension(1) :: IPAR=0
  integer(kind=C_LONG), dimension(21) :: IOUT
  integer(kind=C_LONG) :: j_long

  real(kind=dp), dimension(6) :: ROUT
  real(kind=dp), dimension(1) :: RPAR=0.0d0
  real(kind=c_double), dimension(3) :: detun
  real(kind=dp), dimension(3) :: G,CSD,SND
  real(kind=dp), dimension(8) :: WKG
  real(kind=dp) :: E1,E2,var,norm_last,norm_diff
  real(kind=dp), dimension(3) :: norm
  real(kind=dp), dimension(2,nstates) :: P_A,P_B,P_C
  character*20 NEWNAM1,NEWNAM2,NEWNAM3
  character*4 chnum
  character*30 pfmt
  character (len=*), parameter ::  pul=".ENVG "
  real(kind=dp), parameter :: fs2au=41.3411D+0, ev2au=27.2114D+0
  real(kind=dp), dimension(2,2) :: rho
  
  integer ( c_int )  :: c_ntime,n_fourier
  real ( c_double )  :: c_ti, c_stept
  real(kind=C_DOUBLE),dimension(:),allocatable :: re_rho12_t,im_rho12_t, re_rho23_t, im_rho23_t,G12_t,G23_t
  real(kind=C_DOUBLE),dimension(2**m_fourier) :: sigma_xas,sigma_rixs


  !----------------------------------------------

  !-- output files
  open(unit=42,file='pulses.dat')
  open(unit=43,file='total_population.dat')
  open(unit=49,file='rho_G.dat')
  
  write(pfmt,'(A,I4,A)')'(',nstates+1,'(1X,ES18.9))'

  open(unit=44,file='potA_projection.dat')
  open(unit=45,file='potB_projection.dat')
  open(unit=46,file='potC_projection.dat')

  VPOT(:)=VPOT_in(1:3*n)
  t_total = 0.0d00

  if(dim(1:3).EQ.'.1D')then
     dim_int = 0
  else if(dim(1:5).EQ.'.2DCT')then
     dim_int = 2
  else if(dim(1:3).EQ.'.2D')then
     dim_int = 1
  else
     print *, "wrong dimension!!"
     stop
  end if
  
  call initialize_initguess(dim_int,nd,n,np,xp,xi,sh,shm,VPOT,VMINB,&
       VMINC,ti)

  !-----number of propagation steps
  ntime = 1 + (tf - ti)/dt
  c_ntime=ntime;

  !allocate arrays for cross-section calculation
  allocate(re_rho12_t(ntime))
  allocate(im_rho12_t(ntime))
  allocate(re_rho23_t(ntime))
  allocate(im_rho23_t(ntime))
  allocate(G12_t(ntime))
  allocate(G23_t(ntime))

  print*, 'inside coupled3'
  !     debugs
  !     check initial condition
  call PRPT2('check_initwf.dat ', 9, nd, 0, np, xp, xi, sh, u0_in, v0_in)

  !     check potentials
  call PRTPT('potentialA.dat ', 9,nd,0,np,xp,xi,sh,VPOT(1:n))
  call PRTPT('potentialB.dat ', 9,nd,0,np,xp,xi,sh,VPOT(n+1:2*n))
  call PRTPT('potentialC.dat ', 9,nd,0,np,xp,xi,sh,VPOT(2*n +1:3*n))

  print*, 'potential B minimum',VMINB,'a.u.'
  print*, 'potential C minimum',VMINC,'a.u.'

  write(*,*)
  print*,'gamma A',gamma(1)
  print*,'gamma B',gamma(2)
  print*,'gamma C',gamma(3)
  write(*,*)
  print*,'dipole moment d12',tdipol(1)
  print*,'dipole moment d23',tdipol(2)
  write(*,*)
  print*, 'Pulse one parameters: '
  write(*,*),'E0',e0(1),'omega',omega(1),'t0',t0(1),'phase',sni(1),'shape',kl(1)
  print*, 'Pulse two parameters: '
  print*,'E0',e0(2),'omega',omega(2),'t0',t0(2),'phase',sni(2),'shape',kl(1)
  write(*,*)


  !--------set initial condition
  print *, "psize = ", psize
  do j=1,n
     !---  x_1
     j_long = j
     Y(j_long) = U0_in(j)  !Re part
     j_long = j + n
     Y(j_long) = V0_in(j)  !Im part

     !---  x_2
     j_long = j+2*n
     Y(j_long) = 0.0d+0    !Re part
     j_long = j+3*n
     Y(j_long) = 0.0d+0    !Im part
     
     !---  x_3
     j_long = j+4*n     
     Y(j_long) = 0.0d+0    !Re part
     j_long = j+5*n
     Y(j_long) = 0.0d+0    !Im part
  end do
  !------------------------------

  !------- convert to au
  omega(1) = omega(1)/ev2au
  omega(2) = omega(2)/ev2au
  gamma(1) = gamma(1)/ev2au
  gamma(2) = gamma(2)/ev2au
  gamma(3) = gamma(3)/ev2au


  !------detuning from resonances
  detun(1) = omega(1) - VMINB
  detun(2) = omega(2) - (VMINB - VMINC)


  call initialize_field(e0,t0,td,tp,detun,sni,kl,tdipol,gamma)


  !--------time parameters
  print*,'ti ',ti,'tf ',tf,'step',dt
  print*,'number of points',ntime
  write(*,*)

  print*,'initialize n vector'
  call FNVINITS(1,psize,ierr);call chkierr(ierr)
  print*,'allocating'
  call FCVMALLOC(ti*fs2au, Y, 1, & ! 2 = BDF, 1 = Adams
       1, & ! 1 = functional iteration, 2 = Newton iteration
       1, & ! 1 = scalar absolute tolerance (?), 2 = array absolute tolerance (?), 3 = user-defined
       1.0d-16, &! Relative tolerance
       1.0d-16, &! Absolute tolerance
       IOUT, ROUT, & ! Integer and real optional outputs
       IPAR, RPAR, & ! Integer and real optional parameters
       ierr);call chkierr(ierr)

  ! call FCVMALLOC(ti, Y, 2, & ! 2 = BDF, 1 = Adams
  !      1, & ! 1 = functional iteration, 2 = Newton iteration
  !      1, & ! 1 = scalar absolute tolerance (?), 2 = array absolute tolerance (?), 3 = user-defined
  !      1.0d-6, &! Relative tolerance
  !      1.0d-6, &! Absolute tolerance
  !      IOUT, ROUT, & ! Integer and real optional outputs
  !      IPAR, RPAR, & ! Integer and real optional parameters
  !      ierr);call chkierr(ierr)

  print*,'set pre-conditioner'
  call FCVDIAG(ierr);call chkierr(ierr)

  print*,' set the solver'
  ! call FCVSPGMR(0,& ! 0 = for no preconditioner (problematic!)
  !      2, & ! Modified GS orthogonalization
  !      5000, & ! Maximum size of Krylov subspace
  !      1.0d-12, & ! Tolerance
  !      ierr);call chkierr(ierr)
  call FCVLAPACKDENSE(psize,ierr)

  !--- initial state properties
  t = ti
  !-- check norm at t = ti
  call wp3_norm(Y,n,norm)
  write(*,'(A10,6X,A15,3X,A15,3X,A15,8X,A10,8X,A9,2X,A18)') 't','norm1','norm2','norm3','norm','norm diff','iteration time'
  write(*,'(5ES18.9 )')t,norm(1),norm(2),norm(3),norm(1)+norm(2)+norm(3)
  write(43,'(5ES18.9 )')t,norm(1),norm(2),norm(3),norm(1)+norm(2)+norm(3)
  norm_last = norm(1)+norm(2)+norm(3)

  !-- projections at t = ti
  !-- potential A, x_1
  call projection(Y(1:2*n),eigA,n,nstates,mxdct,lmtreort,P_A)
  !-- potential B, x_2
  call projection(Y(2*n + 1:4*n),eigB,n,nstates,mxdct,lmtreort,P_B)
  !-- potential C, x_3
  call projection(Y(4*n + 1:6*n),eigC,n,nstates,mxdct,lmtreort,P_C)

  write(44,pfmt)t,(P_A(1,k)*P_A(1,k) + P_A(2,k)*P_A(2,k),k=1,nstates)
  write(45,pfmt)t,(P_B(1,k)*P_B(1,k) + P_B(2,k)*P_B(2,k),k=1,nstates)
  write(46,pfmt)t,(P_C(1,k)*P_C(1,k) + P_C(2,k)*P_C(2,k),k=1,nstates)

  ! first position must be zero at t=0 naturally.
  re_rho12_t(1)= 0.0e+0
  im_rho12_t(1)= 0.0e+0
  re_rho23_t(1)= 0.0e+0
  im_rho23_t(1)= 0.0e+0
  G12_t(1)= 0.0e+0
  G23_t(1)= 0.0e+0
  write(49,'(7ES21.9)')t,re_rho12_t(1),im_rho12_t(1),re_rho23_t(1),im_rho23_t(1),G12_t(1),G23_t(1)


  !---------time loop
  do i=2,ntime
     t = ti + (i - 1)*dt
     t = t * fs2au

     !---------Andrey's propagation routines   
     call cpu_time(start_t)
     call FCVODE(t, tout, Y, &
          1,  & ! 1 = overshoot and interpolate to t,  2 = one-step mode
          ierr);call chkierr2(ierr)
     call cpu_time(stop_t)
     !--------------------------------------

     !-- reinitialize problem
     call FCVREINIT(t,Y,1,1.0d-16,1.0d-16,ierr); call  chkierr(ierr)
     !-----------------------

     t = t/fs2au

     if(prteigvc2) then
        !-----------printing wavepackets----------------------------
        WRITE(CHNUM,'(I4)')i
        IF(i.LT.10)THEN
           NEWNAM3 = 'ReImC_000'//CHNUM(4:4)//'.dat'
           NEWNAM2 = 'ReImB_000'//CHNUM(4:4)//'.dat'
           NEWNAM1 = 'ReImA_000'//CHNUM(4:4)//'.dat'
        ELSEIF(i.LT.100)THEN
           NEWNAM3 = 'ReImC_00'//CHNUM(3:4)//'.dat'
           NEWNAM2 = 'ReImB_00'//CHNUM(3:4)//'.dat'
           NEWNAM1 = 'ReImA_00'//CHNUM(3:4)//'.dat'
        ELSEIF(i.LT.1000)THEN
           NEWNAM3 = 'ReImC_0'//CHNUM(2:4)//'.dat'
           NEWNAM2 = 'ReImB_0'//CHNUM(2:4)//'.dat'
           NEWNAM1 = 'ReImA_0'//CHNUM(2:4)//'.dat'
        ELSEIF(i.LT.10000)THEN
           NEWNAM3 = 'ReImC_'//CHNUM(1:4)//'.dat'
           NEWNAM2 = 'ReImB_'//CHNUM(1:4)//'.dat'
           NEWNAM1 = 'ReImA_'//CHNUM(1:4)//'.dat'
        ELSE
           WRITE(*,*)'Too many files to be printed! Stopping'
           STOP
        ENDIF

        CALL PRPT2(NEWNAM1,8,ND,t/fs2au,NP,XP,XI,SH,Y(1:n),Y(n+1:2*n))
        CALL PRPT2(NEWNAM2,8,ND,t/fs2au,NP,XP,XI,SH,Y(2*n+1:3*n),Y(3*n+1:4*n))
        CALL PRPT2(NEWNAM3,8,ND,t/fs2au,NP,XP,XI,SH,Y(4*n+1:5*n),Y(5*n+1:6*n))
        !-----------------------------------------------------------------
     end if

     !--- Print field
     E1 = EF(pul,e0(1),t0(1),td(1),tp(1),detun(1),sni(1),kl(1),t)
     E2 = EF(pul,e0(2),t0(2),td(2),tp(2),detun(2),sni(2),kl(2),t)   

     write(42,*)t,E1,E2
     !---------------


     !-- check norm
     call wp3_norm(Y,n,norm)
     norm_diff = dabs(norm_last - (norm(1)+norm(2)+norm(3)) )
     norm_last = norm(1)+norm(2)+norm(3)
     write(*,'(6ES21.9,5X,F5.2,A2)')t,norm(1),norm(2),norm(3),norm(1)+norm(2)+norm(3),norm_diff,stop_t-start_t,'s'
     write(43,'(5ES21.9)')t,norm(1),norm(2),norm(3),norm(1)+norm(2)+norm(3)
     t_total = t_total + stop_t-start_t

     !-----------wavepacket projection on eigenstates------------------
     
     !-- potential A, x_1
     call projection(Y(1:2*n),eigA,n,nstates,mxdct,lmtreort,P_A)

     !-- potential B, x_2
     call projection(Y(2*n + 1:4*n),eigB,n,nstates,mxdct,lmtreort,P_B)

     !-- potential C, x_3
     call projection(Y(4*n + 1:6*n),eigC,n,nstates,mxdct,lmtreort,P_C)

    
     write(44,pfmt)t,(P_A(1,k)*P_A(1,k) + P_A(2,k)*P_A(2,k),k=1,nstates)
     write(45,pfmt)t,(P_B(1,k)*P_B(1,k) + P_B(2,k)*P_B(2,k),k=1,nstates)
     write(46,pfmt)t,(P_C(1,k)*P_C(1,k) + P_C(2,k)*P_C(2,k),k=1,nstates)

     !-----------------------------------------------------------------

     !------- cross section related quantities
     call wp3_rho(Y,n,rho)
     re_rho12_t(i)= rho(1,1)
     im_rho12_t(i)= rho(1,2)
     re_rho23_t(i)= rho(2,1)
     im_rho23_t(i)= rho(2,2)
     G12_t(i)= tdipol(1) * E1
     G23_t(i)= tdipol(2) * E2

     ! to avoid printing issues
     if(G12_t(i).lt.1.0d-99)then
        G12_t(i) = 0.0d+0
     end if
     
     if(G23_t(i).lt.1.0d-99)then
        G23_t(i) = 0.0d+0
     end if


     write(49,'(7ES21.9)')t,re_rho12_t(i),im_rho12_t(i),re_rho23_t(i),im_rho23_t(i),G12_t(i),G23_t(i)
     !-----------------------------------------------------------------


  end do

  
  write(*,*)
  write(*,'(A13,F10.2,A)') 'total time = ',t_total,'s'
  write(*,*)'propagation finished!'
  write(*,*)
  write(*,*)

  write(*,*) 'Computing Cross-sections'
  c_ti=ti; c_stept=dt
  call f_csection (c_ti,c_stept,c_ntime,re_rho12_t,im_rho12_t,re_rho23_t,im_rho23_t,G12_t,G23_t,detun,m_fourier)
  

  t_ierr = 0

  !------------------------------------------------------
  write(*,*)
  write(*,*)'exiting coupled3'
  write(*,*)
  !------------------------------------------------------
end subroutine coupled3


!---------------------------------
!  wavepacket projection routine
!---------------------------------

subroutine projection(Y,eig,n,nstates,mxdct,lmtreort,P)
  use iso_c_binding
  integer, parameter :: dp=kind(1.0d00)
  integer, intent(in)  :: n,nstates,mxdct,lmtreort
  real(kind=dp), dimension(2*n),intent(in) :: Y
  real(kind=dp), dimension(mxdct,lmtreort),intent(in) :: eig
  real(kind=dp), dimension(2,nstates),intent(out) :: P

  integer(kind=C_LONG) i,k

  P = 0.0d00

  do i=1,n
     do k =1,nstates
        P(1,k) = P(1,k) + Y(i)   * eig(i,k) ! real part
        P(2,k) = P(2,k) + Y(i+n) * eig(i,k) ! imaginary part
     end do
  end do

end subroutine projection


!---------------------------------
!  wavepacket norm function
!---------------------------------
subroutine wp3_norm(Y,n,norm)
  use iso_c_binding
  integer, parameter :: dp=kind(1.0d00)
  integer, intent(in)  :: n
  real(kind=dp), dimension(6*n),intent(in) :: Y

  integer(kind=C_LONG) i,jr,ji
  real(kind=dp), dimension(3), intent(out) :: norm

  norm = 0.0d00

  do i=1,n
     jr = i; ji = n + i;
     norm(1) = norm(1) + Y(jr)*Y(jr) + Y(ji)*Y(ji)

     jr = 2*n + i; ji = 3*n + i;
     norm(2) = norm(2) + Y(jr)*Y(jr) + Y(ji)*Y(ji)

     jr = 4*n + i; ji = 5*n + i;
     norm(3) = norm(3) + Y(jr)*Y(jr) + Y(ji)*Y(ji)
  end do

  do i=1,3
     if(norm(i).LT.1.0d-99)then
        norm(i) = 0.0d0
     end if
  end do

end subroutine wp3_norm


!---------------------------------
!  computes the scalar products < x_1 | x_2 > and < x_2 | x_3 >
!---------------------------------
subroutine wp3_rho(Y,n,rho)
  use iso_c_binding
  integer, parameter :: dp=kind(1.0d00)
  integer, intent(in)  :: n
  real(kind=dp), dimension(6*n),intent(in) :: Y

  integer(kind=C_LONG) i,jr,ji,kr,ki
  real(kind=dp), dimension(2,2), intent(out) :: rho

  rho = 0.0d00

  do i=1,n
     ! rho_12
     jr = i; ji = n + i;
     kr = 2*n + i; ki = 3*n + i;
     rho(1,1) = rho(1,1) + Y(jr)*Y(kr) + Y(ji)*Y(ki) !real part
     rho(1,2) = rho(1,2) + Y(ji)*Y(kr) - Y(jr)*Y(ki) !imaginary part
     
     ! rho_23
     jr = 2*n + i; ji = 3*n + i;
     kr = 4*n + i; ki = 5*n + i;
     rho(2,1) = rho(2,1) + Y(jr)*Y(kr) + Y(ji)*Y(ki) !real part
     rho(2,2) = rho(2,2) + Y(ji)*Y(kr) - Y(jr)*Y(ki) !imaginary part
        
  end do

end subroutine wp3_rho


!---------------------------------
!  function used by sundials to compute YDOT at a given time t (righthand side)
!---------------------------------
subroutine FCVFUN(T, Y, YDOT, IPAR, RPAR, IER)
  use rhside
  integer, parameter :: dp=kind(1.0d00)
  real(kind=dp), dimension(psize), intent(in) :: Y
  real(kind=dp), dimension(psize), intent(out) :: YDOT
  real(kind=dp), dimension(1), intent(in) :: RPAR
  real(kind=dp), intent(in)         :: t
  integer, dimension(1), intent(in) :: IPAR
  integer :: i, j, k
  integer, intent(out) :: IER

  !     call SetRightHandSide(T,Y,YDOT)
  !print *, "Entering compute_rhside at t=",t
  call compute_rhside(T,Y,YDOT)


  IER = 0
end subroutine FCVFUN

!wrapping fortran function for the cross-sections c routine
subroutine f_csection(ti,stept,n,re_wr12, im_wr12, re_wr23, im_wr23,G12,G23,detun,n_fourier)
  use iso_c_binding
  integer ( c_int ),  intent(in) :: n, n_fourier
  integer ( c_int ), parameter :: k=3
  real ( c_double ),  intent(in) :: ti, stept
  real ( c_double ),  intent(in), dimension(n)  :: re_wr12, im_wr12, re_wr23, im_wr23
  real ( c_double ),  intent(in), dimension(n)  ::  G12,  G23
  real ( c_double ),  intent(in), dimension(k)  :: detun
  
  !---- interface to csection.c
  interface
     subroutine csection (  ti,  stept, n, re_wr12,  im_wr12, re_wr23,  im_wr23,  G12,  G23, detun, n_fourier) bind ( c )
       use iso_c_binding
       integer ( c_int ),  VALUE :: n, n_fourier
       real ( c_double ),  VALUE :: ti, stept
       real ( c_double ), dimension(*)  :: re_wr12, im_wr12, re_wr23, im_wr23
       real ( c_double ), dimension(*)  :: G12, G23,detun
     end subroutine csection
  end interface
  !print*,"fortran wrap, array size: ",re_wr12(m)
  call csection (ti,stept,n,re_wr12, im_wr12, re_wr23, im_wr23,G12,G23,detun,n_fourier)
end subroutine f_csection
