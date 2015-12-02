subroutine coupled3(dim,nd,n,np,xp,xi,sh,shm,U0_in,V0_in,VPOT_in,VMINB,VMINC,gamma,&
     tdipol,omega,e0,tp,td,t0,sni,kl,ti,tf,dt,nshot,mxdct)
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

  implicit none
  integer, parameter :: dp=kind(1.0d00)
  !------external arguments (from eSPec)
  integer, intent(in) :: n,nd,nshot,mxdct
  integer, dimension(3), intent(in) :: np,kl
  character (len=*) :: dim
  real(kind=dp) VMINB,VMINC,t,dt, ti, tf,EF,xp,tout,start_t,stop_t
  real(kind=dp), dimension(3) :: omega,e0,tp,td,t0,sni,tdipol,shm,sh,xi,gamma
  real(kind=dp), dimension(mxdct), intent(in) :: U0_in,V0_in
  real(kind=dp), dimension(mxdct), intent(in) :: VPOT_in
  real(kind=dp), dimension(3*n) :: VPOT
  real(kind=dp), dimension(6*n) :: Y, F, HY

  !-----local variables
  integer i,j,k,ntime,ierr,dim_int
  integer(kind=C_LONG), dimension(1) :: IPAR=0
  integer(kind=C_LONG), dimension(21) :: IOUT
  integer(kind=C_LONG) :: j_long

  real(kind=dp), dimension(6) :: ROUT
  real(kind=dp), dimension(1) :: RPAR=0.0d0
  real(kind=dp), dimension(3) :: detun,G,CSD,SND
  real(kind=dp), dimension(8) :: WKG
  real(kind=dp) :: E1,E2,var,norm
  character*20 NEWNAM1,NEWNAM2,NEWNAM3
  character*4 chnum
  real(kind=dp), parameter :: fs2au=41.3411D+0, ev2au=27.2114D+0

  !----------------------------------------------

  VPOT(:)=VPOT_in(1:3*n)

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

  !------detuning from resonances
  detun(1) = omega(1) - VMINB
  detun(2) = omega(2) - (VMINB - VMINC)

  !------- convert to au
  detun(1) = detun(1)/ev2au
  detun(2) = detun(2)/ev2au
  gamma(1) = gamma(1)/ev2au
  gamma(2) = gamma(2)/ev2au
  gamma(3) = gamma(3)/ev2au



  call initialize_field(e0,t0,td,tp,detun,sni,kl,tdipol,gamma)


  !--------time parameters
  print*,'ti ',ti,'tf ',tf,'step',dt
  print*,'number of points',ntime
  write(*,*)

  print*,'initialize n vector'
  call FNVINITS(1,psize,ierr);call chkierr(ierr)
  print*,'allocating'
  call FCVMALLOC(ti, Y, 2, & ! 2 = BDF, 1 = Adams
       2, & ! 1 = functional iteration, 2 = Newton iteration
       1, & ! 1 = scalar absolute tolerance (?), 2 = array absolute tolerance (?), 3 = user-defined
       1.0d-8, &! Relative tolerance
       1.0d-8, &! Absolute tolerance
       IOUT, ROUT, & ! Integer and real optional outputs
       IPAR, RPAR, & ! Integer and real optional parameters
       ierr);call chkierr(ierr)

  print*,'set pre-conditioner'
  call FCVDIAG(ierr);call chkierr(ierr)

  print*,' set the solver'
  call FCVSPGMR(0,& ! 0 = for no preconditioner (problematic!)
       1, & ! Modified GS orthogonalization
       5000, & ! Maximum size of Krylov subspace
       1.0d-8, & ! Tolerance
       ierr);call chkierr(ierr)

  !---------time loop
  do i=2,ntime
     t = ti + (i - 1)*dt
     t = t * fs2au

     !---------Andrey's propagation routines   
     call cpu_time(start_t)
     call FCVODE(t, tout, Y, &
          1,  & ! 1 = overshoot and interpolate to t,  2 = one-step mode
          ierr);call chkierr(ierr)
     call cpu_time(stop_t)
     write(*,*) 'iteration time = ',stop_t-start_t
     !--------------------------------------

     !-----------printing wavepackets
    
     WRITE(CHNUM,'(I4)')i
     t = t/fs2au
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

     CALL PRPT2(NEWNAM1,8,ND,t,NP,XP,XI,SH,Y(1:n),Y(n+1:2*n))
     CALL PRPT2(NEWNAM2,8,ND,T,NP,XP,XI,SH,Y(2*n+1:3*n),Y(3*n+1:4*n))
     CALL PRPT2(NEWNAM3,8,ND,T,NP,XP,XI,SH,Y(4*n+1:5*n),Y(5*n+1:6*n))

  end do
  !--------------------


  !------------------------------------------------------
  write(*,*)
  write(*,*)'exiting coupled3'
  write(*,*)
  !------------------------------------------------------
end subroutine coupled3


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
