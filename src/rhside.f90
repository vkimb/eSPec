!----- module for setting right hand side of wavepacket coupled propagation
!
module rhside
  use iso_c_binding
  implicit none
  integer, parameter :: dp=kind(1.0d00)
  integer :: n,nd,nshot,dim
  integer(kind=C_LONG) psize
  integer, dimension(3) :: np,kl
  real(kind=dp) VMINB,VMINC,xp
  real(kind=dp), allocatable, dimension(:) :: U0,V0,VPOT
  real(kind=dp), dimension(3) :: omega,e0,tp,td,t0,sni,tdipol,gamma,xi,sh,shm,detun

  public :: initialize_initguess, initialize_field, compute_rhside,psize,chkierr,chkierr2
  private
contains

  subroutine initialize_initguess(dim_in,nd_in,n_in,np_in,xp_in,xi_in,sh_in,shm_in,VPOT_in,VMINB_in,&
       VMINC_in,ti_in)
    integer, intent(in) :: n_in,nd_in,dim_in
    integer,  intent(in), dimension(3) :: np_in
    real(kind=dp), intent(in) :: VMINB_in,VMINC_in,xp_in,ti_in
    real(kind=dp), intent(in),  dimension(3) :: xi_in,sh_in,shm_in
    real(kind=dp), intent(in),  dimension(3*n_in) :: VPOT_in
    integer :: ierr
   
    allocate(U0(n_in),stat=ierr);call chkierr(ierr)
    allocate(V0(n_in),stat=ierr);call chkierr(ierr)
    allocate(VPOT(3*n_in),stat=ierr);call chkierr(ierr)

    n = n_in
    psize = 6*n
    np=np_in
    nd = nd_in
    VMINB=VMINB_in
    VMINC=VMINC_in
    xp = xp_in; xi=xi_in; sh=sh_in;
    shm=shm_in; VPOT=VPOT_in; dim = dim_in


  end subroutine initialize_initguess

!---------------------------------------

  subroutine initialize_field(e0_in,t0_in,td_in,tp_in,detun_in,sni_in,kl_in,&
       tdipol_in,gamma_in)
    integer, intent(in), dimension(3) :: kl_in
    real(kind=dp), intent(in), dimension(3) :: e0_in,t0_in,td_in,tp_in,detun_in,sni_in,tdipol_in,gamma_in


    e0=e0_in; t0=t0_in; td=td_in; tp=tp_in; detun=detun_in; 
    sni=sni_in; tdipol=tdipol_in;gamma=gamma_in;kl=kl_in

  end subroutine initialize_field

!----------------------------------------

  subroutine compute_rhside(t,Y,F)
    real(kind=dp), intent(in) :: t
    real(kind=dp), intent(in), dimension(6*n) :: Y
    real(kind=dp), intent(out), dimension(6*n) :: F

    integer :: j
    real(kind=dp) :: E1,E2,var,EF
    real(kind=dp), dimension(8) :: wkg
    real(kind=dp), dimension(2) :: G,CSD,SND
    real(kind=dp), dimension(6*n) :: HY
    character (len=*), parameter ::  pul=".ENVG "
    character (len=6) :: dimc
    real(kind=dp), parameter :: fs2au=41.3411D+0, ev2au=27.2114D+0

    !open(unit=43,file='pulses_deb.dat') 

    if(dim.EQ.0)then
       dimc = '.1D'
    else if(dim.EQ.1)then
       dimc = '.2D'
    else if(dim.EQ.2)then
       dimc = '.2DCT'
    else
       print *, "wrong dimension!!"
       stop
    end if

    ! write(*,*)
    ! write(*,*) 'debug'
    ! write(*,*)t,t/fs2au
    ! write(*,*)
    ! write(*,*)e0(1),t0(1),td(1),tp(1),detun(1),sni(1),kl(1)
    ! write(*,*)
    ! write(*,*)e0(2),t0(2),td(2),tp(2),detun(2),sni(2),kl(2)
    ! read(*,*)

    E1 = EF(pul,e0(1),t0(1),td(1),tp(1),detun(1),sni(1),kl(1),t/fs2au)
    E2 = EF(pul,e0(2),t0(2),td(2),tp(2),detun(2),sni(2),kl(2),t/fs2au)   
    
    !write(43,*)t/fs2au,E1,E2

    G(1) = -0.5D+0 * (E1 * tdipol(1)) ! G12
    G(2) = -0.5D+0 * (E2 * tdipol(2)) ! G23 = G32*   

    !-------- temporal oscilation within RWA
    CSD(1) = dcos(detun(1) * t)
    SND(1) = dsin(detun(1) * t)
    CSD(2) = dcos(detun(2) * t)
    SND(2) = dsin(detun(2) * t)
    
    !----------------------------
    !---------Set right hand side
    !---------apply hamiltonian
    call AU(dimc,n,np,shm,VPOT,Y(1),HY(1),var)
    call AU(dimc,n,np,shm,VPOT,Y(n + 1),HY(n + 1),var)

    call AU(dimc,n,np,shm,VPOT(n+1),Y(2*n + 1),HY(2*n + 1),var)
    call AU(dimc,n,np,shm,VPOT(n+1),Y(3*n + 1),HY(3*n + 1),var)

    call AU(dimc,n,np,shm,VPOT(2*n + 1),Y(4*n + 1),HY(4*n + 1),var)
    call AU(dimc,n,np,shm,VPOT(2*n + 1),Y(5*n + 1),HY(5*n + 1),var)
    !---------------------------
    
    do j=1,n
       !-----------x_1
       wkg(1) = G(1)*(Y(2*n+j)*CSD(1) - Y(3*n+j)*SND(1)) ! exp(i detun1) x_2 = wkg1 + i wkg2
       wkg(2) = G(1)*(Y(3*n+j)*CSD(1) + Y(2*n+j)*SND(1))

       !F(j) =  HY(n+j) + wkg(2) -gamma(1)*Y(n+j)
       !F(n+j) = -(HY(j) + wkg(1) -gamma(1)*Y(j))

       F(j) =  HY(n+j) + wkg(2) -gamma(1)*Y(j)
       F(n+j) = -(HY(j) + wkg(1) +gamma(1)*Y(n+j))


       !-----------x_2            
       wkg(3) = G(1)*(Y(j)*CSD(1) + Y(n+j)*SND(1)) ! exp(-i detun1) x_1 = wkg3 + i wkg4
       wkg(5) = G(2)*(Y(4*n+j)*CSD(2) + Y(5*n+j)*SND(2)) ! exp(-i detun2) x_3 = wkg5 + i wkg6

       wkg(4) = G(1)*(Y(n+j)*CSD(1) - Y(j)*SND(1))
       wkg(6) = G(2)*(Y(5*n+j)*CSD(2) - Y(4*n+j)*SND(2))

       !F(2*n+j) = HY(3*n+j) + wkg(4) + wkg(6) -gamma(2)*Y(3*n+j)
       !F(3*n+j) = -(HY(2*n+j) + wkg(3) + wkg(5) -gamma(2)*Y(2*n+j))

       F(2*n+j) = HY(3*n+j) + wkg(4) + wkg(6) -gamma(2)*Y(2*n+j)
       F(3*n+j) = -(HY(2*n+j) + wkg(3) + wkg(5) +gamma(2)*Y(3*n+j))

       !-----------x_3            
       wkg(7) = G(2)*(Y(2*n+j)*CSD(2) - Y(3*n+j)*SND(2)) ! exp(i detun2) x_2 = wkg7 + i wkg8            
       wkg(8) = G(2)*(Y(3*n+j)*CSD(2) + Y(2*n+j)*SND(2))

       !F(4*n+j) = HY(5*n+j) + wkg(8) -gamma(3)*Y(5*n+j)
       !F(5*n+j) = -(HY(4*n+j) + wkg(7) -gamma(3)*Y(4*n+j))

       F(4*n+j) = HY(5*n+j) + wkg(8) -gamma(3)*Y(4*n+j)
       F(5*n+j) = -(HY(4*n+j) + wkg(7) +gamma(3)*Y(5*n+j))

    end do
    !----------------------------


  end subroutine compute_rhside

!----------------------------------------

  subroutine chkierr(ierr)
    integer, intent(in) :: ierr

    if(ierr.NE.0)then
       write(*,*) 'error in allocation'
       stop
    end if

  end subroutine chkierr

!-----------------------------

  subroutine chkierr2(ierr)
    integer, intent(in) :: ierr

    if(ierr.LT.0)then
       write(*,*) 'error in iteration'
       stop
    end if

  end subroutine chkierr2

end module rhside
