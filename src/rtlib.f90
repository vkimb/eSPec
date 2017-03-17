module rtlib

  implicit none 

  integer, parameter :: dp=kind(1.0d0)
  type rtlib_type
     real(kind=dp) :: chirp_1p, chirp_2p, omega0_1p, omega0_2p
     real(kind=dp) :: resonant_energy_potAB, resonant_energy_potBC
     real(kind=dp) :: phase_1p, phase_2p ! Initial phases in radians, from input
  end type rtlib_type

  type(rtlib_type) :: rt
  
  public :: rt, compute_chirped_omega, ha_to_ev, ev_to_ha
  private
  
contains
  pure elemental real(kind=dp) function ha_to_ev(x)
    implicit none
    real(kind=dp), intent(in) :: x
    ha_to_ev = x*27.2113850560
  end function ha_to_ev

  pure elemental real(kind=dp) function ev_to_ha(x)
    implicit none
    real(kind=dp), intent(in) :: x
    ev_to_ha = x/27.2113850560
  end function ev_to_ha
  
  subroutine compute_chirped_omega(omega,T,T0,pulse_idx)
    real(kind=dp), intent(inout) :: omega
    real(kind=dp), intent(in)    :: T, T0
    integer, intent(in)          :: pulse_idx
    
    select case (pulse_idx)
    case(1)
       omega = rt%omega0_1p + rt%chirp_1p*(T-T0)**2
    case(2)
       omega = rt%omega0_2p + rt%chirp_2p*(T-T0)**2
    case default
       write(*,*) " Unkown pulse number ", pulse_idx
    end select
  end subroutine compute_chirped_omega


end module rtlib
