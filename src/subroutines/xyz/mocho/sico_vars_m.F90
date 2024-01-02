!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ v a r s _ m
!
!> @file
!!
!! Declarations of global variables for SICOPOLIS (for the MOCHO domain).
!!
!! @section Copyright
!!
!! Copyright 2009-2024 Ralf Greve, Eduardo Flandez, Matthias Scheiter
!!
!! @section License
!!
!! This file is part of SICOPOLIS.
!!
!! SICOPOLIS is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! SICOPOLIS is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with SICOPOLIS.  If not, see <http://www.gnu.org/licenses/>.
!<
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------------------------------------------------------------------------------
!> Declarations of global variables for SICOPOLIS (for the MOCHO domain).
!<------------------------------------------------------------------------------
module sico_vars_m

use sico_types_m

implicit none
save

!> smb_corr_in(j,i): Prescribed SMB correction
   real(dp), dimension(0:JMAX,0:IMAX) :: smb_corr_in

#if (SURFACE_FORCING==2)
!> temp_0: Surface temperature at z=0
   real(dp) :: temp_0
!> gamma_t: Lapse rate of the surface temperature
   real(dp) :: gamma_t
!> s_0: Maximum accumulation rate
   real(dp) :: s_0
!> m_0: Melting gradient (change of accumulation rate with elevation)
   real(dp) :: m_0
!> ela: Equilibrium line altitude (ELA)
   real(dp) :: ela
!> dela_dts: Change of ELA with surface temperature
   real(dp) :: dela_dts
#elif (SURFACE_FORCING==3)
!> temp_0: Surface temperature at z=0
   real(dp) :: temp_0
!> gamma_t: Lapse rate of the surface temperature
   real(dp) :: gamma_t
!> s_0: Maximum accumulation rate
   real(dp) :: s_0
!> m_0: Melting gradient (change of accumulation rate with elevation)
   real(dp) :: m_0
!> ela: Equilibrium line altitude (ELA)
   real(dp) :: ela
!> ela_amp: Amplitude for aspect-dependent ELA variation
   real(dp) :: ela_amp
!> dela_dts: Change of ELA with surface temperature
   real(dp) :: dela_dts
!> phi_0: Offset angle for aspect-dependent ELA variation
   real(dp) :: phi_0
!> x_gip: x-coord of summit for aspect-dependent ELA variation
   real(dp) :: x_gip
!> y_gip: y-coord of summit for aspect-dependent ELA variation
   real(dp) :: y_gip
#elif (SURFACE_FORCING==4)
!> temp_0: Surface temperature at z=0
   real(dp) :: temp_0
!> gamma_t: Lapse rate of the surface temperature
   real(dp) :: gamma_t
!> s_0: Maximum accumulation rate
   real(dp) :: s_0
!> m_0: Melting gradient (change of accumulation rate with elevation)
   real(dp) :: m_0
!> ela: Equilibrium line altitude (ELA)
   real(dp) :: ela
!> ela_amp: Amplitude for aspect-dependent ELA variation
   real(dp) :: ela_amp
!> dela_dts: Change of ELA with surface temperature
   real(dp) :: dela_dts
!> phi_0: Offset angle for aspect-dependent ELA variation
   real(dp) :: phi_0
!> x_gip: x-coord of summit for aspect-dependent ELA variation
   real(dp) :: x_gip
!> y_gip: y-coord of summit for aspect-dependent ELA variation
   real(dp) :: y_gip
!> m_1: Melting gradient (change of accumulation rate with elevation) above elevation z_gc
   real(dp) :: m_1
!> z_gc: Gradient change elevation (above, m_1 rather than m_0 is used)
   real(dp) :: z_gc
#elif (SURFACE_FORCING==5)
!> temp_0: Surface temperature at z=0
   real(dp) :: temp_0
!> gamma_t: Lapse rate of the surface temperature
   real(dp) :: gamma_t
!> s_0: Maximum accumulation rate
   real(dp) :: s_0
!> m_0: Melting gradient (change of accumulation rate with elevation)
   real(dp) :: m_0
!> ela: Equilibrium line altitude (ELA)
   real(dp) :: ela
!> ela_amp: Amplitude for aspect-dependent ELA variation
   real(dp) :: ela_amp
!> dela_dts: Change of ELA with surface temperature
   real(dp) :: dela_dts
!> phi_0: Offset angle for aspect-dependent ELA variation
   real(dp) :: phi_0
!> m_1: Melting gradient (change of accumulation rate with elevation) above elevation z_gc
   real(dp) :: m_1
!> z_gc: Gradient change elevation (above, m_1 rather than m_0 is used)
   real(dp) :: z_gc
#elif (SURFACE_FORCING==6)
!> temp_0: Surface temperature at z=0
   real(dp) :: temp_0
!> gamma_t: Lapse rate of the surface temperature
   real(dp) :: gamma_t
!> s_0: Maximum accumulation rate
   real(dp) :: s_0
!> m_0: Melting gradient (change of accumulation rate with elevation)
   real(dp) :: m_0
!> ela: Equilibrium line altitude (ELA)
   real(dp) :: ela
!> ela_amp: Amplitude for aspect-dependent ELA variation
   real(dp) :: ela_amp
!> dela_dts: Change of ELA with surface temperature
   real(dp) :: dela_dts
!> phi_0: Offset angle for aspect-dependent ELA variation
   real(dp) :: phi_0
!> tgt: Topographic gradient threshold for higher ELA
   real(dp) :: tgt
#elif (SURFACE_FORCING==7)
!> temp_0: Surface temperature at z=0
   real(dp) :: temp_0
!> gamma_t: Lapse rate of the surface temperature
   real(dp) :: gamma_t
!> s_0: Maximum accumulation rate
   real(dp) :: s_0
!> m_0: Melting gradient (change of accumulation rate with elevation)
   real(dp) :: m_0
!> ela: Equilibrium line altitude (ELA)
   real(dp) :: ela
!> ela_amp: Amplitude for aspect-dependent ELA variation
   real(dp) :: ela_amp
!> dela_dts: Change of ELA with surface temperature
   real(dp) :: dela_dts
!> phi_0: Offset angle for aspect-dependent ELA variation
   real(dp) :: phi_0
#elif (SURFACE_FORCING==8)
!> temp_0: Surface temperature at z=0
   real(dp) :: temp_0
!> gamma_t: Lapse rate of the surface temperature
   real(dp) :: gamma_t
!> s_0: Maximum accumulation rate
   real(dp) :: s_0
!> m_0: Melting gradient (change of accumulation rate with elevation)
   real(dp) :: m_0
!> ela: Equilibrium line altitude (ELA)
   real(dp) :: ela
!> ela_amp: Amplitude for aspect-dependent ELA variation
   real(dp) :: ela_amp
!> dela_dts: Change of ELA with surface temperature
   real(dp) :: dela_dts
!> phi_0: Offset angle for aspect-dependent ELA variation
   real(dp) :: phi_0
!> x_gip: x-coord of summit for aspect-dependent ELA variation
   real(dp) :: x_gip
!> y_gip: y-coord of summit for aspect-dependent ELA variation
   real(dp) :: y_gip
!> x_gip: x-coord of summit for aspect-dependent ELA variation (Choshuenco)
   real(dp) :: x_gip2
!> y_gip: y-coord of summit for aspect-dependent ELA variation (Choshuenco)
   real(dp) :: y_gip2
!> ela_amp: Amplitude for aspect-dependent ELA variation (Choshuenco)
   real(dp) :: ela_amp2
#endif

end module sico_vars_m
!
