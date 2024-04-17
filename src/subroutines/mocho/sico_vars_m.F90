!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ v a r s _ m
!
!! MOCHO domain: Declarations of global variables for SICOPOLIS.
!!
!!##### Authors
!!
!! Ralf Greve, Eduardo Flandez, Matthias Scheiter
!!
!!##### License
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
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with SICOPOLIS. If not, see <https://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------------------------------------------------------------------------------
!> MOCHO domain: Declarations of global variables for SICOPOLIS.
!-------------------------------------------------------------------------------
module sico_vars_m

use sico_types_m

implicit none
save

real(dp), dimension(0:JMAX,0:IMAX) :: smb_corr_in
   !! Prescribed SMB correction

#if (SURFACE_FORCING==2)

real(dp) :: temp_0
   !! Surface temperature at z=0

real(dp) :: gamma_t
   !! Lapse rate of the surface temperature

real(dp) :: s_0
   !! Maximum accumulation rate

real(dp) :: m_0
   !! Melting gradient (change of accumulation rate with elevation)

real(dp) :: ela
   !! Equilibrium line altitude (ELA)

real(dp) :: dela_dts
   !! Change of ELA with surface temperature

#elif (SURFACE_FORCING==3)

real(dp) :: temp_0
   !! Surface temperature at z=0

real(dp) :: gamma_t
   !! Lapse rate of the surface temperature

real(dp) :: s_0
   !! Maximum accumulation rate

real(dp) :: m_0
   !! Melting gradient (change of accumulation rate with elevation)

real(dp) :: ela
   !! Equilibrium line altitude (ELA)

real(dp) :: ela_amp
   !! Amplitude for aspect-dependent ELA variation

real(dp) :: dela_dts
   !! Change of ELA with surface temperature

real(dp) :: phi_0
   !! Offset angle for aspect-dependent ELA variation

real(dp) :: x_gip
   !! x-coord of summit for aspect-dependent ELA variation

real(dp) :: y_gip
   !! y-coord of summit for aspect-dependent ELA variation

#elif (SURFACE_FORCING==4)

real(dp) :: temp_0
   !! Surface temperature at z=0

real(dp) :: gamma_t
   !! Lapse rate of the surface temperature

real(dp) :: s_0
   !! Maximum accumulation rate

real(dp) :: m_0
   !! Melting gradient (change of accumulation rate with elevation)

real(dp) :: ela
   !! Equilibrium line altitude (ELA)

real(dp) :: ela_amp
   !! Amplitude for aspect-dependent ELA variation

real(dp) :: dela_dts
   !! Change of ELA with surface temperature

real(dp) :: phi_0
   !! Offset angle for aspect-dependent ELA variation

real(dp) :: x_gip
   !! x-coord of summit for aspect-dependent ELA variation

real(dp) :: y_gip
   !! y-coord of summit for aspect-dependent ELA variation

real(dp) :: m_1
   !! Melting gradient (change of accumulation rate with elevation) above elevation z_gc

real(dp) :: z_gc
   !! Gradient change elevation (above, m_1 rather than m_0 is used)

#elif (SURFACE_FORCING==5)

real(dp) :: temp_0
   !! Surface temperature at z=0

real(dp) :: gamma_t
   !! Lapse rate of the surface temperature

real(dp) :: s_0
   !! Maximum accumulation rate

real(dp) :: m_0
   !! Melting gradient (change of accumulation rate with elevation)

real(dp) :: ela
   !! Equilibrium line altitude (ELA)

real(dp) :: ela_amp
   !! Amplitude for aspect-dependent ELA variation

real(dp) :: dela_dts
   !! Change of ELA with surface temperature

real(dp) :: phi_0
   !! Offset angle for aspect-dependent ELA variation

real(dp) :: m_1
   !! Melting gradient (change of accumulation rate with elevation) above elevation z_gc

real(dp) :: z_gc
   !! Gradient change elevation (above, m_1 rather than m_0 is used)

#elif (SURFACE_FORCING==6)

real(dp) :: temp_0
   !! Surface temperature at z=0

real(dp) :: gamma_t
   !! Lapse rate of the surface temperature

real(dp) :: s_0
   !! Maximum accumulation rate

real(dp) :: m_0
   !! Melting gradient (change of accumulation rate with elevation)

real(dp) :: ela
   !! Equilibrium line altitude (ELA)

real(dp) :: ela_amp
   !! Amplitude for aspect-dependent ELA variation

real(dp) :: dela_dts
   !! Change of ELA with surface temperature

real(dp) :: phi_0
   !! Offset angle for aspect-dependent ELA variation

real(dp) :: tgt
   !! Topographic gradient threshold for higher ELA

#elif (SURFACE_FORCING==7)

real(dp) :: temp_0
   !! Surface temperature at z=0

real(dp) :: gamma_t
   !! Lapse rate of the surface temperature

real(dp) :: s_0
   !! Maximum accumulation rate

real(dp) :: m_0
   !! Melting gradient (change of accumulation rate with elevation)

real(dp) :: ela
   !! Equilibrium line altitude (ELA)

real(dp) :: ela_amp
   !! Amplitude for aspect-dependent ELA variation

real(dp) :: dela_dts
   !! Change of ELA with surface temperature

real(dp) :: phi_0
   !! Offset angle for aspect-dependent ELA variation

#elif (SURFACE_FORCING==8)

real(dp) :: temp_0
   !! Surface temperature at z=0

real(dp) :: gamma_t
   !! Lapse rate of the surface temperature

real(dp) :: s_0
   !! Maximum accumulation rate

real(dp) :: m_0
   !! Melting gradient (change of accumulation rate with elevation)

real(dp) :: ela
   !! Equilibrium line altitude (ELA)

real(dp) :: ela_amp
   !! Amplitude for aspect-dependent ELA variation

real(dp) :: dela_dts
   !! Change of ELA with surface temperature

real(dp) :: phi_0
   !! Offset angle for aspect-dependent ELA variation

real(dp) :: x_gip
   !! x-coord of summit for aspect-dependent ELA variation

real(dp) :: y_gip
   !! y-coord of summit for aspect-dependent ELA variation

real(dp) :: x_gip2
   !! x-coord of summit for aspect-dependent ELA variation (Choshuenco)

real(dp) :: y_gip2
   !! y-coord of summit for aspect-dependent ELA variation (Choshuenco)

real(dp) :: ela_amp2
   !! Amplitude for aspect-dependent ELA variation (Choshuenco)

#endif

end module sico_vars_m
!
