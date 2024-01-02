!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ v a r s _ m
!
!> @file
!!
!! Declarations of global variables for SICOPOLIS (for the EISMINT domain).
!!
!! @section Copyright
!!
!! Copyright 2009-2024 Ralf Greve
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
!> Declarations of global variables for SICOPOLIS (for the EISMINT domain).
!<------------------------------------------------------------------------------
module sico_vars_m

use sico_types_m

implicit none
save

#if (!defined(SURFACE_FORCING) || SURFACE_FORCING==1)
!> temp_min: Minimum surface temperature
   real(dp) :: temp_min
!> s_t: Gradient of surface temperature change with horizontal distance
   real(dp) :: s_t
!> b_max: Maximum accumulation rate
   real(dp) :: b_max
!> s_b: Gradient of accumulation rate change with horizontal distance
   real(dp) :: s_b
!> eld: Equilibrium line distance
   real(dp) :: eld
#elif (SURFACE_FORCING==2)
!> temp_0: Surface temperature at z=0
   real(dp) :: temp_0
!> gamma_t: Lapse rate of the surface temperature
   real(dp) :: gamma_t
!> s_0: Maximum accumulation rate
   real(dp) :: s_0
!> m_0: Melting gradient (change of accumulation rate with elevation)
   real(dp) :: m_0
!> ela: Equilibrium line altitude
   real(dp) :: ela
#endif

end module sico_vars_m
!
