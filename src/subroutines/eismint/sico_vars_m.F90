!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ v a r s _ m
!
!! EISMINT domain: Declarations of global variables for SICOPOLIS.
!!
!!##### Authors
!!
!! Ralf Greve
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
!> EISMINT domain: Declarations of global variables for SICOPOLIS.
!-------------------------------------------------------------------------------
module sico_vars_m

use sico_types_m

implicit none
save

#if (!defined(SURFACE_FORCING) || SURFACE_FORCING==1)

real(dp) :: temp_min
   !! Minimum surface temperature

real(dp) :: s_t
   !! Gradient of surface temperature change with horizontal distance

real(dp) :: b_max
   !! Maximum accumulation rate

real(dp) :: s_b
   !! Gradient of accumulation rate change with horizontal distance

real(dp) :: eld
   !! Equilibrium line distance

#elif (SURFACE_FORCING==2)

real(dp) :: temp_0
   !! Surface temperature at z=0

real(dp) :: gamma_t
   !! Lapse rate of the surface temperature

real(dp) :: s_0
   !! Maximum accumulation rate

real(dp) :: m_0
   !! Melting gradient (change of accumulation rate with elevation)

real(dp) :: ela
   !! Equilibrium line altitude

#endif

end module sico_vars_m
!
