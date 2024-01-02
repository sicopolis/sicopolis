!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ v a r s _ m
!
!> @file
!!
!! Declarations of global variables for SICOPOLIS (for the HEINO domain).
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
!> Declarations of global variables for SICOPOLIS (for the HEINO domain).
!<------------------------------------------------------------------------------
module sico_vars_m

use sico_types_m

implicit none
save

!> temp_min: Minimum surface temperature
   real(dp) :: temp_min
!> s_t: Gradient of surface temperature change with horizontal distance
   real(dp) :: s_t
!> x_hat: Coordinate xi (= x) of the centre of the model domain
   real(dp) :: x_hat
!> y_hat: Coordinate eta (= y) of the centre of the model domain
   real(dp) :: y_hat
!> rad: Radius of the model domain
   real(dp) :: rad
!> b_min: Minimum accumulation rate
   real(dp) :: b_min
!> b_max: Maximum accumulation rate
   real(dp) :: b_max

end module sico_vars_m
!
