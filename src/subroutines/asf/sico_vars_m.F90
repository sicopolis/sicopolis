!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ v a r s _ m
!
!! ASF domain: Declarations of global variables for SICOPOLIS.
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
!> ASF domain: Declarations of global variables for SICOPOLIS.
!-------------------------------------------------------------------------------
module sico_vars_m

use sico_types_m

implicit none
save

integer(i4b) :: n_surf
   !! Number of surface points for which time-series data are written 

integer(i4b), parameter :: n_surf_max = 256
   !! Maximum allowed value of n_surf

real(dp), dimension(n_surf_max) :: lambda_surf
   !! Geographical longitude of the prescribed surface points

real(dp), dimension(n_surf_max) :: phi_surf
   !! Geographical latitude of the prescribed surface points

real(dp), dimension(n_surf_max) :: x_surf
   !! Coordinate xi (= x) of the prescribed surface points

real(dp), dimension(n_surf_max) :: y_surf
   !! Coordinate eta (= y) of the prescribed surface points

end module sico_vars_m
!
