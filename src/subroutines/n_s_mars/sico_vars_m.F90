!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ v a r s _ m
!
!! NMARS/SMARS domains: Declarations of global variables for SICOPOLIS.
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
!> NMARS/SMARS domains: Declarations of global variables for SICOPOLIS.
!-------------------------------------------------------------------------------
module sico_vars_m

use sico_types_m

implicit none
save

integer(i4b) :: insol_time_min
   !! Minimum time of the data values for the insolation etc.

integer(i4b) :: insol_time_stp
   !! Time step of the data values for the insolation etc.

integer(i4b) :: insol_time_max
   !! Maximum time of the data values for the insolation etc.

real(dp), dimension(0:100000) :: insol_ma_90
   !! Data values for the mean-annual north- or south-polar insolation

real(dp), dimension(0:100000) :: obl_data
   !! Data values for the obliquity

real(dp), dimension(0:100000) :: ecc_data
   !! Data values for the eccentricity

real(dp), dimension(0:100000) :: ave_data
   !! Data values for the anomaly of vernal equinox
   !! (= 360 deg - Ls of perihelion )

real(dp), dimension(0:100000) :: cp_data
   !! Data values for Laskar's climate parameter
   !!  = eccentricity
   !!    *sin(Laskar's longitude of perihelion from moving equinox),
   !!      ( where Laskar's longitude of perihelion from moving equinox
   !!              = Ls of perihelion - 180 deg )

integer(i4b), dimension(0:JMAX,0:IMAX) :: mask_chasm
   !! Chasma mask.
   !!  0: grounded ice,
   !!  1: ice-free land (normal area),
   !!  7: chasma area

real(dp) :: time_chasm_init
   !! Initial time for active chasma area

real(dp) :: time_chasm_end
   !! Final time for active chasma area

real(dp), dimension(0:JMAX,0:IMAX) :: q_geo_normal
   !! Geothermal heat flux for normal (non-chasma) areas

real(dp) :: RHO_I
   !! Density of ice

real(dp) :: RHO_C
   !! Density of crustal material (dust)

real(dp) :: KAPPA_C
   !! Heat conductivity of crustal material (dust)

real(dp) :: C_C
   !! Specific heat of crustal material (dust)

real(dp) :: rho_inv
   !! Inverse of the density of ice-dust mixture

end module sico_vars_m
!
