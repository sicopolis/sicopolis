!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ v a r s _ m
!
!> @file
!!
!! Declarations of global variables for SICOPOLIS (for the ANT domain).
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
!> Declarations of global variables for SICOPOLIS (for the ANT domain).
!<------------------------------------------------------------------------------
module sico_vars_m

use sico_types_m

implicit none
save

#if (TSURFACE==6 && ACCSURFACE==6 && ABLSURFACE==6)
!> temp_maat_climatol(j,i): Surface-temperature (MAAT) climatology
   real(dp), dimension(0:JMAX,0:IMAX) :: temp_maat_climatol
!> smb_climatol(j,i): SMB climatology
   real(dp), dimension(0:JMAX,0:IMAX) :: smb_climatol
!> temp_maat_anom(j,i): Surface-temperature (MAAT) anomaly
   real(dp), dimension(0:JMAX,0:IMAX) :: temp_maat_anom
!> smb_anom(j,i): SMB anomaly
   real(dp), dimension(0:JMAX,0:IMAX) :: smb_anom
!> dtemp_maat_dz(j,i): Surface-temperature (MAAT) vertical gradient
   real(dp), dimension(0:JMAX,0:IMAX) :: dtemp_maat_dz
!> dsmb_dz(j,i): SMB vertical gradient
   real(dp), dimension(0:JMAX,0:IMAX) :: dsmb_dz
#endif

!> smb_corr_in(j,i): Prescribed SMB correction
   real(dp), dimension(0:JMAX,0:IMAX) :: smb_corr_in

!> flag_initmip_asmb: Flag for use of InitMIP SMB anomaly
logical :: flag_initmip_asmb

!> smb_anom_initmip(j,i): InitMIP SMB anomaly
   real(dp), dimension(0:JMAX,0:IMAX) :: smb_anom_initmip

!> flag_initmip_abmb: Flag for use of InitMIP sub-ice-shelf-melt anomaly
logical :: flag_initmip_abmb

!> ab_anom_initmip(j,i): InitMIP sub-ice-shelf-melt anomaly
   real(dp), dimension(0:JMAX,0:IMAX) :: ab_anom_initmip

!> flag_larmip: Flag for use of LARMIP sub-ice-shelf-melt anomaly
logical :: flag_larmip

!> n_larmip_region(j,i): LARMIP regions for ice shelf basal melting
   integer(i4b), dimension(0:JMAX,0:IMAX) :: n_larmip_region
!> ab_anom_larmip(n): LARMIP sub-ice-shelf-melt anomaly
   real(dp), dimension(0:7) :: ab_anom_larmip

!> n_bm_region(j,i): Regions for ice shelf basal melting
   integer(i4b), dimension(0:JMAX,0:IMAX) :: n_bm_region

#if (FLOATING_ICE_BASAL_MELTING==6)
!> z_tf_bm_present(n): Equidistant depth points of the
!>                     present-day thermal forcing data of the ocean
   real(dp), dimension(0:NZ_TF_BM) :: z_tf_bm_present
!> tf_bm_present(n,j,i): Present-day thermal forcing data of the ocean
   real(dp), dimension(0:NZ_TF_BM,0:JMAX,0:IMAX) :: tf_bm_present
!> z_tf_bm(n): Equidistant depth points of the
!>             thermal forcing data of the ocean
   real(dp), dimension(0:NZ_TF_BM) :: z_tf_bm
!> tf_bm(n,j,i): Thermal forcing data of the ocean
   real(dp), dimension(0:NZ_TF_BM,0:JMAX,0:IMAX) :: tf_bm
#endif

#if (ICE_SHELF_COLLAPSE_MASK==1)
!> H_ref_retreat(j,i): Reference ice thickness for the ice-shelf collapse mask
   real(dp), dimension(0:JMAX,0:IMAX) :: H_ref_retreat
!> r_mask_retreat(j,i): Ice-shelf collapse mask
   real(dp), dimension(0:JMAX,0:IMAX) :: r_mask_retreat
#endif

end module sico_vars_m
!
