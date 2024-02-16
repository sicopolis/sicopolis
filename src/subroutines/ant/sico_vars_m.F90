!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ v a r s _ m
!
!! ANT domain: Declarations of global variables for SICOPOLIS.
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
!> ANT domain: Declarations of global variables for SICOPOLIS.
!-------------------------------------------------------------------------------
module sico_vars_m

use sico_types_m

implicit none
save

#if (TSURFACE==6 && ACCSURFACE==6 && ABLSURFACE==6)

real(dp), dimension(0:JMAX,0:IMAX) :: temp_maat_climatol
   !! Surface-temperature (MAAT) climatology

real(dp), dimension(0:JMAX,0:IMAX) :: smb_climatol
   !! SMB climatology

real(dp), dimension(0:JMAX,0:IMAX) :: temp_maat_anom
   !! Surface-temperature (MAAT) anomaly

real(dp), dimension(0:JMAX,0:IMAX) :: smb_anom
   !! SMB anomaly

real(dp), dimension(0:JMAX,0:IMAX) :: dtemp_maat_dz
   !! Surface-temperature (MAAT) vertical gradient

real(dp), dimension(0:JMAX,0:IMAX) :: dsmb_dz
   !! SMB vertical gradient

#endif

real(dp), dimension(0:JMAX,0:IMAX) :: smb_corr_in
   !! Prescribed SMB correction

logical :: flag_initmip_asmb
   !! Flag for use of InitMIP SMB anomaly

real(dp), dimension(0:JMAX,0:IMAX) :: smb_anom_initmip
   !! InitMIP SMB anomaly

logical :: flag_initmip_abmb
   !! Flag for use of InitMIP sub-ice-shelf-melt anomaly

real(dp), dimension(0:JMAX,0:IMAX) :: ab_anom_initmip
   !! InitMIP sub-ice-shelf-melt anomaly

logical :: flag_larmip
   !! Flag for use of LARMIP sub-ice-shelf-melt anomaly

integer(i4b), dimension(0:JMAX,0:IMAX) :: n_larmip_region
   !! LARMIP regions for ice shelf basal melting

real(dp), dimension(0:7) :: ab_anom_larmip
   !! LARMIP sub-ice-shelf-melt anomaly

integer(i4b), dimension(0:JMAX,0:IMAX) :: n_bm_region
   !! Regions for ice shelf basal melting

#if (FLOATING_ICE_BASAL_MELTING==6)

real(dp), dimension(0:NZ_TF_BM) :: z_tf_bm_present
   !! Equidistant depth points of the
   !! present-day thermal forcing data of the ocean

real(dp), dimension(0:NZ_TF_BM,0:JMAX,0:IMAX) :: tf_bm_present
   !! Present-day thermal forcing data of the ocean

real(dp), dimension(0:NZ_TF_BM) :: z_tf_bm
   !! Equidistant depth points of the
   !! thermal forcing data of the ocean

real(dp), dimension(0:NZ_TF_BM,0:JMAX,0:IMAX) :: tf_bm
   !! Thermal forcing data of the ocean

#endif

#if (ICE_SHELF_COLLAPSE_MASK==1)

real(dp), dimension(0:JMAX,0:IMAX) :: H_ref_retreat
   !! Reference ice thickness for the ice-shelf collapse mask

real(dp), dimension(0:JMAX,0:IMAX) :: r_mask_retreat
   !! Ice-shelf collapse mask

#endif

end module sico_vars_m
!
