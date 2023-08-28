!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ v a r s _ m
!
!> @file
!!
!! Declarations of global variables for SICOPOLIS (for the GRL domain).
!!
!! @section Copyright
!!
!! Copyright 2009-2023 Ralf Greve
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
!> Declarations of global variables for SICOPOLIS (for the GRL domain).
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

#if (defined(INITMIP_SMB_ANOM_FILE))
!> smb_anom_initmip(j,i): InitMIP anomaly of the accumulation-ablation function
!>                        at the ice surface (surface mass balance)
   real(dp), dimension(0:JMAX,0:IMAX) :: smb_anom_initmip
#endif

#if (DISC==2)
!> glann_time_min: Minimum time of the data values for the
!>                 global annual temperature anomaly
   integer(i4b) :: glann_time_min
!> glann_time_stp: Time step of the data values for the
!>                 global annual temperature anomaly
   integer(i4b) :: glann_time_stp
!> glann_time_max: Maximum time of the data values for the
!>                 global annual temperature anomaly
   integer(i4b) :: glann_time_max
!> ndata_glann: Number of data values for the global annual temperature anomaly
   integer(i4b) :: ndata_glann
!> ndata_glann_max: Maximum allowed value of ndata_glann
   integer(i4b) , parameter :: ndata_glann_max = 262143
!> dT_glann_CLIMBER(n): Data values for the global annual temperature anomaly
   real(dp), dimension(0:ndata_glann_max) :: dT_glann_CLIMBER
#endif

#if (RETREAT_MASK==1)
!> H_ref_retreat(j,i): Reference ice thickness for the retreat mask
   real(dp), dimension(0:JMAX,0:IMAX) :: H_ref_retreat
!> r_mask_retreat(j,i): Retreat mask
   real(dp), dimension(0:JMAX,0:IMAX) :: r_mask_retreat
#endif

end module sico_vars_m
!
