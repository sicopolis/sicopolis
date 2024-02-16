!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ v a r s _ m
!
!! GRL domain: Declarations of global variables for SICOPOLIS.
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
!> GRL domain: Declarations of global variables for SICOPOLIS.
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

#if (DISC==2)

integer(i4b) :: glann_time_min
   !! Minimum time of the data values for the
   !! global annual temperature anomaly

integer(i4b) :: glann_time_stp
   !! Time step of the data values for the
   !! global annual temperature anomaly

integer(i4b) :: glann_time_max
   !! Maximum time of the data values for the
   !! global annual temperature anomaly

integer(i4b) :: ndata_glann
   !! Number of data values for the global annual temperature anomaly

integer(i4b) , parameter :: ndata_glann_max = 262143
   !! Maximum allowed value of ndata_glann

real(dp), dimension(0:ndata_glann_max) :: dT_glann_CLIMBER
   !! Data values for the global annual temperature anomaly

#endif

#if (RETREAT_MASK==1)

real(dp), dimension(0:JMAX,0:IMAX) :: H_ref_retreat
   !! Reference ice thickness for the retreat mask

real(dp), dimension(0:JMAX,0:IMAX) :: r_mask_retreat
   !! Retreat mask

#endif

end module sico_vars_m
!
