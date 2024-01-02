!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  b o u n d a r y _ m
!
!> @file
!!
!! Computation of the surface temperature (must be less than 0 degC)
!! and of the accumulation-ablation function.
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
!> Computation of the surface temperature (must be less than 0 degC)
!! and of the accumulation-ablation function.
!<------------------------------------------------------------------------------
module boundary_m

  use sico_types_m
  use sico_variables_m
  use sico_vars_m
  use error_m
#if defined(ALLOW_TAPENADE)
  use globals
#endif
  implicit none

  public

contains

!-------------------------------------------------------------------------------
!> Main routine of boundary_m:
!! Computation of the surface temperature (must be less than 0 degC)
!! and of the accumulation-ablation function.
!<------------------------------------------------------------------------------
subroutine boundary(time, dtime, dxi, deta, &
                    delta_ts, glac_index, z_mar)

#if defined(ALLOW_TAPENADE) /* Tapenade */
  use ctrl_m, only: myceiling, myfloor
  use sico_maths_m, only: my_erfc
#endif /* Tapenade */

  use netcdf
  use nc_check_m

#if ((MARGIN==2) \
      && (MARINE_ICE_FORMATION==2) \
      && (MARINE_ICE_CALVING==9))
  use calving_m
#endif

#if (DISC > 0)
  use discharge_workers_m, only: discharge, dis_perp, dT_glann
#endif

  use mask_update_sea_level_m
  use pdd_m

implicit none

real(dp), intent(in) :: time, dtime, dxi, deta

real(dp), intent(out)   :: delta_ts, glac_index, z_mar

! Further return variables
! (defined as global variables in module sico_variables_m):
!
!    accum(j,i), runoff(j,i), as_perp(j,i), calving(j,i), temp_s(j,i)

integer(i4b) :: i, j, n
integer(i4b) :: i_gr, i_kl
integer(i4b) :: n_year_CE
integer(i4b) :: ios
real(dp), dimension(0:JMAX,0:IMAX) :: z_sl_old
real(dp) :: z_sl_old_mean
real(dp) :: z_sl_min, t1, t2, t3, t4, t5, t6
real(dp) :: time_in_years
real(dp) :: rho_inv
real(dp) :: time_gr, time_kl
real(dp) :: z_sle_present, z_sle_help

real(dp), dimension(0:JMAX,0:IMAX,0:12) :: precip
real(dp), dimension(0:JMAX,0:IMAX,12)   :: temp_mm
real(dp), dimension(12)                 :: temp_mm_help
real(dp), dimension(0:JMAX,0:IMAX)      :: temp_ma, temp_ampl
real(dp), dimension(0:JMAX,0:IMAX)      :: accum_prescribed, &
                                           runoff_prescribed

real(dp) :: temp_jja_help
real(dp) :: theta_ma, theta_ma_offset, c_ma, kappa_ma, gamma_ma, &
            theta_mj, theta_mj_offset, c_mj, kappa_mj, gamma_mj
real(dp) :: sine_factor
real(dp) :: gamma_p, zs_thresh, &
            temp_rain, temp_snow, &
            inv_delta_temp_rain_snow, coeff(0:5), inv_sqrt2_s_stat, &
            precip_fact, frac_solid
real(dp) :: s_stat, beta1, beta2, Pmax, mu, lambda_lti, temp_lti
real(dp) :: r_aux
character(len=256) :: ch_aux
logical, dimension(0:JMAX,0:IMAX) :: check_point
#if !defined(ALLOW_TAPENADE)
logical, save                     :: firstcall = .true.
#endif
#if (TSURFACE==6 && ACCSURFACE==6 && ABLSURFACE==6)
integer(i4b)       :: n_year_CE_aux
#if !defined(ALLOW_TAPENADE)
integer(i4b), save :: n_year_CE_aux_save = -9999
#endif
integer(i4b)       :: n_cnt
real(dp)           :: delta_ts_sum
character(len= 16) :: ch_year_CE
character(len=256) :: filename_with_path
real(dp), dimension(0:IMAX,0:JMAX) :: temp_maat_anom_conv, &
                                      dtemp_maat_dz_conv, &
                                      smb_anom_conv, &
                                      dsmb_dz_conv
integer(i4b)        :: n_smb_anom_unit_length, n_dsmb_dz_unit_length
character(len=64)   :: ch_smb_anom_unit, ch_dsmb_dz_unit
#endif

#if (ACCSURFACE==7 && ABLSURFACE==7)
real(dp) :: target_topo_tau_inv
real(dp) :: smb_no_ice
#endif

#if (RETREAT_MASK==1)
integer(i4b)       :: n_year_CE_rtr
integer(i4b), save :: n_year_CE_rtr_save = -9999
character(len= 16) :: ch_year_CE_rtr
character(len=256) :: filename_rtr
real(dp), dimension(0:IMAX,0:JMAX) :: r_mask_retreat_conv
#endif

integer(i4b) :: ncid
!     ncid:      File ID
integer(i4b) :: ncv
!     ncv:       Variable ID
integer(i4b) :: nc3cor(3)
!     nc3cor(3): Corner of a 3-d array
integer(i4b) :: nc3cnt(3)
!     nc3cnt(3): Count of a 3-d array

real(dp), parameter :: &
          inv_twelve = 1.0_dp/12.0_dp, one_third = 1.0_dp/3.0_dp

character(len=64), parameter :: thisroutine = 'boundary'

#if defined(ALLOW_TAPENADE) /* Tapenade */
integer(i4b) :: i_time_in_years
real(dp)     :: temp_val
#endif /* Tapenade */

time_in_years = time*sec2year
n_year_CE     = floor((time_in_years+YEAR_ZERO)+eps_sp_dp)

rho_inv       = 1.0_dp/RHO

!-------- Initialization of variables --------

z_sl_old      = z_sl
z_sl_old_mean = z_sl_mean

delta_ts   = 0.0_dp
glac_index = 0.0_dp
z_sl       = 0.0_dp
dzsl_dtau  = 0.0_dp
z_mar      = 0.0_dp

!-------- Surface-temperature deviation from present values --------

#if (TSURFACE==1)
delta_ts = DELTA_TS0
!                           ! Steady state with prescribed constant
!                           ! air-temperature deviation
#elif (TSURFACE==3)
delta_ts = SINE_AMPLIT &
           *cos(2.0_dp*pi*time_in_years/SINE_PERIOD) &
           -SINE_AMPLIT
!                           ! Sinusoidal air-temperature forcing
#elif (TSURFACE==4)

!  ------ delta_ts from ice-core record

if (time_in_years < real(grip_time_min,dp)) then
   delta_ts = griptemp(0)
else if (time_in_years < real(grip_time_max,dp)) then

#if !defined(ALLOW_TAPENADE) /* Normal */
   i_kl = floor((time_in_years &
          -real(grip_time_min,dp))/real(grip_time_stp,dp))
#else /* Tapenade */
   call myfloor((time_in_years &
                -real(grip_time_min,dp))/real(grip_time_stp,dp),i_kl)
#endif /* Normal vs. Tapenade */
   i_kl = max(i_kl, 0)

#if !defined(ALLOW_TAPENADE) /* Normal */
   i_gr = ceiling((time_in_years &
          -real(grip_time_min,dp))/real(grip_time_stp,dp))
#else /* Tapenade */
   call myceiling(((time_in_years &
          -real(grip_time_min,dp))/real(grip_time_stp,dp)),i_gr)
#endif /* Normal vs. Tapenade */
   i_gr = min(i_gr, ndata_grip)

   if (i_kl == i_gr) then

      delta_ts = griptemp(i_kl)

   else

      time_kl = (grip_time_min + i_kl*grip_time_stp) *year2sec
      time_gr = (grip_time_min + i_gr*grip_time_stp) *year2sec

      delta_ts = griptemp(i_kl) &
                +(griptemp(i_gr)-griptemp(i_kl)) &
                *(time-time_kl)/(time_gr-time_kl)
                 ! linear interpolation of the ice-core data

   end if

else
   delta_ts  = griptemp(ndata_grip)
end if

delta_ts = delta_ts * GRIP_TEMP_FACT
!                ! modification by constant factor

!-------- Glacial index --------

#elif (TSURFACE==5)

if (time_in_years < real(gi_time_min,dp)) then
   glac_index = glacial_index(0)
else if (time_in_years < real(gi_time_max,dp)) then

#if !defined(ALLOW_TAPENADE)/* Normal */
   i_kl = floor((time_in_years &
          -real(gi_time_min,dp))/real(gi_time_stp,dp))
#else /* Tapenade */
   call myfloor(((time_in_years &
          -real(gi_time_min,dp))/real(gi_time_stp,dp)), i_kl)
#endif /* Normal vs. Tapenade */
   i_kl = max(i_kl, 0)

#if !defined(ALLOW_TAPENADE)/* Normal */
   i_gr = ceiling((time_in_years &
          -real(gi_time_min,dp))/real(gi_time_stp,dp))
#else /* Tapenade */
   call myceiling(((time_in_years &
          -real(gi_time_min,dp))/real(gi_time_stp,dp)), i_gr)
#endif /* Normal vs. Tapenade */
   i_gr = min(i_gr, ndata_gi)

   if (i_kl == i_gr) then

      glac_index = glacial_index(i_kl)

   else

      time_kl = (gi_time_min + i_kl*gi_time_stp) *year2sec
      time_gr = (gi_time_min + i_gr*gi_time_stp) *year2sec

      glac_index = glacial_index(i_kl) &
                +(glacial_index(i_gr)-glacial_index(i_kl)) &
                *(time-time_kl)/(time_gr-time_kl)
                 ! linear interpolation of the glacial-index data

   end if

else
   glac_index  = glacial_index(ndata_gi)
end if

!-------- Reading of surface temperature and SMB
!                                        directly from data --------

#elif (TSURFACE==6 && ACCSURFACE==6 && ABLSURFACE==6)

n_year_CE_aux = n_year_CE

if (n_year_CE_aux < TEMP_SMB_ANOM_TIME_MIN) then
   n_year_CE_aux = TEMP_SMB_ANOM_TIME_MIN
else if (n_year_CE_aux > TEMP_SMB_ANOM_TIME_MAX) then
   n_year_CE_aux = TEMP_SMB_ANOM_TIME_MAX
end if

if ( firstcall.or.(n_year_CE_aux /= n_year_CE_aux_save) ) then

   write(ch_year_CE, '(i0)') n_year_CE_aux

!  ------ Surface-temperature anomaly

   if ( (trim(adjustl(TEMP_ANOM_FILES)) /= 'none') &
        .and. &
        (trim(adjustl(TEMP_ANOM_FILES)) /= 'None') &
        .and. &
        (trim(adjustl(TEMP_ANOM_FILES)) /= 'NONE') ) then

      filename_with_path = trim(TEMP_SMB_ANOM_DIR)//'/'// &
                           trim(TEMP_ANOM_SUBDIR)//'/'// &
                           trim(TEMP_ANOM_FILES)//trim(ch_year_CE)//'.nc'

      ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

      if (ios /= nf90_noerr) then
         errormsg = ' >>> boundary: Error when opening the file' &
                  //                end_of_line &
                  //'               for the surface-temperature anomaly!'
         call error(errormsg)
      end if

      call check( nf90_inq_varid(ncid, 'aST', ncv), thisroutine )
      call check( nf90_get_var(ncid, ncv, temp_maat_anom_conv), thisroutine )

      call check( nf90_close(ncid), thisroutine )

   else

      temp_maat_anom_conv = 0.0_dp

   end if

!  ------ Surface-temperature vertical gradient

   if ( (trim(adjustl(dTEMPdz_FILES)) /= 'none') &
        .and. &
        (trim(adjustl(dTEMPdz_FILES)) /= 'None') &
        .and. &
        (trim(adjustl(dTEMPdz_FILES)) /= 'NONE') ) then

      if ( (trim(adjustl(dTEMPdz_SUBDIR)) /= 'value') &
           .and. &
           (trim(adjustl(dTEMPdz_SUBDIR)) /= 'Value') &
           .and. &
           (trim(adjustl(dTEMPdz_SUBDIR)) /= 'VALUE') ) then  ! read from file

         filename_with_path = trim(TEMP_SMB_ANOM_DIR)//'/'// &
                              trim(dTEMPdz_SUBDIR)//'/'// &
                              trim(dTEMPdz_FILES)//trim(ch_year_CE)//'.nc'

         ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

         if (ios /= nf90_noerr) then
            errormsg = ' >>> boundary: Error when opening the file' &
                     //                end_of_line &
                     //'        for the surface-temperature vertical gradient!'
            call error(errormsg)
         end if

         call check( nf90_inq_varid(ncid, 'dSTdz', ncv), thisroutine )
         call check( nf90_get_var(ncid, ncv, dtemp_maat_dz_conv), thisroutine )

         call check( nf90_close(ncid), thisroutine )

      else   ! use constant value

         ch_aux = trim(adjustl(dTEMPdz_FILES))
         read(ch_aux, *) r_aux
         dtemp_maat_dz_conv = r_aux

      end if

   else

      dtemp_maat_dz_conv = 0.0_dp

   end if

!  ------ SMB anomaly

   if ( (trim(adjustl(SMB_ANOM_FILES)) /= 'none') &
        .and. &
        (trim(adjustl(SMB_ANOM_FILES)) /= 'None') &
        .and. &
        (trim(adjustl(SMB_ANOM_FILES)) /= 'NONE') ) then

      filename_with_path = trim(TEMP_SMB_ANOM_DIR)//'/'// &
                           trim(SMB_ANOM_SUBDIR)//'/'// &
                           trim(SMB_ANOM_FILES)//trim(ch_year_CE)//'.nc'

      ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

      if (ios /= nf90_noerr) then
         errormsg = ' >>> boundary: Error when opening the file' &
                  //                end_of_line &
                  //'               for the SMB anomaly!'
         call error(errormsg)
      end if

      call check( nf90_inq_varid(ncid, 'aSMB', ncv), thisroutine )
      call check( nf90_get_var(ncid, ncv, smb_anom_conv), thisroutine )
      call check( nf90_inquire_attribute(ncid, ncv, 'units', &
                                         len=n_smb_anom_unit_length) )
      call check( nf90_get_att(ncid, ncv, 'units', ch_smb_anom_unit) )

      call check( nf90_close(ncid), thisroutine )

      if (trim(adjustl(ch_smb_anom_unit))=='kg m-2 s-1') then
         smb_anom_conv = smb_anom_conv *rho_inv
                                       ! kg/(m2*s) -> m/s ice equiv.
      else if (trim(adjustl(ch_smb_anom_unit))=='m a-1') then
         smb_anom_conv = smb_anom_conv *sec2year
                                       ! m/a ice equiv. -> m/s ice equiv.
      else
         errormsg = ' >>> boundary: ' &
                       //'Unit of smb_anom_conv could not be determined!'
         call error(errormsg)
      end if

   else

      smb_anom_conv = 0.0_dp

   end if

!  ------ SMB vertical gradient

   if ( (trim(adjustl(dSMBdz_FILES)) /= 'none') &
        .and. &
        (trim(adjustl(dSMBdz_FILES)) /= 'None') &
        .and. &
        (trim(adjustl(dSMBdz_FILES)) /= 'NONE') ) then

      if ( (trim(adjustl(dSMBdz_SUBDIR)) /= 'value') &
           .and. &
           (trim(adjustl(dSMBdz_SUBDIR)) /= 'Value') &
           .and. &
           (trim(adjustl(dSMBdz_SUBDIR)) /= 'VALUE') ) then  ! read from file

         filename_with_path = trim(TEMP_SMB_ANOM_DIR)//'/'// &
                              trim(dSMBdz_SUBDIR)//'/'// &
                              trim(dSMBdz_FILES)//trim(ch_year_CE)//'.nc'

         ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

         if (ios /= nf90_noerr) then
            errormsg = ' >>> boundary: Error when opening the file' &
                     //                end_of_line &
                     //'               for the SMB vertical gradient!'
            call error(errormsg)
         end if

         call check( nf90_inq_varid(ncid, 'dSMBdz', ncv), thisroutine )
         call check( nf90_get_var(ncid, ncv, dsmb_dz_conv), thisroutine )
         call check( nf90_inquire_attribute(ncid, ncv, 'units', &
                                            len=n_dsmb_dz_unit_length) )
         call check( nf90_get_att(ncid, ncv, 'units', ch_dsmb_dz_unit) )

         call check( nf90_close(ncid), thisroutine )

         if ( &
              (trim(adjustl(ch_dsmb_dz_unit))=='kg m-2 s-1 m-1') &
              .or. &
              (trim(adjustl(ch_dsmb_dz_unit))=='kg m-3 s-1') &
            ) then
            dsmb_dz_conv = dsmb_dz_conv *rho_inv
                                      ! [kg/(m2*s)]/m -> [m/s ice equiv.]/m
         else if ( &
                   (trim(adjustl(ch_dsmb_dz_unit))=='m a-1 m-1') &
                   .or. &
                   (trim(adjustl(ch_dsmb_dz_unit))=='a-1') &
                 ) then
            dsmb_dz_conv = dsmb_dz_conv *sec2year
                                      ! [m/a ice equiv.]/m -> [m/s ice equiv.]/m
         else
            errormsg = ' >>> boundary: ' &
                          //'Unit of dsmb_dz_conv could not be determined!'
            call error(errormsg)
         end if

      else   ! use constant value

         ch_aux = trim(adjustl(dSMBdz_FILES))
         read(ch_aux, *) r_aux
         dsmb_dz_conv = r_aux *sec2year
                              ! [m/a ice equiv.]/m -> [m/s ice equiv.]/m
      end if

   else

      dsmb_dz_conv = 0.0_dp

   end if

!  ------ Conversion of all data

   do i=0, IMAX
   do j=0, JMAX
      temp_maat_anom(j,i) = temp_maat_anom_conv(i,j)
      dtemp_maat_dz(j,i)  = dtemp_maat_dz_conv(i,j)
      smb_anom(j,i)       = smb_anom_conv(i,j)
      dsmb_dz(j,i)        = dsmb_dz_conv(i,j)
   end do
   end do

end if

!  ------ Average surface temperature anomaly

n_cnt        = 0
delta_ts_sum = 0.0_dp

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i) /= 2) then   ! no ocean points
      n_cnt        = n_cnt + 1
      delta_ts_sum = delta_ts_sum + temp_maat_anom(j,i)
   end if

end do
end do

if (n_cnt > 0) then
   delta_ts = delta_ts_sum/real(n_cnt,dp)
else
   delta_ts = no_value_neg_2
end if

!  ------ Save value of n_year_CE_aux

n_year_CE_aux_save = n_year_CE_aux

#endif

!-------- Global annual temperature anomaly
!         (for parameterising the sub-ocean temperature anomaly
!                             for the ice discharge parameterisation) --------

#if (DISC==1)

dT_glann = 0.0_dp   ! anomaly set to zero

#elif (DISC==2)

!  ------ dT_glann from CLIMBER simulations

if (time_in_years < real(glann_time_min,dp)) then
   dT_glann = dT_glann_CLIMBER(0)
else if (time_in_years < real(glann_time_max,dp)) then

#if !defined(ALLOW_TAPENADE) /* Normal */
   i_kl = floor((time_in_years &
          -real(glann_time_min,dp))/real(glann_time_stp,dp))
   i_kl = max(i_kl, 0)

   i_gr = ceiling((time_in_years &
          -real(glann_time_min,dp))/real(glann_time_stp,dp))
   i_gr = min(i_gr, ndata_glann)
#else /* Tapenade */
   call myfloor(((time_in_years &
          -real(glann_time_min,dp))/real(glann_time_stp,dp)), i_kl)
   i_kl = max(i_kl, 0)

   call myceiling(((time_in_years &
          -real(glann_time_min,dp))/real(glann_time_stp,dp)), i_gr)
   i_gr = min(i_gr, ndata_glann)
#endif /* Normal vs. Tapenade */

   if (i_kl == i_gr) then

      dT_glann = dT_glann_CLIMBER(i_kl)

   else

      time_kl = (glann_time_min + i_kl*glann_time_stp) *year2sec
      time_gr = (glann_time_min + i_gr*glann_time_stp) *year2sec

      dT_glann = dT_glann_CLIMBER(i_kl) &
                +(dT_glann_CLIMBER(i_gr)-dT_glann_CLIMBER(i_kl)) &
                *(time-time_kl)/(time_gr-time_kl)
                 ! linear interpolation of the data

   end if

else
   dT_glann  = dT_glann_CLIMBER(ndata_glann)
end if

#endif

!-------- Sea level --------

#if (SEA_LEVEL==1)

!  ------ Temporally constant sea level

z_sl = Z_SL0

#elif (SEA_LEVEL==3)

!  ------ Time-dependent sea level from data

if (time_in_years < real(specmap_time_min,dp)) then
   z_sl = specmap_zsl(0)
else if (time_in_years < real(specmap_time_max,dp)) then

#if !defined(ALLOW_TAPENADE) /* Normal */
   i_kl = floor((time_in_years &
          -real(specmap_time_min,dp))/real(specmap_time_stp,dp))
   i_kl = max(i_kl, 0)

   i_gr = ceiling((time_in_years &
          -real(specmap_time_min,dp))/real(specmap_time_stp,dp))
   i_gr = min(i_gr, ndata_specmap)
#else /* Tapenade */
   call myfloor(((time_in_years &
          -real(specmap_time_min,dp))/real(specmap_time_stp,dp)), i_kl)
   i_kl = max(i_kl, 0)

   call myceiling(((time_in_years &
          -real(specmap_time_min,dp))/real(specmap_time_stp,dp)), i_gr)
   i_gr = min(i_gr, ndata_specmap)
#endif /* Normal vs. Tapenade */

   if (i_kl == i_gr) then

      z_sl = specmap_zsl(i_kl)

   else

      time_kl = (specmap_time_min + i_kl*specmap_time_stp) *year2sec
      time_gr = (specmap_time_min + i_gr*specmap_time_stp) *year2sec

      z_sl = specmap_zsl(i_kl) &
                +(specmap_zsl(i_gr)-specmap_zsl(i_kl)) &
                *(time-time_kl)/(time_gr-time_kl)
                 ! linear interpolation of the sea-level data

   end if

else
   z_sl  = specmap_zsl(ndata_specmap)
end if

#else

errormsg = ' >>> boundary: Parameter SEA_LEVEL must be either 1 or 3!'
call error(errormsg)

#endif

!  ------ Mean sea level

z_sl_mean = sum(z_sl*cell_area)/sum(cell_area)

!  ------ Time derivative of the sea level

if ( z_sl_old_mean > -999999.9_dp ) then
   dzsl_dtau = (z_sl-z_sl_old)/dtime
else   ! only dummy value for z_sl_old available
   dzsl_dtau = 0.0_dp
end if

!  ------ Minimum bedrock elevation for extent of marine ice

#if (MARGIN==2)

#if (MARINE_ICE_CALVING==2 || MARINE_ICE_CALVING==3)
z_mar = Z_MAR
#elif (MARINE_ICE_CALVING==4 || MARINE_ICE_CALVING==5)
z_mar = FACT_Z_MAR*z_sl_mean
#elif (MARINE_ICE_CALVING==6 || MARINE_ICE_CALVING==7)
if (z_sl_mean >= -80.0_dp) then
   z_mar = 2.5_dp*z_sl_mean
else
   z_mar = 10.25_dp*(z_sl_mean+80.0_dp)-200.0_dp
end if
z_mar = FACT_Z_MAR*z_mar
#endif

#endif

!  ------ Update of the mask according to the sea level

!    ---- Check all sea and floating-ice points and their direct
!         neighbours

do i=0, IMAX
do j=0, JMAX
   check_point(j,i) = .false.
end do
end do

do i=1, IMAX-1
do j=1, JMAX-1
   if (mask(j,i) >= 2) then
      check_point(j  ,i  ) = .true.
      check_point(j  ,i+1) = .true.
      check_point(j  ,i-1) = .true.
      check_point(j+1,i  ) = .true.
      check_point(j-1,i  ) = .true.
   end if
end do
end do

do i=1, IMAX-1
do j=1, JMAX-1
   if (check_point(j,i)) then
      mask_new(j,i) = mask_update_sea_level(i, j)
   end if
end do
end do

!    ---- Assign new values of the mask

do i=1, IMAX-1
do j=1, JMAX-1
   if (check_point(j,i)) then
      mask(j,i) = mask_new(j,i)
   end if
end do
end do

!-------- Surface air temperatures --------

#if (TSURFACE<=5)

#if (TEMP_PRESENT_PARA==1)   /* Parameterization by Ritz et al. (1997) */

#if (defined(TEMP_PRESENT_OFFSET))
theta_ma_offset = TEMP_PRESENT_OFFSET
theta_mj_offset = TEMP_PRESENT_OFFSET
#else
theta_ma_offset = 0.0_dp
theta_mj_offset = 0.0_dp
#endif

theta_ma = 49.13_dp + theta_ma_offset
gamma_ma = -7.992e-03_dp
c_ma     = -0.7576_dp
kappa_ma =  0.0_dp

theta_mj = 30.38_dp + theta_mj_offset
gamma_mj = -6.277e-03_dp
c_mj     = -0.3262_dp
kappa_mj =  0.0_dp

#elif (TEMP_PRESENT_PARA==2)   /* Parameterization by Fausto et al. (2009) */

#if (defined(TEMP_PRESENT_OFFSET))
theta_ma_offset = TEMP_PRESENT_OFFSET
theta_mj_offset = TEMP_PRESENT_OFFSET
#else
theta_ma_offset = 0.0_dp
theta_mj_offset = 0.0_dp
#endif

theta_ma = 41.83_dp + theta_ma_offset
gamma_ma = -6.309e-03_dp
c_ma     = -0.7189_dp
kappa_ma = -0.0672_dp

theta_mj = 14.70_dp + theta_mj_offset
gamma_mj = -5.426e-03_dp
c_mj     = -0.1585_dp
kappa_mj = -0.0518_dp

#else

errormsg = ' >>> boundary: Parameter TEMP_PRESENT_PARA must be either 1 or 2!'
call error(errormsg)

#endif

#endif

do i=0, IMAX
do j=0, JMAX

!  ------ Present-day mean-annual air temperature

#if (TSURFACE<=5)

   temp_ma_present(j,i) = theta_ma &
                  + gamma_ma*zs(j,i) &
                  + c_ma*phi(j,i)*rad2deg &
#if !defined(ALLOW_TAPENADE) /* Normal */
                  + kappa_ma*(modulo(lambda(j,i)+pi,2.0_dp*pi)-pi)*rad2deg
                              ! west longitudes counted negatively
#else /* Tapenade */
                  + kappa_ma*((lambda(j,i)+pi - &
                              (2.0_dp*pi * int((lambda(j,i)+pi)/(2.0_dp*pi)))) &
                              - pi) &
                            *rad2deg
#endif /* Normal vs. Tapenade */

#elif (TSURFACE==6)

   temp_ma_present(j,i) = temp_maat_climatol(j,i)

#endif

!  ------ Present-day mean-July (summer) air temperature

#if (TSURFACE<=5)

   temp_mj_present(j,i) = theta_mj &
                    + gamma_mj*zs(j,i) &
                    + c_mj*phi(j,i)*rad2deg &
#if !defined(ALLOW_TAPENADE) /* Normal */
                    + kappa_mj*(modulo(lambda(j,i)+pi,2.0_dp*pi)-pi)*rad2deg
                                ! west longitudes counted negatively
#else /* Tapenade */
                    + kappa_mj*((lambda(j,i)+pi - &
                                (2.0_dp*pi * int((lambda(j,i)+pi)/(2.0_dp*pi)))) &
                                - pi) &
                              *rad2deg
#endif /* Normal vs. Tapenade */

#elif (TSURFACE==6)

   temp_mj_present(j,i) = temp_ma_present(j,i)
                   ! not contained in read data, thus set to dummy value

#endif

#if (TSURFACE<=4)

!  ------ Correction with deviation delta_ts

   temp_ma(j,i)   = temp_ma_present(j,i) + delta_ts
   temp_mm(j,i,7) = temp_mj_present(j,i) + delta_ts

#elif (TSURFACE==5)

!  ------ Correction with LGM anomaly and glacial index

   temp_ma(j,i)   = temp_ma_present(j,i) + glac_index*temp_ma_lgm_anom(j,i)
   temp_mm(j,i,7) = temp_mj_present(j,i) + glac_index*temp_mj_lgm_anom(j,i)

#elif (TSURFACE==6)

!  ------ Mean-annual and mean-July (summer) surface temperatures
!                                            from data read above --------

   temp_ma(j,i)   = temp_maat_climatol(j,i) &
                       + temp_maat_anom(j,i) &
                          + dtemp_maat_dz(j,i)*(zs(j,i)-zs_ref(j,i))

   temp_mm(j,i,7) = temp_ma(j,i)
                   ! not contained in read data, thus set to dummy value

#endif

!  ------ Amplitude of the annual cycle

   temp_ampl(j,i) = temp_mm(j,i,7) - temp_ma(j,i)

   if (temp_ampl(j,i) < eps) then
      temp_ampl(j,i) = eps   ! Correction of amplitude, if required
   end if

end do
end do

!  ------ Monthly temperatures

do n=1, 12   ! month counter

   sine_factor = sin((real(n,dp)-4.0_dp)*pi/6.0_dp)

   do i=0, IMAX
   do j=0, JMAX
      temp_mm(j,i,n) = temp_ma(j,i) + sine_factor*temp_ampl(j,i)
   end do
   end do

end do

!  ------ Save mean-annual air temperature

temp_maat = temp_ma

!-------- Accumulation-ablation function as_perp --------

#if (ACCSURFACE<=3)

#if (ELEV_DESERT==1)

gamma_p   = GAMMA_P*1.0e-03_dp   ! Precipitation lapse rate
                                 ! for elevation desertification, in m^(-1)
zs_thresh = ZS_THRESH            ! Elevation threshold, in m

#endif

#endif

#if (ACCSURFACE<=5)

#if (SOLID_PRECIP==1)     /* Marsiat (1994) */

temp_rain =    7.0_dp   ! Threshold monthly mean temperature for
                        ! precipitation = 100% rain, in degC
temp_snow =  -10.0_dp   ! Threshold monthly mean temperature for &
                        ! precipitation = 100% snow, in degC

inv_delta_temp_rain_snow = 1.0_dp/(temp_rain-temp_snow)

#elif (SOLID_PRECIP==2)   /* Bales et al. (2009) */

temp_rain =    7.2_dp   ! Threshold monthly mean temperature for &
                        ! precipitation = 100% rain, in degC
temp_snow =  -11.6_dp   ! Threshold monthly mean temperature for &
                        ! precipitation = 100% snow, in degC

coeff(0) =  5.4714e-01_dp   ! Coefficients
coeff(1) = -9.1603e-02_dp   ! of
coeff(2) = -3.314e-03_dp    ! the
coeff(3) =  4.66e-04_dp     ! fifth-order
coeff(4) =  3.8e-05_dp      ! polynomial
coeff(5) =  6.0e-07_dp      ! fit

#elif (SOLID_PRECIP==3)   /* Huybrechts and de Wolde (1999) */

temp_rain = 2.0_dp      ! Threshold instantaneous temperature for &
                        ! precipitation = 100% rain, in degC
temp_snow = temp_rain   ! Threshold instantaneous temperature for &
                        ! precipitation = 100% snow, in degC

#if (defined(S_STAT_0))
s_stat = S_STAT_0    ! Standard deviation of the air termperature
                     ! (same parameter as in the PDD model)
#else
errormsg = ' >>> boundary: ' &
           // 'Parameters for PDD model not defined in run-specs header!'
call error(errormsg)
#endif

inv_sqrt2_s_stat = 1.0_dp/(sqrt(2.0_dp)*s_stat)

#endif

#endif

#if (ABLSURFACE==1 || ABLSURFACE==2)

#if (defined(S_STAT_0) && defined(BETA1_0) && defined(BETA2_0) \
                       && defined(PMAX_0) && defined(MU_0))
s_stat = S_STAT_0
beta1  = BETA1_0  *(0.001_dp/86400.0_dp)*(RHO_W/RHO)
                           ! (mm WE)/(d*degC) -> (m IE)/(s*degC)
beta2  = BETA2_0  *(0.001_dp/86400.0_dp)*(RHO_W/RHO)
                           ! (mm WE)/(d*degC) -> (m IE)/(s*degC)
Pmax   = PMAX_0
mu     = MU_0     *(1000.0_dp*86400.0_dp)*(RHO/RHO_W)
                           ! (d*degC)/(mm WE) -> (s*degC)/(m IE)
#else
errormsg = ' >>> boundary: ' &
           // 'Parameters for PDD model not defined in run-specs header!'
call error(errormsg)
#endif

#elif (ABLSURFACE==3)

lambda_lti = LAMBDA_LTI  *(0.001_dp*sec2year)*(RHO_W/RHO)
                         ! (mm WE)/(a*degC) -> (m IE)/(s*degC)
temp_lti   = TEMP_LTI
mu         = 0.0_dp      ! no superimposed ice considered

#endif

do i=0, IMAX
do j=0, JMAX

!  ------ Accumulation

#if (ACCSURFACE<=3)

!    ---- Elevation desertification of precipitation

#if (ELEV_DESERT==0)

   precip_fact = 1.0_dp   ! no elevation desertification

#elif (ELEV_DESERT==1)

   if (zs_ref(j,i) < zs_thresh) then
      precip_fact &
         = exp(gamma_p*(max(zs(j,i),zs_thresh)-zs_thresh))
   else
      precip_fact &
         = exp(gamma_p*(max(zs(j,i),zs_thresh)-zs_ref(j,i)))
   end if

#else
   errormsg = ' >>> boundary: Parameter ELEV_DESERT must be either 0 or 1!'
   call error(errormsg)
#endif

   do n=1, 12   ! month counter
      precip(j,i,n) = precip_present(j,i,n)*precip_fact
   end do

#endif

!    ---- Precipitation change related to changing climate

#if (ACCSURFACE==1)
   precip_fact = ACCFACT
#elif (ACCSURFACE==2)
   precip_fact = 1.0_dp + GAMMA_S*delta_ts
#elif (ACCSURFACE==3)
   precip_fact = exp(GAMMA_S*delta_ts)
#endif

#if (ACCSURFACE<=3)

   precip(j,i,0) = 0.0_dp   ! initialisation value for mean annual precip

   do n=1, 12   ! month counter
      precip(j,i,n) = precip(j,i,n)*precip_fact   ! monthly precip
      precip(j,i,0) = precip(j,i,0) + precip(j,i,n)*inv_twelve
                                                  ! mean annual precip
   end do

#elif (ACCSURFACE==5)

   precip(j,i,0) = 0.0_dp   ! initialisation value for mean annual precip

   do n=1, 12   ! month counter

#if (PRECIP_ANOM_INTERPOL==1)
      precip_fact = 1.0_dp-glac_index+glac_index*precip_lgm_anom(j,i,n)
                    ! interpolation with a linear function
#elif (PRECIP_ANOM_INTERPOL==2)
      precip_fact = exp(-glac_index*gamma_precip_lgm_anom(j,i,n))
                    ! interpolation with an exponential function
#endif

      precip(j,i,n) = precip_present(j,i,n)*precip_fact   ! monthly precip
      precip(j,i,0) = precip(j,i,0) + precip(j,i,n)*inv_twelve
                                                      ! mean annual precip
   end do

#elif (ACCSURFACE==6 || ACCSURFACE==7)

   precip(j,i,0) = 0.0_dp   ! only total SMB computed (further below)

   do n=1, 12   ! month counter
      precip(j,i,n) = 0.0_dp   ! only total SMB computed (further below)
   end do

#endif

!    ---- Annual accumulation, snowfall and rainfall rates

#if (ACCSURFACE<=5)

   accum(j,i) = precip(j,i,0)

   snowfall(j,i) = 0.0_dp   ! initialisation value

   do n=1, 12   ! month counter

#if (SOLID_PRECIP==1)     /* Marsiat (1994) */

      if (temp_mm(j,i,n) >= temp_rain) then
         frac_solid = 0.0_dp
      else if (temp_mm(j,i,n) <= temp_snow) then
         frac_solid = 1.0_dp
      else
         frac_solid = (temp_rain-temp_mm(j,i,n))*inv_delta_temp_rain_snow
      end if

#elif (SOLID_PRECIP==2)   /* Bales et al. (2009) */

      if (temp_mm(j,i,n) >= temp_rain) then
         frac_solid = 0.0_dp
      else if (temp_mm(j,i,n) <= temp_snow) then
         frac_solid = 1.0_dp
      else
         frac_solid = coeff(0) + temp_mm(j,i,n) * ( coeff(1) &
                               + temp_mm(j,i,n) * ( coeff(2) &
                               + temp_mm(j,i,n) * ( coeff(3) &
                               + temp_mm(j,i,n) * ( coeff(4) &
                               + temp_mm(j,i,n) *   coeff(5) ) ) ) )
               ! evaluation of 5th-order polynomial by Horner scheme
      end if

#elif (SOLID_PRECIP==3)   /* Huybrechts and de Wolde (1999) */

#if !defined(ALLOW_TAPENADE) /* Normal */
      frac_solid = 1.0_dp &
                   - 0.5_dp*erfc((temp_rain-temp_mm(j,i,n))*inv_sqrt2_s_stat)
#else /* Tapenade */
      call my_erfc((temp_rain-temp_mm(j,i,n))*inv_sqrt2_s_stat, temp_val)
      frac_solid = 1.0_dp - 0.5_dp*temp_val
#endif /* Normal vs. Tapenade */

#endif

      snowfall(j,i) = snowfall(j,i) + precip(j,i,n)*frac_solid*inv_twelve

   end do

   rainfall(j,i) = precip(j,i,0) - snowfall(j,i)

   if (snowfall(j,i) < 0.0_dp) snowfall(j,i) = 0.0_dp   ! correction of
   if (rainfall(j,i) < 0.0_dp) rainfall(j,i) = 0.0_dp   ! negative values

#elif (ACCSURFACE==6 || ACCSURFACE==7)

   accum(j,i)    = 0.0_dp   ! only total SMB computed (further below)
   snowfall(j,i) = 0.0_dp
   rainfall(j,i) = 0.0_dp

#endif

!  ------ Ablation

!    ---- Runoff

#if (ABLSURFACE==1 || ABLSURFACE==2)

!      -- Temperature excess ET

   do n=1, 12   ! month counter
      temp_mm_help(n) = temp_mm(j,i,n)
   end do

   call pdd(temp_mm_help, s_stat, ET(j,i))

!      -- Formation rate of superimposed ice (melt_star), melt rate (melt)
!         and runoff rate (runoff)

#if (ABLSURFACE==1)

   if ((beta1*ET(j,i)) <= (Pmax*snowfall(j,i))) then
      melt_star(j,i) = beta1*ET(j,i)
      melt(j,i)      = 0.0_dp
      runoff(j,i)    = melt(j,i)+rainfall(j,i)
   else
      melt_star(j,i) = Pmax*snowfall(j,i)
      melt(j,i)      = beta2*(ET(j,i)-melt_star(j,i)/beta1)
      runoff(j,i)    = melt(j,i)+rainfall(j,i)
   end if

#elif (ABLSURFACE==2)

   if ( rainfall(j,i) <= (Pmax*snowfall(j,i)) ) then

      if ( (rainfall(j,i)+beta1*ET(j,i)) <= (Pmax*snowfall(j,i)) ) then
         melt_star(j,i) = rainfall(j,i)+beta1*ET(j,i)
         melt(j,i)      = 0.0_dp
         runoff(j,i)    = melt(j,i)
      else
         melt_star(j,i) = Pmax*snowfall(j,i)
         melt(j,i)      = beta2 &
                          *(ET(j,i)-(melt_star(j,i)-rainfall(j,i))/beta1)
         runoff(j,i)    = melt(j,i)
      end if

   else

      melt_star(j,i) = Pmax*snowfall(j,i)
      melt(j,i)      = beta2*ET(j,i)
      runoff(j,i)    = melt(j,i) + rainfall(j,i)-Pmax*snowfall(j,i)

   end if

#endif

#elif (ABLSURFACE==3)

   temp_jja_help  = one_third*(temp_mm(j,i,6)+temp_mm(j,i,7)+temp_mm(j,i,8))

   melt_star(j,i) = 0.0_dp   ! no superimposed ice considered
   melt(j,i)      = lambda_lti*max((temp_jja_help-temp_lti), 0.0_dp)
   runoff(j,i)    = melt(j,i) + rainfall(j,i)

#elif (ABLSURFACE==6 || ABLSURFACE==7)

   melt_star(j,i) = 0.0_dp   ! only total SMB computed (further below)
   melt(j,i)      = 0.0_dp
   runoff(j,i)    = 0.0_dp

#endif

end do
end do

!  ------ SMB = precipitation minus runoff

#if (ACCSURFACE<=5 && ABLSURFACE<=5)

as_perp = accum - runoff

#elif (ACCSURFACE==6 && ABLSURFACE==6)

as_perp = smb_climatol  + smb_anom + dsmb_dz*(zs-zs_ref)

accum  =  max(as_perp, 0.0_dp)
runoff = -min(as_perp, 0.0_dp)

#elif (ACCSURFACE==7 && ABLSURFACE==7)

target_topo_tau_inv = 1.0_dp/target_topo_tau_0

smb_no_ice = -1000.0_dp*sec2year   ! -1000 m/a -> m/s

do i=0, IMAX
do j=0, JMAX

   if ((mask_target(j,i)==0).or.(mask_target(j,i)==3)) then
      as_perp(j,i) = (zs_target(j,i)-zs(j,i))*target_topo_tau_inv
   else
      as_perp(j,i) = smb_no_ice
   end if

end do
end do

accum  =  max(as_perp, 0.0_dp)
runoff = -min(as_perp, 0.0_dp)

#endif

!    ---- Prescribed SMB correction

smb_corr_prescribed = smb_corr_in

if (flag_initmip_asmb) then   ! Correction for ISMIP InitMIP

#if defined(ALLOW_TAPENADE) /* Tapenade */
   call myfloor(time_in_years, i_time_in_years)
#endif /* Tapenade */

   if ((time_in_years > 0.0_dp).and.(time_in_years <= 40.0_dp)) then

#if !defined(ALLOW_TAPENADE) /* Normal */
      smb_corr_prescribed = smb_corr_prescribed &
                              + 0.025_dp*floor(time_in_years) * smb_anom_initmip
#else /* Tapenade */
      smb_corr_prescribed = smb_corr_prescribed &
                              + 0.025_dp*(i_time_in_years) * smb_anom_initmip
#endif /* Normal vs. Tapenade */

   else if (time_in_years > 40.0_dp) then

      smb_corr_prescribed = smb_corr_prescribed + smb_anom_initmip

   end if

end if

as_perp = as_perp + smb_corr_prescribed

accum_prescribed  =  max(smb_corr_prescribed, 0.0_dp)
runoff_prescribed = -min(smb_corr_prescribed, 0.0_dp)

accum  = accum  + accum_prescribed
runoff = runoff + runoff_prescribed

!  ------ Ice-surface temperature (10-m firn temperature) temp_s,
!         including empirical firn-warming correction due to
!         refreezing meltwater when superimposed ice is formed

#if !defined(ALLOW_TAPENADE) /* Normal */

#if (TSURFACE<=5)

where (melt_star >= melt)
   temp_s = temp_ma + mu*(melt_star-melt)
elsewhere
   temp_s = temp_ma
end where

#elif (TSURFACE==6)

temp_s = temp_ma

#endif

where (temp_s > -0.001_dp) temp_s = -0.001_dp
                            ! Cut-off of positive air temperatures

#else /* Tapenade */

#if (TSURFACE<=5)

do i=0, IMAX
do j=0, JMAX
   if (melt_star(j,i) >= melt(j,i)) then
      temp_s(j,i) = temp_ma(j,i) + mu*(melt_star(j,i)-melt(j,i))
   else
      temp_s(j,i) = temp_ma(j,i)
   end if
end do
end do

#elif (TSURFACE==6)

temp_s = temp_ma

#endif

do i=0, IMAX
do j=0, JMAX
   if (temp_s(j,i) > -0.001_dp) then
      temp_s(j,i) = -0.001_dp   ! Cut-off of positive air temperatures
   end if
end do
end do

#endif /* Normal vs. Tapenade */

!-------- Calving --------

calving = 0.0_dp   ! Initialization

#if ((MARGIN==2) \
      && (MARINE_ICE_FORMATION==2) \
      && (MARINE_ICE_CALVING==9))

call calving_underwater_ice()

#endif

#if (DISC>0)   /* Ice discharge parameterization */

call discharge(dxi, deta)
calving = calving + dis_perp

#endif

!-------- Retreat mask due to oceanic forcing --------

#if (RETREAT_MASK==1)

n_year_CE_rtr = n_year_CE

if (n_year_CE_rtr > RETREAT_MASK_TIME_MAX) &
                       n_year_CE_rtr = RETREAT_MASK_TIME_MAX

if (firstcall) r_mask_retreat = 1.0_dp   ! initialization

if (n_year_CE_rtr < RETREAT_MASK_TIME_MIN) then

   r_mask_retreat = 1.0_dp

else if (n_year_CE_rtr /= n_year_CE_rtr_save) then

   write(ch_year_CE_rtr, '(i0)') n_year_CE_rtr

   filename_rtr = trim(RETREAT_MASK_DIR)//'/'// &
                  trim(RETREAT_MASK_FILES)// &
                  trim(ch_year_CE_rtr)//'.nc'

   ios = nf90_open(trim(filename_rtr), NF90_NOWRITE, ncid)

   if (ios /= nf90_noerr) then
      errormsg = ' >>> boundary: Error when opening the file for the' &
               //                end_of_line &
               //'               retreat mask due to oceanic forcing!'
      call error(errormsg)
   end if

   call check( nf90_inq_varid(ncid, 'sftgif', ncv), thisroutine )
   call check( nf90_get_var(ncid, ncv, r_mask_retreat_conv), thisroutine )

   call check( nf90_close(ncid), thisroutine )

   do i=0, IMAX
   do j=0, JMAX
      r_mask_retreat(j,i) = max(min(r_mask_retreat_conv(i,j), 1.0_dp), 0.0_dp)
                                    ! constrain to interval [0,1]
   end do
   end do

end if

n_year_CE_rtr_save = n_year_CE_rtr

#endif

if (firstcall) firstcall = .false.

end subroutine boundary

!-------------------------------------------------------------------------------

end module boundary_m
!
