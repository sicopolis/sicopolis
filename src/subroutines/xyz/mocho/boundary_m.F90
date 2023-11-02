!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  b o u n d a r y _ m
!
!> @file
!!
!! Computation of the surface temperature (must be less than 0 deg C!)
!! and of the accumulation-ablation function.
!!
!! @section Copyright
!!
!! Copyright 2009-2023 Ralf Greve, Eduardo Flandez, Matthias Scheiter
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
!> Computation of the surface temperature (must be less than 0 deg C!)
!! and of the accumulation-ablation function.
!<------------------------------------------------------------------------------
module boundary_m

  use sico_types_m
  use sico_variables_m
  use sico_vars_m
  use error_m

  implicit none

  public

contains

!-------------------------------------------------------------------------------
!> Main routine of boundary_m:
!! Computation of the surface temperature (must be less than 0 deg C!)
!! and of the accumulation-ablation function.
!<------------------------------------------------------------------------------
subroutine boundary(time, dtime, dxi, deta, &
                    delta_ts, glac_index, z_mar)

#if ((MARGIN==2) \
      && (MARINE_ICE_FORMATION==2) \
      && (MARINE_ICE_CALVING==9))
  use calving_m
#endif

  use mask_update_sea_level_m

implicit none

real(dp), intent(in) :: time, dtime, dxi, deta

real(dp), intent(out)   :: delta_ts, glac_index, z_mar

! Further return variables
! (defined as global variables in module sico_variables_m):
!
!    accum(j,i), runoff(j,i), as_perp(j,i), calving(j,i), temp_s(j,i)

integer(i4b) :: i, j
integer(i4b) :: i_gr, i_kl
real(dp), dimension(0:JMAX,0:IMAX) :: z_sl_old
real(dp) :: z_sl_old_mean
real(dp) :: z_sl_min, t1, t2, t3, t4, t5, t6
real(dp) :: time_gr, time_kl, asp_a, asp_b
real(dp) :: ela_now
real(dp), dimension(0:JMAX,0:IMAX) :: ela_dist, asp_dist, grad_dist
real(dp), dimension(0:JMAX,0:IMAX) :: dist, dist2
real(dp), dimension(0:JMAX,0:IMAX) :: accum_prescribed, &
                                      runoff_prescribed
logical, dimension(0:JMAX,0:IMAX) :: check_point
logical, save                     :: firstcall = .true.

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
           *cos(2.0_dp*pi*time/(SINE_PERIOD*year2sec)) &
           -SINE_AMPLIT
!                           ! Sinusoidal air-temperature forcing
#elif (TSURFACE==4)

!  ------ delta_ts from the GRIP record

if (time/year2sec.lt.real(grip_time_min,dp)) then
   delta_ts = griptemp(0)
else if (time/year2sec.lt.real(grip_time_max,dp)) then

   i_kl = floor(((time/year2sec) &
          -real(grip_time_min,dp))/real(grip_time_stp,dp))
   i_kl = max(i_kl, 0)

   i_gr = ceiling(((time/year2sec) &
          -real(grip_time_min,dp))/real(grip_time_stp,dp))
   i_gr = min(i_gr, ndata_grip)

   if (i_kl.eq.i_gr) then

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

#elif (TSURFACE==5)
!     Enter here delta_ts scenario:

#endif

!-------- Sea level --------

#if (SEA_LEVEL==1)

!  ------ Temporally constant sea level

z_sl = Z_SL0

#elif (SEA_LEVEL==3)

!  ------ Time-dependent sea level from data

if (time/year2sec.lt.real(specmap_time_min,dp)) then
   z_sl = specmap_zsl(0)
else if (time/year2sec.lt.real(specmap_time_max,dp)) then

   i_kl = floor(((time/year2sec) &
          -real(specmap_time_min,dp))/real(specmap_time_stp,dp))
   i_kl = max(i_kl, 0)

   i_gr = ceiling(((time/year2sec) &
          -real(specmap_time_min,dp))/real(specmap_time_stp,dp))
   i_gr = min(i_gr, ndata_specmap)

   if (i_kl.eq.i_gr) then

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

!-------- Surface air temperature and 
!                     accumulation-ablation function (SMB) --------

do i=0, IMAX
do j=0, JMAX

   temp_s(j,i) = temp_0 - gamma_t*zs(j,i)
   temp_s(j,i) = temp_s(j,i) + delta_ts
                 ! Correction with temperature deviation delta_ts
   temp_maat(j,i) = temp_s(j,i)
                    ! Save mean-annual air temperature
   if (temp_s(j,i) > -0.001_dp) temp_s(j,i) = -0.001_dp
                               ! Cut-off of positive air temperatures

end do
end do

#if (!defined(SURFACE_FORCING) || SURFACE_FORCING==1)

errormsg = ' >>> boundary: Option SURFACE_FORCING=1 not supported any more!'
call error(errormsg)

#elif (SURFACE_FORCING==2)

do i=0, IMAX
do j=0, JMAX

   accum(j,i) = s_0

   ela_now = ela + dela_dts * delta_ts

   as_perp(j,i) = m_0*(zs(j,i)-ela_now)
   if (as_perp(j,i) > accum(j,i)) as_perp(j,i) = accum(j,i)

   runoff(j,i) = max((accum(j,i)-as_perp(j,i)), 0.0_dp)

end do
end do

#elif (SURFACE_FORCING==3)

do i=0, IMAX
do j=0, JMAX

   accum(j,i) = s_0

   ela_now = ela + dela_dts * delta_ts

   if (xi(i) /= x_gip .or. eta(j) /= y_gip) then
      ela_dist(j,i) = ela_now &
                      + ela_amp &
                        * sin(atan2((eta(j)-y_gip), (xi(i)-x_gip)) + phi_0)
   else
      ela_dist(j,i) = ela_now
   end if

   as_perp(j,i) = m_0*(zs(j,i)-ela_dist(j,i))
   if (as_perp(j,i) > accum(j,i)) as_perp(j,i) = accum(j,i)

   runoff(j,i) = max((accum(j,i)-as_perp(j,i)), 0.0_dp)

end do
end do

#elif (SURFACE_FORCING==4)

do i=0, IMAX
do j=0, JMAX

   accum(j,i) = s_0

   ela_now = ela + dela_dts * delta_ts

   if (xi(i) /= x_gip .or. eta(j) /= y_gip) then
      ela_dist(j,i) = ela_now &
                      + ela_amp &
                        * sin(atan2((eta(j)-y_gip), (xi(i)-x_gip)) + phi_0)
   else
      ela_dist(j,i) = ela_now
   end if

   if (zs(j,i) <= z_gc) then
      as_perp(j,i) = m_0*(zs(j,i)-ela_dist(j,i))
   else   ! (zs(j,i) > z_gc)
      as_perp(j,i) = m_1*(zs(j,i)-ela_dist(j,i))+m_0*(z_gc-ela_dist(j,i))
   end if

   accum(j,i)  = max(accum(j,i), as_perp(j,i))
   runoff(j,i) = max((accum(j,i)-as_perp(j,i)), 0.0_dp)

end do
end do

#elif (SURFACE_FORCING==5)

asp_a = 0.0_dp
asp_b = 0.0_dp

do i=0, IMAX
do j=0, JMAX

   accum(j,i) = s_0
   
   asp_a = zs(j,i+1) - zs(j,i-1)
   asp_b = zs(j-1,i) - zs(j+1,i)

   if (asp_b < 0.0_dp) then
      asp_dist(j,i) = mod(atan(asp_a/asp_b), 2.0_dp*pi)
   else   ! (asp_b >= 0.0_dp)
      asp_dist(j,i) = atan(asp_a/asp_b) + pi
                      !!! singular for asp_b == 0.0_dp !!!
   end if

   ela_now = ela + dela_dts * delta_ts

   ela_dist(j,i) = ela_now + ela_amp * sin(asp_dist(j,i) - phi_0 + pi)

   if (zs(j,i) <= z_gc) then
      as_perp(j,i) = m_0*(zs(j,i)-ela_dist(j,i))
   else   ! (zs(j,i) > z_gc)
      as_perp(j,i) = m_1*(zs(j,i)-ela_dist(j,i))+m_0*(z_gc-ela_dist(j,i))
   end if

   accum(j,i)  = max(accum(j,i), as_perp(j,i))
   runoff(j,i) = max((accum(j,i)-as_perp(j,i)), 0.0_dp)

end do
end do

#elif (SURFACE_FORCING==6)

asp_a = 0.0_dp
asp_b = 0.0_dp

do i=0, IMAX
do j=0, JMAX

   accum(j,i) = s_0
   
   asp_a = zs(j,i+1) - zs(j,i-1)
   asp_b = zs(j-1,i) - zs(j+1,i)

   if (asp_b < 0.0_dp) then
      asp_dist(j,i) = mod(atan(asp_a/asp_b), 2.0_dp*pi)
   else   ! (asp_b >= 0.0_dp)
      asp_dist(j,i) = atan(asp_a/asp_b) + pi
                      !!! singular for asp_b == 0.0_dp !!!
   end if

   grad_dist(j,i) = sqrt( (asp_a/(2.0_dp*dxi))**2 &
                          + (asp_b/(2.0_dp*deta))**2 ) * 45.0_dp

   ela_now = ela + dela_dts * delta_ts

   ela_dist(j,i) = ela_now + ela_amp * sin(asp_dist(j,i) - phi_0 + pi)
   if (grad_dist(j,i) > tgt) ela_dist(j,i) = ela_dist(j,i) + grad_dist(j,i)*2

   as_perp(j,i) = m_0*(zs(j,i)-ela_dist(j,i))
   if (as_perp(j,i) > accum(j,i)) as_perp(j,i) = accum(j,i)

   runoff(j,i) = max((accum(j,i)-as_perp(j,i)), 0.0_dp)

end do
end do

#elif (SURFACE_FORCING==7)

asp_a = 0.0_dp
asp_b = 0.0_dp

do i=0, IMAX
do j=0, JMAX

   accum(j,i) = s_0
   
   asp_a = zs(j,i+1) - zs(j,i-1)
   asp_b = zs(j-1,i) - zs(j+1,i)

   if (asp_b < 0.0_dp) then
      asp_dist(j,i) = mod(atan(asp_a/asp_b), 2.0_dp*pi)
   else   ! (asp_b >= 0.0_dp)
      asp_dist(j,i) = atan(asp_a/asp_b) + pi
                      !!! singular for asp_b == 0.0_dp !!!
   end if

   ela_now = ela + dela_dts * delta_ts

   ela_dist(j,i) = ela_now + ela_amp * sin(asp_dist(j,i) - phi_0 + pi)

   as_perp(j,i) = m_0*(zs(j,i)-ela_dist(j,i))
   if (as_perp(j,i) > accum(j,i)) as_perp(j,i) = accum(j,i)

   runoff(j,i) = max((accum(j,i)-as_perp(j,i)), 0.0_dp)

end do
end do

#elif (SURFACE_FORCING==8)

do i=0, IMAX
do j=0, JMAX
   dist(j,i) = sqrt( (xi(i)-x_gip)**2 + (eta(j)-y_gip)**2 )
   dist2(j,i) = sqrt( (xi(i)-x_gip2)**2 + (eta(j)-y_gip2)**2 )
end do
end do

do i=0, IMAX
do j=0, JMAX

   accum(j,i) = s_0

   ela_now = ela + dela_dts * delta_ts

   if (dist(j,i) < 3000.0_dp) then
      if (xi(i) /= x_gip .or. eta(j) /= y_gip) then
         ela_dist(j,i) = ela_now &
                         + ela_amp &
                           * sin(atan2((eta(j)-y_gip), (xi(i)-x_gip)) + phi_0)
      else
         ela_dist(j,i) = ela_now
      end if
   else
	ela_dist(j,i) = ela_now
   end if

   if (dist2(j,i) < 1500.0_dp) then
      if (xi(i) /= x_gip2 .or. eta(j) /= y_gip2) then
         ela_dist(j,i) = ela_dist(j,i) &
                         + ela_amp2 &
                           * sin(atan2((eta(j)-y_gip2), (xi(i)-x_gip2)) + phi_0)
      end if
   end if

   as_perp(j,i) = m_0*(zs(j,i)-ela_dist(j,i))
   if (as_perp(j,i) > accum(j,i)) as_perp(j,i) = accum(j,i)

   runoff(j,i) = max((accum(j,i)-as_perp(j,i)), 0.0_dp)

end do
end do

#endif

!  ------ Prescribed SMB correction

smb_corr_prescribed = smb_corr_in

as_perp = as_perp + smb_corr_prescribed

accum_prescribed  =  max(smb_corr_prescribed, 0.0_dp)
runoff_prescribed = -min(smb_corr_prescribed, 0.0_dp)

accum  = accum  + accum_prescribed
runoff = runoff + runoff_prescribed

!-------- Calving --------

calving = 0.0_dp   ! Initialization

#if ((MARGIN==2) \
      && (MARINE_ICE_FORMATION==2) \
      && (MARINE_ICE_CALVING==9))

call calving_underwater_ice()

#endif

if (firstcall) firstcall = .false.

end subroutine boundary

!-------------------------------------------------------------------------------

end module boundary_m
!
