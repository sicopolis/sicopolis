!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  b o u n d a r y _ m
!
!! HEINO domain:
!! Computation of the surface temperature (must be less than 0 degC)
!! and of the accumulation-ablation function.
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
!> HEINO domain:
!! Computation of the surface temperature (must be less than 0 degC)
!! and of the accumulation-ablation function.
!-------------------------------------------------------------------------------
module boundary_m

  use sico_types_m
  use sico_variables_m
  use sico_vars_m

  implicit none

  public

contains

!-------------------------------------------------------------------------------
!> Main routine of boundary_m:
!! Computation of the surface temperature (must be less than 0 degC)
!! and of the accumulation-ablation function.
!-------------------------------------------------------------------------------
subroutine boundary(time, dtime, dxi, deta, &
                    delta_ts, glac_index, z_mar)

#if ((MARGIN==2) && (MARINE_ICE_FORMATION==2) && (MARINE_ICE_CALVING==9))
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
real(dp) :: time_gr, time_kl
real(dp), dimension(0:JMAX,0:IMAX) :: dist
real(dp) :: rad_inv
logical, dimension(0:JMAX,0:IMAX) :: check_point

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

if (time*sec2year.lt.real(grip_time_min,dp)) then
   delta_ts = griptemp(0)
else if (time*sec2year.lt.real(grip_time_max,dp)) then

   i_kl = floor(((time*sec2year) &
          -real(grip_time_min,dp))/real(grip_time_stp,dp))
   i_kl = max(i_kl, 0)

   i_gr = ceiling(((time*sec2year) &
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

if (time*sec2year.lt.real(specmap_time_min,dp)) then
   z_sl = specmap_zsl(0)
else if (time*sec2year.lt.real(specmap_time_max,dp)) then

   i_kl = floor(((time*sec2year) &
          -real(specmap_time_min,dp))/real(specmap_time_stp,dp))
   i_kl = max(i_kl, 0)

   i_gr = ceiling(((time*sec2year) &
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

!-------- Surface air temperature and accumulation-ablation function --------

do i=0, IMAX
do j=0, JMAX
   dist(j,i) = sqrt( (xi(i)-x_hat)**2 + (eta(j)-y_hat)**2 )
end do
end do

rad_inv = 1.0_dp/rad

do i=0, IMAX
do j=0, JMAX

   temp_s(j,i) = temp_min + s_t*dist(j,i)**3
   temp_s(j,i) = temp_s(j,i) + delta_ts
                 ! Correction with temperature deviation delta_ts
   temp_maat(j,i) = temp_s(j,i)
                    ! Save mean-annual air temperature
   if (temp_s(j,i) > -0.001_dp) temp_s(j,i) = -0.001_dp
                               ! Cut-off of positive air temperatures

   accum_present(j,i) = b_min + (b_max-b_min)*rad_inv*dist(j,i)
   accum(j,i)         = accum_present(j,i)

   as_perp(j,i) = accum(j,i)

end do
end do

runoff = accum - as_perp

!-------- Calving --------

calving = 0.0_dp   ! Initialization

#if ((MARGIN==2) && (MARINE_ICE_FORMATION==2) && (MARINE_ICE_CALVING==9))

call calving_underwater_ice()

#endif

if (firstcall%boundary) firstcall%boundary = .false.

end subroutine boundary

!-------------------------------------------------------------------------------

end module boundary_m
!
