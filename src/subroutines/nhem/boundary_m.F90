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
real(dp), dimension(0:JMAX,0:IMAX) :: z_sl_old
real(dp) :: z_sl_old_mean
real(dp) :: z_sl_min, t1, t2, t3, t4, t5, t6
real(dp) :: time_gr, time_kl
real(dp) :: z_sle_present, z_sle_help
real(dp), dimension(0:JMAX,0:IMAX,0:12) :: precip
real(dp), dimension(0:JMAX,0:IMAX,12) :: temp_mm
real(dp), dimension(0:JMAX,0:IMAX) :: temp_ma
real(dp), dimension(12) :: temp_mm_help
real(dp) :: temp_jja_help
real(dp) :: gamma_t, temp_diff
real(dp) :: gamma_p, zs_thresh, &
            temp_rain, temp_snow, &
            inv_delta_temp_rain_snow, coeff(0:5), inv_sqrt2_s_stat, &
            precip_fact, frac_solid
real(dp) :: s_stat, beta1, beta2, Pmax, mu, lambda_lti, temp_lti
logical, dimension(0:JMAX,0:IMAX) :: check_point
logical, save                     :: firstcall = .true.

real(dp), parameter :: &
          inv_twelve = 1.0_dp/12.0_dp, one_third = 1.0_dp/3.0_dp

!-------- Initialization of variables --------

z_sl_old      = z_sl
z_sl_old_mean = z_sl_mean

delta_ts   = 0.0_dp
glac_index = 0.0_dp
z_sl       = 0.0_dp
dzsl_dtau  = 0.0_dp
z_mar      = 0.0_dp

!-------- Surface-temperature deviation from present values --------

#if TSURFACE==1
delta_ts = DELTA_TS0
!                           ! Steady state with prescribed constant
!                           ! air-temperature deviation
#elif TSURFACE==3
delta_ts = SINE_AMPLIT &
           *cos(2.0_dp*pi*time/(SINE_PERIOD*year2sec)) &
           -SINE_AMPLIT
!                           ! Sinusoidal air-temperature forcing
#elif TSURFACE==4

!  ------ delta_ts from ice-core record

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

!-------- Glacial index --------

#elif TSURFACE==5

if (time/year2sec < real(gi_time_min,dp)) then
   glac_index = glacial_index(0)
else if (time/year2sec < real(gi_time_max,dp)) then

   i_kl = floor(((time/year2sec) &
          -real(gi_time_min,dp))/real(gi_time_stp,dp))
   i_kl = max(i_kl, 0)

   i_gr = ceiling(((time/year2sec) &
          -real(gi_time_min,dp))/real(gi_time_stp,dp))
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

#if MARGIN==2

#if ( MARINE_ICE_CALVING==2 || MARINE_ICE_CALVING==3 )
z_mar = Z_MAR
#elif ( MARINE_ICE_CALVING==4 || MARINE_ICE_CALVING==5 )
z_mar = FACT_Z_MAR*z_sl_mean
#elif ( MARINE_ICE_CALVING==6 || MARINE_ICE_CALVING==7 )
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

gamma_t = -6.5e-03_dp   ! atmospheric lapse rate

do i=0, IMAX
do j=0, JMAX

#if (TSURFACE <= 4)

!  ------ Correction of present monthly temperatures with elevation changes
!         and temperature deviation delta_ts

   temp_diff = gamma_t*(zs(j,i)-zs_ref(j,i)) + delta_ts

   do n=1, 12   ! month counter
      temp_mm(j,i,n) = temp_mm_present(j,i,n) + temp_diff
   end do

#elif (TSURFACE == 5)

!  ------ Correction of present monthly temperatures with LGM anomaly and
!         glacial index as well as elevation changes

   temp_diff = gamma_t*(zs(j,i)-zs_ref(j,i))

   do n=1, 12   ! month counter
      temp_mm(j,i,n) = temp_mm_present(j,i,n) &
                       + glac_index*temp_mm_lgm_anom(j,i,n) &
                       + temp_diff
   end do

#endif

!  ------ Mean annual air temperature

   temp_ma(j,i) = 0.0_dp   ! initialisation value

   do n=1, 12   ! month counter
      temp_ma(j,i) = temp_ma(j,i) + temp_mm(j,i,n)*inv_twelve
   end do

end do
end do

!  ------ Save mean-annual air temperature

temp_maat = temp_ma

!-------- Accumulation-ablation function as_perp --------

#if (ELEV_DESERT == 1)

gamma_p   = GAMMA_P*1.0e-03_dp   ! Precipitation lapse rate
                                 ! for elevation desertification, in m^(-1)
zs_thresh = ZS_THRESH            ! Elevation threshold, in m

#endif

#if (SOLID_PRECIP == 1)     /* Marsiat (1994) */

temp_rain =    7.0_dp   ! Threshold monthly mean temperature for
                        ! precipitation = 100% rain, in deg C
temp_snow =  -10.0_dp   ! Threshold monthly mean temperature for &
                        ! precipitation = 100% snow, in deg C

inv_delta_temp_rain_snow = 1.0_dp/(temp_rain-temp_snow)

#elif (SOLID_PRECIP == 2)   /* Bales et al. (2009) */

temp_rain =    7.2_dp   ! Threshold monthly mean temperature for &
                        ! precipitation = 100% rain, in deg C
temp_snow =  -11.6_dp   ! Threshold monthly mean temperature for &
                        ! precipitation = 100% snow, in deg C

coeff(0) =  5.4714e-01_dp   ! Coefficients
coeff(1) = -9.1603e-02_dp   ! of
coeff(2) = -3.314e-03_dp    ! the
coeff(3) =  4.66e-04_dp     ! fifth-order
coeff(4) =  3.8e-05_dp      ! polynomial
coeff(5) =  6.0e-07_dp      ! fit

#elif (SOLID_PRECIP == 3)   /* Huybrechts and de Wolde (1999) */

temp_rain = 2.0_dp      ! Threshold instantaneous temperature for &
                        ! precipitation = 100% rain, in deg C
temp_snow = temp_rain   ! Threshold instantaneous temperature for &
                        ! precipitation = 100% snow, in deg C

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

#if (ABLSURFACE==1 || ABLSURFACE==2)

#if (defined(S_STAT_0) && defined(BETA1_0) && defined(BETA2_0) \
                       && defined(PMAX_0) && defined(MU_0))
s_stat = S_STAT_0
beta1  = BETA1_0  *(0.001_dp/86400.0_dp)*(RHO_W/RHO)
                        ! (mm WE)/(d*deg C) --> (m IE)/(s*deg C)
beta2  = BETA2_0  *(0.001_dp/86400.0_dp)*(RHO_W/RHO)
                        ! (mm WE)/(d*deg C) --> (m IE)/(s*deg C)
Pmax   = PMAX_0
mu     = MU_0     *(1000.0_dp*86400.0_dp)*(RHO/RHO_W)
                        ! (d*deg C)/(mm WE) --> (s*deg C)/(m IE)
#else
errormsg = ' >>> boundary: ' &
           // 'Parameters for PDD model not defined in run-specs header!'
call error(errormsg)
#endif

#elif (ABLSURFACE==3)

lambda_lti = LAMBDA_LTI  *(0.001_dp/year2sec)*(RHO_W/RHO)
                         ! (mm WE)/(a*deg C) --> (m IE)/(s*deg C)
temp_lti   = TEMP_LTI
mu         = 0.0_dp      ! no superimposed ice considered

#endif

do i=0, IMAX
do j=0, JMAX

!  ------ Accumulation

#if (ACCSURFACE <= 3)

!    ---- Elevation desertification of precipitation

#if (ELEV_DESERT == 0)

   precip_fact = 1.0_dp   ! no elevation desertification

#elif (ELEV_DESERT == 1)

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

#if ACCSURFACE==1
   precip_fact = ACCFACT
#elif ACCSURFACE==2
   precip_fact = 1.0_dp + GAMMA_S*delta_ts
#elif ACCSURFACE==3
   precip_fact = exp(GAMMA_S*delta_ts)
#endif

#if (ACCSURFACE <= 3)

   precip(j,i,0) = 0.0_dp   ! initialisation value for mean annual precip

   do n=1, 12   ! month counter
      precip(j,i,n) = precip(j,i,n)*precip_fact   ! monthly precip
      precip(j,i,0) = precip(j,i,0) + precip(j,i,n)*inv_twelve
                                              ! mean annual precip
   end do

#elif (ACCSURFACE == 5)

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

#endif

!    ---- Annual accumulation, snowfall and rainfall rates

   accum(j,i) = precip(j,i,0)

   snowfall(j,i) = 0.0_dp   ! initialisation value

   do n=1, 12   ! month counter

#if (SOLID_PRECIP == 1)     /* Marsiat (1994) */

      if (temp_mm(j,i,n) >= temp_rain) then
         frac_solid = 0.0_dp
      else if (temp_mm(j,i,n) <= temp_snow) then
         frac_solid = 1.0_dp
      else
         frac_solid = (temp_rain-temp_mm(j,i,n))*inv_delta_temp_rain_snow
      end if

#elif (SOLID_PRECIP == 2)   /* Bales et al. (2009) */

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

#elif (SOLID_PRECIP == 3)   /* Huybrechts and de Wolde (1999) */

      frac_solid = 1.0_dp &
                   - 0.5_dp*erfc((temp_rain-temp_mm(j,i,n))*inv_sqrt2_s_stat)

#endif

      snowfall(j,i) = snowfall(j,i) + precip(j,i,n)*frac_solid*inv_twelve

   end do

   rainfall(j,i) = precip(j,i,0) - snowfall(j,i)

   if (snowfall(j,i) < 0.0_dp) snowfall(j,i) = 0.0_dp   ! correction of
   if (rainfall(j,i) < 0.0_dp) rainfall(j,i) = 0.0_dp   ! negative values

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

#endif

end do
end do

!  ------ SMB = precipitation minus runoff

as_perp = accum - runoff

!  ------ Ice-surface temperature (10-m firn temperature) temp_s,
!         including empirical firn-warming correction due to
!         refreezing meltwater when superimposed ice is formed

where (melt_star >= melt)
   temp_s = temp_ma + mu*(melt_star-melt)
elsewhere
   temp_s = temp_ma
end where

where (temp_s > -0.001_dp) temp_s = -0.001_dp
                            ! Cut-off of positive air temperatures

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
