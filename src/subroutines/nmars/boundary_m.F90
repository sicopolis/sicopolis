!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  b o u n d a r y _ m
!
!! NMARS domain:
!! Mars Atmosphere-Ice Coupler MAIC-1.5.
!! Computation of the surface temperature (must be less than 0 degC)
!! and of the accumulation-ablation rate for the north polar cap of Mars.
!! Computation of the geothermal heat flux.
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
!> NMARS domain:
!! Mars Atmosphere-Ice Coupler MAIC-1.5.
!! Computation of the surface temperature (must be less than 0 degC)
!! and of the accumulation-ablation rate for the north polar cap of Mars.
!! Computation of the geothermal heat flux.
!-------------------------------------------------------------------------------
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
!! Mars Atmosphere-Ice Coupler MAIC-1.5:
!! Computation of the surface temperature (must be less than 0 degC)
!! and of the accumulation-ablation rate for the north polar cap of Mars.
!! Computation of the geothermal heat flux.
!-------------------------------------------------------------------------------
subroutine boundary(time, dtime, dxi, deta, &
                    delta_ts, glac_index, z_mar)

  use mask_update_sea_level_m

  use output_m, only : borehole

#if (TSURFACE==6)
  use mars_instemp_m
#endif

#if ((MARGIN==2) && (MARINE_ICE_FORMATION==2) && (MARINE_ICE_CALVING==9))
  use calving_m
#endif

implicit none

real(dp), intent(in) :: time, dtime, dxi, deta

real(dp), intent(out)   :: delta_ts, glac_index, z_mar

! Further return variables
! (defined as global variables in module sico_variables_m):
!
!    accum(j,i), runoff(j,i), as_perp(j,i), calving(j,i), temp_s(j,i)

integer(i4b) :: i, j
integer(i4b) :: i_gr, i_kl
integer(i4b) :: ndata_insol
real(dp), dimension(0:JMAX,0:IMAX) :: z_sl_old
real(dp) :: z_sl_old_mean
real(dp) :: z_sl_min, t1, t2, t3, t4, t5, t6
real(dp) :: time_gr, time_kl
real(dp) :: z_sle_present, z_sle_help
real(dp), dimension(0:JMAX,0:IMAX) :: temp_ma
real(dp) :: eld, g_mb, accum_factor
real(dp) :: erosion_chasm
real(dp) :: t_obliq_main, t_obliq_mod, &
            obliq_ampl_max, obliq_ampl_min, obliq_ampl, obliq0, obliq, &
            ecc, ecc0
real(dp) :: insol_ma_90_now, insol_ma_90_present, &
            obl_now, obl_present, ecc_now, ecc_present, &
            ave_now, ave_present, cp_now, cp_present, &
            temp_ma_90, temp_ma_90_present, zs_90
real(dp), dimension(0:JMAX,0:IMAX) :: dist
real(dp) :: ave_data_i_kl, ave_data_i_gr
real(dp) :: q_geo_chasm
logical                           :: flag_in_domain
logical, dimension(0:JMAX,0:IMAX) :: check_point
logical, save                     :: firstcall = .true.

#if (TSURFACE==6)
type (ins) :: temp_now, temp_present
#endif

real(dp), parameter :: &
          time_present  = 0.0_dp, &     ! Present time [s]
          zs_90_present = -2.0e+03_dp   ! Present elevation of the
                                        ! north pole [m]
real(dp), parameter :: &
          sol     = 590.0_dp, &   ! Solar constant [W/m2]
          epsilon = 1.0_dp,   &   ! Emissivity
          sigma   = 5.67e-08_dp   ! Stefan-Boltzmann constant [W/(m2*K4)]

real(dp), parameter :: &
          lambda_H2O = 2.86e+06_dp, &   ! Latent heat [J/kg]
          R_H2O      = 461.5_dp,  &     ! Gas constant [J/(kg*K)]
          temp0_ref  = 173.0_dp         ! Atmopheric reference temperature [K]

real(dp), parameter :: &
          temp_s_min = -125.0_dp   ! Minimum ice-surface temperature 
                                   ! (sublimation temperature of CO2) [C]

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
delta_ts = -SINE_AMPLIT &
           *cos(2.0_dp*pi*time/(SINE_PERIOD*year2sec)) &
           +SINE_AMPLIT
!                           ! Sinusoidal air-temperature forcing
#elif (TSURFACE==4)

!  ------ Obliquity (main cycle and first modulation) and eccentricity

t_obliq_main = 1.25e+05_dp *year2sec
t_obliq_mod  = 1.3e+06_dp  *year2sec

obliq0         = 25.2_dp *deg2rad
obliq_ampl_max = 10.0_dp *deg2rad
obliq_ampl_min =  2.5_dp *deg2rad

obliq_ampl = 0.5_dp*(obliq_ampl_max+obliq_ampl_min) &
            -0.5_dp*(obliq_ampl_max-obliq_ampl_min) &
                  *cos(2.0_dp*pi*time/t_obliq_mod)

obliq      = obliq0 + obliq_ampl*sin(2.0_dp*pi*time/t_obliq_main)

ecc  = 0.0_dp   ! Values not correct, but influence
ecc0 = 0.0_dp   ! on mean-annual north-polar insolation is negligible

!  ------ Mean annual insolation at the north pole

insol_ma_90_now     = (sol/pi)*sin(obliq)/sqrt(1.0_dp-ecc**2)

insol_ma_90_present = (sol/pi)*sin(obliq0)/sqrt(1.0_dp-ecc0**2)
                       ! present value                       

#elif (TSURFACE==5 || TSURFACE==6)

!  ------ Mean annual insolation at the north pole

ndata_insol = (insol_time_max-insol_time_min)/insol_time_stp

if (time*sec2year.lt.real(insol_time_min,dp)) then

   insol_ma_90_now = insol_ma_90(0)
   obl_now         = obl_data(0)
   ecc_now         = ecc_data(0)
   ave_now         = ave_data(0)
   cp_now          = cp_data(0)

else if (time*sec2year.lt.real(insol_time_max,dp)) then

   i_kl = floor(((time*sec2year) &
          -real(insol_time_min,dp))/real(insol_time_stp,dp))
   i_kl = max(i_kl, 0)

   i_gr = ceiling(((time*sec2year) &
          -real(insol_time_min,dp))/real(insol_time_stp,dp))
   i_gr = min(i_gr, ndata_insol)

   if (i_kl.eq.i_gr) then

      insol_ma_90_now = insol_ma_90(i_kl)
      obl_now         = obl_data(i_kl)
      ecc_now         = ecc_data(i_kl)
      ave_now         = ave_data(i_kl)
      cp_now          = cp_data(i_kl)

   else

      time_kl = (insol_time_min + i_kl*insol_time_stp) *year2sec
      time_gr = (insol_time_min + i_gr*insol_time_stp) *year2sec

      insol_ma_90_now = insol_ma_90(i_kl) &
                +(insol_ma_90(i_gr)-insol_ma_90(i_kl)) &
                *(time-time_kl)/(time_gr-time_kl)
      obl_now = obl_data(i_kl) &
                +(obl_data(i_gr)-obl_data(i_kl)) &
                *(time-time_kl)/(time_gr-time_kl)
      ecc_now = ecc_data(i_kl) &
                +(ecc_data(i_gr)-ecc_data(i_kl)) &
                *(time-time_kl)/(time_gr-time_kl)

      if ( abs(ave_data(i_gr)-ave_data(i_kl)) < pi ) then   ! regular case
            ave_data_i_kl = ave_data(i_kl)
            ave_data_i_gr = ave_data(i_gr)
      else
         if ( ave_data(i_gr) > ave_data(i_kl) ) then
            ave_data_i_kl = ave_data(i_kl) + 2.0_dp*pi
            ave_data_i_gr = ave_data(i_gr)
         else
            ave_data_i_kl = ave_data(i_kl)
            ave_data_i_gr = ave_data(i_gr) + 2.0_dp*pi
         end if
      end if

      ave_now = ave_data_i_kl &
                +(ave_data_i_gr-ave_data_i_kl) &
                *(time-time_kl)/(time_gr-time_kl)

      cp_now  = cp_data(i_kl) &
                +(cp_data(i_gr)-cp_data(i_kl)) &
                *(time-time_kl)/(time_gr-time_kl)
                 ! linear interpolation of the data
   end if

else

   insol_ma_90_now = insol_ma_90(ndata_insol)
   obl_now         = obl_data(ndata_insol)
   ecc_now         = ecc_data(ndata_insol)
   ave_now         = ave_data(ndata_insol)
   cp_now          = cp_data(ndata_insol)

end if

!    ---- Present value

if (time_present*sec2year.lt.real(insol_time_min,dp)) then

   insol_ma_90_present = insol_ma_90(0)
   obl_present         = obl_data(0)
   ecc_present         = ecc_data(0)
   ave_present         = ave_data(0)
   cp_present          = cp_data(0)

else if (time_present*sec2year.lt.real(insol_time_max,dp)) then

   i_kl = floor(((time_present*sec2year) &
          -real(insol_time_min,dp))/real(insol_time_stp,dp))
   i_kl = max(i_kl, 0)

   i_gr = ceiling(((time_present*sec2year) &
          -real(insol_time_min,dp))/real(insol_time_stp,dp))
   i_gr = min(i_gr, ndata_insol)

   if (i_kl.eq.i_gr) then

      insol_ma_90_present = insol_ma_90(i_kl)
      obl_present         = obl_data(i_kl)
      ecc_present         = ecc_data(i_kl)
      ave_present         = ave_data(i_kl)
      cp_present          = cp_data(i_kl)

   else

      time_kl = (insol_time_min + i_kl*insol_time_stp) *year2sec
      time_gr = (insol_time_min + i_gr*insol_time_stp) *year2sec

      insol_ma_90_present = insol_ma_90(i_kl) &
                +(insol_ma_90(i_gr)-insol_ma_90(i_kl)) &
                *(time_present-time_kl)/(time_gr-time_kl)
      obl_present = obl_data(i_kl) &
                +(obl_data(i_gr)-obl_data(i_kl)) &
                *(time_present-time_kl)/(time_gr-time_kl)
      ecc_present = ecc_data(i_kl) &
                +(ecc_data(i_gr)-ecc_data(i_kl)) &
                *(time_present-time_kl)/(time_gr-time_kl)

      if ( abs(ave_data(i_gr)-ave_data(i_kl)) < pi ) then   ! regular case
            ave_data_i_kl = ave_data(i_kl)
            ave_data_i_gr = ave_data(i_gr)
      else
         if ( ave_data(i_gr) > ave_data(i_kl) ) then
            ave_data_i_kl = ave_data(i_kl) + 2.0_dp*pi
            ave_data_i_gr = ave_data(i_gr)
         else
            ave_data_i_kl = ave_data(i_kl)
            ave_data_i_gr = ave_data(i_gr) + 2.0_dp*pi
         end if
      end if

      ave_present = ave_data_i_kl &
                +(ave_data_i_gr-ave_data_i_kl) &
                *(time_present-time_kl)/(time_gr-time_kl)

      cp_present  = cp_data(i_kl) &
                +(cp_data(i_gr)-cp_data(i_kl)) &
                *(time_present-time_kl)/(time_gr-time_kl)
                 ! linear interpolation of the data
   end if

else

   insol_ma_90_present = insol_ma_90(ndata_insol)
   obl_present         = obl_data(ndata_insol)
   ecc_present         = ecc_data(ndata_insol)
   ave_present         = ave_data(ndata_insol)
   cp_present          = cp_data(ndata_insol)

end if

#endif

#if (TSURFACE==4 || TSURFACE==5)

!  ------ Mean-annual surface temperature at the north pole

temp_ma_90 = sqrt(sqrt(insol_ma_90_now*(1.0_dp-ALBEDO)/(epsilon*sigma))) &
             -273.15_dp   ! K -> C

!    ---- Present value

temp_ma_90_present &
        = sqrt(sqrt(insol_ma_90_present*(1.0_dp-ALBEDO)/(epsilon*sigma))) &
          -273.15_dp   ! K -> C

!  ------ Surface-temperature deviation at the north pole

delta_ts = temp_ma_90 - temp_ma_90_present

#elif (TSURFACE==6)

!  ------ Mean-annual surface temperature at the north pole

call setinstemp(temp_now, &
                ecc = ecc_now, ave = ave_now*rad2deg, &
                obl = obl_now*rad2deg, sa = ALBEDO, ct = 148.7_dp)

temp_ma_90 = instam(temp_now, 90.0_dp) + (DELTA_TS0) &
                                       - 273.15_dp   ! K -> C

!    ---- Present value

call setinstemp(temp_present, &
                ecc = ecc_present, ave = ave_present*rad2deg, &
                obl = obl_present*rad2deg, sa = ALBEDO, ct = 148.7_dp)

temp_ma_90_present = instam(temp_present, 90.0_dp) + (DELTA_TS0) &
                                                   - 273.15_dp   ! K -> C

!  ------ Surface-temperature deviation at the north pole

delta_ts = temp_ma_90 - temp_ma_90_present

#endif

!-------- Sea level --------

#if (SEA_LEVEL==1)

!  ------ Temporally constant sea level

z_sl = Z_SL0

#else

errormsg = ' >>> boundary: Parameter SEA_LEVEL must be 1!'
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

!-------- Surface air temperature !

call borehole(zs, 0.0_dp, 0.0_dp, dxi, deta, 'grid', &
              zs_90, flag_in_domain)
                             ! zs_90: elevation of the north pole

do i=0, IMAX
do j=0, JMAX

!  ------ Present mean annual air temperature temp_ma

#if (TSURFACE==1 || TSURFACE==2 || TSURFACE==3)

   temp_ma_90_present = TEMP0_MA_90N

   temp_ma(j,i) = temp_ma_90_present &
                  + (GAMMA_MA)*(zs(j,i)-zs_90_present) &
                  + (C_MA*rad2deg)*((90.0_dp*deg2rad)-phi(j,i))

#elif (TSURFACE==4 || TSURFACE==5)

   temp_ma(j,i) = temp_ma_90_present &
                  + (GAMMA_MA)*(zs(j,i)-zs_90_present) &
                  + (C_MA*rad2deg)*((90.0_dp*deg2rad)-phi(j,i))

#endif

!  ------ Correction with deviation delta_ts

#if (TSURFACE==1 || TSURFACE==2 || TSURFACE==3)

   temp_ma(j,i) = temp_ma(j,i) + delta_ts

#elif (TSURFACE==4 || TSURFACE==5)

   temp_ma(j,i) = temp_ma(j,i) + delta_ts &
                               - (GAMMA_MA)*(zs_90-zs_90_present)

#endif

!  ------ Direct computation of the mean annual air temperature temp_ma

#if (TSURFACE==6)

   temp_ma(j,i) = instam(temp_now, phi(j,i)*rad2deg) + (DELTA_TS0) &
                                                        - 273.15_dp   ! K -> C

#endif

end do
end do

!  ------ Save mean-annual air temperature

temp_maat = temp_ma

!-------- Accumulation-ablation function as_perp --------

!  ------ Accumulation

#if (ACCSURFACE==1)
   accum_factor = 1.0_dp+GAMMA_S*delta_ts
#elif (ACCSURFACE==2)
   accum_factor = exp(GAMMA_S*delta_ts)
#elif (ACCSURFACE==3)
   accum_factor = exp(lambda_H2O/(R_H2O*temp0_ref) &
                     -lambda_H2O/(R_H2O*(temp0_ref+delta_ts)))
#endif

do i=0, IMAX
do j=0, JMAX
   accum(j,i) = accum_present(j,i)*accum_factor
   if (accum(j,i).lt.0.0_dp) then
      errormsg = ' >>> boundary: Negative accumulation rate!'
      call error(errormsg)
   end if
end do
end do

!  ------ Ablation

#if (ABLSURFACE==1 || ABLSURFACE==2)

eld  = ELD_0 *1.0e+03_dp   ! km -> m

#if (ACC_UNIT==1)
g_mb = G_0 *(1.0e-06_dp*sec2year)*(RHO_W/RHO_I) &
           *(1.0_dp/(1.0_dp-FRAC_DUST))
           ! [mm/a water equiv.]/km -> [m/s (ice+dust) equiv.]/m
#elif (ACC_UNIT==2)
g_mb = G_0 *(1.0e-06_dp*sec2year)
           ! [mm/a (ice+dust) equiv.]/km -> [m/s (ice+dust) equiv.]/m
#endif

#if (ABLSURFACE==1)
eld  = eld*(1.0_dp+GAMMA_ELD*delta_ts)
g_mb = g_mb*(GAMMA_G*accum_factor)
#elif (ABLSURFACE==2)
eld  = eld*exp(GAMMA_ELD*delta_ts)
g_mb = g_mb*(GAMMA_G*accum_factor)
#endif

if (eld.lt.eps) then
   errormsg = ' >>> boundary: Negative equilibrium line distance eld!'
   call error(errormsg)
end if

do i=0, IMAX
do j=0, JMAX
   dist(j,i) = sqrt( xi(i)**2 + eta(j)**2 )
end do
end do

#endif

#if (CHASM==2)

#if (ACC_UNIT==1)
erosion_chasm = EROSION_CHASM *(1.0e-03_dp*sec2year)*(RHO_W/RHO_I) &
                              *(1.0_dp/(1.0_dp-FRAC_DUST))
           ! [mm/a water equiv.] -> [m/s (ice+dust) equiv.]
#elif (ACC_UNIT==2)
erosion_chasm = EROSION_CHASM  *(1.0e-03_dp*sec2year)
           ! [mm/a (ice+dust) equiv.] -> [m/s (ice+dust) equiv.]
#endif

#endif

do i=0, IMAX
do j=0, JMAX

#if (ABLSURFACE==1 || ABLSURFACE==2)

   as_perp(j,i) = g_mb*(eld-dist(j,i))
   if (as_perp(j,i).gt.accum(j,i)) as_perp(j,i) = accum(j,i)

#endif

   runoff(j,i) = accum(j,i) - as_perp(j,i)

#if (CHASM==2)

!    ---- Correction for additional wind erosion in the chasm area

   if ( (mask_chasm(j,i) == 7) &
        .and.(time >= time_chasm_init) &
        .and.(time <= time_chasm_end) ) then   ! active chasm area
      runoff(j,i)  = runoff(j,i) + erosion_chasm
      as_perp(j,i) = accum(j,i) - runoff(j,i)
   end if

#endif

end do
end do

!-------- Ice-surface temperature (10-m firn temperature) temp_s --------

temp_s = min(temp_ma, -eps)        ! Cut-off of positive air temperatures

temp_s = max(temp_s, temp_s_min)   ! Cut-off of air temperatures below the
                                   ! sublimation temperature of CO2

!-------- Calving --------

calving = 0.0_dp   ! Initialization

#if ((MARGIN==2) && (MARINE_ICE_FORMATION==2) && (MARINE_ICE_CALVING==9))

call calving_underwater_ice()

#endif

!-------- Geothermal heat flux --------

#if (CHASM==1)

q_geo = q_geo_normal

#elif (CHASM==2)

q_geo_chasm = Q_GEO_CHASM *1.0e-03_dp   ! mW/m2 -> W/m2

do i=0, IMAX
do j=0, JMAX

   if ( (mask_chasm(j,i) == 7) &
        .and.(time >= time_chasm_init) &
        .and.(time <= time_chasm_end) ) then   ! active chasm area
      q_geo(j,i) = q_geo_chasm
   else
      q_geo(j,i) = q_geo_normal(j,i)
   end if

end do
end do

#endif

if (firstcall) firstcall = .false.

end subroutine boundary

!-------------------------------------------------------------------------------

end module boundary_m
!
