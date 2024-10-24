!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ b a s _ m e l t _ m
!
!! Computation of the basal melting rate.
!!
!!##### Authors
!!
!! Ralf Greve, Ben Galton-Fenzi, Tatsuru Sato
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
!> Computation of the basal melting rate.
!-------------------------------------------------------------------------------
module calc_bas_melt_m

  use sico_types_m
  use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

  use error_m

  implicit none

  private
  public :: calc_qbm

contains

!-------------------------------------------------------------------------------
!> Computation of the basal melting rate Q_bm.
!! Summation of Q_bm and Q_tld (water drainage rate from the temperate layer).
!-------------------------------------------------------------------------------
subroutine calc_qbm(time, dzeta_c, dzeta_r)

use ice_material_properties_m, only : kappa_val

implicit none

real(dp), intent(in) :: time
real(dp), intent(in) :: dzeta_c, dzeta_r

integer(i4b) :: i, j
integer(i4b) :: n_year_CE
integer(i4b) :: n_ocean, n_float
real(dp) :: time_in_years
real(dp) :: rhow_rho_ratio
real(dp) :: aqbm1, aqbm2, aqbm3a, aqbm3b, aqbm4
real(dp) :: z_abyssal
real(dp) :: frictional_heating
real(dp) :: Q_bm_grounded, Q_bm_marine, Q_bm_floating
real(dp) :: qbm_min, qbm_max
real(dp) :: H_w_now, Q_bm_scaling_factor

!-------- Term abbreviations --------

time_in_years = time*sec2year
n_year_CE     = floor((time_in_years+YEAR_ZERO)+eps_sp_dp)

rhow_rho_ratio = RHO_W/RHO

if (flag_aa_nonzero) then
   aqbm1 = (ea-1.0_dp)/(RHO*L*aa*dzeta_c)
else
   aqbm1 = 1.0_dp/(RHO*L*dzeta_c)
end if

aqbm2  = KAPPA_R/(RHO*L*H_R*dzeta_r)
aqbm3a = G/L
aqbm3b = 1.0_dp/(RHO*L)
aqbm4  = BETA/(RHO*L)

!-------- Computation of Q_bm --------

Q_bm = 0.0_dp   ! initialisation

#if (defined(Z_ABYSS))
z_abyssal = Z_ABYSS
#else
z_abyssal = no_value_neg_1
#endif

do i=1, IMAX-1
do j=1, JMAX-1

   if (mask(j,i)==0) then   ! grounded ice

      if (n_cts(j,i)==-1) then

         frictional_heating = 0.0_dp
         Q_bm(j,i)          = 0.0_dp

      else if (n_cts(j,i)==0) then

#if (DYNAMICS==2 || DYNAMICS==3)
         if (.not.flag_shelfy_stream(j,i)) then
#endif
            frictional_heating &
               = -aqbm3a*H_c(j,i)*0.5_dp &
                  * ( (vx_t(0,j,i)+vx_t(0,j,i-1))*dzs_dxi_g(j,i) &
                     +(vy_t(0,j,i)+vy_t(0,j-1,i))*dzs_deta_g(j,i) )
#if (DYNAMICS==2 || DYNAMICS==3)
         else   ! flag_shelfy_stream(j,i) == .true.

            frictional_heating &
               = aqbm3b &
                    * c_drag(j,i) &
                    * sqrt(vx_b_g(j,i)**2  &
                          +vy_b_g(j,i)**2) &
                                    **(1.0_dp+p_weert_inv(j,i))
         end if
#endif
         Q_bm(j,i) =   aqbm1*(temp_c(1,j,i)-temp_c(0,j,i))/H_c(j,i) &
                            *kappa_val(temp_c(0,j,i)) &
                     - aqbm2*(temp_r(KRMAX,j,i)-temp_r(KRMAX-1,j,i)) &
                     + frictional_heating

      else   ! n_cts(j,i)==1

#if (DYNAMICS==2 || DYNAMICS==3)
         if (.not.flag_shelfy_stream(j,i)) then
#endif
            frictional_heating &
               = -aqbm3a*H(j,i)*0.5_dp &
                 * ( (vx_t(0,j,i)+vx_t(0,j,i-1))*dzs_dxi_g(j,i) &
                    +(vy_t(0,j,i)+vy_t(0,j-1,i))*dzs_deta_g(j,i) )
#if (DYNAMICS==2 || DYNAMICS==3)
         else   ! flag_shelfy_stream(j,i) == .true.

            frictional_heating &
               = aqbm3b &
                    * c_drag(j,i) &
                    * sqrt(vx_b_g(j,i)**2  &
                          +vy_b_g(j,i)**2) &
                                    **(1.0_dp+p_weert_inv(j,i))
         end if
#endif
         Q_bm(j,i) =   aqbm4*kappa_val(temp_t_m(0,j,i)) &
                     - aqbm2*(temp_r(KRMAX,j,i)-temp_r(KRMAX-1,j,i)) &
                     + frictional_heating
           
      end if

#if (MARGIN==1 || MARGIN==2)

#if (!defined(MARINE_ICE_BASAL_MELTING) || MARINE_ICE_BASAL_MELTING==1)

      !%% continue

#elif (MARINE_ICE_BASAL_MELTING==2)

      if ( (zb(j,i) < z_sl(j,i)) &          ! marine ice
           .and. &
           (     (mask(j,i+1)==2) &   ! at least one
             .or.(mask(j,i-1)==2) &   ! nearest neighbour
             .or.(mask(j+1,i)==2) &   ! is
             .or.(mask(j-1,i)==2) &   ! ocean
           ) &
         ) then

         Q_bm(j,i) = QBM_MARINE *sec2year*rhow_rho_ratio
                                ! m/a water equiv. -> m/s ice equiv.

      end if

#elif (MARINE_ICE_BASAL_MELTING==3)

      if ( (zb(j,i) < z_sl(j,i)) &          ! marine ice
           .and. &
           (     (mask(j,i+1)==2) &   ! at least one
             .or.(mask(j,i-1)==2) &   ! nearest neighbour
             .or.(mask(j+1,i)==2) &   ! is
             .or.(mask(j-1,i)==2) &   ! ocean
           ) &
         ) then

         n_ocean = 0
         if (mask(j,i+1)==2) n_ocean = n_ocean+1
         if (mask(j,i-1)==2) n_ocean = n_ocean+1
         if (mask(j+1,i)==2) n_ocean = n_ocean+1
         if (mask(j-1,i)==2) n_ocean = n_ocean+1

         if ( n_ocean > 0 ) then

            Q_bm_grounded = Q_bm(j,i)
            Q_bm_marine   = QBM_MARINE *sec2year*rhow_rho_ratio
                                       ! m/a water equiv. -> m/s ice equiv.

            Q_bm(j,i) = (1.0_dp-0.25_dp*real(n_ocean,dp)) * Q_bm_grounded &
                               +0.25_dp*real(n_ocean,dp)  * Q_bm_marine
                      ! weighed average of grounded ice melting (computed)
                      ! and marine ice melting (prescribed)
         else
            errormsg = ' >>> calc_qbm: Marine ice margin point does not' &
                     //                end_of_line &
                     //'               have an ocean neighbour!'
            call error(errormsg)
         end if

      end if

#else
      errormsg = ' >>> calc_qbm: MARINE_ICE_BASAL_MELTING must be 1, 2 or 3!'
      call error(errormsg)
#endif

#elif (MARGIN==3)

      if (flag_grounding_line_1(j,i)) then
                                ! grounding line (grounded-ice side)
         !%% continue

      else if ( (zb(j,i) < z_sl(j,i)) &          ! marine ice margin
                .and. &
                (     (mask(j,i+1)>=2) &   !  (at least one
                  .or.(mask(j,i-1)>=2) &   !   nearest neighbour
                  .or.(mask(j+1,i)>=2) &   !   is
                  .or.(mask(j-1,i)>=2) &   !   ocean)
                ) &
              ) then

#if (MARINE_ICE_BASAL_MELTING==1)

         !%% continue

#elif (MARINE_ICE_BASAL_MELTING==2)

         Q_bm(j,i) = QBM_MARINE *sec2year*rhow_rho_ratio
                                ! m/a water equiv. -> m/s ice equiv.

#elif (MARINE_ICE_BASAL_MELTING==3)

         n_ocean = 0
         if (mask(j,i+1)>=2) n_ocean = n_ocean+1
         if (mask(j,i-1)>=2) n_ocean = n_ocean+1
         if (mask(j+1,i)>=2) n_ocean = n_ocean+1
         if (mask(j-1,i)>=2) n_ocean = n_ocean+1

         if ( n_ocean > 0 ) then

            Q_bm_grounded = Q_bm(j,i)
            Q_bm_marine   = QBM_MARINE *sec2year*rhow_rho_ratio
                                       ! m/a water equiv. -> m/s ice equiv.

            Q_bm(j,i) = (1.0_dp-0.25_dp*real(n_ocean,dp)) * Q_bm_grounded &
                               +0.25_dp*real(n_ocean,dp)  * Q_bm_marine
                      ! weighed average of grounded ice melting (computed)
                      ! and marine ice melting (prescribed)
         else
            errormsg = ' >>> calc_qbm: Marine ice margin point does not' &
                     //                end_of_line &
                     //'               have a floating ice or ocean neighbour!'
            call error(errormsg)
         end if

#else
         errormsg = ' >>> calc_qbm: MARINE_ICE_BASAL_MELTING must be 1, 2 or 3!'
         call error(errormsg)
#endif

      end if

   else if ( (mask(j,i)==2).or.(mask(j,i)==3) ) then
                                                ! floating ice or ocean

#if (FLOATING_ICE_BASAL_MELTING==1)

      if ( zl(j,i) > z_abyssal ) then   ! floating ice over continental shelf

         Q_bm(j,i) = QBM_FLOAT_1 *sec2year*rhow_rho_ratio
                                 ! m/a water equiv. -> m/s ice equiv.

      else   ! ( zl(j,i) <= z_abyssal ), floating ice over abyssal ocean

         Q_bm(j,i) = QBM_FLOAT_3 *sec2year*rhow_rho_ratio
                                 ! m/a water equiv. -> m/s ice equiv.

      end if

#elif (FLOATING_ICE_BASAL_MELTING==4 || FLOATING_ICE_BASAL_MELTING==5)

      if ( zl(j,i) > z_abyssal ) then   ! floating ice over continental shelf
         call sub_ice_shelf_melting_param_1(time, time_in_years, &
                                            rhow_rho_ratio, &
                                            i, j, Q_bm_floating)
         Q_bm(j,i) = Q_bm_floating
      else   ! ( zl(j,i) <= z_abyssal ), floating ice over abyssal ocean
         Q_bm(j,i) = QBM_FLOAT_3 *sec2year*rhow_rho_ratio
                                 ! m/a water equiv. -> m/s ice equiv.
      end if

#elif (FLOATING_ICE_BASAL_MELTING==6)

      !%% continue
      !%% (will be computed below by subroutine sub_ice_shelf_melting_param_2)

#else
      errormsg = ' >>> calc_qbm: FLOATING_ICE_BASAL_MELTING' &
               //                end_of_line &
               //'               must be 1, 4, 5 or 6!'
      call error(errormsg)
#endif

#else
      errormsg = ' >>> calc_qbm: MARGIN must be 1, 2 or 3!'
      call error(errormsg)
#endif

   end if

end do
end do

#if (FLOATING_ICE_BASAL_MELTING==6)

call sub_ice_shelf_melting_param_2(time, time_in_years, &
                                   rhow_rho_ratio, z_abyssal, &
                                   n_year_CE)

#endif

!-------- Sum of Q_bm and Q_tld --------

Q_b_tot = Q_bm + Q_tld

!-------- Limitation of Q_bm, Q_tld and Q_b_tot --------

#if (defined(QBM_MIN))
   qbm_min = QBM_MIN *sec2year*rhow_rho_ratio
                      ! m/a water equiv. -> m/s ice equiv.
#elif (defined(QB_MIN)) /* obsolete */
   qbm_min = QB_MIN   ! m/s ice equiv.
#else
   errormsg = ' >>> calc_qbm: Limiter QBM_MIN is undefined!'
   call error(errormsg)
#endif

#if (defined(QBM_MAX))
   qbm_max = QBM_MAX *sec2year*rhow_rho_ratio
                      ! m/a water equiv. -> m/s ice equiv.
#elif (defined(QB_MAX)) /* obsolete */
   qbm_max = QB_MAX   ! m/s ice equiv.
#else
   errormsg = ' >>> calc_qbm: Limiter QBM_MAX is undefined!'
   call error(errormsg)
#endif

#if ( defined(MARINE_ICE_BASAL_MELTING) && ( MARINE_ICE_BASAL_MELTING==2 || MARINE_ICE_BASAL_MELTING==3 ) )

if (QBM_MARINE*sec2year*rhow_rho_ratio > qbm_max) then
   errormsg = ' >>> calc_qbm: QBM_MARINE' &
            //                end_of_line &
            //'               (basal melting rate at the ice front)' &
            //                end_of_line &
            //'               is larger than the limiter qbm_max!'
   call error(errormsg)
end if

#endif

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i)==0) then   ! grounded ice

      if (Q_bm(j,i)    < qbm_min) Q_bm(j,i)    = 0.0_dp
      if (Q_bm(j,i)    > qbm_max) Q_bm(j,i)    = qbm_max
      if (Q_tld(j,i)   < qbm_min) Q_tld(j,i)   = 0.0_dp
      if (Q_tld(j,i)   > qbm_max) Q_tld(j,i)   = qbm_max
      if (Q_b_tot(j,i) < qbm_min) Q_b_tot(j,i) = 0.0_dp
      if (Q_b_tot(j,i) > qbm_max) Q_b_tot(j,i) = qbm_max

   else if ( (mask(j,i)==3).and.( zl(j,i) > z_abyssal ) ) then
                                 ! floating ice over continental shelf
#if (defined(H_W_0))
      if (H_W_0 > eps_dp) then
         H_w_now             = max((zb(j,i)-zl(j,i)), 0.0_dp)
         Q_bm_scaling_factor = tanh(H_w_now/H_W_0)
         Q_bm(j,i)           = Q_bm(j,i)*Q_bm_scaling_factor
         Q_b_tot(j,i)        = Q_bm(j,i) + Q_tld(j,i)
      end if
#endif

   end if

end do
end do

end subroutine calc_qbm

!-------------------------------------------------------------------------------
!> Local sub-ice-shelf melting parameterization.
!-------------------------------------------------------------------------------
subroutine sub_ice_shelf_melting_param_1(time, time_in_years, &
                                         rhow_rho_ratio, &
                                         i, j, Q_bm_floating)

#if defined(ALLOW_TAPENADE) /* Tapenade */
  use ctrl_m, only: myfloor
#endif /* Tapenade */

implicit none

integer(i4b), intent(in)    :: i, j
real(dp), intent(in) :: time
real(dp), intent(in) :: time_in_years
real(dp), intent(in) :: rhow_rho_ratio

real(dp), intent(out) :: Q_bm_floating

integer(i4b) :: n
real(dp) :: Toc, Tmb, T_forcing, Omega, alpha, draft0, draft
real(dp) :: lon_d, lat_d, Phi_par

real(dp), parameter :: beta_sw = 7.61e-04_dp

#if defined(ALLOW_TAPENADE) /* Tapenade */
integer(i4b) :: i_time_in_years
#endif /* Tapenade */

#if (FLOATING_ICE_BASAL_MELTING==5)

real(dp), parameter :: c_sw = 3974.0_dp, &
                       g_t  =    5.0e-05_dp

#endif

#if defined(ALLOW_TAPENADE) /* Tapenade */
call myfloor(time_in_years, i_time_in_years)
#endif /* Tapenade */

#if (FLOATING_ICE_BASAL_MELTING==4)

Toc   = TEMP_OCEAN
Omega = OMEGA_QBM *sec2year*rhow_rho_ratio
                  ! m/[a*degC^alpha] water equiv. -> m/[s*degC^alpha] ice equiv.
alpha = ALPHA_QBM

#elif (FLOATING_ICE_BASAL_MELTING==5)

#if (defined(ANT))   /* Antarctic ice sheet */

!-------- Definition of the sectors --------

lon_d = lambda(j,i) *rad2deg
lat_d = phi(j,i)    *rad2deg

#if !defined(ALLOW_TAPENADE) /* Normal */
lon_d = modulo(lon_d+180.0_dp, 360.0_dp)-180.0_dp
                                  ! mapping to interval
                                  ! [-180 deg, +180 deg)
#else /* Tapenade */
lon_d = (lon_d+180.0_sp) &
           - ((360.0_sp * int(lon_d+180.0_sp)/360.0_sp)) - 180.0_sp
#endif /* Normal vs. Tapenade */

if ((lon_d>=-10.0_dp).and.(lon_d<60.0_dp)) then
                                      ! Western East Antarctica
   n_bm_region(j,i) = 1

   Toc    =   0.23_dp
   Omega  =   0.0076_dp
   alpha  =   1.33_dp
   draft0 = 200.0_dp

else if ((lon_d>=60.0_dp).and.(lon_d<80.0_dp)) then
                                      ! Amery/Prydz Bay
   n_bm_region(j,i) = 2

   Toc    =  -1.61_dp
   Omega  =   0.0128_dp
   alpha  =   0.94_dp
   draft0 = 200.0_dp

else if ((lon_d>=80.0_dp).and.(lon_d<130.0_dp)) then
                                      ! Sabrina Coast/
                                      ! Aurora subglacial basin
   n_bm_region(j,i) = 3

   Toc    =  -0.16_dp
   Omega  =   0.0497_dp
   alpha  =   0.85_dp
   draft0 = 200.0_dp

else if ( ((lon_d>=130.0_dp).and.(lon_d<159.0_dp)) &
          .or. &
          ((lon_d>=159.0_dp).and.(lon_d<170.0_dp) &
                            .and.(lat_d>=-72.0_dp)) ) then
                                      ! George V Coast/
                                      ! Wilkes subglacial basin
   n_bm_region(j,i) = 4

   Toc    =  -0.22_dp
   Omega  =   0.0423_dp
   alpha  =   0.11_dp
   draft0 = 200.0_dp

else if ( (lon_d>=159.0_dp) &
          .or. &
          (lon_d<-140.0_dp) &
          .or. &
          ((lon_d>=-140.0_dp).and.(lon_d<-120.0_dp) &
                             .and.(lat_d<-77.0_dp)) ) then
                                      ! Ross Sea
   n_bm_region(j,i) = 5

   Toc    =  -1.8_dp
   Omega  =   0.0181_dp
   alpha  =   0.73_dp
   draft0 = 200.0_dp

else if ( ((lon_d>=-140.0_dp).and.(lon_d<-120.0_dp) &
                             .and.(lat_d>=-77.0_dp)) &
          .or. &
          ((lon_d>=-120.0_dp).and.(lon_d<-90.0_dp)) ) then
                                      ! Amundsen Sea
   n_bm_region(j,i) = 6

   Toc    =   0.42_dp
   Omega  =   0.1023_dp
   alpha  =   0.26_dp
   draft0 = 200.0_dp

else if ( ((lon_d>=-90.0_dp).and.(lon_d<-66.0_dp) &
                            .and.(lat_d>=-74.0_dp)) ) then
                                      ! Bellingshausen Sea
   n_bm_region(j,i) = 7

   Toc    =   0.62_dp
   Omega  =   0.0644_dp
   alpha  =   0.16_dp
   draft0 = 200.0_dp

else
                                      ! Weddell Sea
   n_bm_region(j,i) = 8

   Toc    =  -1.55_dp
   Omega  =   0.0197_dp
   alpha  =   0.26_dp
   draft0 = 200.0_dp

endif

#else   /* not Antarctic ice sheet */

Q_bm_floating = 0.0_dp   ! dummy return value

errormsg = ' >>> sub_ice_shelf_melting_param_1:' &
         //          end_of_line &
         //'         Case FLOATING_ICE_BASAL_MELTING==5' &
         //          end_of_line &
         //'         only defined for Antarctica!'
call error(errormsg)

#endif

#else

!%% continue

#endif

!-------- Computation of the sub-ice-shelf melting rate --------

if (mask(j,i)==2) then   ! ocean
   draft = 0.0_dp
else if (mask(j,i)==3) then   ! floating ice
   draft = max((z_sl(j,i)-zb(j,i)), 0.0_dp)
else
   errormsg = ' >>> sub_ice_shelf_melting_param_1:' &
            //          end_of_line &
            //'         Routine must not be called for mask(j,i) < 2!'
   call error(errormsg)
end if

Tmb = -beta_sw*draft - DELTA_TM_SW   ! temperature at the ice shelf base

T_forcing = max((Toc-Tmb), 0.0_dp)   ! thermal forcing

#if (FLOATING_ICE_BASAL_MELTING==4)

Q_bm_floating  = Omega*T_forcing**alpha

#elif (FLOATING_ICE_BASAL_MELTING==5)

Phi_par = RHO_SW*c_sw*g_t/(RHO*L)

#if !defined(ALLOW_TAPENADE)
Q_bm_floating = Phi_par*Omega*T_forcing*(draft/draft0)**alpha
#else
if (draft.eq.0) then
Q_bm_floating = 0.0_dp 
else
Q_bm_floating = Phi_par*Omega*T_forcing*(draft/draft0)**alpha
end if
#endif
#else

Q_bm_floating = 0.0_dp   ! dummy return value

errormsg = ' >>> sub_ice_shelf_melting_param_1:' &
         //          end_of_line &
         //'         FLOATING_ICE_BASAL_MELTING must be 4 or 5!'
call error(errormsg)

#endif

!  ------ Correction for ISMIP6 InitMIP

#if (defined(ANT)) /* Antarctic ice sheet */

if (flag_initmip_abmb) then

   if ((time_in_years > 0.0_dp).and.(time_in_years <= 40.0_dp)) then

#if !defined(ALLOW_TAPENADE) /* Normal */
      Q_bm_floating = Q_bm_floating &
                      + 0.025_dp*floor(time_in_years)*ab_anom_initmip(j,i)
#else /* Tapenade */
      Q_bm_floating = Q_bm_floating &
                      + 0.025_dp*i_time_in_years*ab_anom_initmip(j,i)
#endif /* Normal vs. Tapenade */

   else if (time_in_years > 40.0_dp) then

      Q_bm_floating = Q_bm_floating + ab_anom_initmip(j,i)

   end if

end if

#endif

!  ------ Correction for ISMIP6 LARMIP

#if (defined(ANT)) /* Antarctic ice sheet */

if (flag_larmip) then
   n             = n_larmip_region(j,i)
   Q_bm_floating = Q_bm_floating + ab_anom_larmip(n)
end if

#endif

end subroutine sub_ice_shelf_melting_param_1

!-------------------------------------------------------------------------------
!> Non-local sub-ice-shelf melting parameterization by ISMIP6.
!-------------------------------------------------------------------------------
subroutine sub_ice_shelf_melting_param_2(time, time_in_years, &
                                         rhow_rho_ratio, z_abyssal, &
                                         n_year_CE)

  use netcdf
  use nc_check_m

#if !defined(ALLOW_TAPENADE) /* Normal */
  use compare_float_m
#endif /* Normal */

#if defined(ALLOW_TAPENADE) /* Tapenade */
  use ctrl_m, only: myfloor
#endif /* Tapenade */

implicit none

real(dp)    , intent(in) :: time
real(dp)    , intent(in) :: time_in_years
real(dp)    , intent(in) :: rhow_rho_ratio, z_abyssal
integer(i4b), intent(in) :: n_year_CE

#if (defined(ANT) && FLOATING_ICE_BASAL_MELTING==6)

integer(i4b) :: i, j, n

integer(i4b)       :: ios
integer(i4b)       :: n_year_CE_bas_melt
character(len= 16) :: ch_year_CE
character(len=256) :: filename_with_path
real(dp), dimension(0:IMAX,0:JMAX,0:NZ_TF_BM) :: tf_bm_aux

character(len=64), parameter :: thisroutine = 'sub_ice_shelf_melting_param_2'

integer(i4b) :: ncid
!     ncid:      File ID
integer(i4b) :: ncv
!     ncv:       Variable ID

integer(i4b) :: ninf, nsup, n_bm_regions
real(dp)     :: dz_inv, tf_bm_no_value_neg
real(dp)     :: cst_bm
real(dp)     :: weigh
real(dp)     :: draft
real(dp)     :: real_n

#if defined(ALLOW_TAPENADE) /* Tapenade */
integer(i4b) :: i_time_in_years
#endif /* Tapenade */

#if (!defined(N_BM_REGIONS) || N_BM_REGIONS<=1)
real(dp) :: gamma0_bm_aux(1)
real(dp) :: delta_tf_bm_aux(1)
real(dp) :: tf_bm_ave(1)
real(dp) :: sum_weigh(1)
#else
real(dp) :: gamma0_bm_aux(N_BM_REGIONS)
real(dp) :: delta_tf_bm_aux(N_BM_REGIONS)
real(dp) :: tf_bm_ave(N_BM_REGIONS)
real(dp) :: sum_weigh(N_BM_REGIONS)
#endif

real(dp), dimension(0:JMAX,0:IMAX) :: gamma0_bm, delta_tf_bm, tf_bm_local

real(dp), parameter :: rhoi_bm  = 918.0_dp
                                  ! Ice density (kg/m3)
real(dp), parameter :: rhosw_bm = 1028.0_dp
                                  ! Sea water density (kg/m3)
real(dp), parameter :: rhofw_bm = 1000.0_dp
                                  ! Fresh water density (kg/m3)
real(dp), parameter :: Lf_bm    = 3.34e+05_dp
                                  ! Latent heat of fusion for ice (J/kg)
real(dp), parameter :: cpw_bm   = 3974.0_dp
                                  ! Specific heat of sea water (J/kg/K)

#if defined(ALLOW_TAPENADE) /* Tapenade */
call myfloor(time_in_years, i_time_in_years)
#endif /* Tapenade */

!-------- Read file with the thermal forcing data of the ocean --------

n_year_CE_bas_melt = n_year_CE

if (n_year_CE_bas_melt < TF_BM_TIME_MIN) then
   n_year_CE_bas_melt = TF_BM_TIME_MIN
else if (n_year_CE_bas_melt > TF_BM_TIME_MAX) then
   n_year_CE_bas_melt = TF_BM_TIME_MAX
end if

if ( firstcall%sub_ice_shelf_melting_param_2 &
     .or.(n_year_CE_bas_melt /= n_year_CE_bas_melt_save) ) then

   write(ch_year_CE, '(i0)') n_year_CE_bas_melt

   if ( (trim(adjustl(TF_BM_FILES)) /= 'none') &
        .and. &
        (trim(adjustl(TF_BM_FILES)) /= 'None') &
        .and. &
        (trim(adjustl(TF_BM_FILES)) /= 'NONE') ) then

!  ------ Read data from file

      filename_with_path = trim(TF_BM_DIR)//'/'// &
                           trim(TF_BM_FILES)//trim(ch_year_CE)//'.nc'

      ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

      if (ios /= nf90_noerr) then
         errormsg = ' >>> sub_ice_shelf_melting_param_2:' &
                  //                end_of_line &
                  //'               Error when opening the file' &
                  //                end_of_line &
                  //'               for the thermal forcing of the ocean!'
         call error(errormsg)
      end if

      call check( nf90_inq_varid(ncid, 'z', ncv) )
      call check( nf90_get_var(ncid, ncv, z_tf_bm) )

      call check( nf90_inq_varid(ncid, 'thermal_forcing', ncv), thisroutine )
      call check( nf90_get_var(ncid, ncv, tf_bm_aux), thisroutine )

      call check( nf90_close(ncid), thisroutine )

!  ------ Ensure positive depth values

      if ( (z_tf_bm(0) < eps_dp).and.(z_tf_bm(NZ_TF_BM) < eps_dp) ) &
         z_tf_bm = -z_tf_bm

!  ------ Swap indices -> SICOPOLIS standard

      do i=0, IMAX
      do j=0, JMAX
      do n=0, NZ_TF_BM
         if (isnan(tf_bm_aux(i,j,n))) then
            tf_bm(n,j,i) = no_value_neg_2
         else
            tf_bm(n,j,i) = tf_bm_aux(i,j,n)
         end if
      end do
      end do
      end do

!  ------ Check consistency of the depth (z_tf_bm) data

#if !defined(ALLOW_TAPENADE) /* Normal */

      if (.not.(approx_equal(z_tf_bm(0), ZMIN_TF_BM, eps_sp_dp))) then
         errormsg = ' >>> sub_ice_shelf_melting_param_2:' &
                  //         end_of_line &
                  //'        Inconsistency between' &
                  //         end_of_line &
                  //'        read z_tf_bm data' &
                  //         end_of_line &
                  //'        and parameter ZMIN_TF_BM!'
         call error(errormsg)
      end if

      if (.not.(approx_equal(z_tf_bm(NZ_TF_BM)-z_tf_bm(0), &
                             NZ_TF_BM*DZ_TF_BM, eps_sp_dp))) then
         errormsg = ' >>> sub_ice_shelf_melting_param_2:' &
                  //         end_of_line &
                  //'        Inconsistency between' &
                  //         end_of_line &
                  //'        read z_tf_bm data' &
                  //         end_of_line &
                  //'        and parameters NZ_TF_BM, DZ_TF_BM!'
         call error(errormsg)
      end if

#endif /* Normal */

   else   ! ( trim(adjustl(TF_BM_FILES)) == 'none' or 'None' or 'NONE' )

!  ------ Use present-day data

      tf_bm_aux = 0.0_dp
      tf_bm     = tf_bm_present
      z_tf_bm   = z_tf_bm_present

   end if

end if

!  ------ Save value of n_year_CE_bas_melt

n_year_CE_bas_melt_save = n_year_CE_bas_melt

!-------- Parameters for the parameterization --------

cst_bm = ((rhosw_bm*cpw_bm)/(rhoi_bm*Lf_bm))**2   ! (1/K2)

#if (!defined(N_BM_REGIONS) || N_BM_REGIONS<=1)
n_bm_regions = 1
#else
n_bm_regions = N_BM_REGIONS
#endif

#if (defined(GAMMA0_BM))
  gamma0_bm_aux = GAMMA0_BM
  gamma0_bm_aux = gamma0_bm_aux *sec2year*(rhofw_bm/rhoi_bm)
                                ! m/a water equiv. -> m/s ice equiv.
#else
  errormsg = ' >>> sub_ice_shelf_melting_param_2: GAMMA0_BM must be defined!'
  call error(errormsg)
#endif

#if (defined(DELTA_TF_BM))
  delta_tf_bm_aux = DELTA_TF_BM
#else
  errormsg = ' >>> sub_ice_shelf_melting_param_2: DELTA_TF_BM must be defined!'
  call error(errormsg)
#endif

do i=0, IMAX
do j=0, JMAX

   n = n_bm_region(j,i)

   if ( (n >= 1).and.(n <= n_bm_regions) ) then
      gamma0_bm(j,i)   = gamma0_bm_aux(n)
      delta_tf_bm(j,i) = delta_tf_bm_aux(n)
   else
      errormsg = ' >>> sub_ice_shelf_melting_param_2: ' &
                    //'Region number out of allowed range!'
      call error(errormsg)
   end if

end do
end do

!-------- Local thermal forcing --------

dz_inv = 1.0_dp/DZ_TF_BM

tf_bm_no_value_neg = 0.999_dp*no_value_neg_2

tf_bm_local = 1.0_dp   ! default value

do i=0, IMAX
do j=0, JMAX

   if ( (mask(j,i)==2).or.(mask(j,i)==3) ) then
                                                ! floating ice or ocean

      if ( zl(j,i) > z_abyssal ) then   ! continental shelf

         if (mask(j,i)==2) then   ! ocean
            draft = 0.0_dp
         else if (mask(j,i)==3) then   ! floating ice
            draft = max((z_sl(j,i)-zb(j,i)), 0.0_dp)
         end if

         real_n = (draft-ZMIN_TF_BM)*dz_inv
         ninf   = floor(real_n)
         nsup   = ninf+1

         if ((ninf >= 0).and.(nsup <= NZ_TF_BM)) then

            if ( &
                 (tf_bm(ninf,j,i) > tf_bm_no_value_neg) &
                 .and. &
                 (tf_bm(nsup,j,i) > tf_bm_no_value_neg) &
               ) &
               tf_bm_local(j,i) = tf_bm(ninf,j,i) &
                                    + ((draft-z_tf_bm(ninf))*dz_inv) &
                                      *(tf_bm(nsup,j,i)-tf_bm(ninf,j,i))

         else if (ninf < 0) then

            ninf = 0

            if (tf_bm(ninf,j,i) > tf_bm_no_value_neg) &
               tf_bm_local(j,i) = tf_bm(ninf,j,i)

         else   ! nsup > NZ_TF_BM

            nsup = NZ_TF_BM

            if (tf_bm(nsup,j,i) > tf_bm_no_value_neg) &
               tf_bm_local(j,i) = tf_bm(nsup,j,i)

         end if

      end if

   end if

end do
end do

!-------- Region-averaged thermal forcing --------

tf_bm_ave = 0.0_dp   ! initialization
sum_weigh = 0.0_dp   ! initialization

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i)==3) then   ! floating ice

      n = n_bm_region(j,i)

      weigh        = cell_area(j,i)
      tf_bm_ave(n) = tf_bm_ave(n) + weigh*tf_bm_local(j,i)
      sum_weigh(n) = sum_weigh(n) + weigh

   end if

end do
end do

do n=1, n_bm_regions

   if (sum_weigh(n) > eps_dp) then
      tf_bm_ave(n) = tf_bm_ave(n)/sum_weigh(n)
   else
      tf_bm_ave(n) = 1.0_dp   ! default value
   end if

end do

!-------- Computation of the sub-ice-shelf melting rate --------

do i=0, IMAX
do j=0, JMAX

   if ( (mask(j,i)==2).or.(mask(j,i)==3) ) then
                                                ! floating ice or ocean

!  ------ Melting rate over continental shelf

      if ( zl(j,i) > z_abyssal ) then

         n = n_bm_region(j,i)

         Q_bm(j,i) = gamma0_bm(j,i)*cst_bm &
                        *(tf_bm_local(j,i)+delta_tf_bm(j,i)) &
                        *abs(tf_bm_ave(n)+delta_tf_bm(j,i))

         if ((mask(j,i)==2).and.(Q_bm(j,i) < 0.0_dp)) Q_bm(j,i) = 0.0_dp
                                     ! avoid negative values for the open ocean

!    ---- Correction for ISMIP6 InitMIP

#if (defined(ANT)) /* Antarctic ice sheet */

         if (flag_initmip_abmb) then

            if ((time_in_years > 0.0_dp).and.(time_in_years <= 40.0_dp)) then

#if !defined(ALLOW_TAPENADE) /* Normal */
               Q_bm(j,i) = Q_bm(j,i) &
                           + 0.025_dp*floor(time_in_years)*ab_anom_initmip(j,i)
#else /* Tapenade */
               Q_bm(j,i) = Q_bm(j,i) &
                           + 0.025_dp*i_time_in_years*ab_anom_initmip(j,i)
#endif /* Normal vs. Tapenade */

            else if (time_in_years > 40.0_dp) then

               Q_bm(j,i) = Q_bm(j,i) + ab_anom_initmip(j,i)

            end if

         end if

#endif

!    ---- Correction for ISMIP6 LARMIP

#if (defined(ANT)) /* Antarctic ice sheet */

         if (flag_larmip) then
            n         = n_larmip_region(j,i)
            Q_bm(j,i) = Q_bm(j,i) + ab_anom_larmip(n)
         end if

#endif

!  ------ Melting rate over abyssal ocean

      else   ! ( zl(j,i) <= z_abyssal )

         Q_bm(j,i) = QBM_FLOAT_3 *sec2year*rhow_rho_ratio
                                 ! m/a water equiv. -> m/s ice equiv.

      end if

   end if

end do
end do

if (firstcall%sub_ice_shelf_melting_param_2) &
    firstcall%sub_ice_shelf_melting_param_2 = .false.

#else   /* not (defined(ANT) && FLOATING_ICE_BASAL_MELTING==6) */

errormsg = ' >>> sub_ice_shelf_melting_param_2:' &
         //          end_of_line &
         //'         Routine only valid for Antarctica' &
         //          end_of_line &
         //'         and FLOATING_ICE_BASAL_MELTING==6!'
call error(errormsg)

#endif

end subroutine sub_ice_shelf_melting_param_2

!-------------------------------------------------------------------------------

end module calc_bas_melt_m
!
