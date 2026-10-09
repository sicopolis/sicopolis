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
z_abyssal = r_no_value_neg_1
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

#elif (FLOATING_ICE_BASAL_MELTING>=6 && FLOATING_ICE_BASAL_MELTING<=8)

      !%% continue
      !%% (will be computed below by subroutine sub_ice_shelf_melting_param_2)

#else
      errormsg = ' >>> calc_qbm: FLOATING_ICE_BASAL_MELTING' &
               //                end_of_line &
               //'               must be 1, 4, 5, 6, 7 or 8!'
      call error(errormsg)
#endif

#else
      errormsg = ' >>> calc_qbm: MARGIN must be 1, 2 or 3!'
      call error(errormsg)
#endif

   end if

end do
end do

#if (FLOATING_ICE_BASAL_MELTING>=6 && FLOATING_ICE_BASAL_MELTING<=8)

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

!  ------ Correction for ISMIP6/7 LARMIP

#if (defined(ANT)) /* Antarctic ice sheet */

if (flag_larmip) then
   n             = n_larmip_region(j,i)
   Q_bm_floating = Q_bm_floating + ab_anom_larmip(n)
end if

#endif

end subroutine sub_ice_shelf_melting_param_1

!-------------------------------------------------------------------------------
!> Local or semi-local sub-ice-shelf melting parameterization by ISMIP6/7
!! (FLOATING_ICE_BASAL_MELTING==6/7),
!! or the local sector-wise power relationships
!! by Galton-Fenzi (2026, Sci Adv, submitted) on the same thermal forcing
!! (FLOATING_ICE_BASAL_MELTING==8).
!!
!! FLOATING_ICE_BASAL_MELTING==8 (added October 2026):
!! sector-wise power relationships by Galton-Fenzi (2026, Sci Adv, submitted),
!!   Q_bm = lambda*Gamma*TF^alpha,
!! fitted to the WAOM4 and RISE multi-model-mean ocean products in eight
!! Antarctic sectors. The thermal forcing TF is the ISMIP6-style field read
!! by sub_ice_shelf_melting_param_2 (shelf water extrapolated to the draft);
!! the sector follows from the grid longitude/latitude. Optional macros:
!!   BM_SECTOR_PRODUCT  1 = WAOM4 (default), 2 = RISE MMM
!!   BM_SECTOR_APPROACH 1 = thermal driving at the ice base, sector melt-weighted
!!                          exponent alpha' (recommended where that driving is known)
!!                      2 = continental-shelf water at the draft, sector exponent
!!                          alpha_cs of Table 1 with its recommendations
!!                      3 = continental-shelf water at the draft, global alpha_cs with
!!                          the coefficient of the sector (default: the ISMIP6-style
!!                          forcing read here, the relationship the paper recommends)
!!                      Every option has a coefficient per sector.
!!   BM_SECTOR_CALIB    1 = coefficient set on the fitted cells (lambda); 2 = on the whole
!!                          base to the net mass loss (lambda_w; default for options 2 and 3,
!!                          which predict melt wherever the shelf driving is positive)
!!   BM_SECTOR_GLOBAL   1 = one global relationship everywhere (default 0)
!!   DELTA_TF_SECTOR    uniform offset of the thermal forcing [degC] (default 0)
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

#if (defined(ANT) && (FLOATING_ICE_BASAL_MELTING>=6 && FLOATING_ICE_BASAL_MELTING<=8))

integer(i4b) :: i, j, n

integer(i4b)       :: ios
integer(i4b)       :: istat, istat1, istat2
integer(i4b)       :: n_year_CE_tf_bm
character(len= 16) :: ch_year_CE
character(len=256) :: filename_with_path
logical            :: flag_param_semilocal
real(dp), dimension(0:IMAX,0:JMAX,0:NZ_TF_BM) :: tf_bm_aux

#if (FLOATING_ICE_BASAL_MELTING==7)
integer(i4b)       :: n_year_CE_so_bm
real(dp), dimension(0:IMAX,0:JMAX,0:NZ_SO_BM) :: so_bm_aux
#endif

character(len=64), parameter :: thisroutine = 'sub_ice_shelf_melting_param_2'

integer(i4b) :: ncid
!     ncid:      File ID
integer(i4b) :: ncv
!     ncv:       Variable ID

integer(i4b) :: ninf, nsup, n_bm_regions
real(dp)     :: dz_tf_inv, tf_bm_no_value_neg
real(dp)     :: slope_mean
real(dp)     :: cst_bm
real(dp)     :: weigh
real(dp)     :: draft
real(dp)     :: real_n
real(dp)     :: tf_bm_help

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

#if (FLOATING_ICE_BASAL_MELTING==7)
real(dp) :: dz_so_inv
real(dp), dimension(0:JMAX,0:IMAX) :: so_bm_local
#endif

#if (FLOATING_ICE_BASAL_MELTING==6)
real(dp), parameter :: rhoi_bm = 918.0_dp
                                  ! Ice density (kg/m3)
#elif (FLOATING_ICE_BASAL_MELTING==7)
real(dp), parameter :: rhoi_bm = 917.0_dp
                                  ! Ice density (kg/m3)
#endif

real(dp), parameter :: rhosw_bm = 1028.0_dp
                                  ! Sea water density (kg/m3)
real(dp), parameter :: rhofw_bm = 1000.0_dp
                                  ! Fresh water density (kg/m3)
real(dp), parameter :: Lf_bm    = 3.34e+05_dp
                                  ! Latent heat of fusion for ice (J/kg)
real(dp), parameter :: cpw_bm   = 3974.0_dp
                                  ! Specific heat of sea water (J/kg/K)

#if (FLOATING_ICE_BASAL_MELTING==7)
real(dp), parameter :: beta_s_bm = 7.86e-04_dp
                                   ! Salt contraction ocefficient (1/PSU)
real(dp), parameter :: coriolis_par_bm = -1.4e-04_dp
                                         ! Coriolis parameter (1/s)
real(dp), parameter :: so_bm_default = 34.7_dp
                                   ! Default value of the oceanic salinity (PSU)
#endif

#if (FLOATING_ICE_BASAL_MELTING==8)
integer(i4b) :: n_sector, i_product, i_approach, i_calib
real(dp) :: lam_sector_use
real(dp)     :: lon_d, lat_d, tf_sector, delta_tf_sector
! >>> BEGIN sector melt coefficients
! (auto-generated by ana/export_fortran_coeffs.m -- do not edit)
! Sector-wise power relationships by Galton-Fenzi (2026, Sci Adv, submitted),
! fitted on the annual mean of WAOM4 (product 1, recommended)
! and on the RISE multi-model mean (product 2):
!    m = lambda * Gamma * TF**alpha     [m/s, ice equivalent]; TF = thermal driving [degC]
!    alpha  = the exponent in use; Gamma = coefficient fitted in log space at that exponent on the
!    sector's cells; lambda = retransformation correction, so that lambda*Gamma reproduces the
!    sector's simulated mass loss. The number that multiplies TF**alpha is lambda*Gamma.
! Options (second index of the sec_* arrays; the named arrays below spell each one out):
!    1 interface-mw  TF = thermal driving at the ice base; sector melt-weighted exponent alpha'
!                    (the exponent that governs the integrated mass loss)
!    2 shelf-sector  TF = continental-shelf thermal driving at the draft (ISMIP6-style thermal
!                    forcing); sector exponent alpha_cs of Table 1 with the rules below
!    3 shelf-global  as 2 with the global alpha_cs and the coefficient of the sector (the relationship
!                    the paper recommends where only the continental-shelf driving is known)
! Two corrections per option: lambda, set on the fitted cells (lambda*Gamma reproduces the melt of the
! sector on the cells where every driving is positive), and lambda_w, set on the whole ice-shelf base
! with the sector's NET mass loss (melt minus refreezing) as the target, for a model that applies
! the relationship on every cell with a positive driving. The continental-shelf driving is positive
! almost everywhere, so with lambda the shelf options over-predict where the cavity refreezes; use
! lambda_w with them unless the refreezing area is masked. The ice-base option predicts nothing
! where TF <= 0 and its lambda_w equals lambda.
! Recommendation of the paper: option 1 wherever the thermal driving at the ice base is known
! (coupled ice-ocean models); otherwise option 3, the global alpha_cs with the coefficient of the
! sector and lambda_w. Option 2 is kept for comparison. Each option carries a coefficient per sector.
! Recommendations of Table 1: a sector fit not identified (sector R^2 < 0.20 x global, or a negative exponent)
! defaults to the linear relationship in a cold sector and to the global exponent in the warm
! Amundsen and Bellingshausen sectors; a sector exponent above 2 is replaced by the global one; each
! sector keeps its own coefficient. Sectors (first index): 0 = global relationship, then
!    1 Dronning Maud Land
!    2 Amery-Prydz Bay
!    3 Sabrina-Aurora
!    4 George V-Wilkes
!    5 Ross Sea
!    6 Amundsen Sea
!    7 Bellingshausen Sea
!    8 Weddell Sea
! ---- WAOM4, option 1 (interface-mw) ----
  real(kind=dp), parameter :: waom4_interface_mw_alpha(0:8) = (/ &
       1.204793_dp     , &   ! 0 Global               (global)
       1.119486_dp     , &   ! 1 Dronning Maud Land   (sector)
       1.273379_dp     , &   ! 2 Amery-Prydz Bay      (sector)
       1.244757_dp     , &   ! 3 Sabrina-Aurora       (sector)
       1.164304_dp     , &   ! 4 George V-Wilkes      (sector)
       1.093666_dp     , &   ! 5 Ross Sea             (sector)
       1.207782_dp     , &   ! 6 Amundsen Sea         (sector)
       1.094358_dp     , &   ! 7 Bellingshausen Sea   (sector)
       1.149885_dp      /)   ! 8 Weddell Sea          (sector)
  real(kind=dp), parameter :: waom4_interface_mw_Gamma(0:8) = (/ &
       6.913375e-07_dp , &   ! 0 Global
       4.198805e-07_dp , &   ! 1 Dronning Maud Land
       8.864637e-07_dp , &   ! 2 Amery-Prydz Bay
       7.616341e-07_dp , &   ! 3 Sabrina-Aurora
       5.960857e-07_dp , &   ! 4 George V-Wilkes
       6.246149e-07_dp , &   ! 5 Ross Sea
       5.099717e-07_dp , &   ! 6 Amundsen Sea
       3.444989e-07_dp , &   ! 7 Bellingshausen Sea
       5.886092e-07_dp  /)   ! 8 Weddell Sea
  real(kind=dp), parameter :: waom4_interface_mw_lambda(0:8) = (/ &
       0.846912_dp     , &   ! 0 Global               lambda*Gamma = 5.8550e-07
       1.132528_dp     , &   ! 1 Dronning Maud Land   lambda*Gamma = 4.7553e-07
       0.824267_dp     , &   ! 2 Amery-Prydz Bay      lambda*Gamma = 7.3068e-07
       1.061908_dp     , &   ! 3 Sabrina-Aurora       lambda*Gamma = 8.0879e-07
       1.129884_dp     , &   ! 4 George V-Wilkes      lambda*Gamma = 6.7351e-07
       0.883052_dp     , &   ! 5 Ross Sea             lambda*Gamma = 5.5157e-07
       1.077939_dp     , &   ! 6 Amundsen Sea         lambda*Gamma = 5.4972e-07
       1.145895_dp     , &   ! 7 Bellingshausen Sea   lambda*Gamma = 3.9476e-07
       0.807999_dp      /)   ! 8 Weddell Sea          lambda*Gamma = 4.7560e-07
  real(kind=dp), parameter :: waom4_interface_mw_lambdaw(0:8) = (/ &
       0.846912_dp     , &   ! 0 Global               whole base: lambda_w*Gamma = 5.8550e-07
       1.132528_dp     , &   ! 1 Dronning Maud Land   whole base: lambda_w*Gamma = 4.7553e-07
       0.824267_dp     , &   ! 2 Amery-Prydz Bay      whole base: lambda_w*Gamma = 7.3068e-07
       1.061908_dp     , &   ! 3 Sabrina-Aurora       whole base: lambda_w*Gamma = 8.0879e-07
       1.129884_dp     , &   ! 4 George V-Wilkes      whole base: lambda_w*Gamma = 6.7351e-07
       0.883052_dp     , &   ! 5 Ross Sea             whole base: lambda_w*Gamma = 5.5157e-07
       1.077939_dp     , &   ! 6 Amundsen Sea         whole base: lambda_w*Gamma = 5.4972e-07
       1.145895_dp     , &   ! 7 Bellingshausen Sea   whole base: lambda_w*Gamma = 3.9476e-07
       0.807999_dp      /)   ! 8 Weddell Sea          whole base: lambda_w*Gamma = 4.7560e-07
! ---- WAOM4, option 2 (shelf-sector) ----
  real(kind=dp), parameter :: waom4_shelf_sector_alpha(0:8) = (/ &
       1.326176_dp     , &   ! 0 Global               (global)
       1.000000_dp     , &   ! 1 Dronning Maud Land   (linear)
       1.295501_dp     , &   ! 2 Amery-Prydz Bay      (sector)
       1.000000_dp     , &   ! 3 Sabrina-Aurora       (linear)
       0.776736_dp     , &   ! 4 George V-Wilkes      (sector)
       0.713885_dp     , &   ! 5 Ross Sea             (sector)
       1.326176_dp     , &   ! 6 Amundsen Sea         (global)
       0.932748_dp     , &   ! 7 Bellingshausen Sea   (sector)
       0.971173_dp      /)   ! 8 Weddell Sea          (sector)
  real(kind=dp), parameter :: waom4_shelf_sector_Gamma(0:8) = (/ &
       2.874823e-08_dp , &   ! 0 Global
       3.412576e-08_dp , &   ! 1 Dronning Maud Land
       2.632288e-08_dp , &   ! 2 Amery-Prydz Bay
       4.116090e-08_dp , &   ! 3 Sabrina-Aurora
       8.014067e-08_dp , &   ! 4 George V-Wilkes
       1.032885e-08_dp , &   ! 5 Ross Sea
       3.766895e-08_dp , &   ! 6 Amundsen Sea
       4.109245e-08_dp , &   ! 7 Bellingshausen Sea
       1.811728e-08_dp  /)   ! 8 Weddell Sea
  real(kind=dp), parameter :: waom4_shelf_sector_lambda(0:8) = (/ &
       1.967524_dp     , &   ! 0 Global               lambda*Gamma = 5.6563e-08
       1.682877_dp     , &   ! 1 Dronning Maud Land   lambda*Gamma = 5.7429e-08
       1.777221_dp     , &   ! 2 Amery-Prydz Bay      lambda*Gamma = 4.6782e-08
       2.004172_dp     , &   ! 3 Sabrina-Aurora       lambda*Gamma = 8.2493e-08
       1.561021_dp     , &   ! 4 George V-Wilkes      lambda*Gamma = 1.2510e-07
       4.661218_dp     , &   ! 5 Ross Sea             lambda*Gamma = 4.8145e-08
       1.626200_dp     , &   ! 6 Amundsen Sea         lambda*Gamma = 6.1257e-08
       1.661934_dp     , &   ! 7 Bellingshausen Sea   lambda*Gamma = 6.8293e-08
       1.799278_dp      /)   ! 8 Weddell Sea          lambda*Gamma = 3.2598e-08
  real(kind=dp), parameter :: waom4_shelf_sector_lambdaw(0:8) = (/ &
       1.696172_dp     , &   ! 0 Global               whole base: lambda_w*Gamma = 4.8762e-08
       1.670544_dp     , &   ! 1 Dronning Maud Land   whole base: lambda_w*Gamma = 5.7009e-08
       1.205560_dp     , &   ! 2 Amery-Prydz Bay      whole base: lambda_w*Gamma = 3.1734e-08
       1.787927_dp     , &   ! 3 Sabrina-Aurora       whole base: lambda_w*Gamma = 7.3593e-08
       1.543120_dp     , &   ! 4 George V-Wilkes      whole base: lambda_w*Gamma = 1.2367e-07
       2.432838_dp     , &   ! 5 Ross Sea             whole base: lambda_w*Gamma = 2.5128e-08
       1.588245_dp     , &   ! 6 Amundsen Sea         whole base: lambda_w*Gamma = 5.9828e-08
       1.594763_dp     , &   ! 7 Bellingshausen Sea   whole base: lambda_w*Gamma = 6.5533e-08
       1.134482_dp      /)   ! 8 Weddell Sea          whole base: lambda_w*Gamma = 2.0554e-08
! ---- WAOM4, option 3 (shelf-global) ----
  real(kind=dp), parameter :: waom4_shelf_global_alpha(0:8) = (/ &
       1.326176_dp     , &   ! 0 Global               (global)
       1.326176_dp     , &   ! 1 Dronning Maud Land   (global)
       1.326176_dp     , &   ! 2 Amery-Prydz Bay      (global)
       1.326176_dp     , &   ! 3 Sabrina-Aurora       (global)
       1.326176_dp     , &   ! 4 George V-Wilkes      (global)
       1.326176_dp     , &   ! 5 Ross Sea             (global)
       1.326176_dp     , &   ! 6 Amundsen Sea         (global)
       1.326176_dp     , &   ! 7 Bellingshausen Sea   (global)
       1.326176_dp      /)   ! 8 Weddell Sea          (global)
  real(kind=dp), parameter :: waom4_shelf_global_Gamma(0:8) = (/ &
       2.874823e-08_dp , &   ! 0 Global
       3.379460e-08_dp , &   ! 1 Dronning Maud Land
       2.691158e-08_dp , &   ! 2 Amery-Prydz Bay
       4.339429e-08_dp , &   ! 3 Sabrina-Aurora
       8.745674e-08_dp , &   ! 4 George V-Wilkes
       2.053243e-08_dp , &   ! 5 Ross Sea
       3.766895e-08_dp , &   ! 6 Amundsen Sea
       5.843699e-08_dp , &   ! 7 Bellingshausen Sea
       2.209787e-08_dp  /)   ! 8 Weddell Sea
  real(kind=dp), parameter :: waom4_shelf_global_lambda(0:8) = (/ &
       1.967524_dp     , &   ! 0 Global               lambda*Gamma = 5.6563e-08
       1.604411_dp     , &   ! 1 Dronning Maud Land   lambda*Gamma = 5.4220e-08
       1.768129_dp     , &   ! 2 Amery-Prydz Bay      lambda*Gamma = 4.7583e-08
       1.878156_dp     , &   ! 3 Sabrina-Aurora       lambda*Gamma = 8.1501e-08
       1.434884_dp     , &   ! 4 George V-Wilkes      lambda*Gamma = 1.2549e-07
       3.450508_dp     , &   ! 5 Ross Sea             lambda*Gamma = 7.0847e-08
       1.626200_dp     , &   ! 6 Amundsen Sea         lambda*Gamma = 6.1257e-08
       1.508579_dp     , &   ! 7 Bellingshausen Sea   lambda*Gamma = 8.8157e-08
       1.495040_dp      /)   ! 8 Weddell Sea          lambda*Gamma = 3.3037e-08
  real(kind=dp), parameter :: waom4_shelf_global_lambdaw(0:8) = (/ &
       1.552875_dp     , &   ! 0 Global               whole base: lambda_w*Gamma = 4.4642e-08
       1.593479_dp     , &   ! 1 Dronning Maud Land   whole base: lambda_w*Gamma = 5.3851e-08
       1.201014_dp     , &   ! 2 Amery-Prydz Bay      whole base: lambda_w*Gamma = 3.2321e-08
       1.674516_dp     , &   ! 3 Sabrina-Aurora       whole base: lambda_w*Gamma = 7.2664e-08
       1.419316_dp     , &   ! 4 George V-Wilkes      whole base: lambda_w*Gamma = 1.2413e-07
       1.953503_dp     , &   ! 5 Ross Sea             whole base: lambda_w*Gamma = 4.0110e-08
       1.588245_dp     , &   ! 6 Amundsen Sea         whole base: lambda_w*Gamma = 5.9828e-08
       1.450563_dp     , &   ! 7 Bellingshausen Sea   whole base: lambda_w*Gamma = 8.4767e-08
       0.990120_dp      /)   ! 8 Weddell Sea          whole base: lambda_w*Gamma = 2.1880e-08
! ---- MMM, option 1 (interface-mw) ----
  real(kind=dp), parameter :: mmm_interface_mw_alpha(0:8) = (/ &
       1.056553_dp     , &   ! 0 Global               (global)
       0.841538_dp     , &   ! 1 Dronning Maud Land   (sector)
       0.982465_dp     , &   ! 2 Amery-Prydz Bay      (sector)
       0.946784_dp     , &   ! 3 Sabrina-Aurora       (sector)
       1.000000_dp     , &   ! 4 George V-Wilkes      (linear)
       1.239062_dp     , &   ! 5 Ross Sea             (sector)
       1.144371_dp     , &   ! 6 Amundsen Sea         (sector)
       1.235183_dp     , &   ! 7 Bellingshausen Sea   (sector)
       0.867543_dp      /)   ! 8 Weddell Sea          (sector)
  real(kind=dp), parameter :: mmm_interface_mw_Gamma(0:8) = (/ &
       1.582271e-07_dp , &   ! 0 Global
       1.111188e-07_dp , &   ! 1 Dronning Maud Land
       1.715262e-07_dp , &   ! 2 Amery-Prydz Bay
       1.381142e-07_dp , &   ! 3 Sabrina-Aurora
       1.699615e-07_dp , &   ! 4 George V-Wilkes
       2.539123e-07_dp , &   ! 5 Ross Sea
       1.773964e-07_dp , &   ! 6 Amundsen Sea
       1.764021e-07_dp , &   ! 7 Bellingshausen Sea
       1.031251e-07_dp  /)   ! 8 Weddell Sea
  real(kind=dp), parameter :: mmm_interface_mw_lambda(0:8) = (/ &
       1.103601_dp     , &   ! 0 Global               lambda*Gamma = 1.7462e-07
       1.097264_dp     , &   ! 1 Dronning Maud Land   lambda*Gamma = 1.2193e-07
       1.082761_dp     , &   ! 2 Amery-Prydz Bay      lambda*Gamma = 1.8572e-07
       1.087065_dp     , &   ! 3 Sabrina-Aurora       lambda*Gamma = 1.5014e-07
       1.041607_dp     , &   ! 4 George V-Wilkes      lambda*Gamma = 1.7703e-07
       0.969604_dp     , &   ! 5 Ross Sea             lambda*Gamma = 2.4619e-07
       1.040994_dp     , &   ! 6 Amundsen Sea         lambda*Gamma = 1.8467e-07
       1.019083_dp     , &   ! 7 Bellingshausen Sea   lambda*Gamma = 1.7977e-07
       1.108276_dp      /)   ! 8 Weddell Sea          lambda*Gamma = 1.1429e-07
  real(kind=dp), parameter :: mmm_interface_mw_lambdaw(0:8) = (/ &
       1.103601_dp     , &   ! 0 Global               whole base: lambda_w*Gamma = 1.7462e-07
       1.097264_dp     , &   ! 1 Dronning Maud Land   whole base: lambda_w*Gamma = 1.2193e-07
       1.082761_dp     , &   ! 2 Amery-Prydz Bay      whole base: lambda_w*Gamma = 1.8572e-07
       1.087065_dp     , &   ! 3 Sabrina-Aurora       whole base: lambda_w*Gamma = 1.5014e-07
       1.041607_dp     , &   ! 4 George V-Wilkes      whole base: lambda_w*Gamma = 1.7703e-07
       0.969604_dp     , &   ! 5 Ross Sea             whole base: lambda_w*Gamma = 2.4619e-07
       1.040994_dp     , &   ! 6 Amundsen Sea         whole base: lambda_w*Gamma = 1.8467e-07
       1.019083_dp     , &   ! 7 Bellingshausen Sea   whole base: lambda_w*Gamma = 1.7977e-07
       1.108276_dp      /)   ! 8 Weddell Sea          whole base: lambda_w*Gamma = 1.1429e-07
! ---- MMM, option 2 (shelf-sector) ----
  real(kind=dp), parameter :: mmm_shelf_sector_alpha(0:8) = (/ &
       1.412234_dp     , &   ! 0 Global               (global)
       0.931626_dp     , &   ! 1 Dronning Maud Land   (sector)
       1.626327_dp     , &   ! 2 Amery-Prydz Bay      (sector)
       1.842830_dp     , &   ! 3 Sabrina-Aurora       (sector)
       1.823145_dp     , &   ! 4 George V-Wilkes      (sector)
       1.000000_dp     , &   ! 5 Ross Sea             (linear)
       1.412234_dp     , &   ! 6 Amundsen Sea         (above2)
       1.412234_dp     , &   ! 7 Bellingshausen Sea   (above2)
       1.000000_dp      /)   ! 8 Weddell Sea          (linear)
  real(kind=dp), parameter :: mmm_shelf_sector_Gamma(0:8) = (/ &
       1.705908e-08_dp , &   ! 0 Global
       3.912644e-08_dp , &   ! 1 Dronning Maud Land
       6.377208e-08_dp , &   ! 2 Amery-Prydz Bay
       5.655324e-08_dp , &   ! 3 Sabrina-Aurora
       6.765658e-08_dp , &   ! 4 George V-Wilkes
       5.805170e-09_dp , &   ! 5 Ross Sea
       6.532956e-08_dp , &   ! 6 Amundsen Sea
       3.740965e-08_dp , &   ! 7 Bellingshausen Sea
       1.280205e-08_dp  /)   ! 8 Weddell Sea
  real(kind=dp), parameter :: mmm_shelf_sector_lambda(0:8) = (/ &
       1.979433_dp     , &   ! 0 Global               lambda*Gamma = 3.3767e-08
       1.162794_dp     , &   ! 1 Dronning Maud Land   lambda*Gamma = 4.5496e-08
       1.132115_dp     , &   ! 2 Amery-Prydz Bay      lambda*Gamma = 7.2197e-08
       1.180511_dp     , &   ! 3 Sabrina-Aurora       lambda*Gamma = 6.6762e-08
       1.117179_dp     , &   ! 4 George V-Wilkes      lambda*Gamma = 7.5584e-08
       2.533372_dp     , &   ! 5 Ross Sea             lambda*Gamma = 1.4707e-08
       1.198968_dp     , &   ! 6 Amundsen Sea         lambda*Gamma = 7.8328e-08
       1.242801_dp     , &   ! 7 Bellingshausen Sea   lambda*Gamma = 4.6493e-08
       1.366851_dp      /)   ! 8 Weddell Sea          lambda*Gamma = 1.7498e-08
  real(kind=dp), parameter :: mmm_shelf_sector_lambdaw(0:8) = (/ &
       1.667697_dp     , &   ! 0 Global               whole base: lambda_w*Gamma = 2.8449e-08
       1.060845_dp     , &   ! 1 Dronning Maud Land   whole base: lambda_w*Gamma = 4.1507e-08
       1.091007_dp     , &   ! 2 Amery-Prydz Bay      whole base: lambda_w*Gamma = 6.9576e-08
       1.054351_dp     , &   ! 3 Sabrina-Aurora       whole base: lambda_w*Gamma = 5.9627e-08
       0.928053_dp     , &   ! 4 George V-Wilkes      whole base: lambda_w*Gamma = 6.2789e-08
       1.805131_dp     , &   ! 5 Ross Sea             whole base: lambda_w*Gamma = 1.0479e-08
       1.139248_dp     , &   ! 6 Amundsen Sea         whole base: lambda_w*Gamma = 7.4427e-08
       1.163527_dp     , &   ! 7 Bellingshausen Sea   whole base: lambda_w*Gamma = 4.3527e-08
       0.949550_dp      /)   ! 8 Weddell Sea          whole base: lambda_w*Gamma = 1.2156e-08
! ---- MMM, option 3 (shelf-global) ----
  real(kind=dp), parameter :: mmm_shelf_global_alpha(0:8) = (/ &
       1.412234_dp     , &   ! 0 Global               (global)
       1.412234_dp     , &   ! 1 Dronning Maud Land   (global)
       1.412234_dp     , &   ! 2 Amery-Prydz Bay      (global)
       1.412234_dp     , &   ! 3 Sabrina-Aurora       (global)
       1.412234_dp     , &   ! 4 George V-Wilkes      (global)
       1.412234_dp     , &   ! 5 Ross Sea             (global)
       1.412234_dp     , &   ! 6 Amundsen Sea         (global)
       1.412234_dp     , &   ! 7 Bellingshausen Sea   (global)
       1.412234_dp      /)   ! 8 Weddell Sea          (global)
  real(kind=dp), parameter :: mmm_shelf_global_Gamma(0:8) = (/ &
       1.705908e-08_dp , &   ! 0 Global
       5.011110e-08_dp , &   ! 1 Dronning Maud Land
       6.070391e-08_dp , &   ! 2 Amery-Prydz Bay
       4.898873e-08_dp , &   ! 3 Sabrina-Aurora
       6.153472e-08_dp , &   ! 4 George V-Wilkes
       6.662865e-09_dp , &   ! 5 Ross Sea
       6.532956e-08_dp , &   ! 6 Amundsen Sea
       3.740965e-08_dp , &   ! 7 Bellingshausen Sea
       1.394415e-08_dp  /)   ! 8 Weddell Sea
  real(kind=dp), parameter :: mmm_shelf_global_lambda(0:8) = (/ &
       1.979433_dp     , &   ! 0 Global               lambda*Gamma = 3.3767e-08
       1.154975_dp     , &   ! 1 Dronning Maud Land   lambda*Gamma = 5.7877e-08
       1.159670_dp     , &   ! 2 Amery-Prydz Bay      lambda*Gamma = 7.0396e-08
       1.201389_dp     , &   ! 3 Sabrina-Aurora       lambda*Gamma = 5.8855e-08
       1.123631_dp     , &   ! 4 George V-Wilkes      lambda*Gamma = 6.9142e-08
       2.510084_dp     , &   ! 5 Ross Sea             lambda*Gamma = 1.6724e-08
       1.198968_dp     , &   ! 6 Amundsen Sea         lambda*Gamma = 7.8328e-08
       1.242801_dp     , &   ! 7 Bellingshausen Sea   lambda*Gamma = 4.6493e-08
       1.311055_dp      /)   ! 8 Weddell Sea          lambda*Gamma = 1.8282e-08
  real(kind=dp), parameter :: mmm_shelf_global_lambdaw(0:8) = (/ &
       1.566711_dp     , &   ! 0 Global               whole base: lambda_w*Gamma = 2.6727e-08
       1.051080_dp     , &   ! 1 Dronning Maud Land   whole base: lambda_w*Gamma = 5.2671e-08
       1.115217_dp     , &   ! 2 Amery-Prydz Bay      whole base: lambda_w*Gamma = 6.7698e-08
       1.096507_dp     , &   ! 3 Sabrina-Aurora       whole base: lambda_w*Gamma = 5.3716e-08
       0.940393_dp     , &   ! 4 George V-Wilkes      whole base: lambda_w*Gamma = 5.7867e-08
       1.791075_dp     , &   ! 5 Ross Sea             whole base: lambda_w*Gamma = 1.1934e-08
       1.139248_dp     , &   ! 6 Amundsen Sea         whole base: lambda_w*Gamma = 7.4427e-08
       1.163527_dp     , &   ! 7 Bellingshausen Sea   whole base: lambda_w*Gamma = 4.3527e-08
       0.919948_dp      /)   ! 8 Weddell Sea          whole base: lambda_w*Gamma = 1.2828e-08
! the arrays the melt code indexes: (sector 0:8, option 1:3, product 1:2)
  real(kind=dp), parameter :: sec_alpha(0:8,3,2) = reshape( (/ &
       waom4_interface_mw_alpha, waom4_shelf_sector_alpha, waom4_shelf_global_alpha, &
       mmm_interface_mw_alpha, mmm_shelf_sector_alpha, mmm_shelf_global_alpha /), (/ 9, 3, 2 /) )
  real(kind=dp), parameter :: sec_Gamma(0:8,3,2) = reshape( (/ &
       waom4_interface_mw_Gamma, waom4_shelf_sector_Gamma, waom4_shelf_global_Gamma, &
       mmm_interface_mw_Gamma, mmm_shelf_sector_Gamma, mmm_shelf_global_Gamma /), (/ 9, 3, 2 /) )
  real(kind=dp), parameter :: sec_lambda(0:8,3,2) = reshape( (/ &
       waom4_interface_mw_lambda, waom4_shelf_sector_lambda, waom4_shelf_global_lambda, &
       mmm_interface_mw_lambda, mmm_shelf_sector_lambda, mmm_shelf_global_lambda /), (/ 9, 3, 2 /) )
  real(kind=dp), parameter :: sec_lambdaw(0:8,3,2) = reshape( (/ &
       waom4_interface_mw_lambdaw, waom4_shelf_sector_lambdaw, waom4_shelf_global_lambdaw, &
       mmm_interface_mw_lambdaw, mmm_shelf_sector_lambdaw, mmm_shelf_global_lambdaw /), (/ 9, 3, 2 /) )
! <<< END sector melt coefficients
#endif

#if defined(ALLOW_TAPENADE) /* Tapenade */
call myfloor(time_in_years, i_time_in_years)
#endif /* Tapenade */

!-------- Determine local vs. semi-local parameterization --------

#if (FLOATING_ICE_BASAL_MELTING==6)
flag_param_semilocal = .true.
                      ! default is semi-local for the ISMIP6 parameterization
#elif (FLOATING_ICE_BASAL_MELTING==7)
flag_param_semilocal = .false.
                      ! default is local for the ISMIP7 parameterization
#elif (FLOATING_ICE_BASAL_MELTING==8)
flag_param_semilocal = .false.
                      ! the sector relationships take the local thermal forcing
#endif

#if (FLOATING_ICE_BASAL_MELTING==6 || FLOATING_ICE_BASAL_MELTING==7)

#if (defined(PARAM_LOCAL_SEMILOCAL))
#if (PARAM_LOCAL_SEMILOCAL==1)
   flag_param_semilocal = .false.  ! local
#elif (PARAM_LOCAL_SEMILOCAL==2)
   flag_param_semilocal = .true.   ! semi-local
#else
   errormsg = ' >>> sub_ice_shelf_melting_param_2:' &
            //         end_of_line &
            //'        PARAM_LOCAL_SEMILOCAL must be either 1 or 2!'
   call error(errormsg)
#endif
#endif

#endif

!-------- Read file with the thermal forcing data of the ocean --------

if (TF_BM_TIME_MIN > TF_BM_TIME_MAX) then
   errormsg = ' >>> sub_ice_shelf_melting_param_2:' &
            //         end_of_line &
            //'        TF_BM_TIME_MIN > TF_BM_TIME_MAX!'
   call error(errormsg)
end if

n_year_CE_tf_bm = n_year_CE

if (n_year_CE_tf_bm < TF_BM_TIME_MIN) then
   n_year_CE_tf_bm = TF_BM_TIME_MIN
else if (n_year_CE_tf_bm > TF_BM_TIME_MAX) then
   n_year_CE_tf_bm = TF_BM_TIME_MAX
end if

if ( firstcall%sub_ice_shelf_melting_param_2 &
     .or.(n_year_CE_tf_bm /= n_year_CE_tf_bm_save) ) then

   write(ch_year_CE, '(i0)') n_year_CE_tf_bm

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

      istat = nf90_inq_varid(ncid, 'z', ncv)
      if (istat /= nf90_noerr) then
         errormsg = ' >>> sub_ice_shelf_melting_param_2:' &
                  //                end_of_line &
                  //'               Error when inquiring the variable' &
                  //                end_of_line &
                  //'               for the vertical coordinate!'
         call error(errormsg)
      end if

      call check( nf90_get_var(ncid, ncv, z_tf_bm) )

      istat1 = nf90_inq_varid(ncid, 'thermal_forcing', ncv)
      if (istat1 /= nf90_noerr) then
         istat2 = nf90_inq_varid(ncid, 'tf', ncv)
         if (istat2 /= nf90_noerr) then
            errormsg = ' >>> sub_ice_shelf_melting_param_2:' &
                     //                end_of_line &
                     //'               Error when inquiring the variable' &
                     //                end_of_line &
                     //'               for the thermal forcing!'
            call error(errormsg)
         end if
      end if

      call check( nf90_get_var(ncid, ncv, tf_bm_aux), thisroutine )

      call check( nf90_close(ncid), thisroutine )

!  ------ Ensure positive depth values

      if ( (z_tf_bm(0) < eps_dp).and.(z_tf_bm(NZ_TF_BM) < eps_dp) ) &
         z_tf_bm = -z_tf_bm

!  ------ Swap indices -> SICOPOLIS standard

      do i=0, IMAX
      do j=0, JMAX
      do n=0, NZ_TF_BM

         tf_bm(n,j,i) = tf_bm_aux(i,j,n)

         if (tf_bm(n,j,i) > r_no_value_pos_1) &
            tf_bm(n,j,i) = r_no_value_neg_2

         if (isnan(tf_bm(n,j,i))) &
            tf_bm(n,j,i) = r_no_value_neg_2
               !%% \!/ 'isnan' checks don't work reliably
               !%%     when using aggressive maths optimization \!/

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

!  ------ Save value of n_year_CE_tf_bm

n_year_CE_tf_bm_save = n_year_CE_tf_bm

!-------- Read file with the salinity data of the ocean --------

#if (FLOATING_ICE_BASAL_MELTING==7)

if (SO_BM_TIME_MIN > SO_BM_TIME_MAX) then
   errormsg = ' >>> sub_ice_shelf_melting_param_2:' &
            //         end_of_line &
            //'        SO_BM_TIME_MIN > SO_BM_TIME_MAX!'
   call error(errormsg)
end if

n_year_CE_so_bm = n_year_CE

if (n_year_CE_so_bm < SO_BM_TIME_MIN) then
   n_year_CE_so_bm = SO_BM_TIME_MIN
else if (n_year_CE_so_bm > SO_BM_TIME_MAX) then
   n_year_CE_so_bm = SO_BM_TIME_MAX
end if

if ( firstcall%sub_ice_shelf_melting_param_2 &
     .or.(n_year_CE_so_bm /= n_year_CE_so_bm_save) ) then

   write(ch_year_CE, '(i0)') n_year_CE_so_bm

   if ( (trim(adjustl(SO_BM_FILES)) /= 'none') &
        .and. &
        (trim(adjustl(SO_BM_FILES)) /= 'None') &
        .and. &
        (trim(adjustl(SO_BM_FILES)) /= 'NONE') ) then

!  ------ Read data from file

      filename_with_path = trim(SO_BM_DIR)//'/'// &
                           trim(SO_BM_FILES)//trim(ch_year_CE)//'.nc'

      ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

      if (ios /= nf90_noerr) then
         errormsg = ' >>> sub_ice_shelf_melting_param_2:' &
                  //                end_of_line &
                  //'               Error when opening the file' &
                  //                end_of_line &
                  //'               for the salinity of the ocean!'
         call error(errormsg)
      end if

      istat = nf90_inq_varid(ncid, 'z', ncv)
      if (istat /= nf90_noerr) then
         errormsg = ' >>> sub_ice_shelf_melting_param_2:' &
                  //                end_of_line &
                  //'               Error when inquiring the variable' &
                  //                end_of_line &
                  //'               for the vertical coordinate!'
         call error(errormsg)
      end if

      call check( nf90_get_var(ncid, ncv, z_so_bm) )

      istat1 = nf90_inq_varid(ncid, 'salinity', ncv)
      if (istat1 /= nf90_noerr) then
         istat2 = nf90_inq_varid(ncid, 'so', ncv)
         if (istat2 /= nf90_noerr) then
            errormsg = ' >>> sub_ice_shelf_melting_param_2:' &
                     //                end_of_line &
                     //'               Error when inquiring the variable' &
                     //                end_of_line &
                     //'               for the salinity!'
            call error(errormsg)
         end if
      end if

      call check( nf90_get_var(ncid, ncv, so_bm_aux), thisroutine )

      call check( nf90_close(ncid), thisroutine )

!  ------ Ensure positive depth values

      if ( (z_so_bm(0) < eps_dp).and.(z_so_bm(NZ_SO_BM) < eps_dp) ) &
         z_so_bm = -z_so_bm

!  ------ Swap indices -> SICOPOLIS standard

      do i=0, IMAX
      do j=0, JMAX
      do n=0, NZ_SO_BM

         so_bm(n,j,i) = so_bm_aux(i,j,n)

         if ( (so_bm(n,j,i) > r_no_value_pos_2) &
              .or. &
              (so_bm(n,j,i) < 0.0_dp) ) then
            so_bm(n,j,i) = so_bm_default
         end if

      end do
      end do
      end do

!  ------ Check consistency of the depth (z_so_bm) data

#if !defined(ALLOW_TAPENADE) /* Normal */

      if (.not.(approx_equal(z_so_bm(0), ZMIN_SO_BM, eps_sp_dp))) then
         errormsg = ' >>> sub_ice_shelf_melting_param_2:' &
                  //         end_of_line &
                  //'        Inconsistency between' &
                  //         end_of_line &
                  //'        read z_so_bm data' &
                  //         end_of_line &
                  //'        and parameter ZMIN_SO_BM!'
         call error(errormsg)
      end if

      if (.not.(approx_equal(z_so_bm(NZ_SO_BM)-z_so_bm(0), &
                             NZ_SO_BM*DZ_SO_BM, eps_sp_dp))) then
         errormsg = ' >>> sub_ice_shelf_melting_param_2:' &
                  //         end_of_line &
                  //'        Inconsistency between' &
                  //         end_of_line &
                  //'        read z_so_bm data' &
                  //         end_of_line &
                  //'        and parameters NZ_SO_BM, DZ_SO_BM!'
         call error(errormsg)
      end if

#endif /* Normal */

   else   ! ( trim(adjustl(SO_BM_FILES)) == 'none' or 'None' or 'NONE' )

!  ------ Use present-day data

      so_bm_aux = 0.0_dp
      so_bm     = so_bm_present
      z_so_bm   = z_so_bm_present

   end if

end if

!  ------ Save value of n_year_CE_so_bm

n_year_CE_so_bm_save = n_year_CE_so_bm

#endif

!-------- Parameters for the parameterization --------

#if (FLOATING_ICE_BASAL_MELTING==6)
slope_mean = 0.0_dp   ! dummy value
#elif (FLOATING_ICE_BASAL_MELTING==7)
#if (defined(SLOPE_MEAN_ICE_DRAFT))
slope_mean = SLOPE_MEAN_ICE_DRAFT
#else
errormsg = ' >>> sub_ice_shelf_melting_param_2:' &
         //         end_of_line &
         //'        SLOPE_MEAN_ICE_DRAFT must be defined!'
#endif
#endif

#if (FLOATING_ICE_BASAL_MELTING==6)
cst_bm = ((rhosw_bm*cpw_bm)/(rhoi_bm*Lf_bm))**2   ! 1/K2
#elif (FLOATING_ICE_BASAL_MELTING==7)
cst_bm = (rhosw_bm/rhoi_bm)*(cpw_bm/Lf_bm)**2     ! 1/K2
cst_bm = cst_bm * (0.5_dp*G/abs(coriolis_par_bm)) ! (m/s)/K2
#endif

#if (!defined(N_BM_REGIONS) || N_BM_REGIONS<=1)
n_bm_regions = 1
#else
n_bm_regions = N_BM_REGIONS
#endif

#if (FLOATING_ICE_BASAL_MELTING==8)

#if (defined(BM_SECTOR_PRODUCT))
i_product = BM_SECTOR_PRODUCT
#else
i_product = 1   ! WAOM4
#endif
#if (defined(BM_SECTOR_DRIVING))
errormsg = ' >>> sub_ice_shelf_melting_param_2: ' &
              //'BM_SECTOR_DRIVING was removed; use BM_SECTOR_APPROACH 1 to 3!'
call error(errormsg)
#endif
#if (defined(BM_SECTOR_APPROACH))
i_approach = BM_SECTOR_APPROACH
#else
i_approach = 3   ! shelf-global: the ISMIP6-style thermal forcing read by this routine, global alpha_cs
#endif
#if (defined(BM_SECTOR_CALIB))
i_calib = BM_SECTOR_CALIB
#else
if (i_approach == 1) then
   i_calib = 1   ! fitted cells: the ice-base option predicts nothing where TF <= 0
else
   i_calib = 2   ! whole base: the shelf options predict melt wherever the shelf driving is positive
end if
#endif
#if (defined(DELTA_TF_SECTOR))
delta_tf_sector = DELTA_TF_SECTOR
#else
delta_tf_sector = 0.0_dp
#endif
if ((i_product < 1).or.(i_product > 2).or.(i_approach < 1).or.(i_approach > 3) &
    .or.(i_calib < 1).or.(i_calib > 2)) then
   errormsg = ' >>> sub_ice_shelf_melting_param_2: ' &
                 //'BM_SECTOR_PRODUCT must be 1 or 2, BM_SECTOR_APPROACH 1 to 3, BM_SECTOR_CALIB 1 or 2!'
   call error(errormsg)
end if

#else   /* FLOATING_ICE_BASAL_MELTING==6 or 7 */

#if (defined(GAMMA0_BM))
#if (FLOATING_ICE_BASAL_MELTING==6)
  gamma0_bm_aux = GAMMA0_BM
  gamma0_bm_aux = gamma0_bm_aux *sec2year*(rhofw_bm/rhoi_bm)
                                ! m/a water equiv. -> m/s ice equiv.
#elif (FLOATING_ICE_BASAL_MELTING==7)
  gamma0_bm_aux = GAMMA0_BM     ! dimensionless parameter
#endif
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

#endif   /* FLOATING_ICE_BASAL_MELTING==6 / 7 / 8 */

!-------- Local thermal forcing --------

dz_tf_inv = 1.0_dp/DZ_TF_BM

tf_bm_no_value_neg = 0.999_dp*r_no_value_neg_2

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

         real_n = (draft-ZMIN_TF_BM)*dz_tf_inv
         ninf   = floor(real_n)
         nsup   = ninf+1

         if ((ninf >= 0).and.(nsup <= NZ_TF_BM)) then

            if ( &
                 (tf_bm(ninf,j,i) > tf_bm_no_value_neg) &
                 .and. &
                 (tf_bm(nsup,j,i) > tf_bm_no_value_neg) &
               ) &
               tf_bm_local(j,i) = tf_bm(ninf,j,i) &
                                    + ((draft-z_tf_bm(ninf))*dz_tf_inv) &
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

if (flag_param_semilocal) then

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

end if

!-------- Local salinity --------

#if (FLOATING_ICE_BASAL_MELTING==7)

dz_so_inv = 1.0_dp/DZ_SO_BM

so_bm_local = so_bm_default   ! initialization by default value

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

         real_n = (draft-ZMIN_SO_BM)*dz_so_inv
         ninf   = floor(real_n)
         nsup   = ninf+1

         if ((ninf >= 0).and.(nsup <= NZ_SO_BM)) then

            so_bm_local(j,i) = so_bm(ninf,j,i) &
                                 + ((draft-z_so_bm(ninf))*dz_so_inv) &
                                   *(so_bm(nsup,j,i)-so_bm(ninf,j,i))

         else if (ninf < 0) then

            ninf = 0
            so_bm_local(j,i) = so_bm(ninf,j,i)

         else   ! nsup > NZ_SO_BM

            nsup = NZ_SO_BM
            so_bm_local(j,i) = so_bm(nsup,j,i)

         end if

      end if

   end if

end do
end do

#endif

!-------- Computation of the sub-ice-shelf melting rate --------

do i=0, IMAX
do j=0, JMAX

   if ( (mask(j,i)==2).or.(mask(j,i)==3) ) then
                                                ! floating ice or ocean

!  ------ Melting rate over continental shelf

      if ( zl(j,i) > z_abyssal ) then

         n = n_bm_region(j,i)

         if (flag_param_semilocal) then
            tf_bm_help = tf_bm_ave(n)
         else
            tf_bm_help = tf_bm_local(j,i)
         end if

#if (FLOATING_ICE_BASAL_MELTING==6)

         Q_bm(j,i) = gamma0_bm(j,i)*cst_bm &
                        *(tf_bm_local(j,i)+delta_tf_bm(j,i)) &
                        *abs(tf_bm_help+delta_tf_bm(j,i))

#elif (FLOATING_ICE_BASAL_MELTING==7)

         Q_bm(j,i) = (gamma0_bm(j,i)*sin(slope_mean)*cst_bm) &
                        *(beta_s_bm*so_bm_local(j,i)) &
                        *(tf_bm_local(j,i)+delta_tf_bm(j,i)) &
                        *abs(tf_bm_help+delta_tf_bm(j,i))

#elif (FLOATING_ICE_BASAL_MELTING==8)

!    ---- Sector-wise power relationship Q_bm = lambda*Gamma*TF^alpha
!         (m/s ice equivalent; TF <= 0 gives no melt, no refreezing predicted)

#if (defined(BM_SECTOR_GLOBAL) && BM_SECTOR_GLOBAL==1)
         n_sector = 0
#else
         lon_d = lambda(j,i) *rad2deg
         lat_d = phi(j,i)    *rad2deg
         n_sector = sector_from_lon_lat(lon_d, lat_d)
#endif

         tf_sector = max(tf_bm_local(j,i)+delta_tf_sector, 0.0_dp)

         if (i_calib == 2) then
            lam_sector_use = sec_lambdaw(n_sector,i_approach,i_product)
         else
            lam_sector_use = sec_lambda(n_sector,i_approach,i_product)
         end if
         Q_bm(j,i) = lam_sector_use &
                     *sec_Gamma(n_sector,i_approach,i_product) &
                     *tf_sector**sec_alpha(n_sector,i_approach,i_product)

#endif

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

!    ---- Correction for ISMIP6/7 LARMIP

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

#else   /* not (defined(ANT) && FLOATING_ICE_BASAL_MELTING==6, 7 or 8) */

errormsg = ' >>> sub_ice_shelf_melting_param_2:' &
         //          end_of_line &
         //'         Routine only valid for Antarctica' &
         //          end_of_line &
         //'         and FLOATING_ICE_BASAL_MELTING==6, 7 or 8!'
call error(errormsg)

#endif

end subroutine sub_ice_shelf_melting_param_2

#if (FLOATING_ICE_BASAL_MELTING==8)
! >>> BEGIN sector helper routines
!-------------------------------------------------------------------------------
!> Sector ids 1..8 from longitude/latitude in degrees, the box logic of
!! classify_sectors.m of the analysis pipeline (first match wins, the rest is
!! the Weddell Sea): 1 Dronning Maud Land, 2 Amery-Prydz Bay, 3 Sabrina-Aurora,
!! 4 George V-Wilkes, 5 Ross Sea, 6 Amundsen Sea, 7 Bellingshausen Sea,
!! 8 Weddell Sea.
!-------------------------------------------------------------------------------
function sector_from_lon_lat(lon_in, lat_d) result(n_sector)

implicit none

real(dp), intent(in) :: lon_in, lat_d
integer(i4b)         :: n_sector

real(dp) :: lon_d

#if !defined(ALLOW_TAPENADE) /* Normal */
lon_d = modulo(lon_in+180.0_dp, 360.0_dp)-180.0_dp
                                  ! mapping to interval
                                  ! [-180 deg, +180 deg)
#else /* Tapenade */
lon_d = (lon_in+180.0_sp) &
           - ((360.0_sp * int(lon_in+180.0_sp)/360.0_sp)) - 180.0_sp
#endif /* Normal vs. Tapenade */

if ((lon_d>=-10.0_dp).and.(lon_d<60.0_dp)) then
   n_sector = 1                          ! Dronning Maud Land
else if ((lon_d>=60.0_dp).and.(lon_d<80.0_dp)) then
   n_sector = 2                          ! Amery/Prydz Bay
else if ((lon_d>=80.0_dp).and.(lon_d<130.0_dp)) then
   n_sector = 3                          ! Sabrina Coast/Aurora subglacial basin
else if ( ((lon_d>=130.0_dp).and.(lon_d<159.0_dp)) &
          .or. &
          ((lon_d>=159.0_dp).and.(lon_d<170.0_dp) &
                            .and.(lat_d>=-72.0_dp)) ) then
   n_sector = 4                          ! George V Coast/Wilkes subglacial basin
else if ( (lon_d>=159.0_dp) &
          .or. &
          (lon_d<-140.0_dp) &
          .or. &
          ((lon_d>=-140.0_dp).and.(lon_d<-120.0_dp) &
                             .and.(lat_d<-77.0_dp)) ) then
   n_sector = 5                          ! Ross Sea
else if ( ((lon_d>=-140.0_dp).and.(lon_d<-120.0_dp) &
                             .and.(lat_d>=-77.0_dp)) &
          .or. &
          ((lon_d>=-120.0_dp).and.(lon_d<-90.0_dp)) ) then
   n_sector = 6                          ! Amundsen Sea
else if ( ((lon_d>=-90.0_dp).and.(lon_d<-66.0_dp) &
                            .and.(lat_d>=-74.0_dp)) ) then
   n_sector = 7                          ! Bellingshausen Sea
else
   n_sector = 8                          ! Weddell Sea
end if

end function sector_from_lon_lat
! <<< END sector helper routines
#endif

!-------------------------------------------------------------------------------

end module calc_bas_melt_m
!
