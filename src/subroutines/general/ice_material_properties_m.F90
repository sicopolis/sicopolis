!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  i c e _ m a t e r i a l _ p r o p e r t i e s _ m
!
!> @file
!!
!! Material properties of ice:
!! Rate factor, heat conductivity, specific heat (heat capacity),
!! creep function, viscosity.
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
!> Material properties of ice:
!! Rate factor, heat conductivity, specific heat (heat capacity),
!! creep function, viscosity.
!<------------------------------------------------------------------------------
module ice_material_properties_m

use sico_types_m
use error_m

implicit none
save

#if !defined(ALLOW_TAPENADE) /* Normal */

!> RF(n): Tabulated values for the rate factor of cold ice
   real(dp), dimension(-256:255), private       :: RF

!> R_T: Coefficient of the water-content dependence in the rate factor
!>      for temperate ice
   real(dp)                     , private       :: R_T

!> KAPPA(n): Tabulated values for the heat conductivity of ice
   real(dp), dimension(-256:255), private       :: KAPPA

!> C(n): Tabulated values for the specific heat of ice
   real(dp), dimension(-256:255), private       :: C

!> n_temp_min: Lower index limit of properly defined values in RF, KAPPA and C
!>             (n_temp_min >= -256).
   integer(i4b), private :: n_temp_min

!> n_temp_max: Upper index limit of properly defined values in RF, KAPPA and C
!>             (n_temp_max <= 255).
   integer(i4b), private :: n_temp_max

!> RHO_I: Density of ice
   real(dp)                     , private       :: RHO_I
                                                ! only for the Martian ice caps
!> RHO_C: Density of crustal material (dust)
   real(dp)                     , private       :: RHO_C
                                                ! only for the Martian ice caps
!> KAPPA_C: Heat conductivity of crustal material (dust)
   real(dp)                     , private       :: KAPPA_C
                                                ! only for the Martian ice caps
!> C_C: Specific heat of crustal material (dust)
   real(dp)                     , private       :: C_C
                                                ! only for the Martian ice caps

private
public :: ice_mat_eqs_pars, &
          ratefac_c, ratefac_t, ratefac_c_t, kappa_val, c_val, &
          viscosity, creep

#else /* Tapenade */

real(dp), dimension(-256:255), public       :: RF_IMQ
real(dp)                     , public       :: R_T_IMQ
real(dp), dimension(-256:255), public       :: KAPPA_IMQ
real(dp), dimension(-256:255), public       :: C_IMQ
integer(i4b), public                        :: n_temp_min_IMQ
integer(i4b), public                        :: n_temp_max_IMQ
real(dp)                     , public       :: RHO_I_IMQ
real(dp)                     , public       :: RHO_C_IMQ
real(dp)                     , public       :: KAPPA_C_IMQ
real(dp)                     , public       :: C_C_IMQ

public :: ice_mat_eqs_pars, &
          ratefac_c, ratefac_t, ratefac_c_t, kappa_val, c_val, &
          viscosity, creep
private :: myfloor_local

#endif /* Normal vs. Tapenade */

contains

!-------------------------------------------------------------------------------
!> Setting of required physical parameters.
!<------------------------------------------------------------------------------
subroutine ice_mat_eqs_pars(RF_table, R_T_val, KAPPA_table, C_table, &
                            n_tmp_min, n_tmp_max, &
                            RHO_I_val, RHO_C_val, KAPPA_C_val, C_C_val)

implicit none

integer(i4b),       intent(in) :: n_tmp_min, n_tmp_max
real(dp), dimension(n_tmp_min:n_tmp_max), &
                    intent(in) :: RF_table, KAPPA_table, C_table
real(dp),           intent(in) :: R_T_val
real(dp), optional, intent(in) :: RHO_I_val, RHO_C_val, KAPPA_C_val, C_C_val

integer(i4b)       :: n
character(len=256) :: errormsgg

!-------- Initialisation --------

#if !defined(ALLOW_TAPENADE) /* Normal */

RF    = 0.0_dp
KAPPA = 0.0_dp
C     = 0.0_dp
n_temp_min = n_tmp_min
n_temp_max = n_tmp_max

#else /* Tapenade */

RF_IMQ    = 0.0_dp
KAPPA_IMQ = 0.0_dp
C_IMQ     = 0.0_dp
n_temp_min_IMQ = n_tmp_min
n_temp_max_IMQ = n_tmp_max

#endif /* Normal vs. Tapenade */

#if !defined(ALLOW_TAPENADE) /* Normal */
if ((n_temp_min <= -256).or.(n_temp_max >= 255)) then
#else /* Tapenade */
if ((n_temp_min_IMQ <= -256).or.(n_temp_max_IMQ >= 255)) then
#endif /* Normal vs. Tapenade */
   errormsgg = ' >>> ice_mat_eqs_pars: ' &
                  //'Temperature indices out of allowed range!'
   call error(errormsgg)
end if

!-------- Assignment --------

#if !defined(ALLOW_TAPENADE) /* Normal */

do n=n_temp_min, n_temp_max
   RF(n)    = RF_table(n)
   KAPPA(n) = KAPPA_table(n)
   C(n)     = C_table(n)
end do

#else /* Tapenade */

do n=n_temp_min_IMQ, n_temp_max_IMQ
   RF_IMQ(n)    = RF_table(n)
   KAPPA_IMQ(n) = KAPPA_table(n)
   C_IMQ(n)     = C_table(n)
end do

#endif /* Normal vs. Tapenade */

#if !defined(ALLOW_TAPENADE) /* Normal */

do n=-256, n_temp_min-1
   RF(n)    = RF(n_temp_min)      ! dummy values
   KAPPA(n) = KAPPA(n_temp_min)   ! dummy values
   C(n)     = C(n_temp_min)       ! dummy values
end do

#else /* Tapenade */

do n=-256, n_temp_min_IMQ-1
   RF_IMQ(n)    = RF_IMQ(n_temp_min_IMQ)      ! dummy values
   KAPPA_IMQ(n) = KAPPA_IMQ(n_temp_min_IMQ)   ! dummy values
   C_IMQ(n)     = C_IMQ(n_temp_min_IMQ)       ! dummy values
end do

#endif /* Normal vs. Tapenade */

#if !defined(ALLOW_TAPENADE) /* Normal */

do n=n_temp_max+1, 255
   RF(n)    = RF(n_temp_max)      ! dummy values
   KAPPA(n) = KAPPA(n_temp_max)   ! dummy values
   C(n)     = C(n_temp_max)       ! dummy values
end do

#else /* Tapenade */

do n=n_temp_max_IMQ+1, 255
   RF_IMQ(n)    = RF_IMQ(n_temp_max_IMQ)      ! dummy values
   KAPPA_IMQ(n) = KAPPA_IMQ(n_temp_max_IMQ)   ! dummy values
   C_IMQ(n)     = C_IMQ(n_temp_max_IMQ)       ! dummy values
end do

#endif /* Normal vs. Tapenade */

#if !defined(ALLOW_TAPENADE) /* Normal */
R_T = R_T_val
#else /* Tapenade */
R_T_IMQ = R_T_val
#endif /* Normal vs. Tapenade */

!-------- Martian stuff --------

#if !defined(ALLOW_TAPENADE) /* Normal */

if ( present(RHO_I_val) ) then
   RHO_I = RHO_I_val
else
   RHO_I = 0.0_dp   ! dummy value
end if

#else /* Tapenade */

if ( present(RHO_I_val) ) then
   RHO_I_IMQ = RHO_I_val
else
   RHO_I_IMQ = 0.0_dp   ! dummy value
end if

#endif /* Normal vs. Tapenade */

#if !defined(ALLOW_TAPENADE) /* Normal */

if ( present(RHO_C_val) ) then
   RHO_C = RHO_C_val
else
   RHO_C = 0.0_dp   ! dummy value
end if

#else /* Tapenade */

if ( present(RHO_C_val) ) then
   RHO_C_IMQ = RHO_C_val
else
   RHO_C_IMQ = 0.0_dp   ! dummy value
end if

#endif /* Normal vs. Tapenade */

#if !defined(ALLOW_TAPENADE) /* Normal */

if ( present(KAPPA_C_val) ) then
   KAPPA_C = KAPPA_C_val
else
   KAPPA_C = 0.0_dp   ! dummy value
end if

#else /* Tapenade */

if ( present(KAPPA_C_val) ) then
   KAPPA_C_IMQ = KAPPA_C_val
else
   KAPPA_C_IMQ = 0.0_dp   ! dummy value
end if

#endif /* Normal vs. Tapenade */

#if !defined(ALLOW_TAPENADE) /* Normal */

if ( present(C_C_val) ) then
   C_C = C_C_val
else
   C_C = 0.0_dp   ! dummy value
end if

#else /* Tapenade */

if ( present(C_C_val) ) then
   C_C_IMQ = C_C_val
else
   C_C_IMQ = 0.0_dp   ! dummy value
end if

#endif /* Normal vs. Tapenade */

end subroutine ice_mat_eqs_pars

!-------------------------------------------------------------------------------
!> Rate factor for cold ice:
!! Linear interpolation of tabulated values in RF(.).
!<------------------------------------------------------------------------------
function ratefac_c(temp_val, temp_m_val)

use sico_variables_m
use sico_vars_m

implicit none
real(dp)             :: ratefac_c
real(dp), intent(in) :: temp_val, temp_m_val

integer(i4b) :: n_temp_1, n_temp_2
real(dp)     :: temp_h_val
#if defined(ALLOW_TAPENADE) /* Tapenade */
integer(i4b) :: itemp_h_val
#endif /* Tapenade */

temp_h_val = temp_val-temp_m_val

#if !defined(ALLOW_TAPENADE) /* Normal */

n_temp_1 = floor(temp_h_val)
n_temp_1 = max(min(n_temp_1, n_temp_max-1), n_temp_min)
n_temp_2 = n_temp_1 + 1

ratefac_c = RF(n_temp_1) &
            + (RF(n_temp_2)-RF(n_temp_1)) &
              * (temp_h_val-real(n_temp_1,dp))   ! Linear interpolation

#else /* Tapenade */

n_temp_1 = myfloor_local(temp_h_val)
n_temp_1 = max(min(n_temp_1, n_temp_max_IMQ-1), n_temp_min_IMQ)
n_temp_2 = n_temp_1 + 1

ratefac_c = RF_IMQ(n_temp_1) &
              +(RF_IMQ(n_temp_2)-RF_IMQ(n_temp_1)) &
               *(temp_h_val-real(n_temp_1,dp))

#endif /* Normal vs. Tapenade */

end function ratefac_c

!-------------------------------------------------------------------------------
!> Rate factor for temperate ice.
!<------------------------------------------------------------------------------
function ratefac_t(omega_val)

use sico_variables_m
use sico_vars_m

implicit none
real(dp)             :: ratefac_t
real(dp), intent(in) :: omega_val

#if !defined(ALLOW_TAPENADE) /* Normal */
ratefac_t = RF(0)*(1.0_dp+R_T*(omega_val))
#else /* Tapenade */
ratefac_t = RF_IMQ(0)*(1.0_dp+R_T_IMQ*(omega_val))
#endif /* Normal vs. Tapenade */

end function ratefac_t

!-------------------------------------------------------------------------------
!> Rate factor for cold and temperate ice:
!! Combination of ratefac_c and ratefac_t (only for the enthalpy method).
!<------------------------------------------------------------------------------
function ratefac_c_t(temp_val, omega_val, temp_m_val)

use sico_variables_m
use sico_vars_m

implicit none
real(dp)             :: ratefac_c_t
real(dp), intent(in) :: temp_val, temp_m_val, omega_val

integer(i4b) :: n_temp_1, n_temp_2
real(dp)     :: temp_h_val
#if defined(ALLOW_TAPENADE) /* Tapenade */
integer(i4b) :: itemp_h_val
#endif /* Tapenade */

temp_h_val = temp_val-temp_m_val

#if !defined(ALLOW_TAPENADE) /* Normal */

n_temp_1 = floor(temp_h_val)
n_temp_1 = max(min(n_temp_1, n_temp_max-1), n_temp_min)
n_temp_2 = n_temp_1 + 1

ratefac_c_t = ( RF(n_temp_1) &
                + (RF(n_temp_2)-RF(n_temp_1)) &
                  * (temp_h_val-real(n_temp_1,dp)) ) &
              *(1.0_dp+R_T*(omega_val))

#else /* Tapenade */

n_temp_1 = myfloor_local(temp_h_val)
n_temp_1 = max(min(n_temp_1, n_temp_max_IMQ-1), n_temp_min_IMQ)
n_temp_2 = n_temp_1 + 1

ratefac_c_t = ( RF_IMQ(n_temp_1) &
                  +(RF_IMQ(n_temp_2)-RF_IMQ(n_temp_1)) &
                   *(temp_h_val-real(n_temp_1,dp)) ) &
              *(1.0_dp+R_T_IMQ*(omega_val))

#endif /* Normal vs. Tapenade */

end function ratefac_c_t

!-------------------------------------------------------------------------------
!> Heat conductivity of ice:
!! Linear interpolation of tabulated values in KAPPA(.).
!<------------------------------------------------------------------------------
function kappa_val(temp_val)

use sico_variables_m
use sico_vars_m

implicit none
real(dp)             :: kappa_val
real(dp), intent(in) :: temp_val

integer(i4b) :: n_temp_1, n_temp_2
real(dp) :: kappa_ice
#if defined(ALLOW_TAPENADE) /* Tapenade */
integer(i4b) :: itemp_val
#endif /* Tapenade */

!-------- Heat conductivity of pure ice --------

#if !defined(ALLOW_TAPENADE) /* Normal */

n_temp_1 = floor(temp_val)
n_temp_1 = max(min(n_temp_1, n_temp_max-1), n_temp_min)
n_temp_2 = n_temp_1 + 1

#if defined(FRAC_DUST)
kappa_ice = KAPPA(n_temp_1) &
            + (KAPPA(n_temp_2)-KAPPA(n_temp_1)) &
              * (temp_val-real(n_temp_1,dp))
#else
kappa_val = KAPPA(n_temp_1) &
            + (KAPPA(n_temp_2)-KAPPA(n_temp_1)) &
              * (temp_val-real(n_temp_1,dp))
#endif

#else /* Tapenade */

n_temp_1 = myfloor_local(temp_val)
n_temp_1 = max(min(n_temp_1, n_temp_max_IMQ-1), n_temp_min_IMQ)
n_temp_2 = n_temp_1 + 1

#if defined(FRAC_DUST)
kappa_ice = KAPPA_IMQ(n_temp_1) &
            + (KAPPA_IMQ(n_temp_2)-KAPPA_IMQ(n_temp_1)) &
              * (temp_val-real(n_temp_1,dp))
#else
kappa_val = KAPPA_IMQ(n_temp_1) &
            + (KAPPA_IMQ(n_temp_2)-KAPPA_IMQ(n_temp_1)) &
              * (temp_val-real(n_temp_1,dp))
#endif

#endif /* Normal vs. Tapenade */

!-------- If dust is present (polar caps of Mars):
!         Heat conductivity of ice-dust mixture --------

#if defined(FRAC_DUST)
#if !defined(ALLOW_TAPENADE) /* Normal */
kappa_val = (1.0_dp-FRAC_DUST)*kappa_ice + FRAC_DUST*KAPPA_C
#else /* Tapenade */
kappa_val = (1.0_dp-FRAC_DUST)*kappa_ice + FRAC_DUST*KAPPA_C_IMQ
#endif /* Normal vs. Tapenade */
#endif

end function kappa_val

!-------------------------------------------------------------------------------
!> Specific heat of ice:
!! Linear interpolation of tabulated values in C(.).
!<------------------------------------------------------------------------------
function c_val(temp_val)

use sico_variables_m
use sico_vars_m

implicit none
real(dp)             :: c_val
real(dp), intent(in) :: temp_val

integer(i4b) :: n_temp_1, n_temp_2
real(dp) :: c_ice
#if defined(ALLOW_TAPENADE) /* Tapenade */
integer(i4b) :: itemp_val
#endif /* Tapenade */

!-------- Specific heat of pure ice --------

#if !defined(ALLOW_TAPENADE) /* Normal */

n_temp_1 = floor(temp_val)
n_temp_1 = max(min(n_temp_1, n_temp_max-1), n_temp_min)
n_temp_2 = n_temp_1 + 1

#if defined(FRAC_DUST)
c_ice = C(n_temp_1) &
        + (C(n_temp_2)-C(n_temp_1)) &
          * (temp_val-real(n_temp_1,dp))
#else
c_val = C(n_temp_1) &
        + (C(n_temp_2)-C(n_temp_1)) &
          * (temp_val-real(n_temp_1,dp))
#endif

#else /* Tapenade */

n_temp_1 = myfloor_local(temp_val)
n_temp_1 = max(min(n_temp_1, n_temp_max_IMQ-1), n_temp_min_IMQ)
n_temp_2 = n_temp_1 + 1

#if defined(FRAC_DUST)
c_ice = C_IMQ(n_temp_1) &
        + (C_IMQ(n_temp_2)-C_IMQ(n_temp_1)) &
          * (temp_val-real(n_temp_1,dp))
#else
c_val = C_IMQ(n_temp_1) &
        + (C_IMQ(n_temp_2)-C_IMQ(n_temp_1)) &
          * (temp_val-real(n_temp_1,dp))
#endif

#endif /* Normal vs. Tapenade */

!-------- If dust is present (polar caps of Mars):
!         Specific heat of ice-dust mixture --------

#if defined(FRAC_DUST)
#if !defined(ALLOW_TAPENADE) /* Normal */
c_val = rho_inv * ( (1.0_dp-FRAC_DUST)*RHO_I*c_ice + FRAC_DUST*RHO_C*C_C )
#else /* Tapenade */
c_val = rho_inv * ( (1.0_dp-FRAC_DUST)*RHO_I_IMQ*c_ice +
FRAC_DUST*RHO_C_IMQ*C_C_IMQ )
#endif /* Normal vs. Tapenade */
#endif

end function c_val

!-------------------------------------------------------------------------------
!> Creep response function for ice.
!<------------------------------------------------------------------------------
function creep(sigma_val)

use sico_variables_m
use sico_vars_m

implicit none

real(dp)             :: creep
real(dp), intent(in) :: sigma_val

#if (FLOW_LAW==4)
real(dp), parameter :: sm_coeff_1 = 8.5112e-15_dp, &   ! s^-1 Pa^-1
                       sm_coeff_2 = 8.1643e-25_dp, &   ! s^-1 Pa^-3
                       sm_coeff_3 = 9.2594e-12_dp      ! Pa^-2
#endif

#if (FIN_VISC==1)

#if (FLOW_LAW==1)

creep = sigma_val*sigma_val
!       Glen's flow law (n=3)

#elif (FLOW_LAW==2)

creep = sigma_val**0.8_dp * GR_SIZE**(-1.4_dp)
!       Goldsby-Kohlstedt flow law (n=1.8, d=1.4)

#elif (FLOW_LAW==3)

creep = sigma_val*sigma_val*sigma_val
!       Durham's flow law (n=4)

#endif

#elif (FIN_VISC==2)

#if (FLOW_LAW==1)

creep = sigma_val*sigma_val + SIGMA_RES*SIGMA_RES
!       Glen's flow (n=3) with additional finite viscosity

#elif (FLOW_LAW==2)

creep = (sigma_val**0.8_dp + SIGMA_RES**0.8_dp) * GR_SIZE**(-1.4_dp)
!       Goldsby-Kohlstedt flow law (n=1.8, d=1.4)
!       with additional finite viscosity

#elif (FLOW_LAW==3)

creep = sigma_val*sigma_val*sigma_val + SIGMA_RES*SIGMA_RES*SIGMA_RES
!       Durham's flow law (n=4) with additional finite viscosity

#endif

#endif

#if (FLOW_LAW==4)

creep = sm_coeff_1 &
        + sm_coeff_2 * (sigma_val*sigma_val) &
          * (1.0_dp + sm_coeff_3 * (sigma_val*sigma_val))
!       Smith-Morland (polynomial) flow law, normalised to a dimensionless
!       rate factor with A(-10C)=1.

#endif

end function creep

!-------------------------------------------------------------------------------
!> Ice viscosity as a function of the effective strain rate and the temperature
!! (in cold ice) or the water content (in temperate ice) or both (for the
!! enthalpy method).
!<------------------------------------------------------------------------------
function viscosity(de_val, temp_val, temp_m_val, omega_val, enh_val, &
                   i_flag_cold_temp)

use sico_variables_m
use sico_vars_m

implicit none

real(dp)                 :: viscosity
real(dp)    , intent(in) :: de_val
real(dp)    , intent(in) :: temp_val, temp_m_val
real(dp)    , intent(in) :: omega_val
real(dp)    , intent(in) :: enh_val
integer(i4b), intent(in) :: i_flag_cold_temp

real(dp) :: ratefac_val
real(dp) :: de_val_m

real(dp) :: n_power_law, inv_n_power_law, n_grain_size

real(dp), parameter :: de_min = 1.0e-30_dp   ! minimum value for the
                                             ! effective strain rate

#if (FLOW_LAW==4)
real(dp), parameter :: sm_coeff_1 = 8.5112e-15_dp, &   ! s^-1 Pa^-1
                       sm_coeff_2 = 8.1643e-25_dp, &   ! s^-1 Pa^-3
                       sm_coeff_3 = 9.2594e-12_dp      ! Pa^-2
#endif

!-------- Rate factor and effective strain rate --------

if (i_flag_cold_temp == 0) then   ! cold ice
   ratefac_val = ratefac_c(temp_val, temp_m_val)
else if (i_flag_cold_temp == 1) then   ! temperate ice
   ratefac_val = ratefac_t(omega_val)
else   ! enthalpy method
   ratefac_val = ratefac_c_t(temp_val, omega_val, temp_m_val)
end if

de_val_m = max(de_val, de_min)

#if (FIN_VISC==1)

#if (FLOW_LAW==1)

!-------- Glen's flow law (n=3) --------

inv_n_power_law = 0.333333333333333_dp   ! 1/3

viscosity = 0.5_dp * de_val_m**(inv_n_power_law-1.0_dp) &
                   * (enh_val*ratefac_val)**(-inv_n_power_law)

#elif (FLOW_LAW==2)

!-------- Goldsby-Kohlstedt flow law (n=1.8, d=1.4) --------

inv_n_power_law = 0.555555555555555_dp   ! 1/1.8
n_grain_size    = 1.4_dp

viscosity = 0.5_dp * de_val_m**(inv_n_power_law-1.0_dp) &
                   * GR_SIZE**(n_grain_size*inv_n_power_law) &
                   * (enh_val*ratefac_val)**(-inv_n_power_law)

#elif (FLOW_LAW==3)

!-------- Durham's flow law (n=4) --------

inv_n_power_law = 0.25_dp   ! 1/4

viscosity = 0.5_dp * de_val_m**(inv_n_power_law-1.0_dp) &
                   * (enh_val*ratefac_val)**(-inv_n_power_law)

#endif

#elif (FIN_VISC==2)

#if (FLOW_LAW==1)

!-------- Glen's flow (n=3) with additional finite viscosity --------

n_power_law = 3.0_dp

viscosity = visc_iter(de_val_m, ratefac_val, enh_val, n_power_law, SIGMA_RES)

#elif (FLOW_LAW==2)

!-------- Goldsby-Kohlstedt flow law (n=1.8, d=1.4)
!                           with additional finite viscosity --------

n_power_law  = 1.8_dp
n_grain_size = 1.4_dp

errormsg = ' >>> viscosity: Computation of the viscosity as a function' &
         //         end_of_line &
         //'        of the effective strain rate not yet implemented' &
         //         end_of_line &
         //'        for grain-size-dependent finite viscosity flow laws!'
call error(errormsg)

#elif (FLOW_LAW==3)

!-------- Durham's flow law (n=4) with additional finite viscosity --------

n_power_law = 4.0_dp

viscosity = visc_iter(de_val_m, ratefac_val, enh_val, n_power_law, SIGMA_RES)

#endif

#endif

#if (FLOW_LAW==4)

!-------- Smith-Morland (polynomial) flow law --------

viscosity = visc_iter_sm(de_val_m, ratefac_val, enh_val, &
                         sm_coeff_1, sm_coeff_2, sm_coeff_3)

#endif

end function viscosity

#if (FIN_VISC==2)

!-------------------------------------------------------------------------------
!> Iterative computation of the viscosity by solving equation (4.28)
!! by Greve and Blatter (Springer, 2009).
!<------------------------------------------------------------------------------
function visc_iter(de_val_m, ratefac_val, enh_val, n_power_law, sigma_res)

implicit none

real(dp)             :: visc_iter
real(dp), intent(in) :: de_val_m
real(dp), intent(in) :: ratefac_val, enh_val
real(dp), intent(in) :: n_power_law, sigma_res

integer(i4b) :: n
integer(i4b) :: max_iters
real(dp)     :: visc_val, res
logical      :: flag_rescheck1, flag_rescheck2

real(dp), parameter :: eps = 1.0e-05_dp   ! convergence parameter

!-------- Determination of the order of magnitude --------

visc_val = 1.0e+10_dp   ! initial guess (very low value)

flag_rescheck1 = .false.
n              = 0
max_iters      = 30

do while ((.not.flag_rescheck1).and.(n <= max_iters))

   n = n+1

   res = fct_visc(de_val_m, ratefac_val, enh_val, visc_val, &
                  n_power_law, sigma_res)

   if (res < 0.0_dp) then
      visc_val = 10.0_dp*visc_val
   else
      flag_rescheck1 = .true.
   end if

end do

!-------- Newton's method --------

if (flag_rescheck1) then
   ! only if order of magnitude could be detected successfully

   flag_rescheck2 = .false.
   n              = 0
   max_iters      = 1000

   do while ((.not.flag_rescheck2).and.(n <= max_iters))

      n = n+1

      visc_val = visc_val &
                 - res &
                   /fct_visc_deriv(de_val_m, ratefac_val, enh_val, visc_val, &
                                   n_power_law, sigma_res)

      res = fct_visc(de_val_m, ratefac_val, enh_val, visc_val, &
                     n_power_law, sigma_res)

      if (abs(res) < eps) then 
         flag_rescheck2 = .true. 
      end if

   end do

end if

visc_iter = visc_val

end function visc_iter

!-------------------------------------------------------------------------------
!> Viscosity polynomial
!! [equation (4.28) by Greve and Blatter (Springer, 2009)].
!<------------------------------------------------------------------------------
function fct_visc(de_val_m, ratefac_val, enh_val, visc_val, &
                  n_power_law, sigma_res)

implicit none

real(dp)             :: fct_visc
real(dp), intent(in) :: de_val_m
real(dp), intent(in) :: ratefac_val, enh_val
real(dp), intent(in) :: visc_val
real(dp), intent(in) :: n_power_law, sigma_res

fct_visc = 2.0_dp**n_power_law &
             *enh_val*ratefac_val &
             *de_val_m**(n_power_law-1.0_dp) &
             *visc_val**n_power_law &
          + 2.0_dp*enh_val*ratefac_val &
             *sigma_res**(n_power_law-1.0_dp) &
             *visc_val &
          - 1.0_dp

end function fct_visc

!-------------------------------------------------------------------------------
!> Derivative of the viscosity polynomial
!! [equation (4.28) by Greve and Blatter (Springer, 2009)].
!<------------------------------------------------------------------------------
function fct_visc_deriv(de_val_m, ratefac_val, enh_val, visc_val, &
                        n_power_law, sigma_res)

implicit none

real(dp)             :: fct_visc_deriv
real(dp), intent(in) :: de_val_m
real(dp), intent(in) :: ratefac_val, enh_val
real(dp), intent(in) :: visc_val
real(dp), intent(in) :: n_power_law, sigma_res

fct_visc_deriv = 2.0_dp**n_power_law*n_power_law &
                   *enh_val*ratefac_val &
                   *de_val_m**(n_power_law-1.0_dp) &
                   *visc_val**(n_power_law-1.0_dp) &
                 + 2.0_dp*enh_val*ratefac_val &
                   *sigma_res**(n_power_law-1.0_dp)
         
end function fct_visc_deriv

#endif /* FIN_VISC==2 */

#if (FLOW_LAW==4)

!-------------------------------------------------------------------------------
!> Iterative computation of the viscosity by solving equation (4.33)
!! [analogous to (4.28)] by Greve and Blatter (Springer, 2009).
!<------------------------------------------------------------------------------
function visc_iter_sm(de_val_m, ratefac_val, enh_val, &
                      sm_coeff_1, sm_coeff_2, sm_coeff_3)

implicit none

real(dp)             :: visc_iter_sm
real(dp), intent(in) :: de_val_m
real(dp), intent(in) :: ratefac_val, enh_val
real(dp), intent(in) :: sm_coeff_1, sm_coeff_2, sm_coeff_3

integer(i4b) :: n
integer(i4b) :: max_iters
real(dp)     :: visc_val, res
logical      :: flag_rescheck1, flag_rescheck2

real(dp), parameter :: eps = 1.0e-05_dp   ! convergence parameter

!-------- Determination of the order of magnitude --------

visc_val = 1.0e+10_dp   ! initial guess (very low value)

flag_rescheck1 = .false.
n              = 0
max_iters      = 30

do while ((.not.flag_rescheck1).and.(n <= max_iters))

   n = n+1

   res = fct_visc_sm(de_val_m, ratefac_val, enh_val, visc_val, &
                     sm_coeff_1, sm_coeff_2, sm_coeff_3)

   if (res < 0.0_dp) then
      visc_val = 10.0_dp*visc_val
   else
      flag_rescheck1 = .true.
   end if

end do

!-------- Newton's method --------

if (flag_rescheck1) then
   ! only if order of magnitude could be detected successfully

   flag_rescheck2 = .false.
   n              = 0
   max_iters      = 1000

   do while ((.not.flag_rescheck2).and.(n <= max_iters))

      n = n+1

      visc_val = visc_val &
                 - res &
                   /fct_visc_sm_deriv(de_val_m, ratefac_val, &
                                      enh_val, visc_val, &
                                      sm_coeff_1, sm_coeff_2, sm_coeff_3)

      res = fct_visc_sm(de_val_m, ratefac_val, enh_val, visc_val, &
                        sm_coeff_1, sm_coeff_2, sm_coeff_3)

      if (abs(res) < eps) then
         flag_rescheck2 = .true. 
      end if

   end do

end if

visc_iter_sm = visc_val
         
end function visc_iter_sm

!-------------------------------------------------------------------------------
!> Viscosity polynomial
!! [equation (4.33) by Greve and Blatter (Springer, 2009)].
!<------------------------------------------------------------------------------
function fct_visc_sm(de_val_m, ratefac_val, enh_val, visc_val, &
                     sm_coeff_1, sm_coeff_2, sm_coeff_3)

implicit none

real(dp)             :: fct_visc_sm
real(dp), intent(in) :: de_val_m
real(dp), intent(in) :: ratefac_val, enh_val
real(dp), intent(in) :: visc_val
real(dp), intent(in) :: sm_coeff_1, sm_coeff_2, sm_coeff_3

real(dp) :: de_visc_factor

de_visc_factor = de_val_m*de_val_m*visc_val*visc_val
                   
fct_visc_sm = 2.0_dp*enh_val*ratefac_val*visc_val &
              * ( sm_coeff_1 &
                  + 4.0_dp*sm_coeff_2*de_visc_factor &
                    * ( 1.0_dp + 4.0_dp*sm_coeff_3*de_visc_factor ) ) &
              - 1.0_dp

end function fct_visc_sm

!-------------------------------------------------------------------------------
!> Derivative of the viscosity polynomial
!! [equation (4.33) by Greve and Blatter (Springer, 2009)].
!<------------------------------------------------------------------------------
function fct_visc_sm_deriv(de_val_m, ratefac_val, enh_val, visc_val, &
                           sm_coeff_1, sm_coeff_2, sm_coeff_3)

implicit none

real(dp)             :: fct_visc_sm_deriv
real(dp), intent(in) :: de_val_m
real(dp), intent(in) :: ratefac_val, enh_val
real(dp), intent(in) :: visc_val
real(dp), intent(in) :: sm_coeff_1, sm_coeff_2, sm_coeff_3

real(dp) :: de_visc_factor

real(dp), parameter :: twenty_over_three = 6.666666666666667_dp
                          
de_visc_factor = de_val_m*de_val_m*visc_val*visc_val
                          
fct_visc_sm_deriv = 2.0_dp*sm_coeff_1*enh_val*ratefac_val &
                   + 24.0_dp*sm_coeff_2*enh_val*ratefac_val*de_visc_factor &
                     * ( 1.0_dp + twenty_over_three*sm_coeff_3*de_visc_factor )
         
end function fct_visc_sm_deriv

#endif /* FLOW_LAW==4 */

!-------------------------------------------------------------------------------

#if defined(ALLOW_TAPENADE) /* Tapenade */

  function myfloor_local(num)
  
  use sico_variables_m
  
  implicit none

  real(dp), intent(in) :: num
  
  integer(i4b) :: inum
  integer(i4b) :: myfloor_local
  
  inum = int(num);
  
  if (abs(num-real(inum,dp)) <= &
    (abs(num+real(inum,dp))*myepsilon_dp) ) then
    myfloor_local = inum;
  else if (num>0) then
    myfloor_local = inum;
  else if (num<0) then
    myfloor_local = inum - 1;
  end if
  
  end function myfloor_local

#endif /* Tapenade */

!-------------------------------------------------------------------------------

end module ice_material_properties_m
!
