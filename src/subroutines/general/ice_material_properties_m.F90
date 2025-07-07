!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  i c e _ m a t e r i a l _ p r o p e r t i e s _ m
!
!! Material properties of ice:
!! Rate factor, heat conductivity, specific heat (heat capacity),
!! creep function, viscosity.
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
!> Material properties of ice:
!! Rate factor, heat conductivity, specific heat (heat capacity),
!! creep function, viscosity.
!-------------------------------------------------------------------------------
module ice_material_properties_m

use sico_types_m
use sico_variables_m, only : RF_imp, R_T_imp, KAPPA_imp, C_imp, &
                             n_temp_min_imp, n_temp_max_imp, &
                             RHO_I_imp, RHO_C_imp, KAPPA_C_imp, C_C_imp
use error_m

implicit none

#if !defined(ALLOW_TAPENADE) /* Normal */

private
public :: ice_mat_eqs_pars, &
          ratefac_c, ratefac_t, ratefac_c_t, kappa_val, c_val, &
          viscosity, creep

#else /* Tapenade */

public :: ice_mat_eqs_pars, &
          ratefac_c, ratefac_t, ratefac_c_t, kappa_val, c_val, &
          viscosity, creep
private :: myfloor_local

#endif /* Normal vs. Tapenade */

contains

!-------------------------------------------------------------------------------
!> Setting of required physical parameters.
!-------------------------------------------------------------------------------
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

!-------- Initialization --------

RF_imp    = 0.0_dp
KAPPA_imp = 0.0_dp
C_imp     = 0.0_dp
n_temp_min_imp = n_tmp_min
n_temp_max_imp = n_tmp_max

if ((n_temp_min_imp <= -256).or.(n_temp_max_imp >= 255)) then
   errormsgg = ' >>> ice_mat_eqs_pars: ' &
                  //'Temperature indices out of allowed range!'
   call error(errormsgg)
end if

!-------- Assignment --------

do n=n_temp_min_imp, n_temp_max_imp
   RF_imp(n)    = RF_table(n)
   KAPPA_imp(n) = KAPPA_table(n)
   C_imp(n)     = C_table(n)
end do

do n=-256, n_temp_min_imp-1
   RF_imp(n)    = RF_imp(n_temp_min_imp)      ! dummy values
   KAPPA_imp(n) = KAPPA_imp(n_temp_min_imp)   ! dummy values
   C_imp(n)     = C_imp(n_temp_min_imp)       ! dummy values
end do

do n=n_temp_max_imp+1, 255
   RF_imp(n)    = RF_imp(n_temp_max_imp)      ! dummy values
   KAPPA_imp(n) = KAPPA_imp(n_temp_max_imp)   ! dummy values
   C_imp(n)     = C_imp(n_temp_max_imp)       ! dummy values
end do

R_T_imp = R_T_val

!-------- Martian stuff --------

if ( present(RHO_I_val) ) then
   RHO_I_imp = RHO_I_val
else
   RHO_I_imp = 0.0_dp   ! dummy value
end if

if ( present(RHO_C_val) ) then
   RHO_C_imp = RHO_C_val
else
   RHO_C_imp = 0.0_dp   ! dummy value
end if

if ( present(KAPPA_C_val) ) then
   KAPPA_C_imp = KAPPA_C_val
else
   KAPPA_C_imp = 0.0_dp   ! dummy value
end if

if ( present(C_C_val) ) then
   C_C_imp = C_C_val
else
   C_C_imp = 0.0_dp   ! dummy value
end if

end subroutine ice_mat_eqs_pars

!-------------------------------------------------------------------------------
!> Rate factor for cold ice:
!! Linear interpolation of tabulated values in RF(.).
!-------------------------------------------------------------------------------
function ratefac_c(temp_val, temp_m_val)

use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

implicit none
real(dp)             :: ratefac_c
real(dp), intent(in) :: temp_val, temp_m_val

integer(i4b) :: n_temp_1, n_temp_2
real(dp)     :: temp_h_val

temp_h_val = temp_val-temp_m_val

#if !defined(ALLOW_TAPENADE) /* Normal */
n_temp_1 = floor(temp_h_val)
#else /* Tapenade */
n_temp_1 = myfloor_local(temp_h_val)
#endif /* Normal vs. Tapenade */

n_temp_1 = max(min(n_temp_1, n_temp_max_imp-1), n_temp_min_imp)
n_temp_2 = n_temp_1 + 1

ratefac_c = RF_imp(n_temp_1) &
              + (RF_imp(n_temp_2)-RF_imp(n_temp_1)) &
                * (temp_h_val-real(n_temp_1,dp))   ! Linear interpolation

end function ratefac_c

!-------------------------------------------------------------------------------
!> Rate factor for temperate ice.
!-------------------------------------------------------------------------------
function ratefac_t(omega_val)

use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

implicit none
real(dp)             :: ratefac_t
real(dp), intent(in) :: omega_val

ratefac_t = RF_imp(0)*(1.0_dp+R_T_imp*(omega_val))

end function ratefac_t

!-------------------------------------------------------------------------------
!> Rate factor for cold and temperate ice:
!! Combination of ratefac_c and ratefac_t (only for the enthalpy method).
!-------------------------------------------------------------------------------
function ratefac_c_t(temp_val, omega_val, temp_m_val)

use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

implicit none
real(dp)             :: ratefac_c_t
real(dp), intent(in) :: temp_val, temp_m_val, omega_val

integer(i4b) :: n_temp_1, n_temp_2
real(dp)     :: temp_h_val

temp_h_val = temp_val-temp_m_val

#if !defined(ALLOW_TAPENADE) /* Normal */
n_temp_1 = floor(temp_h_val)
#else /* Tapenade */
n_temp_1 = myfloor_local(temp_h_val)
#endif /* Normal vs. Tapenade */

n_temp_1 = max(min(n_temp_1, n_temp_max_imp-1), n_temp_min_imp)
n_temp_2 = n_temp_1 + 1

ratefac_c_t = ( RF_imp(n_temp_1) &
                  + (RF_imp(n_temp_2)-RF_imp(n_temp_1)) &
                    * (temp_h_val-real(n_temp_1,dp)) ) &
              *(1.0_dp+R_T_imp*(omega_val))   ! Linear interpolation

end function ratefac_c_t

!-------------------------------------------------------------------------------
!> Heat conductivity of ice:
!! Linear interpolation of tabulated values in KAPPA(.).
!-------------------------------------------------------------------------------
function kappa_val(temp_val)

use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

implicit none
real(dp)             :: kappa_val
real(dp), intent(in) :: temp_val

integer(i4b) :: n_temp_1, n_temp_2
real(dp) :: kappa_ice

!-------- Heat conductivity of pure ice --------

#if !defined(ALLOW_TAPENADE) /* Normal */
n_temp_1 = floor(temp_val)
#else /* Tapenade */
n_temp_1 = myfloor_local(temp_val)
#endif /* Normal vs. Tapenade */

n_temp_1 = max(min(n_temp_1, n_temp_max_imp-1), n_temp_min_imp)
n_temp_2 = n_temp_1 + 1

#if defined(FRAC_DUST)
kappa_ice = KAPPA_imp(n_temp_1) &
            + (KAPPA_imp(n_temp_2)-KAPPA_imp(n_temp_1)) &
              * (temp_val-real(n_temp_1,dp))   ! Linear interpolation
#else
kappa_val = KAPPA_imp(n_temp_1) &
            + (KAPPA_imp(n_temp_2)-KAPPA_imp(n_temp_1)) &
              * (temp_val-real(n_temp_1,dp))   ! Linear interpolation
#endif

!-------- If dust is present (polar caps of Mars):
!         Heat conductivity of ice-dust mixture --------

#if defined(FRAC_DUST)
kappa_val = (1.0_dp-FRAC_DUST)*kappa_ice + FRAC_DUST*KAPPA_C_imp
#endif

end function kappa_val

!-------------------------------------------------------------------------------
!> Specific heat of ice:
!! Linear interpolation of tabulated values in C(.).
!-------------------------------------------------------------------------------
function c_val(temp_val)

use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

implicit none
real(dp)             :: c_val
real(dp), intent(in) :: temp_val

integer(i4b) :: n_temp_1, n_temp_2
real(dp) :: c_ice

!-------- Specific heat of pure ice --------

#if !defined(ALLOW_TAPENADE) /* Normal */
n_temp_1 = floor(temp_val)
#else /* Tapenade */
n_temp_1 = myfloor_local(temp_val)
#endif /* Normal vs. Tapenade */

n_temp_1 = max(min(n_temp_1, n_temp_max_imp-1), n_temp_min_imp)
n_temp_2 = n_temp_1 + 1

#if defined(FRAC_DUST)
c_ice = C_imp(n_temp_1) &
        + (C_imp(n_temp_2)-C_imp(n_temp_1)) &
          * (temp_val-real(n_temp_1,dp))   ! Linear interpolation
#else
c_val = C_imp(n_temp_1) &
        + (C_imp(n_temp_2)-C_imp(n_temp_1)) &
          * (temp_val-real(n_temp_1,dp))   ! Linear interpolation
#endif

!-------- If dust is present (polar caps of Mars):
!         Specific heat of ice-dust mixture --------

#if defined(FRAC_DUST)
c_val = rho_inv &
        * ( (1.0_dp-FRAC_DUST)*RHO_I_imp*c_ice + FRAC_DUST*RHO_C_imp*C_C_imp )
#endif

end function c_val

!-------------------------------------------------------------------------------
!> Creep response function for ice.
!-------------------------------------------------------------------------------
function creep(sigma_val)

use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

implicit none

real(dp)             :: creep
real(dp), intent(in) :: sigma_val

#if (FLOW_LAW==1)
real(dp) :: d_n_power_law
#elif (FLOW_LAW==4)
real(dp) :: sm_coeff_1, sm_coeff_2, sm_coeff_3
#endif

#if (FLOW_LAW==1)

#if (defined(N_POWER_LAW))
d_n_power_law = real(N_POWER_LAW,dp)
#else
d_n_power_law = 3.0_dp   ! default n=3
#endif

#if (FIN_VISC==1)
creep = sigma_val**(d_n_power_law-1.0_dp)
        ! Nye-Glen flow law
#elif (FIN_VISC==2)
creep = sigma_val**(d_n_power_law-1.0_dp) + SIGMA_RES**(d_n_power_law-1.0_dp)
        ! Nye-Glen flow law with additional finite viscosity
#endif

#elif (FLOW_LAW==4)

#if (defined(SM_COEFF1))
sm_coeff_1 = SM_COEFF1
#else
sm_coeff_1 = 8.5112e-15_dp   ! s^-1 Pa^-1
#endif

#if (defined(SM_COEFF2))
sm_coeff_2 = SM_COEFF2
#else
sm_coeff_2 = 8.1643e-25_dp   ! s^-1 Pa^-3
#endif

#if (defined(SM_COEFF3))
sm_coeff_3 = SM_COEFF3
#else
sm_coeff_3 = 9.2594e-12_dp   ! Pa^-2
#endif

creep = sm_coeff_1 &
        + sm_coeff_2 * (sigma_val*sigma_val) &
          * (1.0_dp + sm_coeff_3 * (sigma_val*sigma_val))
        ! Smith-Morland (polynomial) flow law,
        ! normalized to a dimensionless rate factor with A(-10C)=1

#endif

end function creep

!-------------------------------------------------------------------------------
!> Ice viscosity as a function of the effective strain rate and the temperature
!! (in cold ice) or the water content (in temperate ice) or both (for the
!! enthalpy method).
!-------------------------------------------------------------------------------
function viscosity(de_val, temp_val, temp_m_val, omega_val, enh_val, &
                   i_flag_cold_temp)

use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

implicit none

real(dp)                 :: viscosity
real(dp)    , intent(in) :: de_val
real(dp)    , intent(in) :: temp_val, temp_m_val
real(dp)    , intent(in) :: omega_val
real(dp)    , intent(in) :: enh_val
integer(i4b), intent(in) :: i_flag_cold_temp

real(dp) :: ratefac_val
real(dp) :: de_val_m

#if (FLOW_LAW==1)
real(dp) :: d_n_power_law, d_inv_n_power_law
#elif (FLOW_LAW==4)
real(dp) :: sm_coeff_1, sm_coeff_2, sm_coeff_3
#endif

real(dp), parameter :: de_min = 1.0e-30_dp   ! minimum value for the
                                             ! effective strain rate

!-------- Rate factor and effective strain rate --------

if (i_flag_cold_temp == 0) then   ! cold ice
   ratefac_val = ratefac_c(temp_val, temp_m_val)
else if (i_flag_cold_temp == 1) then   ! temperate ice
   ratefac_val = ratefac_t(omega_val)
else   ! enthalpy method
   ratefac_val = ratefac_c_t(temp_val, omega_val, temp_m_val)
end if

de_val_m = max(de_val, de_min)

#if (FLOW_LAW==1)

#if (FIN_VISC==1)

#if (defined(N_POWER_LAW))
d_inv_n_power_law = 1.0_dp/real(N_POWER_LAW,dp)
#else
d_inv_n_power_law = 1.0_dp/3.0_dp   ! default n=3
#endif

viscosity = 0.5_dp * de_val_m**(d_inv_n_power_law-1.0_dp) &
                   * (enh_val*ratefac_val)**(-d_inv_n_power_law)
            ! Nye-Glen flow law

#elif (FIN_VISC==2)

#if (defined(N_POWER_LAW))
d_n_power_law = real(N_POWER_LAW,dp)
#else
d_n_power_law = 3.0_dp   ! default n=3
#endif

viscosity = visc_iter(de_val_m, ratefac_val, enh_val, d_n_power_law, SIGMA_RES)
            ! Nye-Glen flow with additional finite viscosity

#endif

#elif (FLOW_LAW==4)

#if (defined(SM_COEFF1))
sm_coeff_1 = SM_COEFF1
#else
sm_coeff_1 = 8.5112e-15_dp   ! s^-1 Pa^-1
#endif

#if (defined(SM_COEFF2))
sm_coeff_2 = SM_COEFF2
#else
sm_coeff_2 = 8.1643e-25_dp   ! s^-1 Pa^-3
#endif

#if (defined(SM_COEFF3))
sm_coeff_3 = SM_COEFF3
#else
sm_coeff_3 = 9.2594e-12_dp   ! Pa^-2
#endif

viscosity = visc_iter_sm(de_val_m, ratefac_val, enh_val, &
                         sm_coeff_1, sm_coeff_2, sm_coeff_3)
            ! Smith-Morland (polynomial) flow law,
            ! normalized to a dimensionless rate factor with A(-10C)=1

#endif

end function viscosity

#if (FIN_VISC==2)

!-------------------------------------------------------------------------------
!> Iterative computation of the viscosity by solving equation (4.28)
!! by Greve and Blatter (Springer, 2009).
!-------------------------------------------------------------------------------
function visc_iter(de_val_m, ratefac_val, enh_val, d_n_power_law, sigma_resid)

implicit none

real(dp)             :: visc_iter
real(dp), intent(in) :: de_val_m
real(dp), intent(in) :: ratefac_val, enh_val
real(dp), intent(in) :: d_n_power_law, sigma_resid

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
                  d_n_power_law, sigma_resid)

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
                                   d_n_power_law, sigma_resid)

      res = fct_visc(de_val_m, ratefac_val, enh_val, visc_val, &
                     d_n_power_law, sigma_resid)

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
!-------------------------------------------------------------------------------
function fct_visc(de_val_m, ratefac_val, enh_val, visc_val, &
                  d_n_power_law, sigma_resid)

implicit none

real(dp)             :: fct_visc
real(dp), intent(in) :: de_val_m
real(dp), intent(in) :: ratefac_val, enh_val
real(dp), intent(in) :: visc_val
real(dp), intent(in) :: d_n_power_law, sigma_resid

fct_visc = 2.0_dp**d_n_power_law &
             *enh_val*ratefac_val &
             *de_val_m**(d_n_power_law-1.0_dp) &
             *visc_val**d_n_power_law &
          + 2.0_dp*enh_val*ratefac_val &
             *sigma_resid**(d_n_power_law-1.0_dp) &
             *visc_val &
          - 1.0_dp

end function fct_visc

!-------------------------------------------------------------------------------
!> Derivative of the viscosity polynomial
!! [equation (4.28) by Greve and Blatter (Springer, 2009)].
!-------------------------------------------------------------------------------
function fct_visc_deriv(de_val_m, ratefac_val, enh_val, visc_val, &
                        d_n_power_law, sigma_resid)

implicit none

real(dp)             :: fct_visc_deriv
real(dp), intent(in) :: de_val_m
real(dp), intent(in) :: ratefac_val, enh_val
real(dp), intent(in) :: visc_val
real(dp), intent(in) :: d_n_power_law, sigma_resid

fct_visc_deriv = 2.0_dp**d_n_power_law*d_n_power_law &
                   *enh_val*ratefac_val &
                   *de_val_m**(d_n_power_law-1.0_dp) &
                   *visc_val**(d_n_power_law-1.0_dp) &
                 + 2.0_dp*enh_val*ratefac_val &
                   *sigma_resid**(d_n_power_law-1.0_dp)
         
end function fct_visc_deriv

#endif /* FIN_VISC==2 */

#if (FLOW_LAW==4)

!-------------------------------------------------------------------------------
!> Iterative computation of the viscosity by solving equation (4.33)
!! [analogous to (4.28)] by Greve and Blatter (Springer, 2009).
!-------------------------------------------------------------------------------
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
!-------------------------------------------------------------------------------
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
!-------------------------------------------------------------------------------
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
