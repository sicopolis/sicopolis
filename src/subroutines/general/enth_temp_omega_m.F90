!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  e n t h _ t e m p _ o m e g a _ m
!
!! Conversion from temperature (temp) and water content (omega) to enthalpy
!! (enth) and vice versa.
!!
!!##### Authors
!!
!! Ralf Greve, Heinz Blatter
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
!> Conversion from temperature (temp) and water content (omega) to enthalpy
!! (enth) and vice versa.
!-------------------------------------------------------------------------------
module enth_temp_omega_m

use sico_types_m
use sico_variables_m, only : c_int_table, c_int_inv_table, &
                             n_temp_min, n_temp_max, n_enth_min, n_enth_max, &
                             latent_heat, latent_heat_inv
use error_m

implicit none

public :: calc_c_int_table, calc_c_int_inv_table
public :: enth_fct_temp_omega, temp_fct_enth, omega_fct_enth
public :: c_int_val, c_int_inv_val

#if defined(ALLOW_TAPENADE)   /* Tapenade */
private :: myfloor_local_ETO
#endif

contains

!-------------------------------------------------------------------------------
!> Computation of the temperature integral of the specific heat of ice as a
!! table (c_int_table). Further, definition of the latent heat of ice.
!-------------------------------------------------------------------------------
subroutine calc_c_int_table(c_table, n_tmp_min, n_tmp_max, latent_heat_val)

implicit none

integer(i4b),                             intent(in) :: n_tmp_min, n_tmp_max
real(dp), dimension(n_tmp_min:n_tmp_max), intent(in) :: c_table
real(dp), optional,                       intent(in) :: latent_heat_val

integer(i4b)       :: n
real(dp)           :: c_int_zero
character(len=256) :: errormsgg

!-------- Initialization --------

c_int_table = 0.0_dp

n_temp_min = n_tmp_min
n_temp_max = n_tmp_max

if ((n_temp_min <= -256).or.(n_temp_max >= 255)) then
   errormsgg = ' >>> calc_c_int_table: ' &
                  //'Temperature indices out of allowed range!'
   call error(errormsgg)
end if

!-------- Numerical integration with the trapezoidal rule (spacing
!         of data in c_table and c_int_table assumed to be 1 degC) --------

do n=n_temp_min+1, n_temp_max
   c_int_table(n) = c_int_table(n-1) + 0.5_dp*(c_table(n-1)+c_table(n))
                    ! that's the real stuff
end do

do n=n_temp_max+1, 255
   c_int_table(n) = c_int_table(n_temp_max)   ! dummy values
end do

!-------- Shift of the zero level to 0 degC --------

c_int_zero = c_int_table(0)

do n=-256, 255
   c_int_table(n) = c_int_table(n) - c_int_zero
end do

!-------- Further stuff: Latent heat of ice --------

if ( present(latent_heat_val) ) then
   latent_heat = latent_heat_val
else
   latent_heat = 3.35e+05_dp   ! in J/kg
end if

latent_heat_inv = 1.0_dp/latent_heat

end subroutine calc_c_int_table

!-------------------------------------------------------------------------------
!> Computation of the inverse of the temperature integral of the specific heat
!! of ice as a table (c_int_inv_table).
!-------------------------------------------------------------------------------
subroutine calc_c_int_inv_table()

#if defined(ALLOW_TAPENADE) /* Tapenade */
  use ctrl_m, only: myceiling, myfloor
#endif /* Tapenade */

implicit none

integer(i4b)       :: n
integer(i4b)       :: n_temp_1, n_temp_2
real(dp)           :: enth_min, enth_max
real(dp)           :: enth_val, enth_1, enth_2
character(len=256) :: errormsgg

!-------- Initialization --------

c_int_inv_table = 0.0_dp

enth_min = c_int_val(real(n_temp_min,dp))
enth_max = c_int_val(real(n_temp_max,dp))

#if !defined(ALLOW_TAPENADE) /* Normal */
n_enth_min = ceiling(enth_min)
#else /* Tapenade */
call myceiling(enth_min,n_enth_min)
#endif /* Normal vs. Tapenade */

#if !defined(ALLOW_TAPENADE) /* Normal */
n_enth_max = floor(enth_max)
#else /* Tapenade */
call myfloor(enth_max, n_enth_max)
#endif /* Normal vs. Tapenade */

if ((n_enth_min <= -524288).or.(n_enth_max >= 524287)) then
   errormsgg = ' >>> calc_c_int_inv_table: ' &
                  //'Enthalpy indices out of allowed range!'
   call error(errormsgg)
end if

!-------- Linear interpolation between adjacent enthalpy values --------

n_temp_1 = n_temp_min
n_temp_2 = n_temp_min+1

do n=n_enth_min, n_enth_max

   enth_val = real(n,dp)

#if !defined(ALLOW_TAPENADE) /* Normal */

   do

      if ((n_temp_1 > n_temp_max).or.(n_temp_2 > n_temp_max)) then
         errormsgg = ' >>> calc_c_int_inv_table: ' &
                        //'Temperature indices out of allowed range!'
         call error(errormsgg)
      end if

      enth_1 = c_int_val(real(n_temp_1,dp))
      enth_2 = c_int_val(real(n_temp_2,dp))

      if ( (enth_1 <= enth_val).and.(enth_2 >= enth_val) ) exit

      n_temp_1 = n_temp_1+1
      n_temp_2 = n_temp_2+1

   end do

#else /* Tapenade */

   enth_1 = c_int_val(real(n_temp_1,dp))
   enth_2 = c_int_val(real(n_temp_2,dp))

   do while ( .not.((enth_1 <= enth_val).and.(enth_2 >= enth_val)) ) 
   
      if ((n_temp_1 > n_temp_max).or.(n_temp_2 > n_temp_max)) then
         errormsgg = ' >>> calc_c_int_inv_table: ' &
                        //'Temperature indices out of allowed range!'
         call error(errormsgg)
      end if

      enth_1 = c_int_val(real(n_temp_1,dp))
      enth_2 = c_int_val(real(n_temp_2,dp))

      if ( .not.((enth_1 <= enth_val).and.(enth_2 >= enth_val)) ) then
         n_temp_1 = n_temp_1+1
         n_temp_2 = n_temp_2+1
      end if 

   end do

#endif /* Normal vs. Tapenade */

   c_int_inv_table(n) = real(n_temp_1,dp) &
                        + (real(n_temp_2,dp)-real(n_temp_1,dp)) &
                          * (enth_val-enth_1)/(enth_2-enth_1)
                                              ! Linear interpolation
end do

do n=-524288, n_enth_min-1
   c_int_inv_table(n) = c_int_inv_table(n_enth_min)   ! dummy values
end do

do n=n_enth_max+1, 524287
   c_int_inv_table(n) = c_int_inv_table(n_enth_max)   ! dummy values
end do

end subroutine calc_c_int_inv_table

!-------------------------------------------------------------------------------
!> Temperature integral of the specific heat of ice
!! (enthalpy as function of temperature).
!-------------------------------------------------------------------------------
function c_int_val(temp_val)

implicit none

real(dp)             :: c_int_val

real(dp), intent(in) :: temp_val

integer(i4b) :: n_temp_1, n_temp_2

! character(len=256) :: errormsgg

#if !defined(ALLOW_TAPENADE) /* Normal */
n_temp_1 = floor(temp_val)
#else /* Tapenade */
n_temp_1 = myfloor_local_ETO(temp_val)
#endif /* Normal vs. Tapenade */

n_temp_1 = max(min(n_temp_1, n_temp_max-1), n_temp_min)
n_temp_2 = n_temp_1 + 1

! if ((n_temp_1 < n_temp_min-1).or.(n_temp_2 > n_temp_max+1)) then
!    errormsgg = ' >>> c_int_val: Temperature argument out of allowed range!'
!    call error(errormsgg)
! end if
!    *** Commented out after some testing in order to save computing time. ***

c_int_val = c_int_table(n_temp_1) &
            + (c_int_table(n_temp_2)-c_int_table(n_temp_1)) &
              * (temp_val-real(n_temp_1,dp))   ! Linear interpolation

end function c_int_val

!-------------------------------------------------------------------------------
!> Inverse function of c_int_val (temperature as function of enthalpy).
!-------------------------------------------------------------------------------
function c_int_inv_val(enth_val)

implicit none

real(dp)             :: c_int_inv_val

real(dp), intent(in) :: enth_val

integer(i4b) :: n_enth_1, n_enth_2

! character(len=256) :: errormsgg

#if !defined(ALLOW_TAPENADE) /* Normal */
n_enth_1 = floor(enth_val)
#else /* Tapenade */
n_enth_1 = myfloor_local_ETO(enth_val)
#endif /* Normal vs. Tapenade */

n_enth_1 = max(min(n_enth_1, n_enth_max-1), n_enth_min)
n_enth_2 = n_enth_1 + 1

! if ((n_enth_1 < n_enth_min-1).or.(n_enth_2 > n_enth_max+1)) then
!    errormsgg = ' >>> c_int_inv_val: Enthalpy argument out of allowed range!'
!    call error(errormsgg)
! end if
!    *** Commented out after some testing in order to save computing time. ***

c_int_inv_val = c_int_inv_table(n_enth_1) &
                + (c_int_inv_table(n_enth_2)-c_int_inv_table(n_enth_1)) &
                  * (enth_val-real(n_enth_1,dp))   ! Linear interpolation

end function c_int_inv_val

!-------------------------------------------------------------------------------
!> Enthalpy as a function of temperature and water content.
!-------------------------------------------------------------------------------
function enth_fct_temp_omega(temp_val, omega_val)

implicit none

real(dp)              :: enth_fct_temp_omega

real(dp), intent(in)  :: temp_val, omega_val

enth_fct_temp_omega = c_int_val(temp_val) + latent_heat*omega_val

end function enth_fct_temp_omega

!-------------------------------------------------------------------------------
!> Temperature as a function of enthalpy.
!-------------------------------------------------------------------------------
function temp_fct_enth(enth_val, temp_m_val)

implicit none

real(dp)             :: temp_fct_enth

real(dp), intent(in) :: enth_val
real(dp), intent(in) :: temp_m_val

real(dp) :: enth_i

enth_i = c_int_val(temp_m_val)   ! Enthalpy of pure ice at the melting point

if (enth_val < enth_i) then   ! cold ice
   temp_fct_enth = c_int_inv_val(enth_val)
else   ! temperate ice
   temp_fct_enth = temp_m_val
end if

end function temp_fct_enth

!-------------------------------------------------------------------------------
!> Water content as a function of enthalpy.
!-------------------------------------------------------------------------------
function omega_fct_enth(enth_val, temp_m_val)

implicit none

real(dp)             :: omega_fct_enth

real(dp), intent(in) :: enth_val
real(dp), intent(in) :: temp_m_val

real(dp) :: enth_i

enth_i = c_int_val(temp_m_val)   ! Enthalpy of pure ice at the melting point

omega_fct_enth = max((enth_val-enth_i)*latent_heat_inv, 0.0_dp)

end function omega_fct_enth

!-------------------------------------------------------------------------------

#if defined(ALLOW_TAPENADE) /* Tapenade */

  function myfloor_local_ETO(num)

  use sico_variables_m

  implicit none

  integer(i4b)         :: inum
  real(dp), intent(in) :: num
  integer(i4b) :: myfloor_local_ETO

  inum = int(num);

  if (abs(num-real(inum,dp)) <= &
    (abs(num+real(inum,dp))*myepsilon_dp) ) then
    myfloor_local_ETO = inum;
  else if (num>0) then
    myfloor_local_ETO = inum;
  else if (num<0) then
    myfloor_local_ETO = inum - 1;
  end if

  end function myfloor_local_ETO

#endif /* Tapenade */

!-------------------------------------------------------------------------------

end module enth_temp_omega_m
!
