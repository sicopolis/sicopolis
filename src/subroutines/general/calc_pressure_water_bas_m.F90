!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ p r e s s u r e _ w a t e r _ b a s _ m
!
!! Computation of the basal water pressure.
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
!> Computation of the basal water pressure.
!-------------------------------------------------------------------------------
module calc_pressure_water_bas_m

  use sico_types_m
  use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

  use error_m

  implicit none

  public

contains

!-------------------------------------------------------------------------------
!> Main subroutine of calc_pressure_water_bas_m:
!! Computation of the basal water pressure.
!-------------------------------------------------------------------------------
  subroutine calc_pressure_water_bas()

  implicit none

  integer(i4b) :: i, j, ij

!-------- Basal water pressure --------

  do ij=1, (IMAX+1)*(JMAX+1)

     i = n2i(ij)   ! i=0...IMAX
     j = n2j(ij)   ! j=0...JMAX

#if (!defined(BASAL_WATER_PRESSURE) || BASAL_WATER_PRESSURE==0)

     p_b_w(j,i) = 0.0_dp
                   ! zero everywhere

#elif (BASAL_WATER_PRESSURE==1)

     p_b_w(j,i) = RHO_SW*G*(z_sl(j,i)-zb(j,i))
                   ! ocean pressure without cut-off (can become negative)

#elif (BASAL_WATER_PRESSURE==2)

     p_b_w(j,i) = RHO_SW*G*max((z_sl(j,i)-zb(j,i)), 0.0_dp)
                   ! ocean pressure with cut-off

#else

     errormsg = ' >>> calc_pressure_water_bas: ' &
              // 'Parameter BASAL_WATER_PRESSURE must be 0, 1 or 2!'
     call error(errormsg)

#endif

  end do

  end subroutine calc_pressure_water_bas

!-------------------------------------------------------------------------------

end module calc_pressure_water_bas_m
!
