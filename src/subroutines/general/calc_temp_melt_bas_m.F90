!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ t e m p _ m e l t _ b a s _ m
!
!> @file
!!
!! Computation of the melting and basal temperatures.
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
!> Computation of the melting and basal temperatures.
!<------------------------------------------------------------------------------
module calc_temp_melt_bas_m

  use sico_types_m
  use sico_variables_m
  use sico_vars_m

  implicit none

  private
  public :: calc_temp_melt, calc_temp_bas

contains

!-------------------------------------------------------------------------------
!> Computation of the melting temperatures.
!<------------------------------------------------------------------------------
  subroutine calc_temp_melt()

  implicit none
  integer(i4b) :: i, j, kc, kt
  real(dp), dimension(0:KCMAX) :: atm1
  real(dp), dimension(0:KTMAX) :: atm2

!-------- Term abbreviations --------

  atm1 = BETA*(1.0_dp-eaz_c_quotient)
  atm2 = BETA*(1.0_dp-zeta_t)

!-------- Compute the melting temperatures --------

  do i=0, IMAX
  do j=0, JMAX

     do kt=0, KTMAX
        temp_t_m(kt,j,i) = -(BETA*H_c(j,i)+atm2(kt)*H_t(j,i))
     end do

     do kc=0, KCMAX
        temp_c_m(kc,j,i) = -atm1(kc)*H_c(j,i)
     end do

  end do
  end do

  end subroutine calc_temp_melt

!-------------------------------------------------------------------------------
!> Computation of the basal temperatures.
!<------------------------------------------------------------------------------
  subroutine calc_temp_bas()

  implicit none
  integer(i4b) :: i, j

!-------- Computation of the basal temperatures --------

  do i=0, IMAX
  do j=0, JMAX

     if ( (mask(j,i) == 0).or.(mask(j,i) == 3) ) then
                                   ! glaciated land or floating ice

        if (n_cts(j,i) == -1) then   ! cold ice base

           temp_b(j,i)  = temp_c(0,j,i)
           temph_b(j,i) = temp_c(0,j,i) - temp_c_m(0,j,i)
                          ! relative to the pressure melting point

        else   ! n_cts(j,i) == 0 or 1, temperate ice base

           temp_b(j,i)  = temp_t_m(0,j,i)
           temph_b(j,i) = 0.0_dp
                          ! relative to the pressure melting point

        end if

     else   ! mask(j,i) == 1 or 2, ice-free land or sea

        temp_b(j,i)  = temp_c(0,j,i)
        temph_b(j,i) = temp_c(0,j,i) - temp_c_m(0,j,i)
                       ! relative to the pressure melting point

     end if

  end do
  end do

  end subroutine calc_temp_bas

!-------------------------------------------------------------------------------

end module calc_temp_melt_bas_m
!
