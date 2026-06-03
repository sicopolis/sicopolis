!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ t e m p _ a u x _ m
!
!! Computation of the melting temperature, basal temperature, vertical
!! temperature gradient at the base and mean (depth-averaged) temperature.
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
!> Computation of the melting temperature, basal temperature, vertical
!! temperature gradient at the base and mean (depth-averaged) temperature.
!-------------------------------------------------------------------------------
module calc_temp_aux_m

  use sico_types_m
  use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

  implicit none

  public

contains

!-------------------------------------------------------------------------------
!> Computation of the melting temperature.
!-------------------------------------------------------------------------------
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
!> Computation of the basal temperature.
!-------------------------------------------------------------------------------
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
!> Computation of the vertical temperature gradient at the base.
!-------------------------------------------------------------------------------
  subroutine calc_temp_bas_grad(dzeta_c, dzeta_t)

  implicit none

  real(dp), intent(in) :: dzeta_c, dzeta_t

  integer(i4b) :: i, j
  real(dp)     :: fct_c

  if (flag_aa_nonzero) then
     fct_c = (ea-1.0_dp)/aa
  else
     fct_c = 1.0_dp
  end if

  do i=0, IMAX
  do j=0, JMAX

     if ( (mask(j,i) == 0).or.(mask(j,i) == 3) ) then
                                   ! glaciated land or floating ice

        if ( (n_cts(j,i) == -1).or.(n_cts(j,i) == 0) ) then

           dtemp_dz_b(j,i) = fct_c*(temp_c(1,j,i)-temp_c(0,j,i)) &
                                   / (H_c(j,i)*dzeta_c)

        else   ! n_cts(j,i) == 1, temperate ice layer in kt domain

           dtemp_dz_b(j,i) = (temp_t_m(1,j,i)-temp_t_m(0,j,i)) &
                             / (H_t(j,i)*dzeta_t)

        end if

     else   ! mask(j,i) == 1 or 2, ice-free land or sea

        dtemp_dz_b(j,i) = 0.0_dp

     end if

  end do
  end do

  end subroutine calc_temp_bas_grad

!-------------------------------------------------------------------------------
!> Computation of the mean (depth-averaged) temperature.
!-------------------------------------------------------------------------------
  subroutine calc_temp_mean(dzeta_c, dzeta_t)

  implicit none

  real(dp), intent(in) :: dzeta_c, dzeta_t

  integer(i4b) :: i, j, kc, kt
  real(dp), dimension(0:KCMAX) :: ctemp_c
  real(dp), dimension(0:KTMAX) :: ctemp_t

  do i=0, IMAX
  do j=0, JMAX

     if ( (mask(j,i) == 0).or.(mask(j,i) == 3) ) then
                                   ! glaciated land or floating ice

        if (n_cts(j,i) == 1) then
           do kt=0, KTMAX
              ctemp_t(kt) = (H_t(j,i)*dzeta_t) * temp_t_m(kt,j,i)
           end do
        else
           ctemp_t = 0.0_dp   ! not needed
        end if

        do kc=0, KCMAX
           ctemp_c(kc) = (H_c(j,i)*(aa*eaz_c(kc)/(ea-1.0_dp))*dzeta_c) &
                         * temp_c(kc,j,i)
        end do

        temp_mean(j,i) = 0.0_dp

        if (n_cts(j,i) == 1) then
           do kt=0, KTMAX-1
              temp_mean(j,i) = temp_mean(j,i)+0.5_dp*(ctemp_t(kt+1)+ctemp_t(kt))
           end do
        end if

        do kc=0, KCMAX-1
           temp_mean(j,i) = temp_mean(j,i)+0.5_dp*(ctemp_c(kc+1)+ctemp_c(kc))
        end do

        temp_mean(j,i) = temp_mean(j,i)/H(j,i)

     else   ! mask(j,i) == 1 or 2, ice-free land or sea

        ctemp_c = 0.0_dp   ! not needed
        ctemp_t = 0.0_dp   ! not needed

        temp_mean(j,i) = temp_c(0,j,i)

     end if

  end do
  end do

  end subroutine calc_temp_mean

!-------------------------------------------------------------------------------

end module calc_temp_aux_m
!
