!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ e n h a n c e _ m
!
!! Computation of the flow enhancement factor.
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
!> Computation of the flow enhancement factor.
!-------------------------------------------------------------------------------
module calc_enhance_m

  use sico_types_m
  use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

  use error_m

  implicit none

  private
  public :: calc_enhance_1, calc_enhance_2, calc_enhance_3
  public :: calc_enhance_4, calc_enhance_5
  public :: calc_enhance_hybrid_weighted

contains

!-------------------------------------------------------------------------------
!> Computation of the flow enhancement factor.
!! Case ENHMOD==1:
!! constant for grounded ice, constant for floating ice.
!-------------------------------------------------------------------------------
  subroutine calc_enhance_1()

  implicit none

  enh_t = ENH_FACT
  enh_c = ENH_FACT

  call calc_enhance_stream_const()   ! ice streams

#if (MARGIN==3)   /* floating ice */
  call calc_enhance_floating_const()
#endif

#if (defined(NMARS) || defined(SMARS))   /* Martian polar caps */
  call mod_enhance_dust()
#endif

  end subroutine calc_enhance_1

!-------------------------------------------------------------------------------
!> Computation of the flow enhancement factor.
!! Case ENHMOD==2:
!! two different values depending on age for grounded ice,
!! constant for floating ice.
!-------------------------------------------------------------------------------
  subroutine calc_enhance_2()

  implicit none

  integer(i4b) :: i, j, kc, kt
  real(dp)     :: age_trans

  age_trans = AGE_TRANS_0*year2sec

  do i=0, IMAX
  do j=0, JMAX

     do kt=0, KTMAX
        if (age_t(kt,j,i) < age_trans) then
           enh_t(kt,j,i) = ENH_INTG   ! Holocene ice
        else
           enh_t(kt,j,i) = ENH_FACT   ! Pleistocene ice
        end if
     end do

     do kc=0, KCMAX
        if (age_c(kc,j,i) < age_trans) then
           enh_c(kc,j,i) = ENH_INTG   ! Holocene ice
        else
           enh_c(kc,j,i) = ENH_FACT   ! Pleistocene ice
        end if
     end do

  end do
  end do

  call calc_enhance_stream_const()   ! ice streams

#if (MARGIN==3)   /* floating ice */
  call calc_enhance_floating_const()
#endif

#if (defined(NMARS) || defined(SMARS))   /* Martian polar caps */
  call mod_enhance_dust()
#endif

  end subroutine calc_enhance_2

!-------------------------------------------------------------------------------
!> Computation of the flow enhancement factor.
!! Case ENHMOD==3:
!! two different values depending on time of deposition for grounded ice,
!! constant for floating ice.
!-------------------------------------------------------------------------------
  subroutine calc_enhance_3(time)

  implicit none

  real(dp), intent(in) :: time

  integer(i4b) :: i, j, kc, kt
  real(dp)     :: date_trans1, date_trans2, date_trans3

  date_trans1 = DATE_TRANS1_0*year2sec
  date_trans2 = DATE_TRANS2_0*year2sec
  date_trans3 = DATE_TRANS3_0*year2sec

  do i=0, IMAX
  do j=0, JMAX

     do kt=0, KTMAX
        if ( (time-age_t(kt,j,i)) < date_trans1 ) then
           enh_t(kt,j,i) = ENH_FACT   ! pre-Eemian ice
        else if ( ((time-age_t(kt,j,i)) >= date_trans1).and. &
                  ((time-age_t(kt,j,i)) <  date_trans2) ) then
           enh_t(kt,j,i) = ENH_INTG   ! Eemian ice
        else if ( ((time-age_t(kt,j,i)) >= date_trans2).and. &
                  ((time-age_t(kt,j,i)) <  date_trans3) ) then
           enh_t(kt,j,i) = ENH_FACT   ! Weichselian ice
        else
           enh_t(kt,j,i) = ENH_INTG   ! Holocene ice
        end if
     end do

     do kc=0, KCMAX
        if ( (time-age_c(kc,j,i)) < date_trans1 ) then
           enh_c(kc,j,i) = ENH_FACT   ! pre-Eemian ice
        else if ( ((time-age_c(kc,j,i)) >= date_trans1).and. &
                  ((time-age_c(kc,j,i)) <  date_trans2) ) then
           enh_c(kc,j,i) = ENH_INTG   ! Eemian ice
        else if ( ((time-age_c(kc,j,i)) >= date_trans2).and. &
                  ((time-age_c(kc,j,i)) <  date_trans3) ) then
           enh_c(kc,j,i) = ENH_FACT   ! Weichselian ice
        else
           enh_c(kc,j,i) = ENH_INTG   ! Holocene ice
        end if
     end do

  end do
  end do

  call calc_enhance_stream_const()   ! ice streams

#if (MARGIN==3)   /* floating ice */
  call calc_enhance_floating_const()
#endif

#if (defined(NMARS) || defined(SMARS))   /* Martian polar caps */
  call mod_enhance_dust()
#endif

  end subroutine calc_enhance_3

!-------------------------------------------------------------------------------
!> Computation of the flow enhancement factor.
!! Case ENHMOD==4:
!! minimal anisotropic enhancement factor for grounded ice,
!! constant for floating ice.
!-------------------------------------------------------------------------------
  subroutine calc_enhance_4()

  implicit none

  call calc_enhance_aniso()

#if (MARGIN==3)   /* floating ice */
  call calc_enhance_floating_const()
#endif

#if (defined(NMARS) || defined(SMARS))   /* Martian polar caps */
  call mod_enhance_dust()
#endif

  end subroutine calc_enhance_4

!-------------------------------------------------------------------------------
!> Computation of the flow enhancement factor.
!! Case ENHMOD==5:
!! minimal anisotropic enhancement factor for grounded and floating ice.
!-------------------------------------------------------------------------------
  subroutine calc_enhance_5()

  implicit none

  call calc_enhance_aniso()

#if (defined(NMARS) || defined(SMARS))   /* Martian polar caps */
  call mod_enhance_dust()
#endif

  end subroutine calc_enhance_5

!-------------------------------------------------------------------------------
!> Weighted enhancement factor for SIA/SStA hybrid dynamics.
!-------------------------------------------------------------------------------
  subroutine calc_enhance_hybrid_weighted(weigh_ssta_sia)

  implicit none

  real(dp), dimension(0:JMAX,0:IMAX), intent(in) :: weigh_ssta_sia

  integer(i4b) :: i, j, kc, kt

  if (flag_enh_stream) then

     do i=0, IMAX
     do j=0, JMAX

        if (flag_shelfy_stream(j,i)) then   ! shelfy stream

           do kt=0, KTMAX
              enh_t(kt,j,i) = weigh_ssta_sia(j,i)*enh_stream &
                              + (1.0_dp-weigh_ssta_sia(j,i))*enh_t(kt,j,i)
           end do

           do kc=0, KCMAX
              enh_c(kc,j,i) = weigh_ssta_sia(j,i)*enh_stream &
                              + (1.0_dp-weigh_ssta_sia(j,i))*enh_c(kc,j,i)
           end do

        end if

     end do
     end do

  end if

  end subroutine calc_enhance_hybrid_weighted

!-------------------------------------------------------------------------------
!> Minimal anisotropic flow enhancement factor.
!-------------------------------------------------------------------------------
  subroutine calc_enhance_aniso()

  implicit none

  integer(i4b) :: i, j, kc, kt
  real(dp)     :: enh_shear, enh_compr
  real(dp)     :: enh_shear_compr_diff

#if (defined(ENH_SHEAR))
  enh_shear = ENH_SHEAR
#else
  enh_shear = 1.0_dp
#endif

#if (defined(ENH_COMPR))
  enh_compr = ENH_COMPR
#else
  enh_compr = 1.0_dp
#endif

  enh_shear_compr_diff = enh_shear-enh_compr

  do i=0, IMAX
  do j=0, JMAX

     do kt=0, KTMAX
        enh_t(kt,j,i) = enh_compr &
                        + enh_shear_compr_diff*lambda_shear_t(kt,j,i)**2
     end do

     do kc=0, KCMAX
        enh_c(kc,j,i) = enh_compr &
                        + enh_shear_compr_diff*lambda_shear_c(kc,j,i)**2
     end do

  end do
  end do

  end subroutine calc_enhance_aniso

!-------------------------------------------------------------------------------
!> Constant, prescribed flow enhancement factor for ice streams
!! (fast-flowing ice).
!-------------------------------------------------------------------------------
  subroutine calc_enhance_stream_const()

  implicit none

  flag_enh_stream = .false.

#if ((DYNAMICS==2 || DYNAMICS==3) && defined(ENH_STREAM))
  enh_stream = ENH_STREAM
  if (enh_stream >= 0.0_dp) flag_enh_stream = .true.
#else
  enh_stream = no_value_neg_1   ! negative dummy value
#endif

  end subroutine calc_enhance_stream_const

!-------------------------------------------------------------------------------
!> Constant, prescribed flow enhancement factor for floating ice.
!-------------------------------------------------------------------------------
  subroutine calc_enhance_floating_const()

  implicit none

  integer(i4b) :: i, j, kc, kt
  real(dp)     :: enh_shelf

#if (defined(ENH_SHELF))
  enh_shelf = ENH_SHELF
#else
  enh_shelf = 1.0_dp
#endif

  do i=0, IMAX
  do j=0, JMAX

     if ( mask(j,i)==3 ) then   ! floating ice

        do kt=0, KTMAX
           enh_t(kt,j,i) = enh_shelf
        end do

        do kc=0, KCMAX
           enh_c(kc,j,i) = enh_shelf
        end do

     end if

  end do
  end do

  end subroutine calc_enhance_floating_const

!-------------------------------------------------------------------------------
!> Modification of the flow enhancement factor due to dust content
!! (for the Martian polar caps).
!-------------------------------------------------------------------------------
  subroutine mod_enhance_dust()

  implicit none

  real(dp) :: frac_dust
  real(dp) :: enh_mult

#if (defined(FRAC_DUST))
  frac_dust = FRAC_DUST
#else
  frac_dust = 0.0_dp
#endif

#if (FLOW_LAW==1)
  enh_mult = exp(-3.0_dp*2.0_dp*frac_dust)
#elif (FLOW_LAW==2)
  enh_mult = exp(-1.8_dp*2.0_dp*frac_dust)
#elif (FLOW_LAW==3)
  enh_mult = exp(-4.0_dp*2.0_dp*frac_dust)
#elif (FLOW_LAW==4)
  errormsg = ' >>> mod_enhance_dust: FLOW_LAW==4 has not been implemented yet!'
  call error(errormsg)
#endif

  enh_t = enh_t * enh_mult
  enh_c = enh_c * enh_mult

  end subroutine mod_enhance_dust

!-------------------------------------------------------------------------------

end module calc_enhance_m
!
