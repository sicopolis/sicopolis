!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l v i n g _ m
!
!! Calving of grounded or floating ice.
!!
!!##### Authors
!!
!! Ralf Greve, Thorben Dunse
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
!> Calving of grounded or floating ice.
!-------------------------------------------------------------------------------
module calving_m

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
!> Calving of grounded "underwater ice".
!-------------------------------------------------------------------------------
  subroutine calving_underwater_ice()

  implicit none

  real(dp)                           :: rhosw_rho_ratio
  real(dp)                           :: calv_uw_coeff, r1_calv_uw, r2_calv_uw
  real(dp)                           :: H0_flt
  real(dp), dimension(0:JMAX,0:IMAX) :: H_sea, calv_uw_ice
  integer(i4b)                       :: i, j

!-------- Term abbreviations --------

  rhosw_rho_ratio = RHO_SW/RHO

!-------- Setting of parameters --------

#if (defined(CALV_UW_COEFF))
  calv_uw_coeff = CALV_UW_COEFF *sec2year
#else
  errormsg = ' >>> calving_underwater_ice: CALV_UW_COEFF undefined!'
  call error(errormsg)
#endif

#if (defined(R1_CALV_UW))
  r1_calv_uw = R1_CALV_UW
#else
  errormsg = ' >>> calving_underwater_ice: R1_CALV_UW undefined!'
  call error(errormsg)
#endif

#if (defined(R2_CALV_UW))
  r2_calv_uw = R2_CALV_UW
#else
  errormsg = ' >>> calving_underwater_ice: R2_CALV_UW undefined!'
  call error(errormsg)
#endif

#if (defined(H0_FLOAT))
  H0_flt = H0_FLOAT
#else
  H0_flt = 0.0_dp
#endif

!-------- Sea depth --------

  H_sea = max(z_sl - zl, 0.0_dp)   ! sea depth

!-------- Calving of "underwater ice" --------

  do i=0, IMAX
  do j=0, JMAX

     if ( (mask(j,i) == 0) &
          .and. (H(j,i) < rhosw_rho_ratio*H_sea(j,i)+H0_flt) ) then
        calv_uw_ice(j,i) = calv_uw_coeff &
                           * H(j,i)**r1_calv_uw * H_sea(j,i)**r2_calv_uw
     else
        calv_uw_ice(j,i) = 0.0_dp
     end if

  end do
  end do

  calving = calving + calv_uw_ice

  end subroutine calving_underwater_ice

#if (RETREAT_MASK==1 || ICE_SHELF_COLLAPSE_MASK==1)
!-------------------------------------------------------------------------------
!> Adjustment of the newly computed ice thickness distribution due to either
!! the retreat mask due to oceanic forcing or the ice-shelf collapse mask
!! (counted as calving).
!-------------------------------------------------------------------------------
  subroutine calving_retreat_mask(time, dtime, i, j)

  implicit none

  integer(i4b), intent(in) :: i, j
  real(dp)    , intent(in) :: time, dtime

  real(dp) :: H_new_tmp, dHdt_retreat
  real(dp) :: calv_retreat_mask
  real(dp) :: dtime_inv
  real(dp) :: dtime_1year, dtime_1year_inv

  dtime_inv       = 1.0_dp/dtime
  dtime_1year     = year2sec   ! 1 year (in seconds)
  dtime_1year_inv = 1.0_dp/dtime_1year

!-------- Saving computed H_new before any adjustments --------

  H_new_tmp = H_new(j,i)

!-------- Adjustment due to the retreat mask --------

  dHdt_retreat = 0.0_dp   ! initialization

#if (RETREAT_MASK==1)
  if (H_new(j,i) > 0.0_dp) then
#elif (ICE_SHELF_COLLAPSE_MASK==1)
  if ((H_new(j,i) > 0.0_dp).and.(mask(j,i)==3)) then
#endif

    dHdt_retreat = -(1.0_dp-r_mask_retreat(j,i))*H_ref_retreat(j,i) &
                                                 *dtime_1year_inv

    H_new(j,i) = max((H_new(j,i) + dHdt_retreat*dtime), 0.0_dp)

  end if

!-------- Computation of the mass balance adjustment --------

  calv_retreat_mask = (H_new_tmp-H_new(j,i))*dtime_inv
                      ! calving is counted as positive for mass loss

  calving(j,i) = calving(j,i) + calv_retreat_mask

  end subroutine calving_retreat_mask

#endif   /* (RETREAT_MASK==1 || ICE_SHELF_COLLAPSE_MASK==1) */

!-------------------------------------------------------------------------------

#if (ICE_SHELF_CALVING==6)
  subroutine frontal_calving(i, j, dtime)

#if !(defined(FRONTAL_CALVING_RATE))
  errormsg = ' >>> xxx: FRONTAL_CALVING_RATE undefined!'
  call error(errormsg)
#endif  

  implicit none

  integer(i4b), intent(in) :: i, j
  real(dp)    , intent(in) :: dtime, H_lost_rate

  if ((mask(j,i)==3).and.( (mask(j, i-1)==2).or.(mask(j, i+1)==2).or.(mask(j-1, i)==2).or.(mask(j+1, i)==2) ))  then
    H_lost_rate = H_new(j,i)*FRONTAL_CALVING_RATE/DX
  else
    H_lost_rate = 0
  endif

  ! H_new(j,i) = max(H_new(j,i) - H_lost_rate*dtime, 0.0_dp)

  ! calving(j,i) = calving(j,i) + H_lost_rate
  
end subroutine frontal_calving

#endif /* (ICE_SHELF_CALVING==6) */
end module calving_m
!
