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

#if ((MARGIN==3) && (ICE_SHELF_CALVING==5))
  subroutine frontal_calving(dtime)

#if !(defined(FRONTAL_CALVING_RATE))
  errormsg = ' >>> frontal_calving: FRONTAL_CALVING_RATE undefined! Required for ICE_SHELF_CALVING==5'
  call error(errormsg)
#endif 

  implicit none

  real(dp), intent(in) :: dtime
  integer(i4b)         :: i, j, ij
  real(dp)             :: dHdt_retreat, DX_inMeter_inv, F_rate, H_new_tmp, dtime_inv, calv_retreat_mask

  dtime_inv = 1.0_dp/dtime
  DX_inMeter_inv = 1.0_dp/(DX*1000.0_dp)
  F_rate = FRONTAL_CALVING_RATE*DX_inMeter_inv/year2sec
  
  do ij=1, (IMAX+1)*(JMAX+1)
    ! add inner point flag here

    i = n2i(ij)   ! i=0...IMAX
    j = n2j(ij)   ! j=0...JMAX

! #if (ICE_SHELF_CALVING_TYPE==0)
!     if ( ( ((mask(j,i)==0).and.(zl(j,i)<z_sl(j,i))).or.(mask(j,i)==3) ) &   ! grounded or floating ice
! #elif (ICE_SHELF_CALVING_TYPE==1)
!     if ( (mask(j,i)==3) &   ! floating ice
! #elif (ICE_SHELF_CALVING_TYPE==2)
!     if ( ((mask(j,i)==0).and.(zl(j,i)<z_sl(j,i))) &   ! grounded ice
! #else
!     errormsg = 'calc_thk_mask_update_aux3: ICE_SHELF_CALVING_TYPE must be O or 1 or 2.'
!     call error(errormsg)
!     if (.true. &
! #endif
    if (.true. &
      .and. &
         (    (mask(j,i+1)==2)   &   ! with
          .or.(mask(j,i-1)==2)   &   ! one
          .or.(mask(j+1,i)==2)   &   ! neighbouring
          .or.(mask(j-1,i)==2) ) &   ! sea point
    ) then
      ! dHdt_retreat = H(j,i)*F_rate
      ! calving(j,i) = calving(j,i) + dHdt_retreat
      
      H_new_tmp = H_new(j,i)
      dHdt_retreat = H(j,i)*F_rate
      H_new(j,i) = max((H_new(j,i) - dHdt_retreat*dtime), 0.0_dp)
      calv_retreat_mask = (H_new_tmp-H_new(j,i))*dtime_inv
      calving(j,i) = calving(j,i) + calv_retreat_mask
    end if
  end do

  end subroutine frontal_calving


#endif /* ((MARGIN==3) && (ICE_SHELF_CALVING==5)) */


end module calving_m

