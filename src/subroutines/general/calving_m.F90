!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l v i n g _ m
!
!! Calving of grounded or floating ice.
!!
!!##### Authors
!!
!! Ralf Greve, Thorben Dunse, Marius Schaefer
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
!> Subroutine which applies a constant specific mass loss at grid cells
!! which have ocean in at one of its neighbouring grid cells.
!-------------------------------------------------------------------------------
  subroutine constant_calving()

  implicit none

  real(dp)                           :: year_sec_inv
  real(dp)                           :: rhosw_rho_ratio
  real(dp)                           :: calv_const
  real(dp), dimension(0:JMAX,0:IMAX) :: H_sea, calv_ice
  integer(i4b)                       :: i, j
  integer(i4b)                       :: im1, ip1, jm1, jp1, n_ocn

!-------- Term abbreviations --------

  year_sec_inv = 1.0_dp/year2sec

  rhosw_rho_ratio = RHO_SW/RHO

!-------- Setting of parameters --------

#if (defined(KONST_CALV))
  calv_const = KONST_CALV * year_sec_inv
#else
  errormsg = ' >>> calving_constant_ice: KONST_CALV undefined!'
  call error(errormsg)
#endif

!-------- Sea depth --------

  H_sea = max(z_sl - zl, 0.0_dp)   ! sea depth

! ------- At least one ocean neighbours to allow calving? -------

  do i=0, IMAX
  do j=0, JMAX
    ! Define neighbor indices
    im1 = max(i-1,0)
    ip1 = min(i+1,IMAX)
    jm1 = max(j-1,0)
    jp1 = min(j+1,JMAX)
    
    n_ocn=count([mask(j,im1),mask(j,ip1),mask(jm1,i),mask(jp1,i)].gt. 1 )

!-------- Constant calving --------

    if ( (mask(j,i) == 0) .and. (n_ocn .gt. 0) ) then
       !write(*,*)  'CALVING!'
       calv_ice(j,i) = calv_const
    else
       calv_ice(j,i) = 0.0_dp
    end if

  end do
  end do

  calving = calving + calv_ice

  end subroutine constant_calving
  
!-------------------------------------------------------------------------------
!> Subroutine which applies a velocity dependent specific mass loss
!! at grid cells which have ocean at one of its neighbouring grid cells.
!-------------------------------------------------------------------------------
  subroutine velocity_calving()

  implicit none

  real(dp)                           :: year_sec_inv
  real(dp)                           :: scale_factor
  real(dp)                           :: cellwidth
  real(dp), dimension(0:JMAX,0:IMAX) :: calv_ice, vmean
  integer(i4b)                       :: i, j
  integer(i4b)                       :: im1, ip1, jm1, jp1, n_ocn

!-------- Term abbreviations --------

  year_sec_inv = 1.0_dp/year2sec

!-------- Setting of parameters --------

  cellwidth=DX
#if (defined(VELOCITY_CALV))
  scale_factor = VELOCITY_CALV * year_sec_inv
#else
  errormsg = ' >>> calving_constant_ice: VELOCITY_CALV undefined!'
  call error(errormsg)
#endif

! ------- At least one ocean neighbours to allow calving? -------

  do i=0, IMAX
  do j=0, JMAX
    ! Define neighbor indices
    im1 = max(i-1,0)
    ip1 = min(i+1,IMAX)
    jm1 = max(j-1,0)
    jp1 = min(j+1,JMAX)
    
    n_ocn=count([mask(j,im1),mask(j,ip1),mask(jm1,i),mask(jp1,i)].gt. 1 )

!-------- Constant calving --------

    if ( (mask(j,i) == 0) .and. (n_ocn .gt. 0) ) then
       !write(*,*)  'CALVING!'
       vmean(j,i)=sqrt(vx_m(j,i)**2+vy_m(j,i)**2)
       calv_ice(j,i) = scale_factor*H(j,i)*vmean(j,i)/(cellwidth*1000)
    else
       calv_ice(j,i) = 0.0_dp
    end if

  end do
  end do

  calving = calving + calv_ice

  end subroutine velocity_calving

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
  subroutine calving_retreat_mask(time, dtime)

  implicit none

  real(dp), intent(in) :: time, dtime

  integer(i4b) :: i, j
  real(dp), dimension(0:JMAX,0:IMAX) :: H_new_tmp, dHdt_retreat
  real(dp), dimension(0:JMAX,0:IMAX) :: calv_retreat_mask
  real(dp)                           :: dtime_inv
  real(dp)                           :: dtime_1year, dtime_1year_inv

  dtime_inv       = 1.0_dp/dtime
  dtime_1year     = year2sec   ! 1 year (in seconds)
  dtime_1year_inv = 1.0_dp/dtime_1year

!-------- Saving computed H_new before any adjustments --------

  H_new_tmp = H_new

!-------- Adjustment due to the retreat mask --------

  dHdt_retreat = 0.0_dp   ! initialization

  do i=0, IMAX
  do j=0, JMAX

#if (RETREAT_MASK==1)
     if (H_new(j,i) > 0.0_dp) then
#elif (ICE_SHELF_COLLAPSE_MASK==1)
     if ((H_new(j,i) > 0.0_dp).and.(mask(j,i)==3)) then
#endif

        dHdt_retreat(j,i) = -(1.0_dp-r_mask_retreat(j,i))*H_ref_retreat(j,i) &
                                                         *dtime_1year_inv

        H_new(j,i) = max((H_new(j,i) + dHdt_retreat(j,i)*dtime), 0.0_dp)

     end if

  end do
  end do

!-------- Computation of the mass balance adjustment --------

  calv_retreat_mask = (H_new_tmp-H_new)*dtime_inv
                      ! calving is counted as positive for mass loss

  calving = calving + calv_retreat_mask

  end subroutine calving_retreat_mask

#endif   /* (RETREAT_MASK==1 || ICE_SHELF_COLLAPSE_MASK==1) */

!-------------------------------------------------------------------------------

end module calving_m
!
