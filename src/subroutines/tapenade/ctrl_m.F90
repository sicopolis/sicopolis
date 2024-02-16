!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module  :  c t r l _ m
!
!! Declarations of control variables for adjointing.
!!
!!##### Authors
!!
!! Liz Curry-Logan, Shreyas Sunil Gaikwad,
!! Sri Hari Krishna Narayanan
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
!> Declarations of control variables for adjointing.
!-------------------------------------------------------------------------------
module ctrl_m

  use sico_types_m  
  use sico_variables_m
  use sico_vars_m

  implicit none

  public :: cost_final 
#if (defined(ALLOW_GRDCHK) || defined(ALLOW_TAPENADE))
  public :: myceiling, myfloor 
#endif

contains
 
!-------------------------------------------------------------------------------
!> This is the final cost calculation.
!! The cost function structure is defined here. 
!!
!! Currently there are three options - 
!!
!! 1. AGE_COST
!! Currently is a "observed age" - modeled age summed over the entire 
!! domain. The "observed age" is a fake, generated age
!! field performed by the 125 ka run in headers.
!!
!! 2. BEDMACHINE_COST
!! Based on L2 misfit with dataset in -
!! BedMachine v3: Complete Bed Topography and Ocean Bathymetry Mapping of
!! Greenland From Multibeam Echo Sounding Combined With Mass Conservation
!! by Morlighem et. al in 2017.
!!
!! 3. Total volume - Default (Only ALLOW_COST used)
!! This just defines the cost function to be the total volume of the ice sheet.
!! This is currently the default.
!!
!! Other cost functions are certainly possible, and recommended! 
!-------------------------------------------------------------------------------
  subroutine cost_final()
  
  implicit none
  
  integer(i4b) :: i, j, k, kc, kt, ios, KDATA
  
  !-------- Calculate the difference between the modeled and 'observed' ages:
fc = 0.0

#ifdef AGE_COST

#if (CALCMOD!=1)
  KDATA = KCMAX
#else 
  KDATA = KCMAX + KTMAX
print *, '>>> error: CALCMOD == 1 but final cost not properly working for '
print *, '           AGE_COST simulations'
#endif
 
  do k=0, KDATA 
    do j=0, JMAX
      do i=0, IMAX

        ! only counting points that are real in the data: 
        if (  age_data(k,j,i) .ge. -0.5) then

          fc = fc + 1.e-20*((age_data(k,j,i) - age_c(k,j,i)))**2

        end if

      end do
    end do
  end do

  do i=0, IMAX
    do j=0, JMAX
      fc = fc + 1.e-20*(H(j,i) - H_data(j,i))**2
    end do
  end do

#elif defined(BEDMACHINE_COST)
    do i=0, IMAX
      do j=0, JMAX
        fc = fc &
        + (H(j,i) - H_BedMachine_data(j,i))**2/H_unc_BedMachine_data(j,i)**2
      end do
    end do

#else
    do i=0, IMAX
      do j=0, JMAX
        !--- Other cost functions:
        fc = fc + (H_c(j,i) + H_t(j,i))*cell_area(j,i)
      end do
    end do

#endif
  
  !-------- Print to screen just in case something gets
  !         crazy with the file outputting:
  print *, 'Final cost, fc = ', fc
  print *, trim(OUT_PATH)
  end subroutine cost_final

#if (defined(ALLOW_GRDCHK) || defined(ALLOW_TAPENADE))
  subroutine myceiling(num, output)
  
  implicit none
  
  integer(i4b)             :: inum
  real(dp), intent(in) :: num
  integer(i4b), intent(out) :: output
  
  inum = int(num);
  
  inum = int(num);
  if (abs(num-real(inum,dp)) <= &
    (abs(num+real(inum,dp))*myepsilon_dp) ) then
    output = inum;
  else if (num>0) then
    output = inum + 1;
  else if (num<0) then
    output = inum;
  end if
  end subroutine myceiling
  
  subroutine myfloor(num, output)
  
  implicit none
  
  integer(i4b)             :: inum
  real(dp), intent(in) :: num
  integer(i4b), intent(out) :: output
  
  inum = int(num);
  
  if (abs(num-real(inum,dp)) <= &
    (abs(num+real(inum,dp))*myepsilon_dp) ) then
    output = inum;
  else if (num>0) then
    output = inum;
  else if (num<0) then
    output = inum - 1;
  end if
  
  end subroutine myfloor
#endif

end module ctrl_m
!
