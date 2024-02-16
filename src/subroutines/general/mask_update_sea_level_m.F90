!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  m a s k _ u p d a t e _ s e a _ l e v e l _ m
!
!! Update of the ice-land-ocean mask due to changes of the sea level.
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
!> Update of the ice-land-ocean mask due to changes of the sea level.
!-------------------------------------------------------------------------------
module mask_update_sea_level_m

  use sico_types_m
  use sico_variables_m
  use sico_vars_m

  implicit none

  private
  public :: mask_update_sea_level

contains

!-------------------------------------------------------------------------------
!> Main function of mask_update_m:
!! Update of the ice-land-ocean mask due to changes of the sea level.
!-------------------------------------------------------------------------------
  function mask_update_sea_level(i, j)

  implicit none

  integer(i4b), intent(in) :: i, j

  integer(i4b) :: mask_update_sea_level
  real(dp)     :: rhosw_rho_ratio, H_sea

  rhosw_rho_ratio = RHO_SW/RHO

!-------- Previously ice-free land --------

  if (mask(j,i) == 1) then

     if (zl(j,i) >= z_sl(j,i)) then
        mask_update_sea_level = 1   ! now ice-free land (as before)
        return
     else
        mask_update_sea_level = 2   ! now sea point
        return
     end if

!-------- Previously sea point --------

  else if (mask(j,i) == 2) then

     if (zl(j,i) <= z_sl(j,i)) then
        mask_update_sea_level = 2   ! now sea point (as before)
        return
     else
        mask_update_sea_level = 1   ! now ice-free land
        return
     end if

!-------- Previously grounded-ice or floating-ice point --------

  else   ! (mask(j,i) == 0, 3)

     if (zl(j,i) > z_sl(j,i)) then

        mask_update_sea_level = 0   ! now grounded ice
        return

     else

        H_sea = z_sl(j,i)-zl(j,i)   ! sea depth

        if ( H(j,i) < (rhosw_rho_ratio*H_sea) ) then

#if (MARGIN==1 || (MARGIN==2 && MARINE_ICE_FORMATION==1))
           mask_update_sea_level = 2     ! ice becomes floating, therefore
                                         ! now sea point (ice cut off)
#elif (MARGIN==2 && MARINE_ICE_FORMATION==2)
           mask_update_sea_level = 0     ! now "underwater ice"
#elif (MARGIN==3)
           mask_update_sea_level = 3     ! now floating ice
#endif
           return

        else

           mask_update_sea_level = 0     ! now grounded ice
           return

        end if

     end if

  end if

  end function mask_update_sea_level

!-------------------------------------------------------------------------------

end module mask_update_sea_level_m
!
