!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  f l a g _ u p d a t e _ g f _ g l _ c f _ m
!
!! Update of the flags for the land-terminating grounded front,
!! marine-terminating grounded front, grounding line and calving front.
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
!> Update of the flags for the land-terminating grounded front,
!! marine-terminating grounded front, grounding line and calving front.
!-------------------------------------------------------------------------------
module flag_update_gf_gl_cf_m

  use sico_types_m
  use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

  implicit none

  private
  public :: flag_update_gf_gl_cf

contains

!-------------------------------------------------------------------------------
!> Main subroutine of flag_update_gf_gl_cf_m:
!! Update of the flags for the land-terminating grounded front,
!! marine-terminating grounded front, grounding line and calving front.
!-------------------------------------------------------------------------------
  subroutine flag_update_gf_gl_cf()

  implicit none

  integer(i4b) :: i, j

  flag_grounded_front_a_1 = .false.
  flag_grounded_front_a_2 = .false.

  flag_grounded_front_b_1 = .false.
  flag_grounded_front_b_2 = .false.

  flag_grounding_line_1 = .false.
  flag_grounding_line_2 = .false.

  flag_calving_front_1 = .false.
  flag_calving_front_2 = .false.

  do i=1, IMAX-1
  do j=1, JMAX-1

     if ( (mask(j,i)==0) &   ! grounded ice point
          .and. &
            (    (mask(j,i+1)==1)   &   ! with
             .or.(mask(j,i-1)==1)   &   ! one
             .or.(mask(j+1,i)==1)   &   ! neighbouring
             .or.(mask(j-1,i)==1) ) &   ! ice-free land point
        ) &
     flag_grounded_front_a_1(j,i) = .true.

     if ( (mask(j,i)==1) &   ! ice-free land point
          .and. &
            (    (mask(j,i+1)==0)   &   ! with
             .or.(mask(j,i-1)==0)   &   ! one
             .or.(mask(j+1,i)==0)   &   ! neighbouring
             .or.(mask(j-1,i)==0) ) &   ! grounded ice point
        ) &
     flag_grounded_front_a_2(j,i) = .true.

     if ( (mask(j,i)==0) &   ! grounded ice point
          .and. &
            (    (mask(j,i+1)==2)   &   ! with
             .or.(mask(j,i-1)==2)   &   ! one
             .or.(mask(j+1,i)==2)   &   ! neighbouring
             .or.(mask(j-1,i)==2) ) &   ! ocean point
        ) &
     flag_grounded_front_b_1(j,i) = .true.

     if ( (mask(j,i)==2) &   ! ocean point
          .and. &
            (    (mask(j,i+1)==0)   &   ! with
             .or.(mask(j,i-1)==0)   &   ! one
             .or.(mask(j+1,i)==0)   &   ! neighbouring
             .or.(mask(j-1,i)==0) ) &   ! grounded ice point
        ) &
     flag_grounded_front_b_2(j,i) = .true.

#if (MARGIN==3)

     if ( (mask(j,i)==0) &   ! grounded ice point
          .and. &
            (    (mask(j,i+1)==3)   &   ! with
             .or.(mask(j,i-1)==3)   &   ! one
             .or.(mask(j+1,i)==3)   &   ! neighbouring
             .or.(mask(j-1,i)==3) ) &   ! floating ice point
        ) &
     flag_grounding_line_1(j,i) = .true.

     if ( (mask(j,i)==3) &   ! floating ice point
          .and. &
            (    (mask(j,i+1)==0)   &   ! with
             .or.(mask(j,i-1)==0)   &   ! one
             .or.(mask(j+1,i)==0)   &   ! neighbouring
             .or.(mask(j-1,i)==0) ) &   ! grounded ice point
        ) &
     flag_grounding_line_2(j,i) = .true.

     if ( (mask(j,i)==3) &   ! floating ice point
          .and. &
            (    (mask(j,i+1)==2)   &   ! with
             .or.(mask(j,i-1)==2)   &   ! one
             .or.(mask(j+1,i)==2)   &   ! neighbouring
             .or.(mask(j-1,i)==2) ) &   ! ocean point
        ) &
     flag_calving_front_1(j,i) = .true.

     if ( (mask(j,i)==2) &   ! ocean point
          .and. &
            (    (mask(j,i+1)==3)   &   ! with
             .or.(mask(j,i-1)==3)   &   ! one
             .or.(mask(j+1,i)==3)   &   ! neighbouring
             .or.(mask(j-1,i)==3) ) &   ! floating ice point
        ) &
     flag_calving_front_2(j,i) = .true.

#endif   /* MARGIN==3 */

  end do
  end do

  do i=1, IMAX-1

     j=0

     if ( (mask(j,i)==1) &   ! ice-free land point
          .and. (mask(j+1,i)==0) &   ! with one neighbouring
        ) &                               ! grounded ice point
     flag_grounded_front_a_2(j,i) = .true.

     if ( (mask(j,i)==2) &   ! ocean point
          .and. (mask(j+1,i)==0) &   ! with one neighbouring
        ) &                               ! grounded ice point
     flag_grounded_front_b_2(j,i) = .true.

#if (MARGIN==3)

     if ( (mask(j,i)==2) &   ! ocean point
          .and. (mask(j+1,i)==3) &   ! with one neighbouring
        ) &                               ! floating ice point
     flag_calving_front_2(j,i) = .true.

#endif   /* MARGIN==3 */

     j=JMAX

     if ( (mask(j,i)==1) &   ! ice-free land point
          .and. (mask(j-1,i)==0) &   ! with one neighbouring
        ) &                               ! grounded ice point
     flag_grounded_front_a_2(j,i) = .true.

     if ( (mask(j,i)==2) &   ! ocean point
          .and. (mask(j-1,i)==0) &   ! with one neighbouring
        ) &                               ! grounded ice point
     flag_grounded_front_b_2(j,i) = .true.

#if (MARGIN==3)

     if ( (mask(j,i)==2) &   ! ocean point
          .and. (mask(j-1,i)==3) &   ! with one neighbouring
        ) &                               ! floating ice point
     flag_calving_front_2(j,i) = .true.

#endif   /* MARGIN==3 */

  end do

  do j=1, JMAX-1

     i=0

     if ( (mask(j,i)==1) &   ! ice-free land point
          .and. (mask(j,i+1)==0) &   ! with one neighbouring
        ) &                               ! grounded ice point
     flag_grounded_front_a_2(j,i) = .true.

     if ( (mask(j,i)==2) &   ! ocean point
          .and. (mask(j,i+1)==0) &   ! with one neighbouring
        ) &                               ! grounded ice point
     flag_grounded_front_b_2(j,i) = .true.

#if (MARGIN==3)

     if ( (mask(j,i)==2) &   ! ocean point
          .and. (mask(j,i+1)==3) &   ! with one neighbouring
        ) &                               ! floating ice point
     flag_calving_front_2(j,i) = .true.

#endif   /* MARGIN==3 */

     i=IMAX

     if ( (mask(j,i)==1) &   ! ice-free land point
          .and. (mask(j,i-1)==0) &   ! with one neighbouring
        ) &                               ! grounded ice point
     flag_grounded_front_a_2(j,i) = .true.

     if ( (mask(j,i)==2) &   ! ocean point
          .and. (mask(j,i-1)==0) &   ! with one neighbouring
        ) &                               ! grounded ice point
     flag_grounded_front_b_2(j,i) = .true.

#if (MARGIN==3)

     if ( (mask(j,i)==2) &   ! ocean point
          .and. (mask(j,i-1)==3) &   ! with one neighbouring
        ) &                               ! floating ice point
     flag_calving_front_2(j,i) = .true.

#endif   /* MARGIN==3 */

  end do

  end subroutine flag_update_gf_gl_cf

!-------------------------------------------------------------------------------

end module flag_update_gf_gl_cf_m
!
