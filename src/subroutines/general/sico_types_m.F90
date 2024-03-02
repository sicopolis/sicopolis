!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ t y p e s _ m
!
!! Type declarations for SICOPOLIS.
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
!> Type declarations for SICOPOLIS.
!-------------------------------------------------------------------------------
module sico_types_m

implicit none
save

#if !defined(ALLOW_TAPENADE)

integer, parameter :: i4b = selected_int_kind(9)
   !! 4-byte integers
integer, parameter :: sp  = kind(1.0)
   !! Single-precision reals
integer, parameter :: dp  = kind(1.0d0)
   !! Double-precision reals

#else

integer, parameter :: i4b = 4
   !! 4-byte integers
integer, parameter :: sp  = 4
   !! Single-precision reals
integer, parameter :: dp  = 8
   !! Double-precision reals

#endif

type flag_firstcall
   !! First-call flags
   logical :: boundary = .true.
      !! First-call flag for the routine 'boundary'
   logical :: calc_thk_water_bas = .true.
      !! First-call flag for the routine 'calc_thk_water_bas'
   logical :: output1 = .true.
      !! First-call flag for the routine 'output1'
   logical :: output2 = .true.
      !! First-call flag for the routine 'output2'
   logical :: output4 = .true.
      !! First-call flag for the routine 'output4'
   logical :: sub_ice_shelf_melting_param_2 = .true.
      !! First-call flag for the routine 'sub_ice_shelf_melting_param_2'
end type flag_firstcall

end module sico_types_m
!
