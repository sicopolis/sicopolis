!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ t y p e s _ m
!
!> @file
!!
!! Declarations of kind types for SICOPOLIS.
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
!> Declarations of kind types for SICOPOLIS.
!<------------------------------------------------------------------------------
module sico_types_m

implicit none
save

#if !defined(ALLOW_TAPENADE)

! integer, parameter :: i1b = selected_int_kind(2)   !< 1-byte integers
integer, parameter :: i4b = selected_int_kind(9)   !< 4-byte integers
integer, parameter :: sp  = kind(1.0)              !< Single-precision reals
integer, parameter :: dp  = kind(1.0d0)            !< Double-precision reals

#else

! integer, parameter :: i1b = 4
!    ! SHK: found that setting i1b=4 (4-byte integers) is the only way
!    !      to make OpenAD not have a whirl opcode error
integer, parameter :: i4b = 4
integer, parameter :: sp  = 4
integer, parameter :: dp  = 8

#endif

end module sico_types_m
!
