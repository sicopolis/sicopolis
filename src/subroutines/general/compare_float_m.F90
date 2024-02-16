!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c o m p a r e _ f l o a t _ m
!
!! Comparison of single- or double-precision floating-point numbers.
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
!> Comparison of single- or double-precision floating-point numbers.
!-------------------------------------------------------------------------------
module compare_float_m

#if (defined(MODEL_SICOPOLIS))
  use sico_types_m, only: sp, dp
#endif

  implicit none

  private
  public :: approx_equal, approx_equal_integer, approx_integer_multiple

#if (!defined(MODEL_SICOPOLIS))
  integer, parameter :: sp  = kind(1.0)     ! single precision
  integer, parameter :: dp  = kind(1.0d0)   ! double precision
#endif

  interface approx_equal
     module procedure approx_equal_sp
     module procedure approx_equal_dp
  end interface

  interface approx_equal_integer
     module procedure approx_equal_integer_sp
     module procedure approx_equal_integer_dp
  end interface

  interface approx_integer_multiple
     module procedure approx_integer_multiple_sp
     module procedure approx_integer_multiple_dp
  end interface

contains

!-------------------------------------------------------------------------------
!> Check whether single-precision x and y are approximately equal
!! (within a given eps limit).
!-------------------------------------------------------------------------------
  function approx_equal_sp(x, y, eps)

  implicit none

  real(sp), intent(in) :: x, y
  real(sp), intent(in) :: eps

  logical :: approx_equal_sp

  if ( abs(x-y) <= abs(x+y)*eps ) then
     approx_equal_sp = .true.
  else
     approx_equal_sp = .false.
  end if

  end function approx_equal_sp

!-------------------------------------------------------------------------------
!> Check whether double-precision x and y are approximately equal
!! (within a given eps limit).
!-------------------------------------------------------------------------------
  function approx_equal_dp(x, y, eps)

  implicit none

  real(dp), intent(in) :: x, y
  real(dp), intent(in) :: eps

  logical :: approx_equal_dp

  if ( abs(x-y) <= abs(x+y)*eps ) then
     approx_equal_dp = .true.
  else
     approx_equal_dp = .false.
  end if

  end function approx_equal_dp

!-------------------------------------------------------------------------------
!> Check whether single-precision x is approximately equal to an integer
!! (within a given eps limit).
!-------------------------------------------------------------------------------
  function approx_equal_integer_sp(x, eps)

  implicit none

  real(sp), intent(in) :: x
  real(sp), intent(in) :: eps

  logical :: approx_equal_integer_sp

  if (x == 0.0_sp) then
     approx_equal_integer_sp = .true.
  else if (approx_equal_sp(x, real(nint(x),sp), eps)) then
     approx_equal_integer_sp = .true.
  else
     approx_equal_integer_sp = .false.
  end if

  end function approx_equal_integer_sp

!-------------------------------------------------------------------------------
!> Check whether double-precision x is approximately equal to an integer
!! (within a given eps limit).
!-------------------------------------------------------------------------------
  function approx_equal_integer_dp(x, eps)

  implicit none

  real(dp), intent(in) :: x
  real(dp), intent(in) :: eps

  logical :: approx_equal_integer_dp

  if (x == 0.0_dp) then
     approx_equal_integer_dp = .true.
  else if (approx_equal_dp(x, real(nint(x),dp), eps)) then
     approx_equal_integer_dp = .true.
  else
     approx_equal_integer_dp = .false.
  end if

  end function approx_equal_integer_dp

!-------------------------------------------------------------------------------
!> Check whether single-precision x is approximately an integer multiple of
!! single-precision y (within a given eps limit).
!-------------------------------------------------------------------------------
  function approx_integer_multiple_sp(x, y, eps)

  implicit none

  real(sp), intent(in) :: x, y
  real(sp), intent(in) :: eps

  logical :: approx_integer_multiple_sp

  if (y == 0.0_sp) then
     approx_integer_multiple_sp = .false.
  else if (nint(x/y) == 0) then
     approx_integer_multiple_sp = .false.
  else if (approx_equal_integer_sp(x/y, eps)) then
     approx_integer_multiple_sp = .true.
  else
     approx_integer_multiple_sp = .false.
  end if

  end function approx_integer_multiple_sp

!-------------------------------------------------------------------------------
!> Check whether double-precision x is approximately an integer multiple of
!! double-precision y (within a given eps limit).
!-------------------------------------------------------------------------------
  function approx_integer_multiple_dp(x, y, eps)

  implicit none

  real(dp), intent(in) :: x, y
  real(dp), intent(in) :: eps

  logical :: approx_integer_multiple_dp

  if (y == 0.0_dp) then
     approx_integer_multiple_dp = .false.
  else if (nint(x/y) == 0) then
     approx_integer_multiple_dp = .false.
  else if (approx_equal_integer_dp(x/y, eps)) then
     approx_integer_multiple_dp = .true.
  else
     approx_integer_multiple_dp = .false.
  end if

  end function approx_integer_multiple_dp

!-------------------------------------------------------------------------------

end module compare_float_m
!
