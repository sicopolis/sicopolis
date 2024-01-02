!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  e r r o r _ m
!
!> @file
!!
!! Writing of error messages and stopping execution.
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
!> Writing of error messages and stopping execution.
!<------------------------------------------------------------------------------
module error_m

  implicit none

  public

contains

!-------------------------------------------------------------------------------
!> Main routine of error_m: Writing of error messages and stopping execution.
!<------------------------------------------------------------------------------
  subroutine error(error_message)

  implicit none

  character(len=256), intent(in) :: error_message

  write(6, fmt='(/,a,/)') trim(error_message)

#if !defined(ALLOW_TAPENADE) /* Normal */
  stop
#else /* Tapenade */
  !!! continue
  !!! (Tapenade cannot deal with stop statements!)
#endif /* Normal vs. Tapenade */

  end subroutine error

!-------------------------------------------------------------------------------

end module error_m
!
