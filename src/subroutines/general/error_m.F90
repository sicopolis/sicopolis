!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  e r r o r _ m
!
!! Handling of error or warning messages, stopping execution if applicable.
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
!> Handling of error or warning messages, stopping execution if applicable.
!-------------------------------------------------------------------------------
module error_m

  implicit none

  public

contains

!-------------------------------------------------------------------------------
!> Writing of error messages and stopping execution.
!-------------------------------------------------------------------------------
  subroutine error(error_message)

  implicit none

  character(len=256), intent(in) :: error_message

  character(len=256)   :: error_msg
  character, parameter :: end_of_line = char(10)
                          ! End-of-line string

  error_msg = ' ERROR:' &
              // end_of_line &
              // trim(error_message)

  write(6, fmt='(/,a,/)') trim(error_msg)

#if !defined(ALLOW_TAPENADE) /* Normal */
  stop
#else /* Tapenade */
  !%% continue
  !%% (Tapenade cannot deal with stop statements!)
#endif /* Normal vs. Tapenade */

  end subroutine error

!-------------------------------------------------------------------------------
!> Writing of warning messages.
!-------------------------------------------------------------------------------
  subroutine warning(warning_message)

  implicit none

  character(len=256), intent(in) :: warning_message

  character(len=256)   :: warning_msg
  character, parameter :: end_of_line = char(10)
                          ! End-of-line string

  warning_msg = ' Warning:' &
                // end_of_line &
                // trim(warning_message)

  write(6, fmt='(/,a)') trim(warning_msg)

  end subroutine warning

!-------------------------------------------------------------------------------

end module error_m
!
