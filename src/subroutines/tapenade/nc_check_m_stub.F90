!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  n c _ c h e c k _ m _ s t u b
!
!> @file
!!
!! Stub file for NetCDF error capturing.
!!
!! @section Copyright
!!
!! Copyright 2009-2024 Ralf Greve, Shreyas Sunil Gaikwad
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
!> Stub file for NetCDF error capturing.
!<------------------------------------------------------------------------------
module nc_check_m

use sico_types_m
use error_m

use netcdf

contains

!-------------------------------------------------------------------------------
!> NetCDF error capturing.
!<------------------------------------------------------------------------------
subroutine check(status, ch_calling_routine)

implicit none

integer(i4b),                intent(in) :: status
character(len=64), optional, intent(in) :: ch_calling_routine

character(len=64)  :: ch_clrt
character(len=256) :: errormsgg

character, parameter :: ch_end_of_line = char(10)

end subroutine check

end module nc_check_m
!
