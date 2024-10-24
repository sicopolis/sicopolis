!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ e n d _ m
!
!! Ending of SICOPOLIS.
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
!> Ending of SICOPOLIS.
!-------------------------------------------------------------------------------
module sico_end_m

  use sico_types_m
  use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

  implicit none

  public

contains

!-------------------------------------------------------------------------------
!> Main routine of sico_end_m: Ending of SICOPOLIS.
!-------------------------------------------------------------------------------
  subroutine sico_end()

  use netcdf
  use nc_check_m

  implicit none

  integer(i4b) :: n
  integer(i4b) :: ierr

  character(len=64), parameter :: thisroutine = 'sico_end'

  write(unit=6, fmt='(/a/)') ' -------- sico_end --------'

  close(unit=12, status='keep')  ! Close time-series files
  close(unit=14, status='keep')

#if (defined(ASF) && WRITE_SER_FILE_STAKES==1) /* Austfonna */
  close(unit=41, status='keep')
  close(unit=42, status='keep')
  close(unit=43, status='keep')
  close(unit=44, status='keep')
  close(unit=45, status='keep')
  close(unit=46, status='keep')
  close(unit=47, status='keep')
  close(unit=48, status='keep')
  close(unit=49, status='keep')
  close(unit=50, status='keep')
#endif

#if (defined(HEINO))
  close(unit=15, status='keep')
#endif

#if !(defined(ALLOW_GRDCHK) || defined(ALLOW_TAPENADE))

  do n=0, maxval(mask_region)
     call check( nf90_sync(ncid_ser(n)),  thisroutine )
     call check( nf90_close(ncid_ser(n)), thisroutine )
          ! Closing of NetCDF time-series output files
  end do

  if (n_site >= 1) then
     call check( nf90_sync(ncid_site),  thisroutine )
     call check( nf90_close(ncid_site), thisroutine )
          ! Closing of NetCDF time-series output file for the specified sites
          ! (i.e., ice cores)
  end if

#endif

#if (MARGIN==3 || DYNAMICS==2 || DYNAMICS==3)
#if !defined(ALLOW_TAPENADE)
  call lis_finalize(ierr)   ! Finalise execution environment of the
                            ! Library of Iterative Solvers Lis, if required
#else
  call lis_finalize_f(ierr) ! Finalise execution environment of the
                            ! Library of Iterative Solvers Lis, if required
#endif
#endif

  write(unit=6, fmt='(/a/)') '    * * * sicopolis.F90  r e a d y * * *'

  end subroutine sico_end

!-------------------------------------------------------------------------------

end module sico_end_m
!
