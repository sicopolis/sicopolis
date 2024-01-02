!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  m e t r i c _ m
!
!> @file
!!
!! Definition of the components g11 and g22 of the metric tensor of the
!! applied numerical coordinates.
!!
!! @section Copyright
!!
!! Copyright 2009-2024 Ralf Greve, Roland Warner
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
!> Definition of the components g11 and g22 of the metric tensor of the
!! applied numerical coordinates.
!<------------------------------------------------------------------------------
module metric_m

  use sico_types_m
  use error_m

  implicit none

  public :: metric

contains

!-------------------------------------------------------------------------------
!> Main routine of module metric_m:
!! Definition of the components g11 and g22 of the metric tensor of the
!! applied numerical coordinates.
!<------------------------------------------------------------------------------
  subroutine metric()

  use sico_variables_m, only : xi, eta, &
                               sq_g11_g, sq_g22_g, insq_g11_g, insq_g22_g, &
                               sq_g11_sgx, sq_g22_sgx, insq_g11_sgx, &
                               sq_g11_sgy, sq_g22_sgy, insq_g22_sgy, &
                               PHI0, pi, eps_dp, &
                               errormsg

  implicit none

  integer(i4b) :: i, j
  real(dp)     :: K
  real(dp)     :: g11_g(0:JMAX,0:IMAX), g22_g(0:JMAX,0:IMAX), &
                  g11_sgx(0:JMAX,0:IMAX), g11_sgy(0:JMAX,0:IMAX), &
                  g22_sgx(0:JMAX,0:IMAX), g22_sgy(0:JMAX,0:IMAX)

#if (GRID==0)   /* Stereographic projection (distortion neglected) */

!-------- Components g11, g22 on the grid points (_g) and between
!         the grid points (_sg) --------

  g11_g   = 1.0_dp
  g22_g   = 1.0_dp
  g11_sgx = 1.0_dp
  g11_sgy = 1.0_dp
  g22_sgx = 1.0_dp
  g22_sgy = 1.0_dp

#elif (GRID==1)   /* Stereographic projection */

  if (PHI0 > eps_dp) then   ! for northern hemisphere
     K = (cos(0.25_dp*pi-0.5_dp*PHI0))**2
  else if (PHI0 < (-eps_dp)) then   ! for southern hemisphere
     K = (cos(0.25_dp*pi+0.5_dp*PHI0))**2
  else
     errormsg = ' >>> metric: PHI0 must be different from zero!'
     call error(errormsg)
  end if

!-------- Components g11, g22 on the grid points (_g) --------

  do i=0, IMAX
  do j=0, JMAX
     call metric_stereo(xi(i), eta(j), K, g11_g(j,i), g22_g(j,i))
  end do
  end do

!-------- Components g11, g22 between the grid points (_sg) --------

  do i=0, IMAX-1
  do j=0, JMAX
     call metric_stereo(0.5_dp*(xi(i)+xi(i+1)), eta(j), K, &
                        g11_sgx(j,i), g22_sgx(j,i))
  end do
  end do

  do i=0, IMAX
  do j=0, JMAX-1
     call metric_stereo(xi(i), 0.5_dp*(eta(j)+eta(j+1)), K, &
                        g11_sgy(j,i), g22_sgy(j,i))
  end do
  end do

#elif (GRID==2)   /* Geographical coordinates */

!-------- Components g11, g22 on the grid points (_g) --------

  do i=0, IMAX
  do j=0, JMAX
     call metric_geogr(eta(j), g11_g(j,i), g22_g(j,i))
  end do
  end do

!-------- Components g11, g22 between the grid points (_sg) --------

  do i=0, IMAX-1
  do j=0, JMAX
     call metric_geogr(eta(j), g11_sgx(j,i), g22_sgx(j,i))
  end do
  end do

  do i=0, IMAX
  do j=0, JMAX-1
     call metric_geogr(0.5_dp*(eta(j)+eta(j+1)), g11_sgy(j,i), g22_sgy(j,i))
  end do
  end do

#endif

!-------- Square roots (sq_) and inverse square roots (insq_) of
!         g11 and g22 --------

  do i=0, IMAX
  do j=0, JMAX
     sq_g11_g(j,i)   = sqrt(g11_g(j,i))
     sq_g22_g(j,i)   = sqrt(g22_g(j,i))
     insq_g11_g(j,i) = 1.0_dp/sq_g11_g(j,i)
     insq_g22_g(j,i) = 1.0_dp/sq_g22_g(j,i)
  end do
  end do

  do i=0, IMAX-1
  do j=0, JMAX
     sq_g11_sgx(j,i)   = sqrt(g11_sgx(j,i))
     sq_g22_sgx(j,i)   = sqrt(g22_sgx(j,i))
     insq_g11_sgx(j,i) = 1.0_dp/sq_g11_sgx(j,i)
  end do
  end do

  do i=0, IMAX
  do j=0, JMAX-1
     sq_g11_sgy(j,i)   = sqrt(g11_sgy(j,i))
     sq_g22_sgy(j,i)   = sqrt(g22_sgy(j,i))
     insq_g22_sgy(j,i) = 1.0_dp/sq_g22_sgy(j,i)
  end do
  end do

  end subroutine metric

!-------------------------------------------------------------------------------
!> Components g11 and g22 of the metric tensor for the
!! stereographical projection.
!<------------------------------------------------------------------------------
  subroutine metric_stereo(x_val, y_val, K, g11_r, g22_r)

  use sico_variables_m, only : R

  implicit none

  real(dp), intent(in)  :: x_val, y_val
  real(dp), intent(in)  :: K
  real(dp), intent(out) :: g11_r, g22_r

  g11_r = 1.0_dp / ( K**2*(1.0_dp+(x_val**2+y_val**2)/(2.0_dp*R*K)**2)**2 )

  g22_r = g11_r

  end subroutine metric_stereo

!-------------------------------------------------------------------------------
!> Components g11 and g22 of the metric tensor for geographical coordinates.
!<------------------------------------------------------------------------------
  subroutine metric_geogr(phi_val, g11_r, g22_r)

  use sico_variables_m, only : R

  implicit none

  real(dp), intent(in)  :: phi_val
  real(dp), intent(out) :: g11_r, g22_r

  g11_r = R**2*(cos(phi_val))**2

  g22_r = R**2

  end subroutine metric_geogr

!-------------------------------------------------------------------------------

end module metric_m
!
