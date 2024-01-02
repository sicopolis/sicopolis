!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s t e r e o _ p r o j _ m
!
!> @file
!!
!! Computation of the forward or inverse stereographic projection,
!! alternatively for a spherical or an ellipsoidal planet.
!!
!! @section Copyright
!!
!! Copyright 2009-2024 Ralf Greve, Reinhard Calov, Alex Robinson
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
!> Computation of the forward or inverse stereographic projection,
!! alternatively for a spherical or an ellipsoidal planet.
!<------------------------------------------------------------------------------
module stereo_proj_m

  use sico_types_m
  use sico_variables_m
  use error_m

  implicit none

  public :: stereo_forw_ellipsoid, stereo_inv_ellipsoid, &
            stereo_forw_sphere, stereo_inv_sphere

contains

!-------------------------------------------------------------------------------
!> Forward stereographic projection for an ellipsoidal planet.
!<------------------------------------------------------------------------------
  subroutine stereo_forw_ellipsoid(lon_val, lat_val, A, B, &
                                   lon0, lat0, x_val, y_val)

  implicit none

  real(dp), intent(in)  :: lon_val, lat_val, A, B, lon0, lat0
  real(dp), intent(out) :: x_val, y_val
  
  integer(i4b) :: l
  integer(i4b) :: sign_lat0
  
  real(dp) :: latd0, lat0_aux, lat_aux
  real(dp) :: e, mc, t, tc, kp, rho, lat_p
  real(dp) :: sinlat0, coslat0

!-------- Parameters, constraints --------

  latd0 = lat0*rad2deg

  if (lat0 > eps) then   ! for northern hemisphere
     sign_lat0 =  1
     latd0 = min(latd0,  (90.0_dp-eps))
  else if (lat0 < (-eps)) then   ! for southern hemisphere
     sign_lat0 = -1
     latd0 = max(latd0, -(90.0_dp-eps))
  else
     errormsg = ' >>> stereo_forw_ellipsoid: lat0 must be different from zero!'
     call error(errormsg)
  end if

  lat0_aux = (latd0*deg2rad) * sign_lat0

  e = sqrt((A**2-B**2)/(A**2))

!-------- Stereographic coordinates x,y --------

  sinlat0 = sin(lat0_aux)
  coslat0 = cos(lat0_aux)
  
  lat_aux = lat_val * sign_lat0

  mc = coslat0/sqrt(1.0_dp-e*e*sinlat0*sinlat0)
  t = sqrt(((1.0_dp-sin(lat_aux))/(1.0_dp+sin(lat_aux)))* &
          ((1.0_dp+e*sin(lat_aux))/(1.0_dp-e*sin(lat_aux)))**e)
  tc = sqrt(((1.0_dp-sinlat0)/(1.0_dp+sinlat0))* &
           ((1.0_dp+e*sinlat0)/(1.0_dp-e*sinlat0))**e)
  rho = A*mc*t/tc

  x_val =              rho*sin(lon_val-lon0)
  y_val = -sign_lat0 * rho*cos(lon_val-lon0)

  end subroutine stereo_forw_ellipsoid

!-------------------------------------------------------------------------------
!> Inverse stereographic projection for an ellipsoidal planet.
!<------------------------------------------------------------------------------
  subroutine stereo_inv_ellipsoid(x_val, y_val, A, B, &
                                  lon0, lat0, lon_val, lat_val)

  implicit none

  real(dp), intent(in)  :: x_val, y_val, A, B, lon0, lat0
  real(dp), intent(out) :: lon_val, lat_val

  integer(i4b) :: l
  integer(i4b) :: sign_lat0

  real(dp) :: latd0, lat0_aux, lat_aux
  real(dp) :: e, mc, t, tc, kp, rho, lat_p, residual
  real(dp) :: sinlat0, coslat0

  real(dp), parameter :: eps_residual = 1.0e-09_dp

!-------- Parameters, constraints --------

  latd0 = lat0*rad2deg

  if (lat0 > eps) then   ! for northern hemisphere
     sign_lat0 =  1
     latd0 = min(latd0,  (90.0_dp-eps))
  else if (lat0 < (-eps)) then   ! for southern hemisphere
     sign_lat0 = -1
     latd0 = max(latd0, -(90.0_dp-eps))
  else
     errormsg = ' >>> stereo_inv_ellipsoid: lat0 must be different from zero!'
     call error(errormsg)
  end if

  lat0_aux = (latd0*deg2rad) * sign_lat0

  e = sqrt((A**2-B**2)/(A**2))

!-------- Longitude --------

  if ((x_val /= 0.0_dp).or.(y_val /= 0.0_dp)) then
     lon_val = lon0 + sign_lat0*atan2(y_val,x_val) + 0.5_dp*pi
  else
     lon_val = lon0 + 0.5_dp*pi
  end if

!-------- Fix-point iteration for latitude --------

  sinlat0 = sin(lat0_aux)
  coslat0 = cos(lat0_aux)
  
  tc = sqrt(((1.0_dp-sinlat0)/(1.0_dp+sinlat0))* &
           ((1.0_dp+e*sinlat0)/(1.0_dp-e*sinlat0))**e)
  mc = coslat0/sqrt(1.0_dp-e*e*sinlat0*sinlat0)
  rho = sqrt(x_val*x_val+y_val*y_val)
  t = rho*tc/(A*mc)

  lat_p = 0.5_dp*pi-2.0_dp*atan(t)
  l = 0
  residual = 3600.0_dp

  do while(residual >= eps_residual)
     l = l+1
     lat_aux = 0.5_dp*pi-2.0_dp*atan(t*((1.0_dp-e*sin(lat_p))/ &
               (1.0_dp+e*sin(lat_p)))**(0.5_dp*e))
     residual = abs(lat_aux-lat_p)
     lat_p = lat_aux
  end do

  lat_val = lat_aux * sign_lat0

!-------- Constrain longitude to [0, 2*pi) --------

  if (lon_val < 0.0_dp) then
     lon_val = lon_val + 2.0_dp*pi
  else if (lon_val >= (2.0_dp*pi)) then
     lon_val = lon_val - 2.0_dp*pi
  end if

  end subroutine stereo_inv_ellipsoid

!-------------------------------------------------------------------------------
!> Forward stereographic projection for a spherical planet.
!<------------------------------------------------------------------------------
  subroutine stereo_forw_sphere(lon_val, lat_val, R, lon0, lat0, &
                                x_val, y_val)

  implicit none

  real(dp), intent(in)  :: lon_val, lat_val, R, lon0, lat0
  real(dp), intent(out) :: x_val, y_val

  integer(i4b) :: sign_lat0

  real(dp) :: latd0, lat0_aux, lat_aux
  real(dp) :: K

!-------- Parameters, constraints --------

  latd0 = lat0*rad2deg

  if (lat0 > eps) then   ! for northern hemisphere
     sign_lat0 =  1
     latd0 = min(latd0,  90.0_dp)
  else if (lat0 < (-eps)) then   ! for southern hemisphere
     sign_lat0 = -1
     latd0 = max(latd0, -90.0_dp)
  else
     errormsg = ' >>> stereo_forw_sphere: lat0 must be different from zero!'
     call error(errormsg)
  end if

  lat0_aux = (latd0*deg2rad) * sign_lat0

  K = (cos(0.25_dp*pi-0.5_dp*lat0_aux))**2

!-------- Stereographic coordinates x,y --------

  lat_aux = lat_val * sign_lat0

  x_val =              2.0_dp*R*K*tan(0.25_dp*pi-0.5_dp*lat_aux) &
                                  *sin(lon_val-lon0)

  y_val = -sign_lat0 * 2.0_dp*R*K*tan(0.25_dp*pi-0.5_dp*lat_aux) &
                                 *cos(lon_val-lon0)

  end subroutine stereo_forw_sphere

!-------------------------------------------------------------------------------
!> Inverse stereographic projection for a spherical planet.
!<------------------------------------------------------------------------------
  subroutine stereo_inv_sphere(x_val, y_val, R, lon0, lat0, &
                               lon_val, lat_val)

  implicit none

  real(dp), intent(in)  :: x_val, y_val, R, lon0, lat0
  real(dp), intent(out) :: lon_val, lat_val

  integer(i4b) :: sign_lat0

  real(dp) :: latd0, lat0_aux, lat_aux
  real(dp) :: K

!-------- Parameters, constraints --------

  latd0 = lat0*rad2deg

  if (lat0 > eps) then   ! for northern hemisphere
     sign_lat0 =  1
     latd0 = min(latd0,  90.0_dp)
  else if (lat0 < (-eps)) then   ! for southern hemisphere
     sign_lat0 = -1
     latd0 = max(latd0, -90.0_dp)
  else
     errormsg = ' >>> stereo_inv_sphere: lat0 must be different from zero!'
     call error(errormsg)
  end if

  lat0_aux = (latd0*deg2rad) * sign_lat0

  K = (cos(0.25_dp*pi-0.5_dp*lat0_aux))**2

!-------- Longitude --------

  if ((x_val /= 0.0_dp).or.(y_val /= 0.0_dp)) then
     lon_val = lon0 + sign_lat0*atan2(y_val,x_val) + 0.5_dp*pi
  else
     lon_val = lon0 + 0.5_dp*pi
  end if

!-------- Latitude --------

  lat_aux = 0.5_dp*pi &
            -2.0_dp*atan(sqrt(x_val**2+y_val**2)/(2.0_dp*R*K))

  lat_val = lat_aux * sign_lat0

!-------- Constrain longitude to [0, 2*pi) --------

  if (lon_val < 0.0_dp) then
     lon_val = lon_val + 2.0_dp*pi
  else if (lon_val >= (2.0_dp*pi)) then
     lon_val = lon_val - 2.0_dp*pi
  end if

  end subroutine stereo_inv_sphere

!-------------------------------------------------------------------------------

end module stereo_proj_m
!
