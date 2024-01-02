!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  m a r s _ i n s t e m p _ m
!
!> @file
!!
!! Computation of the daily mean surface temperature of Mars based on
!! obliquity, eccentricity and the anomaly of vernal equinox
!! (local insolation temperature scheme = LIT scheme).
!!
!! @section Copyright
!!
!! Copyright 2009-2024 Bjoern Grieger, Oliver J. Stenzel, Ralf Greve
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
!> Computation of the daily mean surface temperature of Mars based on
!! obliquity, eccentricity and the anomaly of vernal equinox
!! (local insolation temperature scheme = LIT scheme).
!<------------------------------------------------------------------------------
module mars_instemp_m

   use sico_types_m

   implicit none

   public

   !> Martian surface temperatures
   type ins
      !> Temperature as a function of time and latitude
      real(dp) :: t(0:360,-90:90)
      !> Annual mean temperature as a function of latitude
      real(dp) :: tam(-90:90)
      !> Annual maximum temperature as a function of latitude
      real(dp) :: tmax(-90:90)
   end type ins

contains

!-------------------------------------------------------------------------------
!> Main subroutine of module mars_instemp_m
!<------------------------------------------------------------------------------
   subroutine setinstemp ( o, ecc, ave, obl, sma, sa, sac, op, ct )

      implicit none

      type(ins)          :: o
      real(dp), optional :: ecc   ! eccentricity
      real(dp), optional :: ave   ! anom. of vernal equinox in degree
      real(dp), optional :: obl   ! obliquity in degree
      real(dp), optional :: sma   ! semi-major axis in AU
      real(dp), optional :: sa    ! surface albedo (planetary mean)
      real(dp), optional :: sac   ! surface albedo (seasonal CO2 ice caps)
      real(dp), optional :: op    ! orbital period in seconds
      real(dp), optional :: ct    ! CO2 condensation temperature in K

      logical      :: co2layer(-90:90)
      integer(i4b) :: iphi, ipsi, year
      integer(i4b), parameter :: yearmax = 2
      real(dp), parameter :: pi = 3.141592653589793_dp, deg2rad = pi/180.0_dp
      real(dp), parameter :: SB = 5.67051e-8_dp
                             ! Wikipedia gives 5.670400E-8 +- 0.000040E-8
      real(dp) :: tptd, seps, sdelta, cdelta
      real(dp) :: albact, albact_co2, delta, du, e, eps, f, &
                  lambda, phi, psi, psi0, r
      real(dp) :: tau0, teq, ufac, usum, w, wg
      real(dp), dimension(-90:90) :: co2, t

      real(dp) :: j0
      real(dp) :: a
      real(dp) :: alb, alb_co2
      real(dp) :: u
      real(dp) :: tco2

      j0   = 1367.6_dp   ! solar constant for Earth in W/m**2

      if ( present(ecc) ) then
         e = ecc
      else
         e = 0.0935_dp
      end if
      if ( present(ave) ) then
         psi0 = ave
      else
         psi0 = 109.13_dp
      end if
      if ( present(obl) ) then
         eps = obl
      else
         eps = 25.19_dp
      end if
      if ( present(sma) ) then
         a = sma
      else
         a = 1.524_dp
      end if
      if ( present(sa) ) then
         alb = sa
      else
         alb = 0.3_dp
      end if
      if ( present(sac) ) then
         alb_co2 = sac
      else
         alb_co2 = alb
      end if
      if ( present(op) ) then
         u = op
      else
         u = 686.95_dp*24._dp*3600._dp
      end if
      if ( present(ct) ) then
         tco2 = ct
      else
         tco2 = 146._dp
      end if

      j0 = j0 / a**2
      f = a * e
      psi0 = psi0 * deg2rad
      eps = eps * deg2rad
      
      seps = sin(eps)

      usum = 0.0_dp
      do ipsi = 0, 360, 1
         psi = ipsi * deg2rad
         r = ( a**2 - f**2 ) / ( a + f*cos(psi) )
         du = r**2
         if (ipsi<360) usum = usum + du
      end do
      ufac = u / usum

      t = 273.15_dp
      co2 = 0.0_dp
      co2layer = .false.
      do year = 1, yearmax
         usum = 0.0_dp
         do ipsi = 0, 360, 1
            psi = ipsi * deg2rad
            r = ( a**2 - f**2 ) / ( a + f*cos(psi) )
            du = ufac * r**2
            if (ipsi<360) usum = usum + du
            lambda = psi - psi0
            delta  = asin( seps * sin(lambda) )
            sdelta = sin(delta)
            cdelta = cos(delta)
            do iphi = -89, 89, 1
               phi = iphi * deg2rad
               tptd = -tan(phi) * tan(delta)
               if ( tptd .le. -1.0_dp ) then
                  tau0 = pi
               else if ( tptd .ge. 1.0_dp ) then
                  tau0 = 0.0_dp
               else
                  tau0 = acos( tptd )
               end if
               w = ( tau0*sin(phi)*sdelta + &
                     cos(phi)*cdelta*sin(tau0) ) &
                   * j0 * a**2 / ( pi * r**2 )
               wg = 0.0_dp ! 0.045 * j0
               if ( t(iphi) < -999._dp ) then
                  albact     = 0.95_dp
                  albact_co2 = 0.95_dp
               else
                  albact     = alb
                  albact_co2 = alb_co2
               end if
               teq = ( ( wg + (1._dp-albact) * w ) / SB )**.25_dp
               if ( teq < tco2 ) then
                  co2layer(iphi) = .true.
                  co2(iphi) = co2(iphi) + SB * tco2**4 * du &
                              - (1._dp-albact_co2) * w * du
                  t(iphi) = tco2
               else
                  if ( co2layer(iphi) ) then
                     co2(iphi) = co2(iphi) + SB * tco2**4 * du &
                                 - (1._dp-albact_co2) * w * du
                     if ( co2(iphi) .le. 0.0_dp ) then
                        co2(iphi) = 0.0_dp
                        co2layer(iphi) = .false.
                        t(iphi) = tco2
                     end if
                  else
                     t(iphi) = teq
                  end if
               end if
               if ( year .eq. yearmax ) o%t(ipsi,iphi) = t(iphi)
            end do
            if ( year .eq. yearmax .and. ipsi .eq. 0 ) then
               o%tam = 0.0_dp
               o%tmax = 0.0_dp
            else
               o%tam = o%tam + t
               do iphi = -89, 89
                  o%tmax(iphi) = max( o%tmax(iphi), t(iphi) )
               end do
            end if
         end do
         o%tam = o%tam / 360._dp
      end do
      o%t(:,-90)  = o%t(:,-89)  + ( o%t(:,-89)  - o%t(:,-88) )  / 2._dp
      o%t(:, 90)  = o%t(:, 89)  + ( o%t(:, 89)  - o%t(:, 88) )  / 2._dp
      o%tam(-90)  = o%tam(-89)  + ( o%tam(-89)  - o%tam(-88) )  / 2._dp
      o%tam( 90)  = o%tam( 89)  + ( o%tam( 89)  - o%tam( 88) )  / 2._dp
      o%tmax(-90) = o%tmax(-89) + ( o%tmax(-89) - o%tmax(-88) ) / 2._dp
      o%tmax( 90) = o%tmax( 89) + ( o%tmax( 89) - o%tmax( 88) ) / 2._dp

   end subroutine

!-------------------------------------------------------------------------------
!> Annual mean temperature at latitude phi
!<------------------------------------------------------------------------------
   function instam ( o, phi )

      implicit none
      real(dp)     :: instam
      type(ins)    :: o
      real(dp)     :: phi

      integer(i4b) :: iphi1, iphi2

      iphi1 = nint( phi - 0.5_dp )
      if ( iphi1 < -90 ) then
         iphi1 = -90
      else if ( iphi1 > 89 ) then
         iphi1 = 89
      end if
      iphi2 = iphi1 + 1
      instam = o%tam(iphi1) + ( o%tam(iphi2) - o%tam(iphi1) ) &
               * ( phi - iphi1 )

   end function instam

!-------------------------------------------------------------------------------
!> Annual maximum temperature at latitude phi
!<------------------------------------------------------------------------------
   function instmax ( o, phi )

      implicit none
      real(dp)     :: instmax
      type(ins)    :: o
      real(dp)     :: phi

      integer(i4b) :: iphi1, iphi2

      iphi1 = nint( phi - 0.5_dp )
      if ( iphi1 < -90 ) then
         iphi1 = -90
      else if ( iphi1 > 89 ) then
         iphi1 = 89
      end if
      iphi2 = iphi1 + 1
      instmax = o%tmax(iphi1) + ( o%tmax(iphi2) - o%tmax(iphi1) ) &
               * ( phi - iphi1 )

   end function instmax

!-------------------------------------------------------------------------------
!> Temperature at orbit position psi and latitude phi
!<------------------------------------------------------------------------------
   function inst ( o, psi, phi )

      implicit none
      real(dp)     :: inst
      type(ins)    :: o
      real(dp)     :: psi, phi

      integer(i4b) :: ipsi1, ipsi2

      ipsi1 = nint( psi - 0.5_dp )
      if ( ipsi1 < 0 ) then
         ipsi1 = 0
      else if ( ipsi1 > 359 ) then
         ipsi1 = 359
      end if
      ipsi2 = ipsi1 + 1
      inst = inst1(o,ipsi1,phi) + ( inst1(o,ipsi2,phi) - inst1(o,ipsi1,phi) ) &
               * ( psi - ipsi1 )

   end function inst

!-------------------------------------------------------------------------------
!> Temperature at orbit position ipsi (integer) and latitude phi
!<------------------------------------------------------------------------------
   function inst1 ( o, ipsi, phi )

      implicit none
      real(dp)     :: inst1
      type(ins)    :: o
      integer(i4b) :: ipsi
      real(dp)     :: phi

      integer(i4b) :: iphi1, iphi2

      iphi1 = nint( phi - 0.5_dp )
      if ( iphi1 < -90 ) then
         iphi1 = -90
      else if ( iphi1 > 89 ) then
         iphi1 = 89
      end if
      iphi2 = iphi1 + 1
      inst1 = o%t(ipsi,iphi1) + ( o%t(ipsi,iphi2) - o%t(ipsi,iphi1) ) &
               * ( phi - iphi1 )

   end function inst1

end module mars_instemp_m
!
