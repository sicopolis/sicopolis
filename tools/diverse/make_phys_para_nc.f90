!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Program   : m a k e _ p h y s _ p a r a _ n c
!
!> @file
!!
!! Generating NetCDF files for the physical parameters required by SICOPOLIS.
!!
!! @section Date
!!
!! 2023-11-09
!!
!! @section Copyright
!!
!! Copyright 2021-2023 Ralf Greve
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

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !
   ! To be compiled with netCDF. Using the SICOPOLIS set-up,
   ! this can be done by executing the following commands:
   !
   !    source ../../runs/sico_configs.sh
   !    $FC make_phys_para_nc.f90 -o make_phys_para_nc.x $FCFLAGS
   !
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!-------------------------------------------------------------------------------
!> Main program:
!! Generating NetCDF files for the physical parameters required by SICOPOLIS.
!<------------------------------------------------------------------------------
program make_phys_para_nc

use netcdf

implicit none

integer, parameter :: i4b = selected_int_kind(9)   ! 4-byte integers
integer, parameter :: dp  = kind(1.0d0)            ! double-precision reals

real(dp) :: RHO, RHO_W, RHO_SW, L, G, NUE, BETA, DELTA_TM_SW, OMEGA_MAX, &
            H_R, RHO_C_R, KAPPA_R, RHO_A, R_T, &
            R, A, B, F_INV, LOND0, LATD0

real(dp), dimension(-190:10) :: RF, KAPPA, C

character(len=256) :: filename_nc, ch_info1, ch_info2, ch_date_time
character(len= 32) :: ch_domain
character(len=  8) :: ch_date
character(len= 10) :: ch_time

integer(i4b) :: n

integer(i4b) :: istat, ncid, ncd
integer(i4b), dimension(100) :: varid

real(dp), dimension(-190:0) :: temp

real(dp), parameter :: temp_C2K = 273.15_dp
real(dp), parameter :: gas_constant = 8.314_dp
real(dp), parameter :: very_large = 1.1111e+11_dp

!-------- Settings --------

ch_domain = 'CALVINGMIP'

! Simulated domain:
!   ANT   - Antarctica
!   ASF   - Austfonna
!   GRL   - Greenland
!   NHEM  - Northern hemisphere
!   SCAND - Scandinavia
!   TIBET - Tibet

!!! !   NMARS - North polar cap of Mars [not yet considered]
!!! !   SMARS - South polar cap of Mars [not yet considered]

!   EISMINT    - EISMINT (Phase 2 SGE and modifications)
!   CALVINGMIP - Calving MIP
!   VIALOV     - 3D Vialov profile
!   HEINO      - ISMIP HEINO

!   MOCHO  - Mocho-Choshuenco ice cap
!   NPI    - North Patagonian ice field
!   SHMARS - Southern hemisphere of Mars

ch_info1 = 'Physical-parameter file for domain ' // trim(ch_domain)
ch_info2 = ''   ! anything to be added can go here

call date_and_time(date=ch_date, time=ch_time)

ch_date_time = trim(ch_date) // 'T' // ch_time(1:6)

filename_nc = 'my_phys_para_' // trim(ch_domain) // '_' &
                              // trim(ch_date_time) // '.nc'

!-------- Values of physical parameters --------

!  ------ General parameters

select case(ch_domain)
   case ('CALVINGMIP')
      RHO = 9.17e+02_dp
      ! Density of ice = 917 kg/m3 
   case default
      RHO = 9.1e+02_dp
      ! Density of ice = 910 kg/m3 
end select

RHO_W = 1.0e+03_dp
! Density of pure water = 1000 kg/m3

RHO_SW = 1.028e+03_dp
! Density of sea water = 1028 kg/m3

L = 3.35e+05_dp
! Latent heat = 3.35e+05 J/kg 

select case(ch_domain)
   case ('SHMARS')
      G = 3.72_dp
      ! Gravity acceleration = 3.72 m/s2
   case default
      G = 9.81_dp
      ! Gravity acceleration = 9.81 m/s2
end select

NUE = 1.0e-06_dp
! Water diffusivity = 1.0e-06 kg/(m*s) 

select case(ch_domain)
   case ('SHMARS')
      BETA = 3.30e-04_dp
      ! Clausius-Clapeyron gradient = 3.30e-04 K/m
   case default
      BETA = 8.7e-04_dp
      ! Clausius-Clapeyron gradient = 8.7e-04 K/m
end select

DELTA_TM_SW = 1.85_dp
! Melting point depression of sea water
! due to its average salinity = 1.85 degC

OMEGA_MAX = 1.0e-02_dp
! Threshold value for the water content = 1!

H_R = 2.0e+03_dp
! Thickness of the modelled lithosphere layer = 2000 m

RHO_C_R = 2.0e+06_dp
! Density times specific heat of the lithosphere = 2.0e+06 J/(m3*K) 

KAPPA_R = 3.0_dp
! Heat conductivity of the lithosphere = 3 W/(m*K) 

RHO_A = 3.3e+03_dp
! Density of the asthenosphere = 3300 kg/m3 

select case(ch_domain)
   case ('EISMINT', 'CALVINGMIP', 'VIALOV', 'HEINO')
      R_T = 0.0_dp
      ! Coefficient of the water-content dependence in the rate factor
      ! for temperate ice = 0
   case default
      R_T = 1.8125e+02_dp
      ! Coefficient of the water-content dependence in the rate factor
      ! for temperate ice = 181.25
end select

!  ------ Stereographic projection

select case(ch_domain)
   case ('SHMARS')
      R = 3.3762e+06_dp
      ! Radius of Mars = 3376200 m
   case default
      R = 6.371e+06_dp
      ! Radius of the Earth = 6371000 m
end select

select case(ch_domain)
   case ('SHMARS')
      A = 3.3762e+06_dp
      ! Semi-major axis of Mars = 3376200 m (sphere assumed)
   case default
      A = 6.378137e+06_dp
      ! Semi-major axis of the Earth = 6378137 m
end select

select case(ch_domain)
   case ('SHMARS')
      F_INV = very_large
      ! Inverse flattening of Mars = inf (sphere assumed)
   case default
      F_INV = 298.257223563_dp
      ! Inverse flattening of the Earth = 298.257223563
end select

select case(ch_domain)
   case ('ANT')
      LATD0 = -71.0_dp
      LOND0 =   0.0_dp
   case ('ASF')
      LATD0 = 79.7208_dp
      LOND0 = 15.0_dp   ! central median of UTM33
   case ('GRL')
      LATD0 =  70.0_dp
      LOND0 = -45.0_dp
   case ('NHEM')
      LATD0 =  71.0_dp
      LOND0 = -44.0_dp
   case ('SCAND')
      LATD0 =  71.0_dp
      LOND0 =  75.0_dp
   case ('TIBET')
      LATD0 =  35.0_dp
      LOND0 =  90.0_dp
   case ('MOCHO')
      LATD0 = -40.0_dp
      LOND0 = -72.5_dp
   case ('NPI')
      LATD0 = -47.0_dp
      LOND0 = -75.0_dp   ! central meridian of UTM S18 (72W-78W)
   case ('SHMARS')
      LATD0 = -60.0_dp
      LOND0 =   0.0_dp
   case ('EISMINT', 'CALVINGMIP', 'VIALOV', 'HEINO')
      LATD0 =  71.0_dp   ! for simple, "academic" geometry
      LOND0 =   0.0_dp   ! for simple, "academic" geometry
   case default
      write(6,'(/a)') 'Warning: LATD0 and LOND0 undefined!'
      LATD0 = very_large   ! undefined
      LOND0 = very_large   ! undefined
end select

! LATD0: standard parallel (degrees N)
!        \!/ In the ASCII files, this is called PHI0 (in rad) \!/
! LOND0: Reference longitude (degrees E)
!        \!/ In the ASCII files, this is called LAMBDA0 (in rad) \!/

!  ------ Rate factor RF(T)
!         [RF] = 1/(s*Pa3), [T] = degC

do n=-190, 0
   temp(n) = temp_C2K + real(n,dp)
end do

select case(ch_domain)

   case ('VIALOV')

!    ---- RF(T) = const

      RF = 3.16887646e-24_dp

   case ('CALVINGMIP')

!    ---- RF(T) = const

      RF = 9.30920838e-26_dp

   case ('EISMINT', 'HEINO')

!    ---- RF(T), set-up of EISMINT Phase 2 SGE

      do n=-190, -11
         RF(n) = 3.61e-13_dp * exp( -6.0e+04_dp/(gas_constant*temp(n)))
      end do

      do n=-10, 0
         RF(n) = 1.73e+03_dp * exp(-13.9e+04_dp/(gas_constant*temp(n)))
      end do

      RF(1:10) = RF(0)   ! Dummies

   case default

!    ---- RF(T) by Cuffey and Paterson, "The Physics of Glaciers",
!                                      4th ed. 2010, Sect. 3.4.6)

      do n=-190, -11
         RF(n) = 3.5e-25_dp &
                    * exp( ( -6.0e+04_dp/gas_constant) &
                           * (1.0_dp/temp(n)-1.0_dp/(temp_C2K-10.0_dp)) )
      end do

      do n=-10, 0
         RF(n) = 3.5e-25_dp &
                    * exp( (-11.5e+04_dp/gas_constant) &
                           * (1.0_dp/temp(n)-1.0_dp/(temp_C2K-10.0_dp)) )
      end do

      RF(1:10) = RF(0)   ! Dummies

end select

!  ------ Heat conductivity KAPPA(T)
!         [KAPPA] = W/(m*K), [T] = degC

select case(ch_domain)

   case ('EISMINT', 'CALVINGMIP', 'VIALOV', 'HEINO')

!    ---- KAPPA(T) = const, by EISMINT Phase 2 SGE

      KAPPA = 2.1_dp

   case default

!    ---- KAPPA(T) by Greve, Weis and Hutter, 1998,
!                     Paleoclimates 2 (2-3), 133-161

      do n=-190, 0
         KAPPA(n) = 9.828_dp*exp(-0.0057_dp*temp(n))
      end do

      KAPPA(1:10) = KAPPA(0)   ! Dummies

end select

!  ------ Specific heat C(T)
!         [C] = J/(kg*K), [T] = degC

select case(ch_domain)

   case ('EISMINT', 'CALVINGMIP', 'VIALOV', 'HEINO')

!    ---- C(T) = const, by EISMINT Phase 2 SGE

      C = 2.009e+03_dp

   case default

!    ---- C(T) by Greve, Weis and Hutter, 1998,
!                 Paleoclimates 2 (2-3), 133-161

      do n=-190, 0
         C(n) = 146.3_dp + 7.253_dp*temp(n)
      end do

      C(1:10) = C(0)   ! Dummies

end select

!-------- Write physical parameters on file --------

!  ------ Create NetCDF file

istat = nf90_create(trim(filename_nc), NF90_NOCLOBBER, ncid)

!  ------ Global attributes

istat = nf90_put_att(ncid, NF90_GLOBAL, 'Info1' , trim(ch_info1))
istat = nf90_put_att(ncid, NF90_GLOBAL, 'Info2' , trim(ch_info2))
istat = nf90_put_att(ncid, NF90_GLOBAL, 'Date' , trim(ch_date_time))

!  ------ Define dimension

istat = nf90_def_dim(ncid, 'tables', size(RF), ncd)

!  ------ Define variables

varid = -9999   ! initialization

istat = nf90_def_var(ncid, 'RHO', NF90_DOUBLE, varid(1))

istat = nf90_put_att(ncid, varid(1), 'units', 'kg m-3')
istat = nf90_put_att(ncid, varid(1), 'standard_name', 'ice_density')
istat = nf90_put_att(ncid, varid(1), 'long_name', 'Density of ice')

istat = nf90_def_var(ncid, 'RHO_W', NF90_DOUBLE, varid(2))

istat = nf90_put_att(ncid, varid(2), 'units', 'kg m-3')
istat = nf90_put_att(ncid, varid(2), 'standard_name', 'water_density')
istat = nf90_put_att(ncid, varid(2), 'long_name', 'Density of pure water')

istat = nf90_def_var(ncid, 'RHO_SW', NF90_DOUBLE, varid(3))

istat = nf90_put_att(ncid, varid(3), 'units', 'kg m-3')
istat = nf90_put_att(ncid, varid(3), 'standard_name', 'sea_water_density')
istat = nf90_put_att(ncid, varid(3), 'long_name', 'Density of sea water')

istat = nf90_def_var(ncid, 'L', NF90_DOUBLE, varid(4))

istat = nf90_put_att(ncid, varid(4), 'units', 'J kg-1')
istat = nf90_put_att(ncid, varid(4), 'standard_name', 'ice_latent_heat')
istat = nf90_put_att(ncid, varid(4), 'long_name', 'Latent heat of ice')

istat = nf90_def_var(ncid, 'G', NF90_DOUBLE, varid(5))

istat = nf90_put_att(ncid, varid(5), 'units', 'm s-2')
istat = nf90_put_att(ncid, varid(5), 'standard_name', 'gravity_acceleration')
istat = nf90_put_att(ncid, varid(5), 'long_name', 'Gravity acceleration')

istat = nf90_def_var(ncid, 'NUE', NF90_DOUBLE, varid(6))

istat = nf90_put_att(ncid, varid(6), 'units', 'kg m-1 s-1')
istat = nf90_put_att(ncid, varid(6), 'standard_name', 'ice_water_diffusivity')
istat = nf90_put_att(ncid, varid(6), 'long_name', 'Water diffusivity in ice')

istat = nf90_def_var(ncid, 'BETA', NF90_DOUBLE, varid(7))

istat = nf90_put_att(ncid, varid(7), 'units', 'K m-1')
istat = nf90_put_att(ncid, varid(7), &
   'standard_name', &
   'clausius_clapeyron_gradient')
istat = nf90_put_att(ncid, varid(7), &
   'long_name', &
   'Clausius-Clapeyron gradient')

istat = nf90_def_var(ncid, 'DELTA_TM_SW', NF90_DOUBLE, varid(8))

istat = nf90_put_att(ncid, varid(8), 'units', 'degC')
istat = nf90_put_att(ncid, varid(8), &
   'standard_name', &
   'sea_water_melting_point_depression_due_to_salinity')
istat = nf90_put_att(ncid, varid(8), &
   'long_name', &
   'Melting point depression of sea water due to salinity')

istat = nf90_def_var(ncid, 'OMEGA_MAX', NF90_DOUBLE, varid(9))

istat = nf90_put_att(ncid, varid(9), 'units', '-')
istat = nf90_put_att(ncid, varid(9), &
   'standard_name', &
   'water_content_threshold_value')
istat = nf90_put_att(ncid, varid(9), &
   'long_name', &
   'Threshold value for the water content')

istat = nf90_def_var(ncid, 'H_R', NF90_DOUBLE, varid(10))

istat = nf90_put_att(ncid, varid(10), 'units', 'm')
istat = nf90_put_att(ncid, varid(10), &
   'standard_name', &
   'thermal_lithosphere_layer_thickness')
istat = nf90_put_att(ncid, varid(10), &
   'long_name', &
   'Thickness of the thermal lithosphere layer')

istat = nf90_def_var(ncid, 'RHO_C_R', NF90_DOUBLE, varid(11))

istat = nf90_put_att(ncid, varid(11), 'units', 'J m-3 K-1')
istat = nf90_put_att(ncid, varid(11), &
   'standard_name', &
   'lithosphere_density_times_specific_heat')
istat = nf90_put_att(ncid, varid(11), &
   'long_name', &
   'Density times specific heat of the lithosphere')

istat = nf90_def_var(ncid, 'KAPPA_R', NF90_DOUBLE, varid(12))

istat = nf90_put_att(ncid, varid(12), 'units', 'W m-1 K-1')
istat = nf90_put_att(ncid, varid(12), &
   'standard_name', &
   'lithosphere_heat_conductivity')
istat = nf90_put_att(ncid, varid(12), &
   'long_name', &
   'Heat conductivity of the lithosphere')

istat = nf90_def_var(ncid, 'RHO_A', NF90_DOUBLE, varid(13))

istat = nf90_put_att(ncid, varid(13), 'units', 'kg m-3')
istat = nf90_put_att(ncid, varid(13), 'standard_name', 'asthenosphere_density')
istat = nf90_put_att(ncid, varid(13), 'long_name', 'Density of the asthenosphere')

istat = nf90_def_var(ncid, 'R_T', NF90_DOUBLE, varid(14))

istat = nf90_put_att(ncid, varid(14), 'units', '-')
istat = nf90_put_att(ncid, varid(14), &
   'standard_name', &
   'ice_rate_factor_water_content_dependence_coefficient')
istat = nf90_put_att(ncid, varid(14), &
   'long_name', &
   'Coefficient of the water-content dependence ' &
   // 'of the rate factor of ice')

istat = nf90_def_var(ncid, 'R', NF90_DOUBLE, varid(21))

istat = nf90_put_att(ncid, varid(21), 'units', 'm')
istat = nf90_put_att(ncid, varid(21), &
   'standard_name', &
   'earth_radius')
istat = nf90_put_att(ncid, varid(21), &
   'long_name', &
   'Radius of the Earth')

istat = nf90_def_var(ncid, 'A', NF90_DOUBLE, varid(22))

istat = nf90_put_att(ncid, varid(22), 'units', 'm')
istat = nf90_put_att(ncid, varid(22), &
   'standard_name', &
   'earth_semi_major_axis')
istat = nf90_put_att(ncid, varid(22), &
   'long_name', &
   'Semi-major axis of the Earth')

istat = nf90_def_var(ncid, 'F_INV', NF90_DOUBLE, varid(23))

istat = nf90_put_att(ncid, varid(23), 'units', '-')
istat = nf90_put_att(ncid, varid(23), &
   'standard_name', &
   'earth_inverse_flattening')
istat = nf90_put_att(ncid, varid(23), &
   'long_name', &
   'Inverse flattening of the Earth')

istat = nf90_def_var(ncid, 'LATD0', NF90_DOUBLE, varid(24))

istat = nf90_put_att(ncid, varid(24), 'units', 'degrees_N')
istat = nf90_put_att(ncid, varid(24), 'standard_name', 'standard_parallel')
istat = nf90_put_att(ncid, varid(24), 'long_name', 'Standard parallel')

istat = nf90_def_var(ncid, 'LOND0', NF90_DOUBLE, varid(25))

istat = nf90_put_att(ncid, varid(25), 'units', 'degrees_E')
istat = nf90_put_att(ncid, varid(25), 'standard_name', 'central_meridian')
istat = nf90_put_att(ncid, varid(25), 'long_name', 'Central meridian')

istat = nf90_def_var(ncid, 'RF', NF90_DOUBLE, ncd, varid(51))

istat = nf90_put_att(ncid, varid(51), 'units', 's-1 Pa-3')
istat = nf90_put_att(ncid, varid(51), 'standard_name', 'ice_rate_factor')
istat = nf90_put_att(ncid, varid(51), 'long_name', 'Rate factor of ice')

istat = nf90_def_var(ncid, 'KAPPA', NF90_DOUBLE, ncd, varid(52))

istat = nf90_put_att(ncid, varid(52), 'units', 'W m-1 K-1')
istat = nf90_put_att(ncid, varid(52), 'standard_name', 'ice_heat_conductivity')
istat = nf90_put_att(ncid, varid(52), 'long_name', 'Heat conductivity of ice')

istat = nf90_def_var(ncid, 'C', NF90_DOUBLE, ncd, varid(53))

istat = nf90_put_att(ncid, varid(53), 'units', 'J kg-1 K-1')
istat = nf90_put_att(ncid, varid(53), 'standard_name', 'ice_specific_heat')
istat = nf90_put_att(ncid, varid(53), 'long_name', 'Specific heat of ice')

istat = nf90_enddef(ncid)

!  ------ Write variables on file

istat = nf90_put_var(ncid, varid( 1), RHO)
istat = nf90_put_var(ncid, varid( 2), RHO_W)
istat = nf90_put_var(ncid, varid( 3), RHO_SW)
istat = nf90_put_var(ncid, varid( 4), L)
istat = nf90_put_var(ncid, varid( 5), G)
istat = nf90_put_var(ncid, varid( 6), NUE)
istat = nf90_put_var(ncid, varid( 7), BETA)
istat = nf90_put_var(ncid, varid( 8), DELTA_TM_SW)
istat = nf90_put_var(ncid, varid( 9), OMEGA_MAX)
istat = nf90_put_var(ncid, varid(10), H_R)
istat = nf90_put_var(ncid, varid(11), RHO_C_R)
istat = nf90_put_var(ncid, varid(12), KAPPA_R)
istat = nf90_put_var(ncid, varid(13), RHO_A)
istat = nf90_put_var(ncid, varid(14), R_T)

istat = nf90_put_var(ncid, varid(21), R)
istat = nf90_put_var(ncid, varid(22), A)
istat = nf90_put_var(ncid, varid(23), F_INV)
istat = nf90_put_var(ncid, varid(24), LATD0)
istat = nf90_put_var(ncid, varid(25), LOND0)

istat = nf90_put_var(ncid, varid(51), RF)
istat = nf90_put_var(ncid, varid(52), KAPPA)
istat = nf90_put_var(ncid, varid(53), C)

!  ------ Close file

istat = nf90_sync(ncid)
istat = nf90_close(ncid)

!-------- End of program --------

write(6,'(/a/)') 'Done.'

end program make_phys_para_nc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                        End of make_phys_para_nc.F90
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
