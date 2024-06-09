!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ g i a _ m
!
!! Computation of the glacial isostatic adjustment of the lithosphere surface.
!!
!!##### Authors
!!
!! Ralf Greve, Sascha Knell
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
!> Computation of the glacial isostatic adjustment of the lithosphere surface.
!-------------------------------------------------------------------------------
module calc_gia_m

  use sico_types_m
  use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

  use error_m

  implicit none

  private
  public :: calc_gia

contains

!-------------------------------------------------------------------------------
!> Main subroutine of calc_gia_m:
!! Computation of the glacial isostatic adjustment of the lithosphere surface.
!-------------------------------------------------------------------------------
subroutine calc_gia(time, dtime, dxi, deta, itercount, iter_wss)

#if defined(ALLOW_TAPENADE) /* Tapenade */
#if (THK_EVOL==2)
  use ctrl_m, only: myfloor, myceiling
#endif
#endif /* Tapenade */

implicit none

integer(i4b), intent(in) :: itercount, iter_wss
real(dp),     intent(in) :: time
real(dp),     intent(in) :: dtime, dxi, deta

integer(i4b) :: i, j
integer(i4b) :: n1, n2
real(dp) :: time_in_years
real(dp) :: time1, time2
real(dp) :: tldt_inv(0:JMAX,0:IMAX)
real(dp) :: time_ratio_1(0:JMAX,0:IMAX), time_ratio_2(0:JMAX,0:IMAX)
real(dp) :: load_ice_water(0:JMAX,0:IMAX)
real(dp) :: dtime_inv
real(dp) :: rho_g, rhosw_g, rhoa_g_inv

!-------- Term abbreviations --------

time_in_years = time*sec2year

do i=0, IMAX
do j=0, JMAX
   tldt_inv(j,i) = 1.0_dp/(time_lag_asth(j,i)+dtime)
   time_ratio_1(j,i) = tldt_inv(j,i) * time_lag_asth(j,i)
   time_ratio_2(j,i) = tldt_inv(j,i) * dtime
end do
end do

rho_g      = RHO*G
rhosw_g    = RHO_SW*G
rhoa_g_inv = 1.0_dp/(RHO_A*G)

dtime_inv = 1.0_dp/dtime

!-------- Load due to ice and sea water --------

#if (REBOUND==0)

load_ice_water = 0.0_dp   ! not needed, thus not computed

#elif (REBOUND==1 || REBOUND==2)

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i) <= 1) then   ! grounded ice or ice-free land

      load_ice_water(j,i) = rho_g * H(j,i)

   else   ! (mask(j,i) >= 2, floating ice or ocean)

      load_ice_water(j,i) = rhosw_g * z_sl(j,i)
                   ! Water load relative to the present sea-level stand (0 m)
                   ! -> can be positive or negative

   end if

end do
end do

#endif

!-------- Steady-state displacement of the lithosphere
!                              (wss, positive downward) --------

#if (REBOUND==0)

wss = 0.0_dp

#elif (REBOUND==1)

wss = FRAC_LLRA * ( load_ice_water * rhoa_g_inv )

#elif (REBOUND==2)

if ((mod(itercount, iter_wss)==0).or.(itercount==1)) then
   write(6, fmt='(10x,a)') 'Computation of wss (EL)'
   call calc_el(load_ice_water, dxi, deta)
end if

#endif

!-------- Computation of zl_new and its time derivative --------

#if (REBOUND==0)

do i=0, IMAX
do j=0, JMAX
   zl_new(j,i)   = zl(j,i)
   dzl_dtau(j,i) = (zl_new(j,i)-zl(j,i))*dtime_inv
end do
end do

#elif (REBOUND==1 || REBOUND==2)

do i=0, IMAX
do j=0, JMAX
   zl_new(j,i)   = time_ratio_1(j,i)*zl(j,i) &
                      + time_ratio_2(j,i)*(zl0(j,i)-wss(j,i))
   dzl_dtau(j,i) = (zl_new(j,i)-zl(j,i))*dtime_inv
end do
end do

#endif

!  ------ Adjustment due to prescribed target topography

#if (THK_EVOL==2)

if (time_in_years < real(target_topo_tau0_time_min,dp)) then

   target_topo_tau = target_topo_tau0(0)

else if (time_in_years < real(target_topo_tau0_time_max,dp)) then

#if !defined(ALLOW_TAPENADE) /* Normal */
   n1 = floor((time_in_years &
           -real(target_topo_tau0_time_min,dp)) &
              /real(target_topo_tau0_time_stp,dp))
#else /* Tapenade */
   call myfloor((time_in_years &
           -real(target_topo_tau0_time_min,dp)) &
              /real(target_topo_tau0_time_stp,dp),n1)
#endif /* Normal vs. Tapenade */

   n1 = max(n1, 0)

#if !defined(ALLOW_TAPENADE) /* Normal */
   n2 = ceiling((time_in_years &
           -real(target_topo_tau0_time_min,dp)) &
              /real(target_topo_tau0_time_stp,dp))
#else /* Tapenade */
   call myceiling(((time_in_years &
            -real(target_topo_tau0_time_min,dp)) &
               /real(target_topo_tau0_time_stp,dp)),n2)
#endif /* Normal vs. Tapenade */

   n2 = min(n2, ndata_target_topo_tau0)

   if (n1 == n2) then

      target_topo_tau = target_topo_tau0(n1)

   else

      time1 = (target_topo_tau0_time_min + n1*target_topo_tau0_time_stp) &
              *year2sec
      time2 = (target_topo_tau0_time_min + n2*target_topo_tau0_time_stp) &
              *year2sec

      target_topo_tau = target_topo_tau0(n1) &
                +(target_topo_tau0(n2)-target_topo_tau0(n1)) &
                *(time-time1)/(time2-time1)
                 ! linear interpolation of the relaxation-time data

   end if

else

   target_topo_tau = target_topo_tau0(ndata_target_topo_tau0)

end if

if (target_topo_tau*sec2year > no_value_pos_1) then
           ! relaxation time target_topo_tau interpreted as infinity

#if (!defined(ALLOW_GRDCHK) && !defined(ALLOW_TAPENADE)) /* Normal */
   target_topo_tau = huge(1.0_dp)
#endif

   !%% zl_new = zl_new

else if (target_topo_tau*sec2year < epsi) then
           ! relaxation time target_topo_tau interpreted as zero

   target_topo_tau = 0.0_dp
   zl_new          = zl_target

else

   zl_new =   ( target_topo_tau*zl_new + dtime*zl_target ) &
            / ( target_topo_tau        + dtime )

end if

dzl_dtau = (zl_new-zl)*dtime_inv

#elif (THK_EVOL==3)

target_topo_tau = target_topo_tau_0

if (target_topo_tau*sec2year > no_value_pos_1) then
           ! relaxation time target_topo_tau interpreted as infinity

#if (!defined(ALLOW_GRDCHK) && !defined(ALLOW_TAPENADE)) /* Normal */
   target_topo_tau = huge(1.0_dp)
#endif

   !%% zl_new = zl_new

else if (target_topo_tau*sec2year < epsi) then
           ! relaxation time target_topo_tau interpreted as zero

   target_topo_tau = 0.0_dp
   zl_new          = zl_target

else

   zl_new =   ( target_topo_tau*zl_new + dtime*zl_target ) &
            / ( target_topo_tau        + dtime )

end if

dzl_dtau = (zl_new-zl)*dtime_inv

#endif

!======== Generation of an isostatically relaxed
!         lithosphere surface topography. Not to be used regularly! ========

#if (defined(EXEC_MAKE_ZL0))

call make_zl0()

errormsg = ' >>> calc_gia: Routine make_zl0 successfully completed,' &
         //         end_of_line &
         //'        isostatically relaxed lithosphere topography' &
         //         end_of_line &
         //'        written on file (in output directory of this run).' &
         //         end_of_line &
         //'        Execution of SICOPOLIS stopped.'
call error(errormsg)   ! actually not an error,
                       ! just a regular stop with an info message

#endif

end subroutine calc_gia

!-------------------------------------------------------------------------------
!> Computation of the isostatic steady-state displacement of the lithosphere
!! for the elastic-lithosphere (EL) model.
!-------------------------------------------------------------------------------
subroutine calc_el(load_ice_water, dxi, deta)

!$ use omp_lib

#if defined(ALLOW_TAPENADE) /* Tapenade */
  use ctrl_m, only: myfloor
#endif /* Tapenade */

implicit none

real(dp), intent(in), dimension(0:JMAX,0:IMAX) :: load_ice_water
real(dp), intent(in)                           :: dxi, deta

integer(i4b) :: i, j, ij, ir, jr, il, jl, n
integer(i4b) :: ir_max, jr_max, min_imax_jmax
integer(i4b) :: il_begin, il_end, jl_begin, jl_end
real(dp)                              :: rhoa_g_inv
real(dp)                              :: dxi_inv, deta_inv
real(dp)                              :: kei_r_incr_inv
real(dp), dimension(0:JMAX,0:IMAX)    :: l_r, l_r_inv, fac_wss
real(dp)                              :: ra_max

#if !defined(ALLOW_TAPENADE) /* Normal */
real(dp), allocatable, dimension(:,:) :: f_0
#else /* Tapenade */
real(dp), dimension(-JMAX:2*JMAX,-IMAX:2*IMAX) :: f_0
#endif /* Normal vs. Tapenade */

real(dp), parameter :: r_infl = 8.0_dp
                       ! Radius of non-locality influence of the elastic
                       ! lithosphere plate (in units of l_r)

!-------- Initialisations --------

rhoa_g_inv = 1.0_dp/(RHO_A*G)

dxi_inv   = 1.0_dp/dxi
deta_inv  = 1.0_dp/deta

kei_r_incr_inv = 1.0_dp/kei_r_incr

do i=0, IMAX
do j=0, JMAX
   l_r(j,i)     = (flex_rig_lith(j,i)*rhoa_g_inv)**0.25_dp
   l_r_inv(j,i) = 1.0_dp/l_r(j,i)
   fac_wss(j,i) = l_r(j,i)*l_r(j,i)/(2.0_dp*pi*flex_rig_lith(j,i))
end do
end do

ra_max  = r_infl*maxval(l_r)   ! Radius of non-locality influence (in m)

#if !defined(ALLOW_TAPENADE) /* Normal */

ir_max = floor(ra_max*dxi_inv)
jr_max = floor(ra_max*deta_inv)
         ! Radius of non-locality influence (in grid points)

#else /* Tapenade */

call myfloor(ra_max*dxi_inv, ir_max)
call myfloor(ra_max*deta_inv, jr_max)
         ! Radius of non-locality influence (in grid points)

if (ir_max.gt.IMAX) then
  write(*,*) "ERROR: unhandled case ir_max greater than IMAX for f_0 dimensions"
end if
if (jr_max.gt.JMAX) then
  write(*,*) "ERROR: unhandled case jr_max greater than JMAX for f_0 dimensions"
end if

#endif /* Normal vs. Tapenade */

min_imax_jmax = min(IMAX, JMAX)
ir_max        = min(ir_max, min_imax_jmax)
jr_max        = min(jr_max, min_imax_jmax)
                ! ir_max and jr_max constrained by size of domain

il_begin = -ir_max
il_end   = ir_max+IMAX
jl_begin = -jr_max
jl_end   = jr_max+JMAX

!-------- Ice/water load --------

#if !defined(ALLOW_TAPENADE) /* Normal */
allocate(f_0(jl_begin:jl_end, il_begin:il_end))
#endif /* Normal */

do il=il_begin, il_end
do jl=jl_begin, jl_end

   i = min(max(il, 0), IMAX)
   j = min(max(jl, 0), JMAX)

   if ((mask(j,i)==0).or.(mask(j,i)>=2)) then
      f_0(jl,il) = load_ice_water(j,i) * cell_area(j,i)
   else   ! (mask(j,i)==1)
      f_0(jl,il) = 0.0_dp
   end if

end do
end do

!-------- Steady-state displacement of the lithosphere --------

!$omp parallel do private(ij,i,j,ir,jr,n)
do ij=1, (IMAX+1)*(JMAX+1)

   i = n2i(ij)
   j = n2j(ij)

   wss(j,i) = 0.0_dp

   do ir=-ir_max, ir_max
   do jr=-jr_max, jr_max

#if (GRID==0)
      n = nint( dist_dxdy(jr,ir)*l_r_inv(j,i)*kei_r_incr_inv )
#elif (GRID==1)
      n = nint( sq_g11_g(j,i)*dist_dxdy(jr,ir)*l_r_inv(j,i)*kei_r_incr_inv )
              ! The correction with g11 only is OK because g11=g22 for this case
#elif (GRID==2)
      n = nint( dist_dxdy(jr,ir)*l_r_inv(j,i)*kei_r_incr_inv )
              ! The distortion due to g11 and g22 is already accounted for in
              ! the variable dist_dxdy(jr,ir)
#endif

      wss(j,i) = wss(j,i) - fac_wss(j,i)*f_0(jr+j,ir+i)*kei(n)

   end do
   end do

end do
!$omp end parallel do

#if !defined(ALLOW_TAPENADE) /* Normal */
deallocate (f_0)
#endif /* Normal */

end subroutine calc_el

!-------------------------------------------------------------------------------
!> Generation of an isostatically relaxed lithosphere surface topography
!! for either the rigid lithosphere, the local lithosphere or the elastic
!! lithosphere model (depending on the setting of the parameter REBOUND).
!! This routine is not to be used regularly, and it is only executed if the
!! parameter EXEC_MAKE_ZL0 is defined in the header file.
!-------------------------------------------------------------------------------
subroutine make_zl0()

#if defined(ALLOW_TAPENADE) /* Tapenade */
  use ctrl_m, only: myceiling 
#endif /* Tapenade */

  use netcdf
  use nc_check_m

implicit none

integer(i4b)                       :: i, j, m, n, i_f, j_f, n_filter
integer(i4b)                       :: ios
real(dp), dimension(0:JMAX,0:IMAX) :: zl0_raw, zl0_smoothed
real(dp)                           :: dx
real(dp)                           :: filter_width, sigma_filter
real(dp)                           :: dist, weigh, sum_weigh
character(len=  8)                 :: ch_resolution
character(len= 64)                 :: ch_model
character(len=256)                 :: filename, filename_with_path

real(sp), dimension(0:IMAX,0:JMAX) :: zl0_conv

integer(i4b) :: ncid, ncv
integer(i4b) :: ncd, nc1d, nc2d(2)
integer(i4b) :: nc1cor_i(1), nc1cor_j(1), nc2cor_ij(2)
integer(i4b) :: nc1cnt_i(1), nc1cnt_j(1), nc2cnt_ij(2)
character(len=64), parameter :: thisroutine = 'make_zl0'

!-------- Checking of the grid --------

#if (GRID==0 || GRID==1)

dx = real(DX, dp)   ! horizontal grid spacing (in km)

#else

dx = 0.0_dp   ! dummy value

errormsg = ' >>> make_zl0: Routine works only for GRID 0 or 1!'
call error(errormsg)

#endif

!-------- Computation of the raw (unsmoothed) topography --------

#if (REBOUND==0)

zl0_raw = zl   ! rigid lithosphere

#elif (REBOUND==1 || REBOUND==2)

zl0_raw = zl + wss   ! local or elastic lithosphere

#else

errormsg = ' >>> make_zl0: REBOUND must be either 0, 1 or 2!'
call error(errormsg)

#endif

!-------- Smoothing of the topography (Gaussian filter) --------

#if (defined(ZL0_FILTER_WIDTH))
  filter_width = real(ZL0_FILTER_WIDTH, dp)
                      ! filter width (half span of filtered area), in km
#else
  filter_width = dx   ! default value: one grid spacing, in km
#endif

sigma_filter = filter_width/dx   ! half span of filtered area,
                                 ! in grid points

#if !defined(ALLOW_TAPENADE) /* Normal */
n_filter = ceiling(2.0_dp*sigma_filter)
#else /* Tapenade */
call myceiling(2.0_dp*sigma_filter, n_filter)
#endif /* Normal vs. Tapenade */

n_filter = max(n_filter, 5)

if (sigma_filter > eps_sp_dp) then

   zl0_smoothed = 0.0_dp

   do i=0, IMAX
   do j=0, JMAX

      sum_weigh = 0.0_dp

      do m=-n_filter, n_filter
      do n=-n_filter, n_filter

         i_f = i+m
         j_f = j+n

         if (i_f <    0) i_f =    0
         if (i_f > IMAX) i_f = IMAX

         if (j_f <    0) j_f =    0
         if (j_f > JMAX) j_f = JMAX

         dist      = sqrt(real(m,dp)**2+real(n,dp)**2)
         weigh     = exp(-(dist/sigma_filter)**2)
         sum_weigh = sum_weigh + weigh

         zl0_smoothed(j,i) = zl0_smoothed(j,i) + weigh*zl0_raw(j_f,i_f)

      end do
      end do

      zl0_smoothed(j,i) = zl0_smoothed(j,i)/sum_weigh

   end do
   end do

else

   zl0_smoothed = zl0_raw   ! no smoothing

end if

!-------- Writing on ASCII file --------

if ( abs(dx-nint(dx)) < eps ) then
   write(ch_resolution, fmt='(i8)') nint(dx)
else
   write(ch_resolution, fmt='(f8.2)') dx
end if

#if !defined(ALLOW_TAPENADE) /* Normal */
ch_resolution = adjustl(ch_resolution)
#endif /* Normal */

#if (REBOUND==0)
ch_model = 'rigid_lithosphere'
#elif (REBOUND==1)
ch_model = 'local_lithosphere'
#elif (REBOUND==2)
ch_model = 'elastic_lithosphere'
#endif

filename           = trim(ch_domain_short)//'_'//trim(ch_resolution)
filename           = trim(filename)//'_zl0_'//trim(ch_model)//'.dat'
filename_with_path = trim(OUT_PATH)//'/'//trim(filename)

open(23, iostat=ios, file=trim(filename_with_path), recl=rcl1, status='replace')

if (ios /= 0) then
   errormsg = ' >>> make_zl0: Error when opening the new zl0 ASCII file!'
   call error(errormsg)
end if

write(23,'(a)') &
'%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
write(23,'(a)') '% '//trim(ch_domain_long)//':'
write(23,'(a)') &
'% Relaxed lithosphere surface topography without ice load, in m AMSL.'
write(23,'(a)') &
'% Horizontal resolution '//trim(ch_resolution)//' km.'
write(23,'(a,i3,a,i3,a,i3,a,i3,a)') &
'% ', JMAX+1, ' records [j = ', &
JMAX, ' (-1) 0] with ', IMAX+1, &
' values [i = 0 (1) ', IMAX, '] in each record.'
write(23,'(a)') &
'%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

do j=JMAX, 0, -1
   do i=0, IMAX-1
      zl0_smoothed(j,i) = max(zl0_smoothed(j,i), no_value_neg_2)
      write(23,'(f8.1)', advance='no') zl0_smoothed(j,i)
   end do
   i=IMAX
      zl0_smoothed(j,i) = max(zl0_smoothed(j,i), no_value_neg_2)
      write(23,'(f8.1)') zl0_smoothed(j,i)
end do

close(23, status='keep')

!-------- Writing on NetCDF file --------

!  ------ Open NetCDF file

filename           = trim(ch_domain_short)//'_'//trim(ch_resolution)
filename           = trim(filename)//'_zl0_'//trim(ch_model)//'.nc'
filename_with_path = trim(OUT_PATH)//'/'//trim(filename)

ios = nf90_create(trim(filename_with_path), NF90_NOCLOBBER, ncid)

if (ios /= nf90_noerr) then
   errormsg = ' >>> make_zl0: Error when opening the new zl0 NetCDF file!'
   call error(errormsg)
end if

!  ------ Definition of the dimensions

call check( nf90_def_dim(ncid, 'x', IMAX+1, ncd), thisroutine )
call check( nf90_def_dim(ncid, 'y', JMAX+1, ncd), thisroutine )

!  ------ Definition of the variables

!    ---- x (= xi)

call check( nf90_inq_dimid(ncid, 'x', nc1d), thisroutine )
call check( nf90_def_var(ncid, 'x', NF90_DOUBLE, nc1d, ncv), thisroutine )
call check( nf90_put_att(ncid, ncv, 'units', 'm'), thisroutine )
call check( nf90_put_att(ncid, ncv, &
                         'standard_name', 'projection_x_coordinate'), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, &
                         'long_name', 'x-coordinate of the grid point i'), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'axis', 'x'), thisroutine )

!    ---- y (= eta)

call check( nf90_inq_dimid(ncid, 'y', nc1d), thisroutine )
call check( nf90_def_var(ncid, 'y', NF90_DOUBLE, nc1d, ncv), thisroutine )
call check( nf90_put_att(ncid, ncv, 'units', 'm'), thisroutine )
call check( nf90_put_att(ncid, ncv, &
                         'standard_name', 'projection_y_coordinate'), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, &
                         'long_name', 'y-coordinate of the grid point j'), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'axis', 'y'), thisroutine )

!    ---- zl0

call check( nf90_inq_dimid(ncid, 'x', nc2d(1)), thisroutine )
call check( nf90_inq_dimid(ncid, 'y', nc2d(2)), thisroutine )
call check( nf90_def_var(ncid, 'zl0', NF90_FLOAT, nc2d, ncv), thisroutine )
call check( nf90_put_att(ncid, ncv, 'units', 'm'), thisroutine )
call check( nf90_put_att(ncid, ncv, &
               'standard_name', &
               'isostatically_relaxed_bedrock_altitude'), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, &
               'long_name', &
               'Topography of the isostatically relaxed lithosphere surface'), &
            thisroutine )

!    ---- End of definitions

call check( nf90_enddef(ncid), thisroutine )

!  ------ Writing variables on file

nc1cor_i = [ 1 ]
nc1cnt_i = [ IMAX+1 ]

nc1cor_j = [ 1 ]
nc1cnt_j = [ JMAX+1 ]

nc2cor_ij = [ 1, 1 ]
nc2cnt_ij = [ IMAX+1, JMAX+1 ]

call check( nf90_inq_varid(ncid, 'x', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, xi, &
                         start=nc1cor_i, count=nc1cnt_i), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'y', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, eta, &
                         start=nc1cor_j, count=nc1cnt_j), &
            thisroutine )

do i=0, IMAX
do j=0, JMAX
   zl0_conv(i,j) = real(zl0_smoothed(j,i),sp)
end do
end do

call check( nf90_inq_varid(ncid, 'zl0', ncv), &
            thisroutine )
call check( nf90_put_var(ncid, ncv, zl0_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

!  ------ Closing of NetCDF file

call check( nf90_sync(ncid), thisroutine )

call check( nf90_close(ncid), thisroutine )

end subroutine make_zl0

!-------------------------------------------------------------------------------

end module calc_gia_m
!
