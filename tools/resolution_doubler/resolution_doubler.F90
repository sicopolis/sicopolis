!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Program   : r e s o l u t i o n _ d o u b l e r . F 9 0
!
!! Doubling the horizontal resolution of a NetCDF time-slice output file
!! produced by SICOPOLIS.
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

!-------- Inclusion of specification header --------

#include RUN_SPECS_HEADER

!-------------------------------------------------------------------------------
!> Declarations of kind types.
!-------------------------------------------------------------------------------
module resolution_doubler_types

integer, parameter :: i4b = selected_int_kind(9)   ! 4-byte integers
integer, parameter :: sp  = kind(1.0)              ! single-precision reals
integer, parameter :: dp  = kind(1.0d0)            ! double-precision reals

end module resolution_doubler_types

!-------------------------------------------------------------------------------
!> Declarations of global variables.
!-------------------------------------------------------------------------------
module resolution_doubler_vars

use resolution_doubler_types

integer(i4b) :: mapping_erg
integer(i4b), dimension(0:IMAX,0:JMAX) :: mask_erg, mask_old_erg, &
                                          n_cts_erg, &
                                          flag_shelfy_stream_x_erg, &
                                          flag_shelfy_stream_y_erg, &
                                          flag_shelfy_stream_erg, &
                                          flag_grounding_line_1_erg, &
                                          flag_grounding_line_2_erg, &
                                          flag_calving_front_1_erg, &
                                          flag_calving_front_2_erg, &
                                          flag_grounded_front_a_1_erg, &
                                          flag_grounded_front_a_2_erg, &
                                          flag_grounded_front_b_1_erg, &
                                          flag_grounded_front_b_2_erg
integer(i4b), dimension(0:IMAX,0:JMAX) :: kc_cts_erg
real(dp) :: mapping_semi_major_axis_erg, &
            mapping_inv_flattening_erg, &
            mapping_radius_of_sphere_erg, &
            mapping_latitude_origin_erg, &
            mapping_standard_parallel_erg, &
            mapping_reference_longitude_erg, &
            mapping_false_E_erg, &
            mapping_false_N_erg
real(dp) :: year2sec_erg, time_erg, &
            delta_ts_erg, glac_index_erg, z_sl_mean_erg, &
            V_tot_erg, V_af_erg, A_grounded_erg, A_floating_erg, &
            xi_erg(0:IMAX), eta_erg(0:JMAX), &
            sigma_level_c_erg(0:KCMAX), sigma_level_t_erg(0:KTMAX), &
            sigma_level_r_erg(0:KRMAX)
real(sp) :: H_R_erg
real(sp), dimension(0:IMAX,0:JMAX) :: lambda_erg, phi_erg, &
            lon_erg, lat_erg, &
            temp_s_erg, prec_erg, &
            snowfall_erg, rainfall_erg, pdd_erg, & 
            as_perp_erg, as_perp_apl_erg, smb_corr_erg, &
            z_sl_erg, &
            q_geo_erg, &
            zs_erg, zm_erg, zb_erg, zl_erg, zl0_erg, wss_erg, &
            H_cold_erg, H_temp_erg, H_erg, &
            Q_bm_erg, Q_tld_erg, calving_erg, &
            am_perp_erg, &
            qx_erg, qy_erg, &
            vx_m_sia_erg, vy_m_sia_erg, vx_m_ssa_erg, vy_m_ssa_erg, &
            dzs_dtau_erg, dzm_dtau_erg, dzb_dtau_erg, dzl_dtau_erg, &
            dH_c_dtau_erg, dH_t_dtau_erg, dH_dtau_erg, &
            vx_b_g_erg, vy_b_g_erg, vz_b_erg, vh_b_erg, &
            vx_s_g_erg, vy_s_g_erg, vz_s_erg, vh_s_erg, &
            vx_m_g_erg, vy_m_g_erg,           vh_m_erg, &
            temp_b_erg, temph_b_erg, &
            tau_dr_erg, tau_b_erg, &
            p_b_w_erg, q_w_erg, q_w_x_erg, q_w_y_erg, H_w_erg, &
            q_gl_g_erg, &
            ratio_sl_sia_x_erg, ratio_sl_sia_y_erg, &
            vis_ave_g_erg, vis_int_g_erg
real(sp), dimension(0:IMAX,0:JMAX) :: r_kc_cts_erg
real(sp), dimension(0:IMAX,0:JMAX,0:KCMAX) :: vx_c_erg, vy_c_erg, vz_c_erg, &
                                              temp_c_erg, age_c_erg, &
                                              enth_c_erg, omega_c_erg, &
                                              enh_c_erg, strain_heating_c_erg
real(sp), dimension(0:IMAX,0:JMAX,0:KTMAX) :: vx_t_erg, vy_t_erg, vz_t_erg, &
                                              omega_t_erg, age_t_erg, &
                                              enth_t_erg, &
                                              enh_t_erg, strain_heating_t_erg
real(sp), dimension(0:IMAX,0:JMAX,0:KRMAX) :: temp_r_erg
character(len=64) :: mapping_grid_mapping_name_erg, mapping_ellipsoid_erg

logical :: flag_wss

#if (DISC>0)   /* Ice discharge parameterisation */
integer(i4b), dimension(0:IMAX,0:JMAX) :: mask_mar_erg
real(sp),     dimension(0:IMAX,0:JMAX) :: dis_perp_erg, &
                                          cst_dist_erg, cos_grad_tc_erg
#endif

integer(i4b) :: mapping_dbl
integer(i4b), dimension(0:2*IMAX,0:2*JMAX) :: mask_dbl, mask_old_dbl, &
                                              n_cts_dbl, &
                                              flag_shelfy_stream_x_dbl, &
                                              flag_shelfy_stream_y_dbl, &
                                              flag_shelfy_stream_dbl, &
                                              flag_grounding_line_1_dbl, &
                                              flag_grounding_line_2_dbl, &
                                              flag_calving_front_1_dbl, &
                                              flag_calving_front_2_dbl, &
                                              flag_grounded_front_a_1_dbl, &
                                              flag_grounded_front_a_2_dbl, &
                                              flag_grounded_front_b_1_dbl, &
                                              flag_grounded_front_b_2_dbl
integer(i4b), dimension(0:2*IMAX,0:2*JMAX) :: kc_cts_dbl
real(dp) :: mapping_semi_major_axis_dbl, &
            mapping_inv_flattening_dbl, &
            mapping_radius_of_sphere_dbl, &
            mapping_latitude_origin_dbl, &
            mapping_standard_parallel_dbl, &
            mapping_reference_longitude_dbl, &
            mapping_false_E_dbl, &
            mapping_false_N_dbl
real(dp) :: year2sec_dbl, time_dbl, &
            delta_ts_dbl, glac_index_dbl, z_sl_mean_dbl, &
            V_tot_dbl, V_af_dbl, A_grounded_dbl, A_floating_dbl, &
            xi_dbl(0:2*IMAX), eta_dbl(0:2*JMAX), &
            sigma_level_c_dbl(0:KCMAX), sigma_level_t_dbl(0:KTMAX), &
            sigma_level_r_dbl(0:KRMAX)
real(sp) :: H_R_dbl
real(sp), dimension(0:2*IMAX,0:2*JMAX) :: lambda_dbl, phi_dbl, &
            lon_dbl, lat_dbl, &
            temp_s_dbl, prec_dbl, &
            snowfall_dbl, rainfall_dbl, pdd_dbl, & 
            as_perp_dbl, as_perp_apl_dbl, smb_corr_dbl, &
            z_sl_dbl, &
            q_geo_dbl, &
            zs_dbl, zm_dbl, zb_dbl, zl_dbl, zl0_dbl, wss_dbl, &
            H_cold_dbl, H_temp_dbl, H_dbl, &
            Q_bm_dbl, Q_tld_dbl, calving_dbl, &
            am_perp_dbl, &
            qx_dbl, qy_dbl, &
            vx_m_sia_dbl, vy_m_sia_dbl, vx_m_ssa_dbl, vy_m_ssa_dbl, &
            dzs_dtau_dbl, dzm_dtau_dbl, dzb_dtau_dbl, dzl_dtau_dbl, &
            dH_c_dtau_dbl, dH_t_dtau_dbl, dH_dtau_dbl, &
            vx_b_g_dbl, vy_b_g_dbl, vz_b_dbl, vh_b_dbl, &
            vx_s_g_dbl, vy_s_g_dbl, vz_s_dbl, vh_s_dbl, &
            vx_m_g_dbl, vy_m_g_dbl,           vh_m_dbl, &
            temp_b_dbl, temph_b_dbl, &
            tau_dr_dbl, tau_b_dbl, &
            p_b_w_dbl, q_w_dbl, q_w_x_dbl, q_w_y_dbl, H_w_dbl, &
            q_gl_g_dbl, &
            ratio_sl_sia_x_dbl, ratio_sl_sia_y_dbl, &
            vis_ave_g_dbl, vis_int_g_dbl
real(sp), dimension(0:2*IMAX,0:2*JMAX) :: r_kc_cts_dbl
real(sp), dimension(0:2*IMAX,0:2*JMAX,0:KCMAX) :: vx_c_dbl, vy_c_dbl, &
                                                  vz_c_dbl, &
                                                  temp_c_dbl, age_c_dbl, &
                                                  enth_c_dbl, omega_c_dbl, &
                                                  enh_c_dbl, &
                                                  strain_heating_c_dbl
real(sp), dimension(0:2*IMAX,0:2*JMAX,0:KTMAX) :: vx_t_dbl, vy_t_dbl, &
                                                  vz_t_dbl, &
                                                  omega_t_dbl, age_t_dbl, &
                                                  enth_t_dbl, &
                                                  enh_t_dbl, &
                                                  strain_heating_t_dbl
real(sp), dimension(0:2*IMAX,0:2*JMAX,0:KRMAX) :: temp_r_dbl
character(len=64) :: mapping_grid_mapping_name_dbl, mapping_ellipsoid_dbl

#if (DISC>0)   /* Ice discharge parameterisation */
integer(i4b), dimension(0:2*IMAX,0:2*JMAX) :: mask_mar_dbl
real(sp),     dimension(0:2*IMAX,0:2*JMAX) :: dis_perp_dbl, &
                                              cst_dist_dbl, cos_grad_tc_dbl
#endif

real(dp), parameter :: no_value_neg_dp = -9999.0_dp
real(dp), parameter :: eps_dp = 1.0e-05_dp

end module resolution_doubler_vars

!-------------------------------------------------------------------------------
!> Main program:
!! Doubling the horizontal resolution of a NetCDF time-slice output file
!! produced by SICOPOLIS.
!-------------------------------------------------------------------------------
program resolution_doubler

use resolution_doubler_types
use resolution_doubler_vars

implicit none

integer(i4b) :: i, j, n1, n2
integer(i4b) :: forcing_flag
character(len=256) :: run_name
character(len=  4) :: ergnum

n1 = len('sico_specs_')+1
n2 = len(trim(RUN_SPECS_HEADER))-len('.h')
run_name = trim(RUN_SPECS_HEADER)
run_name = run_name(n1:n2)

call read_nc(run_name, ergnum, forcing_flag)

call double_res_interpol

call write_nc_double(run_name, ergnum, forcing_flag)

contains

!-------------------------------------------------------------------------------
!> Reading of data of time-slice files *.nc (NetCDF format).
!-------------------------------------------------------------------------------
subroutine read_nc(run_name, ergnum, forcing_flag)

use resolution_doubler_types
use resolution_doubler_vars
use netcdf

implicit none

character(len=256), intent(in)  :: run_name

character(len=  4), intent(out) :: ergnum
integer(i4b),       intent(out) :: forcing_flag

character(len=256) :: filename, filename_with_path
character(len=256) :: ch_msg
character, parameter :: end_of_line = char(10)
                        ! End-of-line string

integer(i4b) :: ios
integer(i4b) :: istat, istat1, istat2
logical      :: flag_z_sl_xy_array

integer(i4b) :: ncid, ncv, ncv_test1, ncv_test2
!     ncid:      ID of the NetCDF file
!     ncv:       Variable ID

real(dp) :: r_aux

character(len=64) :: ch_aux

year2sec_erg   = 0.0_dp
time_erg       = 0.0_dp
delta_ts_erg   = 0.0_dp
glac_index_erg = 0.0_dp
z_sl_mean_erg  = 0.0_dp
V_tot_erg      = 0.0_dp
V_af_erg       = 0.0_dp
A_grounded_erg = 0.0_dp
A_floating_erg = 0.0_dp

!-------- Enter name of time-slice file --------

write(6,'(1x,a)') ' '
write(6,'(1x,a,a)') 'Name of run: ', trim(run_name)
write(6,'(1x,a)',advance='no') &
'Number of time-slice file (with leading zeros, 4 digits) > '
read (5,'(a)') ergnum

!-------- Name of time-slice file --------

filename = trim(run_name)//trim(ergnum)//'.nc'

!-------- Reading of data from time-slice file --------

filename_with_path = trim(OUT_PATH)//'/'//trim(filename)

write (6,'(/a)') ' Now reading '//trim(filename_with_path)//' ...'

ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

if (ios /= nf90_noerr) then
   ch_msg = ' >>> read_nc: Error when opening the nc file!'
   call error(ch_msg)
end if

mapping_erg = 1                     ! initial value
mapping_grid_mapping_name_erg = 'xxx'   ! initial value
mapping_ellipsoid_erg         = 'xxx'   ! initial value

mapping_semi_major_axis_erg     = no_value_neg_dp   ! initial value
mapping_inv_flattening_erg      = no_value_neg_dp   ! initial value
mapping_radius_of_sphere_erg    = no_value_neg_dp   ! initial value
mapping_latitude_origin_erg     = no_value_neg_dp   ! initial value
mapping_standard_parallel_erg   = no_value_neg_dp   ! initial value
mapping_reference_longitude_erg = no_value_neg_dp   ! initial value
mapping_false_E_erg             = no_value_neg_dp   ! initial value
mapping_false_N_erg             = no_value_neg_dp   ! initial value

if ( nf90_inq_varid(ncid, 'mapping', ncv_test1) == nf90_noerr ) then
   ncv = ncv_test1
   call check( nf90_get_var(ncid, ncv, mapping_erg) )
else if ( nf90_inq_varid(ncid, 'crs', ncv_test2) == nf90_noerr ) then
   ncv = ncv_test2
   call check( nf90_get_var(ncid, ncv, mapping_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''mapping'' ' &
          //               end_of_line &
          //'              not available in read nc file.'
   call warning(ch_msg)
end if

istat = nf90_get_att(ncid, ncv, 'grid_mapping_name', ch_aux)
if (istat == nf90_noerr) mapping_grid_mapping_name_erg = trim(ch_aux)

istat = nf90_get_att(ncid, ncv, 'ellipsoid', ch_aux)
if (istat == nf90_noerr) mapping_ellipsoid_erg = trim(ch_aux)

istat = nf90_get_att(ncid, ncv, 'semi_major_axis', r_aux)
if (istat == nf90_noerr) then
   mapping_semi_major_axis_erg = r_aux
   istat1 = nf90_get_att(ncid, ncv, 'inverse_flattening', r_aux)
   if (istat1 == nf90_noerr) mapping_inv_flattening_erg = r_aux
else
   istat2 = nf90_get_att(ncid, ncv, 'radius_of_sphere', r_aux)
   if (istat2 == nf90_noerr) mapping_radius_of_sphere_erg = r_aux
end if

#if (GRID==0 || GRID==1)

istat = nf90_get_att(ncid, ncv, 'latitude_of_projection_origin', r_aux)
if (istat == nf90_noerr) mapping_latitude_origin_erg = r_aux

istat = nf90_get_att(ncid, ncv, 'standard_parallel', r_aux)
if (istat == nf90_noerr) mapping_standard_parallel_erg = r_aux

istat = nf90_get_att(ncid, ncv, 'straight_vertical_longitude_from_pole', r_aux)
if (istat == nf90_noerr) mapping_reference_longitude_erg = r_aux

istat = nf90_get_att(ncid, ncv, 'false_easting', r_aux)
if (istat == nf90_noerr) mapping_false_E_erg = r_aux

istat = nf90_get_att(ncid, ncv, 'false_northing', r_aux)
if (istat == nf90_noerr) mapping_false_N_erg = r_aux

#endif

istat = nf90_inq_varid(ncid, 'year2sec', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, year2sec_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''year2sec'' ' &
          //               end_of_line &
          //'              not available in read nc file.'
   call warning(ch_msg)
   year2sec_erg = 3.1556925445e+07_dp   ! default value
                     ! IUPAC-IUGS year for epoch 2000.0
                     ! (Holden et al., 2011, PAC, doi:10.1351/PAC-REC-09-01-22)
end if

istat = nf90_inq_varid(ncid, 'time', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, time_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''time'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

if ( nf90_inq_varid(ncid, 'delta_ts', ncv) == nf90_noerr ) then
   forcing_flag = 1
   call check( nf90_get_var(ncid, ncv, delta_ts_erg) )
else if ( nf90_inq_varid(ncid, 'glac_index', ncv) == nf90_noerr ) then
   forcing_flag = 2
   call check( nf90_get_var(ncid, ncv, glac_index_erg) )
else
   ch_msg = ' >>> read_nc: Neither variable ''delta_ts'' nor ''glac_index'' ' &
          //               end_of_line &
          //'              available in read nc file.'
   call warning(ch_msg)
   forcing_flag = 0
   delta_ts_erg   = 0.0_dp
   glac_index_erg = 0.0_dp
end if

istat1 = nf90_inq_varid(ncid, 'z_sl_mean', ncv)
if (istat1 == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, z_sl_mean_erg) )
   flag_z_sl_xy_array = .true.
else
   istat2 = nf90_inq_varid(ncid, 'z_sl', ncv)
   if (istat2 == nf90_noerr) then
      call check( nf90_get_var(ncid, ncv, z_sl_mean_erg) )
      flag_z_sl_xy_array = .false.
   else
      ch_msg = ' >>> read_nc: Variable ''z_sl(_mean)'' ' &
             //               end_of_line &
             //'              not available in read nc file!'
      call error(ch_msg)
   end if
end if

istat = nf90_inq_varid(ncid, 'V_tot', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, V_tot_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''V_tot'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'V_af', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, V_af_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''V_af'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'A_grounded', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, A_grounded_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''A_grounded'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'A_floating', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, A_floating_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''A_floating'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'x', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, xi_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''x'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'y', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, eta_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''y'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'sigma_level_c', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, sigma_level_c_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''sigma_level_c'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'sigma_level_t', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, sigma_level_t_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''sigma_level_t'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'sigma_level_r', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, sigma_level_r_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''sigma_level_r'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'lambda', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, lambda_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''lambda'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'phi', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, phi_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''phi'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'lon', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, lon_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''lon'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'lat', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, lat_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''lat'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'temp_s', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, temp_s_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''temp_s'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'prec', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, prec_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''prec'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'snowfall', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, snowfall_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''snowfall'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'rainfall', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, rainfall_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''rainfall'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'pdd', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, pdd_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''pdd'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'as_perp', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, as_perp_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''as_perp'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'as_perp_apl', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, as_perp_apl_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''as_perp_apl'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'smb_corr', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, smb_corr_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''smb_corr'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

if (flag_z_sl_xy_array) then
   istat = nf90_inq_varid(ncid, 'z_sl', ncv)
   if (istat == nf90_noerr) then
      call check( nf90_get_var(ncid, ncv, z_sl_erg) )
   else
      ch_msg = ' >>> read_nc: Variable ''z_sl'' ' &
             //               end_of_line &
             //'              not available in read nc file!'
      call error(ch_msg)
   end if
else
   z_sl_erg = z_sl_mean_erg
end if

#if (DISC>0)   /* Ice discharge parameterisation */

istat = nf90_inq_varid(ncid, 'dis_perp', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, dis_perp_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''dis_perp'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'cst_dist', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, cst_dist_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''cst_dist'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'cos_grad_tc', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, cos_grad_tc_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''cos_grad_tc'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'mask_mar', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, mask_mar_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''mask_mar'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

#endif

istat = nf90_inq_varid(ncid, 'q_geo', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, q_geo_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''q_geo'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

if ( nf90_inq_varid(ncid, 'mask', ncv) == nf90_noerr ) then
   call check( nf90_get_var(ncid, ncv, mask_erg) )
else if ( nf90_inq_varid(ncid, 'maske', ncv) == nf90_noerr ) then
   call check( nf90_get_var(ncid, ncv, mask_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''mask'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

if ( nf90_inq_varid(ncid, 'mask_old', ncv) == nf90_noerr ) then
   call check( nf90_get_var(ncid, ncv, mask_old_erg) )
else if ( nf90_inq_varid(ncid, 'maske_old', ncv) == nf90_noerr ) then
   call check( nf90_get_var(ncid, ncv, mask_old_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''mask_old'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'n_cts', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, n_cts_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''n_cts'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'kc_cts', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, kc_cts_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''kc_cts'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'zs', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, zs_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''zs'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'zm', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, zm_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''zm'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'zb', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, zb_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''zb'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'zl', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, zl_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''zl'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'zl0', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, zl0_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''zl0'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'wss', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, wss_erg) )
   flag_wss = .true.
else
   ch_msg = ' >>> read_nc: Variable ''wss'' ' &
          //               end_of_line &
          //'              not available in read nc file.'
   call warning(ch_msg)
   wss_erg = 0.0_sp
   flag_wss = .false.
end if

istat = nf90_inq_varid(ncid, 'H_cold', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, H_cold_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''H_cold'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'H_temp', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, H_temp_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''H_temp'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'H', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, H_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''H'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'H_R', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, H_R_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''H_R'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'Q_bm', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, Q_bm_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''Q_bm'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'Q_tld', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, Q_tld_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''Q_tld'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat1 = nf90_inq_varid(ncid, 'calving', ncv)
if (istat1 == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, calving_erg) )
else
   istat2 = nf90_inq_varid(ncid, 'q_cf_g', ncv)
   if (istat2 == nf90_noerr) then
      call check( nf90_get_var(ncid, ncv, calving_erg) )
   end if
end if

istat = nf90_inq_varid(ncid, 'am_perp', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, am_perp_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''am_perp'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'qx', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, qx_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''qx'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'qy', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, qy_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''qy'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'vx_m_sia', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, vx_m_sia_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''vx_m_sia'' ' &
          //               end_of_line &
          //'              not available in read nc file.'
   call warning(ch_msg)
   vx_m_sia_erg = 0.0_sp
end if

istat = nf90_inq_varid(ncid, 'vy_m_sia', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, vy_m_sia_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''vy_m_sia'' ' &
          //               end_of_line &
          //'              not available in read nc file.'
   call warning(ch_msg)
   vy_m_sia_erg = 0.0_sp
end if

istat = nf90_inq_varid(ncid, 'vx_m_ssa', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, vx_m_ssa_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''vx_m_ssa'' ' &
          //               end_of_line &
          //'              not available in read nc file.'
   call warning(ch_msg)
   vx_m_ssa_erg = 0.0_sp
end if

istat = nf90_inq_varid(ncid, 'vy_m_ssa', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, vy_m_ssa_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''vy_m_ssa'' ' &
          //               end_of_line &
          //'              not available in read nc file.'
   call warning(ch_msg)
   vy_m_ssa_erg = 0.0_sp
end if

istat = nf90_inq_varid(ncid, 'dzs_dt', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, dzs_dtau_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''dzs_dt'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'dzm_dt', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, dzm_dtau_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''dzm_dt'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'dzb_dt', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, dzb_dtau_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''dzb_dt'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'dzl_dt', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, dzl_dtau_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''dzl_dt'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'dH_c_dt', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, dH_c_dtau_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''dH_c_dt'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'dH_t_dt', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, dH_t_dtau_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''dH_t_dt'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'dH_dt', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, dH_dtau_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''dH_dt'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'vx_b_g', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, vx_b_g_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''vx_b_g'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'vy_b_g', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, vy_b_g_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''vy_b_g'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'vz_b', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, vz_b_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''vz_b'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'vh_b', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, vh_b_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''vh_b'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'vx_s_g', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, vx_s_g_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''vx_s_g'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'vy_s_g', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, vy_s_g_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''vy_s_g'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'vz_s', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, vz_s_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''vz_s'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'vh_s', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, vh_s_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''vh_s'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'vx_m_g', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, vx_m_g_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''vx_m_g'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'vy_m_g', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, vy_m_g_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''vy_m_g'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'vh_m', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, vh_m_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''vh_m'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'temp_b', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, temp_b_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''temp_b'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'temph_b', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, temph_b_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''temph_b'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'tau_dr', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, tau_dr_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''tau_dr'' ' &
          //               end_of_line &
          //'              not available in read nc file.'
   call warning(ch_msg)
   tau_dr_erg = 0.0_sp
end if

istat = nf90_inq_varid(ncid, 'tau_b', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, tau_b_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''tau_b'' ' &
          //               end_of_line &
          //'              not available in read nc file.'
   call warning(ch_msg)
   tau_b_erg = 0.0_sp
end if

istat = nf90_inq_varid(ncid, 'p_b_w', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, p_b_w_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''p_b_w'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'q_w', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, q_w_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''q_w'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'q_w_x', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, q_w_x_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''q_w_x'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'q_w_y', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, q_w_y_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''q_w_y'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'H_w', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, H_w_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''H_w'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'q_gl_g', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, q_gl_g_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''q_gl_g'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'ratio_sl_sia_x', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, ratio_sl_sia_x_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''ratio_sl_sia_x'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'ratio_sl_sia_y', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, ratio_sl_sia_y_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''ratio_sl_sia_y'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'flag_shelfy_stream_x', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, flag_shelfy_stream_x_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''flag_shelfy_stream_x'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'flag_shelfy_stream_y', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, flag_shelfy_stream_y_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''flag_shelfy_stream_y'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'flag_shelfy_stream', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, flag_shelfy_stream_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''flag_shelfy_stream'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'flag_grounding_line_1', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, flag_grounding_line_1_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''flag_grounding_line_1'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'flag_grounding_line_2', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, flag_grounding_line_2_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''flag_grounding_line_2'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'flag_calving_front_1', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, flag_calving_front_1_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''flag_calving_front_1'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'flag_calving_front_2', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, flag_calving_front_2_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''flag_calving_front_2'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'flag_grounded_front_a_1', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, flag_grounded_front_a_1_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''flag_calving_front_a_1'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'flag_grounded_front_a_2', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, flag_grounded_front_a_2_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''flag_calving_front_a_2'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'flag_grounded_front_b_1', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, flag_grounded_front_b_1_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''flag_calving_front_b_1'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'flag_grounded_front_b_2', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, flag_grounded_front_b_2_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''flag_calving_front_b_2'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'vis_ave_g', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, vis_ave_g_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''vis_ave_g'' ' &
          //               end_of_line &
          //'              not available in read nc file.'
   call warning(ch_msg)
   vis_ave_g_erg = 0.0_sp
end if

istat = nf90_inq_varid(ncid, 'vis_int_g', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, vis_int_g_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''vis_int_g'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

#if (OUTPUT==3 || ERGDAT==1)

istat = nf90_inq_varid(ncid, 'vx_c', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, vx_c_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''vx_c'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'vy_c', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, vy_c_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''vy_c'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'vz_c', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, vz_c_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''vz_c'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'vx_t', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, vx_t_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''vx_t'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'vy_t', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, vy_t_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''vy_t'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'vz_t', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, vz_t_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''vz_t'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'temp_c', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, temp_c_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''temp_c'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'omega_t', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, omega_t_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''omega_t'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'temp_r', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, temp_r_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''temp_r'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'enth_c', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, enth_c_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''enth_c'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'enth_t', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, enth_t_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''enth_t'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'omega_c', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, omega_c_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''omega_c'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'enh_c', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, enh_c_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''enh_c'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'enh_t', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, enh_t_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''enh_t'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'strain_heating_c', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, strain_heating_c_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''strain_heating_c'' ' &
          //               end_of_line &
          //'              not available in read nc file.'
   call warning(ch_msg)
   strain_heating_c_erg = 0.0_sp
end if

istat = nf90_inq_varid(ncid, 'strain_heating_t', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, strain_heating_t_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''strain_heating_t'' ' &
          //               end_of_line &
          //'              not available in read nc file.'
   call warning(ch_msg)
   strain_heating_t_erg = 0.0_sp
end if

istat = nf90_inq_varid(ncid, 'age_c', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, age_c_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''age_c'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

istat = nf90_inq_varid(ncid, 'age_t', ncv)
if (istat == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, age_t_erg) )
else
   ch_msg = ' >>> read_nc: Variable ''age_t'' ' &
          //               end_of_line &
          //'              not available in read nc file!'
   call error(ch_msg)
end if

#else

ch_msg = ' >>> read_nc: Resolution doubling requires 3D fields!'
call error(ch_msg)

#endif

call check( nf90_close(ncid) )

end subroutine read_nc

!-------------------------------------------------------------------------------
!> Interpolation of data to double resolution.
!-------------------------------------------------------------------------------
subroutine double_res_interpol

use resolution_doubler_types
use resolution_doubler_vars

implicit none

integer(i4b) :: i, j, i1, i2, j1, j2, ii, jj, kc, kt, kr

write (6,'(/a)') ' Now interpolating data to double resolution ...'

!-------- Interpolation for scalars --------

mapping_dbl                     = mapping_erg
mapping_grid_mapping_name_dbl   = trim(mapping_grid_mapping_name_erg)
mapping_ellipsoid_dbl           = trim(mapping_ellipsoid_erg)
mapping_semi_major_axis_dbl     = mapping_semi_major_axis_erg
mapping_inv_flattening_dbl      = mapping_inv_flattening_erg
mapping_radius_of_sphere_dbl    = mapping_radius_of_sphere_erg
mapping_latitude_origin_dbl     = mapping_latitude_origin_erg
mapping_standard_parallel_dbl   = mapping_standard_parallel_erg
mapping_reference_longitude_dbl = mapping_reference_longitude_erg
mapping_false_E_dbl             = mapping_false_E_erg
mapping_false_N_dbl             = mapping_false_N_erg

year2sec_dbl   = year2sec_erg
time_dbl       = time_erg
delta_ts_dbl   = delta_ts_erg
glac_index_dbl = glac_index_erg
z_sl_mean_dbl  = z_sl_mean_erg
V_tot_dbl      = V_tot_erg
V_af_dbl       = V_af_erg
A_grounded_dbl = A_grounded_erg
A_floating_dbl = A_floating_erg
H_R_dbl        = H_R_erg

!-------- Interpolation for 1D fields --------

do ii = 0, 2*IMAX, 2
   i  = ii/2
   xi_dbl(ii) = xi_erg(i)
end do

do ii = 1, 2*IMAX-1, 2
   i1 = (ii-1)/2
   i2 = (ii+1)/2
   xi_dbl(ii) = 0.5*(xi_erg(i1)+xi_erg(i2))
end do

do jj = 0, 2*JMAX, 2
   j  = jj/2
   eta_dbl(jj) = eta_erg(j)
end do

do jj = 1, 2*JMAX-1, 2
   j1 = (jj-1)/2
   j2 = (jj+1)/2
   eta_dbl(jj) = 0.5*(eta_erg(j1)+eta_erg(j2))
end do

sigma_level_c_dbl = sigma_level_c_erg
sigma_level_t_dbl = sigma_level_t_erg
sigma_level_r_dbl = sigma_level_r_erg

!-------- Interpolation for 2D and 3D fields --------

r_kc_cts_erg = real(kc_cts_erg,sp)

do ii = 0, 2*IMAX, 2
do jj = 0, 2*JMAX, 2
   i  = ii/2
   j  = jj/2
   lambda_dbl(ii,jj)    = lambda_erg(i,j)
   phi_dbl(ii,jj)       = phi_erg(i,j)
   lon_dbl(ii,jj)       = lon_erg(i,j)
   lat_dbl(ii,jj)       = lat_erg(i,j)
   temp_s_dbl(ii,jj)    = temp_s_erg(i,j)
   prec_dbl(ii,jj)      = prec_erg(i,j)
   snowfall_dbl(ii,jj)  = snowfall_erg(i,j)
   rainfall_dbl(ii,jj)  = rainfall_erg(i,j)
   pdd_dbl(ii,jj)       = pdd_erg(i,j)
   as_perp_dbl(ii,jj)   = as_perp_erg(i,j)
   as_perp_apl_dbl(ii,jj) = as_perp_apl_erg(i,j)
   smb_corr_dbl(ii,jj)  = smb_corr_erg(i,j)
   z_sl_dbl(ii,jj)      = z_sl_erg(i,j)
#if (DISC>0)   /* Ice discharge parameterisation */
   dis_perp_dbl(ii,jj)    = dis_perp_erg(i,j)
   cst_dist_dbl(ii,jj)    = cst_dist_erg(i,j)
   cos_grad_tc_dbl(ii,jj) = cos_grad_tc_erg(i,j)
#endif
   q_geo_dbl(ii,jj)     = q_geo_erg(i,j)
   zs_dbl(ii,jj)        = zs_erg(i,j)
   zm_dbl(ii,jj)        = zm_erg(i,j)
   zb_dbl(ii,jj)        = zb_erg(i,j)
   zl_dbl(ii,jj)        = zl_erg(i,j)
   zl0_dbl(ii,jj)       = zl0_erg(i,j)
   wss_dbl(ii,jj)       = wss_erg(i,j)
   H_cold_dbl(ii,jj)    = H_cold_erg(i,j)
   H_temp_dbl(ii,jj)    = H_temp_erg(i,j)
   H_dbl(ii,jj)         = H_erg(i,j)
   r_kc_cts_dbl(ii,jj)  = r_kc_cts_erg(i,j)
   vz_c_dbl(ii,jj,:)    = vz_c_erg(i,j,:)
   vz_t_dbl(ii,jj,:)    = vz_t_erg(i,j,:)
   temp_c_dbl(ii,jj,:)  = temp_c_erg(i,j,:)
   age_c_dbl(ii,jj,:)   = age_c_erg(i,j,:)
   omega_t_dbl(ii,jj,:) = omega_t_erg(i,j,:)
   age_t_dbl(ii,jj,:)   = age_t_erg(i,j,:)
   temp_r_dbl(ii,jj,:)  = temp_r_erg(i,j,:)
   enth_c_dbl(ii,jj,:)  = enth_c_erg(i,j,:)
   enth_t_dbl(ii,jj,:)  = enth_t_erg(i,j,:)
   omega_c_dbl(ii,jj,:) = omega_c_erg(i,j,:)
   enh_c_dbl(ii,jj,:)   = enh_c_erg(i,j,:)
   enh_t_dbl(ii,jj,:)   = enh_t_erg(i,j,:)
   strain_heating_c_dbl(ii,jj,:) = strain_heating_c_erg(i,j,:)
   strain_heating_t_dbl(ii,jj,:) = strain_heating_t_erg(i,j,:)
   vx_b_g_dbl(ii,jj)    = vx_b_g_erg(i,j)
   vy_b_g_dbl(ii,jj)    = vy_b_g_erg(i,j)
   vz_b_dbl(ii,jj)      = vz_b_erg(i,j)
   vh_b_dbl(ii,jj)      = vh_b_erg(i,j)
   vx_s_g_dbl(ii,jj)    = vx_s_g_erg(i,j)
   vy_s_g_dbl(ii,jj)    = vy_s_g_erg(i,j)
   vz_s_dbl(ii,jj)      = vz_s_erg(i,j)
   vh_s_dbl(ii,jj)      = vh_s_erg(i,j)
   vx_m_g_dbl(ii,jj)    = vx_m_g_erg(i,j)
   vy_m_g_dbl(ii,jj)    = vy_m_g_erg(i,j)
   vh_m_dbl(ii,jj)      = vh_m_erg(i,j)
   temp_b_dbl(ii,jj)    = temp_b_erg(i,j)
   temph_b_dbl(ii,jj)   = temph_b_erg(i,j)
   Q_bm_dbl(ii,jj)      = Q_bm_erg(i,j)
   Q_tld_dbl(ii,jj)     = Q_tld_erg(i,j)
   calving_dbl(ii,jj)   = calving_erg(i,j)
   am_perp_dbl(ii,jj)   = am_perp_erg(i,j)
   dzs_dtau_dbl(ii,jj)  = dzs_dtau_erg(i,j)
   dzm_dtau_dbl(ii,jj)  = dzm_dtau_erg(i,j)
   dzb_dtau_dbl(ii,jj)  = dzb_dtau_erg(i,j)
   dzl_dtau_dbl(ii,jj)  = dzl_dtau_erg(i,j)
   dH_c_dtau_dbl(ii,jj) = dH_c_dtau_erg(i,j)
   dH_t_dtau_dbl(ii,jj) = dH_t_dtau_erg(i,j)
   dH_dtau_dbl(ii,jj)   = dH_dtau_erg(i,j)
   tau_dr_dbl(ii,jj)    = tau_dr_erg(i,j)
   tau_b_dbl(ii,jj)     = tau_b_erg(i,j)
   p_b_w_dbl(ii,jj)     = p_b_w_erg(i,j)
   q_w_dbl(ii,jj)       = q_w_erg(i,j)
   q_w_x_dbl(ii,jj)     = q_w_x_erg(i,j)
   q_w_y_dbl(ii,jj)     = q_w_y_erg(i,j)
   H_w_dbl(ii,jj)       = H_w_erg(i,j)
   q_gl_g_dbl(ii,jj)    = q_gl_g_erg(i,j)
   ratio_sl_sia_x_dbl(ii,jj) = ratio_sl_sia_x_erg(i,j)
   ratio_sl_sia_y_dbl(ii,jj) = ratio_sl_sia_y_erg(i,j)
   vis_ave_g_dbl(ii,jj) = vis_ave_g_erg(i,j)
   vis_int_g_dbl(ii,jj) = vis_int_g_erg(i,j)

end do
end do

do ii = 1, 2*IMAX-1, 2
do jj = 0, 2*JMAX, 2
   i1 = (ii-1)/2
   i2 = (ii+1)/2
   j  = jj/2
   lambda_dbl(ii,jj)    = 0.5*(lambda_erg(i1,j)+lambda_erg(i2,j))
   phi_dbl(ii,jj)       = 0.5*(phi_erg(i1,j)+phi_erg(i2,j))
   lon_dbl(ii,jj)       = 0.5*(lon_erg(i1,j)+lon_erg(i2,j))
   lat_dbl(ii,jj)       = 0.5*(lat_erg(i1,j)+lat_erg(i2,j))
   temp_s_dbl(ii,jj)    = 0.5*(temp_s_erg(i1,j)+temp_s_erg(i2,j))
   prec_dbl(ii,jj)      = 0.5*(prec_erg(i1,j)+prec_erg(i2,j))
   snowfall_dbl(ii,jj)  = 0.5*(snowfall_erg(i1,j)+snowfall_erg(i2,j))
   rainfall_dbl(ii,jj)  = 0.5*(rainfall_erg(i1,j)+rainfall_erg(i2,j))
   pdd_dbl(ii,jj)       = 0.5*(pdd_erg(i1,j)+pdd_erg(i2,j))
   as_perp_dbl(ii,jj)   = 0.5*(as_perp_erg(i1,j)+as_perp_erg(i2,j))
   as_perp_apl_dbl(ii,jj) = 0.5*(as_perp_apl_erg(i1,j)+as_perp_apl_erg(i2,j))
   smb_corr_dbl(ii,jj)  = 0.5*(smb_corr_erg(i1,j)+smb_corr_erg(i2,j))
   z_sl_dbl(ii,jj)      = 0.5*(z_sl_erg(i1,j)+z_sl_erg(i2,j))
#if (DISC>0)   /* Ice discharge parameterisation */
   dis_perp_dbl(ii,jj)    = 0.5*(dis_perp_erg(i1,j)+dis_perp_erg(i2,j))
   cst_dist_dbl(ii,jj)    = 0.5*(cst_dist_erg(i1,j)+cst_dist_erg(i2,j))
   cos_grad_tc_dbl(ii,jj) = 0.5*(cos_grad_tc_erg(i1,j)+cos_grad_tc_erg(i2,j))
#endif
   q_geo_dbl(ii,jj)     = 0.5*(q_geo_erg(i1,j)+q_geo_erg(i2,j))
   zs_dbl(ii,jj)        = 0.5*(zs_erg(i1,j)+zs_erg(i2,j))
   zm_dbl(ii,jj)        = 0.5*(zm_erg(i1,j)+zm_erg(i2,j))
   zb_dbl(ii,jj)        = 0.5*(zb_erg(i1,j)+zb_erg(i2,j))
   zl_dbl(ii,jj)        = 0.5*(zl_erg(i1,j)+zl_erg(i2,j))
   zl0_dbl(ii,jj)       = 0.5*(zl0_erg(i1,j)+zl0_erg(i2,j))
   wss_dbl(ii,jj)       = 0.5*(wss_erg(i1,j)+wss_erg(i2,j))
   H_cold_dbl(ii,jj)    = 0.5*(H_cold_erg(i1,j)+H_cold_erg(i2,j))
   H_temp_dbl(ii,jj)    = 0.5*(H_temp_erg(i1,j)+H_temp_erg(i2,j))
   H_dbl(ii,jj)         = 0.5*(H_erg(i1,j)+H_erg(i2,j))
   r_kc_cts_dbl(ii,jj)  = 0.5*(r_kc_cts_erg(i1,j)+r_kc_cts_erg(i2,j))
   vz_c_dbl(ii,jj,:)    = 0.5*(vz_c_erg(i1,j,:)+vz_c_erg(i2,j,:))
   vz_t_dbl(ii,jj,:)    = 0.5*(vz_t_erg(i1,j,:)+vz_t_erg(i2,j,:))
   temp_c_dbl(ii,jj,:)  = 0.5*(temp_c_erg(i1,j,:)+temp_c_erg(i2,j,:))
   age_c_dbl(ii,jj,:)   = 0.5*(age_c_erg(i1,j,:)+age_c_erg(i2,j,:))
   omega_t_dbl(ii,jj,:) = 0.5*(omega_t_erg(i1,j,:)+omega_t_erg(i2,j,:))
   age_t_dbl(ii,jj,:)   = 0.5*(age_t_erg(i1,j,:)+age_t_erg(i2,j,:))
   temp_r_dbl(ii,jj,:)  = 0.5*(temp_r_erg(i1,j,:)+temp_r_erg(i2,j,:))
   enth_c_dbl(ii,jj,:)  = 0.5*(enth_c_erg(i1,j,:)+enth_c_erg(i2,j,:))
   enth_t_dbl(ii,jj,:)  = 0.5*(enth_t_erg(i1,j,:)+enth_t_erg(i2,j,:))
   omega_c_dbl(ii,jj,:) = 0.5*(omega_c_erg(i1,j,:)+omega_c_erg(i2,j,:))
   enh_c_dbl(ii,jj,:)   = 0.5*(enh_c_erg(i1,j,:)+enh_c_erg(i2,j,:))
   enh_t_dbl(ii,jj,:)   = 0.5*(enh_t_erg(i1,j,:)+enh_t_erg(i2,j,:))
   strain_heating_c_dbl(ii,jj,:) = 0.5*( strain_heating_c_erg(i1,j,:) &
                                        +strain_heating_c_erg(i2,j,:))
   strain_heating_t_dbl(ii,jj,:) = 0.5*( strain_heating_t_erg(i1,j,:) &
                                        +strain_heating_t_erg(i2,j,:))
   vx_b_g_dbl(ii,jj)    = 0.5*(vx_b_g_erg(i1,j)+vx_b_g_erg(i2,j))
   vy_b_g_dbl(ii,jj)    = 0.5*(vy_b_g_erg(i1,j)+vy_b_g_erg(i2,j))
   vz_b_dbl(ii,jj)      = 0.5*(vz_b_erg(i1,j)+vz_b_erg(i2,j))
   vh_b_dbl(ii,jj)      = 0.5*(vh_b_erg(i1,j)+vh_b_erg(i2,j))
   vx_s_g_dbl(ii,jj)    = 0.5*(vx_s_g_erg(i1,j)+vx_s_g_erg(i2,j))
   vy_s_g_dbl(ii,jj)    = 0.5*(vy_s_g_erg(i1,j)+vy_s_g_erg(i2,j))
   vz_s_dbl(ii,jj)      = 0.5*(vz_s_erg(i1,j)+vz_s_erg(i2,j))
   vh_s_dbl(ii,jj)      = 0.5*(vh_s_erg(i1,j)+vh_s_erg(i2,j))
   vx_m_g_dbl(ii,jj)    = 0.5*(vx_m_g_erg(i1,j)+vx_m_g_erg(i2,j))
   vy_m_g_dbl(ii,jj)    = 0.5*(vy_m_g_erg(i1,j)+vy_m_g_erg(i2,j))
   vh_m_dbl(ii,jj)      = 0.5*(vh_m_erg(i1,j)+vh_m_erg(i2,j))
   temp_b_dbl(ii,jj)    = 0.5*(temp_b_erg(i1,j)+temp_b_erg(i2,j))
   temph_b_dbl(ii,jj)   = 0.5*(temph_b_erg(i1,j)+temph_b_erg(i2,j))
   Q_bm_dbl(ii,jj)      = 0.5*(Q_bm_erg(i1,j)+Q_bm_erg(i2,j))
   Q_tld_dbl(ii,jj)     = 0.5*(Q_tld_erg(i1,j)+Q_tld_erg(i2,j))
   calving_dbl(ii,jj)   = 0.5*(calving_erg(i1,j)+calving_erg(i2,j))
   am_perp_dbl(ii,jj)   = 0.5*(am_perp_erg(i1,j)+am_perp_erg(i2,j))
   dzs_dtau_dbl(ii,jj)  = 0.5*(dzs_dtau_erg(i1,j)+dzs_dtau_erg(i2,j))
   dzm_dtau_dbl(ii,jj)  = 0.5*(dzm_dtau_erg(i1,j)+dzm_dtau_erg(i2,j))
   dzb_dtau_dbl(ii,jj)  = 0.5*(dzb_dtau_erg(i1,j)+dzb_dtau_erg(i2,j))
   dzl_dtau_dbl(ii,jj)  = 0.5*(dzl_dtau_erg(i1,j)+dzl_dtau_erg(i2,j))
   dH_c_dtau_dbl(ii,jj) = 0.5*(dH_c_dtau_erg(i1,j)+dH_c_dtau_erg(i2,j))
   dH_t_dtau_dbl(ii,jj) = 0.5*(dH_t_dtau_erg(i1,j)+dH_t_dtau_erg(i2,j))
   dH_dtau_dbl(ii,jj)   = 0.5*(dH_dtau_erg(i1,j)+dH_dtau_erg(i2,j))
   tau_dr_dbl(ii,jj)    = 0.5*(tau_dr_erg(i1,j)+tau_dr_erg(i2,j))
   tau_b_dbl(ii,jj)     = 0.5*(tau_b_erg(i1,j)+tau_b_erg(i2,j))
   p_b_w_dbl(ii,jj)     = 0.5*(p_b_w_erg(i1,j)+p_b_w_erg(i2,j))
   q_w_dbl(ii,jj)       = 0.5*(q_w_erg(i1,j)+q_w_erg(i2,j))
   q_w_x_dbl(ii,jj)     = 0.5*(q_w_x_erg(i1,j)+q_w_x_erg(i2,j))
   q_w_y_dbl(ii,jj)     = 0.5*(q_w_y_erg(i1,j)+q_w_y_erg(i2,j))
   H_w_dbl(ii,jj)       = 0.5*(H_w_erg(i1,j)+H_w_erg(i2,j))
   q_gl_g_dbl(ii,jj)    = 0.5*(q_gl_g_erg(i1,j)+q_gl_g_erg(i2,j))
   ratio_sl_sia_x_dbl(ii,jj) = 0.5*( ratio_sl_sia_x_erg(i1,j) &
                                    +ratio_sl_sia_x_erg(i2,j))
   ratio_sl_sia_y_dbl(ii,jj) = 0.5*( ratio_sl_sia_y_erg(i1,j) &
                                    +ratio_sl_sia_y_erg(i2,j))
   vis_ave_g_dbl(ii,jj) = 0.5*(vis_ave_g_erg(i1,j)+vis_ave_g_erg(i2,j))
   vis_int_g_dbl(ii,jj) = 0.5*(vis_int_g_erg(i1,j)+vis_int_g_erg(i2,j))
end do
end do

do ii = 0, 2*IMAX, 2
do jj = 1, 2*JMAX-1, 2
   i  = ii/2
   j1 = (jj-1)/2
   j2 = (jj+1)/2
   lambda_dbl(ii,jj)    = 0.5*(lambda_erg(i,j1)+lambda_erg(i,j2))
   phi_dbl(ii,jj)       = 0.5*(phi_erg(i,j1)+phi_erg(i,j2))
   lon_dbl(ii,jj)       = 0.5*(lon_erg(i,j1)+lon_erg(i,j2))
   lat_dbl(ii,jj)       = 0.5*(lat_erg(i,j1)+lat_erg(i,j2))
   temp_s_dbl(ii,jj)    = 0.5*(temp_s_erg(i,j1)+temp_s_erg(i,j2))
   prec_dbl(ii,jj)      = 0.5*(prec_erg(i,j1)+prec_erg(i,j2))
   snowfall_dbl(ii,jj)  = 0.5*(snowfall_erg(i,j1)+snowfall_erg(i,j2))
   rainfall_dbl(ii,jj)  = 0.5*(rainfall_erg(i,j1)+rainfall_erg(i,j2))
   pdd_dbl(ii,jj)       = 0.5*(pdd_erg(i,j1)+pdd_erg(i,j2))
   as_perp_dbl(ii,jj)   = 0.5*(as_perp_erg(i,j1)+as_perp_erg(i,j2))
   as_perp_apl_dbl(ii,jj) = 0.5*(as_perp_apl_erg(i,j1)+as_perp_apl_erg(i,j2))
   smb_corr_dbl(ii,jj)  = 0.5*(smb_corr_erg(i,j1)+smb_corr_erg(i,j2))
   z_sl_dbl(ii,jj)      = 0.5*(z_sl_erg(i,j1)+z_sl_erg(i,j2))
#if (DISC>0)   /* Ice discharge parameterisation */
   dis_perp_dbl(ii,jj)    = 0.5*(dis_perp_erg(i,j1)+dis_perp_erg(i,j2))
   cst_dist_dbl(ii,jj)    = 0.5*(cst_dist_erg(i,j1)+cst_dist_erg(i,j2))
   cos_grad_tc_dbl(ii,jj) = 0.5*(cos_grad_tc_erg(i,j1)+cos_grad_tc_erg(i,j2))
#endif
   q_geo_dbl(ii,jj)     = 0.5*(q_geo_erg(i,j1)+q_geo_erg(i,j2))
   zs_dbl(ii,jj)        = 0.5*(zs_erg(i,j1)+zs_erg(i,j2))
   zm_dbl(ii,jj)        = 0.5*(zm_erg(i,j1)+zm_erg(i,j2))
   zb_dbl(ii,jj)        = 0.5*(zb_erg(i,j1)+zb_erg(i,j2))
   zl_dbl(ii,jj)        = 0.5*(zl_erg(i,j1)+zl_erg(i,j2))
   zl0_dbl(ii,jj)       = 0.5*(zl0_erg(i,j1)+zl0_erg(i,j2))
   wss_dbl(ii,jj)       = 0.5*(wss_erg(i,j1)+wss_erg(i,j2))
   H_cold_dbl(ii,jj)    = 0.5*(H_cold_erg(i,j1)+H_cold_erg(i,j2))
   H_temp_dbl(ii,jj)    = 0.5*(H_temp_erg(i,j1)+H_temp_erg(i,j2))
   H_dbl(ii,jj)         = 0.5*(H_erg(i,j1)+H_erg(i,j2))
   r_kc_cts_dbl(ii,jj)  = 0.5*(r_kc_cts_erg(i,j1)+r_kc_cts_erg(i,j2))
   vz_c_dbl(ii,jj,:)    = 0.5*(vz_c_erg(i,j1,:)+vz_c_erg(i,j2,:))
   vz_t_dbl(ii,jj,:)    = 0.5*(vz_t_erg(i,j1,:)+vz_t_erg(i,j2,:))
   temp_c_dbl(ii,jj,:)  = 0.5*(temp_c_erg(i,j1,:)+temp_c_erg(i,j2,:))
   age_c_dbl(ii,jj,:)   = 0.5*(age_c_erg(i,j1,:)+age_c_erg(i,j2,:))
   omega_t_dbl(ii,jj,:) = 0.5*(omega_t_erg(i,j1,:)+omega_t_erg(i,j2,:))
   age_t_dbl(ii,jj,:)   = 0.5*(age_t_erg(i,j1,:)+age_t_erg(i,j2,:))
   temp_r_dbl(ii,jj,:)  = 0.5*(temp_r_erg(i,j1,:)+temp_r_erg(i,j2,:))
   enth_c_dbl(ii,jj,:)  = 0.5*(enth_c_erg(i,j1,:)+enth_c_erg(i,j2,:))
   enth_t_dbl(ii,jj,:)  = 0.5*(enth_t_erg(i,j1,:)+enth_t_erg(i,j2,:))
   omega_c_dbl(ii,jj,:) = 0.5*(omega_c_erg(i,j1,:)+omega_c_erg(i,j2,:))
   enh_c_dbl(ii,jj,:)   = 0.5*(enh_c_erg(i,j1,:)+enh_c_erg(i,j2,:))
   enh_t_dbl(ii,jj,:)   = 0.5*(enh_t_erg(i,j1,:)+enh_t_erg(i,j2,:))
   strain_heating_c_dbl(ii,jj,:) = 0.5*( strain_heating_c_erg(i,j1,:) &
                                        +strain_heating_c_erg(i,j2,:))
   strain_heating_t_dbl(ii,jj,:) = 0.5*( strain_heating_t_erg(i,j1,:) &
                                        +strain_heating_t_erg(i,j2,:))
   vx_b_g_dbl(ii,jj)    = 0.5*(vx_b_g_erg(i,j1)+vx_b_g_erg(i,j2))
   vy_b_g_dbl(ii,jj)    = 0.5*(vy_b_g_erg(i,j1)+vy_b_g_erg(i,j2))
   vz_b_dbl(ii,jj)      = 0.5*(vz_b_erg(i,j1)+vz_b_erg(i,j2))
   vh_b_dbl(ii,jj)      = 0.5*(vh_b_erg(i,j1)+vh_b_erg(i,j2))
   vx_s_g_dbl(ii,jj)    = 0.5*(vx_s_g_erg(i,j1)+vx_s_g_erg(i,j2))
   vy_s_g_dbl(ii,jj)    = 0.5*(vy_s_g_erg(i,j1)+vy_s_g_erg(i,j2))
   vz_s_dbl(ii,jj)      = 0.5*(vz_s_erg(i,j1)+vz_s_erg(i,j2))
   vh_s_dbl(ii,jj)      = 0.5*(vh_s_erg(i,j1)+vh_s_erg(i,j2))
   vx_m_g_dbl(ii,jj)    = 0.5*(vx_m_g_erg(i,j1)+vx_m_g_erg(i,j2))
   vy_m_g_dbl(ii,jj)    = 0.5*(vy_m_g_erg(i,j1)+vy_m_g_erg(i,j2))
   vh_m_dbl(ii,jj)      = 0.5*(vh_m_erg(i,j1)+vh_m_erg(i,j2))
   temp_b_dbl(ii,jj)    = 0.5*(temp_b_erg(i,j1)+temp_b_erg(i,j2))
   temph_b_dbl(ii,jj)   = 0.5*(temph_b_erg(i,j1)+temph_b_erg(i,j2))
   Q_bm_dbl(ii,jj)      = 0.5*(Q_bm_erg(i,j1)+Q_bm_erg(i,j2))
   Q_tld_dbl(ii,jj)     = 0.5*(Q_tld_erg(i,j1)+Q_tld_erg(i,j2))
   calving_dbl(ii,jj)   = 0.5*(calving_erg(i,j1)+calving_erg(i,j2))
   am_perp_dbl(ii,jj)   = 0.5*(am_perp_erg(i,j1)+am_perp_erg(i,j2))
   dzs_dtau_dbl(ii,jj)  = 0.5*(dzs_dtau_erg(i,j1)+dzs_dtau_erg(i,j2))
   dzm_dtau_dbl(ii,jj)  = 0.5*(dzm_dtau_erg(i,j1)+dzm_dtau_erg(i,j2))
   dzb_dtau_dbl(ii,jj)  = 0.5*(dzb_dtau_erg(i,j1)+dzb_dtau_erg(i,j2))
   dzl_dtau_dbl(ii,jj)  = 0.5*(dzl_dtau_erg(i,j1)+dzl_dtau_erg(i,j2))
   dH_c_dtau_dbl(ii,jj) = 0.5*(dH_c_dtau_erg(i,j1)+dH_c_dtau_erg(i,j2))
   dH_t_dtau_dbl(ii,jj) = 0.5*(dH_t_dtau_erg(i,j1)+dH_t_dtau_erg(i,j2))
   dH_dtau_dbl(ii,jj)   = 0.5*(dH_dtau_erg(i,j1)+dH_dtau_erg(i,j2))
   tau_dr_dbl(ii,jj)    = 0.5*(tau_dr_erg(i,j1)+tau_dr_erg(i,j2))
   tau_b_dbl(ii,jj)     = 0.5*(tau_b_erg(i,j1)+tau_b_erg(i,j2))
   p_b_w_dbl(ii,jj)     = 0.5*(p_b_w_erg(i,j1)+p_b_w_erg(i,j2))
   q_w_dbl(ii,jj)       = 0.5*(q_w_erg(i,j1)+q_w_erg(i,j2))
   q_w_x_dbl(ii,jj)     = 0.5*(q_w_x_erg(i,j1)+q_w_x_erg(i,j2))
   q_w_y_dbl(ii,jj)     = 0.5*(q_w_y_erg(i,j1)+q_w_y_erg(i,j2))
   H_w_dbl(ii,jj)       = 0.5*(H_w_erg(i,j1)+H_w_erg(i,j2))
   q_gl_g_dbl(ii,jj)    = 0.5*(q_gl_g_erg(i,j1)+q_gl_g_erg(i,j2))
   ratio_sl_sia_x_dbl(ii,jj) = 0.5*( ratio_sl_sia_x_erg(i,j1) &
                                    +ratio_sl_sia_x_erg(i,j2))
   ratio_sl_sia_y_dbl(ii,jj) = 0.5*( ratio_sl_sia_y_erg(i,j1) &
                                    +ratio_sl_sia_y_erg(i,j2))
   vis_ave_g_dbl(ii,jj) = 0.5*(vis_ave_g_erg(i,j1)+vis_ave_g_erg(i,j2))
   vis_int_g_dbl(ii,jj) = 0.5*(vis_int_g_erg(i,j1)+vis_int_g_erg(i,j2))
end do
end do

do ii = 1, 2*IMAX-1, 2
do jj = 1, 2*JMAX-1, 2
   i1 = (ii-1)/2
   i2 = (ii+1)/2
   j1 = (jj-1)/2
   j2 = (jj+1)/2
   lambda_dbl(ii,jj)    = 0.25*( lambda_erg(i1,j1)+lambda_erg(i2,j1) &
                                +lambda_erg(i1,j2)+lambda_erg(i2,j2) )
   phi_dbl(ii,jj)       = 0.25*( phi_erg(i1,j1)+phi_erg(i2,j1) &
                                +phi_erg(i1,j2)+phi_erg(i2,j2) )
   lon_dbl(ii,jj)       = 0.25*( lon_erg(i1,j1)+lon_erg(i2,j1) &
                                +lon_erg(i1,j2)+lon_erg(i2,j2) )
   lat_dbl(ii,jj)       = 0.25*( lat_erg(i1,j1)+lat_erg(i2,j1) &
                                +lat_erg(i1,j2)+lat_erg(i2,j2) )
   temp_s_dbl(ii,jj)    = 0.25*( temp_s_erg(i1,j1)+temp_s_erg(i2,j1) &
                                +temp_s_erg(i1,j2)+temp_s_erg(i2,j2) )
   prec_dbl(ii,jj)      = 0.25*( prec_erg(i1,j1)+prec_erg(i2,j1) &
                                +prec_erg(i1,j2)+prec_erg(i2,j2) )
   snowfall_dbl(ii,jj)  = 0.25*( snowfall_erg(i1,j1)+snowfall_erg(i2,j1) &
                                +snowfall_erg(i1,j2)+snowfall_erg(i2,j2) )
   rainfall_dbl(ii,jj)  = 0.25*( rainfall_erg(i1,j1)+rainfall_erg(i2,j1) &
                                +rainfall_erg(i1,j2)+rainfall_erg(i2,j2) )
   pdd_dbl(ii,jj)       = 0.25*( pdd_erg(i1,j1)+pdd_erg(i2,j1) &
                                +pdd_erg(i1,j2)+pdd_erg(i2,j2) )
   as_perp_dbl(ii,jj)   = 0.25*( as_perp_erg(i1,j1)+as_perp_erg(i2,j1) &
                                +as_perp_erg(i1,j2)+as_perp_erg(i2,j2) )
   as_perp_apl_dbl(ii,jj) = 0.25*( as_perp_apl_erg(i1,j1) &
                                  +as_perp_apl_erg(i2,j1) &
                                  +as_perp_apl_erg(i1,j2) &
                                  +as_perp_apl_erg(i2,j2) )
   smb_corr_dbl(ii,jj)  = 0.25*( smb_corr_erg(i1,j1)+smb_corr_erg(i2,j1) &
                                +smb_corr_erg(i1,j2)+smb_corr_erg(i2,j2) )
   z_sl_dbl(ii,jj)      = 0.25*( z_sl_erg(i1,j1)+z_sl_erg(i2,j1) &
                                +z_sl_erg(i1,j2)+z_sl_erg(i2,j2) )
#if (DISC>0)   /* Ice discharge parameterisation */
   dis_perp_dbl(ii,jj)    = 0.25*( dis_perp_erg(i1,j1) &
                                  +dis_perp_erg(i2,j1) &
                                  +dis_perp_erg(i1,j2) &
                                  +dis_perp_erg(i2,j2) )
   cst_dist_dbl(ii,jj)    = 0.25*( cst_dist_erg(i1,j1) &
                                  +cst_dist_erg(i2,j1) &
                                  +cst_dist_erg(i1,j2) &
                                  +cst_dist_erg(i2,j2) )
   cos_grad_tc_dbl(ii,jj) = 0.25*( cos_grad_tc_erg(i1,j1) &
                                  +cos_grad_tc_erg(i2,j1) &
                                  +cos_grad_tc_erg(i1,j2) &
                                  +cos_grad_tc_erg(i2,j2) )
#endif
   q_geo_dbl(ii,jj)     = 0.25*( q_geo_erg(i1,j1)+q_geo_erg(i2,j1) &
                                +q_geo_erg(i1,j2)+q_geo_erg(i2,j2) )
   zs_dbl(ii,jj)        = 0.25*( zs_erg(i1,j1)+zs_erg(i2,j1) &
                                +zs_erg(i1,j2)+zs_erg(i2,j2) )
   zm_dbl(ii,jj)        = 0.25*( zm_erg(i1,j1)+zm_erg(i2,j1) &
                                +zm_erg(i1,j2)+zm_erg(i2,j2) )
   zb_dbl(ii,jj)        = 0.25*( zb_erg(i1,j1)+zb_erg(i2,j1) &
                                +zb_erg(i1,j2)+zb_erg(i2,j2) )
   zl_dbl(ii,jj)        = 0.25*( zl_erg(i1,j1)+zl_erg(i2,j1) &
                                +zl_erg(i1,j2)+zl_erg(i2,j2) )
   zl0_dbl(ii,jj)       = 0.25*( zl0_erg(i1,j1)+zl0_erg(i2,j1) &
                                +zl0_erg(i1,j2)+zl0_erg(i2,j2) )
   wss_dbl(ii,jj)       = 0.25*( wss_erg(i1,j1)+wss_erg(i2,j1) &
                                +wss_erg(i1,j2)+wss_erg(i2,j2) )
   H_cold_dbl(ii,jj)    = 0.25*( H_cold_erg(i1,j1)+H_cold_erg(i2,j1) &
                                +H_cold_erg(i1,j2)+H_cold_erg(i2,j2) )
   H_temp_dbl(ii,jj)    = 0.25*( H_temp_erg(i1,j1)+H_temp_erg(i2,j1) &
                                +H_temp_erg(i1,j2)+H_temp_erg(i2,j2) )
   H_dbl(ii,jj)         = 0.25*( H_erg(i1,j1)+H_erg(i2,j1) &
                                +H_erg(i1,j2)+H_erg(i2,j2) )
   r_kc_cts_dbl(ii,jj)  = 0.25*( r_kc_cts_erg(i1,j1)+r_kc_cts_erg(i2,j1) &
                                +r_kc_cts_erg(i1,j2)+r_kc_cts_erg(i2,j2) )
   vz_c_dbl(ii,jj,:)    = 0.25*( vz_c_erg(i1,j1,:)+vz_c_erg(i2,j1,:) &
                                +vz_c_erg(i1,j2,:)+vz_c_erg(i2,j2,:) )
   vz_t_dbl(ii,jj,:)    = 0.25*( vz_t_erg(i1,j1,:)+vz_t_erg(i2,j1,:) &
                                +vz_t_erg(i1,j2,:)+vz_t_erg(i2,j2,:) )
   temp_c_dbl(ii,jj,:)  = 0.25*( temp_c_erg(i1,j1,:)+temp_c_erg(i2,j1,:) &
                                +temp_c_erg(i1,j2,:)+temp_c_erg(i2,j2,:) )
   age_c_dbl(ii,jj,:)   = 0.25*( age_c_erg(i1,j1,:)+age_c_erg(i2,j1,:) &
                                +age_c_erg(i1,j2,:)+age_c_erg(i2,j2,:) )
   omega_t_dbl(ii,jj,:) = 0.25*( omega_t_erg(i1,j1,:)+omega_t_erg(i2,j1,:) &
                                +omega_t_erg(i1,j2,:)+omega_t_erg(i2,j2,:) )
   age_t_dbl(ii,jj,:)   = 0.25*( age_t_erg(i1,j1,:)+age_t_erg(i2,j1,:) &
                                +age_t_erg(i1,j2,:)+age_t_erg(i2,j2,:) )
   temp_r_dbl(ii,jj,:)  = 0.25*( temp_r_erg(i1,j1,:)+temp_r_erg(i2,j1,:) &
                                +temp_r_erg(i1,j2,:)+temp_r_erg(i2,j2,:) )
   enth_c_dbl(ii,jj,:)  = 0.25*( enth_c_erg(i1,j1,:)+enth_c_erg(i2,j1,:) &
                                +enth_c_erg(i1,j2,:)+enth_c_erg(i2,j2,:) )
   enth_t_dbl(ii,jj,:)  = 0.25*( enth_t_erg(i1,j1,:)+enth_t_erg(i2,j1,:) &
                                +enth_t_erg(i1,j2,:)+enth_t_erg(i2,j2,:) )
   omega_c_dbl(ii,jj,:) = 0.25*( omega_c_erg(i1,j1,:)+omega_c_erg(i2,j1,:) &
                                +omega_c_erg(i1,j2,:)+omega_c_erg(i2,j2,:) )
   enh_c_dbl(ii,jj,:)   = 0.25*( enh_c_erg(i1,j1,:)+enh_c_erg(i2,j1,:) &
                                +enh_c_erg(i1,j2,:)+enh_c_erg(i2,j2,:) )
   enh_t_dbl(ii,jj,:)   = 0.25*( enh_t_erg(i1,j1,:)+enh_t_erg(i2,j1,:) &
                                +enh_t_erg(i1,j2,:)+enh_t_erg(i2,j2,:) )
   strain_heating_c_dbl(ii,jj,:) = 0.25*( strain_heating_c_erg(i1,j1,:) &
                                         +strain_heating_c_erg(i2,j1,:) &
                                         +strain_heating_c_erg(i1,j2,:) &
                                         +strain_heating_c_erg(i2,j2,:) )
   strain_heating_t_dbl(ii,jj,:) = 0.25*( strain_heating_t_erg(i1,j1,:) &
                                         +strain_heating_t_erg(i2,j1,:) &
                                         +strain_heating_t_erg(i1,j2,:) &
                                         +strain_heating_t_erg(i2,j2,:) )
   vx_b_g_dbl(ii,jj)    = 0.25*( vx_b_g_erg(i1,j1)+vx_b_g_erg(i2,j1) &
                                +vx_b_g_erg(i1,j2)+vx_b_g_erg(i2,j2) )
   vy_b_g_dbl(ii,jj)    = 0.25*( vy_b_g_erg(i1,j1)+vy_b_g_erg(i2,j1) &
                                +vy_b_g_erg(i1,j2)+vy_b_g_erg(i2,j2) )
   vz_b_dbl(ii,jj)      = 0.25*( vz_b_erg(i1,j1)+vz_b_erg(i2,j1) &
                                +vz_b_erg(i1,j2)+vz_b_erg(i2,j2) )
   vh_b_dbl(ii,jj)      = 0.25*( vh_b_erg(i1,j1)+vh_b_erg(i2,j1) &
                                +vh_b_erg(i1,j2)+vh_b_erg(i2,j2) )
   vx_s_g_dbl(ii,jj)    = 0.25*( vx_s_g_erg(i1,j1)+vx_s_g_erg(i2,j1) &
                                +vx_s_g_erg(i1,j2)+vx_s_g_erg(i2,j2) )
   vy_s_g_dbl(ii,jj)    = 0.25*( vy_s_g_erg(i1,j1)+vy_s_g_erg(i2,j1) &
                                +vy_s_g_erg(i1,j2)+vy_s_g_erg(i2,j2) )
   vz_s_dbl(ii,jj)      = 0.25*( vz_s_erg(i1,j1)+vz_s_erg(i2,j1) &
                                +vz_s_erg(i1,j2)+vz_s_erg(i2,j2) )
   vh_s_dbl(ii,jj)      = 0.25*( vh_s_erg(i1,j1)+vh_s_erg(i2,j1) &
                                +vh_s_erg(i1,j2)+vh_s_erg(i2,j2) )
   vx_m_g_dbl(ii,jj)    = 0.25*( vx_m_g_erg(i1,j1)+vx_m_g_erg(i2,j1) &
                                +vx_m_g_erg(i1,j2)+vx_m_g_erg(i2,j2) )
   vy_m_g_dbl(ii,jj)    = 0.25*( vy_m_g_erg(i1,j1)+vy_m_g_erg(i2,j1) &
                                +vy_m_g_erg(i1,j2)+vy_m_g_erg(i2,j2) )
   vh_m_dbl(ii,jj)      = 0.25*( vh_m_erg(i1,j1)+vh_m_erg(i2,j1) &
                                +vh_m_erg(i1,j2)+vh_m_erg(i2,j2) )
   temp_b_dbl(ii,jj)    = 0.25*( temp_b_erg(i1,j1)+temp_b_erg(i2,j1) &
                                +temp_b_erg(i1,j2)+temp_b_erg(i2,j2) )
   temph_b_dbl(ii,jj)   = 0.25*( temph_b_erg(i1,j1)+temph_b_erg(i2,j1) &
                                +temph_b_erg(i1,j2)+temph_b_erg(i2,j2) )
   Q_bm_dbl(ii,jj)      = 0.25*( Q_bm_erg(i1,j1)+Q_bm_erg(i2,j1) &
                                +Q_bm_erg(i1,j2)+Q_bm_erg(i2,j2) )
   Q_tld_dbl(ii,jj)     = 0.25*( Q_tld_erg(i1,j1)+Q_tld_erg(i2,j1) &
                                +Q_tld_erg(i1,j2)+Q_tld_erg(i2,j2) )
   calving_dbl(ii,jj)   = 0.25*( calving_erg(i1,j1)+calving_erg(i2,j1) &
                                +calving_erg(i1,j2)+calving_erg(i2,j2) )
   am_perp_dbl(ii,jj)   = 0.25*( am_perp_erg(i1,j1)+am_perp_erg(i2,j1) &
                                +am_perp_erg(i1,j2)+am_perp_erg(i2,j2) )
   dzs_dtau_dbl(ii,jj)  = 0.25*( dzs_dtau_erg(i1,j1)+dzs_dtau_erg(i2,j1) &
                                +dzs_dtau_erg(i1,j2)+dzs_dtau_erg(i2,j2) )
   dzm_dtau_dbl(ii,jj)  = 0.25*( dzm_dtau_erg(i1,j1)+dzm_dtau_erg(i2,j1) &
                                +dzm_dtau_erg(i1,j2)+dzm_dtau_erg(i2,j2) )
   dzb_dtau_dbl(ii,jj)  = 0.25*( dzb_dtau_erg(i1,j1)+dzb_dtau_erg(i2,j1) &
                                +dzb_dtau_erg(i1,j2)+dzb_dtau_erg(i2,j2) )
   dzl_dtau_dbl(ii,jj)  = 0.25*( dzl_dtau_erg(i1,j1)+dzl_dtau_erg(i2,j1) &
                                +dzl_dtau_erg(i1,j2)+dzl_dtau_erg(i2,j2) )
   dH_c_dtau_dbl(ii,jj) = 0.25*( dH_c_dtau_erg(i1,j1)+dH_c_dtau_erg(i2,j1) &
                                +dH_c_dtau_erg(i1,j2)+dH_c_dtau_erg(i2,j2) )
   dH_t_dtau_dbl(ii,jj) = 0.25*( dH_t_dtau_erg(i1,j1)+dH_t_dtau_erg(i2,j1) &
                                +dH_t_dtau_erg(i1,j2)+dH_t_dtau_erg(i2,j2) )
   dH_dtau_dbl(ii,jj)   = 0.25*( dH_dtau_erg(i1,j1)+dH_dtau_erg(i2,j1) &
                                +dH_dtau_erg(i1,j2)+dH_dtau_erg(i2,j2) )
   tau_dr_dbl(ii,jj)    = 0.25*( tau_dr_erg(i1,j1)+tau_dr_erg(i2,j1) &
                                +tau_dr_erg(i1,j2)+tau_dr_erg(i2,j2) )
   tau_b_dbl(ii,jj)     = 0.25*( tau_b_erg(i1,j1)+tau_b_erg(i2,j1) &
                                +tau_b_erg(i1,j2)+tau_b_erg(i2,j2) )
   p_b_w_dbl(ii,jj)     = 0.25*( p_b_w_erg(i1,j1)+p_b_w_erg(i2,j1) &
                                +p_b_w_erg(i1,j2)+p_b_w_erg(i2,j2) )
   q_w_dbl(ii,jj)       = 0.25*( q_w_erg(i1,j1)+q_w_erg(i2,j1) &
                                +q_w_erg(i1,j2)+q_w_erg(i2,j2) )
   q_w_x_dbl(ii,jj)     = 0.25*( q_w_x_erg(i1,j1)+q_w_x_erg(i2,j1) &
                                +q_w_x_erg(i1,j2)+q_w_x_erg(i2,j2) )
   q_w_y_dbl(ii,jj)     = 0.25*( q_w_y_erg(i1,j1)+q_w_y_erg(i2,j1) &
                                +q_w_y_erg(i1,j2)+q_w_y_erg(i2,j2) )
   H_w_dbl(ii,jj)       = 0.25*( H_w_erg(i1,j1)+H_w_erg(i2,j1) &
                                +H_w_erg(i1,j2)+H_w_erg(i2,j2) )
   q_gl_g_dbl(ii,jj)    = 0.25*( q_gl_g_erg(i1,j1)+q_gl_g_erg(i2,j1) &
                                +q_gl_g_erg(i1,j2)+q_gl_g_erg(i2,j2) )
   ratio_sl_sia_x_dbl(ii,jj) = 0.25*( ratio_sl_sia_x_erg(i1,j1) &
                                     +ratio_sl_sia_x_erg(i2,j1) &
                                     +ratio_sl_sia_x_erg(i1,j2) &
                                     +ratio_sl_sia_x_erg(i2,j2) )
   ratio_sl_sia_y_dbl(ii,jj) = 0.25*( ratio_sl_sia_y_erg(i1,j1) &
                                     +ratio_sl_sia_y_erg(i2,j1) &
                                     +ratio_sl_sia_y_erg(i1,j2) &
                                     +ratio_sl_sia_y_erg(i2,j2) )
   vis_ave_g_dbl(ii,jj) = 0.25*( vis_ave_g_erg(i1,j1)+vis_ave_g_erg(i2,j1) &
                                +vis_ave_g_erg(i1,j2)+vis_ave_g_erg(i2,j2) )
   vis_int_g_dbl(ii,jj) = 0.25*( vis_int_g_erg(i1,j1)+vis_int_g_erg(i2,j1) &
                                +vis_int_g_erg(i1,j2)+vis_int_g_erg(i2,j2) )
end do
end do

kc_cts_dbl = ceiling(r_kc_cts_dbl,i4b)

!-------- Interpolation for 2D and 3D fields
!                                  (staggered grid in x-direction) --------

qx_dbl(2*IMAX,:)       = 0.0_sp   ! outside domain -> undefined
vx_m_sia_dbl(2*IMAX,:) = 0.0_sp   ! outside domain -> undefined
vx_m_ssa_dbl(2*IMAX,:) = 0.0_sp   ! outside domain -> undefined
vx_c_dbl(2*IMAX,:,:)   = 0.0_sp   ! outside domain -> undefined
vx_t_dbl(2*IMAX,:,:)   = 0.0_sp   ! outside domain -> undefined

do jj = 0, 2*JMAX, 2

   j  = jj/2

   ii = 0
   i  = 0
   qx_dbl(ii,jj)     = qx_erg(i,j)
   vx_m_sia_dbl(ii,jj) = vx_m_sia_erg(i,j)
   vx_m_ssa_dbl(ii,jj) = vx_m_ssa_erg(i,j)
   vx_c_dbl(ii,jj,:) = vx_c_erg(i,j,:)
   vx_t_dbl(ii,jj,:) = vx_t_erg(i,j,:)

   ii = 2*IMAX-1
   i  = IMAX-1
   qx_dbl(ii,jj)     = qx_erg(i,j)
   vx_m_sia_dbl(ii,jj) = vx_m_sia_erg(i,j)
   vx_m_ssa_dbl(ii,jj) = vx_m_ssa_erg(i,j)
   vx_c_dbl(ii,jj,:) = vx_c_erg(i,j,:)
   vx_t_dbl(ii,jj,:) = vx_t_erg(i,j,:)

   do ii = 1, 2*IMAX-3, 2
      i1 = (ii-1)/2
      i2 = (ii+1)/2
      qx_dbl(ii,jj)     = 0.75*qx_erg(i1,j)+0.25*qx_erg(i2,j)
      vx_m_sia_dbl(ii,jj) = 0.75*vx_m_sia_erg(i1,j)+0.25*vx_m_sia_erg(i2,j)
      vx_m_ssa_dbl(ii,jj) = 0.75*vx_m_ssa_erg(i1,j)+0.25*vx_m_ssa_erg(i2,j)
      vx_c_dbl(ii,jj,:) = 0.75*vx_c_erg(i1,j,:)+0.25*vx_c_erg(i2,j,:)
      vx_t_dbl(ii,jj,:) = 0.75*vx_t_erg(i1,j,:)+0.25*vx_t_erg(i2,j,:)
   end do

   do ii = 2, 2*IMAX-2, 2
      i1 = (ii-2)/2
      i2 = ii/2
      qx_dbl(ii,jj)     = 0.25*qx_erg(i1,j)+0.75*qx_erg(i2,j)
      vx_m_sia_dbl(ii,jj) = 0.25*vx_m_sia_erg(i1,j)+0.75*vx_m_sia_erg(i2,j)
      vx_m_ssa_dbl(ii,jj) = 0.25*vx_m_ssa_erg(i1,j)+0.75*vx_m_ssa_erg(i2,j)
      vx_c_dbl(ii,jj,:) = 0.25*vx_c_erg(i1,j,:)+0.75*vx_c_erg(i2,j,:)
      vx_t_dbl(ii,jj,:) = 0.25*vx_t_erg(i1,j,:)+0.75*vx_t_erg(i2,j,:)
   end do

end do

do jj = 1, 2*JMAX-1, 2

   j1 = (jj-1)/2
   j2 = (jj+1)/2

   ii = 0
   i  = 0
   qx_dbl(ii,jj)     = 0.5*(qx_erg(i,j1)+qx_erg(i,j2))
   vx_m_sia_dbl(ii,jj) = 0.5*(vx_m_sia_erg(i,j1)+vx_m_sia_erg(i,j2))
   vx_m_ssa_dbl(ii,jj) = 0.5*(vx_m_ssa_erg(i,j1)+vx_m_ssa_erg(i,j2))
   vx_c_dbl(ii,jj,:) = 0.5*(vx_c_erg(i,j1,:)+vx_c_erg(i,j2,:))
   vx_t_dbl(ii,jj,:) = 0.5*(vx_t_erg(i,j1,:)+vx_t_erg(i,j2,:))

   ii = 2*IMAX-1
   i  = IMAX-1
   qx_dbl(ii,jj)     = 0.5*(qx_erg(i,j1)+qx_erg(i,j2))
   vx_m_sia_dbl(ii,jj) = 0.5*(vx_m_sia_erg(i,j1)+vx_m_sia_erg(i,j2))
   vx_m_ssa_dbl(ii,jj) = 0.5*(vx_m_ssa_erg(i,j1)+vx_m_ssa_erg(i,j2))
   vx_c_dbl(ii,jj,:) = 0.5*(vx_c_erg(i,j1,:)+vx_c_erg(i,j2,:))
   vx_t_dbl(ii,jj,:) = 0.5*(vx_t_erg(i,j1,:)+vx_t_erg(i,j2,:))

   do ii = 1, 2*IMAX-3, 2
      i1 = (ii-1)/2
      i2 = (ii+1)/2
      qx_dbl(ii,jj)     = 0.5*((0.75*qx_erg(i1,j1)+0.25*qx_erg(i2,j1)) &
                              +(0.75*qx_erg(i1,j2)+0.25*qx_erg(i2,j2)))
      vx_m_sia_dbl(ii,jj) = 0.5*((0.75*vx_m_sia_erg(i1,j1) &
                                       +0.25*vx_m_sia_erg(i2,j1)) &
                                +(0.75*vx_m_sia_erg(i1,j2) &
                                       +0.25*vx_m_sia_erg(i2,j2)))
      vx_m_ssa_dbl(ii,jj) = 0.5*((0.75*vx_m_ssa_erg(i1,j1) &
                                       +0.25*vx_m_ssa_erg(i2,j1)) &
                                +(0.75*vx_m_ssa_erg(i1,j2) &
                                       +0.25*vx_m_ssa_erg(i2,j2)))
      vx_c_dbl(ii,jj,:) = 0.5*((0.75*vx_c_erg(i1,j1,:)+0.25*vx_c_erg(i2,j1,:)) &
                              +(0.75*vx_c_erg(i1,j2,:)+0.25*vx_c_erg(i2,j2,:)))
      vx_t_dbl(ii,jj,:) = 0.5*((0.75*vx_t_erg(i1,j1,:)+0.25*vx_t_erg(i2,j1,:)) &
                              +(0.75*vx_t_erg(i1,j2,:)+0.25*vx_t_erg(i2,j2,:)))
   end do

   do ii = 2, 2*IMAX-2, 2
      i1 = (ii-2)/2
      i2 = ii/2
      qx_dbl(ii,jj)     = 0.5*((0.25*qx_erg(i1,j1)+0.75*qx_erg(i2,j1)) &
                              +(0.25*qx_erg(i1,j2)+0.75*qx_erg(i2,j2)))
      vx_m_sia_dbl(ii,jj) = 0.5*((0.25*vx_m_sia_erg(i1,j1) &
                                       +0.75*vx_m_sia_erg(i2,j1)) &
                                +(0.25*vx_m_sia_erg(i1,j2) &
                                       +0.75*vx_m_sia_erg(i2,j2)))
      vx_m_ssa_dbl(ii,jj) = 0.5*((0.25*vx_m_ssa_erg(i1,j1) &
                                       +0.75*vx_m_ssa_erg(i2,j1)) &
                                +(0.25*vx_m_ssa_erg(i1,j2) &
                                       +0.75*vx_m_ssa_erg(i2,j2)))
      vx_c_dbl(ii,jj,:) = 0.5*((0.25*vx_c_erg(i1,j1,:)+0.75*vx_c_erg(i2,j1,:)) &
                              +(0.25*vx_c_erg(i1,j2,:)+0.75*vx_c_erg(i2,j2,:)))
      vx_t_dbl(ii,jj,:) = 0.5*((0.25*vx_t_erg(i1,j1,:)+0.75*vx_t_erg(i2,j1,:)) &
                              +(0.25*vx_t_erg(i1,j2,:)+0.75*vx_t_erg(i2,j2,:)))
   end do

end do

!-------- Interpolation for 2D and 3D fields
!                                   (staggered grid in y-direction) --------

qy_dbl(:,2*JMAX)       = 0.0_sp   ! outside domain -> undefined
vy_m_sia_dbl(:,2*JMAX) = 0.0_sp   ! outside domain -> undefined
vy_m_ssa_dbl(:,2*JMAX) = 0.0_sp   ! outside domain -> undefined
vy_c_dbl(:,2*JMAX,:)   = 0.0_sp   ! outside domain -> undefined
vy_t_dbl(:,2*JMAX,:)   = 0.0_sp   ! outside domain -> undefined

do ii = 0, 2*IMAX, 2

   i  = ii/2

   jj = 0
   j  = 0
   qy_dbl(ii,jj)     = qy_erg(i,j)
   vy_m_sia_dbl(ii,jj) = vy_m_sia_erg(i,j)
   vy_m_ssa_dbl(ii,jj) = vy_m_ssa_erg(i,j)
   vy_c_dbl(ii,jj,:) = vy_c_erg(i,j,:)
   vy_t_dbl(ii,jj,:) = vy_t_erg(i,j,:)

   jj = 2*JMAX-1
   j  = JMAX-1
   qy_dbl(ii,jj)     = qy_erg(i,j)
   vy_m_sia_dbl(ii,jj) = vy_m_sia_erg(i,j)
   vy_m_ssa_dbl(ii,jj) = vy_m_ssa_erg(i,j)
   vy_c_dbl(ii,jj,:) = vy_c_erg(i,j,:)
   vy_t_dbl(ii,jj,:) = vy_t_erg(i,j,:)

   do jj = 1, 2*JMAX-3, 2
      j1 = (jj-1)/2
      j2 = (jj+1)/2
      qy_dbl(ii,jj)     = 0.75*qy_erg(i,j1)+0.25*qy_erg(i,j2)
      vy_m_sia_dbl(ii,jj) = 0.75*vy_m_sia_erg(i,j1)+0.25*vy_m_sia_erg(i,j2)
      vy_m_ssa_dbl(ii,jj) = 0.75*vy_m_ssa_erg(i,j1)+0.25*vy_m_ssa_erg(i,j2)
      vy_c_dbl(ii,jj,:) = 0.75*vy_c_erg(i,j1,:)+0.25*vy_c_erg(i,j2,:)
      vy_t_dbl(ii,jj,:) = 0.75*vy_t_erg(i,j1,:)+0.25*vy_t_erg(i,j2,:)
   end do

   do jj = 2, 2*JMAX-2, 2
      j1 = (jj-2)/2
      j2 = jj/2
      qy_dbl(ii,jj)     = 0.25*qy_erg(i,j1)+0.75*qy_erg(i,j2)
      vy_m_sia_dbl(ii,jj) = 0.25*vy_m_sia_erg(i,j1)+0.75*vy_m_sia_erg(i,j2)
      vy_m_ssa_dbl(ii,jj) = 0.25*vy_m_ssa_erg(i,j1)+0.75*vy_m_ssa_erg(i,j2)
      vy_c_dbl(ii,jj,:) = 0.25*vy_c_erg(i,j1,:)+0.75*vy_c_erg(i,j2,:)
      vy_t_dbl(ii,jj,:) = 0.25*vy_t_erg(i,j1,:)+0.75*vy_t_erg(i,j2,:)
   end do

end do

do ii = 1, 2*IMAX-1, 2

   i1 = (ii-1)/2
   i2 = (ii+1)/2

   jj = 0
   j  = 0
   qy_dbl(ii,jj)     = 0.5*(qy_erg(i1,j)+qy_erg(i2,j))
   vy_m_sia_dbl(ii,jj) = 0.5*(vy_m_sia_erg(i1,j)+vy_m_sia_erg(i2,j))
   vy_m_ssa_dbl(ii,jj) = 0.5*(vy_m_ssa_erg(i1,j)+vy_m_ssa_erg(i2,j))
   vy_c_dbl(ii,jj,:) = 0.5*(vy_c_erg(i1,j,:)+vy_c_erg(i2,j,:))
   vy_t_dbl(ii,jj,:) = 0.5*(vy_t_erg(i1,j,:)+vy_t_erg(i2,j,:))

   jj = 2*JMAX-1
   j  = JMAX-1
   qy_dbl(ii,jj)     = 0.5*(qy_erg(i1,j)+qy_erg(i2,j))
   vy_m_sia_dbl(ii,jj) = 0.5*(vy_m_sia_erg(i1,j)+vy_m_sia_erg(i2,j))
   vy_m_ssa_dbl(ii,jj) = 0.5*(vy_m_ssa_erg(i1,j)+vy_m_ssa_erg(i2,j))
   vy_c_dbl(ii,jj,:) = 0.5*(vy_c_erg(i1,j,:)+vy_c_erg(i2,j,:))
   vy_t_dbl(ii,jj,:) = 0.5*(vy_t_erg(i1,j,:)+vy_t_erg(i2,j,:))

   do jj = 1, 2*JMAX-3, 2
      j1 = (jj-1)/2
      j2 = (jj+1)/2
      qy_dbl(ii,jj)     = 0.5*((0.75*qy_erg(i1,j1)+0.25*qy_erg(i1,j2)) &
                              +(0.75*qy_erg(i2,j1)+0.25*qy_erg(i2,j2)))
      vy_m_sia_dbl(ii,jj) = 0.5*((0.75*vy_m_sia_erg(i1,j1) &
                                       +0.25*vy_m_sia_erg(i1,j2)) &
                                +(0.75*vy_m_sia_erg(i2,j1) &
                                       +0.25*vy_m_sia_erg(i2,j2)))
      vy_m_ssa_dbl(ii,jj) = 0.5*((0.75*vy_m_ssa_erg(i1,j1) &
                                       +0.25*vy_m_ssa_erg(i1,j2)) &
                                +(0.75*vy_m_ssa_erg(i2,j1) &
                                       +0.25*vy_m_ssa_erg(i2,j2)))
      vy_c_dbl(ii,jj,:) = 0.5*((0.75*vy_c_erg(i1,j1,:)+0.25*vy_c_erg(i1,j2,:)) &
                              +(0.75*vy_c_erg(i2,j1,:)+0.25*vy_c_erg(i2,j2,:)))
      vy_t_dbl(ii,jj,:) = 0.5*((0.75*vy_t_erg(i1,j1,:)+0.25*vy_t_erg(i1,j2,:)) &
                              +(0.75*vy_t_erg(i2,j1,:)+0.25*vy_t_erg(i2,j2,:)))
   end do

   do jj = 2, 2*JMAX-2, 2
      j1 = (jj-2)/2
      j2 = jj/2
      qy_dbl(ii,jj)     = 0.5*((0.25*qy_erg(i1,j1)+0.75*qy_erg(i1,j2)) &
                              +(0.25*qy_erg(i2,j1)+0.75*qy_erg(i2,j2)))
      vy_m_sia_dbl(ii,jj) = 0.5*((0.25*vy_m_sia_erg(i1,j1) &
                                       +0.75*vy_m_sia_erg(i1,j2)) &
                                +(0.25*vy_m_sia_erg(i2,j1) &
                                       +0.75*vy_m_sia_erg(i2,j2)))
      vy_m_ssa_dbl(ii,jj) = 0.5*((0.25*vy_m_ssa_erg(i1,j1) &
                                       +0.75*vy_m_ssa_erg(i1,j2)) &
                                +(0.25*vy_m_ssa_erg(i2,j1) &
                                       +0.75*vy_m_ssa_erg(i2,j2)))
      vy_c_dbl(ii,jj,:) = 0.5*((0.25*vy_c_erg(i1,j1,:)+0.75*vy_c_erg(i1,j2,:)) &
                              +(0.25*vy_c_erg(i2,j1,:)+0.75*vy_c_erg(i2,j2,:)))
      vy_t_dbl(ii,jj,:) = 0.5*((0.25*vy_t_erg(i1,j1,:)+0.75*vy_t_erg(i1,j2,:)) &
                              +(0.25*vy_t_erg(i2,j1,:)+0.75*vy_t_erg(i2,j2,:)))
   end do

end do

!-------- Masks --------

do ii = 0, 2*IMAX
do jj = 0, 2*JMAX

   if ((mod(ii,2) == 0).and.(mod(jj,2) == 0)) then
      ! ii even, jj even

      i  = ii/2
      j  = jj/2

      mask_dbl(ii,jj)     = mask_erg(i,j)
      mask_old_dbl(ii,jj) = mask_old_erg(i,j)
      n_cts_dbl(ii,jj)     = n_cts_erg(i,j)

#if (DISC>0)   /* Ice discharge parameterisation */
      mask_mar_dbl(ii,jj)  = mask_mar_erg(i,j)
#endif

   else if ((mod(ii,2) /= 0).and.(mod(jj,2) == 0)) then
      ! ii odd, jj even

      i1 = (ii-1)/2
      i2 = (ii+1)/2
      j  = jj/2

      if (    (mask_erg(i1,j)==3) &
          .or.(mask_erg(i2,j)==3) ) &
      then
         mask_dbl(ii,jj) = 3
      else if (    (mask_erg(i1,j)==0) &
               .or.(mask_erg(i2,j)==0) ) &
      then
         mask_dbl(ii,jj) = 0
      else if (    (mask_erg(i1,j)==1) &
               .or.(mask_erg(i2,j)==1) ) &
      then
         mask_dbl(ii,jj) = 1
      else
         mask_dbl(ii,jj) = 2
      end if

      if (    (mask_old_erg(i1,j)==3) &
          .or.(mask_old_erg(i2,j)==3) ) &
      then
         mask_old_dbl(ii,jj) = 3
      else if (    (mask_old_erg(i1,j)==0) &
               .or.(mask_old_erg(i2,j)==0) ) &
      then
         mask_old_dbl(ii,jj) = 0
      else if (    (mask_old_erg(i1,j)==1) &
               .or.(mask_old_erg(i2,j)==1) ) &
      then
         mask_old_dbl(ii,jj) = 1
      else
         mask_old_dbl(ii,jj) = 2
      end if

      if (    (n_cts_erg(i1,j)==1) &
          .or.(n_cts_erg(i2,j)==1) ) &
      then
         n_cts_dbl(ii,jj) = 1
      else if (    (n_cts_erg(i1,j)==0) &
               .or.(n_cts_erg(i2,j)==0) ) &
      then
         n_cts_dbl(ii,jj) = 0
      else
         n_cts_dbl(ii,jj) = -1
      end if

#if (DISC>0)   /* Ice discharge parameterisation */

      if ( (mask_mar_erg(i1,j)+mask_mar_erg(i2,j)) >= 1 ) then
         mask_mar_dbl(ii,jj) = 1
      else
         mask_mar_dbl(ii,jj) = 0
      end if

#endif

   else if ((mod(ii,2) == 0).and.(mod(jj,2) /= 0)) then
      ! ii even, jj odd

      i  = ii/2
      j1 = (jj-1)/2
      j2 = (jj+1)/2

      if (    (mask_erg(i,j1)==3) &
          .or.(mask_erg(i,j2)==3) ) &
      then
         mask_dbl(ii,jj) = 3
      else if (    (mask_erg(i,j1)==0) &
               .or.(mask_erg(i,j2)==0) ) &
      then
         mask_dbl(ii,jj) = 0
      else if (    (mask_erg(i,j1)==1) &
               .or.(mask_erg(i,j2)==1) ) &
      then
         mask_dbl(ii,jj) = 1
      else
         mask_dbl(ii,jj) = 2
      end if

      if (    (mask_old_erg(i,j1)==3) &
          .or.(mask_old_erg(i,j2)==3) ) &
      then
         mask_old_dbl(ii,jj) = 3
      else if (    (mask_old_erg(i,j1)==0) &
               .or.(mask_old_erg(i,j2)==0) ) &
      then
         mask_old_dbl(ii,jj) = 0
      else if (    (mask_old_erg(i,j1)==1) &
               .or.(mask_old_erg(i,j2)==1) ) &
      then
         mask_old_dbl(ii,jj) = 1
      else
         mask_old_dbl(ii,jj) = 2
      end if

      if (    (n_cts_erg(i,j1)==1) &
          .or.(n_cts_erg(i,j2)==1) ) &
      then
         n_cts_dbl(ii,jj) = 1
      else if (    (n_cts_erg(i,j1)==0) &
               .or.(n_cts_erg(i,j2)==0) ) &
      then
         n_cts_dbl(ii,jj) = 0
      else
         n_cts_dbl(ii,jj) = -1
      end if

#if (DISC>0)   /* Ice discharge parameterisation */

      if ( (mask_mar_erg(i,j1)+mask_mar_erg(i,j2)) >= 1 ) then
         mask_mar_dbl(ii,jj) = 1
      else
         mask_mar_dbl(ii,jj) = 0
      end if

#endif

   else   ! (mod(ii,2) /= 0).and.(mod(jj,2) /= 0)
      ! ii odd, jj odd

      i1 = (ii-1)/2
      i2 = (ii+1)/2
      j1 = (jj-1)/2
      j2 = (jj+1)/2

      if (    (mask_erg(i1,j1)==3) &
          .or.(mask_erg(i2,j1)==3) &
          .or.(mask_erg(i1,j2)==3) &
          .or.(mask_erg(i2,j2)==3) ) &
      then
         mask_dbl(ii,jj) = 3
      else if (    (mask_erg(i1,j1)==0) &
               .or.(mask_erg(i2,j1)==0) &
               .or.(mask_erg(i1,j2)==0) &
               .or.(mask_erg(i2,j2)==0) ) &
      then
         mask_dbl(ii,jj) = 0
      else if (    (mask_erg(i1,j1)==1) &
               .or.(mask_erg(i2,j1)==1) &
               .or.(mask_erg(i1,j2)==1) &
               .or.(mask_erg(i2,j2)==1) ) &
      then
         mask_dbl(ii,jj) = 1
      else
         mask_dbl(ii,jj) = 2
      end if

      if (    (mask_old_erg(i1,j1)==3) &
          .or.(mask_old_erg(i2,j1)==3) &
          .or.(mask_old_erg(i1,j2)==3) &
          .or.(mask_old_erg(i2,j2)==3) ) &
      then
         mask_old_dbl(ii,jj) = 3
      else if (    (mask_old_erg(i1,j1)==0) &
               .or.(mask_old_erg(i2,j1)==0) &
               .or.(mask_old_erg(i1,j2)==0) &
               .or.(mask_old_erg(i2,j2)==0) ) &
      then
         mask_old_dbl(ii,jj) = 0
      else if (    (mask_old_erg(i1,j1)==1) &
               .or.(mask_old_erg(i2,j1)==1) &
               .or.(mask_old_erg(i1,j2)==1) &
               .or.(mask_old_erg(i2,j2)==1) ) &
      then
         mask_old_dbl(ii,jj) = 1
      else
         mask_old_dbl(ii,jj) = 2
      end if

      if (    (n_cts_erg(i1,j1)==1) &
          .or.(n_cts_erg(i2,j1)==1) &
          .or.(n_cts_erg(i1,j2)==1) &
          .or.(n_cts_erg(i2,j2)==1) ) &
      then
         n_cts_dbl(ii,jj) = 1
      else if (    (n_cts_erg(i1,j1)==0) &
               .or.(n_cts_erg(i2,j1)==0) &
               .or.(n_cts_erg(i1,j2)==0) &
               .or.(n_cts_erg(i2,j2)==0) ) &
      then
         n_cts_dbl(ii,jj) = 0
      else
         n_cts_dbl(ii,jj) = -1
      end if

#if (DISC>0)   /* Ice discharge parameterisation */

      if ( (mask_mar_erg(i1,j1)+mask_mar_erg(i2,j1) &
           +mask_mar_erg(i1,j2)+mask_mar_erg(i2,j2)) >= 2 ) then
         mask_mar_dbl(ii,jj) = 1
      else
         mask_mar_dbl(ii,jj) = 0
      end if

#endif

   end if

end do
end do

!-------- Flags --------

flag_shelfy_stream_x_dbl    = 0   ! all
flag_shelfy_stream_y_dbl    = 0   ! set
flag_shelfy_stream_dbl      = 0   ! to
flag_grounding_line_1_dbl   = 0   ! 0
flag_grounding_line_2_dbl   = 0   ! (false),
flag_calving_front_1_dbl    = 0   ! will
flag_calving_front_2_dbl    = 0   ! be
flag_grounded_front_a_1_dbl = 0   ! re-set
flag_grounded_front_a_2_dbl = 0   ! by
flag_grounded_front_b_1_dbl = 0   ! SICOPOLIS
flag_grounded_front_b_2_dbl = 0   ! anyway

end subroutine double_res_interpol

!-------------------------------------------------------------------------------
!> Writing of double resolution data on NetCDF file.
!-------------------------------------------------------------------------------
subroutine write_nc_double(run_name, ergnum, forcing_flag)

use resolution_doubler_types
use resolution_doubler_vars
use netcdf

implicit none

integer(i4b),       intent(in) :: forcing_flag
character(len=256), intent(in) :: run_name
character(len=  4), intent(in) :: ergnum

integer(i4b) :: ios
integer(i4b) :: cmode
integer(i4b) :: n_deflate_level
logical      :: flag_shuffle
integer(i4b) :: ncid, ncv
!     ncid:      ID of the output file
!     ncv:       Variable ID
integer(i4b) :: ncd, nc1d, nc2d(2), nc3d(3)
!     ncd:       Dimension ID
!     nc1d:      Dimension of a 1D array
!     nc2d:      Vector with the dimensions of a 2D array
!     nc3d:      Vector with the dimensions of a 3D array
integer(i4b) :: nc2flag(2), nc3flag(3), nc4flag(4)
!     nc2flag:   Vector with the 2 possible values of a flag variable
!     nc3flag:   Vector with the 3 possible values of a flag variable
!     nc4flag:   Vector with the 4 possible values of a flag variable
integer(i4b) :: nc1cor_i(1), nc1cor_j(1), &
                nc1cor_kc(1), nc1cor_kt(1), nc1cor_kr(1), &
                nc2cor_ij(2), &
                nc3cor_ijkc(3), nc3cor_ijkt(3), nc3cor_ijkr(3)
!     nc1cor(1): Corner of a 1D array
!     nc2cor(2): Corner of a 2D array
!     nc3cor(3): Corner of a 3D array
integer(i4b) :: nc1cnt_i(1), nc1cnt_j(1), &
                nc1cnt_kc(1), nc1cnt_kt(1), nc1cnt_kr(1), &
                nc2cnt_ij(2), &
                nc3cnt_ijkc(3), nc3cnt_ijkt(3), nc3cnt_ijkr(3)
!     nc1cnt(1): Count of a 1D array
!     nc2cnt(2): Count of a 2D array
!     nc3cnt(3): Count of a 3D array
character(len= 16) :: ch_date, ch_time, ch_zone
character(len= 16), allocatable :: coord_id(:)
character(len=256) :: filename_with_path, buffer
character(len=256) :: ch_msg
character, parameter :: end_of_line = char(10)
                        ! End-of-line string

nc1cor_i = (/ 1 /)
nc1cnt_i = (/ 2*IMAX+1 /)

nc1cor_j = (/ 1 /)
nc1cnt_j = (/ 2*JMAX+1 /)

nc1cor_kc = (/ 1 /)
nc1cnt_kc = (/ KCMAX+1 /)

nc1cor_kt = (/ 1 /)
nc1cnt_kt = (/ KTMAX+1 /)

nc1cor_kr = (/ 1 /)
nc1cnt_kr = (/ KRMAX+1 /)

nc2cor_ij = (/ 1, 1 /)
nc2cnt_ij = (/ 2*IMAX+1, 2*JMAX+1 /)

nc3cor_ijkc = (/ 1, 1, 1 /)
nc3cnt_ijkc = (/ 2*IMAX+1, 2*JMAX+1, KCMAX+1 /)

nc3cor_ijkt = (/ 1, 1, 1 /)
nc3cnt_ijkt = (/ 2*IMAX+1, 2*JMAX+1, KTMAX+1 /)

nc3cor_ijkr = (/ 1, 1, 1 /)
nc3cnt_ijkr = (/ 2*IMAX+1, 2*JMAX+1, KRMAX+1 /)

if (allocated(coord_id)) deallocate(coord_id); allocate(coord_id(5))
coord_id(1) = 'x'; coord_id(2) = 'y'
coord_id(3) = 'zeta_c'; coord_id(4) = 'zeta_t'; coord_id(5) = 'zeta_r'

#if (NETCDF4_ENABLED==1)
  cmode           = ior(NF90_NETCDF4, NF90_CLASSIC_MODEL)
  n_deflate_level = 1
  flag_shuffle    = .true.
#else
  cmode           = NF90_NOCLOBBER
  n_deflate_level = 0
  flag_shuffle    = .false.
#endif

!-------- NetCDF initialization --------

!  ------ Open NetCDF file

write (6,'(/a)') ' Now opening new NetCDF file ...'

filename_with_path = trim(OUT_PATH)//'/'//trim(run_name)//'_dbl_'// &
                     trim(ergnum)//'.nc'

ios = nf90_create(trim(filename_with_path), cmode, ncid)

if (ios /= nf90_noerr) then
   ch_msg = ' >>> write_nc_double:' &
          //         end_of_line &
          //'        Error when opening the resolution-doubled nc file!'
   call error(ch_msg)
end if

!  ------ Global attributes

buffer = 'Time-slice output no. '//trim(ergnum)//' '// &
         'of simulation '//trim(run_name)//' (doubled horizontal resolution)'
call check( nf90_put_att(ncid, NF90_GLOBAL, 'title', trim(buffer)) )

call set_ch_institution(buffer)
call check( nf90_put_att(ncid, NF90_GLOBAL, 'institution', trim(buffer)) )

buffer = 'SICOPOLIS'
call check( nf90_put_att(ncid, NF90_GLOBAL, 'source', trim(buffer)) )

call date_and_time(ch_date, ch_time, ch_zone)
buffer = ch_date(1:4)//'-'//ch_date(5:6)//'-'//ch_date(7:8)//' '// &
         ch_time(1:2)//':'//ch_time(3:4)//':'//ch_time(5:6)//' '// &
         ch_zone(1:3)//':'//ch_zone(4:5)//' - Data produced'
call check( nf90_put_att(ncid, NF90_GLOBAL, 'history', trim(buffer)) )

call set_ch_website(buffer)
call check( nf90_put_att(ncid, NF90_GLOBAL, 'references', trim(buffer)) )

!  ------ Definition of the dimensions (doubled resolution!)

call check( nf90_def_dim(ncid, trim(coord_id(1)), 2*IMAX+1, ncd) )
call check( nf90_def_dim(ncid, trim(coord_id(2)), 2*JMAX+1, ncd) )
call check( nf90_def_dim(ncid, trim(coord_id(3)), KCMAX+1,  ncd) )
call check( nf90_def_dim(ncid, trim(coord_id(4)), KTMAX+1,  ncd) )
call check( nf90_def_dim(ncid, trim(coord_id(5)), KRMAX+1,  ncd) )

!  ------ Definition of the variables

!    ---- mapping

call check( nf90_def_var(ncid, 'mapping', NF90_BYTE, ncv) )

call check( nf90_put_att(ncid, ncv, 'grid_mapping_name', &
                                     trim(mapping_grid_mapping_name_dbl)) )

if (trim(mapping_ellipsoid_dbl) /= 'xxx' ) &
   call check( nf90_put_att(ncid, ncv, 'ellipsoid', &
                                        trim(mapping_ellipsoid_dbl)) )

if (mapping_semi_major_axis_dbl > (no_value_neg_dp+eps_dp) ) &
   call check( nf90_put_att(ncid, ncv, 'semi_major_axis', &
                                        mapping_semi_major_axis_dbl) )

if (mapping_inv_flattening_dbl > (no_value_neg_dp+eps_dp) ) &
   call check( nf90_put_att(ncid, ncv, 'inverse_flattening', &
                                        mapping_inv_flattening_dbl) )

if (mapping_radius_of_sphere_dbl > (no_value_neg_dp+eps_dp) ) &
   call check( nf90_put_att(ncid, ncv, 'radius_of_sphere', &
                                        mapping_radius_of_sphere_dbl) )

call check( nf90_put_att(ncid, ncv, 'latitude_of_projection_origin', &
                                     mapping_latitude_origin_dbl) )
call check( nf90_put_att(ncid, ncv, 'standard_parallel', &
                                     mapping_standard_parallel_dbl) )
call check( nf90_put_att(ncid, ncv, 'straight_vertical_longitude_from_pole', &
                                     mapping_reference_longitude_dbl) )
call check( nf90_put_att(ncid, ncv, 'false_easting', mapping_false_E_dbl) )
call check( nf90_put_att(ncid, ncv, 'false_northing', mapping_false_N_dbl) )

!    ---- year2sec

call check( nf90_def_var(ncid, 'year2sec', NF90_DOUBLE, ncv) )

buffer = 's a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'seconds_per_year'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = '1 year (1 a) in seconds'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

!    ---- time

call check( nf90_def_var(ncid, 'time', NF90_DOUBLE, ncv) )

buffer = 'a'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'time'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Time'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

if (forcing_flag == 1) then

!    ---- delta_ts

   call check( nf90_def_var(ncid, 'delta_ts', NF90_DOUBLE, ncv) )

   buffer = 'degC'
   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
   buffer = 'surface_temperature_anomaly'
   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
   buffer = 'Surface temperature anomaly'
   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

else if (forcing_flag == 2) then

!    ---- glac_index

   call check( nf90_def_var(ncid, 'glac_index', NF90_DOUBLE, ncv) )

   buffer = '1'
   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
   buffer = 'glacial_index'
   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
   buffer = 'Glacial index'
   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

end if

!    ---- z_sl_mean

call check( nf90_def_var(ncid, 'z_sl_mean', NF90_DOUBLE, ncv) )

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'global_average_sea_level_change'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Mean sea level'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

!    ---- V_tot

call check( nf90_def_var(ncid, 'V_tot', NF90_DOUBLE, ncv) )

buffer = 'm3'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_volume'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Ice volume'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

!    ---- V_af

call check( nf90_def_var(ncid, 'V_af', NF90_DOUBLE, ncv) )

buffer = 'm3'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_volume_not_displacing_sea_water'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Ice volume above flotation'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

!    ---- A_grounded

call check( nf90_def_var(ncid, 'A_grounded', NF90_DOUBLE, ncv) )

buffer = 'm2'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'grounded_land_ice_area'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Area covered by grounded ice'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

!    ---- A_floating

call check( nf90_def_var(ncid, 'A_floating', NF90_DOUBLE, ncv) )

buffer = 'm2'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'floating_ice_shelf_area'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Area covered by floating ice'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

!    ---- x (= xi)

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc1d) )

call check( nf90_def_var(ncid, 'x', NF90_DOUBLE, nc1d, ncv) )

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'projection_x_coordinate'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'x-coordinate of the grid point i'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'axis', 'x') )

!    ---- y (= eta)

call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc1d) )

call check( nf90_def_var(ncid, 'y', NF90_DOUBLE, nc1d, ncv) )

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'projection_y_coordinate'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'y-coordinate of the grid point j'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'axis', 'y') )

!    ---- sigma_level_c

call check( nf90_inq_dimid(ncid, trim(coord_id(3)), nc1d) )

call check( nf90_def_var(ncid, 'sigma_level_c', NF90_DOUBLE, nc1d, ncv) )

buffer = 'up'
call check( nf90_put_att(ncid, ncv, 'positive', trim(buffer)) )
buffer = 'land_ice_kc_layer_sigma_coordinate'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'sigma-coordinate of the grid point kc'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

!    ---- sigma_level_t

call check( nf90_inq_dimid(ncid, trim(coord_id(4)), nc1d) )

call check( nf90_def_var(ncid, 'sigma_level_t', NF90_DOUBLE, nc1d, ncv) )

buffer = 'up'
call check( nf90_put_att(ncid, ncv, 'positive', trim(buffer)) )
buffer = 'land_ice_kt_layer_sigma_coordinate'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'sigma-coordinate of the grid point kt'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

!    ---- sigma_level_r

call check( nf90_inq_dimid(ncid, trim(coord_id(5)), nc1d) )

call check( nf90_def_var(ncid, 'sigma_level_r', NF90_DOUBLE, nc1d, ncv) )

buffer = 'up'
call check( nf90_put_att(ncid, ncv, 'positive', trim(buffer)) )
buffer = 'lithosphere_layer_sigma_coordinate'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'sigma-coordinate of the grid point kr'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

!    ---- lon

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'lon', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'lon', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'degrees_E'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'longitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Geographical longitude'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

!    ---- lat

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'lat', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'lat', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'degrees_N'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'latitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Geographical latitude'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

!    ---- lambda

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'lambda', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'lambda', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'rad'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'longitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Geographical longitude'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- phi

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'phi', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'phi', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'rad'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'latitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Geographical latitude'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- temp_s

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'temp_s', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'temp_s', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'degC'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'surface_temperature'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Temperature at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- prec

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'prec', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'prec', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_precipitation'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Annual precipitation at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- snowfall

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'snowfall', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'snowfall', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_snowfall'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Annual snowfall at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- rainfall

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'rainfall', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'rainfall', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_rainfall'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Annual rainfall at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- pdd

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'pdd', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'pdd', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'degC a'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_positive_degree_days'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Positive degree days at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- as_perp

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'as_perp', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'as_perp', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_surface_mass_balance'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Mass balance at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- as_perp_apl

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'as_perp_apl', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'as_perp_apl', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'applied_land_ice_surface_mass_balance'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Applied mass balance at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- smb_corr

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'smb_corr', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'smb_corr', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_surface_mass_balance_diagnosed_correction'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Diagnosed correction of the mass balance at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- z_sl

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'z_sl', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'z_sl', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'sea_level_change'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Sea level'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

#if (DISC>0)   /* Ice discharge parameterisation */

!    ---- dis_perp

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'dis_perp', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'dis_perp', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'ice_discharge'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Ice discharge'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- cst_dist

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'cst_dist', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'cst_dist', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'km'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'coastal_distance'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Coastal distance'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- cos_grad_tc

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'cos_grad_tc', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'cos_grad_tc', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = '1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'cos_alpha'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Cosine of angle between surface gradient and cst dist gradient'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- mask_mar

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'mask_mar', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'mask_mar', NF90_BYTE, nc2d, ncv) )
#endif

buffer = 'marginal_ring_mask'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Marginal ring mask'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
nc2flag = (/ 0, 1 /)
call check( nf90_put_att(ncid, ncv, 'flag_values', nc2flag) )
buffer = 'no_ring '// &
         'ring'
call check( nf90_put_att(ncid, ncv, 'flag_meanings', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

#endif

!    ---- q_geo

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'q_geo', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'q_geo', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'W m-2'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'upward_geothermal_heat_flux_at_ground_level'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Geothermal heat flux'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- mask

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'mask', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'mask', NF90_BYTE, nc2d, ncv) )
#endif

buffer = 'ice_land_sea_mask'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Ice-land-sea mask'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
nc4flag = (/ 0, 1, 2, 3 /)
call check( nf90_put_att(ncid, ncv, 'flag_values', nc4flag) )
buffer = 'glaciated_land '// &
         'ice_free_land '// &
         'sea '// &
         'floating_ice'
call check( nf90_put_att(ncid, ncv, 'flag_meanings', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- mask_old

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'mask_old', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'mask_old', NF90_BYTE, nc2d, ncv) )
#endif

buffer = 'ice_land_sea_mask_old'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Ice-land-sea mask (old)'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
nc4flag = (/ 0, 1, 2, 3 /)
call check( nf90_put_att(ncid, ncv, 'flag_values', nc4flag) )
buffer = 'glaciated_land '// &
         'ice_free_land '// &
         'sea '// &
         'floating_ice'
call check( nf90_put_att(ncid, ncv, 'flag_meanings', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- n_cts

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'n_cts', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'n_cts', NF90_BYTE, nc2d, ncv) )
#endif

buffer = 'polythermal_condition_mask'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Mask for polythermal conditions'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
nc3flag = (/ -1, 0, 1 /)
call check( nf90_put_att(ncid, ncv, 'flag_values', nc3flag) )
buffer = 'cold_base '// &
         'temperate_base_with_cold_ice_above '// &
         'temperate_base_with_temperate_ice_above'
call check( nf90_put_att(ncid, ncv, 'flag_meanings', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- kc_cts

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'kc_cts', NF90_INT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'kc_cts', NF90_INT, nc2d, ncv) )
#endif

buffer = 'CTS_position_grid_index'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Grid index of the CTS position'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- zs

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'zs', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'zs', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'surface_altitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Topography of the free surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- zm

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'zm', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'zm', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'zm_interface_altitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Topography of the z=zm interface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- zb

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'zb', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'zb', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'ice_base_altitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Topography of the ice base'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- zl

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'zl', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'zl', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'bedrock_altitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Topography of the lithosphere surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- zl0

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'zl0', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'zl0', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'isostatically_relaxed_bedrock_altitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Topography of the isostatically relaxed lithosphere surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- wss

if (flag_wss) then

   call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
   call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
   call check( nf90_def_var(ncid, 'wss', NF90_FLOAT, nc2d, ncv, &
               deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
   call check( nf90_def_var(ncid, 'wss', NF90_FLOAT, nc2d, ncv) )
#endif

   buffer = 'm'
   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
   buffer = 'positive downward'
   call check( nf90_put_att(ncid, ncv, 'note', trim(buffer)) )
   buffer = 'isostatic_steady_state_bedrock_displacement'
   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
   buffer = 'Isostatic steady-state displacement of the lithosphere'
   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
   call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

end if

!    ---- H_cold

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'H_cold', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'H_cold', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_cold_layer_thickness'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Thickness of the cold ice layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- H_temp

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'H_temp', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'H_temp', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_temperate_layer_thickness'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Thickness of the temperate ice layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- H

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'H', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'H', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_thickness'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Ice thickness'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- H_R

call check( nf90_def_var(ncid, 'H_R', NF90_FLOAT, ncv) )

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'lithosphere_layer_thickness'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Thickness of the lithosphere layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )

!    ---- Q_bm

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'Q_bm', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'Q_bm', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_basal_melt_rate'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Basal melting rate'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- Q_tld

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'Q_tld', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'Q_tld', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_temperate_layer_water_drainage'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Water drainage from the temperate layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- calving

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'calving', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'calving', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_volume_flux_due_to_calving'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Calving flux'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- am_perp

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'am_perp', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'am_perp', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_volume_flux_across_zm_interface'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Volume flux across the z=zm interface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- qx

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'qx', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'qx', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm2 a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_vertical_integral_x_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Horizontal volume flux qx'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
buffer = 'Staggered grid variable, defined at (j,i+1/2)'
call check( nf90_put_att(ncid, ncv, 'comment', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- qy

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'qy', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'qy', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm2 a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_vertical_integral_y_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Horizontal volume flux qy'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
buffer = 'Staggered grid variable, defined at (j+1/2,i)'
call check( nf90_put_att(ncid, ncv, 'comment', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- vx_m_sia

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vx_m_sia', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'vx_m_sia', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_vertical_mean_x_velocity_sia'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Vertical mean of horizontal velocity vx (SIA)'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
buffer = 'Staggered grid variable, defined at (j,i+1/2)'
call check( nf90_put_att(ncid, ncv, 'comment', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- vy_m_sia

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vy_m_sia', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'vy_m_sia', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_vertical_mean_y_velocity_sia'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Vertical mean of horizontal velocity vy (SIA)'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
buffer = 'Staggered grid variable, defined at (j+1/2,i)'
call check( nf90_put_att(ncid, ncv, 'comment', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- vx_m_ssa

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vx_m_ssa', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'vx_m_ssa', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_vertical_mean_x_velocity_ssa'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Vertical mean of horizontal velocity vx (SSA)'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
buffer = 'Staggered grid variable, defined at (j,i+1/2)'
call check( nf90_put_att(ncid, ncv, 'comment', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- vy_m_ssa

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vy_m_ssa', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'vy_m_ssa', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_vertical_mean_y_velocity_ssa'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Vertical mean of horizontal velocity vy (SSA)'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
buffer = 'Staggered grid variable, defined at (j+1/2,i)'
call check( nf90_put_att(ncid, ncv, 'comment', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- dzs_dt (= dzs_dtau)

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'dzs_dt', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'dzs_dt', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'tendency_of_surface_altitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Rate of change of the topography of the free surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- dzm_dt (= dzm_dtau)

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'dzm_dt', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'dzm_dt', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'tendency_of_zm_interface_altitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Rate of change of the topography of the z=zm interface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- dzb_dt (= dzb_dtau)

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'dzb_dt', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'dzb_dt', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'tendency_of_ice_base_altitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Rate of change of the topography of the ice base'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- dzl_dt (= dzl_dtau)

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'dzl_dt', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'dzl_dt', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'tendency_of_bedrock_altitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Rate of change of the topography of the lithosphere surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- dH_c_dt (= dH_c_dtau)

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'dH_c_dt', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'dH_c_dt', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'tendency_of_land_ice_kc_layer_thickness'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Rate of change of the thickness of the upper (kc) ice layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- dH_t_dt (= dH_t_dtau)

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'dH_t_dt', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'dH_t_dt', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'tendency_of_land_ice_kt_layer_thickness'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Rate of change of the thickness of the lower (kt) ice layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- dH_dt (= dH_dtau)

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'dH_dt', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'dH_dt', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'tendency_of_land_ice_thickness'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Rate of change of the ice thickness'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- vx_b_g

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vx_b_g', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'vx_b_g', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_base_x_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Horizontal velocity vx at the ice base'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- vy_b_g

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vy_b_g', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'vy_b_g', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_base_y_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Horizontal velocity vy at the ice base'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- vz_b

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vz_b', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'vz_b', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_base_z_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Vertical velocity vz at the ice base'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- vh_b

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vh_b', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'vh_b', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_base_horizontal_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Horizontal velocity vh at the ice base'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- vx_s_g

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vx_s_g', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'vx_s_g', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_surface_x_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Horizontal velocity vx at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- vy_s_g

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vy_s_g', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'vy_s_g', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_surface_y_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Horizontal velocity vy at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- vz_s

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vz_s', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'vz_s', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_surface_z_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Vertical velocity vz at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- vh_s

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vh_s', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'vh_s', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_surface_horizontal_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Horizontal velocity vh at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- vx_m_g

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vx_m_g', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'vx_m_g', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_vertical_mean_x_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Vertical mean of horizontal velocity vx'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- vy_m_g

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vy_m_g', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'vy_m_g', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_vertical_mean_y_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Vertical mean of horizontal velocity vy'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- vh_m

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vh_m', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'vh_m', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_vertical_mean_horizontal_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Vertical mean of horizontal velocity vh'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- temp_b

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'temp_b', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'temp_b', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'degC'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'basal_temperature'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Temperature at the ice base'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- temph_b

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'temph_b', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'temph_b', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'degC'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'basal_temperature_rel_to_pmp'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Temperature at the ice base relative to the pressure melting point'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- tau_dr

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'tau_dr', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'tau_dr', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'Pa'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'magnitude_of_land_ice_driving_stress'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Driving stress'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- tau_b

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'tau_b', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'tau_b', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'Pa'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'magnitude_of_land_ice_basal_drag'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Basal drag'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- p_b_w

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'p_b_w', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'p_b_w', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'Pa'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'basal_water_pressure'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Basal water pressure'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- q_w

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'q_w', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'q_w', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm2 a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'basal_water_flux'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Volume flux of basal water'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- q_w_x

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'q_w_x', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'q_w_x', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm2 a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'basal_water_flux_x'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Volume flux of basal water in x-direction'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- q_w_y

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'q_w_y', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'q_w_y', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm2 a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'basal_water_flux_y'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Volume flux of basal water in y-direction'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- H_w

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'H_w', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'H_w', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'water_column_thickness'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Thickness of the water column under the ice base'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- q_gl_g

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'q_gl_g', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'q_gl_g', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'm2 a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_volume_flux_across_gl'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Horizontal volume flux across the grounding line'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- ratio_sl_sia_x

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'ratio_sl_sia_x', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'ratio_sl_sia_x', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_x_slip_ratio_sia'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Ratio of SIA basal to surface velocity (slip ratio) ' &
         // 'in x-direction, at (i+1/2,j)'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- ratio_sl_sia_y

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'ratio_sl_sia_y', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'ratio_sl_sia_y', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_y_slip_ratio_sia'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Ratio of SIA basal to surface velocity (slip ratio) ' &
         // 'in y-direction, at (i+1/2,j)'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- flag_shelfy_stream_x

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'flag_shelfy_stream_x', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'flag_shelfy_stream_x', NF90_BYTE, nc2d, ncv) )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_x_shelfy_stream_flag'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Shelfy stream flag in x-direction, at (i+1/2,j)'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- flag_shelfy_stream_y

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'flag_shelfy_stream_y', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'flag_shelfy_stream_y', NF90_BYTE, nc2d, ncv) )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_y_shelfy_stream_flag'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Shelfy stream flag in y-direction, at (i,j+1/2)'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- flag_shelfy_stream

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'flag_shelfy_stream', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'flag_shelfy_stream', NF90_BYTE, nc2d, ncv) )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_shelfy_stream_flag'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Shelfy stream flag'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- flag_grounding_line_1

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'flag_grounding_line_1', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'flag_grounding_line_1', NF90_BYTE, nc2d, ncv) )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_grounding_line_1_flag'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Grounding line flag 1'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- flag_grounding_line_2

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'flag_grounding_line_2', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'flag_grounding_line_2', NF90_BYTE, nc2d, ncv) )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_grounding_line_2_flag'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Grounding line flag 2'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- flag_calving_front_1

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'flag_calving_front_1', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'flag_calving_front_1', NF90_BYTE, nc2d, ncv) )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_calving_front_1_flag'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Calving front flag 1'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- flag_calving_front_2

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'flag_calving_front_2', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'flag_calving_front_2', NF90_BYTE, nc2d, ncv) )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_calving_front_2_flag'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Calving front flag 2'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- flag_grounded_front_a_1

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'flag_grounded_front_a_1', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'flag_grounded_front_a_1', NF90_BYTE, nc2d, ncv) )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_grounded_front_a_1_flag'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Grounded front flag a 1'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- flag_grounded_front_a_2

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'flag_grounded_front_a_2', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'flag_grounded_front_a_2', NF90_BYTE, nc2d, ncv) )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_grounded_front_a_2_flag'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Grounded front flag a 2'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- flag_grounded_front_b_1

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'flag_grounded_front_b_1', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'flag_grounded_front_b_1', NF90_BYTE, nc2d, ncv) )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_grounded_front_b_1_flag'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Grounded front flag b 1'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- flag_grounded_front_b_2

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'flag_grounded_front_b_2', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'flag_grounded_front_b_2', NF90_BYTE, nc2d, ncv) )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_grounded_front_b_2_flag'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Grounded front flag b 2'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- vis_ave_g

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vis_ave_g', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'vis_ave_g', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'Pa s'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_vertical_average_viscosity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Depth-averaged viscosity (SIA/SSA)'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- vis_int_g

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vis_int_g', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'vis_int_g', NF90_FLOAT, nc2d, ncv) )
#endif

buffer = 'Pa s m'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_vertical_integral_viscosity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Depth-integrated viscosity (SIA/SSA)'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- vx_c

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(3)), nc3d(3)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vx_c', NF90_FLOAT, nc3d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'vx_c', NF90_FLOAT, nc3d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_kc_layer_x_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Horizontal velocity vx in the upper (kc) ice layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
buffer = 'Staggered grid variable, defined at (kc,j,i+1/2)'
call check( nf90_put_att(ncid, ncv, 'comment', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- vy_c

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(3)), nc3d(3)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vy_c', NF90_FLOAT, nc3d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'vy_c', NF90_FLOAT, nc3d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_kc_layer_y_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Horizontal velocity vy in the upper (kc) ice layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
buffer = 'Staggered grid variable, defined at (kc,j+1/2,i)'
call check( nf90_put_att(ncid, ncv, 'comment', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- vz_c

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(3)), nc3d(3)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vz_c', NF90_FLOAT, nc3d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'vz_c', NF90_FLOAT, nc3d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_kc_layer_z_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Vertical velocity vz in the upper (kc) ice layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
buffer = 'Staggered grid variable, defined at (kc+1/2,j,i)'
call check( nf90_put_att(ncid, ncv, 'comment', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- vx_t

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(4)), nc3d(3)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vx_t', NF90_FLOAT, nc3d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'vx_t', NF90_FLOAT, nc3d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_kt_layer_x_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Horizontal velocity vx in the lower (kt) ice layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
buffer = 'Staggered grid variable, defined at (kt,j,i+1/2)'
call check( nf90_put_att(ncid, ncv, 'comment', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- vy_t

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(4)), nc3d(3)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vy_t', NF90_FLOAT, nc3d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'vy_t', NF90_FLOAT, nc3d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_kt_layer_y_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Horizontal velocity vy in the lower (kt) ice layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
buffer = 'Staggered grid variable, defined at (kt,j+1/2,i)'
call check( nf90_put_att(ncid, ncv, 'comment', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- vz_t

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(4)), nc3d(3)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vz_t', NF90_FLOAT, nc3d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'vz_t', NF90_FLOAT, nc3d, ncv) )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_kt_layer_z_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Vertical velocity vz in the lower (kt) ice layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
buffer = 'Staggered grid variable, defined at (kt+1/2,j,i)'
call check( nf90_put_att(ncid, ncv, 'comment', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- temp_c

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(3)), nc3d(3)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'temp_c', NF90_FLOAT, nc3d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'temp_c', NF90_FLOAT, nc3d, ncv) )
#endif

buffer = 'degC'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_kc_layer_temperature'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Temperature in the upper (kc) ice layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- omega_t

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(4)), nc3d(3)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'omega_t', NF90_FLOAT, nc3d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'omega_t', NF90_FLOAT, nc3d, ncv) )
#endif

buffer = '1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_kt_layer_water_content'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Water content in the lower (kt) ice layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- temp_r

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(5)), nc3d(3)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'temp_r', NF90_FLOAT, nc3d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'temp_r', NF90_FLOAT, nc3d, ncv) )
#endif

buffer = 'degC'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'lithosphere_layer_temperature'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Temperature in the lithosphere layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- enth_c

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(3)), nc3d(3)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'enth_c', NF90_FLOAT, nc3d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'enth_c', NF90_FLOAT, nc3d, ncv) )
#endif

buffer = 'J kg-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_kc_layer_enthalpy'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Enthalpy in the upper (kc) ice layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- enth_t

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(4)), nc3d(3)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'enth_t', NF90_FLOAT, nc3d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'enth_t', NF90_FLOAT, nc3d, ncv) )
#endif

buffer = 'J kg-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_kt_layer_enthalpy'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Enthalpy in the lower (kt) ice layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- omega_c

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(3)), nc3d(3)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'omega_c', NF90_FLOAT, nc3d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'omega_c', NF90_FLOAT, nc3d, ncv) )
#endif

buffer = '1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_kc_layer_water_content'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Water content in the upper (kc) ice layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- enh_c

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(3)), nc3d(3)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'enh_c', NF90_FLOAT, nc3d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'enh_c', NF90_FLOAT, nc3d, ncv) )
#endif

buffer = '1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_kc_layer_flow_enhancement_factor'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Flow enhancement factor in the upper (kc) ice layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- enh_t

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(4)), nc3d(3)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'enh_t', NF90_FLOAT, nc3d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'enh_t', NF90_FLOAT, nc3d, ncv) )
#endif

buffer = '1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_kt_layer_flow_enhancement_factor'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Flow enhancement factor in the lower (kt) ice layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- strain_heating_c

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(3)), nc3d(3)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'strain_heating_c', NF90_FLOAT, nc3d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'strain_heating_c', NF90_FLOAT, nc3d, ncv) )
#endif

buffer = 'W kg-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_kc_layer_strain_heating'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Strain heating in the upper (kc) ice layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- strain_heating_t

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(4)), nc3d(3)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'strain_heating_t', NF90_FLOAT, nc3d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'strain_heating_t', NF90_FLOAT, nc3d, ncv) )
#endif

buffer = 'W kg-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_kt_layer_strain_heating'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Strain heating in the lower (kt) ice layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- age_c

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(3)), nc3d(3)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'age_c', NF90_FLOAT, nc3d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'age_c', NF90_FLOAT, nc3d, ncv) )
#endif

buffer = 'a'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_kc_layer_age'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Age in the upper (kc) ice layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- age_t

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)) )
call check( nf90_inq_dimid(ncid, trim(coord_id(4)), nc3d(3)) )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'age_t', NF90_FLOAT, nc3d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle) )
#else
call check( nf90_def_var(ncid, 'age_t', NF90_FLOAT, nc3d, ncv) )
#endif

buffer = 'a'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)) )
buffer = 'land_ice_kt_layer_age'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)) )
buffer = 'Age in the lower (kt) ice layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)) )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping') )

!    ---- End of the definitions

call check( nf90_enddef(ncid) )

!-------- Writing of data on NetCDF file --------

write (6,'(/a)') ' Now writing data ...'

call check( nf90_inq_varid(ncid, 'mapping', ncv) )
call check( nf90_put_var(ncid, ncv, mapping_dbl) )

call check( nf90_inq_varid(ncid, 'year2sec', ncv) )
call check( nf90_put_var(ncid, ncv, year2sec_dbl) )

call check( nf90_inq_varid(ncid, 'time', ncv) )
call check( nf90_put_var(ncid, ncv, time_dbl) )

if (forcing_flag == 1) then

   call check( nf90_inq_varid(ncid, 'delta_ts', ncv) )
   call check( nf90_put_var(ncid, ncv, delta_ts_dbl) )

else if (forcing_flag == 2) then

   call check( nf90_inq_varid(ncid, 'glac_index', ncv) )
   call check( nf90_put_var(ncid, ncv, glac_index_dbl) )

end if

call check( nf90_inq_varid(ncid, 'z_sl_mean', ncv) )
call check( nf90_put_var(ncid, ncv, z_sl_mean_dbl) )

call check( nf90_inq_varid(ncid, 'V_tot', ncv) )
call check( nf90_put_var(ncid, ncv, V_tot_dbl) )

call check( nf90_inq_varid(ncid, 'V_af', ncv) )
call check( nf90_put_var(ncid, ncv, V_af_dbl) )

call check( nf90_inq_varid(ncid, 'A_grounded', ncv) )
call check( nf90_put_var(ncid, ncv, A_grounded_dbl) )

call check( nf90_inq_varid(ncid, 'A_floating', ncv) )
call check( nf90_put_var(ncid, ncv, A_floating_dbl) )

call check( nf90_inq_varid(ncid, 'x', ncv) )
call check( nf90_put_var(ncid, ncv, xi_dbl, &
                         start=nc1cor_i, count=nc1cnt_i) )

call check( nf90_inq_varid(ncid, 'y', ncv) )
call check( nf90_put_var(ncid, ncv, eta_dbl, &
                         start=nc1cor_j, count=nc1cnt_j) )

call check( nf90_inq_varid(ncid, 'sigma_level_c', ncv) )
call check( nf90_put_var(ncid, ncv, sigma_level_c_dbl, &
                         start=nc1cor_kc, count=nc1cnt_kc) )

call check( nf90_inq_varid(ncid, 'sigma_level_t', ncv) )
call check( nf90_put_var(ncid, ncv, sigma_level_t_dbl, &
                         start=nc1cor_kt, count=nc1cnt_kt) )

call check( nf90_inq_varid(ncid, 'sigma_level_r', ncv) )
call check( nf90_put_var(ncid, ncv, sigma_level_r_dbl, &
                         start=nc1cor_kr, count=nc1cnt_kr) )

call check( nf90_inq_varid(ncid, 'lon', ncv) )
call check( nf90_put_var(ncid, ncv, lon_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'lat', ncv) )
call check( nf90_put_var(ncid, ncv, lat_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'lambda', ncv) )
call check( nf90_put_var(ncid, ncv, lambda_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'phi', ncv) )
call check( nf90_put_var(ncid, ncv, phi_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'temp_s', ncv) )
call check( nf90_put_var(ncid, ncv, temp_s_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'prec', ncv) )
call check( nf90_put_var(ncid, ncv, prec_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'snowfall', ncv) )
call check( nf90_put_var(ncid, ncv, snowfall_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'rainfall', ncv) )
call check( nf90_put_var(ncid, ncv, rainfall_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'pdd', ncv) )
call check( nf90_put_var(ncid, ncv, pdd_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'as_perp', ncv) )
call check( nf90_put_var(ncid, ncv, as_perp_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'as_perp_apl', ncv) )
call check( nf90_put_var(ncid, ncv, as_perp_apl_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'smb_corr', ncv) )
call check( nf90_put_var(ncid, ncv, smb_corr_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'z_sl', ncv) )
call check( nf90_put_var(ncid, ncv, z_sl_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

#if (DISC>0)   /* Ice discharge parameterisation */

call check( nf90_inq_varid(ncid, 'dis_perp', ncv) )
call check( nf90_put_var(ncid, ncv, dis_perp_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'cst_dist', ncv) )
call check( nf90_put_var(ncid, ncv, cst_dist_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'cos_grad_tc', ncv) )
call check( nf90_put_var(ncid, ncv, cos_grad_tc_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'mask_mar', ncv) )
call check( nf90_put_var(ncid, ncv, mask_mar_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

#endif

call check( nf90_inq_varid(ncid, 'q_geo', ncv) )
call check( nf90_put_var(ncid, ncv, q_geo_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'mask', ncv) )
call check( nf90_put_var(ncid, ncv, mask_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'mask_old', ncv) )
call check( nf90_put_var(ncid, ncv, mask_old_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'n_cts', ncv) )
call check( nf90_put_var(ncid, ncv, n_cts_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'kc_cts', ncv) )
call check( nf90_put_var(ncid, ncv, kc_cts_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'zs', ncv) )
call check( nf90_put_var(ncid, ncv, zs_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'zm', ncv) )
call check( nf90_put_var(ncid, ncv, zm_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'zb', ncv) )
call check( nf90_put_var(ncid, ncv, zb_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'zl', ncv) )
call check( nf90_put_var(ncid, ncv, zl_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'zl0', ncv) )
call check( nf90_put_var(ncid, ncv, zl0_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

if (flag_wss) then
   call check( nf90_inq_varid(ncid, 'wss', ncv) )
   call check( nf90_put_var(ncid, ncv, wss_dbl, &
        start=nc2cor_ij, count=nc2cnt_ij) )
end if

call check( nf90_inq_varid(ncid, 'H_cold', ncv) )
call check( nf90_put_var(ncid, ncv, H_cold_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'H_temp', ncv) )
call check( nf90_put_var(ncid, ncv, H_temp_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'H', ncv) )
call check( nf90_put_var(ncid, ncv, H_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'H_R', ncv) )
call check( nf90_put_var(ncid, ncv, H_R_dbl) )

call check( nf90_inq_varid(ncid, 'Q_bm', ncv) )
call check( nf90_put_var(ncid, ncv, Q_bm_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'Q_tld', ncv) )
call check( nf90_put_var(ncid, ncv, Q_tld_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'calving', ncv) )
call check( nf90_put_var(ncid, ncv, calving_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'am_perp', ncv) )
call check( nf90_put_var(ncid, ncv, am_perp_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'qx', ncv) )
call check( nf90_put_var(ncid, ncv, qx_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'qy', ncv) )
call check( nf90_put_var(ncid, ncv, qy_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'vx_m_sia', ncv) )
call check( nf90_put_var(ncid, ncv, vx_m_sia_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'vy_m_sia', ncv) )
call check( nf90_put_var(ncid, ncv, vy_m_sia_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'vx_m_ssa', ncv) )
call check( nf90_put_var(ncid, ncv, vx_m_ssa_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'vy_m_ssa', ncv) )
call check( nf90_put_var(ncid, ncv, vy_m_ssa_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'dzs_dt', ncv) )
call check( nf90_put_var(ncid, ncv, dzs_dtau_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'dzm_dt', ncv) )
call check( nf90_put_var(ncid, ncv, dzm_dtau_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'dzb_dt', ncv) )
call check( nf90_put_var(ncid, ncv, dzb_dtau_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'dzl_dt', ncv) )
call check( nf90_put_var(ncid, ncv, dzl_dtau_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'dH_c_dt', ncv) )
call check( nf90_put_var(ncid, ncv, dH_c_dtau_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'dH_t_dt', ncv) )
call check( nf90_put_var(ncid, ncv, dH_t_dtau_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'dH_dt', ncv) )
call check( nf90_put_var(ncid, ncv, dH_dtau_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'vx_b_g', ncv) )
call check( nf90_put_var(ncid, ncv, vx_b_g_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'vy_b_g', ncv) )
call check( nf90_put_var(ncid, ncv, vy_b_g_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'vz_b', ncv) )
call check( nf90_put_var(ncid, ncv, vz_b_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'vh_b', ncv) )
call check( nf90_put_var(ncid, ncv, vh_b_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'vx_s_g', ncv) )
call check( nf90_put_var(ncid, ncv, vx_s_g_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'vy_s_g', ncv) )
call check( nf90_put_var(ncid, ncv, vy_s_g_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'vz_s', ncv) )
call check( nf90_put_var(ncid, ncv, vz_s_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'vh_s', ncv) )
call check( nf90_put_var(ncid, ncv, vh_s_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'vx_m_g', ncv) )
call check( nf90_put_var(ncid, ncv, vx_m_g_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'vy_m_g', ncv) )
call check( nf90_put_var(ncid, ncv, vy_m_g_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'vh_m', ncv) )
call check( nf90_put_var(ncid, ncv, vh_m_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'temp_b', ncv) )
call check( nf90_put_var(ncid, ncv, temp_b_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'temph_b', ncv) )
call check( nf90_put_var(ncid, ncv, temph_b_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'tau_dr', ncv) )
call check( nf90_put_var(ncid, ncv, tau_dr_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'tau_b', ncv) )
call check( nf90_put_var(ncid, ncv, tau_b_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'p_b_w', ncv) )
call check( nf90_put_var(ncid, ncv, p_b_w_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'q_w', ncv) )
call check( nf90_put_var(ncid, ncv, q_w_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'q_w_x', ncv) )
call check( nf90_put_var(ncid, ncv, q_w_x_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'q_w_y', ncv) )
call check( nf90_put_var(ncid, ncv, q_w_y_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'H_w', ncv) )
call check( nf90_put_var(ncid, ncv, H_w_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'q_gl_g', ncv) )
call check( nf90_put_var(ncid, ncv, q_gl_g_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'ratio_sl_sia_x', ncv) )
call check( nf90_put_var(ncid, ncv, ratio_sl_sia_x_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'ratio_sl_sia_y', ncv) )
call check( nf90_put_var(ncid, ncv, ratio_sl_sia_y_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'flag_shelfy_stream_x', ncv) )
call check( nf90_put_var(ncid, ncv, flag_shelfy_stream_x_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'flag_shelfy_stream_y', ncv) )
call check( nf90_put_var(ncid, ncv, flag_shelfy_stream_y_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'flag_shelfy_stream', ncv) )
call check( nf90_put_var(ncid, ncv, flag_shelfy_stream_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'flag_grounding_line_1', ncv) )
call check( nf90_put_var(ncid, ncv, flag_grounding_line_1_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'flag_grounding_line_2', ncv) )
call check( nf90_put_var(ncid, ncv, flag_grounding_line_2_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'flag_calving_front_1', ncv) )
call check( nf90_put_var(ncid, ncv, flag_calving_front_1_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'flag_calving_front_2', ncv) )
call check( nf90_put_var(ncid, ncv, flag_calving_front_2_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'flag_grounded_front_a_1', ncv) )
call check( nf90_put_var(ncid, ncv, flag_grounded_front_a_1_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'flag_grounded_front_a_2', ncv) )
call check( nf90_put_var(ncid, ncv, flag_grounded_front_a_2_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'flag_grounded_front_b_1', ncv) )
call check( nf90_put_var(ncid, ncv, flag_grounded_front_b_1_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'flag_grounded_front_b_2', ncv) )
call check( nf90_put_var(ncid, ncv, flag_grounded_front_b_2_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'vis_ave_g', ncv) )
call check( nf90_put_var(ncid, ncv, vis_ave_g_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'vis_int_g', ncv) )
call check( nf90_put_var(ncid, ncv, vis_int_g_dbl, &
                         start=nc2cor_ij, count=nc2cnt_ij) )

call check( nf90_inq_varid(ncid, 'vx_c', ncv) )
call check( nf90_put_var(ncid, ncv, vx_c_dbl, &
                         start=nc3cor_ijkc, count=nc3cnt_ijkc) )

call check( nf90_inq_varid(ncid, 'vy_c', ncv) )
call check( nf90_put_var(ncid, ncv, vy_c_dbl, &
                         start=nc3cor_ijkc, count=nc3cnt_ijkc) )

call check( nf90_inq_varid(ncid, 'vz_c', ncv) )
call check( nf90_put_var(ncid, ncv, vz_c_dbl, &
                         start=nc3cor_ijkc, count=nc3cnt_ijkc) )

call check( nf90_inq_varid(ncid, 'vx_t', ncv) )
call check( nf90_put_var(ncid, ncv, vx_t_dbl, &
                         start=nc3cor_ijkt, count=nc3cnt_ijkt) )

call check( nf90_inq_varid(ncid, 'vy_t', ncv) )
call check( nf90_put_var(ncid, ncv, vy_t_dbl, &
                         start=nc3cor_ijkt, count=nc3cnt_ijkt) )

call check( nf90_inq_varid(ncid, 'vz_t', ncv) )
call check( nf90_put_var(ncid, ncv, vz_t_dbl, &
                         start=nc3cor_ijkt, count=nc3cnt_ijkt) )

call check( nf90_inq_varid(ncid, 'temp_c', ncv) )
call check( nf90_put_var(ncid, ncv, temp_c_dbl, &
                         start=nc3cor_ijkc, count=nc3cnt_ijkc) )

call check( nf90_inq_varid(ncid, 'omega_t', ncv) )
call check( nf90_put_var(ncid, ncv, omega_t_dbl, &
                         start=nc3cor_ijkt, count=nc3cnt_ijkt) )

call check( nf90_inq_varid(ncid, 'temp_r', ncv) )
call check( nf90_put_var(ncid, ncv, temp_r_dbl, &
                         start=nc3cor_ijkr, count=nc3cnt_ijkr) )

call check( nf90_inq_varid(ncid, 'enth_c', ncv) )
call check( nf90_put_var(ncid, ncv, enth_c_dbl, &
                         start=nc3cor_ijkc, count=nc3cnt_ijkc) )

call check( nf90_inq_varid(ncid, 'enth_t', ncv) )
call check( nf90_put_var(ncid, ncv, enth_t_dbl, &
                         start=nc3cor_ijkt, count=nc3cnt_ijkt) )

call check( nf90_inq_varid(ncid, 'omega_c', ncv) )
call check( nf90_put_var(ncid, ncv, omega_c_dbl, &
                         start=nc3cor_ijkc, count=nc3cnt_ijkc) )

call check( nf90_inq_varid(ncid, 'enh_c', ncv) )
call check( nf90_put_var(ncid, ncv, enh_c_dbl, &
                         start=nc3cor_ijkc, count=nc3cnt_ijkc) )

call check( nf90_inq_varid(ncid, 'enh_t', ncv) )
call check( nf90_put_var(ncid, ncv, enh_t_dbl, &
                         start=nc3cor_ijkt, count=nc3cnt_ijkt) )

call check( nf90_inq_varid(ncid, 'strain_heating_c', ncv) )
call check( nf90_put_var(ncid, ncv, strain_heating_c_dbl, &
                         start=nc3cor_ijkc, count=nc3cnt_ijkc) )

call check( nf90_inq_varid(ncid, 'strain_heating_t', ncv) )
call check( nf90_put_var(ncid, ncv, strain_heating_t_dbl, &
                         start=nc3cor_ijkt, count=nc3cnt_ijkt) )

call check( nf90_inq_varid(ncid, 'age_c', ncv) )
call check( nf90_put_var(ncid, ncv, age_c_dbl, &
                         start=nc3cor_ijkc, count=nc3cnt_ijkc) )

call check( nf90_inq_varid(ncid, 'age_t', ncv) )
call check( nf90_put_var(ncid, ncv, age_t_dbl, &
                         start=nc3cor_ijkt, count=nc3cnt_ijkt) )

!-------- Closing of NetCDF file --------

write (6,'(/a)') ' Now closing new NetCDF file ...'

call check( nf90_sync(ncid) )

call check( nf90_close(ncid) )

write (6,'(/a/)') ' Done!'

end subroutine write_nc_double

!-------------------------------------------------------------------------------
!> Set the value of the institution string ch_institution.
!-------------------------------------------------------------------------------
  subroutine set_ch_institution(ch_institution)

  use resolution_doubler_types

  implicit none

  character(len=*), intent(out) :: ch_institution

  integer(i4b)       :: istat
  character(len=256) :: ch_institution_default, ch_value

  ch_institution_default = 'Institute of Low Temperature Science, '// &
                           'Hokkaido University'

  call get_environment_variable(name='SICO_INSTITUTION', value=ch_value, &
                                status=istat, trim_name=.true.)

  if (istat /= 0) then

     ch_institution = ch_institution_default

  else

     if (     (trim(ch_value)=='default') &
          .or.(trim(ch_value)=='Default') &
          .or.(trim(ch_value)=='DEFAULT') ) then

        ch_institution = ch_institution_default

     else

        ch_institution = trim(ch_value)

     end if

  end if

  end subroutine set_ch_institution

!-------------------------------------------------------------------------------
!> Set the value of the website string ch_website.
!-------------------------------------------------------------------------------
  subroutine set_ch_website(ch_website)

  implicit none

  character(len=*), intent(out) :: ch_website

  integer(i4b)       :: istat
  character(len=256) :: ch_website_default, ch_value

  ch_website_default = ' '

  call get_environment_variable(name='SICO_WEBSITE', value=ch_value, &
                                status=istat, trim_name=.true.)

  if (istat /= 0) then
    ch_website = ch_website_default
  else
     ch_website = trim(ch_value)
  end if

  end subroutine set_ch_website

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

  stop

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
!> NetCDF error capturing.
!-------------------------------------------------------------------------------
  subroutine check(status)

  use resolution_doubler_types
  use resolution_doubler_vars
  use netcdf

  implicit none

  integer(i4b), intent(in) :: status

  if (status /= nf90_noerr) then
     write(6,'(1x,a)') trim(nf90_strerror(status))
     stop ' >>> check: Stopped due to NetCDF error!'
  end if

  end subroutine check  

!-------- End of program --------

end program resolution_doubler

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                        End of resolution_doubler.F90
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
