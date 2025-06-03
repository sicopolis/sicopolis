!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  r e a d _ m
!
!! Reading of several input data.
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
!> Reading of several input data.
!-------------------------------------------------------------------------------
module read_m

  use sico_types_m
  use error_m

  implicit none

  private
  public :: read_tms_nc, read_target_topo_nc, &
            read_scalar_input, read_2d_input, read_phys_para, read_kei

#if (defined(ALLOW_NORMAL) || defined(ALLOW_GRDCHK) || defined(ALLOW_TAPENADE))
  public :: read_age_data, read_BedMachine_data
#endif

contains

!-------------------------------------------------------------------------------
!> Reading of time-slice files in NetCDF format.
!-------------------------------------------------------------------------------
  subroutine read_tms_nc(filename, &
                         opt_mask, opt_n_cts, opt_kc_cts, &
                         opt_H_cold, opt_H_temp, opt_H, &
                         opt_temp_r, opt_omega_t, opt_age_t, &
                         opt_temp_c, opt_age_c, opt_omega_c, &
                         opt_flag_temp_age_only)

  use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

  use netcdf
  use nc_check_m

  implicit none

  character(len=256), intent(in) :: filename

  logical, optional,  intent(in)    :: opt_flag_temp_age_only
  integer(i4b), optional, dimension(0:JMAX,0:IMAX), &
                      intent(inout) :: opt_mask, opt_n_cts
  integer(i4b), optional, dimension(0:JMAX,0:IMAX), &
                      intent(inout) :: opt_kc_cts
  real(dp), optional, dimension(0:JMAX,0:IMAX), &
                      intent(inout) :: opt_H_cold, opt_H_temp, opt_H
  real(dp), optional, dimension(0:KRMAX,0:JMAX,0:IMAX), &
                      intent(inout) :: opt_temp_r
  real(dp), optional, dimension(0:KTMAX,0:JMAX,0:IMAX), &
                      intent(inout) :: opt_omega_t, opt_age_t
  real(dp), optional, dimension(0:KCMAX,0:JMAX,0:IMAX), &
                      intent(inout) :: opt_temp_c, opt_age_c, opt_omega_c

  integer(i4b)       :: i, j, kc, kt, kr, n
  integer(i4b)       :: ios, istat, istat2
  character(len=256) :: anfdat_path
  character(len=256) :: filename_with_path
  logical            :: flag_temp_age_only
  logical            :: flag_z_sl_xy_array

  integer(i4b) :: ncid, ncv
  !     ncid:      ID of the output file
  !     ncv:       Variable ID

  integer(i4b), dimension(0:IMAX,0:JMAX) :: mask_conv, mask_old_conv, &
                                            n_cts_conv, &
                                            flag_shelfy_stream_x_conv, &
                                            flag_shelfy_stream_y_conv, &
                                            flag_shelfy_stream_conv, &
                                            flag_grounding_line_1_conv, &
                                            flag_grounding_line_2_conv, &
                                            flag_calving_front_1_conv, &
                                            flag_calving_front_2_conv, &
                                            flag_grounded_front_a_1_conv, &
                                            flag_grounded_front_a_2_conv, &
                                            flag_grounded_front_b_1_conv, &
                                            flag_grounded_front_b_2_conv
  integer(i4b), dimension(0:IMAX,0:JMAX) :: kc_cts_conv

  real(dp) :: year2sec_conv, time_conv, &
              delta_ts_conv, glac_index_conv, z_sl_mean_conv, &
              xi_conv(0:IMAX), eta_conv(0:JMAX), &
              sigma_level_c_conv(0:KCMAX), sigma_level_t_conv(0:KTMAX), &
              sigma_level_r_conv(0:KRMAX)

  real(sp) :: H_R_conv

  real(sp), dimension(0:IMAX,0:JMAX) :: lambda_conv, phi_conv, &
              lon_conv, lat_conv, &
              cell_area_conv, &
              temp_maat_conv, temp_s_conv, accum_conv, &
              snowfall_conv, rainfall_conv, pdd_conv, & 
              as_perp_conv, as_perp_apl_conv, smb_corr_conv, &
              z_sl_conv, &
              q_geo_conv, &
              zs_conv, zm_conv, zb_conv, zl_conv, zl0_conv, wss_conv, &
              H_cold_conv, H_temp_conv, H_conv, &
              Q_bm_conv, Q_tld_conv, &
              am_perp_conv, &
              qx_conv, qy_conv, &
              vx_m_sia_conv, vy_m_sia_conv, vx_m_ssa_conv, vy_m_ssa_conv, &
              dzs_dtau_conv, dzm_dtau_conv, dzb_dtau_conv, dzl_dtau_conv, &
              dH_c_dtau_conv, dH_t_dtau_conv, dH_dtau_conv, &
              vx_b_g_conv, vy_b_g_conv, vz_b_conv, vh_b_conv, &
              vx_s_g_conv, vy_s_g_conv, vz_s_conv, vh_s_conv, &
              vx_m_g_conv, vy_m_g_conv,            vh_m_conv, &
              temp_b_conv, temph_b_conv, &
              tau_dr_conv, tau_b_conv, &
              p_b_w_conv, q_w_conv, q_w_x_conv, q_w_y_conv, H_w_conv, &
              q_gl_g_conv, &
              ratio_sl_sia_x_conv, ratio_sl_sia_y_conv, ratio_sl_sia_conv, &
              vis_ave_g_conv, vis_int_g_conv
  real(sp), dimension(0:IMAX,0:JMAX,0:KCMAX) :: vx_c_conv, vy_c_conv, vz_c_conv, &
                                                temp_c_conv, age_c_conv, &
                                                enth_c_conv, omega_c_conv, &
                                                enh_c_conv, strain_heating_c_conv
  real(sp), dimension(0:IMAX,0:JMAX,0:KTMAX) :: vx_t_conv, vy_t_conv, vz_t_conv, &
                                                omega_t_conv, age_t_conv, &
                                                enth_t_conv, &
                                                enh_t_conv, strain_heating_t_conv
  real(sp), dimension(0:IMAX,0:JMAX,0:KRMAX) :: temp_r_conv

#if (DISC>0)   /* Ice discharge parameterization */
  integer(i4b), dimension(0:IMAX,0:JMAX) :: mask_mar_conv
  real(sp),     dimension(0:IMAX,0:JMAX) :: dis_perp_conv, &
                                            cst_dist_conv, cos_grad_tc_conv
#endif

  real(dp) :: visc_min, visc_max, visc_init
  logical :: flag_ratio_sl_sia, flag_vis_ave_g

  real(dp), parameter :: one_year = 1.0_dp

!-------- Read data from time-slice file of previous simulation --------

  if (present(opt_flag_temp_age_only)) then
     flag_temp_age_only = opt_flag_temp_age_only
  else
     flag_temp_age_only = .false.
  end if

#if (defined(ANF_DAT_PATH))
  anfdat_path = trim(ANF_DAT_PATH)
#else
  anfdat_path = 'dummy'
  errormsg = ' >>> read_tms_nc: ANF_DAT_PATH must be defined!'
  call error(errormsg)
#endif

  filename_with_path = trim(anfdat_path)//'/'//trim(filename)

  ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

  if (ios /= nf90_noerr) then
     errormsg = ' >>> read_tms_nc: Error when opening the' &
              //                   end_of_line &
              //'                  time-slice NetCDF file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'year2sec', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, year2sec_conv) )
  else
     year2sec_conv = 0.0_dp
  end if

  istat = nf90_inq_varid(ncid, 'time', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, time_conv) )
  else
     time_conv = 0.0_dp
  end if

  istat = nf90_inq_varid(ncid, 'delta_ts', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, delta_ts_conv) )
  else
     delta_ts_conv = 0.0_dp
  end if

  istat = nf90_inq_varid(ncid, 'glac_index', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, glac_index_conv) )
  else
     glac_index_conv = 0.0_dp
  end if

  istat = nf90_inq_varid(ncid, 'z_sl_mean', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, z_sl_mean_conv) )
     flag_z_sl_xy_array = .true.
  else
     istat2 = nf90_inq_varid(ncid, 'z_sl', ncv)
     if (istat2 == nf90_noerr) then
        call check( nf90_get_var(ncid, ncv, z_sl_mean_conv) )
        flag_z_sl_xy_array = .false.
     else
        errormsg = ' >>> read_tms_nc: Variable ''z_sl(_mean)'' ' &
                 //                   end_of_line &
                 //'                  not available in read nc file!'
        call error(errormsg)
     end if
  end if

  istat = nf90_inq_varid(ncid, 'x', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, xi_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''x'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'y', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, eta_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''y'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'sigma_level_c', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, sigma_level_c_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''sigma_level_c'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'sigma_level_t', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, sigma_level_t_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''sigma_level_t'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'sigma_level_r', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, sigma_level_r_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''sigma_level_r'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'lon', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, lon_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''lon'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'lat', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, lat_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''lat'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'lambda', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, lambda_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''lambda'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'phi', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, phi_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''phi'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'cell_area', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, cell_area_conv) )
  else
     warningmsg = ' >>> read_tms_nc: Variable ''cell_area'' ' &
                //                   end_of_line &
                //'                  not available in read nc file.'
     call warning(warningmsg)
     cell_area_conv = 0.0_sp
  end if

  istat = nf90_inq_varid(ncid, 'temp_maat', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, temp_maat_conv) )
  else
     warningmsg = ' >>> read_tms_nc: Variable ''temp_maat'' ' &
                //                   end_of_line &
                //'                  not available in read nc file.'
     call warning(warningmsg)
     temp_maat_conv = 0.0_sp
  end if

  istat = nf90_inq_varid(ncid, 'temp_s', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, temp_s_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''temp_s'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'prec', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, accum_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''prec'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'snowfall', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, snowfall_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''snowfall'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'rainfall', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, rainfall_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''rainfall'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'pdd', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, pdd_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''pdd'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'as_perp', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, as_perp_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''as_perp'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'as_perp_apl', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, as_perp_apl_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''as_perp_apl'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'smb_corr', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, smb_corr_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''smb_corr'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  if (flag_z_sl_xy_array) then
     istat = nf90_inq_varid(ncid, 'z_sl', ncv)
     if (istat == nf90_noerr) then
        call check( nf90_get_var(ncid, ncv, z_sl_conv) )
     else
        errormsg = ' >>> read_tms_nc: Variable ''z_sl'' ' &
                 //                   end_of_line &
                 //'                  not available in read nc file!'
        call error(errormsg)
     end if
  else
     z_sl_conv = z_sl_mean_conv
  end if

#if (DISC>0)   /* Ice discharge parameterization */

  istat = nf90_inq_varid(ncid, 'dis_perp', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, dis_perp_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''dis_perp'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'cst_dist', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, cst_dist_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''cst_dist'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'cos_grad_tc', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, cos_grad_tc_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''cos_grad_tc'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'mask_mar', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, mask_mar_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''mask_mar'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

#endif

  istat = nf90_inq_varid(ncid, 'q_geo', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, q_geo_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''q_geo'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'mask', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, mask_conv) )
  else
     istat2 = nf90_inq_varid(ncid, 'maske', ncv)
     if (istat2 == nf90_noerr) then
        call check( nf90_get_var(ncid, ncv, mask_conv) )
     else
        errormsg = ' >>> read_tms_nc: Variable ''mask'' ' &
                 //                   end_of_line &
                 //'                  not available in read nc file!'
        call error(errormsg)
     end if
  end if

  istat = nf90_inq_varid(ncid, 'mask_old', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, mask_old_conv) )
  else
     istat2 = nf90_inq_varid(ncid, 'maske_old', ncv)
     if (istat2 == nf90_noerr) then
        call check( nf90_get_var(ncid, ncv, mask_old_conv) )
     else
        errormsg = ' >>> read_tms_nc: Variable ''mask_old'' ' &
                 //                   end_of_line &
                 //'                  not available in read nc file!'
        call error(errormsg)
     end if
  end if

  istat = nf90_inq_varid(ncid, 'n_cts', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, n_cts_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''n_cts'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'kc_cts', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, kc_cts_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''kc_cts'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'zs', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, zs_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''zs'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'zm', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, zm_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''zm'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'zb', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, zb_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''zb'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'zl', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, zl_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''zl'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'zl0', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, zl0_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''zl0'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'wss', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, wss_conv) )
  else
#if ((REBOUND==2) && (ANF_DAT==3) && !defined(LEGACY_RESTART))
     errormsg = ' >>> read_tms_nc: Variable ''wss'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
#else
     warningmsg = ' >>> read_tms_nc: Variable ''wss'' ' &
                //                   end_of_line &
                //'                  not available in read nc file.'
     call warning(warningmsg)
     wss_conv = 0.0_sp
#endif
  end if

  istat = nf90_inq_varid(ncid, 'H_cold', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, H_cold_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''H_cold'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'H_temp', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, H_temp_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''H_temp'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'H', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, H_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''H'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'H_R', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, H_R_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''H_R'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'Q_bm', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, Q_bm_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''Q_bm'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'Q_tld', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, Q_tld_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''Q_tld'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'am_perp', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, am_perp_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''am_perp'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'qx', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, qx_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''qx'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'qy', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, qy_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''qy'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'vx_m_sia', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, vx_m_sia_conv) )
  else
     warningmsg = ' >>> read_tms_nc: Variable ''vx_m_sia'' ' &
                //                   end_of_line &
                //'                  not available in read nc file.'
     call warning(warningmsg)
     vx_m_sia_conv = 0.0_sp
  end if

  istat = nf90_inq_varid(ncid, 'vy_m_sia', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, vy_m_sia_conv) )
  else
     warningmsg = ' >>> read_tms_nc: Variable ''vy_m_sia'' ' &
                //                   end_of_line &
                //'                  not available in read nc file.'
     call warning(warningmsg)
     vy_m_sia_conv = 0.0_sp
  end if

  istat = nf90_inq_varid(ncid, 'vx_m_ssa', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, vx_m_ssa_conv) )
  else
     warningmsg = ' >>> read_tms_nc: Variable ''vx_m_ssa'' ' &
                //                   end_of_line &
                //'                  not available in read nc file.'
     call warning(warningmsg)
     vx_m_ssa_conv = 0.0_sp
  end if

  istat = nf90_inq_varid(ncid, 'vy_m_ssa', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, vy_m_ssa_conv) )
  else
     warningmsg = ' >>> read_tms_nc: Variable ''vy_m_ssa'' ' &
                //                   end_of_line &
                //'                  not available in read nc file.'
     call warning(warningmsg)
     vy_m_ssa_conv = 0.0_sp
  end if

  istat = nf90_inq_varid(ncid, 'dzs_dt', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, dzs_dtau_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''dzs_dt'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'dzm_dt', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, dzm_dtau_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''dzm_dt'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'dzb_dt', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, dzb_dtau_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''dzb_dt'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'dzl_dt', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, dzl_dtau_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''dzl_dt'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'dH_c_dt', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, dH_c_dtau_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''dH_c_dt'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'dH_t_dt', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, dH_t_dtau_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''dH_t_dt'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'dH_dt', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, dH_dtau_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''dH_dt'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'vx_b_g', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, vx_b_g_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''vx_b_g'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'vy_b_g', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, vy_b_g_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''vy_b_g'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'vz_b', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, vz_b_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''vz_b'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'vh_b', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, vh_b_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''vh_b'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'vx_s_g', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, vx_s_g_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''vx_s_g'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'vy_s_g', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, vy_s_g_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''vy_s_g'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'vz_s', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, vz_s_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''vz_s'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'vh_s', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, vh_s_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''vh_s'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'vx_m_g', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, vx_m_g_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''vx_m_g'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'vy_m_g', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, vy_m_g_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''vy_m_g'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'vh_m', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, vh_m_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''vh_m'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'temp_b', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, temp_b_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''temp_b'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'temph_b', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, temph_b_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''temph_b'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'tau_dr', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, tau_dr_conv) )
  else
     warningmsg = ' >>> read_tms_nc: Variable ''tau_dr'' ' &
                //                   end_of_line &
                //'                  not available in read nc file.'
     call warning(warningmsg)
     tau_dr_conv = 0.0_sp
  end if

  istat = nf90_inq_varid(ncid, 'tau_b', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, tau_b_conv) )
  else
     warningmsg = ' >>> read_tms_nc: Variable ''tau_b'' ' &
                //                   end_of_line &
                //'                  not available in read nc file.'
     call warning(warningmsg)
     tau_b_conv = 0.0_sp
  end if

  istat = nf90_inq_varid(ncid, 'p_b_w', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, p_b_w_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''p_b_w'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'q_w', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, q_w_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''q_w'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'q_w_x', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, q_w_x_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''q_w_x'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'q_w_y', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, q_w_y_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''q_w_y'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'H_w', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, H_w_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''H_w'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'q_gl_g', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, q_gl_g_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''q_gl_g'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'ratio_sl_sia_x', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, ratio_sl_sia_x_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''ratio_sl_sia_x'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'ratio_sl_sia_y', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, ratio_sl_sia_y_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''ratio_sl_sia_y'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'ratio_sl_sia', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, ratio_sl_sia_conv) )
     flag_ratio_sl_sia = .true.
  else
     warningmsg = ' >>> read_tms_nc: Variable ''ratio_sl_sia'' ' &
                //                   end_of_line &
                //'                  not available in read nc file.'
     call warning(warningmsg)
     ratio_sl_sia_conv = 0.0_sp
     flag_ratio_sl_sia = .false.
  end if

  istat = nf90_inq_varid(ncid, 'flag_shelfy_stream_x', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, flag_shelfy_stream_x_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''flag_shelfy_stream_x'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'flag_shelfy_stream_y', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, flag_shelfy_stream_y_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''flag_shelfy_stream_y'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'flag_shelfy_stream', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, flag_shelfy_stream_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''flag_shelfy_stream'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'flag_grounding_line_1', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, flag_grounding_line_1_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''flag_grounding_line_1'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'flag_grounding_line_2', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, flag_grounding_line_2_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''flag_grounding_line_2'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'flag_calving_front_1', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, flag_calving_front_1_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''flag_calving_front_1'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'flag_calving_front_2', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, flag_calving_front_2_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''flag_calving_front_2'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'flag_grounded_front_a_1', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, flag_grounded_front_a_1_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''flag_calving_front_a_1'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'flag_grounded_front_a_2', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, flag_grounded_front_a_2_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''flag_calving_front_a_2'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'flag_grounded_front_b_1', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, flag_grounded_front_b_1_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''flag_calving_front_b_1'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'flag_grounded_front_b_2', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, flag_grounded_front_b_2_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''flag_calving_front_b_2'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'vis_ave_g', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, vis_ave_g_conv) )
     flag_vis_ave_g = .true.
  else
     warningmsg = ' >>> read_tms_nc: Variable ''vis_ave_g'' ' &
                //                   end_of_line &
                //'                  not available in read nc file.'
     call warning(warningmsg)
     vis_ave_g_conv = 0.0_sp
     flag_vis_ave_g = .false.
  end if

  istat = nf90_inq_varid(ncid, 'vis_int_g', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, vis_int_g_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''vis_int_g'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'vx_c', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, vx_c_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''vx_c'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'vy_c', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, vy_c_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''vy_c'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'vz_c', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, vz_c_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''vz_c'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'vx_t', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, vx_t_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''vx_t'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'vy_t', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, vy_t_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''vy_t'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'vz_t', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, vz_t_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''vz_t'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'temp_c', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, temp_c_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''temp_c'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'omega_t', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, omega_t_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''omega_t'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'temp_r', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, temp_r_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''temp_r'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'enth_c', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, enth_c_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''enth_c'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'enth_t', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, enth_t_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''enth_t'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'omega_c', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, omega_c_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''omega_c'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'enh_c', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, enh_c_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''enh_c'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'enh_t', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, enh_t_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''enh_t'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'strain_heating_c', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, strain_heating_c_conv) )
  else
     warningmsg = ' >>> read_tms_nc: Variable ''strain_heating_c'' ' &
                //                   end_of_line &
                //'                  not available in read nc file.'
     call warning(warningmsg)
     strain_heating_c_conv = 0.0_sp
  end if

  istat = nf90_inq_varid(ncid, 'strain_heating_t', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, strain_heating_t_conv) )
  else
     warningmsg = ' >>> read_tms_nc: Variable ''strain_heating_t'' ' &
                //                   end_of_line &
                //'                  not available in read nc file.'
     call warning(warningmsg)
     strain_heating_t_conv = 0.0_sp
  end if

  istat = nf90_inq_varid(ncid, 'age_c', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, age_c_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''age_c'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'age_t', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, age_t_conv) )
  else
     errormsg = ' >>> read_tms_nc: Variable ''age_t'' ' &
              //                   end_of_line &
              //'                  not available in read nc file!'
     call error(errormsg)
  end if

  call check( nf90_close(ncid) )

!-------- Convert data to real*8 and years to seconds --------

  if (.not.flag_temp_age_only) then

     ! year2sec = year2sec_conv   ! set in sico_init
     ! time     = time_conv       ! set in sico_init

     delta_ts   = delta_ts_conv
     glac_index = glac_index_conv
     z_sl_mean  = z_sl_mean_conv

     do i=0, IMAX
        xi(i) = xi_conv(i)
     end do

     do j=0, JMAX
        eta(j) = eta_conv(j)
     end do

     H_R = real(H_R_conv,dp)

     do i=0, IMAX
     do j=0, JMAX

        mask(j,i)     = mask_conv(i,j)
        mask_old(j,i) = mask_old_conv(i,j)
        n_cts(j,i)    = n_cts_conv(i,j)
        kc_cts(j,i)   = kc_cts_conv(i,j)

        temp_maat(j,i)   = real(temp_maat_conv(i,j),dp)
        temp_s(j,i)      = real(temp_s_conv(i,j),dp)
        accum(j,i)       = real(accum_conv(i,j),dp)
        snowfall(j,i)    = real(snowfall_conv(i,j),dp)
        rainfall(j,i)    = real(rainfall_conv(i,j),dp)
        ET(j,i)          = real(pdd_conv(i,j),dp)/one_year
        as_perp(j,i)     = real(as_perp_conv(i,j),dp)
        as_perp_apl(j,i) = real(as_perp_apl_conv(i,j),dp)
        smb_corr(j,i)    = real(smb_corr_conv(i,j),dp)

        z_sl(j,i) = real(z_sl_conv(i,j),dp)

#if (DISC>0)   /* Ice discharge parameterization */
        dis_perp(j,i)    = real(dis_perp_conv(i,j),dp)
        cst_dist(j,i)    = real(cst_dist_conv(i,j),dp)
        cos_grad_tc(j,i) = real(cos_grad_tc_conv(i,j),dp)
        mask_mar(j,i)    = mask_mar_conv(i,j)
#endif

        q_geo(j,i) = real(q_geo_conv(i,j),dp)

        zs(j,i)   = real(zs_conv(i,j),dp)
        zm(j,i)   = real(zm_conv(i,j),dp)
        zb(j,i)   = real(zb_conv(i,j),dp)
        zl(j,i)   = real(zl_conv(i,j),dp)
        zl0(j,i)  = real(zl0_conv(i,j),dp)
        wss(j,i)  = real(wss_conv(i,j),dp)
        H(j,i)    = real(H_conv(i,j),dp)
#if (CALCMOD==1)
        H_c(j,i)  = real(H_cold_conv(i,j),dp)
        H_t(j,i)  = real(H_temp_conv(i,j),dp)
#elif (CALCMOD==0 || CALCMOD==2 || CALCMOD==3 || CALCMOD==-1)
        H_c(j,i)  = H(j,i)
        H_t(j,i)  = 0.0_dp
#endif
        Q_bm(j,i)    = real(Q_bm_conv(i,j),dp)*sec2year
        Q_tld(j,i)   = real(Q_tld_conv(i,j),dp)*sec2year
        am_perp(j,i) = real(am_perp_conv(i,j),dp)*sec2year
        qx(j,i)      = real(qx_conv(i,j),dp)*sec2year
        qy(j,i)      = real(qy_conv(i,j),dp)*sec2year
        vx_m_sia(j,i) = real(vx_m_sia_conv(i,j),dp)*sec2year
        vy_m_sia(j,i) = real(vy_m_sia_conv(i,j),dp)*sec2year
        vx_m_ssa(j,i) = real(vx_m_ssa_conv(i,j),dp)*sec2year
        vy_m_ssa(j,i) = real(vy_m_ssa_conv(i,j),dp)*sec2year
        dzs_dtau(j,i)  = real(dzs_dtau_conv(i,j),dp)*sec2year
        dzm_dtau(j,i)  = real(dzm_dtau_conv(i,j),dp)*sec2year
        dzb_dtau(j,i)  = real(dzb_dtau_conv(i,j),dp)*sec2year
        dzl_dtau(j,i)  = real(dzl_dtau_conv(i,j),dp)*sec2year
        dH_dtau(j,i)   = real(dH_dtau_conv(i,j),dp)*sec2year
        dH_c_dtau(j,i) = real(dH_c_dtau_conv(i,j),dp)*sec2year
        dH_t_dtau(j,i) = real(dH_t_dtau_conv(i,j),dp)*sec2year
        vx_b_g(j,i)  = real(vx_b_g_conv(i,j),dp)*sec2year
        vy_b_g(j,i)  = real(vy_b_g_conv(i,j),dp)*sec2year
        vz_b(j,i)    = real(vz_b_conv(i,j),dp)*sec2year
        vx_s_g(j,i)  = real(vx_s_g_conv(i,j),dp)*sec2year
        vy_s_g(j,i)  = real(vy_s_g_conv(i,j),dp)*sec2year
        vz_s(j,i)    = real(vz_s_conv(i,j),dp)*sec2year
        temp_b(j,i)  = real(temp_b_conv(i,j),dp)
        temph_b(j,i) = real(temph_b_conv(i,j),dp)
        p_b_w(j,i)   = real(p_b_w_conv(i,j),dp)
        q_w(j,i)     = real(q_w_conv(i,j),dp)*sec2year
        q_w_x(j,i)   = real(q_w_x_conv(i,j),dp)*sec2year
        q_w_y(j,i)   = real(q_w_y_conv(i,j),dp)*sec2year
        H_w(j,i)     = real(H_w_conv(i,j),dp)
        ratio_sl_sia_x(j,i) = real(ratio_sl_sia_x_conv(i,j),dp)
        ratio_sl_sia_y(j,i) = real(ratio_sl_sia_y_conv(i,j),dp)
        ratio_sl_sia(j,i)   = real(ratio_sl_sia_conv(i,j),dp)

        if (flag_shelfy_stream_x_conv(i,j) == 1) then
           flag_shelfy_stream_x(j,i) = .true.
        else
           flag_shelfy_stream_x(j,i) = .false.
        end if

        if (flag_shelfy_stream_y_conv(i,j) == 1) then
           flag_shelfy_stream_y(j,i) = .true.
        else
           flag_shelfy_stream_y(j,i) = .false.
        end if

        if (flag_shelfy_stream_conv(i,j) == 1) then
           flag_shelfy_stream(j,i) = .true.
        else
           flag_shelfy_stream(j,i) = .false.
        end if

        if (flag_grounding_line_1_conv(i,j) == 1) then
           flag_grounding_line_1(j,i) = .true.
        else
           flag_grounding_line_1(j,i) = .false.
        end if

        if (flag_grounding_line_2_conv(i,j) == 1) then
           flag_grounding_line_2(j,i) = .true.
        else
           flag_grounding_line_2(j,i) = .false.
        end if

        if (flag_calving_front_1_conv(i,j) == 1) then
           flag_calving_front_1(j,i) = .true.
        else
           flag_calving_front_1(j,i) = .false.
        end if

        if (flag_calving_front_2_conv(i,j) == 1) then
           flag_calving_front_2(j,i) = .true.
        else
           flag_calving_front_2(j,i) = .false.
        end if

        if (flag_grounded_front_a_1_conv(i,j) == 1) then
           flag_grounded_front_a_1(j,i) = .true.
        else
           flag_grounded_front_a_1(j,i) = .false.
        end if

        if (flag_grounded_front_a_2_conv(i,j) == 1) then
           flag_grounded_front_a_2(j,i) = .true.
        else
           flag_grounded_front_a_2(j,i) = .false.
        end if

        if (flag_grounded_front_b_1_conv(i,j) == 1) then
           flag_grounded_front_b_1(j,i) = .true.
        else
           flag_grounded_front_b_1(j,i) = .false.
        end if

        if (flag_grounded_front_b_2_conv(i,j) == 1) then
           flag_grounded_front_b_2(j,i) = .true.
        else
           flag_grounded_front_b_2(j,i) = .false.
        end if

        vis_ave_g(j,i)  = real(vis_ave_g_conv(i,j),dp)
        vis_int_g(j,i)  = real(vis_int_g_conv(i,j),dp)

        do kr=0, KRMAX
           temp_r(kr,j,i) = real(temp_r_conv(i,j,kr),dp)
        end do

        do kt=0, KTMAX
           vx_t(kt,j,i)    = real(vx_t_conv(i,j,kt),dp)*sec2year
           vy_t(kt,j,i)    = real(vy_t_conv(i,j,kt),dp)*sec2year
           vz_t(kt,j,i)    = real(vz_t_conv(i,j,kt),dp)*sec2year
           omega_t(kt,j,i) = real(omega_t_conv(i,j,kt),dp)
           age_t(kt,j,i)   = real(age_t_conv(i,j,kt),dp)*year2sec
           enth_t(kt,j,i)  = real(enth_t_conv(i,j,kt),dp)
           enh_t(kt,j,i)   = real(enh_t_conv(i,j,kt),dp)
           strain_heating_t(kt,j,i) = real(strain_heating_t_conv(i,j,kt),dp)
        end do

        do kc=0, KCMAX
           vx_c(kc,j,i)    = real(vx_c_conv(i,j,kc),dp)*sec2year
           vy_c(kc,j,i)    = real(vy_c_conv(i,j,kc),dp)*sec2year
           vz_c(kc,j,i)    = real(vz_c_conv(i,j,kc),dp)*sec2year
           temp_c(kc,j,i)  = real(temp_c_conv(i,j,kc),dp)
           age_c(kc,j,i)   = real(age_c_conv(i,j,kc),dp)*year2sec
           enth_c(kc,j,i)  = real(enth_c_conv(i,j,kc),dp)
           omega_c(kc,j,i) = real(omega_c_conv(i,j,kc),dp)
           enh_c(kc,j,i)   = real(enh_c_conv(i,j,kc),dp)
           strain_heating_c(kc,j,i) = real(strain_heating_c_conv(i,j,kc),dp)
        end do

     end do
     end do

     if (.not.flag_ratio_sl_sia) then
              ! reconstruct ratio_sl_sia from ratio_sl_sia_x/y

        ratio_sl_sia = 0.0_dp

        do i=1, IMAX-1
        do j=1, JMAX-1

           if (mask(j,i) == 0) &   ! grounded ice
              ratio_sl_sia(j,i) &
                 = 0.25_dp &
                      * (   ratio_sl_sia_x(j,i-1) + ratio_sl_sia_x(j,i) &
                          + ratio_sl_sia_y(j-1,i) + ratio_sl_sia_y(j,i) )
        end do
        end do

     end if

     if (.not.flag_vis_ave_g) then   ! reconstruct vis_ave_g from vis_int_g

#if (defined(VISC_MIN) && defined(VISC_MAX))
        visc_min = VISC_MIN
        visc_max = VISC_MAX
#else
        visc_min = 1.0e+10_dp   ! Pa s
        visc_max = 1.0e+25_dp   ! Pa s
#endif

#if (defined(VISC_INIT_SSA))
        visc_init = VISC_INIT_SSA
#else
        visc_init = 1.0e+15_dp   ! Pa s
#endif

        do i=0, IMAX
        do j=0, JMAX

           if ((mask(j,i)==0).or.(mask(j,i)==3)) then
                  ! grounded or floating ice

              vis_ave_g(j,i) = vis_int_g(j,i)/max(H(j,i), eps_dp)

              vis_ave_g(j,i) = max(min(vis_ave_g(j,i), visc_max), visc_min)

           else   ! (mask(j,i)==1).or.(mask(j,i)==2),
                  ! ice-free land or ocean

              vis_ave_g(j,i) = visc_init   ! dummy value

           end if

        end do
        end do

     end if   ! (.not.flag_vis_ave_g)

  else   ! flag_temp_age_only == .true.

     if ( present(opt_mask) &
          .and.present(opt_n_cts) &
          .and.present(opt_kc_cts) &
          .and.present(opt_H_cold) &
          .and.present(opt_H_temp) &
          .and.present(opt_H) &
          .and.present(opt_temp_r) &
          .and.present(opt_omega_t) &
          .and.present(opt_age_t) &
          .and.present(opt_temp_c) &
          .and.present(opt_age_c) &
          .and.present(opt_omega_c) &
        ) then

        do i=0, IMAX
        do j=0, JMAX

           opt_mask(j,i)  = mask_conv(i,j)
           opt_n_cts(j,i)  = n_cts_conv(i,j)
           opt_kc_cts(j,i) = kc_cts_conv(i,j)

           opt_H_cold(j,i) = real(H_cold_conv(i,j),dp)
           opt_H_temp(j,i) = real(H_temp_conv(i,j),dp)
           opt_H(j,i)      = real(H_conv(i,j),dp)

           do kr=0, KRMAX
              opt_temp_r(kr,j,i) = real(temp_r_conv(i,j,kr),dp)
           end do

           do kt=0, KTMAX
              opt_omega_t(kt,j,i) = real(omega_t_conv(i,j,kt),dp)
              opt_age_t(kt,j,i)   = real(age_t_conv(i,j,kt),dp)*year2sec
           end do

           do kc=0, KCMAX
              opt_temp_c(kc,j,i)  = real(temp_c_conv(i,j,kc),dp)
              opt_age_c(kc,j,i)   = real(age_c_conv(i,j,kc),dp)*year2sec
              opt_omega_c(kc,j,i) = real(omega_c_conv(i,j,kc),dp)
           end do

        end do
        end do

     else

        errormsg = ' >>> read_tms_nc: optional argument(s) missing!'
        call error(errormsg)

     end if

  end if   ! flag_temp_age_only == .false./.true.

  end subroutine read_tms_nc

!-------------------------------------------------------------------------------
!> Reading of the target-topography file (in NetCDF format).
!-------------------------------------------------------------------------------
  subroutine read_target_topo_nc(target_topo_dat_name)

  use sico_variables_m, only : mask_target, &
                               zs_target, zb_target, zl_target, &
                               H_target, errormsg, end_of_line

  use netcdf
  use nc_check_m

  implicit none

  character(len=256), intent(in) :: target_topo_dat_name

! Return variables
! (defined as global variables in module sico_variables_m):
!
!    mask_target, zs_target, zb_target, zl_target, H_target

  integer(i4b)                           :: ios, istat, istat2
  integer(i4b), dimension(0:IMAX,0:JMAX) :: mask_conv
  real(sp), dimension(0:IMAX,0:JMAX)     :: zs_conv, zb_conv, zl_conv, H_conv
  character(len=256)                     :: target_topo_dat_path
  character(len=256)                     :: filename_with_path

  integer(i4b) :: ncid, ncv
  !     ncid:      ID of the output file
  !     ncv:       Variable ID

  character(len=64), parameter :: thisroutine = 'read_target_topo_nc'

#if defined(ALLOW_TAPENADE) /* Tapenade */
  integer(i4b) :: i, j
#endif /* Tapenade */

!-------- Read target-topography data
!         (from a time-slice file of a previous simulation) --------

#if (defined(TARGET_TOPO_PATH))
  target_topo_dat_path = trim(TARGET_TOPO_PATH)
#else
  target_topo_dat_path = 'dummy'
  errormsg = ' >>> read_target_topo_nc: TARGET_TOPO_PATH must be defined!'
  call error(errormsg)
#endif

  filename_with_path &
     = trim(target_topo_dat_path)//'/'//trim(target_topo_dat_name)

  ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

  if (ios /= nf90_noerr) then
     errormsg = ' >>> read_target_topo_nc: Error when opening the' &
              //                           end_of_line &
              //'                          target-topography NetCDF file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'mask', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, mask_conv), thisroutine )
  else
     istat2 = nf90_inq_varid(ncid, 'maske', ncv)
     if (istat2 == nf90_noerr) then
        call check( nf90_get_var(ncid, ncv, mask_conv), thisroutine )
     else
        errormsg = ' >>> read_target_topo_nc: Variable ''mask'' ' &
                 //                           end_of_line &
                 //'                          not available in read nc file!'
        call error(errormsg)
     end if
  end if

  istat = nf90_inq_varid(ncid, 'zs', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, zs_conv), thisroutine )
  else
     errormsg = ' >>> read_target_topo_nc: Variable ''zs'' ' &
              //                           end_of_line &
              //'                          not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'zb', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, zb_conv), thisroutine )
  else
     errormsg = ' >>> read_target_topo_nc: Variable ''zb'' ' &
              //                           end_of_line &
              //'                          not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'zl', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, zl_conv), thisroutine )
  else
     errormsg = ' >>> read_target_topo_nc: Variable ''zl'' ' &
              //                           end_of_line &
              //'                          not available in read nc file!'
     call error(errormsg)
  end if

  istat = nf90_inq_varid(ncid, 'H', ncv)
  if (istat == nf90_noerr) then
     call check( nf90_get_var(ncid, ncv, H_conv), thisroutine )
  else
     errormsg = ' >>> read_target_topo_nc: Variable ''H'' ' &
              //                           end_of_line &
              //'                          not available in read nc file!'
     call error(errormsg)
  end if

  call check( nf90_close(ncid), thisroutine )

!-------- Convert data to double precision --------

#if !defined(ALLOW_TAPENADE) /* Normal */
  mask_target = transpose(mask_conv)
  zs_target    = real(transpose(zs_conv),dp)
  zb_target    = real(transpose(zb_conv),dp)
  zl_target    = real(transpose(zl_conv),dp)
  H_target     = real(transpose(H_conv) ,dp)
#else /* Tapenade */
  do i=0, IMAX
  do j=0, JMAX
     mask_target(j,i) = mask_conv(i,j)
     zs_target(j,i)    = real(zs_conv(i,j),dp)
     zb_target(j,i)    = real(zb_conv(i,j),dp)
     zl_target(j,i)    = real(zl_conv(i,j),dp)
     H_target(j,i)     = real(H_conv(i,j) ,dp)
  end do
  end do
#endif /* Normal vs. Tapenade */

  end subroutine read_target_topo_nc

!-------------------------------------------------------------------------------
!> Reading of scalar, time-dependent input files in NetCDF or ASCII format.
!-------------------------------------------------------------------------------
  subroutine read_scalar_input(filename_with_path, ch_var_name, ndata_max, &
                               n_time_min, n_time_stp, n_time_max, ndata, &
                               scalar_data)

  use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

  use netcdf
  use nc_check_m

  implicit none

  character(len=256), intent(in) :: filename_with_path
  character(len=  *), intent(in) :: ch_var_name
  integer(i4b)      , intent(in) :: ndata_max

  integer(i4b), intent(out) :: n_time_min, n_time_stp, n_time_max
  integer(i4b), intent(out) :: ndata
  real(dp)    , intent(out), dimension(0:ndata_max) :: scalar_data

  integer(i4b)       :: n
  integer(i4b)       :: ios
  character(len=256) :: filename_aux
  character(len=  3) :: ch_nc_test
  character          :: ch_dummy
  real(dp)           :: d_dummy
  logical            :: flag_nc

  integer(i4b) :: ncid, ncv
  !     ncid:      ID of the output file
  !     ncv:       Variable ID

!-------- Determining file type --------

  filename_aux = adjustr(filename_with_path)
  n            = len(filename_aux)
  ch_nc_test   = filename_aux(n-2:n)
  filename_aux = adjustl(filename_aux)

  if (ch_nc_test == '.nc') then
     flag_nc = .true.   ! NetCDF file
  else
     flag_nc = .false.  ! ASCII file
  end if

!-------- Reading file --------

  if (flag_nc) then   ! NetCDF file

     errormsg = ' >>> read_scalar_input: Reading NetCDF not yet implemented!'
     call error(errormsg)

     ios = nf90_open(trim(filename_aux), NF90_NOWRITE, ncid)

     ! ...

     call check( nf90_close(ncid) )

  else   ! ASCII file

     open(21, iostat=ios, file=trim(filename_with_path), status='old')

     if (ios /= 0) then
        errormsg = ' >>> read_scalar_input: Error when opening the ' &
                         // trim(adjustl(ch_var_name)) // ' ASCII file!'
        call error(errormsg)
     end if

     read(21, fmt=*) ch_dummy, n_time_min, n_time_stp, n_time_max

     if (ch_dummy /= '#') then
        errormsg = ' >>> read_scalar_input: ' &
                 //         'n_time_min, n_time_stp, n_time_max' &
                 //         end_of_line &
                 //'        not defined in ' &
                 //         trim(adjustl(ch_var_name)) // ' ASCII file!'
        call error(errormsg)
     end if

     ndata = (n_time_max-n_time_min)/n_time_stp

     if (ndata > ndata_max) then
        errormsg = ' >>> read_scalar_input: ' &
                 //         'ndata <= ndata_max required!' &
                 //         end_of_line &
                 //'        Increase corresponding ndata_max value for ' &
                 //         trim(adjustl(ch_var_name)) &
                 //         end_of_line &
                 //'        in sico_variables_m or sico_vars_m!'
        call error(errormsg)
     end if

     scalar_data = 0.0_dp   ! initialization

     do n=0, ndata
        read(21, fmt=*) d_dummy, scalar_data(n)
     end do

     close(21, status='keep')

  end if

  end subroutine read_scalar_input

!-------------------------------------------------------------------------------
!> Reading of 2D input files in NetCDF or ASCII format.
!-------------------------------------------------------------------------------
  subroutine read_2d_input(filename_with_path, ch_var_name, n_var_type, &
                           n_ascii_header, field2d_r)

  use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

  use netcdf
  use nc_check_m

  implicit none

  character(len=256), intent(in)                  :: filename_with_path
  character(len=  *), intent(in)                  :: ch_var_name
  integer(i4b), intent(in)                        :: n_var_type
  integer(i4b), intent(in)                        :: n_ascii_header
  real(dp), dimension(0:JMAX,0:IMAX), intent(out) :: field2d_r

  integer(i4b)       :: i, j
  integer(i4b)       :: n
  integer(i4b)       :: ios, ierr
  character(len=256) :: ch_var_name_aux
  character(len=256) :: filename_aux
  character(len=  3) :: ch_nc_test
  character(len=  3) :: ch_month
  character          :: ch_dummy
  logical            :: flag_nc
  logical            :: flag_nc_mm

  integer(i4b), dimension(0:IMAX,0:JMAX)    :: mask_aux_conv
  integer(i4b), dimension(0:IMAX,0:JMAX)    :: n_aux_conv
  real(dp)    , dimension(0:IMAX,0:JMAX)    :: r_aux_conv
  real(dp)    , dimension(0:IMAX,0:JMAX,12) :: r_aux_conv_mm

  integer(i4b) :: ncid, ncv
  !     ncid:      ID of the output file
  !     ncv:       Variable ID

  character(len=  8) :: ch_imax
  character(len=128) :: fmt4

  write(ch_imax, fmt='(i8)') IMAX
  write(fmt4,    fmt='(a)')  '('//trim(adjustl(ch_imax))//'(i1),i1)'

  ch_var_name_aux = trim(adjustl(ch_var_name))

  flag_nc_mm = .false.

!-------- Determining file type --------

  filename_aux = adjustr(filename_with_path)
  n            = len(filename_aux)
  ch_nc_test   = filename_aux(n-2:n)
  filename_aux = adjustl(filename_aux)

  if (ch_nc_test == '.nc') then
     flag_nc = .true.   ! NetCDF file
  else
     flag_nc = .false.  ! ASCII file
  end if

!-------- Reading file --------

  if (flag_nc) then   ! NetCDF file

     ios = nf90_open(trim(filename_aux), NF90_NOWRITE, ncid)

     if (ios /= nf90_noerr) then
        errormsg = ' >>> read_2d_input: Error when opening the ' &
                         // trim(ch_var_name_aux) // ' NetCDF file!'
        call error(errormsg)
     end if

     ierr = nf90_inq_varid(ncid, trim(ch_var_name_aux), ncv)

     if ((ierr /= nf90_noerr).and.(n_var_type==0)) then

        ch_var_name_aux = adjustr(ch_var_name_aux)
        n               = len(ch_var_name_aux)
        ch_month        = ch_var_name_aux(n-2:n)
        ch_var_name_aux = ch_var_name_aux(1:n-4)
        ch_var_name_aux = adjustl(ch_var_name_aux)
        flag_nc_mm      = .true.

        ierr = nf90_inq_varid(ncid, trim(ch_var_name_aux), ncv)
                  ! Trying again without the month code ('_jan', '_feb' etc.)
     end if

     if (ierr /= nf90_noerr) then
        errormsg = ' >>> read_2d_input: NetCDF variable ' &
                         // trim(ch_var_name_aux) // ' not present!'
        call error(errormsg)
     end if

     if ((n_var_type==0).and.(flag_nc_mm)) then
        call check( nf90_get_var(ncid, ncv, r_aux_conv_mm) )
     else if ((n_var_type==0).and.(.not.flag_nc_mm)) then
        call check( nf90_get_var(ncid, ncv, r_aux_conv) )
     else if (n_var_type==1) then
        call check( nf90_get_var(ncid, ncv, r_aux_conv) )
     else if (n_var_type==2) then
        call check( nf90_get_var(ncid, ncv, n_aux_conv) )
     else if (n_var_type==3) then
        call check( nf90_get_var(ncid, ncv, mask_aux_conv) )
     else
        errormsg = ' >>> read_2d_input: n_var_type must be between 0 and 3!'
        call error(errormsg)
     end if

     call check( nf90_close(ncid) )

  else   ! ASCII file

     if ((n_var_type==0).or.(n_var_type==1)) then
        open(21, iostat=ios, file=trim(filename_aux), recl=rcl1, status='old')
     else if ((n_var_type==2).or.(n_var_type==3)) then
        open(21, iostat=ios, file=trim(filename_aux), recl=rcl2, status='old')
     else
        errormsg = ' >>> read_2d_input: n_var_type must be between 0 and 3!'
        call error(errormsg)
     end if

     if (ios /= 0) then
        errormsg = ' >>> read_2d_input: Error when opening the ' &
                         // trim(ch_var_name_aux) // ' ASCII file!'
        call error(errormsg)
     end if

     do n=1, n_ascii_header; read(21, fmt='(a)') ch_dummy; end do

     do j=JMAX, 0, -1

        if ((n_var_type==0).or.(n_var_type==1)) then
           read(21, fmt=*) (r_aux_conv(i,j), i=0,IMAX)
        else if (n_var_type==2) then
           read(21, fmt=*) (n_aux_conv(i,j), i=0,IMAX)
        else if (n_var_type==3) then
           read(21, fmt=trim(fmt4)) (mask_aux_conv(i,j), i=0,IMAX)
        else
           errormsg = ' >>> read_2d_input: n_var_type must be between 0 and 3!'
           call error(errormsg)
        end if

     end do

     close(21, status='keep')

  end if

!-------- Converting read field --------

  do i=0, IMAX
  do j=0, JMAX

     if ((n_var_type==0).and.(flag_nc_mm)) then
        if (ch_month == 'jan') then
           field2d_r(j,i) = r_aux_conv_mm(i,j, 1)
        else if (ch_month == 'feb') then
           field2d_r(j,i) = r_aux_conv_mm(i,j, 2)
        else if (ch_month == 'mar') then
           field2d_r(j,i) = r_aux_conv_mm(i,j, 3)
        else if (ch_month == 'apr') then
           field2d_r(j,i) = r_aux_conv_mm(i,j, 4)
        else if (ch_month == 'may') then
           field2d_r(j,i) = r_aux_conv_mm(i,j, 5)
        else if (ch_month == 'jun') then
           field2d_r(j,i) = r_aux_conv_mm(i,j, 6)
        else if (ch_month == 'jul') then
           field2d_r(j,i) = r_aux_conv_mm(i,j, 7)
        else if (ch_month == 'aug') then
           field2d_r(j,i) = r_aux_conv_mm(i,j, 8)
        else if (ch_month == 'sep') then
           field2d_r(j,i) = r_aux_conv_mm(i,j, 9)
        else if (ch_month == 'oct') then
           field2d_r(j,i) = r_aux_conv_mm(i,j,10)
        else if (ch_month == 'nov') then
           field2d_r(j,i) = r_aux_conv_mm(i,j,11)
        else if (ch_month == 'dec') then
           field2d_r(j,i) = r_aux_conv_mm(i,j,12)
        else
           errormsg = ' >>> read_2d_input: Month code ' // ch_month &
                    //         end_of_line &
                    //'        invalid when trying to convert variable ' &
                    //         end_of_line &
                    //'        ' // trim(ch_var_name_aux) // '!'
           call error(errormsg)
        end if
     else if ((n_var_type==0).and.(.not.flag_nc_mm)) then
        field2d_r(j,i) = r_aux_conv(i,j)
     else if (n_var_type==1) then
        field2d_r(j,i) = r_aux_conv(i,j)
     else if (n_var_type==2) then
        field2d_r(j,i) = real(n_aux_conv(i,j),dp)
     else if (n_var_type==3) then
        field2d_r(j,i) = real(mask_aux_conv(i,j),dp)
     else
        errormsg = ' >>> read_2d_input: n_var_type must be between 0 and 3!'
        call error(errormsg)
     end if

  end do
  end do

  end subroutine read_2d_input

!-------------------------------------------------------------------------------
!> Reading of the tabulated kei function.
!-------------------------------------------------------------------------------
  subroutine read_kei()

  use sico_variables_m, only : kei, kei_r_max, kei_r_incr, n_data_kei, errormsg

  implicit none

  integer(i4b)       :: n, n_data_kei_max
  integer(i4b)       :: ios
  real(dp)           :: r_val, d_dummy
  character          :: ch_dummy
  character(len=256) :: filename_with_path

  n_data_kei_max = 10000

  kei = 0.0_dp

!-------- Reading of tabulated values --------

  filename_with_path = trim(IN_PATH)//'/general/kei.dat'

#if !defined(ALLOW_TAPENADE) /* Normal */
  open(unit=11, iostat=ios, file=trim(filename_with_path), status='old')
#else /* Tapenade */
  open(unit=11, iostat=ios, file=trim(filename_with_path))
#endif /* Normal vs. Tapenade */

  if (ios /= 0) then
     errormsg = ' >>> read_kei: Error when opening the kei file!'
     call error(errormsg)
  end if

  read(unit=11, fmt='(a)') ch_dummy
  read(unit=11, fmt='(15x,f5.0)') kei_r_max
  read(unit=11, fmt='(15x,f5.0)') kei_r_incr
  read(unit=11, fmt='(a)') ch_dummy

  n_data_kei = nint(kei_r_max/kei_r_incr)

  if (n_data_kei > n_data_kei_max) then
     errormsg = ' >>> read_kei: Array kei too small!'
     call error(errormsg)
  end if

  read(unit=11, fmt='(a)') ch_dummy
  read(unit=11, fmt='(a)') ch_dummy
  read(unit=11, fmt='(a)') ch_dummy

  do n=-n_data_kei, n_data_kei
     read(unit=11, fmt=*) d_dummy, kei(n)
  end do

  close(unit=11, status='keep')

  end subroutine read_kei

!-------------------------------------------------------------------------------
!> Reading of physical parameters
!! (tabulated values of the rate factor, heat conductivity and specific heat).
!-------------------------------------------------------------------------------
  subroutine read_phys_para()

  use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

  use netcdf
  use nc_check_m

  implicit none

  integer(i4b), parameter :: n_unit=31
  integer(i4b) :: ios
  integer(i4b) :: n
  character(len=256) :: filename_with_path
  character(len=256) :: filename_aux
  character(len=  3) :: ch_nc_test
  logical            :: flag_nc

  integer(i4b) :: ncid, ncv
  !     ncid:      ID of the output file
  !     ncv:       Variable ID

  character(len=64), parameter :: thisroutine = 'read_phys_para'

!-------- Determining file type --------

#if (defined(RF_KAPPA_C_FILE))
  filename_with_path = trim(IN_PATH)//'/general/'// &
                       trim(RF_KAPPA_C_FILE)
#endif

  filename_aux = adjustr(filename_with_path)
  n            = len(filename_aux)
  ch_nc_test   = filename_aux(n-2:n)
  filename_aux = adjustl(filename_aux)

  if (ch_nc_test == '.nc') then
     flag_nc = .true.   ! NetCDF file
  else
     flag_nc = .false.  ! ASCII file
  end if

!-------- Reading of parameters from file --------

!  ------ Opening file

  if (flag_nc) then   ! NetCDF file

     ios = nf90_open(trim(filename_aux), NF90_NOWRITE, ncid)

     if (ios /= nf90_noerr) then
        errormsg = ' >>> read_phys_para: ' &
                         // 'Error when opening the RF_KAPPA_C NetCDF file!'
        call error(errormsg)
     end if

  else   ! ASCII file

     open(n_unit, iostat=ios, file=trim(filename_aux))

     if (ios /= 0) then
        errormsg = ' >>> read_phys_para: ' &
                         // 'Error when opening the RF_KAPPA_C ASCII file!'
        call error(errormsg)
     end if

  end if

!  ------ Reading the rate factor, heat conductivity and specific heat

  if (flag_nc) then
     call check( nf90_inq_varid(ncid, 'RF', ncv), thisroutine )
     call check( nf90_get_var(ncid, ncv, RF), thisroutine )
  else
     do n=10, -190, -1
        call read_phys_para_value(n_unit, 'RF(.)', RF(n))
     end do
  end if

  if (flag_nc) then
     call check( nf90_inq_varid(ncid, 'KAPPA', ncv), thisroutine )
     call check( nf90_get_var(ncid, ncv, KAPPA), thisroutine )
  else
     do n=10, -190, -1
        call read_phys_para_value(n_unit, 'KAPPA(.)', KAPPA(n))
     end do
  end if

  if (flag_nc) then
     call check( nf90_inq_varid(ncid, 'C', ncv), thisroutine )
     call check( nf90_get_var(ncid, ncv, C), thisroutine )
  else
     do n=10, -190, -1
        call read_phys_para_value(n_unit, 'C(.)', C(n))
     end do
  end if

!  ------ Closing file

  if (flag_nc) then   ! NetCDF file
     call check( nf90_close(ncid), thisroutine )
  else
     close(n_unit, status='keep')
  end if

  end subroutine read_phys_para

!-------------------------------------------------------------------------------
!> Reading of a value of a physical parameter from the input file.
!-------------------------------------------------------------------------------
  subroutine read_phys_para_value(n_unit, ch_para, d_para)

  use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

  implicit none

  integer(i4b)    , intent(in) :: n_unit
  character(len=*), intent(in) :: ch_para

  real(dp), intent(out) :: d_para

  integer            :: lch, n_equals
  character(len=256) :: ch_line
  character(len= 64) :: ch_para_1a, ch_para_1b, ch_para_2a, ch_para_2b
  character(len=  8) :: ch_equals
  character          :: ch1, ch2

  lch        = len(ch_para)
  ch_para_1a = ch_para

  ch1 = '%'

  do while (ch1 == '%')

     read(unit=n_unit, fmt='(a)') ch_line

     ch1 = ch_line(1:1)

     if (ch1 /= '%') then   ! no comment line

#if !defined(ALLOW_TAPENADE) /* Normal */
        n_equals = index(ch_line, '=')
#else /* Tapenade: cannot differentiate through index function */
        n_equals = 15   ! assuming the equals sign to be
                        ! always in row 15 of the input file
#endif /* Normal vs. Tapenade */

        ch_para_2a = ch_line(1:n_equals-1)
        ch2        = ch_line(n_equals:n_equals)

        if (ch2 == '=') then
           read(ch_line(n_equals+1:), fmt=*) d_para
        else
           write(ch_equals, '(i0)') n_equals
           errormsg = ' >>> read_phys_para_value:' &
                    //         end_of_line &
                    //'        Trying to read ' // trim(ch_para_1a) &
                    //         ' from the RF_KAPPA_C file,' &
                    //         end_of_line &
                    //'        but equals sign (=) is not in row ' &
                    //         trim(ch_equals)//'!'
           call error(errormsg)
        end if

     end if

  end do

  ch_para_1a = adjustl(ch_para_1a)
  ch_para_2a = adjustl(ch_para_2a)

  if (ch_para(lch-2:lch) == '(.)') then
     ch_para_1b = ch_para_1a(1:lch-3)
     ch_para_2b = ch_para_2a(1:lch-3)
  else
     ch_para_1b = ch_para_1a
     ch_para_2b = ch_para_2a
  end if

  if (trim(ch_para_1b) /= trim(ch_para_2b)) then
     errormsg = ' >>> read_phys_para_value:' &
              //         end_of_line &
              //'        Trying to read '//trim(ch_para_1a)//', ' &
              //         end_of_line &
              //'        but found ' &
              //         trim(ch_para_2a)//' in the RF_KAPPA_C file!'
     call error(errormsg)
  end if

  end subroutine read_phys_para_value

#if (defined(ALLOW_NORMAL) || defined(ALLOW_GRDCHK) || defined(ALLOW_TAPENADE))
!-------------------------------------------------------------------------------
!> Reading in of adjoint related data (only set up for ages right now).
!-------------------------------------------------------------------------------
  subroutine read_age_data()

  use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

  implicit none

    integer(i4b) :: ios, i, j, k, kc, kt, KDATA
    logical      :: debug = .true.
! path stores path with filename
! assume you are in src directory when defining path
    character(len=256) :: path 

#if (defined(AGE_COST)) /* Only designed for AGE_COST right now */
! Warning: code below is not robust against
! CALCMOD.NE.3,2
#if (CALCMOD==3)
    KDATA = KCMAX
#else
    KDATA = KCMAX + KTMAX
#endif

#if (defined(AGE_DATA_PATH))
  path = trim(AGE_DATA_PATH)
#else
  errormsg = ' >>> read_tms_nc: AGE_DATA_PATH must be defined!'
  call error(errormsg)
#endif

open(unit=1092, file=path, status='old')

! Some issues need to be sorted with age_unc data first

!#if (defined(AGE_UNC_PATH))
!  path = trim(AGE_UNC_PATH)
!#else
!  errormsg = ' >>> read_tms_nc: AGE_UNC_PATH must be defined!'
!  call error(errormsg)
!#endif

!open(unit=1093, file=path, status='old')

#if (defined(H_DATA_PATH))
  path = trim(H_DATA_PATH)
#else
  errormsg = ' >>> read_tms_nc: H_DATA_PATH must be defined!'
  call error(errormsg)
#endif

open(unit=1094, file=path, status='old')

! Age layer data only has H_data, but not H_unc
! Will have to consider using BedMachine data for H_unc

!#if (defined(H_UNC_PATH))
!  path = trim(H_UNC_PATH)
!#else
!  errormsg = ' >>> read_tms_nc: H_UNC_PATH must be defined!'
!  call error(errormsg)
!#endif

!open(unit=1095, file=path, status='old')

    do i=0,IMAX
       do j=0,JMAX
          read(1094, *) H_data(j,i)
          !read(1095, *) H_unc(j,i)
       end do
    end do

    do k=0,KDATA
       do i=0,IMAX
          do j=0,JMAX
             read(1092, *) age_data(k,j,i)
             !read(1093, *) age_unc(k,j,i)
          end do
       end do
    end do

    ! ages a -> s
    do k=0,KDATA
       do j=0,JMAX
          do i=0,IMAX
                age_data(k,j,i) = age_data(k,j,i) * year2sec
                !age_unc(k,j,i)  = age_unc(k,j,i)  * year2sec
          end do
       end do
    end do
    close(unit=1092)
    !close(unit=1093)
    close(unit=1094)
    !close(unit=1095)
#endif

  end subroutine read_age_data

  subroutine read_BedMachine_data()

  use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

  implicit none

    integer(i4b) :: i, j
! path stores path with filename
! assume you are in src directory when defining path
    character(len=256) :: path 

#if (defined(BEDMACHINE_COST)) 

#if (defined(BEDMACHINE_PATH))
  path = trim(BEDMACHINE_PATH)
#else
  errormsg = ' >>> read_tms_nc: BEDMACHINE_PATH must be defined!'
  call error(errormsg)
#endif

open(unit=92, file=path, status='old')

#if (defined(BEDMACHINE_UNC_PATH))
  path = trim(BEDMACHINE_UNC_PATH)
#else
  errormsg = ' >>> read_tms_nc: BEDMACHINE_UNC_PATH must be defined!'
  call error(errormsg)
#endif

open(unit=93, file=path, status='old')

       do i=0,IMAX
          do j=0,JMAX
             read(92, *) H_BedMachine_data(j,i)
             read(93, *) H_unc_BedMachine_data(j,i)
          end do
       end do
    close(unit=92)
    close(unit=93)

#endif
  end subroutine read_BedMachine_data
#endif /* inclusion of Tapenade only routine */

!-------------------------------------------------------------------------------

end module read_m
!
