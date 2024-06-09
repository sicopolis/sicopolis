!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  o u t p u t _ m
!
!! Writing of output data on files.
!!
!!##### Authors
!!
!! Ralf Greve, Reinhard Calov, Thomas Goelles, Thorben Dunse
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
!> Writing of output data on files.
!-------------------------------------------------------------------------------
module output_m

  use sico_types_m
  use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

  use error_m

  implicit none

  private
  public :: output1, output2, output4, borehole
#if (defined(ASF) && WRITE_SER_FILE_STAKES==1) /* Austfonna */
  public :: output5
#endif

contains

!-------------------------------------------------------------------------------
!> Writing of time-slice files in NetCDF format.
!-------------------------------------------------------------------------------
subroutine output1(time, delta_ts, glac_index, &
                   flag_3d_output, ndat2d, ndat3d, &
                   opt_flag_compute_flux_vars_only)

#if (CALCMOD==1 || CALCMOD==0 || CALCMOD==-1)
  use enth_temp_omega_m, only : enth_fct_temp_omega
#endif

  use netcdf
  use nc_check_m

#if (DISC>0)
  use discharge_workers_m, only: dis_perp, cst_dist, cos_grad_tc, mask_mar
#endif

implicit none

real(dp),           intent(in) :: time, delta_ts, glac_index
logical,            intent(in) :: flag_3d_output

logical, optional,  intent(in) :: opt_flag_compute_flux_vars_only

integer(i4b),    intent(inout) :: ndat2d, ndat3d

integer(i4b) :: i, j, kc, kt, kr
integer(i4b) :: ios
integer(i4b) :: ndat
real(dp), dimension(0:JMAX,0:IMAX) :: H_cold, H_temp
real(dp), dimension(0:JMAX,0:IMAX) :: vx_m_g, vy_m_g
real(dp) :: V_tot, V_grounded, V_floating, V_gr_redu, V_af
real(dp) :: A_grounded, A_floating
real(dp) :: lond0, latd0
real(dp) :: rhosw_rho_ratio
character(len=256) :: filename, filename_with_path
character(len= 16) :: ch_date, ch_time, ch_zone
character(len=  4) :: ch_ndat
logical :: flag_compute_flux_vars_only

real(dp)         , parameter :: one_year = 1.0_dp
character(len=64), parameter :: thisroutine = 'output1'

#if (OUTPUT_FLUX_VARS==2)   /* averaging of flux variables */

integer(i4b), save :: n_flx_ave_cnt = 0
real(dp)           :: r_n_flx_ave_cnt_inv

real(dp), dimension(0:JMAX,0:IMAX), save :: accum_sum         = 0.0_dp, &
                                            snowfall_sum      = 0.0_dp, &
                                            rainfall_sum      = 0.0_dp, &
                                            as_perp_sum       = 0.0_dp, &
                                            as_perp_apl_sum   = 0.0_dp, &
                                            smb_corr_sum      = 0.0_dp, &
                                            mb_source_apl_sum = 0.0_dp, &
                                            runoff_sum        = 0.0_dp, &
                                            runoff_apl_sum    = 0.0_dp, &
                                            Q_b_tot_sum       = 0.0_dp, &
                                            Q_b_apl_sum       = 0.0_dp, &
                                            calving_sum       = 0.0_dp, &
                                            calving_apl_sum   = 0.0_dp, &
#if (DISC>0)   /* Ice discharge parameterisation */
                                            dis_perp_sum      = 0.0_dp, &
#endif
                                            q_geo_sum         = 0.0_dp, &
                                            Q_bm_sum          = 0.0_dp, &
                                            Q_tld_sum         = 0.0_dp, &
                                            am_perp_sum       = 0.0_dp, &
                                            dzs_dtau_sum      = 0.0_dp, &
                                            dzm_dtau_sum      = 0.0_dp, &
                                            dzb_dtau_sum      = 0.0_dp, &
                                            dzl_dtau_sum      = 0.0_dp, &
                                            dH_c_dtau_sum     = 0.0_dp, &
                                            dH_t_dtau_sum     = 0.0_dp, &
                                            dH_dtau_sum       = 0.0_dp, &
                                            q_gl_g_sum        = 0.0_dp

#if (defined(CLIMATOLOGY_EXTRACTION_HACK))
    !%% Climatology extraction hack (must not be used routinely)!
        real(dp), dimension(0:JMAX,0:IMAX), save :: temp_maat_sum = 0.0_dp, &
                                                    zs_sum        = 0.0_dp
#endif

#endif

real(dp), dimension(0:JMAX,0:IMAX) :: accum_flx         , &
                                      snowfall_flx      , &
                                      rainfall_flx      , &
                                      as_perp_flx       , &
                                      as_perp_apl_flx   , &
                                      smb_corr_flx      , &
                                      mb_source_apl_flx , &
                                      runoff_flx        , &
                                      runoff_apl_flx    , &
                                      Q_b_tot_flx       , &
                                      Q_b_apl_flx       , &
                                      calving_flx       , &
                                      calving_apl_flx   , &
#if (DISC>0)   /* Ice discharge parameterisation */
                                      dis_perp_flx      , &
#endif
                                      q_geo_flx         , &
                                      Q_bm_flx          , &
                                      Q_tld_flx         , &
                                      am_perp_flx       , &
                                      dzs_dtau_flx      , &
                                      dzm_dtau_flx      , &
                                      dzb_dtau_flx      , &
                                      dzl_dtau_flx      , &
                                      dH_c_dtau_flx     , &
                                      dH_t_dtau_flx     , &
                                      dH_dtau_flx       , &
                                      q_gl_g_flx

#if (defined(CLIMATOLOGY_EXTRACTION_HACK))
    !%% Climatology extraction hack (must not be used routinely)!
        real(dp), dimension(0:JMAX,0:IMAX) :: temp_maat_flx , &
                                              zs_flx
#endif

integer(i4b), dimension(0:IMAX,0:JMAX) :: mask_conv, mask_old_conv, &
                                          mask_ablation_type_conv, &
                                          n_cts_conv
integer(i4b), dimension(0:IMAX,0:JMAX) :: mask_region_conv
integer(i4b), dimension(0:IMAX,0:JMAX) :: kc_cts_conv
integer(i4b), dimension(0:IMAX,0:JMAX) :: mask_mar_conv
integer(i4b), dimension(0:IMAX,0:JMAX) :: flag_shelfy_stream_x_conv, &
                                          flag_shelfy_stream_y_conv, &
                                          flag_shelfy_stream_conv
integer(i4b), dimension(0:IMAX,0:JMAX) :: flag_grounding_line_1_conv, &
                                          flag_grounding_line_2_conv, &
                                          flag_calving_front_1_conv, &
                                          flag_calving_front_2_conv, &
                                          flag_grounded_front_a_1_conv, &
                                          flag_grounded_front_a_2_conv, &
                                          flag_grounded_front_b_1_conv, &
                                          flag_grounded_front_b_2_conv

real(dp) :: year2sec_conv, time_conv, &
            delta_ts_conv, glac_index_conv, z_sl_mean_conv, &
            V_tot_conv, V_af_conv, A_grounded_conv, A_floating_conv, &
            xi_conv(0:IMAX), eta_conv(0:JMAX), &
            sigma_level_c_conv(0:KCMAX), sigma_level_t_conv(0:KTMAX), &
            sigma_level_r_conv(0:KRMAX)

real(sp) :: H_R_conv

real(sp), dimension(0:IMAX,0:JMAX) :: lambda_conv, phi_conv, &
            lond_conv, latd_conv, &
            cell_area_conv, &
            temp_maat_conv, temp_s_conv, accum_conv, &
            snowfall_conv, rainfall_conv, pdd_conv, &
            as_perp_conv, as_perp_apl_conv, smb_corr_conv, &
            mb_source_apl_conv, runoff_conv, runoff_apl_conv, &
            z_sl_conv, &
            Q_b_tot_conv, Q_b_apl_conv, &
            calving_conv, calving_apl_conv, &
            q_geo_conv, &
            zs_conv, zm_conv, zb_conv, zl_conv, zl0_conv, &
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
            c_slide_conv, &
            temp_b_conv, temph_b_conv, &
            tau_dr_conv, tau_b_conv, &
            p_b_w_conv, q_w_conv, q_w_x_conv, q_w_y_conv, H_w_conv, &
            q_gl_g_conv, &
            cst_dist_conv, cos_grad_tc_conv, dis_perp_conv, &
            ratio_sl_x_conv, ratio_sl_y_conv, ratio_sl_conv, &
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

integer(i4b) :: cmode
integer(i4b) :: n_deflate_level
logical      :: flag_shuffle
integer(i4b) :: ncid, ncv
!     ncid:      ID of the output file
!     ncv:       Variable ID
integer(i4b) :: ncd, nc1d, nc2d(2), nc3d(3)
!     ncd:       Dimension ID
!     nc1d:      Dimension of a 1-d array
!     nc2d:      Vector with the dimensions of a 2-d array
!     nc3d:      Vector with the dimensions of a 3-d array
integer(i4b) :: nc2flag(2), nc3flag(3), nc4flag(4), nc5flag(5)
!     nc2flag:   Vector with the 2 possible values of a flag variable
!     nc3flag:   Vector with the 3 possible values of a flag variable
!     nc4flag:   Vector with the 4 possible values of a flag variable
!     nc5flag:   Vector with the 5 possible values of a flag variable
integer(i4b) :: nc1cor_i(1), nc1cor_j(1), &
                nc1cor_kc(1), nc1cor_kt(1), nc1cor_kr(1), &
                nc2cor_ij(2), &
                nc3cor_ijkc(3), nc3cor_ijkt(3), nc3cor_ijkr(3)
!     nc1cor(1): Corner of a 1-d array
!     nc2cor(2): Corner of a 2-d array
!     nc3cor(3): Corner of a 3-d array
integer(i4b) :: nc1cnt_i(1), nc1cnt_j(1), &
                nc1cnt_kc(1), nc1cnt_kt(1), nc1cnt_kr(1), &
                nc2cnt_ij(2), &
                nc3cnt_ijkc(3), nc3cnt_ijkt(3), nc3cnt_ijkr(3)
!     nc1cnt(1): Count of a 1-d array
!     nc2cnt(2): Count of a 2-d array
!     nc3cnt(3): Count of a 3-d array
character(len=256) :: buffer
character(len= 16), parameter :: filename_extension = '.nc'
character(len= 16), allocatable :: coord_id(:)

if (present(opt_flag_compute_flux_vars_only)) then
   flag_compute_flux_vars_only = opt_flag_compute_flux_vars_only
else
   flag_compute_flux_vars_only = .false.
end if

nc1cor_i = [ 1 ]
nc1cnt_i = [ IMAX+1 ]

nc1cor_j = [ 1 ]
nc1cnt_j = [ JMAX+1 ]

nc1cor_kc = [ 1 ]
nc1cnt_kc = [ KCMAX+1 ]

nc1cor_kt = [ 1 ]
nc1cnt_kt = [ KTMAX+1 ]

nc1cor_kr = [ 1 ]
nc1cnt_kr = [ KRMAX+1 ]

nc2cor_ij = [ 1, 1 ]
nc2cnt_ij = [ IMAX+1, JMAX+1 ]

nc3cor_ijkc = [ 1, 1, 1 ]
nc3cnt_ijkc = [ IMAX+1, JMAX+1, KCMAX+1 ]

nc3cor_ijkt = [ 1, 1, 1 ]
nc3cnt_ijkt = [ IMAX+1, JMAX+1, KTMAX+1 ]

nc3cor_ijkr = [ 1, 1, 1 ]
nc3cnt_ijkr = [ IMAX+1, JMAX+1, KRMAX+1 ]

!-------- Create consecutively numbered file names --------

if (.not.flag_compute_flux_vars_only) then

if (flag_3d_output) then
   ndat = ndat3d
else
   ndat = ndat2d
end if

if (ndat > 9999) then 
   errormsg = ' >>> output1: Too many time-slice files!'
   call error(errormsg)
endif

write(ch_ndat, '(i0.4)') ndat

if (flag_3d_output) then
   filename = trim(run_name)//ch_ndat//trim(filename_extension)
else
   filename = trim(run_name)//'_2d_'//ch_ndat//trim(filename_extension)
end if

filename_with_path = trim(OUT_PATH)//'/'//trim(filename)

end if   ! (.not.flag_compute_flux_vars_only)

!-------- File initialization --------

if (.not.flag_compute_flux_vars_only) then

if (allocated(coord_id)) deallocate(coord_id); allocate(coord_id(5))
coord_id(1) = 'x'; coord_id(2) = 'y'
coord_id(3) = 'zeta_c'; coord_id(4) = 'zeta_t'; coord_id(5) = 'zeta_r'

!  ------ Open NetCDF file

call set_cmode(cmode, n_deflate_level, flag_shuffle)

ios = nf90_create(trim(filename_with_path), cmode, ncid)

if (ios /= nf90_noerr) then
   errormsg = ' >>> output1: Error when opening a' &
            //               end_of_line &
            //'              NetCDF time-slice file!'
   call error(errormsg)
end if

!  ------ Global attributes

buffer = 'Time-slice output no. '//ch_ndat//' of simulation ' &
                                 //trim(run_name)
call check( nf90_put_att(ncid, NF90_GLOBAL, 'title', trim(buffer)), &
            thisroutine )

buffer = 'SICOPOLIS Version '//VERSION
call check( nf90_put_att(ncid, NF90_GLOBAL, 'source', trim(buffer)), &
            thisroutine )

call date_and_time(ch_date, ch_time, ch_zone)
buffer = ch_date(1:4)//'-'//ch_date(5:6)//'-'//ch_date(7:8)//' '// &
         ch_time(1:2)//':'//ch_time(3:4)//':'//ch_time(5:6)//' '// &
         ch_zone(1:3)//':'//ch_zone(4:5)//' - Data produced'
call check( nf90_put_att(ncid, NF90_GLOBAL, 'history', trim(buffer)), &
            thisroutine )

call set_ch_website(buffer)
call check( nf90_put_att(ncid, NF90_GLOBAL, 'references', trim(buffer)), &
            thisroutine )

call set_ch_institution(buffer)
call check( nf90_put_att(ncid, NF90_GLOBAL, 'institution', trim(buffer)), &
            thisroutine )

!  ------ Definition of the dimensions

call check( nf90_def_dim(ncid, trim(coord_id(1)), IMAX+1,  ncd), thisroutine )
call check( nf90_def_dim(ncid, trim(coord_id(2)), JMAX+1,  ncd), thisroutine )
call check( nf90_def_dim(ncid, trim(coord_id(3)), KCMAX+1, ncd), thisroutine )
call check( nf90_def_dim(ncid, trim(coord_id(4)), KTMAX+1, ncd), thisroutine )
call check( nf90_def_dim(ncid, trim(coord_id(5)), KRMAX+1, ncd), thisroutine )

!  ------ Definition of the variables

!    ---- mapping

call check( nf90_def_var(ncid, 'mapping', NF90_BYTE, ncv), thisroutine ) 

#if (GRID==0 || GRID==1)
buffer = 'polar_stereographic'
call check( nf90_put_att(ncid, ncv, 'grid_mapping_name', trim(buffer)), &
            thisroutine )
#elif (GRID==2)
buffer = 'latitude_longitude'
call check( nf90_put_att(ncid, ncv, 'grid_mapping_name', trim(buffer)), &
            thisroutine )
#endif

if (F_INV > 1.0e+10_dp) then   ! interpreted as infinity, thus no flattening

   call check( nf90_put_att(ncid, ncv, 'radius_of_sphere', R), &
               thisroutine )

else   ! finite inverse flattening

   call check( nf90_put_att(ncid, ncv, 'semi_major_axis', A), &
               thisroutine )

   call check( nf90_put_att(ncid, ncv, 'inverse_flattening', F_INV), &
               thisroutine )

end if

#if (GRID==0 || GRID==1)

lond0 = LAMBDA0*rad2deg
lond0 = modulo(lond0+180.0_dp, 360.0_dp)-180.0_dp
lond0 = nint(lond0*1.0e+04_dp)*1.0e-04_dp

latd0 = PHI0*rad2deg
if (latd0 >  90.0_dp) latd0 =  90.0_dp
if (latd0 < -90.0_dp) latd0 = -90.0_dp
latd0 = nint(latd0*1.0e+04_dp)*1.0e-04_dp
       ! reference longitude and standard parallel in deg rounded to 4 digits

if (latd0 >= 0.0_dp) then
   call check( nf90_put_att(ncid, ncv, &
                            'latitude_of_projection_origin',  90.0_dp), &
                            thisroutine )
else
   call check( nf90_put_att(ncid, ncv, &
                            'latitude_of_projection_origin', -90.0_dp), &
                            thisroutine )
end if

call check( nf90_put_att(ncid, ncv, &
                         'standard_parallel', latd0), &
                         thisroutine )

call check( nf90_put_att(ncid, ncv, &
                         'straight_vertical_longitude_from_pole', lond0), &
                         thisroutine )

call check( nf90_put_att(ncid, ncv, 'false_easting', 0.0_dp), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'false_northing', 0.0_dp), &
            thisroutine )

#endif

!    ---- year2sec

call check( nf90_def_var(ncid, 'year2sec', NF90_DOUBLE, ncv), &
            thisroutine )

buffer = 's a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'seconds_per_year'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = '1 year (1 a) in seconds'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )

!    ---- time

call check( nf90_def_var(ncid, 'time', NF90_DOUBLE, ncv), &
            thisroutine )

buffer = 'a'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'time'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Time'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )

if ((forcing_flag == 1).or.(forcing_flag == 3)) then

!    ---- delta_ts

   call check( nf90_def_var(ncid, 'delta_ts', NF90_DOUBLE, ncv), &
               thisroutine )

   buffer = 'degC'
   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
               thisroutine )
   buffer = 'surface_temperature_anomaly'
   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
               thisroutine )
   buffer = 'Surface temperature anomaly'
   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
               thisroutine )

else if (forcing_flag == 2) then

!    ---- glac_index

   call check( nf90_def_var(ncid, 'glac_index', NF90_DOUBLE, ncv), &
               thisroutine )

   buffer = '1'
   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
               thisroutine )
   buffer = 'glacial_index'
   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
               thisroutine )
   buffer = 'Glacial index'
   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
               thisroutine )

end if

!    ---- z_sl_mean

call check( nf90_def_var(ncid, 'z_sl_mean', NF90_DOUBLE, ncv), &
            thisroutine )

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'global_average_sea_level_change'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Mean sea level'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )

!    ---- V_tot

call check( nf90_def_var(ncid, 'V_tot', NF90_DOUBLE, ncv), &
            thisroutine )

buffer = 'm3'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_volume'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Ice volume'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )

!    ---- V_af

call check( nf90_def_var(ncid, 'V_af', NF90_DOUBLE, ncv), &
            thisroutine )

buffer = 'm3'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_volume_not_displacing_sea_water'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Ice volume above flotation'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )

!    ---- A_grounded

call check( nf90_def_var(ncid, 'A_grounded', NF90_DOUBLE, ncv), &
            thisroutine )

buffer = 'm2'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'grounded_land_ice_area'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Area covered by grounded ice'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )

!    ---- A_floating

call check( nf90_def_var(ncid, 'A_floating', NF90_DOUBLE, ncv), &
            thisroutine )

buffer = 'm2'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'floating_ice_shelf_area'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Area covered by floating ice'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )

!    ---- x (= xi)

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc1d), &
            thisroutine )

call check( nf90_def_var(ncid, 'x', NF90_DOUBLE, nc1d, ncv), &
            thisroutine )

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'projection_x_coordinate'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'x-coordinate of the grid point i'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'axis', 'x'), &
            thisroutine )

!    ---- y (= eta)

call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc1d), &
            thisroutine )

call check( nf90_def_var(ncid, 'y', NF90_DOUBLE, nc1d, ncv), &
            thisroutine )

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'projection_y_coordinate'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'y-coordinate of the grid point j'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'axis', 'y'), &
            thisroutine )

!    ---- sigma_level_c

call check( nf90_inq_dimid(ncid, trim(coord_id(3)), nc1d), &
            thisroutine )

call check( nf90_def_var(ncid, 'sigma_level_c', NF90_DOUBLE, nc1d, ncv), &
            thisroutine )

buffer = 'up'
call check( nf90_put_att(ncid, ncv, 'positive', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_kc_layer_sigma_coordinate'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'sigma-coordinate of the grid point kc'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )

!    ---- sigma_level_t

call check( nf90_inq_dimid(ncid, trim(coord_id(4)), nc1d), &
            thisroutine )

call check( nf90_def_var(ncid, 'sigma_level_t', NF90_DOUBLE, nc1d, ncv), &
            thisroutine )

buffer = 'up'
call check( nf90_put_att(ncid, ncv, 'positive', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_kt_layer_sigma_coordinate'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'sigma-coordinate of the grid point kt'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )

!    ---- sigma_level_r

call check( nf90_inq_dimid(ncid, trim(coord_id(5)), nc1d), &
            thisroutine )

call check( nf90_def_var(ncid, 'sigma_level_r', NF90_DOUBLE, nc1d, ncv), &
            thisroutine )

buffer = 'up'
call check( nf90_put_att(ncid, ncv, 'positive', trim(buffer)), &
            thisroutine )
buffer = 'lithosphere_layer_sigma_coordinate'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'sigma-coordinate of the grid point kr'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )

!    ---- lon

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'lon', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'lon', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'degrees_E'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'longitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Geographical longitude'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )

!    ---- lat

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'lat', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'lat', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'degrees_N'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'latitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Geographical latitude'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )

!    ---- lambda

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'lambda', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'lambda', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'rad'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'longitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Geographical longitude'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- phi

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'phi', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'phi', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'rad'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'latitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Geographical latitude'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- cell_area

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'cell_area', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'cell_area', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm2'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'cell_area'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Area of grid cell'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- temp_maat

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'temp_maat', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'temp_maat', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'degC'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'mean_annual_air_temperature'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Mean annual air temperature'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- temp_s

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'temp_s', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'temp_s', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'degC'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'surface_temperature'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Temperature at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- accum

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'prec', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'prec', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_precipitation'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Annual precipitation at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- snowfall

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'snowfall', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'snowfall', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_snowfall'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Annual snowfall at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- rainfall

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'rainfall', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'rainfall', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_rainfall'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Annual rainfall at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- pdd

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'pdd', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'pdd', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'degC a'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_positive_degree_days'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Positive degree days at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- as_perp

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'as_perp', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'as_perp', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_surface_mass_balance'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Mass balance at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- as_perp_apl

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'as_perp_apl', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'as_perp_apl', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'applied_land_ice_surface_mass_balance'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Applied mass balance at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- smb_corr

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'smb_corr', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'smb_corr', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_surface_mass_balance_diagnosed_correction'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Diagnosed correction of the mass balance at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- mb_source_apl

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'mb_source_apl', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'mb_source_apl', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'applied_land_ice_mass_balance'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Applied mass balance'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- runoff

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'runoff', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'runoff', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_surface_runoff'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Runoff rate at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- runoff_apl

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'runoff_apl', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'runoff_apl', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'applied_land_ice_surface_runoff'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Applied runoff rate at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- z_sl

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'z_sl', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'z_sl', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'sea_level_change'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Sea level'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- Q_b_tot

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'Q_b_tot', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'Q_b_tot', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'total_basal_melt'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Total basal melt'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- Q_b_apl

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'Q_b_apl', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'Q_b_apl', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'applied_total_basal_melt'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Applied total basal melt'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- calving

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'calving', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'calving', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_volume_flux_due_to_calving'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Calving flux'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- calving_apl

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'calving_apl', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'calving_apl', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'applied_land_ice_volume_flux_due_to_calving'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Applied calving flux'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

#if (DISC>0)   /* Ice discharge parameterisation */

!    ---- dis_perp

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'dis_perp', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'dis_perp', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'ice_discharge'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Ice discharge'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- cst_dist

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'cst_dist', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'cst_dist', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'km'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'coastal_distance'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Coastal distance'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- cos_grad_tc

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'cos_grad_tc', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'cos_grad_tc', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = '1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'cos_alpha'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Cosine of angle between surface gradient and cst dist gradient'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- mask_mar

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'mask_mar', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'mask_mar', NF90_BYTE, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'marginal_ring_mask'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Marginal ring mask'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
nc2flag = [ 0, 1 ]
call check( nf90_put_att(ncid, ncv, 'flag_values', nc2flag), &
            thisroutine )
buffer = 'no_ring '// &
         'ring'
call check( nf90_put_att(ncid, ncv, 'flag_meanings', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

#endif

!    ---- q_geo

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'q_geo', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'q_geo', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'W m-2'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'upward_geothermal_heat_flux_at_ground_level'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Geothermal heat flux'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- mask

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'mask', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'mask', NF90_BYTE, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'ice_land_sea_mask'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Ice-land-sea mask'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
nc4flag = [ 0, 1, 2, 3 ]
call check( nf90_put_att(ncid, ncv, 'flag_values', nc4flag), &
            thisroutine )
buffer = 'glaciated_land '// &
         'ice_free_land '// &
         'sea '// &
         'floating_ice'
call check( nf90_put_att(ncid, ncv, 'flag_meanings', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- mask_old

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'mask_old', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'mask_old', NF90_BYTE, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'ice_land_sea_mask_old'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Ice-land-sea mask (old)'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
nc4flag = [ 0, 1, 2, 3 ]
call check( nf90_put_att(ncid, ncv, 'flag_values', nc4flag), &
            thisroutine )
buffer = 'glaciated_land '// &
         'ice_free_land '// &
         'sea '// &
         'floating_ice'
call check( nf90_put_att(ncid, ncv, 'flag_meanings', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- mask_ablation_type

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'mask_ablation_type', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'mask_ablation_type', NF90_BYTE, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'mask_indicating_ablation_type'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Mask indicating ablation type'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
nc5flag = [ -2, -1, 1, 3, 9 ]
call check( nf90_put_att(ncid, ncv, 'flag_values', nc5flag), &
            thisroutine )
buffer = 'hidden_(ocean) '// &
         'hidden_(land) '// &
         'visible_(grounded_ice) '// &
         'visible_(floating_ice) '// &
         'visible_(misaccounted)'
call check( nf90_put_att(ncid, ncv, 'flag_meanings', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- mask_region

if (maxval(mask_region) > 0) then

   call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
               thisroutine )

#if (NETCDF4_ENABLED==1)
   call check( nf90_def_var(ncid, 'mask_region', NF90_INT, nc2d, ncv, &
               deflate_level=n_deflate_level, shuffle=flag_shuffle), &
               thisroutine )
#else
   call check( nf90_def_var(ncid, 'mask_region', NF90_INT, nc2d, ncv), &
               thisroutine )
#endif

   buffer = 'mask_region'
   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
               thisroutine )
   buffer = 'Region mask'
   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
               thisroutine )
   buffer = '0, 1, 2, ...'
   call check( nf90_put_att(ncid, ncv, 'flag_values', trim(buffer)), &
               thisroutine )
   buffer = 'undefined '// &
            'region_1 '// &
            'region_2 '// &
            '...'
   call check( nf90_put_att(ncid, ncv, 'flag_meanings', trim(buffer)), &
               thisroutine )
   call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
               thisroutine )

end if

!    ---- n_cts

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'n_cts', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'n_cts', NF90_BYTE, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'polythermal_condition_mask'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Mask for polythermal conditions'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
nc3flag = [ -1, 0, 1 ]
call check( nf90_put_att(ncid, ncv, 'flag_values', nc3flag), &
            thisroutine )
buffer = 'cold_base '// &
         'temperate_base_with_cold_ice_above '// &
         'temperate_base_with_temperate_ice_above'
call check( nf90_put_att(ncid, ncv, 'flag_meanings', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- kc_cts

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'kc_cts', NF90_INT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'kc_cts', NF90_INT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'CTS_position_grid_index'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Grid index of the CTS position'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- zs

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'zs', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'zs', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'surface_altitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Topography of the free surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- zm

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'zm', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'zm', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'zm_interface_altitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Topography of the z=zm interface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- zb

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'zb', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'zb', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'ice_base_altitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Topography of the ice base'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- zl

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'zl', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'zl', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'bedrock_altitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Topography of the lithosphere surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- zl0

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'zl0', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'zl0', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'isostatically_relaxed_bedrock_altitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Topography of the isostatically relaxed lithosphere surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- H_cold

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'H_cold', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'H_cold', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_cold_layer_thickness'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Thickness of the cold ice layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- H_temp

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'H_temp', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'H_temp', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_temperate_layer_thickness'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Thickness of the temperate ice layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- H

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'H', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'H', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_thickness'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Ice thickness'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- H_R

call check( nf90_def_var(ncid, 'H_R', NF90_FLOAT, ncv), &
            thisroutine )

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'lithosphere_layer_thickness'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Thickness of the lithosphere layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )

!    ---- Q_bm

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'Q_bm', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'Q_bm', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_basal_melt_rate'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Basal melting rate'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- Q_tld

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'Q_tld', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'Q_tld', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_temperate_layer_water_drainage'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Water drainage from the temperate layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- am_perp

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'am_perp', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'am_perp', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_volume_flux_across_zm_interface'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Volume flux across the z=zm interface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- qx

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'qx', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'qx', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm2 a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_vertical_integral_x_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Horizontal volume flux qx'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
buffer = 'Staggered grid variable, defined at (j,i+1/2)'
call check( nf90_put_att(ncid, ncv, 'comment', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- qy

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'qy', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'qy', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm2 a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_vertical_integral_y_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Horizontal volume flux qy'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
buffer = 'Staggered grid variable, defined at (j+1/2,i)'
call check( nf90_put_att(ncid, ncv, 'comment', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- vx_m_sia

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vx_m_sia', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'vx_m_sia', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_vertical_mean_x_velocity_sia'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Vertical mean of horizontal velocity vx (SIA)'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
buffer = 'Staggered grid variable, defined at (j,i+1/2)'
call check( nf90_put_att(ncid, ncv, 'comment', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- vy_m_sia

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vy_m_sia', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'vy_m_sia', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_vertical_mean_y_velocity_sia'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Vertical mean of horizontal velocity vy (SIA)'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
buffer = 'Staggered grid variable, defined at (j+1/2,i)'
call check( nf90_put_att(ncid, ncv, 'comment', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- vx_m_ssa

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vx_m_ssa', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'vx_m_ssa', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_vertical_mean_x_velocity_ssa'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Vertical mean of horizontal velocity vx (SSA)'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
buffer = 'Staggered grid variable, defined at (j,i+1/2)'
call check( nf90_put_att(ncid, ncv, 'comment', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- vy_m_ssa

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vy_m_ssa', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'vy_m_ssa', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_vertical_mean_y_velocity_ssa'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Vertical mean of horizontal velocity vy (SSA)'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
buffer = 'Staggered grid variable, defined at (j+1/2,i)'
call check( nf90_put_att(ncid, ncv, 'comment', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- dzs_dt (= dzs_dtau)

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'dzs_dt', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'dzs_dt', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'tendency_of_surface_altitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Rate of change of the topography of the free surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- dzm_dt (= dzm_dtau)

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'dzm_dt', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'dzm_dt', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'tendency_of_zm_interface_altitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Rate of change of the topography of the z=zm interface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- dzb_dt (= dzb_dtau)

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'dzb_dt', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'dzb_dt', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'tendency_of_ice_base_altitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Rate of change of the topography of the ice base'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- dzl_dt (= dzl_dtau)

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'dzl_dt', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'dzl_dt', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'tendency_of_bedrock_altitude'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Rate of change of the topography of the lithosphere surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- dH_c_dt (= dH_c_dtau)

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'dH_c_dt', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'dH_c_dt', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'tendency_of_land_ice_kc_layer_thickness'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Rate of change of the thickness of the upper (kc) ice layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- dH_t_dt (= dH_t_dtau)

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'dH_t_dt', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'dH_t_dt', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'tendency_of_land_ice_kt_layer_thickness'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Rate of change of the thickness of the lower (kt) ice layer'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- dH_dt (= dH_dtau)

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'dH_dt', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'dH_dt', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'tendency_of_land_ice_thickness'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Rate of change of the ice thickness'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- vx_b_g

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vx_b_g', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'vx_b_g', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_base_x_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Horizontal velocity vx at the ice base'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- vy_b_g

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vy_b_g', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'vy_b_g', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_base_y_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Horizontal velocity vy at the ice base'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- vz_b

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vz_b', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'vz_b', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_base_z_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Vertical velocity vz at the ice base'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- vh_b

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vh_b', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'vh_b', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_base_horizontal_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Horizontal velocity vh at the ice base'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- vx_s_g

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vx_s_g', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'vx_s_g', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_surface_x_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Horizontal velocity vx at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- vy_s_g

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vy_s_g', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'vy_s_g', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_surface_y_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Horizontal velocity vy at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- vz_s

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vz_s', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'vz_s', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_surface_z_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Vertical velocity vz at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- vh_s

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vh_s', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'vh_s', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_surface_horizontal_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Horizontal velocity vh at the ice surface'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- vx_m_g

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vx_m_g', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'vx_m_g', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_vertical_mean_x_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Vertical mean of horizontal velocity vx'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- vy_m_g

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vy_m_g', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'vy_m_g', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_vertical_mean_y_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Vertical mean of horizontal velocity vy'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- vh_m

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vh_m', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'vh_m', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_vertical_mean_horizontal_velocity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Vertical mean of horizontal velocity vh'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- c_slide

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'c_slide', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'c_slide', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm a-1 Pa-(p-q)'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_base_sliding_coefficient'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Basal sliding coefficient'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- temp_b

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'temp_b', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'temp_b', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'degC'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'basal_temperature'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Temperature at the ice base'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- temph_b

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'temph_b', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'temph_b', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'degC'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'basal_temperature_rel_to_pmp'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Temperature at the ice base relative to the pressure melting point'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- tau_dr

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'tau_dr', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'tau_dr', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'Pa'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'magnitude_of_land_ice_driving_stress'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Driving stress'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- tau_b

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'tau_b', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'tau_b', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'Pa'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'magnitude_of_land_ice_basal_drag'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Basal drag'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- p_b_w

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'p_b_w', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'p_b_w', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'Pa'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'basal_water_pressure'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Basal water pressure'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- q_w

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'q_w', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'q_w', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm2 a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'basal_water_flux'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Volume flux of basal water'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- q_w_x

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'q_w_x', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'q_w_x', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm2 a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'basal_water_flux_x'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Volume flux of basal water in x-direction'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- q_w_y

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'q_w_y', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'q_w_y', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm2 a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'basal_water_flux_y'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Volume flux of basal water in y-direction'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- H_w

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'H_w', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'H_w', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'water_column_thickness'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Thickness of the water column under the ice base'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- q_gl_g

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'q_gl_g', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'q_gl_g', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'm2 a-1'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_volume_flux_across_gl'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Horizontal volume flux across the grounding line'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- ratio_sl_x

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'ratio_sl_x', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'ratio_sl_x', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_x_slip_ratio'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Ratio of basal to surface velocity (slip ratio) in x-direction, ' &
         // 'at (i+1/2,j)'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- ratio_sl_y

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'ratio_sl_y', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'ratio_sl_y', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_y_slip_ratio'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Ratio of basal to surface velocity (slip ratio) in y-direction, ' &
         // 'at (i+1/2,j)'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- ratio_sl

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'ratio_sl', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'ratio_sl', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_slip_ratio'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Ratio of basal to surface velocity (slip ratio)'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- flag_shelfy_stream_x

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'flag_shelfy_stream_x', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'flag_shelfy_stream_x', NF90_BYTE, nc2d, ncv), &
            thisroutine )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
        thisroutine )
buffer = 'land_ice_x_shelfy_stream_flag'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Shelfy stream flag in x-direction, at (i+1/2,j)'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
        thisroutine )

!    ---- flag_shelfy_stream_y

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'flag_shelfy_stream_y', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'flag_shelfy_stream_y', NF90_BYTE, nc2d, ncv), &
            thisroutine )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_y_shelfy_stream_flag'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Shelfy stream flag in y-direction, at (i,j+1/2)'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- flag_shelfy_stream

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'flag_shelfy_stream', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'flag_shelfy_stream', NF90_BYTE, nc2d, ncv), &
            thisroutine )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_shelfy_stream_flag'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Shelfy stream flag'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- flag_grounding_line_1

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'flag_grounding_line_1', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'flag_grounding_line_1', NF90_BYTE, nc2d, ncv), &
            thisroutine )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_grounding_line_1_flag'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Grounding line flag 1'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- flag_grounding_line_2

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'flag_grounding_line_2', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'flag_grounding_line_2', NF90_BYTE, nc2d, ncv), &
            thisroutine )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_grounding_line_2_flag'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Grounding line flag 2'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- flag_calving_front_1

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'flag_calving_front_1', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'flag_calving_front_1', NF90_BYTE, nc2d, ncv), &
            thisroutine )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_calving_front_1_flag'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Calving front flag 1'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- flag_calving_front_2

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'flag_calving_front_2', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'flag_calving_front_2', NF90_BYTE, nc2d, ncv), &
            thisroutine )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_calving_front_2_flag'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Calving front flag 2'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- flag_grounded_front_a_1

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'flag_grounded_front_a_1', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'flag_grounded_front_a_1', NF90_BYTE, nc2d, ncv), &
            thisroutine )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_grounded_front_a_1_flag'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Grounded front flag a 1'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- flag_grounded_front_a_2

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'flag_grounded_front_a_2', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'flag_grounded_front_a_2', NF90_BYTE, nc2d, ncv), &
            thisroutine )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_grounded_front_a_2_flag'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Grounded front flag a 2'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- flag_grounded_front_b_1

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'flag_grounded_front_b_1', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'flag_grounded_front_b_1', NF90_BYTE, nc2d, ncv), &
            thisroutine )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_grounded_front_b_1_flag'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Grounded front flag b 1'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- flag_grounded_front_b_2

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'flag_grounded_front_b_2', NF90_BYTE, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'flag_grounded_front_b_2', NF90_BYTE, nc2d, ncv), &
            thisroutine )
#endif

buffer = '-'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_grounded_front_b_2_flag'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Grounded front flag b 2'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- vis_ave_g

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vis_ave_g', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'vis_ave_g', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'Pa s'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_vertical_average_viscosity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Depth-averaged viscosity (SIA/SSA)'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

!    ---- vis_int_g

call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc2d(1)), &
            thisroutine )
call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc2d(2)), &
            thisroutine )

#if (NETCDF4_ENABLED==1)
call check( nf90_def_var(ncid, 'vis_int_g', NF90_FLOAT, nc2d, ncv, &
            deflate_level=n_deflate_level, shuffle=flag_shuffle), &
            thisroutine )
#else
call check( nf90_def_var(ncid, 'vis_int_g', NF90_FLOAT, nc2d, ncv), &
            thisroutine )
#endif

buffer = 'Pa s m'
call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
buffer = 'land_ice_vertical_integral_viscosity'
call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
            thisroutine )
buffer = 'Depth-integrated viscosity (SIA/SSA)'
call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
            thisroutine )
call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
            thisroutine )

if (flag_3d_output) then

!    ---- vx_c

   call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(3)), nc3d(3)), &
               thisroutine )

#if (NETCDF4_ENABLED==1)
   call check( nf90_def_var(ncid, 'vx_c', NF90_FLOAT, nc3d, ncv, &
               deflate_level=n_deflate_level, shuffle=flag_shuffle), &
               thisroutine )
#else
   call check( nf90_def_var(ncid, 'vx_c', NF90_FLOAT, nc3d, ncv), &
               thisroutine )
#endif

   buffer = 'm a-1'
   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
               thisroutine )
   buffer = 'land_ice_kc_layer_x_velocity'
   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
               thisroutine )
   buffer = 'Horizontal velocity vx in the upper (kc) ice layer'
   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
               thisroutine )
   buffer = 'Staggered grid variable, defined at (kc,j,i+1/2)'
   call check( nf90_put_att(ncid, ncv, 'comment', trim(buffer)), &
               thisroutine )
   call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
               thisroutine )

!    ---- vy_c

   call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(3)), nc3d(3)), &
               thisroutine )

#if (NETCDF4_ENABLED==1)
   call check( nf90_def_var(ncid, 'vy_c', NF90_FLOAT, nc3d, ncv, &
               deflate_level=n_deflate_level, shuffle=flag_shuffle), &
               thisroutine )
#else
   call check( nf90_def_var(ncid, 'vy_c', NF90_FLOAT, nc3d, ncv), &
               thisroutine )
#endif

   buffer = 'm a-1'
   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
               thisroutine )
   buffer = 'land_ice_kc_layer_y_velocity'
   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
               thisroutine )
   buffer = 'Horizontal velocity vy in the upper (kc) ice layer'
   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
               thisroutine )
   buffer = 'Staggered grid variable, defined at (kc,j+1/2,i)'
   call check( nf90_put_att(ncid, ncv, 'comment', trim(buffer)), &
               thisroutine )
   call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
               thisroutine )

!    ---- vz_c

   call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(3)), nc3d(3)), &
               thisroutine )

#if (NETCDF4_ENABLED==1)
   call check( nf90_def_var(ncid, 'vz_c', NF90_FLOAT, nc3d, ncv, &
               deflate_level=n_deflate_level, shuffle=flag_shuffle), &
               thisroutine )
#else
   call check( nf90_def_var(ncid, 'vz_c', NF90_FLOAT, nc3d, ncv), &
               thisroutine )
#endif

   buffer = 'm a-1'
   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
               thisroutine )
   buffer = 'land_ice_kc_layer_z_velocity'
   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
               thisroutine )
   buffer = 'Vertical velocity vz in the upper (kc) ice layer'
   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
               thisroutine )
   buffer = 'Staggered grid variable, defined at (kc+1/2,j,i)'
   call check( nf90_put_att(ncid, ncv, 'comment', trim(buffer)), &
               thisroutine )
   call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
               thisroutine )

!    ---- vx_t

   call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(4)), nc3d(3)), &
               thisroutine )

#if (NETCDF4_ENABLED==1)
   call check( nf90_def_var(ncid, 'vx_t', NF90_FLOAT, nc3d, ncv, &
               deflate_level=n_deflate_level, shuffle=flag_shuffle), &
               thisroutine )
#else
   call check( nf90_def_var(ncid, 'vx_t', NF90_FLOAT, nc3d, ncv), &
               thisroutine )
#endif

   buffer = 'm a-1'
   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
               thisroutine )
   buffer = 'land_ice_kt_layer_x_velocity'
   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
               thisroutine )
   buffer = 'Horizontal velocity vx in the lower (kt) ice layer'
   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
               thisroutine )
   buffer = 'Staggered grid variable, defined at (kt,j,i+1/2)'
   call check( nf90_put_att(ncid, ncv, 'comment', trim(buffer)), &
               thisroutine )
   call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
               thisroutine )

!    ---- vy_t

   call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(4)), nc3d(3)), &
               thisroutine )

#if (NETCDF4_ENABLED==1)
   call check( nf90_def_var(ncid, 'vy_t', NF90_FLOAT, nc3d, ncv, &
               deflate_level=n_deflate_level, shuffle=flag_shuffle), &
               thisroutine )
#else
   call check( nf90_def_var(ncid, 'vy_t', NF90_FLOAT, nc3d, ncv), &
               thisroutine )
#endif

   buffer = 'm a-1'
   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
            thisroutine )
   buffer = 'land_ice_kt_layer_y_velocity'
   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
               thisroutine )
   buffer = 'Horizontal velocity vy in the lower (kt) ice layer'
   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
               thisroutine )
   buffer = 'Staggered grid variable, defined at (kt,j+1/2,i)'
   call check( nf90_put_att(ncid, ncv, 'comment', trim(buffer)), &
               thisroutine )
   call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
               thisroutine )

!    ---- vz_t

   call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(4)), nc3d(3)), &
               thisroutine )

#if (NETCDF4_ENABLED==1)
   call check( nf90_def_var(ncid, 'vz_t', NF90_FLOAT, nc3d, ncv, &
               deflate_level=n_deflate_level, shuffle=flag_shuffle), &
               thisroutine )
#else
   call check( nf90_def_var(ncid, 'vz_t', NF90_FLOAT, nc3d, ncv), &
               thisroutine )
#endif

   buffer = 'm a-1'
   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
               thisroutine )
   buffer = 'land_ice_kt_layer_z_velocity'
   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
               thisroutine )
   buffer = 'Vertical velocity vz in the lower (kt) ice layer'
   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
               thisroutine )
   buffer = 'Staggered grid variable, defined at (kt+1/2,j,i)'
   call check( nf90_put_att(ncid, ncv, 'comment', trim(buffer)), &
               thisroutine )
   call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
               thisroutine )

!    ---- temp_c

   call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(3)), nc3d(3)), &
               thisroutine )

#if (NETCDF4_ENABLED==1)
   call check( nf90_def_var(ncid, 'temp_c', NF90_FLOAT, nc3d, ncv, &
               deflate_level=n_deflate_level, shuffle=flag_shuffle), &
               thisroutine )
#else
   call check( nf90_def_var(ncid, 'temp_c', NF90_FLOAT, nc3d, ncv), &
               thisroutine )
#endif

   buffer = 'degC'
   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
               thisroutine )
   buffer = 'land_ice_kc_layer_temperature'
   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
               thisroutine )
   buffer = 'Temperature in the upper (kc) ice layer'
   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
               thisroutine )
   call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
               thisroutine )

!    ---- omega_t

   call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(4)), nc3d(3)), &
               thisroutine )

#if (NETCDF4_ENABLED==1)
   call check( nf90_def_var(ncid, 'omega_t', NF90_FLOAT, nc3d, ncv, &
               deflate_level=n_deflate_level, shuffle=flag_shuffle), &
               thisroutine )
#else
   call check( nf90_def_var(ncid, 'omega_t', NF90_FLOAT, nc3d, ncv), &
               thisroutine )
#endif

   buffer = '1'
   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
               thisroutine )
   buffer = 'land_ice_kt_layer_water_content'
   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
               thisroutine )
   buffer = 'Water content in the lower (kt) ice layer'
   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
               thisroutine )
   call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
               thisroutine )

!    ---- temp_r

   call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(5)), nc3d(3)), &
               thisroutine )

#if (NETCDF4_ENABLED==1)
   call check( nf90_def_var(ncid, 'temp_r', NF90_FLOAT, nc3d, ncv, &
               deflate_level=n_deflate_level, shuffle=flag_shuffle), &
               thisroutine )
#else
   call check( nf90_def_var(ncid, 'temp_r', NF90_FLOAT, nc3d, ncv), &
               thisroutine )
#endif

   buffer = 'degC'
   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
               thisroutine )
   buffer = 'lithosphere_layer_temperature'
   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
               thisroutine )
   buffer = 'Temperature in the lithosphere layer'
   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
               thisroutine )
   call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
               thisroutine )

!    ---- enth_c

   call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(3)), nc3d(3)), &
               thisroutine )

#if (NETCDF4_ENABLED==1)
   call check( nf90_def_var(ncid, 'enth_c', NF90_FLOAT, nc3d, ncv, &
               deflate_level=n_deflate_level, shuffle=flag_shuffle), &
               thisroutine )
#else
   call check( nf90_def_var(ncid, 'enth_c', NF90_FLOAT, nc3d, ncv), &
               thisroutine )
#endif

   buffer = 'J kg-1'
   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
               thisroutine )
   buffer = 'land_ice_kc_layer_enthalpy'
   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
               thisroutine )
   buffer = 'Enthalpy in the upper (kc) ice layer'
   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
               thisroutine )
   call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
               thisroutine )

!    ---- enth_t

   call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(4)), nc3d(3)), &
               thisroutine )

#if (NETCDF4_ENABLED==1)
   call check( nf90_def_var(ncid, 'enth_t', NF90_FLOAT, nc3d, ncv, &
               deflate_level=n_deflate_level, shuffle=flag_shuffle), &
               thisroutine )
#else
   call check( nf90_def_var(ncid, 'enth_t', NF90_FLOAT, nc3d, ncv), &
               thisroutine )
#endif

   buffer = 'J kg-1'
   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
               thisroutine )
   buffer = 'land_ice_kt_layer_enthalpy'
   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
               thisroutine )
   buffer = 'Enthalpy in the lower (kt) ice layer'
   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
               thisroutine )
   call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
               thisroutine )

!    ---- omega_c

   call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(3)), nc3d(3)), &
               thisroutine )

#if (NETCDF4_ENABLED==1)
   call check( nf90_def_var(ncid, 'omega_c', NF90_FLOAT, nc3d, ncv, &
               deflate_level=n_deflate_level, shuffle=flag_shuffle), &
               thisroutine )
#else
   call check( nf90_def_var(ncid, 'omega_c', NF90_FLOAT, nc3d, ncv), &
               thisroutine )
#endif

   buffer = '1'
   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
               thisroutine )
   buffer = 'land_ice_kc_layer_water_content'
   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
               thisroutine )
   buffer = 'Water content in the upper (kc) ice layer'
   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
               thisroutine )
   call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
               thisroutine )

!    ---- enh_c

   call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(3)), nc3d(3)), &
               thisroutine )

#if (NETCDF4_ENABLED==1)
   call check( nf90_def_var(ncid, 'enh_c', NF90_FLOAT, nc3d, ncv, &
               deflate_level=n_deflate_level, shuffle=flag_shuffle), &
               thisroutine )
#else
   call check( nf90_def_var(ncid, 'enh_c', NF90_FLOAT, nc3d, ncv), &
               thisroutine )
#endif

   buffer = '1'
   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
               thisroutine )
   buffer = 'land_ice_kc_layer_flow_enhancement_factor'
   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
               thisroutine )
   buffer = 'Flow enhancement factor in the upper (kc) ice layer'
   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
               thisroutine )
   call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
               thisroutine )

!    ---- enh_t

   call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(4)), nc3d(3)), &
               thisroutine )

#if (NETCDF4_ENABLED==1)
   call check( nf90_def_var(ncid, 'enh_t', NF90_FLOAT, nc3d, ncv, &
               deflate_level=n_deflate_level, shuffle=flag_shuffle), &
               thisroutine )
#else
   call check( nf90_def_var(ncid, 'enh_t', NF90_FLOAT, nc3d, ncv), &
               thisroutine )
#endif

   buffer = '1'
   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
               thisroutine )
   buffer = 'land_ice_kt_layer_flow_enhancement_factor'
   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
               thisroutine )
   buffer = 'Flow enhancement factor in the lower (kt) ice layer'
   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
               thisroutine )
   call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
               thisroutine )

!    ---- strain_heating_c

   call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(3)), nc3d(3)), &
               thisroutine )

#if (NETCDF4_ENABLED==1)
   call check( nf90_def_var(ncid, 'strain_heating_c', NF90_FLOAT, nc3d, ncv, &
               deflate_level=n_deflate_level, shuffle=flag_shuffle), &
               thisroutine )
#else
   call check( nf90_def_var(ncid, 'strain_heating_c', NF90_FLOAT, nc3d, ncv), &
               thisroutine )
#endif

   buffer = 'W kg-1'
   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
               thisroutine )
   buffer = 'land_ice_kc_layer_strain_heating'
   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
               thisroutine )
   buffer = 'Strain heating in the upper (kc) ice layer'
   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
               thisroutine )
   call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
               thisroutine )

!    ---- strain_heating_t

   call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(4)), nc3d(3)), &
               thisroutine )

#if (NETCDF4_ENABLED==1)
   call check( nf90_def_var(ncid, 'strain_heating_t', NF90_FLOAT, nc3d, ncv, &
               deflate_level=n_deflate_level, shuffle=flag_shuffle), &
               thisroutine )
#else
   call check( nf90_def_var(ncid, 'strain_heating_t', NF90_FLOAT, nc3d, ncv), &
               thisroutine )
#endif

   buffer = 'W kg-1'
   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
               thisroutine )
   buffer = 'land_ice_kt_layer_strain_heating'
   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
               thisroutine )
   buffer = 'Strain heating in the lower (kt) ice layer'
   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
               thisroutine )
   call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
               thisroutine )

!    ---- age_c

   call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(3)), nc3d(3)), &
               thisroutine )

#if (NETCDF4_ENABLED==1)
   call check( nf90_def_var(ncid, 'age_c', NF90_FLOAT, nc3d, ncv, &
               deflate_level=n_deflate_level, shuffle=flag_shuffle), &
               thisroutine )
#else
   call check( nf90_def_var(ncid, 'age_c', NF90_FLOAT, nc3d, ncv), &
               thisroutine )
#endif

   buffer = 'a'
   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
               thisroutine )
   buffer = 'land_ice_kc_layer_age'
   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
               thisroutine )
   buffer = 'Age in the upper (kc) ice layer'
   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
               thisroutine )
   call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
               thisroutine )

!    ---- age_t

   call check( nf90_inq_dimid(ncid, trim(coord_id(1)), nc3d(1)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(2)), nc3d(2)), &
               thisroutine )
   call check( nf90_inq_dimid(ncid, trim(coord_id(4)), nc3d(3)), &
               thisroutine )

#if (NETCDF4_ENABLED==1)
   call check( nf90_def_var(ncid, 'age_t', NF90_FLOAT, nc3d, ncv, &
               deflate_level=n_deflate_level, shuffle=flag_shuffle), &
               thisroutine )
#else
   call check( nf90_def_var(ncid, 'age_t', NF90_FLOAT, nc3d, ncv), &
               thisroutine )
#endif

   buffer = 'a'
   call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
               thisroutine )
   buffer = 'land_ice_kt_layer_age'
   call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
               thisroutine )
   buffer = 'Age in the lower (kt) ice layer'
   call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
               thisroutine )
   call check( nf90_put_att(ncid, ncv, 'grid_mapping', 'mapping'), &
               thisroutine )
                               
end if

!    ---- End of the definitions

call check( nf90_enddef(ncid), thisroutine )

end if   ! (.not.flag_compute_flux_vars_only)

!-------- Thickness of the cold and temperate layers --------

if (.not.flag_compute_flux_vars_only) then

H_cold = 0.0_dp
H_temp = 0.0_dp

#if (CALCMOD==1)
do i=1, IMAX-1
do j=1, JMAX-1
   H_temp(j,i) = H_t(j,i)
end do
end do
#elif (CALCMOD==0 || CALCMOD==2 || CALCMOD==3 || CALCMOD==-1)
do i=1, IMAX-1
do j=1, JMAX-1
   H_temp(j,i) = H_c(j,i)*eaz_c_quotient(kc_cts(j,i))
end do
end do
#endif

H_cold = H - H_temp

end if   ! (.not.flag_compute_flux_vars_only)

!-------- Enthalpies for the non-enthalpy methods (POLY, COLD, ISOT) --------

if (.not.flag_compute_flux_vars_only) then

#if (CALCMOD==1)

do i=0, IMAX
do j=0, JMAX

   do kc=0, KCMAX
      enth_c(kc,j,i) = enth_fct_temp_omega(temp_c(kc,j,i), 0.0_dp)
   end do

   if ( (mask(j,i)==0).and.(n_cts(j,i)==1) ) then
      do kt=0, KTMAX
         enth_t(kt,j,i) = enth_fct_temp_omega(temp_t_m(kt,j,i), omega_t(kt,j,i))
      end do
   else
      do kt=0, KTMAX
         enth_t(kt,j,i) = enth_c(0,j,i)
      end do
   end if

end do
end do

#elif (CALCMOD==0 || CALCMOD==-1)

do i=0, IMAX
do j=0, JMAX

   do kc=0, KCMAX
      enth_c(kc,j,i) = enth_fct_temp_omega(temp_c(kc,j,i), 0.0_dp)
   end do

   do kt=0, KTMAX
      enth_t(kt,j,i) = enth_c(0,j,i)
   end do

end do
end do

#endif

end if   ! (.not.flag_compute_flux_vars_only)

!-------- Vertical mean of horizontal velocities --------

if (.not.flag_compute_flux_vars_only) then

vx_m_g = 0.0_dp
vy_m_g = 0.0_dp

do i=1, IMAX-1
do j=1, JMAX-1

   if ( (mask(j,i)==0).or.(mask(j,i)==3) ) then
      vx_m_g(j,i) = 0.5_dp*(vx_m(j,i)+vx_m(j,i-1))
      vy_m_g(j,i) = 0.5_dp*(vy_m(j,i)+vy_m(j-1,i))
   end if

end do
end do

end if   ! (.not.flag_compute_flux_vars_only)

!-------- Computation of scalar volumes and areas --------

if (.not.flag_compute_flux_vars_only) then

rhosw_rho_ratio = RHO_SW/RHO

V_grounded = 0.0_dp
V_floating = 0.0_dp
V_gr_redu  = 0.0_dp
A_grounded = 0.0_dp
A_floating = 0.0_dp

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i)==0) then   ! grounded ice

      V_grounded = V_grounded + H(j,i)*cell_area(j,i)
      A_grounded = A_grounded + cell_area(j,i)
      V_gr_redu  = V_gr_redu &
                   + rhosw_rho_ratio  &
                        *max((z_sl(j,i)-zl(j,i)),0.0_dp)*cell_area(j,i)

   else if (mask(j,i)==3) then   ! floating ice

      V_floating = V_floating + H(j,i)*cell_area(j,i)
      A_floating = A_floating + cell_area(j,i)

   end if

end do
end do

V_tot   = V_grounded + V_floating
V_af    = V_grounded - V_gr_redu

end if   ! (.not.flag_compute_flux_vars_only)

!-------- Flux variables --------

#if (!defined(OUTPUT_FLUX_VARS) || OUTPUT_FLUX_VARS==1)
       ! snapshots of flux variables

accum_flx         = accum
snowfall_flx      = snowfall
rainfall_flx      = rainfall
as_perp_flx       = as_perp
as_perp_apl_flx   = as_perp_apl
smb_corr_flx      = smb_corr
mb_source_apl_flx = mb_source_apl
runoff_flx        = runoff
runoff_apl_flx    = runoff_apl
Q_b_tot_flx       = Q_b_tot
Q_b_apl_flx       = Q_b_apl
calving_flx       = calving
calving_apl_flx   = calving_apl
#if (DISC>0)   /* Ice discharge parameterisation */
dis_perp_flx      = dis_perp
#endif
q_geo_flx         = q_geo
Q_bm_flx          = Q_bm
Q_tld_flx         = Q_tld
am_perp_flx       = am_perp
dzs_dtau_flx      = dzs_dtau
dzm_dtau_flx      = dzm_dtau
dzb_dtau_flx      = dzb_dtau
dzl_dtau_flx      = dzl_dtau
dH_c_dtau_flx     = dH_c_dtau
dH_t_dtau_flx     = dH_t_dtau
dH_dtau_flx       = dH_dtau
q_gl_g_flx        = q_gl_g

#if (defined(CLIMATOLOGY_EXTRACTION_HACK))
    !%% Climatology extraction hack (must not be used routinely)!
        temp_maat_flx = temp_maat
        zs_flx        = zs
#endif

#elif (OUTPUT_FLUX_VARS==2)

if ( .not.((OUTPUT==3).and.(flag_3d_output)) ) then
   ! averaging of flux variables

   if (n_flx_ave_cnt==0) then

      accum_sum         = 0.0_dp
      snowfall_sum      = 0.0_dp
      rainfall_sum      = 0.0_dp
      as_perp_sum       = 0.0_dp
      as_perp_apl_sum   = 0.0_dp
      smb_corr_sum      = 0.0_dp
      mb_source_apl_sum = 0.0_dp
      runoff_sum        = 0.0_dp
      runoff_apl_sum    = 0.0_dp
      Q_b_tot_sum       = 0.0_dp
      Q_b_apl_sum       = 0.0_dp
      calving_sum       = 0.0_dp
      calving_apl_sum   = 0.0_dp
#if (DISC>0)   /* Ice discharge parameterisation */
      dis_perp_sum      = 0.0_dp
#endif
      q_geo_sum         = 0.0_dp
      Q_bm_sum          = 0.0_dp
      Q_tld_sum         = 0.0_dp
      am_perp_sum       = 0.0_dp
      dzs_dtau_sum      = 0.0_dp
      dzm_dtau_sum      = 0.0_dp
      dzb_dtau_sum      = 0.0_dp
      dzl_dtau_sum      = 0.0_dp
      dH_c_dtau_sum     = 0.0_dp
      dH_t_dtau_sum     = 0.0_dp
      dH_dtau_sum       = 0.0_dp
      q_gl_g_sum        = 0.0_dp

#if (defined(CLIMATOLOGY_EXTRACTION_HACK))
    !%% Climatology extraction hack (must not be used routinely)!
        temp_maat_sum = 0.0_dp
        zs_sum        = 0.0_dp
#endif

   end if

   n_flx_ave_cnt = n_flx_ave_cnt + 1

   accum_sum         = accum_sum         + accum
   snowfall_sum      = snowfall_sum      + snowfall
   rainfall_sum      = rainfall_sum      + rainfall
   as_perp_sum       = as_perp_sum       + as_perp
   as_perp_apl_sum   = as_perp_apl_sum   + as_perp_apl
   smb_corr_sum      = smb_corr_sum      + smb_corr
   mb_source_apl_sum = mb_source_apl_sum + mb_source_apl
   runoff_sum        = runoff_sum        + runoff
   runoff_apl_sum    = runoff_apl_sum    + runoff_apl
   Q_b_tot_sum       = Q_b_tot_sum       + Q_b_tot
   Q_b_apl_sum       = Q_b_apl_sum       + Q_b_apl
   calving_sum       = calving_sum       + calving
   calving_apl_sum   = calving_apl_sum   + calving_apl
#if (DISC>0)   /* Ice discharge parameterisation */
   dis_perp_sum      = dis_perp_sum      +  dis_perp
#endif
   q_geo_sum         = q_geo_sum         + q_geo
   Q_bm_sum          = Q_bm_sum          + Q_bm
   Q_tld_sum         = Q_tld_sum         + Q_tld
   am_perp_sum       = am_perp_sum       + am_perp
   dzs_dtau_sum      = dzs_dtau_sum      + dzs_dtau
   dzm_dtau_sum      = dzm_dtau_sum      + dzm_dtau
   dzb_dtau_sum      = dzb_dtau_sum      + dzb_dtau
   dzl_dtau_sum      = dzl_dtau_sum      + dzl_dtau
   dH_c_dtau_sum     = dH_c_dtau_sum     + dH_c_dtau
   dH_t_dtau_sum     = dH_t_dtau_sum     + dH_t_dtau
   dH_dtau_sum       = dH_dtau_sum       + dH_dtau
   q_gl_g_sum        = q_gl_g_sum        + q_gl_g
      !   \!/ constant time step dtime assumed
      !       (otherwise, weighting with changing dtime would be required)

#if (defined(CLIMATOLOGY_EXTRACTION_HACK))
    !%% Climatology extraction hack (must not be used routinely)!
        temp_maat_sum = temp_maat_sum + temp_maat
        zs_sum        = zs_sum        + zs
#endif

   if (.not.flag_compute_flux_vars_only) then

      r_n_flx_ave_cnt_inv = 1.0_dp/real(n_flx_ave_cnt,dp)

      accum_flx         = accum_sum         * r_n_flx_ave_cnt_inv
      snowfall_flx      = snowfall_sum      * r_n_flx_ave_cnt_inv
      rainfall_flx      = rainfall_sum      * r_n_flx_ave_cnt_inv
      as_perp_flx       = as_perp_sum       * r_n_flx_ave_cnt_inv
      as_perp_apl_flx   = as_perp_apl_sum   * r_n_flx_ave_cnt_inv
      smb_corr_flx      = smb_corr_sum      * r_n_flx_ave_cnt_inv
      mb_source_apl_flx = mb_source_apl_sum * r_n_flx_ave_cnt_inv
      runoff_flx        = runoff_sum        * r_n_flx_ave_cnt_inv
      runoff_apl_flx    = runoff_apl_sum    * r_n_flx_ave_cnt_inv
      Q_b_tot_flx       = Q_b_tot_sum       * r_n_flx_ave_cnt_inv
      Q_b_apl_flx       = Q_b_apl_sum       * r_n_flx_ave_cnt_inv
      calving_flx       = calving_sum       * r_n_flx_ave_cnt_inv
      calving_apl_flx   = calving_apl_sum   * r_n_flx_ave_cnt_inv
#if (DISC>0)   /* Ice discharge parameterisation */
      dis_perp_flx      = dis_perp_sum      * r_n_flx_ave_cnt_inv
#endif
      q_geo_flx         = q_geo_sum         * r_n_flx_ave_cnt_inv
      Q_bm_flx          = Q_bm_sum          * r_n_flx_ave_cnt_inv
      Q_tld_flx         = Q_tld_sum         * r_n_flx_ave_cnt_inv
      am_perp_flx       = am_perp_sum       * r_n_flx_ave_cnt_inv
      dzs_dtau_flx      = dzs_dtau_sum      * r_n_flx_ave_cnt_inv
      dzm_dtau_flx      = dzm_dtau_sum      * r_n_flx_ave_cnt_inv
      dzb_dtau_flx      = dzb_dtau_sum      * r_n_flx_ave_cnt_inv
      dzl_dtau_flx      = dzl_dtau_sum      * r_n_flx_ave_cnt_inv
      dH_c_dtau_flx     = dH_c_dtau_sum     * r_n_flx_ave_cnt_inv
      dH_t_dtau_flx     = dH_t_dtau_sum     * r_n_flx_ave_cnt_inv
      dH_dtau_flx       = dH_dtau_sum       * r_n_flx_ave_cnt_inv
      q_gl_g_flx        = q_gl_g_sum        * r_n_flx_ave_cnt_inv

#if (defined(CLIMATOLOGY_EXTRACTION_HACK))
    !%% Climatology extraction hack (must not be used routinely)!
        temp_maat_flx = temp_maat_sum * r_n_flx_ave_cnt_inv
        zs_flx        = zs_sum        * r_n_flx_ave_cnt_inv
#endif

      n_flx_ave_cnt = 0

   end if   ! (.not.flag_compute_flux_vars_only)

else   ! (OUTPUT==3).and.(flag_3d_output)
       !    -> only snapshots of flux variables despite OUTPUT_FLUX_VARS==2

   accum_flx         = accum
   snowfall_flx      = snowfall
   rainfall_flx      = rainfall
   as_perp_flx       = as_perp
   as_perp_apl_flx   = as_perp_apl
   smb_corr_flx      = smb_corr
   mb_source_apl_flx = mb_source_apl
   runoff_flx        = runoff
   runoff_apl_flx    = runoff_apl
   Q_b_tot_flx       = Q_b_tot
   Q_b_apl_flx       = Q_b_apl
   calving_flx       = calving
   calving_apl_flx   = calving_apl
#if (DISC>0)   /* Ice discharge parameterisation */
   dis_perp_flx      = dis_perp
#endif
   q_geo_flx         = q_geo
   Q_bm_flx          = Q_bm
   Q_tld_flx         = Q_tld
   am_perp_flx       = am_perp
   dzs_dtau_flx      = dzs_dtau
   dzm_dtau_flx      = dzm_dtau
   dzb_dtau_flx      = dzb_dtau
   dzl_dtau_flx      = dzl_dtau
   dH_c_dtau_flx     = dH_c_dtau
   dH_t_dtau_flx     = dH_t_dtau
   dH_dtau_flx       = dH_dtau
   q_gl_g_flx        = q_gl_g

#if (defined(CLIMATOLOGY_EXTRACTION_HACK))
    !%% Climatology extraction hack (must not be used routinely)!
        temp_maat_flx = temp_maat
        zs_flx        = zs
#endif

end if

#else

errormsg = ' >>> output1: Parameter OUTPUT_FLUX_VARS must be either 1 or 2!'
call error(errormsg)

#endif

!-------- Convert data to real*4 and seconds to years --------

if (.not.flag_compute_flux_vars_only) then

year2sec_conv = year2sec

#if (!defined(OUT_TIMES) || OUT_TIMES==1)
time_conv = time*sec2year
#elif (OUT_TIMES==2)
time_conv = (time+year_zero)*sec2year
#else
errormsg = ' >>> output1: OUT_TIMES must be either 1 or 2!'
call error(errormsg)
#endif

delta_ts_conv   = delta_ts
glac_index_conv = glac_index
z_sl_mean_conv  = z_sl_mean
V_tot_conv      = V_tot
V_af_conv       = V_af
A_grounded_conv = A_grounded
A_floating_conv = A_floating

do i=0, IMAX
   xi_conv(i) = xi(i)
end do

do j=0, JMAX
   eta_conv(j) = eta(j)
end do

do kc=0, KCMAX
   sigma_level_c_conv(kc) = eaz_c_quotient(kc)
end do

do kt=0, KTMAX
   sigma_level_t_conv(kt) = zeta_t(kt)
end do

do kr=0, KRMAX
   sigma_level_r_conv(kr) = real(kr,dp)/real(KRMAX,dp)
end do

H_R_conv = real(H_R,sp)

do i=0, IMAX
do j=0, JMAX

   mask_conv(i,j)              = mask(j,i)
   mask_old_conv(i,j)          = mask_old(j,i)
   mask_ablation_type_conv(i,j) = mask_ablation_type(j,i)
   mask_region_conv(i,j)        = mask_region(j,i)
   n_cts_conv(i,j)              = n_cts(j,i)
   kc_cts_conv(i,j)             = kc_cts(j,i)

   lambda_conv(i,j)    = real(lambda(j,i),sp) 
   phi_conv(i,j)       = real(phi(j,i),sp) 
   lond_conv(i,j)      = real(lambda(j,i)*rad2deg,sp)   ! longitude in deg
   lond_conv(i,j)      = modulo(lond_conv(i,j)+180.0_sp, 360.0_sp)-180.0_sp
                                    ! mapping to interval [-180 deg, +180 deg)
   latd_conv(i,j)      = real(phi(j,i)   *rad2deg,sp)   ! latitute  in deg
   if (latd_conv(i,j) >  90.0_sp) latd_conv(i,j) =  90.0_sp
   if (latd_conv(i,j) < -90.0_sp) latd_conv(i,j) = -90.0_sp
                                 ! constraining to interval [-90 deg, +90 deg]
   cell_area_conv(i,j)   = real(cell_area(j,i),sp)
   temp_maat_conv(i,j)   = real(temp_maat(j,i),sp)
   temp_s_conv(i,j)      = real(temp_s(j,i),sp)
   accum_conv(i,j)       = real(accum_flx(j,i)*year2sec,sp)
   snowfall_conv(i,j)    = real(snowfall_flx(j,i)*year2sec,sp)
   rainfall_conv(i,j)    = real(rainfall_flx(j,i)*year2sec,sp)
   pdd_conv(i,j)         = real(ET(j,i)*one_year,sp)
   as_perp_conv(i,j)     = real(as_perp_flx(j,i)*year2sec,sp)
   as_perp_apl_conv(i,j) = real(as_perp_apl_flx(j,i)*year2sec,sp)
   smb_corr_conv(i,j)    = real(smb_corr_flx(j,i)*year2sec,sp)
   mb_source_apl_conv(i,j) = real(mb_source_apl_flx(j,i)*year2sec,sp)
   runoff_conv(i,j)        = real(runoff_flx(j,i)*year2sec,sp)
   runoff_apl_conv(i,j)    = real(runoff_apl_flx(j,i)*year2sec,sp)
   z_sl_conv(i,j)          = real(z_sl(j,i),sp)
   Q_b_tot_conv(i,j)       = real(Q_b_tot_flx(j,i)*year2sec,sp)
   Q_b_apl_conv(i,j)       = real(Q_b_apl_flx(j,i)*year2sec,sp)
   calving_conv(i,j)       = real(calving_flx(j,i)*year2sec,sp)
   calving_apl_conv(i,j)   = real(calving_apl_flx(j,i)*year2sec,sp)

#if (DISC>0)   /* Ice discharge parameterisation */
   dis_perp_conv(i,j)  = real(dis_perp_flx(j,i)*year2sec,sp)
   cst_dist_conv(i,j)  = real(cst_dist(j,i)*0.001_dp,sp)
   cos_grad_tc_conv(i,j) = real(cos_grad_tc(j,i),sp)
   mask_mar_conv(i,j)  = mask_mar(j,i)
#endif

   q_geo_conv(i,j)     = real(q_geo_flx(j,i),sp)
   zs_conv(i,j)        = real(zs(j,i),sp)
   zm_conv(i,j)        = real(zm(j,i),sp)
   zb_conv(i,j)        = real(zb(j,i),sp)
   zl_conv(i,j)        = real(zl(j,i),sp)
   zl0_conv(i,j)       = real(zl0(j,i),sp)
   H_cold_conv(i,j)    = real(H_cold(j,i),sp)
   H_temp_conv(i,j)    = real(H_temp(j,i),sp)
   H_conv(i,j)         = real(H(j,i),sp)
   Q_bm_conv(i,j)      = real(Q_bm_flx(j,i)*year2sec,sp)
   Q_tld_conv(i,j)     = real(Q_tld_flx(j,i)*year2sec,sp)
   am_perp_conv(i,j)   = real(am_perp_flx(j,i)*year2sec,sp)
   qx_conv(i,j)        = real(qx(j,i)*year2sec,sp)
   qy_conv(i,j)        = real(qy(j,i)*year2sec,sp)
   vx_m_sia_conv(i,j)  = real(vx_m_sia(j,i)*year2sec,sp)
   vy_m_sia_conv(i,j)  = real(vy_m_sia(j,i)*year2sec,sp)
   vx_m_ssa_conv(i,j)  = real(vx_m_ssa(j,i)*year2sec,sp)
   vy_m_ssa_conv(i,j)  = real(vy_m_ssa(j,i)*year2sec,sp)
   dzs_dtau_conv(i,j)  = real(dzs_dtau_flx(j,i)*year2sec,sp)
   dzm_dtau_conv(i,j)  = real(dzm_dtau_flx(j,i)*year2sec,sp)
   dzb_dtau_conv(i,j)  = real(dzb_dtau_flx(j,i)*year2sec,sp)
   dzl_dtau_conv(i,j)  = real(dzl_dtau_flx(j,i)*year2sec,sp)
   dH_c_dtau_conv(i,j) = real(dH_c_dtau_flx(j,i)*year2sec,sp)
   dH_t_dtau_conv(i,j) = real(dH_t_dtau_flx(j,i)*year2sec,sp)
   dH_dtau_conv(i,j)   = real(dH_dtau_flx(j,i)*year2sec,sp)
   vx_b_g_conv(i,j)    = real(vx_b_g(j,i)*year2sec,sp)
   vy_b_g_conv(i,j)    = real(vy_b_g(j,i)*year2sec,sp)
   vz_b_conv(i,j)      = real(vz_b(j,i)*year2sec,sp)
   vh_b_conv(i,j)      = sqrt( vx_b_g_conv(i,j)**2 + vy_b_g_conv(i,j)**2 )
   vx_s_g_conv(i,j)    = real(vx_s_g(j,i)*year2sec,sp)
   vy_s_g_conv(i,j)    = real(vy_s_g(j,i)*year2sec,sp)
   vz_s_conv(i,j)      = real(vz_s(j,i)*year2sec,sp)
   vh_s_conv(i,j)      = sqrt( vx_s_g_conv(i,j)**2 + vy_s_g_conv(i,j)**2 )
   vx_m_g_conv(i,j)    = real(vx_m_g(j,i)*year2sec,sp)
   vy_m_g_conv(i,j)    = real(vy_m_g(j,i)*year2sec,sp)
   vh_m_conv(i,j)      = sqrt( vx_m_g_conv(i,j)**2 + vy_m_g_conv(i,j)**2 )
   c_slide_conv(i,j)   = real(c_slide(j,i)*year2sec,sp)
   temp_b_conv(i,j)    = real(temp_b(j,i),sp)
   temph_b_conv(i,j)   = real(temph_b(j,i),sp)
   tau_dr_conv(i,j)    = real(tau_dr(j,i),sp)
   tau_b_conv(i,j)     = real(tau_b(j,i),sp)
   p_b_w_conv(i,j)     = real(p_b_w(j,i),sp)
   q_w_conv(i,j)       = real(q_w(j,i)*year2sec,sp)
   q_w_x_conv(i,j)     = real(q_w_x(j,i)*year2sec,sp)
   q_w_y_conv(i,j)     = real(q_w_y(j,i)*year2sec,sp)
   H_w_conv(i,j)       = real(H_w(j,i),sp)
   q_gl_g_conv(i,j)    = real(q_gl_g_flx(j,i)*year2sec,sp)
   ratio_sl_x_conv(i,j) = real(ratio_sl_x(j,i),sp)
   ratio_sl_y_conv(i,j) = real(ratio_sl_y(j,i),sp)
   ratio_sl_conv(i,j)   = real(ratio_sl(j,i),sp)

#if (defined(CLIMATOLOGY_EXTRACTION_HACK))
    !%% Climatology extraction hack (must not be used routinely)!
        temp_maat_conv(i,j) = real(temp_maat_flx(j,i),sp)
        zs_conv(i,j)        = real(zs_flx(j,i),sp)
#endif

   if (flag_shelfy_stream_x(j,i)) then
      flag_shelfy_stream_x_conv(i,j) = 1
   else
      flag_shelfy_stream_x_conv(i,j) = 0
   end if

   if (flag_shelfy_stream_y(j,i)) then
      flag_shelfy_stream_y_conv(i,j) = 1
   else
      flag_shelfy_stream_y_conv(i,j) = 0
   end if

   if (flag_shelfy_stream(j,i)) then
      flag_shelfy_stream_conv(i,j) = 1
   else
      flag_shelfy_stream_conv(i,j) = 0
   end if

   if (flag_grounding_line_1(j,i)) then
      flag_grounding_line_1_conv(i,j) = 1
   else
      flag_grounding_line_1_conv(i,j) = 0
   end if

   if (flag_grounding_line_2(j,i)) then
      flag_grounding_line_2_conv(i,j) = 1
   else
      flag_grounding_line_2_conv(i,j) = 0
   end if

   if (flag_calving_front_1(j,i)) then
      flag_calving_front_1_conv(i,j) = 1
   else
      flag_calving_front_1_conv(i,j) = 0
   end if

   if (flag_calving_front_2(j,i)) then
      flag_calving_front_2_conv(i,j) = 1
   else
      flag_calving_front_2_conv(i,j) = 0
   end if

   if (flag_grounded_front_a_1(j,i)) then
      flag_grounded_front_a_1_conv(i,j) = 1
   else
      flag_grounded_front_a_1_conv(i,j) = 0
   end if

   if (flag_grounded_front_a_2(j,i)) then
      flag_grounded_front_a_2_conv(i,j) = 1
   else
      flag_grounded_front_a_2_conv(i,j) = 0
   end if

   if (flag_grounded_front_b_1(j,i)) then
      flag_grounded_front_b_1_conv(i,j) = 1
   else
      flag_grounded_front_b_1_conv(i,j) = 0
   end if

   if (flag_grounded_front_b_2(j,i)) then
      flag_grounded_front_b_2_conv(i,j) = 1
   else
      flag_grounded_front_b_2_conv(i,j) = 0
   end if

   vis_ave_g_conv(i,j) = real(vis_ave_g(j,i),sp)
   vis_int_g_conv(i,j) = real(vis_int_g(j,i),sp)

   do kr=0, KRMAX
      temp_r_conv(i,j,kr) = real(temp_r(kr,j,i),sp)
   end do

   do kt=0, KTMAX
      vx_t_conv(i,j,kt)    = real(vx_t(kt,j,i)*year2sec,sp)
      vy_t_conv(i,j,kt)    = real(vy_t(kt,j,i)*year2sec,sp)
      vz_t_conv(i,j,kt)    = real(vz_t(kt,j,i)*year2sec,sp)
      omega_t_conv(i,j,kt) = real(omega_t(kt,j,i),sp)
      age_t_conv(i,j,kt)   = real(age_t(kt,j,i)*sec2year,sp)
      enth_t_conv(i,j,kt)  = real(enth_t(kt,j,i),sp)
      enh_t_conv(i,j,kt)   = real(enh_t(kt,j,i),sp)
      strain_heating_t_conv(i,j,kt) = real(strain_heating_t(kt,j,i),sp)
   end do

   do kc=0, KCMAX
      vx_c_conv(i,j,kc)    = real(vx_c(kc,j,i)*year2sec,sp)
      vy_c_conv(i,j,kc)    = real(vy_c(kc,j,i)*year2sec,sp)
      vz_c_conv(i,j,kc)    = real(vz_c(kc,j,i)*year2sec,sp)
      temp_c_conv(i,j,kc)  = real(temp_c(kc,j,i),sp)
      age_c_conv(i,j,kc)   = real(age_c(kc,j,i)*sec2year,sp)
      enth_c_conv(i,j,kc)  = real(enth_c(kc,j,i),sp)
      omega_c_conv(i,j,kc) = real(omega_c(kc,j,i),sp)
      enh_c_conv(i,j,kc)   = real(enh_c(kc,j,i),sp)
      strain_heating_c_conv(i,j,kc) = real(strain_heating_c(kc,j,i),sp)
   end do

end do
end do

end if   ! (.not.flag_compute_flux_vars_only)

!-------- Write data on file --------

if (.not.flag_compute_flux_vars_only) then

call check( nf90_inq_varid(ncid, 'mapping', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, 0), thisroutine )

call check( nf90_inq_varid(ncid, 'year2sec', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, year2sec_conv), thisroutine )

call check( nf90_inq_varid(ncid, 'time', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, time_conv), thisroutine )

if ((forcing_flag == 1).or.(forcing_flag == 3)) then

   call check( nf90_inq_varid(ncid, 'delta_ts', ncv), thisroutine )
   call check( nf90_put_var(ncid, ncv, delta_ts_conv), thisroutine )

else if (forcing_flag == 2) then

   call check( nf90_inq_varid(ncid, 'glac_index', ncv), thisroutine )
   call check( nf90_put_var(ncid, ncv, glac_index_conv), thisroutine )

end if

call check( nf90_inq_varid(ncid, 'z_sl_mean', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, z_sl_mean_conv), thisroutine )

call check( nf90_inq_varid(ncid, 'V_tot', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, V_tot_conv), thisroutine )

call check( nf90_inq_varid(ncid, 'V_af', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, V_af_conv), thisroutine )

call check( nf90_inq_varid(ncid, 'A_grounded', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, A_grounded_conv), thisroutine )

call check( nf90_inq_varid(ncid, 'A_floating', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, A_floating_conv), thisroutine )

call check( nf90_inq_varid(ncid, 'x', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, xi_conv, &
                         start=nc1cor_i, count=nc1cnt_i), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'y', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, eta_conv, &
                         start=nc1cor_j, count=nc1cnt_j), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'sigma_level_c', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, sigma_level_c_conv, &
                         start=nc1cor_kc, count=nc1cnt_kc), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'sigma_level_t', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, sigma_level_t_conv, &
                         start=nc1cor_kt, count=nc1cnt_kt), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'sigma_level_r', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, sigma_level_r_conv, &
                         start=nc1cor_kr, count=nc1cnt_kr), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'lon', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, lond_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'lat', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, latd_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'lambda', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, lambda_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'phi', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, phi_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'cell_area', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, cell_area_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'temp_maat', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, temp_maat_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'temp_s', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, temp_s_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'prec', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, accum_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'snowfall', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, snowfall_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'rainfall', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, rainfall_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'pdd', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, pdd_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'as_perp', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, as_perp_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'as_perp_apl', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, as_perp_apl_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'smb_corr', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, smb_corr_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'mb_source_apl', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, mb_source_apl_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'runoff', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, runoff_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'runoff_apl', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, runoff_apl_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'z_sl', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, z_sl_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'Q_b_tot', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, Q_b_tot_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'Q_b_apl', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, Q_b_apl_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'calving', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, calving_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'calving_apl', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, calving_apl_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

#if (DISC>0)   /* Ice discharge parameterisation */

call check( nf90_inq_varid(ncid, 'dis_perp', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, dis_perp_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'cst_dist', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, cst_dist_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'cos_grad_tc', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, cos_grad_tc_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'mask_mar', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, mask_mar_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

#endif

call check( nf90_inq_varid(ncid, 'q_geo', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, q_geo_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'mask', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, mask_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'mask_old', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, mask_old_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'mask_ablation_type', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, mask_ablation_type_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

if (maxval(mask_region) > 0) then
   call check( nf90_inq_varid(ncid, 'mask_region', ncv), thisroutine )
   call check( nf90_put_var(ncid, ncv, mask_region_conv, &
                            start=nc2cor_ij, count=nc2cnt_ij), &
               thisroutine )
end if

call check( nf90_inq_varid(ncid, 'n_cts', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, n_cts_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'kc_cts', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, kc_cts_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'zs', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, zs_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'zm', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, zm_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'zb', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, zb_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'zl', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, zl_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'zl0', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, zl0_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'H_cold', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, H_cold_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'H_temp', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, H_temp_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'H', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, H_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'H_R', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, H_R_conv), thisroutine )

call check( nf90_inq_varid(ncid, 'Q_bm', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, Q_bm_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'Q_tld', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, Q_tld_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'am_perp', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, am_perp_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'qx', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, qx_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'qy', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, qy_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'vx_m_sia', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, vx_m_sia_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'vy_m_sia', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, vy_m_sia_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'vx_m_ssa', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, vx_m_ssa_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'vy_m_ssa', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, vy_m_ssa_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'dzs_dt', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, dzs_dtau_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'dzm_dt', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, dzm_dtau_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'dzb_dt', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, dzb_dtau_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'dzl_dt', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, dzl_dtau_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'dH_c_dt', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, dH_c_dtau_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'dH_t_dt', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, dH_t_dtau_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'dH_dt', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, dH_dtau_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'vx_b_g', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, vx_b_g_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'vy_b_g', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, vy_b_g_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'vz_b', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, vz_b_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'vh_b', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, vh_b_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'vx_s_g', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, vx_s_g_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'vy_s_g', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, vy_s_g_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'vz_s', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, vz_s_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'vh_s', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, vh_s_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'vx_m_g', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, vx_m_g_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'vy_m_g', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, vy_m_g_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'vh_m', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, vh_m_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'c_slide', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, c_slide_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'temp_b', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, temp_b_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'temph_b', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, temph_b_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'tau_dr', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, tau_dr_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'tau_b', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, tau_b_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'p_b_w', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, p_b_w_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'q_w', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, q_w_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'q_w_x', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, q_w_x_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'q_w_y', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, q_w_y_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'H_w', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, H_w_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'q_gl_g', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, q_gl_g_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'ratio_sl_x', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, ratio_sl_x_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'ratio_sl_y', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, ratio_sl_y_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'ratio_sl', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, ratio_sl_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'flag_shelfy_stream_x', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, flag_shelfy_stream_x_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'flag_shelfy_stream_y', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, flag_shelfy_stream_y_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'flag_shelfy_stream', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, flag_shelfy_stream_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'flag_grounding_line_1', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, flag_grounding_line_1_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'flag_grounding_line_2', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, flag_grounding_line_2_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'flag_calving_front_1', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, flag_calving_front_1_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'flag_calving_front_2', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, flag_calving_front_2_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'flag_grounded_front_a_1', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, flag_grounded_front_a_1_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'flag_grounded_front_a_2', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, flag_grounded_front_a_2_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'flag_grounded_front_b_1', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, flag_grounded_front_b_1_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'flag_grounded_front_b_2', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, flag_grounded_front_b_2_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'vis_ave_g', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, vis_ave_g_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

call check( nf90_inq_varid(ncid, 'vis_int_g', ncv), thisroutine )
call check( nf90_put_var(ncid, ncv, vis_int_g_conv, &
                         start=nc2cor_ij, count=nc2cnt_ij), &
            thisroutine )

if (flag_3d_output) then

   call check( nf90_inq_varid(ncid, 'vx_c', ncv), thisroutine )
   call check( nf90_put_var(ncid, ncv, vx_c_conv, &
                            start=nc3cor_ijkc, count=nc3cnt_ijkc), &
               thisroutine )

   call check( nf90_inq_varid(ncid, 'vy_c', ncv), thisroutine )
   call check( nf90_put_var(ncid, ncv, vy_c_conv, &
                            start=nc3cor_ijkc, count=nc3cnt_ijkc), &
               thisroutine )

   call check( nf90_inq_varid(ncid, 'vz_c', ncv), thisroutine )
   call check( nf90_put_var(ncid, ncv, vz_c_conv, &
                            start=nc3cor_ijkc, count=nc3cnt_ijkc), &
               thisroutine )

   call check( nf90_inq_varid(ncid, 'vx_t', ncv), thisroutine )
   call check( nf90_put_var(ncid, ncv, vx_t_conv, &
                            start=nc3cor_ijkt, count=nc3cnt_ijkt), &
               thisroutine )

   call check( nf90_inq_varid(ncid, 'vy_t', ncv), thisroutine )
   call check( nf90_put_var(ncid, ncv, vy_t_conv, &
                            start=nc3cor_ijkt, count=nc3cnt_ijkt), &
               thisroutine )

   call check( nf90_inq_varid(ncid, 'vz_t', ncv), thisroutine )
   call check( nf90_put_var(ncid, ncv, vz_t_conv, &
                            start=nc3cor_ijkt, count=nc3cnt_ijkt), &
               thisroutine )

   call check( nf90_inq_varid(ncid, 'temp_c', ncv), thisroutine )
   call check( nf90_put_var(ncid, ncv, temp_c_conv, &
                            start=nc3cor_ijkc, count=nc3cnt_ijkc), &
               thisroutine )

   call check( nf90_inq_varid(ncid, 'omega_t', ncv), thisroutine )
   call check( nf90_put_var(ncid, ncv, omega_t_conv, &
                            start=nc3cor_ijkt, count=nc3cnt_ijkt), &
               thisroutine )

   call check( nf90_inq_varid(ncid, 'temp_r', ncv), thisroutine )
   call check( nf90_put_var(ncid, ncv, temp_r_conv, &
                            start=nc3cor_ijkr, count=nc3cnt_ijkr), &
               thisroutine )

   call check( nf90_inq_varid(ncid, 'enth_c', ncv), thisroutine )
   call check( nf90_put_var(ncid, ncv, enth_c_conv, &
                            start=nc3cor_ijkc, count=nc3cnt_ijkc), &
               thisroutine )

   call check( nf90_inq_varid(ncid, 'enth_t', ncv), thisroutine )
   call check( nf90_put_var(ncid, ncv, enth_t_conv, &
                            start=nc3cor_ijkt, count=nc3cnt_ijkt), &
               thisroutine )

   call check( nf90_inq_varid(ncid, 'omega_c', ncv), thisroutine )
   call check( nf90_put_var(ncid, ncv, omega_c_conv, &
                            start=nc3cor_ijkc, count=nc3cnt_ijkc), &
               thisroutine )

   call check( nf90_inq_varid(ncid, 'enh_c', ncv), thisroutine )
   call check( nf90_put_var(ncid, ncv, enh_c_conv, &
                            start=nc3cor_ijkc, count=nc3cnt_ijkc), &
               thisroutine )

   call check( nf90_inq_varid(ncid, 'enh_t', ncv), thisroutine )
   call check( nf90_put_var(ncid, ncv, enh_t_conv, &
                            start=nc3cor_ijkt, count=nc3cnt_ijkt), &
               thisroutine )

   call check( nf90_inq_varid(ncid, 'strain_heating_c', ncv), thisroutine )
   call check( nf90_put_var(ncid, ncv, strain_heating_c_conv, &
                            start=nc3cor_ijkc, count=nc3cnt_ijkc), &
               thisroutine )

   call check( nf90_inq_varid(ncid, 'strain_heating_t', ncv), thisroutine )
   call check( nf90_put_var(ncid, ncv, strain_heating_t_conv, &
                            start=nc3cor_ijkt, count=nc3cnt_ijkt), &
               thisroutine )

   call check( nf90_inq_varid(ncid, 'age_c', ncv), thisroutine )
   call check( nf90_put_var(ncid, ncv, age_c_conv, &
                            start=nc3cor_ijkc, count=nc3cnt_ijkc), &
               thisroutine )

   call check( nf90_inq_varid(ncid, 'age_t', ncv), thisroutine )
   call check( nf90_put_var(ncid, ncv, age_t_conv, &
                            start=nc3cor_ijkt, count=nc3cnt_ijkt), &
               thisroutine )

end if

end if   ! (.not.flag_compute_flux_vars_only)

!-------- Close file --------

if (.not.flag_compute_flux_vars_only) then

call check( nf90_sync(ncid),  thisroutine )
call check( nf90_close(ncid), thisroutine )

deallocate(coord_id)

end if   ! (.not.flag_compute_flux_vars_only)

!-------- Increase file counter --------

if (.not.flag_compute_flux_vars_only) then

ndat = ndat+1

if (flag_3d_output) then
   ndat3d = ndat
else
   ndat2d = ndat
end if

end if   ! (.not.flag_compute_flux_vars_only)

if (firstcall%output1) firstcall%output1 = .false.

end subroutine output1

!-------------------------------------------------------------------------------
!> Writing of time-series data on file in ASCII and NetCDF format.
!-------------------------------------------------------------------------------
subroutine output2(time, dxi, deta, delta_ts, glac_index, &
                   opt_flag_compute_flux_vars_only)

  use netcdf
  use nc_check_m

#if (DISC>0)
  use discharge_workers_m, only: dT_glann, dT_sub
#endif

implicit none

real(dp), intent(in) :: time, dxi, deta, delta_ts, glac_index

logical, optional, intent(in) :: opt_flag_compute_flux_vars_only

integer(i4b) :: i, j, n
integer(i4b) :: ios
integer(i4b) :: n_base, n_tempbase
real(dp) :: time_val, &
            V_tot, V_grounded, V_floating, &
            A_tot, A_grounded, A_floating, &
            V_af, V_sle, V_temp, A_temp, &
            H_max, H_t_max, zs_max, vs_max, Tbh_max, &
            dV_dt, Q_s, precip_tot, runoff_tot, &
            Q_b, Q_temp, bmb_tot, bmb_gr_tot, bmb_fl_tot, &
            calv_tot, mbp, mb_resid, mb_mis, disc_lsc, disc_ssc
real(dp) :: x_pos, y_pos
real(dp), dimension(0:JMAX,0:IMAX) :: H_cold, H_temp
real(dp) :: Tbh_help
real(dp) :: H_ave_sed, Tbh_ave_sed, Atb_sed
real(dp) :: sum_area_sed
logical :: flag_compute_flux_vars_only
logical, dimension(0:JMAX,0:IMAX) :: flag_region

integer(i4b)       :: cmode
integer(i4b)       :: n_deflate_level
logical            :: flag_shuffle
integer(i4b), dimension(0:99), save :: ncid
integer(i4b)       :: ncd, ncv, nc1d
integer(i4b)       :: nc1cor(1), nc1cnt(1)
integer(i4b)       :: n_sync
real(dp), save     :: time_add_offset_val
character(len= 16) :: ch_date, ch_time, ch_zone
character(len=256) :: filename, filename_with_path, buffer
character(len=  2) :: ch2_aux
logical, save      :: grads_nc_tweaks

integer(i4b), save :: counter = 0

character(len=64), parameter :: thisroutine = 'output2'

#if (OUTPUT_FLUX_VARS==2)   /* averaging of flux variables */

integer(i4b), dimension(0:99), save :: n_flx_ave_cnt = 0
real(dp)                            :: r_n_flx_ave_cnt_inv

real(dp), dimension(0:99), save :: Q_s_sum         = 0.0_dp, &
                                   precip_tot_sum  = 0.0_dp, &
                                   runoff_tot_sum  = 0.0_dp, &
                                   bmb_tot_sum     = 0.0_dp, &
#if (MARGIN==3)
                                   bmb_gr_tot_sum  = 0.0_dp, &
                                   bmb_fl_tot_sum  = 0.0_dp, &
#endif
                                   Q_b_sum         = 0.0_dp, &
                                   Q_temp_sum      = 0.0_dp, &
                                   calv_tot_sum    = 0.0_dp, &
#if (DISC>0)
                                   disc_lsc_sum    = 0.0_dp, &
                                   disc_ssc_sum    = 0.0_dp, &
#endif
                                   dV_dt_sum       = 0.0_dp, &
                                   mb_resid_sum    = 0.0_dp, &
                                   mb_mis_sum      = 0.0_dp, &
                                   mbp_sum         = 0.0_dp

#endif

real(dp) :: Q_s_flx        , &
            precip_tot_flx , &
            runoff_tot_flx , &
            bmb_tot_flx    , &
#if (MARGIN==3)
            bmb_gr_tot_flx , &
            bmb_fl_tot_flx , &
#endif
            Q_b_flx        , &
            Q_temp_flx     , &
            calv_tot_flx   , &
#if (DISC>0)
            disc_lsc_flx   , &
            disc_ssc_flx   , &
#endif
            dV_dt_flx      , &
            mb_resid_flx   , &
            mb_mis_flx     , &
            mbp_flx

character(len=128), parameter :: &
                    fmt1 = '(1pe13.6,2(1pe13.4),/,' // &
                           '13x,6(1pe13.4),/,'      // &
                           '26x,3(1pe13.4),/,'      // &
                           '26x,5(1pe13.4),/)'

character(len=128), parameter :: &
                    fmt2 = '(1pe13.6,2(1pe13.4),/,' // &
                           '13x,3(1pe13.4),/)'

character(len=128), parameter :: &
                    fmt3 = '(1pe13.6,3(1pe13.4),/)'

if (present(opt_flag_compute_flux_vars_only)) then
   flag_compute_flux_vars_only = opt_flag_compute_flux_vars_only
else
   flag_compute_flux_vars_only = .false.
end if

if (.not.flag_compute_flux_vars_only) &
   counter = counter + 1

!-------- Thickness of the cold and temperate layers --------

H_cold = 0.0_dp
H_temp = 0.0_dp

#if (CALCMOD==1)
do i=1, IMAX-1
do j=1, JMAX-1
   H_temp(j,i) = H_t(j,i)
end do
end do
#elif (CALCMOD==0 || CALCMOD==2 || CALCMOD==3 || CALCMOD==-1)
do i=1, IMAX-1
do j=1, JMAX-1
   H_temp(j,i) = H_c(j,i)*eaz_c_quotient(kc_cts(j,i))
end do
end do
#endif

H_cold = H - H_temp

!-------- Begin loop over regions --------

if (maxval(mask_region) > 99) then
   errormsg = ' >>> output2: Not more than 99 regions allowed!'
   call error(errormsg)
end if

do n=0, maxval(mask_region)   ! n=0: entire ice sheet, n>0: defined regions

!-------- Computing the scalar output variables --------

   if (n==0) then

      call scalar_variables(time, &
                            H_cold, H_temp, &
                            time_val, &
                            V_tot, V_grounded, V_floating, &
                            A_tot, A_grounded, A_floating, &
                            V_af, V_sle, V_temp, A_temp, &
                            H_max, H_t_max, zs_max, vs_max, Tbh_max, &
                            dV_dt, Q_s, precip_tot, runoff_tot, &
                            Q_b, Q_temp, bmb_tot, bmb_gr_tot, bmb_fl_tot, &
                            calv_tot, disc_lsc, disc_ssc, &
                            mbp, mb_resid, mb_mis)

   else

      flag_region = .false.
      flag_region = (mask_region==n)

      call scalar_variables(time, &
                            H_cold, H_temp, &
                            time_val, &
                            V_tot, V_grounded, V_floating, &
                            A_tot, A_grounded, A_floating, &
                            V_af, V_sle, V_temp, A_temp, &
                            H_max, H_t_max, zs_max, vs_max, Tbh_max, &
                            dV_dt, Q_s, precip_tot, runoff_tot, &
                            Q_b, Q_temp, bmb_tot, bmb_gr_tot, bmb_fl_tot, &
                            calv_tot, disc_lsc, disc_ssc, &
                            mbp, mb_resid, mb_mis, &
                            opt_flag_region=flag_region)

   end if

!-------- Flux variables --------

#if (!defined(OUTPUT_FLUX_VARS) || OUTPUT_FLUX_VARS==1)
       ! snapshots of flux variables

   Q_s_flx         = Q_s
   precip_tot_flx  = precip_tot
   runoff_tot_flx  = runoff_tot
   bmb_tot_flx     = bmb_tot
#if (MARGIN==3)
   bmb_gr_tot_flx  = bmb_gr_tot
   bmb_fl_tot_flx  = bmb_fl_tot
#endif
   Q_b_flx         = Q_b
   Q_temp_flx      = Q_temp
   calv_tot_flx    = calv_tot
#if (DISC>0)
   disc_lsc_flx    = disc_lsc
   disc_ssc_flx    = disc_ssc
#endif
   dV_dt_flx       = dV_dt
   mb_resid_flx    = mb_resid
   mb_mis_flx      = mb_mis
   mbp_flx         = mbp

#elif (OUTPUT_FLUX_VARS==2)
       ! averaging of flux variables

   if (n_flx_ave_cnt(n)==0) then

      Q_s_sum(n)         = 0.0_dp
      precip_tot_sum(n)  = 0.0_dp
      runoff_tot_sum(n)  = 0.0_dp
      bmb_tot_sum(n)     = 0.0_dp
#if (MARGIN==3)
      bmb_gr_tot_sum(n)  = 0.0_dp
      bmb_fl_tot_sum(n)  = 0.0_dp
#endif
      Q_b_sum(n)         = 0.0_dp
      Q_temp_sum(n)      = 0.0_dp
      calv_tot_sum(n)    = 0.0_dp
#if (DISC>0)
      disc_lsc_sum(n)    = 0.0_dp
      disc_ssc_sum(n)    = 0.0_dp
#endif
      dV_dt_sum(n)       = 0.0_dp
      mb_resid_sum(n)    = 0.0_dp
      mb_mis_sum(n)      = 0.0_dp
      mbp_sum(n)         = 0.0_dp

   end if

   n_flx_ave_cnt(n) = n_flx_ave_cnt(n) + 1

   Q_s_sum(n)         = Q_s_sum(n)        + Q_s
   precip_tot_sum(n)  = precip_tot_sum(n) + precip_tot
   runoff_tot_sum(n)  = runoff_tot_sum(n) + runoff_tot
   bmb_tot_sum(n)     = bmb_tot_sum(n)    + bmb_tot
#if (MARGIN==3)
   bmb_gr_tot_sum(n)  = bmb_gr_tot_sum(n) + bmb_gr_tot
   bmb_fl_tot_sum(n)  = bmb_fl_tot_sum(n) + bmb_fl_tot
#endif
   Q_b_sum(n)         = Q_b_sum(n)        + Q_b
   Q_temp_sum(n)      = Q_temp_sum(n)     + Q_temp
   calv_tot_sum(n)    = calv_tot_sum(n)   + calv_tot
#if (DISC>0)
   disc_lsc_sum(n)    = disc_lsc_sum(n)   + disc_lsc
   disc_ssc_sum(n)    = disc_ssc_sum(n)   + disc_ssc
#endif
   dV_dt_sum(n)       = dV_dt_sum(n)      + dV_dt
   mb_resid_sum(n)    = mb_resid_sum(n)   + mb_resid
   mb_mis_sum(n)      = mb_mis_sum(n)     + mb_mis
   mbp_sum(n)         = mbp_sum(n)        + mbp
      !   \!/ constant time step dtime assumed
      !       (otherwise, weighting with changing dtime would be required)

   if (.not.flag_compute_flux_vars_only) then

      r_n_flx_ave_cnt_inv = 1.0_dp/real(n_flx_ave_cnt(n),dp)

      Q_s_flx         = Q_s_sum(n)        * r_n_flx_ave_cnt_inv
      precip_tot_flx  = precip_tot_sum(n) * r_n_flx_ave_cnt_inv
      runoff_tot_flx  = runoff_tot_sum(n) * r_n_flx_ave_cnt_inv
      bmb_tot_flx     = bmb_tot_sum(n)    * r_n_flx_ave_cnt_inv
#if (MARGIN==3)
      bmb_gr_tot_flx  = bmb_gr_tot_sum(n) * r_n_flx_ave_cnt_inv
      bmb_fl_tot_flx  = bmb_fl_tot_sum(n) * r_n_flx_ave_cnt_inv
#endif
      Q_b_flx         = Q_b_sum(n)        * r_n_flx_ave_cnt_inv
      Q_temp_flx      = Q_temp_sum(n)     * r_n_flx_ave_cnt_inv
      calv_tot_flx    = calv_tot_sum(n)   * r_n_flx_ave_cnt_inv
#if (DISC>0)
      disc_lsc_flx    = disc_lsc_sum(n)   * r_n_flx_ave_cnt_inv
      disc_ssc_flx    = disc_ssc_sum(n)   * r_n_flx_ave_cnt_inv
#endif
      dV_dt_flx       = dV_dt_sum(n)      * r_n_flx_ave_cnt_inv
      mb_resid_flx    = mb_resid_sum(n)   * r_n_flx_ave_cnt_inv
      mb_mis_flx      = mb_mis_sum(n)     * r_n_flx_ave_cnt_inv
      mbp_flx         = mbp_sum(n)        * r_n_flx_ave_cnt_inv

      n_flx_ave_cnt(n) = 0

   end if   ! (.not.flag_compute_flux_vars_only)

#else

   errormsg = ' >>> output2: Parameter OUTPUT_FLUX_VARS must be either 1 or 2!'
   call error(errormsg)

#endif

!-------- Writing of data on ASCII time-series file --------

   if ((n==0).and.(.not.flag_compute_flux_vars_only)) then

      if ((forcing_flag == 1).or.(forcing_flag == 3)) then

         write(unit=12, fmt=trim(fmt1)) time_val, delta_ts, z_sl_mean, &
            V_tot, V_grounded, V_floating, A_tot, A_grounded, A_floating, &
            V_sle, V_temp, A_temp, &
            H_max, H_t_max, zs_max, vs_max, Tbh_max

      else if (forcing_flag == 2) then

         write(unit=12, fmt=trim(fmt1)) time_val, glac_index, z_sl_mean, &
            V_tot, V_grounded, V_floating, A_tot, A_grounded, A_floating, &
            V_sle, V_temp, A_temp, &
            H_max, H_t_max, zs_max, vs_max, Tbh_max

      end if

   end if   ! ((n==0).and.(.not.flag_compute_flux_vars_only))

!-------- Special output (ASCII) --------

!  ------ Time-series data for the sediment region of HEINO

#if (defined(HEINO))

!    ---- Average ice thickness, average basal temperature rel. to pmp,
!         temperate basal area

   if ((n==0).and.(.not.flag_compute_flux_vars_only)) then

      H_ave_sed    = 0.0_dp
      Tbh_ave_sed  = 0.0_dp
      Atb_sed      = 0.0_dp
      sum_area_sed = 0.0_dp

      do i=1, IMAX-1
      do j=1, JMAX-1

         if (n_slide_region(j,i)==3) then   ! sediment region

            sum_area_sed = sum_area_sed + cell_area(j,i)

            H_ave_sed = H_ave_sed + cell_area(j,i)*H(j,i)

            if (n_cts(j,i) /= -1) then   ! temperate base
               Tbh_help    = 0.0_dp
            else   ! cold base
               Tbh_help    = min((temp_c(0,j,i)-temp_c_m(0,j,i)), 0.0_dp)
            end if
            Tbh_ave_sed = Tbh_ave_sed + cell_area(j,i)*Tbh_help

            if (n_cts(j,i) /= -1) Atb_sed = Atb_sed + cell_area(j,i)

         end if

      end do
      end do

      if (sum_area_sed > eps) then
         H_ave_sed   = H_ave_sed   / sum_area_sed
         Tbh_ave_sed = Tbh_ave_sed / sum_area_sed
      else
         errormsg = ' >>> output2: No sediment area found!'
         call error(errormsg)
      end if

!    ---- Writing of data on time-series file

      if ((forcing_flag == 1).or.(forcing_flag == 3)) then
         write(unit=15, fmt=trim(fmt2)) time_val, delta_ts, z_sl_mean, &
                                        H_ave_sed, Tbh_ave_sed, Atb_sed
      else if (forcing_flag == 2) then
         write(unit=15, fmt=trim(fmt2)) time_val, glac_index, z_sl_mean, &
                                        H_ave_sed, Tbh_ave_sed, Atb_sed
      end if

   end if   ! ((n==0).and.(.not.flag_compute_flux_vars_only))

#endif

!-------- Extended time-series file in NetCDF format --------

   if (firstcall%output2) then

!  ------ Open NetCDF file

      if (n==0) then
         filename = trim(run_name)//'_ser.nc'
      else
         write(ch2_aux, '(i0.2)') n
         filename = trim(run_name)//'_ser_region'//ch2_aux//'.nc'
      end if

      filename_with_path = trim(OUT_PATH)//'/'//trim(filename)

      !%% call set_cmode(cmode, n_deflate_level, flag_shuffle)

      cmode           = NF90_NOCLOBBER
      n_deflate_level = 0
      flag_shuffle    = .false.

      ios = nf90_create(trim(filename_with_path), cmode, ncid(n))

      if (ios /= nf90_noerr) then
         errormsg = ' >>> output2: Error when opening the' &
                  //               end_of_line &
                  //'              NetCDF time-series file!'
         call error(errormsg)
      end if

      ncid_ser(n) = ncid(n)

!  ------ Global attributes

      buffer = 'Time-series output of simulation '//trim(run_name)
      if (n>0) buffer = trim(buffer)//' (region '//ch2_aux//')'
      call check( nf90_put_att(ncid(n), NF90_GLOBAL, 'title', trim(buffer)), &
                  thisroutine )

      buffer = 'SICOPOLIS Version '//VERSION
      call check( nf90_put_att(ncid(n), NF90_GLOBAL, 'source', trim(buffer)), &
                  thisroutine )

      call date_and_time(ch_date, ch_time, ch_zone)
      buffer = ch_date(1:4)//'-'//ch_date(5:6)//'-'//ch_date(7:8)//' '// &
               ch_time(1:2)//':'//ch_time(3:4)//':'//ch_time(5:6)//' '// &
               ch_zone(1:3)//':'//ch_zone(4:5)//' - Data produced'
      call check( nf90_put_att(ncid(n), NF90_GLOBAL, 'history', trim(buffer)), &
                  thisroutine )

      call set_ch_website(buffer)
      call check( nf90_put_att(ncid(n), NF90_GLOBAL, 'references', trim(buffer)), &
                  thisroutine )

      call set_ch_institution(buffer)
      call check( nf90_put_att(ncid(n), NF90_GLOBAL, 'institution', trim(buffer)), &
                  thisroutine )

!  ------ Definition of the dimensions

      call set_grads_nc_tweaks(grads_nc_tweaks)

      if (grads_nc_tweaks) then
         call check( nf90_def_dim(ncid(n), 'x', 1, ncd), thisroutine )
         call check( nf90_def_dim(ncid(n), 'y', 1, ncd), thisroutine )
      end if

      call check( nf90_def_dim(ncid(n), 't', NF90_UNLIMITED, ncd), thisroutine )

!  ------ Definition of the variables

      if (grads_nc_tweaks) then

!    ---- x

         call check( nf90_inq_dimid(ncid(n), 'x', nc1d), thisroutine )

         call check( nf90_def_var(ncid(n), 'x', NF90_FLOAT, nc1d, ncv), &
                     thisroutine )

         buffer = 'm'
         call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                     thisroutine )
         buffer = 'x'
         call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                     thisroutine )
         buffer = 'Dummy x-coordinate for one point'
         call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                     thisroutine )
         call check( nf90_put_att(ncid(n), ncv, 'axis', 'x'), thisroutine )

!    ---- y

         call check( nf90_inq_dimid(ncid(n), 'y', nc1d), thisroutine )

         call check( nf90_def_var(ncid(n), 'y', NF90_FLOAT, nc1d, ncv), &
                     thisroutine )

         buffer = 'm'
         call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                     thisroutine )
         buffer = 'y'
         call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                     thisroutine )
         buffer = 'Dummy y-coordinate for one point'
         call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                     thisroutine )
         call check( nf90_put_att(ncid(n), ncv, 'axis', 'y'), thisroutine )

      end if

!    ---- year2sec

      call check( nf90_def_var(ncid(n), 'year2sec', NF90_DOUBLE, ncv), &
               thisroutine )

      buffer = 's a-1'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'seconds_per_year'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = '1 year (1 a) in seconds'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- Time

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 't', NF90_DOUBLE, nc1d, ncv), &
                  thisroutine )

      buffer = 'a'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'time'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Time'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )
      call check( nf90_put_att(ncid(n), ncv, 'axis', 't'), thisroutine )

      if (grads_nc_tweaks) then

!    ---- Time offset

         call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

         call check( nf90_def_var(ncid(n), 't_add_offset', NF90_DOUBLE, nc1d, ncv), &
                     thisroutine )

         buffer = 'a'
         call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                     thisroutine )
         buffer = 'time_add_offset'
         call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                     thisroutine )
         buffer = 'Time offset'
         call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                     thisroutine )

         time_add_offset_val = min(time_val, 0.0_dp)

      end if

      if ((forcing_flag == 1).or.(forcing_flag == 3)) then

!    ---- delta_ts

         call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

         call check( nf90_def_var(ncid(n), 'delta_ts', NF90_FLOAT, nc1d, ncv), &
                     thisroutine )

         buffer = 'degC'
         call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                     thisroutine )
         buffer = 'surface_temperature_anomaly'
         call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                     thisroutine )
         buffer = 'Surface temperature anomaly'
         call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                     thisroutine )

      else if (forcing_flag == 2) then

!    ---- glac_index

         call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

         call check( nf90_def_var(ncid(n), 'glac_index', NF90_FLOAT, nc1d, ncv), &
                     thisroutine )

         buffer = '1'
         call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                     thisroutine )
         buffer = 'glacial_index'
         call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                     thisroutine )
         buffer = 'Glacial index'
         call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                     thisroutine )

      end if

!    ---- z_sl_mean

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'z_sl_mean', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'global_average_sea_level_change'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Mean sea level'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- V_tot

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'V_tot', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm3'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'land_ice_volume'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Ice volume'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- V_grounded

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'V_grounded', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm3'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'grounded_land_ice_volume'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Volume of grounded ice'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- V_floating

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'V_floating', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm3'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'floating_ice_shelf_volume'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Volume of floating ice'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- A_tot

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'A_tot', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm2'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'land_ice_area'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Ice area'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- A_grounded

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'A_grounded', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm2'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'grounded_land_ice_area'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Area covered by grounded ice'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- A_floating

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'A_floating', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm2'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'floating_ice_shelf_area'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Area covered by floating ice'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- V_af

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'V_af', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm3'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'land_ice_volume_not_displacing_sea_water'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Ice volume above flotation'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- V_sle

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'V_sle', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm SLE'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'land_ice_volume_sle'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Ice volume in SLE'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- V_temp

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'V_temp', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm3'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'temperate_land_ice_volume'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Volume of temperate ice'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- A_temp

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'A_temp', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm2'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'temperate_land_ice_area'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Area covered by temperate ice'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- Q_s

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'Q_s', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm3 ice equiv. a-1'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'tendency_of_land_ice_volume_due_to_surface_mass_balance'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Total surface mass balance'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- precip_tot

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'precip_tot', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm3 ice equiv. a-1'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'precipitation_rate'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Precipitation rate'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- runoff_tot

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'runoff_tot', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm3 ice equiv. a-1'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'runoff_rate'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Runoff rate'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- bmb_tot

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'bmb_tot', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm3 ice equiv. a-1'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'tendency_of_land_ice_volume_due_to_basal_mass_balance'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Total basal mass balance'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

#if (MARGIN==3)

!    ---- bmb_gr_tot

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'bmb_gr_tot', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm3 ice equiv. a-1'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'tendency_of_land_ice_volume_due_to_basal_mass_balance_for_grounded_ice'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Total basal mass balance for grounded ice'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- bmb_fl_tot

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'bmb_fl_tot', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm3 ice equiv. a-1'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'tendency_of_land_ice_volume_due_to_basal_mass_balance_for_floating_ice'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Total basal mass balance for floating ice'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

#endif

!    ---- Q_b

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'Q_b', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm3 ice equiv. a-1'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'basal_melting_rate'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Basal melting rate'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- Q_temp

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'Q_temp', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm3 ice equiv. a-1'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'temperate_layer_drainage_rate'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Drainage rate from the temperate ice layer'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- calv_tot

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'calv_tot', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm3 ice equiv. a-1'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'tendency_of_land_ice_volume_due_to_calving'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Total calving rate'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

#if (DISC>0)

!    ---- disc_lsc

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'disc_lsc', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm3 ice equiv. a-1'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'large_scale_ice_lost_into_the_ocean'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Large scale ice lost into the ocean'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- disc_ssc

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'disc_ssc', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm3 ice equiv. a-1'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'small_scale_ice_lost_into_the_ocean'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Small scale ice lost into the ocean'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- dT_glann

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'dT_glann', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'degC'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'global_annual_temperature_anomaly'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Global annual temperature anomaly'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- dT_sub

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'dT_sub', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'degC'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'subsurface_ocean_temperature_anomaly'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Subsurface ocean temperature anomaly'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

#endif

!    ---- dV_dt

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'dV_dt', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm3 a-1'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'tendency_of_land_ice_volume'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Rate of ice volume change'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

      if (n==0) then

!    ---- mb_resid

         call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

         call check( nf90_def_var(ncid(n), 'mb_resid', NF90_FLOAT, nc1d, ncv), &
                     thisroutine )

         buffer = 'm3 ice equiv. a-1'
         call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                     thisroutine )
         buffer = 'total_mass_balance_residual'
         call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                     thisroutine )
         buffer = 'Residual of the total mass balance'
         call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                     thisroutine )

!    ---- mb_mis

         call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

         call check( nf90_def_var(ncid(n), 'mb_mis', NF90_FLOAT, nc1d, ncv), &
                     thisroutine )

         buffer = 'm3 ice equiv. a-1'
         call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                     thisroutine )
         buffer = 'total_mass_balance_misacconted'
         call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                     thisroutine )
         buffer = 'Total mass balance misacconted'
         call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                     thisroutine )

      end if

!    ---- mbp

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'mbp', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = '1'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'mass_balance_partition'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Mass balance partition'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- H_max

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'H_max', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'maximum_ice_thickness'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Maximum ice thickness'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- H_t_max

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'H_t_max', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'maximum_thickness_of_temperate_ice'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Maximum thickness of temperate ice'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- zs_max

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'zs_max', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'maximum_surface_elevation'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Maximum surface elevation'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- vs_max

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'vs_max', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm a-1'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'maximum_surface_speed'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Maximum surface speed'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- Tbh_max

      call check( nf90_inq_dimid(ncid(n), 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid(n), 'Tbh_max', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'degC'
      call check( nf90_put_att(ncid(n), ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'maximum_basal_temperature_relative_to_pmp'
      call check( nf90_put_att(ncid(n), ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Maximum basal temperature relative to pmp'
      call check( nf90_put_att(ncid(n), ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- End of the definitions

      call check( nf90_enddef(ncid(n)), thisroutine )

   end if

!  ------ Writing of data on NetCDF file

   if (firstcall%output2) then

      if (grads_nc_tweaks) then

         nc1cor = [ 1 ]

         call check( nf90_inq_varid(ncid(n), 'x', ncv), thisroutine )
         call check( nf90_put_var(ncid(n), ncv, 0.0_sp, &
                                             start=nc1cor), thisroutine )

         call check( nf90_inq_varid(ncid(n), 'y', ncv), thisroutine )
         call check( nf90_put_var(ncid(n), ncv, 0.0_sp, &
                                             start=nc1cor), thisroutine )

      end if

      call check( nf90_inq_varid(ncid(n), 'year2sec', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, year2sec), thisroutine )

      call check( nf90_sync(ncid(n)), thisroutine )

   end if

   if (.not.flag_compute_flux_vars_only) then

      nc1cor(1) = counter
      ! nc1cnt(1) = 1   ! (not needed, and causes troubles for whatever reason)

      call check( nf90_inq_varid(ncid(n), 't', ncv), thisroutine )

      if (.not.grads_nc_tweaks) then
         call check( nf90_put_var(ncid(n), ncv, time_val, &
                                  start=nc1cor), thisroutine )
      else
         call check( nf90_put_var(ncid(n), ncv, time_val-time_add_offset_val, &
                                  start=nc1cor), thisroutine )
                               ! this makes sure that all times are >=0
                               ! (GrADS doesn't like negative numbers)
         call check( nf90_inq_varid(ncid(n), 't_add_offset', ncv), thisroutine )
         call check( nf90_put_var(ncid(n), ncv, time_add_offset_val, &
                                  start=nc1cor), thisroutine )
      end if

      if ((forcing_flag == 1).or.(forcing_flag == 3)) then

         call check( nf90_inq_varid(ncid(n), 'delta_ts', ncv), thisroutine )
         call check( nf90_put_var(ncid(n), ncv, real(delta_ts,sp), &
                                  start=nc1cor), thisroutine )

      else if (forcing_flag == 2) then

         call check( nf90_inq_varid(ncid(n), 'glac_index', ncv), thisroutine )
         call check( nf90_put_var(ncid(n), ncv, real(glac_index,sp), &
                                  start=nc1cor), thisroutine )

      end if

      call check( nf90_inq_varid(ncid(n), 'z_sl_mean', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(z_sl_mean,sp), &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid(n), 'V_tot', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(V_tot,sp), &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid(n), 'V_grounded', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(V_grounded,sp), &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid(n), 'V_floating', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(V_floating,sp), &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid(n), 'A_tot', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(A_tot,sp), &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid(n), 'A_grounded', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(A_grounded,sp), &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid(n), 'A_floating', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(A_floating,sp), &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid(n), 'V_af', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(V_af,sp), &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid(n), 'V_sle', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(V_sle,sp), &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid(n), 'V_temp', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(V_temp,sp), &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid(n), 'A_temp', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(A_temp,sp), &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid(n), 'Q_s', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(Q_s_flx,sp), &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid(n), 'precip_tot', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(precip_tot_flx,sp), &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid(n), 'runoff_tot', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(runoff_tot_flx,sp), &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid(n), 'bmb_tot', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(bmb_tot_flx,sp), &
                               start=nc1cor), thisroutine )

#if (MARGIN==3)

      call check( nf90_inq_varid(ncid(n), 'bmb_gr_tot', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(bmb_gr_tot_flx,sp), &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid(n), 'bmb_fl_tot', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(bmb_fl_tot_flx,sp), &
                               start=nc1cor), thisroutine )
#endif

      call check( nf90_inq_varid(ncid(n), 'Q_b', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(Q_b_flx,sp), &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid(n), 'Q_temp', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(Q_temp_flx,sp), &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid(n), 'calv_tot', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(calv_tot_flx,sp), &
                               start=nc1cor), thisroutine )

#if (DISC>0)
      call check( nf90_inq_varid(ncid(n), 'disc_lsc', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(disc_lsc_flx,sp), &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid(n), 'disc_ssc', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(disc_ssc_flx,sp), &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid(n), 'dT_glann', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(dT_glann,sp), &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid(n), 'dT_sub', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(dT_sub,sp), &
                               start=nc1cor), thisroutine )
#endif

      call check( nf90_inq_varid(ncid(n), 'dV_dt', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(dV_dt_flx,sp), &
                               start=nc1cor), thisroutine )

      if (n==0) then

         call check( nf90_inq_varid(ncid(n), 'mb_resid', ncv), thisroutine )
         call check( nf90_put_var(ncid(n), ncv, real(mb_resid_flx,sp), &
                                  start=nc1cor), thisroutine )

         call check( nf90_inq_varid(ncid(n), 'mb_mis', ncv), thisroutine )
         call check( nf90_put_var(ncid(n), ncv, real(mb_mis_flx,sp), &
                                  start=nc1cor), thisroutine )

      end if

      call check( nf90_inq_varid(ncid(n), 'mbp', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(mbp_flx,sp), &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid(n), 'H_max', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(H_max,sp), &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid(n), 'H_t_max', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(H_t_max,sp), &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid(n), 'zs_max', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(zs_max,sp), &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid(n), 'vs_max', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(vs_max,sp), &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid(n), 'Tbh_max', ncv), thisroutine )
      call check( nf90_put_var(ncid(n), ncv, real(Tbh_max,sp), &
                               start=nc1cor), thisroutine )

!  ------ Syncing NetCDF file (every 100th time)

      n_sync = 100

      if ( mod((counter-1), n_sync) == 0 ) &
         call check( nf90_sync(ncid(n)),  thisroutine )

   end if   ! (.not.flag_compute_flux_vars_only)

!-------- End loop over regions --------

end do   ! n=0, maxval(mask_region)

if (firstcall%output2) firstcall%output2 = .false.

end subroutine output2

!-------------------------------------------------------------------------------
!> Computation of the scalar output variables.
!-------------------------------------------------------------------------------
subroutine scalar_variables(time, &
                            H_cold, H_temp, &
                            time_val, &
                            V_tot, V_grounded, V_floating, &
                            A_tot, A_grounded, A_floating, &
                            V_af, V_sle, V_temp, A_temp, &
                            H_max, H_t_max, zs_max, vs_max, Tbh_max, &
                            dV_dt, Q_s, precip_tot, runoff_tot, &
                            Q_b, Q_temp, bmb_tot, bmb_gr_tot, bmb_fl_tot, &
                            calv_tot, disc_lsc, disc_ssc, &
                            mbp, mb_resid, mb_mis, &
                            opt_flag_region)

implicit none

real(dp), intent(in) :: time
real(dp), dimension(0:JMAX,0:IMAX), intent(in) :: H_cold, H_temp
logical, dimension(0:JMAX,0:IMAX), optional, intent(in) :: opt_flag_region

real(dp), intent(out) :: time_val, &
                         V_tot, V_grounded, V_floating, &
                         A_tot, A_grounded, A_floating, &
                         V_af, V_sle, V_temp, A_temp, &
                         H_max, H_t_max, zs_max, vs_max, Tbh_max, &
                         dV_dt, Q_s, precip_tot, runoff_tot, &
                         Q_b, Q_temp, bmb_tot, bmb_gr_tot, bmb_fl_tot, &
                         calv_tot, mbp, mb_resid, &
                         mb_mis, disc_lsc, disc_ssc

integer(i4b) :: i, j

real(dp) :: V_gr_redu, A_surf, rhosw_rho_ratio
real(dp) :: vs_help, Tbh_help
real(dp) :: MB, LMT, GIMB, SIMB, LMH, OMH, PAT, PAH, LQH

logical, dimension(0:JMAX,0:IMAX) :: flag_region

if (present(opt_flag_region)) then
   flag_region = opt_flag_region
else
   flag_region = .true.
end if

!-------- Computing scalar variables --------

#if (defined(NMARS) || defined(SMARS))   /* Martian ice sheet */
rhosw_rho_ratio = 0.0_dp   ! dummy value
#else   /* terrestrial ice sheet */
rhosw_rho_ratio = RHO_SW/RHO
#endif

V_grounded = 0.0_dp
V_gr_redu  = 0.0_dp
V_floating = 0.0_dp
A_grounded = 0.0_dp
A_floating = 0.0_dp
V_temp     = 0.0_dp
A_temp     = 0.0_dp
H_max      = 0.0_dp
H_t_max    = 0.0_dp
zs_max     = no_value_neg_2
vs_max     = 0.0_dp
Tbh_max    = no_value_neg_2

do i=0, IMAX
do j=0, JMAX

   if (flag_inner_point(j,i).and.flag_region(j,i)) then

      if (mask(j,i)==0) then   ! grounded ice

         if (zs(j,i)       > zs_max)  zs_max  = zs(j,i)
         if (H(j,i)        > H_max  ) H_max   = H(j,i)
         if (H_temp(j,i)   > H_t_max) H_t_max = H_temp(j,i)

         V_grounded = V_grounded + H(j,i)     *cell_area(j,i)
         V_temp     = V_temp     + H_temp(j,i)*cell_area(j,i)

#if (!defined(NMARS) && !defined(SMARS))   /* terrestrial ice sheet */
         V_gr_redu = V_gr_redu &
                     + rhosw_rho_ratio &
                          *max((z_sl(j,i)-zl(j,i)),0.0_dp)*cell_area(j,i)
#endif

         A_grounded = A_grounded + cell_area(j,i)

         if (n_cts(j,i) /= -1) A_temp = A_temp + cell_area(j,i)

         vs_help = sqrt(0.25_dp &
                          * ( (vx_c(KCMAX,j,i)+vx_c(KCMAX,j,i-1))**2 &
                             +(vy_c(KCMAX,j,i)+vy_c(KCMAX,j-1,i))**2 ) )
         if (vs_help > vs_max) vs_max = vs_help

         if (n_cts(j,i) >= 0) then   ! temperate base
            Tbh_max = 0.0_dp
         else   ! cold base
            Tbh_help = min((temp_c(0,j,i)-temp_c_m(0,j,i)), 0.0_dp)
            if (Tbh_help > Tbh_max) Tbh_max = Tbh_help
         end if

      else if (mask(j,i)==3) then   ! floating ice
                                    ! (basal temperature assumed to be below
                                    ! the pressure melting point for pure ice)

         if (zs(j,i)    > zs_max) zs_max = zs(j,i)
         if (H(j,i)     > H_max)  H_max  = H(j,i)

         V_floating = V_floating + H(j,i)*cell_area(j,i)
         A_floating = A_floating + cell_area(j,i)

         vs_help = sqrt(0.25_dp &
                          * ( (vx_c(KCMAX,j,i)+vx_c(KCMAX,j,i-1))**2 &
                             +(vy_c(KCMAX,j,i)+vy_c(KCMAX,j-1,i))**2 ) )
         if (vs_help > vs_max) vs_max = vs_help

         Tbh_help = min((temp_c(0,j,i)-temp_c_m(0,j,i)), 0.0_dp)
         if (Tbh_help > Tbh_max) Tbh_max = Tbh_help

      end if

   end if

end do
end do

!  ------ Unit conversion

#if (defined(NMARS) || defined(SMARS))   /* Martian ice sheet */

A_surf = 1.44371391e+14_dp   ! surface area of Mars, in m2
            ! Source: https://solarsystem.nasa.gov/planets/mars/by-the-numbers/
            !         (accessed on 2019-11-23)

V_af   = V_grounded
V_sle  = V_af*(RHO_I/RHO_W)*(1.0_dp-FRAC_DUST)/A_surf
                            ! m3 (ice+dust) equiv./m2 -> m water equiv.

#elif (defined(SHMARS))   /* Martian southern-hemisphere ice sheet */

A_surf = 1.44371391e+14_dp   ! surface area of Mars, in m2
            ! Source: https://solarsystem.nasa.gov/planets/mars/by-the-numbers/
            !         (accessed on 2019-11-23)

V_af   = V_grounded
V_sle  = V_af*(RHO/RHO_W)/A_surf   ! m3 ice equiv./m2 -> m water equiv.

#else   /* terrestrial ice sheet */

A_surf = 3.61132e+14_dp   ! global ocean area, in m2

V_af   = V_grounded - V_gr_redu
V_sle  = V_af*(RHO/RHO_W)/A_surf   ! m3 ice equiv./m2 -> m water equiv.

#endif

#if (!defined(OUT_TIMES) || OUT_TIMES==1)
time_val = time *sec2year   ! s -> a
#elif (OUT_TIMES==2)
time_val = (time+year_zero) *sec2year   ! s -> a
#else
errormsg = ' >>> scalar_variables: OUT_TIMES must be either 1 or 2!'
call error(errormsg)
#endif

V_tot = V_grounded + V_floating   ! in m3
A_tot = A_grounded + A_floating   ! in m2

vs_max = vs_max *year2sec   ! m/s -> m/a

!-------- Computing further scalar variables (mass balance components) --------

Q_s    = 0.0_dp
Q_b    = 0.0_dp
Q_temp = 0.0_dp

do i=0, IMAX
do j=0, JMAX

   if (flag_inner_point(j,i).and.flag_region(j,i)) then

      Q_s = Q_s + as_perp_apl(j,i) * cell_area(j,i)

      if (     (mask(j,i)==0).or.(mask_old(j,i)==0) &
           .or.(mask(j,i)==3).or.(mask_old(j,i)==3) &
         ) &   ! grounded or floating ice before or after the time step
         Q_b = Q_b + Q_bm(j,i) * cell_area(j,i)
                                 !!% Also *_apl required (or delete?)

      if ( (mask(j,i)==0).or.(mask_old(j,i)==0) &
         ) &   ! grounded ice before or after the time step
         Q_temp = Q_temp + Q_tld(j,i) * cell_area(j,i)
                                        !!% Also *_apl required (or delete?)

   end if

end do
end do

Q_s    = Q_s    *year2sec   ! m3/s ice equiv. -> m3/a ice equiv.
Q_b    = Q_b    *year2sec   ! m3/s ice equiv. -> m3/a ice equiv.
Q_temp = Q_temp *year2sec   ! m3/s ice equiv. -> m3/a ice equiv.

dV_dt      = 0.0_dp
precip_tot = 0.0_dp

do i=0, IMAX
do j=0, JMAX

   if (flag_inner_point(j,i).and.flag_region(j,i)) then

      if (     (mask(j,i)==0).or.(mask_old(j,i)==0) &
           .or.(mask(j,i)==3).or.(mask_old(j,i)==3) ) then
               ! grounded or floating ice before or after the time step
         dV_dt = dV_dt + (dzs_dtau(j,i)-dzb_dtau(j,i))*cell_area(j,i)
                 !!% change to more direct, V_tot-based computation?
      end if

      precip_tot = precip_tot + accum_apl(j,i)*cell_area(j,i)

   end if

end do
end do

precip_tot = precip_tot *year2sec   ! m3/s ice equiv. -> m3/a ice equiv.
dV_dt      = dV_dt      *year2sec   ! m3/s ice equiv. -> m3/a ice equiv.

! MB:   total mass balance as computed in subroutine apply_smb
! LMT:  total lost on land at the top of the ice sheet including ice shelf
! GIMB: total basal lost grounded ice
! SIMB: total basal lost ice shelf
! LMH:  total hidden lost through runoff on land
! OMH:  total hidden lost through flow into the ocean
! PAT:  total lost through ice discharge/calving parameterizations
! PAH:  total hidden lost through ice discharge/calving parameterization
! LQH:  total hidden lost through melt at the base on land
! mb_mis: misaccounted

MB     = 0.0_dp
LMT    = 0.0_dp
GIMB   = 0.0_dp
SIMB   = 0.0_dp
LMH    = 0.0_dp
LQH    = 0.0_dp
OMH    = 0.0_dp
PAT    = 0.0_dp
PAH    = 0.0_dp
mb_mis = 0.0_dp

do i=0, IMAX
do j=0, JMAX

   if (flag_inner_point(j,i).and.flag_region(j,i)) then

      if ( mask_ablation_type(j,i) /= 0 ) then
                     ! glaciated land and ocean (including hidden melt points)

         ! Quantify what types of melt occurred
         select case ( mask_ablation_type(j,i) )
            case( 3 )
               LMT  = LMT + runoff_apl(j,i)  * cell_area(j,i)
               PAT  = PAT + calving_apl(j,i) * cell_area(j,i)
                               ! could cause problems for Greenland
               SIMB = SIMB + Q_b_apl(j,i)    * cell_area(j,i)
            case( 1 )
               LMT  = LMT + runoff_apl(j,i)  * cell_area(j,i)
               PAT  = PAT + calving_apl(j,i) * cell_area(j,i)  ! ok
               GIMB = GIMB + Q_b_apl(j,i)    * cell_area(j,i)
            case( 9 )
               mb_mis = mb_mis + mb_source_apl(j,i) * cell_area(j,i)
            case( -1 )
               LMH = LMH + runoff_apl(j,i)  * cell_area(j,i)
               PAH = PAH + calving_apl(j,i) * cell_area(j,i)
               LQH = LQH + Q_b_apl(j,i)     * cell_area(j,i)
            case( -2 )
               OMH = OMH + calving_apl(j,i) * cell_area(j,i) ! only one contribution
         end select

      end if

      ! Actual ice mass balance (from top melt, bottom melt and calving)
      MB = MB + mb_source_apl(j,i)*cell_area(j,i)

   end if

end do
end do

! Runoff on land (excluding basal melt)
runoff_tot = LMT + LMH

! Ice discharge (excluding basal melt)
calv_tot   = OMH + PAT + PAH

! Ice discharge from ice flow, large scale
disc_lsc   = OMH

! Ice discharge from parameterization, small scale
disc_ssc   = PAT + PAH

! Basal mass balance
bmb_tot    = -GIMB-SIMB-LQH

! Grounded ice
bmb_gr_tot = -GIMB-LQH

! Ice shelf
bmb_fl_tot = -SIMB   ! hidden is counted as large scale calving

MB         = MB         * year2sec
                        ! m3/s ice equiv. -> m3/a ice equiv.
mb_mis     = mb_mis     * year2sec
                        ! m3/s ice equiv. -> m3/a ice equiv.
disc_lsc   = disc_lsc   * year2sec
                        ! m3/s ice equiv. -> m3/a ice equiv.
disc_ssc   = disc_ssc   * year2sec
                        ! m3/s ice equiv. -> m3/a ice equiv.
bmb_tot    = bmb_tot    * year2sec
                        ! m3/s ice equiv. -> m3/a ice equiv.
bmb_fl_tot = bmb_fl_tot * year2sec
                        ! m3/s ice equiv. -> m3/a ice equiv.
bmb_gr_tot = bmb_gr_tot * year2sec
                        ! m3/s ice equiv. -> m3/a ice equiv.
runoff_tot = runoff_tot * year2sec
                        ! m3/s ice equiv. -> m3/a ice equiv.
calv_tot   = calv_tot   * year2sec
                        ! m3/s ice equiv. -> m3/a ice equiv.

if (precip_tot /= 0.0_dp) then
   mbp = calv_tot/precip_tot
else
   mbp = 0.0_dp
end if

mb_resid = Q_s + bmb_tot - calv_tot - dV_dt
!!% (previously) mb_resid = MB - dV_dt

end subroutine scalar_variables

!-------------------------------------------------------------------------------
!> Writing of time-series data of the specified sites (i.e., ice cores) on file
!! in ASCII and NetCDF format.
!-------------------------------------------------------------------------------
subroutine output4(time, dxi, deta, delta_ts, glac_index)

  use netcdf
  use nc_check_m

implicit none

real(dp), intent(in) :: time, dxi, deta, delta_ts, glac_index

integer(i4b)                        :: i, j, n
integer(i4b)                        :: ios
real(dp), dimension(:), allocatable :: r_n_site
real(dp)                            :: time_val
real(dp), dimension(0:JMAX,0:IMAX)  :: field
real(dp), dimension(:), allocatable :: H_site, temp_b_site, &
                                       vx_b_site, vy_b_site, vh_b_site, &
                                       vx_s_site, vy_s_site, vh_s_site, &
                                       Rx_b_site, Ry_b_site, R_b_site, &
                                       bmb_site
real(dp), dimension(:), allocatable :: arg1

integer(i4b)       :: cmode
integer(i4b)       :: n_deflate_level
logical            :: flag_shuffle
integer(i4b), save :: ncid
integer(i4b)       :: ncd, ncv, nc1d, nc2d(2)
integer(i4b)       :: nc1cor(1), nc1cnt(1), nc2cor(2), nc2cnt(2)
integer(i4b)       :: n_sync
real(dp), save     :: time_add_offset_val
character(len= 16) :: ch_date, ch_time, ch_zone
character(len=256) :: filename, filename_with_path, buffer
logical            :: flag_in_domain
logical, save      :: grads_nc_tweaks

integer(i4b), save :: counter = 0

character(len=64), parameter :: thisroutine = 'output4'

counter = counter + 1

if (n_site >= 1) then

   allocate(r_n_site(n_site), H_site(n_site), temp_b_site(n_site), &
            vx_b_site(n_site), vy_b_site(n_site), vh_b_site(n_site), &
            vx_s_site(n_site), vy_s_site(n_site), vh_s_site(n_site), &
            Rx_b_site(n_site), Ry_b_site(n_site), R_b_site(n_site), &
            bmb_site(n_site))
   allocate(arg1(n_site))

!-------- Determination of site data --------

   do n=1, n_site

!  ------ Site number

      r_n_site(n) = real(n,dp)

!  ------ Ice thickness

      field = H_c + H_t
      call borehole(field, x_site(n), y_site(n), dxi, deta, 'grid', &
                    H_site(n), flag_in_domain)

      if ( (firstcall%output4).and.(.not.flag_in_domain) ) then
         write(6, fmt='(a,i0,a)') &
            ' >>> output4: WARNING: Borehole position ', n, ' out of domain!'
      end if

!  ------ Basal velocity

      do i=0, IMAX
      do j=0, JMAX
         field(j,i) = vx_t(0,j,i)
      end do
      end do

      call borehole(field, x_site(n), y_site(n), dxi, deta, 'sg_x', &
                    vx_b_site(n), flag_in_domain)

      do i=0, IMAX
      do j=0, JMAX
         field(j,i) = vy_t(0,j,i)
      end do
      end do

      call borehole(field, x_site(n), y_site(n), dxi, deta, 'sg_y', &
                    vy_b_site(n), flag_in_domain)

      vh_b_site(n) = sqrt(vx_b_site(n)**2+vy_b_site(n)**2)

!  ------ Surface velocity

      do i=0, IMAX
      do j=0, JMAX
         field(j,i) = vx_c(KCMAX,j,i)
      end do
      end do

      call borehole(field, x_site(n), y_site(n), dxi, deta, 'sg_x', &
                    vx_s_site(n), flag_in_domain)

      do i=0, IMAX
      do j=0, JMAX
         field(j,i) = vy_c(KCMAX,j,i)
      end do
      end do

      call borehole(field, x_site(n), y_site(n), dxi, deta, 'sg_y', &
                    vy_s_site(n), flag_in_domain)

      vh_s_site(n) = sqrt(vx_s_site(n)**2+vy_s_site(n)**2)

!  ------ Basal temperature

      do i=0, IMAX
      do j=0, JMAX
         field(j,i) = temp_r(KRMAX,j,i)
      end do
      end do

      call borehole(field, x_site(n), y_site(n), dxi, deta, 'grid', &
                    temp_b_site(n), flag_in_domain)

!  ------ Basal frictional heating

      do i=0, IMAX
      do j=0, JMAX
         field(j,i) = vx_t(0,j,i)*txz_t(0,j,i)
      end do
      end do

      call borehole(field, x_site(n), y_site(n), dxi, deta, 'sg_x', &
                    Rx_b_site(n), flag_in_domain)

      do i=0, IMAX
      do j=0, JMAX
         field(j,i) = vy_t(0,j,i)*tyz_t(0,j,i)
      end do
      end do

      call borehole(field, x_site(n), y_site(n), dxi, deta, 'sg_y', &
                    Ry_b_site(n), flag_in_domain)

      R_b_site(n) = Rx_b_site(n) + Ry_b_site(n)

!  ------ Basal mass balance

      field = -(Q_bm+Q_tld)   ! positive for supply, negative for loss
      call borehole(field, x_site(n), y_site(n), dxi, deta, 'grid', &
                    bmb_site(n), flag_in_domain)

   end do

!  ------ Conversion

#if (!defined(OUT_TIMES) || OUT_TIMES==1)
   time_val = time *sec2year   ! s -> a
#elif (OUT_TIMES==2)
   time_val = (time+year_zero) *sec2year   ! s -> a
#else
   errormsg = ' >>> output4: OUT_TIMES must be either 1 or 2!'
   call error(errormsg)
#endif

   vh_b_site = vh_b_site *year2sec   ! m/s -> m/a
   vh_s_site = vh_s_site *year2sec   ! m/s -> m/a

   bmb_site  = bmb_site  *year2sec   ! m ice equiv./s -> m ice equiv./a

!-------- Writing of data on file --------

   if ((forcing_flag == 1).or.(forcing_flag == 3)) then
      write(unit=14, fmt='(1pe13.6,2(1pe13.4))') time_val, delta_ts, z_sl_mean
   else if (forcing_flag == 2) then
      write(unit=14, fmt='(1pe13.6,2(1pe13.4))') time_val, glac_index, z_sl_mean
   end if

   n=1
   write(unit=14, fmt='(13x,1pe13.4)', advance='no') H_site(n)
   do n=2, n_site-1
      write(unit=14, fmt='(1pe13.4)', advance='no') H_site(n)
   end do
   n=n_site
   write(unit=14, fmt='(1pe13.4)') H_site(n)

   n=1
   write(unit=14, fmt='(13x,1pe13.4)', advance='no') vh_s_site(n)
   do n=2, n_site-1
      write(unit=14, fmt='(1pe13.4)', advance='no') vh_s_site(n)
   end do
   n=n_site
   write(unit=14, fmt='(1pe13.4)') vh_s_site(n)

   n=1
   write(unit=14, fmt='(13x,1pe13.4)', advance='no') temp_b_site(n)
   do n=2, n_site-1
      write(unit=14, fmt='(1pe13.4)', advance='no') temp_b_site(n)
   end do
   n=n_site
   write(unit=14, fmt='(1pe13.4,/)') temp_b_site(n)

!-------- Extended time-series file in NetCDF format --------

   if (firstcall%output4) then

!  ------ Open NetCDF file

      filename           = trim(run_name)//'_site.nc'
      filename_with_path = trim(OUT_PATH)//'/'//trim(filename)

      !%% call set_cmode(cmode, n_deflate_level, flag_shuffle)

      cmode           = NF90_NOCLOBBER
      n_deflate_level = 0
      flag_shuffle    = .false.

      ios = nf90_create(trim(filename_with_path), cmode, ncid)

      if (ios /= nf90_noerr) then
         errormsg = ' >>> output4: Error when opening the' &
                  //               end_of_line &
                  //'              NetCDF site time-series file!'
         call error(errormsg)
      end if

      ncid_site = ncid

!  ------ Global attributes

      buffer = 'Time-series output for the specified sites of simulation '// &
               trim(run_name)
      call check( nf90_put_att(ncid, NF90_GLOBAL, 'title', trim(buffer)), &
                  thisroutine )

      buffer = 'SICOPOLIS Version '//VERSION
      call check( nf90_put_att(ncid, NF90_GLOBAL, 'source', trim(buffer)), &
                  thisroutine )

      call date_and_time(ch_date, ch_time, ch_zone)
      buffer = ch_date(1:4)//'-'//ch_date(5:6)//'-'//ch_date(7:8)//' '// &
               ch_time(1:2)//':'//ch_time(3:4)//':'//ch_time(5:6)//' '// &
               ch_zone(1:3)//':'//ch_zone(4:5)//' - Data produced'
      call check( nf90_put_att(ncid, NF90_GLOBAL, 'history', trim(buffer)), &
                  thisroutine )

      call set_ch_website(buffer)
      call check( nf90_put_att(ncid, NF90_GLOBAL, 'references', trim(buffer)), &
                  thisroutine )

      call set_ch_institution(buffer)
      call check( nf90_put_att(ncid, NF90_GLOBAL, 'institution', &
                               trim(buffer)), &
                  thisroutine )

!  ------ Definition of the dimensions

      call set_grads_nc_tweaks(grads_nc_tweaks)

      if (grads_nc_tweaks) then
         call check( nf90_def_dim(ncid, 'x', 1, ncd), thisroutine )
         call check( nf90_def_dim(ncid, 'y', 1, ncd), thisroutine )
      end if

      call check( nf90_def_dim(ncid, 'n', n_site,  ncd), thisroutine )
      call check( nf90_def_dim(ncid, 't', NF90_UNLIMITED, ncd), thisroutine )

!  ------ Definition of the variables

      if (grads_nc_tweaks) then

!    ---- x

         call check( nf90_inq_dimid(ncid, 'x', nc1d), thisroutine )

         call check( nf90_def_var(ncid, 'x', NF90_FLOAT, nc1d, ncv), &
                     thisroutine )

         buffer = 'm'
         call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
                     thisroutine )
         buffer = 'x'
         call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
                     thisroutine )
         buffer = 'Dummy x-coordinate for one point'
         call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
                     thisroutine )
         call check( nf90_put_att(ncid, ncv, 'axis', 'x'), thisroutine )

!    ---- y

         call check( nf90_inq_dimid(ncid, 'y', nc1d), thisroutine )

         call check( nf90_def_var(ncid, 'y', NF90_FLOAT, nc1d, ncv), &
                     thisroutine )

         buffer = 'm'
         call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
                     thisroutine )
         buffer = 'y'
         call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
                     thisroutine )
         buffer = 'Dummy y-coordinate for one point'
         call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
                     thisroutine )
         call check( nf90_put_att(ncid, ncv, 'axis', 'y'), thisroutine )

      end if

!    ---- year2sec

      call check( nf90_def_var(ncid, 'year2sec', NF90_DOUBLE, ncv), &
               thisroutine )

      buffer = 's a-1'
      call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'seconds_per_year'
      call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = '1 year (1 a) in seconds'
      call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- Time

      call check( nf90_inq_dimid(ncid, 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid, 't', NF90_DOUBLE, nc1d, ncv), &
                  thisroutine )

      buffer = 'a'
      call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'time'
      call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Time'
      call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
                  thisroutine )
      call check( nf90_put_att(ncid, ncv, 'axis', 't'), thisroutine )

      if (grads_nc_tweaks) then

!    ---- Time offset

         call check( nf90_inq_dimid(ncid, 't', nc1d), thisroutine )

         call check( nf90_def_var(ncid, 't_add_offset', &
                                  NF90_DOUBLE, nc1d, ncv), &
                     thisroutine )

         buffer = 'a'
         call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
                     thisroutine )
         buffer = 'time_add_offset'
         call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
                     thisroutine )
         buffer = 'Time offset'
         call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
                     thisroutine )

         time_add_offset_val = min(time_val, 0.0_dp)

      end if

!    ---- Site number

      call check( nf90_inq_dimid(ncid, 'n', nc1d), thisroutine )

      call check( nf90_def_var(ncid, 'n', NF90_DOUBLE, nc1d, ncv), &
                  thisroutine )

      buffer = '1'
      call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'site_number'
      call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Site number'
      call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
                  thisroutine )
      buffer = trim(ch_site(1))
      do n=2, n_site; buffer = trim(buffer)//', '//trim(ch_site(n)); end do
      call check( nf90_put_att(ncid, ncv, 'site_names', trim(buffer)), &
                  thisroutine )
      call check( nf90_put_att(ncid, ncv, 'axis', 'n'), thisroutine )

      if ((forcing_flag == 1).or.(forcing_flag == 3)) then

!    ---- delta_ts

         call check( nf90_inq_dimid(ncid, 't', nc1d), thisroutine )

         call check( nf90_def_var(ncid, 'delta_ts', NF90_FLOAT, nc1d, ncv), &
                     thisroutine )

         buffer = 'degC'
         call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
                     thisroutine )
         buffer = 'surface_temperature_anomaly'
         call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
                     thisroutine )
         buffer = 'Surface temperature anomaly'
         call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
                     thisroutine )

      else if (forcing_flag == 2) then

!    ---- glac_index

         call check( nf90_inq_dimid(ncid, 't', nc1d), thisroutine )

         call check( nf90_def_var(ncid, 'glac_index', NF90_FLOAT, nc1d, ncv), &
                     thisroutine )

         buffer = '1'
         call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
                     thisroutine )
         buffer = 'glacial_index'
         call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
                     thisroutine )
         buffer = 'Glacial index'
         call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
                     thisroutine )

      end if

!    ---- z_sl_mean

      call check( nf90_inq_dimid(ncid, 't', nc1d), thisroutine )

      call check( nf90_def_var(ncid, 'z_sl_mean', NF90_FLOAT, nc1d, ncv), &
                  thisroutine )

      buffer = 'm'
      call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'global_average_sea_level_change'
      call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Mean sea level'
      call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- H_site

      call check( nf90_inq_dimid(ncid, 'n', nc2d(1)), thisroutine )
      call check( nf90_inq_dimid(ncid, 't', nc2d(2)), thisroutine )

      call check( nf90_def_var(ncid, 'H_site', NF90_FLOAT, nc2d, ncv), &
                  thisroutine )

      buffer = 'm'
      call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'land_ice_thickness'
      call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Ice thickness'
      call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- vh_b_site

      call check( nf90_inq_dimid(ncid, 'n', nc2d(1)), thisroutine )
      call check( nf90_inq_dimid(ncid, 't', nc2d(2)), thisroutine )

      call check( nf90_def_var(ncid, 'vh_b_site', NF90_FLOAT, nc2d, ncv), &
                  thisroutine )

      buffer = 'm a-1'
      call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'land_ice_basal_horizontal_velocity'
      call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Horizontal velocity at the ice base'
      call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- vh_s_site

      call check( nf90_inq_dimid(ncid, 'n', nc2d(1)), thisroutine )
      call check( nf90_inq_dimid(ncid, 't', nc2d(2)), thisroutine )

      call check( nf90_def_var(ncid, 'vh_s_site', NF90_FLOAT, nc2d, ncv), &
                  thisroutine )

      buffer = 'm a-1'
      call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'land_ice_surface_horizontal_velocity'
      call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Horizontal velocity at the ice surface'
      call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- temp_b_site

      call check( nf90_inq_dimid(ncid, 'n', nc2d(1)), thisroutine )
      call check( nf90_inq_dimid(ncid, 't', nc2d(2)), thisroutine )

      call check( nf90_def_var(ncid, 'temp_b_site', NF90_FLOAT, nc2d, ncv), &
                  thisroutine )

      buffer = 'degC'
      call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'basal_temperature'
      call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Temperature at the ice base'
      call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- R_b_site

      call check( nf90_inq_dimid(ncid, 'n', nc2d(1)), thisroutine )
      call check( nf90_inq_dimid(ncid, 't', nc2d(2)), thisroutine )

      call check( nf90_def_var(ncid, 'R_b_site', NF90_FLOAT, nc2d, ncv), &
                  thisroutine )

      buffer = 'W m-2'
      call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'basal_frictional_heating'
      call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Basal frictional heating'
      call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- bmb_site

      call check( nf90_inq_dimid(ncid, 'n', nc2d(1)), thisroutine )
      call check( nf90_inq_dimid(ncid, 't', nc2d(2)), thisroutine )

      call check( nf90_def_var(ncid, 'bmb_site', NF90_FLOAT, nc2d, ncv), &
                  thisroutine )

      buffer = 'm ice equiv. a-1'
      call check( nf90_put_att(ncid, ncv, 'units', trim(buffer)), &
                  thisroutine )
      buffer = 'land_ice_basal_mass_balance'
      call check( nf90_put_att(ncid, ncv, 'standard_name', trim(buffer)), &
                  thisroutine )
      buffer = 'Basal mass balance'
      call check( nf90_put_att(ncid, ncv, 'long_name', trim(buffer)), &
                  thisroutine )

!    ---- End of the definitions

      call check( nf90_enddef(ncid), thisroutine )

   end if

!  ------ Writing of data on NetCDF file

   if (firstcall%output4) then

      if (grads_nc_tweaks) then

         nc1cor = [ 1 ]

         call check( nf90_inq_varid(ncid, 'x', ncv), thisroutine )
         call check( nf90_put_var(ncid, ncv, 0.0_sp, &
                                             start=nc1cor), thisroutine )

         call check( nf90_inq_varid(ncid, 'y', ncv), thisroutine )
         call check( nf90_put_var(ncid, ncv, 0.0_sp, &
                                             start=nc1cor), thisroutine )

      end if

      nc1cor = [ 1 ]

      call check( nf90_inq_varid(ncid, 'n', ncv), thisroutine )
      call check( nf90_put_var(ncid, ncv, r_n_site, &
                               start=nc1cor), thisroutine )

      call check( nf90_inq_varid(ncid, 'year2sec', ncv), thisroutine )
      call check( nf90_put_var(ncid, ncv, year2sec), thisroutine )

      call check( nf90_sync(ncid), thisroutine )

   end if

   nc1cor(1) = counter
   ! nc1cnt(1) = 1   ! (not needed, and causes troubles for whatever reason)

   nc2cor(1) = 1
   nc2cor(2) = counter
   nc2cnt(1) = n_site
   nc2cnt(2) = 1

   call check( nf90_inq_varid(ncid, 't', ncv), thisroutine )

   if (.not.grads_nc_tweaks) then
      call check( nf90_put_var(ncid, ncv, time_val, &
                               start=nc1cor), thisroutine )
   else
      call check( nf90_put_var(ncid, ncv, time_val-time_add_offset_val, &
                               start=nc1cor), thisroutine )
                            ! this makes sure that all times are >=0
                            ! (GrADS doesn't like negative numbers)
      call check( nf90_inq_varid(ncid, 't_add_offset', ncv), thisroutine )
      call check( nf90_put_var(ncid, ncv, time_add_offset_val, &
                               start=nc1cor), thisroutine )
   end if

   if ((forcing_flag == 1).or.(forcing_flag == 3)) then

      call check( nf90_inq_varid(ncid, 'delta_ts', ncv), thisroutine )
      call check( nf90_put_var(ncid, ncv, delta_ts, &
                               start=nc1cor), thisroutine )

   else if (forcing_flag == 2) then

      call check( nf90_inq_varid(ncid, 'glac_index', ncv), thisroutine )
      call check( nf90_put_var(ncid, ncv, glac_index, &
                               start=nc1cor), thisroutine )

   end if

   call check( nf90_inq_varid(ncid, 'z_sl_mean', ncv), thisroutine )
   call check( nf90_put_var(ncid, ncv, z_sl_mean, &
                            start=nc1cor), thisroutine )

   call check( nf90_inq_varid(ncid, 'H_site', ncv), thisroutine )
   arg1 = H_site
   call check( nf90_put_var(ncid, ncv, arg1, &
                            start=nc2cor, count=nc2cnt), thisroutine )

   call check( nf90_inq_varid(ncid, 'vh_b_site', ncv), thisroutine )
   arg1 = vh_b_site
   call check( nf90_put_var(ncid, ncv, arg1, &
                            start=nc2cor, count=nc2cnt), thisroutine )

   call check( nf90_inq_varid(ncid, 'vh_s_site', ncv), thisroutine )
   arg1 = vh_s_site
   call check( nf90_put_var(ncid, ncv, arg1, &
                            start=nc2cor, count=nc2cnt), thisroutine )

   call check( nf90_inq_varid(ncid, 'temp_b_site', ncv), thisroutine )
   arg1 = temp_b_site
   call check( nf90_put_var(ncid, ncv, arg1, &
                            start=nc2cor, count=nc2cnt), thisroutine )

   call check( nf90_inq_varid(ncid, 'R_b_site', ncv), thisroutine )
   arg1 = R_b_site
   call check( nf90_put_var(ncid, ncv, arg1, &
                            start=nc2cor, count=nc2cnt), thisroutine )

   call check( nf90_inq_varid(ncid, 'bmb_site', ncv), thisroutine )
   arg1 = bmb_site
   call check( nf90_put_var(ncid, ncv, arg1, &
                            start=nc2cor, count=nc2cnt), thisroutine )

!  ------ Syncing NetCDF file (every 100th time)

   n_sync = 100

   if ( mod((counter-1), n_sync) == 0 ) &
      call check( nf90_sync(ncid),  thisroutine )

   deallocate(r_n_site, H_site, &
              vx_b_site, vy_b_site, vh_b_site, &
              vx_s_site, vy_s_site, vh_s_site, &
              temp_b_site, &
              Rx_b_site, Ry_b_site, R_b_site, &
              bmb_site)
   deallocate(arg1)

!%% else   ! (n_site == 0 -> do nothing)
!%%
!%%    continue

end if

if (firstcall%output4) firstcall%output4 = .false.

end subroutine output4

#if (defined(ASF) && WRITE_SER_FILE_STAKES==1) /* Austfonna */

!-------------------------------------------------------------------------------
!> Writing of time-series data for all defined surface points on file
!! in ASCII format. Modification of Tolly's output7 by Thorben Dunse.
!-------------------------------------------------------------------------------
subroutine output5(time, dxi, deta, delta_ts, glac_index)

implicit none

real(dp), intent(in) :: time, dxi, deta, delta_ts, glac_index

integer(i4b) :: n, k
real(dp) :: time_val
real(dp), dimension(0:JMAX,0:IMAX) :: field
real(dp), dimension(:), allocatable :: zl_surf, zs_surf, &
                  accum_surf, as_perp_surf, &
		  snowfall_surf, rainfall_surf, runoff_surf, &
		  vx_surf, vy_surf, vz_surf, &
                  vx_base, vy_base, vz_base, &
		  temp_base_pmp
logical :: flag_in_domain

allocate(zl_surf(n_surf), zs_surf(n_surf), &
        accum_surf(n_surf), &
	as_perp_surf(n_surf), snowfall_surf(n_surf), &
	rainfall_surf(n_surf), runoff_surf(n_surf), &
	vx_surf(n_surf), vy_surf(n_surf), vz_surf(n_surf), &
	vx_base(n_surf), vy_base(n_surf), vz_base(n_surf), &
 	temp_base_pmp(n_surf))

!-------- Determination of site data --------

do n=1, n_surf

!  ------ Bedrock elevation

   field = zl
   call borehole(field, x_surf(n), y_surf(n), dxi, deta, 'grid', &
                 zl_surf(n), flag_in_domain)

!  ------ Surface elevation

   field = zs
   call borehole(field, x_surf(n), y_surf(n), dxi, deta, 'grid', &
                 zs_surf(n), flag_in_domain)

!  ------ Accumulation

   field = accum
   call borehole(field, x_surf(n), y_surf(n), dxi, deta, 'grid', &
                 accum_surf(n), flag_in_domain)

!  ------ Surface mass balance

   field = as_perp
   call borehole(field, x_surf(n), y_surf(n), dxi, deta, 'grid', &
                 as_perp_surf(n), flag_in_domain)

!  ------ Snowfall

   field = snowfall
   call borehole(field, x_surf(n), y_surf(n), dxi, deta, 'grid', &
                 snowfall_surf(n), flag_in_domain)

!  ------ Rainfall

   field = rainfall
   call borehole(field, x_surf(n), y_surf(n), dxi, deta, 'grid', &
                 rainfall_surf(n), flag_in_domain)

!  ------ Runoff

   field = runoff
   call borehole(field, x_surf(n), y_surf(n), dxi, deta, 'grid', &
                 runoff_surf(n), flag_in_domain)

! ------ Surface velocities

   field = vx_s_g
   call borehole(field, x_surf(n), y_surf(n), dxi, deta, 'grid', &
                 vx_surf(n), flag_in_domain)

   field = vy_s_g
   call borehole(field, x_surf(n), y_surf(n), dxi, deta, 'grid', &
                 vy_surf(n), flag_in_domain)

   field = vz_s
   call borehole(field, x_surf(n), y_surf(n), dxi, deta, 'grid', &
                 vz_surf(n), flag_in_domain)

! ------ Basal velocities

   field = vx_b_g
   call borehole(field, x_surf(n), y_surf(n), dxi, deta, 'grid', &
                 vx_base(n), flag_in_domain)

   field = vy_b_g
   call borehole(field, x_surf(n), y_surf(n), dxi, deta, 'grid', &
                 vy_base(n), flag_in_domain)

   field = vz_b
   call borehole(field, x_surf(n), y_surf(n), dxi, deta, 'grid', &
                 vz_base(n), flag_in_domain)

! ------ Basal temperature relative to pressure melting point

   field = temph_b
   call borehole(field, x_surf(n), y_surf(n), dxi, deta, 'grid', &
                 temp_base_pmp(n), flag_in_domain)

end do

!  ------ Conversion

#if (!defined(OUT_TIMES) || OUT_TIMES==1)
time_val = time *sec2year   ! s -> a
#elif (OUT_TIMES==2)
time_val = (time+year_zero) *sec2year   ! s -> a
#else
errormsg = ' >>> output5: OUT_TIMES must be either 1 or 2!'
call error(errormsg)
#endif

do n=1, n_surf
   accum_surf(n)    = accum_surf(n)    *year2sec   ! m/s -> m/a
   as_perp_surf(n)  = as_perp_surf(n)  *year2sec   ! m/s -> m/a
   snowfall_surf(n) = snowfall_surf(n) *year2sec   ! m/s -> m/a
   rainfall_surf(n) = rainfall_surf(n) *year2sec   ! m/s -> m/a
   runoff_surf(n)   = runoff_surf(n)   *year2sec   ! m/s -> m/a
   vx_surf(n)       = vx_surf(n)       *year2sec   ! m/s -> m/a
   vy_surf(n)       = vy_surf(n)       *year2sec   ! m/s -> m/a
   vz_surf(n)       = vz_surf(n)       *year2sec   ! m/s -> m/a
   vx_base(n)       = vx_base(n)       *year2sec   ! m/s -> m/a
   vy_base(n)       = vy_base(n)       *year2sec   ! m/s -> m/a
   vz_base(n)       = vz_base(n)       *year2sec   ! m/s -> m/a
end do

!-------- Writing of data on file --------

if ((forcing_flag == 1).or.(forcing_flag == 3)) then
   write(41,'(1pe13.6,2(1pe13.4))') time_val, delta_ts, z_sl_mean
   write(42,'(1pe13.6,2(1pe13.4))') time_val, delta_ts, z_sl_mean
   write(43,'(1pe13.6,2(1pe13.4))') time_val, delta_ts, z_sl_mean
   write(44,'(1pe13.6,2(1pe13.4))') time_val, delta_ts, z_sl_mean
   write(45,'(1pe13.6,2(1pe13.4))') time_val, delta_ts, z_sl_mean
   write(46,'(1pe13.6,2(1pe13.4))') time_val, delta_ts, z_sl_mean
   write(47,'(1pe13.6,2(1pe13.4))') time_val, delta_ts, z_sl_mean
   write(48,'(1pe13.6,2(1pe13.4))') time_val, delta_ts, z_sl_mean
   write(49,'(1pe13.6,2(1pe13.4))') time_val, delta_ts, z_sl_mean
   write(50,'(1pe13.6,2(1pe13.4))') time_val, delta_ts, z_sl_mean
   write(51,'(1pe13.6,2(1pe13.4))') time_val, delta_ts, z_sl_mean
   write(52,'(1pe13.6,2(1pe13.4))') time_val, delta_ts, z_sl_mean
   write(53,'(1pe13.6,2(1pe13.4))') time_val, delta_ts, z_sl_mean
   write(54,'(1pe13.6,2(1pe13.4))') time_val, delta_ts, z_sl_mean
else if (forcing_flag == 2) then
   write(41,'(1pe13.6,2(1pe13.4))') time_val, glac_index, z_sl_mean
   write(42,'(1pe13.6,2(1pe13.4))') time_val, glac_index, z_sl_mean
   write(43,'(1pe13.6,2(1pe13.4))') time_val, glac_index, z_sl_mean
   write(44,'(1pe13.6,2(1pe13.4))') time_val, glac_index, z_sl_mean
   write(45,'(1pe13.6,2(1pe13.4))') time_val, glac_index, z_sl_mean
   write(46,'(1pe13.6,2(1pe13.4))') time_val, glac_index, z_sl_mean
   write(47,'(1pe13.6,2(1pe13.4))') time_val, glac_index, z_sl_mean
   write(48,'(1pe13.6,2(1pe13.4))') time_val, glac_index, z_sl_mean
   write(49,'(1pe13.6,2(1pe13.4))') time_val, glac_index, z_sl_mean
   write(50,'(1pe13.6,2(1pe13.4))') time_val, glac_index, z_sl_mean
   write(51,'(1pe13.6,2(1pe13.4))') time_val, glac_index, z_sl_mean
   write(52,'(1pe13.6,2(1pe13.4))') time_val, glac_index, z_sl_mean
   write(53,'(1pe13.6,2(1pe13.4))') time_val, glac_index, z_sl_mean
   write(54,'(1pe13.6,2(1pe13.4))') time_val, glac_index, z_sl_mean
end if

do n=1, n_surf-1
   write(41,'(1pe13.4)',advance='no') zl_surf(n)
   write(42,'(1pe13.4)',advance='no') zs_surf(n)
   write(43,'(1pe13.4)',advance='no') accum_surf(n)
   write(44,'(1pe13.4)',advance='no') as_perp_surf(n)
   write(45,'(1pe13.4)',advance='no') snowfall_surf(n)
   write(46,'(1pe13.4)',advance='no') rainfall_surf(n)
   write(47,'(1pe13.4)',advance='no') runoff_surf(n)
   write(48,'(1pe13.4)',advance='no') vx_surf(n)
   write(49,'(1pe13.4)',advance='no') vy_surf(n)
   write(50,'(1pe13.4)',advance='no') vz_surf(n)
   write(51,'(1pe13.4)',advance='no') vx_base(n)
   write(52,'(1pe13.4)',advance='no') vy_base(n)
   write(53,'(1pe13.4)',advance='no') vz_base(n)
   write(54,'(1pe13.4)',advance='no') temp_base_pmp(n) 
end do

n=n_surf
   write(41,'(1pe13.4)') zl_surf(n)
   write(42,'(1pe13.4)') zs_surf(n)
   write(43,'(1pe13.4)') accum_surf(n)
   write(44,'(1pe13.4)') as_perp_surf(n)
   write(45,'(1pe13.4)') snowfall_surf(n)
   write(46,'(1pe13.4)') rainfall_surf(n)
   write(47,'(1pe13.4)') runoff_surf(n)
   write(48,'(1pe13.4)') vx_surf(n)
   write(49,'(1pe13.4)') vy_surf(n)
   write(50,'(1pe13.4)') vz_surf(n)
   write(51,'(1pe13.4)') vx_base(n)
   write(52,'(1pe13.4)') vy_base(n)
   write(53,'(1pe13.4)') vz_base(n)
   write(54,'(1pe13.4)') temp_base_pmp(n) 

deallocate(zl_surf, zs_surf, accum_surf, as_perp_surf, &
           snowfall_surf, rainfall_surf, runoff_surf, &
	   vx_surf, vy_surf, vz_surf, &
	   vx_base, vy_base, vz_base, temp_base_pmp)

end subroutine output5

#endif /* (defined(ASF) && WRITE_SER_FILE_STAKES==1) Austfonna*/

!-------------------------------------------------------------------------------
!> Computation of an arbitrary field quantity for a given borehole
!! position x_pos, y_pos by weighed averaging of the corresponding
!! gridded 2-d field.
!-------------------------------------------------------------------------------
  subroutine borehole(field, x_pos, y_pos, dxi, deta, ch_grid, &
                      field_val, flag_in_domain)

#if !defined(ALLOW_TAPENADE) /* Normal */
  use sico_maths_m, only : bilinint
#else /* Tapenade */
  use sico_maths_m
#endif /* Normal vs. Tapenade */

  implicit none

  real(dp), dimension(0:JMAX,0:IMAX), intent(in) :: field
  real(dp)                          , intent(in) :: x_pos, y_pos, dxi, deta
  character (len=*)                 , intent(in) :: ch_grid

  real(dp), intent(out) :: field_val
  logical , intent(out) :: flag_in_domain

  integer(i4b) :: i1, i2, j1, j2
  real(dp)     :: real_i, real_j

!-------- Neighbour points --------

  real_i = (x_pos-xi(0)) /dxi
  real_j = (y_pos-eta(0))/deta

  if (ch_grid=='sg_x') real_i = real_i - 0.5_dp
  if (ch_grid=='sg_y') real_j = real_j - 0.5_dp

  if (real_i < 0.5_dp*real(IMAX,dp)) then
     i1 = floor(real_i)
     i2 = i1 + 1
  else
     i2 = ceiling(real_i)
     i1 = i2 - 1
  end if

  if (real_j < 0.5_dp*real(JMAX,dp)) then
     j1 = floor(real_j)
     j2 = j1 + 1
  else
     j2 = ceiling(real_j)
     j1 = j2 - 1
  end if

  flag_in_domain = .true.

  if (ch_grid=='grid') then
                    ! field(j,i) defined on grid points

     if ((i1 < 0).or.(j1 < 0).or.(i2 > IMAX).or.(j2 > JMAX)) then
        flag_in_domain = .false.
     end if

  else if (ch_grid=='sg_x') then
                    ! field(j,i) defined on staggered grid in x direction

     if ((i1 < 0).or.(j1 < 0).or.(i2 > IMAX-1).or.(j2 > JMAX)) then
        flag_in_domain = .false.
     end if

  else if (ch_grid=='sg_y') then
                    ! field(j,i) defined on staggered grid in y direction

     if ((i1 < 0).or.(j1 < 0).or.(i2 > IMAX).or.(j2 > JMAX-1)) then
        flag_in_domain = .false.
     end if

  else

     errormsg = ' >>> borehole: Parameter ch_grid has undefined value!'
     call error(errormsg)

  end if

!-------- Weighing of the four adjacent grid values --------

  if (flag_in_domain) then

     field_val = bilinint(real(i1,dp), real(i2,dp), real(j1,dp), real(j2,dp), &
                      field(j1,i1), field(j2,i1), field(j1,i2), field(j2,i2), &
                      real_i, real_j)

  else

     field_val = 0.0_dp

  end if

  end subroutine borehole

!-------------------------------------------------------------------------------
!> Set the creation mode and compression type for NetCDF files.
!-------------------------------------------------------------------------------
  subroutine set_cmode(cmode, n_deflate_level, flag_shuffle)

  use netcdf

  implicit none

  integer(i4b), intent(out) :: cmode
  integer(i4b), intent(out) :: n_deflate_level
  logical,      intent(out) :: flag_shuffle

#if (NETCDF4_ENABLED==1)

  cmode           = ior(NF90_NETCDF4, NF90_CLASSIC_MODEL)
  n_deflate_level = 1
  flag_shuffle    = .true.

#else

  cmode           = NF90_NOCLOBBER
  n_deflate_level = 0
  flag_shuffle    = .false.

#endif

  end subroutine set_cmode

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
!> Set the value of the institution string ch_institution.
!-------------------------------------------------------------------------------
  subroutine set_ch_institution(ch_institution)

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
!> Set the value of the auxiliary variable grads_nc_tweaks.
!-------------------------------------------------------------------------------
  subroutine set_grads_nc_tweaks(grads_nc_tweaks)

  implicit none

  logical, intent(out) :: grads_nc_tweaks

  character(len=16) :: ch_value

  grads_nc_tweaks = .false.

!-------- Try environment variable --------

  call get_environment_variable('SICO_GRADS_NC_TWEAKS', ch_value)

  if ( (trim(ch_value)=='true') &
       .or.(trim(ch_value)=='True').or.(trim(ch_value)=='TRUE') ) &
     grads_nc_tweaks = .true.

  if ( (trim(ch_value)=='yes') &   ! obsolete, but still supported
       .or.(trim(ch_value)=='Yes').or.(trim(ch_value)=='YES') &
       .or.(trim(ch_value)=='y').or.(trim(ch_value)=='Y') ) &
     grads_nc_tweaks = .true.

!-------- Try preprocessor switch --------

#if (defined(GRADS_NC_TWEAKS))
#if (GRADS_NC_TWEAKS==1)
  grads_nc_tweaks = .true.
#else
  grads_nc_tweaks = .false.
#endif
#endif

  end subroutine set_grads_nc_tweaks

!-------------------------------------------------------------------------------

end module output_m
!
