!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ i n i t _ m
!
!! Initializations for SICOPOLIS.
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
!> Initializations for SICOPOLIS.
!-------------------------------------------------------------------------------
module sico_init_m

  use sico_types_m
  use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

  use error_m

  implicit none

  public

contains

!-------------------------------------------------------------------------------
!> Main routine of sico_init_m: Initializations for SICOPOLIS.
!-------------------------------------------------------------------------------
subroutine sico_init(delta_ts, glac_index, &
               mean_accum, &
               dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser, &
               time, time_init, time_end, time_output, &
               dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
               z_mar, &
               ndat2d, ndat3d, n_output)

#if !defined(ALLOW_TAPENADE) /* Normal */
  use compare_float_m
#endif /* Normal */

  use ice_material_properties_m, only : ice_mat_eqs_pars
  use enth_temp_omega_m, only : calc_c_int_table, calc_c_int_inv_table, &
                                enth_fct_temp_omega

  use netcdf
  use nc_check_m

#if (defined(GRL) && DISC>0)
#if (!defined(EXEC_MAKE_C_DIS_0))
  use discharge_workers_m, only: disc_param, disc_fields
#else /* defined(EXEC_MAKE_C_DIS_0) */
  use discharge_workers_m, only: disc_param, disc_fields, calc_c_dis_0
#endif
#endif

#if (GRID==0 || GRID==1)
  use stereo_proj_m
#endif

  use read_m, only : read_target_topo_nc, &
                     read_scalar_input, read_2d_input, read_kei, read_phys_para

#if (defined(ALLOW_NORMAL) || defined(ALLOW_GRDCHK) || defined(ALLOW_TAPENADE))
  use read_m, only : read_age_data, read_BedMachine_data
#endif

  use boundary_m
  use init_temp_water_age_m
  use calc_enhance_m
  use flag_update_gf_gl_cf_m
  use calc_vxy_m
  use calc_vz_m
  use calc_dxyz_m
  use calc_temp_melt_bas_m

#if !defined(ALLOW_TAPENADE) /* Normal */
  use output_m
#endif /* Normal */

implicit none

integer(i4b), intent(out) :: ndat2d, ndat3d
integer(i4b), intent(out) :: n_output
real(dp),     intent(out) :: delta_ts, glac_index
real(dp),     intent(out) :: mean_accum
real(dp),     intent(out) :: dtime, dtime_temp, dtime_wss
real(dp),     intent(out) :: dtime_out, dtime_ser
real(dp),     intent(out) :: time, time_init, time_end, time_output(100)
real(dp),     intent(out) :: dxi, deta, dzeta_c, dzeta_t, dzeta_r
real(dp),     intent(out) :: z_mar

integer(i4b) :: i, j, kc, kt, kr, m, n, ir, jr, n1, n2
integer(i4b) :: ios, ios1, ios2, ios3, ios4
integer(i4b) :: istat, ierr
integer(i4b) :: n_q_geo_mod
real(dp) :: dtime0, dtime_temp0, dtime_wss0, dtime_out0, dtime_ser0
real(dp) :: time_init0, time_end0
#if (OUTPUT==2 || OUTPUT==3)
real(dp) :: time_output0(N_OUTPUT)
#endif
real(dp) :: d_dummy

#if (defined(ANT))
real(dp) :: larmip_qbm_anom_aux(5)
#endif

character(len=256) :: anfdatname, target_topo_dat_name
character(len=256) :: filename_with_path
character(len=256) :: shell_command
character(len=256) :: ch_revision
character(len= 64) :: ch_var_name
character(len=  3) :: ch_month(12)
character          :: ch_dummy
logical            :: flag_precip_monthly_mean
logical            :: flag_init_output, flag_3d_output

integer(i4b), dimension(0:JMAX,0:IMAX) :: mask_present
real(dp),     dimension(0:JMAX,0:IMAX) :: zs_present

#if (TSURFACE<=5)
logical :: flag_temp_zs_ref_file
#endif

#if (ACCSURFACE<=5)
logical :: flag_precip_zs_ref_file
#endif

real(dp), dimension(0:JMAX,0:IMAX) :: field2d_aux
real(dp), dimension(0:IMAX,0:JMAX) :: field2d_tra_aux

integer(i4b) :: n_slide_regions
#if (!defined(N_SLIDE_REGIONS) || N_SLIDE_REGIONS<=1)
integer(i4b) :: p_weert_aux(1)
integer(i4b) :: q_weert_aux(1)
real(dp) :: c_slide_aux(1)
real(dp) :: gamma_slide_aux(1)
#else
integer(i4b) :: p_weert_aux(N_SLIDE_REGIONS)
integer(i4b) :: q_weert_aux(N_SLIDE_REGIONS)
real(dp) :: c_slide_aux(N_SLIDE_REGIONS)
real(dp) :: gamma_slide_aux(N_SLIDE_REGIONS)
#endif

#if (!defined(N_BM_REGIONS) || N_BM_REGIONS<=1)
real(dp) :: gamma0_bm_aux(1)
real(dp) :: delta_tf_bm_aux(1)
#else
real(dp) :: gamma0_bm_aux(N_BM_REGIONS)
real(dp) :: delta_tf_bm_aux(N_BM_REGIONS)
#endif

#if (TSURFACE==6 && ACCSURFACE==6 && ABLSURFACE==6)
real(dp), dimension(0:IMAX,0:JMAX) :: temp_maat_climatol_conv, &
                                      smb_climatol_conv, &
                                      zs_ref_conv
integer(i4b)        :: n_st_unit_length, n_smb_unit_length
character(len=64)   :: ch_st_unit, ch_smb_unit
real(dp), parameter :: temp_C_to_K = 273.15_dp
#endif

#if (defined(ANT) || defined(GRL))
character(len=64) :: ch_initmip_smb_anom_file
#endif
#if (defined(ANT))
character(len=64) :: ch_initmip_bmb_anom_file
character(len=64) :: ch_larmip_regions_file
#endif

#if (FLOATING_ICE_BASAL_MELTING==6)
real(dp), dimension(0:IMAX,0:JMAX,0:NZ_TF_BM) :: tf_bm_present_aux
#endif

#if (defined(ANT) && ICE_SHELF_COLLAPSE_MASK==1)
real(dp), dimension(0:IMAX,0:JMAX) :: H_ref_retreat_conv
#endif

#if (defined(GRL) && RETREAT_MASK==1)
real(dp), dimension(0:IMAX,0:JMAX) :: H_ref_retreat_conv
#endif

integer(i4b) :: dimid, ncid, ncv
!   dimid:       Dimension ID
!    ncid:       File ID
!     ncv:       Variable ID

real(dp), parameter :: inv_twelve = 1.0_dp/12.0_dp

character(len=64), parameter :: thisroutine = 'sico_init'

character(len=64), parameter :: fmt1 = '(a)', &
                                fmt2 = '(a,i0)', &
                                fmt3 = '(a,es13.5)', &
                                fmt4 = '(a,es20.12)'

write(unit=6, fmt='(a)') ' '
write(unit=6, fmt='(a)') ' -------- sico_init --------'
write(unit=6, fmt='(a)') ' '

!-------- Name of the computational domain --------

#if (defined(ANT))
ch_domain_long  = 'Antarctica'
ch_domain_short = 'ant'

#elif (defined(GRL))
ch_domain_long  = 'Greenland'
ch_domain_short = 'grl'

#elif (defined(NHEM))
ch_domain_long  = 'Entire northern hemisphere'
ch_domain_short = 'nhem'

#elif (defined(LCIS))
ch_domain_long  = 'Laurentide and Cordilleran ice sheets'
ch_domain_short = 'lcis'

#elif (defined(SCAND))
ch_domain_long  = 'Fennoscandian and Eurasian ice sheets'
ch_domain_short = 'scand'

#elif (defined(TIBET))
ch_domain_long  = 'Tibetan ice sheet'
ch_domain_short = 'tibet'

#elif (defined(ASF))
ch_domain_long  = 'Austfonna'
ch_domain_short = 'asf'

#elif (defined(NPI))
ch_domain_long  = 'Northern Patagonian ice field'
ch_domain_short = 'npi'

#elif (defined(MOCHO))
ch_domain_long  = 'Mocho-Choshuenco ice cap'
ch_domain_short = 'mocho'

#elif (defined(EISMINT))
ch_domain_long  = 'EISMINT'
ch_domain_short = 'eismint'

#elif (defined(HEINO))
ch_domain_long  = 'ISMIP HEINO'
ch_domain_short = 'heino'

#elif (defined(NMARS))
ch_domain_long  = 'North polar cap of Mars'
ch_domain_short = 'nmars'

#elif (defined(SMARS))
ch_domain_long  = 'South polar cap of Mars'
ch_domain_short = 'smars'

#else
ch_domain_long  = 'Unspecified domain'
ch_domain_short = 'xyz'

#endif

!-------- Some initial values --------

n_output = 0

delta_ts   = 0.0_dp
glac_index = 0.0_dp
z_mar      = 0.0_dp

dtime      = 0.0_dp
dtime_temp = 0.0_dp
dtime_wss  = 0.0_dp
dtime_out  = 0.0_dp
dtime_ser  = 0.0_dp

time        = 0.0_dp
time_init   = 0.0_dp
time_end    = 0.0_dp
time_output = 0.0_dp

!-------- Initialization of the Library of Iterative Solvers Lis,
!                                                     if required --------

#if (MARGIN==3 || DYNAMICS==2)
#if !defined(ALLOW_TAPENADE) /* Normal */
  call lis_initialize(ierr)
#else /* Tapenade */
  call lis_init_f(ierr)
#endif /* Normal vs. Tapenade */
#endif

!-------- Physical parameters --------

#if (defined(YEAR_SEC))
year2sec = YEAR_SEC
#else
year2sec = 3.1556925445e+07_dp
              ! IUPAC-IUGS year for epoch 2000.0
              ! (Holden et al., 2011, PAC, doi:10.1351/PAC-REC-09-01-22)
#endif

sec2year = 1.0_dp/year2sec

#if (defined(PARAM_RHO))
RHO = real(PARAM_RHO,dp)
#else
errormsg = ' >>> sico_init: ' &
           // 'Parameter PARAM_RHO not defined in run-specs header!'
call error(errormsg)
#endif

#if (defined(PARAM_RHO_W))
RHO_W = real(PARAM_RHO_W,dp)
#else
errormsg = ' >>> sico_init: ' &
           // 'Parameter PARAM_RHO_W not defined in run-specs header!'
call error(errormsg)
#endif

#if (defined(PARAM_RHO_SW))
RHO_SW = real(PARAM_RHO_SW,dp)
#else
errormsg = ' >>> sico_init: ' &
           // 'Parameter PARAM_RHO_SW not defined in run-specs header!'
call error(errormsg)
#endif

#if (defined(PARAM_L))
L = real(PARAM_L,dp)
#else
errormsg = ' >>> sico_init: ' &
           // 'Parameter PARAM_L not defined in run-specs header!'
call error(errormsg)
#endif

#if (defined(PARAM_G))
G = real(PARAM_G,dp)
#else
errormsg = ' >>> sico_init: ' &
           // 'Parameter PARAM_G not defined in run-specs header!'
call error(errormsg)
#endif

#if (defined(PARAM_NUE))
NUE = real(PARAM_NUE,dp)
#else
errormsg = ' >>> sico_init: ' &
           // 'Parameter PARAM_NUE not defined in run-specs header!'
call error(errormsg)
#endif

#if (defined(PARAM_BETA))
BETA = real(PARAM_BETA,dp)
#else
errormsg = ' >>> sico_init: ' &
           // 'Parameter PARAM_BETA not defined in run-specs header!'
call error(errormsg)
#endif

#if (defined(PARAM_DELTA_TM_SW))
DELTA_TM_SW = real(PARAM_DELTA_TM_SW,dp)
#else
errormsg = ' >>> sico_init: ' &
           // 'Parameter PARAM_DELTA_TM_SW not defined in run-specs header!'
call error(errormsg)
#endif

#if (defined(PARAM_OMEGA_MAX))
OMEGA_MAX = real(PARAM_OMEGA_MAX,dp)
#else
errormsg = ' >>> sico_init: ' &
           // 'Parameter PARAM_OMEGA_MAX not defined in run-specs header!'
call error(errormsg)
#endif

#if (defined(PARAM_H_R))
H_R = real(PARAM_H_R,dp)
#else
errormsg = ' >>> sico_init: ' &
           // 'Parameter PARAM_H_R not defined in run-specs header!'
call error(errormsg)
#endif

#if (defined(PARAM_RHO_C_R))
RHO_C_R = real(PARAM_RHO_C_R,dp)
#else
errormsg = ' >>> sico_init: ' &
           // 'Parameter PARAM_RHO_C_R not defined in run-specs header!'
call error(errormsg)
#endif

#if (defined(PARAM_KAPPA_R))
KAPPA_R = real(PARAM_KAPPA_R,dp)
#else
errormsg = ' >>> sico_init: ' &
           // 'Parameter PARAM_KAPPA_R not defined in run-specs header!'
call error(errormsg)
#endif

#if (defined(PARAM_RHO_A))
RHO_A = real(PARAM_RHO_A,dp)
#else
errormsg = ' >>> sico_init: ' &
           // 'Parameter PARAM_RHO_A not defined in run-specs header!'
call error(errormsg)
#endif

#if (defined(PARAM_R_T))
R_T = real(PARAM_R_T,dp)
#else
errormsg = ' >>> sico_init: ' &
           // 'Parameter PARAM_R_T not defined in run-specs header!'
call error(errormsg)
#endif

#if (!defined(RF_KAPPA_C_FILE))
errormsg = ' >>> sico_init: ' &
           // 'File RF_KAPPA_C_FILE not defined in run-specs header!'
call error(errormsg)
#endif

call read_phys_para()
     ! read tabulated values of the
     ! rate factor, heat conductivity and specific heat

call ice_mat_eqs_pars(RF, R_T, KAPPA, C, -190, 10)

!  ------ Some auxiliary quantities required for the enthalpy method

call calc_c_int_table(C, -190, 10, L)
call calc_c_int_inv_table()

!-------- Check whether the dynamics and thermodynamics modes are defined

#if (!defined(DYNAMICS))
errormsg = ' >>> sico_init: DYNAMICS not defined in the header file!'
call error(errormsg)
#endif

#if (!defined(CALCMOD))
errormsg = ' >>> sico_init: CALCMOD not defined in the header file!'
call error(errormsg)
#endif

#if (defined(ENTHMOD))
errormsg = ' >>> sico_init: ENTHMOD must not be defined any more.' &
         //         end_of_line &
         //'        Please update your header file!'
call error(errormsg)
#endif

!-------- Compatibility check of the horizontal resolution with the
!         number of grid points --------

#if !defined(ALLOW_TAPENADE) /* Normal */

#if (!defined(CHECK_RES_IMAX_JMAX) || CHECK_RES_IMAX_JMAX==1)

#if (GRID==0 || GRID==1)

#if (defined(ANT)) /* Antarctic ice sheet */

if (approx_equal(DX, 64.0_dp, eps_sp_dp)) then

   if ( (IMAX /= 95).or.(JMAX /= 95) ) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else if (approx_equal(DX, 40.0_dp, eps_sp_dp)) then

   if ( (IMAX /= 152).or.(JMAX /= 152) ) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else if (approx_equal(DX, 32.0_dp, eps_sp_dp)) then

   if ( (IMAX /= 190).or.(JMAX /= 190) ) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else if (approx_equal(DX, 20.0_dp, eps_sp_dp)) then

   if ( (IMAX /= 304).or.(JMAX /= 304) ) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else if (approx_equal(DX, 16.0_dp, eps_sp_dp)) then

   if ( (IMAX /= 380).or.(JMAX /= 380) ) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else if (approx_equal(DX, 10.0_dp, eps_sp_dp)) then

   if ( (IMAX /= 608).or.(JMAX /= 608) ) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else if (approx_equal(DX, 8.0_dp, eps_sp_dp)) then

   if ( (IMAX /= 760).or.(JMAX /= 760) ) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else
   errormsg = ' >>> sico_init: DX wrong!'
   call error(errormsg)
end if

#elif (defined(GRL)) /* Greenland ice sheet */

if (approx_equal(DX, 40.0_dp, eps_sp_dp)) then

   if ( (IMAX /= 42).or.(JMAX /= 72) ) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else if (approx_equal(DX, 20.0_dp, eps_sp_dp)) then

   if ( (IMAX /= 84).or.(JMAX /= 144) ) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else if (approx_equal(DX, 16.0_dp, eps_sp_dp)) then

   if ( (IMAX /= 105).or.(JMAX /= 180) ) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else if (approx_equal(DX, 10.0_dp, eps_sp_dp)) then

   if ( (IMAX /= 168).or.(JMAX /= 288) ) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else if (approx_equal(DX, 8.0_dp, eps_sp_dp)) then

   if ( (IMAX /= 210).or.(JMAX /= 360) ) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else if (approx_equal(DX, 5.0_dp, eps_sp_dp)) then

   if ( (IMAX /= 336).or.(JMAX /= 576) ) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else if (approx_equal(DX, 4.0_dp, eps_sp_dp)) then

   if ( (IMAX /= 420).or.(JMAX /= 720) ) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else if (approx_equal(DX, 2.0_dp, eps_sp_dp)) then

   if ( (IMAX /= 840).or.(JMAX /= 1440) ) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else
   errormsg = ' >>> sico_init: DX wrong!'
   call error(errormsg)
end if

#elif (defined(NHEM)) /* Entire northern hemisphere */

if (approx_equal(DX, 80.0_dp, eps_sp_dp)) then

   if ((IMAX /= 156).or.(JMAX /= 156)) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else if (approx_equal(DX, 40.0_dp, eps_sp_dp)) then

   if ((IMAX /= 312).or.(JMAX /= 312)) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else if (approx_equal(DX, 20.0_dp, eps_sp_dp)) then

   if ((IMAX /= 624).or.(JMAX /= 624)) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else
   errormsg = ' >>> sico_init: DX wrong!'
   call error(errormsg)
end if

#elif (defined(LCIS)) /* Laurentide and Cordilleran ice sheets */

if (approx_equal(DX, 80.0_dp, eps_sp_dp)) then

   if ((IMAX /= 105).or.(JMAX /= 78)) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else if (approx_equal(DX, 40.0_dp, eps_sp_dp)) then

   if ((IMAX /= 211).or.(JMAX /= 157)) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else if (approx_equal(DX, 20.0_dp, eps_sp_dp)) then

   if ((IMAX /= 423).or.(JMAX /= 315)) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else
   errormsg = ' >>> sico_init: DX wrong!'
   call error(errormsg)
end if

#elif (defined(SCAND)) /* Fennoscandian and Eurasian ice sheets */

if (approx_equal(DX, 40.0_dp, eps_sp_dp)) then

   if ((IMAX /= 150).or.(JMAX /= 70)) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else if (approx_equal(DX, 20.0_dp, eps_sp_dp)) then

   if ((IMAX /= 300).or.(JMAX /= 140)) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else if (approx_equal(DX, 10.0_dp, eps_sp_dp)) then

   if ((IMAX /= 600).or.(JMAX /= 280)) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else
   errormsg = ' >>> sico_init: DX wrong!'
   call error(errormsg)
end if

#elif (defined(ASF)) /* Austfonna */

if (approx_equal(DX, 4.0_dp, eps_sp_dp)) then

   if ((IMAX /= 34).or.(JMAX /= 33)) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else if (approx_equal(DX, 2.0_dp, eps_sp_dp)) then

   if ((IMAX /= 68).or.(JMAX /= 66)) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else if (approx_equal(DX, 1.0_dp, eps_sp_dp)) then

   if ((IMAX /= 136).or.(JMAX /= 132)) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else
   errormsg = ' >>> sico_init: DX wrong!'
   call error(errormsg)
end if

#elif (defined(NPI)) /* Northern Patagonian ice field */

if (approx_equal(DX, 0.45_dp, eps_sp_dp)) then

   if ((IMAX /= 182).or.(JMAX /= 314)) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

   if ( (.not.(approx_equal(X0,   66.082d0, eps_sp_dp))) &
        .or. &
        (.not.(approx_equal(Y0, 4288.950d0, eps_sp_dp))) ) then
      errormsg = ' >>> sico_init: X0 and/or Y0 wrong!'
      call error(errormsg)
   end if

else if (approx_equal(DX, 0.9_dp, eps_sp_dp)) then

   if ((IMAX /=  91).or.(JMAX /= 157)) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

   if ( (.not.(approx_equal(X0,   66.082d0, eps_sp_dp))) &
        .or. &
        (.not.(approx_equal(Y0, 4288.950d0, eps_sp_dp))) ) then
      errormsg = ' >>> sico_init: X0 and/or Y0 wrong!'
      call error(errormsg)
   end if

else
   errormsg = ' >>> sico_init: DX wrong!'
   call error(errormsg)
end if

#endif /* Different computational domains */

#elif (GRID==2)

#if (defined(TIBET)) /* Tibetan ice sheet */

if (      (approx_equal(DLAMBDA, 1.0_dp/3.0_dp, eps_sp_dp)) &
     .and.(approx_equal(DPHI   , 1.0_dp/3.0_dp, eps_sp_dp)) ) then

   if ((IMAX /= 135).or.(JMAX /= 51)) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else if (      (approx_equal(DLAMBDA, 1.0_dp/6.0_dp, eps_sp_dp)) &
          .and.(approx_equal(DPHI   , 1.0_dp/6.0_dp, eps_sp_dp)) ) then

   if ((IMAX /= 270).or.(JMAX /= 102)) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else if (      (approx_equal(DLAMBDA, 1.0_dp/12.0_dp, eps_sp_dp)) &
          .and.(approx_equal(DPHI   , 1.0_dp/12.0_dp, eps_sp_dp)) ) then

   if ((IMAX /= 540).or.(JMAX /= 204)) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else

   errormsg = ' >>> sico_init: DLAMBDA / DPHI wrong!'
   call error(errormsg)

end if

#endif

#endif

#else /* CHECK_RES_IMAX_JMAX==0 */

write(6, fmt='(a)') ' >>> sico_init: CHECK_RES_IMAX_JMAX==0'
write(6, fmt='(a)') '      -> compatibility check between horizontal resolution'
write(6, fmt='(a)') '         and number of grid points not performed.'
write(6, fmt='(a)') ' '

#endif /* CHECK_RES_IMAX_JMAX */

#else /* Tapenade */

print *, ' >>> sico_init: grid compatibility check not performed'
print *, '          in adjoint applications; check manually.' 

#endif /* Normal vs. Tapenade */

!-------- Compatibility check of the thermodynamics mode
!         (cold vs. polythermal vs. enthalpy method)
!         and the number of grid points in the lower (kt) ice domain --------

#if (CALCMOD==0 || CALCMOD==2 || CALCMOD==3 || CALCMOD==-1)

if (KTMAX > 2) then
   write(6, fmt='(a)') ' >>> sico_init: For options CALCMOD==0, 2, 3 or -1,'
   write(6, fmt='(a)') '                the separate kt domain is redundant.'
   write(6, fmt='(a)') '                Therefore, consider setting KTMAX to 2.'
   write(6, fmt='(a)') ' '
end if

#endif

!-------- Compatibility check of surface-temperature
!                         and precipitation switches --------

#if (TSURFACE == 5 && ACCSURFACE != 5)
errormsg = ' >>> sico_init: ' &
              //'Options TSURFACE==5 and ACCSURFACE==5 must be used together!'
call error(errormsg)
#endif

#if (TSURFACE != 5 && ACCSURFACE == 5)
errormsg = ' >>> sico_init: ' &
              //'Options TSURFACE==5 and ACCSURFACE==5 must be used together!'
call error(errormsg)
#endif

#if (TSURFACE == 6)
#if (ACCSURFACE != 6 || ABLSURFACE != 6)
errormsg = ' >>> sico_init: Options TSURFACE==6, ACCSURFACE==6, ABLSURFACE==6' &
         //         end_of_line &
         //'        must be used together!'
call error(errormsg)
#endif
#endif

#if (ACCSURFACE == 6)
#if (TSURFACE != 6 || ABLSURFACE != 6)
errormsg = ' >>> sico_init: Options TSURFACE==6, ACCSURFACE==6, ABLSURFACE==6' &
         //         end_of_line &
         //'        must be used together!'
call error(errormsg)
#endif
#endif

#if (ABLSURFACE == 6)
#if (TSURFACE != 6 || ACCSURFACE != 6)
errormsg = ' >>> sico_init: Options TSURFACE==6, ACCSURFACE==6, ABLSURFACE==6' &
         //         end_of_line &
         //'        must be used together!'
call error(errormsg)
#endif
#endif

#if (ACCSURFACE == 7)
#if (ABLSURFACE != 7)
errormsg = ' >>> sico_init: Options ACCSURFACE==7, ABLSURFACE==7' &
         //         end_of_line &
         //'        must be used together!'
call error(errormsg)
#endif
#endif

#if (ABLSURFACE == 7)
#if (ACCSURFACE != 7)
errormsg = ' >>> sico_init: Options ACCSURFACE==7, ABLSURFACE==7' &
         //         end_of_line &
         //'        must be used together!'
call error(errormsg)
#endif
#endif

!-------- Compatibility check of discretization schemes for the horizontal and
!         vertical advection terms in the temperature and age equations --------

#if (ADV_HOR==1)
errormsg = ' >>> sico_init: ' &
              //'Option ADV_HOR==1 (central differences) not defined!'
call error(errormsg)
#endif

!-------- Check whether for the shallow shelf
!               or shelfy stream approximation
!                  the chosen grid is Cartesian coordinates
!                             without distortion correction (GRID==0) --------

#if ((MARGIN==3 && DYNAMICS==1) || DYNAMICS==2)   /* requires SSA or SStA */
#if (GRID != 0)
write(6, fmt='(a)') ' >>> sico_init: WARNING:'
write(6, fmt='(a)') '                Distortion correction for GRID.ne.0'
write(6, fmt='(a)') '                not yet implemented'
write(6, fmt='(a)') '                for the shallow shelf approximation (SSA)'
write(6, fmt='(a)') '                or the shelfy stream approximation (SStA).'
write(6, fmt='(a)') ' '
#endif
#endif

!-------- Setting of forcing flag --------

#if (TSURFACE <= 4)

forcing_flag = 1   ! forcing by delta_ts

#elif (TSURFACE == 5)

forcing_flag = 2   ! forcing by glac_index

#elif (TSURFACE == 6)

forcing_flag = 3   ! forcing by time-dependent surface temperature
                   ! and SMB data

#endif

!-------- Initialization of numerical time steps --------

dtime0      = DTIME0
dtime_temp0 = DTIME_TEMP0
#if (REBOUND==2)
dtime_wss0  = DTIME_WSS0
#endif

!-------- Further initializations --------

#if (defined(PLANET_R))
R = real(PLANET_R,dp)   ! mean radius of the planet
#else
errormsg = ' >>> sico_init: ' &
           // 'Parameter PLANET_R not defined in run-specs header!'
call error(errormsg)
#endif

#if (defined(PLANET_A))
A = real(PLANET_A,dp)   ! semi-major axis of the planet
#else
errormsg = ' >>> sico_init: ' &
           // 'Parameter PLANET_A not defined in run-specs header!'
call error(errormsg)
#endif

#if (defined(PLANET_F_INV))
F_INV = real(PLANET_F_INV,dp)  ! inverse flattening of the planet
#else
errormsg = ' >>> sico_init: ' &
           // 'Parameter PLANET_F_INV not defined in run-specs header!'
call error(errormsg)
#endif

if (F_INV > 1.0e+10_dp) then   ! interpreted as infinity -> no flattening
   B = A
else   ! finite inverse flattening
   B = A - A/F_INV
end if

#if (GRID==0 || GRID==1)

#if (defined(STEREO_PROJ_LATD0))
PHI0 = real(STEREO_PROJ_LATD0,dp) *deg2rad   ! deg -> rad
              ! central meridian of the stereographic projection
#else
errormsg = ' >>> sico_init: ' &
           // 'Parameter STEREO_PROJ_LATD0 not defined in run-specs header!'
call error(errormsg)
#endif

#if (defined(STEREO_PROJ_LOND0))
LAMBDA0 = real(STEREO_PROJ_LOND0,dp) *deg2rad   ! deg -> rad
              ! standard parallel of the stereographic projection
#else
errormsg = ' >>> sico_init: ' &
           // 'Parameter STEREO_PROJ_LOND0 not defined in run-specs header!'
call error(errormsg)
#endif

#else

PHI0    = 0.0_dp   ! dummy value
LAMBDA0 = 0.0_dp   ! dummy value

#endif

dzeta_c = 1.0_dp/real(KCMAX,dp)
dzeta_t = 1.0_dp/real(KTMAX,dp)
dzeta_r = 1.0_dp/real(KRMAX,dp)

ndat2d = 1
ndat3d = 1

flag_calc_temp = .true.

#if (defined(ANT) || defined(GRL))
flag_initmip_asmb = .false.
#endif
#if (defined(ANT))
flag_initmip_abmb = .false.
flag_larmip       = .false.
#endif

#if (ACCSURFACE==2 || ACCSURFACE==3)
gamma_s = GAMMA_S
#endif

#if (ABLSURFACE==1 || ABLSURFACE==2)

#if (defined(S_STAT_0) && defined(BETA1_0) && defined(BETA2_0) && defined(PMAX_0) && defined(MU_0))
s_stat = S_STAT_0
beta1  = BETA1_0 *(0.001_dp/86400.0_dp)*(RHO_W/RHO)
                          ! (mm WE)/(d*degC) -> (m IE)/(s*degC)
beta2  = BETA2_0 *(0.001_dp/86400.0_dp)*(RHO_W/RHO)
                          ! (mm WE)/(d*degC) -> (m IE)/(s*degC)
Pmax   = PMAX_0
mu     = MU_0    *(1000.0_dp*86400.0_dp)*(RHO/RHO_W)
                          ! (d*degC)/(mm WE) -> (s*degC)/(m IE)
#else
errormsg = ' >>> sico_init: ' &
           // 'Parameters for PDD model not defined in run-specs header!'
call error(errormsg)
#endif

#elif (ABLSURFACE==3)

lambda_lti = LAMBDA_LTI *(0.001_dp*sec2year)*(RHO_W/RHO)
                        ! (mm WE)/(a*degC) -> (m IE)/(s*degC)
temp_lti   = TEMP_LTI

#endif

!-------- General abbreviations --------

!  ------ kc domain

if (DEFORM >= eps) then

   flag_aa_nonzero = .true.   ! non-equidistant grid

   aa = DEFORM
   ea = exp(aa)

   kc=0
   zeta_c(kc)         = 0.0_dp
   eaz_c(kc)          = 1.0_dp
   eaz_c_quotient(kc) = 0.0_dp

   do kc=1, KCMAX-1
      zeta_c(kc) = kc*dzeta_c
      eaz_c(kc)  = exp(aa*zeta_c(kc))
      eaz_c_quotient(kc) = (eaz_c(kc)-1.0_dp)/(ea-1.0_dp)
   end do

   kc=KCMAX
   zeta_c(kc)         = 1.0_dp
   eaz_c(kc)          = exp(aa)
   eaz_c_quotient(kc) = 1.0_dp

else

   flag_aa_nonzero = .false.   ! equidistant grid

   aa = 0.0_dp
   ea = 1.0_dp

   kc=0
   zeta_c(kc)         = 0.0_dp
   eaz_c(kc)          = 1.0_dp
   eaz_c_quotient(kc) = 0.0_dp

   do kc=1, KCMAX-1
      zeta_c(kc) = kc*dzeta_c
      eaz_c(kc)  = 1.0_dp
      eaz_c_quotient(kc) = zeta_c(kc)
   end do

   kc=KCMAX
   zeta_c(kc)         = 1.0_dp
   eaz_c(kc)          = 1.0_dp
   eaz_c_quotient(kc) = 1.0_dp

end if

!  ------ kt domain

kt=0
zeta_t(kt) = 0.0_dp

do kt=1, KTMAX-1
   zeta_t(kt) = kt*dzeta_t
end do

kt=KTMAX
zeta_t(kt) = 1.0_dp

!  ------ kr domain

kr=0
zeta_r(kr) = 0.0_dp

do kr=1, KRMAX-1
   zeta_r(kr) = kr*dzeta_r
end do

kr=KRMAX
zeta_r(kr) = 1.0_dp

!-------- Reshaping of a 2-d array (with indices i, j)
!                                  to a vector (with index n) --------

n=1

do i=0, IMAX
do j=0, JMAX
   n2i(n)    = i
   n2j(n)    = j
   ij2n(j,i) = n
   n=n+1
end do
end do

!-------- Specification of current simulation --------

n1 = len('sico_specs_')+1
n2 = len(trim(RUN_SPECS_HEADER))-len('.h')
run_name = trim(RUN_SPECS_HEADER)
run_name = run_name(n1:n2)

anfdatname = trim(ANFDATNAME)

#if (defined(YEAR_ZERO))
year_zero  = YEAR_ZERO
#else
year_zero  = 2000.0_dp   ! default value 2000 CE
#endif

time_init0 = TIME_INIT0
time_end0  = TIME_END0
dtime_ser0 = DTIME_SER0

#if (OUTPUT==1 || OUTPUT==3)
dtime_out0 = DTIME_OUT0
#endif

#if (OUTPUT==2 || OUTPUT==3)

#if (N_OUTPUT<=100)
n_output = N_OUTPUT
#else
errormsg = ' >>> sico_init: N_OUTPUT > 100 not allowed!'
call error(errormsg)
#endif

#if (defined(TIME_OUT0))

time_output0 = TIME_OUT0

#else   /* !defined(TIME_OUT0), legacy mode */

#if (N_OUTPUT>= 1)
time_output0( 1) = TIME_OUT0_01
#endif
#if (N_OUTPUT>= 2)
time_output0( 2) = TIME_OUT0_02
#endif
#if (N_OUTPUT>= 3)
time_output0( 3) = TIME_OUT0_03
#endif
#if (N_OUTPUT>= 4)
time_output0( 4) = TIME_OUT0_04
#endif
#if (N_OUTPUT>= 5)
time_output0( 5) = TIME_OUT0_05
#endif
#if (N_OUTPUT>= 6)
time_output0( 6) = TIME_OUT0_06
#endif
#if (N_OUTPUT>= 7)
time_output0( 7) = TIME_OUT0_07
#endif
#if (N_OUTPUT>= 8)
time_output0( 8) = TIME_OUT0_08
#endif
#if (N_OUTPUT>= 9)
time_output0( 9) = TIME_OUT0_09
#endif
#if (N_OUTPUT>=10)
time_output0(10) = TIME_OUT0_10
#endif
#if (N_OUTPUT>=11)
time_output0(11) = TIME_OUT0_11
#endif
#if (N_OUTPUT>=12)
time_output0(12) = TIME_OUT0_12
#endif
#if (N_OUTPUT>=13)
time_output0(13) = TIME_OUT0_13
#endif
#if (N_OUTPUT>=14)
time_output0(14) = TIME_OUT0_14
#endif
#if (N_OUTPUT>=15)
time_output0(15) = TIME_OUT0_15
#endif
#if (N_OUTPUT>=16)
time_output0(16) = TIME_OUT0_16
#endif
#if (N_OUTPUT>=17)
time_output0(17) = TIME_OUT0_17
#endif
#if (N_OUTPUT>=18)
time_output0(18) = TIME_OUT0_18
#endif
#if (N_OUTPUT>=19)
time_output0(19) = TIME_OUT0_19
#endif
#if (N_OUTPUT>=20)
time_output0(20) = TIME_OUT0_20
#endif

#endif

#endif

call set_flag_grads_nc_tweaks()

!-------- Maximum ice extent yes/no --------

#if (!defined(MASK_MAXEXTENT_FILE) || THK_EVOL==0)

flag_mask_maxextent = .false.   ! no maximum ice extent specified

#else

if ( (trim(adjustl(MASK_MAXEXTENT_FILE)) == 'none') &
     .or. &
     (trim(adjustl(MASK_MAXEXTENT_FILE)) == 'None') &
     .or. &
     (trim(adjustl(MASK_MAXEXTENT_FILE)) == 'NONE') ) then

   flag_mask_maxextent = .false.   ! no maximum ice extent specified

else

   flag_mask_maxextent = .true.   ! maximum ice extent specified

end if

#endif

!-------- Type of the geothermal heat flux (GHF) --------

#if (!defined(Q_GEO_FILE))

n_q_geo_mod = 1   ! spatially constant GHF

#else

if ( (trim(adjustl(Q_GEO_FILE)) == 'none') &
     .or. &
     (trim(adjustl(Q_GEO_FILE)) == 'None') &
     .or. &
     (trim(adjustl(Q_GEO_FILE)) == 'NONE') ) then

   n_q_geo_mod = 1   ! spatially constant GHF

else

   n_q_geo_mod = 2   ! spatially varying GHF

end if

#endif

!-------- Write log file --------

shell_command = 'if [ ! -d'
shell_command = trim(shell_command)//' '//OUT_PATH
shell_command = trim(shell_command)//' '//'] ; then mkdir'
shell_command = trim(shell_command)//' '//OUT_PATH
shell_command = trim(shell_command)//' '//'; fi'
call system(trim(shell_command))
     ! Check whether directory OUT_PATH exists. If not, it is created.

filename_with_path = trim(OUT_PATH)//'/'//trim(run_name)//'.log'

#if !defined(ALLOW_TAPENADE) /* Normal */
open(10, iostat=ios, file=trim(filename_with_path), status='new')
#else /* Tapenade */
open(10, iostat=ios, file=trim(filename_with_path))
#endif /* Normal vs. Tapenade */

if (ios /= 0) then
   errormsg = ' >>> sico_init: Error when opening the log file!'
   call error(errormsg)
end if

write(10, fmt=trim(fmt1)) 'Computational domain:'
write(10, fmt=trim(fmt1)) '   '//trim(ch_domain_long)
write(10, fmt=trim(fmt1)) ' '

#if (defined(PARAM_RHO))
write(10, fmt=trim(fmt3)) 'RHO         =', PARAM_RHO
#endif

#if (defined(PARAM_RHO_W))
write(10, fmt=trim(fmt3)) 'RHO_W       =', PARAM_RHO_W
#endif

#if (defined(PARAM_RHO_SW))
write(10, fmt=trim(fmt3)) 'RHO_SW      =', PARAM_RHO_SW
#endif

#if (defined(PARAM_L))
write(10, fmt=trim(fmt3)) 'L           =', PARAM_L
#endif

#if (defined(PARAM_G))
write(10, fmt=trim(fmt3)) 'G           =', PARAM_G
#endif

#if (defined(PARAM_NUE))
write(10, fmt=trim(fmt3)) 'NUE         =', PARAM_NUE
#endif

#if (defined(PARAM_BETA))
write(10, fmt=trim(fmt3)) 'BETA        =', PARAM_BETA
#endif

#if (defined(PARAM_DELTA_TM_SW))
write(10, fmt=trim(fmt3)) 'DELTA_TM_SW =', PARAM_DELTA_TM_SW
#endif

#if (defined(PARAM_OMEGA_MAX))
write(10, fmt=trim(fmt3)) 'OMEGA_MAX   =', PARAM_OMEGA_MAX
#endif

#if (defined(PARAM_H_R))
write(10, fmt=trim(fmt3)) 'H_R         =', PARAM_H_R
#endif

#if (defined(PARAM_RHO_C_R))
write(10, fmt=trim(fmt3)) 'RHO_C_R     =', PARAM_RHO_C_R
#endif

#if (defined(PARAM_KAPPA_R))
write(10, fmt=trim(fmt3)) 'KAPPA_R     =', PARAM_KAPPA_R
#endif

#if (defined(PARAM_RHO_A))
write(10, fmt=trim(fmt3)) 'RHO_A       =', PARAM_RHO_A
#endif

#if (defined(PARAM_R_T))
write(10, fmt=trim(fmt3)) 'R_T         =', PARAM_R_T
#endif

write(10, fmt=trim(fmt1)) ' '

#if (defined(RF_KAPPA_C_FILE))
write(10, fmt=trim(fmt1)) 'RF_KAPPA_C_FILE = ' &
                          // trim(adjustl(RF_KAPPA_C_FILE))
#endif
write(10, fmt=trim(fmt1)) ' '

write(10, fmt=trim(fmt2)) 'GRID = ', GRID
write(10, fmt=trim(fmt1)) ' '

#if (defined(PLANET_R))
write(10, fmt=trim(fmt4)) 'R =', PLANET_R
#endif

#if (defined(PLANET_A))
write(10, fmt=trim(fmt4)) 'A =', PLANET_A
#endif

#if (defined(PLANET_F_INV))
write(10, fmt=trim(fmt4)) 'F_INV =', PLANET_F_INV
#endif

#if (GRID==0 || GRID==1)

#if (defined(STEREO_PROJ_LATD0))
write(10, fmt=trim(fmt3)) 'LATD0 =', STEREO_PROJ_LATD0
#endif

#if (defined(STEREO_PROJ_LOND0))
write(10, fmt=trim(fmt3)) 'LOND0 =', STEREO_PROJ_LOND0
#endif

#endif

write(10, fmt=trim(fmt1)) ' '

write(10, fmt=trim(fmt2)) 'imax  = ', IMAX
write(10, fmt=trim(fmt2)) 'jmax  = ', JMAX
write(10, fmt=trim(fmt2)) 'kcmax = ', KCMAX
write(10, fmt=trim(fmt2)) 'ktmax = ', KTMAX
write(10, fmt=trim(fmt2)) 'krmax = ', KRMAX
write(10, fmt=trim(fmt1)) ' '

write(10, fmt=trim(fmt3)) 'a =', aa
write(10, fmt=trim(fmt1)) ' '

#if (GRID==0 || GRID==1)
write(10, fmt=trim(fmt3)) 'x0 =', X0
write(10, fmt=trim(fmt3)) 'y0 =', Y0
write(10, fmt=trim(fmt3)) 'dx =', DX
#elif (GRID==2)
write(10, fmt=trim(fmt3)) 'lambda0 =', LAMBDA_0
write(10, fmt=trim(fmt3)) 'phi0    =', PHI_0
write(10, fmt=trim(fmt3)) 'dlambda =', DLAMBDA
write(10, fmt=trim(fmt3)) 'dphi    =', DPHI
#endif
write(10, fmt=trim(fmt1)) ' '

#if (defined(CHECK_RES_IMAX_JMAX))
write(10, fmt=trim(fmt2)) 'CHECK_RES_IMAX_JMAX = ', CHECK_RES_IMAX_JMAX
write(10, fmt=trim(fmt1)) ' '
#endif

write(10, fmt=trim(fmt3)) 'year_zero  =', year_zero
write(10, fmt=trim(fmt3)) 'time_init  =', time_init0
write(10, fmt=trim(fmt3)) 'time_end   =', time_end0
write(10, fmt=trim(fmt3)) 'dtime      =', dtime0
write(10, fmt=trim(fmt3)) 'dtime_temp =', dtime_temp0
#if (REBOUND==2)
write(10, fmt=trim(fmt3)) 'dtime_wss  =', dtime_wss0
#endif
#if (defined(GRL) && DISC>0)
write(10, fmt=trim(fmt3)) 'dtime_mar_coa =', DTIME_MAR_COA0
#endif
write(10, fmt=trim(fmt1)) ' '

write(10, fmt=trim(fmt2)) 'DYNAMICS = ', DYNAMICS
#if (DYNAMICS==2 && defined(HYB_MODE))
write(10, fmt=trim(fmt2)) 'HYB_MODE = ', HYB_MODE
#endif
#if ((DYNAMICS==1 && MARGIN==3) || DYNAMICS==2)
#if (defined(LIS_OPTS))
write(10, fmt=trim(fmt1)) 'lis_opts = '//LIS_OPTS
#endif
#if (defined(TOL_ITER_SSA))
write(10, fmt=trim(fmt3)) 'tol_iter_ssa =', TOL_ITER_SSA
#endif
#if (defined(N_ITER_SSA))
write(10, fmt=trim(fmt2)) 'n_iter_ssa = ', N_ITER_SSA
#endif
#if (defined(N_ITER_SSA_MIN))
write(10, fmt=trim(fmt2)) 'n_iter_ssa_min = ', N_ITER_SSA_MIN
#endif
#if (defined(ITER_INIT_SSA))
write(10, fmt=trim(fmt2)) 'iter_init_ssa = ', ITER_INIT_SSA
#endif
#if (defined(VISC_INIT_SSA))
write(10, fmt=trim(fmt3)) 'visc_init_ssa =', VISC_INIT_SSA
#endif
#if (defined(N_VISC_SMOOTH))
write(10, fmt=trim(fmt2)) 'n_visc_smooth = ', N_VISC_SMOOTH
#endif
#if (defined(VISC_SMOOTH_DIFF))
write(10, fmt=trim(fmt3)) 'visc_smooth_diff =', VISC_SMOOTH_DIFF
#endif
#if (defined(RELAX_FACT_SSA))
write(10, fmt=trim(fmt3)) 'relax_fact_ssa =', RELAX_FACT_SSA
#endif
#endif
#if (DYNAMICS==2 && HYB_MODE==0 && defined(RATIO_SL_THRESH))
write(10, fmt=trim(fmt3)) 'ratio_sl_thresh =', RATIO_SL_THRESH
#endif
#if (DYNAMICS==2 && HYB_MODE==0 && defined(SSTA_SIA_WEIGH_FCT))
write(10, fmt=trim(fmt2)) 'SSTA_SIA_WEIGH_FCT = ', SSTA_SIA_WEIGH_FCT
#endif
#if (DYNAMICS==2 && HYB_MODE==1 && defined(HYB_REF_SPEED))
write(10, fmt=trim(fmt3)) 'hyb_ref_speed =', HYB_REF_SPEED
#endif
write(10, fmt=trim(fmt1)) ' '

write(10, fmt=trim(fmt2)) 'CALCMOD = ', CALCMOD
#if (CALCMOD==-1 && defined(TEMP_CONST))
write(10, fmt=trim(fmt3)) 'TEMP_CONST =', TEMP_CONST
#endif
#if (CALCMOD==-1 && defined(AGE_CONST))
write(10, fmt=trim(fmt3)) 'AGE_CONST =', AGE_CONST
#endif
#if (CALCMOD==1 && defined(CTS_MELTING_FREEZING))
write(10, fmt=trim(fmt2)) 'CTS_MELTING_FREEZING = ', CTS_MELTING_FREEZING
#endif
write(10, fmt=trim(fmt1)) ' '

write(10, fmt=trim(fmt2)) 'MARGIN = ', MARGIN
#if (MARGIN==2)
write(10, fmt=trim(fmt2)) 'MARINE_ICE_FORMATION = ', MARINE_ICE_FORMATION
write(10, fmt=trim(fmt2)) 'MARINE_ICE_CALVING   = ', MARINE_ICE_CALVING
#if (MARINE_ICE_CALVING==2 || MARINE_ICE_CALVING==3)
write(10, fmt=trim(fmt3)) 'z_mar =', Z_MAR
#elif (MARINE_ICE_CALVING==4 || MARINE_ICE_CALVING==5 || MARINE_ICE_CALVING==6 || MARINE_ICE_CALVING==7)
write(10, fmt=trim(fmt3)) 'fact_z_mar =', FACT_Z_MAR
#elif (MARINE_ICE_FORMATION==2 && MARINE_ICE_CALVING==9)
write(10, fmt=trim(fmt3)) 'calv_uw_coeff =', CALV_UW_COEFF
write(10, fmt=trim(fmt3)) 'r1_calv_uw    =', R1_CALV_UW
write(10, fmt=trim(fmt3)) 'r2_calv_uw    =', R2_CALV_UW
#endif
#elif (MARGIN==3)
write(10, fmt=trim(fmt2)) 'ICE_SHELF_CALVING = ', ICE_SHELF_CALVING
#if (ICE_SHELF_CALVING==2)
write(10, fmt=trim(fmt3)) 'H_calv =', H_CALV
#endif
#endif
write(10, fmt=trim(fmt1)) ' '

write(10, fmt=trim(fmt2)) 'FLOW_LAW = ', FLOW_LAW
write(10, fmt=trim(fmt2)) 'FIN_VISC = ', FIN_VISC
#if (FLOW_LAW==2)
write(10, fmt=trim(fmt3)) 'gr_size =', GR_SIZE
#endif
#if (FIN_VISC==2)
write(10, fmt=trim(fmt3)) 'sigma_res =', SIGMA_RES
#endif
write(10, fmt=trim(fmt1)) ' '

write(10, fmt=trim(fmt2)) 'ENHMOD = ', ENHMOD
#if (ENHMOD==1 || ENHMOD==2 || ENHMOD==3)
write(10, fmt=trim(fmt3)) 'enh_fact =', ENH_FACT
#endif
#if (ENHMOD==2 || ENHMOD==3)
write(10, fmt=trim(fmt3)) 'enh_intg =', ENH_INTG
#endif
#if (ENHMOD==2)
write(10, fmt=trim(fmt3)) 'age_trans =', AGE_TRANS_0
#endif
#if (ENHMOD==3)
write(10, fmt=trim(fmt3)) 'date_trans1 =', DATE_TRANS1_0
write(10, fmt=trim(fmt3)) 'date_trans2 =', DATE_TRANS2_0
write(10, fmt=trim(fmt3)) 'date_trans3 =', DATE_TRANS3_0
#endif
#if (ENHMOD==4 || ENHMOD==5)
write(10, fmt=trim(fmt3)) 'enh_compr =', ENH_COMPR
write(10, fmt=trim(fmt3)) 'enh_shear =', ENH_SHEAR
#endif
#if (DYNAMICS==2 && defined(ENH_STREAM))
if (ENH_STREAM >= 0.0_dp) &
   write(10, fmt=trim(fmt3)) 'enh_stream =', ENH_STREAM
#endif
#if ((ENHMOD==1 || ENHMOD==2 || ENHMOD==3 || ENHMOD==4) && MARGIN==3)
write(10, fmt=trim(fmt3)) 'enh_shelf =', ENH_SHELF
#endif
write(10, fmt=trim(fmt1)) ' '

write(10, fmt=trim(fmt2)) 'ANF_DAT = ', ANF_DAT
write(10, fmt=trim(fmt1)) 'zs_present file = '//ZS_PRESENT_FILE
#if (ANF_DAT==1)
#if (defined(ZB_PRESENT_FILE))
write(10, fmt=trim(fmt1)) 'zb_present file = '//ZB_PRESENT_FILE
#endif
write(10, fmt=trim(fmt1)) 'zl_present file = '//ZL_PRESENT_FILE
#endif
write(10, fmt=trim(fmt1)) 'zl0 file = '//ZL0_FILE
write(10, fmt=trim(fmt1)) 'mask_present file = '//MASK_PRESENT_FILE
#if (defined(MASK_REGION_FILE))
if ( (trim(adjustl(MASK_REGION_FILE)) /= 'none') &
     .and. &
     (trim(adjustl(MASK_REGION_FILE)) /= 'None') &
     .and. &
     (trim(adjustl(MASK_REGION_FILE)) /= 'NONE') ) then
   write(10, fmt=trim(fmt1)) 'mask_region file = '//MASK_REGION_FILE
   write(10, fmt=trim(fmt1)) ' '
end if
#endif
#if (ANF_DAT==1)
write(10, fmt=trim(fmt2)) 'TEMP_INIT = ', TEMP_INIT
#if (TEMP_INIT==1 && defined(TEMP_INIT_VAL))
write(10, fmt=trim(fmt3)) 'temp_init_val =', TEMP_INIT_VAL
#endif
#endif
#if (ANF_DAT==3 || (ANF_DAT==1 && TEMP_INIT==5))
write(10, fmt=trim(fmt1)) 'Initial-value file = '//ANFDATNAME
write(10, fmt=trim(fmt1)) 'Path to initial-value file = '//ANF_DAT_PATH
#endif
#if (ANF_DAT==3 && defined(RESTART))
write(10, fmt=trim(fmt2)) 'RESTART = ', RESTART
#endif
write(10, fmt=trim(fmt1)) ' '

#if (defined(THK_EVOL))
write(10, fmt=trim(fmt2)) 'THK_EVOL = ', THK_EVOL
#else
errormsg = ' >>> sico_init: Define THK_EVOL in header file!'
call error(errormsg)
#endif
#if (defined(CALCTHK))
write(10, fmt=trim(fmt2)) 'CALCTHK = ', CALCTHK
#else
errormsg = ' >>> sico_init: Define CALCTHK in header file!'
call error(errormsg)
#endif
#if (defined(OCEAN_CONNECTIVITY))
write(10, fmt=trim(fmt2)) 'OCEAN_CONNECTIVITY = ', OCEAN_CONNECTIVITY
#endif
#if (defined(H_ISOL_MAX))
write(10, fmt=trim(fmt3)) 'H_isol_max =', H_ISOL_MAX
#endif

#if (THK_EVOL==2)
write(10, fmt=trim(fmt1)) 'Target-topography relaxation-time file = ' &
                          //TARGET_TOPO_TAU0_FILE
write(10, fmt=trim(fmt1)) 'Target-topography file = '//TARGET_TOPO_DAT_NAME
write(10, fmt=trim(fmt1)) 'Path to target-topography file = '//TARGET_TOPO_PATH
#endif

#if (THK_EVOL==3)
write(10, fmt=trim(fmt3)) 'target_topo_tau_0 =', TARGET_TOPO_TAU0
write(10, fmt=trim(fmt1)) 'Target-topography file = '//TARGET_TOPO_DAT_NAME
write(10, fmt=trim(fmt1)) 'Path to target-topography file = '//TARGET_TOPO_PATH
#endif

#if (defined(MASK_MAXEXTENT_FILE))
if (flag_mask_maxextent) &
   write(10, fmt=trim(fmt1)) 'Maximum ice extent mask file = ' &
                             // trim(adjustl(MASK_MAXEXTENT_FILE))
#endif

#if (CALCTHK==2)
write(10, fmt=trim(fmt3))  'ovi_weight =', OVI_WEIGHT
write(10, fmt=trim(fmt3))  'omega_sor =', OMEGA_SOR
#if (ITER_MAX_SOR>0)
write(10, fmt=trim(fmt2)) 'iter_max_sor = ', ITER_MAX_SOR
#endif
#endif
write(10, fmt=trim(fmt1)) ' '

write(10, fmt=trim(fmt2)) 'ADV_HOR  = ', ADV_HOR
write(10, fmt=trim(fmt2)) 'ADV_VERT = ', ADV_VERT
write(10, fmt=trim(fmt2)) 'TOPOGRAD = ', TOPOGRAD
#if (MARGIN==3 && defined(GL_SURF_GRAD))
write(10, fmt=trim(fmt2)) 'GL_SURF_GRAD = ', GL_SURF_GRAD
#endif
write(10, fmt=trim(fmt1)) ' '

write(10, fmt=trim(fmt2)) 'TSURFACE = ', TSURFACE

#if (TSURFACE<=5)

flag_temp_zs_ref_file = .false.

#if (defined(ANT) || defined(GRL)) /* Antarctica or Greenland */

write(10, fmt=trim(fmt2)) 'TEMP_PRESENT_PARA = ', TEMP_PRESENT_PARA
#if (defined(TEMP_PRESENT_OFFSET))
write(10, fmt=trim(fmt3))  'TEMP_PRESENT_OFFSET =', TEMP_PRESENT_OFFSET
#endif

#else /* other than Antarctica or Greenland */

write(10, fmt=trim(fmt1)) 'temp_present file = '//TEMP_PRESENT_FILE
#if (defined(TOPO_LAPSE_RATE))
write(10, fmt=trim(fmt3)) 'topo_lapse_rate =', TOPO_LAPSE_RATE
#endif

#if (defined(TEMP_ZS_REF_FILE))
if ( (trim(adjustl(TEMP_ZS_REF_FILE)) /= 'none') &
     .and. &
     (trim(adjustl(TEMP_ZS_REF_FILE)) /= 'None') &
     .and. &
     (trim(adjustl(TEMP_ZS_REF_FILE)) /= 'NONE') ) then
   flag_temp_zs_ref_file = .true.
   write(10, fmt=trim(fmt1)) 'temp_zs_ref_file = '//TEMP_ZS_REF_FILE
end if
#endif

#endif

#endif

#if (TSURFACE==1)
write(10, fmt=trim(fmt3)) 'delta_ts0 =', DELTA_TS0
#elif (TSURFACE==3)
write(10, fmt=trim(fmt3)) 'sine_amplit =', SINE_AMPLIT
write(10, fmt=trim(fmt3)) 'sine_period =', SINE_PERIOD
#elif (TSURFACE==4)
write(10, fmt=trim(fmt1)) 'GRIP file = '//GRIP_TEMP_FILE
write(10, fmt=trim(fmt3)) 'grip_temp_fact =', GRIP_TEMP_FACT
#elif (TSURFACE==5)
write(10, fmt=trim(fmt1)) 'Glacial-index file = '//GLAC_IND_FILE
write(10, fmt=trim(fmt1)) 'temp_anom file  = '//TEMP_ANOM_FILE
write(10, fmt=trim(fmt3)) 'temp_anom fact  = ', TEMP_ANOM_FACT
#endif
write(10, fmt=trim(fmt1)) ' '

write(10, fmt=trim(fmt2)) 'ACCSURFACE = ', ACCSURFACE

#if (ACCSURFACE<=5)

#if (!defined(PRECIP_PRESENT_FILE))
errormsg = ' >>> sico_init: PRECIP_PRESENT_FILE not defined in the header file!'
call error(errormsg)
#endif

#if (!defined(PRECIP_MA_PRESENT_FILE))
errormsg = ' >>> sico_init: ' &
              //'PRECIP_MA_PRESENT_FILE not defined in the header file!'
call error(errormsg)
#endif

#if (defined(PRECIP_PRESENT_FILE) && defined(PRECIP_MA_PRESENT_FILE))
if ( (trim(adjustl(PRECIP_PRESENT_FILE)) /= 'none') &
     .and. &
     (trim(adjustl(PRECIP_PRESENT_FILE)) /= 'None') &
     .and. &
     (trim(adjustl(PRECIP_PRESENT_FILE)) /= 'NONE') ) then
   write(10, fmt=trim(fmt1)) 'precip_present_file = '//PRECIP_PRESENT_FILE
   flag_precip_monthly_mean = .true.
else if ( (trim(adjustl(PRECIP_MA_PRESENT_FILE)) /= 'none') &
          .and. &
          (trim(adjustl(PRECIP_MA_PRESENT_FILE)) /= 'None') &
          .and. &
          (trim(adjustl(PRECIP_MA_PRESENT_FILE)) /= 'NONE') ) then
   write(10, fmt=trim(fmt1)) 'precip_ma_present_file = '//PRECIP_MA_PRESENT_FILE
   flag_precip_monthly_mean = .false.
else
   errormsg = ' >>> sico_init: Neither PRECIP_PRESENT_FILE' &
            //         end_of_line &
            //'        nor PRECIP_MA_PRESENT_FILE' &
            //         end_of_line &
            //'        specified in the header file!'
   call error(errormsg)
end if
#endif

flag_precip_zs_ref_file = .false.
#if (defined(PRECIP_ZS_REF_FILE))
if ( (trim(adjustl(PRECIP_ZS_REF_FILE)) /= 'none') &
     .and. &
     (trim(adjustl(PRECIP_ZS_REF_FILE)) /= 'None') &
     .and. &
     (trim(adjustl(PRECIP_ZS_REF_FILE)) /= 'NONE') ) then
   flag_precip_zs_ref_file = .true.
   write(10, fmt=trim(fmt1)) 'precip_zs_ref_file = '//PRECIP_ZS_REF_FILE
end if
#endif

#endif

#if (ACCSURFACE==1)
write(10, fmt=trim(fmt3)) 'accfact =', ACCFACT
#elif (ACCSURFACE==2 || ACCSURFACE==3)
write(10, fmt=trim(fmt3)) 'gamma_s =', GAMMA_S
#endif
#if (ACCSURFACE<=5)
write(10, fmt=trim(fmt2)) 'ELEV_DESERT = ', ELEV_DESERT
#if (ELEV_DESERT == 1)
write(10, fmt=trim(fmt3)) 'gamma_p   =', GAMMA_P
write(10, fmt=trim(fmt3)) 'zs_thresh =', ZS_THRESH
#endif
#endif
#if (ACCSURFACE==5)
write(10, fmt=trim(fmt1)) 'precip_anom file = '//PRECIP_ANOM_FILE
write(10, fmt=trim(fmt3)) 'precip_anom fact = ', PRECIP_ANOM_FACT
write(10, fmt=trim(fmt2)) 'PRECIP_ANOM_INTERPOL = ', PRECIP_ANOM_INTERPOL
#endif
#if (ACCSURFACE<=5)
write(10, fmt=trim(fmt2)) 'SOLID_PRECIP = ', SOLID_PRECIP
#endif
write(10, fmt=trim(fmt1)) ' '

write(10, fmt=trim(fmt2)) 'ABLSURFACE = ', ABLSURFACE

#if (ABLSURFACE==1 || ABLSURFACE==2)
#if (defined(S_STAT_0) && defined(BETA1_0) && defined(BETA2_0) && defined(PMAX_0) && defined(MU_0))
write(10, fmt=trim(fmt3)) 's_stat =', S_STAT_0
write(10, fmt=trim(fmt3)) 'beta1  =', BETA1_0
write(10, fmt=trim(fmt3)) 'beta2  =', BETA2_0
write(10, fmt=trim(fmt3)) 'Pmax   =', PMAX_0
write(10, fmt=trim(fmt3)) 'mu     =', MU_0
#else
errormsg = ' >>> sico_init: ' &
           // 'Parameters for PDD model not defined in run-specs header!'
call error(errormsg)
#endif
#elif (ABLSURFACE==3)
write(10, fmt=trim(fmt3)) 'lambda_lti =', LAMBDA_LTI
write(10, fmt=trim(fmt3)) 'temp_lti   =', TEMP_LTI
#endif

#if (defined(MB_ACCOUNT))
write(10, fmt=trim(fmt2)) 'MB_ACCOUNT = ', MB_ACCOUNT
#endif
write(10, fmt=trim(fmt1)) ' '

#if (TSURFACE==6 && ACCSURFACE==6 && ABLSURFACE==6)
write(10, fmt=trim(fmt1)) 'TEMP_SMB_CLIMATOLOGY_FILE = ' &
                             //TEMP_SMB_CLIMATOLOGY_FILE
write(10, fmt=trim(fmt1)) 'TEMP_SMB_ANOM_DIR = '//TEMP_SMB_ANOM_DIR
write(10, fmt=trim(fmt1)) 'TEMP_ANOM_SUBDIR  = '//TEMP_ANOM_SUBDIR
write(10, fmt=trim(fmt1)) 'TEMP_ANOM_FILES   = '//TEMP_ANOM_FILES
write(10, fmt=trim(fmt1)) 'dTEMPdz_SUBDIR    = '//dTEMPdz_SUBDIR
write(10, fmt=trim(fmt1)) 'dTEMPdz_FILES     = '//dTEMPdz_FILES
write(10, fmt=trim(fmt1)) 'SMB_ANOM_SUBDIR   = '//SMB_ANOM_SUBDIR
write(10, fmt=trim(fmt1)) 'SMB_ANOM_FILES    = '//SMB_ANOM_FILES
write(10, fmt=trim(fmt1)) 'dSMBdz_SUBDIR     = '//dSMBdz_SUBDIR
write(10, fmt=trim(fmt1)) 'dSMBdz_FILES      = '//dSMBdz_FILES
write(10, fmt=trim(fmt2)) 'TEMP_SMB_ANOM_TIME_MIN = ', TEMP_SMB_ANOM_TIME_MIN
write(10, fmt=trim(fmt2)) 'TEMP_SMB_ANOM_TIME_MAX = ', TEMP_SMB_ANOM_TIME_MAX
write(10, fmt=trim(fmt1)) ' '
#endif

#if (ACCSURFACE==7 && ABLSURFACE==7)
write(10, fmt=trim(fmt3)) 'target_topo_tau_0 =', TARGET_TOPO_TAU0
write(10, fmt=trim(fmt1)) 'Target-topography file = '//TARGET_TOPO_DAT_NAME
write(10, fmt=trim(fmt1)) 'Path to target-topography file = '//TARGET_TOPO_PATH
write(10, fmt=trim(fmt1)) ' '
#endif

#if (defined(SMB_CORR_FILE))
if ( (trim(adjustl(SMB_CORR_FILE)) /= 'none') &
     .and. &
     (trim(adjustl(SMB_CORR_FILE)) /= 'None') &
     .and. &
     (trim(adjustl(SMB_CORR_FILE)) /= 'NONE') ) then
   write(10, fmt=trim(fmt1)) 'smb_corr_file = '//SMB_CORR_FILE
   write(10, fmt=trim(fmt1)) ' '
end if
#endif

#if ((defined(ANT) || defined(GRL)) && defined(INITMIP_SMB_ANOM_FILE))
if ( (trim(adjustl(INITMIP_SMB_ANOM_FILE)) /= 'none') &
     .and. &
     (trim(adjustl(INITMIP_SMB_ANOM_FILE)) /= 'None') &
     .and. &
     (trim(adjustl(INITMIP_SMB_ANOM_FILE)) /= 'NONE') ) then
   flag_initmip_asmb = .true.
   ch_initmip_smb_anom_file = trim(adjustl(INITMIP_SMB_ANOM_FILE))
   write(10, fmt=trim(fmt1)) 'initmip_smb_anom file = ' &
                                // trim(ch_initmip_smb_anom_file)
   write(10, fmt=trim(fmt1)) ' '
end if
#endif

#if (defined(ANT) && defined(ICE_SHELF_COLLAPSE_MASK))
write(10, fmt=trim(fmt2)) 'ICE_SHELF_COLLAPSE_MASK = ', ICE_SHELF_COLLAPSE_MASK
#if (ICE_SHELF_COLLAPSE_MASK==1)
write(10, fmt=trim(fmt1)) 'ICE_SHELF_COLLAPSE_MASK_DIR   = ' &
                             //ICE_SHELF_COLLAPSE_MASK_DIR
write(10, fmt=trim(fmt1)) 'ICE_SHELF_COLLAPSE_MASK_FILES = ' &
                             //ICE_SHELF_COLLAPSE_MASK_FILES
write(10, fmt=trim(fmt1)) 'ICE_SHELF_COLLAPSE_MASK_H_REF_FILE = ' &
                             //ICE_SHELF_COLLAPSE_MASK_H_REF_FILE
write(10, fmt=trim(fmt2)) 'ICE_SHELF_COLLAPSE_MASK_TIME_MIN = ', &
                             ICE_SHELF_COLLAPSE_MASK_TIME_MIN
write(10, fmt=trim(fmt2)) 'ICE_SHELF_COLLAPSE_MASK_TIME_MAX = ', &
                             ICE_SHELF_COLLAPSE_MASK_TIME_MAX
#endif
write(10, fmt=trim(fmt1)) ' '
#endif

#if (defined(GRL) && defined(RETREAT_MASK))
write(10, fmt=trim(fmt2)) 'RETREAT_MASK = ', RETREAT_MASK
#if (RETREAT_MASK==1)
write(10, fmt=trim(fmt1)) 'RETREAT_MASK_DIR   = '//RETREAT_MASK_DIR
write(10, fmt=trim(fmt1)) 'RETREAT_MASK_FILES = '//RETREAT_MASK_FILES
write(10, fmt=trim(fmt1)) 'RETREAT_MASK_H_REF_FILE = '//RETREAT_MASK_H_REF_FILE
write(10, fmt=trim(fmt2)) 'RETREAT_MASK_TIME_MIN = ', RETREAT_MASK_TIME_MIN
write(10, fmt=trim(fmt2)) 'RETREAT_MASK_TIME_MAX = ', RETREAT_MASK_TIME_MAX
#endif
write(10, fmt=trim(fmt1)) ' '
#endif

#if (defined(GRL) && defined(DISC))
write(10, fmt=trim(fmt2)) 'DISC = ', DISC
#if (DISC>0)
write(10, fmt=trim(fmt3)) 'c_dis_0   =', C_DIS_0
write(10, fmt=trim(fmt3)) 'c_dis_fac =', C_DIS_FAC
write(10, fmt=trim(fmt3)) 'm_H       =', M_H
write(10, fmt=trim(fmt3)) 'm_D       =', M_D
write(10, fmt=trim(fmt3)) 'r_mar_eff =', R_MAR_EFF
#if (defined(S_DIS))
write(10, fmt=trim(fmt3)) 's_dis     =', S_DIS
#endif
#if (defined(ALPHA_SUB))
write(10, fmt=trim(fmt3)) 'alpha_sub =', ALPHA_SUB
#endif
#if (defined(ALPHA_O))
write(10, fmt=trim(fmt3)) 'alpha_o   =', ALPHA_O
#endif
#endif
write(10, fmt=trim(fmt1)) ' '
#endif

write(10, fmt=trim(fmt2)) 'SEA_LEVEL = ', SEA_LEVEL
#if (SEA_LEVEL==1)
write(10, fmt=trim(fmt3)) 'z_sl0 =', Z_SL0
#elif (SEA_LEVEL==3)
write(10, fmt=trim(fmt1)) 'sea-level file = '//SEA_LEVEL_FILE
#endif
write(10, fmt=trim(fmt1)) ' '

#if (defined(BASAL_HYDROLOGY))
write(10, fmt=trim(fmt2)) 'BASAL_HYDROLOGY = ', BASAL_HYDROLOGY
#if (BASAL_HYDROLOGY==1 && defined(MELT_DRAIN))
write(10, fmt=trim(fmt2)) 'MELT_DRAIN = ', MELT_DRAIN
#endif
#endif

write(10, fmt=trim(fmt2)) 'SLIDE_LAW = ', SLIDE_LAW

#if (SLIDE_LAW==1)

#if (defined(N_SLIDE_REGIONS))
write(10, fmt=trim(fmt2)) 'N_SLIDE_REGIONS = ', N_SLIDE_REGIONS
#if (N_SLIDE_REGIONS>1)
write(10, fmt=trim(fmt1)) 'SLIDE_REGIONS_FILE = '//SLIDE_REGIONS_FILE
#endif
#endif

#if (!defined(N_SLIDE_REGIONS) || N_SLIDE_REGIONS<=1)
n_slide_regions = 1
#else
n_slide_regions = N_SLIDE_REGIONS
#endif

#if (defined(BASAL_WATER_PRESSURE))
write(10, fmt=trim(fmt2)) 'BASAL_WATER_PRESSURE = ', BASAL_WATER_PRESSURE
#endif

c_slide_aux = C_SLIDE
gamma_slide_aux = GAMMA_SLIDE
p_weert_aux = P_WEERT
q_weert_aux = Q_WEERT

write(10, fmt=trim(fmt3)) 'c_slide =', c_slide_aux(1)
#if (N_SLIDE_REGIONS>1)
do n=2, n_slide_regions
   write(10, fmt=trim(fmt3)) '         ', c_slide_aux(n)
end do
#endif

write(10, fmt=trim(fmt3)) 'gamma_slide =', gamma_slide_aux(1)
#if (N_SLIDE_REGIONS>1)
do n=2, n_slide_regions
   write(10, fmt=trim(fmt3)) '             ', gamma_slide_aux(n)
end do
#endif

write(10, fmt=trim(fmt2)) 'p_weert = ', p_weert_aux(1)
#if (N_SLIDE_REGIONS>1)
do n=2, n_slide_regions
   write(10, fmt=trim(fmt2)) '          ', p_weert_aux(n)
end do
#endif

write(10, fmt=trim(fmt2)) 'q_weert = ', q_weert_aux(1)
#if (N_SLIDE_REGIONS>1)
do n=2, n_slide_regions
   write(10, fmt=trim(fmt2)) '          ', q_weert_aux(n)
end do
#endif

#if (defined(C_SLIDE_FILTER_WIDTH))
write(10, fmt=trim(fmt3)) 'c_slide_filter_width =', C_SLIDE_FILTER_WIDTH
#endif
#if (defined(TIME_RAMP_UP_SLIDE))
write(10, fmt=trim(fmt3)) 'time_ramp_up_slide =', TIME_RAMP_UP_SLIDE
#endif
write(10, fmt=trim(fmt3)) 'red_pres_limit_fact =', RED_PRES_LIMIT_FACT
#if (BASAL_HYDROLOGY==1 && defined(HYDRO_SLIDE_SAT_FCT) && defined(C_HW_SLIDE) && defined(HW0_SLIDE))
write(10, fmt=trim(fmt2)) 'HYDRO_SLIDE_SAT_FCT = ', HYDRO_SLIDE_SAT_FCT
write(10, fmt=trim(fmt3)) 'c_Hw_slide =', C_HW_SLIDE
write(10, fmt=trim(fmt3)) 'Hw0_slide  =', HW0_SLIDE
#endif

#endif

write(10, fmt=trim(fmt1)) ' '

if (n_q_geo_mod==1) then
   write(10, fmt=trim(fmt3)) 'q_geo =', Q_GEO
else if (n_q_geo_mod==2) then
   write(10, fmt=trim(fmt1)) 'q_geo file = '//Q_GEO_FILE
end if
write(10, fmt=trim(fmt2)) 'Q_LITHO = ', Q_LITHO
write(10, fmt=trim(fmt1)) ' '

#if (defined(MARINE_ICE_BASAL_MELTING))
write(10, fmt=trim(fmt2)) 'MARINE_ICE_BASAL_MELTING = ', MARINE_ICE_BASAL_MELTING
#if (MARINE_ICE_BASAL_MELTING==2 || MARINE_ICE_BASAL_MELTING==3)
write(10, fmt=trim(fmt3)) 'qbm_marine =', QBM_MARINE
#endif
write(10, fmt=trim(fmt1)) ' '
#endif

#if (MARGIN==3)

write(10, fmt=trim(fmt2)) 'FLOATING_ICE_BASAL_MELTING = ', FLOATING_ICE_BASAL_MELTING
#if (FLOATING_ICE_BASAL_MELTING==1)
write(10, fmt=trim(fmt3)) 'qbm_float_1 =', QBM_FLOAT_1
#endif
write(10, fmt=trim(fmt3)) 'qbm_float_3 =', QBM_FLOAT_3
write(10, fmt=trim(fmt3)) 'z_abyss =', Z_ABYSS
#if (FLOATING_ICE_BASAL_MELTING==4)
write(10, fmt=trim(fmt3)) 'temp_ocean =', TEMP_OCEAN
write(10, fmt=trim(fmt3)) 'Omega_qbm  =', OMEGA_QBM
write(10, fmt=trim(fmt3)) 'alpha_qbm  =', ALPHA_QBM
#endif
write(10, fmt=trim(fmt3)) 'H_w_0 =', H_W_0
write(10, fmt=trim(fmt1)) ' '

#if (FLOATING_ICE_BASAL_MELTING==6)
write(10, fmt=trim(fmt2)) 'n_bm_regions = ', N_BM_REGIONS
write(10, fmt=trim(fmt1)) 'bm_regions_file = '//BM_REGIONS_FILE
gamma0_bm_aux   = GAMMA0_BM
delta_tf_bm_aux = DELTA_TF_BM
write(10, fmt=trim(fmt3)) 'gamma0_bm =', gamma0_bm_aux(1)
do n=2, N_BM_REGIONS
   write(10, fmt=trim(fmt3)) '           ', gamma0_bm_aux(n)
end do
write(10, fmt=trim(fmt3)) 'delta_tf_bm =', delta_tf_bm_aux(1)
do n=2, N_BM_REGIONS
   write(10, fmt=trim(fmt3)) '             ', delta_tf_bm_aux(n)
end do
write(10, fmt=trim(fmt1)) 'tf_bm_present_file = '//TF_BM_PRESENT_FILE
write(10, fmt=trim(fmt1)) 'tf_bm_dir   = '//TF_BM_DIR
write(10, fmt=trim(fmt1)) 'tf_bm_files = '//TF_BM_FILES
write(10, fmt=trim(fmt2)) 'tf_bm_time_min = ', TF_BM_TIME_MIN
write(10, fmt=trim(fmt2)) 'tf_bm_time_max = ', TF_BM_TIME_MAX
write(10, fmt=trim(fmt3)) 'zmin_tf_bm =',  ZMIN_TF_BM
write(10, fmt=trim(fmt2)) 'nz_tf_bm   = ', NZ_TF_BM
write(10, fmt=trim(fmt3)) 'dz_tf_bm   =',  DZ_TF_BM
write(10, fmt=trim(fmt1)) ' '
#endif

#if (defined(ANT) && defined(INITMIP_BMB_ANOM_FILE))
if ( (trim(adjustl(INITMIP_BMB_ANOM_FILE)) /= 'none') &
     .and. &
     (trim(adjustl(INITMIP_BMB_ANOM_FILE)) /= 'None') &
     .and. &
     (trim(adjustl(INITMIP_BMB_ANOM_FILE)) /= 'NONE') ) then
   flag_initmip_abmb = .true.
   ch_initmip_bmb_anom_file = trim(adjustl(INITMIP_BMB_ANOM_FILE))
   write(10, fmt=trim(fmt1)) 'initmip_bmb_anom file = ' &
                                // trim(ch_initmip_bmb_anom_file)
   write(10, fmt=trim(fmt1)) ' '
end if
#endif

#if (defined(ANT) && defined(LARMIP_REGIONS_FILE))
if ( (trim(adjustl(LARMIP_REGIONS_FILE)) /= 'none') &
     .and. &
     (trim(adjustl(LARMIP_REGIONS_FILE)) /= 'None') &
     .and. &
     (trim(adjustl(LARMIP_REGIONS_FILE)) /= 'NONE') ) then
   flag_larmip = .true.
   ch_larmip_regions_file = trim(adjustl(LARMIP_REGIONS_FILE))
   write(10, fmt=trim(fmt1)) 'larmip_regions_file = ' &
                                // trim(ch_larmip_regions_file)
   larmip_qbm_anom_aux = LARMIP_QBM_ANOM
   write(10, fmt=trim(fmt3)) 'larmip_qbm_anom_1 =', larmip_qbm_anom_aux(1)
   write(10, fmt=trim(fmt3)) 'larmip_qbm_anom_2 =', larmip_qbm_anom_aux(2)
   write(10, fmt=trim(fmt3)) 'larmip_qbm_anom_3 =', larmip_qbm_anom_aux(3)
   write(10, fmt=trim(fmt3)) 'larmip_qbm_anom_4 =', larmip_qbm_anom_aux(4)
   write(10, fmt=trim(fmt3)) 'larmip_qbm_anom_5 =', larmip_qbm_anom_aux(5)
   write(10, fmt=trim(fmt1)) ' '
end if
#endif

#endif /* MARGIN==3 */

write(10, fmt=trim(fmt2)) 'REBOUND = ', REBOUND
#if (REBOUND==1)
write(10, fmt=trim(fmt3)) 'frac_llra =', FRAC_LLRA
#endif
#if (REBOUND==1 || REBOUND==2)
write(10, fmt=trim(fmt2)) 'TIME_LAG_MOD = ', TIME_LAG_MOD
#if (TIME_LAG_MOD==1)
write(10, fmt=trim(fmt3)) 'time_lag =', TIME_LAG
#elif (TIME_LAG_MOD==2)
write(10, fmt=trim(fmt1)) 'time_lag_file = '//TIME_LAG_FILE
#else
errormsg = ' >>> sico_init: TIME_LAG_MOD must be either 1 or 2!'
call error(errormsg)
#endif
#endif
#if (REBOUND==2)
write(10, fmt=trim(fmt2)) 'FLEX_RIG_MOD = ', FLEX_RIG_MOD
#if (FLEX_RIG_MOD==1)
write(10, fmt=trim(fmt3)) 'flex_rig =', FLEX_RIG
#elif (FLEX_RIG_MOD==2)
write(10, fmt=trim(fmt1)) 'flex_rig_file = '//FLEX_RIG_FILE
#else
errormsg = ' >>> sico_init: FLEX_RIG_MOD must be either 1 or 2!'
call error(errormsg)
#endif
#endif
write(10, fmt=trim(fmt1)) ' '

write(10, fmt=trim(fmt3)) 'numdiff_H_t =', NUMDIFF_H_T
write(10, fmt=trim(fmt3)) 'tau_cts     =', TAU_CTS
write(10, fmt=trim(fmt3)) 'vh_max      =', VH_MAX
write(10, fmt=trim(fmt3)) 'hd_min      =', HD_MIN
write(10, fmt=trim(fmt3)) 'hd_max      =', HD_MAX
#if (defined(VISC_MIN) && defined(VISC_MAX))
write(10, fmt=trim(fmt3)) 'visc_min    =', VISC_MIN
write(10, fmt=trim(fmt3)) 'visc_max    =', VISC_MAX
#endif
write(10, fmt=trim(fmt3)) 'qbm_min     =', QBM_MIN
write(10, fmt=trim(fmt3)) 'qbm_max     =', QBM_MAX
write(10, fmt=trim(fmt3)) 'age_min     =', AGE_MIN
write(10, fmt=trim(fmt3)) 'age_max     =', AGE_MAX
write(10, fmt=trim(fmt3)) 'mean_accum  =', MEAN_ACCUM
#if (ADV_VERT==1)
write(10, fmt=trim(fmt3)) 'age_diff    =', AGEDIFF
#endif
write(10, fmt=trim(fmt1)) ' '

#if (defined(NETCDF4_ENABLED))
write(10, fmt=trim(fmt2)) 'NETCDF4_ENABLED = ', NETCDF4_ENABLED
#endif
#if (defined(OUT_TIMES))
write(10, fmt=trim(fmt2)) 'OUT_TIMES = ', OUT_TIMES
#endif
#if (defined(OUTPUT_INIT))
write(10, fmt=trim(fmt2)) 'OUTPUT_INIT = ', OUTPUT_INIT
#endif
write(10, fmt=trim(fmt2)) 'OUTPUT = ', OUTPUT
#if (OUTPUT==1 || OUTPUT==3)
write(10, fmt=trim(fmt3))  'dtime_out =' , dtime_out0
#endif
write(10, fmt=trim(fmt3))  'dtime_ser =' , dtime_ser0
#if (OUTPUT==1 || OUTPUT==2)
write(10, fmt=trim(fmt2)) 'ERGDAT = ', ERGDAT
#endif
#if (defined(OUTPUT_FLUX_VARS))
write(10, fmt=trim(fmt2)) 'OUTPUT_FLUX_VARS = ', OUTPUT_FLUX_VARS
#endif
#if (defined(CLIMATOLOGY_EXTRACTION_HACK))
    !%% Climatology extraction hack (must not be used routinely)!
write(10, fmt=trim(fmt1)) '!!! CLIMATOLOGY_EXTRACTION_HACK defined !!!'
#endif
#if (OUTPUT==2 || OUTPUT==3)
write(10, fmt=trim(fmt2)) 'n_output = ', n_output
do n=1, n_output
   if (n==1) then
      write(10, fmt=trim(fmt3))  'time_output =' , time_output0(n)
   else
      write(10, fmt=trim(fmt3))  '             ' , time_output0(n)
   end if
end do
#endif
write(10, fmt=trim(fmt1)) ' '

#if (defined(ASF) && defined(WRITE_SER_FILE_STAKES)) /* Austfonna */
write(10, fmt=trim(fmt2)) 'WRITE_SER_FILE_STAKES = ', WRITE_SER_FILE_STAKES
write(10, fmt=trim(fmt1)) ' '
#endif

call get_environment_variable(name='REPO_REVISION', value=ch_revision, &
                              status=istat, trim_name=.true.)
write(10, fmt=trim(fmt1)) 'Program version and date: '//VERSION//' / '//DATE
write(10, fmt=trim(fmt1)) 'Git revision identifier : ' // trim(ch_revision)

close(10, status='keep')

#if (defined(CLIMATOLOGY_EXTRACTION_HACK))
    !%% Climatology extraction hack (must not be used routinely)!
#if (OUTPUT_FLUX_VARS==2)
write(6,'(/1x,a)') &
   '>>> sico_init: CLIMATOLOGY_EXTRACTION_HACK defined!'
write(6, '(1x,a)') &
   '               Only for extracting a climatology,'
write(6, '(1x,a/)') &
   '               must not be used routinely!'
#else
errormsg = ' >>> sico_init: CLIMATOLOGY_EXTRACTION_HACK requires' &
         //                 end_of_line &
         //'                OUTPUT_FLUX_VARS==2!'
call error(errormsg)
#endif
#endif

!-------- Conversion of time quantities --------

year_zero  = year_zero*year2sec     ! a -> s
time_init  = time_init0*year2sec    ! a -> s
time_end   = time_end0*year2sec     ! a -> s
dtime      = dtime0*year2sec        ! a -> s
dtime_temp = dtime_temp0*year2sec   ! a -> s
#if (REBOUND==2)
dtime_wss  = dtime_wss0*year2sec    ! a -> s
#endif
dtime_ser  = dtime_ser0*year2sec    ! a -> s
#if (OUTPUT==1 || OUTPUT==3)
dtime_out  = dtime_out0*year2sec    ! a -> s
#endif
#if (OUTPUT==2 || OUTPUT==3)
do n=1, n_output
   time_output(n) = time_output0(n)*year2sec  ! a -> s
end do
#endif

#if !defined(ALLOW_TAPENADE) /* Normal */

if (.not.approx_integer_multiple(dtime_temp, dtime, eps_sp_dp)) then
   errormsg = ' >>> sico_init: dtime_temp must be a multiple of dtime!'
   call error(errormsg)
end if

#if (REBOUND==2)
if (.not.approx_integer_multiple(dtime_wss, dtime, eps_sp_dp)) then
   errormsg = ' >>> sico_init: dtime_wss must be a multiple of dtime!'
   call error(errormsg)
end if
#endif

if (.not.approx_integer_multiple(dtime_ser, dtime, eps_sp_dp)) then
   errormsg = ' >>> sico_init: dtime_ser must be a multiple of dtime!'
   call error(errormsg)
end if

#if (OUTPUT==1 || OUTPUT==3)
if (.not.approx_integer_multiple(dtime_out, dtime, eps_sp_dp)) then
   errormsg = ' >>> sico_init: dtime_out must be a multiple of dtime!'
   call error(errormsg)
end if
#endif

#else /* Tapenade */

print *, ' >>> sico_init: not checking that time steps are '
print *, '                multiples of each other in adjoint mode;'
print *, '                check manually.'

#endif /* Normal vs. Tapenade */

#if (THK_EVOL==2)

filename_with_path = trim(IN_PATH)//'/general/'//trim(TARGET_TOPO_TAU0_FILE)

call read_scalar_input(filename_with_path, &
                       'target_topo_tau0', ndata_target_topo_tau0_max, &
                       target_topo_tau0_time_min, target_topo_tau0_time_stp, &
                       target_topo_tau0_time_max, &
                       ndata_target_topo_tau0, target_topo_tau0)

target_topo_tau0 = target_topo_tau0 *year2sec   ! a -> s

#endif

#if (THK_EVOL==3)
target_topo_tau_0 = TARGET_TOPO_TAU0 *year2sec   ! a -> s
#endif

#if (ACCSURFACE==7 && ABLSURFACE==7)
target_topo_tau_0 = TARGET_TOPO_TAU0 *year2sec   ! a -> s
#endif

time = time_init

!-------- Reading of present-day ice-land-ocean mask mask
!                                and surface elevation --------

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(MASK_PRESENT_FILE)

call read_2d_input(filename_with_path, &
                   ch_var_name='mask', n_var_type=3, n_ascii_header=6, &
                   field2d_r=field2d_aux)

mask_present = nint(field2d_aux)

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(ZS_PRESENT_FILE)

call read_2d_input(filename_with_path, &
                   ch_var_name='zs', n_var_type=1, n_ascii_header=6, &
                   field2d_r=field2d_aux)

zs_present = field2d_aux

do i=0, IMAX
do j=0, JMAX
   if (mask_present(j,i) >= 2) zs_present(j,i) = max(zs_present(j,i), 0.0_dp)
                 ! resetting negative elevations (bathymetry data)
                 ! to the present-day sea surface
end do
end do

!-------- Reading of present-day
!         monthly mean or mean annual precipitation rate --------

#if (ACCSURFACE<=5)

#if (defined(PRECIP_PRESENT_FILE) && defined(PRECIP_MA_PRESENT_FILE))

if (flag_precip_monthly_mean) then
   filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                        trim(PRECIP_PRESENT_FILE)
else
   filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                        trim(PRECIP_MA_PRESENT_FILE)
end if

#endif

if (flag_precip_monthly_mean) then

   ch_month = [ 'jan', 'feb', 'mar', 'apr', 'may', 'jun', &
                'jul', 'aug', 'sep', 'oct', 'nov', 'dec' ]

   do n=1, 12   ! month counter

      ch_var_name = 'precip_present_' // trim(ch_month(n))

      call read_2d_input(filename_with_path, &
                         ch_var_name=trim(ch_var_name), &
                         n_var_type=0, n_ascii_header=6+3*n+(JMAX+1)*(n-1), &
                         field2d_r=field2d_aux)

      precip_present(:,:,n) = field2d_aux *(1.0e-03_dp*sec2year)*(RHO_W/RHO)
                                           ! mm/a water equiv. -> m/s ice equiv.

   end do

else

   call read_2d_input(filename_with_path, &
                      ch_var_name='precip_ma_present', &
                      n_var_type=1, n_ascii_header=6, &
                      field2d_r=field2d_aux)

   precip_ma_present = field2d_aux *(1.0e-03_dp*sec2year)*(RHO_W/RHO)
                                    ! mm/a water equiv. -> m/s ice equiv.

end if

!  ------ Computation of the still undefined present-day
!         mean annual or monthly mean precipitation rate

if (flag_precip_monthly_mean) then

   precip_ma_present = 0.0_dp   ! initialization

   do n=1, 12   ! month counter
      do i=0, IMAX
      do j=0, JMAX
         precip_ma_present(j,i) = precip_ma_present(j,i) &
                                  + precip_present(j,i,n)*inv_twelve
      end do
      end do
   end do

else

   do n=1, 12   ! month counter
      do i=0, IMAX
      do j=0, JMAX
         precip_present(j,i,n) = precip_ma_present(j,i)
                ! monthly means assumed to be equal
                ! to the mean annual precipitation rate
      end do
      end do
   end do

end if

#endif

!-------- Reading of LGM monthly-mean precipitation-rate anomalies --------

#if (ACCSURFACE==5)

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(PRECIP_ANOM_FILE)

ch_month = [ 'jan', 'feb', 'mar', 'apr', 'may', 'jun', &
             'jul', 'aug', 'sep', 'oct', 'nov', 'dec' ]

do n=1, 12   ! month counter

   ch_var_name = 'precip_lgm_anom_' // trim(ch_month(n))

   call read_2d_input(filename_with_path, &
                      ch_var_name=trim(ch_var_name), &
                      n_var_type=0, n_ascii_header=6+3*n+(JMAX+1)*(n-1), &
                      field2d_r=field2d_aux)

   precip_lgm_anom(:,:,n) = field2d_aux

end do

precip_lgm_anom = precip_lgm_anom * PRECIP_ANOM_FACT

do i=0, IMAX
do j=0, JMAX

#if (PRECIP_ANOM_INTERPOL==1)
   do n=1, 12   ! month counter
      gamma_precip_lgm_anom(j,i,n) = 0.0_dp   ! dummy values
   end do
#elif (PRECIP_ANOM_INTERPOL==2)
   do n=1, 12   ! month counter
      gamma_precip_lgm_anom(j,i,n) = -log(precip_lgm_anom(j,i,n))
   end do
#else
   errormsg = ' >>> sico_init: Wrong value of switch PRECIP_ANOM_INTERPOL!'
   call error(errormsg)
#endif

end do
end do

#endif

!-------- Mean accumulation --------

mean_accum = MEAN_ACCUM*(1.0e-03_dp*sec2year)*(RHO_W/RHO)
                       ! mm/a water equiv. -> m/s ice equiv.

!-------- Read file defining the regions for the sliding laws --------

#if (!defined(N_SLIDE_REGIONS) || N_SLIDE_REGIONS<=1)

n_slide_region = 1

#else

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(SLIDE_REGIONS_FILE)

call read_2d_input(filename_with_path, &
                   ch_var_name='n_basin', n_var_type=2, n_ascii_header=6, &
                   field2d_r=field2d_aux)

n_slide_region = nint(field2d_aux)

#endif

!-------- Ice shelf basal melting --------

n_bm_region = 0   ! initialization

#if (FLOATING_ICE_BASAL_MELTING==6)

!  ------ Read file defining the regions for ice shelf basal melting

#if (!defined(N_BM_REGIONS) || N_BM_REGIONS<=1)

n_bm_region = 1

#else

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(BM_REGIONS_FILE)

call read_2d_input(filename_with_path, &
                   ch_var_name='n_bm_region', n_var_type=2, n_ascii_header=6, &
                   field2d_r=field2d_aux)

n_bm_region = nint(field2d_aux)

#endif

!  ------ Read file with the present-day thermal forcing data of the ocean

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(TF_BM_PRESENT_FILE)

ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

if (ios /= nf90_noerr) then
   errormsg = ' >>> sico_init: Error when opening the file' &
            //                 end_of_line &
            //'                for the present-day thermal forcing data' &
            //                 end_of_line &
            //'                of the ocean!'
   call error(errormsg)
end if

call check( nf90_inq_varid(ncid, 'z', ncv) )
call check( nf90_get_var(ncid, ncv, z_tf_bm_present) )

call check( nf90_inq_varid(ncid, 'thermal_forcing', ncv) )
call check( nf90_get_var(ncid, ncv, tf_bm_present_aux) )

call check( nf90_close(ncid) )

if ( (z_tf_bm_present(0) < eps_dp).and.(z_tf_bm_present(NZ_TF_BM) < eps_dp) ) &
   z_tf_bm_present = -z_tf_bm_present   ! ensure positive depth values

do i=0, IMAX
do j=0, JMAX
do n=0, NZ_TF_BM
   if (isnan(tf_bm_present_aux(i,j,n))) then
      tf_bm_present(n,j,i) = no_value_neg_2
   else
      tf_bm_present(n,j,i) = tf_bm_present_aux(i,j,n)
                             ! swap indices -> SICOPOLIS standard
   end if
end do
end do
end do

if (.not.(approx_equal(z_tf_bm_present(0), ZMIN_TF_BM, eps_sp_dp))) then
   errormsg = ' >>> sico_init: Inconsistency between' &
            //         end_of_line &
            //'        read z_tf_bm data' &
            //         end_of_line &
            //'        and parameter ZMIN_TF_BM!'
   call error(errormsg)
end if

if (.not.(approx_equal(z_tf_bm_present(NZ_TF_BM)-z_tf_bm_present(0), &
                       NZ_TF_BM*DZ_TF_BM, eps_sp_dp))) then
   errormsg = ' >>> sico_init: Inconsistency between' &
            //         end_of_line &
            //'        read z_tf_bm data' &
            //         end_of_line &
            //'        and parameters NZ_TF_BM, DZ_TF_BM!'
   call error(errormsg)
end if

#endif   /* (FLOATING_ICE_BASAL_MELTING==6) */

!-------- Reading of the prescribed target topography --------

#if ( (THK_EVOL==2 || THK_EVOL==3) || (ACCSURFACE==7 && ABLSURFACE==7) )

target_topo_dat_name = trim(TARGET_TOPO_DAT_NAME)

call read_target_topo_nc(target_topo_dat_name)

#endif

!-------- Reading of the maximum ice extent mask --------

mask_maxextent = 1   ! default (no constraint)

#if (defined(MASK_MAXEXTENT_FILE))

if (flag_mask_maxextent) then

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(MASK_MAXEXTENT_FILE)

call read_2d_input(filename_with_path, &
                   ch_var_name='mask_maxextent', &
                   n_var_type=3, n_ascii_header=6, &
                   field2d_r=field2d_aux)

mask_maxextent = nint(field2d_aux)

end if

#endif

!-------- Reading of present-day monthly mean surface temperature --------

#if (!defined(ANT) && !defined(GRL)) /* other than Antarctica or Greenland */

#if (TSURFACE<=5)

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(TEMP_PRESENT_FILE)

ch_month = [ 'jan', 'feb', 'mar', 'apr', 'may', 'jun', &
             'jul', 'aug', 'sep', 'oct', 'nov', 'dec' ]

do n=1, 12   ! month counter

   ch_var_name = 'temp_present_' // trim(ch_month(n))

   call read_2d_input(filename_with_path, &
                      ch_var_name=trim(ch_var_name), &
                      n_var_type=0, n_ascii_header=6+3*n+(JMAX+1)*(n-1), &
                      field2d_r=field2d_aux)

   temp_present(:,:,n) = field2d_aux

end do

#endif

#endif

!-------- Reading of LGM monthly-mean surface-temperature anomalies --------

#if (TSURFACE==5)

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(TEMP_ANOM_FILE)

ch_month = [ 'jan', 'feb', 'mar', 'apr', 'may', 'jun', &
             'jul', 'aug', 'sep', 'oct', 'nov', 'dec' ]

do n=1, 12   ! month counter

   ch_var_name = 'temp_lgm_anom_' // trim(ch_month(n))

   call read_2d_input(filename_with_path, &
                      ch_var_name=trim(ch_var_name), &
                      n_var_type=0, n_ascii_header=6+3*n+(JMAX+1)*(n-1), &
                      field2d_r=field2d_aux)

   temp_lgm_anom(:,:,n) = field2d_aux

end do

temp_lgm_anom = temp_lgm_anom * TEMP_ANOM_FACT

#endif

!-------- Reference elevations
!         (for surface-temperature and precipitation data) --------

#if (TSURFACE<=5)

if (flag_temp_zs_ref_file) then

#if (defined(TEMP_ZS_REF_FILE))
   filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                        trim(TEMP_ZS_REF_FILE)
#endif

   call read_2d_input(filename_with_path, &
                      ch_var_name='zs', n_var_type=1, n_ascii_header=6, &
                      field2d_r=field2d_aux)

   zs_ref_temp = max(field2d_aux, 0.0_dp)
            ! resetting negative elevations (bathymetry data)
            ! to the present-day sea surface


else
   zs_ref_temp = zs_present
end if

#endif

#if (ACCSURFACE<=5)

if (flag_precip_zs_ref_file) then

#if (defined(PRECIP_ZS_REF_FILE))
   filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                        trim(PRECIP_ZS_REF_FILE)
#endif

   call read_2d_input(filename_with_path, &
                      ch_var_name='zs', n_var_type=1, n_ascii_header=6, &
                      field2d_r=field2d_aux)

   zs_ref_precip = max(field2d_aux, 0.0_dp)
            ! resetting negative elevations (bathymetry data)
            ! to the present-day sea surface

else
   zs_ref_precip = zs_present
end if

#endif

!-------- Read data for delta_ts --------

#if (TSURFACE==4)

filename_with_path = trim(IN_PATH)//'/general/'//trim(GRIP_TEMP_FILE)

call read_scalar_input(filename_with_path, &
                       'delta_ts', ndata_grip_max, &
                       grip_time_min, grip_time_stp, grip_time_max, &
                       ndata_grip, griptemp)

#endif

!-------- Read data for the glacial index --------

#if (TSURFACE==5)

filename_with_path = trim(IN_PATH)//'/general/'//trim(GLAC_IND_FILE)

call read_scalar_input(filename_with_path, &
                       'gi', ndata_gi_max, &
                       gi_time_min, gi_time_stp, gi_time_max, &
                       ndata_gi, glacial_index)

#endif

!-------- Reading of the surface-temperature and SMB climatology --------

#if (TSURFACE==6 && ACCSURFACE==6 && ABLSURFACE==6)

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)// &
                                   '/'//trim(TEMP_SMB_CLIMATOLOGY_FILE)

ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

if (ios /= nf90_noerr) then
   errormsg = ' >>> sico_init: Error when opening the file' &
            //                 end_of_line &
            //'                for the surface-temperature and SMB climatology!'
   call error(errormsg)
end if

call check( nf90_inq_varid(ncid, 'ST_clim', ncv) )
call check( nf90_get_var(ncid, ncv, temp_maat_climatol_conv) )
call check( nf90_inquire_attribute(ncid, ncv, 'units', len=n_st_unit_length) )
call check( nf90_get_att(ncid, ncv, 'units', ch_st_unit) )

call check( nf90_inq_varid(ncid, 'SMB_clim', ncv) )
call check( nf90_get_var(ncid, ncv, smb_climatol_conv) )
call check( nf90_inquire_attribute(ncid, ncv, 'units', len=n_smb_unit_length) )
call check( nf90_get_att(ncid, ncv, 'units', ch_smb_unit) )

ios = nf90_inq_varid(ncid, 'surf_elev_ref', ncv)

if (ios == nf90_noerr) then
   call check( nf90_get_var(ncid, ncv, zs_ref_conv) )
else
   do i=0, IMAX
   do j=0, JMAX
      zs_ref_conv(i,j) = zs_present(j,i)
                         ! use previously read present-day topography
   end do
   end do
end if

call check( nf90_close(ncid) )

do i=0, IMAX
do j=0, JMAX

   if (trim(adjustl(ch_st_unit))=='K') then
      temp_maat_climatol(j,i) = temp_maat_climatol_conv(i,j) - temp_C_to_K
                                       ! K -> degC
   else if (trim(adjustl(ch_st_unit))=='degC') then
      temp_maat_climatol(j,i) = temp_maat_climatol_conv(i,j)
                                       ! degC
   else
      errormsg = ' >>> sico_init: Unit of ST_clim could not be determined!'
      call error(errormsg)
   end if

   if (trim(adjustl(ch_smb_unit))=='kg m-2 s-1') then
      smb_climatol(j,i) = smb_climatol_conv(i,j) /RHO
                                       ! kg/(m2*s) -> m/s ice equiv.
   else if (trim(adjustl(ch_smb_unit))=='m a-1') then
      smb_climatol(j,i) = smb_climatol_conv(i,j) *sec2year
                                       ! m/a ice equiv. -> m/s ice equiv.
   else
      errormsg = ' >>> sico_init: Unit of SMB_clim could not be determined!'
      call error(errormsg)
   end if

   zs_ref_climatol(j,i) = zs_ref_conv(i,j)

end do
end do

#endif

!-------- Prescribed surface mass balance correction --------

smb_corr_in = 0.0_dp

#if (defined(SMB_CORR_FILE))

if ( (trim(adjustl(SMB_CORR_FILE)) /= 'none') &
     .and. &
     (trim(adjustl(SMB_CORR_FILE)) /= 'None') &
     .and. &
     (trim(adjustl(SMB_CORR_FILE)) /= 'NONE') ) then

   filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                        trim(SMB_CORR_FILE)

   call read_2d_input(filename_with_path, &
                      ch_var_name='DSMB', n_var_type=1, n_ascii_header=6, &
                      field2d_r=field2d_aux)

   smb_corr_in = field2d_aux *sec2year
                             ! m/a ice equiv. -> m/s ice equiv.

end if

#endif

!-------- Reading of ISMIP6 SMB and BMB anomaly data --------

!  ------ Antarctica or Greenland: SMB (InitMIP)

#if (defined(ANT) || defined(GRL))

if (flag_initmip_asmb) then

   filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)// &
                                      '/'//trim(ch_initmip_smb_anom_file)

   ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

   if (ios /= nf90_noerr) then
      errormsg = ' >>> sico_init: Error when opening the file' &
               //                 end_of_line &
               //'                for the ISMIP6 SMB anomaly data!'
      call error(errormsg)
   end if

#if (defined(GRL))
   call check( nf90_inq_varid(ncid, 'DSMB', ncv) )
#else
   call check( nf90_inq_varid(ncid, 'asmb', ncv) )
#endif

   call check( nf90_get_var(ncid, ncv, field2d_tra_aux) )

   call check( nf90_close(ncid) )

   do i=0, IMAX
   do j=0, JMAX
      smb_anom_initmip(j,i) = field2d_tra_aux(i,j) *sec2year
                                   ! m/a ice equiv. -> m/s ice equiv.
   end do
   end do

else
   smb_anom_initmip = 0.0_dp
end if

#endif

!  ------ Antarctica only: BMB (InitMIP)

#if (defined(ANT))

if (flag_initmip_abmb) then

   filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)// &
                                      '/'//trim(ch_initmip_bmb_anom_file)

   ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

   if (ios /= nf90_noerr) then
      errormsg = ' >>> sico_init: Error when opening the file' &
               //                 end_of_line &
               //'                for the ISMIP6 BMB anomaly data!'
      call error(errormsg)
   end if

   call check( nf90_inq_varid(ncid, 'abmb', ncv) )
   call check( nf90_get_var(ncid, ncv, field2d_tra_aux) )

   call check( nf90_close(ncid) )

   do i=0, IMAX
   do j=0, JMAX
      ab_anom_initmip(j,i) = field2d_tra_aux(i,j) *sec2year
                                   ! m/a ice equiv. -> m/s ice equiv.
   end do
   end do

else
   ab_anom_initmip = 0.0_dp
end if

#endif

!  ------ Antarctica only: BMB (LARMIP)

#if (defined(ANT))

if (flag_larmip) then

   filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)// &
                                      '/'//trim(ch_larmip_regions_file)

   ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

   if (ios /= nf90_noerr) then
      errormsg = ' >>> sico_init: Error when opening the file' &
               //                 end_of_line &
               //'                for the LARMIP BMB regions!'
      call error(errormsg)
   end if

   call check( nf90_inq_varid(ncid, 'regions', ncv) )
   call check( nf90_get_var(ncid, ncv, field2d_tra_aux) )

   call check( nf90_close(ncid) )

   do i=0, IMAX
   do j=0, JMAX
      n_larmip_region(j,i) = nint(field2d_tra_aux(i,j))
   end do
   end do

   ab_anom_larmip      = 0.0_dp
   ab_anom_larmip(1:5) = larmip_qbm_anom_aux *sec2year
                                   ! m/a ice equiv. -> m/s ice equiv.
else
   n_larmip_region = 0
   ab_anom_larmip  = 0.0_dp
end if

#endif

!-------- Greenland only:
!         Reading of the global annual temperature anomaly
!         (for parameterizing the sub-ocean temperature anomaly
!         for the ice discharge parameterization)

#if (defined(GRL) && DISC==2)

filename_with_path = trim(IN_PATH)//'/general/dTg_paleo.dat'

call read_scalar_input(filename_with_path, &
                       'dT_glann', ndata_glann_max, &
                       glann_time_min, glann_time_stp, glann_time_max, &
                       ndata_glann, dT_glann_CLIMBER)

#endif

!-------- Read data for z_sl --------

#if (SEA_LEVEL==3)

filename_with_path = trim(IN_PATH)//'/general/'//trim(SEA_LEVEL_FILE)

call read_scalar_input(filename_with_path, &
                       'z_sl', ndata_specmap_max, &
                       specmap_time_min, specmap_time_stp, specmap_time_max, &
                       ndata_specmap, specmap_zsl)

#endif

!-------- Determination of the geothermal heat flux --------

if (n_q_geo_mod==1) then

!  ------ Constant value

   do i=0, IMAX
   do j=0, JMAX
      q_geo(j,i) = Q_GEO *1.0e-03_dp   ! mW/m2 -> W/m2
   end do
   end do

else if (n_q_geo_mod==2) then

!  ------ Read data from file

   filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                        trim(Q_GEO_FILE)

   call read_2d_input(filename_with_path, &
                      ch_var_name='GHF', n_var_type=1, n_ascii_header=6, &
                      field2d_r=field2d_aux)

   do i=0, IMAX
   do j=0, JMAX
      q_geo(j,i) = field2d_aux(j,i) *1.0e-03_dp   ! mW/m2 -> W/m2
   end do
   end do

end if

!-------- Reading of tabulated kei function--------

#if (REBOUND==0 || REBOUND==1)

kei = 0.0_dp; n_data_kei = 0; kei_r_max = 0.0_dp; kei_r_incr = 0.0_dp
                                                    ! dummy values
#elif (REBOUND==2)

call read_kei()

#endif

!-------- Determination of the time lag
!                              of the relaxing asthenosphere --------

#if (REBOUND==1 || REBOUND==2)

#if (TIME_LAG_MOD==1)

time_lag_asth = TIME_LAG*year2sec   ! a -> s

#elif (TIME_LAG_MOD==2)

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(TIME_LAG_FILE)

call read_2d_input(filename_with_path, &
                   ch_var_name='time_lag_asth', &
                   n_var_type=1, n_ascii_header=6, &
                   field2d_r=field2d_aux)

time_lag_asth = field2d_aux *year2sec   ! a -> s

#endif

#elif (REBOUND==0)

time_lag_asth = 0.0_dp   ! dummy values

#endif

!-------- Determination of the flexural rigidity of the lithosphere --------

#if (REBOUND==2)

#if (FLEX_RIG_MOD==1)

flex_rig_lith = FLEX_RIG

#elif (FLEX_RIG_MOD==2)

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(FLEX_RIG_FILE)

call read_2d_input(filename_with_path, &
                   ch_var_name='flex_rig_lith', &
                   n_var_type=1, n_ascii_header=6, &
                   field2d_r=field2d_aux)

flex_rig_lith = field2d_aux

#endif

#elif (REBOUND==0 || REBOUND==1)

flex_rig_lith = 0.0_dp   ! dummy values

#endif

!-------- Antarctica only:
!         Reading of the reference ice thickness for the
!                                  ice-shelf collapse masks --------

#if (defined(ANT) && ICE_SHELF_COLLAPSE_MASK==1)

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)// &
                                   '/'//trim(ICE_SHELF_COLLAPSE_MASK_H_REF_FILE)

ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

if (ios /= nf90_noerr) then
   errormsg = ' >>> sico_init: Error when opening the file' &
            //                 end_of_line &
            //'                for the reference ice thickness' &
            //                 end_of_line &
            //'                for the ice-shelf collapse masks!'
   call error(errormsg)
end if

call check( nf90_inq_varid(ncid, 'H', ncv) )
call check( nf90_get_var(ncid, ncv, H_ref_retreat_conv) )

call check( nf90_close(ncid) )

do i=0, IMAX
do j=0, JMAX
   H_ref_retreat(j,i) = max(H_ref_retreat_conv(i,j), 0.0_dp)
end do
end do

#endif

!-------- Greenland only:
!         Reading of the reference ice thickness for the retreat masks --------

#if (defined(GRL) && RETREAT_MASK==1)

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)// &
                                   '/'//trim(RETREAT_MASK_H_REF_FILE)

ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

if (ios /= nf90_noerr) then
   errormsg = ' >>> sico_init: Error when opening the file' &
            //                 end_of_line &
            //'                for the reference ice thickness' &
            //                 end_of_line &
            //'                for the retreat masks!'
   call error(errormsg)
end if

call check( nf90_inq_varid(ncid, 'H', ncv) )
call check( nf90_get_var(ncid, ncv, H_ref_retreat_conv) )

call check( nf90_close(ncid) )

do i=0, IMAX
do j=0, JMAX
   H_ref_retreat(j,i) = max(H_ref_retreat_conv(i,j), 0.0_dp)
end do
end do

#endif

!-------- Definition of initial values --------

!  ------ Present topography

#if (ANF_DAT==1)

call topography1(dxi, deta)

#if (defined(GRL) && DISC>0) /* Ice discharge parameterization for Greenland */
call disc_param(dtime)
call disc_fields()
#endif

z_sl      = -1.11e+11_dp   ! dummy values for initial call
z_sl_mean = -1.11e+11_dp   ! of subroutine boundary

call boundary(time_init, dtime, dxi, deta, &
              delta_ts, glac_index, z_mar)

do i=0, IMAX
do j=0, JMAX

   if ((mask(j,i)==0).or.(mask(j,i)==3)) then
                 ! grounded or floating ice
      as_perp_apl(j,i) = as_perp(j,i)
   else          ! mask==1 or 2, ice-free land or sea
      as_perp_apl(j,i) = 0.0_dp
   end if

end do
end do

smb_corr = 0.0_dp

Q_bm    = 0.0_dp
Q_tld   = 0.0_dp
Q_b_tot = 0.0_dp 

p_b_w   = 0.0_dp
H_w     = 0.0_dp

#if (TEMP_INIT==1)
  call init_temp_water_age_1_1()
#elif (TEMP_INIT==2)
  call init_temp_water_age_1_2()
#elif (TEMP_INIT==3)
  call init_temp_water_age_1_3()
#elif (TEMP_INIT==4)
  call init_temp_water_age_1_4()
#elif (TEMP_INIT==5)
  call init_temp_water_age_1_5(anfdatname)
#else
  errormsg = ' >>> sico_init: TEMP_INIT must be between 1 and 5!'
  call error(errormsg)
#endif

#if (defined(ALLOW_NORMAL) || defined(ALLOW_GRDCHK) || defined(ALLOW_TAPENADE))
#if defined(AGE_COST)
  call read_age_data()
#endif
#if defined(BEDMACHINE_COST)
  call read_BedMachine_data()
#endif
#endif

#if (ENHMOD==1)
   call calc_enhance_1()
#elif (ENHMOD==2)
   call calc_enhance_2()
#elif (ENHMOD==3)
   call calc_enhance_3(time_init)
#elif (ENHMOD==4)
   call calc_enhance_4()
#elif (ENHMOD==5)
   call calc_enhance_5()
#else
   errormsg = ' >>> sico_init: Parameter ENHMOD must be between 1 and 5!'
   call error(errormsg)
#endif

vx_m_sia = 0.0_dp
vy_m_sia = 0.0_dp
vx_m_ssa = 0.0_dp
vy_m_ssa = 0.0_dp

#if (defined(VISC_INIT_SSA))
  vis_ave_g = VISC_INIT_SSA
#else
  vis_ave_g = 1.0e+15_dp   ! Pa s
#endif

vis_int_g = 0.0_dp

!  ------ Ice-free, relaxed bedrock

#elif (ANF_DAT==2)

call topography2(dxi, deta)

#if (defined(GRL) && DISC>0) /* Ice discharge parameterization for Greenland */
call disc_param(dtime)
call disc_fields()
#endif

z_sl      = -1.11e+11_dp   ! dummy values for initial call
z_sl_mean = -1.11e+11_dp   ! of subroutine boundary

call boundary(time_init, dtime, dxi, deta, &
              delta_ts, glac_index, z_mar)

as_perp_apl = 0.0_dp

smb_corr = 0.0_dp

Q_bm    = 0.0_dp
Q_tld   = 0.0_dp
Q_b_tot = 0.0_dp 

p_b_w   = 0.0_dp
H_w     = 0.0_dp

call init_temp_water_age_2()

#if (ENHMOD==1)
   call calc_enhance_1()
#elif (ENHMOD==2)
   call calc_enhance_2()
#elif (ENHMOD==3)
   call calc_enhance_3(time_init)
#elif (ENHMOD==4)
   call calc_enhance_4()
#elif (ENHMOD==5)
   call calc_enhance_5()
#else
   errormsg = ' >>> sico_init: Parameter ENHMOD must be between 1 and 5!'
   call error(errormsg)
#endif

vx_m_sia = 0.0_dp
vy_m_sia = 0.0_dp
vx_m_ssa = 0.0_dp
vy_m_ssa = 0.0_dp

#if (defined(VISC_MIN) && defined(VISC_MAX))
  vis_ave_g = VISC_MAX
#else
  vis_ave_g = 1.0e+25_dp   ! Pa s
#endif

vis_int_g = 0.0_dp

!  ------ Read initial state from output of previous simulation

#elif (ANF_DAT==3)

call topography3(dxi, deta, anfdatname)

#if (defined(GRL) && DISC>0) /* Ice discharge parameterization for Greenland */
call disc_param(dtime)
call disc_fields()
#endif

#if (!(ANF_DAT==3 && RESTART==1))

call boundary(time_init, dtime, dxi, deta, &
              delta_ts, glac_index, z_mar)

do i=0, IMAX
do j=0, JMAX

   if ((mask(j,i)==0).or.(mask(j,i)==3)) then
                 ! grounded or floating ice
      as_perp_apl(j,i) = as_perp(j,i)
   else          ! mask==1 or 2, ice-free land or sea
      as_perp_apl(j,i) = 0.0_dp
   end if

end do
end do

smb_corr = 0.0_dp

#endif /* !(ANF_DAT==3 && RESTART==1) */

Q_b_tot = Q_bm + Q_tld

#if (!(ANF_DAT==3 && RESTART==1))

#if (ENHMOD==1)
   call calc_enhance_1()
#elif (ENHMOD==2)
   call calc_enhance_2()
#elif (ENHMOD==3)
   call calc_enhance_3(time_init)
#elif (ENHMOD==4)
   !%% call calc_enhance_4()
   !%%    (for anisotropic flow enhancement factor,
   !%%     use values read from output of previous simulation)
#elif (ENHMOD==5)
   !%% call calc_enhance_5()
   !%%    (for anisotropic flow enhancement factor,
   !%%     use values read from output of previous simulation)
#else
   errormsg = ' >>> sico_init: Parameter ENHMOD must be between 1 and 5!'
   call error(errormsg)
#endif

#endif /* !(ANF_DAT==3 && RESTART==1) */

#endif

!-------- Inner-point and staggered-grid flags --------

flag_inner_point            = .true.
flag_inner_point(0   ,:   ) = .false.
flag_inner_point(JMAX,:   ) = .false.
flag_inner_point(:   ,0   ) = .false.
flag_inner_point(:   ,IMAX) = .false.

flag_inner_inner_point                = flag_inner_point
flag_inner_inner_point(1     ,:     ) = .false.
flag_inner_inner_point(JMAX-1,:     ) = .false.
flag_inner_inner_point(:     ,1     ) = .false.
flag_inner_inner_point(:     ,IMAX-1) = .false.

flag_sg_x         = .true.
flag_sg_x(:,IMAX) = .false.

flag_sg_y         = .true.
flag_sg_y(JMAX,:) = .false.

flag_sg_x_inner_y         = flag_sg_x
flag_sg_x_inner_y(0   ,:) = .false.
flag_sg_x_inner_y(JMAX,:) = .false.

flag_sg_y_inner_x         = flag_sg_y
flag_sg_y_inner_x(:,0   ) = .false.
flag_sg_y_inner_x(:,IMAX) = .false.

!-------- Distance between grid points with delta_i=ir, delta_j=jr --------

#if (GRID==0 || GRID==1)

do ir=-IMAX, IMAX
do jr=-JMAX, JMAX
   dist_dxdy(jr,ir) = sqrt( (real(ir,dp)*dxi)**2 + (real(jr,dp)*deta)**2 )
                  ! distortion due to stereographic projection not accounted for
end do
end do

#elif (GRID==2)

do ir=-IMAX, IMAX
do jr=-JMAX, JMAX

   dist_dxdy(jr,ir) = sqrt( (sq_g11_g(JMAX/2,IMAX/2)*real(ir,dp)*dxi)**2 &
                          + (sq_g22_g(JMAX/2,IMAX/2)*real(jr,dp)*deta)**2 )

                      ! This uses the metric tensor in the center of the domain
		      ! for the entire domain; quite DIRTY TRICK!

end do
end do

#endif

!-------- Initial velocities --------

call calc_temp_melt()
call flag_update_gf_gl_cf()
call calc_vxy_b_init()
call calc_dzs_dxy_aux(dxi, deta)

#if (!(ANF_DAT==3 && RESTART==1))

#if (DYNAMICS==1 || DYNAMICS==2)

call calc_vxy_b_sia(time)
call calc_vxy_sia(dzeta_c, dzeta_t)

#if (MARGIN==3 || DYNAMICS==2)
call calc_vxy_ssa(dxi, deta, dzeta_c, dzeta_t)
#endif

call calc_vz_grounded(dxi, deta, dzeta_c, dzeta_t)

#if (MARGIN==3)
call calc_vz_floating(dxi, deta, dzeta_c)
#endif

#elif (DYNAMICS==0)

call calc_vxy_static()
call calc_vz_static()

#else
errormsg = ' >>> sico_init: DYNAMICS must be either 0, 1 or 2!'
call error(errormsg)
#endif

#endif /* !(ANF_DAT==3 && RESTART==1) */

call calc_dxyz(dxi, deta, dzeta_c, dzeta_t)

!-------- Initial enthalpies --------

#if (CALCMOD==0 || CALCMOD==-1)

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

#elif (CALCMOD==1)

do i=0, IMAX
do j=0, JMAX

   do kc=0, KCMAX
      enth_c(kc,j,i) = enth_fct_temp_omega(temp_c(kc,j,i), 0.0_dp)
   end do

   if ( (mask(j,i) == 0).and.(n_cts(j,i) == 1) ) then
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

#elif (CALCMOD==2 || CALCMOD==3)

do i=0, IMAX
do j=0, JMAX

   do kc=0, KCMAX
      enth_c(kc,j,i) = enth_fct_temp_omega(temp_c(kc,j,i), omega_c(kc,j,i))
   end do

   do kt=0, KTMAX
      enth_t(kt,j,i) = enth_c(0,j,i)
   end do

end do
end do

#else

errormsg = ' >>> sico_init: Parameter CALCMOD must be either -1, 0, 1, 2 or 3!'
call error(errormsg)

#endif

!-------- Initialize time-series files --------

!  ------ Time-series file for the entire ice sheet

filename_with_path = trim(OUT_PATH)//'/'//trim(run_name)//'.ser'

#if !defined(ALLOW_TAPENADE) /* Normal */
open(12, iostat=ios, file=trim(filename_with_path), status='new')
#else /* Tapenade */
open(12, iostat=ios, file=trim(filename_with_path))
#endif /* Normal vs. Tapenade */

if (ios /= 0) then
   errormsg = ' >>> sico_init: Error when opening the ser file!'
   call error(errormsg)
end if

if ((forcing_flag == 1).or.(forcing_flag == 3)) then

   write(12,1102)
   write(12,1103)

   1102 format('         t(a)   D_Ts(degC) z_sl_mean(m)',/, &
               '                    V(m^3)     V_g(m^3)     V_f(m^3)', &
               '       A(m^2)     A_g(m^2)     A_f(m^2)',/, &
               '                               V_sle(m)     V_t(m^3)', &
               '     A_t(m^2)',/, &
               '                               H_max(m)   H_t_max(m)', &
               '    zs_max(m)  vs_max(m/a)   Tbh_max(C)')
   1103 format('----------------------------------------------------', &
               '---------------------------------------')

else if (forcing_flag == 2) then

   write(12,1112)
   write(12,1113)

   1112 format('         t(a)  glac_ind(1) z_sl_mean(m)',/, &
               '                    V(m^3)     V_g(m^3)     V_f(m^3)', &
               '       A(m^2)     A_g(m^2)     A_f(m^2)',/, &
               '                               V_sle(m)     V_t(m^3)', &
               '     A_t(m^2)',/, &
               '                               H_max(m)   H_t_max(m)', &
               '    zs_max(m)  vs_max(m/a)   Tbh_max(C)')
   1113 format('----------------------------------------------------', &
               '---------------------------------------')

end if

!  ------ Time-series file for specified sites (i.e., ice cores)

#if (defined(ANT))
n_site = 6   ! Vostok, Dome A, Dome C, Dome F, Kohnen, Byrd
#elif (defined(GRL))
n_site = 7   ! GRIP, GISP2, Dye3, Camp Century (CC),
             ! NorthGRIP (NGRIP), NEEM, EastGRIP (EGRIP)
#else
n_site = 0   ! No sites defined
#endif

if (n_site > n_site_max) then
   errormsg = ' >>> sico_init: n_site <= n_site_max required!' &
            //         end_of_line &
            //'        Increase value of n_site_max in sico_variables_m!'
   call error(errormsg)
end if

#if (defined(ANT))

ch_site(1)     = 'Vostok'
phi_site(1)    = -78.467_dp *deg2rad    ! Geographical position of Vostok,
lambda_site(1) = 106.8_dp   *deg2rad    ! conversion deg -> rad
 
ch_site(2)     = 'Dome A'
phi_site(2)    = -80.367_dp *deg2rad    ! Geographical position of Dome A,
lambda_site(2) =  77.35_dp  *deg2rad    ! conversion deg -> rad

ch_site(3)     = 'Dome C'
phi_site(3)    = -75.1_dp *deg2rad      ! Geographical position of Dome C,
lambda_site(3) = 123.4_dp *deg2rad      ! conversion deg -> rad

ch_site(4)     = 'Dome F'
phi_site(4)    = -77.317_dp *deg2rad    ! Geographical position of Dome F,
lambda_site(4) =  39.7_dp   *deg2rad    ! conversion deg -> rad

ch_site(5)     = 'Kohnen'
phi_site(5)    = -75.0_dp   *deg2rad    ! Geographical position of Kohnen,
lambda_site(5) =   0.067_dp *deg2rad    ! conversion deg -> rad

ch_site(6)     = 'Byrd'
phi_site(6)    =  -80.017_dp *deg2rad   ! Geographical position of Byrd,
lambda_site(6) = -119.517_dp *deg2rad   ! conversion deg -> rad

#elif (defined(GRL))

ch_site(1)     = 'GRIP'
phi_site(1)    =  72.58722_dp *deg2rad   ! Geographical position of GRIP,
lambda_site(1) = -37.64222_dp *deg2rad   ! conversion deg -> rad
 
ch_site(2)     = 'GISP2'
phi_site(2)    =  72.58833_dp *deg2rad   ! Geographical position of GISP2
lambda_site(2) = -38.45750_dp *deg2rad   ! conversion deg -> rad

ch_site(3)     = 'Dye3'
phi_site(3)    =  65.15139_dp *deg2rad   ! Geographical position of Dye3,
lambda_site(3) = -43.81722_dp *deg2rad   ! conversion deg -> rad

ch_site(4)     = 'Camp Century'
phi_site(4)    =  77.17970_dp *deg2rad   ! Geographical position of CC,
lambda_site(4) = -61.10975_dp *deg2rad   ! conversion deg -> rad

ch_site(5)     = 'NGRIP'
phi_site(5)    =  75.09694_dp *deg2rad   ! Geographical position of NGRIP,
lambda_site(5) = -42.31956_dp *deg2rad   ! conversion deg -> rad

ch_site(6)     = 'NEEM'
phi_site(6)    =  77.5_dp     *deg2rad   ! Geographical position of NEEM,
lambda_site(6) = -50.9_dp     *deg2rad   ! conversion deg -> rad

ch_site(7)     = 'EGRIP'
phi_site(7)    =  75.6299_dp  *deg2rad   ! Geographical position of EGRIP,
lambda_site(7) = -35.9867_dp  *deg2rad   ! conversion deg -> rad

#endif

if (n_site > 0) then

#if (GRID==0 || GRID==1)   /* Stereographic projection */

   do n=1, n_site

      if (F_INV > 1.0e+10_dp) then   ! interpreted as infinity,
                                     ! thus no flattening (spherical planet)

         call stereo_forw_sphere(lambda_site(n), phi_site(n), &
                                 R, LAMBDA0, PHI0, x_site(n), y_site(n))

      else   ! finite inverse flattening (ellipsoidal planet)

         call stereo_forw_ellipsoid(lambda_site(n), phi_site(n), &
                                    A, B, LAMBDA0, PHI0, x_site(n), y_site(n))

      end if

   end do

#elif (GRID==2)   /* Geographical coordinates */

   x_site = lambda_site
   y_site = phi_site

#endif

end if

filename_with_path = trim(OUT_PATH)//'/'//trim(run_name)//'.site'

#if !defined(ALLOW_TAPENADE) /* Normal */
open(14, iostat=ios, file=trim(filename_with_path), status='new')
#else /* Tapenade */
open(14, iostat=ios, file=trim(filename_with_path))
#endif /* Normal vs. Tapenade */

if (ios /= 0) then
   errormsg = ' >>> sico_init: Error when opening the site file!'
   call error(errormsg)
end if

if (n_site == 0) then

   write(14,'(1x,a)') '-----------------'
   write(14,'(1x,a)') 'No sites defined.'
   write(14,'(1x,a)') '-----------------'

else

#if (defined(ANT))

   if ((forcing_flag == 1).or.(forcing_flag == 3)) then

      write(14,1106)
      write(14,1107)

      1106 format('         t(a)      D_Ts(C) z_sl_mean(m)',/, &
                  '                   H_Vo(m)      H_DA(m)      H_DC(m)', &
                  '      H_DF(m)      H_Ko(m)      H_By(m)',/, &
                  '                 v_Vo(m/a)    v_DA(m/a)    v_DC(m/a)', &
                  '    v_DF(m/a)    v_Ko(m/a)    v_By(m/a)',/, &
                  '                   T_Vo(C)      T_DA(C)      T_DC(C)', &
                  '      T_DF(C)      T_Ko(C)      T_By(C)')
      1107 format('----------------------------------------------------', &
                  '---------------------------------------')

   else if (forcing_flag == 2) then

      write(14,1116)
      write(14,1117)

      1116 format('         t(a)  glac_ind(1) z_sl_mean(m)',/, &
                  '                   H_Vo(m)      H_DA(m)      H_DC(m)', &
                  '      H_DF(m)      H_Ko(m)      H_By(m)',/, &
                  '                 v_Vo(m/a)    v_DA(m/a)    v_DC(m/a)', &
                  '    v_DF(m/a)    v_Ko(m/a)    v_By(m/a)',/, &
                  '                   T_Vo(C)      T_DA(C)      T_DC(C)', &
                  '      T_DF(C)      T_Ko(C)      T_By(C)')
      1117 format('----------------------------------------------------', &
                  '---------------------------------------')

   end if

#elif (defined(GRL))

   if ((forcing_flag == 1).or.(forcing_flag == 3)) then

      write(14,1106)
      write(14,1107)

      1106 format('         t(a)      D_Ts(C) z_sl_mean(m)',/, &
                  '                   H_GR(m)      H_G2(m)      H_D3(m)', &
                  '      H_CC(m)      H_NG(m)      H_NE(m)      H_EG(m)',/, &
                  '                 v_GR(m/a)    v_G2(m/a)    v_D3(m/a)', &
                  '    v_CC(m/a)    v_NG(m/a)    v_NE(m/a)    v_EG(m/a)',/, &
                  '                   T_GR(C)      T_G2(C)      T_D3(C)', &
                  '      T_CC(C)      T_NG(C)      T_NE(C)      T_EG(C)')
      1107 format('----------------------------------------------------', &
                  '----------------------------------------------------')

   else if (forcing_flag == 2) then

      write(14,1116)
      write(14,1117)

      1116 format('         t(a)  glac_ind(1) z_sl_mean(m)',/, &
                  '                   H_GR(m)      H_G2(m)      H_D3(m)', &
                  '      H_CC(m)      H_NG(m)      H_NE(m)      H_EG(m)',/, &
                  '                 v_GR(m/a)    v_G2(m/a)    v_D3(m/a)', &
                  '    v_CC(m/a)    v_NG(m/a)    v_NE(m/a)    v_EG(m/a)',/, &
                  '                   T_GR(C)      T_G2(C)      T_D3(C)', &
                  '      T_CC(C)      T_NG(C)      T_NE(C)      T_EG(C)')
      1117 format('----------------------------------------------------', &
                  '----------------------------------------------------')

   end if

#endif

end if

!  ------ Time-series file for mass balance stakes etc.

#if (defined(ASF) && WRITE_SER_FILE_STAKES==1) /* Austfonna */

n_surf = 163   ! 19 mass balance stakes + 18 cores (Pinglot) 
               ! 10 points on flowlines (Duvebreen & B3)
	       ! 116 points along ELA

if (n_surf > n_surf_max) then
   errormsg = ' >>> sico_init: n_surf <= n_surf_max required!' &
            //         end_of_line &
            //'        Increase value of n_surf_max in sico_vars_m!'
   call error(errormsg)
end if

!    ---- Mass balance stakes

phi_surf(1)    =  79.8322_dp *deg2rad    ! Geographical position
lambda_surf(1) =  24.0043_dp *deg2rad    ! at 79.8322N, 24.0043E

phi_surf(2)    =  79.3613_dp *deg2rad    ! Geographical position
lambda_surf(2) =  23.5622_dp *deg2rad    ! at 79.3613N, 23.5622E

phi_surf(3)    =  79.4497_dp *deg2rad    ! Geographical position
lambda_surf(3) =  23.6620_dp *deg2rad    ! at 79.4497N, 23.6620E

phi_surf(4)    =  79.5388_dp *deg2rad    ! Geographical position
lambda_surf(4) =  23.7644_dp *deg2rad    ! at 79.5388N, 23.7644E

phi_surf(5)    =  79.6421_dp *deg2rad    ! Geographical position
lambda_surf(5) =  23.8834_dp *deg2rad    ! at 79.6421N, 23.8834E

phi_surf(6)    =  79.7302_dp *deg2rad    ! Geographical position
lambda_surf(6) =  23.9872_dp *deg2rad    ! at 79.7302N, 23.9872E

phi_surf(7)    =  79.5829_dp *deg2rad    ! Geographical position
lambda_surf(7) =  24.6716_dp *deg2rad    ! at 79.5829N, 24.6716E

phi_surf(8)    =  79.7171_dp *deg2rad    ! Geographical position
lambda_surf(8) =  22.1832_dp *deg2rad    ! at 79.7171N, 22.1832E

phi_surf(9)    =  79.7335_dp *deg2rad    ! Geographical position
lambda_surf(9) =  22.4169_dp *deg2rad    ! at 79.7335N, 22.4169E

phi_surf(10)    =  79.7519_dp *deg2rad   ! Geographical position
lambda_surf(10) =  22.6407_dp *deg2rad   ! at 79.7519N, 22.6407E

phi_surf(11)    =  79.7670_dp *deg2rad   ! Geographical position
lambda_surf(11) =  22.8271_dp *deg2rad   ! at 79.7670N, 22.8271E

phi_surf(12)    =  79.7827_dp *deg2rad   ! Geographical position
lambda_surf(12) =  23.1220_dp *deg2rad   ! at 79.7827N, 23.1220E

phi_surf(13)    =  79.5884_dp *deg2rad   ! Geographical position
lambda_surf(13) =  25.5511_dp *deg2rad   ! at 79.5884N, 25.5511E

phi_surf(14)    =  79.6363_dp *deg2rad   ! Geographical position
lambda_surf(14) =  25.3582_dp *deg2rad   ! at 79.6363N, 25.3582E

phi_surf(15)    =  80.0638_dp *deg2rad   ! Geographical position
lambda_surf(15) =  24.2605_dp *deg2rad   ! at 80.0638N, 24.2605E

phi_surf(16)    =  79.9426_dp *deg2rad   ! Geographical position
lambda_surf(16) =  24.2433_dp *deg2rad   ! at 79.9426N, 24.2433E

phi_surf(17)    =  79.8498_dp *deg2rad   ! Geographical position
lambda_surf(17) =  26.6449_dp *deg2rad   ! at 79.8498N, 26.6449E

phi_surf(18)    =  79.8499_dp *deg2rad   ! Geographical position
lambda_surf(18) =  26.1354_dp *deg2rad   ! at 79.8499N, 26.1354E

phi_surf(19)    =  79.8499_dp *deg2rad   ! Geographical position
lambda_surf(19) =  25.7261_dp *deg2rad   ! at 79.8499N, 25.7261E

!    ---- Pinglot's shallow cores

phi_surf(20)    =  79.833333_dp *deg2rad    ! Geographical position
lambda_surf(20) =  24.935833_dp *deg2rad    ! at 79.833333N, 24.935833E

phi_surf(21)    =  79.783333_dp *deg2rad    ! Geographical position
lambda_surf(21) =  25.400000_dp *deg2rad    ! at 79.783333N, 25.400000E

phi_surf(22)    =  79.750000_dp *deg2rad    ! Geographical position
lambda_surf(22) =  23.866667_dp *deg2rad    ! at 79.750000N, 23.866667E

phi_surf(23)    =  79.615000_dp *deg2rad    ! Geographical position
lambda_surf(23) =  23.490556_dp *deg2rad    ! at 79.615000N, 23.490556E

phi_surf(24)    =  79.797778_dp *deg2rad    ! Geographical position
lambda_surf(24) =  23.997500_dp *deg2rad    ! at 79.797778N, 23.997500E

phi_surf(25)    =  79.765000_dp *deg2rad    ! Geographical position
lambda_surf(25) =  24.809722_dp *deg2rad    ! at 79.765000N, 24.809722E

phi_surf(26)    =  79.874722_dp *deg2rad    ! Geographical position
lambda_surf(26) =  23.541667_dp *deg2rad    ! at 79.874722N, 23.541667E

phi_surf(27)    =  79.697778_dp *deg2rad    ! Geographical position
lambda_surf(27) =  25.096111_dp *deg2rad    ! at 79.697778N, 25.096111E

phi_surf(28)    =  79.897500_dp *deg2rad    ! Geographical position
lambda_surf(28) =  23.230278_dp *deg2rad    ! at 79.897500N, 23.230278E

phi_surf(29)    =  79.874444_dp *deg2rad    ! Geographical position
lambda_surf(29) =  24.046111_dp *deg2rad    ! at 79.874444N, 24.046111E

phi_surf(30)    =  79.962500_dp *deg2rad    ! Geographical position
lambda_surf(30) =  24.169722_dp *deg2rad    ! at 79.962500N, 24.169722E

phi_surf(31)    =  79.664444_dp *deg2rad    ! Geographical position
lambda_surf(31) =  25.235833_dp *deg2rad    ! at 79.664444N, 25.235833E

phi_surf(32)    =  79.681111_dp *deg2rad    ! Geographical position
lambda_surf(32) =  23.713056_dp *deg2rad    ! at 79.681111N, 23.713056E

phi_surf(33)    =  79.554167_dp *deg2rad    ! Geographical position
lambda_surf(33) =  23.796944_dp *deg2rad    ! at 79.554167N, 23.796944E

phi_surf(34)    =  79.511667_dp *deg2rad    ! Geographical position
lambda_surf(34) =  24.032778_dp *deg2rad    ! at 79.511667N, 24.032778E

phi_surf(35)    =  79.552222_dp *deg2rad    ! Geographical position
lambda_surf(35) =  22.799167_dp *deg2rad    ! at 79.552222N, 22.799167E

phi_surf(36)    =  79.847778_dp *deg2rad    ! Geographical position
lambda_surf(36) =  25.788611_dp *deg2rad    ! at 79.847778N, 25.788611E

phi_surf(37)    =  79.830000_dp *deg2rad    ! Geographical position
lambda_surf(37) =  24.001389_dp *deg2rad    ! at 79.830000N, 24.001389E

!    ---- Flowline points

phi_surf(38)    =  80.1427268586056_dp *deg2rad    ! Geographical position of
lambda_surf(38) =  23.9534492294493_dp *deg2rad    ! Duve-1

phi_surf(39)    =  80.1124108950185_dp *deg2rad    ! Geographical position of
lambda_surf(39) =  24.0629399381155_dp *deg2rad    ! Duve-2

phi_surf(40)    =  80.0765637664780_dp *deg2rad    ! Geographical position of
lambda_surf(40) =  24.0481674197519_dp *deg2rad    ! Duve-3

phi_surf(41)    =  80.0409891299823_dp *deg2rad    ! Geographical position of
lambda_surf(41) =  24.0074110976581_dp *deg2rad    ! Duve-4

phi_surf(42)    =  80.0049393359201_dp *deg2rad    ! Geographical position of
lambda_surf(42) =  23.9894145095442_dp *deg2rad    ! Duve-5

phi_surf(43)    =  79.4994665039616_dp *deg2rad    ! Geographical position of
lambda_surf(43) =  25.4790616341716_dp *deg2rad    ! B3-1

phi_surf(44)    =  79.4973763443781_dp *deg2rad    ! Geographical position of
lambda_surf(44) =  25.2826485444194_dp *deg2rad    ! B3-2

phi_surf(45)    =  79.5028080484427_dp *deg2rad    ! Geographical position of
lambda_surf(45) =  25.0844021770897_dp *deg2rad    ! B3-3

phi_surf(46)    =  79.5131051861579_dp *deg2rad    ! Geographical position of
lambda_surf(46) =  24.8934874838598_dp *deg2rad    ! B3-4

phi_surf(47)    =  79.5275754154375_dp *deg2rad    ! Geographical position of
lambda_surf(47) =  24.7125320718015_dp *deg2rad    ! B3-5

!    ---- Basin control points on ELA (N:450m, S:300m)

phi_surf(48)    =  79.6232572730302_dp *deg2rad    ! Geographical position of
lambda_surf(48) =  22.4297425686265_dp *deg2rad    ! Eton-1

phi_surf(49)    =  79.6355048663362_dp *deg2rad    ! Geographical position of
lambda_surf(49) =  22.5023513660534_dp *deg2rad    ! Eton-2

phi_surf(50)    =  79.6477359930900_dp *deg2rad    ! Geographical position of
lambda_surf(50) =  22.5751300038166_dp *deg2rad    ! Eton-3

phi_surf(51)    =  79.6599505942585_dp *deg2rad    ! Geographical position of
lambda_surf(51) =  22.6480788556811_dp *deg2rad    ! Eton-4

phi_surf(52)    =  79.6730674725108_dp *deg2rad    ! Geographical position of
lambda_surf(52) =  22.7116449352996_dp *deg2rad    ! Eton-5

phi_surf(53)    =  79.6907455504277_dp *deg2rad    ! Geographical position of
lambda_surf(53) =  22.7278148586532_dp *deg2rad    ! Eton-6

phi_surf(54)    =  79.7084227767215_dp *deg2rad    ! Geographical position of
lambda_surf(54) =  22.7440404597164_dp *deg2rad    ! Eton-7

phi_surf(55)    =  79.7260991471427_dp *deg2rad    ! Geographical position of
lambda_surf(55) =  22.7603220234687_dp *deg2rad    ! Eton-8

phi_surf(56)    =  79.7437746574126_dp *deg2rad    ! Geographical position of
lambda_surf(56) =  22.7766598368173_dp *deg2rad    ! Eton-9

phi_surf(57)    =  79.7615003936967_dp *deg2rad    ! Geographical position of
lambda_surf(57) =  22.7895141723757_dp *deg2rad    ! Eton-10

phi_surf(58)    =  79.7794141201101_dp *deg2rad    ! Geographical position of
lambda_surf(58) =  22.7893415392149_dp *deg2rad    ! B-16s-1

phi_surf(59)    =  79.7973278451020_dp *deg2rad    ! Geographical position of
lambda_surf(59) =  22.7891690597211_dp *deg2rad    ! B-16s-2

phi_surf(60)    =  79.8152415686728_dp *deg2rad    ! Geographical position of
lambda_surf(60) =  22.7889967333372_dp *deg2rad    ! B-16s-3

phi_surf(61)    =  79.8331552908230_dp *deg2rad    ! Geographical position of
lambda_surf(61) =  22.7888245595023_dp *deg2rad    ! B-16s-4

phi_surf(62)    =  79.8504448969531_dp *deg2rad    ! Geographical position of
lambda_surf(62) =  22.8027142916594_dp *deg2rad    ! B-16n-1

phi_surf(63)    =  79.8662041154283_dp *deg2rad    ! Geographical position of
lambda_surf(63) =  22.8510765245997_dp *deg2rad    ! B-16n-2

phi_surf(64)    =  79.8819561232071_dp *deg2rad    ! Geographical position of
lambda_surf(64) =  22.8995882814793_dp *deg2rad    ! B-16n-3

phi_surf(65)    =  79.8977008864609_dp *deg2rad    ! Geographical position of
lambda_surf(65) =  22.9482501953580_dp *deg2rad    ! B-16n-4

phi_surf(66)    =  79.9134383711667_dp *deg2rad    ! Geographical position of
lambda_surf(66) =  22.9970629023954_dp *deg2rad    ! B-16n-5

phi_surf(67)    =  79.9291685431056_dp *deg2rad    ! Geographical position of
lambda_surf(67) =  23.0460270418662_dp *deg2rad    ! B-16n-6

phi_surf(68)    =  79.9448913678619_dp *deg2rad    ! Geographical position of
lambda_surf(68) =  23.0951432561750_dp *deg2rad    ! B-16n-7

phi_surf(69)    =  79.9606068108212_dp *deg2rad    ! Geographical position of
lambda_surf(69) =  23.1444121908719_dp *deg2rad    ! B-16n-8

phi_surf(70)    =  79.9741572381786_dp *deg2rad    ! Geographical position of
lambda_surf(70) =  23.2092211687201_dp *deg2rad    ! Botnevika-1

phi_surf(71)    =  79.9859141894524_dp *deg2rad    ! Geographical position of
lambda_surf(71) =  23.2868821248161_dp *deg2rad    ! Botnevika-2

phi_surf(72)    =  79.9976529206869_dp *deg2rad    ! Geographical position of
lambda_surf(72) =  23.3647236505600_dp *deg2rad    ! Botnevika-3

phi_surf(73)    =  80.0093733670701_dp *deg2rad    ! Geographical position of
lambda_surf(73) =  23.4427461021207_dp *deg2rad    ! Botnevika-4

phi_surf(74)    =  80.0201320622880_dp *deg2rad    ! Geographical position of
lambda_surf(74) =  23.5253067161782_dp *deg2rad    ! Botnevika-5

phi_surf(75)    =  80.0308022109253_dp *deg2rad    ! Geographical position of
lambda_surf(75) =  23.6083570802514_dp *deg2rad    ! Botnevika-6

phi_surf(76)    =  80.0414516357850_dp *deg2rad    ! Geographical position of
lambda_surf(76) =  23.6915833394057_dp *deg2rad    ! Botnevika-7

phi_surf(77)    =  80.0520802696857_dp *deg2rad    ! Geographical position of
lambda_surf(77) =  23.7749857156754_dp *deg2rad    ! Botnevika-8

phi_surf(78)    =  80.0547633949370_dp *deg2rad    ! Geographical position of
lambda_surf(78) =  23.8736736708044_dp *deg2rad    ! Duvebreen-1

phi_surf(79)    =  80.0548013447126_dp *deg2rad    ! Geographical position of
lambda_surf(79) =  23.9773687987851_dp *deg2rad    ! Duvebreen-2

phi_surf(80)    =  80.0548073397268_dp *deg2rad    ! Geographical position of
lambda_surf(80) =  24.0810636270044_dp *deg2rad    ! Duvebreen-3

phi_surf(81)    =  80.0547813803758_dp *deg2rad    ! Geographical position of
lambda_surf(81) =  24.1847574925018_dp *deg2rad    ! Duvebreen-4

phi_surf(82)    =  80.0646160588300_dp *deg2rad    ! Geographical position of
lambda_surf(82) =  24.2700368789878_dp *deg2rad    ! Duvebreen-5

phi_surf(83)    =  80.0750999374003_dp *deg2rad    ! Geographical position of
lambda_surf(83) =  24.3542380951582_dp *deg2rad    ! Duvebreen-6

phi_surf(84)    =  80.0846920877530_dp *deg2rad    ! Geographical position of
lambda_surf(84) =  24.4407004402100_dp *deg2rad    ! Duvebreen-7

phi_surf(85)    =  80.0875193831616_dp *deg2rad    ! Geographical position of
lambda_surf(85) =  24.5434121380084_dp *deg2rad    ! Schweigaardbreen-1

phi_surf(86)    =  80.0903153574351_dp *deg2rad    ! Geographical position of
lambda_surf(86) =  24.6461808494348_dp *deg2rad    ! Schweigaardbreen-2

phi_surf(87)    =  80.0924166470023_dp *deg2rad    ! Geographical position of
lambda_surf(87) =  24.7486469216956_dp *deg2rad    ! Schweigaardbreen-3

phi_surf(88)    =  80.0864319373603_dp *deg2rad    ! Geographical position of
lambda_surf(88) =  24.8467147281595_dp *deg2rad    ! Schweigaardbreen-4

phi_surf(89)    =  80.0804188683848_dp *deg2rad    ! Geographical position of
lambda_surf(89) =  24.9446644540260_dp *deg2rad    ! Schweigaardbreen-5

phi_surf(90)    =  80.0743774931913_dp *deg2rad    ! Geographical position of
lambda_surf(90) =  25.0424957604751_dp *deg2rad    ! Schweigaardbreen-6

phi_surf(91)    =  80.0713340422000_dp *deg2rad    ! Geographical position of
lambda_surf(91) =  25.1439126047994_dp *deg2rad    ! Nilsenbreen B-12-1

phi_surf(92)    =  80.0700730909331_dp *deg2rad    ! Geographical position of
lambda_surf(92) =  25.2475056357563_dp *deg2rad    ! Nilsenbreen B-12-2

phi_surf(93)    =  80.0687803205250_dp *deg2rad    ! Geographical position of
lambda_surf(93) =  25.3510715226335_dp *deg2rad    ! Nilsenbreen B-12-3

phi_surf(94)    =  80.0647501708291_dp *deg2rad    ! Geographical position of
lambda_surf(94) =  25.4519066363393_dp *deg2rad    ! Sexebreen B-11-1

phi_surf(95)    =  80.0595181102431_dp *deg2rad    ! Geographical position of
lambda_surf(95) =  25.5506489732496_dp *deg2rad    ! Leighbreen-1

phi_surf(96)    =  80.0494857323914_dp *deg2rad    ! Geographical position of
lambda_surf(96) =  25.6365356440635_dp *deg2rad    ! Leighbreen-2

phi_surf(97)    =  80.0394316265850_dp *deg2rad    ! Geographical position of
lambda_surf(97) =  25.7222505219501_dp *deg2rad    ! Leighbreen-3

phi_surf(98)    =  80.0293558606091_dp *deg2rad    ! Geographical position of
lambda_surf(98) =  25.8077937609009_dp *deg2rad    ! Leighbreen-4

phi_surf(99)    =  80.0192585021221_dp *deg2rad    ! Geographical position of
lambda_surf(99) =  25.8931655175225_dp *deg2rad    ! Leighbreen-5

phi_surf(100)    =  80.0091396186553_dp *deg2rad    ! Geographical position of
lambda_surf(100) =  25.9783659510134_dp *deg2rad    ! Leighbreen-6

phi_surf(101)    =  79.9989992776120_dp *deg2rad    ! Geographical position of
lambda_surf(101) =  26.0633952231408_dp *deg2rad    ! Leighbreen-7

phi_surf(102)    =  79.9888375462661_dp *deg2rad    ! Geographical position of
lambda_surf(102) =  26.1482534982178_dp *deg2rad    ! Leighbreen-8

phi_surf(103)    =  79.9786544917617_dp *deg2rad    ! Geographical position of
lambda_surf(103) =  26.2329409430807_dp *deg2rad    ! Leighbreen-9

phi_surf(104)    =  79.9683923353960_dp *deg2rad    ! Geographical position of
lambda_surf(104) =  26.3172101192864_dp *deg2rad    ! Leighbreen-10

phi_surf(105)    =  80.0241705082505_dp *deg2rad    ! Geographical position of
lambda_surf(105) =  26.7558248932553_dp *deg2rad    ! Worsleybreen-1 (B9-1)

phi_surf(106)    =  80.0069243536208_dp *deg2rad    ! Geographical position of
lambda_surf(106) =  26.7836310921011_dp *deg2rad    ! Worsleybreen-2 (B9-2)

phi_surf(107)    =  79.9896760170551_dp *deg2rad    ! Geographical position of
lambda_surf(107) =  26.8113433337043_dp *deg2rad    ! Worsleybreen-3 (B9-3)

phi_surf(108)    =  79.9723667157507_dp *deg2rad    ! Geographical position of
lambda_surf(108) =  26.8350524380302_dp *deg2rad    ! Worsleybreen-4 (B9-4)

phi_surf(109)    =  79.9545472297622_dp *deg2rad    ! Geographical position of
lambda_surf(109) =  26.8248911276131_dp *deg2rad    ! B8-1

phi_surf(110)    =  79.9367274171506_dp *deg2rad    ! Geographical position of
lambda_surf(110) =  26.8147665774914_dp *deg2rad    ! B8-2

phi_surf(111)    =  79.9189072796258_dp *deg2rad    ! Geographical position of
lambda_surf(111) =  26.8046785944172_dp *deg2rad    ! B8-3

phi_surf(112)    =  79.9009446914988_dp *deg2rad    ! Geographical position of
lambda_surf(112) =  26.7957185084455_dp *deg2rad    ! B7-1

phi_surf(113)    =  79.8843576455373_dp *deg2rad    ! Geographical position of
lambda_surf(113) =  26.7616970403497_dp *deg2rad    ! B7-2

phi_surf(114)    =  79.8676428266616_dp *deg2rad    ! Geographical position of
lambda_surf(114) =  26.7251472990965_dp *deg2rad    ! B7-3

phi_surf(115)    =  79.8509238637717_dp *deg2rad    ! Geographical position of
lambda_surf(115) =  26.6887177159393_dp *deg2rad    ! B7-4

phi_surf(116)    =  79.8342007771708_dp *deg2rad    ! Geographical position of
lambda_surf(116) =  26.6524077251556_dp *deg2rad    ! B7-5

phi_surf(117)    =  79.8189961177120_dp *deg2rad    ! Geographical position of
lambda_surf(117) =  26.6017802396904_dp *deg2rad    ! B6-1

phi_surf(118)    =  79.8054200039019_dp *deg2rad    ! Geographical position of
lambda_surf(118) =  26.5357666498664_dp *deg2rad    ! B6-2

phi_surf(119)    =  79.7918304753589_dp *deg2rad    ! Geographical position of
lambda_surf(119) =  26.4699273874801_dp *deg2rad    ! B6-3

phi_surf(120)    =  79.7782275858515_dp *deg2rad    ! Geographical position of
lambda_surf(120) =  26.4042619219016_dp *deg2rad    ! B6-4

phi_surf(121)    =  79.7646113889145_dp *deg2rad    ! Geographical position of
lambda_surf(121) =  26.3387697236600_dp *deg2rad    ! B6-5

phi_surf(122)    =  79.7518386380187_dp *deg2rad    ! Geographical position of
lambda_surf(122) =  26.2683717557144_dp *deg2rad    ! B5-1

phi_surf(123)    =  79.7395107596368_dp *deg2rad    ! Geographical position of
lambda_surf(123) =  26.1954158840248_dp *deg2rad    ! B5-2

phi_surf(124)    =  79.7271664326874_dp *deg2rad    ! Geographical position of
lambda_surf(124) =  26.1226336416600_dp *deg2rad    ! B5-3

phi_surf(125)    =  79.7148057168060_dp *deg2rad    ! Geographical position of
lambda_surf(125) =  26.0500246274899_dp *deg2rad    ! B5-4

phi_surf(126)    =  79.7024286714212_dp *deg2rad    ! Geographical position of
lambda_surf(126) =  25.9775884402940_dp *deg2rad    ! B5-5

phi_surf(127)    =  79.6900353557545_dp *deg2rad    ! Geographical position of
lambda_surf(127) =  25.9053246787703_dp *deg2rad    ! B5-6

phi_surf(128)    =  79.6776258288211_dp *deg2rad    ! Geographical position of
lambda_surf(128) =  25.8332329415456_dp *deg2rad    ! B5-7

phi_surf(129)    =  79.6652001494302_dp *deg2rad    ! Geographical position of
lambda_surf(129) =  25.7613128271851_dp *deg2rad    ! B5-8

phi_surf(130)    =  79.6527583761852_dp *deg2rad    ! Geographical position of
lambda_surf(130) =  25.6895639342015_dp *deg2rad    ! B5-9

phi_surf(131)    =  79.6403005674845_dp *deg2rad    ! Geographical position of
lambda_surf(131) =  25.6179858610658_dp *deg2rad    ! B5-10

phi_surf(132)    =  79.6272788783125_dp *deg2rad    ! Geographical position of
lambda_surf(132) =  25.5497696493382_dp *deg2rad    ! B4-1

phi_surf(133)    =  79.6138476738577_dp *deg2rad    ! Geographical position of
lambda_surf(133) =  25.4840259325117_dp *deg2rad    ! B4-2

phi_surf(134)    =  79.6004029370116_dp *deg2rad    ! Geographical position of
lambda_surf(134) =  25.4184506246986_dp *deg2rad    ! B4-3

phi_surf(135)    =  79.5869447205062_dp *deg2rad    ! Geographical position of
lambda_surf(135) =  25.3530432378053_dp *deg2rad    ! B4-4

phi_surf(136)    =  79.5734730768545_dp *deg2rad    ! Geographical position of
lambda_surf(136) =  25.2878032846200_dp *deg2rad    ! B4-5

phi_surf(137)    =  79.5599880583521_dp *deg2rad    ! Geographical position of
lambda_surf(137) =  25.2227302788170_dp *deg2rad    ! B4-6

phi_surf(138)    =  79.5464897170775_dp *deg2rad    ! Geographical position of
lambda_surf(138) =  25.1578237349623_dp *deg2rad    ! B4-7

phi_surf(139)    =  79.5340825476013_dp *deg2rad    ! Geographical position of
lambda_surf(139) =  25.0873713598923_dp *deg2rad    ! B3-1

phi_surf(140)    =  79.5231871974923_dp *deg2rad    ! Geographical position of
lambda_surf(140) =  25.0091720580033_dp *deg2rad    ! B3-2

phi_surf(141)    =  79.5122726145574_dp *deg2rad    ! Geographical position of
lambda_surf(141) =  24.9311335486110_dp *deg2rad    ! B3-3

phi_surf(142)    =  79.5013388593293_dp *deg2rad    ! Geographical position of
lambda_surf(142) =  24.8532556096146_dp *deg2rad    ! B3-4

phi_surf(143)    =  79.4881304535468_dp *deg2rad    ! Geographical position of
lambda_surf(143) =  24.7885573077964_dp *deg2rad    ! B3-5

phi_surf(144)    =  79.4734132097634_dp *deg2rad    ! Geographical position of
lambda_surf(144) =  24.7326565135170_dp *deg2rad    ! B3-6

phi_surf(145)    =  79.4586860312332_dp *deg2rad    ! Geographical position of
lambda_surf(145) =  24.6769105936574_dp *deg2rad    ! B3-7

phi_surf(146)    =  79.4439489597131_dp *deg2rad    ! Geographical position of
lambda_surf(146) =  24.6213190006049_dp *deg2rad    ! B3-8

phi_surf(147)    =  79.4321693404700_dp *deg2rad    ! Geographical position of
lambda_surf(147) =  24.5500779464491_dp *deg2rad    ! B2-1

phi_surf(148)    =  79.4223453273505_dp *deg2rad    ! Geographical position of
lambda_surf(148) =  24.4684716320257_dp *deg2rad    ! B2-2

phi_surf(149)    =  79.4125002037095_dp *deg2rad    ! Geographical position of
lambda_surf(149) =  24.3870150299917_dp *deg2rad    ! B2-3

phi_surf(150)    =  79.4026340289842_dp *deg2rad    ! Geographical position of
lambda_surf(150) =  24.3057080421768_dp *deg2rad    ! B2-4

phi_surf(151)    =  79.3927468625203_dp *deg2rad    ! Geographical position of
lambda_surf(151) =  24.2245505685362_dp *deg2rad    ! B2-5

phi_surf(152)    =  79.3909641358607_dp *deg2rad    ! Geographical position of
lambda_surf(152) =  24.1356247611452_dp *deg2rad    ! Brasvellbreen-1

phi_surf(153)    =  79.3950618239069_dp *deg2rad    ! Geographical position of
lambda_surf(153) =  24.0409163942958_dp *deg2rad    ! Brasvellbreen-2

phi_surf(154)    =  79.3991312122811_dp *deg2rad    ! Geographical position of
lambda_surf(154) =  23.9461351693152_dp *deg2rad    ! Brasvellbreen-3

phi_surf(155)    =  79.4031722671433_dp *deg2rad    ! Geographical position of
lambda_surf(155) =  23.8512815066396_dp *deg2rad    ! Brasvellbreen-4

phi_surf(156)    =  79.4071849548373_dp *deg2rad    ! Geographical position of
lambda_surf(156) =  23.7563558291274_dp *deg2rad    ! Brasvellbreen-5

phi_surf(157)    =  79.4111692418918_dp *deg2rad    ! Geographical position of
lambda_surf(157) =  23.6613585620463_dp *deg2rad    ! Brasvellbreen-6

phi_surf(158)    =  79.4127149901435_dp *deg2rad    ! Geographical position of
lambda_surf(158) =  23.5647431017868_dp *deg2rad    ! Brasvellbreen-7

phi_surf(159)    =  79.4129320057492_dp *deg2rad    ! Geographical position of
lambda_surf(159) =  23.4672773246991_dp *deg2rad    ! Brasvellbreen-8

phi_surf(160)    =  79.4131190508990_dp *deg2rad    ! Geographical position of
lambda_surf(160) =  23.3698071241014_dp *deg2rad    ! Brasvellbreen-9

phi_surf(161)    =  79.4132761235192_dp *deg2rad    ! Geographical position of
lambda_surf(161) =  23.2723330506382_dp *deg2rad    ! Brasvellbreen-10

phi_surf(162)    =  79.4134032217989_dp *deg2rad    ! Geographical position of
lambda_surf(162) =  23.1748556552727_dp *deg2rad    ! Brasvellbreen-11

phi_surf(163)    =  79.4135003441905_dp *deg2rad    ! Geographical position of
lambda_surf(163) =  23.0773754892604_dp *deg2rad    ! Brasvellbreen-12

#if (GRID==0 || GRID==1) /* Stereographic projection */

do n=1, n_surf

   if (F_INV > 1.0e+10_dp) then   ! interpreted as infinity, thus no flattening
                                  ! (spherical planet)

      call stereo_forw_sphere(lambda_surf(n), phi_surf(n), &
                              R, LAMBDA0, PHI0, x_surf(n), y_surf(n))

   else   ! finite inverse flattening (ellipsoidal planet)

      call stereo_forw_ellipsoid(lambda_surf(n), phi_surf(n), &
                                 A, B, LAMBDA0, PHI0, x_surf(n), y_surf(n))

   end if

end do

#elif (GRID==2) /* Geographical coordinates */

x_surf = lambda_surf
y_surf = phi_surf

#endif

#endif /* (defined(ASF) && WRITE_SER_FILE_STAKES==1) Austfonna */

!-------- Output of the initial state --------

#if !defined(ALLOW_TAPENADE) /* Normal */

#if (defined(OUTPUT_INIT))

#if (OUTPUT_INIT==0)
   flag_init_output = .false.
#elif (OUTPUT_INIT==1)
   flag_init_output = .true.
#else
   errormsg = ' >>> sico_init: OUTPUT_INIT must be either 0 or 1!'
   call error(errormsg)
#endif

#else
   flag_init_output = .true.   ! default for undefined parameter OUTPUT_INIT
#endif

#if (OUTPUT==1)

#if (ERGDAT==0)
   flag_3d_output = .false.
#elif (ERGDAT==1)
   flag_3d_output = .true.
#else
   errormsg = ' >>> sico_init: ERGDAT must be either 0 or 1!'
   call error(errormsg)
#endif

   if (flag_init_output) &
      call output1(time_init, delta_ts, glac_index, &
                   flag_3d_output, ndat2d, ndat3d)

#elif (OUTPUT==2)

if (time_output(1) <= time_init+eps) then

#if (ERGDAT==0)
   flag_3d_output = .false.
#elif (ERGDAT==1)
   flag_3d_output = .true.
#else
   errormsg = ' >>> sico_init: ERGDAT must be either 0 or 1!'
   call error(errormsg)
#endif

   call output1(time_init, delta_ts, glac_index, &
                flag_3d_output, ndat2d, ndat3d)

end if

#elif (OUTPUT==3)

   flag_3d_output = .false.

   if (flag_init_output) &
      call output1(time_init, delta_ts, glac_index, &
                   flag_3d_output, ndat2d, ndat3d)

if (time_output(1) <= time_init+eps) then

   flag_3d_output = .true.

   call output1(time_init, delta_ts, glac_index, &
                flag_3d_output, ndat2d, ndat3d)

end if

#else
   errormsg = ' >>> sico_init: OUTPUT must be either 1, 2 or 3!'
   call error(errormsg)
#endif

if (flag_init_output) then

   call output2(time_init, dxi, deta, delta_ts, glac_index)
   call output4(time_init, dxi, deta, delta_ts, glac_index)

#if (defined(ASF) && WRITE_SER_FILE_STAKES==1) /* Austfonna */
   call output5(time, dxi, deta, delta_ts, glac_index)
#endif

end if

#else /* Tapenade */

print *, ' >>> sico_init: not producing initial, typical outputs'
print *, '                in adjoint mode.'

#endif /* Normal vs. Tapenade */

#if (defined(GRL) && defined(EXEC_MAKE_C_DIS_0))

#if (DISC>0)

call calc_c_dis_0(dxi, deta)

errormsg = ' >>> sico_init: Routine calc_c_dis_0 successfully completed,' &
         //         end_of_line &
         //'        c_dis_0 written on file out_run_name.dat' &
         //         end_of_line &
         //'        (in directory specified by OUT_PATH).' &
         //         end_of_line &
         //'        Execution of SICOPOLIS stopped.'
call error(errormsg)   ! actually not an error,
                       ! just a regular stop with an info message

#else
  errormsg = ' >>> sico_init: EXEC_MAKE_C_DIS_0 requires DISC>0!'
  call error(errormsg)
#endif

#endif

end subroutine sico_init

!-------------------------------------------------------------------------------
!> Definition of the initial surface and bedrock topography
!! (including gradients) and of the horizontal grid spacings dxi, deta.
!! For present-day initial topography.
!-------------------------------------------------------------------------------
subroutine topography1(dxi, deta)

  use read_m, only : read_2d_input

#if (GRID==0 || GRID==1)
  use stereo_proj_m
#endif

  use metric_m
  use topograd_m

implicit none

real(dp), intent(out) :: dxi, deta

integer(i4b) :: i, j, n
real(dp)     :: xi0, eta0
real(dp)     :: H_ice, freeboard_ratio

character(len=256) :: filename_with_path

real(dp), dimension(0:JMAX,0:IMAX) :: field2d_aux

!-------- Read topography --------

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(ZS_PRESENT_FILE)

call read_2d_input(filename_with_path, &
                   ch_var_name='zs', n_var_type=1, n_ascii_header=6, &
                   field2d_r=field2d_aux)

zs = field2d_aux

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(ZL_PRESENT_FILE)

call read_2d_input(filename_with_path, &
                   ch_var_name='zl', n_var_type=1, n_ascii_header=6, &
                   field2d_r=field2d_aux)

zl = field2d_aux

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(ZL0_FILE)

call read_2d_input(filename_with_path, &
                   ch_var_name='zl0', n_var_type=1, n_ascii_header=6, &
                   field2d_r=field2d_aux)

zl0 = field2d_aux

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(MASK_PRESENT_FILE)

call read_2d_input(filename_with_path, &
                   ch_var_name='mask', n_var_type=3, n_ascii_header=6, &
                   field2d_r=field2d_aux)

mask = nint(field2d_aux)

#if (defined(ZB_PRESENT_FILE))

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(ZB_PRESENT_FILE)

call read_2d_input(filename_with_path, &
                   ch_var_name='zb', n_var_type=1, n_ascii_header=6, &
                   field2d_r=field2d_aux)

zb = field2d_aux

#else

write(6, fmt='(a)') ' >>> topography1: ZB_PRESENT_FILE not defined,'
write(6, fmt='(a)') '                  thus zb = zl assumed.'

zb = zl

#endif

!-------- Further stuff --------

#if (GRID==0 || GRID==1)

dxi  = DX *1000.0_dp   ! km -> m
deta = DX *1000.0_dp   ! km -> m

xi0  = X0 *1000.0_dp   ! km -> m
eta0 = Y0 *1000.0_dp   ! km -> m

#elif (GRID==2)

dxi  = DLAMBDA *deg2rad
deta = DPHI    *deg2rad

xi0  = LAMBDA_0 *deg2rad
eta0 = PHI_0    *deg2rad

#endif

freeboard_ratio = (RHO_SW-RHO)/RHO_SW

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i) <= 1) then

      zb(j,i) = zl(j,i)   ! ensure consistency

   else if (mask(j,i) == 2) then

#if (MARGIN==1 || MARGIN==2)
      zs(j,i) = zl(j,i)   ! ensure
      zb(j,i) = zl(j,i)   ! consistency
#elif (MARGIN==3)
      zs(j,i) = 0.0_dp    ! present-day
      zb(j,i) = 0.0_dp    ! sea level
#endif

   else if (mask(j,i) == 3) then

#if (MARGIN==1 || (MARGIN==2 && MARINE_ICE_FORMATION==1))
      mask(j,i) = 2   ! floating ice cut off
      zs(j,i) = zl(j,i)
      zb(j,i) = zl(j,i)
#elif (MARGIN==2 && MARINE_ICE_FORMATION==2)
      mask(j,i) = 0   ! floating ice becomes "underwater ice"
      H_ice   = zs(j,i)-zb(j,i)   ! ice thickness
      zs(j,i) = zl(j,i)+H_ice
      zb(j,i) = zl(j,i)
#elif (MARGIN==3)
      H_ice = zs(j,i)-zb(j,i)   ! ice thickness
      zs(j,i) = freeboard_ratio*H_ice   ! ensure properly
      zb(j,i) = zs(j,i)-H_ice           ! floating ice
#endif

   end if

   xi(i)  = xi0  + real(i,dp)*dxi
   eta(j) = eta0 + real(j,dp)*deta

   zm(j,i) = zb(j,i)
   n_cts(j,i) = -1
   kc_cts(j,i) = 0

   H(j,i)   = zs(j,i)-zm(j,i)
   H_c(j,i) = H(j,i)
   H_t(j,i) = 0.0_dp

   dzs_dtau(j,i)  = 0.0_dp
   dzm_dtau(j,i)  = 0.0_dp
   dzb_dtau(j,i)  = 0.0_dp
   dzl_dtau(j,i)  = 0.0_dp
   dH_dtau(j,i)   = 0.0_dp
   dH_c_dtau(j,i) = 0.0_dp
   dH_t_dtau(j,i) = 0.0_dp

end do
end do

mask_old = mask

!-------- Geographic coordinates, metric tensor,
!                                 gradients of the topography --------

do i=0, IMAX
do j=0, JMAX

#if (GRID==0 || GRID==1)   /* Stereographic projection */

   if (F_INV > 1.0e+10_dp) then   ! interpreted as infinity, thus no flattening
                                  ! (spherical planet)

      call stereo_inv_sphere(xi(i), eta(j), R, &
                             LAMBDA0, PHI0, lambda(j,i), phi(j,i))

   else   ! finite inverse flattening (ellipsoidal planet)

      call stereo_inv_ellipsoid(xi(i), eta(j), A, B, &
                                LAMBDA0, PHI0, lambda(j,i), phi(j,i))

   end if

#elif (GRID==2)   /* Geographic coordinates */

   lambda(j,i) = xi(i)
   phi(j,i)    = eta(j)

#endif

end do
end do

call metric()

#if (TOPOGRAD==0)
call topograd_1(dxi, deta, 1)
#elif (TOPOGRAD==1)
call topograd_2(dxi, deta, 1)
#endif

!-------- Corresponding area of grid cells --------

do i=0, IMAX
do j=0, JMAX
   cell_area(j,i) = sq_g11_g(j,i)*sq_g22_g(j,i)*dxi*deta
end do
end do

!-------- Region mask --------

mask_region = -1

#if (defined(MASK_REGION_FILE))

if ( (trim(adjustl(MASK_REGION_FILE)) /= 'none') &
     .and. &
     (trim(adjustl(MASK_REGION_FILE)) /= 'None') &
     .and. &
     (trim(adjustl(MASK_REGION_FILE)) /= 'NONE') ) then
                                      ! read mask_region from file

   filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                        trim(MASK_REGION_FILE)

   call read_2d_input(filename_with_path, &
                      ch_var_name='mask_region', &
                      n_var_type=2, n_ascii_header=6, &
                      field2d_r=field2d_aux)

   mask_region = nint(field2d_aux)

end if

#endif

if (mask_region(0,0) == -1) then

#if (defined(ANT)) /* Antarctic ice sheet */

   do i=0, IMAX
   do j=0, JMAX

      ! Set default values for mask_region

      if ( (phi(j,i) > -77.5_dp*deg2rad) &
           .and. &
           (lambda(j,i) > 277.5_dp*deg2rad) &
           .and. &
           (lambda(j,i) < 307.5_dp*deg2rad) ) then

         mask_region(j,i) = 3   ! AP (Antarctic Peninsula)

      else if ( (lambda(j,i) > 195.0_dp*deg2rad) &
                .and. &
                (lambda(j,i) < 315.0_dp*deg2rad) ) then

         mask_region(j,i) = 2   ! WAIS (West-Antarctic ice sheet)

      else

         mask_region(j,i) = 1   ! EAIS (East-Antarctic ice sheet)

      end if

   end do
   end do

#else /* not Antarctic ice sheet */

   mask_region = 0   ! regions undefined

#endif

end if

end subroutine topography1

!-------------------------------------------------------------------------------
!> Definition of the initial surface and bedrock topography
!! (including gradients) and of the horizontal grid spacings dxi, deta.
!! For ice-free initial topography with relaxed lithosphere.
!-------------------------------------------------------------------------------
subroutine topography2(dxi, deta)

  use read_m, only : read_2d_input

#if (GRID==0 || GRID==1)
  use stereo_proj_m
#endif

  use metric_m
  use topograd_m

implicit none

real(dp), intent(out) :: dxi, deta

integer(i4b) :: i, j, n
real(dp)     :: xi0, eta0

character(len=256) :: filename_with_path

real(dp), dimension(0:JMAX,0:IMAX) :: field2d_aux

!-------- Read topography --------

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(ZL0_FILE)

call read_2d_input(filename_with_path, &
                   ch_var_name='zl0', n_var_type=1, n_ascii_header=6, &
                   field2d_r=field2d_aux)

zl0 = field2d_aux

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(MASK_PRESENT_FILE)

call read_2d_input(filename_with_path, &
                   ch_var_name='mask', n_var_type=3, n_ascii_header=6, &
                   field2d_r=field2d_aux)

mask = nint(field2d_aux)

!-------- Further stuff --------

#if (GRID==0 || GRID==1)

dxi  = DX *1000.0_dp   ! km -> m
deta = DX *1000.0_dp   ! km -> m

xi0  = X0 *1000.0_dp   ! km -> m
eta0 = Y0 *1000.0_dp   ! km -> m

#elif (GRID==2)

dxi  = DLAMBDA *deg2rad
deta = DPHI    *deg2rad

xi0  = LAMBDA_0 *deg2rad
eta0 = PHI_0    *deg2rad

#endif

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i) <= 1) then
      mask(j,i) = 1
      zs(j,i) = zl0(j,i)
      zb(j,i) = zl0(j,i)
      zl(j,i) = zl0(j,i)
   else   ! (mask(j,i) >= 2)
      mask(j,i) = 2
#if (MARGIN==1 || MARGIN==2)
      zs(j,i) = zl0(j,i)
      zb(j,i) = zl0(j,i)
#elif (MARGIN==3)
      zs(j,i) = 0.0_dp   ! present-day
      zb(j,i) = 0.0_dp   ! sea level
#endif
      zl(j,i) = zl0(j,i)
   end if

   xi(i)  = xi0  + real(i,dp)*dxi
   eta(j) = eta0 + real(j,dp)*deta

   zm(j,i) = zb(j,i)
   n_cts(j,i) = -1
   kc_cts(j,i) = 0

   H(j,i)   = 0.0_dp
   H_c(j,i) = 0.0_dp
   H_t(j,i) = 0.0_dp

   dzs_dtau(j,i)  = 0.0_dp
   dzm_dtau(j,i)  = 0.0_dp
   dzb_dtau(j,i)  = 0.0_dp
   dzl_dtau(j,i)  = 0.0_dp
   dH_dtau(j,i)   = 0.0_dp
   dH_c_dtau(j,i) = 0.0_dp
   dH_t_dtau(j,i) = 0.0_dp

end do
end do

mask_old = mask

!-------- Geographic coordinates, metric tensor,
!                                 gradients of the topography --------

do i=0, IMAX
do j=0, JMAX

#if (GRID==0 || GRID==1)   /* Stereographic projection */

   if (F_INV > 1.0e+10_dp) then   ! interpreted as infinity, thus no flattening
                                  ! (spherical planet)

      call stereo_inv_sphere(xi(i), eta(j), R, &
                             LAMBDA0, PHI0, lambda(j,i), phi(j,i))

   else   ! finite inverse flattening (ellipsoidal planet)

      call stereo_inv_ellipsoid(xi(i), eta(j), A, B, &
                                LAMBDA0, PHI0, lambda(j,i), phi(j,i))

   end if

#elif (GRID==2)   /* Geographic coordinates */

   lambda(j,i) = xi(i)
   phi(j,i)    = eta(j)

#endif

end do
end do

call metric()

#if (TOPOGRAD==0)
call topograd_1(dxi, deta, 1)
#elif (TOPOGRAD==1)
call topograd_2(dxi, deta, 1)
#endif

!-------- Corresponding area of grid cells --------

do i=0, IMAX
do j=0, JMAX
   cell_area(j,i) = sq_g11_g(j,i)*sq_g22_g(j,i)*dxi*deta
end do
end do

!-------- Region mask --------

mask_region = -1

#if (defined(MASK_REGION_FILE))

if ( (trim(adjustl(MASK_REGION_FILE)) /= 'none') &
     .and. &
     (trim(adjustl(MASK_REGION_FILE)) /= 'None') &
     .and. &
     (trim(adjustl(MASK_REGION_FILE)) /= 'NONE') ) then
                                      ! read mask_region from file

   filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                        trim(MASK_REGION_FILE)

   call read_2d_input(filename_with_path, &
                      ch_var_name='mask_region', &
                      n_var_type=2, n_ascii_header=6, &
                      field2d_r=field2d_aux)

   mask_region = nint(field2d_aux)

end if

#endif

if (mask_region(0,0) == -1) then

#if (defined(ANT)) /* Antarctic ice sheet */

   do i=0, IMAX
   do j=0, JMAX

      ! Set default values for mask_region

      if ( (phi(j,i) > -77.5_dp*deg2rad) &
           .and. &
           (lambda(j,i) > 277.5_dp*deg2rad) &
           .and. &
           (lambda(j,i) < 307.5_dp*deg2rad) ) then

         mask_region(j,i) = 3   ! AP (Antarctic Peninsula)

      else if ( (lambda(j,i) > 195.0_dp*deg2rad) &
                .and. &
                (lambda(j,i) < 315.0_dp*deg2rad) ) then

         mask_region(j,i) = 2   ! WAIS (West-Antarctic ice sheet)

      else

         mask_region(j,i) = 1   ! EAIS (East-Antarctic ice sheet)

      end if

   end do
   end do

#else /* not Antarctic ice sheet */

   mask_region = 0   ! regions undefined

#endif

end if

end subroutine topography2

!-------------------------------------------------------------------------------
!> Definition of the initial surface and bedrock topography
!! (including gradients) and of the horizontal grid spacings dxi, deta.
!! For initial topography from previous simulation.
!-------------------------------------------------------------------------------
subroutine topography3(dxi, deta, anfdatname)

  use read_m, only : read_tms_nc, read_2d_input

#if (GRID==0 || GRID==1)
  use stereo_proj_m
#endif

  use metric_m
  use topograd_m

implicit none

character(len=256), intent(in) :: anfdatname

real(dp),          intent(out) :: dxi, deta

integer(i4b) :: i, j, n

character(len=256) :: filename_with_path

real(dp), dimension(0:JMAX,0:IMAX) :: field2d_aux

!-------- Read data from time-slice file of previous simulation --------

call read_tms_nc(anfdatname)

!-------- Read topography of the relaxed bedrock --------

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(ZL0_FILE)

call read_2d_input(filename_with_path, &
                   ch_var_name='zl0', n_var_type=1, n_ascii_header=6, &
                   field2d_r=field2d_aux)

zl0 = field2d_aux

!-------- Further stuff --------

#if (GRID==0 || GRID==1)

dxi  = DX *1000.0_dp   ! km -> m
deta = DX *1000.0_dp   ! km -> m

#elif (GRID==2)

dxi  = DLAMBDA *deg2rad
deta = DPHI    *deg2rad

#endif

!-------- Geographic coordinates, metric tensor,
!                                 gradients of the topography --------

do i=0, IMAX
do j=0, JMAX

#if (GRID==0 || GRID==1)   /* Stereographic projection */

   if (F_INV > 1.0e+10_dp) then   ! interpreted as infinity, thus no flattening
                                  ! (spherical planet)

      call stereo_inv_sphere(xi(i), eta(j), R, &
                             LAMBDA0, PHI0, lambda(j,i), phi(j,i))

   else   ! finite inverse flattening (ellipsoidal planet)

      call stereo_inv_ellipsoid(xi(i), eta(j), A, B, &
                                LAMBDA0, PHI0, lambda(j,i), phi(j,i))

   end if

#elif (GRID==2)   /* Geographic coordinates */

   lambda(j,i) = xi(i)
   phi(j,i)    = eta(j)

#endif

end do
end do

call metric()

#if (TOPOGRAD==0)
call topograd_1(dxi, deta, 1)
#elif (TOPOGRAD==1)
call topograd_2(dxi, deta, 1)
#endif

!-------- Corresponding area of grid cells --------

do i=0, IMAX
do j=0, JMAX
   cell_area(j,i) = sq_g11_g(j,i)*sq_g22_g(j,i)*dxi*deta
end do
end do

!-------- Region mask --------

mask_region = -1

#if (defined(MASK_REGION_FILE))

if ( (trim(adjustl(MASK_REGION_FILE)) /= 'none') &
     .and. &
     (trim(adjustl(MASK_REGION_FILE)) /= 'None') &
     .and. &
     (trim(adjustl(MASK_REGION_FILE)) /= 'NONE') ) then
                                      ! read mask_region from file

   filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                        trim(MASK_REGION_FILE)

   call read_2d_input(filename_with_path, &
                      ch_var_name='mask_region', &
                      n_var_type=2, n_ascii_header=6, &
                      field2d_r=field2d_aux)

   mask_region = nint(field2d_aux)

end if

#endif

if (mask_region(0,0) == -1) then

#if (defined(ANT)) /* Antarctic ice sheet */

   do i=0, IMAX
   do j=0, JMAX

      ! Set default values for mask_region

      if ( (phi(j,i) > -77.5_dp*deg2rad) &
           .and. &
           (lambda(j,i) > 277.5_dp*deg2rad) &
           .and. &
           (lambda(j,i) < 307.5_dp*deg2rad) ) then

         mask_region(j,i) = 3   ! AP (Antarctic Peninsula)

      else if ( (lambda(j,i) > 195.0_dp*deg2rad) &
                .and. &
                (lambda(j,i) < 315.0_dp*deg2rad) ) then

         mask_region(j,i) = 2   ! WAIS (West-Antarctic ice sheet)

      else

         mask_region(j,i) = 1   ! EAIS (East-Antarctic ice sheet)

      end if

   end do
   end do

#else /* not Antarctic ice sheet */

   mask_region = 0   ! regions undefined

#endif

end if

end subroutine topography3

!-------------------------------------------------------------------------------
!> Set the value of the auxiliary variable flag_grads_nc_tweaks.
!-------------------------------------------------------------------------------
  subroutine set_flag_grads_nc_tweaks()

  implicit none

  character(len=16) :: ch_value

  flag_grads_nc_tweaks = .false.   ! default

!-------- Try environment variable --------

  call get_environment_variable('SICO_GRADS_NC_TWEAKS', ch_value)

  if ( (trim(ch_value)=='true') &
       .or.(trim(ch_value)=='True').or.(trim(ch_value)=='TRUE') ) &
     flag_grads_nc_tweaks = .true.

  if ( (trim(ch_value)=='yes') &   ! obsolete, but still supported
       .or.(trim(ch_value)=='Yes').or.(trim(ch_value)=='YES') &
       .or.(trim(ch_value)=='y').or.(trim(ch_value)=='Y') ) &
     flag_grads_nc_tweaks = .true.

!-------- Try preprocessor switch --------

#if (defined(GRADS_NC_TWEAKS))
#if (GRADS_NC_TWEAKS==1)
  flag_grads_nc_tweaks = .true.
#else
  flag_grads_nc_tweaks = .false.
#endif
#endif

  end subroutine set_flag_grads_nc_tweaks

!-------------------------------------------------------------------------------

end module sico_init_m
!
