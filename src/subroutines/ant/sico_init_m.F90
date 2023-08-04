!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ i n i t _ m
!
!> @file
!!
!! Initialisations for SICOPOLIS.
!!
!! @section Copyright
!!
!! Copyright 2009-2023 Ralf Greve
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
!> Initialisations for SICOPOLIS.
!<------------------------------------------------------------------------------
module sico_init_m

  use sico_types_m
  use sico_variables_m
  use sico_vars_m
  use error_m
#if defined(ALLOW_TAPENADE)
  use globals
#endif
  implicit none

  public

contains

!-------------------------------------------------------------------------------
!> Main routine of sico_init_m: Initialisations for SICOPOLIS.
!<------------------------------------------------------------------------------
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

#if (GRID==0 || GRID==1)
  use stereo_proj_m
#endif

  use read_m, only : read_target_topo_nc, &
                     read_2d_input, read_kei, read_phys_para

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

integer(i4b),       intent(out) :: ndat2d, ndat3d
integer(i4b),       intent(out) :: n_output
real(dp),           intent(out) :: delta_ts, glac_index
real(dp),           intent(out) :: mean_accum
real(dp),           intent(out) :: dtime, dtime_temp, dtime_wss, &
                                   dtime_out, dtime_ser
real(dp),           intent(out) :: time, time_init, time_end, time_output(100)
real(dp),           intent(out) :: dxi, deta, dzeta_c, dzeta_t, dzeta_r
real(dp),           intent(out) :: z_mar

integer(i4b)       :: i, j, kc, kt, kr, m, n, ir, jr, n1, n2
integer(i4b)       :: ios, ios1, ios2, ios3, ios4
integer(i4b)       :: ierr
integer(i4b)       :: n_q_geo_mod
integer(i4b), dimension(0:JMAX,0:IMAX) :: mask_ref
real(dp)           :: dtime0, dtime_temp0, dtime_wss0, dtime_out0, dtime_ser0
real(dp)           :: time_init0, time_end0
#if (OUTPUT==2 || OUTPUT==3)
real(dp)           :: time_output0(N_OUTPUT)
#endif
real(dp)           :: d_dummy
real(dp)           :: larmip_qbm_anom_aux(5)
character(len=256) :: anfdatname, target_topo_dat_name
character(len=256) :: filename_with_path
character(len=256) :: shell_command
character          :: ch_dummy
logical            :: flag_init_output, flag_3d_output

real(dp), dimension(0:JMAX,0:IMAX) :: field2d_aux

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

#if (defined(SMB_CORR_FILE))
real(dp), dimension(0:IMAX,0:JMAX) :: smb_corr_in_conv
#endif

#if (defined(INITMIP_SMB_ANOM_FILE))
real(dp), dimension(0:IMAX,0:JMAX) :: smb_anom_initmip_conv
#endif

#if (defined(INITMIP_BMB_ANOM_FILE))
real(dp), dimension(0:IMAX,0:JMAX) :: ab_anom_initmip_conv
#endif

#if (defined(LARMIP_REGIONS_FILE))
real(dp), dimension(0:IMAX,0:JMAX) :: regions_larmip_conv
#endif

#if (FLOATING_ICE_BASAL_MELTING==6)
real(dp), dimension(0:IMAX,0:JMAX,0:NZ_TF_BM) :: tf_bm_present_aux
#endif

#if (ICE_SHELF_COLLAPSE_MASK==1)
real(dp), dimension(0:IMAX,0:JMAX) :: H_ref_retreat_conv
#endif

integer(i4b) :: dimid, ncid, ncv
!   dimid:       Dimension ID
!    ncid:       File ID
!     ncv:       Variable ID

character(len=64), parameter :: thisroutine = 'sico_init'

character(len=64), parameter :: fmt1 = '(a)', &
                                fmt2 = '(a,i0)', &
                                fmt3 = '(a,es12.4)'

write(unit=6, fmt='(a)') ' '
write(unit=6, fmt='(a)') ' -------- sico_init --------'
write(unit=6, fmt='(a)') ' '

!-------- Name of the computational domain --------

#if (defined(ANT))
ch_domain_long  = 'Antarctica'
ch_domain_short = 'ant'

#elif (defined(ASF))
ch_domain_long  = 'Austfonna'
ch_domain_short = 'asf'

#elif (defined(EISMINT))
ch_domain_long  = 'EISMINT'
ch_domain_short = 'eismint'

#elif (defined(GRL))
ch_domain_long  = 'Greenland'
ch_domain_short = 'grl'

#elif (defined(NHEM))
ch_domain_long  = 'Northern hemisphere'
ch_domain_short = 'nhem'

#elif (defined(SCAND))
ch_domain_long  = 'Scandinavia and Eurasia'
ch_domain_short = 'scand'

#elif (defined(TIBET))
ch_domain_long  = 'Tibet'
ch_domain_short = 'tibet'

#elif (defined(NMARS))
ch_domain_long  = 'North polar cap of Mars'
ch_domain_short = 'nmars'

#elif (defined(SMARS))
ch_domain_long  = 'South polar cap of Mars'
ch_domain_short = 'smars'

#elif (defined(XYZ))
ch_domain_long  = 'XYZ'
ch_domain_short = 'xyz'
#if (defined(HEINO))
ch_domain_long  = trim(ch_domain_long)//'/ISMIP HEINO'
#endif

#else
errormsg = ' >>> sico_init: No valid domain specified!'
call error(errormsg)
#endif

!-------- Some initial values --------

n_output = 0

dtime       = 0.0_dp
dtime_temp  = 0.0_dp
dtime_wss   = 0.0_dp
dtime_out   = 0.0_dp
dtime_ser   = 0.0_dp

time        = 0.0_dp
time_init   = 0.0_dp
time_end    = 0.0_dp
time_output = 0.0_dp

!-------- Initialisation of the Library of Iterative Solvers Lis,
!                                                     if required --------

#if (CALCTHK==3 || CALCTHK==6 || MARGIN==3 || DYNAMICS==2)
#if !defined(ALLOW_TAPENADE) /* Normal */
  call lis_initialize(ierr)
#else /* Tapenade */
  call lis_init_f(ierr)
#endif /* Normal vs. Tapenade */
#endif

!-------- Read physical parameters --------

#if (defined(YEAR_SEC))
year2sec = YEAR_SEC
#else
year2sec = 3.1556925445e+07_dp
              ! IUPAC-IUGS year for epoch 2000.0
              ! (Holden et al., 2011, PAC, doi:10.1351/PAC-REC-09-01-22)
#endif

sec2year = 1.0_dp/year2sec

call read_phys_para()

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

#elif (GRID==2)
errormsg = ' >>> sico_init: GRID==2 not allowed for this application!'
call error(errormsg)

#endif

#else /* CHECK_RES_IMAX_JMAX==0 */

write(6, fmt='(a)') ' >>> sico_init: CHECK_RES_IMAX_JMAX==0'
write(6, fmt='(a)') '      -> compatibility check between horizontal resolution'
write(6, fmt='(a)') '         and number of grid points not performed.'
write(6, fmt='(a)') ' '

#endif /* CHECK_RES_IMAX_JMAX */

#else /* Tapenade */

print *, ' >>> sico_init: not using compare_float for adjoint applications!'
print *, '                double-check resolutions are integer multiples of'
print *, '                domain size.'

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
                   ! and precipitation data

#endif

!-------- Initialization of numerical time steps --------

dtime0      = DTIME0
dtime_temp0 = DTIME_TEMP0
#if (REBOUND==2)
dtime_wss0  = DTIME_WSS0
#endif

!-------- Further initializations --------

dzeta_c = 1.0_dp/real(KCMAX,dp)
dzeta_t = 1.0_dp/real(KTMAX,dp)
dzeta_r = 1.0_dp/real(KRMAX,dp)

ndat2d = 1
ndat3d = 1

flag_calc_temp = .true.

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

write(10, fmt=trim(fmt1)) 'Physical-parameter file = '//PHYS_PARA_FILE
write(10, fmt=trim(fmt1)) ' '

write(10, fmt=trim(fmt2)) 'GRID = ', GRID
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
errormsg = ' >>> sico_init: GRID==2 not allowed for this application!'
call error(errormsg)
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
#elif (MARINE_ICE_CALVING==4 || MARINE_ICE_CALVING==5 \
       || MARINE_ICE_CALVING==6 || MARINE_ICE_CALVING==7)
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
if ( trim(adjustl(MASK_REGION_FILE)) /= 'none' ) then
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
write(10, fmt=trim(fmt3)) 'time_target_topo_init  =', TIME_TARGET_TOPO_INIT0
write(10, fmt=trim(fmt3)) 'time_target_topo_final =', TIME_TARGET_TOPO_FINAL0
write(10, fmt=trim(fmt3)) 'target_topo_tau_0 =', TARGET_TOPO_TAU0
write(10, fmt=trim(fmt1)) 'Target-topography file = '//TARGET_TOPO_DAT_NAME
write(10, fmt=trim(fmt1)) 'Path to target-topography file = '//TARGET_TOPO_PATH
#endif

#if (THK_EVOL==3)
write(10, fmt=trim(fmt3)) 'target_topo_tau_0 =', TARGET_TOPO_TAU0
write(10, fmt=trim(fmt1)) 'Target-topography file = '//TARGET_TOPO_DAT_NAME
write(10, fmt=trim(fmt1)) 'Path to target-topography file = '//TARGET_TOPO_PATH
#endif

#if (THK_EVOL==4)
write(10, fmt=trim(fmt1)) 'Maximum ice extent mask file = '//MASK_MAXEXTENT_FILE
#endif

#if (CALCTHK==2 || CALCTHK==3 || CALCTHK==5 || CALCTHK==6)
write(10, fmt=trim(fmt3))  'ovi_weight =', OVI_WEIGHT
#if (CALCTHK==2 || CALCTHK==5)
write(10, fmt=trim(fmt3))  'omega_sor =', OMEGA_SOR
#if (ITER_MAX_SOR>0)
write(10, fmt=trim(fmt2)) 'iter_max_sor = ', ITER_MAX_SOR
#endif
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
write(10, fmt=trim(fmt2)) 'TEMP_PRESENT_PARA = ', TEMP_PRESENT_PARA
#if (defined(TEMP_PRESENT_OFFSET))
write(10, fmt=trim(fmt3))  'TEMP_PRESENT_OFFSET =', TEMP_PRESENT_OFFSET
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
write(10, fmt=trim(fmt1)) 'temp_ma_anom file = '//TEMP_MA_ANOM_FILE
write(10, fmt=trim(fmt3)) 'temp_ma_anom fact = ', TEMP_MA_ANOM_FACT
write(10, fmt=trim(fmt1)) 'temp_mj_anom file = '//TEMP_MJ_ANOM_FILE
write(10, fmt=trim(fmt3)) 'temp_mj_anom fact = ', TEMP_MJ_ANOM_FACT
#endif
write(10, fmt=trim(fmt1)) ' '

write(10, fmt=trim(fmt2)) 'ACCSURFACE = ', ACCSURFACE

#if (ACCSURFACE<=5)
write(10, fmt=trim(fmt1)) 'precip_present_file = '//PRECIP_PRESENT_FILE
#endif

#if (ACCSURFACE==1)
write(10, fmt=trim(fmt3)) 'accfact =', ACCFACT
#elif (ACCSURFACE==2 || ACCSURFACE==3)
write(10, fmt=trim(fmt3)) 'gamma_s =', GAMMA_S
#endif
#if (ACCSURFACE<=3)
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

#if (ABLSURFACE==3)
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
if ( trim(adjustl(SMB_CORR_FILE)) /= 'none' ) then
   write(10, fmt=trim(fmt1)) 'smb_corr_file = '//SMB_CORR_FILE
   write(10, fmt=trim(fmt1)) ' '
end if
#endif

#if (defined(INITMIP_SMB_ANOM_FILE))
if ( trim(adjustl(INITMIP_SMB_ANOM_FILE)) /= 'none' ) then
   write(10, fmt=trim(fmt1)) 'initmip_smb_anom file = '//INITMIP_SMB_ANOM_FILE
end if
write(10, fmt=trim(fmt1)) ' '
#endif

#if (defined(ICE_SHELF_COLLAPSE_MASK))
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
#if (SLIDE_LAW==2 || SLIDE_LAW==3)
write(10, fmt=trim(fmt3)) 'red_pres_limit_fact =', RED_PRES_LIMIT_FACT
#endif
#if (BASAL_HYDROLOGY==1 \
       && defined(HYDRO_SLIDE_SAT_FCT) \
       && defined(C_HW_SLIDE) && defined(HW0_SLIDE))
write(10, fmt=trim(fmt2)) 'HYDRO_SLIDE_SAT_FCT = ', HYDRO_SLIDE_SAT_FCT
write(10, fmt=trim(fmt3)) 'c_Hw_slide =', C_HW_SLIDE
write(10, fmt=trim(fmt3)) 'Hw0_slide  =', HW0_SLIDE
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
#endif

#if (defined(INITMIP_BMB_ANOM_FILE))
if ( trim(adjustl(INITMIP_BMB_ANOM_FILE)) /= 'none' ) then
   write(10, fmt=trim(fmt1)) 'initmip_bmb_anom file = '//INITMIP_BMB_ANOM_FILE
end if
#endif

#if (defined(LARMIP_REGIONS_FILE))
write(10, fmt=trim(fmt1)) 'larmip_regions_file = '//LARMIP_REGIONS_FILE
if ( trim(adjustl(LARMIP_REGIONS_FILE)) /= 'none' ) then
   larmip_qbm_anom_aux = LARMIP_QBM_ANOM
   write(10, fmt=trim(fmt3)) 'larmip_qbm_anom_1 =', larmip_qbm_anom_aux(1)
   write(10, fmt=trim(fmt3)) 'larmip_qbm_anom_2 =', larmip_qbm_anom_aux(2)
   write(10, fmt=trim(fmt3)) 'larmip_qbm_anom_3 =', larmip_qbm_anom_aux(3)
   write(10, fmt=trim(fmt3)) 'larmip_qbm_anom_4 =', larmip_qbm_anom_aux(4)
   write(10, fmt=trim(fmt3)) 'larmip_qbm_anom_5 =', larmip_qbm_anom_aux(5)
end if
#endif
write(10, fmt=trim(fmt1)) ' '

#endif

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
    !!! Climatology extraction hack (must not be used routinely) !!!
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

write(10, fmt=trim(fmt1)) 'Program version and date: '//VERSION//' / '//DATE

close(10, status='keep')

#if (defined(CLIMATOLOGY_EXTRACTION_HACK))
    !!! Climatology extraction hack (must not be used routinely) !!!
#if (OUTPUT_FLUX_VARS==2)
write(6,'(/1x,a)') &
   '>>> sico_init: CLIMATOLOGY_EXTRACTION_HACK defined!'
write(6, '(1x,a)') &
   '               Only for extracting a climatology,'
write(6, '(1x,a/)') &
   '               must not be used routinely!!!'
#else
errormsg = ' >>> sico_init: CLIMATOLOGY_EXTRACTION_HACK requires' &
         //                 end_of_line &
         //'                OUTPUT_FLUX_VARS==2!'
call error(errormsg)
#endif
#endif

!-------- Conversion of time quantities --------

year_zero  = year_zero*year2sec     ! a --> s
time_init  = time_init0*year2sec    ! a --> s
time_end   = time_end0*year2sec     ! a --> s
dtime      = dtime0*year2sec        ! a --> s
dtime_temp = dtime_temp0*year2sec   ! a --> s
#if (REBOUND==2)
dtime_wss  = dtime_wss0*year2sec    ! a --> s
#endif
dtime_ser  = dtime_ser0*year2sec    ! a --> s
#if (OUTPUT==1 || OUTPUT==3)
dtime_out  = dtime_out0*year2sec    ! a --> s
#endif
#if (OUTPUT==2 || OUTPUT==3)
do n=1, n_output
   time_output(n) = time_output0(n)*year2sec  ! a --> s
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

print *, ' >>> sico_init: compare_float not used in adjoint applications!'
print *, '                double check your time step sizes are multples' 
print *, '                of each other.' 

#endif /* Normal vs. Tapenade */

#if (THK_EVOL==2)
time_target_topo_init  = TIME_TARGET_TOPO_INIT0 *year2sec   ! a --> s
time_target_topo_final = TIME_TARGET_TOPO_FINAL0*year2sec   ! a --> s
target_topo_tau_0 = TARGET_TOPO_TAU0 *year2sec   ! a --> s
#endif

#if (THK_EVOL==3)
target_topo_tau_0 = TARGET_TOPO_TAU0 *year2sec   ! a --> s
#endif

#if (ACCSURFACE==7 && ABLSURFACE==7)
target_topo_tau_0 = TARGET_TOPO_TAU0 *year2sec   ! a --> s
#endif

time = time_init

!-------- Reading of present-day mean-annual precipitation rate --------

#if (ACCSURFACE<=5)

#if (GRID==0 || GRID==1)

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(PRECIP_PRESENT_FILE)

call read_2d_input(filename_with_path, &
                   ch_var_name='precip_ma_present', &
                   n_var_type=1, n_ascii_header=6, &
                   field2d_r=field2d_aux)

precip_ma_present = field2d_aux *(1.0e-03_dp/year2sec)*(RHO_W/RHO)
                                 ! mm/a water equiv. -> m/s ice equiv.

!  ------ Present monthly precipitation rates
!         (assumed to be equal to the mean annual precipitation rate)

do n=1, 12   ! month counter
do i=0, IMAX
do j=0, JMAX
   precip_present(j,i,n) = precip_ma_present(j,i)
end do
end do
end do

#elif (GRID==2)

errormsg = ' >>> sico_init: GRID==2 not allowed for this application!'
call error(errormsg)

#endif

#endif

!-------- Read reference topography (present-day topography) --------

#if (GRID==0 || GRID==1)

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(ZS_PRESENT_FILE)

call read_2d_input(filename_with_path, &
                   ch_var_name='zs', n_var_type=1, n_ascii_header=6, &
                   field2d_r=field2d_aux)

zs_ref = field2d_aux

#elif (GRID==2)

errormsg = ' >>> sico_init: GRID==2 not allowed for this application!'
call error(errormsg)

#endif

#if (GRID==0 || GRID==1)

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(MASK_PRESENT_FILE)

call read_2d_input(filename_with_path, &
                   ch_var_name='mask', n_var_type=3, n_ascii_header=6, &
                   field2d_r=field2d_aux)

mask_ref = nint(field2d_aux)

do i=0, IMAX
do j=0, JMAX
   if (mask_ref(j,i) >= 2) zs_ref(j,i) = 0.0_dp
                 ! resetting elevations over the ocean
                 ! to the present-day sea surface
end do
end do

#elif (GRID==2)

errormsg = ' >>> sico_init: GRID==2 not allowed for this application!'
call error(errormsg)

#endif

!-------- Reading of LGM mean-annual precipitation-rate anomaly --------

#if (ACCSURFACE==5)

#if (GRID==0 || GRID==1)

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(PRECIP_ANOM_FILE)

call read_2d_input(filename_with_path, &
                   ch_var_name='precip_ma_lgm_anom', &
                   n_var_type=1, n_ascii_header=6, &
                   field2d_r=field2d_aux)

precip_ma_lgm_anom = field2d_aux

precip_ma_lgm_anom = precip_ma_lgm_anom * PRECIP_ANOM_FACT

#endif

!  ------ LGM monthly precipitation-rate anomalies (assumed to be
!             equal to the mean annual precipitation-rate anomaly)

do n=1, 12   ! month counter
do i=0, IMAX
do j=0, JMAX

   precip_lgm_anom(j,i,n) = max(precip_ma_lgm_anom(j,i), eps)
                            ! positive values ensured

#if (PRECIP_ANOM_INTERPOL==1)
   gamma_precip_lgm_anom(j,i,n) = 0.0_dp   ! dummy values
#elif (PRECIP_ANOM_INTERPOL==2)
   gamma_precip_lgm_anom(j,i,n) = -log(precip_lgm_anom(j,i,n))
#else
   errormsg = ' >>> sico_init: Wrong value of switch PRECIP_ANOM_INTERPOL!'
   call error(errormsg)
#endif

end do
end do
end do

#endif

!-------- Mean accumulation --------

mean_accum = MEAN_ACCUM*(1.0e-03_dp/year2sec)*(RHO_W/RHO)
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

#if ( (THK_EVOL==2 || THK_EVOL==3) \
      || \
      (ACCSURFACE==7 && ABLSURFACE==7) )

target_topo_dat_name = trim(TARGET_TOPO_DAT_NAME)

call read_target_topo_nc(target_topo_dat_name)

#endif

!-------- Reading of the maximum ice extent mask --------

#if (THK_EVOL==4)

#if (GRID==0 || GRID==1)

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(MASK_MAXEXTENT_FILE)

call read_2d_input(filename_with_path, &
                   ch_var_name='mask_maxextent', &
                   n_var_type=3, n_ascii_header=6, &
                   field2d_r=field2d_aux)

mask_maxextent = nint(field2d_aux)

#elif (GRID==2)

errormsg = ' >>> sico_init: GRID==2 not allowed for this application!'
call error(errormsg)

#endif

#endif

!-------- Reading of LGM mean-annual and mean-January (summer)
!                                surface-temperature anomalies --------

#if (TSURFACE==5)

#if (GRID==0 || GRID==1)

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(TEMP_MA_ANOM_FILE)

call read_2d_input(filename_with_path, &
                   ch_var_name='temp_ma_lgm_anom', &
                   n_var_type=1, n_ascii_header=6, &
                   field2d_r=field2d_aux)

temp_ma_lgm_anom = field2d_aux

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(TEMP_MJ_ANOM_FILE)

call read_2d_input(filename_with_path, &
                   ch_var_name='temp_mj_lgm_anom', &
                   n_var_type=1, n_ascii_header=6, &
                   field2d_r=field2d_aux)

temp_mj_lgm_anom = field2d_aux

temp_ma_lgm_anom = temp_ma_lgm_anom * TEMP_MA_ANOM_FACT
temp_mj_lgm_anom = temp_mj_lgm_anom * TEMP_MJ_ANOM_FACT

#endif

#endif

!-------- Read data for delta_ts --------

#if (TSURFACE==4)

filename_with_path = trim(IN_PATH)//'/general/'//trim(GRIP_TEMP_FILE)

open(21, iostat=ios, file=trim(filename_with_path), status='old')

if (ios /= 0) then
   errormsg = ' >>> sico_init: Error when opening the data file for delta_ts!'
   call error(errormsg)
end if
#if !defined(ALLOW_TAPENADE)
read(21, fmt=*) ch_dummy, grip_time_min, grip_time_stp, grip_time_max

if (ch_dummy /= '#') then
   errormsg = ' >>> sico_init: grip_time_min, grip_time_stp, grip_time_max' &
            //         end_of_line &
            //'        not defined in data file for delta_ts!'
   call error(errormsg)
end if

ndata_grip = (grip_time_max-grip_time_min)/grip_time_stp

allocate(griptemp(0:ndata_grip))
#else 
read(21, fmt=*)
#endif
do n=0, ndata_grip
   read(21, fmt=*) d_dummy, griptemp(n)
end do

close(21, status='keep')

#endif

!-------- Read data for the glacial index --------

#if (TSURFACE==5)

filename_with_path = trim(IN_PATH)//'/general/'//trim(GLAC_IND_FILE)

open(21, iostat=ios, file=trim(filename_with_path), status='old')

if (ios /= 0) then
   errormsg = ' >>> sico_init: Error when opening the glacial-index file!'
   call error(errormsg)
end if

read(21, fmt=*) ch_dummy, gi_time_min, gi_time_stp, gi_time_max

if (ch_dummy /= '#') then
   errormsg = ' >>> sico_init: gi_time_min, gi_time_stp, gi_time_max' &
            //         end_of_line &
            //'        not defined in glacial-index file!'
   call error(errormsg)
end if

ndata_gi = (gi_time_max-gi_time_min)/gi_time_stp

allocate(glacial_index(0:ndata_gi))

do n=0, ndata_gi
   read(21, fmt=*) d_dummy, glacial_index(n)
end do

close(21, status='keep')

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
      zs_ref_conv(i,j) = zs_ref(j,i)
                         ! use previously read reference topography
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
      smb_climatol(j,i) = smb_climatol_conv(i,j) /year2sec
                                       ! m/a ice equiv. -> m/s ice equiv.
   else
      errormsg = ' >>> sico_init: Unit of SMB_clim could not be determined!'
      call error(errormsg)
   end if

   zs_ref(j,i) = zs_ref_conv(i,j)

end do
end do

#endif

!-------- Prescribed surface mass balance correction --------

smb_corr_in = 0.0_dp

#if (defined(SMB_CORR_FILE))

if (trim(adjustl(SMB_CORR_FILE)) /= 'none') then

   filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                        trim(SMB_CORR_FILE)

   call read_2d_input(filename_with_path, &
                      ch_var_name='DSMB', n_var_type=1, n_ascii_header=6, &
                      field2d_r=field2d_aux)

   smb_corr_in = field2d_aux /year2sec
                             ! m/a ice equiv. -> m/s ice equiv.

end if

#endif

!-------- Reading of ISMIP6 SMB and BMB anomaly data --------

!  ------ SMB (InitMIP)

#if (defined(INITMIP_SMB_ANOM_FILE))

if ( trim(adjustl(INITMIP_SMB_ANOM_FILE)) /= 'none' ) then

   filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)// &
                                      '/'//trim(INITMIP_SMB_ANOM_FILE)

   ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

   if (ios /= nf90_noerr) then
      errormsg = ' >>> sico_init: Error when opening the file' &
               //                 end_of_line &
               //'                for the ISMIP6 SMB anomaly data!'
      call error(errormsg)
   end if

   call check( nf90_inq_varid(ncid, 'asmb', ncv) )
   call check( nf90_get_var(ncid, ncv, smb_anom_initmip_conv) )

   call check( nf90_close(ncid) )

   do i=0, IMAX
   do j=0, JMAX
      smb_anom_initmip(j,i) = smb_anom_initmip_conv(i,j) /year2sec
                                   ! m/a ice equiv. -> m/s ice equiv.
   end do
   end do

else
   smb_anom_initmip = 0.0_dp
end if

#endif

!  ------ BMB (InitMIP)

#if (defined(INITMIP_BMB_ANOM_FILE))

if ( trim(adjustl(INITMIP_BMB_ANOM_FILE)) /= 'none' ) then

   filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)// &
                                      '/'//trim(INITMIP_BMB_ANOM_FILE)

   ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

   if (ios /= nf90_noerr) then
      errormsg = ' >>> sico_init: Error when opening the file' &
               //                 end_of_line &
               //'                for the ISMIP6 BMB anomaly data!'
      call error(errormsg)
   end if

   call check( nf90_inq_varid(ncid, 'abmb', ncv) )
   call check( nf90_get_var(ncid, ncv, ab_anom_initmip_conv) )

   call check( nf90_close(ncid) )

   do i=0, IMAX
   do j=0, JMAX
      ab_anom_initmip(j,i) = ab_anom_initmip_conv(i,j) /year2sec
                                   ! m/a ice equiv. -> m/s ice equiv.
   end do
   end do

else
   ab_anom_initmip = 0.0_dp
end if

#endif

!  ------ BMB (LARMIP)

#if (defined(LARMIP_REGIONS_FILE))

if ( trim(adjustl(LARMIP_REGIONS_FILE)) /= 'none' ) then

   filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)// &
                                      '/'//trim(LARMIP_REGIONS_FILE)

   ios = nf90_open(trim(filename_with_path), NF90_NOWRITE, ncid)

   if (ios /= nf90_noerr) then
      errormsg = ' >>> sico_init: Error when opening the file' &
               //                 end_of_line &
               //'                for the LARMIP BMB regions!'
      call error(errormsg)
   end if

   call check( nf90_inq_varid(ncid, 'regions', ncv) )
   call check( nf90_get_var(ncid, ncv, regions_larmip_conv) )

   call check( nf90_close(ncid) )

   do i=0, IMAX
   do j=0, JMAX
      n_larmip_region(j,i) = nint(regions_larmip_conv(i,j))
   end do
   end do

   ab_anom_larmip      = 0.0_dp
   ab_anom_larmip(1:5) = LARMIP_QBM_ANOM / year2sec
                                   ! m/a ice equiv. -> m/s ice equiv.

else
   n_larmip_region = 0
   ab_anom_larmip  = 0.0_dp
end if

#endif

!-------- Read data for z_sl --------

#if (SEA_LEVEL==3)

filename_with_path = trim(IN_PATH)//'/general/'//trim(SEA_LEVEL_FILE)

open(21, iostat=ios, file=trim(filename_with_path), status='old')

if (ios /= 0) then
   errormsg = ' >>> sico_init: Error when opening the data file for z_sl!'
   call error(errormsg)
end if
#if !defined(ALLOW_TAPENADE)
read(21, fmt=*) ch_dummy, specmap_time_min, specmap_time_stp, specmap_time_max

if (ch_dummy /= '#') then
   errormsg = ' >>> sico_init:' &
            //         end_of_line &
            //'        specmap_time_min, specmap_time_stp, specmap_time_max' &
            //         end_of_line &
            //'        not defined in data file for z_sl!'
   call error(errormsg)
end if

ndata_specmap = (specmap_time_max-specmap_time_min)/specmap_time_stp

allocate(specmap_zsl(0:ndata_specmap))
#else 
read(21, fmt=*)
#endif
do n=0, ndata_specmap
   read(21, fmt=*) d_dummy, specmap_zsl(n)
end do

close(21, status='keep')

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

!-------- Reading of the reference ice thickness for the
!                                           ice-shelf collapse masks --------

#if (ICE_SHELF_COLLAPSE_MASK==1)

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

!-------- Definition of initial values --------

!  ------ Present topography

#if (ANF_DAT==1)

call topography1(dxi, deta)

z_sl      = -1.11e+11_dp   ! dummy values for initial call
z_sl_mean = -1.11e+11_dp   ! of subroutine boundary

call boundary(time_init, dtime, dxi, deta, &
              delta_ts, glac_index, z_mar)

#if !defined(ALLOW_TAPENADE) /* Normal */

where ((mask==0).or.(mask==3))
                 ! grounded or floating ice
   as_perp_apl = as_perp
elsewhere        ! mask==1 or 2, ice-free land or sea
   as_perp_apl = 0.0_dp
end where

#else /* Tapenade */

do j=0,JMAX
do i=0,IMAX
  if ((mask(j,i)==0).or.(mask(j,i)==3)) then
    as_perp_apl(j,i) = as_perp(j,i)
  else
    as_perp_apl(j,i) = 0.0_dp
  end if
end do
end do

#endif /* Normal vs. Tapenade */

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

call boundary(time_init, dtime, dxi, deta, &
              delta_ts, glac_index, z_mar)

#if !defined(ALLOW_TAPENADE) /* Normal */

where ((mask==0).or.(mask==3))
                 ! grounded or floating ice
   as_perp_apl = as_perp
elsewhere        ! mask==1 or 2, ice-free land or sea
   as_perp_apl = 0.0_dp
end where

#else /* Tapenade */

do j=0,JMAX
do i=0,IMAX
  if ((mask(j,i)==0).or.(mask(j,i)==3)) then
    as_perp_apl(j,i) = as_perp(j,i)
  else
    as_perp_apl(j,i) = 0.0_dp
  end if
end do
end do

#endif /* Normal vs. Tapenade */

smb_corr = 0.0_dp

Q_b_tot = Q_bm + Q_tld

#if (ENHMOD==1)
   call calc_enhance_1()
#elif (ENHMOD==2)
   call calc_enhance_2()
#elif (ENHMOD==3)
   call calc_enhance_3(time_init)
#elif (ENHMOD==4)
   !!! call calc_enhance_4()
   !!!    (for anisotropic flow enhancement factor,
   !!!     use values read from output of previous simulation)
#elif (ENHMOD==5)
   !!! call calc_enhance_5()
   !!!    (for anisotropic flow enhancement factor,
   !!!     use values read from output of previous simulation)
#else
   errormsg = ' >>> sico_init: Parameter ENHMOD must be between 1 and 5!'
   call error(errormsg)
#endif

#endif

!-------- Inner-point flag --------

flag_inner_point = .true.

flag_inner_point(0,:)    = .false.
flag_inner_point(JMAX,:) = .false.

flag_inner_point(:,0)    = .false.
flag_inner_point(:,IMAX) = .false.

!-------- Distance between grid points with delta_i=ir, delta_j=jr --------

#if (GRID==0 || GRID==1)

do ir=-IMAX, IMAX
do jr=-JMAX, JMAX
   dist_dxdy(jr,ir) = sqrt( (real(ir,dp)*dxi)**2 + (real(jr,dp)*deta)**2 )
                  ! distortion due to stereographic projection not accounted for
end do
end do

#endif

!-------- Initial velocities --------

call calc_temp_melt()
call flag_update_gf_gl_cf()
call calc_vxy_b_init()
call calc_dzs_dxy_aux(dxi, deta)

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

!  ------ Time-series file for the ice sheet on the whole

#if !defined(ALLOW_TAPENADE) /* Normal */

filename_with_path = trim(OUT_PATH)//'/'//trim(run_name)//'.ser'

open(12, iostat=ios, file=trim(filename_with_path), status='new')

if (ios /= 0) then
   errormsg = ' >>> sico_init: Error when opening the ser file!'
   call error(errormsg)
end if

if ((forcing_flag == 1).or.(forcing_flag == 3)) then

   write(12,1102)
   write(12,1103)

   1102 format('         t(a)  D_Ts(deg C) z_sl_mean(m)',/, &
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

#else /* Tapenade */

print *, ' >>> sico_init: not writing to the ser file in '
print *, '                adjoint applications!' 

#endif /* Normal vs. Tapenade */

!  ------ Time-series file for deep boreholes

#if !defined(ALLOW_TAPENADE) /* Normal */

n_core = 6   ! Vostok, Dome A, Dome C, Dome F, Kohnen, Byrd

allocate(lambda_core(n_core), phi_core(n_core), &
         x_core(n_core), y_core(n_core), ch_core(n_core))

ch_core(1)     = 'Vostok'
phi_core(1)    = -78.467_dp *deg2rad    ! Geographical position of Vostok,
lambda_core(1) = 106.8_dp   *deg2rad    ! conversion deg -> rad
 
ch_core(2)     = 'Dome A'
phi_core(2)    = -80.367_dp *deg2rad    ! Geographical position of Dome A,
lambda_core(2) =  77.35_dp  *deg2rad    ! conversion deg -> rad

ch_core(3)     = 'Dome C'
phi_core(3)    = -75.1_dp *deg2rad      ! Geographical position of Dome C,
lambda_core(3) = 123.4_dp *deg2rad      ! conversion deg -> rad

ch_core(4)     = 'Dome F'
phi_core(4)    = -77.317_dp *deg2rad    ! Geographical position of Dome F,
lambda_core(4) =  39.7_dp   *deg2rad    ! conversion deg -> rad

ch_core(5)     = 'Kohnen'
phi_core(5)    = -75.0_dp   *deg2rad    ! Geographical position of Kohnen,
lambda_core(5) =   0.067_dp *deg2rad    ! conversion deg -> rad

ch_core(6)     = 'Byrd'
phi_core(6)    =  -80.017_dp *deg2rad   ! Geographical position of Byrd,
lambda_core(6) = -119.517_dp *deg2rad   ! conversion deg -> rad

#if (GRID==0 || GRID==1)   /* Stereographic projection */

do n=1, n_core

   if (F_INV > 1.0e+10_dp) then   ! interpreted as infinity, thus no flattening
                                  ! (spherical planet)

      call stereo_forw_sphere(lambda_core(n), phi_core(n), &
                              R, LAMBDA0, PHI0, x_core(n), y_core(n))

   else   ! finite inverse flattening (ellipsoidal planet)

      call stereo_forw_ellipsoid(lambda_core(n), phi_core(n), &
                                 A, B, LAMBDA0, PHI0, x_core(n), y_core(n))

   end if

end do

#elif (GRID==2)   /* Geographical coordinates */

x_core = lambda_core
y_core = phi_core

#endif

#endif /* Normal (Tapenade: No core data for adjoint) */

filename_with_path = trim(OUT_PATH)//'/'//trim(run_name)//'.core'

#if !defined(ALLOW_TAPENADE) /* Normal */
open(14, iostat=ios, file=trim(filename_with_path), status='new')
#else /* Tapenade */
open(14, iostat=ios, file=trim(filename_with_path))
#endif /* Normal vs. Tapenade */

if (ios /= 0) then
   errormsg = ' >>> sico_init: Error when opening the core file!'
   call error(errormsg)
end if

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
end if

#else /* Tapenade */

print *, ' >>> sico_init: not producing initial, typical outputs'
print *, '                in adjoint mode.'

#endif /* Normal vs. Tapenade */

end subroutine sico_init

!-------------------------------------------------------------------------------
!> Definition of the initial surface and bedrock topography
!! (including gradients) and of the horizontal grid spacings dxi, deta.
!! For present-day initial topography.
!<------------------------------------------------------------------------------
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

dxi  = DX *1000.0_dp   ! km -> m
deta = DX *1000.0_dp   ! km -> m

xi0  = X0 *1000.0_dp   ! km -> m
eta0 = Y0 *1000.0_dp   ! km -> m

freeboard_ratio = (RHO_SW-RHO)/RHO_SW

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i) <= 1) then

      zb(j,i) = zl(j,i)   ! ensure consistency

      if (zs(j,i) < zb(j,i)) then
         zs(j,i) = zb(j,i)
         mask(j,i) = 1
      end if

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
      if (H_ice < 0.0_dp) then
         H_ice = 0.0_dp
         mask(j,i) = 2
      end if
      zs(j,i) = zl(j,i)+H_ice
      zb(j,i) = zl(j,i)
#elif (MARGIN==3)
      H_ice = zs(j,i)-zb(j,i)   ! ice thickness
      if (H_ice < 0.0_dp) then
         H_ice = 0.0_dp
         mask(j,i) = 2
      end if
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

if ( trim(adjustl(MASK_REGION_FILE)) /= 'none' ) then
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

if (mask_region(0,0) == -1) then   ! set default values for mask_region

   do i=0, IMAX
   do j=0, JMAX

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

end if

end subroutine topography1

!-------------------------------------------------------------------------------
!> Definition of the initial surface and bedrock topography
!! (including gradients) and of the horizontal grid spacings dxi, deta.
!! For ice-free initial topography with relaxed lithosphere.
!<------------------------------------------------------------------------------
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

dxi  = DX *1000.0_dp   ! km -> m
deta = DX *1000.0_dp   ! km -> m

xi0  = X0 *1000.0_dp   ! km -> m
eta0 = Y0 *1000.0_dp   ! km -> m

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

if ( trim(adjustl(MASK_REGION_FILE)) /= 'none' ) then
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

if (mask_region(0,0) == -1) then   ! set default values for mask_region

   do i=0, IMAX
   do j=0, JMAX

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

end if

end subroutine topography2

!-------------------------------------------------------------------------------
!> Definition of the initial surface and bedrock topography
!! (including gradients) and of the horizontal grid spacings dxi, deta.
!! For initial topography from previous simulation.
!<------------------------------------------------------------------------------
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

dxi  = DX *1000.0_dp   ! km -> m
deta = DX *1000.0_dp   ! km -> m

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

if ( trim(adjustl(MASK_REGION_FILE)) /= 'none' ) then
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

if (mask_region(0,0) == -1) then   ! set default values for mask_region

   do i=0, IMAX
   do j=0, JMAX

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

end if

end subroutine topography3

!-------------------------------------------------------------------------------

end module sico_init_m
!
