!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ i n i t _ m
!
!! SCAND domain: Initializations for SICOPOLIS.
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
!> SCAND domain: Initializations for SICOPOLIS.
!-------------------------------------------------------------------------------
module sico_init_m

  use sico_types_m
  use sico_variables_m
  use sico_vars_m
  use error_m

  implicit none

  public

contains

!-------------------------------------------------------------------------------
!> Main routine of sico_init_m: Initialisations for SICOPOLIS.
!-------------------------------------------------------------------------------
subroutine sico_init(delta_ts, glac_index, &
               mean_accum, &
               dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser, &
               time, time_init, time_end, time_output, &
               dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
               z_mar, &
               ndat2d, ndat3d, n_output)

  use compare_float_m
  use ice_material_properties_m, only : ice_mat_eqs_pars
  use enth_temp_omega_m, only : calc_c_int_table, calc_c_int_inv_table, &
                                enth_fct_temp_omega

  use read_m, only : read_scalar_input, read_2d_input, read_kei, read_phys_para

  use boundary_m
  use init_temp_water_age_m
  use calc_enhance_m
  use flag_update_gf_gl_cf_m
  use calc_vxy_m
  use calc_vz_m
  use calc_dxyz_m
  use calc_temp_melt_bas_m

  use output_m

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
integer(i4b)       :: ios
integer(i4b)       :: istat, ierr
integer(i4b)       :: n_q_geo_mod
integer(i4b), dimension(0:JMAX,0:IMAX) :: mask_ref
real(dp)           :: dtime0, dtime_temp0, dtime_wss0, dtime_out0, dtime_ser0
real(dp)           :: time_init0, time_end0
#if (OUTPUT==2 || OUTPUT==3)
real(dp)           :: time_output0(N_OUTPUT)
#endif
real(dp)           :: d_dummy
character(len=256) :: anfdatname
character(len=256) :: filename_with_path
character(len=256) :: shell_command
character(len=256) :: ch_revision
character(len= 64) :: ch_var_name
character(len=  3) :: ch_month(12)
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

#if (MARGIN==3 || DYNAMICS==2)
  call lis_initialize(ierr)
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

#if (!defined(CHECK_RES_IMAX_JMAX) || CHECK_RES_IMAX_JMAX==1)

#if (GRID==0 || GRID==1)

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

#elif (GRID==2)

   errormsg = ' >>> sico_init: ' &
                 //'GRID==2 not allowed for the Scandinavia application!'
   call error(errormsg)

#endif

#else /* CHECK_RES_IMAX_JMAX==0 */

write(6, fmt='(a)') ' >>> sico_init: CHECK_RES_IMAX_JMAX==0'
write(6, fmt='(a)') '      -> compatibility check between horizontal resolution'
write(6, fmt='(a)') '         and number of grid points not performed.'
write(6, fmt='(a)') ' '

#endif /* CHECK_RES_IMAX_JMAX */

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

!-------- Compatibility check of surface-temperature and precipitation
!         determination by interpolation between present and LGM values
!         with a glacial index --------

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

open(10, iostat=ios, file=trim(filename_with_path), status='new')

if (ios /= 0) then
   errormsg = ' >>> sico_init: Error when opening the log file!'
   call error(errormsg)
end if

write(10, fmt=trim(fmt1)) 'Computational domain:'
write(10, fmt=trim(fmt1)) trim(ch_domain_long)
write(10, fmt=trim(fmt1)) ' '

write(10, fmt=trim(fmt1)) 'Physical-parameter file = ' &
                          // trim(adjustl(PHYS_PARA_FILE))
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
write(10, fmt=trim(fmt3)) 'x0      =', X0
write(10, fmt=trim(fmt3)) 'y0      =', Y0
write(10, fmt=trim(fmt3)) 'dx      =', DX
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

write(10, fmt=trim(fmt2)) 'ANF_DAT = ', ANF_DAT
write(10, fmt=trim(fmt1)) 'zs_present file   = '//ZS_PRESENT_FILE
#if (ANF_DAT==1)
#if (defined(ZB_PRESENT_FILE))
write(10, fmt=trim(fmt1)) 'zb_present file   = '//ZB_PRESENT_FILE
#endif
write(10, fmt=trim(fmt1)) 'zl_present file   = '//ZL_PRESENT_FILE
#endif
write(10, fmt=trim(fmt1)) 'zl0 file          = '//ZL0_FILE
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
write(10, fmt=trim(fmt1)) ' '

write(10, fmt=trim(fmt1)) 'Physical-parameter file = '//PHYS_PARA_FILE
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

#if (CALCTHK==2)
write(10, fmt=trim(fmt3))  'ovi_weight   =', OVI_WEIGHT
write(10, fmt=trim(fmt3))  'omega_sor    =', OMEGA_SOR
#if (ITER_MAX_SOR>0)
write(10, fmt=trim(fmt2)) 'iter_max_sor = ', ITER_MAX_SOR
#endif
#endif

write(10, fmt=trim(fmt1)) ' '

write(10, fmt=trim(fmt1)) 'temp_mm_present file = '//TEMP_MM_PRESENT_FILE
#if (TSURFACE==1)
write(10, fmt=trim(fmt3)) 'delta_ts0      =', DELTA_TS0
#elif (TSURFACE==3)
write(10, fmt=trim(fmt3)) 'sine_amplit    =', SINE_AMPLIT
write(10, fmt=trim(fmt3)) 'sine_period    =', SINE_PERIOD
#elif (TSURFACE==4)
write(10, fmt=trim(fmt1)) 'GRIP file      = '//GRIP_TEMP_FILE
write(10, fmt=trim(fmt3)) 'grip_temp_fact =', GRIP_TEMP_FACT
#elif (TSURFACE==5)
write(10, fmt=trim(fmt1)) 'Glacial-index file = '//GLAC_IND_FILE
write(10, fmt=trim(fmt1)) 'temp_mm_anom file  = '//TEMP_MM_ANOM_FILE
write(10, fmt=trim(fmt3)) 'temp_mm_anom fact  = ', TEMP_MM_ANOM_FACT
#endif

write(10, fmt=trim(fmt1)) 'precip_mm_present file = '//PRECIP_MM_PRESENT_FILE
#if (ACCSURFACE==1)
write(10, fmt=trim(fmt3)) 'accfact        =', ACCFACT
#elif (ACCSURFACE==2 || ACCSURFACE==3)
write(10, fmt=trim(fmt3)) 'gamma_s        =', GAMMA_S
#elif (ACCSURFACE==5)
write(10, fmt=trim(fmt1)) 'precip_mm_anom file    = '//PRECIP_MM_ANOM_FILE
write(10, fmt=trim(fmt3)) 'precip_mm_anom fact    = ', PRECIP_MM_ANOM_FACT
#endif
#if (ACCSURFACE <= 3)
write(10, fmt=trim(fmt2)) 'ELEV_DESERT = ', ELEV_DESERT
#if (ELEV_DESERT == 1)
write(10, fmt=trim(fmt3)) 'gamma_p     =', GAMMA_P
write(10, fmt=trim(fmt3)) 'zs_thresh   =', ZS_THRESH
#endif
#endif

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

write(10, fmt=trim(fmt2)) 'SEA_LEVEL  = ', SEA_LEVEL
#if (SEA_LEVEL==1)
write(10, fmt=trim(fmt3)) 'z_sl0          =', Z_SL0
#elif (SEA_LEVEL==3)
write(10, fmt=trim(fmt1)) 'sea-level file = '//SEA_LEVEL_FILE
#endif
write(10, fmt=trim(fmt1)) ' '

#if (MARGIN==2)
#if (MARINE_ICE_CALVING==2 || MARINE_ICE_CALVING==3)
write(10, fmt=trim(fmt3)) 'z_mar          =', Z_MAR
write(10, fmt=trim(fmt1)) ' '
#elif (MARINE_ICE_CALVING==4 || MARINE_ICE_CALVING==5 || MARINE_ICE_CALVING==6 || MARINE_ICE_CALVING==7)
write(10, fmt=trim(fmt3)) 'fact_z_mar     =', FACT_Z_MAR
write(10, fmt=trim(fmt1)) ' '
#elif (MARINE_ICE_FORMATION==2 && MARINE_ICE_CALVING==9)
write(10, fmt=trim(fmt3)) 'calv_uw_coeff  =', CALV_UW_COEFF
write(10, fmt=trim(fmt3)) 'r1_calv_uw     =', R1_CALV_UW
write(10, fmt=trim(fmt3)) 'r2_calv_uw     =', R2_CALV_UW
write(10, fmt=trim(fmt1)) ' '
#endif
#elif (MARGIN==3)
#if (ICE_SHELF_CALVING==2)
write(10, fmt=trim(fmt3)) 'H_calv          =', H_CALV
write(10, fmt=trim(fmt1)) ' '
#endif
#endif

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
#if (BASAL_HYDROLOGY==1 && defined(HYDRO_SLIDE_SAT_FCT) && defined(C_HW_SLIDE) && defined(HW0_SLIDE))
write(10, fmt=trim(fmt2)) 'HYDRO_SLIDE_SAT_FCT = ', HYDRO_SLIDE_SAT_FCT
write(10, fmt=trim(fmt3)) 'c_Hw_slide  =', C_HW_SLIDE
write(10, fmt=trim(fmt3)) 'Hw0_slide   =', HW0_SLIDE
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
write(10, fmt=trim(fmt3)) 'qbm_marine               =', QBM_MARINE
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
#endif

write(10, fmt=trim(fmt2)) 'REBOUND       = ', REBOUND
#if (REBOUND==1)
write(10, fmt=trim(fmt3)) 'frac_llra     =', FRAC_LLRA
#endif
#if (REBOUND==1 || REBOUND==2)
write(10, fmt=trim(fmt2)) 'TIME_LAG_MOD  = ', TIME_LAG_MOD
#if (TIME_LAG_MOD==1)
write(10, fmt=trim(fmt3)) 'time_lag      =', TIME_LAG
#elif (TIME_LAG_MOD==2)
write(10, fmt=trim(fmt1)) 'time_lag_file = '//TIME_LAG_FILE
#else
errormsg = ' >>> sico_init: TIME_LAG_MOD must be either 1 or 2!'
call error(errormsg)
#endif
#endif
#if (REBOUND==2)
write(10, fmt=trim(fmt2)) 'FLEX_RIG_MOD  = ', FLEX_RIG_MOD
#if (FLEX_RIG_MOD==1)
write(10, fmt=trim(fmt3)) 'flex_rig      =', FLEX_RIG
#elif (FLEX_RIG_MOD==2)
write(10, fmt=trim(fmt1)) 'flex_rig_file = '//FLEX_RIG_FILE
#else
errormsg = ' >>> sico_init: FLEX_RIG_MOD must be either 1 or 2!'
call error(errormsg)
#endif
#endif
write(10, fmt=trim(fmt1)) ' '

#if (FLOW_LAW==2)
write(10, fmt=trim(fmt3)) 'gr_size   =', GR_SIZE
write(10, fmt=trim(fmt1)) ' '
#endif
#if (FIN_VISC==2)
write(10, fmt=trim(fmt3)) 'sigma_res =', SIGMA_RES
write(10, fmt=trim(fmt1)) ' '
#endif

write(10, fmt=trim(fmt2)) 'ENHMOD = ', ENHMOD
#if (ENHMOD==1 || ENHMOD==2 || ENHMOD==3)
write(10, fmt=trim(fmt3)) 'enh_fact    =', ENH_FACT
#endif
#if (ENHMOD==2 || ENHMOD==3)
write(10, fmt=trim(fmt3)) 'enh_intg    =', ENH_INTG
#endif
#if (ENHMOD==2)
write(10, fmt=trim(fmt3)) 'age_trans   =', AGE_TRANS_0
#endif
#if (ENHMOD==3)
write(10, fmt=trim(fmt3)) 'date_trans1 =', DATE_TRANS1_0
write(10, fmt=trim(fmt3)) 'date_trans2 =', DATE_TRANS2_0
write(10, fmt=trim(fmt3)) 'date_trans3 =', DATE_TRANS3_0
#endif
#if (ENHMOD==4 || ENHMOD==5)
write(10, fmt=trim(fmt3)) 'enh_compr   =', ENH_COMPR
write(10, fmt=trim(fmt3)) 'enh_shear   =', ENH_SHEAR
#endif
#if (DYNAMICS==2 && defined(ENH_STREAM))
if (ENH_STREAM >= 0.0_dp) &
   write(10, fmt=trim(fmt3)) 'enh_stream =', ENH_STREAM
#endif
#if ((ENHMOD==1 || ENHMOD==2 || ENHMOD==3 || ENHMOD==4) && MARGIN==3)
write(10, fmt=trim(fmt3)) 'enh_shelf   =', ENH_SHELF
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

write(10, fmt=trim(fmt2)) 'CALCMOD    = ', CALCMOD
#if (CALCMOD==-1 && defined(TEMP_CONST))
write(10, fmt=trim(fmt3)) 'TEMP_CONST =', TEMP_CONST
#endif
#if (CALCMOD==-1 && defined(AGE_CONST))
write(10, fmt=trim(fmt3)) 'AGE_CONST  =', AGE_CONST
#endif
#if (CALCMOD==1 && defined(CTS_MELTING_FREEZING))
write(10, fmt=trim(fmt2)) 'CTS_MELTING_FREEZING = ', CTS_MELTING_FREEZING
#endif
write(10, fmt=trim(fmt2)) 'FLOW_LAW   = ', FLOW_LAW
write(10, fmt=trim(fmt2)) 'FIN_VISC   = ', FIN_VISC
write(10, fmt=trim(fmt2)) 'MARGIN     = ', MARGIN
#if (MARGIN==2)
write(10, fmt=trim(fmt2)) 'MARINE_ICE_FORMATION = ', MARINE_ICE_FORMATION
write(10, fmt=trim(fmt2)) 'MARINE_ICE_CALVING   = ', MARINE_ICE_CALVING
#elif (MARGIN==3)
write(10, fmt=trim(fmt2)) 'ICE_SHELF_CALVING = ', ICE_SHELF_CALVING
#endif
write(10, fmt=trim(fmt2)) 'ADV_HOR    = ', ADV_HOR
write(10, fmt=trim(fmt2)) 'ADV_VERT   = ', ADV_VERT
write(10, fmt=trim(fmt2)) 'TOPOGRAD   = ', TOPOGRAD
#if (MARGIN==3 && defined(GL_SURF_GRAD))
write(10, fmt=trim(fmt2)) 'GL_SURF_GRAD = ', GL_SURF_GRAD
#endif
write(10, fmt=trim(fmt2)) 'TSURFACE   = ', TSURFACE
write(10, fmt=trim(fmt2)) 'ACCSURFACE = ', ACCSURFACE
#if (ACCSURFACE==5)
write(10, fmt=trim(fmt2)) 'PRECIP_ANOM_INTERPOL = ', PRECIP_ANOM_INTERPOL
#endif
write(10, fmt=trim(fmt2)) 'SOLID_PRECIP = ', SOLID_PRECIP
write(10, fmt=trim(fmt2)) 'ABLSURFACE = ', ABLSURFACE
#if (defined(MB_ACCOUNT))
write(10, fmt=trim(fmt2)) 'MB_ACCOUNT = ', MB_ACCOUNT
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
write(10, fmt=trim(fmt2)) 'OUT_TIMES   = ', OUT_TIMES
#endif
#if (defined(OUTPUT_INIT))
write(10, fmt=trim(fmt2)) 'OUTPUT_INIT = ', OUTPUT_INIT
#endif
write(10, fmt=trim(fmt2)) 'OUTPUT      = ', OUTPUT
#if (OUTPUT==1 || OUTPUT==3)
write(10, fmt=trim(fmt3))  'dtime_out   =' , dtime_out0
#endif
write(10, fmt=trim(fmt3))  'dtime_ser   =' , dtime_ser0
#if (OUTPUT==1 || OUTPUT==2)
write(10, fmt=trim(fmt2)) 'ERGDAT      = ', ERGDAT
#endif
#if (defined(OUTPUT_FLUX_VARS))
write(10, fmt=trim(fmt2)) 'OUTPUT_FLUX_VARS = ', OUTPUT_FLUX_VARS
#endif
#if (OUTPUT==2 || OUTPUT==3)
write(10, fmt=trim(fmt2)) 'n_output    = ', n_output
do n=1, n_output
   if (n==1) then
      write(10, fmt=trim(fmt3))  'time_output =' , time_output0(n)
   else
      write(10, fmt=trim(fmt3))  '             ' , time_output0(n)
   end if
end do
#endif
write(10, fmt=trim(fmt1)) ' '

call get_environment_variable(name='REPO_REVISION', value=ch_revision, &
                              status=istat, trim_name=.true.)
write(10, fmt=trim(fmt1)) 'Program version and date: '//VERSION//' / '//DATE
write(10, fmt=trim(fmt1)) 'Git revision identifier : ' // trim(ch_revision)

close(10, status='keep')

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

time = time_init

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

!-------- Reading of measurements for present monthly-mean precipitation --------

#if (GRID==0 || GRID==1)

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(PRECIP_MM_PRESENT_FILE)

ch_month = [ 'jan', 'feb', 'mar', 'apr', 'may', 'jun', &
             'jul', 'aug', 'sep', 'oct', 'nov', 'dec' ]

do n=1, 12   ! month counter

   ch_var_name = 'precip_present_' // trim(ch_month(n))

   call read_2d_input(filename_with_path, &
                      ch_var_name=trim(ch_var_name), &
                      n_var_type=1, n_ascii_header=6+3*n+(JMAX+1)*(n-1), &
                      field2d_r=field2d_aux)

   precip_present(:,:,n) = field2d_aux *(1.0e-03_dp*sec2year)*(RHO_W/RHO)
                                        ! mm/a water equiv. -> m/s ice equiv.

end do

#elif (GRID==2)

errormsg = ' >>> sico_init: ' &
              //'GRID==2 not allowed for the Scandinavia application!'
call error(errormsg)

#endif

!-------- Reading of LGM monthly-mean precipitation-rate anomalies --------

#if (ACCSURFACE==5)

#if (GRID==0 || GRID==1)

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(PRECIP_MM_ANOM_FILE)

ch_month = [ 'jan', 'feb', 'mar', 'apr', 'may', 'jun', &
             'jul', 'aug', 'sep', 'oct', 'nov', 'dec' ]

do n=1, 12   ! month counter

   ch_var_name = 'precip_lgm_anom_' // trim(ch_month(n))

   call read_2d_input(filename_with_path, &
                      ch_var_name=trim(ch_var_name), &
                      n_var_type=1, n_ascii_header=6+3*n+(JMAX+1)*(n-1), &
                      field2d_r=field2d_aux)

   precip_lgm_anom(:,:,n) = field2d_aux

end do

precip_lgm_anom = precip_lgm_anom * PRECIP_MM_ANOM_FACT

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

#endif

!-------- Mean accumulation --------

mean_accum = MEAN_ACCUM*(1.0e-03_dp*sec2year)*(RHO_W/RHO)
                       ! mm/a water equiv. -> m/s ice equiv.

!-------- Reading of present topography mask --------

#if (GRID==0 || GRID==1)

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(MASK_PRESENT_FILE)

call read_2d_input(filename_with_path, &
                   ch_var_name='mask', n_var_type=3, n_ascii_header=6, &
                   field2d_r=field2d_aux)

mask_ref = nint(field2d_aux)

#elif (GRID==2)

errormsg = ' >>> sico_init: ' &
              //'GRID==2 not allowed for the Scandinavia application!'
call error(errormsg)

#endif

!-------- Reading of data for present monthly-mean surface temperature --------

#if (GRID==0 || GRID==1)

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(TEMP_MM_PRESENT_FILE)

ch_month = [ 'jan', 'feb', 'mar', 'apr', 'may', 'jun', &
             'jul', 'aug', 'sep', 'oct', 'nov', 'dec' ]

do n=1, 12   ! month counter

   ch_var_name = 'temp_present_' // trim(ch_month(n))

   call read_2d_input(filename_with_path, &
                      ch_var_name=trim(ch_var_name), &
                      n_var_type=1, n_ascii_header=6+3*n+(JMAX+1)*(n-1), &
                      field2d_r=field2d_aux)

   temp_mm_present(:,:,n) = field2d_aux

end do

#elif (GRID==2)

errormsg = ' >>> sico_init: ' &
              //'GRID==2 not allowed for the Scandinavia application!'
call error(errormsg)

#endif

!-------- Reading of LGM monthly-mean surface-temperature anomalies --------

#if (TSURFACE==5)

#if (GRID==0 || GRID==1)

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(TEMP_MM_ANOM_FILE)

ch_month = [ 'jan', 'feb', 'mar', 'apr', 'may', 'jun', &
             'jul', 'aug', 'sep', 'oct', 'nov', 'dec' ]

do n=1, 12   ! month counter

   ch_var_name = 'temp_lgm_anom_' // trim(ch_month(n))

   call read_2d_input(filename_with_path, &
                      ch_var_name=trim(ch_var_name), &
                      n_var_type=1, n_ascii_header=6+3*n+(JMAX+1)*(n-1), &
                      field2d_r=field2d_aux)

   temp_mm_lgm_anom(:,:,n) = field2d_aux

end do

temp_mm_lgm_anom = temp_mm_lgm_anom * TEMP_MM_ANOM_FACT

#endif

#endif

!-------- Present reference elevation
!         (for precipitation and surface-temperature data) --------

#if (GRID==0 || GRID==1)

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(ZS_PRESENT_FILE)

call read_2d_input(filename_with_path, &
                   ch_var_name='zs', n_var_type=1, n_ascii_header=6, &
                   field2d_r=field2d_aux)

zs_ref = field2d_aux

do i=0, IMAX
do j=0, JMAX
   if (mask_ref(j,i) >= 2) zs_ref(j,i) = 0.0_dp
                 ! resetting elevations over the ocean
                 ! to the present-day sea surface
end do
end do

#elif (GRID==2)

errormsg = ' >>> sico_init: ' &
              //'GRID==2 not allowed for the Scandinavia application!'
call error(errormsg)

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

!-------- Definition of initial values --------

!  ------ Present topography

#if (ANF_DAT==1)

errormsg = ' >>> sico_init: ANF_DAT==1 not allowed for scand application!'
call error(errormsg)

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

where ((mask==0).or.(mask==3))
                 ! grounded or floating ice
   as_perp_apl = as_perp
elsewhere        ! mask==1 or 2, ice-free land or sea
   as_perp_apl = 0.0_dp
end where

smb_corr = 0.0_dp

Q_b_tot = Q_bm + Q_tld

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

filename_with_path = trim(OUT_PATH)//'/'//trim(run_name)//'.ser'

open(12, iostat=ios, file=trim(filename_with_path), status='new')

if (ios /= 0) then
   errormsg = ' >>> sico_init: Error when opening the ser file!'
   call error(errormsg)
end if

if (forcing_flag == 1) then

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

!  ------ Time-series file for deep boreholes

n_core = 0   ! No boreholes defined

if (n_core > n_core_max) then
   errormsg = ' >>> sico_init: n_core <= n_core_max required!' &
            //         end_of_line &
            //'        Increase value of n_core_max in sico_variables_m!'
   call error(errormsg)
end if

filename_with_path = trim(OUT_PATH)//'/'//trim(run_name)//'.core'

open(14, iostat=ios, file=trim(filename_with_path), status='new')

write(14,'(1x,a)') '---------------------'
write(14,'(1x,a)') 'No boreholes defined.'
write(14,'(1x,a)') '---------------------'

!-------- Output of the initial state --------

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

end subroutine sico_init

!-------------------------------------------------------------------------------
!> Definition of the initial surface and bedrock topography
!! (including gradients) and of the horizontal grid spacings dxi, deta.
!! For ice-free initial topography with relaxed lithosphere
!! (not defined for this domain, thus routine stops execution of SICOPOLIS).
!-------------------------------------------------------------------------------
subroutine topography1(dxi, deta)

implicit none

real(dp), intent(out) :: dxi, deta

! integer(i4b) :: i, j
! real(dp)     :: xi0, eta0

dxi=0.0_dp; deta=0.0_dp   ! dummy values

errormsg = ' >>> topography1: ANF_DAT==1 not defined for this domain!'
call error(errormsg)

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

if (mask_region(0,0) == -1) mask_region = 0   ! regions undefined

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

if (mask_region(0,0) == -1) mask_region = 0   ! regions undefined

end subroutine topography3

!-------------------------------------------------------------------------------

end module sico_init_m
!
