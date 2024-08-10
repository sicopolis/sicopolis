!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ i n i t _ m
!
!! HEINO domain: Initializations for SICOPOLIS.
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
!> HEINO domain: Initializations for SICOPOLIS.
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
!> Main routine of sico_init_m: Initializations for SICOPOLIS.
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

  use read_m, only : read_scalar_input, read_2d_input, read_phys_para

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
integer(i4b)       :: ios, ios1, ios2, ios3, ios4
integer(i4b)       :: istat, ierr
integer(i4b)       :: n_q_geo_mod
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

dtime       = 0.0_dp
dtime_temp  = 0.0_dp
dtime_wss   = 0.0_dp
dtime_out   = 0.0_dp
dtime_ser   = 0.0_dp

time        = 0.0_dp
time_init   = 0.0_dp
time_end    = 0.0_dp
time_output = 0.0_dp

!-------- Initialization of the Library of Iterative Solvers Lis,
!                                                     if required --------

#if (MARGIN==3 || DYNAMICS==2)
  call lis_initialize(ierr)
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

temp_min  = TEMP_MIN
s_t       = S_T       *1.0e-09_dp            ! K/km3 -> K/m3
x_hat     = X_HAT     *1.0e+03_dp            ! km -> m
y_hat     = Y_HAT     *1.0e+03_dp            ! km -> m
rad       = RAD       *1.0e+03_dp            ! km -> m
b_min     = B_MIN     *sec2year              ! m/a -> m/s
b_max     = B_MAX     *sec2year              ! m/a -> m/s

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

#if (GRID==0)

if (approx_equal(DX, 50.0_dp, eps_sp_dp)) then

   if ((IMAX /= 80).or.(JMAX /= 80)) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else if (approx_equal(DX, 25.0_dp, eps_sp_dp)) then

   if ((IMAX /= 160).or.(JMAX /= 160)) then
      errormsg = ' >>> sico_init: IMAX and/or JMAX wrong!'
      call error(errormsg)
   end if

else

   errormsg = ' >>> sico_init: DX wrong!'
   call error(errormsg)

end if

#elif (GRID==1)

errormsg = ' >>> sico_init: GRID==1 not allowed for this application!'
call error(errormsg)

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

forcing_flag = 1   ! forcing by delta_ts

!-------- Initialization of numerical time steps --------

dtime0      = DTIME0
dtime_temp0 = DTIME_TEMP0

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
write(10, fmt=trim(fmt3)) 'x0      =', X0
write(10, fmt=trim(fmt3)) 'y0      =', Y0
write(10, fmt=trim(fmt3)) 'dx      =', DX
#elif (GRID==2)
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
write(10, fmt=trim(fmt1)) ' '

write(10, fmt=trim(fmt2)) 'ANF_DAT = ', ANF_DAT
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

write(10, fmt=trim(fmt3)) 'temp_min =', TEMP_MIN
write(10, fmt=trim(fmt3)) 's_t      =', S_T
write(10, fmt=trim(fmt3)) 'x_hat    =', X_HAT
write(10, fmt=trim(fmt3)) 'y_hat    =', Y_HAT
write(10, fmt=trim(fmt3)) 'rad      =', RAD
write(10, fmt=trim(fmt3)) 'b_min    =', B_MIN
write(10, fmt=trim(fmt3)) 'b_max    =', B_MAX
write(10, fmt=trim(fmt1)) ' '

#if (TSURFACE==1)
write(10, fmt=trim(fmt3)) 'delta_ts0      =', DELTA_TS0
#elif (TSURFACE==3)
write(10, fmt=trim(fmt3)) 'sine_amplit    =', SINE_AMPLIT
write(10, fmt=trim(fmt3)) 'sine_period    =', SINE_PERIOD
#elif (TSURFACE==4)
write(10, fmt=trim(fmt1)) 'GRIP file      = '//GRIP_TEMP_FILE
write(10, fmt=trim(fmt3)) 'grip_temp_fact =', GRIP_TEMP_FACT
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
write(10, fmt=trim(fmt3)) 'c_Hw_slide  =', C_HW_SLIDE
write(10, fmt=trim(fmt3)) 'Hw0_slide   =', HW0_SLIDE
#endif

#endif

write(10, fmt=trim(fmt1)) ' '

if (n_q_geo_mod==1) then
   write(10, fmt=trim(fmt3)) 'q_geo =', Q_GEO
end if
write(10, fmt=trim(fmt2)) 'Q_LITHO = ', Q_LITHO
write(10, fmt=trim(fmt1)) ' '

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

!-------- Mean accumulation --------

mean_accum = MEAN_ACCUM*(1.0e-03_dp*sec2year)*(RHO_W/RHO)
!                      ! mm/a water equiv. -> m/s ice equiv.

!-------- Read file defining the regions for the sliding laws --------

#if (!defined(N_SLIDE_REGIONS) || N_SLIDE_REGIONS<=1)

n_slide_region = 1

#else

filename_with_path = trim(IN_PATH)//'/'//trim(ch_domain_short)//'/'// &
                     trim(SLIDE_REGIONS_FILE)

call read_2d_input(filename_with_path, &
                   ch_var_name='n_basin', n_var_type=3, n_ascii_header=6, &
                   field2d_r=field2d_aux)

n_slide_region = nint(field2d_aux)

#endif

!-------- Read data for delta_ts --------

#if (TSURFACE==4)

filename_with_path = trim(IN_PATH)//'/general/'//trim(GRIP_TEMP_FILE)

call read_scalar_input(filename_with_path, &
                       'delta_ts', ndata_grip_max, &
                       grip_time_min, grip_time_stp, grip_time_max, &
                       ndata_grip, griptemp)

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

   errormsg = ' >>> sico_init: ' &
                 //'Option Q_GEO_FILE not available for this application!'
   call error(errormsg)

end if

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

!  ------ Initial topography with a thin ice layer everywhere

#if (ANF_DAT==1)

call topography1(dxi, deta)

z_sl      = -1.11e+11_dp   ! dummy values for initial call
z_sl_mean = -1.11e+11_dp   ! of subroutine boundary

call boundary(time_init, dtime, dxi, deta, &
              delta_ts, glac_index, z_mar)

do i=0, IMAX
do j=0, JMAX

   if ((mask(j,i)==0).or.(mask(j,i)==3)) then
                    ! grounded or floating ice
      as_perp_apl(j,i) = as_perp(j,i)
   else             ! mask==1 or 2, ice-free land or sea
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

do i=0, IMAX
do j=0, JMAX

   if ((mask(j,i)==0).or.(mask(j,i)==3)) then
                    ! grounded or floating ice
      as_perp_apl(j,i) = as_perp(j,i)
   else             ! mask==1 or 2, ice-free land or sea
      as_perp_apl(j,i) = 0.0_dp
   end if

end do
end do

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

!  ------ Time-series file for the sediment area ("Hudson Bay, Hudson Strait")

filename_with_path = trim(OUT_PATH)//'/'//trim(run_name)//'.sed'

open(15, iostat=ios, file=trim(filename_with_path), status='new')

if (ios /= 0) then
   errormsg = ' >>> sico_init: Error when opening the sed file!'
   call error(errormsg)
end if

write(15,1108)
write(15,1109)

1108 format('         t(a)   D_Ts(degC) z_sl_mean(m)',/, &
            '                  H_ave(m)   Tbh_ave(C)     Atb(m^2)')
1109 format('----------------------------------------------------')

!  ------ Time-series file for specified sites (i.e., ice cores)

n_site = 7   ! Points P1 - P7

if (n_site > n_site_max) then
   errormsg = ' >>> sico_init: n_site <= n_site_max required!' &
            //         end_of_line &
            //'        Increase value of n_site_max in sico_variables_m!'
   call error(errormsg)
end if

ch_site(1)     = 'P1'
lambda_site(1) =    0.0_dp  ! dummy
phi_site(1)    =    0.0_dp  ! dummy
x_site(1)      = 3900.0_dp *1.0e+03_dp    ! Point P1,
y_site(1)      = 2000.0_dp *1.0e+03_dp    ! conversion km -> m

ch_site(2)     = 'P2'
lambda_site(2) =    0.0_dp  ! dummy
phi_site(2)    =    0.0_dp  ! dummy
x_site(2)      = 3800.0_dp *1.0e+03_dp    ! Point P2,
y_site(2)      = 2000.0_dp *1.0e+03_dp    ! conversion km -> m

ch_site(3)     = 'P3'
lambda_site(3) =    0.0_dp  ! dummy
phi_site(3)    =    0.0_dp  ! dummy
x_site(3)      = 3700.0_dp *1.0e+03_dp    ! Point P3,
y_site(3)      = 2000.0_dp *1.0e+03_dp    ! conversion km -> m

ch_site(4)     = 'P4'
lambda_site(4) =    0.0_dp  ! dummy
phi_site(4)    =    0.0_dp  ! dummy
x_site(4)      = 3500.0_dp *1.0e+03_dp    ! Point P4,
y_site(4)      = 2000.0_dp *1.0e+03_dp    ! conversion km -> m

ch_site(5)     = 'P5'
lambda_site(5) =    0.0_dp  ! dummy
phi_site(5)    =    0.0_dp  ! dummy
x_site(5)      = 3200.0_dp *1.0e+03_dp    ! Point P5,
y_site(5)      = 2000.0_dp *1.0e+03_dp    ! conversion km -> m

ch_site(6)     = 'P6'
lambda_site(6) =    0.0_dp  ! dummy
phi_site(6)    =    0.0_dp  ! dummy
x_site(6)      = 2900.0_dp *1.0e+03_dp    ! Point P6,
y_site(6)      = 2000.0_dp *1.0e+03_dp    ! conversion km -> m

ch_site(7)     = 'P7'
lambda_site(7) =    0.0_dp  ! dummy
phi_site(7)    =    0.0_dp  ! dummy
x_site(7)      = 2600.0_dp *1.0e+03_dp    ! Point P7,
y_site(7)      = 2000.0_dp *1.0e+03_dp    ! conversion km -> m

filename_with_path = trim(OUT_PATH)//'/'//trim(run_name)//'.site'

open(14, iostat=ios, file=trim(filename_with_path), status='new')

if (ios /= 0) then
   errormsg = ' >>> sico_init: Error when opening the site file!'
   call error(errormsg)
end if

if (forcing_flag == 1) then

   write(14,1106)
   write(14,1107)

   1106 format('         t(a)      D_Ts(C) z_sl_mean(m)',/, &
               '                   H_P1(m)      H_P2(m)      H_P3(m)      H_P4(m)', &
               '      H_P5(m)      H_P6(m)      H_P7(m)',/, &
               '                 v_P1(m/a)    v_P2(m/a)    v_P3(m/a)    v_P4(m/a)', &
               '    v_P5(m/a)    v_P6(m/a)    v_P7(m/a)',/, &
               '                   T_P1(C)      T_P2(C)      T_P3(C)      T_P4(C)', &
               '      T_P5(C)      T_P6(C)      T_P7(C)')
   1107 format('-----------------------------------------------------------------', &
               '---------------------------------------')

end if

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
!! For an initial topography with a very thin ice layer everywhere on the
!! land area.
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

character(len=256) :: filename_with_path

real(dp), dimension(0:JMAX,0:IMAX) :: field2d_aux

real(dp), parameter :: H_init = 1.0e-03_dp   ! initial ice thickness

!-------- Set topography --------

zl0 = 0.0_dp

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

errormsg = ' >>> topography1: GRID==2 not allowed for this application!'
call error(errormsg)

#endif

do i=0, IMAX
do j=0, JMAX

   zs(j,i) = zl0(j,i)
   zb(j,i) = zl0(j,i)
   zl(j,i) = zl0(j,i)

   if (mask(j,i) <= 1) then
      mask(j,i) = 0
      zs(j,i) = zs(j,i) + H_init
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

if (mask_region(0,0) == -1) mask_region = 0   ! regions undefined

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

!-------- Set topography --------

zl0 = 0.0_dp

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

errormsg = ' >>> topography2: GRID==2 not allowed for this application!'
call error(errormsg)

#endif

do i=0, IMAX
do j=0, JMAX

   zs(j,i) = zl0(j,i)
   zb(j,i) = zl0(j,i)
   zl(j,i) = zl0(j,i)

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

integer(i4b) :: i, j

character(len=256) :: filename_with_path

real(dp), dimension(0:JMAX,0:IMAX) :: field2d_aux

!-------- Read data from time-slice file of previous simulation --------

call read_tms_nc(anfdatname)

!-------- Set topography of the relaxed bedrock --------

zl0 = 0.0_dp

!-------- Further stuff --------

#if (GRID==0 || GRID==1)

dxi  = DX *1000.0_dp   ! km -> m
deta = DX *1000.0_dp   ! km -> m

#elif (GRID==2)

errormsg = ' >>> topography3: GRID==2 not allowed for this application!'
call error(errormsg)

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

if (mask_region(0,0) == -1) mask_region = 0   ! regions undefined

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
