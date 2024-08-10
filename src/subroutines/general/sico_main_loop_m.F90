!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ m a i n _ l o o p _ m
!
!! Main loop of SICOPOLIS.
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
!> Main loop of SICOPOLIS.
!-------------------------------------------------------------------------------
module sico_main_loop_m
  
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
!> Main routine of sico_main_loop_m: Main loop of SICOPOLIS.
!-------------------------------------------------------------------------------
!@ begin tapenade_extract @
  subroutine sico_main_loop(delta_ts, glac_index, &
                      mean_accum, &
                      dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser, &
                      time, time_init, time_end, time_output, &
                      dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
                      z_mar, &
                      ndat2d, ndat3d, n_output)
!@ end tapenade_extract @
 
  use boundary_m
  
#if (CALCMOD==0 || CALCMOD==1 || CALCMOD==-1)
  use calc_temp_m
#elif (CALCMOD==2 || CALCMOD==3)
  use calc_temp_enth_m
#endif
  
  use calc_enhance_m
  use flag_update_gf_gl_cf_m
  use calc_vxy_m
  use calc_vz_m
  use calc_dxyz_m
  use calc_gia_m
  use calc_thk_m
  use calc_temp_melt_bas_m
  use calc_bas_melt_m
  use calc_thk_water_bas_m
  use calc_pressure_water_bas_m
  use output_m
 
  implicit none

#if !defined(ALLOW_TAPENADE) /* Normal */
  integer(i4b),       intent(in)    :: n_output
  real(dp),           intent(in)    :: mean_accum
  real(dp),           intent(in)    :: dtime, dtime_temp, dtime_wss, &
                                       dtime_out, dtime_ser
  real(dp),           intent(in)    :: time_init, time_end, time_output(100)
  real(dp),           intent(in)    :: dxi, deta, dzeta_c, dzeta_t, dzeta_r
#else
  integer(i4b),       intent(inout)    :: n_output
  real(dp),           intent(inout)    :: mean_accum
  real(dp),           intent(inout)    :: dtime, dtime_temp, dtime_wss, &
                                          dtime_out, dtime_ser
  real(dp),           intent(inout)    :: time_init, time_end, time_output(100)
  real(dp),           intent(inout)    :: dxi, deta, dzeta_c, dzeta_t, dzeta_r
#endif  
  integer(i4b),       intent(inout) :: ndat2d, ndat3d
  real(dp),           intent(inout) :: delta_ts, glac_index
  real(dp),           intent(inout) :: time
  real(dp),           intent(inout) :: z_mar
  
  integer(i4b) :: i, j, kc, kt, kr, n
  integer(i4b) :: itercount, itercount_max
  integer(i4b) :: iter_temp, iter_wss, iter_ser, iter_out, iter_output(100)
  real(dp)     :: dtime_temp_inv
  logical      :: flag_3d_output, flag_output1, flag_output2
  
  !-------- Begin of main loop (time integration) --------
  
  write(unit=6, fmt='(/a/)') ' -------- sico_main_loop --------'
  
  itercount_max = nint((time_end-time_init)/dtime)
  
  iter_temp = nint(dtime_temp/dtime)
#if (REBOUND==2)
  iter_wss = nint(dtime_wss/dtime)
#endif
  iter_ser  = nint(dtime_ser/dtime)
#if (OUTPUT==1 || OUTPUT==3)
  iter_out  = nint(dtime_out/dtime)
#endif
#if (OUTPUT==2 || OUTPUT==3)
  do n=1, n_output
     iter_output(n) = nint((time_output(n)-time_init)/dtime)
  end do
#endif

  !$NO AD BINOMIAL-CKP itercount_max+1 20 1
  main_loop : do itercount=1, itercount_max
  
  write(unit=6, fmt='(2x,i0)') itercount
  
  !-------- Update of time --------
  
  time      = time_init + real(itercount,dp)*dtime
  
  !-------- Save old mask --------
  
  mask_old = mask
  
  !-------- Boundary conditions --------
  
  call boundary(time, dtime, dxi, deta, delta_ts, glac_index, z_mar)
  
  !-------- Temperature, water content, age, flow enhancement factor --------
  
  if ( mod(itercount, iter_temp) == 0 ) then
     flag_calc_temp = .true.
  else
     flag_calc_temp = .false.
  end if

  if ( flag_calc_temp ) then

       write(unit=6, fmt='(10x,a)') 'Computation of T'
  
  !  ------ Temperature, water content, age
  
#if (CALCMOD==1)
     call calc_temp_poly(dxi, deta, dzeta_c, dzeta_t, dzeta_r, dtime_temp)
                    ! polythermal scheme (POLY)
  
#elif (CALCMOD==0)
     call calc_temp_cold(dxi, deta, dzeta_c, dzeta_t, dzeta_r, dtime_temp)
                    ! cold-ice scheme (COLD)
  
#elif (CALCMOD==2 || CALCMOD==3)
     call calc_temp_enth(dxi, deta, dzeta_c, dzeta_t, dzeta_r, dtime_temp)
                    ! enthalpy scheme (ENTC or ENTM)
  
#elif (CALCMOD==-1)
     call calc_temp_const()   ! isothermal scheme (ISOT)
  
#else
     errormsg = ' >>> sico_main_loop: ' &
              //'Parameter CALCMOD must be between -1 and 3!'
     call error(errormsg)
#endif
  
  !  ------ Time derivative of H_t (further time derivatives are
  !         computed in subroutine calc_thk_xxx)
  
     dtime_temp_inv = 1.0_dp/dtime_temp
  
     dH_t_dtau      = (H_t_new - H_t)*dtime_temp_inv
  
  !  ------ New values -> old values
  
     n_cts   = n_cts_new
     kc_cts  = kc_cts_new
     zm      = zm_new
     H_c     = H_c_new
     H_t     = H_t_new
     temp_c  = temp_c_new
     age_c   = age_c_new
     omega_t = omega_t_new
     age_t   = age_t_new
     temp_r  = temp_r_new
  
#if (CALCMOD==2 || CALCMOD==3)
     enth_c  = enth_c_new
     enth_t  = enth_t_new
     omega_c = omega_c_new
#endif
  
  !  ------ Flow enhancement factor
  
#if (ENHMOD==1)
     call calc_enhance_1()
#elif (ENHMOD==2)
     call calc_enhance_2()
#elif (ENHMOD==3)
     call calc_enhance_3(time)
#elif (ENHMOD==4)
     call calc_enhance_4()
#elif (ENHMOD==5)
     call calc_enhance_5()
#else
     errormsg = ' >>> sico_main_loop: ' &
              //'Parameter ENHMOD must be between 1 and 5!'
     call error(errormsg)
#endif
  
  end if
  !     End of computation of temperature, water content, age and
  !     enhancement factor (only if flag_calc_temp == .true.)
  
  !-------- Velocity --------

  call flag_update_gf_gl_cf()
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
  errormsg = ' >>> sico_main_loop: ' &
           //'Parameter DYNAMICS must be either 0, 1 or 2!'
  call error(errormsg)
#endif
  
  call calc_dxyz(dxi, deta, dzeta_c, dzeta_t)
  
  !-------- Glacial isostatic adjustment and ice topography --------
  
  call calc_gia(time, dtime, dxi, deta, itercount, iter_wss)
  
  call calc_thk_init()
  
#if ((MARGIN==3 || DYNAMICS==2) && (CALCTHK==1 || CALCTHK==2))
    errormsg = ' >>> sico_main_loop:' &
             //           end_of_line &
             //'          Non-SIA dynamics combined with' &
             //           end_of_line &
             //'          SIA ice thickness solver!'
    call error(errormsg)
#endif
  
#if (CALCTHK==1)
    call calc_thk_sia_expl(time, dtime, dxi, deta, z_mar)
#elif (CALCTHK==2)
    call calc_thk_sia_impl(time, dtime, dxi, deta, z_mar, mean_accum)
#elif (CALCTHK==4)
    call calc_thk_expl(time, dtime, dxi, deta, z_mar)
#else
    errormsg = ' >>> sico_main_loop: ' &
             //'Parameter CALCTHK must be 1, 2 or 4!'
    call error(errormsg)
#endif
  
#if (MARGIN==3)       /* coupled SIA/SSA or SIA/SStA/SSA dynamics */
    call calc_thk_mask_update(time, dtime, dxi, deta, z_mar, 3)
#elif (DYNAMICS==2)   /* hybrid SIA/SStA dynamics */
    call calc_thk_mask_update(time, dtime, dxi, deta, z_mar, 2)
#else                 /* SIA-only dynamics */
#if (CALCTHK==1 || CALCTHK==2)
    call calc_thk_mask_update(time, dtime, dxi, deta, z_mar, 1)
#elif (CALCTHK==4)
    call calc_thk_mask_update(time, dtime, dxi, deta, z_mar, 2)
#endif
#endif

  call account_mb_source(dtime, z_mar)

  call flag_update_gf_gl_cf()

  !  ------ New values -> old values
  
  zs  = zs_new
  zm  = zm_new
  zb  = zb_new
  zl  = zl_new
  H   = H_new
  H_c = H_c_new
  H_t = H_t_new
  
  !-------- Melting temperature --------

  call calc_temp_melt()

  !-------- Basal temperature --------

  call calc_temp_bas()

  !-------- Basal melting rate --------

  call calc_qbm(time, dzeta_c, dzeta_r)

  !-------- Effective thickness of subglacial water  --------

  call calc_thk_water_bas()

  !-------- Basal water pressure --------

  call calc_pressure_water_bas()

  !-------- Data output --------

#if !(defined(ALLOW_GRDCHK) || defined(ALLOW_TAPENADE))

  !  ------ Time-slice data
  
#if (OUTPUT==1)

  flag_output1 = .false.

  if ( mod(itercount, iter_out) == 0 ) then
  
#if (ERGDAT==0)
     flag_3d_output = .false.
#elif (ERGDAT==1)
     flag_3d_output = .true.
#endif
  
     call output1(time, delta_ts, glac_index, &
                  flag_3d_output, ndat2d, ndat3d)

     flag_output1 = .true.

  end if
  
#if (OUTPUT_FLUX_VARS==2)   /* averaging of flux variables */

  if (.not.flag_output1) then
  
#if (ERGDAT==0)
     flag_3d_output = .false.
#elif (ERGDAT==1)
     flag_3d_output = .true.
#endif

     call output1(time, delta_ts, glac_index, &
                  flag_3d_output, ndat2d, ndat3d, &
                  opt_flag_compute_flux_vars_only=.true.)

  end if

#endif

#elif (OUTPUT==2)
  
  flag_output1 = .false.

  do n=1, n_output
  
     if (itercount == iter_output(n)) then
  
#if (ERGDAT==0)
        flag_3d_output = .false.
#elif (ERGDAT==1)
        flag_3d_output = .true.
#endif
  
        call output1(time, delta_ts, glac_index, &
                     flag_3d_output, ndat2d, ndat3d)

        flag_output1 = .true.

     end if
  
  end do

#if (OUTPUT_FLUX_VARS==2)   /* averaging of flux variables */

  if (.not.flag_output1) then

#if (ERGDAT==0)
     flag_3d_output = .false.
#elif (ERGDAT==1)
     flag_3d_output = .true.
#endif

     call output1(time, delta_ts, glac_index, &
                  flag_3d_output, ndat2d, ndat3d, &
                  opt_flag_compute_flux_vars_only=.true.)

  end if

#endif

#elif (OUTPUT==3)

  flag_output1 = .false.
  
  if ( mod(itercount, iter_out) == 0 ) then
  
     flag_3d_output = .false.
  
     call output1(time, delta_ts, glac_index, &
                  flag_3d_output, ndat2d, ndat3d)

     flag_output1 = .true.
  
  end if

#if (OUTPUT_FLUX_VARS==2)   /* averaging of flux variables */

  if (.not.flag_output1) then

     flag_3d_output = .false.

     call output1(time, delta_ts, glac_index, &
                  flag_3d_output, ndat2d, ndat3d, &
                  opt_flag_compute_flux_vars_only=.true.)

  end if

#endif

  do n=1, n_output
  
     if (itercount == iter_output(n)) then
  
        flag_3d_output = .true.
  
        call output1(time, delta_ts, glac_index, &
                     flag_3d_output, ndat2d, ndat3d)
  
     end if
  
  end do
  
#endif

  !  ------ Time-series data

  flag_output2 = .false.

  if ( mod(itercount, iter_ser) == 0 ) then

     call output2(time, dxi, deta, delta_ts, glac_index)
     flag_output2 = .true.

     call output4(time, dxi, deta, delta_ts, glac_index)

#if (defined(ASF) && WRITE_SER_FILE_STAKES==1) /* Austfonna */
     call output5(time, dxi, deta, delta_ts, glac_index)
#endif

  end if

#if (OUTPUT_FLUX_VARS==2)   /* averaging of flux variables */

  if (.not.flag_output2) &
     call output2(time, dxi, deta, delta_ts, glac_index, &
                  opt_flag_compute_flux_vars_only=.true.)

#endif

#endif   /* !ALLOW_GRDCHK & !ALLOW_TAPENADE */

  end do main_loop   ! End of main loop (time integration)

  end subroutine sico_main_loop

!-------------------------------------------------------------------------------
  
end module sico_main_loop_m
!
