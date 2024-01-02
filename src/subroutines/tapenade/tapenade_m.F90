!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  t a p e n a d e _ m
!
!> @file
!!
!! Catch-all module for tapenade-related subroutines.
!!
!! @section Copyright
!!
!! Copyright 2017-2024 Shreyas Sunil Gaikwad, Laurent Hascoet,
!!                     Sri Hari Krishna Narayanan, Liz Curry-Logan,
!!                     Patrick Heimbach, Ralf Greve
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
!> Catch-all module for tapenade-related subroutines.
!<------------------------------------------------------------------------------
module tapenade_m

  implicit none

  private
#ifdef ALLOW_TAPENADE
  public :: adjoint_master
#endif 

#if   (defined (ALLOW_GRDCHK))
  private :: deldirs
  public :: grdchk_main
#endif

contains

!-------------------------------------------------------------------------------
!> Adjoint master is the main tool by which sicopolis.F90 invokes the
!! adjoint/tlm code. Its job is to figure out what mode of the adjoint code is
!! being invoked and run the appropriate subroutine. 
!<------------------------------------------------------------------------------
#ifdef ALLOW_TAPENADE
  subroutine adjoint_master

use sico_variables_m_diff
#if (defined(GRL) && DISC>0)
  use discharge_workers_m_diff
#endif
  use ice_material_properties_m_diff
  use enth_temp_omega_m_diff
  use sico_init_m_diff
  use globals_diff
  USE CTRL_M_DIFF
  USE SICO_TYPES_M
  USE SICO_VARS_M
  USE SICO_MAIN_LOOP_M_DIFF
  USE SICO_END_M_DIFF

  implicit none
  integer(i4b)                               :: ndat2d, ndat3d
  integer(i4b)                               :: n_output
  real(dp)                                   :: delta_ts, glac_index
  real(dp)                                   :: mean_accum
  real(dp)                                   :: dtime, dtime_temp, &
                                                dtime_wss, dtime_out, dtime_ser
  real(dp)                                   :: time, time_init, time_end
  real(dp), dimension(100)                   :: time_output
  real(dp)                                   :: dxi, deta, dzeta_c, &
                                                dzeta_t, dzeta_r
  real(dp)                                   :: z_mar
  integer(i4b), parameter                    :: points = 5
  integer(i4b), dimension(points)            :: ipoints, jpoints
  integer(i4b)                               :: i, j, p
   !-------- Test points along spines of the ice sheets
   do p = 1, points
#if (defined(GRL))
      ipoints(p) = int(real(IMAX/2))
      jpoints(p) = int(real(JMAX/5)) + (p-1) * points
#elif (defined(ANT))
      ipoints(p) = int(real(IMAX/3)) + int(real((.85-.33)*IMAX/points)) * (p - 1)
      jpoints(p) = int(real(JMAX/2))
#endif
   end do

!@ python_automated_tlm IO begin @

	   open(9999, file='ForwardVals_q_geo_v5_grl20_ss25ka_limited.dat',&
	       form="FORMATTED", status="REPLACE")
	


   !-------- Loop over points
   do p = 1, points !@ python_automated_tlm limited_or_block_or_full @
     i = ipoints(p)
     j = jpoints(p)

  CALL SICO_INIT_D(delta_ts, glac_index, mean_accum, dtime, dtime_temp, &
&            dtime_wss, dtime_out, dtime_ser, time, time_init, time_end&
&            , time_output, dxi, deta, dzeta_c, dzeta_t, dzeta_r, z_mar&
&            , ndat2d, ndat3d, n_output)

!@ python_automated_tlm dep_vard @

		            q_geod = 0.0
		            q_geod(j,i) = 1.0
		
!-------- Main loop --------
  CALL SICO_MAIN_LOOP_D(delta_ts, glac_index, mean_accum, dtime, &
&                 dtime_temp, dtime_wss, dtime_out, dtime_ser, time, &
&                 time_init, time_end, time_output, dxi, deta, dzeta_c, &
&                 dzeta_t, dzeta_r, z_mar, ndat2d, ndat3d, n_output)
  CALL COST_FINAL_D()
     
  CALL SICO_END()

!@ python_automated_tlm IO write @
          write(9999, fmt='(f40.20)') fcd

   end do ! (close loop over points)


!@ python_automated_tlm IO end @
   close(unit=9999)

  end subroutine adjoint_master
#endif

#ifdef ALLOW_GRDCHK
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   Subroutine :  g r d c h k _ m a i n
!   Purpose    :  Gradient check top level routine
!                 Compares the gradients calculated by the adjoint model
!                 to the gradients calculated by finite difference
!                 approximations
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine grdchk_main
   
   use sico_types_m
   use sico_variables_m
   use sico_vars_m
   
   use sico_init_m
   use sico_main_loop_m
   use sico_end_m
   
   use ctrl_m
   
   implicit none
   
   integer(i4b)       :: ndat2d, ndat3d
   integer(i4b)       :: n_output
   real(dp)           :: delta_ts, glac_index
   real(dp)           :: mean_accum
   real(dp)           :: dtime, dtime_temp, dtime_wss, &
                                      dtime_out, dtime_ser
   real(dp)           :: time, time_init, time_end, time_output(100)
   real(dp)           :: dxi, deta, dzeta_c, dzeta_t, dzeta_r
   real(dp)           :: z_mar
   
   !-------- Variable declarations needed for this routine specifically
   real(dp)                          :: orig_val, perturb_val = 1.e-3
   real(dp),     dimension(3)        :: fc_collected
   real(dp),     dimension(3)        :: direction
   real(dp)                          :: gfd0,gfd, perturbation
   integer(i4b), parameter           :: points = 5
   integer(i4b), dimension(points)   :: ipoints, jpoints
   integer(i4b)                      :: i, j, p, d
   character(len=100)                :: fname
   
   !-------- This array holds the direction of perturbations to follow:
   direction(1) = 0
   direction(2) = 1
   direction(3) = -1

#if (!defined(GRL) && !defined(ANT))
   print *, ">>> Adjoint only available for GRL and ANT right now; kill code." 
#endif

   !-------- Test points along spines of the ice sheets
   do p = 1, points
#if (defined(GRL))
      ipoints(p) = int(real(IMAX/2))
      jpoints(p) = int(real(JMAX/5)) + (p-1) * points
#elif (defined(ANT))
      ipoints(p) = int(real(IMAX/3)) + int(real((.85-.33)*IMAX/points)) * (p - 1) 
      jpoints(p) = int(real(JMAX/2)) 
#endif
   end do

   !-------- Initialize output files 
   open(99, file='GradientVals_'//trim(HEADER)//'.dat',&
       form="FORMATTED", status="REPLACE")
   open(98, file='CostVals_'//trim(HEADER)//'.dat',&
       form="FORMATTED", status="REPLACE")

!@ python_automated_grdchk IO begin @

   
   !-------- Loop over points
   do p = 1, points !@ python_automated_grdchk limited_or_full @
     i = ipoints(p)
     j = jpoints(p)

          !-------- Loop over perturbation direction (0, +, -)
          do d = 1, 3 

          !-------- Let yourself know where you are:
          print *, ' point (p, i, j), direction (d) [ ', p , ', ', i, ', ', j, ', ', d, ' ] '

          !-------- One complete forward run 
            call deldirs
        
            call sico_init(delta_ts, glac_index, &
                 mean_accum, &
                 dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser, &
                 time, time_init, time_end, time_output, &
                 dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
                 z_mar, &
                 ndat2d, ndat3d, n_output)

            perturbation = 1 + direction(d) * perturb_val 


          !-------- Controls to be perturbed (add your own here and below in
          !         subroutine print_output()
          !         store original value that will be perturbed
          !         and then perturb it (first in +dir then -dir) 

                
            !@ python_automated_grdchk @

 
            ! -- H_c
            !orig_val = H(j,i)
            !H(j,i) = orig_val * perturbation 

            ! -- mean annual temp 
            !orig_val = temp_ma_present(j,i)
            !temp_ma_present(j,i) = orig_val * perturbation

            ! -- mean annual temp 
            !orig_val = q_geo(j,i)
            !q_geo(j,i) = orig_val * perturbation 

            ! -- mean annual temp 
            !orig_val = vx_c(24,j,i)
            !vx_c(24,j,i) = orig_val * perturbation

            ! -- mean annual temp 
            !orig_val = c_slide(j,i)
            !c_slide(j,i) = orig_val * perturbation

            ! -- mean annual temp 
            !orig_val = c_drag(j,i)
            !c_drag(j,i) = orig_val * perturbation

            ! -- precip_present
            !orig_val = precip_present(j,i,1)
            !precip_present(j,i,1) = orig_val * perturbation

            ! -- experimental field acc_fact
            !orig_val = acc_fact(j,i)
            !acc_fact(j,i) = orig_val * perturbation

            ! -- sanity check
            write(6,fmt='(a,f40.20)') "orig_val = ", orig_val
            write(6,fmt='(a,f40.20)') "pert_val = ", orig_val*perturbation

            call sico_main_loop(delta_ts, glac_index, &
                 mean_accum, &
                 dtime, dtime_temp, dtime_wss, dtime_out, dtime_ser, &
                 time, time_init, time_end, time_output, &
                 dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
                 z_mar, &
                 ndat2d, ndat3d, n_output)
          
            call cost_final()
            call sico_end
       
            ! store cost
            fc_collected(d) = fc
        
            ! --------- calculate simple 2-sided finite difference due to
            !           perturbation: fc(+) - fc(-) / 2*espsilon
            if (orig_val .ne. 0) then
                gfd = (fc_collected(2) - fc_collected(3))/(2.d0 * perturb_val * orig_val)
            else
                gfd = 0.0
            end if          
          end do ! (close perturb loop)

          ! -- sanity check
          write(6, fmt='(a,f40.20)')   "Finite difference is = ", gfd
        
          ! --------- write these values to output file
          write(99, fmt='(f40.20)') gfd
          write(98, fmt='(f40.20)') fc_collected(1)
          write(98, fmt='(f40.20)') fc_collected(2)
          write(98, fmt='(f40.20)') fc_collected(3)
          write(98, fmt='(a)') '----------------------------------'
          
          !@ python_automated_grdchk IO write @
    
   end do ! (close loop over points)
  
   close(unit=99)
   close(unit=98)
   !@ python_automated_grdchk IO end @   
   end subroutine grdchk_main

!!-------------------------------------------------------------------------------
!!> Checks to see if output dir exists. If so, deletes it.
!!<------------------------------------------------------------------------------
  subroutine deldirs

  implicit none

  character(len=256) :: shell_command

  !-------- deleting directories

    shell_command = 'if [ -d'
    shell_command = trim(shell_command)//' '//OUT_PATH
    shell_command = trim(shell_command)//' '//'] ; then rm -rf'
    shell_command = trim(shell_command)//' '//OUT_PATH
    shell_command = trim(shell_command)//' '//'; fi'

    call system(trim(shell_command))

  end subroutine deldirs

#endif /* Only ALLOW_GRDCHK */

end module tapenade_m 
!
