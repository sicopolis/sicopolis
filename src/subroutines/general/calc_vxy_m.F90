!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ v x y _ m
!
!! Computation of the horizontal velocity vx, vy.
!!
!!##### Authors
!!
!! Ralf Greve, Tatsuru Sato, Thomas Goelles, Jorge Bernales, Felix Grandadam
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
!> Computation of the horizontal velocity vx, vy.
!-------------------------------------------------------------------------------
module calc_vxy_m

  use sico_types_m
  use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

  use error_m

  implicit none

  private
  public :: calc_vxy_b_init, calc_dzs_dxy_aux, &
            calc_vxy_b_sia, calc_vxy_sia, calc_vxy_static, calc_vxy_ssa

contains

!-------------------------------------------------------------------------------
!> Initializations for the basal horizontal velocity vx_b, vy_b.
!-------------------------------------------------------------------------------
subroutine calc_vxy_b_init()

implicit none

integer(i4b) :: i, j, m, n
integer(i4b) :: n_slide_regions
integer(i4b) :: i_f, j_f, n_filter

#if (!defined(N_SLIDE_REGIONS) || N_SLIDE_REGIONS<=1)
integer(i4b) :: p_weert_aux(1)
integer(i4b) :: q_weert_aux(1)
real(dp) :: c_slide_aux(1)
real(dp) :: gamma_slide_aux(1)
real(dp) :: gamma_slide_inv_aux(1)
#else
integer(i4b) :: p_weert_aux(N_SLIDE_REGIONS)
integer(i4b) :: q_weert_aux(N_SLIDE_REGIONS)
real(dp) :: c_slide_aux(N_SLIDE_REGIONS)
real(dp) :: gamma_slide_aux(N_SLIDE_REGIONS)
real(dp) :: gamma_slide_inv_aux(N_SLIDE_REGIONS)
#endif

real(dp) :: dx
real(dp) :: filter_width, sigma_filter
real(dp) :: dist, weigh, sum_weigh
real(dp), dimension(0:JMAX,0:IMAX) :: c_slide_smoothed

!-------- Sliding-law coefficients --------

#if (!defined(N_SLIDE_REGIONS) || N_SLIDE_REGIONS<=1)
n_slide_regions = 1
#else
n_slide_regions = N_SLIDE_REGIONS
#endif

#if (SLIDE_LAW==0)

p_weert_aux = 1
q_weert_aux = 0
c_slide_aux = 0.0_dp   ! no-slip
gamma_slide_aux = 1.0_dp

#elif (SLIDE_LAW==1)

p_weert_aux = P_WEERT
q_weert_aux = Q_WEERT
c_slide_aux = C_SLIDE
gamma_slide_aux = GAMMA_SLIDE

#else

errormsg = ' >>> calc_vxy_b_init: SLIDE_LAW must be 0 or 1!' &
         //         end_of_line &
         //'        Change obsolete SLIDE_LAW = 2 or 3' &
         //         end_of_line &
         //'        to SLIDE_LAW = 1, BASAL_WATER_PRESSURE = 2.'
call error(errormsg)

#endif

do n=1, n_slide_regions
   gamma_slide_inv_aux(n) = 1.0_dp/max(gamma_slide_aux(n), eps)
end do

do i=0, IMAX
do j=0, JMAX
   if ( (n_slide_region(j,i) >= 1) &
        .and. &
        (n_slide_region(j,i) <= n_slide_regions) ) then
      p_weert(j,i)         = p_weert_aux(n_slide_region(j,i))
      q_weert(j,i)         = q_weert_aux(n_slide_region(j,i))
      c_slide_init(j,i)    = c_slide_aux(n_slide_region(j,i))*sec2year
      gamma_slide_inv(j,i) = gamma_slide_inv_aux(n_slide_region(j,i))
      sub_melt_flag(j,i)   = (gamma_slide_aux(n_slide_region(j,i)) >= eps)
   else
      errormsg = ' >>> calc_vxy_b_init: ' &
                    //'Region number out of allowed range!'
      call error(errormsg)
   end if
end do
end do

do i=0, IMAX
do j=0, JMAX
   p_weert_inv(j,i) = 1.0_dp/max(real(p_weert(j,i),dp), eps)
end do
end do

!  ------ Smoothing the coefficient c_slide_init by a Gaussian filter
!         (only meaningful if the exponents p_weert and q_weert are constant!)

#if (defined(C_SLIDE_FILTER_WIDTH))

filter_width = real(C_SLIDE_FILTER_WIDTH, dp)
                   ! filter width (half span of filtered area), in km

if (filter_width > eps_sp_dp) then

#if (GRID==0 || GRID==1)
   dx = real(DX, dp)   ! horizontal grid spacing (in km)
#else
   dx = 0.0_dp   ! dummy value
   errormsg = ' >>> calc_vxy_b_init: ' &
           // 'Smoothing of c_slide_init only implemented for GRID 0 or 1!'
   call error(errormsg)
#endif

   sigma_filter = filter_width/dx   ! half span of filtered area,
                                    ! in grid points

   n_filter = ceiling(2.0_dp*sigma_filter)
   n_filter = max(n_filter, 5)

   c_slide_smoothed = 0.0_dp

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

         c_slide_smoothed(j,i) = c_slide_smoothed(j,i) &
                                    + weigh*c_slide_init(j_f,i_f)

      end do
      end do

      c_slide_smoothed(j,i) = c_slide_smoothed(j,i)/sum_weigh

   end do
   end do

   c_slide_init = c_slide_smoothed

end if

#endif

end subroutine calc_vxy_b_init

!-------------------------------------------------------------------------------
!> Computation of the auxiliary surface gradients dzs_dx_aux, dzs_dy_aux
!! (optional one-sided gradients at the grounding line).
!-------------------------------------------------------------------------------
subroutine calc_dzs_dxy_aux(dxi, deta)

implicit none

real(dp), intent(in) :: dxi, deta

integer(i4b) :: i, j
real(dp)     :: inv_dx, inv_dy
real(dp)     :: rhosw_rho_ratio
real(dp)     :: H_mid, zl_mid, zs_mid, z_sl_mid

inv_dx          = 1.0_dp/dxi
inv_dy          = 1.0_dp/deta
rhosw_rho_ratio = RHO_SW/RHO

dzs_dx_aux = dzs_dxi
dzs_dy_aux = dzs_deta

#if (MARGIN==3)

#if (!defined(GL_SURF_GRAD) || GL_SURF_GRAD==1)

!%% continue

#elif (GL_SURF_GRAD==2)

do i=0, IMAX-1
do j=1, JMAX-1
   ! inner point on the staggered grid in x-direction

   if ( (flag_grounding_line_1(j,i).and.flag_grounding_line_2(j,i+1)) &
        .or. &
        (flag_grounding_line_2(j,i).and.flag_grounding_line_1(j,i+1)) &
      ) then
        ! one neighbour is floating ice and the other is grounded ice
        ! (grounding line)

      H_mid    = 0.5_dp*(H(j,i)   +H(j,i+1))
      zl_mid   = 0.5_dp*(zl(j,i)  +zl(j,i+1))
      zs_mid   = 0.5_dp*(zs(j,i)  +zs(j,i+1))
      z_sl_mid = 0.5_dp*(z_sl(j,i)+z_sl(j,i+1))

      if (H_mid < (z_sl_mid-zl_mid)*rhosw_rho_ratio) then
         ! floating condition is satisfied

         if ( &
              (flag_grounding_line_1(j,i).and.flag_grounding_line_2(j,i+1)) &
              .and. &
              (i+2 <= IMAX) &
            ) then

            if ((mask(j,i+2) == 3).or.(mask(j,i+2) == 2)) &
               dzs_dx_aux(j,i) = (0.5_dp*(zs(j,i+1)+zs(j,i+2))-zs_mid)*inv_dx
                                 ! one-sided gradient into floating ice

         else if ( &
              (flag_grounding_line_2(j,i).and.flag_grounding_line_1(j,i+1)) &
              .and. &
              (i-1 >= 0) &
            ) then

            if ((mask(j,i-1) == 3).or.(mask(j,i-1) == 2)) &
               dzs_dx_aux(j,i) = (zs_mid-0.5_dp*(zs(j,i)+zs(j,i-1)))*inv_dx
                                 ! one-sided gradient into floating ice

         end if

      else   ! H_mid >= (z_sl_mid-zl_mid)*rhosw_rho_ratio,
             ! floating condition is not satisfied

         if ( &
              (flag_grounding_line_1(j,i).and.flag_grounding_line_2(j,i+1)) &
              .and. &
              (i-1 >= 0) &
            ) then

            if ((mask(j,i-1) == 0).or.(mask(j,i-1) == 1)) &
               dzs_dx_aux(j,i) = (zs_mid-0.5_dp*(zs(j,i)+zs(j,i-1)))*inv_dx
                                 ! one-sided gradient into grounded ice

         else if ( &
              (flag_grounding_line_2(j,i).and.flag_grounding_line_1(j,i+1)) &
              .and. &
              (i+2 <= IMAX) &
            ) then

            if ((mask(j,i+2) == 0).or.(mask(j,i+2) == 1)) &
               dzs_dx_aux(j,i) = (0.5_dp*(zs(j,i+1)+zs(j,i+2))-zs_mid)*inv_dx
                                 ! one-sided gradient into grounded ice

         end if

      end if


   end if

end do
end do

do i=1, IMAX-1
do j=0, JMAX-1
   ! inner point on the staggered grid in y-direction

   if ( (flag_grounding_line_1(j,i).and.flag_grounding_line_2(j+1,i)) &
        .or. &
        (flag_grounding_line_2(j,i).and.flag_grounding_line_1(j+1,i)) &
      ) then
        ! one neighbour is floating ice and the other is grounded ice
        ! (grounding line)

      H_mid    = 0.5_dp*(H(j,i)   +H(j+1,i))
      zl_mid   = 0.5_dp*(zl(j,i)  +zl(j+1,i))
      zs_mid   = 0.5_dp*(zs(j,i)  +zs(j+1,i))
      z_sl_mid = 0.5_dp*(z_sl(j,i)+z_sl(j+1,i))

      if (H_mid < (z_sl_mid-zl_mid)*rhosw_rho_ratio) then
         ! floating condition is satisfied

         if ( &
              (flag_grounding_line_1(j,i).and.flag_grounding_line_2(j+1,i)) &
              .and. &
              (j+2 <= JMAX) &
            ) then

            if ((mask(j+2,i) == 3).or.(mask(j+2,i) == 2)) &
               dzs_dy_aux(j,i) = (0.5_dp*(zs(j+1,i)+zs(j+2,i))-zs_mid)*inv_dy
                                 ! one-sided gradient into floating ice

         else if ( &
              (flag_grounding_line_2(j,i).and.flag_grounding_line_1(j+1,i)) &
              .and. &
              (j-1 >= 0) &
            ) then

            if ((mask(j-1,i) == 3).or.(mask(j-1,i) == 2)) &
               dzs_dy_aux(j,i) = (zs_mid-0.5_dp*(zs(j,i)+zs(j-1,i)))*inv_dy
                                 ! one-sided gradient into floating ice

         end if

      else   ! H_mid >= (z_sl_mid-zl_mid)*rhosw_rho_ratio,
             ! floating condition is not satisfied

         if ( &
              (flag_grounding_line_1(j,i).and.flag_grounding_line_2(j+1,i)) &
              .and. &
              (j-1 >= 0) &
            ) then

            if ((mask(j-1,i) == 0).or.(mask(j-1,i) == 1)) &
               dzs_dy_aux(j,i) = (zs_mid-0.5_dp*(zs(j,i)+zs(j-1,i)))*inv_dy
                                 ! one-sided gradient into grounded ice

         else if ( &
              (flag_grounding_line_2(j,i).and.flag_grounding_line_1(j+1,i)) &
              .and. &
              (j+2 <= JMAX) &
            ) then

            if ((mask(j+2,i) == 0).or.(mask(j+2,i) == 1)) &
               dzs_dy_aux(j,i) = (0.5_dp*(zs(j+1,i)+zs(j+2,i))-zs_mid)*inv_dy
                                 ! one-sided gradient into grounded ice

         end if

      end if

   end if

end do
end do

#else

errormsg = ' >>> calc_dzs_dxy_aux: ' &
              //'GL_SURF_GRAD must be either 1 or 2!'
call error(errormsg)

#endif

#endif

end subroutine calc_dzs_dxy_aux

!-------------------------------------------------------------------------------
!> Computation of the basal horizontal velocity vx_b, vy_b in the shallow ice
!! approximation.
!-------------------------------------------------------------------------------
subroutine calc_vxy_b_sia(time)

use calc_pressure_water_bas_m

implicit none

real(dp), intent(in) :: time

integer(i4b) :: i, j
real(dp), dimension(0:JMAX,0:IMAX) :: p_b_red_lim
real(dp) :: cvxy1, cvxy1a, cvxy1b, ctau1, ctau1a, ctau1b
real(dp) :: temp_diff
real(dp) :: c_Hw_slide, Hw0_slide, Hw0_slide_inv, ratio_Hw_slide
real(dp) :: vh_max, vh_max_inv
real(dp) :: time_in_years
real(dp) :: ramp_up_factor
real(dp) :: zs_grad_sq

time_in_years = time*sec2year

!  ------ Ramping up basal sliding

ramp_up_factor = 1.0_dp

c_slide = c_slide_init

#if (defined(TIME_RAMP_UP_SLIDE))

if ( (TIME_RAMP_UP_SLIDE > eps_dp) &
     .and. &
     ((time_in_years-(TIME_INIT0)) < (TIME_RAMP_UP_SLIDE)) ) then

   ramp_up_factor = (time_in_years-(TIME_INIT0))/(TIME_RAMP_UP_SLIDE)

   ramp_up_factor = max(min(ramp_up_factor, 1.0_dp), 0.0_dp)
                                  ! constrain to interval [0,1]

   ramp_up_factor = ramp_up_factor*ramp_up_factor*ramp_up_factor &
                                  *(10.0_dp + ramp_up_factor &
                                              *(-15.0_dp+6.0_dp*ramp_up_factor))
                                  ! make transition smooth (quintic function)

   c_slide = c_slide * ramp_up_factor

end if

#endif

!  ------ Coefficients for the contribution of the basal water layer

#if (BASAL_HYDROLOGY==1)

#if (!defined(HYDRO_SLIDE_SAT_FCT))
  errormsg = ' >>> calc_vxy_b_sia: ' &
                //'HYDRO_SLIDE_SAT_FCT must be defined!'
  call error(errormsg)
#endif

#if (defined(C_HW_SLIDE))
  c_Hw_slide = C_HW_SLIDE
#else
  errormsg = ' >>> calc_vxy_b_sia: C_HW_SLIDE must be defined!'
  call error(errormsg)
#endif

#if (defined(HW0_SLIDE))
  Hw0_slide = HW0_SLIDE
#else
  errormsg = ' >>> calc_vxy_b_sia: HW0_SLIDE must be defined!'
  call error(errormsg)
#endif

  Hw0_slide_inv = 1.0_dp/max(Hw0_slide, eps_dp)

#else   /* BASAL_HYDROLOGY==0 */

  c_Hw_slide    = 0.0_dp   ! dummy value
  Hw0_slide     = 1.0_dp   ! dummy value
  Hw0_slide_inv = 1.0_dp   ! dummy value

#endif

!-------- Computation of basal stresses --------

!  ------ Basal pressure p_b, basal water pressure p_b_w,
!         reduced pressure p_b_red

call calc_pressure_water_bas()   ! compute p_b_w

do i=0, IMAX
do j=0, JMAX

   if ((mask(j,i) == 0).or.(mask(j,i) == 3)) then
                     ! grounded or floating ice

      p_b(j,i) = max(RHO*G*H(j,i), 0.0_dp)

      if (mask(j,i) == 0) then
         p_b_red(j,i) = max(p_b(j,i)-p_b_w(j,i), 0.0_dp)
      else   ! (mask(j,i) == 3)
         p_b_red(j,i) = 0.0_dp
      end if

      p_b_red_lim(j,i) = max(p_b_red(j,i), RED_PRES_LIMIT_FACT*p_b(j,i))
                         ! in order to avoid very small values, which may lead
                         ! to huge sliding velocities in the SIA

   else   ! ice-free land or ocean

      p_b(j,i)         = 0.0_dp
      p_b_red(j,i)     = 0.0_dp
      p_b_red_lim(j,i) = 0.0_dp

   end if

end do
end do

!  ------ Absolute value of the absolute value (magnitude)
!         of the driving stress (tau_dr) and the basal shear stress (tau_b)

do i=0, IMAX
do j=0, JMAX

   if ((mask(j,i) == 0).or.(mask(j,i) == 3)) then
                     ! grounded or floating ice

      zs_grad_sq = dzs_dxi_g(j,i)**2+dzs_deta_g(j,i)**2

#if !defined(ALLOW_TAPENADE) /* Normal */

      tau_dr(j,i) = p_b(j,i)*sqrt(zs_grad_sq)

#else /* Tapenade: guarding against non-differentiable sqrt(0) */

      if (zs_grad_sq > 0) then
         tau_dr(j,i) = p_b(j,i)*sqrt(zs_grad_sq)
      else
         tau_dr(j,i) = 0.0_dp
      end if

#endif /* Normal vs. Tapenade */

      tau_b(j,i) = tau_dr(j,i)
                   ! for hybrid dynamics or floating ice,
                   ! this will be corrected later in calc_vxy_ssa

   else   ! ice-free land or ocean

      tau_dr(j,i) = 0.0_dp
      tau_b(j,i)  = 0.0_dp

   end if

end do
end do

!-------- Computation of d_help_b (defined on the grid points (i,j)) --------

do i=0, IMAX
do j=0, JMAX

   if ((mask(j,i) == 0).or.flag_grounding_line_2(j,i)) then
                     ! grounded ice, or floating ice at the grounding line

!  ------ Abbreviations

#if (SLIDE_LAW==0)

      cvxy1 = 0.0_dp   ! No-slip
      ctau1 = 1.0_dp/eps_dp

#elif (SLIDE_LAW==1)

      cvxy1 = c_slide(j,i) &
              * ( (tau_b(j,i)+eps_dp)**(p_weert(j,i)-1) &
                  /(p_b_red_lim(j,i)+eps_dp)**q_weert(j,i) ) &
              * p_b(j,i)
      ctau1 = 1.0_dp/(c_slide(j,i)+eps_dp)**p_weert_inv(j,i) &
              * (p_b_red(j,i)+eps_dp)**(q_weert(j,i)*p_weert_inv(j,i))
              ! Basal sliding at pressure melting

#else

      errormsg = ' >>> calc_vxy_b_sia: SLIDE_LAW must be 0 or 1!' &
               //         end_of_line &
               //'        Change obsolete SLIDE_LAW = 2 or 3' &
               //         end_of_line &
               //'        to SLIDE_LAW = 1, BASAL_WATER_PRESSURE = 2.'
      call error(errormsg)

#endif

!  ------ d_help_b, c_drag

      if (n_cts(j,i) == -1) then   ! cold ice base

         if (sub_melt_flag(j,i)) then
            temp_diff = max((temp_c_m(0,j,i)-temp_c(0,j,i)), 0.0_dp)
            cvxy1a    = exp(-gamma_slide_inv(j,i)*temp_diff)  ! sub-melt sliding
            ctau1a    = 1.0_dp/(cvxy1a+eps_dp)**p_weert_inv(j,i)
         else
            cvxy1a    = 0.0_dp   ! no sub-melt sliding
            ctau1a    = 1.0_dp/eps_dp**p_weert_inv(j,i)   ! dummy value
         end if

         d_help_b(j,i) = cvxy1*cvxy1a
         c_drag(j,i)   = ctau1*ctau1a

      else if (n_cts(j,i) == 0) then   ! temperate ice base

         d_help_b(j,i) = cvxy1   ! basal sliding
         c_drag(j,i)   = ctau1   ! (pressure-melting conditions)

      else   ! n_cts(j,i) == 1, temperate ice layer

         d_help_b(j,i) = cvxy1   ! basal sliding
         c_drag(j,i)   = ctau1   ! (pressure-melting conditions)

      end if

!    ---- Contribution of the basal water layer

#if (BASAL_HYDROLOGY==1)

#if (HYDRO_SLIDE_SAT_FCT==0)

      ratio_Hw_slide = max(H_w(j,i)*Hw0_slide_inv, 0.0_dp)
                        ! constrain to interval [0,infty)

      cvxy1b = 1.0_dp + c_Hw_slide*(1.0_dp-exp(-ratio_Hw_slide))
                        ! exponential saturation function
                        ! by Kleiner and Humbert (2014, J. Glaciol. 60)

#elif (HYDRO_SLIDE_SAT_FCT==1)

      ratio_Hw_slide = max(min(H_w(j,i)*Hw0_slide_inv, 1.0_dp), 0.0_dp)
                        ! constrain to interval [0,1]

      cvxy1b = 1.0_dp + c_Hw_slide*ratio_Hw_slide
                        ! linear saturation function

#elif (HYDRO_SLIDE_SAT_FCT==2)

      ratio_Hw_slide = max(min(H_w(j,i)*Hw0_slide_inv, 1.0_dp), 0.0_dp)
                        ! constrain to interval [0,1]

      cvxy1b = 1.0_dp + c_Hw_slide &
                           *ratio_Hw_slide*ratio_Hw_slide &
                              *(3.0_dp-2.0_dp*ratio_Hw_slide)
                        ! cubic S-shape saturation function

#elif (HYDRO_SLIDE_SAT_FCT==3)

      ratio_Hw_slide = max(min(H_w(j,i)*Hw0_slide_inv, 1.0_dp), 0.0_dp)
                        ! constrain to interval [0,1]

      cvxy1b = 1.0_dp + c_Hw_slide &
                           *ratio_Hw_slide*ratio_Hw_slide*ratio_Hw_slide &
                              *(10.0_dp + ratio_Hw_slide &
                                 *(-15.0_dp+6.0_dp*ratio_Hw_slide))
                        ! quintic S-shape saturation function

#else
      errormsg = ' >>> calc_vxy_b_sia: ' &
                    //'HYDRO_SLIDE_SAT_FCT must be 0, 1, 2 or 3!'
      call error(errormsg)
#endif

      ctau1b = 1.0_dp/(cvxy1b+eps_dp)**p_weert_inv(j,i)

      d_help_b(j,i) = d_help_b(j,i) *cvxy1b
      c_drag(j,i)   = c_drag(j,i)   *ctau1b

#endif

   else   ! mask(j,i) == 1, 2 or 3 away from the grounding line

      d_help_b(j,i) = 0.0_dp
      c_drag(j,i)   = 0.0_dp

   end if

end do
end do

!-------- Computation of vx_b (defined at (i+1/2,j)) --------

do i=0, IMAX-1
do j=1, JMAX-1
   vx_b(j,i) = -0.5_dp*(d_help_b(j,i)+d_help_b(j,i+1))*dzs_dx_aux(j,i)
end do
end do

!-------- Computation of vy_b (defined at (i,j+1/2)) --------

do i=1, IMAX-1
do j=0, JMAX-1
   vy_b(j,i) = -0.5_dp*(d_help_b(j,i)+d_help_b(j+1,i))*dzs_dy_aux(j,i)
end do
end do

!-------- Computation of vx_b_g and vy_b_g (defined at (i,j)) --------

do i=0, IMAX
do j=0, JMAX
   vx_b_g(j,i) = -d_help_b(j,i)*dzs_dxi_g(j,i)
   vy_b_g(j,i) = -d_help_b(j,i)*dzs_deta_g(j,i)
end do
end do

!-------- Limitation of computed vx_b, vy_b, vx_b_g, vy_b_g to the interval
!         [-VH_MAX, VH_MAX] --------

vh_max     = max(VH_MAX, eps_dp)*sec2year
vh_max_inv = 1.0_dp/vh_max

do i=0, IMAX
do j=0, JMAX

   call velocity_limiter_gradual(vx_b(j,i), vh_max, vh_max_inv)
   call velocity_limiter_gradual(vy_b(j,i), vh_max, vh_max_inv)

   call velocity_limiter_gradual(vx_b_g(j,i), vh_max, vh_max_inv)
   call velocity_limiter_gradual(vy_b_g(j,i), vh_max, vh_max_inv)

end do
end do

!-------- Discard basal velocities for HYB_MODE==[1,2] --------

#if ( (DYNAMICS==2) && (HYB_MODE==1 || HYB_MODE==2) )

d_help_b = 0.0_dp
vx_b     = 0.0_dp
vy_b     = 0.0_dp
vx_b_g   = 0.0_dp
vy_b_g   = 0.0_dp

! c_slide and c_drag are not reset because they will be used in the
! computation of the SStA velocity components

!-------- Discard basal velocities for DYNAMICS==3 (DIVA) --------

#elif (DYNAMICS==3)

d_help_b = 0.0_dp
vx_b     = 0.0_dp
vy_b     = 0.0_dp
vx_b_g   = 0.0_dp
vy_b_g   = 0.0_dp

#endif

end subroutine calc_vxy_b_sia

!-------------------------------------------------------------------------------
!> Computation of the shear stresses txz, tyz, the effective shear stress
!! sigma, the depth-averaged fluidity flui_ave_sia, the horizontal
!! velocity vx, vy and the horizontal volume flux qx, qy in the shallow ice
!! approximation.
!-------------------------------------------------------------------------------
subroutine calc_vxy_sia(dzeta_c, dzeta_t)

!$ use omp_lib

use ice_material_properties_m, only : ratefac_c, ratefac_t, ratefac_c_t, creep

implicit none

real(dp), intent(in) :: dzeta_c, dzeta_t

integer(i4b) :: i, j, ij, kc, kt
real(dp) :: avxy3(0:KCMAX), aqxy1(0:KCMAX)
real(dp) :: ctxyz1(0:KCMAX,0:JMAX,0:IMAX), &
            ctxyz2(0:KTMAX,0:JMAX,0:IMAX)
real(dp) :: flui_t(0:KTMAX), flui_c(0:KCMAX)
real(dp) :: cflui0(0:KTMAX), cflui1(0:KCMAX)
real(dp) :: cvxy2(0:KTMAX), cvxy3(0:KCMAX)
real(dp) :: cqxy0(0:KTMAX), cqxy1(0:KCMAX)
real(dp) :: vh_max, vh_max_inv
real(dp) :: flui_min, flui_max, flui_init
real(dp) :: ratio_sl_threshold

!$omp parallel default(shared) private(ij, i, j, kc, kt) &
!$omp private(avxy3, aqxy1) &
!$omp private(flui_t, flui_c, cflui0, cflui1) &
!$omp private(cvxy2, cvxy3, cqxy0, cqxy1) &
!$omp private(vh_max, vh_max_inv) &
!$omp private(flui_min, flui_max, flui_init) &
!$omp private(ratio_sl_threshold)

!-------- Term abbreviations --------

do kc=0, KCMAX
   if (flag_aa_nonzero) then
      avxy3(kc) = aa*eaz_c(kc)/(ea-1.0_dp)*dzeta_c !== aqxy1 no ?
      aqxy1(kc) = aa/(ea-1.0_dp)*eaz_c(kc)*dzeta_c
   else
      avxy3(kc) = dzeta_c
      aqxy1(kc) = dzeta_c
   end if
end do

!-------- Computation of stresses --------

!  ------ Term abbreviations

!$omp do
do ij=1, (IMAX+1)*(JMAX+1)

   i = n2i(ij)   ! i=0...IMAX
   j = n2j(ij)   ! j=0...JMAX

   if ((mask(j,i) == 0).or.flag_grounding_line_2(j,i)) then
                     ! grounded ice, or floating ice at the grounding line

      do kc=0, KCMAX
         ctxyz1(kc,j,i) = RHO*G*H_c(j,i)*(1.0_dp-eaz_c_quotient(kc))
      end do

      if (n_cts(j,i) == 1) then   ! temperate layer

         do kt=0, KTMAX
            ctxyz2(kt,j,i) = RHO*G*H_t(j,i)*(1.0_dp-zeta_t(kt))
         end do

      else   ! cold base (-1), temperate base (0)

         do kt=0, KTMAX
            ctxyz2(kt,j,i) = 0.0_dp
         end do

      end if

   else   ! mask(j,i) == 1, 2 or 3 away from the grounding line

      do kc=0, KCMAX
         ctxyz1(kc,j,i) = 0.0_dp
      end do

      do kt=0, KTMAX
         ctxyz2(kt,j,i) = 0.0_dp
      end do

   end if

end do
!$omp end do

!  ------ Stresses txz, tyz, sigma

!$omp do
do ij=1, (IMAX+1)*(JMAX+1)

   i = n2i(ij)   ! i=0...IMAX
   j = n2j(ij)   ! j=0...JMAX

!    ---- Shear stress txz (defined at (i+1/2,j,kc/t))

   if (flag_sg_x(j,i)) then

      do kc=0, KCMAX
         txz_c(kc,j,i) = -0.5_dp*(ctxyz1(kc,j,i)+ctxyz1(kc,j,i+1)) &
                         *dzs_dx_aux(j,i)
      end do

      do kt=0, KTMAX
         txz_t(kt,j,i) = txz_c(0,j,i) &
                         -0.5_dp*(ctxyz2(kt,j,i)+ctxyz2(kt,j,i+1)) &
                         *dzs_dx_aux(j,i)
      end do

   end if

!    ---- Shear stress tyz (defined at (i,j+1/2,kc/t))

   if (flag_sg_y(j,i)) then

      do kc=0, KCMAX
         tyz_c(kc,j,i) = -0.5_dp*(ctxyz1(kc,j,i)+ctxyz1(kc,j+1,i)) &
                         *dzs_dy_aux(j,i)
      end do

      do kt=0, KTMAX
         tyz_t(kt,j,i) = tyz_c(0,j,i) &
                         -0.5_dp*(ctxyz2(kt,j,i)+ctxyz2(kt,j+1,i)) &
                         *dzs_dy_aux(j,i)
      end do

   end if

!    ---- Effective shear stress sigma (defined at (i,j,kc/t))

   do kc=0, KCMAX

#if !defined(ALLOW_TAPENADE) /* Normal */

      sigma_c(kc,j,i) = ctxyz1(kc,j,i) &
                        *sqrt(dzs_dxi_g(j,i)**2+dzs_deta_g(j,i)**2)

#else /* Tapenade: guarding against non-differentiable sqrt(0) */

      if ( (dzs_dxi_g(j,i)**2+dzs_deta_g(j,i)**2) > 0 ) then
         sigma_c(kc,j,i) = ctxyz1(kc,j,i) &
                           *sqrt(dzs_dxi_g(j,i)**2+dzs_deta_g(j,i)**2)
      else
         sigma_c(kc,j,i) = 0.0_dp 
      end if

#endif /* Normal vs. Tapenade */

   end do

   do kt=0, KTMAX

#if !defined(ALLOW_TAPENADE) /* Normal */

      sigma_t(kt,j,i) = sigma_c(0,j,i) &
                        + ctxyz2(kt,j,i) &
                          *sqrt(dzs_dxi_g(j,i)**2+dzs_deta_g(j,i)**2)

#else /* Tapenade: guarding against non-differentiable sqrt(0) */

      if ( (dzs_dxi_g(j,i)**2+dzs_deta_g(j,i)**2) > 0 ) then
         sigma_t(kt,j,i) = sigma_c(0,j,i) &
                           + ctxyz2(kt,j,i) &
                             *sqrt(dzs_dxi_g(j,i)**2+dzs_deta_g(j,i)**2)
      else
         sigma_t(kt,j,i) = sigma_c(0,j,i) 
      end if 

#endif /* Normal vs. Tapenade */

   end do

end do
!$omp end do

!-------- Computation of the depth-averaged fluidity
!                 (defined on the grid points (i,j)) --------

#if (defined(VISC_MIN) && defined(VISC_MAX))
  flui_min = 1.0_dp/VISC_MAX
  flui_max = 1.0_dp/VISC_MIN
#else
  flui_min = 1.0e-25_dp   ! 1/(Pa s)
  flui_max = 1.0e-10_dp   ! 1/(Pa s)
#endif

#if (defined(VISC_INIT_SSA))
  flui_init = 1.0_dp/VISC_INIT_SSA
#else
  flui_init = 1.0e-15_dp   ! 1/(Pa s)
#endif

!$omp do
do ij=1, (IMAX+1)*(JMAX+1)

   i = n2i(ij)   ! i=0...IMAX
   j = n2j(ij)   ! j=0...JMAX

   if ((mask(j,i) == 0).or.flag_grounding_line_2(j,i)) then
                     ! grounded ice, or floating ice at the grounding line

!  ------ Fluidity, abbreviations

      do kt=0, KTMAX
         flui_t(kt) = 2.0_dp &
                      *enh_t(kt,j,i) &
                      *ratefac_t(omega_t(kt,j,i)) &
                      *creep(sigma_t(kt,j,i))
         cflui0(kt) = H_t(j,i)*flui_t(kt)*dzeta_t
      end do

      do kc=0, KCMAX
         flui_c(kc) = 2.0_dp &
                      *enh_c(kc,j,i) &
#if (CALCMOD==0 || CALCMOD==1 || CALCMOD==-1)
                      *ratefac_c(temp_c(kc,j,i), temp_c_m(kc,j,i)) &
#elif (CALCMOD==2 || CALCMOD==3)
                      *ratefac_c_t(temp_c(kc,j,i), omega_c(kc,j,i), &
                                                   temp_c_m(kc,j,i)) &
#endif
                      *creep(sigma_c(kc,j,i))
         cflui1(kc) = aqxy1(kc)*H_c(j,i)*flui_c(kc)
      end do

!  ------ Depth average

      flui_ave_sia(j,i) = 0.0_dp

      if (n_cts(j,i) == 1) then

         do kt=0, KTMAX-1
            flui_ave_sia(j,i) = flui_ave_sia(j,i)+0.5_dp*(cflui0(kt+1)+cflui0(kt))
         end do

      end if

      do kc=0, KCMAX-1
         flui_ave_sia(j,i) = flui_ave_sia(j,i)+0.5_dp*(cflui1(kc+1)+cflui1(kc))
      end do

      flui_ave_sia(j,i) = flui_ave_sia(j,i)/max(H(j,i), eps_dp)

      flui_ave_sia(j,i) = max(min(flui_ave_sia(j,i), flui_max), flui_min)

   else   ! mask(j,i) == 1, 2 or 3 away from the grounding line

      flui_ave_sia(j,i) = flui_init   ! dummy value

   end if

end do
!$omp end do

!-------- Computation of d_help_c/t
!         (defined on the grid points (i,j,kc/t)) --------

!$omp do
do ij=1, (IMAX+1)*(JMAX+1)

   i = n2i(ij)   ! i=0...IMAX
   j = n2j(ij)   ! j=0...JMAX

   if ((mask(j,i) == 0).or.flag_grounding_line_2(j,i)) then
                     ! grounded ice, or floating ice at the grounding line

!  ------ Abbreviations

      do kt=0, KTMAX
         cvxy2(kt) = 2.0_dp*H_t(j,i) &
                     *enh_t(kt,j,i) &
                     *ratefac_t(omega_t(kt,j,i)) &
                     *creep(sigma_t(kt,j,i)) &
                     *(ctxyz1(0,j,i)+ctxyz2(kt,j,i)) &
                     *dzeta_t
      end do

      do kc=0, KCMAX
         cvxy3(kc) = 2.0_dp*avxy3(kc)*H_c(j,i) &
                     *enh_c(kc,j,i) &
#if (CALCMOD==0 || CALCMOD==1 || CALCMOD==-1)
                     *ratefac_c(temp_c(kc,j,i), temp_c_m(kc,j,i)) &
#elif (CALCMOD==2 || CALCMOD==3)
                     *ratefac_c_t(temp_c(kc,j,i), omega_c(kc,j,i), &
                                                  temp_c_m(kc,j,i)) &
#endif
                     *creep(sigma_c(kc,j,i)) &
                     *ctxyz1(kc,j,i)
      end do

!  ------ d_help_c, d_help_t

      if (n_cts(j,i) == -1) then   ! cold ice base

         do kt=0, KTMAX
            d_help_t(kt,j,i) = d_help_b(j,i)
         end do

         d_help_c(0,j,i) = d_help_t(KTMAX,j,i)

         do kc=0, KCMAX-1
            d_help_c(kc+1,j,i) = d_help_c(kc,j,i) &
                                +0.5_dp*(cvxy3(kc+1)+cvxy3(kc))
         end do

      else if (n_cts(j,i) == 0) then   ! temperate ice base

         do kt=0, KTMAX
            d_help_t(kt,j,i) = d_help_b(j,i)
         end do

         d_help_c(0,j,i) = d_help_t(KTMAX,j,i)

         do kc=0, KCMAX-1
            d_help_c(kc+1,j,i) = d_help_c(kc,j,i) &
                                +0.5_dp*(cvxy3(kc+1)+cvxy3(kc))
         end do

      else   ! n_cts(j,i) == 1, temperate ice layer

         d_help_t(0,j,i) = d_help_b(j,i)

         do kt=0, KTMAX-1
            d_help_t(kt+1,j,i) = d_help_t(kt,j,i) &
                                +0.5_dp*(cvxy2(kt+1)+cvxy2(kt))
         end do

         d_help_c(0,j,i) = d_help_t(KTMAX,j,i)

         do kc=0, KCMAX-1
            d_help_c(kc+1,j,i) = d_help_c(kc,j,i) &
                                +0.5_dp*(cvxy3(kc+1)+cvxy3(kc))
         end do

      end if

   else   ! mask(j,i) == 1, 2 or 3 away from the grounding line

      do kt=0, KTMAX
         d_help_t(kt,j,i) = 0.0_dp
      end do

      do kc=0, KCMAX
         d_help_c(kc,j,i) = 0.0_dp
      end do

   end if

end do
!$omp end do

!-------- Computation of vx_c/t (defined at (i+1/2,j,kc/t))
!                 and of vy_c/t (defined at (i,j+1/2,kc/t)) --------

!$omp do
do ij=1, (IMAX+1)*(JMAX+1)

   i = n2i(ij)   ! i=0...IMAX
   j = n2j(ij)   ! j=0...JMAX

   if (flag_sg_x_inner_y(j,i)) then

      do kt=0, KTMAX
         vx_t(kt,j,i) = -0.5_dp*(d_help_t(kt,j,i)+d_help_t(kt,j,i+1)) &
                        *dzs_dx_aux(j,i)
      end do

      do kc=0, KCMAX
         vx_c(kc,j,i) = -0.5_dp*(d_help_c(kc,j,i)+d_help_c(kc,j,i+1)) &
                        *dzs_dx_aux(j,i)
      end do

   end if

   if (flag_sg_y_inner_x(j,i)) then

      do kt=0, KTMAX
         vy_t(kt,j,i) = -0.5_dp*(d_help_t(kt,j,i)+d_help_t(kt,j+1,i)) &
                        *dzs_dy_aux(j,i)
      end do

      do kc=0, KCMAX
         vy_c(kc,j,i) = -0.5_dp*(d_help_c(kc,j,i)+d_help_c(kc,j+1,i)) &
                        *dzs_dy_aux(j,i)
      end do

   end if

end do
!$omp end do

!-------- Computation of the surface velocities vx_s_g and vy_s_g
!                                              (defined at (i,j)) --------

!$omp do
do ij=1, (IMAX+1)*(JMAX+1)

   i = n2i(ij)   ! i=0...IMAX
   j = n2j(ij)   ! j=0...JMAX

   vx_s_g(j,i) = -d_help_c(KCMAX,j,i)*dzs_dxi_g(j,i)
   vy_s_g(j,i) = -d_help_c(KCMAX,j,i)*dzs_deta_g(j,i)

end do
!$omp end do

!-------- Limitation of computed vx_c/t, vy_c/t, vx_s_g, vy_s_g
!         to the interval [-VH_MAX, VH_MAX] --------

vh_max     = max(VH_MAX, eps_dp)*sec2year
vh_max_inv = 1.0_dp/vh_max

!$omp do
do ij=1, (IMAX+1)*(JMAX+1)

   i = n2i(ij)   ! i=0...IMAX
   j = n2j(ij)   ! j=0...JMAX

   call velocity_limiter_gradual(vx_s_g(j,i), vh_max, vh_max_inv)
   call velocity_limiter_gradual(vy_s_g(j,i), vh_max, vh_max_inv)

   do kt=0, KTMAX
      call velocity_limiter_gradual(vx_t(kt,j,i), vh_max, vh_max_inv)
      call velocity_limiter_gradual(vy_t(kt,j,i), vh_max, vh_max_inv)
   end do

   do kc=0, KCMAX
      call velocity_limiter_gradual(vx_c(kc,j,i), vh_max, vh_max_inv)
      call velocity_limiter_gradual(vy_c(kc,j,i), vh_max, vh_max_inv)
   end do

end do
!$omp end do

!-------- Computation of h_diff
!         (defined on the grid points (i,j)) --------

!$omp do
do ij=1, (IMAX+1)*(JMAX+1)

   i = n2i(ij)   ! i=0...IMAX
   j = n2j(ij)   ! j=0...JMAX

   if ((mask(j,i) == 0).or.flag_grounding_line_2(j,i)) then
                     ! grounded ice, or floating ice at the grounding line

!  ------ Abbreviations

      do kt=0, KTMAX
         cqxy0(kt) = H_t(j,i)*d_help_t(kt,j,i)*dzeta_t
      end do

      do kc=0, KCMAX
         cqxy1(kc) = aqxy1(kc)*H_c(j,i)*d_help_c(kc,j,i)
      end do

!  ------ h_diff

      h_diff(j,i) = 0.0_dp

      if (n_cts(j,i) == 1) then

         do kt=0, KTMAX-1
            h_diff(j,i) = h_diff(j,i)+0.5_dp*(cqxy0(kt+1)+cqxy0(kt))
         end do

      end if

      do kc=0, KCMAX-1
         h_diff(j,i) = h_diff(j,i)+0.5_dp*(cqxy1(kc+1)+cqxy1(kc))
      end do

!  ------ Limitation of h_diff

      if (h_diff(j,i) < HD_MIN) h_diff(j,i) = 0.0_dp
      if (h_diff(j,i) > HD_MAX) h_diff(j,i) = HD_MAX

   else   ! mask(j,i) == 1, 2 or 3 away from the grounding line

      h_diff(j,i) = 0.0_dp

   end if

end do
!$omp end do

!-------- Computation of the horizontal volume flux
!                            and the depth-averaged velocity --------

!$omp do
do ij=1, (IMAX+1)*(JMAX+1)

   i = n2i(ij)   ! i=0...IMAX
   j = n2j(ij)   ! j=0...JMAX

   if (flag_sg_x(j,i)) then

      qx(j,i) = -0.5_dp*(h_diff(j,i)+h_diff(j,i+1))*dzs_dx_aux(j,i)

      if ( (mask(j,i)==0).or.(mask(j,i+1)==0) ) then
                               ! at least one neighbour point is grounded ice

         vx_m(j,i) = qx(j,i) / ( 0.5_dp*(H(j,i)+H(j,i+1)) )

         call velocity_limiter_gradual(vx_m(j,i), vh_max, vh_max_inv)

         ratio_sl_x(j,i) = abs(vx_t(0,j,i)) / max(abs(vx_c(KCMAX,j,i)), eps_dp)

      else 

         vx_m(j,i)       = 0.0_dp
         ratio_sl_x(j,i) = 0.0_dp

      end if

   end if

   if (flag_sg_y(j,i)) then

      qy(j,i) = -0.5_dp*(h_diff(j,i)+h_diff(j+1,i))*dzs_dy_aux(j,i)

      if ( (mask(j,i)==0).or.(mask(j+1,i)==0) ) then
                               ! at least one neighbour point is grounded ice

         vy_m(j,i) = qy(j,i) / ( 0.5_dp*(H(j,i)+H(j+1,i)) )

         call velocity_limiter_gradual(vy_m(j,i), vh_max, vh_max_inv)

         ratio_sl_y(j,i) = abs(vy_t(0,j,i)) / max(abs(vy_c(KCMAX,j,i)), eps_dp)

      else 

         vy_m(j,i)       = 0.0_dp
         ratio_sl_y(j,i) = 0.0_dp

      end if

   end if

   if (flag_inner_point(j,i).and.(mask(j,i) == 0)) then
                                       ! inner point, grounded ice
      ratio_sl(j,i) = 0.25_dp &
                        * (   ratio_sl_x(j,i-1) + ratio_sl_x(j,i) &
                            + ratio_sl_y(j-1,i) + ratio_sl_y(j,i) )
   else
      ratio_sl(j,i) = 0.0_dp
   end if

end do
!$omp end do

!-------- Detection of shelfy stream points --------

!$omp do
do ij=1, (IMAX+1)*(JMAX+1)

   i = n2i(ij)   ! i=0...IMAX
   j = n2j(ij)   ! j=0...JMAX

   flag_shelfy_stream_x(j,i) = .false.
   flag_shelfy_stream_y(j,i) = .false.
   flag_shelfy_stream(j,i)   = .false.

end do
!$omp end do

#if (DYNAMICS==0 || DYNAMICS==1)

ratio_sl_threshold = 1.11e+11_dp   ! dummy value

#elif (DYNAMICS==2)

#if (HYB_MODE==0)
#if (defined(RATIO_SL_THRESH))
ratio_sl_threshold = RATIO_SL_THRESH
#else
ratio_sl_threshold = 0.5_dp   ! default value
#endif
#else
ratio_sl_threshold = 1.11e+11_dp   ! dummy value
#endif

!$omp do
do ij=1, (IMAX+1)*(JMAX+1)

   i = n2i(ij)   ! i=0...IMAX
   j = n2j(ij)   ! j=0...JMAX

   if (flag_sg_x(j,i)) then

#if (HYB_MODE==0)
      if (ratio_sl_x(j,i) > ratio_sl_threshold) &
         flag_shelfy_stream_x(j,i) = .true.
#elif (HYB_MODE==1 || HYB_MODE==2)
      if ( (mask(j,i)==0).or.(mask(j,i+1)==0) ) &
         flag_shelfy_stream_x(j,i) = .true.
#else
      errormsg = ' >>> calc_vxy_sia: HYB_MODE must be 0, 1 or 2!'
      call error(errormsg)
#endif

   end if

   if (flag_sg_y(j,i)) then

#if (HYB_MODE==0)
      if (ratio_sl_y(j,i) > ratio_sl_threshold) &
         flag_shelfy_stream_y(j,i) = .true.
#elif (HYB_MODE==1 || HYB_MODE==2)
      if ( (mask(j,i)==0).or.(mask(j+1,i)==0) ) &
         flag_shelfy_stream_y(j,i) = .true.
#else
      errormsg = ' >>> calc_vxy_sia: HYB_MODE must be 0, 1 or 2!'
      call error(errormsg)
#endif

   end if

   if (flag_inner_point(j,i)) then

      if (mask(j,i) == 0) then   ! grounded ice

         if (     flag_shelfy_stream_x(j,i-1)   &   ! at least
              .or.flag_shelfy_stream_x(j,i)     &   ! one neighbour
              .or.flag_shelfy_stream_y(j-1,i)   &   ! on the staggered grid
              .or.flag_shelfy_stream_y(j,i)   ) &   ! is a shelfy stream point
                  flag_shelfy_stream(j,i) = .true.

      end if

   end if

end do
!$omp end do

#elif (DYNAMICS==3)

ratio_sl_threshold = 1.11e+11_dp   ! dummy value

!$omp do
do ij=1, (IMAX+1)*(JMAX+1)

   i = n2i(ij)   ! i=0...IMAX
   j = n2j(ij)   ! j=0...JMAX

   if (flag_sg_x(j,i)) then
      if ( (mask(j,i)==0).or.(mask(j,i+1)==0) ) &
         flag_shelfy_stream_x(j,i) = .true.
              ! set to true the shelfy stream flag
              ! for the whole grounded ice
              ! (re-used to compute diva in the calc_vxy_ssa)
   end if

   if (flag_sg_y(j,i)) then
      if ( (mask(j,i)==0).or.(mask(j+1,i)==0) ) &
         flag_shelfy_stream_y(j,i) = .true.
              ! set to true the shelfy stream flag
              ! for the whole grounded ice
              ! (re-used to compute diva in the calc_vxy_ssa)
   end if

   if (flag_inner_point(j,i)) then

      if (mask(j,i) == 0) then   ! grounded ice

         if (     flag_shelfy_stream_x(j,i-1)   &   ! at least
              .or.flag_shelfy_stream_x(j,i)     &   ! one neighbour
              .or.flag_shelfy_stream_y(j-1,i)   &   ! on the staggered grid
              .or.flag_shelfy_stream_y(j,i)   ) &   ! is a shelfy stream point
                  flag_shelfy_stream(j,i) = .true.

      end if

   end if

end do
!$omp end do

#else

errormsg = ' >>> calc_vxy_sia: DYNAMICS must be between 0 and 3!'
call error(errormsg)

#endif

!-------- Loop over (i,j) --------

!$omp do
do ij=1, (IMAX+1)*(JMAX+1)

   i = n2i(ij)   ! i=0...IMAX
   j = n2j(ij)   ! j=0...JMAX

!-------- Save mean (depth-averaged) horizontal velocities from SIA --------

   vx_m_sia(j,i) = vx_m(j,i)
   vy_m_sia(j,i) = vy_m(j,i)

!-------- Initialisation of the variable q_gl_g
!         (volume flux across the grounding line, to be
!         computed in the routine calc_vxy_ssa
!         if ice shelves are present)

   q_gl_g(j,i) = 0.0_dp

end do
!$omp end do

!$omp end parallel

end subroutine calc_vxy_sia

!-------------------------------------------------------------------------------
!> Computation of the horizontal velocity vx, vy, the horizontal volume flux
!! qx, qy etc. for static ice.
!-------------------------------------------------------------------------------
subroutine calc_vxy_static()

use calc_pressure_water_bas_m

implicit none

real(dp) :: flui_init

#if (defined(VISC_INIT_SSA))
  flui_init = 1.0_dp/VISC_INIT_SSA
#else
  flui_init = 1.0e-15_dp   ! 1/(Pa s)
#endif

c_slide = 0.0_dp
p_weert = 0
q_weert = 0

call calc_pressure_water_bas()   ! compute p_b_w

d_help_b = 0.0_dp
c_drag   = 0.0_dp

vx_b   = 0.0_dp
vy_b   = 0.0_dp
vx_b_g = 0.0_dp
vy_b_g = 0.0_dp

txz_c = 0.0_dp
txz_t = 0.0_dp

tyz_c = 0.0_dp
tyz_t = 0.0_dp

sigma_c = 0.0_dp
sigma_t = 0.0_dp

flui_ave_sia = flui_init
de_ssa       = 0.0_dp
vis_int_g    = 0.0_dp

d_help_t = 0.0_dp
d_help_c = 0.0_dp

vx_c = 0.0_dp
vy_c = 0.0_dp

vx_t = 0.0_dp
vy_t = 0.0_dp

vx_s_g = 0.0_dp
vy_s_g = 0.0_dp

h_diff = 0.0_dp

qx = 0.0_dp
qy = 0.0_dp

vx_m = 0.0_dp
vy_m = 0.0_dp

ratio_sl_x = 0.0_dp
ratio_sl_y = 0.0_dp

flag_shelfy_stream_x = .false.
flag_shelfy_stream_y = .false.
flag_shelfy_stream   = .false.

vx_m_sia = 0.0_dp
vy_m_sia = 0.0_dp

q_gl_g = 0.0_dp

end subroutine calc_vxy_static

!-------------------------------------------------------------------------------
!> Computation of the horizontal velocity vx, vy, the horizontal volume flux
!! qx, qy and the flux across the grounding line q_gl_g in the shallow shelf
!! approximation (SSA) or the shelfy stream approximation (SStA).
!-------------------------------------------------------------------------------
subroutine calc_vxy_ssa(dxi, deta, dzeta_c, dzeta_t)

#if (DYNAMICS==2 || DYNAMICS==3)
  use calc_enhance_m, only : calc_enhance_hybrid_weighted
#endif

implicit none

real(dp), intent(in) :: dxi, deta, dzeta_c, dzeta_t

integer(i4b) :: i, j, kc, kt, m
integer(i4b) :: iter_ssa_min, iter_ssa_max
real(dp), dimension(0:JMAX,0:IMAX) :: vx_m_prev, vy_m_prev
real(dp) :: tol_ssa, rel_ssa
real(dp) :: res_vxy_m_ssa_1, res_vxy_m_ssa_2, res_vxy_m_ssa
real(dp) :: dxi_inv, deta_inv
real(dp) :: visc_init
real(dp) :: vh_max, vh_max_inv
real(dp) :: ratio_sl_threshold, ratio_help
real(dp), dimension(0:JMAX,0:IMAX) :: weigh_ssta_sia_x, weigh_ssta_sia_y
real(dp), dimension(0:JMAX,0:IMAX) :: weigh_ssta_sia
real(dp) :: qx_gl_g, qy_gl_g
logical, dimension(0:JMAX,0:IMAX) :: flag_calc_vxy_ssa_x, flag_calc_vxy_ssa_y
real(dp) :: v_ref, v_ref_sq_inv
real(dp) :: v_b_sq
real(dp) :: pi_inv


pi_inv   = 1.0_dp/pi

#if (MARGIN==3 || DYNAMICS==2 || DYNAMICS==3)

!-------- Parameters for the relaxation scheme --------

#if (defined(TOL_ITER_SSA))
   tol_ssa = TOL_ITER_SSA   ! tolerance of iterations 
#else
   tol_ssa = 0.025_dp       ! default value
#endif

#if (defined(N_ITER_SSA))
   iter_ssa_max = max(N_ITER_SSA, 1)   ! max. number of iterations
#else
   iter_ssa_max = 25                   ! default value
#endif

#if (defined(N_ITER_SSA_MIN))
   iter_ssa_min = max(N_ITER_SSA_MIN, 1)   ! min. number of iterations
#else
   iter_ssa_min = 2                        ! default value
#endif

if (iter_ssa_min > iter_ssa_max) then
   errormsg = ' >>> calc_vxy_ssa: N_ITER_SSA_MIN > N_ITER_SSA_MAX not allowed!'
   call error(errormsg)
end if

#if (defined(RELAX_FACT_SSA))
   rel_ssa = RELAX_FACT_SSA   ! relaxation factor
#else
   rel_ssa = 0.7_dp           ! default value
#endif

#if (defined(VISC_INIT_SSA))
  visc_init = VISC_INIT_SSA
#else
  visc_init = 1.0e+15_dp   ! Pa s
#endif

vh_max     = max(VH_MAX, eps_dp)*sec2year
vh_max_inv = 1.0_dp/vh_max

write(6,'(10x,a)') 'calc_vxy_ssa:'

!-------- Iterations --------

res_vxy_m_ssa = 1.11e+11_dp   ! initial, very large value of the residual

m=0

do while ( (m < iter_ssa_min) &
           .or. &
           ((m < iter_ssa_max).and.(res_vxy_m_ssa > tol_ssa)) )

   m = m+1

   write(6,'(13x,a,i0,a)', advance='no') 'Iter ', m, ': '

!  ------ Save velocities from previous iteration

   vx_m_prev = vx_m_ssa
   vy_m_prev = vy_m_ssa

!  ------ Depth-integrated viscosity vis_int_g

   if (m > 1) then

      call calc_vis_ssa(dxi, deta, dzeta_c, dzeta_t)

   else   ! (m == 1, first iteration)

#if (!defined(ITER_INIT_SSA) || ITER_INIT_SSA==1)
      vis_int_g = H*visc_init
                  ! constant viscosity times ice thickness
#elif (ITER_INIT_SSA==2)
      vis_int_g = H*vis_ave_g
                  ! previous depth-averaged viscosity times
                  ! ice thickness
#else
      errormsg = ' >>> calc_vxy_ssa: ITER_INIT_SSA must be 1 or 2!'
      call error(errormsg)
#endif

   end if

#if (DYNAMICS==3)
   call calc_F_int_DIVA( dzeta_c, dzeta_t)
#endif

!  ------ Horizontal velocity vx_m_ssa, vy_m_ssa

   flag_calc_vxy_ssa_x = .false.   ! initialization
   flag_calc_vxy_ssa_y = .false.   ! initialization

   call calc_vxy_ssa_matrix(dxi, deta, &
                            flag_calc_vxy_ssa_x, flag_calc_vxy_ssa_y, &
                            dzeta_c, dzeta_t)

   do i=0, IMAX
   do j=0, JMAX
      call velocity_limiter_gradual(vx_m_ssa(j,i), vh_max, vh_max_inv)
      call velocity_limiter_gradual(vy_m_ssa(j,i), vh_max, vh_max_inv)
   end do
   end do

!  ------ Relaxation scheme

   if (m > 1) then
      vx_m_ssa = rel_ssa*vx_m_ssa + (1.0_dp-rel_ssa)*vx_m_prev
      vy_m_ssa = rel_ssa*vy_m_ssa + (1.0_dp-rel_ssa)*vy_m_prev
   end if

!  ------ Residual

   res_vxy_m_ssa_1 &
      = sqrt( sum((vx_m_ssa-vx_m_prev)*(vx_m_ssa-vx_m_prev), &
                  mask=flag_calc_vxy_ssa_x) &
             +sum((vy_m_ssa-vy_m_prev)*(vy_m_ssa-vy_m_prev), &
                  mask=flag_calc_vxy_ssa_y) )
   res_vxy_m_ssa_2 &
      = sqrt( sum((vx_m_ssa+vx_m_prev)*(vx_m_ssa+vx_m_prev), &
                  mask=flag_calc_vxy_ssa_x) &
             +sum((vy_m_ssa+vy_m_prev)*(vy_m_ssa+vy_m_prev), &
                  mask=flag_calc_vxy_ssa_y) )

   res_vxy_m_ssa = 2.0_dp*res_vxy_m_ssa_1/max(res_vxy_m_ssa_2, eps_dp)

   write(6,'(a,es9.2)') 'res =', res_vxy_m_ssa

end do

!  ------ Depth-integrated viscosity vis_int_g

call calc_vis_ssa(dxi, deta, dzeta_c, dzeta_t)

!-------- 3-D velocities, basal velocities and volume flux --------

#if (DYNAMICS==0 || DYNAMICS==1 || DYNAMICS==3)

ratio_sl_threshold = 1.11e+11_dp   ! dummy value
ratio_help         = 1.11e+11_dp   ! dummy value
v_ref              = 1.11e+11_dp   ! dummy value
v_ref_sq_inv       = 1.11e+11_dp   ! dummy value

#elif (DYNAMICS==2)

#if (HYB_MODE==0)
#if (defined(RATIO_SL_THRESH))
ratio_sl_threshold = RATIO_SL_THRESH
#else
ratio_sl_threshold = 0.5_dp   ! default value
#endif
ratio_help = 1.0_dp/(1.0_dp-ratio_sl_threshold)
#else
ratio_sl_threshold = 1.11e+11_dp   ! dummy value
ratio_help         = 1.11e+11_dp   ! dummy value
#endif

#if (HYB_MODE==1)
#if (defined(HYB_REF_SPEED))
v_ref = real(HYB_REF_SPEED,dp)*sec2year
#else
v_ref = 30.0_dp*sec2year   ! default value
#endif
v_ref_sq_inv = 1.0_dp/(v_ref*v_ref)
#else
v_ref        = 1.11e+11_dp   ! dummy value
v_ref_sq_inv = 1.11e+11_dp   ! dummy value
#endif

#else

errormsg = ' >>> calc_vxy_ssa: DYNAMICS must be between 0 and 3!'
call error(errormsg)

#endif

weigh_ssta_sia_x = 0.0_dp
weigh_ssta_sia_y = 0.0_dp


!  ------ x-component

do i=0, IMAX-1
do j=0, JMAX

   if (flag_shelfy_stream_x(j,i)) then   ! shelfy stream

#if (DYNAMICS==2 && HYB_MODE==0)
             ! Sum of weighted sliding SIA and weighted SStA
             ! (by R. Greve)

      weigh_ssta_sia_x(j,i) = (ratio_sl_x(j,i)-ratio_sl_threshold)*ratio_help

      weigh_ssta_sia_x(j,i) = max(min(weigh_ssta_sia_x(j,i), 1.0_dp), 0.0_dp)
                              ! constrain to interval [0,1]

#if (SSTA_SIA_WEIGH_FCT==0)

      ! stick to the linear function set above

#elif (SSTA_SIA_WEIGH_FCT==1)

      weigh_ssta_sia_x(j,i) = weigh_ssta_sia_x(j,i)*weigh_ssta_sia_x(j,i) &
                                   *(3.0_dp-2.0_dp*weigh_ssta_sia_x(j,i))
                              ! make transition smooth (cubic function)

#elif (SSTA_SIA_WEIGH_FCT==2)

      weigh_ssta_sia_x(j,i) = weigh_ssta_sia_x(j,i)*weigh_ssta_sia_x(j,i) &
                                                   *weigh_ssta_sia_x(j,i) &
                                   *(10.0_dp + weigh_ssta_sia_x(j,i) &
                                       *(-15.0_dp+6.0_dp*weigh_ssta_sia_x(j,i)))
                              ! make transition even smoother (quintic function)

#else
      errormsg = ' >>> calc_vxy_ssa: SSTA_SIA_WEIGH_FCT must be 0, 1 or 2!'
      call error(errormsg)
#endif

      do kt=0, KTMAX
         vx_t(kt,j,i) = weigh_ssta_sia_x(j,i)*vx_m_ssa(j,i) &
                        + (1.0_dp-weigh_ssta_sia_x(j,i))*vx_t(kt,j,i)
      end do

      do kc=0, KCMAX
         vx_c(kc,j,i) = weigh_ssta_sia_x(j,i)*vx_m_ssa(j,i) &
                        + (1.0_dp-weigh_ssta_sia_x(j,i))*vx_c(kc,j,i)
      end do

      vx_b(j,i) = vx_t(0,j,i)

      vx_m(j,i) = weigh_ssta_sia_x(j,i)*vx_m_ssa(j,i) &
                  + (1.0_dp-weigh_ssta_sia_x(j,i))*vx_m_sia(j,i)

#elif (DYNAMICS==2 && HYB_MODE==1)
               ! Sum of weighted non-sliding SIA and full SStA
               ! (by J. Bernales)

      weigh_ssta_sia_x(j,i) = (2.0_dp*pi_inv) &
                                 * atan( (vx_m_ssa(j,i)*vx_m_ssa(j,i)) &
                                                            *v_ref_sq_inv )

      do kt=0, KTMAX
         vx_t(kt,j,i) = vx_m_ssa(j,i) &
                           + (1.0_dp-weigh_ssta_sia_x(j,i))*vx_t(kt,j,i)
         vx_t(kt,j,i) = max(vx_t(kt,j,i), -vh_max)
         vx_t(kt,j,i) = min(vx_t(kt,j,i),  vh_max)
      end do

      do kc=0, KCMAX
         vx_c(kc,j,i) = vx_m_ssa(j,i) &
                           + (1.0_dp-weigh_ssta_sia_x(j,i))*vx_c(kc,j,i)
         vx_c(kc,j,i) = max(vx_c(kc,j,i), -vh_max)
         vx_c(kc,j,i) = min(vx_c(kc,j,i),  vh_max)
      end do

      vx_b(j,i) = vx_t(0,j,i)

      vx_m(j,i) = vx_m_ssa(j,i) &
                     + (1.0_dp-weigh_ssta_sia_x(j,i))*vx_m_sia(j,i)
      vx_m(j,i) = max(vx_m(j,i), -vh_max)
      vx_m(j,i) = min(vx_m(j,i),  vh_max)

#elif (DYNAMICS==2 && HYB_MODE==2)   /* Pure SStA (no SIA) */

      do kt=0, KTMAX
         vx_t(kt,j,i) = vx_m_ssa(j,i)
      end do

      do kc=0, KCMAX
         vx_c(kc,j,i) = vx_m_ssa(j,i)
      end do

      vx_b(j,i) = vx_t(0,j,i)

      vx_m(j,i) = vx_m_ssa(j,i)

#elif (DYNAMICS==3)   /* DIVA */ 

!  ------ x-component (for DIVA)

      vx_m(j,i) = vx_m_ssa(j,i)   ! we used the ssa solver to compute vxy_m

      tau_bx(j,i) = 0.5_dp*(beta_eff(j,i) + beta_eff(j,i+1)) * vx_m(j,i)
                    ! shear basal drag in the x-direction (on the staggered grid)

!    ---- Basal velocity
      if (C_SLIDE == 0) then
         vx_b(j,i) = 0.0_dp
      else
         vx_b(j,i) = vx_m(j,i) - (tau_bx(j,i) * F_2_x(j,i))
      endif
!    ---- 3D velocity

      do kt=0, KTMAX
         vx_t(kt,j,i) = vx_b(j,i) + tau_bx(j,i) * 0.5_dp*(F_1_t_g(kt,j,i)+F_1_t_g(kt,j,i+1))
      end do

      do kc=0, KCMAX
         vx_c(kc,j,i) = vx_b(j,i) + tau_bx(j,i) * 0.5_dp*(F_1_c_g(kc,j,i)+F_1_c_g(kc,j,i+1))
      end do

#endif

      qx(j,i) = vx_m(j,i) * 0.5_dp*(H(j,i)+H(j,i+1))

   else if (flag_calc_vxy_ssa_x(j,i)) then   ! floating ice

      do kt=0, KTMAX
         vx_t(kt,j,i) = vx_m_ssa(j,i)
      end do

      do kc=0, KCMAX
         vx_c(kc,j,i) = vx_m_ssa(j,i)
      end do

      vx_b(j,i) = vx_m_ssa(j,i)

      vx_m(j,i) = vx_m_ssa(j,i)

      qx(j,i) = vx_m(j,i) * 0.5_dp*(H(j,i)+H(j,i+1))

#if (DYNAMICS==3)
      tau_bx(j,i) = 0.0_dp
#endif

!  else
!     In all other cases, the depth-averaged velocities vx_m_ssa(j,i) computed
!     by the SSA/SStA solver are discarded.

   end if

end do
end do

!  ------ y-component

do i=0, IMAX
do j=0, JMAX-1

   if (flag_shelfy_stream_y(j,i)) then   ! shelfy stream

#if (DYNAMICS==2 && HYB_MODE==0)
             ! Sum of weighted sliding SIA and weighted SStA
             ! (by R. Greve)

      weigh_ssta_sia_y(j,i) = (ratio_sl_y(j,i)-ratio_sl_threshold)*ratio_help

      weigh_ssta_sia_y(j,i) = max(min(weigh_ssta_sia_y(j,i), 1.0_dp), 0.0_dp)
                              ! constrain to interval [0,1]

#if (SSTA_SIA_WEIGH_FCT==0)

      ! stick to the linear function set above

#elif (SSTA_SIA_WEIGH_FCT==1)

      weigh_ssta_sia_y(j,i) = weigh_ssta_sia_y(j,i)*weigh_ssta_sia_y(j,i) &
                                   *(3.0_dp-2.0_dp*weigh_ssta_sia_y(j,i))
                              ! make transition smooth (cubic function)

#elif (SSTA_SIA_WEIGH_FCT==2)

      weigh_ssta_sia_y(j,i) = weigh_ssta_sia_y(j,i)*weigh_ssta_sia_y(j,i) &
                                                   *weigh_ssta_sia_y(j,i) &
                                   *(10.0_dp + weigh_ssta_sia_y(j,i) &
                                       *(-15.0_dp+6.0_dp*weigh_ssta_sia_y(j,i)))
                              ! make transition even smoother (quintic function)

#else
      errormsg = ' >>> calc_vxy_ssa: SSTA_SIA_WEIGH_FCT must be 0, 1 or 2!'
      call error(errormsg)
#endif

      do kt=0, KTMAX
         vy_t(kt,j,i) = weigh_ssta_sia_y(j,i)*vy_m_ssa(j,i) &
                        + (1.0_dp-weigh_ssta_sia_y(j,i))*vy_t(kt,j,i)
      end do

      do kc=0, KCMAX
         vy_c(kc,j,i) = weigh_ssta_sia_y(j,i)*vy_m_ssa(j,i) &
                        + (1.0_dp-weigh_ssta_sia_y(j,i))*vy_c(kc,j,i)
      end do

      vy_b(j,i) = vy_t(0,j,i)

      vy_m(j,i) = weigh_ssta_sia_y(j,i)*vy_m_ssa(j,i) &
                  + (1.0_dp-weigh_ssta_sia_y(j,i))*vy_m_sia(j,i)

#elif (DYNAMICS==2 && HYB_MODE==1)
               ! Sum of weighted non-sliding SIA and full SStA
               ! (by J. Bernales)

      weigh_ssta_sia_y(j,i) = (2.0_dp*pi_inv) &
                                 * atan( (vy_m_ssa(j,i)*vy_m_ssa(j,i)) &
                                                            *v_ref_sq_inv )

      do kt=0, KTMAX
         vy_t(kt,j,i) = vy_m_ssa(j,i) &
                           + (1.0_dp-weigh_ssta_sia_y(j,i))*vy_t(kt,j,i)
         vy_t(kt,j,i) = max(vy_t(kt,j,i), -vh_max)
         vy_t(kt,j,i) = min(vy_t(kt,j,i),  vh_max)
      end do

      do kc=0, KCMAX
         vy_c(kc,j,i) = vy_m_ssa(j,i) &
                           + (1.0_dp-weigh_ssta_sia_y(j,i))*vy_c(kc,j,i)
         vy_c(kc,j,i) = max(vy_c(kc,j,i), -vh_max)
         vy_c(kc,j,i) = min(vy_c(kc,j,i),  vh_max)
      end do

      vy_b(j,i) = vy_t(0,j,i)

      vy_m(j,i) = vy_m_ssa(j,i) &
                     + (1.0_dp-weigh_ssta_sia_y(j,i))*vy_m_sia(j,i)
      vy_m(j,i) = max(vy_m(j,i), -vh_max)
      vy_m(j,i) = min(vy_m(j,i),  vh_max)

#elif (DYNAMICS==2 && HYB_MODE==2)   /* Pure SStA (no SIA) */

      do kt=0, KTMAX
         vy_t(kt,j,i) = vy_m_ssa(j,i)
      end do

      do kc=0, KCMAX
         vy_c(kc,j,i) = vy_m_ssa(j,i)
      end do

      vy_b(j,i) = vy_t(0,j,i)

      vy_m(j,i) = vy_m_ssa(j,i)

#elif (DYNAMICS==3)   /* DIVA */

!  ------ y-component (for DIVA)

      vy_m(j,i) = vy_m_ssa(j,i)   ! we used the ssa solver to compute vxy_m

      tau_by(j,i) =  0.5_dp*(beta_eff(j,i) + beta_eff(j+1,i)) * vy_m(j,i)
                    ! shear basal drag in the y-direction (on the staggered grid)


!    ---- Basal velocity
      if (C_SLIDE == 0) then
         vy_b(j,i) = 0.0_dp
      else
         vy_b(j,i) = vy_m(j,i) - (tau_by(j,i) * F_2_y(j,i))
      endif
!    ---- 3D velocity

      do kt=0, KTMAX
         vy_t(kt,j,i) = vy_b(j,i) + tau_by(j,i) * 0.5_dp*(F_1_t_g(kt,j,i)+F_1_t_g(kt,j+1,i))
      end do

      do kc=0, KCMAX
         vy_c(kc,j,i) = vy_b(j,i) + tau_by(j,i) * 0.5_dp*(F_1_c_g(kc,j,i)+F_1_c_g(kc,j+1,i))
      end do

#endif

      qy(j,i) = vy_m(j,i) * 0.5_dp*(H(j,i)+H(j+1,i))

   else if (flag_calc_vxy_ssa_y(j,i)) then   ! floating ice

      do kt=0, KTMAX
         vy_t(kt,j,i) = vy_m_ssa(j,i)
      end do

      do kc=0, KCMAX
         vy_c(kc,j,i) = vy_m_ssa(j,i)
      end do

      vy_b(j,i) = vy_m_ssa(j,i)

      vy_m(j,i) = vy_m_ssa(j,i)

      qy(j,i) = vy_m(j,i) * 0.5_dp*(H(j,i)+H(j+1,i))

#if (DYNAMICS==3)
      tau_by(j,i) = 0.0_dp
#endif

!  else
!     In all other cases, the depth-averaged velocities vy_m_ssa(j,i) computed
!     by the SSA/SStA solver are discarded.

   end if

end do
end do

!-------- Weighted flow enhancement factor --------

weigh_ssta_sia = 0.0_dp

#if (DYNAMICS==2 || DYNAMICS==3)

#if (DYNAMICS==2 && HYB_MODE==0)
             ! Sum of weighted sliding SIA and weighted SStA
             ! (by R. Greve)

do i=0, IMAX
do j=0, JMAX

   if (flag_shelfy_stream(j,i)) then   ! shelfy stream

      weigh_ssta_sia(j,i) = (ratio_sl(j,i)-ratio_sl_threshold)*ratio_help

      weigh_ssta_sia(j,i) = max(min(weigh_ssta_sia(j,i), 1.0_dp), 0.0_dp)
                              ! constrain to interval [0,1]

#if (SSTA_SIA_WEIGH_FCT==0)

      ! stick to the linear function set above

#elif (SSTA_SIA_WEIGH_FCT==1)

      weigh_ssta_sia(j,i) = weigh_ssta_sia(j,i)*weigh_ssta_sia(j,i) &
                                 *(3.0_dp-2.0_dp*weigh_ssta_sia(j,i))
                            ! make transition smooth (cubic function)

#elif (SSTA_SIA_WEIGH_FCT==2)

      weigh_ssta_sia(j,i) = weigh_ssta_sia(j,i)*weigh_ssta_sia(j,i) &
                                               *weigh_ssta_sia(j,i) &
                                 *(10.0_dp + weigh_ssta_sia(j,i) &
                                     *(-15.0_dp+6.0_dp*weigh_ssta_sia(j,i)))
                            ! make transition even smoother (quintic function)

#else
      errormsg = ' >>> calc_vxy_ssa: SSTA_SIA_WEIGH_FCT must be 0, 1 or 2!'
      call error(errormsg)
#endif

   end if

end do
end do

#elif (DYNAMICS==2 && HYB_MODE==1)
               ! Sum of weighted non-sliding SIA and full SStA
               ! (by J. Bernales)

do i=1, IMAX-1
do j=1, JMAX-1

   if (flag_shelfy_stream(j,i)) then   ! shelfy stream

      weigh_ssta_sia(j,i) = 0.25_dp * &
                            (   weigh_ssta_sia_x(j,i-1) + weigh_ssta_sia_x(j,i) &
                              + weigh_ssta_sia_y(j-1,i) + weigh_ssta_sia_y(j,i) )

      weigh_ssta_sia(j,i) = max(min(weigh_ssta_sia(j,i), 1.0_dp), 0.0_dp)
                            ! constrain to interval [0,1] (just in case :shrug: )

   end if

end do
end do

#elif (DYNAMICS==2 && HYB_MODE==2)   /* Pure SStA (no SIA) */

do i=1, IMAX-1
do j=1, JMAX-1
   if (flag_shelfy_stream(j,i)) weigh_ssta_sia(j,i) = 1.0_dp   ! shelfy stream
end do
end do

#elif (DYNAMICS==3)   /* DIVA */

do i=1, IMAX-1
do j=1, JMAX-1
   if (flag_shelfy_stream(j,i)) weigh_ssta_sia(j,i) = 1.0_dp   ! shelfy stream
                                !%% \!/ This needs to be reconsidered. \!/
end do
end do

#endif

if (flag_calc_temp) call calc_enhance_hybrid_weighted(weigh_ssta_sia)

#endif   /* DYNAMICS==2 or 3 */

!-------- Surface and basal velocities vx_s_g vy_s_g, vx_b_g vy_b_g
!                                                (defined at (i,j)) --------

do i=1, IMAX-1
do j=1, JMAX-1

   if (flag_shelfy_stream(j,i)) then   ! shelfy stream

      vx_s_g(j,i) = 0.5_dp*(vx_c(KCMAX,j,i-1)+vx_c(KCMAX,j,i))
      vx_b_g(j,i) = 0.5_dp*(vx_b(      j,i-1)+vx_b(      j,i))

      vy_s_g(j,i) = 0.5_dp*(vy_c(KCMAX,j-1,i)+vy_c(KCMAX,j,i))
      vy_b_g(j,i) = 0.5_dp*(vy_b(      j-1,i)+vy_b(      j,i))

   else if (mask(j,i)==3) then   ! floating ice

      vx_s_g(j,i) = 0.5_dp*(vx_m(j,i-1)+vx_m(j,i))
      vx_b_g(j,i) = vx_s_g(j,i)

      vy_s_g(j,i) = 0.5_dp*(vy_m(j-1,i)+vy_m(j,i))
      vy_b_g(j,i) = vy_s_g(j,i)

   end if
#if (DYNAMICS==3)
   !save the norm of the basal velocity to use in calc_vxy_ssa_matrix :
   vb_sq_prev(j,i) = (vx_b_g(j,i))**2  &
                   + (vy_b_g(j,i))**2
#endif
end do
end do

!-------- Magnitude of the basal shear stress (drag) --------

do i=1, IMAX-1
do j=1, JMAX-1

   if (mask(j,i) == 3) then   ! floating ice

      tau_b(j,i) = 0.0_dp

   else if (flag_shelfy_stream(j,i)) then   ! shelfy stream
! this computation seems weird as vxy_b_g is defined on main grid already
!      v_b_sq =   (0.5_dp*(vx_b_g(j,i)+vx_b_g(j,i-1)))**2  &
!               + (0.5_dp*(vy_b_g(j,i)+vy_b_g(j-1,i)))**2
! i would rather have done :
      v_b_sq =   (vx_b_g(j,i))**2  &
               + (vy_b_g(j,i))**2

#if !defined(ALLOW_TAPENADE) /* Normal */

      tau_b(j,i) = c_drag(j,i) * sqrt(v_b_sq)**p_weert_inv(j,i)

#else /* Tapenade: guarding against non-differentiable sqrt(0) */

      if (v_b_sq > 0) then
         tau_b(j,i) = c_drag(j,i) * sqrt(v_b_sq)**p_weert_inv(j,i)
      else
         tau_b(j,i) = 0.0_dp
      end if

#endif /* Normal vs. Tapenade */

   end if

end do
end do

!-------- Computation of the flux across the grounding line q_gl_g

do i=1, IMAX-1
do j=1, JMAX-1

   if ( flag_grounding_line_1(j,i) ) then   ! grounding line

      qx_gl_g = 0.5_dp*(qx(j,i)+qx(j,i-1))
      qy_gl_g = 0.5_dp*(qy(j,i)+qy(j-1,i))

#if !defined(ALLOW_TAPENADE) /* Normal */

      q_gl_g(j,i) = sqrt(qx_gl_g*qx_gl_g+qy_gl_g*qy_gl_g)

#else /* Tapenade: guarding against non-differentiable sqrt(0) */

      if ( (qx_gl_g*qx_gl_g+qy_gl_g*qy_gl_g) > 0 ) then
         q_gl_g(j,i) = sqrt(qx_gl_g*qx_gl_g+qy_gl_g*qy_gl_g)
      else
         q_gl_g(j,i) = 0.0_dp 
      end if

#endif /* Normal vs. Tapenade */

   end if

end do
end do

#else

errormsg = ' >>> calc_vxy_ssa: Only to be called for' &
         //                    end_of_line &
         //'                   MARGIN==3, DYNAMICS==2 or DYNAMICS==3!'
call error(errormsg)

#endif

end subroutine calc_vxy_ssa

!-------------------------------------------------------------------------------
!> Solution of the system of linear equations for the horizontal velocities
!! vx_m_ssa, vy_m_ssa in the SSA/SStA.
!-------------------------------------------------------------------------------
subroutine calc_vxy_ssa_matrix(dxi, deta, &
                               flag_calc_vxy_ssa_x, flag_calc_vxy_ssa_y, &
                               dzeta_c, dzeta_t)

#if (MARGIN==3 || DYNAMICS==2 || DYNAMICS==3)
#if defined(ALLOW_TAPENADE) /* Tapenade */
use sico_maths_m
#endif /* Tapenade */
#endif

implicit none

real(dp), intent(in) :: dxi, deta
real(dp), intent(in) :: dzeta_c, dzeta_t
logical, dimension(0:JMAX,0:IMAX), intent(inout) :: &
                            flag_calc_vxy_ssa_x, flag_calc_vxy_ssa_y

integer(i4b) :: i, j, k, n
integer(i4b) :: i1, j1
real(dp) :: inv_dxi, inv_deta, inv_dxi_deta, inv_dxi2, inv_deta2
real(dp) :: factor_rhs_1, factor_rhs_2, factor_rhs_3a, factor_rhs_3b
real(dp) :: rhosw_rho_ratio
real(dp) :: H_mid, zl_mid, z_sl_mid
real(dp), dimension(0:JMAX,0:IMAX) :: vis_int_sgxy
character(len=256) :: ch_solver_set_option

#if (MARGIN==3 || DYNAMICS==2 || DYNAMICS==3)

#if !defined(ALLOW_TAPENADE) /* Normal */

LIS_INTEGER :: ierr
LIS_INTEGER :: nc, nr
LIS_INTEGER :: lin_iter
LIS_MATRIX  :: lgs_a
LIS_VECTOR  :: lgs_b, lgs_x
LIS_SOLVER  :: solver

LIS_INTEGER, parameter                 :: nmax   =  2*(IMAX+1)*(JMAX+1)
LIS_INTEGER, parameter                 :: n_sprs = 20*(IMAX+1)*(JMAX+1)
LIS_INTEGER, allocatable, dimension(:) :: lgs_a_ptr, lgs_a_index
LIS_SCALAR,  allocatable, dimension(:) :: lgs_a_value, lgs_b_value, lgs_x_value

#else /* Tapenade */

integer(i4b) :: ierr
integer(i4b) :: nc, nr
integer(i4b) :: lin_iter

integer(i4b), parameter                 :: nmax   =  2*(IMAX+1)*(JMAX+1)
integer(i4b), parameter                 :: n_sprs = 20*(IMAX+1)*(JMAX+1)
integer(i4b), dimension(nmax+1)         :: lgs_a_ptr
real(dp), dimension(n_sprs)             :: lgs_a_value
real(dp), dimension(n_sprs)             :: lgs_a_index
integer(i4b), dimension(n_sprs)         :: lgs_a_index_pass
real(dp), dimension(nmax)               :: lgs_b_value, lgs_x_value

#endif /* Normal vs. Tapenade */

!-------- Abbreviations --------

inv_dxi      = 1.0_dp/dxi
inv_deta     = 1.0_dp/deta
inv_dxi_deta = 1.0_dp/(dxi*deta)
inv_dxi2     = 1.0_dp/(dxi*dxi)
inv_deta2    = 1.0_dp/(deta*deta)

factor_rhs_1  = RHO*G
factor_rhs_2  = 0.5_dp*RHO*G*(RHO_SW-RHO)/RHO_SW
factor_rhs_3a = 0.5_dp*RHO*G
factor_rhs_3b = 0.5_dp*RHO_SW*G

rhosw_rho_ratio = RHO_SW/RHO

!-------- Depth-integrated viscosity on the staggered grid
!                                       [at (i+1/2,j+1/2)] --------

vis_int_sgxy = 0.0_dp   ! initialisation

do i=0, IMAX-1
do j=0, JMAX-1

   k=0

   if ((mask(j,i)==0).or.(mask(j,i)==3)) then
      k = k+1                              ! floating or grounded ice
      vis_int_sgxy(j,i) = vis_int_sgxy(j,i) + vis_int_g(j,i)
   end if

   if ((mask(j,i+1)==0).or.(mask(j,i+1)==3)) then
      k = k+1                                  ! floating or grounded ice
      vis_int_sgxy(j,i) = vis_int_sgxy(j,i) + vis_int_g(j,i+1)
   end if

   if ((mask(j+1,i)==0).or.(mask(j+1,i)==3)) then
      k = k+1                                  ! floating or grounded ice
      vis_int_sgxy(j,i) = vis_int_sgxy(j,i) + vis_int_g(j+1,i)
   end if

   if ((mask(j+1,i+1)==0).or.(mask(j+1,i+1)==3)) then
      k = k+1                                      ! floating or grounded ice
      vis_int_sgxy(j,i) = vis_int_sgxy(j,i) + vis_int_g(j+1,i+1)
   end if

   if (k>0) vis_int_sgxy(j,i) = vis_int_sgxy(j,i)/real(k,dp)

end do
end do

!-------- Basal drag parameter (for shelfy stream) --------

beta_drag = 0.0_dp   ! initialisation
beta_eff = 0.0_dp    ! initialisation

do i=1, IMAX-1
do j=1, JMAX-1

   if (flag_shelfy_stream(j,i)) then

# if (DYNAMICS==2)

      beta_drag(j,i) = c_drag(j,i) &
                     / sqrt( (   (0.5_dp*(vx_m_ssa(j,i)+vx_m_ssa(j,i-1)))**2  &
                               + (0.5_dp*(vy_m_ssa(j,i)+vy_m_ssa(j-1,i)))**2 ) &
                             + eps_dp**2 ) &
                                     **(1.0_dp-p_weert_inv(j,i))

#elif (DYNAMICS==3)
      beta_drag(j,i) = c_drag(j,i) &
                     / ((sqrt(vb_sq_prev(j,i)+ eps_dp**2))**(1.0_dp-p_weert_inv(j,i)))
      ! because tau_b = beta_drag * v_b = c_drag * |v_b|**-(1-1/p) * v_b

      beta_eff(j,i) = beta_drag(j,i)/(1.0_dp + beta_drag(j,i)*F_2_g(j,i))
#endif

   end if

end do
end do
!-------- Assembly of the system of linear equations
!                         (matrix storage: compressed sparse row CSR) --------

#if !defined(ALLOW_TAPENADE) /* Normal */
allocate(lgs_a_value(n_sprs), lgs_a_index(n_sprs), lgs_a_ptr(nmax+1))
allocate(lgs_b_value(nmax), lgs_x_value(nmax))
#endif /* Normal */

lgs_a_value = 0.0_dp
#if !defined(ALLOW_TAPENADE) /* Normal */
lgs_a_index = 0
#else /* Tapenade */
lgs_a_index = 0.0_dp
#endif /* Normal vs. Tapenade */
lgs_a_ptr   = 0

lgs_b_value = 0.0_dp
lgs_x_value = 0.0_dp

lgs_a_ptr(1) = 1

k = 0

do n=1, nmax-1, 2

   i = n2i((n+1)/2)
   j = n2j((n+1)/2)

!  ------ Equations for vx_m_ssa (at (i+1/2,j))

   nr = n   ! row counter

   if ( (i /= IMAX).and.(j /= 0).and.(j /= JMAX) ) then
      ! inner point on the staggered grid in x-direction

      H_mid    = 0.5_dp*(H(j,i)   +H(j,i+1))
      zl_mid   = 0.5_dp*(zl(j,i)  +zl(j,i+1))
      z_sl_mid = 0.5_dp*(z_sl(j,i)+z_sl(j,i+1))

      if ( &
           ( (mask(j,i)==3).and.(mask(j,i+1)==3) ) &
           .or. &
           ( flag_grounding_line_1(j,i).and.flag_grounding_line_2(j,i+1) &
             .and.(H_mid < (z_sl_mid-zl_mid)*rhosw_rho_ratio) ) &
           .or. &
           ( flag_grounding_line_2(j,i).and.flag_grounding_line_1(j,i+1) &
             .and.(H_mid < (z_sl_mid-zl_mid)*rhosw_rho_ratio) ) &
         ) then
           ! both neighbours are floating ice
           !   or
           ! one neighbour is floating ice and the other is grounded ice
           ! (grounding line)
           ! and floating conditions are satisfied;
           ! discretization of the x-component of the PDE

         flag_calc_vxy_ssa_x(j,i) = .true.

         flag_shelfy_stream_x(j,i) = .false.
                                   ! make sure not to treat as shelfy stream

         nc = 2*ij2n(j,i-1)-1
                  ! smallest nc (column counter), for vx_m(j,i-1)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 4.0_dp*inv_dxi2*vis_int_g(j,i)
         lgs_a_index(k) = nc

         nc = 2*ij2n(j-1,i)-1
                  ! next nc (column counter), for vx_m(j-1,i)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = inv_deta2*vis_int_sgxy(j-1,i)
         lgs_a_index(k) = nc

         nc = 2*ij2n(j-1,i)
                  ! next nc (column counter), for vy_m(j-1,i)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = inv_dxi_deta &
                                 *(2.0_dp*vis_int_g(j,i)+vis_int_sgxy(j-1,i))
         lgs_a_index(k) = nc

         nc = 2*ij2n(j,i)-1
                  ! next nc (column counter), for vx_m(j,i)
         if (nc /= nr) then   ! (diagonal element)
            errormsg = ' >>> calc_vxy_ssa_matrix: ' &
                          //'Check for diagonal element failed!'
            call error(errormsg)
         end if
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = -4.0_dp*inv_dxi2 &
                                 *(vis_int_g(j,i+1)+vis_int_g(j,i)) &
                          -inv_deta2 &
                                 *(vis_int_sgxy(j,i)+vis_int_sgxy(j-1,i))
         lgs_a_index(k) = nc

         nc = 2*ij2n(j,i)
                  ! next nc (column counter), for vy_m(j,i)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = -inv_dxi_deta &
                                 *(2.0_dp*vis_int_g(j,i)+vis_int_sgxy(j,i))
         lgs_a_index(k) = nc

         nc = 2*ij2n(j+1,i)-1
                  ! next nc (column counter), for vx_m(j+1,i)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = inv_deta2*vis_int_sgxy(j,i)
         lgs_a_index(k) = nc

         nc = 2*ij2n(j-1,i+1)
                  ! next nc (column counter), for vy_m(j-1,i+1)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = -inv_dxi_deta &
                                 *(2.0_dp*vis_int_g(j,i+1)+vis_int_sgxy(j-1,i))
         lgs_a_index(k) = nc

         nc = 2*ij2n(j,i+1)-1
                  ! next nc (column counter), for vx_m(j,i+1)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 4.0_dp*inv_dxi2*vis_int_g(j,i+1)
         lgs_a_index(k) = nc

         nc = 2*ij2n(j,i+1)
                  ! largest nc (column counter), for vy_m(j,i+1)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = inv_dxi_deta &
                                 *(2.0_dp*vis_int_g(j,i+1)+vis_int_sgxy(j,i))
         lgs_a_index(k) = nc

         lgs_b_value(nr) = factor_rhs_1*H_mid*dzs_dx_aux(j,i)

         lgs_x_value(nr) = vx_m_ssa(j,i)
 
      else if (flag_shelfy_stream_x(j,i)) then
           ! shelfy stream (as determined by routine calc_vxy_sia)

         flag_calc_vxy_ssa_x(j,i) = .true.

#if (!defined(BC_SSA_LTGF) || BC_SSA_LTGF==1)

         if ( &
              ( ( flag_grounded_front_b_1(j,i) &
                       .and.flag_grounded_front_b_2(j,i+1) ) &
                .or. &
                ( flag_grounded_front_b_2(j,i) &
                       .and.flag_grounded_front_b_1(j,i+1) ) ) &
              .and. &
              ( zl_mid < z_sl_mid ) &
            ) then
            ! one neighbour is grounded ice and the other is ocean
            ! (ocean-terminating grounded front)

#elif (BC_SSA_LTGF==2)

         if ( &
              ( flag_grounded_front_b_1(j,i) &
                     .and.flag_grounded_front_b_2(j,i+1) ) &
              .or. &
              ( flag_grounded_front_b_2(j,i) &
                     .and.flag_grounded_front_b_1(j,i+1) ) &
            ) then
            ! one neighbour is grounded ice and the other is ocean
            ! (ocean-terminating grounded front)

#else
         errormsg = ' >>> calc_vxy_ssa_matrix: ' &
                       //'BC_SSA_LTGF must be either 1 or 2!'
         call error(errormsg)
#endif

            if (flag_grounded_front_b_1(j,i)) then
               i1 = i     ! grounded ice marker
            else   ! flag_grounded_front_b_1(j,i+1)==.true.
               i1 = i+1   ! grounded ice marker
            end if

            if (.not.( flag_grounded_front_b_2(j,i1-1) &
                       .and. &
                       flag_grounded_front_b_2(j,i1+1) ) ) then
               ! discretization of the x-component of the BC

               nc = 2*ij2n(j,i1-1)-1
                        ! smallest nc (column counter), for vx_m(j,i1-1)
               k  = k+1
               lgs_a_value(k) = -4.0_dp*inv_dxi*vis_int_g(j,i1)
               lgs_a_index(k) = nc

               nc = 2*ij2n(j-1,i1)
                        ! next nc (column counter), for vy_m(j-1,i1)
               k  = k+1
               lgs_a_value(k) = -2.0_dp*inv_deta*vis_int_g(j,i1)
               lgs_a_index(k) = nc

               nc = 2*ij2n(j,i1)-1
                        ! next nc (column counter), for vx_m(j,i1)
               k  = k+1
               lgs_a_value(k) = 4.0_dp*inv_dxi*vis_int_g(j,i1)
               lgs_a_index(k) = nc

               nc = 2*ij2n(j,i1)
                        ! largest nc (column counter), for vy_m(j,i1)
               k  = k+1
               lgs_a_value(k) = 2.0_dp*inv_deta*vis_int_g(j,i1)
               lgs_a_index(k) = nc

               lgs_b_value(nr) = factor_rhs_3a &
                                    *H(j,i1)**2 &
                               - factor_rhs_3b &
                                    *(max((z_sl(j,i1)-zb(j,i1)), 0.0_dp))**2

               lgs_x_value(nr) = vx_m_ssa(j,i)

            else   !      (flag_grounded_front_b_2(j,i1-1)==.true.)
                   ! .and.(flag_grounded_front_b_2(j,i1+1)==.true.);
                   ! velocity assumed to be zero

               k  = k+1
               lgs_a_value(k) = 1.0_dp   ! diagonal element only
               lgs_a_index(k) = nr

               lgs_b_value(nr) = 0.0_dp

               lgs_x_value(nr) = 0.0_dp

            end if

#if (BC_SSA_LTGF==2)

         else if ( &
              ( flag_grounded_front_a_1(j,i) &
                     .and.flag_grounded_front_a_2(j,i+1) ) &
              .or. &
              ( flag_grounded_front_a_2(j,i) &
                     .and.flag_grounded_front_a_1(j,i+1) ) &
            ) then
            ! one neighbour is grounded ice and the other is ice-free land
            ! (land-terminating grounded front)

            if (flag_grounded_front_a_1(j,i)) then
               i1 = i     ! grounded ice marker
            else   ! flag_grounded_front_a_1(j,i+1)==.true.
               i1 = i+1   ! grounded ice marker
            end if

            if (.not.( flag_grounded_front_a_2(j,i1-1) &
                       .and. &
                       flag_grounded_front_a_2(j,i1+1) ) ) then
               ! discretization of the x-component of the BC

               nc = 2*ij2n(j,i1-1)-1
                        ! smallest nc (column counter), for vx_m(j,i1-1)
               k  = k+1
               lgs_a_value(k) = -4.0_dp*inv_dxi*vis_int_g(j,i1)
               lgs_a_index(k) = nc

               nc = 2*ij2n(j-1,i1)
                        ! next nc (column counter), for vy_m(j-1,i1)
               k  = k+1
               lgs_a_value(k) = -2.0_dp*inv_deta*vis_int_g(j,i1)
               lgs_a_index(k) = nc

               nc = 2*ij2n(j,i1)-1
                        ! next nc (column counter), for vx_m(j,i1)
               k  = k+1
               lgs_a_value(k) = 4.0_dp*inv_dxi*vis_int_g(j,i1)
               lgs_a_index(k) = nc

               nc = 2*ij2n(j,i1)
                        ! largest nc (column counter), for vy_m(j,i1)
               k  = k+1
               lgs_a_value(k) = 2.0_dp*inv_deta*vis_int_g(j,i1)
               lgs_a_index(k) = nc

               lgs_b_value(nr) = factor_rhs_3a &
                                    *H(j,i1)**2

               lgs_x_value(nr) = vx_m_ssa(j,i)

            else   !      (flag_grounded_front_a_2(j,i1-1)==.true.)
                   ! .and.(flag_grounded_front_a_2(j,i1+1)==.true.);
                   ! velocity assumed to be zero

               k  = k+1
               lgs_a_value(k) = 1.0_dp   ! diagonal element only
               lgs_a_index(k) = nr

               lgs_b_value(nr) = 0.0_dp

               lgs_x_value(nr) = 0.0_dp

            end if

#endif

         else
            ! inner shelfy stream
            !   or
            ! one neighbour is floating ice and the other is grounded ice
            ! (grounding line)
            ! and floating conditions are not satisfied
#if (!defined(BC_SSA_LTGF) || BC_SSA_LTGF==1)
            !   or
            ! land-terminating grounded front
#endif

            nc = 2*ij2n(j,i-1)-1
                     ! smallest nc (column counter), for vx_m(j,i-1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = 4.0_dp*inv_dxi2*vis_int_g(j,i)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j-1,i)-1
                     ! next nc (column counter), for vx_m(j-1,i)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_deta2*vis_int_sgxy(j-1,i)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j-1,i)
                     ! next nc (column counter), for vy_m(j-1,i)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_dxi_deta &
                                    *(2.0_dp*vis_int_g(j,i)+vis_int_sgxy(j-1,i))
            lgs_a_index(k) = nc

            nc = 2*ij2n(j,i)-1
                     ! next nc (column counter), for vx_m(j,i)
            if (nc /= nr) then   ! (diagonal element)
               errormsg = ' >>> calc_vxy_ssa_matrix: ' &
                             //'Check for diagonal element failed!'
               call error(errormsg)
            end if
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
#if (DYNAMICS==2)
            lgs_a_value(k) = -4.0_dp*inv_dxi2 &
                                    *(vis_int_g(j,i+1)+vis_int_g(j,i)) &
                             -inv_deta2 &
                                    *(vis_int_sgxy(j,i)+vis_int_sgxy(j-1,i)) &
                             -0.5_dp*(beta_drag(j,i+1)+beta_drag(j,i))
#elif (DYNAMICS==3)
            lgs_a_value(k) = -4.0_dp*inv_dxi2 &
                                    *(vis_int_g(j,i+1)+vis_int_g(j,i)) &
                             -inv_deta2 &
                                    *(vis_int_sgxy(j,i)+vis_int_sgxy(j-1,i)) &
                             -0.5_dp*(beta_eff(j,i+1)+beta_eff(j,i))
#endif
            lgs_a_index(k) = nc

            nc = 2*ij2n(j,i)
                     ! next nc (column counter), for vy_m(j,i)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = -inv_dxi_deta &
                                    *(2.0_dp*vis_int_g(j,i)+vis_int_sgxy(j,i))
            lgs_a_index(k) = nc

            nc = 2*ij2n(j+1,i)-1
                     ! next nc (column counter), for vx_m(j+1,i)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_deta2*vis_int_sgxy(j,i)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j-1,i+1)
                     ! next nc (column counter), for vy_m(j-1,i+1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = -inv_dxi_deta &
                                  *(2.0_dp*vis_int_g(j,i+1)+vis_int_sgxy(j-1,i))
            lgs_a_index(k) = nc

            nc = 2*ij2n(j,i+1)-1
                     ! next nc (column counter), for vx_m(j,i+1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = 4.0_dp*inv_dxi2*vis_int_g(j,i+1)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j,i+1)
                     ! largest nc (column counter), for vy_m(j,i+1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_dxi_deta &
                                    *(2.0_dp*vis_int_g(j,i+1)+vis_int_sgxy(j,i))
            lgs_a_index(k) = nc

            lgs_b_value(nr) = factor_rhs_1*H_mid*dzs_dx_aux(j,i)

            lgs_x_value(nr) = vx_m_ssa(j,i)

         end if

      else if ( &
                ( flag_grounding_line_1(j,i).and.flag_grounding_line_2(j,i+1) &
                  .and.(H_mid >= (z_sl_mid-zl_mid)*rhosw_rho_ratio) ) &
                .or. &
                ( flag_grounding_line_2(j,i).and.flag_grounding_line_1(j,i+1) &
                  .and.(H_mid >= (z_sl_mid-zl_mid)*rhosw_rho_ratio) ) &
              ) then 
              ! one neighbour is floating ice and the other is grounded ice
              ! (grounding line)
              ! and floating conditions are not satisfied;
              ! velocity taken from the SIA solution for grounded ice

         flag_calc_vxy_ssa_x(j,i) = .false.

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_dp   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = vx_m_sia(j,i)

         lgs_x_value(nr) = vx_m_sia(j,i)

      else if ( &
                ( (mask(j,i)==3).and.(mask(j,i+1)==1) ) &
                .or. &
                ( (mask(j,i)==1).and.(mask(j,i+1)==3) ) &
              ) then
              ! one neighbour is floating ice and the other is ice-free land;
              ! velocity assumed to be zero

         flag_calc_vxy_ssa_x(j,i) = .true.

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_dp   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = 0.0_dp

         lgs_x_value(nr) = 0.0_dp

      else if ( &
                ( flag_calving_front_1(j,i).and.flag_calving_front_2(j,i+1) ) &
                .or. &
                ( flag_calving_front_2(j,i).and.flag_calving_front_1(j,i+1) ) &
              ) then
              ! one neighbour is floating ice and the other is ocean
              ! (calving front)

         flag_calc_vxy_ssa_x(j,i) = .true.

         if (flag_calving_front_1(j,i)) then
            i1 = i     ! floating ice marker
         else   ! flag_calving_front_1(j,i+1)==.true.
            i1 = i+1   ! floating ice marker
         end if

         if (.not.( flag_calving_front_2(j,i1-1) &
                    .and. &
                    flag_calving_front_2(j,i1+1) ) ) then
            ! discretization of the x-component of the BC

            nc = 2*ij2n(j,i1-1)-1
                     ! smallest nc (column counter), for vx_m(j,i1-1)
            k  = k+1
            lgs_a_value(k) = -4.0_dp*inv_dxi*vis_int_g(j,i1)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j-1,i1)
                     ! next nc (column counter), for vy_m(j-1,i1)
            k  = k+1
            lgs_a_value(k) = -2.0_dp*inv_deta*vis_int_g(j,i1)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j,i1)-1
                     ! next nc (column counter), for vx_m(j,i1)
            k  = k+1
            lgs_a_value(k) = 4.0_dp*inv_dxi*vis_int_g(j,i1)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j,i1)
                     ! largest nc (column counter), for vy_m(j,i1)
            k  = k+1
            lgs_a_value(k) = 2.0_dp*inv_deta*vis_int_g(j,i1)
            lgs_a_index(k) = nc

            lgs_b_value(nr) = factor_rhs_2 &
                                 *H(j,i1)**2

            lgs_x_value(nr) = vx_m_ssa(j,i)

         else   !      (flag_calving_front_2(j,i1-1)==.true.)
                ! .and.(flag_calving_front_2(j,i1+1)==.true.);
                ! velocity assumed to be zero

            k  = k+1
            lgs_a_value(k) = 1.0_dp   ! diagonal element only
            lgs_a_index(k) = nr

            lgs_b_value(nr) = 0.0_dp

            lgs_x_value(nr) = 0.0_dp

         end if

      else if ( (mask(j,i)==0).or.(mask(j,i+1)==0) ) then
           ! neither neighbour is floating ice, but at least one neighbour is
           ! grounded ice; velocity taken from the SIA solution for grounded ice

         flag_calc_vxy_ssa_x(j,i) = .false.

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_dp   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = vx_m_sia(j,i)

         lgs_x_value(nr) = vx_m_sia(j,i)

      else   ! neither neighbour is floating or grounded ice,
             ! velocity assumed to be zero

         flag_calc_vxy_ssa_x(j,i) = .false.

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_dp   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = 0.0_dp

         lgs_x_value(nr) = 0.0_dp

      end if

   else   ! boundary condition, velocity assumed to be zero

      flag_calc_vxy_ssa_x(j,i) = .false.

      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
      lgs_a_value(k) = 1.0_dp   ! diagonal element only
      lgs_a_index(k) = nr

      lgs_b_value(nr) = 0.0_dp

      lgs_x_value(nr) = 0.0_dp

   end if

   lgs_a_ptr(nr+1) = k+1   ! row is completed, store index to next row

!  ------ Equations for vy_m_ssa (at (i,j+1/2))

   nr = n+1   ! row counter

   if ( (j /= JMAX).and.(i /= 0).and.(i /= IMAX) ) then
      ! inner point on the staggered grid in y-direction

      H_mid    = 0.5_dp*(H(j,i)   +H(j+1,i))
      zl_mid   = 0.5_dp*(zl(j,i)  +zl(j+1,i))
      z_sl_mid = 0.5_dp*(z_sl(j,i)+z_sl(j+1,i))
   
      if ( &
           ( (mask(j,i)==3).and.(mask(j+1,i)==3) ) &
           .or. &
           ( flag_grounding_line_1(j,i).and.flag_grounding_line_2(j+1,i) &
             .and.(H_mid < (z_sl_mid-zl_mid)*rhosw_rho_ratio) ) &
           .or. &
           ( flag_grounding_line_2(j,i).and.flag_grounding_line_1(j+1,i) &
             .and.(H_mid < (z_sl_mid-zl_mid)*rhosw_rho_ratio) ) &
         ) then
           ! both neighbours are floating ice
           !   or
           ! one neighbour is floating ice and the other is grounded ice
           ! (grounding line)
           ! and floating conditions are satisfied;
           ! discretization of the y-component of the PDE

         flag_calc_vxy_ssa_y(j,i) = .true.

         flag_shelfy_stream_y(j,i) = .false.
                                   ! make sure not to treat as shelfy stream

         nc = 2*ij2n(j,i-1)-1
                  ! smallest nc (column counter), for vx_m(j,i-1)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = inv_dxi_deta &
                                 *(2.0_dp*vis_int_g(j,i)+vis_int_sgxy(j,i-1))
         lgs_a_index(k) = nc

         nc = 2*ij2n(j,i-1)
                  ! next nc (column counter), for vy_m(j,i-1)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = inv_dxi2*vis_int_sgxy(j,i-1)
         lgs_a_index(k) = nc

         nc = 2*ij2n(j+1,i-1)-1
                  ! next nc (column counter), for vx_m(j+1,i-1)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = -inv_dxi_deta &
                                 *(2.0_dp*vis_int_g(j+1,i)+vis_int_sgxy(j,i-1))
         lgs_a_index(k) = nc

         nc = 2*ij2n(j-1,i)
                  ! next nc (column counter), for vy_m(j-1,i)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 4.0_dp*inv_deta2*vis_int_g(j,i)
         lgs_a_index(k) = nc

         nc = 2*ij2n(j,i)-1
                  ! next nc (column counter), for vx_m(j,i)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = -inv_dxi_deta &
                                 *(2.0_dp*vis_int_g(j,i)+vis_int_sgxy(j,i))
         lgs_a_index(k) = nc

         nc = 2*ij2n(j,i)
                  ! next nc (column counter), for vy_m(j,i)
         if (nc /= nr) then   ! (diagonal element)
            errormsg = ' >>> calc_vxy_ssa_matrix: ' &
                          //'Check for diagonal element failed!'
            call error(errormsg)
         end if
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = -4.0_dp*inv_deta2 &
                                 *(vis_int_g(j+1,i)+vis_int_g(j,i)) &
                          -inv_dxi2 &
                                 *(vis_int_sgxy(j,i)+vis_int_sgxy(j,i-1))
         lgs_a_index(k) = nc

         nc = 2*ij2n(j+1,i)-1
                  ! next nc (column counter), for vx_m(j+1,i)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = inv_dxi_deta &
                                 *(2.0_dp*vis_int_g(j+1,i)+vis_int_sgxy(j,i))
         lgs_a_index(k) = nc

         nc = 2*ij2n(j+1,i)
                  ! next nc (column counter), for vy_m(j+1,i)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 4.0_dp*inv_deta2*vis_int_g(j+1,i)
         lgs_a_index(k) = nc

         nc = 2*ij2n(j,i+1)
                  ! largest nc (column counter), for vy_m(j,i+1)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = inv_dxi2*vis_int_sgxy(j,i)
         lgs_a_index(k) = nc

         lgs_b_value(nr) = factor_rhs_1*H_mid*dzs_dy_aux(j,i)

         lgs_x_value(nr) = vy_m_ssa(j,i)
 
      else if (flag_shelfy_stream_y(j,i)) then
           ! shelfy stream (as determined by routine calc_vxy_sia)

         flag_calc_vxy_ssa_y(j,i) = .true.

#if (!defined(BC_SSA_LTGF) || BC_SSA_LTGF==1)

         if ( &
              ( ( flag_grounded_front_b_1(j,i) &
                       .and.flag_grounded_front_b_2(j+1,i) ) &
                .or. &
                ( flag_grounded_front_b_2(j,i) &
                       .and.flag_grounded_front_b_1(j+1,i) ) ) &
              .and. &
              ( zl_mid < z_sl_mid ) &
            ) then
            ! one neighbour is grounded ice and the other is ocean
            ! (ocean-terminating grounded front)

#elif (BC_SSA_LTGF==2)

         if ( &
              ( flag_grounded_front_b_1(j,i) &
                     .and.flag_grounded_front_b_2(j+1,i) ) &
              .or. &
              ( flag_grounded_front_b_2(j,i) &
                     .and.flag_grounded_front_b_1(j+1,i) ) &
            ) then
            ! one neighbour is grounded ice and the other is ocean
            ! (ocean-terminating grounded front)

#else
         errormsg = ' >>> calc_vxy_ssa_matrix: ' &
                       //'BC_SSA_LTGF must be either 1 or 2!'
         call error(errormsg)
#endif

            if (flag_grounded_front_b_1(j,i)) then
               j1 = j     ! grounded ice marker
            else   ! flag_grounded_front_b_1(j+1,i)==.true.
               j1 = j+1   ! grounded ice marker
            end if

            if (.not.( flag_grounded_front_b_2(j1-1,i) &
                       .and. &
                       flag_grounded_front_b_2(j1+1,i) ) ) then
               ! discretization of the y-component of the BC

               nc = 2*ij2n(j1,i-1)-1
                        ! smallest nc (column counter), for vx_m(j1,i-1)
               k  = k+1
               lgs_a_value(k) = -2.0_dp*inv_dxi*vis_int_g(j1,i)
               lgs_a_index(k) = nc

               nc = 2*ij2n(j1-1,i)
                        ! next nc (column counter), for vy_m(j1-1,i)
               k  = k+1
               lgs_a_value(k) = -4.0_dp*inv_deta*vis_int_g(j1,i)
               lgs_a_index(k) = nc

               nc = 2*ij2n(j1,i)-1
                        ! next nc (column counter), for vx_m(j1,i)
               k  = k+1
               lgs_a_value(k) = 2.0_dp*inv_dxi*vis_int_g(j1,i)
               lgs_a_index(k) = nc

               nc = 2*ij2n(j1,i)
                        ! largest nc (column counter), for vy_m(j1,i)
               k  = k+1
               lgs_a_value(k) = 4.0_dp*inv_deta*vis_int_g(j1,i)
               lgs_a_index(k) = nc

               lgs_b_value(nr) = factor_rhs_3a &
                                    *H(j1,i)**2 &
                               - factor_rhs_3b &
                                    *(max((z_sl(j1,i)-zb(j1,i)), 0.0_dp))**2

               lgs_x_value(nr) = vy_m_ssa(j,i)

            else   !      (flag_grounded_front_b_2(j1-1,i)==.true.)
                   ! .and.(flag_grounded_front_b_2(j1+1,i)==.true.);
                   ! velocity assumed to be zero

               k  = k+1
               lgs_a_value(k) = 1.0_dp   ! diagonal element only
               lgs_a_index(k) = nr

               lgs_b_value(nr) = 0.0_dp

               lgs_x_value(nr) = 0.0_dp

            end if

#if (BC_SSA_LTGF==2)

         else if ( &
              ( flag_grounded_front_a_1(j,i) &
                     .and.flag_grounded_front_a_2(j+1,i) ) &
              .or. &
              ( flag_grounded_front_a_2(j,i) &
                     .and.flag_grounded_front_a_1(j+1,i) ) &
            ) then
            ! one neighbour is grounded ice and the other is ice-free land
            ! (land-terminating grounded front)

            if (flag_grounded_front_a_1(j,i)) then
               j1 = j     ! grounded ice marker
            else   ! flag_grounded_front_a_1(j+1,i)==.true.
               j1 = j+1   ! grounded ice marker
            end if

            if (.not.( flag_grounded_front_a_2(j1-1,i) &
                       .and. &
                       flag_grounded_front_a_2(j1+1,i) ) ) then
               ! discretization of the y-component of the BC

               nc = 2*ij2n(j1,i-1)-1
                        ! smallest nc (column counter), for vx_m(j1,i-1)
               k  = k+1
               lgs_a_value(k) = -2.0_dp*inv_dxi*vis_int_g(j1,i)
               lgs_a_index(k) = nc

               nc = 2*ij2n(j1-1,i)
                        ! next nc (column counter), for vy_m(j1-1,i)
               k  = k+1
               lgs_a_value(k) = -4.0_dp*inv_deta*vis_int_g(j1,i)
               lgs_a_index(k) = nc

               nc = 2*ij2n(j1,i)-1
                        ! next nc (column counter), for vx_m(j1,i)
               k  = k+1
               lgs_a_value(k) = 2.0_dp*inv_dxi*vis_int_g(j1,i)
               lgs_a_index(k) = nc

               nc = 2*ij2n(j1,i)
                        ! largest nc (column counter), for vy_m(j1,i)
               k  = k+1
               lgs_a_value(k) = 4.0_dp*inv_deta*vis_int_g(j1,i)
               lgs_a_index(k) = nc

               lgs_b_value(nr) = factor_rhs_3a &
                                    *H(j1,i)**2

               lgs_x_value(nr) = vy_m_ssa(j,i)

            else   !      (flag_grounded_front_a_2(j1-1,i)==.true.)
                   ! .and.(flag_grounded_front_a_2(j1+1,i)==.true.);
                   ! velocity assumed to be zero

               k  = k+1
               lgs_a_value(k) = 1.0_dp   ! diagonal element only
               lgs_a_index(k) = nr

               lgs_b_value(nr) = 0.0_dp

               lgs_x_value(nr) = 0.0_dp

            end if

#endif

         else
            ! inner shelfy stream
            !   or
            ! one neighbour is floating ice and the other is grounded ice
            ! (grounding line)
            ! and floating conditions are not satisfied
#if (!defined(BC_SSA_LTGF) || BC_SSA_LTGF==1)
            !   or
            ! land-terminating grounded front
#endif

            nc = 2*ij2n(j,i-1)-1
                     ! smallest nc (column counter), for vx_m(j,i-1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_dxi_deta &
                                    *(2.0_dp*vis_int_g(j,i)+vis_int_sgxy(j,i-1))
            lgs_a_index(k) = nc

            nc = 2*ij2n(j,i-1)
                     ! next nc (column counter), for vy_m(j,i-1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_dxi2*vis_int_sgxy(j,i-1)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j+1,i-1)-1
                     ! next nc (column counter), for vx_m(j+1,i-1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = -inv_dxi_deta &
                                  *(2.0_dp*vis_int_g(j+1,i)+vis_int_sgxy(j,i-1))
            lgs_a_index(k) = nc

            nc = 2*ij2n(j-1,i)
                     ! next nc (column counter), for vy_m(j-1,i)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = 4.0_dp*inv_deta2*vis_int_g(j,i)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j,i)-1
                     ! next nc (column counter), for vx_m(j,i)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = -inv_dxi_deta &
                                    *(2.0_dp*vis_int_g(j,i)+vis_int_sgxy(j,i))
            lgs_a_index(k) = nc

            nc = 2*ij2n(j,i)
                     ! next nc (column counter), for vy_m(j,i)
            if (nc /= nr) then   ! (diagonal element)
               errormsg = ' >>> calc_vxy_ssa_matrix: ' &
                             //'Check for diagonal element failed!'
               call error(errormsg)
            end if
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'

#if (DYNAMICS==2)
            lgs_a_value(k) = -4.0_dp*inv_deta2 &
                                    *(vis_int_g(j+1,i)+vis_int_g(j,i)) &
                             -inv_dxi2 &
                                    *(vis_int_sgxy(j,i)+vis_int_sgxy(j,i-1)) &
                             -0.5_dp*(beta_drag(j+1,i)+beta_drag(j,i))
#elif (DYNAMICS==3)
            lgs_a_value(k) = -4.0_dp*inv_deta2 &
                                    *(vis_int_g(j+1,i)+vis_int_g(j,i)) &
                             -inv_dxi2 &
                                    *(vis_int_sgxy(j,i)+vis_int_sgxy(j,i-1)) &
                             -0.5_dp*(beta_eff(j+1,i)+beta_eff(j,i))
#endif
            lgs_a_index(k) = nc

            nc = 2*ij2n(j+1,i)-1
                     ! next nc (column counter), for vx_m(j+1,i)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_dxi_deta &
                                    *(2.0_dp*vis_int_g(j+1,i)+vis_int_sgxy(j,i))
            lgs_a_index(k) = nc

            nc = 2*ij2n(j+1,i)
                     ! next nc (column counter), for vy_m(j+1,i)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = 4.0_dp*inv_deta2*vis_int_g(j+1,i)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j,i+1)
                     ! largest nc (column counter), for vy_m(j,i+1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_dxi2*vis_int_sgxy(j,i)
            lgs_a_index(k) = nc

            lgs_b_value(nr) = factor_rhs_1*H_mid*dzs_dy_aux(j,i)

            lgs_x_value(nr) = vy_m_ssa(j,i)

         end if

      else if ( &
                ( flag_grounding_line_1(j,i).and.flag_grounding_line_2(j+1,i) &
                  .and.(H_mid >= (z_sl_mid-zl_mid)*rhosw_rho_ratio) ) &
                .or. &
                ( flag_grounding_line_2(j,i).and.flag_grounding_line_1(j+1,i) &
                  .and.(H_mid >= (z_sl_mid-zl_mid)*rhosw_rho_ratio) ) &
              ) then
              ! one neighbour is floating ice and the other is grounded ice
              ! (grounding line)
              ! and floating conditions are not satisfied;
              ! velocity taken from the SIA solution for grounded ice

         flag_calc_vxy_ssa_y(j,i) = .false.

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_dp   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = vy_m_sia(j,i)

         lgs_x_value(nr) = vy_m_sia(j,i)

      else if ( &
                ( (mask(j,i)==3).and.(mask(j+1,i)==1) ) &
                .or. &
                ( (mask(j,i)==1).and.(mask(j+1,i)==3) ) &
              ) then
           ! one neighbour is floating ice and the other is ice-free land;
           ! velocity assumed to be zero

         flag_calc_vxy_ssa_y(j,i) = .true.

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_dp   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = 0.0_dp

         lgs_x_value(nr) = 0.0_dp

      else if ( &
                ( flag_calving_front_1(j,i).and.flag_calving_front_2(j+1,i) ) &
                .or. &
                ( flag_calving_front_2(j,i).and.flag_calving_front_1(j+1,i) ) &
              ) then
              ! one neighbour is floating ice and the other is ocean
              ! (calving front)

         flag_calc_vxy_ssa_y(j,i) = .true.

         if (flag_calving_front_1(j,i)) then
            j1 = j     ! floating ice marker
         else   ! flag_calving_front_1(j+1,i)==.true.
            j1 = j+1   ! floating ice marker
         end if

         if (.not.( flag_calving_front_2(j1-1,i) &
                    .and. &
                    flag_calving_front_2(j1+1,i) ) ) then
            ! discretization of the y-component of the BC

            nc = 2*ij2n(j1,i-1)-1
                     ! smallest nc (column counter), for vx_m(j1,i-1)
            k  = k+1
            lgs_a_value(k) = -2.0_dp*inv_dxi*vis_int_g(j1,i)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j1-1,i)
                     ! next nc (column counter), for vy_m(j1-1,i)
            k  = k+1
            lgs_a_value(k) = -4.0_dp*inv_deta*vis_int_g(j1,i)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j1,i)-1
                     ! next nc (column counter), for vx_m(j1,i)
            k  = k+1
            lgs_a_value(k) = 2.0_dp*inv_dxi*vis_int_g(j1,i)
            lgs_a_index(k) = nc

            nc = 2*ij2n(j1,i)
                     ! largest nc (column counter), for vy_m(j1,i)
            k  = k+1
            lgs_a_value(k) = 4.0_dp*inv_deta*vis_int_g(j1,i)
            lgs_a_index(k) = nc

            lgs_b_value(nr) = factor_rhs_2 &
                                 *H(j1,i)**2

            lgs_x_value(nr) = vy_m_ssa(j,i)

         else   !      (flag_calving_front_2(j1-1,i)==.true.)
                ! .and.(flag_calving_front_2(j1+1,i)==.true.);
                ! velocity assumed to be zero

            k  = k+1
            lgs_a_value(k) = 1.0_dp   ! diagonal element only
            lgs_a_index(k) = nr

            lgs_b_value(nr) = 0.0_dp

            lgs_x_value(nr) = 0.0_dp

         end if

      else if ( (mask(j,i)==0).or.(mask(j+1,i)==0) ) then
           ! neither neighbour is floating ice, but at least one neighbour is
           ! grounded ice; velocity taken from the SIA solution for grounded ice

         flag_calc_vxy_ssa_y(j,i) = .false.

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_dp   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = vy_m_sia(j,i)

         lgs_x_value(nr) = vy_m_sia(j,i)

      else   ! neither neighbour is floating or grounded ice,
             ! velocity assumed to be zero

         flag_calc_vxy_ssa_y(j,i) = .false.

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_dp   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = 0.0_dp

         lgs_x_value(nr) = 0.0_dp

      end if

   else   ! boundary condition, velocity assumed to be zero

      flag_calc_vxy_ssa_y(j,i) = .false.

      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
      lgs_a_value(k) = 1.0_dp   ! diagonal element only
      lgs_a_index(k) = nr

      lgs_b_value(nr) = 0.0_dp

      lgs_x_value(nr) = 0.0_dp

   end if

   lgs_a_ptr(nr+1) = k+1   ! row is completed, store index to next row

end do

!-------- Settings for Lis --------

#if !defined(ALLOW_TAPENADE) /* Normal */

call lis_matrix_create(LIS_COMM_WORLD, lgs_a, ierr)
call lis_vector_create(LIS_COMM_WORLD, lgs_b, ierr)
call lis_vector_create(LIS_COMM_WORLD, lgs_x, ierr)

call lis_matrix_set_size(lgs_a, 0, nmax, ierr)
call lis_vector_set_size(lgs_b, 0, nmax, ierr)
call lis_vector_set_size(lgs_x, 0, nmax, ierr)

do nr=1, nmax

   do nc=lgs_a_ptr(nr), lgs_a_ptr(nr+1)-1
      call lis_matrix_set_value(LIS_INS_VALUE, nr, lgs_a_index(nc), &
                                               lgs_a_value(nc), lgs_a, ierr)
   end do

   call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_b_value(nr), lgs_b, ierr)
   call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_x_value(nr), lgs_x, ierr)

end do

call lis_matrix_set_type(lgs_a, LIS_MATRIX_CSR, ierr)
call lis_matrix_assemble(lgs_a, ierr)

!-------- Solution of the system of linear equations with Lis --------

call lis_solver_create(solver, ierr)

#if (defined(LIS_OPTS))
    ch_solver_set_option = trim(LIS_OPTS)
#else
    ch_solver_set_option = '-i bicgsafe -p jacobi '// &
                           '-maxiter 100 -tol 1.0e-4 -initx_zeros false'
#endif

call lis_solver_set_option(trim(ch_solver_set_option), solver, ierr)
call CHKERR(ierr)

call lis_solve(lgs_a, lgs_b, lgs_x, solver, ierr)
call CHKERR(ierr)

call lis_solver_get_iter(solver, lin_iter, ierr)

write(6,'(a,i0,a)', advance='no') 'lin_iter = ', lin_iter, ', '

!%% call lis_solver_get_time(solver,solver_time,ierr)
!%% print *, 'calc_vxy_ssa_matrix: time (s) = ', solver_time

lgs_x_value = 0.0_dp
call lis_vector_gather(lgs_x, lgs_x_value, ierr)
call lis_matrix_destroy(lgs_a, ierr)
call lis_vector_destroy(lgs_b, ierr)
call lis_vector_destroy(lgs_x, ierr)
call lis_solver_destroy(solver, ierr)

#else /* Tapenade */

lgs_a_index_pass = int(lgs_a_index)
call sico_lis_solver(nmax, n_sprs, &
                           lgs_a_ptr, lgs_a_index_pass, &
                           lgs_a_value, lgs_b_value, lgs_x_value)

#endif /* Normal vs. Tapenade */

do n=1, nmax-1, 2

   i = n2i((n+1)/2)
   j = n2j((n+1)/2)

   nr = n
   vx_m_ssa(j,i) = lgs_x_value(nr)

   nr = n+1
   vy_m_ssa(j,i) = lgs_x_value(nr)

end do

#if !defined(ALLOW_TAPENADE) /* Normal */
deallocate(lgs_a_value, lgs_a_index, lgs_a_ptr)
deallocate(lgs_b_value, lgs_x_value)
#endif /* Normal */

#else

errormsg = ' >>> calc_vxy_ssa_matrix: Only to be called for' &
         //                           end_of_line &
         //'                          MARGIN==3, DYNAMICS==2 or DYNAMICS==3!'
call error(errormsg)

#endif

end subroutine calc_vxy_ssa_matrix

!-------------------------------------------------------------------------------
!> Computation of the depth-integrated viscosity vis_int_g in the SSA/SStA.
!-------------------------------------------------------------------------------
subroutine calc_vis_ssa(dxi, deta, dzeta_c, dzeta_t)

use ice_material_properties_m, only : viscosity

implicit none

real(dp), intent(in) :: dxi, deta, dzeta_c, dzeta_t

integer(i4b) :: i, j, kc, kt, m
integer(i4b) :: m_smooth, m_smooth_abs
real(dp) :: visc_min, visc_max, visc_init
real(dp) :: dxi_inv, deta_inv
real(dp) :: dvx_dxi, dvx_deta, dvy_dxi, dvy_deta
real(dp) :: enh_val
real(dp) :: diff_vis
real(dp) :: aqxy1(0:KCMAX)
real(dp) :: cvis0(0:KTMAX), cvis1(0:KCMAX)
real(dp), dimension(0:JMAX,0:IMAX) :: vis_ave_g_smooth

! ---- variables created for diva :

real(dp) :: dvx_dzeta_c, dvy_dzeta_c, dvx_dzeta_t, dvy_dzeta_t
real(dp) :: de_ssa_squared, H_t_ratio, inv_H
real(dp) :: H_dzeta_c_dz(0:KCMAX), H_dz_c_dzeta(0:KCMAX)
real(dp) :: de_tmp
real(dp) :: tau_by_g, tau_bx_g

#if (MARGIN==3 || DYNAMICS==2 || DYNAMICS==3)

!-------- Parameters, term abbreviations --------

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

dxi_inv  = 1.0_dp/dxi
deta_inv = 1.0_dp/deta

!-------- Parameters for the sigma transformation computation --------

do kc=0, KCMAX
   if (flag_aa_nonzero) then
      aqxy1(kc) = aa/(ea-1.0_dp)*eaz_c(kc)*dzeta_c
      H_dzeta_c_dz(kc) = (ea-1.0_dp)/(aa * eaz_c(kc)) ! dzeta_c_dz(kc) times H_c(j,i), for optimisation

      H_dz_c_dzeta(kc) = (aa * eaz_c(kc))/(ea-1.0_dp) ! dzeta_c_dz(kc) divided by H_c(j,i)
   else
      aqxy1(kc) = dzeta_c
      H_dzeta_c_dz(kc) = 1.0_dp ! for linear sigma transformation, dzeta_c_dz=H_c(j,i)

      H_dz_c_dzeta(kc) = 1.0_dp ! for linear sigma transformation, dz_c_dzeta= 1 * H_c(j,i)
   end if
end do

!-------- Computation of the depth-integrated viscosity --------

do i=0, IMAX
do j=0, JMAX

   if ((mask(j,i)==0).and.(.not.flag_shelfy_stream(j,i))) then
                                                   ! grounded ice, but
                                                   ! not shelfy stream
      de_ssa(j,i) = 0.0_dp   ! dummy value

#if (DYNAMICS==3)   /* DIVA */
      de_c_diva(:,j,i) = 0.0_dp   ! dummy values
      de_t_diva(:,j,i) = 0.0_dp   ! dummy values
#endif

      vis_ave_g(j,i) = 1.0_dp/flui_ave_sia(j,i)
      vis_int_g(j,i) = H(j,i) * vis_ave_g(j,i)

   else if ((mask(j,i)==1).or.(mask(j,i)==2)) then
                                                   ! ice-free land or ocean
      de_ssa(j,i) = 0.0_dp   ! dummy value

#if (DYNAMICS==3)   /* DIVA */
      de_c_diva(:,j,i) = 0.0_dp   ! dummy values
      de_t_diva(:,j,i) = 0.0_dp   ! dummy values
#endif

      vis_ave_g(j,i) = visc_init   ! dummy value
      vis_int_g(j,i) = 0.0_dp      ! dummy value

   else   ! (mask(j,i)==3).or.(flag_shelfy_stream(j,i)),
          ! floating ice or shelfy stream; 
          ! must not be at the margin of the computational domain

!  ------ Effective strain rate

      dvx_dxi  = (vx_m_ssa(j,i)-vx_m_ssa(j,i-1))*dxi_inv 
      dvy_deta = (vy_m_ssa(j,i)-vy_m_ssa(j-1,i))*deta_inv

      dvx_deta = 0.25_dp*deta_inv &
                    *( (vx_m_ssa(j+1,i)+vx_m_ssa(j+1,i-1)) &
                      -(vx_m_ssa(j-1,i)+vx_m_ssa(j-1,i-1)) )
      dvy_dxi  = 0.25_dp*dxi_inv &
                    *( (vy_m_ssa(j,i+1)+vy_m_ssa(j-1,i+1)) &
                      -(vy_m_ssa(j,i-1)+vy_m_ssa(j-1,i-1)) )

#if !defined(ALLOW_TAPENADE) /* Normal */

      de_ssa(j,i) = sqrt( dvx_dxi*dvx_dxi &
                        + dvy_deta*dvy_deta &
                        + dvx_dxi*dvy_deta &
                        + 0.25_dp*(dvx_deta+dvy_dxi)*(dvx_deta+dvy_dxi) )

#else /* Tapenade: guarding against non-differentiable sqrt(0) */

      if ( ( dvx_dxi*dvx_dxi &
           + dvy_deta*dvy_deta &
           + dvx_dxi*dvy_deta &
           + 0.25_dp*(dvx_deta+dvy_dxi)*(dvx_deta+dvy_dxi) ) > 0 ) then
         de_ssa(j,i) = sqrt( dvx_dxi*dvx_dxi &
                           + dvy_deta*dvy_deta &
                           + dvx_dxi*dvy_deta &
                           + 0.25_dp*(dvx_deta+dvy_dxi)*(dvx_deta+dvy_dxi) )
      else
         de_ssa(j,i) = 0.0_dp 
      end if

#endif /* Normal vs. Tapenade */

#if (DYNAMICS==3)   /* DIVA */ 

      if (flag_shelfy_stream(j,i)) then   ! (flag name is missleading, just convenient to re-use) ! Compute de_c_diva and de_t_diva here

         !needed 2D values :
         de_ssa_squared = de_ssa(j,i)*de_ssa(j,i)
         inv_H = 1.0_dp/H(j,i)
         tau_bx_g = 0.5_dp * (tau_bx(j,i) + tau_bx(j,i-1))
         tau_by_g = 0.5_dp * (tau_by(j,i) + tau_by(j-1,i))

         ! compute de_c_diva : 

         do kc=0, KCMAX ! from eq 36 Lipscomb-2019
            dvx_dzeta_c = flui_c_diva(kc,j,i) * inv_H * tau_bx_g * (H_c(j,i)) * (1.0_dp - eaz_c_quotient(kc))
            dvy_dzeta_c = flui_c_diva(kc,j,i) * inv_H * tau_by_g * (H_c(j,i)) * (1.0_dp - eaz_c_quotient(kc))

            de_c_diva(kc,j,i) = sqrt(de_ssa_squared + 0.25_dp*(dvx_dzeta_c*dvx_dzeta_c) + 0.25_dp*(dvy_dzeta_c*dvy_dzeta_c))
         end do   

         ! compute de_t_diva : 
         
         if (n_cts(j,i) == 1) then !ensure that there is a physically existant temperate base
            H_t_ratio=H_t(j,i)*inv_H

            do kt=0, KTMAX ! from eq 36 Lipscomb-2019
               dvx_dzeta_t = flui_t_diva(kt,j,i) * tau_bx_g  * ( 1.0_dp - zeta_t(kt) * H_t_ratio )
               dvy_dzeta_t = flui_t_diva(kt,j,i) * tau_by_g  * ( 1.0_dp - zeta_t(kt) * H_t_ratio )

               de_t_diva(kt,j,i) = sqrt(de_ssa_squared + 0.25_dp*(dvx_dzeta_t*dvx_dzeta_t) + 0.25_dp*(dvy_dzeta_t*dvy_dzeta_t))
            end do

         else ! no physical temperate layer
            do kt=0, KTMAX
               de_t_diva(kt,j,i)= de_c_diva(0,j,i) !temperate "slice" collapsed to the bottom
            end do
         endif
   

      else if (mask(j,i)==3) then   ! floating ice

         de_c_diva(:,j,i) = de_ssa(j,i)
         de_t_diva(:,j,i) = de_ssa(j,i)

      else
         errormsg = ' >>> calc_vis_ssa:' &
         //                    end_of_line &
         //'                   Either mask(j,i)==3' &
         //                    end_of_line &
         //'                   or flag_shelfy_stream(j,i) must be true here!'
         call error(errormsg)
      end if

#endif

!  ------ Term abbreviations

#if (DYNAMICS==2 || DYNAMICS==3)
      if (.not.flag_shelfy_stream(j,i)) then ! for DYN = 3, flag_shelfy_stream is True everywhere on grounded ice
#endif

         do kc=0, KCMAX

            enh_val = enh_c(kc,j,i)

            cvis1(kc) = aqxy1(kc)*H_c(j,i) &
                           *viscosity(de_ssa(j,i), &
                              temp_c(kc,j,i), temp_c_m(kc,j,i), &
                              0.0_dp, enh_val, 0)

#if (DYNAMICS==3)
            flui_c_diva(kc,j,i)=1.0_dp/min(max(viscosity(de_ssa(j,i), &
                              temp_c(kc,j,i), temp_c_m(kc,j,i), &
                              0.0_dp, enh_val, 0),visc_min),visc_max)
#endif
         end do
         ! Ice shelves (floating ice) are assumed to consist of cold ice only

#if (DYNAMICS==2 || DYNAMICS==3)
      else   ! flag_shelfy_stream(j,i) == .true.

#if (DYNAMICS==2)
         de_tmp=de_ssa(j,i)
#endif

#if (CALCMOD==-1 || CALCMOD==0)

         do kc=0, KCMAX
#if (DYNAMICS==3)
            de_tmp=de_c_diva(kc,j,i)
#endif

            if (flag_enh_stream) then
               enh_val = enh_stream
            else
               enh_val = enh_c(kc,j,i)
            end if

            cvis1(kc) = aqxy1(kc)*H_c(j,i) &
                           *viscosity(de_tmp, &
                              temp_c(kc,j,i), temp_c_m(kc,j,i), &
                              0.0_dp, enh_val, 0)

#if (DYNAMICS==3)
            flui_c_diva(kc,j,i)=1.0_dp/min(max(viscosity(de_tmp, &
                             temp_c(kc,j,i), temp_c_m(kc,j,i), &
                             0.0_dp, enh_val, 0),visc_min),visc_max)
#endif
         end do

#elif (CALCMOD==1)

         do kt=0, KTMAX
#if (DYNAMICS==3)
            de_tmp=de_t_diva(kt,j,i)
#endif

            if (flag_enh_stream) then
               enh_val = enh_stream
            else
               enh_val = enh_t(kt,j,i)
            end if

            cvis0(kt) = dzeta_t*H_t(j,i) &
                           *viscosity(de_tmp, &
                              temp_t_m(kt,j,i), temp_t_m(kt,j,i), &
                              omega_t(kt,j,i), enh_val, 1)

#if (DYNAMICS==3)
            flui_t_diva(kt,j,i)=1.0_dp/min(max(viscosity(de_tmp, &
                                 temp_t_m(kt,j,i), temp_t_m(kt,j,i), &
                                 omega_t(kt,j,i), enh_val, 1),visc_min),visc_max)
#endif
         end do

         do kc=0, KCMAX
#if (DYNAMICS==3)
            de_tmp=de_c_diva(kc,j,i)
#endif

            if (flag_enh_stream) then
               enh_val = enh_stream
            else
               enh_val = enh_c(kc,j,i)
            end if

            cvis1(kc) = aqxy1(kc)*H_c(j,i) &
                           *viscosity(de_tmp, &
                              temp_c(kc,j,i), temp_c_m(kc,j,i), &
                              0.0_dp, enh_val, 0)

#if (DYNAMICS==3)
            flui_c_diva(kc,j,i)=1.0_dp/min(max(viscosity(de_tmp, &
                                 temp_c(kc,j,i), temp_c_m(kc,j,i), &
                                 0.0_dp, enh_val, 0),visc_min),visc_max)
#endif
                              
         end do

#elif (CALCMOD==2 || CALCMOD==3)

         do kc=0, KCMAX
#if (DYNAMICS==3)
            de_tmp=de_c_diva(kc,j,i)
#endif

            if (flag_enh_stream) then
               enh_val = enh_stream
            else
               enh_val = enh_c(kc,j,i)
            end if

            cvis1(kc) = aqxy1(kc)*H_c(j,i) &
                         *viscosity(de_tmp, &
                           temp_c(kc,j,i), temp_c_m(kc,j,i), &
                           omega_c(kc,j,i), enh_val, 2)

#if (DYNAMICS==3)
            flui_c_diva(kc,j,i)=1.0_dp/min(max(viscosity(de_tmp, &
                              temp_c(kc,j,i), temp_c_m(kc,j,i), &
                              omega_c(kc,j,i), enh_val, 2),visc_min),visc_max)
#endif
         end do

#else
         errormsg = ' >>> calc_vis_ssa: CALCMOD must be -1, 0, 1, 2 or 3!'
         call error(errormsg)
#endif

      end if

#endif   /* DYNAMICS==2 or 3 */

!  ------ Depth-integrated viscosity

      vis_int_g(j,i) = 0.0_dp

#if (CALCMOD==1)
      do kt=0, KTMAX-1
         vis_int_g(j,i) = vis_int_g(j,i)+0.5_dp*(cvis0(kt+1)+cvis0(kt))
      end do
#endif

      do kc=0, KCMAX-1
         vis_int_g(j,i) = vis_int_g(j,i)+0.5_dp*(cvis1(kc+1)+cvis1(kc))
      end do

!  ------ Depth-averaged viscosity

      vis_ave_g(j,i) = vis_int_g(j,i)/max(H(j,i), eps_dp)

      vis_ave_g(j,i) = max(min(vis_ave_g(j,i), visc_max), visc_min)

   end if
  
end do
end do

!-------- Smoothing of the depth-averaged viscosity --------

#if (defined(N_VISC_SMOOTH))
  m_smooth = N_VISC_SMOOTH
#else
  m_smooth = 0
#endif

if (m_smooth /= 0) then

   m_smooth_abs = abs(m_smooth)

#if (defined(VISC_SMOOTH_DIFF))
   diff_vis = VISC_SMOOTH_DIFF
#else
   diff_vis = 0.0_dp
#endif

   if (m_smooth < 0) vis_ave_g = log(vis_ave_g)   ! logarithmic smoothing

   vis_ave_g_smooth = vis_ave_g

   do m=1, m_smooth_abs

      do i=1, IMAX-1
      do j=1, JMAX-1
         vis_ave_g_smooth(j,i) = (1.0_dp-4.0_dp*diff_vis)*vis_ave_g(j,i) &
                                    + diff_vis &
                                       *( (vis_ave_g(j,i+1)+vis_ave_g(j,i-1)) &
                                         +(vis_ave_g(j+1,i)+vis_ave_g(j-1,i)) )
      end do
      end do

      vis_ave_g = vis_ave_g_smooth

   end do

   if (m_smooth < 0) vis_ave_g = exp(vis_ave_g)   ! logarithmic smoothing

end if

!-------- Final depth-integrated viscosity --------

vis_int_g = vis_ave_g*H

#else

errormsg = ' >>> calc_vis_ssa: Only to be called for' &
         //                    end_of_line &
         //'                   MARGIN==3, DYNAMICS==2 or DYNAMICS==3!'
call error(errormsg)

#endif

end subroutine calc_vis_ssa

#if (DYNAMICS==3)

!-------------------------------------------------------------------------------
!> For DYNAMICS==3 (DIVA), subroutine to compute F_2 and F_1 from Lispcomb 2019,
!! Mathematical integral used to compute beta_drag and the 3D velocities,
!! F_2 is 2D, F_1 are 3D, on the main grid.
!-------------------------------------------------------------------------------
subroutine calc_F_int_DIVA( dzeta_c, dzeta_t)

implicit none

real(dp), intent(in) :: dzeta_c, dzeta_t

integer(i4b) ::i, j, kc, kt
real(dp) :: H_dz_c_dzeta(0:KCMAX)
real(dp), dimension(0:KCMAX) :: f_2_pre_int_c, f_1_pre_int_c
real(dp), dimension(0:KTMAX) :: f_2_pre_int_t, f_1_pre_int_t
real(dp) :: H_c_ratio
real(dp) :: H_t_ratio
real(dp) :: enh_val
real(dp) :: inv_H
real(dp) :: var_mid



!-------- Parameters for the sigma transformation computation --------

do kc=0, KCMAX
   if (flag_aa_nonzero) then
      H_dz_c_dzeta(kc) = (aa * eaz_c(kc))/(ea-1.0_dp) ! dzeta_c_dz(kc) divided by H_c(j,i)
   else
      H_dz_c_dzeta(kc) = 1.0_dp ! for linear sigma transformation, dz_c_dzeta= 1 * H_c(j,i)
   end if
end do

do i=0, IMAX
do j=0, JMAX

   if (flag_shelfy_stream(j,i)) then

      ! initialisation to ensure values are always defined
      f_2_pre_int_t(:) = 0.0_dp
      f_1_pre_int_t(:) = 0.0_dp

      f_2_pre_int_c(:) = 0.0_dp
      f_1_pre_int_c(:) = 0.0_dp

      H_t_ratio = 0.0_dp

      inv_H = 1.0_dp/H(j,i)
      H_c_ratio = H_c(j,i) * inv_H
      if ( n_cts(j,i)==1 ) then
         H_t_ratio = H_t(j,i) * inv_H
      end if


      if ( n_cts(j,i)==1 ) then

         do kt=0, KTMAX

            f_1_pre_int_t(kt) = flui_t_diva(kt,j,i) &
                                 * (1.0_dp - H_t_ratio * zeta_t(kt)) &
                                 * H_t(j,i) ! * sigma transformation

            f_2_pre_int_t(kt) = f_1_pre_int_t(kt) &
                                 * (1.0_dp - H_t_ratio * zeta_t(kt))

         end do

      end if


      do kc=0, KCMAX

         f_1_pre_int_c(kc) = flui_c_diva(kc,j,i) &
                              * ( H_c_ratio * (1.0_dp - eaz_c_quotient(kc)) ) &
                              * H_dz_c_dzeta(kc) * H_c(j,i)! * sigma transformation

         f_2_pre_int_c(kc) = f_1_pre_int_c(kc) &
                              * ( H_c_ratio * (1.0_dp - eaz_c_quotient(kc)) )

      end do

      !-------- Integration --------
      F_2_g(j,i) = 0.0_dp
      F_1_t_g(0,j,i) = 0.0_dp

      do kt=1, KTMAX
         var_mid = 0.5_dp * (f_2_pre_int_t(kt) + f_2_pre_int_t(kt-1))
         F_2_g(j,i)     = F_2_g(j,i) + var_mid * dzeta_t

         var_mid = 0.5_dp * (f_1_pre_int_t(kt) + f_1_pre_int_t(kt-1))
         F_1_t_g(kt,j,i) = F_1_t_g(kt-1,j,i) + var_mid * dzeta_t
      end do

      F_1_c_g(0,j,i) = F_1_t_g(KTMAX,j,i)

      do kc=1, KCMAX
         var_mid = 0.5_dp * (f_2_pre_int_c(kc) + f_2_pre_int_c(kc-1))
         F_2_g(j,i)     = F_2_g(j,i) + var_mid * dzeta_c

         var_mid = 0.5_dp * (f_1_pre_int_c(kc) + f_1_pre_int_c(kc-1))
         F_1_c_g(kc,j,i) = F_1_c_g(kc-1,j,i) + var_mid * dzeta_c
      end do

      if (F_2_g(j,i) < 0.0_dp) then
         errormsg = ' >>> calc_F_int_DIVA: F_2_g is negative!'
         call error(errormsg)
      end if

      !-------- Stagger on x and y grid ---------
      !if (flag_shelfy_stream(j,i) .and. flag_shelfy_stream(j,i+1)) then
      F_2_x(j,i) = 0.5_dp*(F_2_g(j,i) + F_2_g(j,i+1))
      F_2_y(j,i) = 0.5_dp*(F_2_g(j,i) + F_2_g(j+1,i))

   else !not shelfy stream :
      F_2_g(j,i)     = 0.0_dp
      F_2_x(j,i)     = 0.0_dp
      F_2_y(j,i)     = 0.0_dp
      F_1_c_g(:,j,i) = 0.0_dp
   end if
end do
end do

end subroutine calc_F_int_DIVA

#endif   /* DYNAMICS==3 */

!-------------------------------------------------------------------------------
!> Gradual limitation of computed horizontal velocities to the interval
!! [-vel_max, vel_max].
!-------------------------------------------------------------------------------
subroutine velocity_limiter_gradual(velocity, vel_max, vel_max_inv)

implicit none

real(dp), intent(in)    :: vel_max, vel_max_inv
real(dp), intent(inout) :: velocity

real(dp) :: vel_abs, vel_sign, vel_scaled, vel_scaled_lim

vel_abs = abs(velocity)

if (vel_abs >= 1.1_dp*vel_max) then

   vel_sign = sign(1.0_dp, velocity)
   velocity = vel_sign * vel_max

else if (vel_abs > 0.9_dp*vel_max) then

   ! gradual limitation between 0.9*vel_max and 1.1*vel_max

   vel_sign = sign(1.0_dp, velocity)

   vel_scaled     = (vel_abs-0.9_dp*vel_max)*(10.0_dp*vel_max_inv)
                       ! between 0 and 2
   vel_scaled_lim = vel_scaled &
                       *(1.0_dp-0.25_dp*vel_scaled*vel_scaled &
                                       *(1.0_dp-0.25_dp*vel_scaled))
                       ! between 0 and 1

   velocity = vel_sign * vel_max * (0.9_dp + 0.1_dp*vel_scaled_lim)

end if

end subroutine velocity_limiter_gradual

!-------------------------------------------------------------------------------

end module calc_vxy_m
!
