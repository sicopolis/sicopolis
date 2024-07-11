!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ v a r i a b l e s _ m
!
!! Declarations of global variables for SICOPOLIS.
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
!> Declarations of global variables for SICOPOLIS.
!-------------------------------------------------------------------------------
module sico_variables_m

use sico_types_m

implicit none
save

!-------- Field quantities --------

integer(i4b), dimension(0:JMAX,0:IMAX) :: mask
   !! Ice-land-ocean mask:
   !!  0: grounded ice,
   !!  1: ice-free land,
   !!  2: ocean,
   !!  3: floating ice

integer(i4b), dimension(0:JMAX,0:IMAX) :: mask_old
   !! Old value of mask (at the previous time step)

integer(i4b), dimension(0:JMAX,0:IMAX) :: mask_new
   !! New value of mask computed during an integration step

integer(i4b), dimension(0:JMAX,0:IMAX) :: n_cts
   !! Mask for thermal conditions:
   !!  -1: cold ice base,
   !!   0: temperate ice base with cold ice above,
   !!   1: temperate ice base with temperate ice layer above (only for POLY)

integer(i4b), dimension(0:JMAX,0:IMAX) :: n_cts_new
   !! New value of quantity computed during an integration step

integer(i4b), dimension(0:JMAX,0:IMAX) :: kc_cts
   !! Position kc of the CTS (for COLD, ENTC, ENTM)

integer(i4b), dimension(0:JMAX,0:IMAX) :: kc_cts_new
   !! New value of quantity computed during an integration step

integer(i4b), dimension(0:JMAX,0:IMAX) :: mask_region
   !! Region mask:
   !!  0: undefined,
   !!  1: EAIS,
   !!  2: WAIS,
   !!  3: AP

logical :: flag_calc_temp
   !! Flag for computation of the temperature, water content,
   !! age and flow enhancement factor during an integration step:
   !!   .true.: temperature etc. computed,
   !!  .false.: temperature etc. not computed

logical, dimension(0:JMAX,0:IMAX) :: flag_inner_point
   !! Inner-point flag:
   !!   .true.: inner point,
   !!  .false.: otherwise (margin point)

logical, dimension(0:JMAX,0:IMAX) :: flag_inner_inner_point
   !! Inner-point flag:
   !!   .true.: inner point (outermost two grid points excepted),
   !!  .false.: otherwise (margin point, or margin-neighbour point)

logical, dimension(0:JMAX,0:IMAX) :: flag_sg_x
   !! Flag for staggered grid in x-direction:
   !!   .true.: staggered-grid point in x-direction,
   !!  .false.: otherwise

logical, dimension(0:JMAX,0:IMAX) :: flag_sg_y
   !! Flag for staggered grid in y-direction:
   !!   .true.: staggered-grid point in y-direction,
   !!  .false.: otherwise

logical, dimension(0:JMAX,0:IMAX) :: flag_sg_x_inner_y
   !! Flag for staggered grid in x-direction and inner point in y-direction:
   !!   .true.: staggered-grid point in x-direction
   !!           and inner point in y-direction,
   !!  .false.: otherwise

logical, dimension(0:JMAX,0:IMAX) :: flag_sg_y_inner_x
   !! Flag for staggered grid in y-direction and inner point in x-direction:
   !!   .true.: staggered-grid point in y-direction
   !!           and inner point in x-direction,
   !!  .false.: otherwise

logical, dimension(0:JMAX,0:IMAX) :: flag_grounding_line_1
   !! Grounding line flag:
   !!   .true.: grounding line point
   !!           (grounded ice point with at least
   !!            one floating ice neighbour),
   !!  .false.: otherwise

logical, dimension(0:JMAX,0:IMAX) :: flag_grounding_line_2
   !! Grounding line flag:
   !!   .true.: grounding line point
   !!           (floating ice point with at least
   !!            one grounded ice neighbour),
   !!  .false.: otherwise

logical, dimension(0:JMAX,0:IMAX) :: flag_calving_front_1
   !! Calving front flag:
   !!   .true.: calving front point
   !!           (floating ice point with at least
   !!            one ocean neighbour),
   !!  .false.: otherwise

logical, dimension(0:JMAX,0:IMAX) :: flag_calving_front_2
   !! Calving front flag:
   !!   .true.: calving front point
   !!           (ocean point with at least
   !!            one floating ice neighbour),
   !!  .false.: otherwise

logical, dimension(0:JMAX,0:IMAX) :: flag_grounded_front_a_1
   !! Land-terminating grounded front flag:
   !!   .true.: grounded front point
   !!           (grounded ice point with at least
   !!            one ice-free land neighbour),
   !!  .false.: otherwise

logical, dimension(0:JMAX,0:IMAX) :: flag_grounded_front_a_2
   !! Land-terminating grounded front flag:
   !!   .true.: grounded front point
   !!           (ice-free land point with at least
   !!            one grounded ice neighbour),
   !!  .false.: otherwise

logical, dimension(0:JMAX,0:IMAX) :: flag_grounded_front_b_1
   !! Marine-terminating grounded front flag:
   !!   .true.: grounded front point
   !!           (grounded ice point with at least
   !!            one ocean neighbour),
   !!  .false.: otherwise

logical, dimension(0:JMAX,0:IMAX) :: flag_grounded_front_b_2
   !! Marine-terminating grounded front flag:
   !!   .true.: grounded front point
   !!           (ocean point with at least
   !!            one grounded ice neighbour),
   !!  .false.: otherwise

logical, dimension(0:JMAX,0:IMAX) :: flag_shelfy_stream_x
   !! Shelfy stream flag in x-direction, at (i+1/2,j):
   !!   .true.: shelfy stream point,
   !!  .false.: otherwise

logical, dimension(0:JMAX,0:IMAX) :: flag_shelfy_stream_y
   !! Shelfy stream flag in y-direction, at (i,j+1/2):
   !!   .true.: shelfy stream point,
   !!  .false.: otherwise

logical, dimension(0:JMAX,0:IMAX) :: flag_shelfy_stream
   !! Shelfy stream flag on the main grid:
   !!   .true.: grounded ice,
   !!           and at least one neighbour on the
   !!           staggered grid is a shelfy stream point,
   !!  .false.: otherwise

real(dp), dimension(0:IMAX) :: xi
   !! Coordinate xi (= x) of grid point i

real(dp), dimension(0:JMAX) :: eta
   !! Coordinate eta (= y) of grid point j

real(dp), dimension(0:KCMAX) :: zeta_c
   !! Sigma coordinate zeta_c of grid point kc

real(dp), dimension(0:KTMAX) :: zeta_t
   !! Sigma coordinate zeta_t of grid point kt

real(dp), dimension(0:KRMAX) :: zeta_r
   !! Sigma coordinate zeta_r of grid point kr

real(dp) :: aa
   !! Exponential stretch parameter of the non-equidistant vertical grid
   !! in the upper (kc) ice domain

logical :: flag_aa_nonzero
   !! Flag for the exponential stretch parameter aa:
   !!   .true.: aa greater than zero (non-equidistant grid),
   !!  .false.: aa equal to zero (equidistant grid)

real(dp) :: ea
   !! Abbreviation for exp(aa)

real(dp), dimension(0:KCMAX) :: eaz_c
   !! Abbreviation for exp(aa*zeta(kc))

real(dp), dimension(0:KCMAX) :: eaz_c_quotient
   !! Abbreviation for (eaz_c(kc)-1.0)/(ea-1.0)

real(dp), dimension(0:JMAX,0:IMAX) :: lambda
   !! Geographic longitude of grid point (i,j)

real(dp), dimension(0:JMAX,0:IMAX) :: phi
   !! Geographic latitude of grid point (i,j)

real(dp), dimension(0:JMAX,0:IMAX) :: cell_area
   !! Area of grid cell associated with grid point (i,j)

real(dp), dimension(0:JMAX,0:IMAX) :: sq_g11_g
   !! Square root of the coefficient g11 of the metric tensor
   !! on grid point (i,j)

real(dp), dimension(0:JMAX,0:IMAX) :: sq_g22_g
   !! Square root of the coefficient g22 of the metric tensor
   !! on grid point (i,j)

real(dp), dimension(0:JMAX,0:IMAX) :: insq_g11_g
   !! Inverse square root of g11 on grid point (i,j)

real(dp), dimension(0:JMAX,0:IMAX) :: insq_g22_g
   !! Inverse square root of g22 on grid point (i,j)

real(dp), dimension(0:JMAX,0:IMAX) :: sq_g11_sgx
   !! Square root of g11, at (i+1/2,j)

real(dp), dimension(0:JMAX,0:IMAX) :: sq_g11_sgy
   !! Square root of g11, at (i,j+1/2)

real(dp), dimension(0:JMAX,0:IMAX) :: sq_g22_sgx
   !! Square root of g22, at (i+1/2,j)

real(dp), dimension(0:JMAX,0:IMAX) :: sq_g22_sgy
   !! Square root of g22, at (i,j+1/2)

real(dp), dimension(0:JMAX,0:IMAX) :: insq_g11_sgx
   !! Inverse square root of g11, at (i+1/2,j)

real(dp), dimension(0:JMAX,0:IMAX) :: insq_g22_sgy
   !! Inverse square root of g22, at (i,j+1/2)

real(dp), dimension(0:JMAX,0:IMAX) :: zs
   !! Coordinate z of the surface topography

real(dp), dimension(0:JMAX,0:IMAX) :: zm
   !! Coordinate z of the bottom of the upper (kc) ice domain
   !!                      = top of the lower (kt) ice domain
   !! (position of the CTS for POLY,
   !!  equal to zb for ISOT, COLD, ENTC, ENTM)

real(dp), dimension(0:JMAX,0:IMAX) :: zb
   !! Coordinate z of the ice base

real(dp), dimension(0:JMAX,0:IMAX) :: zl
   !! Coordinate z of the lithosphere surface

real(dp), dimension(0:JMAX,0:IMAX) :: zl0
   !! zl for isostatically relaxed ice-free conditions

real(dp), dimension(0:JMAX,0:IMAX) :: wss
   !! Isostatic steady-state displacement of the lithosphere

real(dp), dimension(0:JMAX,0:IMAX) :: flex_rig_lith
   !! Flexural rigidity of the lithosphere

real(dp), dimension(0:JMAX,0:IMAX) :: time_lag_asth
   !! Time lag of the relaxing asthenosphere

real(dp), dimension(0:JMAX,0:IMAX) :: H
   !! Ice thickness (= H_c + H_t)

real(dp), dimension(0:JMAX,0:IMAX) :: H_c
   !! Thickness of ice in the upper (kc) domain
   !! (thickness of the cold-ice layer for POLY,
   !!  entire ice thickness for ISOT, COLD, ENTC, ENTM)

real(dp), dimension(0:JMAX,0:IMAX) :: H_t
   !! Thickness of ice in the lower (kt) domain
   !! (thickness of the temperate layer for POLY,
   !!  redundant and thus set to zero for ISOT, COLD, ENTC, ENTM)

real(dp), dimension(0:JMAX,0:IMAX) :: dzs_dxi
   !! Derivative of zs by xi (at i+1/2,j)

real(dp), dimension(0:JMAX,0:IMAX) :: dzs_dx_aux
   !! Auxiliary variable for dzs_dxi

real(dp), dimension(0:JMAX,0:IMAX) :: dzm_dxi
   !! Derivative of zm by xi (at i+1/2,j)

real(dp), dimension(0:JMAX,0:IMAX) :: dzb_dxi
   !! Derivative of zb by xi (at i+1/2,j)

real(dp), dimension(0:JMAX,0:IMAX) :: dH_c_dxi
   !! Derivative of H_c by xi (at i+1/2,j)

real(dp), dimension(0:JMAX,0:IMAX) :: dH_t_dxi
   !! Derivative of H_t by xi (at i+1/2,j)

real(dp), dimension(0:JMAX,0:IMAX) :: dzs_deta
   !! Derivative of zs by eta (at i,j+1/2)

real(dp), dimension(0:JMAX,0:IMAX) :: dzs_dy_aux
   !! Auxiliary variable for dzs_deta

real(dp), dimension(0:JMAX,0:IMAX) :: dzm_deta
   !! Derivative of zm by eta (at i,j+1/2)

real(dp), dimension(0:JMAX,0:IMAX) :: dzb_deta
   !! Derivative of zb by eta (at i,j+1/2)

real(dp), dimension(0:JMAX,0:IMAX) :: dH_c_deta
   !! Derivative of H_c by eta (at i,j+1/2)

real(dp), dimension(0:JMAX,0:IMAX) :: dH_t_deta
   !! Derivative of H_t by eta (at i,j+1/2)

real(dp), dimension(0:JMAX,0:IMAX) :: dzs_dxi_g
   !! Staggered-grid quantity interpolated to grid point (i,j)

real(dp), dimension(0:JMAX,0:IMAX) :: dzm_dxi_g
   !! Staggered-grid quantity interpolated to grid point (i,j)

real(dp), dimension(0:JMAX,0:IMAX) :: dzb_dxi_g
   !! Staggered-grid quantity interpolated to grid point (i,j)

real(dp), dimension(0:JMAX,0:IMAX) :: dH_c_dxi_g
   !! Staggered-grid quantity interpolated to grid point (i,j)

real(dp), dimension(0:JMAX,0:IMAX) :: dH_t_dxi_g
   !! Staggered-grid quantity interpolated to grid point (i,j)

real(dp), dimension(0:JMAX,0:IMAX) :: dzs_deta_g
   !! Staggered-grid quantity interpolated to grid point (i,j)

real(dp), dimension(0:JMAX,0:IMAX) :: dzm_deta_g
   !! Staggered-grid quantity interpolated to grid point (i,j)

real(dp), dimension(0:JMAX,0:IMAX) :: dzb_deta_g
   !! Staggered-grid quantity interpolated to grid point (i,j)

real(dp), dimension(0:JMAX,0:IMAX) :: dH_c_deta_g
   !! Staggered-grid quantity interpolated to grid point (i,j)

real(dp), dimension(0:JMAX,0:IMAX) :: dH_t_deta_g
   !! Staggered-grid quantity interpolated to grid point (i,j)

real(dp), dimension(0:JMAX,0:IMAX) :: dzs_dtau
   !! Derivative of zs by tau (time)

real(dp), dimension(0:JMAX,0:IMAX) :: dzm_dtau
   !! Derivative of zm by tau (time)

real(dp), dimension(0:JMAX,0:IMAX) :: dzb_dtau
   !! Derivative of zb by tau (time)

real(dp), dimension(0:JMAX,0:IMAX) :: dzl_dtau
   !! Derivative of zl by tau (time)

real(dp), dimension(0:JMAX,0:IMAX) :: dH_dtau
   !! Derivative of H by tau (time)

real(dp), dimension(0:JMAX,0:IMAX) :: dH_c_dtau
   !! Derivative of H_c by tau (time)

real(dp), dimension(0:JMAX,0:IMAX) :: dH_t_dtau
   !! Derivative of H_t by tau (time)

integer(i4b), dimension(0:JMAX,0:IMAX) :: n_slide_region
   !! Regions for the basal sliding laws

integer(i4b), dimension(0:JMAX,0:IMAX) :: n_bm_region
   !! Regions for ice shelf basal melting

integer(i4b), dimension(0:JMAX,0:IMAX) :: p_weert
   !! Weertman exponent for the basal shear stress

integer(i4b), dimension(0:JMAX,0:IMAX) :: q_weert
   !! Weertman exponent for the basal pressure

real(dp), dimension(0:JMAX,0:IMAX) :: p_weert_inv
   !! Inverse of p_weert

real(dp), dimension(0:JMAX,0:IMAX) :: c_slide
   !! Basal sliding coefficient

real(dp), dimension(0:JMAX,0:IMAX) :: c_slide_init
   !! Initial basal sliding coefficient

real(dp), dimension(0:JMAX,0:IMAX) :: gamma_slide_inv
   !! Inverse of the sub-melt-sliding parameter

logical, dimension(0:JMAX,0:IMAX) :: sub_melt_flag
   !! Flag for presence of sub-melt sliding

real(dp), dimension(0:JMAX,0:IMAX) :: d_help_b
   !! Auxiliary quantity for the computation of vx_b and vy_b

real(dp), dimension(0:JMAX,0:IMAX) :: c_drag
   !! Auxiliary quantity for the computation of the basal drag

real(dp), dimension(0:JMAX,0:IMAX) :: p_b
   !! Basal pressure

real(dp), dimension(0:JMAX,0:IMAX) :: p_b_w
   !! Basal water pressure

real(dp), dimension(0:JMAX,0:IMAX) :: p_b_red
   !! Reduced basal pressure

real(dp), dimension(0:JMAX,0:IMAX) :: tau_dr
   !! Driving stress

real(dp), dimension(0:JMAX,0:IMAX) :: tau_b
   !! Basal shear stress (drag)

real(dp), dimension(0:JMAX,0:IMAX) :: vx_b
   !! Velocity in x-direction at the ice base, at (i+1/2,j)

real(dp), dimension(0:JMAX,0:IMAX) :: vy_b
   !! Velocity in y-direction at the ice base, at (i,j+1/2)

real(dp), dimension(0:JMAX,0:IMAX) :: vx_m
   !! Mean (depth-averaged) velocity in x-direction, at (i+1/2,j)

real(dp), dimension(0:JMAX,0:IMAX) :: vy_m
   !! Mean (depth-averaged) velocity in y-direction, at (i,j+1/2)

real(dp), dimension(0:JMAX,0:IMAX) :: vx_m_sia
   !! Mean (depth-averaged) SIA velocity in x-direction, at (i+1/2,j)

real(dp), dimension(0:JMAX,0:IMAX) :: vy_m_sia
   !! Mean (depth-averaged) SIA velocity in y-direction, at (i,j+1/2)

real(dp), dimension(0:JMAX,0:IMAX) :: vx_m_ssa
   !! Mean (depth-averaged) SSA velocity in x-direction, at (i+1/2,j)

real(dp), dimension(0:JMAX,0:IMAX) :: vy_m_ssa
   !! Mean (depth-averaged) SSA velocity in y-direction, at (i,j+1/2)

real(dp), dimension(0:JMAX,0:IMAX) :: ratio_sl_x
   !! Ratio of basal to surface velocity (slip ratio) in x-direction,
   !! at (i+1/2,j)

real(dp), dimension(0:JMAX,0:IMAX) :: ratio_sl_y
   !! Ratio of basal to surface velocity (slip ratio) in y-direction,
   !! at (i,j+1/2)

real(dp), dimension(0:JMAX,0:IMAX) :: ratio_sl
   !! Ratio of basal to surface velocity (slip ratio) on the main grid

real(dp), dimension(0:JMAX,0:IMAX) :: vx_b_g
   !! Staggered-grid quantity interpolated to grid point (i,j)

real(dp), dimension(0:JMAX,0:IMAX) :: vy_b_g
   !! Staggered-grid quantity interpolated to grid point (i,j)

real(dp), dimension(0:JMAX,0:IMAX) :: vz_b
   !! Velocity in z-direction at the ice base

real(dp), dimension(0:JMAX,0:IMAX) :: vz_m
   !! Velocity in z-direction at the position z=zm (interface between
   !! the upper (kc) and the lower (kt) domain)

real(dp), dimension(0:JMAX,0:IMAX) :: vx_s_g
   !! Velocity in x-direction at the ice surface

real(dp), dimension(0:JMAX,0:IMAX) :: vy_s_g
   !! Velocity in x-direction at the ice surface

real(dp), dimension(0:JMAX,0:IMAX) :: vz_s
   !! Velocity in z-direction at the ice surface

real(dp), dimension(0:JMAX,0:IMAX) :: flui_ave_sia
   !! Depth-averaged fluidity of the SIA

real(dp), dimension(0:JMAX,0:IMAX) :: h_diff
   !! Diffusivity of the SIA evolution equation of the ice surface

real(dp), dimension(0:JMAX,0:IMAX) :: qx
   !! Volume flux in x-direction (depth-integrated vx, at (i+1/2,j))

real(dp), dimension(0:JMAX,0:IMAX) :: qy
   !! Volume flux in y-direction (depth-integrated vy, at (i,j+1/2))

real(dp), dimension(0:JMAX,0:IMAX) :: q_gl_g
   !! Volume flux across the grounding line

real(dp), dimension(0:JMAX,0:IMAX) :: q_geo
   !! Geothermal heat flux

real(dp), dimension(0:JMAX,0:IMAX) :: temp_b
   !! Basal temperature

real(dp), dimension(0:JMAX,0:IMAX) :: temph_b
   !! Basal temperature relative to the pressure melting point

real(dp), dimension(0:JMAX,0:IMAX) :: Q_bm
   !! Basal melting rate

real(dp), dimension(0:JMAX,0:IMAX) :: Q_tld
   !! Water drainage rate from the temperate layer

real(dp), dimension(0:JMAX,0:IMAX) :: Q_b_tot
   !! Sum of Q_bm and Q_tld

real(dp), dimension(0:JMAX,0:IMAX) :: q_w
   !! Scalar volume flux of the basal water

real(dp), dimension(0:JMAX,0:IMAX) :: q_w_x
   !! Scalar volume flux of the basal water in x-direction

real(dp), dimension(0:JMAX,0:IMAX) :: q_w_y
   !! Scalar volume flux of the basal water in y-direction

real(dp), dimension(0:JMAX,0:IMAX) :: H_w
   !! Thickness of the water column under the ice base

real(dp), dimension(0:JMAX,0:IMAX) :: accum
   !! Accumulation rate at the ice surface
   !! (including liquid precipitation = rainfall)

real(dp), dimension(0:JMAX,0:IMAX) :: snowfall
   !! Snowfall rate at the ice surface

real(dp), dimension(0:JMAX,0:IMAX) :: rainfall
   !! Rainfall rate at the ice surface

real(dp), dimension(0:JMAX,0:IMAX) :: ET
   !! Temperature excess at the ice surface
   !! (positive degree days divided by time)

real(dp), dimension(0:JMAX,0:IMAX) :: melt
   !! Melting rate at the ice surface

real(dp), dimension(0:JMAX,0:IMAX) :: melt_star
   !! Superimposed ice formation rate at the ice surface

real(dp), dimension(0:JMAX,0:IMAX) :: runoff
   !! Runoff rate at the ice surface

real(dp), dimension(0:JMAX,0:IMAX) :: as_perp
   !! Accumulation-ablation function at the ice surface (SMB)

real(dp), dimension(0:JMAX,0:IMAX) :: as_perp_apl
   !! Applied accumulation-ablation function (SMB)

real(dp), dimension(0:JMAX,0:IMAX) :: smb_corr_in
   !! Prescribed SMB correction read from file

real(dp), dimension(0:JMAX,0:IMAX) :: smb_corr_prescribed
   !! Prescribed SMB correction

real(dp), dimension(0:JMAX,0:IMAX) :: smb_corr
   !! Diagnosed SMB correction

real(dp), dimension(0:JMAX,0:IMAX) :: mb_source_apl
   !! Applied mass balance source (SMB, BMB, calving)

real(dp), dimension(0:JMAX,0:IMAX) :: accum_apl
   !! Applied accumulation rate at the ice surface
   !!                 (including liquid precipitation = rainfall!)

real(dp), dimension(0:JMAX,0:IMAX) :: runoff_apl
   !! Applied runoff rate at the ice surface

real(dp), dimension(0:JMAX,0:IMAX) :: Q_b_apl
   !! Applied basal melting rate
   !!         + water drainage rate from the temperate layer

real(dp), dimension(0:JMAX,0:IMAX) :: calving
   !! Calving rate of grounded and floating ice

real(dp), dimension(0:JMAX,0:IMAX) :: calving_apl
   !! Applied calving rate of grounded and floating ice

integer(i4b), dimension(0:JMAX,0:IMAX) :: mask_ablation_type
   !! Mask indicating ablation type:
   !!   2: visible (ocean, for later developments),
   !!   1: visible (grounded ice),
   !!  -1: hidden on land,
   !!  -2: hidden in ocean

real(dp), dimension(0:JMAX,0:IMAX) :: temp_maat
   !! Mean annual air temperature

real(dp), dimension(0:JMAX,0:IMAX) :: temp_s
   !! Ice surface temperature

real(dp), dimension(0:JMAX,0:IMAX) :: z_sl
   !! Sea level

real(dp), dimension(0:JMAX,0:IMAX) :: dzsl_dtau
   !! Derivative of zsl by tau (time)

real(dp) :: z_sl_mean
   !! Mean sea level

real(dp), dimension(0:JMAX,0:IMAX) :: am_perp
   !! Ice volume flux across the z=zm interface

real(dp), dimension(0:JMAX,0:IMAX) :: am_perp_st
   !! Steady-state part of am_perp (without contribution of dzm_dtau)

real(dp), dimension(0:JMAX,0:IMAX) :: zs_new
   !! New value of quantity computed during an integration step

real(dp), dimension(0:JMAX,0:IMAX) :: zm_new
   !! New value of quantity computed during an integration step

real(dp), dimension(0:JMAX,0:IMAX) :: zb_new
   !! New value of quantity computed during an integration step

real(dp), dimension(0:JMAX,0:IMAX) :: zl_new
   !! New value of quantity computed during an integration step

real(dp), dimension(0:JMAX,0:IMAX) :: H_new
   !! New value of quantity computed during an integration step

real(dp), dimension(0:JMAX,0:IMAX) :: H_c_new
   !! New value of quantity computed during an integration step

real(dp), dimension(0:JMAX,0:IMAX) :: H_t_new
   !! New value of quantity computed during an integration step

real(dp), dimension(0:JMAX,0:IMAX) :: H_new_flow
   !! New value of the ice thickness due to glacial flow only
   !! (no source term considered)

real(dp), dimension(0:JMAX,0:IMAX) :: mb_source
   !! Source term for the ice-thickness equation

logical :: flag_thk_solver_explicit
   !! Flag for the type of solver of the ice-thickness equation:
   !!   .true.: explicit solver,
   !!  .false.: implicit solver

#if (TSURFACE<=5)
real(dp), dimension(0:JMAX,0:IMAX) :: zs_ref_temp
   !! Reference elevation for the present-day surface temperature
#endif

#if (ACCSURFACE<=5)
real(dp), dimension(0:JMAX,0:IMAX) :: zs_ref_precip
   !! Reference elevation for the present-day precipitation
#endif

#if (TSURFACE==6 && ACCSURFACE==6 && ABLSURFACE==6)
real(dp), dimension(0:JMAX,0:IMAX) :: zs_ref_climatol
   !! Reference elevation for the present-day climatology
#endif

real(dp), dimension(0:JMAX,0:IMAX) :: accum_present
   !! Present-day accumulation rate at the ice surface
   !! (for EISMINT, ISMIP HEINO and the north and south polar caps of Mars)

real(dp), dimension(0:JMAX,0:IMAX) :: precip_ma_present
   !! Present-day mean annual precipitation rate at the ice surface

real(dp), dimension(0:JMAX,0:IMAX) :: temp_ma_present
   !! Present-day mean annual surface temperature

real(dp), dimension(0:JMAX,0:IMAX) :: temp_mj_present
   !! Present-day mean summer (northern hemisphere: July,
   !! southern hemisphere: January) surface temperature

real(dp), dimension(-JMAX:JMAX,-IMAX:IMAX) :: dist_dxdy
   !! Distance between grid points with delta_i=ir, delta_j=jr

real(dp), dimension(0:JMAX,0:IMAX,12) :: precip_present
   !! Present-day mean monthly precipitation rate
   !!                        at the ice surface

real(dp), dimension(0:JMAX,0:IMAX,12) :: precip_lgm_anom
   !! LGM anomaly (ratio LGM/present) of the mean monthly precipitation rate
   !! at the ice surface

real(dp), dimension(0:JMAX,0:IMAX,12) :: gamma_precip_lgm_anom
   !! Negative natural logarithm of precip_lgm_anom(j,i,n)

real(dp), dimension(0:JMAX,0:IMAX,12) :: temp_present
   !! Present-day mean monthly surface temperature

real(dp), dimension(0:JMAX,0:IMAX,12) :: temp_lgm_anom
   !! LGM anomaly (difference LGM - present) of the mean monthly
   !! surface temperature

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: d_help_c
   !! Auxiliary quantity for the computation of vx, vy und zs

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: vx_c
   !! Velocity in x-direction in the upper (kc) ice domain
   !! (at (i+1/2,j,kc))

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: vy_c
   !! Velocity in y-direction in the upper (kc) ice domain
   !! (at (i,j+1/2,kc))

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: vz_c
   !! Velocity in z-direction in the upper (kc) ice domain
   !! (at (i,j,kc+1/2))

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: temp_c
   !! Temperature in the upper (kc) ice domain

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: temp_c_new
   !! New value of quantity computed during an integration step

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: temp_c_m
   !! Melting temperature in the upper (kc) ice domain

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: age_c
   !! Age in the upper (kc) ice domain

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: age_c_new
   !! New value of quantity computed during an integration step

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: txz_c
   !! Shear stress txz in the upper (kc) ice domain (at (i+1/2,j,kc))

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: tyz_c
   !! Shear stress tyz in the upper (kc) ice domain (at (i,j+1/2,kc))

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: sigma_c
   !! Effective stress in the upper (kc) ice domain

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: enh_c
   !! Flow enhancement factor in the upper (kc) ice domain

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: strain_heating_c
   !! Strain heating in the upper (kc) ice domain

real(dp), dimension(0:JMAX,0:IMAX) :: de_ssa
   !! Effective strain rate of the SSA

real(dp), dimension(0:JMAX,0:IMAX) :: vis_ave_g
   !! Depth-averaged viscosity of the SIA/SSA

real(dp), dimension(0:JMAX,0:IMAX) :: vis_int_g
   !! Depth-integrated viscosity of the SIA/SSA

real(dp), dimension(0:JMAX,0:IMAX) :: vx_g
   !! Velocity in x-direction of the SSA

real(dp), dimension(0:JMAX,0:IMAX) :: vy_g
   !! Velocity in y-direction of the SSA

real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: d_help_t
   !! Auxiliary quantity for the computation of vx, vy und zs

real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: vx_t
   !! Velocity in x-direction in the lower (kt) ice domain (at (i+1/2,j,kt))

real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: vy_t
   !! Velocity in y-direction in the lower (kt) ice domain (at (i,j+1/2,kt))

real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: vz_t
   !! Velocity in z-direction in the lower (kt) ice domain (at (i,j,kt+1/2))

real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: omega_t
   !! Water content in the lower (kt) ice domain

real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: omega_t_new
   !! New value of quantity computed during an integration step

real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: temp_t_m
   !! Melting temperature in the lower (kt) ice domain

real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: age_t
   !! Age in the lower (kt) ice domain

real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: age_t_new
   !! New value of quantity computed during an integration step

real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: txz_t
   !! Shear stress txz in the lower (kt) ice domain (at (i+1/2,j,kt))

real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: tyz_t
   !! Shear stress tyz in the lower (kt) ice domain (at (i,j+1/2,kt))

real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: sigma_t
   !! Effective stress in the lower (kt) ice domain

real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: enh_t
   !! Flow enhancement factor in the lower (kt) ice domain

real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: strain_heating_t
   !! Strain heating in the lower (kt) ice domain

real(dp) :: enh_stream
   !! Flow enhancement factor for ice streams

logical :: flag_enh_stream
   !! Flag for definition of flow enhancement factor for ice streams

real(dp), dimension(0:KRMAX,0:JMAX,0:IMAX) :: temp_r
   !! Temperature in the bedrock

real(dp), dimension(0:KRMAX,0:JMAX,0:IMAX) :: temp_r_new
   !! New value of quantity computed during an integration step

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: enth_c
   !! Enthalpy in the upper (kc) ice domain

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: enth_c_new
   !! New value of quantity computed during an integration step

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: omega_c
   !! Water content in the upper (kc) ice domain

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: omega_c_new
   !! New value of quantity computed during an integration step

real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: enth_t
   !! Enthalpy in the lower (kt) ice domain

real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: enth_t_new
   !! New value of quantity computed during an integration step

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: dxx_c
   !! Strain rate dxx in the upper (kc) ice domain

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: dyy_c
   !! Strain rate dyy in the upper (kc) ice domain

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: dxy_c
   !! Strain rate dxy in the upper (kc) ice domain

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: dxz_c
   !! Strain rate dxz in the upper (kc) ice domain

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: dyz_c
   !! Strain rate dyz in the upper (kc) ice domain

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: de_c
   !! Full effective strain rate in the upper (kc) ice domain

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: lambda_shear_c
   !! Shear fraction in the upper (kc) ice domain

real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: dxx_t
   !! Strain rate dxx in the lower (kt) ice domain

real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: dyy_t
   !! Strain rate dyy in the lower (kt) ice domain

real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: dxy_t
   !! Strain rate dxy in the lower (kt) ice domain

real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: dxz_t
   !! Strain rate dxz in the lower (kt) ice domain

real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: dyz_t
   !! Strain rate dyz in the lower (kt) ice domain

real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: de_t
   !! Full effective strain rate in the lower (kt) ice domain

real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: lambda_shear_t
   !! Shear fraction in the lower (kt) ice domain

!-------- Physical parameters --------

real(dp) :: RHO
   !! Density of ice

real(dp) :: RHO_W
   !! Density of pure water

real(dp) :: RHO_SW
   !! Density of sea water

real(dp) :: L
   !! Latent heat of ice

real(dp) :: G
   !! Acceleration due to gravity

real(dp) :: NUE
   !! Water diffusivity in ice

real(dp) :: BETA
   !! Clausius-Clapeyron gradient of ice

real(dp) :: DELTA_TM_SW
   !! Melting point depression of sea water due to its average salinity

real(dp) :: OMEGA_MAX
   !! Threshold value for the water content of temperate ice

real(dp) :: H_R
   !! Thickness of the modelled lithosphere layer

real(dp) :: RHO_C_R
   !! Density times specific heat of the lithosphere

real(dp) :: KAPPA_R
   !! Heat conductivity of the lithosphere

real(dp) :: RHO_A
   !! Density of the asthenosphere

real(dp) :: R_T
   !! Coefficient of the water-content dependence in the rate factor
   !! for temperate ice

real(dp) :: R
   !! Mean radius of the planet

real(dp) :: A
   !! Semi-major axis of the planet

real(dp) :: B
   !! Semi-minor axis of the planet

real(dp) :: F_INV
   !! Inverse flattening of the planet (>1e10 interpreted as infinity)

real(dp) :: LAMBDA0
   !! Reference longitude (central meridian) of the stereographic projection

real(dp) :: PHI0
   !! Standard parallel of the stereographic projection

real(dp), dimension(-190:10) :: RF
   !! Tabulated values for the rate factor of cold ice

real(dp), dimension(-190:10) :: KAPPA
   !! Tabulated values for the heat conductivity of ice

real(dp), dimension(-190:10) :: C
   !! Tabulated values for the specific heat of ice

!-------- Auxiliary variables for the module ice_material_properties_m --------

real(dp), dimension(-256:255) :: RF_imp
   !! Tabulated values for the rate factor of cold ice

real(dp) :: R_T_imp
   !! Coefficient of the water-content dependence in the rate factor
   !! for temperate ice

real(dp), dimension(-256:255) :: KAPPA_imp
   !! Tabulated values for the heat conductivity of ice

real(dp), dimension(-256:255) :: C_imp
   !! Tabulated values for the specific heat of ice

integer(i4b) :: n_temp_min_imp
   !! Lower index limit of properly defined values in
   !! RF_imp, KAPPA_imp and C_imp (n_temp_min_imp >= -256).

integer(i4b) :: n_temp_max_imp
   !! Upper index limit of properly defined values in
   !! RF_imp, KAPPA_imp and C_imp (n_temp_max_imp <= 255).

real(dp) :: RHO_I_imp
   !! Density of ice
   !! (only for the Martian ice caps)

real(dp) :: RHO_C_imp
   !! Density of crustal material (dust)
   !! (only for the Martian ice caps)

real(dp) :: KAPPA_C_imp
   !! Heat conductivity of crustal material (dust)
   !! (only for the Martian ice caps)

real(dp) :: C_C_imp
   !! Specific heat of crustal material (dust)
   !! (only for the Martian ice caps)

!-------- Auxiliary variables for the module enth_temp_omega_m
!         (converting temperature and water content
!                                 to enthalpy and vice versa) --------

real(dp), dimension(-256:255) :: c_int_table
   !! Temperature integral of the specific heat of ice;
   !! index is temperature in degC

real(dp), dimension(-524288:524287) :: c_int_inv_table
   !! Inverse of the temperature integral of the specific heat of ice;
   !! index is enthalpy in J/kg (zero for 0 degC)

integer(i4b) :: n_temp_min
   !! Lower index limit of properly defined values in c_int_table
   !! (n_temp_min >= -256)

integer(i4b) :: n_temp_max
   !! Upper index limit of properly defined values in c_int_table
   !! (n_temp_max <= 255)

integer(i4b) :: n_enth_min
   !! Lower index limit of properly defined values in c_int_inv_table
   !! (n_enth_min >= -524288)

integer(i4b) :: n_enth_max
   !! Upper index limit of properly defined values in c_int_inv_table
   !! (n_enth_max <= 524287)

real(dp) :: latent_heat
   !! Latent heat of ice

real(dp) :: latent_heat_inv
   !! Inverse of the latent heat of ice

!-------- Temperature-to-precipitation conversion --------

#if (ACCSURFACE==2 || ACCSURFACE==3)

real(dp), dimension(0:JMAX,0:IMAX) :: gamma_s
   !! Temperature-to-precipitation conversion coefficient

#endif

!-------- PDD and LTI parameters --------

#if (ABLSURFACE==1 || ABLSURFACE==2)

real(dp), dimension(0:JMAX,0:IMAX) :: beta1
   !! PDD factor for snow melt

real(dp), dimension(0:JMAX,0:IMAX) :: beta2
   !! PDD factor for ice melt

real(dp), dimension(0:JMAX,0:IMAX) :: s_stat
   !! Standard deviation of air-temperature fluctuations

real(dp), dimension(0:JMAX,0:IMAX) :: Pmax
   !! Saturation factor for the formation of superimposed ice

real(dp), dimension(0:JMAX,0:IMAX) :: mu
   !! Firn-warming correction

#elif (ABLSURFACE==3)

real(dp), dimension(0:JMAX,0:IMAX) :: lambda_lti
   !! Melting coefficient for the LTI method

real(dp), dimension(0:JMAX,0:IMAX) :: temp_lti
   !! Threshold summer temperature for the LTI method

#endif

!-------- ISMIP6-like climate forcing --------

#if (TSURFACE==6 && ACCSURFACE==6 && ABLSURFACE==6)

real(dp), dimension(0:JMAX,0:IMAX) :: temp_maat_climatol
   !! Surface-temperature (MAAT) climatology

real(dp), dimension(0:JMAX,0:IMAX) :: smb_climatol
   !! SMB climatology

real(dp), dimension(0:JMAX,0:IMAX) :: temp_maat_anom
   !! Surface-temperature (MAAT) anomaly

real(dp), dimension(0:JMAX,0:IMAX) :: smb_anom
   !! SMB anomaly

real(dp), dimension(0:JMAX,0:IMAX) :: dtemp_maat_dz
   !! Surface-temperature (MAAT) vertical gradient

real(dp), dimension(0:JMAX,0:IMAX) :: dsmb_dz
   !! SMB vertical gradient

#endif

!-------- ISMIP6 InitMIP --------

#if (defined(ANT) || defined(GRL)) /* Antarctica or Greenland */

logical :: flag_initmip_asmb
   !! Flag for use of InitMIP SMB anomaly

real(dp), dimension(0:JMAX,0:IMAX) :: smb_anom_initmip
   !! InitMIP SMB anomaly

#endif

#if (defined(ANT)) /* Antarctica */

logical :: flag_initmip_abmb
   !! Flag for use of InitMIP sub-ice-shelf-melt anomaly

real(dp), dimension(0:JMAX,0:IMAX) :: ab_anom_initmip
   !! InitMIP sub-ice-shelf-melt anomaly

logical :: flag_larmip
   !! Flag for use of LARMIP sub-ice-shelf-melt anomaly

integer(i4b), dimension(0:JMAX,0:IMAX) :: n_larmip_region
   !! LARMIP regions for ice shelf basal melting

real(dp), dimension(0:7) :: ab_anom_larmip
   !! LARMIP sub-ice-shelf-melt anomaly

#endif

!-------- ISMIP6-like oceanic forcing --------

#if (FLOATING_ICE_BASAL_MELTING==6)

real(dp), dimension(0:NZ_TF_BM) :: z_tf_bm_present
   !! Equidistant depth points of the
   !! present-day thermal forcing data of the ocean

real(dp), dimension(0:NZ_TF_BM,0:JMAX,0:IMAX) :: tf_bm_present
   !! Present-day thermal forcing data of the ocean

real(dp), dimension(0:NZ_TF_BM) :: z_tf_bm
   !! Equidistant depth points of the
   !! thermal forcing data of the ocean

real(dp), dimension(0:NZ_TF_BM,0:JMAX,0:IMAX) :: tf_bm
   !! Thermal forcing data of the ocean

#endif

!-------- ISMIP6-like prescribed ice-shelf collapse
!                                or grounded-ice retreat --------

#if (defined(ANT) && ICE_SHELF_COLLAPSE_MASK==1) /* Antarctica */

real(dp), dimension(0:JMAX,0:IMAX) :: H_ref_retreat
   !! Reference ice thickness for the ice-shelf collapse mask

real(dp), dimension(0:JMAX,0:IMAX) :: r_mask_retreat
   !! Ice-shelf collapse mask

#endif

#if (defined(GRL) && RETREAT_MASK==1) /* Greenland */

real(dp), dimension(0:JMAX,0:IMAX) :: H_ref_retreat
   !! Reference ice thickness for the retreat mask

real(dp), dimension(0:JMAX,0:IMAX) :: r_mask_retreat
   !! Retreat mask

#endif

!-------- Ice discharge parameterization for Greenland --------

#if (defined(GRL) && DISC>0) /* Greenland */

integer(i4b) :: disc
integer(i4b) :: n_discharge_call
integer(i4b) :: iter_mar_coa
real(dp)     :: c_dis_0
real(dp)     :: s_dis
real(dp)     :: c_dis_fac
real(dp)     :: T_sub_PD
real(dp)     :: alpha_sub
real(dp)     :: alpha_o
real(dp)     :: m_H
real(dp)     :: m_D
real(dp)     :: r_mar_eff
real(dp)     :: T_sea_freeze
real(dp)     :: dT_glann
real(dp)     :: dT_sub

integer(i4b), dimension(0:JMAX,0:IMAX) :: mask_mar
real(dp),     dimension(0:JMAX,0:IMAX) :: c_dis
real(dp),     dimension(0:JMAX,0:IMAX) :: cst_dist
real(dp),     dimension(0:JMAX,0:IMAX) :: cos_grad_tc
real(dp),     dimension(0:JMAX,0:IMAX) :: dis_perp

#if (DISC==2)

integer(i4b) :: glann_time_min
   !! Minimum time of the data values for the
   !! global annual temperature anomaly

integer(i4b) :: glann_time_stp
   !! Time step of the data values for the
   !! global annual temperature anomaly

integer(i4b) :: glann_time_max
   !! Maximum time of the data values for the
   !! global annual temperature anomaly

integer(i4b) :: ndata_glann
   !! Number of data values for the global annual temperature anomaly

integer(i4b) , parameter :: ndata_glann_max = 262143
   !! Maximum allowed value of ndata_glann

real(dp), dimension(0:ndata_glann_max) :: dT_glann_CLIMBER
   !! Data values for the global annual temperature anomaly

#endif

#endif

!-------- Austfonna: Additional output for prescribed sites --------

#if (defined(ASF) && WRITE_SER_FILE_STAKES==1) /* Austfonna */

integer(i4b) :: n_surf
   !! Number of surface points for which time-series data are written 

integer(i4b), parameter :: n_surf_max = 256
   !! Maximum allowed value of n_surf

real(dp), dimension(n_surf_max) :: lambda_surf
   !! Geographical longitude of the prescribed surface points

real(dp), dimension(n_surf_max) :: phi_surf
   !! Geographical latitude of the prescribed surface points

real(dp), dimension(n_surf_max) :: x_surf
   !! Coordinate xi (= x) of the prescribed surface points

real(dp), dimension(n_surf_max) :: y_surf
   !! Coordinate eta (= y) of the prescribed surface points

#endif

!-------- Mathematical constants --------

real(dp), parameter :: pi = 3.141592653589793_dp
   !! Constant pi

real(dp), parameter :: deg2rad = pi/180.0_dp
   !! pi divided by 180 (-> deg to rad)

real(dp), parameter :: rad2deg = 180.0_dp/pi
   !! 180 divided by pi (-> rad to deg)

real(dp), parameter :: euler = 2.718281828459045_dp
   !! Euler number

real(dp), parameter :: eps = 1.0e-05_dp
   !! Small number

real(dp), parameter :: epsi = 1.0e-12_dp
   !! Very small number

#if (!defined(ALLOW_GRDCHK) && !defined(ALLOW_TAPENADE)) /* Normal */

real(sp), parameter :: eps_sp = epsilon(1.0_sp)
   !! Small number to single-precision accuracy

real(dp), parameter :: eps_sp_dp = eps_sp
   !! Small number to single-precision accuracy in double precision

real(dp), parameter :: eps_dp = epsilon(1.0_dp)
   !! Small number to double-precision accuracy

#else /* Tapenade */

real(sp), parameter :: myepsilon_sp  = 1.19209290E-07

real(sp), parameter :: eps_sp = myepsilon_sp

real(dp), parameter :: myepsilon_sp_dp  = 1.1920928955078125E-007

real(dp), parameter :: eps_sp_dp = myepsilon_sp_dp

real(dp), parameter :: myepsilon_dp  = 2.2204460492503131E-016

real(dp), parameter :: eps_dp = myepsilon_dp

#endif /* Normal vs. Tapenade */

!-------- Further quantities --------

real(dp) :: year2sec
   !! 1 year (1 a) in seconds

real(dp) :: sec2year
   !! 1 second in years

real(dp) :: year_zero
   !! SICOPOLIS year zero in astronomical year numbering
   !! [ = signed year CE (AD) ]

character(len=256) :: run_name
   !! Name of simulation

character(len=64) :: ch_domain_long
   !! Long name of the computational domain

character(len=16) :: ch_domain_short
   !! Short name of the computational domain

integer(i4b) :: forcing_flag
   !! Flag for the forcing type:
   !!  1: forcing by a spatially constant surface temperature
   !!     anomaly (delta_ts),
   !!  2: forcing by a glacial index (glac_index),
   !!  3: forcing by time-dependent surface temperature
   !!     and precipitation data

integer(i4b) :: n_site
   !! Number of positions to be considered in the time-series file
   !! for specified sites (i.e., ice cores)

integer(i4b), parameter :: n_site_max = 256
   !! Maximum allowed value of n_site

real(dp), dimension(n_site_max) :: lambda_site
   !! Geographical longitude of the prescribed sites

real(dp), dimension(n_site_max) :: phi_site
   !! Geographical latitude of the prescribed sites

real(dp), dimension(n_site_max) :: x_site
   !! Coordinate xi (= x) of the prescribed sites

real(dp), dimension(n_site_max) :: y_site
   !! Coordinate eta (= y) of the prescribed sites

character(len=16), dimension(n_site_max) :: ch_site
   !! Names of the prescribed sites

integer(i4b) :: grip_time_min
   !! Minimum time of the data values for the surface temperature anomaly

integer(i4b) :: grip_time_stp
   !! Time step of the data values for the surface temperature anomaly

integer(i4b) :: grip_time_max
   !! Maximum time of the data values for the surface temperature anomaly

integer(i4b) :: ndata_grip
   !! Number of data values for the surface temperature anomaly

integer(i4b), parameter :: ndata_grip_max = 262143
   !! Maximum allowed value of ndata_grip

real(dp), dimension(0:ndata_grip_max):: griptemp
   !! Data values for the surface temperature anomaly

integer(i4b) :: gi_time_min
   !! Minimum time of the data values for the glacial index

integer(i4b) :: gi_time_stp
   !! Time step of the data values for the glacial index

integer(i4b) :: gi_time_max
   !! Maximum time of the data values for the glacial index

integer(i4b) :: ndata_gi
   !! Number of data values for the glacial index

integer(i4b), parameter :: ndata_gi_max = 262143
   !! Maximum allowed value of ndata_gi

real(dp), dimension(0:ndata_gi_max) :: glacial_index
   !! Data values for the glacial index

integer(i4b) :: specmap_time_min
   !! Minimum time of the data values for the sea level

integer(i4b) :: specmap_time_stp
   !! Time step of the data values for the sea level

integer(i4b) :: specmap_time_max
   !! Maximum time of the data values for the sea level

integer(i4b) :: ndata_specmap
   !! Number of data values for the sea level

integer(i4b), parameter :: ndata_specmap_max = 16383
   !! Maximum allowed value of ndata_specmap

real(dp), dimension(0:ndata_specmap_max) :: specmap_zsl
   !! Data values for the sea level

integer(i4b) :: target_topo_tau0_time_min
   !! Minimum time of the data values
   !! for the relaxation time for the topography nudging

integer(i4b) :: target_topo_tau0_time_stp
   !! Time step of the data values
   !! for the relaxation time for the topography nudging

integer(i4b) :: target_topo_tau0_time_max
   !! Maximum time of the data values
   !! for the relaxation time for the topography nudging

integer(i4b) :: ndata_target_topo_tau0
   !! Number of data values
   !! for the relaxation time for the topography nudging

integer(i4b), parameter :: ndata_target_topo_tau0_max = 16383
   !! Maximum allowed value of ndata_target_topo_tau0

real(dp), dimension(0:ndata_target_topo_tau0_max) :: target_topo_tau0
   !! Data values for the relaxation time for the topography nudging

real(dp) :: target_topo_tau_0
   !! Constant relaxation time for the topography nudging

real(dp) :: target_topo_tau
   !! Relaxation time for the topography nudging

integer(i4b), dimension(0:JMAX,0:IMAX) :: mask_target
   !! Target topography (ice-land-ocean mask)

real(dp), dimension(0:JMAX,0:IMAX) :: zs_target
   !! Target topography (ice surface)

real(dp), dimension(0:JMAX,0:IMAX) :: zb_target
   !! Target topography (ice base)

real(dp), dimension(0:JMAX,0:IMAX) :: zl_target
   !! Target topography (lithosphere surface)

real(dp), dimension(0:JMAX,0:IMAX) :: H_target
   !! Target topography (ice thickness)

logical :: flag_mask_maxextent
   !! Flag for maximum ice extent:
   !!   .true.: specified,
   !!  .false.: not specified

integer(i4b), dimension(0:JMAX,0:IMAX) :: mask_maxextent
   !! Maximum ice extent mask:
   !!  0: not allowed to glaciate,
   !!  1: allowed to glaciate

type(flag_firstcall) :: firstcall
   !! First-call flags for several routines

integer(i4b) :: n_year_CE_surf_clim_save = -9999
   !! Time (in years CE) for ISMIP6-type surface climate forcing data

integer(i4b) :: n_year_CE_bas_melt_save = -9999
   !! Time (in years CE) for the thermal forcing data of the ocean
   !! for the ISMIP6 non-local sub-ice-shelf melting parameterization

integer(i4b) :: n_year_CE_isc_save = -9999
   !! Time (in years CE) for ISMIP6-type ice-shelf collapse masks

integer(i4b) :: n_year_CE_rtr_save = -9999
   !! Time (in years CE) for ISMIP6-type ice-sheet retreat masks

integer(i4b), dimension(0:99) :: ncid_ser
   !! IDs of the NetCDF time-series output files

integer(i4b) :: ncid_site
   !! ID of the NetCDF time-series output file for the specified sites
   !! (i.e., ice cores)

real(dp), dimension(-10000:10000) :: kei
   !! Tabulated values of the kei function (Kelvin function of zero order)

integer(i4b):: n_data_kei
   !! Number of tabulated values of the kei function

real(dp) :: kei_r_max
   !! Maximum value of the argument r of the tabulated kei function

real(dp) :: kei_r_incr
   !! Increment of the argument r of the tabulated kei function

integer(i4b), parameter :: rcl1 = 3*8*(IMAX+1)
   !! Maximum length of record for input files
   !! (with factor 3 safety margin)

integer(i4b), parameter :: rcl2 = 3  *(IMAX+1)
   !! Maximum length of record for input mask files
   !! (with factor 3 safety margin)

integer(i4b), dimension(0:JMAX,0:IMAX) :: ij2n
   !! Conversion from 2d index pair (i,j) to linear index n

integer(i4b), dimension((IMAX+1)*(JMAX+1)) :: n2i
   !! Conversion from linear index n to 2d index i

integer(i4b), dimension((IMAX+1)*(JMAX+1)) :: n2j
   !! Conversion from linear index n to 2d index j

logical :: flag_grads_nc_tweaks
   !! Flag for optimizing NetCDF output for viewing with GrADS:
   !!   .true.: optimized output,
   !!  .false.: normal NetCDF output (default)

real(dp), parameter :: no_value_pos_1 =  1.11e+11_dp
   !! Positive no-value parameter

real(dp), parameter :: no_value_pos_2 =  9.999e+03_dp
   !! Positive no-value parameter

real(dp), parameter :: no_value_neg_1 = -1.11e+11_dp
   !! Negative no-value parameter

real(dp), parameter :: no_value_neg_2 = -9.999e+03_dp
   !! Negative no-value parameter

character(len=256) :: errormsg
   !! Error-message string

character, parameter :: end_of_line = char(10)
   !! End-of-line string

#if (defined(ALLOW_NORMAL) || defined(ALLOW_GRDCHK) || defined(ALLOW_TAPENADE)) /* Tapenade */

real(dp) :: fc
   !! Scalar cost function

#if (defined(AGE_COST))

! Note: for the age cost, CALCMOD!=1 is recommended because
! the gridded ages of the GRL ice sheet are only 25 z-levels.

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: age_data
   !! Array of gridded ages to be used in adjoint mode

real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: age_unc
   !! Array of gridded uncertainty in ages to be used in adjoint mode

real(dp), dimension(0:JMAX,0:IMAX) :: H_data
   !! Array of gridded ice thickness to be used in adjoint mode

real(dp), dimension(0:JMAX,0:IMAX) :: H_unc
   !! Array of gridded uncertainty in ice thickness to be used
   !! in adjoint mode

#endif /* No age cost used */

real(dp), dimension(0:JMAX,0:IMAX) :: acc_fact

#if (defined(BEDMACHINE_COST)) 

real(dp), dimension(0:JMAX,0:IMAX) :: H_BedMachine_data

real(dp), dimension(0:JMAX,0:IMAX) :: H_unc_BedMachine_data

#endif

#endif /* Tapenade */

end module sico_variables_m
!
