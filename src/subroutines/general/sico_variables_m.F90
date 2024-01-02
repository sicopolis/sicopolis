!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ v a r i a b l e s _ m
!
!> @file
!!
!! Declarations of global variables for SICOPOLIS.
!!
!! @section Copyright
!!
!! Copyright 2009-2024 Ralf Greve
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
!> Declarations of global variables for SICOPOLIS.
!<------------------------------------------------------------------------------
module sico_variables_m

use sico_types_m

implicit none
save

!-------- Field quantities --------

!> mask(j,i): Ice-land-ocean mask.
!>             0: grounded ice,
!>             1: ice-free land,
!>             2: ocean,
!>             3: floating ice
   integer(i4b), dimension(0:JMAX,0:IMAX) :: mask
!> mask_old(j,i): Old value of mask (at the previous time step)
   integer(i4b), dimension(0:JMAX,0:IMAX) :: mask_old
!> mask_new(j,i): New value of mask computed during an integration step
   integer(i4b), dimension(0:JMAX,0:IMAX) :: mask_new
!> n_cts(j,i): Mask for thermal conditions.
!>             -1: cold ice base,
!>              0: temperate ice base with cold ice above,
!>              1: temperate ice base with temperate ice layer above
!>                 (only for POLY)
   integer(i4b), dimension(0:JMAX,0:IMAX) :: n_cts
!> (.)_new: New value of quantity (.) computed during an integration step
   integer(i4b), dimension(0:JMAX,0:IMAX) :: n_cts_new
!> kc_cts(j,i): Position kc of the CTS (for COLD, ENTC, ENTM)
   integer(i4b), dimension(0:JMAX,0:IMAX) :: kc_cts
!> (.)_new: New value of quantity (.) computed during an integration step
   integer(i4b), dimension(0:JMAX,0:IMAX) :: kc_cts_new
!> mask_region(j,i): Region mask.
!>             0: undefined,
!>             1: EAIS,
!>             2: WAIS,
!>             3: AP
   integer(i4b), dimension(0:JMAX,0:IMAX) :: mask_region
!> flag_calc_temp: Flag for computation of the temperature, water content,
!>                 age and flow enhancement factor during an integration step.
!>                  .true.: temperature etc. computed
!>                 .false.: temperature etc. not computed
   logical :: flag_calc_temp
!> flag_inner_point(j,i): Inner-point flag.
!>                              .true.: inner point,
!>                             .false.: margin point
   logical, dimension(0:JMAX,0:IMAX) :: flag_inner_point
!> flag_grounding_line_1(j,i): Grounding line flag.
!>                              .true.: grounding line point
!>                                      (grounded ice point with at least
!>                                      one floating ice neighbour),
!>                             .false.: otherwise
   logical, dimension(0:JMAX,0:IMAX) :: flag_grounding_line_1
!> flag_grounding_line_2(j,i): Grounding line flag.
!>                              .true.: grounding line point
!>                                      (floating ice point with at least
!>                                      one grounded ice neighbour),
!>                             .false.: otherwise
   logical, dimension(0:JMAX,0:IMAX) :: flag_grounding_line_2
!> flag_calving_front_1(j,i): Calving front flag.
!>                              .true.: calving front point
!>                                      (floating ice point with at least
!>                                      one ocean neighbour),
!>                             .false.: otherwise
   logical, dimension(0:JMAX,0:IMAX) :: flag_calving_front_1
!> flag_calving_front_2(j,i): Calving front flag.
!>                              .true.: calving front point
!>                                      (ocean point with at least
!>                                      one floating ice neighbour),
!>                             .false.: otherwise
   logical, dimension(0:JMAX,0:IMAX) :: flag_calving_front_2
!> flag_grounded_front_a_1(j,i): Land-terminating grounded front flag.
!>                                .true.: grounded front point
!>                                        (grounded ice point with at least
!>                                        one ice-free land neighbour),
!>                               .false.: otherwise
   logical, dimension(0:JMAX,0:IMAX) :: flag_grounded_front_a_1
!> flag_grounded_front_a_2(j,i): Land-terminating grounded front flag.
!>                                .true.: grounded front point
!>                                        (ice-free land point with at least
!>                                        one grounded ice neighbour),
!>                               .false.: otherwise
   logical, dimension(0:JMAX,0:IMAX) :: flag_grounded_front_a_2
!> flag_grounded_front_b_1(j,i): Marine-terminating grounded front flag.
!>                                .true.: grounded front point
!>                                        (grounded ice point with at least
!>                                        one ocean neighbour),
!>                               .false.: otherwise
   logical, dimension(0:JMAX,0:IMAX) :: flag_grounded_front_b_1
!> flag_grounded_front_b_2(j,i): Marine-terminating grounded front flag.
!>                                .true.: grounded front point
!>                                        (ocean point with at least
!>                                        one grounded ice neighbour),
!>                               .false.: otherwise
   logical, dimension(0:JMAX,0:IMAX) :: flag_grounded_front_b_2
!> flag_shelfy_stream_x(j,i): Shelfy stream flag in x-direction, at (i+1/2,j).
!>                             .true.: shelfy stream point
!>                            .false.: otherwise
   logical, dimension(0:JMAX,0:IMAX) :: flag_shelfy_stream_x
!> flag_shelfy_stream_y(j,i): Shelfy stream flag in y-direction, at (i,j+1/2).
!>                             .true.: shelfy stream point
!>                            .false.: otherwise
   logical, dimension(0:JMAX,0:IMAX) :: flag_shelfy_stream_y
!> flag_shelfy_stream(j,i):   Shelfy stream flag on the main grid.
!>                             .true.: grounded ice,
!>                                     and at least one neighbour on the
!>                                     staggered grid is a shelfy stream point
!>                            .false.: otherwise
   logical, dimension(0:JMAX,0:IMAX) :: flag_shelfy_stream
!> xi(i): Coordinate xi (= x) of grid point i
   real(dp), dimension(0:IMAX) :: xi
!> eta(j): Coordinate eta (= y) of grid point j
   real(dp), dimension(0:JMAX) :: eta
!> zeta_c(kc): Sigma coordinate zeta_c of grid point kc
   real(dp), dimension(0:KCMAX) :: zeta_c
!> zeta_t(kt): Sigma coordinate zeta_t of grid point kt
   real(dp), dimension(0:KTMAX) :: zeta_t
!> zeta_r(kr): Sigma coordinate zeta_r of grid point kr
   real(dp), dimension(0:KRMAX) :: zeta_r
!> aa: Exponential stretch parameter of the non-equidistant vertical grid
!>     in the upper (kc) ice domain
   real(dp) :: aa
!> flag_aa_nonzero: Flag for the exponential stretch parameter aa.
!>                    .true.: aa greater than zero (non-equidistant grid)
!>                   .false.: aa equal to zero (equidistant grid)
   logical :: flag_aa_nonzero
!> ea: Abbreviation for exp(aa)
   real(dp) :: ea
!> eaz_c(kc): Abbreviation for exp(aa*zeta(kc))
   real(dp), dimension(0:KCMAX) :: eaz_c
!> eaz_c_quotient(kc): Abbreviation for (eaz_c(kc)-1.0)/(ea-1.0)
   real(dp), dimension(0:KCMAX) :: eaz_c_quotient

!> lambda(j,i): Geographic longitude of grid point (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: lambda
!> phi(j,i): Geographic latitude of grid point (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: phi
!> cell_area(j,i): Area of grid cell associated with grid point (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: cell_area
!> sq_g11_g(j,i): Square root of the coefficient g11 of the metric tensor
!>                on grid point (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: sq_g11_g
!> sq_g22_g(j,i): Square root of the coefficient g22 of the metric tensor
!>                on grid point (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: sq_g22_g
!> insq_g11_g(j,i): Inverse square root of g11 on grid point (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: insq_g11_g
!> insq_g22_g(j,i): Inverse square root of g22 on grid point (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: insq_g22_g
!> sq_g11_sgx(j,i): Square root of g11, at (i+1/2,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: sq_g11_sgx
!> sq_g11_sgy(j,i): Square root of g11, at (i,j+1/2)
   real(dp), dimension(0:JMAX,0:IMAX) :: sq_g11_sgy
!> sq_g22_sgx(j,i): Square root of g22, at (i+1/2,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: sq_g22_sgx
!> sq_g22_sgy(j,i): Square root of g22, at (i,j+1/2)
   real(dp), dimension(0:JMAX,0:IMAX) :: sq_g22_sgy
!> insq_g11_sgx(j,i): Inverse square root of g11, at (i+1/2,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: insq_g11_sgx
!> insq_g22_sgy(j,i): Inverse square root of g22, at (i,j+1/2)
   real(dp), dimension(0:JMAX,0:IMAX) :: insq_g22_sgy
!> zs(j,i): Coordinate z of the surface topography
   real(dp), dimension(0:JMAX,0:IMAX) :: zs
!> zm(j,i): Coordinate z of the bottom of the upper (kc) ice domain
!>                               = top of the lower (kt) ice domain
!>          (position of the CTS for POLY,
!>          equal to zb for ISOT, COLD, ENTC, ENTM)
   real(dp), dimension(0:JMAX,0:IMAX) :: zm
!> zb(j,i): Coordinate z of the ice base
   real(dp), dimension(0:JMAX,0:IMAX) :: zb
!> zl(j,i): Coordinate z of the lithosphere surface
   real(dp), dimension(0:JMAX,0:IMAX) :: zl
!> zl0(j,i): zl for isostatically relaxed ice-free conditions
   real(dp), dimension(0:JMAX,0:IMAX) :: zl0
!> wss(j,i): Isostatic steady-state displacement of the lithosphere
   real(dp), dimension(0:JMAX,0:IMAX) :: wss
!> flex_rig_lith(j,i): Flexural rigidity of the lithosphere
   real(dp), dimension(0:JMAX,0:IMAX) :: flex_rig_lith
!> time_lag_asth(j,i): Time lag of the relaxing asthenosphere
   real(dp), dimension(0:JMAX,0:IMAX) :: time_lag_asth
!> H(j,i): Ice thickness (= H_c + H_t)
   real(dp), dimension(0:JMAX,0:IMAX) :: H
!> H_c(j,i): Thickness of ice in the upper (kc) domain
!>           (thickness of the cold-ice layer for POLY,
!>           entire ice thickness for ISOT, COLD, ENTC, ENTM)
   real(dp), dimension(0:JMAX,0:IMAX) :: H_c
!> H_t(j,i): Thickness of ice in the lower (kt) domain
!>           (thickness of the temperate layer for POLY,
!>           redundant and thus set to zero for ISOT, COLD, ENTC, ENTM)
   real(dp), dimension(0:JMAX,0:IMAX) :: H_t
!> dzs_dxi(j,i): Derivative of zs by xi (at i+1/2,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: dzs_dxi
!> dzm_dxi(j,i): Derivative of zm by xi (at i+1/2,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: dzm_dxi
!> dzb_dxi(j,i): Derivative of zb by xi (at i+1/2,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: dzb_dxi
!> dH_c_dxi(j,i): Derivative of H_c by xi (at i+1/2,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: dH_c_dxi
!> dH_t_dxi(j,i): Derivative of H_t by xi (at i+1/2,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: dH_t_dxi
!> dzs_deta(j,i): Derivative of zs by eta (at i,j+1/2)
   real(dp), dimension(0:JMAX,0:IMAX) :: dzs_deta
!> dzm_deta(j,i): Derivative of zm by eta (at i,j+1/2)
   real(dp), dimension(0:JMAX,0:IMAX) :: dzm_deta
!> dzb_deta(j,i): Derivative of zb by eta (at i,j+1/2)
   real(dp), dimension(0:JMAX,0:IMAX) :: dzb_deta
!> dH_c_deta(j,i): Derivative of H_c by eta (at i,j+1/2)
   real(dp), dimension(0:JMAX,0:IMAX) :: dH_c_deta
!> dH_t_deta(j,i): Derivative of H_t by eta (at i,j+1/2)
   real(dp), dimension(0:JMAX,0:IMAX) :: dH_t_deta
!> (.)_g(j,i): Staggered-grid quantity (.) interpolated to grid point (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: dzs_dxi_g
!> (.)_g(j,i): Staggered-grid quantity (.) interpolated to grid point (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: dzm_dxi_g
!> (.)_g(j,i): Staggered-grid quantity (.) interpolated to grid point (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: dzb_dxi_g
!> (.)_g(j,i): Staggered-grid quantity (.) interpolated to grid point (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: dH_c_dxi_g
!> (.)_g(j,i): Staggered-grid quantity (.) interpolated to grid point (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: dH_t_dxi_g
!> (.)_g(j,i): Staggered-grid quantity (.) interpolated to grid point (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: dzs_deta_g
!> (.)_g(j,i): Staggered-grid quantity (.) interpolated to grid point (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: dzm_deta_g
!> (.)_g(j,i): Staggered-grid quantity (.) interpolated to grid point (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: dzb_deta_g
!> (.)_g(j,i): Staggered-grid quantity (.) interpolated to grid point (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: dH_c_deta_g
!> (.)_g(j,i): Staggered-grid quantity (.) interpolated to grid point (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: dH_t_deta_g
!> dzs_dtau(j,i): Derivative of zs by tau (time)
   real(dp), dimension(0:JMAX,0:IMAX) :: dzs_dtau
!> dzm_dtau(j,i): Derivative of zm by tau (time)
   real(dp), dimension(0:JMAX,0:IMAX) :: dzm_dtau
!> dzb_dtau(j,i): Derivative of zb by tau (time)
   real(dp), dimension(0:JMAX,0:IMAX) :: dzb_dtau
!> dzl_dtau(j,i): Derivative of zl by tau (time)
   real(dp), dimension(0:JMAX,0:IMAX) :: dzl_dtau
!> dH_dtau(j,i): Derivative of H by tau (time)
   real(dp), dimension(0:JMAX,0:IMAX) :: dH_dtau
!> dH_c_dtau(j,i): Derivative of H_c by tau (time)
   real(dp), dimension(0:JMAX,0:IMAX) :: dH_c_dtau
!> dH_t_dtau(j,i): Derivative of H_t by tau (time)
   real(dp), dimension(0:JMAX,0:IMAX) :: dH_t_dtau
!> n_slide_region(j,i): Regions for the basal sliding laws.
   integer(i4b), dimension(0:JMAX,0:IMAX) :: n_slide_region
!> p_weert(j,i): Weertman exponent for the basal shear stress
   integer(i4b), dimension(0:JMAX,0:IMAX) :: p_weert
!> q_weert(j,i): Weertman exponent for the basal pressure
   integer(i4b), dimension(0:JMAX,0:IMAX) :: q_weert
!> p_weert_inv(j,i): Inverse of p_weert
   real(dp), dimension(0:JMAX,0:IMAX) :: p_weert_inv
!> c_slide(j,i): Basal sliding coefficient
   real(dp), dimension(0:JMAX,0:IMAX) :: c_slide
!> c_slide_init(j,i): Initial basal sliding coefficient
   real(dp), dimension(0:JMAX,0:IMAX) :: c_slide_init
!> gamma_slide_inv(j,i): Inverse of the sub-melt-sliding parameter
   real(dp), dimension(0:JMAX,0:IMAX) :: gamma_slide_inv
!> sub_melt_flag(j,i): Flag for presence of sub-melt sliding
   logical, dimension(0:JMAX,0:IMAX) :: sub_melt_flag
!> d_help_b(j,i): Auxiliary quantity for the computation of vx_b and vy_b
   real(dp), dimension(0:JMAX,0:IMAX) :: d_help_b
!> c_drag(j,i): Auxiliary quantity for the computation of the basal drag
   real(dp), dimension(0:JMAX,0:IMAX) :: c_drag
!> p_b(j,i): Basal pressure
   real(dp), dimension(0:JMAX,0:IMAX) :: p_b
!> p_b_w(j,i): Basal water pressure
   real(dp), dimension(0:JMAX,0:IMAX) :: p_b_w
!> p_b_red(j,i): Reduced basal pressure
   real(dp), dimension(0:JMAX,0:IMAX) :: p_b_red
!> tau_dr(j,i): Driving stress
   real(dp), dimension(0:JMAX,0:IMAX) :: tau_dr
!> tau_b(j,i): Basal shear stress (drag)
   real(dp), dimension(0:JMAX,0:IMAX) :: tau_b
!> vx_b(j,i): Velocity in x-direction at the ice base, at (i+1/2,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: vx_b
!> vy_b(j,i): Velocity in y-direction at the ice base, at (i,j+1/2)
   real(dp), dimension(0:JMAX,0:IMAX) :: vy_b
!> vx_m(j,i): Mean (depth-averaged) velocity in x-direction, at (i+1/2,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: vx_m
!> vy_m(j,i): Mean (depth-averaged) velocity in y-direction, at (i,j+1/2)
   real(dp), dimension(0:JMAX,0:IMAX) :: vy_m
!> vx_m_sia(j,i): Mean (depth-averaged) SIA velocity in x-direction,
!>                at (i+1/2,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: vx_m_sia
!> vy_m_sia(j,i): Mean (depth-averaged) SIA velocity in y-direction,
!>                at (i,j+1/2)
   real(dp), dimension(0:JMAX,0:IMAX) :: vy_m_sia

!> vx_m_ssa(j,i): Mean (depth-averaged) SSA velocity in x-direction,
!>                at (i+1/2,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: vx_m_ssa
!> vy_m_ssa(j,i): Mean (depth-averaged) SSA velocity in y-direction,
!>                at (i,j+1/2)
   real(dp), dimension(0:JMAX,0:IMAX) :: vy_m_ssa
!> ratio_sl_x(j,i): Ratio of basal to surface velocity (slip ratio)
!>                  in x-direction, at (i+1/2,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: ratio_sl_x
!> ratio_sl_y(j,i): Ratio of basal to surface velocity (slip ratio)
!>                  in y-direction, at (i,j+1/2)
   real(dp), dimension(0:JMAX,0:IMAX) :: ratio_sl_y
!> ratio_sl(j,i): Ratio of basal to surface velocity (slip ratio)
!>                on the main grid
   real(dp), dimension(0:JMAX,0:IMAX) :: ratio_sl
!> (.)_g(j,i): Staggered-grid quantity (.) interpolated to grid point (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: vx_b_g
!> (.)_g(j,i): Staggered-grid quantity (.) interpolated to grid point (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: vy_b_g
!> vz_b(j,i): Velocity in z-direction at the ice base, at (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: vz_b
!> vz_m(j,i): Velocity in z-direction at the position z=zm (interface between
!>            the upper (kc) and the lower (kt) domain), at (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: vz_m
!> vx_s_g(j,i): Velocity in x-direction at the ice surface, at (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: vx_s_g
!> vy_s_g(j,i): Velocity in x-direction at the ice surface, at (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: vy_s_g
!> vz_s(j,i): Velocity in z-direction at the ice surface, at (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: vz_s
!> flui_ave_sia(j,i): Depth-averaged fluidity of the SIA
   real(dp), dimension(0:JMAX,0:IMAX) :: flui_ave_sia
!> h_diff(j,i): Diffusivity of the SIA evolution equation of the ice surface
   real(dp), dimension(0:JMAX,0:IMAX) :: h_diff
!> qx(j,i): Volume flux in x-direction (depth-integrated vx, at (i+1/2,j))
   real(dp), dimension(0:JMAX,0:IMAX) :: qx
!> qy(j,i): Volume flux in y-direction (depth-integrated vy, at (i,j+1/2))
   real(dp), dimension(0:JMAX,0:IMAX) :: qy
!> q_gl_g(j,i): Volume flux across the grounding line, at (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: q_gl_g
!> q_geo(j,i): Geothermal heat flux
   real(dp), dimension(0:JMAX,0:IMAX) :: q_geo
!> temp_b(j,i): Basal temperature
   real(dp), dimension(0:JMAX,0:IMAX) :: temp_b
!> temph_b(j,i): Basal temperature relative to the pressure melting point
   real(dp), dimension(0:JMAX,0:IMAX) :: temph_b
!> Q_bm(j,i): Basal melting rate
   real(dp), dimension(0:JMAX,0:IMAX) :: Q_bm
!> Q_tld(j,i): Water drainage rate from the temperate layer
   real(dp), dimension(0:JMAX,0:IMAX) :: Q_tld
!> Q_b_tot(j,i): Sum of Q_bm and Q_tld
   real(dp), dimension(0:JMAX,0:IMAX) :: Q_b_tot
!> q_w(j,i): Scalar volume flux of the basal water
   real(dp), dimension(0:JMAX,0:IMAX) :: q_w
!> q_w_x(j,i): Scalar volume flux of the basal water in x-direction
   real(dp), dimension(0:JMAX,0:IMAX) :: q_w_x
!> q_w_y(j,i): Scalar volume flux of the basal water in y-direction
   real(dp), dimension(0:JMAX,0:IMAX) :: q_w_y
!> H_w(j,i): Thickness of the water column under the ice base
   real(dp), dimension(0:JMAX,0:IMAX) :: H_w
!> accum(j,i): Accumulation rate at the ice surface
!>             (including liquid precipitation = rainfall!)
   real(dp), dimension(0:JMAX,0:IMAX) :: accum
!> snowfall(j,i): Snowfall rate at the ice surface
   real(dp), dimension(0:JMAX,0:IMAX) :: snowfall
!> rainfall(j,i): Rainfall rate at the ice surface
   real(dp), dimension(0:JMAX,0:IMAX) :: rainfall
!> ET(j,i): Temperature excess at the ice surface
!>          (positive degree days divided by time)
   real(dp), dimension(0:JMAX,0:IMAX) :: ET
!> melt(j,i): Melting rate at the ice surface
   real(dp), dimension(0:JMAX,0:IMAX) :: melt
!> melt_star(j,i): Superimposed ice formation rate at the ice surface
   real(dp), dimension(0:JMAX,0:IMAX) :: melt_star
!> runoff(j,i): Runoff rate at the ice surface
   real(dp), dimension(0:JMAX,0:IMAX) :: runoff
!> as_perp(j,i): Accumulation-ablation function at the ice surface (SMB)
   real(dp), dimension(0:JMAX,0:IMAX) :: as_perp
!> as_perp_apl(j,i): Applied accumulation-ablation function (SMB)
   real(dp), dimension(0:JMAX,0:IMAX) :: as_perp_apl
!> smb_corr_prescribed(j,i): Prescribed SMB correction
   real(dp), dimension(0:JMAX,0:IMAX) :: smb_corr_prescribed
!> smb_corr(j,i): Diagnosed SMB correction
   real(dp), dimension(0:JMAX,0:IMAX) :: smb_corr
!> mb_source_apl(j,i): Applied mass balance source (SMB, BMB, calving)
   real(dp), dimension(0:JMAX,0:IMAX) :: mb_source_apl
!> accum_apl(j,i): Applied accumulation rate at the ice surface
!>                 (including liquid precipitation = rainfall!)
   real(dp), dimension(0:JMAX,0:IMAX) :: accum_apl
!> runoff_apl(j,i): Applied runoff rate at the ice surface
   real(dp), dimension(0:JMAX,0:IMAX) :: runoff_apl
!> Q_b_apl(j,i): Applied basal melting rate
!>                       + water drainage rate from the temperate layer
   real(dp), dimension(0:JMAX,0:IMAX) :: Q_b_apl
!> calving(j,i): Calving rate of grounded and floating ice
   real(dp), dimension(0:JMAX,0:IMAX) :: calving
!> calving_apl(j,i): Applied calving rate of grounded and floating ice
   real(dp), dimension(0:JMAX,0:IMAX) :: calving_apl
!> mask_ablation_type(j,i): Mask indicating ablation type.
!>             2: visible (ocean, for later developments),
!>             1: visible (grounded ice),
!>            -1: hidden on land,
!>            -2: hidden in ocean
   integer(i4b), dimension(0:JMAX,0:IMAX) :: mask_ablation_type
!> temp_maat(j,i): Mean annual air temperature
   real(dp), dimension(0:JMAX,0:IMAX) :: temp_maat
!> temp_s(j,i): Ice surface temperature
   real(dp), dimension(0:JMAX,0:IMAX) :: temp_s
!> z_sl(j,i): Sea level
   real(dp), dimension(0:JMAX,0:IMAX) :: z_sl
!> dzsl_dtau(j,i): Derivative of zsl by tau (time)
   real(dp), dimension(0:JMAX,0:IMAX) :: dzsl_dtau
!> z_sl_mean: Mean sea level
   real(dp) :: z_sl_mean
!> am_perp(j,i): Ice volume flux across the z=zm interface
   real(dp), dimension(0:JMAX,0:IMAX) :: am_perp
!> am_perp_st(j,i): Steady-state part of am_perp
!>                  (without contribution of dzm_dtau)
   real(dp), dimension(0:JMAX,0:IMAX) :: am_perp_st
!> (.)_new: New value of quantity (.) computed during an integration step
   real(dp), dimension(0:JMAX,0:IMAX) :: zs_new
!> (.)_new: New value of quantity (.) computed during an integration step
   real(dp), dimension(0:JMAX,0:IMAX) :: zm_new
!> (.)_new: New value of quantity (.) computed during an integration step
   real(dp), dimension(0:JMAX,0:IMAX) :: zb_new
!> (.)_new: New value of quantity (.) computed during an integration step
   real(dp), dimension(0:JMAX,0:IMAX) :: zl_new
!> (.)_new: New value of quantity (.) computed during an integration step
   real(dp), dimension(0:JMAX,0:IMAX) :: H_new
!> (.)_new: New value of quantity (.) computed during an integration step
   real(dp), dimension(0:JMAX,0:IMAX) :: H_c_new
!> (.)_new: New value of quantity (.) computed during an integration step
   real(dp), dimension(0:JMAX,0:IMAX) :: H_t_new

!> zs_ref(j,i): Reference elevation for the present-day climatology
   real(dp), dimension(0:JMAX,0:IMAX) :: zs_ref

!> accum_present(j,i): Present-day accumulation rate at the ice surface
!>                     (for EISMINT, ISMIP HEINO and the north and south
!>                     polar caps of Mars)
   real(dp), dimension(0:JMAX,0:IMAX) :: accum_present
!> precip_ma_present(j,i): Present-day mean annual precipitation rate
!>                         at the ice surface
   real(dp), dimension(0:JMAX,0:IMAX) :: precip_ma_present
!> precip_ma_lgm_anom(j,i): LGM anomaly (ratio LGM/present) of the mean annual
!>                          precipitation rate at the ice surface
   real(dp), dimension(0:JMAX,0:IMAX) :: precip_ma_lgm_anom
!> temp_ma_present(j,i): Present-day mean annual surface temperature
   real(dp), dimension(0:JMAX,0:IMAX) :: temp_ma_present
!> temp_mj_present(j,i): Present-day mean summer (northern hemisphere: July,
!>                       southern hemisphere: January) surface temperature
   real(dp), dimension(0:JMAX,0:IMAX) :: temp_mj_present
!> temp_ma_lgm_anom(j,i): LGM anomaly (difference LGM - present) of the mean
!>                        annual surface temperature
   real(dp), dimension(0:JMAX,0:IMAX) :: temp_ma_lgm_anom
!> temp_mj_lgm_anom(j,i): LGM anomaly (difference LGM - present) of the mean
!>                        summer (northern hemisphere: July, southern
!>                        hemisphere: January) surface temperature
   real(dp), dimension(0:JMAX,0:IMAX) :: temp_mj_lgm_anom

!> dist_dxdy(jr,ir): Distance between grid points with delta_i=ir, delta_j=jr
   real(dp), dimension(-JMAX:JMAX,-IMAX:IMAX) :: dist_dxdy

!> precip_present(j,i,n): Present-day mean monthly precipitation rate
!>                        at the ice surface
   real(dp), dimension(0:JMAX,0:IMAX,12) :: precip_present
!> precip_lgm_anom(j,i,n): LGM anomaly (ratio LGM/present) of the mean monthly
!>                         precipitation rate at the ice surface
   real(dp), dimension(0:JMAX,0:IMAX,12) :: precip_lgm_anom
!> gamma_precip_lgm_anom(j,i,n): negative natural logarithm of
!>                               precip_lgm_anom(j,i,n)
   real(dp), dimension(0:JMAX,0:IMAX,12) :: gamma_precip_lgm_anom
!> temp_mm_present(j,i,n): Present-day mean monthly surface temperature
   real(dp), dimension(0:JMAX,0:IMAX,12) :: temp_mm_present
!> temp_mm_lgm_anom(j,i,n): LGM anomaly (difference LGM - present) of the mean
!>                          monthly surface temperature
   real(dp), dimension(0:JMAX,0:IMAX,12) :: temp_mm_lgm_anom

!> d_help_c(kc,j,i): Auxiliary quantity for the computation of vx, vy und zs
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: d_help_c
!> vx_c(kc,j,i): Velocity in x-direction in the upper (kc) ice domain
!>               (at (i+1/2,j,kc))
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: vx_c
!> vy_c(kc,j,i): Velocity in y-direction in the upper (kc) ice domain
!>               (at (i,j+1/2,kc))
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: vy_c
!> vz_c(kc,j,i): Velocity in z-direction in the upper (kc) ice domain
!>               (at (i,j,kc+1/2))
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: vz_c
!> temp_c(kc,j,i): Temperature in the upper (kc) ice domain
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: temp_c
!> (.)_new: New value of quantity (.) computed during an integration step
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: temp_c_new
!> temp_c_m(kc,j,i): Melting temperature in the upper (kc) ice domain
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: temp_c_m
!> age_c(kc,j,i): Age in the upper (kc) ice domain
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: age_c
!> (.)_new: New value of quantity (.) computed during an integration step
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: age_c_new
!> txz_c(kc,j,i): Shear stress txz in the upper (kc) ice domain
!>                (at (i+1/2,j,kc))
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: txz_c
!> tyz_c(kc,j,i): Shear stress tyz in the upper (kc) ice domain
!>                (at (i,j+1/2,kc))
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: tyz_c
!> sigma_c(kc,j,i): Effective stress in the upper (kc) ice domain
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: sigma_c
!> enh_c(kc,j,i): Flow enhancement factor in the upper (kc) ice domain
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: enh_c
!> strain_heating_c(kc,j,i): Strain heating in the upper (kc) ice domain
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: strain_heating_c

!> de_ssa(j,i): Effective strain rate of the SSA, at (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: de_ssa
!> vis_ave_g(j,i): Depth-averaged viscosity of the SIA/SSA, at (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: vis_ave_g
!> vis_int_g(j,i): Depth-integrated viscosity of the SIA/SSA, at (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: vis_int_g
!> vx_g(j,i): Velocity in x-direction of the SSA, at (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: vx_g
!> vy_g(j,i): Velocity in y-direction of the SSA, at (i,j)
   real(dp), dimension(0:JMAX,0:IMAX) :: vy_g

!> d_help_t(kt,j,i): Auxiliary quantity for the computation of vx, vy und zs
   real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: d_help_t
!> vx_t(kt,j,i): Velocity in x-direction in the lower (kt) ice domain
!>               (at (i+1/2,j,kt))
   real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: vx_t
!> vy_t(kt,j,i): Velocity in y-direction in the lower (kt) ice domain
!>               (at (i,j+1/2,kt))
   real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: vy_t
!> vz_t(kt,j,i): Velocity in z-direction in the lower (kt) ice domain
!>               (at (i,j,kt+1/2))
   real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: vz_t
!> omega_t(kt,j,i): Water content in the lower (kt) ice domain
   real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: omega_t
!> (.)_new: New value of quantity (.) computed during an integration step
   real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: omega_t_new
!> temp_t_m(kt,j,i): Melting temperature in the lower (kt) ice domain
   real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: temp_t_m
!> age_t(kt,j,i): Age in the lower (kt) ice domain
   real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: age_t
!> (.)_new: New value of quantity (.) computed during an integration step
   real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: age_t_new
!> txz_t(kt,j,i): Shear stress txz in the lower (kt) ice domain
!>                (at (i+1/2,j,kt))
   real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: txz_t
!> tyz_t(kt,j,i): Shear stress tyz in the lower (kt) ice domain
!>                (at (i,j+1/2,kt))
   real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: tyz_t
!> sigma_t(kt,j,i): Effective stress in the lower (kt) ice domain
   real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: sigma_t
!> enh_t(kt,j,i): Flow enhancement factor in the lower (kt) ice domain
   real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: enh_t
!> strain_heating_t(kt,j,i): Strain heating in the lower (kt) ice domain
   real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: strain_heating_t

!> temp_r(kr,j,i): Temperature in the bedrock
   real(dp), dimension(0:KRMAX,0:JMAX,0:IMAX) :: temp_r
!> (.)_new: New value of quantity (.) computed during an integration step
   real(dp), dimension(0:KRMAX,0:JMAX,0:IMAX) :: temp_r_new

!> enth_c(kc,j,i): Enthalpy in the upper (kc) ice domain
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: enth_c
!> (.)_new: New value of quantity (.) computed during an integration step
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: enth_c_new
!> omega_c(kc,j,i): Water content in the upper (kc) ice domain
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: omega_c
!> (.)_new: New value of quantity (.) computed during an integration step
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: omega_c_new

!> enth_t(kt,j,i): Enthalpy in the lower (kt) ice domain
   real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: enth_t
!> (.)_new: New value of quantity (.) computed during an integration step
   real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: enth_t_new

!> dxx_c(kc,j,i): Strain rate dxx in the upper (kc) ice domain
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: dxx_c
!> dyy_c(kc,j,i): Strain rate dyy in the upper (kc) ice domain
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: dyy_c
!> dxy_c(kc,j,i): Strain rate dxy in the upper (kc) ice domain
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: dxy_c
!> dxz_c(kc,j,i): Strain rate dxz in the upper (kc) ice domain
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: dxz_c
!> dyz_c(kc,j,i): Strain rate dyz in the upper (kc) ice domain
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: dyz_c
!> de_c(kc,j,i): Full effective strain rate in the upper (kc) ice domain
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: de_c
!> lambda_shear_c(kc,j,i): Shear fraction in the upper (kc) ice domain
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: lambda_shear_c

!> dxx_t(kt,j,i): Strain rate dxx in the lower (kt) ice domain
   real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: dxx_t
!> dyy_t(kt,j,i): Strain rate dyy in the lower (kt) ice domain
   real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: dyy_t
!> dxy_t(kt,j,i): Strain rate dxy in the lower (kt) ice domain
   real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: dxy_t
!> dxz_t(kt,j,i): Strain rate dxz in the lower (kt) ice domain
   real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: dxz_t
!> dyz_t(kt,j,i): Strain rate dyz in the lower (kt) ice domain
   real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: dyz_t
!> de_t(kt,j,i): Full effective strain rate in the lower (kt) ice domain
   real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: de_t
!> lambda_shear_t(kt,j,i): Shear fraction in the lower (kt) ice domain
   real(dp), dimension(0:KTMAX,0:JMAX,0:IMAX) :: lambda_shear_t

!-------- Physical parameters --------

!> RHO: Density of ice
   real(dp) :: RHO
!> RHO_W: Density of pure water
   real(dp) :: RHO_W
!> RHO_SW: Density of sea water
   real(dp) :: RHO_SW
!> L: Latent heat of ice
   real(dp) :: L
!> G: Acceleration due to gravity
   real(dp) :: G
!> NUE: Water diffusivity in ice
   real(dp) :: NUE
!> BETA: Clausius-Clapeyron gradient of ice
   real(dp) :: BETA
!> DELTA_TM_SW: Melting point depression of sea water due to its
!>              average salinity
   real(dp) :: DELTA_TM_SW
!> OMEGA_MAX: Threshold value for the water content of temperate ice
   real(dp) :: OMEGA_MAX
!> H_R: Thickness of the modelled lithosphere layer
   real(dp) :: H_R
!> RHO_C_R: Density times specific heat of the lithosphere
   real(dp) :: RHO_C_R
!> KAPPA_R: Heat conductivity of the lithosphere
   real(dp) :: KAPPA_R
!> RHO_A: Density of the asthenosphere
   real(dp) :: RHO_A
!> R_T: Coefficient of the water-content dependence in the rate factor
!>      for temperate ice
   real(dp) :: R_T

!> R: Mean radius of the planet
   real(dp) :: R
!> A: Semi-major axis of the planet
   real(dp) :: A
!> B: Semi-minor axis of the planet
   real(dp) :: B
!> F_INV: Inverse flattening of the planet (>1e10 interpreted as infinity)
   real(dp) :: F_INV
!> LAMBDA0: Reference longitude (central meridian) of the stereographic
!>          projection
   real(dp) :: LAMBDA0
!> PHI0: Standard parallel of the stereographic projection
   real(dp) :: PHI0

!> RF(n): Tabulated values for the rate factor of cold ice
   real(dp), dimension(-190:10) :: RF
!> KAPPA(n): Tabulated values for the heat conductivity of ice
   real(dp), dimension(-190:10) :: KAPPA
!> C(n): Tabulated values for the specific heat of ice
   real(dp), dimension(-190:10) :: C

!-------- Mathematical constants -------- 

!> pi: Constant pi
   real(dp), parameter :: pi = 3.141592653589793_dp

!> deg2rad: pi divided by 180 (-> deg to rad)
   real(dp), parameter :: deg2rad = pi/180.0_dp
!> rad2deg: 180 divided by pi (-> rad to deg)
   real(dp), parameter :: rad2deg = 180.0_dp/pi

!> euler: Euler number
   real(dp), parameter :: euler = 2.718281828459045_dp

!> eps: Small number
   real(dp), parameter :: eps = 1.0e-05_dp
!> epsi: Very small number
   real(dp), parameter :: epsi = 1.0e-12_dp

#if (!defined(ALLOW_GRDCHK) && !defined(ALLOW_TAPENADE)) /* Normal */

!> eps_sp: Small number to single-precision accuracy
   real(sp), parameter :: eps_sp = epsilon(1.0_sp)
!> eps_sp_dp: Small number to single-precision accuracy in double precision
   real(dp), parameter :: eps_sp_dp = eps_sp
!> eps_dp: Small number to double-precision accuracy
   real(dp), parameter :: eps_dp = epsilon(1.0_dp)

#else /* Tapenade */

   ! myepsilon_sp was computed using the code below. 4 is the value of sp
   ! real(4) :: y = 1.0
   ! print *, EPSILON(y)
   real(sp), parameter :: myepsilon_sp  = 1.19209290E-07
   real(sp), parameter :: eps_sp = myepsilon_sp
   ! myepsilon_sp_dp was computed using the code below. 4 is the value of sp,
   ! 8 is the value of dp
   ! real(4) :: y = 1.0
   ! print *, REAL(EPSILON(y),8)
   real(dp), parameter :: myepsilon_sp_dp  = 1.1920928955078125E-007
   real(dp), parameter :: eps_sp_dp = myepsilon_sp_dp
   ! myepsilon_dp was computed using the code below. 8 is the value of dp
   ! real(8) :: y = 1.0
   ! print *, EPSILON(y)
   real(dp), parameter :: myepsilon_dp  = 2.2204460492503131E-016
   real(dp), parameter :: eps_dp = myepsilon_dp

#endif /* Normal vs. Tapenade */

!-------- Further quantities -------- 

!> year2sec: 1 year (1 a) in seconds
   real(dp) :: year2sec

!> sec2year: 1 second in years
   real(dp) :: sec2year

!> year_zero: SICOPOLIS year zero in astronomical year numbering
!>            [ = signed year CE (AD) ]
   real(dp) :: year_zero

!> run_name: Name of simulation
   character(len=256) :: run_name

!> ch_domain_long: Long name of the computational domain
   character(len=64) :: ch_domain_long
!> ch_domain_short: Short name of the computational domain
   character(len=16) :: ch_domain_short

!> forcing_flag: Flag for the forcing type.
!>               1: forcing by a spatially constant surface temperature
!>                  anomaly (delta_ts),
!>               2: forcing by a glacial index (glac_index),
!>               3: forcing by time-dependent surface temperature
!>                  and precipitation data.
   integer(i4b) :: forcing_flag

!> n_core: Number of positions to be considered in the time-series file
!>         for deep boreholes
   integer(i4b) :: n_core
!> n_core_max: Maximum allowed value of n_core
   integer(i4b), parameter :: n_core_max = 256
!> lambda_core(n): Geographical longitude of the prescribed borehole positions
   real(dp), dimension(n_core_max) :: lambda_core
!> phi_core(n): Geographical latitude of the prescribed borehole positions
   real(dp), dimension(n_core_max) :: phi_core
!> x_core(n): Coordinate xi (= x) of the prescribed borehole positions
   real(dp), dimension(n_core_max) :: x_core
!> y_core(n): Coordinate eta (= y) of the prescribed borehole positions
   real(dp), dimension(n_core_max) :: y_core
!> ch_core(n): Names of the prescribed borehole positions
   character(len=16), dimension(n_core_max) :: ch_core

!> grip_time_min: Minimum time of the data values for the
!>                surface temperature anomaly
   integer(i4b) :: grip_time_min
!> grip_time_stp: Time step of the data values for the
!>                surface temperature anomaly
   integer(i4b) :: grip_time_stp
!> grip_time_max: Maximum time of the data values for the
!>                surface temperature anomaly
   integer(i4b) :: grip_time_max
!> ndata_grip: Number of data values for the surface temperature anomaly
   integer(i4b) :: ndata_grip
!> ndata_grip_max: Maximum allowed value of ndata_grip
   integer(i4b), parameter :: ndata_grip_max = 262143
!> griptemp(n): Data values for the surface temperature anomaly
   real(dp), dimension(0:ndata_grip_max):: griptemp

!> gi_time_min: Minimum time of the data values for the glacial index
   integer(i4b) :: gi_time_min
!> gi_time_stp: Time step of the data values for the glacial index
   integer(i4b) :: gi_time_stp
!> gi_time_max: Maximum time of the data values for the glacial index
   integer(i4b) :: gi_time_max
!> ndata_gi: Number of data values for the glacial index
   integer(i4b) :: ndata_gi
!> ndata_gi_max: Maximum allowed value of ndata_gi
   integer(i4b), parameter :: ndata_gi_max = 262143
!> glacial_index(n): Data values for the glacial index
   real(dp), dimension(0:ndata_gi_max) :: glacial_index

!> specmap_time_min: Minimum time of the data values for the sea level
   integer(i4b) :: specmap_time_min
!> specmap_time_stp: Time step of the data values for the sea level
   integer(i4b) :: specmap_time_stp
!> specmap_time_max: Maximum time of the data values for the sea level
   integer(i4b) :: specmap_time_max
!> ndata_specmap: Number of data values for the sea level
   integer(i4b) :: ndata_specmap
!> ndata_specmap_max: Maximum allowed value of ndata_specmap
   integer(i4b), parameter :: ndata_specmap_max = 16383
!> specmap_zsl(n): Data values for the sea level
   real(dp), dimension(0:ndata_specmap_max) :: specmap_zsl

!> target_topo_tau0_time_min: Minimum time of the data values
!>                            for the relaxation time for the topography nudging
   integer(i4b) :: target_topo_tau0_time_min
!> target_topo_tau0_time_stp: Time step of the data values
!>                            for the relaxation time for the topography nudging
   integer(i4b) :: target_topo_tau0_time_stp
!> target_topo_tau0_time_max: Maximum time of the data values
!>                            for the relaxation time for the topography nudging
   integer(i4b) :: target_topo_tau0_time_max
!> ndata_target_topo_tau0: Number of data values
!>                         for the relaxation time for the topography nudging
   integer(i4b) :: ndata_target_topo_tau0
!> ndata_target_topo_tau0_max: Maximum allowed value of ndata_target_topo_tau0
   integer(i4b), parameter :: ndata_target_topo_tau0_max = 16383
!> target_topo_tau0(n): Data values for the relaxation time
!>                      for the topography nudging
   real(dp), dimension(0:ndata_target_topo_tau0_max) :: target_topo_tau0
!> target_topo_tau_0: Constant relaxation time for the topography nudging
   real(dp) :: target_topo_tau_0
!> target_topo_tau: Relaxation time for the topography nudging
   real(dp) :: target_topo_tau
!> mask_target(j,i): Target topography (ice-land-ocean mask)
   integer(i4b), dimension(0:JMAX,0:IMAX) :: mask_target
!> zs_target(j,i): Target topography (ice surface)
   real(dp), dimension(0:JMAX,0:IMAX) :: zs_target
!> zb_target(j,i): Target topography (ice base)
   real(dp), dimension(0:JMAX,0:IMAX) :: zb_target
!> zl_target(j,i): Target topography (lithosphere surface)
   real(dp), dimension(0:JMAX,0:IMAX) :: zl_target
!> H_target(j,i): Target topography (ice thickness)
   real(dp), dimension(0:JMAX,0:IMAX) :: H_target

!> flag_mask_maxextent: Flag for maximum ice extent.
!>                       .true.: specified
!>                      .false.: not specified
   logical :: flag_mask_maxextent
!> mask_maxextent(j,i): Maximum ice extent mask.
!>                       0: not allowed to glaciate,
!>                       1: allowed to glaciate.
   integer(i4b), dimension(0:JMAX,0:IMAX) :: mask_maxextent

!> ncid_ser: IDs of the NetCDF time-series output files
   integer(i4b), dimension(0:99) :: ncid_ser
!> ncid_core: ID of the NetCDF time-series output file for the deep ice cores
   integer(i4b) :: ncid_core

!> kei(n): Tabulated values of the kei function (Kelvin function of zero order)
   real(dp), dimension(-10000:10000) :: kei
!> n_data_kei: Number of tabulated values of the kei function
   integer(i4b):: n_data_kei
!> kei_r_max: Maximum value of the argument r of the tabulated kei function
   real(dp) :: kei_r_max
!> kei_r_incr: Increment of the argument r of the tabulated kei function
   real(dp) :: kei_r_incr

!> rcl1: Maximum length of record for input files
!>       (with factor 3 safety margin)
   integer(i4b), parameter :: rcl1 = 3*8*(IMAX+1)
!> rcl2: Maximum length of record for input mask files
!>       (with factor 3 safety margin)
   integer(i4b), parameter :: rcl2 = 3  *(IMAX+1)

!> ij2n: Conversion from 2d index pair (i,j) to linear index n
   integer(i4b), dimension(0:JMAX,0:IMAX) :: ij2n
!> n2i: Conversion from linear index n to 2d index i
   integer(i4b), dimension((IMAX+1)*(JMAX+1)) :: n2i
!> n2j: Conversion from linear index n to 2d index j
   integer(i4b), dimension((IMAX+1)*(JMAX+1)) :: n2j

!> no_value_pos_1: Positive no-value parameter
   real(dp), parameter :: no_value_pos_1 =  1.11e+11_dp
!> no_value_pos_2: Positive no-value parameter
   real(dp), parameter :: no_value_pos_2 =  9.999e+03_dp
!> no_value_neg_1: Negative no-value parameter
   real(dp), parameter :: no_value_neg_1 = -1.11e+11_dp
!> no_value_neg_2: Negative no-value parameter
   real(dp), parameter :: no_value_neg_2 = -9.999e+03_dp

!> errormsg: Error-message string
   character(len=256) :: errormsg
!> end_of_line: End-of-line string
   character, parameter :: end_of_line = char(10)

#if (defined(ALLOW_NORMAL) \
    || defined(ALLOW_GRDCHK) || defined(ALLOW_TAPENADE)) /* Tapenade */

!> fc: scalar cost function
   real(dp) :: fc

#if (defined(AGE_COST))

! Note: for the age cost, CALCMOD!=1 is recommended because
! the gridded ages of the GRL ice sheet are only 25 z-levels.

!> age_data: array of gridded ages to be used in adjoint mode
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: age_data
!> age_unc: array of gridded uncertainty in ages to be used in adjoint mode
   real(dp), dimension(0:KCMAX,0:JMAX,0:IMAX) :: age_unc
!> H_data: array of gridded ice thickness to be used in adjoint mode
   real(dp), dimension(0:JMAX,0:IMAX) :: H_data
!> H_unc: array of gridded uncertainty in ice thickness to be used
!>        in adjoint mode
   real(dp), dimension(0:JMAX,0:IMAX) :: H_unc

#endif /* No age cost used */

   real(dp), dimension(0:JMAX,0:IMAX) :: acc_fact

#if (defined(BEDMACHINE_COST)) 

   real(dp), dimension(0:JMAX,0:IMAX) :: H_BedMachine_data
   real(dp), dimension(0:JMAX,0:IMAX) :: H_unc_BedMachine_data

#endif

#endif /* Tapenade */

end module sico_variables_m
!
