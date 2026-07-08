!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  f r o n t a l _ a b l a t i o n _ m
!
!! Frontal ablation (frontal melting, calving).
!!
!!##### Authors
!!
!! Ralf Greve, Thorben Dunse, Nicolas Sartore
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
!> Frontal ablation (frontal melting, calving).
!-------------------------------------------------------------------------------
module frontal_ablation_m

  use sico_types_m
  use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

  use error_m

  implicit none

  public

contains

#if (FRONTAL_MELTING==1)
!-------------------------------------------------------------------------------
!> Frontal melting (at grounded, vertical fronts).
!-------------------------------------------------------------------------------
  subroutine frontal_melting_grounded(my_zb, my_z_sl, dxi, deta)

  implicit none

  real(dp), dimension(0:JMAX,0:IMAX), intent(in) :: my_zb, my_z_sl
  real(dp),                           intent(in) :: dxi, deta

  integer(i4b) :: i, j, ij, n
  real(dp)     :: a_fm, b_fm, alpha_fm, beta_fm
  real(dp)     :: lambda_a_fm, lambda_b_fm

  real(dp), dimension(0:JMAX,0:IMAX) :: H_submerged
  real(dp), dimension(0:JMAX,0:IMAX) :: frontal_area_submerged
  real(dp), dimension(0:JMAX,0:IMAX) :: sgd_normalized

  real(dp), dimension(N_SGD_REGIONS) :: frontal_area_submerged_region

  real(dp) :: frontal_melting_horizontal, area_ratio

  character(len=8) :: ch_i, ch_j, ch_n

  a_fm     = 3.0e-04_dp
  b_fm     = 0.15_dp
  alpha_fm = 0.39_dp
  beta_fm  = 1.18_dp

#if (defined(LAMBDA_A_FRONT_MELT))
  lambda_a_fm = real(LAMBDA_A_FRONT_MELT,dp)
#else
  lambda_a_fm = 1.0_dp
#endif

#if (defined(LAMBDA_B_FRONT_MELT))
  lambda_b_fm = real(LAMBDA_B_FRONT_MELT,dp)
#else
  lambda_b_fm = 1.0_dp
#endif

  a_fm = a_fm * lambda_a_fm   ! scaling of the coefficients
  b_fm = b_fm * lambda_b_fm   ! of the parameterization

  do ij=1, (IMAX+1)*(JMAX+1)

     i = n2i(ij)   ! i=0...IMAX
     j = n2j(ij)   ! j=0...JMAX

     H_submerged(j,i) = my_z_sl(j,i)-my_zb(j,i)
                        ! submerged ice depth (if positive)

     frontal_melting(j,i)        = 0.0_dp
     frontal_area_submerged(j,i) = 0.0_dp
     sgd_normalized(j,i)         = 0.0_dp

  end do

  frontal_area_submerged_region = 0.0_dp

  do ij=1, (IMAX+1)*(JMAX+1)

     i = n2i(ij)   ! i=0...IMAX
     j = n2j(ij)   ! j=0...JMAX

     if ( flag_inner_point(j,i) &
          .and. &
          flag_grounded_front_b_1(j,i) &
          .and. &
          H_submerged(j,i) > 0.0_dp ) then
                             ! inner point, marine-terminating grounded front

        if (flag_grounded_front_b_2(j,i+1)) then
           frontal_area_submerged(j,i) &
                = frontal_area_submerged(j,i) &
                       + H_submerged(j,i)*(deta*sq_g22_sgx(j,i))
        end if

        if (flag_grounded_front_b_2(j,i-1)) then
           frontal_area_submerged(j,i) &
                = frontal_area_submerged(j,i) &
                       + H_submerged(j,i)*(deta*sq_g22_sgx(j,i-1))
        end if

        if (flag_grounded_front_b_2(j+1,i)) then
           frontal_area_submerged(j,i) &
                = frontal_area_submerged(j,i) &
                       + H_submerged(j,i)*(dxi*sq_g11_sgy(j,i))
        end if

        if (flag_grounded_front_b_2(j-1,i)) then
           frontal_area_submerged(j,i) &
                = frontal_area_submerged(j,i) &
                       + H_submerged(j,i)*(dxi*sq_g11_sgy(j-1,i))
        end if

        if (frontal_area_submerged(j,i) < eps_dp) then
           write (ch_i, '(i0)') i; ch_i = adjustl(ch_i)
           write (ch_j, '(i0)') j; ch_j = adjustl(ch_j)
           errormsg = ' >>> frontal_melting_grounded: ' &
                    //         end_of_line &
                    //'        Non-zero area ''frontal_area_submerged(j,i)''' &
                    //         end_of_line &
                    //'        could not be determined for (i,j) =' &
                    //       ' ('//trim(ch_i)//','//trim(ch_j)//')!'
           call error(errormsg)
        end if

        n = n_sgd_region(j,i)

        frontal_area_submerged_region(n) = frontal_area_submerged_region(n) &
                                              + frontal_area_submerged(j,i)

     end if

  end do

  do ij=1, (IMAX+1)*(JMAX+1)

     i = n2i(ij)   ! i=0...IMAX
     j = n2j(ij)   ! j=0...JMAX

     if ( flag_inner_point(j,i) &
          .and. &
          flag_grounded_front_b_1(j,i) &
          .and. &
          H_submerged(j,i) > 0.0_dp ) then
                             ! inner point, marine-terminating grounded front

        n = n_sgd_region(j,i)

        if (frontal_area_submerged_region(n) >= eps_dp) then

           sgd_normalized(j,i) = (sgd(j,i)/frontal_area_submerged_region(n)) &
                                 * day2sec   ! m/s -> m/d

        else

           write (ch_n, '(i0)') n; ch_n = adjustl(ch_n)
           errormsg = ' >>> frontal_melting_grounded: ' &
                    //         end_of_line &
                    //'        Non-zero area' &
                    //       ' ''frontal_area_submerged_region(n)''' &
                    //         end_of_line &
                    //'        could not be determined for n = ' &
                    //         trim(ch_n) // ' !'
           call error(errormsg)

        end if

        frontal_melting_horizontal &
             = ( a_fm * H_submerged(j,i) &
                      * sgd_normalized(j,i)**alpha_fm + b_fm ) &
               * tf(j,i)**beta_fm   ! m/d

        area_ratio = frontal_area_submerged(j,i)/cell_area(j,i)

        frontal_melting(j,i) &
             = (frontal_melting_horizontal*area_ratio) &
               * sec2day  !    m/d = m3/d/(m2 vertical area)
                          ! -> m/s = m3/s/(m2 horizontal area)
 
     end if

  end do

  end subroutine frontal_melting_grounded

#endif   /* (FRONTAL_MELTING==1) */

!-------------------------------------------------------------------------------
!> Detection of detached icebergs (floating ice), which are flagged
!! and interpreted as unresolved calving events to be corrected.
!-------------------------------------------------------------------------------
  subroutine detect_icebergs(my_mask)

  implicit none

  integer(i4b), dimension(0:JMAX,0:IMAX), intent(in) :: my_mask

  ! Local arrays for tracking
  logical :: visited(0:JMAX,0:IMAX)

  ! Queue for breadth-first search (BFS) flood-fill
  ! Max possible size in worst case is the total number of grid points
  integer(i4b) :: queue_i((IMAX+1)*(JMAX+1)), queue_j((IMAX+1)*(JMAX+1))
  integer(i4b) :: q_head, q_tail

  ! Temporary component storage
  integer(i4b) :: comp_i((IMAX+1)*(JMAX+1)), comp_j((IMAX+1)*(JMAX+1))
  integer(i4b) :: comp_count

  ! Loop and neighbour indices
  integer(i4b) :: i, j, ij, ni, nj, k
  integer(i4b) :: cur_i, cur_j
  logical :: touches_grounded

  ! Direction offsets for 4-way connectivity (up, down, left, right)
  integer(i4b), dimension(4), parameter :: di = (/ 1, -1,  0,  0 /)
  integer(i4b), dimension(4), parameter :: dj = (/ 0,  0,  1, -1 /)

  ! Initialize outputs and tracking arrays
  flag_iceberg = .false.
  visited      = .false.

  ! Loop through the entire grid to find unvisited floating ice
  do ij=1, (IMAX+1)*(JMAX+1)

     i = n2i(ij)   ! i=0...IMAX
     j = n2j(ij)   ! j=0...JMAX

     if (my_mask(j,i)==3 .and. .not. visited(j,i)) then

        ! Start a new flood-fill component
        comp_count = 0
        touches_grounded = .false.

        ! Initialize queue
        q_head = 1
        q_tail = 1
        queue_i(q_tail) = i
        queue_j(q_tail) = j
        visited(j,i)    = .true.

        ! BFS loop
        do while (q_head <= q_tail)
           ! Pop from queue
           cur_i = queue_i(q_head)
           cur_j = queue_j(q_head)
           q_head = q_head + 1

           ! Add to the current ice component list
           comp_count = comp_count + 1
           comp_i(comp_count) = cur_i
           comp_j(comp_count) = cur_j

           ! Check 4 neighbours
           do k=1, 4

              ni = cur_i + di(k)
              nj = cur_j + dj(k)

              ! Ensure neighbour is within domain boundaries
              if (ni >= 0 .and. ni <= IMAX .and. nj >= 0 .and. nj <= JMAX) then

                 if (my_mask(nj,ni) == 0) then
                    ! If it hits grounded ice, this whole cluster is a shelf,
                    ! not an iceberg
                    touches_grounded = .true.
                 end if

                 if (my_mask(nj,ni) == 3 .and. .not. visited(nj,ni)) then
                    ! Queue the connected floating ice point
                    q_tail = q_tail + 1
                    queue_i(q_tail) = ni
                    queue_j(q_tail) = nj
                    visited(nj,ni)  = .true.
                 end if

              end if

           end do

        end do

        ! If the components never touched grounded ice, it's a detached iceberg
        if (.not. touches_grounded) then
           do k=1, comp_count
              flag_iceberg(comp_j(k), comp_i(k)) = .true.
           end do
        end if

     end if

  end do

  end subroutine detect_icebergs

!-------------------------------------------------------------------------------
!> Calving of grounded "underwater ice".
!-------------------------------------------------------------------------------
  subroutine calving_underwater_ice()

  implicit none

  integer(i4b) :: i, j, ij
  real(dp)     :: rhosw_rho_ratio
  real(dp)     :: calv_uw_coeff, r1_calv_uw, r2_calv_uw
  real(dp)     :: H0_flt

  real(dp), dimension(0:JMAX,0:IMAX) :: H_water, calv_uw_ice

!-------- Term abbreviations --------

  rhosw_rho_ratio = RHO_SW/RHO

!-------- Setting of parameters --------

#if (defined(CALV_UW_COEFF))
  calv_uw_coeff = CALV_UW_COEFF *sec2year
#else
  errormsg = ' >>> calving_underwater_ice: CALV_UW_COEFF undefined!'
  call error(errormsg)
#endif

#if (defined(R1_CALV_UW))
  r1_calv_uw = R1_CALV_UW
#else
  errormsg = ' >>> calving_underwater_ice: R1_CALV_UW undefined!'
  call error(errormsg)
#endif

#if (defined(R2_CALV_UW))
  r2_calv_uw = R2_CALV_UW
#else
  errormsg = ' >>> calving_underwater_ice: R2_CALV_UW undefined!'
  call error(errormsg)
#endif

#if (defined(H0_FLOAT))
  H0_flt = H0_FLOAT
#else
  H0_flt = 0.0_dp
#endif

!-------- Calving of "underwater ice" --------

  do ij=1, (IMAX+1)*(JMAX+1)

     i = n2i(ij)   ! i=0...IMAX
     j = n2j(ij)   ! j=0...JMAX

     calv_uw_ice(j,i) = 0.0_dp

     H_water(j,i) = max(z_sl(j,i)-zl(j,i), 0.0_dp)   ! water depth

     if ( (mask(j,i) == 0) &
          .and. (H(j,i) < rhosw_rho_ratio*H_water(j,i)+H0_flt) ) then
        calv_uw_ice(j,i) = calv_uw_coeff &
                           * H(j,i)**r1_calv_uw * H_water(j,i)**r2_calv_uw
     end if

     calving(j,i) = calving(j,i) + calv_uw_ice(j,i)

  end do

  end subroutine calving_underwater_ice

#if (ICE_SHELF_COLLAPSE_MASK==1)
!-------------------------------------------------------------------------------
!> Adjustment of the newly computed ice thickness distribution
!! due to the ice-shelf collapse mask (counted as calving).
!-------------------------------------------------------------------------------
  subroutine calving_retreat_mask(time, dtime, i, j)

  implicit none

  integer(i4b), intent(in) :: i, j
  real(dp)    , intent(in) :: time, dtime

  real(dp) :: H_new_tmp, dHdt_retreat
  real(dp) :: calv_retreat_mask
  real(dp) :: dtime_inv
  real(dp) :: dtime_1year, dtime_1year_inv

  dtime_inv       = 1.0_dp/dtime
  dtime_1year     = year2sec   ! 1 year (in seconds)
  dtime_1year_inv = 1.0_dp/dtime_1year

!-------- Saving computed H_new before any adjustments --------

  H_new_tmp = H_new(j,i)

!-------- Adjustment due to the retreat mask --------

  dHdt_retreat = 0.0_dp   ! initialization

  if ((H_new(j,i) > 0.0_dp).and.(mask(j,i)==3)) then

     dHdt_retreat = -(1.0_dp-r_mask_retreat(j,i))*H_ref_retreat(j,i) &
                                                 *dtime_1year_inv

     H_new(j,i) = max((H_new(j,i) + dHdt_retreat*dtime), 0.0_dp)

  end if

!-------- Computation of the mass balance adjustment --------

  calv_retreat_mask = (H_new_tmp-H_new(j,i))*dtime_inv
                      ! calving is counted as positive for mass loss

  calving(j,i) = calving(j,i) + calv_retreat_mask

  end subroutine calving_retreat_mask

#endif   /* (ICE_SHELF_COLLAPSE_MASK==1) */

!-------------------------------------------------------------------------------

end module frontal_ablation_m
!
