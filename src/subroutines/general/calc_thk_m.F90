!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ t h k _ m
!
!! Computation of the ice thickness.
!!
!!##### Authors
!!
!! Ralf Greve, Reinhard Calov, Tatsuru Sato
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
!> Computation of the ice thickness.
!-------------------------------------------------------------------------------
module calc_thk_m

  use sico_types_m
  use sico_variables_m
  use sico_vars_m
  use error_m

  implicit none

  real(dp), parameter :: eps_H = eps_sp_dp
     !! Small number for the ice-thickness computations
           ! If computations are in double precision,
           ! single precision is a good choice for this value

  public

contains

!-------------------------------------------------------------------------------
!> Initialisations for the ice thickness computation.
!-------------------------------------------------------------------------------
subroutine calc_thk_init()

implicit none

#if defined(ALLOW_TAPENADE) /* Tapenade */
integer(i4b) :: i, j
#endif /* Tapenade */

#if (defined(CALCZS))
  errormsg = ' >>> calc_thk_init: Replace CALCZS by CALCTHK in header file!'
  call error(errormsg)
#elif (!defined(CALCTHK))
  errormsg = ' >>> calc_thk_init: Define CALCTHK in header file!'
  call error(errormsg)
#endif

!-------- Computation/initialisation of the ice base topography
!                                       and its time derivative --------

#if (MARGIN==1 || MARGIN==2)   /* only grounded ice */

zb_new   = zl_new
dzb_dtau = dzl_dtau

#elif (MARGIN==3)   /* grounded and floating ice */

#if !defined(ALLOW_TAPENADE) /* Normal */

where (mask <= 1)   ! grounded ice or ice-free land
   zb_new   = zl_new
   dzb_dtau = dzl_dtau
elsewhere   ! (mask >= 2; ocean or floating ice)
   zb_new   = zb       ! initialisation,
   dzb_dtau = 0.0_dp   ! will be overwritten later
end where

#else /* Tapenade */

do i=0, IMAX
do j=0, JMAX
   if (mask(j,i) <= 1) then    ! grounded ice or ice-free land
      zb_new(j,i)   = zl_new(j,i)
      dzb_dtau(j,i) = dzl_dtau(j,i)
   else
      zb_new(j,i)   = zb(j,i)       ! initialisation,
      dzb_dtau(j,i) = 0.0_dp        ! will be overwritten later
   endif
end do
end do

#endif /* Normal vs. Tapenade */

#else

errormsg = ' >>> calc_thk_init: MARGIN must be either 1, 2 or 3!'
call error(errormsg)

#endif

!-------- Initialisation of the ice thickness
!                               and surface topography --------

zs_new = zs   ! initialisation,
H_new  = H    ! will be overwritten later

!-------- Solver type --------

#if (CALCTHK==1 || CALCTHK==4)   /* explicit solver */

  flag_thk_solver_explicit = .true.

#elif (CALCTHK==2)   /* implicit solver */

  flag_thk_solver_explicit = .false.

#else

  errormsg = ' >>> calc_thk_init: CALCTHK must be 1, 2 or 4!'
  call error(errormsg)

#endif

!-------- Source term for the ice thickness equation --------

mb_source = as_perp - Q_b_tot - calving

end subroutine calc_thk_init

!-------------------------------------------------------------------------------
!> Explicit solver for the diffusive SIA ice surface equation.
!-------------------------------------------------------------------------------
subroutine calc_thk_sia_expl(time, dtime, dxi, deta, z_mar)

#if (RETREAT_MASK==1 || ICE_SHELF_COLLAPSE_MASK==1)
  use calving_m
#endif

implicit none

real(dp), intent(in) :: time, dtime, dxi, deta
real(dp), intent(in) :: z_mar

integer(i4b)                       :: i, j
real(dp)                           :: azs2, azs3
real(dp), dimension(0:JMAX,0:IMAX) :: czs2, czs3

!-------- Abbreviations --------

azs2 = dtime/(dxi*dxi)
azs3 = dtime/(deta*deta)

czs2 = 0.0_dp
czs3 = 0.0_dp

do i=0, IMAX-1
do j=0, JMAX
   czs2(j,i) = azs2*0.5_dp*(h_diff(j,i)+h_diff(j,i+1)) &
               *(sq_g22_sgx(j,i)*insq_g11_sgx(j,i))
end do
end do

do i=0, IMAX
do j=0, JMAX-1
   czs3(j,i) = azs3*0.5_dp*(h_diff(j,i)+h_diff(j+1,i)) &
               *(sq_g11_sgy(j,i)*insq_g22_sgy(j,i))
end do
end do

!-------- Solution of the explicit scheme --------

do i=0, IMAX
do j=0, JMAX

   if (flag_inner_point(j,i)) then   ! inner point

      zs_new(j,i) = zs(j,i) &
                    + dtime*(mb_source(j,i)+dzb_dtau(j,i)) &
                    + ( ( czs2(j,i)  *(zs(j,i+1)-zs(j,i)  ) &
                         -czs2(j,i-1)*(zs(j,i)  -zs(j,i-1)) ) &
                       +( czs3(j,i)  *(zs(j+1,i)-zs(j,i)  ) &
                         -czs3(j-1,i)*(zs(j,i)  -zs(j-1,i)) ) ) &
                      *(insq_g11_g(j,i)*insq_g22_g(j,i))

   else
      zs_new(j,i) = zb_new(j,i)   ! zero-thickness boundary condition
   end if

end do
end do

!-------- Ice thickness --------

H_new = zs_new - zb_new

!-------- Applying the source term --------

call apply_mb_source(dtime, z_mar)

!-------- Adjusting the ice thickness, if needed --------

call thk_adjust(time, dtime)

!-------- Calving due to prescribed retreat mask --------

#if (RETREAT_MASK==1 || ICE_SHELF_COLLAPSE_MASK==1)
call calving_retreat_mask(time, dtime)
#endif

end subroutine calc_thk_sia_expl

!-------------------------------------------------------------------------------
!> Over-implicit solver for the diffusive SIA ice surface equation.
!-------------------------------------------------------------------------------
subroutine calc_thk_sia_impl(time, dtime, dxi, deta, z_mar, mean_accum)

use sico_maths_m

#if (RETREAT_MASK==1 || ICE_SHELF_COLLAPSE_MASK==1)
  use calving_m
#endif

implicit none

real(dp), intent(in) :: time, dtime, dxi, deta
real(dp), intent(in) :: z_mar
real(dp), intent(in) :: mean_accum

integer(i4b)                       :: i, j
integer(i4b)                       :: k, nnz
real(dp)                           :: azs2, azs3
real(dp), dimension(0:JMAX,0:IMAX) :: czs2, czs3

integer(i4b)                            :: ierr
integer(i4b)                            :: iter
integer(i4b)                            :: nc, nr
integer(i4b), parameter                 :: nmax   =    (IMAX+1)*(JMAX+1)
integer(i4b), parameter                 :: n_sprs = 10*(IMAX+1)*(JMAX+1)

#if !defined(ALLOW_TAPENADE) /* Normal */
integer(i4b), allocatable, dimension(:) :: lgs_a_ptr, lgs_a_index
integer(i4b), allocatable, dimension(:) :: lgs_a_diag_index
integer(i4b), allocatable, dimension(:) :: lgs_a_index_trim
real(dp),     allocatable, dimension(:) :: lgs_a_value, lgs_b_value, lgs_x_value
real(dp),     allocatable, dimension(:) :: lgs_a_value_trim
#else /* Tapenade */
real(dp),             dimension(n_sprs) :: lgs_a_index
integer(i4b),         dimension(nmax+1) :: lgs_a_ptr
integer(i4b),           dimension(nmax) :: lgs_a_diag_index
integer(i4b),         dimension(n_sprs) :: lgs_a_index_trim
real(dp),             dimension(n_sprs) :: lgs_a_value
real(dp),               dimension(nmax) :: lgs_b_value
real(dp),               dimension(nmax) :: lgs_x_value
real(dp),             dimension(n_sprs) :: lgs_a_value_trim
#endif /* Normal vs. Tapenade */

real(dp)                                :: eps_sor

!-------- Abbreviations --------

azs2 = dtime/(dxi*dxi)
azs3 = dtime/(deta*deta)

czs2 = 0.0_dp
czs3 = 0.0_dp

do i=0, IMAX-1
do j=0, JMAX
  czs2(j,i) = azs2*0.5_dp*(h_diff(j,i)+h_diff(j,i+1)) &
              *(sq_g22_sgx(j,i)*insq_g11_sgx(j,i))
end do
end do

do i=0, IMAX
do j=0, JMAX-1
   czs3(j,i) = azs3*0.5_dp*(h_diff(j,i)+h_diff(j+1,i)) &
               *(sq_g11_sgy(j,i)*insq_g22_sgy(j,i))
end do
end do

!-------- Assembly of the system of linear equations
!                     (matrix storage: compressed sparse row CSR) --------

#if !defined(ALLOW_TAPENADE) /* Normal */
allocate(lgs_a_value(n_sprs), lgs_a_index(n_sprs), lgs_a_ptr(nmax+1))
allocate(lgs_a_diag_index(nmax), lgs_b_value(nmax), lgs_x_value(nmax))
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

do nr=1, nmax   ! loop over rows

   i = n2i(nr)
   j = n2j(nr)

   if (flag_inner_point(j,i)) then   ! inner point

      nc = ij2n(j,i-1)   ! smallest nc (column counter), for zs(j,i-1)
      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_thk_sia_impl: n_sprs too small!'
      lgs_a_value(k) = -czs2(j,i-1)*OVI_WEIGHT &
                        *(insq_g11_g(j,i)*insq_g22_g(j,i))
      lgs_a_index(k) = nc

      nc = ij2n(j-1,i)   ! next nc (column counter), for zs(j-1,i)
      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_thk_sia_impl: n_sprs too small!'
      lgs_a_value(k) = -czs3(j-1,i)*OVI_WEIGHT &
                        *(insq_g11_g(j,i)*insq_g22_g(j,i))
      lgs_a_index(k) = nc

      nc = ij2n(j,i)     ! next nc (column counter), for zs(j,i)
      ! if (nc /= nr) &                     (diagonal element)
      !    stop ' >>> calc_thk_sia_impl: Check for diagonal element failed!'
      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_thk_sia_impl: n_sprs too small!'
      lgs_a_value(k) = 1.0_dp &
                       + ((czs2(j,i)+czs2(j,i-1))+(czs3(j,i)+czs3(j-1,i))) &
                         *OVI_WEIGHT &
                         *(insq_g11_g(j,i)*insq_g22_g(j,i))
      lgs_a_diag_index(nr) = k
      lgs_a_index(k) = nc

      nc = ij2n(j+1,i)   ! next nc (column counter), for zs(j+1,i)
      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_thk_sia_impl: n_sprs too small!'
      lgs_a_value(k) = -czs3(j,i)*OVI_WEIGHT &
                        *(insq_g11_g(j,i)*insq_g22_g(j,i))
      lgs_a_index(k) = nc

      nc = ij2n(j,i+1)   ! largest nc (column counter), for zs(j,i+1)
      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_thk_sia_impl: n_sprs too small!'
      lgs_a_value(k) = -czs2(j,i)*OVI_WEIGHT &
                        *(insq_g11_g(j,i)*insq_g22_g(j,i))
      lgs_a_index(k) = nc

      lgs_b_value(nr) = zs(j,i) &
                          + dtime*(mb_source(j,i)+dzb_dtau(j,i)) &
                          + ( ( czs2(j,i)*(zs(j,i+1)-zs(j,i)) &
                               -czs2(j,i-1)*(zs(j,i)-zs(j,i-1)) ) &
                             +( czs3(j,i)*(zs(j+1,i)-zs(j,i)) &
                               -czs3(j-1,i)*(zs(j,i)-zs(j-1,i)) ) ) &
                            *(1.0_dp-OVI_WEIGHT) &
                            *(insq_g11_g(j,i)*insq_g22_g(j,i))
                                                          ! right-hand side

   else   !  zero-thickness boundary condition

      k = k+1
      ! if (k > n_sprs) stop ' >>> calc_thk_sia_impl: n_sprs too small!'
      lgs_a_value(k)       = 1.0_dp   ! diagonal element only
      lgs_a_diag_index(nr) = k
      lgs_a_index(k)       = nr
      lgs_b_value(nr)      = zb_new(j,i)

   end if

   lgs_x_value(nr) = zs(j,i)   ! old surface topography,
                               ! initial guess for solution vector

   lgs_a_ptr(nr+1) = k+1   ! row is completed, store index to next row

end do

nnz = k   ! number of non-zero elements of the matrix

!-------- Solution of the system of linear equations
!                                   (built-in SOR solver) --------

#if !defined(ALLOW_TAPENADE) /* Normal */
allocate(lgs_a_value_trim(nnz), lgs_a_index_trim(nnz))
#endif /* Normal */

do k=1, nnz   ! relocate matrix to trimmed arrays
   lgs_a_value_trim(k) = lgs_a_value(k)
#if !defined(ALLOW_TAPENADE) /* Normal */
   lgs_a_index_trim(k) = lgs_a_index(k)
#else /* Tapenade */
   lgs_a_index_trim(k) = int(lgs_a_index(k))
#endif /* Normal vs. Tapenade */
end do

#if !defined(ALLOW_TAPENADE) /* Normal */
deallocate(lgs_a_value, lgs_a_index)
#endif /* Normal */

eps_sor = 1.0e-05_dp*mean_accum*dtime   ! convergence parameter

call sor_sprs(lgs_a_value_trim, &
              lgs_a_index_trim, lgs_a_diag_index, lgs_a_ptr, &
              lgs_b_value, &
              nnz, nmax, &
              OMEGA_SOR, eps_sor, lgs_x_value, ierr)

do nr=1, nmax
   i = n2i(nr)
   j = n2j(nr)
   zs_new(j,i) = lgs_x_value(nr)
end do

#if !defined(ALLOW_TAPENADE) /* Normal */
deallocate(lgs_a_value_trim, lgs_a_index_trim, lgs_a_ptr)
deallocate(lgs_a_diag_index, lgs_b_value, lgs_x_value)
#endif /* Normal */

!-------- Ice thickness --------

H_new = zs_new - zb_new

!-------- Applying the source term --------

call apply_mb_source(dtime, z_mar)

!-------- Adjusting the ice thickness, if needed --------

call thk_adjust(time, dtime)

!-------- Calving due to prescribed retreat mask --------

#if (RETREAT_MASK==1 || ICE_SHELF_COLLAPSE_MASK==1)
call calving_retreat_mask(time, dtime)
#endif

end subroutine calc_thk_sia_impl

!-------------------------------------------------------------------------------
!> Explicit solver for the general ice thickness equation.
!-------------------------------------------------------------------------------
subroutine calc_thk_expl(time, dtime, dxi, deta, z_mar)

#if (RETREAT_MASK==1 || ICE_SHELF_COLLAPSE_MASK==1)
  use calving_m
#endif

implicit none

real(dp), intent(in) :: time, dtime, dxi, deta
real(dp), intent(in) :: z_mar

integer(i4b)                       :: i, j
real(dp), dimension(0:JMAX,0:IMAX) :: dt_darea
real(dp), dimension(0:JMAX,0:IMAX) :: vx_m_1, vx_m_2, vy_m_1, vy_m_2
real(dp), dimension(0:JMAX,0:IMAX) :: upH_x_1, upH_x_2, upH_y_1, upH_y_2
real(dp), dimension(0:JMAX,0:IMAX) :: sq_g22_x_1, sq_g22_x_2
real(dp), dimension(0:JMAX,0:IMAX) :: sq_g11_y_1, sq_g11_y_2

!-------- Abbreviations --------

dt_darea = dtime/cell_area

do i=0, IMAX
do j=0, JMAX

   if (flag_inner_point(j,i)) then

      vx_m_1(j,i) = vx_m(j,i-1)
      vx_m_2(j,i) = vx_m(j,i)
      vy_m_1(j,i) = vy_m(j-1,i)
      vy_m_2(j,i) = vy_m(j,i)

      if (vx_m_1(j,i) >= 0.0_dp) then
         upH_x_1(j,i) = H(j,i-1)
      else
         upH_x_1(j,i) = H(j,i)
      end if

      if (vx_m_2(j,i) >= 0.0_dp) then
         upH_x_2(j,i) = H(j,i)
      else
         upH_x_2(j,i) = H(j,i+1)
      end if

      if (vy_m_1(j,i) >= 0.0_dp) then
         upH_y_1(j,i) = H(j-1,i)
      else
         upH_y_1(j,i) = H(j,i)
      end if

      if (vy_m_2(j,i) >= 0.0_dp) then
         upH_y_2(j,i) = H(j,i)
      else
         upH_y_2(j,i) = H(j+1,i)
      end if

      sq_g22_x_1(j,i) = sq_g22_sgx(j,i-1)
      sq_g22_x_2(j,i) = sq_g22_sgx(j,i)
      sq_g11_y_1(j,i) = sq_g11_sgy(j-1,i)
      sq_g11_y_2(j,i) = sq_g11_sgy(j,i)

   else   ! .not.(flag_inner_point(j,i))

      vx_m_1(j,i) = 0.0_dp
      vx_m_2(j,i) = 0.0_dp
      vy_m_1(j,i) = 0.0_dp
      vy_m_2(j,i) = 0.0_dp

      upH_x_1(j,i) = 0.0_dp
      upH_x_2(j,i) = 0.0_dp
      upH_y_1(j,i) = 0.0_dp
      upH_y_2(j,i) = 0.0_dp

      sq_g22_x_1(j,i) = 0.0_dp
      sq_g22_x_2(j,i) = 0.0_dp
      sq_g11_y_1(j,i) = 0.0_dp
      sq_g11_y_2(j,i) = 0.0_dp

   end if

end do
end do

!-------- Solution of the explicit scheme --------

#if !defined(ALLOW_TAPENADE) /* Normal */

where (flag_inner_point)   ! inner point

   H_new = H + dtime*mb_source &
             - dt_darea &
               * (  ( vx_m_2*upH_x_2*sq_g22_x_2*deta   &
                     -vx_m_1*upH_x_1*sq_g22_x_1*deta ) &
                  + ( vy_m_2*upH_y_2*sq_g11_y_2*dxi    &
                     -vy_m_1*upH_y_1*sq_g11_y_1*dxi  ) )

elsewhere
   H_new = 0.0_dp   ! zero-thickness boundary condition
end where

#else /* Tapenade */

do i=0, IMAX
do j=0, JMAX

   if (flag_inner_point(j,i)) then   ! inner point

      H_new(j,i) = H(j,i) + dtime*mb_source(j,i) &
                - dt_darea(j,i) &
                  * (  ( vx_m_2(j,i)*upH_x_2(j,i)*sq_g22_x_2(j,i)*deta   &
                        -vx_m_1(j,i)*upH_x_1(j,i)*sq_g22_x_1(j,i)*deta ) &
                     + ( vy_m_2(j,i)*upH_y_2(j,i)*sq_g11_y_2(j,i)*dxi    &
                        -vy_m_1(j,i)*upH_y_1(j,i)*sq_g11_y_1(j,i)*dxi  ) )

   else
      H_new(j,i) = 0.0_dp   ! zero-thickness boundary condition
   end if

end do
end do

#endif /* Normal vs. Tapenade */

!-------- Applying the source term --------

call apply_mb_source(dtime, z_mar)

!-------- Adjusting the ice thickness, if needed --------

call thk_adjust(time, dtime)

!-------- Calving due to prescribed retreat mask --------

#if (RETREAT_MASK==1 || ICE_SHELF_COLLAPSE_MASK==1)
call calving_retreat_mask(time, dtime)
#endif

end subroutine calc_thk_expl

!-------------------------------------------------------------------------------
!> Ice thickness evolution due to the source term (surface mass balance,
!! basal mass balance and calving).
!-------------------------------------------------------------------------------
subroutine apply_mb_source(dtime, z_mar)

  implicit none

  real(dp), intent(in) :: dtime
  real(dp), intent(in) :: z_mar

  integer(i4b) :: i, j

!-------- Compute new ice thickness H_new_flow due to glacial flow only
!         (no source term considered) --------

  H_new_flow = 0.0_dp

  do i=0, IMAX
  do j=0, JMAX
     if (flag_inner_point(j,i)) then
        H_new_flow(j,i) = H_new(j,i) - dtime*mb_source(j,i)
                          ! new ice thickness due to glacial flow only
     end if
  end do
  end do

!-------- Apply source term --------

  H_new = 0.0_dp

#if (MB_ACCOUNT==0)

  do i=0, IMAX
  do j=0, JMAX
     if (flag_inner_point(j,i)) then
        H_new(j,i) = max(H_new_flow(j,i) + dtime*mb_source(j,i), 0.0_dp)
     end if
  end do
  end do

#elif (MB_ACCOUNT==1)

  do i=2, IMAX-2   ! outermost two grid points excepted,
  do j=2, JMAX-2   ! required for accurate mass balance accounting
     H_new(j,i) = max(H_new_flow(j,i) + dtime*mb_source(j,i), 0.0_dp)
  end do
  end do

#else

  errormsg = ' >>> apply_mb_source: MB_ACCOUNT must be either 0 or 1!'
  call error(errormsg)

#endif

!-------- Reset new ice thickness if needed --------

#if (MARGIN==1)

  do i=0, IMAX
  do j=0, JMAX
     if (mask(j,i) >= 2 .and. H_new(j,i) > eps_H) then
        H_new(j,i) = 0.0_dp
     end if
  end do
  end do

#else   /* MARGIN==2 or 3  */

  do i=0, IMAX
  do j=0, JMAX
     if (zl(j,i) <= Z_ABYSS) then
        H_new(j,i) = 0.0_dp
     end if
  end do
  end do

#endif

end subroutine apply_mb_source

!-------------------------------------------------------------------------------
!> Adjustment of the newly computed ice thickness distribution.
!-------------------------------------------------------------------------------
  subroutine thk_adjust(time, dtime)

#if defined(ALLOW_TAPENADE) /* Tapenade */
#if (THK_EVOL==2)
  use ctrl_m, only: myfloor, myceiling
#endif
#endif /* Tapenade */

  implicit none

  real(dp), intent(in) :: time, dtime

  integer(i4b)                       :: i, j
  integer(i4b)                       :: n1, n2
  real(dp), dimension(0:JMAX,0:IMAX) :: H_new_tmp
  real(dp)                           :: time_in_years
  real(dp)                           :: time1, time2
  real(dp)                           :: dtime_inv

  time_in_years = time*sec2year

!-------- Saving computed H_new before any adjustments --------

  H_new_tmp = H_new

!-------- Correct negative thickness values --------

  do i=0, IMAX  
  do j=0, JMAX 
     if (H_new(j,i) < 0.0_dp) then
        H_new(j,i) = 0.0_dp
     end if
  end do
  end do

!-------- Further adjustments --------

#if (THK_EVOL==0)

!  ------ No evolution of the ice thickness

  H_new = H   ! newly computed ice thickness is discarded

#elif (THK_EVOL==1)

!  ------ Ice thickness evolves freely

  !%% continue

#elif (THK_EVOL==2)

!  ------ Nudging towards prescribed target topography
!                                    with varying relaxation time

  if (time_in_years < real(target_topo_tau0_time_min,dp)) then

     target_topo_tau = target_topo_tau0(0)

  else if (time_in_years < real(target_topo_tau0_time_max,dp)) then

#if !defined(ALLOW_TAPENADE) /* Normal */
     n1 = floor((time_in_years &
             -real(target_topo_tau0_time_min,dp)) &
                /real(target_topo_tau0_time_stp,dp))
#else /* Tapenade */
     call myfloor((time_in_years &
             -real(target_topo_tau0_time_min,dp)) &
                /real(target_topo_tau0_time_stp,dp),n1)
#endif /* Normal vs. Tapenade */

     n1 = max(n1, 0)

#if !defined(ALLOW_TAPENADE) /* Normal */
     n2 = ceiling((time_in_years &
             -real(target_topo_tau0_time_min,dp)) &
                /real(target_topo_tau0_time_stp,dp))
#else /* Tapenade */
     call myceiling(((time_in_years &
              -real(target_topo_tau0_time_min,dp)) &
                 /real(target_topo_tau0_time_stp,dp)),n2)
#endif /* Normal vs. Tapenade */

     n2 = min(n2, ndata_target_topo_tau0)

     if (n1 == n2) then

        target_topo_tau = target_topo_tau0(n1)

     else

        time1 = (target_topo_tau0_time_min + n1*target_topo_tau0_time_stp) &
                *year2sec
        time2 = (target_topo_tau0_time_min + n2*target_topo_tau0_time_stp) &
                *year2sec

        target_topo_tau = target_topo_tau0(n1) &
                  +(target_topo_tau0(n2)-target_topo_tau0(n1)) &
                  *(time-time1)/(time2-time1)
                   ! linear interpolation of the relaxation-time data

     end if

  else

     target_topo_tau = target_topo_tau0(ndata_target_topo_tau0)

  end if

  if (target_topo_tau*sec2year > no_value_pos_1) then
             ! relaxation time target_topo_tau interpreted as infinity

#if (!defined(ALLOW_GRDCHK) && !defined(ALLOW_TAPENADE)) /* Normal */
     target_topo_tau = huge(1.0_dp)
#endif

     !%% H_new = H_new

  else if (target_topo_tau*sec2year < epsi) then
             ! relaxation time target_topo_tau interpreted as zero

     target_topo_tau = 0.0_dp
     H_new           = H_target

  else

     H_new =   ( target_topo_tau*H_new + dtime*H_target ) &
             / ( target_topo_tau       + dtime )

  end if

#elif (THK_EVOL==3)

!  ------ Nudging towards prescribed target topography
!                                    with constant relaxation time

  target_topo_tau = target_topo_tau_0

  if (target_topo_tau*sec2year > no_value_pos_1) then
             ! relaxation time target_topo_tau interpreted as infinity

#if (!defined(ALLOW_GRDCHK) && !defined(ALLOW_TAPENADE)) /* Normal */
     target_topo_tau = huge(1.0_dp)
#endif

     !%% H_new = H_new

  else if (target_topo_tau*sec2year < epsi) then
             ! relaxation time target_topo_tau interpreted as zero

     target_topo_tau = 0.0_dp
     H_new           = H_target

  else

     H_new =   ( target_topo_tau*H_new + dtime*H_target ) &
             / ( target_topo_tau       + dtime )

  end if

#else

  errormsg = ' >>> thk_adjust: THK_EVOL must be between 0 and 3!'
  call error(errormsg)

#endif

!  ------ Maximum ice extent constrained by prescribed mask

  if (flag_mask_maxextent) then

     do i=0, IMAX
     do j=0, JMAX
        if (mask_maxextent(j,i) == 0) then
           H_new(j,i) = 0.0_dp   ! not allowed to glaciate
        end if
     end do
     end do

  end if

!-------- Computation of the mass balance adjustment --------

  dtime_inv = 1.0_dp/dtime

  smb_corr = (H_new-H_new_tmp)*dtime_inv

  as_perp = as_perp + smb_corr
  accum   = accum   + max(smb_corr, 0.0_dp)
  runoff  = runoff  - min(smb_corr, 0.0_dp)
                    ! runoff is counted as positive for mass loss

  end subroutine thk_adjust

!-------------------------------------------------------------------------------
!> Update of the ice-land-ocean mask etc.
!-------------------------------------------------------------------------------
subroutine calc_thk_mask_update(time, dtime, dxi, deta, z_mar, &
                                n_calc_thk_mask_update_aux)

use topograd_m

implicit none

real(dp),     intent(in) :: time
real(dp),     intent(in) :: dtime, dxi, deta
real(dp),     intent(in) :: z_mar
integer(i4b), intent(in) :: n_calc_thk_mask_update_aux

integer(i4b)                       :: i, j
real(dp), dimension(0:JMAX,0:IMAX) :: H_new_tmp
real(dp)                           :: dtime_inv

!-------- Saving computed H_new before any modifications --------

H_new_tmp = H_new

!-------- Update of the mask --------

dtime_inv = 1.0_dp/dtime

if (n_calc_thk_mask_update_aux == 1) then
   call calc_thk_mask_update_aux1(time, dtime, dtime_inv, z_mar)
else if (n_calc_thk_mask_update_aux == 2) then
   call calc_thk_mask_update_aux2(time, dtime, dtime_inv, z_mar)
else if (n_calc_thk_mask_update_aux == 3) then
   call calc_thk_mask_update_aux3(time, dtime, dtime_inv, z_mar)
else
   errormsg = ' >>> calc_thk_mask_update:' &
            //         end_of_line &
            //'        n_calc_thk_mask_update_aux has no valid value!'
   call error(errormsg)
end if

!  ------ Enforce connectivity of the ocean

#if (!defined(OCEAN_CONNECTIVITY) || OCEAN_CONNECTIVITY==1)
   call ocean_connect()
#elif (OCEAN_CONNECTIVITY==0)
   !%% continue
#else
   errormsg = ' >>> calc_thk_mask_update:' &
            //         end_of_line &
            //'        OCEAN_CONNECTIVITY must be either 0 or 1!'
   call error(errormsg)
#endif

!-------- Time derivatives --------

dzs_dtau  = (zs_new-zs)*dtime_inv
dzb_dtau  = (zb_new-zb)*dtime_inv
dzm_dtau  = dH_t_dtau+dzb_dtau
dH_dtau   = (H_new-H)*dtime_inv
dH_c_dtau = dzs_dtau-dzm_dtau

!-------- New gradients --------

#if (TOPOGRAD==0)
call topograd_1(dxi, deta, 2)
#elif (TOPOGRAD==1)
call topograd_2(dxi, deta, 2)
#endif

!-------- Compute the volume flux across the CTS, am_perp --------

#if (CALCMOD==1)

do i=1, IMAX-1
do j=1, JMAX-1

   if ( ((mask(j,i) == 0).or.(mask(j,i) == 3)) &
        .and.(n_cts(j,i) == 1)) then

      am_perp_st(j,i) = &
              (   0.5_dp*(vx_c(0,j,i)+vx_c(0,j,i-1))*dzm_dxi_g(j,i) &
                + 0.5_dp*(vy_c(0,j,i)+vy_c(0,j-1,i))*dzm_deta_g(j,i) ) &
                - 0.5_dp*(vz_c(0,j,i)+vz_t(KTMAX-1,j,i))
      am_perp(j,i) = am_perp_st(j,i) + dzm_dtau(j,i)

   else
      am_perp_st(j,i) = 0.0_dp
      am_perp(j,i)    = 0.0_dp
   end if

end do
end do

#endif

end subroutine calc_thk_mask_update

!-------------------------------------------------------------------------------
!> Update of the ice-land-ocean mask for SIA-only dynamics of ice sheets
!! without ice shelves.
!-------------------------------------------------------------------------------
subroutine calc_thk_mask_update_aux1(time, dtime, dtime_inv, z_mar)

implicit none

real(dp), intent(in) :: time, dtime, dtime_inv, z_mar

integer(i4b) :: i, j, kc, kt
real(dp)     :: rhosw_rho_ratio, rho_rhosw_ratio
real(dp)     :: H_inv

real(dp), dimension(0:JMAX,0:IMAX) :: H_sea_new, H_balance

!-------- Term abbreviations --------

rhosw_rho_ratio = RHO_SW/RHO
rho_rhosw_ratio = RHO/RHO_SW

!-------- Update of the mask --------

zs_new    = zb_new + H_new   ! ice surface topography
H_sea_new = z_sl - zl_new    ! sea depth
H_balance = H_sea_new*rhosw_rho_ratio

#if (THK_EVOL>=1)

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i) <= 1) then   ! grounded ice or ice-free land

      if (H_new(j,i) >= eps_H) then
         mask(j,i) = 0   ! grounded ice
      else
         mask(j,i) = 1   ! ice-free land
      end if

#if (MARGIN==2)

   else   ! (mask(j,i) == 2); sea

      if (H_new(j,i) >= eps_H) then

         if ( H_new(j,i) < (rhosw_rho_ratio*H_sea_new(j,i)) ) then

#if (MARINE_ICE_FORMATION==1)
            mask(j,i) = 2   ! floating ice cut off -> sea
#elif (MARINE_ICE_FORMATION==2)
            mask(j,i) = 0   ! "underwater ice"
#endif
         else
            mask(j,i) = 0   ! grounded ice
         end if

#if (MARINE_ICE_CALVING==2 \
      || MARINE_ICE_CALVING==4 \
      || MARINE_ICE_CALVING==6)
         if (zl0(j,i) < z_mar) mask(j,i) = 2   ! sea
#elif (MARINE_ICE_CALVING==3 \
        || MARINE_ICE_CALVING==5 \
        || MARINE_ICE_CALVING==7)
         if (zl_new(j,i) < z_mar) mask(j,i) = 2   ! sea
#endif

      else
         mask(j,i) = 2   ! sea
      end if

#endif

   end if

end do
end do

#endif

!  ------ Adjustment due to prescribed target topography

#if (THK_EVOL==2 || THK_EVOL==3)
if (target_topo_tau*sec2year < epsi) mask = mask_target
#endif

!-------- Correction of zs_new and zb_new for ice-free land and sea --------

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i) == 1) then   ! ice-free land

      zs_new(j,i) = zb_new(j,i)   ! this prevents zs_new(j,i)
      H_new(j,i)  = 0.0_dp        ! from being below zb_new(j,i)

   else if (mask(j,i) == 2) then   ! sea

      zs_new(j,i) = zb_new(j,i)   ! this prevents zs_new(j,i)
      H_new(j,i)  = 0.0_dp        ! from being below zb_new(j,i)

   else if (mask(j,i) == 3) then   ! floating ice

      errormsg = ' >>> calc_thk_mask_update_aux1:' &
               //           end_of_line &
               //'          mask(j,i)==3 not allowed for' &
               //           end_of_line &
               //'          SIA-only dynamics!'
      call error(errormsg)

   end if

end do
end do

!-------- Limit thickness of isolated ice points --------

call limit_thickness_isolated_ice()

!-------- Computation of further quantities --------

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i) == 0) then   ! grounded ice

      if (n_cts(j,i) == 1) then
         if (H(j,i) > 0.0_dp) then
            H_inv        = 1.0_dp/H(j,i)
            H_c_new(j,i) = H_c(j,i) * H_new(j,i)*H_inv
            H_t_new(j,i) = H_t(j,i) * H_new(j,i)*H_inv
         else
            H_c_new(j,i) = 0.99_dp * H_new(j,i)   ! this case should not occur,
            H_t_new(j,i) = 0.01_dp * H_new(j,i)   ! just for safety
         end if
         zm_new(j,i) = zb_new(j,i)+H_t_new(j,i)
      else
         H_c_new(j,i) = H_new(j,i)
         H_t_new(j,i) = 0.0_dp
         zm_new(j,i)  = zb_new(j,i)
      end if

   else   ! mask(j,i) == 1 or 2, ice-free land or sea

      H_c_new(j,i) = 0.0_dp
      H_t_new(j,i) = 0.0_dp
      H_new(j,i)   = 0.0_dp
      zs_new(j,i)  = zb_new(j,i)
      zm_new(j,i)  = zb_new(j,i)

   end if

end do
end do

end subroutine calc_thk_mask_update_aux1

!-------------------------------------------------------------------------------
!> Update of the ice-land-ocean mask for hybrid SIA/SStA dynamics of ice sheets
!! without ice shelves.
!-------------------------------------------------------------------------------
subroutine calc_thk_mask_update_aux2(time, dtime, dtime_inv, z_mar)

implicit none

real(dp), intent(in) :: time, dtime, dtime_inv, z_mar

integer(i4b) :: i, j, kc, kt
real(dp)     :: rhosw_rho_ratio, rho_rhosw_ratio
real(dp)     :: H_inv

real(dp), dimension(0:JMAX,0:IMAX) :: H_sea_new, H_balance

!-------- Term abbreviations --------

rhosw_rho_ratio = RHO_SW/RHO
rho_rhosw_ratio = RHO/RHO_SW

!-------- Update of the mask --------

zs_new    = zb_new + H_new   ! ice surface topography
H_sea_new = z_sl - zl_new    ! sea depth
H_balance = H_sea_new*rhosw_rho_ratio

#if (THK_EVOL>=1)

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i) <= 1) then   ! grounded ice or ice-free land

      if (H_new(j,i) >= eps_H) then
         mask(j,i) = 0   ! grounded ice
      else
         mask(j,i) = 1   ! ice-free land
      end if

#if (MARGIN==2)

   else   ! (mask(j,i) == 2); sea

      if (H_new(j,i) >= eps_H) then

         if ( H_new(j,i) < (rhosw_rho_ratio*H_sea_new(j,i)) ) then

#if (MARINE_ICE_FORMATION==1)
            mask(j,i) = 2   ! floating ice cut off -> sea
#elif (MARINE_ICE_FORMATION==2)
            mask(j,i) = 0   ! "underwater ice"
#endif
         else
            mask(j,i) = 0   ! grounded ice
         end if

#if (MARINE_ICE_CALVING==2 \
      || MARINE_ICE_CALVING==4 \
      || MARINE_ICE_CALVING==6)
         if (zl0(j,i) < z_mar) mask(j,i) = 2   ! sea
#elif (MARINE_ICE_CALVING==3 \
        || MARINE_ICE_CALVING==5 \
        || MARINE_ICE_CALVING==7)
         if (zl_new(j,i) < z_mar) mask(j,i) = 2   ! sea
#endif

      else
         mask(j,i) = 2   ! sea
      end if

#endif

   end if

end do
end do

#endif

!  ------ Adjustment due to prescribed target topography

#if (THK_EVOL==2 || THK_EVOL==3)
if (target_topo_tau*sec2year < epsi) mask = mask_target
#endif

!-------- Correction of zs_new and zb_new for ice-free land and sea --------

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i) == 1) then   ! ice-free land

      zs_new(j,i) = zb_new(j,i)   ! this prevents zs_new(j,i)
      H_new(j,i)  = 0.0_dp        ! from being below zb_new(j,i)

   else if (mask(j,i) == 2) then   ! sea

      zs_new(j,i) = zb_new(j,i)   ! this prevents zs_new(j,i)
      H_new(j,i)  = 0.0_dp        ! from being below zb_new(j,i)

   else if (mask(j,i) == 3) then   ! floating ice

      errormsg = ' >>> calc_thk_mask_update_aux2:' &
               //           end_of_line &
               //'          mask(j,i)==3 not allowed for' &
               //           end_of_line &
               //'          hybrid SIA/SStA dynamics of ice sheets' &
               //           end_of_line &
               //'          without ice shelves!'
      call error(errormsg)

   end if

end do
end do

!-------- Limit thickness of isolated ice points --------

call limit_thickness_isolated_ice()

!-------- Computation of further quantities --------

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i) == 0) then   ! grounded ice

      if (n_cts(j,i) == 1) then
         if (H(j,i) > 0.0_dp) then
            H_inv        = 1.0_dp/H(j,i)
            H_c_new(j,i) = H_c(j,i) * H_new(j,i)*H_inv
            H_t_new(j,i) = H_t(j,i) * H_new(j,i)*H_inv
         else
            H_c_new(j,i) = 0.99_dp * H_new(j,i)   ! this case should not occur,
            H_t_new(j,i) = 0.01_dp * H_new(j,i)   ! just for safety
         end if
         zm_new(j,i) = zb_new(j,i)+H_t_new(j,i)
      else
         H_c_new(j,i) = H_new(j,i)
         H_t_new(j,i) = 0.0_dp
         zm_new(j,i)  = zb_new(j,i)
      end if

   else   ! mask(j,i) == 1 or 2, ice-free land or sea

      H_c_new(j,i) = 0.0_dp
      H_t_new(j,i) = 0.0_dp
      H_new(j,i)   = 0.0_dp
      zs_new(j,i)  = zb_new(j,i)
      zm_new(j,i)  = zb_new(j,i)

   end if

end do
end do

end subroutine calc_thk_mask_update_aux2

!-------------------------------------------------------------------------------
!> Update of the ice-land-ocean mask for coupled SIA/SSA or
!! SIA/SStA/SSA dynamics of ice sheets with ice shelves.
!-------------------------------------------------------------------------------
subroutine calc_thk_mask_update_aux3(time, dtime, dtime_inv, z_mar)

implicit none

real(dp), intent(in) :: time, dtime, dtime_inv, z_mar

integer(i4b) :: i, j, kc, kt
real(dp)     :: rhosw_rho_ratio, rho_rhosw_ratio
real(dp)     :: H_inv
logical      :: flag_calving_event

real(dp), dimension(0:JMAX,0:IMAX) :: H_sea_new, H_balance

!-------- Term abbreviations --------

rhosw_rho_ratio = RHO_SW/RHO
rho_rhosw_ratio = RHO/RHO_SW

!-------- Update of the mask --------

zs_new    = zb_new + H_new   ! ice surface topography
H_sea_new = z_sl - zl_new    ! sea depth
H_balance = H_sea_new*rhosw_rho_ratio

#if (THK_EVOL>=1)

do i=1, IMAX-1
do j=1, JMAX-1

   ! grounding_line migration check

   if ( ( mask(j,i) <= 1 &
         .and.    (mask(j,i+1)>1.or.mask(j,i-1)>1 &
               .or.mask(j+1,i)>1.or.mask(j-1,i)>1) ) &
        .or. &
        ( mask(j,i)>=2 &
         .and.    (mask(j,i+1)==0.or.mask(j,i-1)==0  &
               .or.mask(j+1,i)==0.or.mask(j-1,i)==0) ) ) then

      if (H_new(j,i) >= eps_H) then

         if (H_new(j,i)<H_balance(j,i).and.zl_new(j,i)<z_sl(j,i)) then
            mask(j,i)     = 3
            zb_new(j,i)   = z_sl(j,i)-rho_rhosw_ratio*H_new(j,i)
            dzb_dtau(j,i) = dtime_inv*(zb_new(j,i)-zb(j,i))
            zs_new(j,i)   = zb_new(j,i)+H_new(j,i)
         else
            mask(j,i)     = 0
            zb_new(j,i)   = zl_new(j,i)
            dzb_dtau(j,i) = dzl_dtau(j,i)
            zs_new(j,i)   = zb_new(j,i)+H_new(j,i)
         end if

      else   ! if (H_new(j,i) <= eps_H) then

         if (zl_new(j,i)>z_sl(j,i)) then
            mask(j,i)     = 1
            zb_new(j,i)   = zl_new(j,i)
            dzb_dtau(j,i) = dtime_inv*(zb_new(j,i)-zb(j,i))
            zs_new(j,i)   = zb_new(j,i)
         else
            mask(j,i)     = 2
            zb_new(j,i)   = z_sl(j,i)
            dzb_dtau(j,i) = 0.0_dp
            zs_new(j,i)   = z_sl(j,i)
         end if

      end if

   else if (mask(j,i) <= 1) then   ! grounded ice or ice-free land

      if (H_new(j,i) >= eps_H) then   ! can change mask 0 or 1
         mask(j,i)     = 0
         zb_new(j,i)   = zl_new(j,i)
         dzb_dtau(j,i) = dzl_dtau(j,i)
         zs_new(j,i)   = zb_new(j,i)+H_new(j,i)
      else
         mask(j,i)     = 1
         zb_new(j,i)   = zl_new(j,i)
         dzb_dtau(j,i) = dzl_dtau(j,i)
         zs_new(j,i)   = zl_new(j,i)
      end if

   else   ! if (mask(j,i)==2.or.mask(j,i)==3) then

      if (H_new(j,i) > eps_H) then

         if (H_new(j,i)<H_balance(j,i).and.zl_new(j,i)<z_sl(j,i)) then
            mask(j,i)     = 3
            zb_new(j,i)   = z_sl(j,i)-rho_rhosw_ratio*H_new(j,i)
            dzb_dtau(j,i) = dtime_inv*(zb_new(j,i)-zb(j,i))
            zs_new(j,i)   = zb_new(j,i)+H_new(j,i)
         else
            mask(j,i)     = 0
            zb_new(j,i)   = zl_new(j,i)
            dzb_dtau(j,i) = dzl_dtau(j,i)
            zs_new(j,i)   = zb_new(j,i)+H_new(j,i)
         end if

      else   ! if (H_new(j,i) <= eps_H) then

         if (zl_new(j,i)>z_sl(j,i)) then
            mask(j,i)     = 1
            zb_new(j,i)   = zl_new(j,i)
            dzb_dtau(j,i) = dtime_inv*(zb_new(j,i)-zb(j,i))
            zs_new(j,i)   = zb_new(j,i)
         else
            mask(j,i)     = 2
            zb_new(j,i)   = z_sl(j,i)
            dzb_dtau(j,i) = 0.0_dp
            zs_new(j,i)   = z_sl(j,i)
         end if

      end if

   end if

end do
end do

#if (ICE_SHELF_CALVING==2)

#if !defined(ALLOW_TAPENADE) /* Normal */

do

   flag_calving_front_1 = .false.
   flag_calving_event   = .false.

   do i=1, IMAX-1
   do j=1, JMAX-1

      if ( (mask(j,i)==3) &   ! floating ice
           .and. &
             (    (mask(j,i+1)==2)   &   ! with
              .or.(mask(j,i-1)==2)   &   ! one
              .or.(mask(j+1,i)==2)   &   ! neighbouring
              .or.(mask(j-1,i)==2) ) &   ! sea point
         ) &
         flag_calving_front_1(j,i) = .true.   ! preliminary detection
                                              ! of the calving front

      if ( flag_calving_front_1(j,i).and.(H_new(j,i) < H_CALV) ) then
         flag_calving_event = .true.  ! calving event,
         mask(j,i)          = 2   ! floating ice point changes to sea point
      end if

   end do
   end do

   if (.not.flag_calving_event) exit

end do

#else /* Tapenade */

   flag_calving_event   = .true.

do while (flag_calving_event)

   flag_calving_front_1 = .false.
   flag_calving_event   = .false.

   do i=1, IMAX-1
   do j=1, JMAX-1

      if ( (mask(j,i)==3) &   ! floating ice
           .and. &
             (    (mask(j,i+1)==2)   &   ! with
              .or.(mask(j,i-1)==2)   &   ! one
              .or.(mask(j+1,i)==2)   &   ! neighbouring
              .or.(mask(j-1,i)==2) ) &   ! sea point
         ) &
         flag_calving_front_1(j,i) = .true.       ! preliminary detection
                                                  ! of the calving front

      if ( (flag_calving_front_1(j,i)).and.(H_new(j,i) < H_CALV) ) then
         flag_calving_event = .true.  ! calving event,
         mask(j,i)          = 2   ! floating ice point changes to sea point
      end if

   end do
   end do

end do

#endif /* Normal vs. Tapenade */

#elif (ICE_SHELF_CALVING==3)

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i)==3) mask(j,i) = 2
      ! float-kill: all floating ice points changed to sea points

end do
end do

#endif

#endif

!  ------ Adjustment due to prescribed target topography

#if (THK_EVOL==2 || THK_EVOL==3)
if (target_topo_tau*sec2year < epsi) mask = mask_target
#endif

!-------- Correction of zs_new and zb_new
!         for ice-free land, sea and floating ice --------

do i=0, IMAX
do j=0, JMAX

   if (mask(j,i) == 1) then   ! ice-free land

      zs_new(j,i) = zb_new(j,i)   ! this prevents zs_new(j,i)
      H_new(j,i)  = 0.0_dp        ! from being below zb_new(j,i)

   else if (mask(j,i) == 2) then   ! sea

      zs_new(j,i) = z_sl(j,i)   ! both zs_new(j,i) and zb_new(j,i)
      zb_new(j,i) = z_sl(j,i)   ! set to the current sea level stand
      H_new(j,i)  = 0.0_dp

   else if (mask(j,i) == 3) then   ! floating ice

      H_new(j,i) = zs_new(j,i)-zb_new(j,i)   ! ice thickness

      zb_new(j,i) = z_sl(j,i) - rho_rhosw_ratio*H_new(j,i)  ! floating condition
      zs_new(j,i) = zb_new(j,i) + H_new(j,i)

   end if

end do
end do

!-------- Limit thickness of isolated ice points --------

call limit_thickness_isolated_ice()

!-------- Computation of further quantities --------

do i=0, IMAX
do j=0, JMAX

   if ( (mask(j,i) == 0).or.(mask(j,i) == 3) ) then
                                        ! grounded or floating ice

      if (n_cts(j,i) == 1) then
         if (H(j,i) > 0.0_dp) then
            H_inv        = 1.0_dp/H(j,i)
            H_c_new(j,i) = H_c(j,i) * H_new(j,i)*H_inv
            H_t_new(j,i) = H_t(j,i) * H_new(j,i)*H_inv
         else
            H_c_new(j,i) = 0.99_dp * H_new(j,i)   ! this case should not occur,
            H_t_new(j,i) = 0.01_dp * H_new(j,i)   ! just for safety
         end if
         zm_new(j,i) = zb_new(j,i)+H_t_new(j,i)
      else
         H_c_new(j,i) = H_new(j,i)
         H_t_new(j,i) = 0.0_dp
         zm_new(j,i)  = zb_new(j,i)
      end if

   else   ! mask(j,i) == 1 or 2, ice-free land or sea

      H_c_new(j,i) = 0.0_dp
      H_t_new(j,i) = 0.0_dp
      H_new(j,i)   = 0.0_dp
      zs_new(j,i)  = zb_new(j,i)
      zm_new(j,i)  = zb_new(j,i)

   end if

end do
end do

end subroutine calc_thk_mask_update_aux3

!-------------------------------------------------------------------------------
!> Limit thickness of isolated ice points.
!-------------------------------------------------------------------------------
subroutine limit_thickness_isolated_ice()

implicit none

integer(i4b) :: i, j
real(dp)     :: rho_rhosw_ratio
real(dp)     :: H_isol_max
logical      :: flag_kill_isolated_ice

rho_rhosw_ratio = RHO/RHO_SW

#if (defined(H_ISOL_MAX))
  H_isol_max = H_ISOL_MAX
#else
  H_isol_max = 1.0e+03_dp   ! default value of 1000 m
#endif

if (H_isol_max < eps_H) then
   H_isol_max = 0.0_dp
   flag_kill_isolated_ice = .true.
else
   flag_kill_isolated_ice = .false.
end if

do i=1, IMAX-1
do j=1, JMAX-1

   if (mask(j,i) == 0) then   ! grounded ice

      if ( ((mask(j,i+1) == 1).or.(mask(j,i+1) == 2)) &
           .and. ((mask(j,i-1) == 1).or.(mask(j,i-1) == 2)) &
           .and. ((mask(j+1,i) == 1).or.(mask(j+1,i) == 2)) &
           .and. ((mask(j-1,i) == 1).or.(mask(j-1,i) == 2)) &
         ) then   ! all nearest neighbours ice-free

         H_new(j,i)  = min(H_new(j,i), H_isol_max)
         zs_new(j,i) = zb_new(j,i)+H_new(j,i)

         if (flag_kill_isolated_ice) then
            if (zb_new(j,i) >= z_sl(j,i)) then
               mask(j,i) = 1
            else
               mask(j,i) = 2
            end if
         end if

      end if

   else if (mask(j,i) == 3) then   ! floating ice

      if ( ((mask(j,i+1) == 1).or.(mask(j,i+1) == 2)) &
           .and. ((mask(j,i-1) == 1).or.(mask(j,i-1) == 2)) &
           .and. ((mask(j+1,i) == 1).or.(mask(j+1,i) == 2)) &
           .and. ((mask(j-1,i) == 1).or.(mask(j-1,i) == 2)) &
         ) then   ! all nearest neighbours ice-free

         H_new(j,i)  = min(H_new(j,i), H_isol_max)
         zb_new(j,i) = z_sl(j,i)-rho_rhosw_ratio*H_new(j,i)
         zs_new(j,i) = zb_new(j,i)+H_new(j,i)

         if (flag_kill_isolated_ice) then
            mask(j,i) = 2
         end if

      end if

   end if

end do
end do

end subroutine limit_thickness_isolated_ice

!-------------------------------------------------------------------------------
!> Enforce connectivity of the ocean.
!-------------------------------------------------------------------------------
subroutine ocean_connect()

implicit none

integer(i4b)                           :: i, j
integer(i4b), dimension(0:JMAX,0:IMAX) :: mask_connect
integer(i4b), dimension(0:JMAX,0:IMAX) :: mask_connect_save, mask_connect_diff
logical                                :: flag_change

!-------- Determine connected area allowed to be ocean --------

mask_connect = 0

mask_connect(0:1         , :          ) = 1
mask_connect(JMAX-1:JMAX , :          ) = 1
mask_connect(:           , 0:1        ) = 1
mask_connect(:           , IMAX-1:IMAX) = 1

flag_change = .true.

do while (flag_change)

   mask_connect_save = mask_connect

   do i=1, IMAX-1
   do j=1, JMAX-1

      if (mask_connect_save(j,i) == 1) then
         if (mask(j  ,i+1) >= 2) mask_connect(j  ,i+1) = 1
         if (mask(j  ,i-1) >= 2) mask_connect(j  ,i-1) = 1
         if (mask(j+1,i  ) >= 2) mask_connect(j+1,i  ) = 1
         if (mask(j-1,i  ) >= 2) mask_connect(j-1,i  ) = 1
         if (mask(j+1,i+1) >= 2) mask_connect(j+1,i+1) = 1
         if (mask(j+1,i-1) >= 2) mask_connect(j+1,i-1) = 1
         if (mask(j-1,i+1) >= 2) mask_connect(j-1,i+1) = 1
         if (mask(j-1,i-1) >= 2) mask_connect(j-1,i-1) = 1
      end if

   end do
   end do

#if !defined(ALLOW_TAPENADE) /* Normal */

   mask_connect_diff = abs(mask_connect-mask_connect_save)

   if (maxval(mask_connect_diff) > 0) then
      flag_change = .true.
   else
      flag_change = .false.
   end if

#else /* Tapenade */

   flag_change = .false.

   do i=0, IMAX
   do j=0, JMAX
      mask_connect_diff(j,i) = mask_connect(j,i)-mask_connect_save(j,i)
      if (mask_connect_diff(j,i) /= 0) flag_change = .true.
   end do
   end do

#endif /* Normal vs. Tapenade */

end do

!-------- Reset disconnected "ocean islands" to ice-free land --------

#if !defined(ALLOW_TAPENADE) /* Normal */

where ((mask == 2).and.(mask_connect == 0)) mask = 1

#else /* Tapenade */

do i=0, IMAX
do j=0, JMAX
   if ((mask(j,i) == 2).and.(mask_connect(j,i) == 0)) &
      mask(j,i) = 1
end do
end do

#endif /* Normal vs. Tapenade */

end subroutine ocean_connect

!-------------------------------------------------------------------------------
!> Determination of the several components of the mass balance:
!! Accumulation (precipitation), runoff, net SMB, calving, basal melt. 
!-------------------------------------------------------------------------------
subroutine account_mb_source(dtime, z_mar)

  implicit none

  real(dp), intent(in) :: dtime
  real(dp), intent(in) :: z_mar

  integer(i4b) :: i, j
  real(dp)     :: dtime_inv

  dtime_inv = 1.0_dp/dtime

!-------- Compute mass balance components --------

  mb_source_apl      = 0.0_dp
  as_perp_apl        = 0.0_dp
  accum_apl          = 0.0_dp
  runoff_apl         = 0.0_dp
  Q_b_apl            = 0.0_dp
  calving_apl        = 0.0_dp
  mask_ablation_type = 0

  do i=0, IMAX
  do j=0, JMAX

     if (flag_inner_point(j,i)) then   ! inner point

        if (H_new(j,i) >= eps_H) then   ! glaciated point

           mb_source_apl(j,i) = (H_new(j,i)-H_new_flow(j,i))*dtime_inv
                                ! applied MB, based on ice thickness change

           ! Store melt quantities here for later accounting
           accum_apl(j,i)   = accum(j,i)
           runoff_apl(j,i)  = runoff(j,i)
           Q_b_apl(j,i)     = Q_b_tot(j,i)
           calving_apl(j,i) = calving(j,i)
           as_perp_apl(j,i) = accum_apl(j,i) - runoff_apl(j,i)

           ! Store melting mask value. This became an ice point,
           ! should be either grounded or floating ice.
           if (mask(j,i)==0) then
              mask_ablation_type(j,i) = 1 ! grounded ice
           else if (mask(j,i)==3) then
              mask_ablation_type(j,i) = 3 ! floating ice
           else
              mask_ablation_type(j,i) = 9 ! misaccounted (if appears)
           end if

        else if (H_new(j,i) < eps_H .and. H_new_flow(j,i) >= eps_H) then
                                           ! ice disappeared

           mb_source_apl(j,i) = (H_new(j,i)-H_new_flow(j,i))*dtime_inv
                                ! applied MB, based on ice thickness change

           accum_apl(j,i) = accum(j,i)

           if (mask(j,i)==2) then
              ! all mass that flows into the ocean counted as
              ! (large-scale) calving
              calving_apl(j,i) = -(mb_source_apl(j,i)-accum_apl(j,i))

           else if (runoff(j,i) > 0.0_dp .or. calving(j,i) > 0.0_dp &
                                         .or. Q_b_tot(j,i) > 0.0_dp) then
              ! land or shelf where ice disappered,
              ! three melting contingents shared proportionally
              runoff_apl(j,i)  = -(mb_source_apl(j,i)-accum_apl(j,i)) &
                                       *runoff(j,i) &
                                       /(runoff(j,i)+calving(j,i) &
                                                    +Q_b_tot(j,i))
              calving_apl(j,i) = -(mb_source_apl(j,i)-accum_apl(j,i)) &
                                       *calving(j,i) &
                                       /(runoff(j,i)+calving(j,i) &
                                                    +Q_b_tot(j,i))
              Q_b_apl(j,i)     = -(mb_source_apl(j,i)-accum_apl(j,i)) &
                                       *Q_b_tot(j,i) &
                                       /(runoff(j,i)+calving(j,i) &
                                                    +Q_b_tot(j,i))

           !%% else runoff(j,i)=0.0_dp .and. calving(j,i) = 0.0_dp &
           !%%                         .and. Q_b_tot(j,i) = 0.0_dp

           end if

           as_perp_apl(j,i) = accum_apl(j,i) - runoff_apl(j,i)

           ! Store melting mask value.
           ! This grid point is ice free now, therefore it
           ! can be either ocean or land (using actual mask value).

           if (mask(j,i) == 2) then   ! ocean point (hidden melt)
              mask_ablation_type(j,i) = -2
           else                            ! land (hidden melt)
              mask_ablation_type(j,i) = -1
           end if

        end if

     end if

  end do
  end do

end subroutine account_mb_source

!-------------------------------------------------------------------------------

end module calc_thk_m
!
