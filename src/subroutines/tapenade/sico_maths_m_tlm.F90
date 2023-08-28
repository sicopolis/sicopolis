!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ m a t h s _ m _ t l m
!
!> @file
!!
!! Stub file for solvers in tangent linear (forward) mode of SICOPOLIS-AD v2.
!!
!! @section Copyright
!!
!! Copyright 2009-2023 Shreyas Sunil Gaikwad, Liz Curry-Logan, Ralf Greve
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
!> Stub file for solvers in tangent linear (forward) mode of SICOPOLIS-AD v2.
!<------------------------------------------------------------------------------
module sico_maths_m_diff

  use sico_types_m

  implicit none
  public

  interface tri_sle
      module procedure tri_sle_stub
  end interface

  interface tri_sle_d
      module procedure tri_sle_stub_d
  end interface

  interface my_erfc
      module procedure my_erfc_stub
  end interface

  interface my_erfc_d
      module procedure my_erfc_stub_d
  end interface

  interface sor_sprs
      module procedure sor_sprs_stub
  end interface

  interface sor_sprs_d
      module procedure sor_sprs_stub_d
  end interface

#if defined(BUILD_LIS) && (MARGIN==3 || DYNAMICS==2)
  interface sico_lis_solver
      module procedure sico_lis_solver_stub
  end interface

  interface sico_lis_solver_d
      module procedure sico_lis_solver_stub_d
  end interface
#endif

contains
!-------------------------------------------------------------------------------
!> Differentiation of sor_sprs_stub in forward (tangent) mode:
!! variations of useful results: lgs_x_value,
!! with respect to varying inputs: lgs_b_value lgs_x_value lgs_a_value.
!<------------------------------------------------------------------------------
  subroutine sor_sprs_stub_d(lgs_a_value, lgs_a_valued, lgs_a_index, &
                             lgs_a_diag_index, lgs_a_ptr, &
                             lgs_b_value, lgs_b_valued, &
                             nnz, nmax, &
                             omega, eps_sor, lgs_x_value, lgs_x_valued, ierr)

  implicit none

  integer(i4b), intent(in) :: nnz, nmax
  real(dp), intent(in) :: omega, eps_sor
  integer(i4b), dimension(nmax+1), intent(in) :: lgs_a_ptr
  integer(i4b), dimension(nnz), intent(in) :: lgs_a_index
  integer(i4b), dimension(nmax), intent(in) :: lgs_a_diag_index
  real(dp), dimension(nnz), intent(in) :: lgs_a_value
  real(dp), dimension(nnz), intent(in) :: lgs_a_valued
  real(dp), dimension(nmax), intent(in) :: lgs_b_value
  real(dp), dimension(nmax), intent(in) :: lgs_b_valued
  integer(i4b), intent(out) :: ierr
  real(dp), dimension(nmax), intent(inout) :: lgs_x_value
  real(dp), dimension(nmax), intent(inout) :: lgs_x_valued
  integer(i4b) :: nr, k
  real(dp), dimension(nmax) :: lgs_rhs_value

  call sor_sprs(lgs_a_value, lgs_a_index, lgs_a_diag_index, &
                lgs_a_ptr, lgs_b_value, nnz, nmax, omega, eps_sor, &
                lgs_x_value, ierr)

  do nr=1, nmax

     lgs_rhs_value(nr) = lgs_b_valued(nr)

     do k=lgs_a_ptr(nr), lgs_a_ptr(nr+1)-1
        lgs_rhs_value(nr) = lgs_rhs_value(nr) &
                            - lgs_a_valued(k)*lgs_x_value(lgs_a_index(k))
     end do

  end do

  call sor_sprs(lgs_a_value, lgs_a_index, lgs_a_diag_index, &
                lgs_a_ptr, lgs_rhs_value, nnz, nmax, omega, eps_sor/10000.0, &
                lgs_x_valued, ierr)

  end subroutine sor_sprs_stub_d

!-------------------------------------------------------------------------------
!> SOR solver for a system of linear equations lgs_a*lgs_x=lgs_b
!! [matrix storage: compressed sparse row CSR,
!! represented by arrays lgs_a_value(values), lgs_a_index (indices)
!! and lgs_a_ptr (pointers)].
!<------------------------------------------------------------------------------
  subroutine sor_sprs_stub(lgs_a_value, lgs_a_index, lgs_a_diag_index, &
                           lgs_a_ptr, lgs_b_value, &
                           nnz, nmax, &
                           omega, eps_sor, lgs_x_value, ierr)

  implicit none

  integer(i4b), intent(in) :: nnz, nmax
  real(dp),     intent(in) :: omega, eps_sor

  integer(i4b), dimension(nmax+1), intent(in) :: lgs_a_ptr
  integer(i4b), dimension(nnz),    intent(in) :: lgs_a_index
  integer(i4b), dimension(nmax),   intent(in) :: lgs_a_diag_index
  real(dp),     dimension(nnz),    intent(in) :: lgs_a_value
  real(dp),     dimension(nmax),   intent(in) :: lgs_b_value

  integer(i4b), intent(out) :: ierr

  real(dp), dimension(nmax), intent(inout) :: lgs_x_value

  integer(i4b) :: iter
  integer(i4b) :: iter_max
  integer(i4b) :: nr, k
  real(dp), allocatable, dimension(:) :: lgs_x_value_prev
  real(dp)     :: b_nr
  logical      :: flag_convergence

#if (ITER_MAX_SOR > 0)
  iter_max = ITER_MAX_SOR
#else
  iter_max = 1000   ! default value
#endif

  allocate(lgs_x_value_prev(nmax))

  iter_loop : do iter=1, iter_max

     lgs_x_value_prev = lgs_x_value

     do nr=1, nmax

        b_nr = 0.0_dp 

        do k=lgs_a_ptr(nr), lgs_a_ptr(nr+1)-1
           b_nr = b_nr + lgs_a_value(k)*lgs_x_value(lgs_a_index(k))
        end do

        lgs_x_value(nr) = lgs_x_value(nr) &
                          -omega*(b_nr-lgs_b_value(nr)) &
                                /lgs_a_value(lgs_a_diag_index(nr))

     end do

     flag_convergence = .true.
     do nr=1, nmax
        if (abs(lgs_x_value(nr)-lgs_x_value_prev(nr)) > eps_sor) then
           flag_convergence = .false.
           exit
        end if
     end do

     if (flag_convergence) then
        write(6,'(10x,a,i0)') 'sor_sprs: iter = ', iter
        ierr = 0   ! convergence criterion fulfilled
        deallocate(lgs_x_value_prev)
        return
     end if

  end do iter_loop

  write(6,'(10x,a,i0)') 'sor_sprs: iter = ', iter
  ierr = -1   ! convergence criterion not fulfilled
  deallocate(lgs_x_value_prev)

  end subroutine sor_sprs_stub

!-------------------------------------------------------------------------------
!> Differentiation of tri_sle_stub in forward (tangent) mode:
!! variations of useful results: x b,
!! with respect to varying inputs: x a0 a1 a2 b.
!<------------------------------------------------------------------------------
  subroutine tri_sle_stub_d(a0, a0d, a1, a1d, a2, a2d, x, xd, b, bd, nrows)

  implicit none

  integer(i4b), intent(in) :: nrows
  real(dp), dimension(0:nrows), intent(in) :: a0, a2
  real(dp), dimension(0:nrows), intent(in) :: a0d, a2d
  real(dp), dimension(0:nrows), intent(inout) :: a1, b
  real(dp), dimension(0:nrows), intent(inout) :: a1d, bd
  real(dp), dimension(0:nrows), intent(out) :: x
  real(dp), dimension(0:nrows), intent(out) :: xd
  integer(i4b) :: n
  real(dp), dimension(0:nrows) :: rhs, a1copy
  integer(i4b) :: n1

  do n1 = 0, nrows
     a1copy(n1) = a1(n1)
  end do

  call tri_sle(a0, a1, a2, x, b, nrows)

  rhs(0) = bd(0) - a1d(0)*x(0) - a2d(0)*x(1)
  rhs(nrows) = bd(nrows) - a0d(nrows)*x(nrows-1) - a1d(nrows)*x(nrows)

  do n=1,nrows-1
     rhs(n) = bd(n) - a0d(n)*x(n-1) - a1d(n)*x(n) - a2d(n)*x(n+1)
  end do

  call tri_sle(a0, a1copy, a2, xd, rhs, nrows)

  end subroutine tri_sle_stub_d

!-------------------------------------------------------------------------------
!> Solution of a system of linear equations Ax=b with tridiagonal matrix A.
!<------------------------------------------------------------------------------
  subroutine tri_sle_stub(a0, a1, a2, x, b, nrows)

  implicit none

  integer(i4b),                 intent(in)  :: nrows
  real(dp), dimension(0:nrows), intent(in)  :: a0, a1, a2, b
  real(dp), dimension(0:nrows), intent(out) :: x

     ! a0: a0(j) is element A_(j,j-1) of matrix A
     ! a1: a1(j) is element A_(j,j)   of matrix A
     ! a2: a2(j) is element A_(j,j+1) of matrix A
     ! b: inhomogeneity vector
     ! nrows: size of matrix A (indices run from 0 (!!!) to nrows)
     ! x: solution vector

  integer(i4b) :: n
  real(dp), dimension(0:nrows) :: a0_aux, a1_aux, a2_aux, b_aux, x_aux

!--------  Define local variables --------

  a0_aux = a0
  a1_aux = a1
  a2_aux = a2
  b_aux  = b
  x_aux  = 0.0_dp   ! initialization

!--------  Generate an upper triangular matrix --------

  do n=1, nrows
     a1_aux(n) = a1_aux(n) - a0_aux(n)/a1_aux(n-1)*a2_aux(n-1)
  end do

  do n=1, nrows
     b_aux(n) = b_aux(n) - a0_aux(n)/a1_aux(n-1)*b_aux(n-1)
     ! a0_aux(n) = 0.0_dp , not needed in the following, therefore not set
  end do

!-------- Iterative solution of the new system --------

  x_aux(0) = b_aux(nrows)/a1_aux(nrows)

  do n=1, nrows
     x_aux(n) = b_aux(nrows-n)/a1_aux(nrows-n) &
                -a2_aux(nrows-n)/a1_aux(nrows-n)*x_aux(n-1)
  end do

  do n=0, nrows
     x(n) = x_aux(nrows-n)
  end do

  !  WARNING: Subroutine does not check for elements of the main
  !           diagonal becoming zero. In this case it crashes even
  !           though the system may be solvable. Otherwise ok.

  end subroutine tri_sle_stub

!-------------------------------------------------------------------------------
!> Bilinear interpolation.
!<------------------------------------------------------------------------------
  function bilinint(x1, x2, y1, y2, z11, z12, z21, z22, x, y)
  
  implicit none

  real(dp), intent(in) :: x1, x2, y1, y2, z11, z12, z21, z22, x, y

  real(dp) :: t, u
  real(dp) :: bilinint

  real(dp), parameter :: I = 1.0_dp

  t = (x-x1)/(x2-x1)
  u = (y-y1)/(y2-y1)

  bilinint = (I-t)*(I-u)*z11 + (I-t)*u*z12 + t*(I-u)*z21 + t*u*z22

  end function bilinint

!-------------------------------------------------------------------------------
!> Differentiation of my_erfc_stub in forward (tangent) mode:
!! variations of useful results: retval,
!! with respect to varying inputs: x.
!<------------------------------------------------------------------------------
  subroutine my_erfc_stub_d(x, xd, retval, retvald)

  implicit none

  real(dp), intent(in) :: x
  real(dp), intent(in) :: xd
  real(dp), intent(out) :: retval
  real(dp), intent(out) :: retvald
  real(dp) :: t, z
  real(dp) :: td, zd

  intrinsic abs
  intrinsic exp

  real(dp) :: arg1
  real(dp) :: arg1d
  real(dp) :: temp
  real(dp) :: temp0
  real(dp) :: temp1
  real(dp) :: temp2

  if (x >= 0.) then
     zd = xd
     z = x
  else
     zd = -xd
     z = -x
  end if

  temp = 1.0/(0.5_dp*z+1.0_dp)
  td = -(temp*0.5_dp*zd/(0.5_dp*z+1.0_dp))
  t = temp
  temp = t*(0.17087277_dp*t-0.82215223_dp) + 1.48851587_dp
  temp0 = t*temp - 1.13520398_dp
  temp1 = t*(t*temp0+0.27886807_dp) - 0.18628806_dp
  temp2 = t*(t*temp1+0.09678418_dp) + 0.37409196_dp
  arg1d = (t*temp2+t*(temp2+t*(t*temp1+t*(temp1+t*(t*temp0+t*(temp0+t* &
          (temp+t*(0.17087277_dp*t+t*0.17087277_dp-0.82215223_dp)))+ &
          0.27886807_dp))+0.09678418_dp))+1.00002368_dp)*td - 2*z*zd
  arg1 = t*(t*temp2+1.00002368_dp) - z*z - 1.26551223_dp
  temp2 = exp(arg1)
  retvald = temp2*td + t*exp(arg1)*arg1d
  retval = t*temp2

  if (x < 0.0_dp) then
     retvald = -retvald
     retval = 2.0_dp - retval
  end if

  end subroutine my_erfc_stub_d

!-------------------------------------------------------------------------------
!> Computation of the complementary error function erfc(x) = 1-erf(x)
!! with a fractional error everywhere less than 1.2 x 10^(-7)
!! (formula by Press et al., 'Numerical Recipes in Fortran 77').
!<------------------------------------------------------------------------------
  subroutine my_erfc_stub(x, retval)

  implicit none

  real(dp), intent(in) :: x
  real(dp), intent(out) :: retval
  real(dp) :: t, z

  intrinsic abs
  intrinsic exp

  real(dp) :: arg1

  if (x >= 0.) then
     z = x
  else
     z = -x
  end if

  t = 1.0_dp/(1.0_dp+0.5_dp*z)
  arg1 = -(z*z) - 1.26551223_dp + t*(1.00002368_dp+t*(0.37409196_dp+t* &
         (0.09678418_dp+t*(-0.18628806_dp+t*(0.27886807_dp+t*(- &
         1.13520398_dp+t*(1.48851587_dp+t*(-0.82215223_dp+t*0.17087277_dp)) &
         ))))))
  retval = t*exp(arg1)
  if (x < 0.0_dp) retval = 2.0_dp - retval

  end subroutine my_erfc_stub

#if defined(BUILD_LIS) && (MARGIN==3 || DYNAMICS==2)
!-------------------------------------------------------------------------------
!> Template for Tapenade to help differentiate through the LIS solver.
!<------------------------------------------------------------------------------
#include "lisf.h"
  subroutine sico_lis_solver_stub(nmax, nnz, lgs_a_ptr, lgs_a_index, &
                                  lgs_a_value, lgs_b_value, lgs_x_value)

  implicit none

  integer(i4b) :: ierr
  integer(i4b) :: iter
  integer(i4b) :: nc, nr

  integer(i4b), intent(in) :: nmax
  integer(i4b), intent(in) :: nnz

  integer(i4b), dimension(nmax+1), intent(in) :: lgs_a_ptr
  integer(i4b), dimension(nnz),    intent(in) :: lgs_a_index

  LIS_MATRIX :: lgs_a
  LIS_VECTOR :: lgs_b
  LIS_VECTOR :: lgs_x
  LIS_SOLVER :: solver

  real(dp), dimension(nnz),  intent(in)    :: lgs_a_value
  real(dp), dimension(nmax), intent(in)    :: lgs_b_value
  real(dp), dimension(nmax), intent(inout) :: lgs_x_value

  character(len=256) :: ch_solver_set_option

!  ------ Settings for Lis

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

!  ------ Solution with Lis

  call lis_solver_create(solver, ierr)

#if (defined(LIS_OPTS))
  ch_solver_set_option = trim(LIS_OPTS)
#else
  ch_solver_set_option = '-i bicgsafe -p jacobi '// &
                         '-maxiter 2000 -tol 1.0e-18 -initx_zeros false'
#endif

  call lis_solver_set_option(trim(ch_solver_set_option), solver, ierr)
  call CHKERR(ierr)

  call lis_solve(lgs_a, lgs_b, lgs_x, solver, ierr)
  call CHKERR(ierr)

  call lis_solver_get_iter(solver, iter, ierr)
  write(6,'(10x,a,i0)') 'calc_thk_sia_impl: iter = ', iter

  lgs_x_value = 0.0_dp
  call lis_vector_gather(lgs_x, lgs_x_value, ierr)
  call lis_matrix_destroy(lgs_a, ierr)
  call lis_vector_destroy(lgs_b, ierr)
  call lis_vector_destroy(lgs_x, ierr)
  call lis_solver_destroy(solver, ierr)

  end subroutine sico_lis_solver_stub

!-------------------------------------------------------------------------------
!> Differentiation of sico_lis_solver_stub in forward (tangent) mode.
!<------------------------------------------------------------------------------
  subroutine sico_lis_solver_stub_d(nmax, nnz, &
                                    lgs_a_ptr, lgs_a_index, &
                                    lgs_a_value, lgs_a_valued, lgs_b_value, &
                                    lgs_b_valued, lgs_x_value, lgs_x_valued)

  implicit none

  integer(i4b) :: ierr
  integer(i4b) :: nc, nr

  integer(i4b), intent(in) :: nmax
  integer(i4b), intent(in) :: nnz

  integer(i4b), dimension(nmax+1), intent(in) :: lgs_a_ptr
  integer(i4b), dimension(nnz),    intent(in) :: lgs_a_index

  real(dp), dimension(nnz),  intent(in)    :: lgs_a_value
  real(dp), dimension(nmax), intent(in)    :: lgs_b_value
  real(dp), dimension(nmax), intent(inout) :: lgs_x_value

  real(dp), dimension(nnz),  intent(inout) :: lgs_a_valued
  real(dp), dimension(nmax), intent(inout) :: lgs_b_valued
  real(dp), dimension(nmax), intent(inout) :: lgs_x_valued

  real(dp), dimension(nmax) :: rhs_value

  LIS_SCALAR :: alpha = -1.0
  LIS_MATRIX :: lgs_ad
  LIS_VECTOR :: lgs_bd
  LIS_VECTOR :: lgs_x
  LIS_VECTOR :: rhs
 
  call lis_init_f(ierr)
  call lis_matrix_create(LIS_COMM_WORLD, lgs_ad, ierr)
  call lis_vector_create(LIS_COMM_WORLD, lgs_bd, ierr)
  call lis_vector_create(LIS_COMM_WORLD, lgs_x, ierr)
  call lis_vector_create(LIS_COMM_WORLD, rhs, ierr)

  call lis_matrix_set_size(lgs_ad, 0, nmax, ierr)
  call lis_vector_set_size(lgs_bd, 0, nmax, ierr)
  call lis_vector_set_size(lgs_x, 0, nmax, ierr)
  call lis_vector_set_size(rhs, 0, nmax, ierr)

  do nr=1, nmax

     do nc=lgs_a_ptr(nr), lgs_a_ptr(nr+1)-1
        call lis_matrix_set_value(LIS_INS_VALUE, nr, lgs_a_index(nc), &
                                  lgs_a_valued(nc), lgs_ad, ierr)
     end do

     call lis_vector_set_value(LIS_INS_VALUE, &
                               nr, lgs_b_valued(nr), lgs_bd, ierr)

  end do

  call lis_matrix_set_type(lgs_ad, LIS_MATRIX_CSR, ierr)
  call lis_matrix_assemble(lgs_ad, ierr)

  lgs_x_value = 0.0
  call sico_lis_solver(nmax, nnz, &
                       lgs_a_ptr, lgs_a_index, &
                       lgs_a_value, lgs_b_value, lgs_x_value)

  do nr=1, nmax
     call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_x_value(nr), lgs_x, ierr)
  end do

  rhs_value = 0.0
  call lis_matvec(lgs_ad, lgs_x, rhs, ierr) 
  call lis_vector_xpay(lgs_bd, alpha, rhs, ierr)
  call lis_vector_gather(rhs, rhs_value, ierr)

  lgs_x_valued = 0.0
  call sico_lis_solver(nmax, nnz, &
                       lgs_a_ptr, lgs_a_index, &
                       lgs_a_value, rhs_value, lgs_x_valued)

  call lis_matrix_destroy(lgs_ad, ierr)
  call lis_vector_destroy(lgs_bd, ierr)
  call lis_vector_destroy(lgs_x, ierr)
  call lis_vector_destroy(rhs, ierr)
  call lis_finalize_f(ierr)

  end subroutine sico_lis_solver_stub_d

#endif

!-------------------------------------------------------------------------------

end module sico_maths_m_diff
!
