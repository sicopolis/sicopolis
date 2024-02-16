!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ m a t h s _ m _ g r a d
!
!! Stub file for solvers in adjoint (reverse) mode of SICOPOLIS-AD v2.
!!
!!##### Authors
!!
!! Shreyas Sunil Gaikwad, Liz Curry-Logan, Ralf Greve
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
!> Stub file for solvers in adjoint (reverse) mode of SICOPOLIS-AD v2.
!-------------------------------------------------------------------------------
module sico_maths_m_diff

  use sico_types_m

  implicit none
  public

  interface sor_sprs
      module procedure sor_sprs_stub
  end interface

  interface sor_sprs_b
      module procedure sor_sprs_stub_b
  end interface

  interface tri_sle
      module procedure tri_sle_stub
  end interface

  interface tri_sle_b
      module procedure tri_sle_stub_b
  end interface

  interface my_erfc
      module procedure my_erfc_stub
  end interface

  interface my_erfc_b
      module procedure my_erfc_stub_b
  end interface

#if defined(BUILD_LIS) && (MARGIN==3 || DYNAMICS==2)
  interface sico_lis_solver
      module procedure sico_lis_solver_stub
  end interface

  interface sico_lis_solver_b
      module procedure sico_lis_solver_stub_b
  end interface
#endif

contains

!-------------------------------------------------------------------------------
!> Subroutine to transpose matrix lgs_a for the SOR solver.
!! [matrix storage: compressed sparse row CSR,
!! represented by arrays lgs_a_value(values), lgs_a_index (indices)
!! and lgs_a_ptr (pointers)].
!-------------------------------------------------------------------------------
subroutine transpose_csr(a_value, a_index, a_diag_index, a_ptr, &
                         nnz, nmax, &
                         b_value, b_index, b_diag_index, b_ptr)

  implicit none

  integer(i4b),                    intent(in)  :: nnz, nmax
  integer(i4b), dimension(nmax+1), intent(in)  :: a_ptr
  integer(i4b), dimension(nnz),    intent(in)  :: a_index
  integer(i4b), dimension(nmax),   intent(in)  :: a_diag_index
  real(dp),     dimension(nnz),    intent(in)  :: a_value
  integer(i4b), dimension(nmax+1), intent(out) :: b_ptr
  integer(i4b), dimension(nnz),    intent(out) :: b_index
  integer(i4b), dimension(nmax),   intent(out) :: b_diag_index
  real(dp),     dimension(nnz),    intent(out) :: b_value

  integer(i4b) :: nr, k
  integer(i4b), dimension(nmax) :: b_ptr_plus

  b_ptr = 0

  do nr=1, nmax
     do k=a_ptr(nr), a_ptr(nr+1)-1
        b_ptr(a_index(k)+1) = b_ptr(a_index(k)+1) + 1
     end do
  end do

  b_ptr(1) = 1

  do nr=1, nmax
    b_ptr(nr+1) = b_ptr(nr) + b_ptr(nr+1)
  end do

  b_ptr_plus = 0

  do nr=1, nmax
     do k=a_ptr(nr), a_ptr(nr+1)-1
        b_index(b_ptr(a_index(k)) + b_ptr_plus(a_index(k))) = nr
        if (a_index(k) == nr) then
           b_diag_index(nr) = b_ptr(a_index(k)) + b_ptr_plus(a_index(k))
        end if
        b_value(b_ptr(a_index(k)) + b_ptr_plus(a_index(k))) = a_value(k)
        b_ptr_plus(a_index(k)) = b_ptr_plus(a_index(k)) + 1
     end do
  end do

  end subroutine transpose_csr

!-------------------------------------------------------------------------------
!> Differentiation of sor_sprs_stub in reverse (adjoint) mode:
!! gradient of useful results: lgs_x_value,
!! with respect to varying inputs: lgs_b_value lgs_x_value lgs_a_value.
!-------------------------------------------------------------------------------
  subroutine sor_sprs_stub_b(lgs_a_value, lgs_a_valueb, lgs_a_index, &
                             lgs_a_diag_index, lgs_a_ptr, &
                             lgs_b_value, lgs_b_valueb, &
                             nnz, nmax, &
                             omega, eps_sor, lgs_x_value, lgs_x_valueb, ierr)

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
  real(dp), dimension(nnz),  intent(inout) :: lgs_a_valueb
  real(dp), dimension(nmax), intent(inout) :: lgs_b_valueb
  real(dp), dimension(nmax), intent(inout) :: lgs_x_valueb

  integer(i4b), dimension(nmax+1) :: lgs_aT_ptr
  integer(i4b), dimension(nnz)    :: lgs_aT_index
  integer(i4b), dimension(nmax)   :: lgs_aT_diag_index
  real(dp),     dimension(nnz)    :: lgs_aT_value
  real(dp),     dimension(nmax)   :: incrbb

  integer(i4b) :: iter
  integer(i4b) :: nr, k
  real(dp), allocatable, dimension(:) :: lgs_x_value_prev
  real(dp)     :: b_nr
  logical      :: flag_convergence
  logical      :: inif, infor

  incrbb = 0.0
  call transpose_csr(lgs_a_value, lgs_a_index, lgs_a_diag_index, lgs_a_ptr, &
                     nnz, nmax, &
                     lgs_aT_value, lgs_aT_index, lgs_aT_diag_index, lgs_aT_ptr)

  call sor_sprs(lgs_aT_value, lgs_aT_index, lgs_aT_diag_index, lgs_aT_ptr, &
                lgs_x_valueb, &
                nnz, nmax, omega, eps_sor, incrbb, ierr)

  do nr=1,nmax
     lgs_b_valueb(nr) = lgs_b_valueb(nr) + incrbb(nr)
  end do

  call sor_sprs(lgs_a_value, lgs_a_index, lgs_a_diag_index, lgs_a_ptr, &
                lgs_b_value, &
                nnz, nmax, omega, eps_sor, lgs_x_value, ierr)

  do nr=1, nmax
     inif = .false. 
     infor = .false. 
     if (nr == lgs_a_diag_index(nr)) then
        inif = .true.
        lgs_a_valueb(lgs_a_diag_index(nr)) &
           = lgs_a_valueb(lgs_a_diag_index(nr)) &
             - lgs_x_value(lgs_a_index(nr)) * incrbb(nr)
     end if
     do k=lgs_a_ptr(nr), lgs_a_ptr(nr+1)-1
        infor = .true.
        lgs_a_valueb(k) = lgs_a_valueb(k) &
                          - lgs_x_value(lgs_a_index(k)) *incrbb(nr)
     end do
  end do

  lgs_x_valueb = 0.0

  end subroutine sor_sprs_stub_b

!-------------------------------------------------------------------------------
!> SOR solver for a system of linear equations lgs_a*lgs_x=lgs_b
!! [matrix storage: compressed sparse row CSR,
!! represented by arrays lgs_a_value(values), lgs_a_index (indices)
!! and lgs_a_ptr (pointers)].
!-------------------------------------------------------------------------------
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
!> Differentiation of tri_sle_stub in reverse (adjoint) mode:
!! gradient of useful results: x,
!! with respect to varying inputs: x a0 a1 a2 b.
!-------------------------------------------------------------------------------
  subroutine tri_sle_stub_b(a0, a0b, a1, a1b, a2, a2b, x, xb, b, bb, nrows)

  implicit none

  integer(i4b), intent(in) :: nrows
  real(dp), dimension(0:nrows), intent(in) :: a0, a1, a2, b
  real(dp), dimension(0:nrows) :: a0b, a1b, a2b, bb
  real(dp), dimension(0:nrows) :: x, x_copy
  real(dp), dimension(0:nrows) :: xb
  integer(i4b) :: n

  real(dp), dimension(0:nrows) :: a0T, a1T, a2T
  real(dp), dimension(0:nrows) :: incrbb
  integer(i4b) :: i

  a0T(0) = 0.0   
  a0T(1:nrows) = a2(0:nrows-1)
  a1T(0:nrows) = a1(0:nrows)
  a2T(0:nrows-1) = a0(1:nrows)
  a2T(nrows) = 0.0

  call tri_sle(a0T, a1T, a2T, incrbb, xb, nrows)

  do i=0,nrows
     bb(i) = incrbb(i)
  end do

  ! The results go wrong when x is calculated here, so calculate x_copy.

  call tri_sle(a0, a1, a2, x_copy, b, nrows)

  do i=1,nrows
     a0b(i) = - x_copy(i-1)*incrbb(i)
  end do
  do i=0,nrows
     a1b(i) = - x_copy(i)*incrbb(i)
  end do
  do i=0,nrows
     a2b(i) = - x_copy(i+1)*incrbb(i)
  end do
  xb = 0.0

  end subroutine tri_sle_stub_b

!-------------------------------------------------------------------------------
!> Solution of a system of linear equations Ax=b with tridiagonal matrix A.
!-------------------------------------------------------------------------------
  subroutine tri_sle_stub(a0, a1, a2, x, b, nrows)

  implicit none

  integer(i4b),                 intent(in)  :: nrows
  real(dp), dimension(0:nrows), intent(in)  :: a0, a1, a2, b
  real(dp), dimension(0:nrows), intent(out) :: x

     ! a0: a0(j) is element A_(j,j-1) of matrix A
     ! a1: a1(j) is element A_(j,j)   of matrix A
     ! a2: a2(j) is element A_(j,j+1) of matrix A
     ! b: inhomogeneity vector
     ! nrows: size of matrix A (indices run from 0 (!) to nrows)
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
!-------------------------------------------------------------------------------
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
!> Differentiation of my_erfc_stub in reverse (adjoint) mode:
!! gradient of useful results: retval,
!! with respect to varying inputs: x.
!-------------------------------------------------------------------------------
  subroutine my_erfc_stub_b(x, xb, retval, retvalb)

  implicit none

  real(dp), intent(in) :: x

  real(dp) :: xb
  real(dp) :: retval
  real(dp) :: retvalb
  real(dp) :: t, z
  real(dp) :: tb, zb

  intrinsic abs
  intrinsic exp

  real(dp) :: temp
  real(dp) :: temp0
  real(dp) :: temp1
  real(dp) :: temp2
  real(dp) :: temp3
  real(dp) :: temp4
  real(dp) :: tempb
  real(dp) :: tempb0
  real(dp) :: tempb1
  real(dp) :: tempb2
  integer*4 :: branch

  if (x >= 0.) then
     z = x
     call PUSHCONTROL1B(0)
  else
     z = -x
     call PUSHCONTROL1B(1)
  end if

  t = 1.0_dp/(1.0_dp+0.5_dp*z)

  if (x < 0.0_dp) retvalb = -retvalb

  temp = t*(0.17087277_dp*t-0.82215223_dp) + 1.48851587_dp
  temp0 = t*temp - 1.13520398_dp
  temp1 = t*(t*temp0+0.27886807_dp) - 0.18628806_dp
  temp2 = t*(t*temp1+0.09678418_dp) + 0.37409196_dp
  temp3 = t*temp2 + 1.00002368_dp
  temp4 = t*temp3 - z*z - 1.26551223_dp
  tempb = exp(temp4)*t*retvalb
  tempb0 = t**2*tempb
  tempb1 = t**2*tempb0
  tempb2 = t**2*tempb1

  tb = exp(temp4)*retvalb + (temp3+temp2*t)*tempb + (t*temp1+temp1*t+ &
       0.09678418_dp)*tempb0 + (t*temp0+temp0*t+0.27886807_dp)*tempb1 + ( &
       temp+(0.17087277_dp*t-0.82215223_dp)*t+0.17087277_dp*t**2)*tempb2

  zb = -(2*z*tempb) - 0.5_dp*tb/(0.5_dp*z+1.0_dp)**2

  call POPCONTROL1B(branch)

  if (branch == 0) then
     xb = zb
  else
     xb = -zb
  end if

  end subroutine my_erfc_stub_b

!-------------------------------------------------------------------------------
!> Computation of the complementary error function erfc(x) = 1-erf(x)
!! with a fractional error everywhere less than 1.2 x 10^(-7)
!! (formula by Press et al., 'Numerical Recipes in Fortran 77').
!-------------------------------------------------------------------------------
  subroutine my_erfc_stub(x, retval)

  implicit none

  real(dp), intent(in) :: x
  real(dp), intent(out) :: retval
  real(dp) :: t, z

  intrinsic abs
  intrinsic exp

  if (x >= 0.) then
     z = x
  else
     z = -x
  end if

  t = 1.0_dp/(1.0_dp+0.5_dp*z)

  retval = t*exp(-(z*z)-1.26551223_dp+t*(1.00002368_dp+t*( &
           0.37409196_dp+t*(0.09678418_dp+t*(-0.18628806_dp+t*(0.27886807_dp+ &
           t*(-1.13520398_dp+t*(1.48851587_dp+t*(-0.82215223_dp+t* &
           0.17087277_dp)))))))))

  if (x < 0.0_dp) retval = 2.0_dp - retval

  end subroutine my_erfc_stub

#if defined(BUILD_LIS) && (MARGIN==3 || DYNAMICS==2)
!-------------------------------------------------------------------------------
!> Subroutine to transpose matrix lgs_a for the LIS solver.
!! [matrix storage: compressed sparse row CSR,
!! represented by arrays lgs_a_value(values), lgs_a_index (indices)
!! and lgs_a_ptr (pointers)].
!-------------------------------------------------------------------------------
  subroutine transpose_csr_no_diagonal(a_value, a_index, a_ptr, &
                                       nnz, nmax, &
                                       b_value, b_index, b_ptr)

  implicit none

  integer(i4b),                    intent(in)  :: nnz, nmax
  integer(i4b), dimension(nmax+1), intent(in)  :: a_ptr
  integer(i4b), dimension(nnz),    intent(in)  :: a_index
  real(dp),     dimension(nnz),    intent(in)  :: a_value
  integer(i4b), dimension(nmax+1), intent(out) :: b_ptr
  integer(i4b), dimension(nnz),    intent(out) :: b_index
  real(dp),     dimension(nnz),    intent(out) :: b_value

  integer(i4b) :: nr, k
  integer(i4b), dimension(nmax) :: b_ptr_plus

  b_ptr = 0

  do nr=1, nmax
     do k=a_ptr(nr), a_ptr(nr+1)-1
        b_ptr(a_index(k)+1) = b_ptr(a_index(k)+1) + 1
     end do
  end do

  b_ptr(1) = 1

  do nr=1, nmax
     b_ptr(nr+1) = b_ptr(nr) + b_ptr(nr+1)
  end do

  b_ptr_plus = 0

  do nr=1, nmax
     do k=a_ptr(nr), a_ptr(nr+1)-1
        b_index(b_ptr(a_index(k)) + b_ptr_plus(a_index(k))) = nr
        b_value(b_ptr(a_index(k)) + b_ptr_plus(a_index(k))) = a_value(k)
        b_ptr_plus(a_index(k)) = b_ptr_plus(a_index(k)) + 1
     end do
  end do

  end subroutine transpose_csr_no_diagonal

!-------------------------------------------------------------------------------
!> LIS solver for a system of linear equations lgs_a*lgs_x=lgs_b
!! [matrix storage: compressed sparse row CSR,
!! represented by arrays lgs_a_value(values), lgs_a_index (indices)
!! and lgs_a_ptr (pointers)].
!-------------------------------------------------------------------------------
#include "lisf.h"
  subroutine sico_lis_solver_stub(nmax, nnz, lgs_a_ptr, lgs_a_index, &
                                  lgs_a_value, lgs_b_value, lgs_x_value)

  implicit none

  integer(i4b) :: ierr
  integer(i4b) :: iter
  integer(i4b) :: nc, nr

  integer(i4b), intent(in) :: nmax
  integer(i4b), intent(in) :: nnz

  integer(i4b), dimension(nmax+1),  intent(in) :: lgs_a_ptr
  integer(i4b), dimension(nnz),     intent(in) :: lgs_a_index

  LIS_MATRIX :: lgs_a
  LIS_VECTOR :: lgs_b
  LIS_VECTOR :: lgs_x
  LIS_SOLVER :: solver

  real(dp), dimension(nnz),  intent(in)    :: lgs_a_value
  real(dp), dimension(nmax), intent(in)    :: lgs_b_value
  real(dp), dimension(nmax), intent(inout) :: lgs_x_value

  character(len=256) :: ch_solver_set_option

!  ------ Settings for Lis

  call lis_init_f(ierr)
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

  ch_solver_set_option = trim(LIS_OPTS)

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
  call lis_finalize_f(ierr)

  end subroutine sico_lis_solver_stub

!-------------------------------------------------------------------------------
!> Differentiation of sico_lis_solver_stub in reverse (adjoint) mode:
!! gradient of useful results: lgs_x_value,
!! with respect to varying inputs: lgs_b_value lgs_x_value lgs_a_value.
!-------------------------------------------------------------------------------
  subroutine sico_lis_solver_stub_b(nmax, nnz, lgs_a_ptr, lgs_a_index, &
                                    lgs_a_value, lgs_a_valueb, &
                                    lgs_b_value, lgs_b_valueb, lgs_x_value, &
                                    lgs_x_valueb)

  implicit none

  integer(i4b) :: nc, nr

  integer(i4b), intent(in) :: nmax
  integer(i4b), intent(in) :: nnz

  integer(i4b), dimension(nmax+1), intent(in) :: lgs_a_ptr
  integer(i4b), dimension(nnz),    intent(in) :: lgs_a_index

  real(dp), dimension(nnz),  intent(in)    :: lgs_a_value
  real(dp), dimension(nmax), intent(in)    :: lgs_b_value
  real(dp), dimension(nmax), intent(inout) :: lgs_x_value

  real(dp), dimension(nnz),  intent(inout) :: lgs_a_valueb
  real(dp), dimension(nmax), intent(inout) :: lgs_b_valueb
  real(dp), dimension(nmax), intent(inout) :: lgs_x_valueb

  integer(i4b), dimension(nmax+1) :: lgs_aT_ptr
  integer(i4b), dimension(nnz)    :: lgs_aT_index
  real(dp),     dimension(nnz)    :: lgs_aT_value
  real(dp),     dimension(nmax)   :: incrbb

  lgs_aT_value = 0.
  lgs_aT_index = 0
  lgs_aT_value = 0.
  call transpose_csr_no_diagonal(lgs_a_value, lgs_a_index, lgs_a_ptr, &
                                 nnz, nmax, &
                                 lgs_aT_value, lgs_aT_index, lgs_aT_ptr)

  incrbb = 0.
  call sico_lis_solver(nmax, nnz, &
                       lgs_aT_ptr, lgs_aT_index, &
                       lgs_aT_value, lgs_x_valueb, incrbb)

  do nr=1,nmax
     lgs_b_valueb(nr) = lgs_b_valueb(nr) + incrbb(nr)
  end do

  call sico_lis_solver(nmax, nnz, &
                       lgs_a_ptr, lgs_a_index, &
                       lgs_a_value, lgs_b_value, lgs_x_value)

  do nr=1, nmax
     do nc=lgs_a_ptr(nr), lgs_a_ptr(nr+1)-1
        lgs_a_valueb(nc) = lgs_a_valueb(nc) &
                           - lgs_x_value(lgs_a_index(nc))*incrbb(nr)
     end do
  end do
  lgs_x_valueb = 0.0

  end subroutine sico_lis_solver_stub_b

#endif

!-------------------------------------------------------------------------------

end module sico_maths_m_diff
!
