!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  s i c o _ m a t h s _ m _ s t u b
!
!! Stub file for the mathematical tools in the sico_maths_m module.
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
!> Stub file for the mathematical tools in the sico_maths_m module.
!-------------------------------------------------------------------------------
module sico_maths_m

  use sico_types_m

  public :: sor_sprs, tri_sle, my_erfc
#if (MARGIN==3 || DYNAMICS==2)
  public :: sico_lis_solver
#endif

  interface sor_sprs
     module procedure sor_sprs_stub
  end interface

  interface tri_sle
     module procedure tri_sle_stub
  end interface

  interface bilinint 
     module procedure bilinint_stub
  end interface

  interface my_erfc
     module procedure my_erfc_stub
  end interface

#if (MARGIN==3 || DYNAMICS==2)
  interface sico_lis_solver
     module procedure sico_lis_solver_stub
  end interface
#endif

contains

!-------------------------------------------------------------------------------
!> SOR solver for a system of linear equations lgs_a*lgs_x=lgs_b
!! [matrix storage: compressed sparse row CSR,
!! represented by arrays lgs_a_value(values), lgs_a_index (indices)
!! and lgs_a_ptr (pointers)].
!-------------------------------------------------------------------------------
  subroutine sor_sprs_stub(lgs_a_value, lgs_a_index, lgs_a_diag_index, &
                           lgs_a_ptr, &
                           lgs_b_value, &
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
  real(dp), dimension(nmax) :: lgs_x_value_prev
  real(dp) :: temp1, temp2
  logical  :: isnanflag1, isnanflag2, isnanflag3
  real(dp) :: b_nr
  logical  :: flag_convergence

#if (ITER_MAX_SOR > 0)
  iter_max = ITER_MAX_SOR
#else
  iter_max = 1000   ! default value
#endif

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
        return
     end if

  end do iter_loop

  write(6,'(10x,a,i0)') 'sor_sprs: iter = ', iter
  ierr = -1   ! convergence criterion not fulfilled

  end subroutine sor_sprs_stub

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
  subroutine bilinint_stub(x1, x2, y1, y2, z11, z12, z21, z22, x, y, retval)

  implicit none

  real(dp), intent(in) :: x1, x2, y1, y2, z11, z12, z21, z22, x, y

  real(dp) :: t, u
  real(dp), intent(out) :: retval

  real(dp), parameter :: I = 1.0_dp

  t = (x-x1)/(x2-x1)
  u = (y-y1)/(y2-y1)

  retval = (I-t)*(I-u)*z11 + (I-t)*u*z12 + t*(I-u)*z21 + t*u*z22

  end subroutine bilinint_stub

!-------------------------------------------------------------------------------
!> Computation of the complementary error function erfc(x) = 1-erf(x)
!! with a fractional error everywhere less than 1.2 x 10^(-7)
!! (formula by Press et al., 'Numerical Recipes in Fortran 77').
!-------------------------------------------------------------------------------
  subroutine my_erfc_stub(x, retval)

  implicit none

  real(dp), intent(in)  :: x
  real(dp), intent(out) :: retval

  real(dp) :: t, z

  z = abs(x)
  t = 1.0_dp/(1.0_dp+0.5_dp*z)

  retval = t * exp( -z*z     -1.26551223_dp &
                     + t  * (  1.00002368_dp &
                     + t  * (  0.37409196_dp &
                     + t  * (  0.09678418_dp &
                     + t  * ( -0.18628806_dp &
                     + t  * (  0.27886807_dp &
                     + t  * ( -1.13520398_dp &
                     + t  * (  1.48851587_dp &
                     + t  * ( -0.82215223_dp &
                     + t  *    0.17087277_dp ) ) ) ) ) ) ) ) )

  if (x < 0.0_dp) retval = 2.0_dp-retval

  end subroutine my_erfc_stub

#if (MARGIN==3 || DYNAMICS==2)
!-------------------------------------------------------------------------------
!> A stub or dummy subroutine for sico_lis_solver.
!! Note that there are no actual calls to the LIS library since it is treated
!! as a black box, and we only need this stub to be a placeholder.
!-------------------------------------------------------------------------------
  subroutine sico_lis_solver_stub(nmax, nnz, &
                                  lgs_a_ptr, lgs_a_index, &
                                  lgs_a_value, lgs_b_value, lgs_x_value)

  implicit none

  integer(i4b) :: ierr
  integer(i4b) :: iter
  integer(i4b) :: nc, nr

  integer(i4b), intent(in) :: nmax
  integer(i4b), intent(in) :: nnz

  integer(i4b), dimension(nmax+1), intent(in) :: lgs_a_ptr
  integer(i4b), dimension(nnz),    intent(in) :: lgs_a_index
  real(dp),     dimension(nnz),    intent(in) :: lgs_a_value
  real(dp),     dimension(nmax),   intent(in) :: lgs_b_value

#if 0
  real(dp), dimension(nmax), intent(in)    :: lgs_x_value
#else
  real(dp), dimension(nmax), intent(inout) :: lgs_x_value
#endif

  character(len=256) :: ch_solver_set_option

  intrinsic sum

  integer(i4b) :: k
  real(dp)     :: b_nr

  do nr=1, nmax

     b_nr = 0.0_dp

     do k=lgs_a_ptr(nr), lgs_a_ptr(nr+1)-1
        b_nr = b_nr + lgs_a_value(k)*lgs_b_value(lgs_a_index(k))
     end do
     b_nr = b_nr + sum(lgs_x_value)

     lgs_x_value(nr) = lgs_x_value(nr)&
                       -(b_nr-lgs_b_value(nr))
     lgs_x_value(nr) = lgs_x_value(nr) &
                       + sum(lgs_a_value) + sum(lgs_b_value) + sum(lgs_a_ptr)

  end do

  end subroutine sico_lis_solver_stub

#endif

!-------------------------------------------------------------------------------

end module sico_maths_m
!
