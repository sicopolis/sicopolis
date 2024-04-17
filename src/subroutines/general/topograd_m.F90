!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  t o p o g r a d _ m
!
!! Calculation of topography gradients on the staggered grid and on the grid
!! points (including length rescaling with the corresponding components of the
!! metric tensor).
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
!> Calculation of topography gradients on the staggered grid and on the grid
!! points (including length rescaling with the corresponding components of the
!! metric tensor).
!-------------------------------------------------------------------------------
module topograd_m

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
!> Calculation of topography gradients on the staggered grid and on the grid
!! points (the latter by second-order discretization).
!-------------------------------------------------------------------------------
  subroutine topograd_1(dxi, deta, n_switch)

  implicit none

  integer(i4b), intent(in) :: n_switch
  real(dp),     intent(in) :: dxi, deta

  integer(i4b)                       :: i, j
  real(dp)                           :: dxi_inv, deta_inv
  real(dp), dimension(0:JMAX,0:IMAX) :: zs_aux, zm_aux, zb_aux

  dxi_inv  = 1.0_dp/dxi
  deta_inv = 1.0_dp/deta

!-------- Distinguish between old and new topography data --------

  if (n_switch == 1) then
     zs_aux = zs
     zm_aux = zm
     zb_aux = zb
  else if (n_switch == 2) then
     zs_aux = zs_new
     zm_aux = zm_new
     zb_aux = zb_new
  else
     errormsg = ' >>> topograd_1: Wrong value for n_switch!'
     call error(errormsg)
  end if

!-------- Topography gradients on the staggered grid --------

!  ------ x-derivatives

  do i=0, IMAX-1
  do j=0, JMAX
     dzs_dxi(j,i)  = (zs_aux(j,i+1)-zs_aux(j,i))*dxi_inv &
                     *insq_g11_sgx(j,i)
     dzm_dxi(j,i)  = (zm_aux(j,i+1)-zm_aux(j,i))*dxi_inv &
                     *insq_g11_sgx(j,i)
     dzb_dxi(j,i)  = (zb_aux(j,i+1)-zb_aux(j,i))*dxi_inv &
                     *insq_g11_sgx(j,i)
     dH_c_dxi(j,i) = dzs_dxi(j,i)-dzm_dxi(j,i)
     dH_t_dxi(j,i) = dzm_dxi(j,i)-dzb_dxi(j,i)
  end do
  end do

!  ------ y-derivatives

  do i=0, IMAX
  do j=0, JMAX-1
     dzs_deta(j,i)  = (zs_aux(j+1,i)-zs_aux(j,i))*deta_inv &
                      *insq_g22_sgy(j,i)
     dzm_deta(j,i)  = (zm_aux(j+1,i)-zm_aux(j,i))*deta_inv &
                      *insq_g22_sgy(j,i)
     dzb_deta(j,i)  = (zb_aux(j+1,i)-zb_aux(j,i))*deta_inv &
                      *insq_g22_sgy(j,i)
     dH_c_deta(j,i) = dzs_deta(j,i)-dzm_deta(j,i)
     dH_t_deta(j,i) = dzm_deta(j,i)-dzb_deta(j,i)
  end do
  end do

!-------- Topography gradients on the grid points --------

!  ------ x-derivatives

  do i=1, IMAX-1
  do j=0, JMAX
     dzs_dxi_g(j,i)  = (zs_aux(j,i+1)-zs_aux(j,i-1))*0.5_dp*dxi_inv &
                       *insq_g11_g(j,i)
     dzm_dxi_g(j,i)  = (zm_aux(j,i+1)-zm_aux(j,i-1))*0.5_dp*dxi_inv &
                       *insq_g11_g(j,i)
     dzb_dxi_g(j,i)  = (zb_aux(j,i+1)-zb_aux(j,i-1))*0.5_dp*dxi_inv &
                       *insq_g11_g(j,i)
     dH_c_dxi_g(j,i) = dzs_dxi_g(j,i)-dzm_dxi_g(j,i)
     dH_t_dxi_g(j,i) = dzm_dxi_g(j,i)-dzb_dxi_g(j,i)
  end do
  end do

  do j=0, JMAX
     dzs_dxi_g(j,0)     = (zs_aux(j,1)-zs_aux(j,0))*dxi_inv &
                          *insq_g11_g(j,0)
     dzm_dxi_g(j,0)     = (zm_aux(j,1)-zm_aux(j,0))*dxi_inv &
                          *insq_g11_g(j,0)
     dzb_dxi_g(j,0)     = (zb_aux(j,1)-zb_aux(j,0))*dxi_inv &
                          *insq_g11_g(j,0)
     dH_c_dxi_g(j,0)    = dzs_dxi_g(j,0)-dzm_dxi_g(j,0)
     dH_t_dxi_g(j,0)    = dzm_dxi_g(j,0)-dzb_dxi_g(j,0)
     dzs_dxi_g(j,IMAX)  = (zs_aux(j,IMAX)-zs_aux(j,IMAX-1)) &
                          *dxi_inv &
                          *insq_g11_g(j,IMAX)
     dzm_dxi_g(j,IMAX)  = (zm_aux(j,IMAX)-zm_aux(j,IMAX-1)) &
                          *dxi_inv &
                          *insq_g11_g(j,IMAX)
     dzb_dxi_g(j,IMAX)  = (zb_aux(j,IMAX)-zb_aux(j,IMAX-1)) &
                          *dxi_inv &
                          *insq_g11_g(j,IMAX)
     dH_c_dxi_g(j,IMAX) = dzs_dxi_g(j,IMAX)-dzm_dxi_g(j,IMAX)
     dH_t_dxi_g(j,IMAX) = dzm_dxi_g(j,IMAX)-dzb_dxi_g(j,IMAX)
  end do

!  ------ y-derivatives

  do i=0, IMAX
  do j=1, JMAX-1
     dzs_deta_g(j,i)  = (zs_aux(j+1,i)-zs_aux(j-1,i)) &
                        *0.5_dp*deta_inv &
                        *insq_g22_g(j,i)
     dzm_deta_g(j,i)  = (zm_aux(j+1,i)-zm_aux(j-1,i)) &
                        *0.5_dp*deta_inv &
                        *insq_g22_g(j,i)
     dzb_deta_g(j,i)  = (zb_aux(j+1,i)-zb_aux(j-1,i)) &
                        *0.5_dp*deta_inv &
                        *insq_g22_g(j,i)
     dH_c_deta_g(j,i) = dzs_deta_g(j,i)-dzm_deta_g(j,i)
     dH_t_deta_g(j,i) = dzm_deta_g(j,i)-dzb_deta_g(j,i)
  end do
  end do

  do i=0, IMAX
     dzs_deta_g(0,i)     = (zs_aux(1,i)-zs_aux(0,i))*deta_inv &
                           *insq_g22_g(0,i)
     dzm_deta_g(0,i)     = (zm_aux(1,i)-zm_aux(0,i))*deta_inv &
                           *insq_g22_g(0,i)
     dzb_deta_g(0,i)     = (zb_aux(1,i)-zb_aux(0,i))*deta_inv &
                           *insq_g22_g(0,i)
     dH_c_deta_g(0,i)    = dzs_deta_g(0,i)-dzm_deta_g(0,i)
     dH_t_deta_g(0,i)    = dzm_deta_g(0,i)-dzb_deta_g(0,i)
     dzs_deta_g(JMAX,i)  = (zs_aux(JMAX,i)-zs_aux(JMAX-1,i)) &
                           *deta_inv &
                           *insq_g22_g(JMAX,i)
     dzm_deta_g(JMAX,i)  = (zm_aux(JMAX,i)-zm_aux(JMAX-1,i)) &
                           *deta_inv &
                           *insq_g22_g(JMAX,i)
     dzb_deta_g(JMAX,i)  = (zb_aux(JMAX,i)-zb_aux(JMAX-1,i)) &
                           *deta_inv &
                           *insq_g22_g(JMAX,i)
     dH_c_deta_g(JMAX,i) = dzs_deta_g(JMAX,i)-dzm_deta_g(JMAX,i)
     dH_t_deta_g(JMAX,i) = dzm_deta_g(JMAX,i)-dzb_deta_g(JMAX,i)
  end do

  end subroutine topograd_1

!-------------------------------------------------------------------------------
!> Calculation of topography gradients on the staggered grid and on the grid
!! points (the latter by fourth-order discretization).
!-------------------------------------------------------------------------------
  subroutine topograd_2(dxi, deta, n_switch)

  implicit none

  integer(i4b), intent(in) :: n_switch
  real(dp),     intent(in) :: dxi, deta

  integer(i4b)                       :: i, j
  real(dp)                           :: dxi_inv, deta_inv, dxi12_inv, deta12_inv
  real(dp), dimension(0:JMAX,0:IMAX) :: zs_aux, zm_aux, zb_aux

  dxi_inv  = 1.0_dp/dxi
  deta_inv = 1.0_dp/deta

  dxi12_inv  = 1.0_dp/(12.0_dp*dxi)
  deta12_inv = 1.0_dp/(12.0_dp*deta)

!-------- Distinguish between old and new topography data --------

  if (n_switch == 1) then
     zs_aux = zs
     zm_aux = zm
     zb_aux = zb
  else if (n_switch == 2) then
     zs_aux = zs_new
     zm_aux = zm_new
     zb_aux = zb_new
  else
     errormsg = ' >>> topograd_2: Wrong value for n_switch!'
     call error(errormsg)
  end if

!-------- Topography gradients on the staggered grid --------

!  ------ x-derivatives

  do i=0, IMAX-1
  do j=0, JMAX
     dzs_dxi(j,i)  = (zs_aux(j,i+1)-zs_aux(j,i))*dxi_inv &
                     *insq_g11_sgx(j,i)
     dzm_dxi(j,i)  = (zm_aux(j,i+1)-zm_aux(j,i))*dxi_inv &
                     *insq_g11_sgx(j,i)
     dzb_dxi(j,i)  = (zb_aux(j,i+1)-zb_aux(j,i))*dxi_inv &
                     *insq_g11_sgx(j,i)
     dH_c_dxi(j,i) = dzs_dxi(j,i)-dzm_dxi(j,i)
     dH_t_dxi(j,i) = dzm_dxi(j,i)-dzb_dxi(j,i)
  end do
  end do

!  ------ y-derivatives

  do i=0, IMAX
  do j=0, JMAX-1
     dzs_deta(j,i)  = (zs_aux(j+1,i)-zs_aux(j,i))*deta_inv &
                      *insq_g22_sgy(j,i)
     dzm_deta(j,i)  = (zm_aux(j+1,i)-zm_aux(j,i))*deta_inv &
                      *insq_g22_sgy(j,i)
     dzb_deta(j,i)  = (zb_aux(j+1,i)-zb_aux(j,i))*deta_inv &
                      *insq_g22_sgy(j,i)
     dH_c_deta(j,i) = dzs_deta(j,i)-dzm_deta(j,i)
     dH_t_deta(j,i) = dzm_deta(j,i)-dzb_deta(j,i)
  end do
  end do

!-------- Topography gradients on the grid points --------

!  ------ x-derivatives

  do i=2, IMAX-2
  do j=0, JMAX
     dzs_dxi_g(j,i) &
        = ( (8.0_dp*zs_aux(j,i+1)-zs_aux(j,i+2)) &
           -(8.0_dp*zs_aux(j,i-1)-zs_aux(j,i-2)) ) &
          *dxi12_inv &
          *insq_g11_g(j,i)
     dzm_dxi_g(j,i) &
        = ( (8.0_dp*zm_aux(j,i+1)-zm_aux(j,i+2)) &
           -(8.0_dp*zm_aux(j,i-1)-zm_aux(j,i-2)) ) &
          *dxi12_inv &
          *insq_g11_g(j,i)
     dzb_dxi_g(j,i) &
        = ( (8.0_dp*zb_aux(j,i+1)-zb_aux(j,i+2)) &
           -(8.0_dp*zb_aux(j,i-1)-zb_aux(j,i-2)) ) &
          *dxi12_inv &
          *insq_g11_g(j,i)
     dH_c_dxi_g(j,i) = dzs_dxi_g(j,i)-dzm_dxi_g(j,i)
     dH_t_dxi_g(j,i) = dzm_dxi_g(j,i)-dzb_dxi_g(j,i)
  end do
  end do

  do j=0, JMAX

     dzs_dxi_g(j,0)       = (zs_aux(j,1)-zs_aux(j,0))*dxi_inv &
                            *insq_g11_g(j,0)
     dzm_dxi_g(j,0)       = (zm_aux(j,1)-zm_aux(j,0))*dxi_inv &
                            *insq_g11_g(j,0)
     dzb_dxi_g(j,0)       = (zb_aux(j,1)-zb_aux(j,0))*dxi_inv &
                            *insq_g11_g(j,0)
     dH_c_dxi_g(j,0)      = dzs_dxi_g(j,0)-dzm_dxi_g(j,0)
     dH_t_dxi_g(j,0)      = dzm_dxi_g(j,0)-dzb_dxi_g(j,0)

     dzs_dxi_g(j,1)       = (zs_aux(j,2)-zs_aux(j,0)) &
                            *0.5_dp*dxi_inv &
                            *insq_g11_g(j,1)
     dzm_dxi_g(j,1)       = (zm_aux(j,2)-zm_aux(j,0)) &
                            *0.5_dp*dxi_inv &
                            *insq_g11_g(j,1)
     dzb_dxi_g(j,1)       = (zb_aux(j,2)-zb_aux(j,0)) &
                            *0.5_dp*dxi_inv &
                            *insq_g11_g(j,1)
     dH_c_dxi_g(j,1)      = dzs_dxi_g(j,1)-dzm_dxi_g(j,1)
     dH_t_dxi_g(j,1)      = dzm_dxi_g(j,1)-dzb_dxi_g(j,1)

     dzs_dxi_g(j,IMAX-1)  = (zs_aux(j,IMAX)-zs_aux(j,IMAX-2)) &
                            *0.5_dp*dxi_inv &
                            *insq_g11_g(j,IMAX-1)
     dzm_dxi_g(j,IMAX-1)  = (zm_aux(j,IMAX)-zm_aux(j,IMAX-2)) &
                            *0.5_dp*dxi_inv &
                            *insq_g11_g(j,IMAX-1)
     dzb_dxi_g(j,IMAX-1)  = (zb_aux(j,IMAX)-zb_aux(j,IMAX-2)) &
                            *0.5_dp*dxi_inv &
                            *insq_g11_g(j,IMAX-1)
     dH_c_dxi_g(j,IMAX-1) = dzs_dxi_g(j,IMAX-1) &
                            -dzm_dxi_g(j,IMAX-1)
     dH_t_dxi_g(j,IMAX-1) = dzm_dxi_g(j,IMAX-1) &
                            -dzb_dxi_g(j,IMAX-1)

     dzs_dxi_g(j,IMAX)    = (zs_aux(j,IMAX)-zs_aux(j,IMAX-1)) &
                            *dxi_inv &
                            *insq_g11_g(j,IMAX)
     dzm_dxi_g(j,IMAX)    = (zm_aux(j,IMAX)-zm_aux(j,IMAX-1)) &
                            *dxi_inv &
                            *insq_g11_g(j,IMAX)
     dzb_dxi_g(j,IMAX)    = (zb_aux(j,IMAX)-zb_aux(j,IMAX-1)) &
                            *dxi_inv &
                            *insq_g11_g(j,IMAX)
     dH_c_dxi_g(j,IMAX)   = dzs_dxi_g(j,IMAX)-dzm_dxi_g(j,IMAX)
     dH_t_dxi_g(j,IMAX)   = dzm_dxi_g(j,IMAX)-dzb_dxi_g(j,IMAX)

  end do

!  ------ y-derivatives

  do i=0, IMAX
  do j=2, JMAX-2
     dzs_deta_g(j,i) &
        = ( (8.0_dp*zs_aux(j+1,i)-zs_aux(j+2,i)) &
           -(8.0_dp*zs_aux(j-1,i)-zs_aux(j-2,i)) ) &
          *deta12_inv &
          *insq_g22_g(j,i)
     dzm_deta_g(j,i) &
        = ( (8.0_dp*zm_aux(j+1,i)-zm_aux(j+2,i)) &
           -(8.0_dp*zm_aux(j-1,i)-zm_aux(j-2,i)) ) &
          *deta12_inv &
          *insq_g22_g(j,i)
     dzb_deta_g(j,i) &
        = ( (8.0_dp*zb_aux(j+1,i)-zb_aux(j+2,i)) &
           -(8.0_dp*zb_aux(j-1,i)-zb_aux(j-2,i)) ) &
          *deta12_inv &
          *insq_g22_g(j,i)
     dH_c_deta_g(j,i) = dzs_deta_g(j,i)-dzm_deta_g(j,i)
     dH_t_deta_g(j,i) = dzm_deta_g(j,i)-dzb_deta_g(j,i)
  end do
  end do

  do i=0, IMAX

     dzs_deta_g(0,i)       = (zs_aux(1,i)-zs_aux(0,i))*deta_inv &
                             *insq_g22_g(0,i)
     dzm_deta_g(0,i)       = (zm_aux(1,i)-zm_aux(0,i))*deta_inv &
                             *insq_g22_g(0,i)
     dzb_deta_g(0,i)       = (zb_aux(1,i)-zb_aux(0,i))*deta_inv &
                             *insq_g22_g(0,i)
     dH_c_deta_g(0,i)      = dzs_deta_g(0,i)-dzm_deta_g(0,i)
     dH_t_deta_g(0,i)      = dzm_deta_g(0,i)-dzb_deta_g(0,i)

     dzs_deta_g(1,i)       = (zs_aux(2,i)-zs_aux(0,i)) &
                             *0.5_dp*deta_inv &
                             *insq_g22_g(1,i)
     dzm_deta_g(1,i)       = (zm_aux(2,i)-zm_aux(0,i)) &
                             *0.5_dp*deta_inv &
                             *insq_g22_g(1,i)
     dzb_deta_g(1,i)       = (zb_aux(2,i)-zb_aux(0,i)) &
                             *0.5_dp*deta_inv &
                             *insq_g22_g(1,i)
     dH_c_deta_g(1,i)      = dzs_deta_g(1,i)-dzm_deta_g(1,i)
     dH_t_deta_g(1,i)      = dzm_deta_g(1,i)-dzb_deta_g(1,i)

     dzs_deta_g(JMAX-1,i)  = (zs_aux(JMAX,i)-zs_aux(JMAX-2,i)) &
                             *0.5_dp*deta_inv &
                             *insq_g22_g(JMAX-1,i)
     dzm_deta_g(JMAX-1,i)  = (zm_aux(JMAX,i)-zm_aux(JMAX-2,i)) &
                             *0.5_dp*deta_inv &
                             *insq_g22_g(JMAX-1,i)
     dzb_deta_g(JMAX-1,i)  = (zb_aux(JMAX,i)-zb_aux(JMAX-2,i)) &
                             *0.5_dp*deta_inv &
                             *insq_g22_g(JMAX-1,i)
     dH_c_deta_g(JMAX-1,i) = dzs_deta_g(JMAX-1,i) &
                             -dzm_deta_g(JMAX-1,i)
     dH_t_deta_g(JMAX-1,i) = dzm_deta_g(JMAX-1,i) &
                             -dzb_deta_g(JMAX-1,i)

     dzs_deta_g(JMAX,i)    = (zs_aux(JMAX,i)-zs_aux(JMAX-1,i)) &
                             *deta_inv &
                             *insq_g22_g(JMAX,i)
     dzm_deta_g(JMAX,i)    = (zm_aux(JMAX,i)-zm_aux(JMAX-1,i)) &
                             *deta_inv &
                             *insq_g22_g(JMAX,i)
     dzb_deta_g(JMAX,i)    = (zb_aux(JMAX,i)-zb_aux(JMAX-1,i)) &
                             *deta_inv &
                             *insq_g22_g(JMAX,i)
     dH_c_deta_g(JMAX,i)   = dzs_deta_g(JMAX,i)-dzm_deta_g(JMAX,i)
     dH_t_deta_g(JMAX,i)   = dzm_deta_g(JMAX,i)-dzb_deta_g(JMAX,i)

  end do

  end subroutine topograd_2

!-------------------------------------------------------------------------------

end module topograd_m
!
