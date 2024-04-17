!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ d x y z _ m
!
!! Computation of all components of the strain-rate tensor, the full
!! effective strain rate and the shear fraction.
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
!> Computation of all components of the strain-rate tensor, the full
!! effective strain rate and the shear fraction.
!-------------------------------------------------------------------------------
module calc_dxyz_m

  use sico_types_m
  use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

  use error_m

  implicit none

  private
  public :: calc_dxyz

contains

!-------------------------------------------------------------------------------
!> Main subroutine of calc_dxyz_m:
!! Computation of all components of the strain-rate tensor, the full
!! effective strain rate and the shear fraction.
!-------------------------------------------------------------------------------
  subroutine calc_dxyz(dxi, deta, dzeta_c, dzeta_t)

  implicit none

  real(dp), intent(in) :: dxi, deta, dzeta_c, dzeta_t

  integer(i4b)                       :: i, j, kc, kt
  real(dp)                           :: dxi_inv, deta_inv
  real(dp), dimension(0:JMAX,0:IMAX) :: fact_x, fact_y
  real(dp), dimension(0:KCMAX)       :: fact_z_c
  real(dp)                           :: fact_z_t
  real(dp)                           :: H_c_inv, H_t_inv
  real(dp), dimension(0:KCMAX)       :: lxy_c, lyx_c, lxz_c, lzx_c, lyz_c, lzy_c
  real(dp), dimension(0:KTMAX)       :: lxy_t, lyx_t, lxz_t, lzx_t, lyz_t, lzy_t
  real(dp), dimension(0:KCMAX)       :: shear_c_squared
  real(dp), dimension(0:KTMAX)       :: shear_t_squared
  real(dp)                           :: abs_v_ssa_inv, nx, ny
  real(dp)                           :: shear_x_help, shear_y_help
  real(dp)                           :: lambda_shear_help

!-------- Term abbreviations --------

  dxi_inv  = 1.0_dp/dxi
  deta_inv = 1.0_dp/deta

  fact_x   = dxi_inv *insq_g11_g
  fact_y   = deta_inv*insq_g22_g

  do kc=0, KCMAX
     if (flag_aa_nonzero) then
        fact_z_c(kc)  = (ea-1.0_dp)/(aa*eaz_c(kc)*dzeta_c)
     else
        fact_z_c(kc)  = 1.0_dp/dzeta_c
     end if
  end do

  fact_z_t  = 1.0_dp/dzeta_t

!-------- Initialisation --------

  dxx_c          = 0.0_dp
  dyy_c          = 0.0_dp
  dxy_c          = 0.0_dp
  dxz_c          = 0.0_dp
  dyz_c          = 0.0_dp
  de_c           = 0.0_dp
  lambda_shear_c = 0.0_dp

  dxx_t          = 0.0_dp
  dyy_t          = 0.0_dp
  dxy_t          = 0.0_dp
  dxz_t          = 0.0_dp
  dyz_t          = 0.0_dp
  de_t           = 0.0_dp
  lambda_shear_t = 0.0_dp

!-------- Computation --------

  do i=1, IMAX-1
  do j=1, JMAX-1

     if ((mask(j,i) == 0).or.(mask(j,i) == 3)) then
                                                 ! grounded or floating ice

        H_c_inv = 1.0_dp/(abs(H_c(j,i))+eps_dp)

        kc=0

           dxx_c(kc,j,i) = (vx_c(kc,j,i)-vx_c(kc,j,i-1))*fact_x(j,i)
           dyy_c(kc,j,i) = (vy_c(kc,j,i)-vy_c(kc,j-1,i))*fact_y(j,i)

           lxy_c(kc) = (  (vx_c(kc,j+1,i)+vx_c(kc,j+1,i-1)) &
                        - (vx_c(kc,j-1,i)+vx_c(kc,j-1,i-1)) ) &
                       *0.25_dp*fact_y(j,i)

           lyx_c(kc) = (  (vy_c(kc,j,i+1)+vy_c(kc,j-1,i+1)) &
                        - (vy_c(kc,j,i-1)+vy_c(kc,j-1,i-1)) ) &
                       *0.25_dp*fact_x(j,i)

           lzx_c(kc) = (vz_m(j,i+1)-vz_m(j,i-1))*0.5_dp*fact_x(j,i)
           lzy_c(kc) = (vz_m(j+1,i)-vz_m(j-1,i))*0.5_dp*fact_y(j,i)

           lxz_c(kc) = (  (vx_c(kc+1,j,i)+vx_c(kc+1,j,i-1)) &
                        - (vx_c(kc  ,j,i)+vx_c(kc  ,j,i-1)) ) &
                       *0.5_dp*fact_z_c(kc)*H_c_inv

           lyz_c(kc) = (  (vy_c(kc+1,j,i)+vy_c(kc+1,j-1,i)) &
                        - (vy_c(kc  ,j,i)+vy_c(kc  ,j-1,i)) ) &
                       *0.5_dp*fact_z_c(kc)*H_c_inv

        ! end kc=0

        do kc=1, KCMAX-1

           dxx_c(kc,j,i) = (vx_c(kc,j,i)-vx_c(kc,j,i-1))*fact_x(j,i)
           dyy_c(kc,j,i) = (vy_c(kc,j,i)-vy_c(kc,j-1,i))*fact_y(j,i)

           lxy_c(kc) = (  (vx_c(kc,j+1,i)+vx_c(kc,j+1,i-1)) &
                        - (vx_c(kc,j-1,i)+vx_c(kc,j-1,i-1)) ) &
                       *0.25_dp*fact_y(j,i)

           lyx_c(kc) = (  (vy_c(kc,j,i+1)+vy_c(kc,j-1,i+1)) &
                        - (vy_c(kc,j,i-1)+vy_c(kc,j-1,i-1)) ) &
                       *0.25_dp*fact_x(j,i)

           lzx_c(kc) = (  (vz_c(kc,j,i+1)+vz_c(kc-1,j,i+1)) &
                        - (vz_c(kc,j,i-1)+vz_c(kc-1,j,i-1)) ) &
                       *0.25_dp*fact_x(j,i)

           lzy_c(kc) = (  (vz_c(kc,j+1,i)+vz_c(kc-1,j+1,i)) &
                        - (vz_c(kc,j-1,i)+vz_c(kc-1,j-1,i)) ) &
                       *0.25_dp*fact_y(j,i)

           lxz_c(kc) = (  (vx_c(kc+1,j,i)+vx_c(kc+1,j,i-1)) &
                        - (vx_c(kc-1,j,i)+vx_c(kc-1,j,i-1)) ) &
                       *0.25_dp*fact_z_c(kc)*H_c_inv

           lyz_c(kc) = (  (vy_c(kc+1,j,i)+vy_c(kc+1,j-1,i)) &
                        - (vy_c(kc-1,j,i)+vy_c(kc-1,j-1,i)) ) &
                       *0.25_dp*fact_z_c(kc)*H_c_inv

        end do

        kc=KCMAX

           dxx_c(kc,j,i) = (vx_c(kc,j,i)-vx_c(kc,j,i-1))*fact_x(j,i)
           dyy_c(kc,j,i) = (vy_c(kc,j,i)-vy_c(kc,j-1,i))*fact_y(j,i)

           lxy_c(kc) = (  (vx_c(kc,j+1,i)+vx_c(kc,j+1,i-1)) &
                        - (vx_c(kc,j-1,i)+vx_c(kc,j-1,i-1)) ) &
                       *0.25_dp*fact_y(j,i)

           lyx_c(kc) = (  (vy_c(kc,j,i+1)+vy_c(kc,j-1,i+1)) &
                        - (vy_c(kc,j,i-1)+vy_c(kc,j-1,i-1)) ) &
                       *0.25_dp*fact_x(j,i)

           lzx_c(kc) = (vz_s(j,i+1)-vz_s(j,i-1))*0.5_dp*fact_x(j,i)
           lzy_c(kc) = (vz_s(j+1,i)-vz_s(j-1,i))*0.5_dp*fact_y(j,i)

           lxz_c(kc) = (  (vx_c(kc  ,j,i)+vx_c(kc  ,j,i-1)) &
                        - (vx_c(kc-1,j,i)+vx_c(kc-1,j,i-1)) ) &
                       *0.5_dp*fact_z_c(kc)*H_c_inv

           lyz_c(kc) = (  (vy_c(kc  ,j,i)+vy_c(kc  ,j-1,i)) &
                        - (vy_c(kc-1,j,i)+vy_c(kc-1,j-1,i)) ) &
                       *0.5_dp*fact_z_c(kc)*H_c_inv

        ! end kc=KCMAX

        dxy_c(:,j,i) = 0.5_dp*(lxy_c+lyx_c)
        dxz_c(:,j,i) = 0.5_dp*(lxz_c+lzx_c)
        dyz_c(:,j,i) = 0.5_dp*(lyz_c+lzy_c)

        do kc=0, KCMAX

           shear_c_squared(kc) =   dxz_c(kc,j,i)*dxz_c(kc,j,i) &
                                 + dyz_c(kc,j,i)*dyz_c(kc,j,i)

           de_c(kc,j,i) =  sqrt( ( ( dxx_c(kc,j,i)*dxx_c(kc,j,i) &
                                    +dyy_c(kc,j,i)*dyy_c(kc,j,i) ) &
                                 + ( dxx_c(kc,j,i)*dyy_c(kc,j,i) &
                                    +dxy_c(kc,j,i)*dxy_c(kc,j,i) ) &
                                 + shear_c_squared(kc) ) &
                               + eps_dp*eps_dp )
           
           lambda_shear_c(kc,j,i) = sqrt(shear_c_squared(kc)) &
                                                 /(de_c(kc,j,i)+eps_dp)

        end do

#if (CALCMOD==1)

        if ((n_cts(j,i) == -1).or.(n_cts(j,i) == 0)) then
                          ! cold ice base, temperate ice base

           dxx_t(:,j,i)          = dxx_c(0,j,i)
           dyy_t(:,j,i)          = dyy_c(0,j,i)
           dxy_t(:,j,i)          = dxy_c(0,j,i)
           dxz_t(:,j,i)          = dxz_c(0,j,i)
           dyz_t(:,j,i)          = dyz_c(0,j,i)
           de_t(:,j,i)           = de_c(0,j,i)
           lambda_shear_t(:,j,i) = lambda_shear_c(0,j,i)

        else   ! n_cts(j,i) == 1, temperate ice layer

           H_t_inv = 1.0_dp/(abs(H_t(j,i))+eps_dp)

           kt=0

              dxx_t(kt,j,i) = (vx_t(kt,j,i)-vx_t(kt,j,i-1))*fact_x(j,i)
              dyy_t(kt,j,i) = (vy_t(kt,j,i)-vy_t(kt,j-1,i))*fact_y(j,i)

              lxy_t(kt) = (  (vx_t(kt,j+1,i)+vx_t(kt,j+1,i-1)) &
                           - (vx_t(kt,j-1,i)+vx_t(kt,j-1,i-1)) ) &
                          *0.25_dp*fact_y(j,i)

              lyx_t(kt) = (  (vy_t(kt,j,i+1)+vy_t(kt,j-1,i+1)) &
                           - (vy_t(kt,j,i-1)+vy_t(kt,j-1,i-1)) ) &
                          *0.25_dp*fact_x(j,i)

              lzx_t(kt) = (vz_b(j,i+1)-vz_b(j,i-1))*0.5_dp*fact_x(j,i)
              lzy_t(kt) = (vz_b(j+1,i)-vz_b(j-1,i))*0.5_dp*fact_y(j,i)

              lxz_t(kt) = (  (vx_t(kt+1,j,i)+vx_t(kt+1,j,i-1)) &
                           - (vx_t(kt  ,j,i)+vx_t(kt  ,j,i-1)) ) &
                          *0.5_dp*fact_z_t*H_t_inv

              lyz_t(kt) = (  (vy_t(kt+1,j,i)+vy_t(kt+1,j-1,i)) &
                           - (vy_t(kt  ,j,i)+vy_t(kt  ,j-1,i)) ) &
                          *0.5_dp*fact_z_t*H_t_inv

           ! end kt=0

           do kt=1, KTMAX-1

              dxx_t(kt,j,i) = (vx_t(kt,j,i)-vx_t(kt,j,i-1))*fact_x(j,i)
              dyy_t(kt,j,i) = (vy_t(kt,j,i)-vy_t(kt,j-1,i))*fact_y(j,i)

              lxy_t(kt) = (  (vx_t(kt,j+1,i)+vx_t(kt,j+1,i-1)) &
                           - (vx_t(kt,j-1,i)+vx_t(kt,j-1,i-1)) ) &
                          *0.25_dp*fact_y(j,i)

              lyx_t(kt) = (  (vy_t(kt,j,i+1)+vy_t(kt,j-1,i+1)) &
                           - (vy_t(kt,j,i-1)+vy_t(kt,j-1,i-1)) ) &
                          *0.25_dp*fact_x(j,i)

              lzx_t(kt) = (  (vz_t(kt,j,i+1)+vz_t(kt-1,j,i+1)) &
                           - (vz_t(kt,j,i-1)+vz_t(kt-1,j,i-1)) ) &
                          *0.25_dp*fact_x(j,i)

              lzy_t(kt) = (  (vz_t(kt,j+1,i)+vz_t(kt-1,j+1,i)) &
                           - (vz_t(kt,j-1,i)+vz_t(kt-1,j-1,i)) ) &
                          *0.25_dp*fact_y(j,i)

              lxz_t(kt) = (  (vx_t(kt+1,j,i)+vx_t(kt+1,j,i-1)) &
                           - (vx_t(kt-1,j,i)+vx_t(kt-1,j,i-1)) ) &
                          *0.25_dp*fact_z_t*H_t_inv

              lyz_t(kt) = (  (vy_t(kt+1,j,i)+vy_t(kt+1,j-1,i)) &
                           - (vy_t(kt-1,j,i)+vy_t(kt-1,j-1,i)) ) &
                          *0.25_dp*fact_z_t*H_t_inv

           end do

           kt=KTMAX

              dxx_t(kt,j,i) = (vx_t(kt,j,i)-vx_t(kt,j,i-1))*fact_x(j,i)
              dyy_t(kt,j,i) = (vy_t(kt,j,i)-vy_t(kt,j-1,i))*fact_y(j,i)

              lxy_t(kt) = (  (vx_t(kt,j+1,i)+vx_t(kt,j+1,i-1)) &
                           - (vx_t(kt,j-1,i)+vx_t(kt,j-1,i-1)) ) &
                          *0.25_dp*fact_y(j,i)

              lyx_t(kt) = (  (vy_t(kt,j,i+1)+vy_t(kt,j-1,i+1)) &
                           - (vy_t(kt,j,i-1)+vy_t(kt,j-1,i-1)) ) &
                          *0.25_dp*fact_x(j,i)

              lzx_t(kt) = (vz_m(j,i+1)-vz_m(j,i-1))*0.5_dp*fact_x(j,i)
              lzy_t(kt) = (vz_m(j+1,i)-vz_m(j-1,i))*0.5_dp*fact_y(j,i)

              lxz_t(kt) = (  (vx_t(kt  ,j,i)+vx_t(kt  ,j,i-1)) &
                           - (vx_t(kt-1,j,i)+vx_t(kt-1,j,i-1)) ) &
                          *0.5_dp*fact_z_t*H_t_inv

              lyz_t(kt) = (  (vy_t(kt  ,j,i)+vy_t(kt  ,j-1,i)) &
                           - (vy_t(kt-1,j,i)+vy_t(kt-1,j-1,i)) ) &
                          *0.5_dp*fact_z_t*H_t_inv

           ! end kt=KTMAX

           dxy_t(:,j,i) = 0.5_dp*(lxy_t+lyx_t)
           dxz_t(:,j,i) = 0.5_dp*(lxz_t+lzx_t)
           dyz_t(:,j,i) = 0.5_dp*(lyz_t+lzy_t)

           do kt=0, KTMAX

              shear_t_squared(kt) =   dxz_t(kt,j,i)*dxz_t(kt,j,i) &
                                    + dyz_t(kt,j,i)*dyz_t(kt,j,i)

              de_t(kt,j,i)  = sqrt( ( ( dxx_t(kt,j,i)*dxx_t(kt,j,i) &
                                       +dyy_t(kt,j,i)*dyy_t(kt,j,i) ) &
                                    + ( dxx_t(kt,j,i)*dyy_t(kt,j,i) &
                                       +dxy_t(kt,j,i)*dxy_t(kt,j,i) ) &
                                    + shear_t_squared(kt) ) &
                                  + eps_dp*eps_dp )

              lambda_shear_t(kt,j,i) = sqrt(shear_t_squared(kt)) &
                                                 /(de_t(kt,j,i)+eps_dp)

           end do

        end if

#elif (CALCMOD==0 || CALCMOD==2 || CALCMOD==3 || CALCMOD==-1)

        dxx_t(:,j,i)          = dxx_c(0,j,i)
        dyy_t(:,j,i)          = dyy_c(0,j,i)
        dxy_t(:,j,i)          = dxy_c(0,j,i)
        dxz_t(:,j,i)          = dxz_c(0,j,i)
        dyz_t(:,j,i)          = dyz_c(0,j,i)
        de_t(:,j,i)           = de_c(0,j,i)
        lambda_shear_t(:,j,i) = lambda_shear_c(0,j,i)

#else
        errormsg = ' >>> calc_dxyz: Parameter CALCMOD must be -1, 0, 1, 2 or 3!'
        call error(errormsg)
#endif

!  ------ Modification of the shear fraction for floating ice (ice shelves)

        if (mask(j,i) == 3) then   ! floating ice

           abs_v_ssa_inv = 1.0_dp / &
                           sqrt( ( 0.25_dp*(vx_m(j,i)+vx_m(j,i-1))**2 &
                                  +0.25_dp*(vy_m(j,i)+vy_m(j-1,i))**2 ) &
                                 + eps_dp**2 )

           nx = -0.5_dp*(vy_m(j,i)+vy_m(j-1,i)) * abs_v_ssa_inv
           ny =  0.5_dp*(vx_m(j,i)+vx_m(j,i-1)) * abs_v_ssa_inv

           shear_x_help = ( dxx_c(KCMAX,j,i)*nx + dxy_c(KCMAX,j,i)*ny ) &
                          - (        dxx_c(KCMAX,j,i)*nx**3      &
                             +2.0_dp*dxy_c(KCMAX,j,i)*(nx**2*ny) &
                             +       dyy_c(KCMAX,j,i)*(nx*ny**2) )
                       ! strain rate for ice shelves independent of depth,
                       ! thus surface values used here

           shear_y_help = ( dyy_c(KCMAX,j,i)*ny + dxy_c(KCMAX,j,i)*nx ) &
                          - (        dyy_c(KCMAX,j,i)*ny**3      &
                             +2.0_dp*dxy_c(KCMAX,j,i)*(ny**2*nx) &
                             +       dxx_c(KCMAX,j,i)*(ny*nx**2) )
                       ! strain rate for ice shelves independent of depth,
                       ! thus surface values used here

           lambda_shear_help = sqrt( ( shear_x_help**2 + shear_y_help**2 ) &
                                     + eps_dp**2 ) &
                               / (de_ssa(j,i)+eps_dp)


           lambda_shear_c(:,j,i) = lambda_shear_help
           lambda_shear_t(:,j,i) = lambda_shear_help

        end if

!  ------ Constrain the shear fraction to reasonable [0,1] interval

#if !defined(ALLOW_TAPENADE)
        lambda_shear_c(:,j,i) = min(max(lambda_shear_c(:,j,i), 0.0_dp), 1.0_dp)
        lambda_shear_t(:,j,i) = min(max(lambda_shear_t(:,j,i), 0.0_dp), 1.0_dp)
#else
        lambda_shear_c(:,j,i) = max(lambda_shear_c(:,j,i), 0.0_dp)
        lambda_shear_c(:,j,i) = min(lambda_shear_c(:,j,i), 1.0_dp)
        lambda_shear_t(:,j,i) = max(lambda_shear_t(:,j,i), 0.0_dp)
        lambda_shear_t(:,j,i) = min(lambda_shear_t(:,j,i), 1.0_dp)
#endif

     else   ! mask(j,i) == 1 or 2; ice-free land or ocean

        dxx_c(:,j,i)          = 0.0_dp
        dyy_c(:,j,i)          = 0.0_dp
        dxy_c(:,j,i)          = 0.0_dp
        dxz_c(:,j,i)          = 0.0_dp
        dyz_c(:,j,i)          = 0.0_dp
        de_c(:,j,i)           = 0.0_dp
        lambda_shear_c(:,j,i) = 0.0_dp

        dxx_t(:,j,i)          = 0.0_dp
        dyy_t(:,j,i)          = 0.0_dp
        dxy_t(:,j,i)          = 0.0_dp
        dxz_t(:,j,i)          = 0.0_dp
        dyz_t(:,j,i)          = 0.0_dp
        de_t(:,j,i)           = 0.0_dp
        lambda_shear_t(:,j,i) = 0.0_dp

     end if

  end do
  end do

  end subroutine calc_dxyz

!-------------------------------------------------------------------------------

end module calc_dxyz_m
!
