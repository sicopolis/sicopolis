!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ v z _ m
!
!! Computation of the vertical velocity vz.
!!
!!##### Authors
!!
!! Ralf Greve, Tatsuru Sato
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
!> Computation of the vertical velocity vz.
!-------------------------------------------------------------------------------
module calc_vz_m

  use sico_types_m
  use sico_variables_m
  use sico_vars_m

  implicit none

  private
  public :: calc_vz_grounded, calc_vz_floating, calc_vz_static

contains

!-------------------------------------------------------------------------------
!> Computation of the vertical velocity vz for grounded ice.
!-------------------------------------------------------------------------------
subroutine calc_vz_grounded(dxi, deta, dzeta_c, dzeta_t)

implicit none

real(dp), intent(in) :: dxi, deta, dzeta_c, dzeta_t

integer(i4b) :: i, j, kc, kt
real(dp)     :: avz3(0:KCMAX)
real(dp)     :: cvz00, cvz0(0:KTMAX), cvz1(0:KTMAX), cvz2(0:KTMAX), &
                cvz3(0:KCMAX), cvz4(0:KCMAX), cvz5(0:KCMAX)
real(dp)     :: dxi_inv, deta_inv

!-------- Term abbreviations --------

dxi_inv  = 1.0_dp/dxi
deta_inv = 1.0_dp/deta

do kc=0, KCMAX
   if (flag_aa_nonzero) then
      avz3(kc)  = aa*eaz_c(kc)/(ea-1.0_dp)*dzeta_c
   else
      avz3(kc)  = dzeta_c
   end if
end do

do i=1, IMAX-1
do j=1, JMAX-1

   if (mask(j,i)==0) then   ! grounded ice

!-------- Abbreviations --------

      cvz00 = ( 0.5_dp*(vx_t(0,j,i)+vx_t(0,j,i-1))*dzb_dxi_g(j,i) &
               +0.5_dp*(vy_t(0,j,i)+vy_t(0,j-1,i))*dzb_deta_g(j,i) ) &
              +dzb_dtau(j,i)-Q_b_tot(j,i)

      kt=0
      cvz0(kt) = H_t(j,i)*(insq_g11_g(j,i)*insq_g22_g(j,i)) &
                 *( (vx_t(kt,j,i)  *sq_g22_sgx(j,i) &
                    -vx_t(kt,j,i-1)*sq_g22_sgx(j,i-1))*dxi_inv &
                   +(vy_t(kt,j,i)  *sq_g11_sgy(j,i) &
                    -vy_t(kt,j-1,i)*sq_g11_sgy(j-1,i))*deta_inv ) &
                 *dzeta_t
      cvz1(kt) = (dzb_dxi_g(j,i)+zeta_t(kt)*dH_t_dxi_g(j,i)) &
                 *( (vx_t(kt+1,j,i)+vx_t(kt+1,j,i-1)) &
                   -(vx_t(kt,j,i)  +vx_t(kt,j,i-1)  ) )*0.5_dp
      cvz2(kt) = (dzb_deta_g(j,i)+zeta_t(kt)*dH_t_deta_g(j,i)) &
                 *( (vy_t(kt+1,j,i)+vy_t(kt+1,j-1,i)) &
                   -(vy_t(kt,j,i)  +vy_t(kt,j-1,i)  ) )*0.5_dp

      do kt=1, KTMAX-1
         cvz0(kt) = H_t(j,i)*(insq_g11_g(j,i)*insq_g22_g(j,i)) &
                    *( (vx_t(kt,j,i)  *sq_g22_sgx(j,i) &
                       -vx_t(kt,j,i-1)*sq_g22_sgx(j,i-1))*dxi_inv &
                      +(vy_t(kt,j,i)  *sq_g11_sgy(j,i) &
                       -vy_t(kt,j-1,i)*sq_g11_sgy(j-1,i))*deta_inv ) &
                    *dzeta_t
         cvz1(kt) = (dzb_dxi_g(j,i)+zeta_t(kt)*dH_t_dxi_g(j,i)) &
                    *( (vx_t(kt+1,j,i)+vx_t(kt+1,j,i-1)) &
                      -(vx_t(kt-1,j,i)+vx_t(kt-1,j,i-1)) )*0.25_dp
         cvz2(kt) = (dzb_deta_g(j,i)+zeta_t(kt)*dH_t_deta_g(j,i)) &
                    *( (vy_t(kt+1,j,i)+vy_t(kt+1,j-1,i)) &
                      -(vy_t(kt-1,j,i)+vy_t(kt-1,j-1,i)) )*0.25_dp
      end do

      kt=KTMAX
      cvz0(kt) = H_t(j,i)*(insq_g11_g(j,i)*insq_g22_g(j,i)) &
                 *( (vx_t(kt,j,i)  *sq_g22_sgx(j,i) &
                    -vx_t(kt,j,i-1)*sq_g22_sgx(j,i-1))*dxi_inv &
                   +(vy_t(kt,j,i)  *sq_g11_sgy(j,i) &
                    -vy_t(kt,j-1,i)*sq_g11_sgy(j-1,i))*deta_inv ) &
                 *dzeta_t
      cvz1(kt) = (dzb_dxi_g(j,i)+zeta_t(kt)*dH_t_dxi_g(j,i)) &
                 *( (vx_t(kt,j,i)  +vx_t(kt,j,i-1)  ) &
                   -(vx_t(kt-1,j,i)+vx_t(kt-1,j,i-1)) )*0.5_dp
      cvz2(kt) = (dzb_deta_g(j,i)+zeta_t(kt)*dH_t_deta_g(j,i)) &
                 *( (vy_t(kt,j,i)  +vy_t(kt,j-1,i)  ) &
                   -(vy_t(kt-1,j,i)+vy_t(kt-1,j-1,i)) )*0.5_dp

      kc=0
      cvz3(kc) = avz3(kc)*H_c(j,i) &
                 *(insq_g11_g(j,i)*insq_g22_g(j,i)) &
                 *( (vx_c(kc,j,i)  *sq_g22_sgx(j,i) &
                    -vx_c(kc,j,i-1)*sq_g22_sgx(j,i-1))*dxi_inv &
                   +(vy_c(kc,j,i)  *sq_g11_sgy(j,i) &
                    -vy_c(kc,j-1,i)*sq_g11_sgy(j-1,i))*deta_inv )
      cvz4(kc) = (dzm_dxi_g(j,i) &
                 +eaz_c_quotient(kc)*dH_c_dxi_g(j,i)) &
                 *( (vx_c(kc+1,j,i)+vx_c(kc+1,j,i-1)) &
                   -(vx_c(kc,j,i)  +vx_c(kc,j,i-1)  ) )*0.5_dp
      cvz5(kc) = (dzm_deta_g(j,i) &
                 +eaz_c_quotient(kc)*dH_c_deta_g(j,i)) &
                 *( (vy_c(kc+1,j,i)+vy_c(kc+1,j-1,i)) &
                   -(vy_c(kc,j,i)  +vy_c(kc,j-1,i)  ) )*0.5_dp

      do kc=1, KCMAX-1
         cvz3(kc) = avz3(kc)*H_c(j,i) &
                    *(insq_g11_g(j,i)*insq_g22_g(j,i)) &
                    *( (vx_c(kc,j,i)  *sq_g22_sgx(j,i) &
                       -vx_c(kc,j,i-1)*sq_g22_sgx(j,i-1))*dxi_inv &
                      +(vy_c(kc,j,i)  *sq_g11_sgy(j,i) &
                       -vy_c(kc,j-1,i)*sq_g11_sgy(j-1,i))*deta_inv )
         cvz4(kc) = (dzm_dxi_g(j,i) &
                    +eaz_c_quotient(kc)*dH_c_dxi_g(j,i)) &
                    *( (vx_c(kc+1,j,i)+vx_c(kc+1,j,i-1)) &
                      -(vx_c(kc-1,j,i)+vx_c(kc-1,j,i-1)) )*0.25_dp
         cvz5(kc) = (dzm_deta_g(j,i) &
                    +eaz_c_quotient(kc)*dH_c_deta_g(j,i)) &
                    *( (vy_c(kc+1,j,i)+vy_c(kc+1,j-1,i)) &
                      -(vy_c(kc-1,j,i)+vy_c(kc-1,j-1,i)) )*0.25_dp
      end do

      kc=KCMAX
      cvz3(kc) = avz3(kc)*H_c(j,i) &
                 *(insq_g11_g(j,i)*insq_g22_g(j,i)) &
                 *( (vx_c(kc,j,i)  *sq_g22_sgx(j,i) &
                    -vx_c(kc,j,i-1)*sq_g22_sgx(j,i-1))*dxi_inv &
                   +(vy_c(kc,j,i)  *sq_g11_sgy(j,i) &
                    -vy_c(kc,j-1,i)*sq_g11_sgy(j-1,i))*deta_inv )
      cvz4(kc) = (dzm_dxi_g(j,i) &
                 +eaz_c_quotient(kc)*dH_c_dxi_g(j,i)) &
                 *( (vx_c(kc,j,i)  +vx_c(kc,j,i-1)  ) &
                   -(vx_c(kc-1,j,i)+vx_c(kc-1,j,i-1)) )*0.5_dp
      cvz5(kc) = (dzm_deta_g(j,i) &
                 +eaz_c_quotient(kc)*dH_c_deta_g(j,i)) &
                 *( (vy_c(kc,j,i)  +vy_c(kc,j-1,i)  ) &
                   -(vy_c(kc-1,j,i)+vy_c(kc-1,j-1,i)) )*0.5_dp

!-------- Computation of vz_b --------

      vz_b(j,i) = cvz00

!-------- Computation of vz --------

      if ((n_cts(j,i) == -1).or.(n_cts(j,i) == 0)) then
                        ! cold ice base, temperate ice base

         do kt=0, KTMAX-1
            vz_t(kt,j,i) = vz_b(j,i)
         end do

         vz_m(j,i) = vz_b(j,i)

         vz_c(0,j,i) = vz_m(j,i) &
                       +0.5_dp*(-cvz3(0)+cvz4(0)+cvz5(0))

         do kc=1, KCMAX-1
            vz_c(kc,j,i) = vz_c(kc-1,j,i) &
                           +(-cvz3(kc)+cvz4(kc)+cvz5(kc))
         end do

      else   ! n_cts(j,i) == 1, temperate ice layer

         vz_t(0,j,i) = vz_b(j,i) &
                       +0.5_dp*(-cvz0(0)+cvz1(0)+cvz2(0))

         do kt=1, KTMAX-1
            vz_t(kt,j,i) = vz_t(kt-1,j,i) &
                           +(-cvz0(kt)+cvz1(kt)+cvz2(kt))
         end do

         vz_m(j,i) = vz_t(KTMAX-1,j,i) &
                     +0.5_dp*(-cvz0(KTMAX)+cvz1(KTMAX)+cvz2(KTMAX))

         vz_c(0,j,i) = vz_m(j,i) &
                       +0.5_dp*(-cvz3(0)+cvz4(0)+cvz5(0))

         do kc=1, KCMAX-1
            vz_c(kc,j,i) = vz_c(kc-1,j,i) &
                           +(-cvz3(kc)+cvz4(kc)+cvz5(kc))
         end do

      end if

!-------- Computation of vz_s --------

      kc=KCMAX

      vz_s(j,i) = vz_c(kc-1,j,i) &
                  +0.5_dp*(-cvz3(kc)+cvz4(kc)+cvz5(kc))

   else   ! mask(j,i) /= 0 (not grounded ice)

      vz_b(j,i) = 0.0_dp

      do kt=0, KTMAX-1
         vz_t(kt,j,i) = 0.0_dp
      end do

      vz_m(j,i) = 0.0_dp

      do kc=0, KCMAX-1
         vz_c(kc,j,i) = 0.0_dp
      end do

      vz_s(j,i) = 0.0_dp

   end if

end do
end do

end subroutine calc_vz_grounded

!-------------------------------------------------------------------------------
!> Computation of the vertical velocity vz for floating ice.
!-------------------------------------------------------------------------------
subroutine calc_vz_floating(dxi, deta, dzeta_c)

implicit none

real(dp), intent(in) :: dxi, deta, dzeta_c

integer(i4b) :: i, j, kt, kc
real(dp), dimension(0:JMAX,0:IMAX) :: vz_sl
real(dp), dimension(0:KCMAX) :: zeta_c_sgz, eaz_c_sgz, eaz_c_quotient_sgz
real(dp) :: dxi_inv, deta_inv
real(dp) :: dvx_dxi, dvy_deta

dxi_inv  = 1.0_dp/dxi
deta_inv = 1.0_dp/deta

!-------- Abbreviations --------

do kc=0, KCMAX-1

   zeta_c_sgz(kc) = (kc+0.5_dp)*dzeta_c

   if (flag_aa_nonzero) then
      eaz_c_sgz(kc)          = exp(aa*zeta_c_sgz(kc))
      eaz_c_quotient_sgz(kc) = (eaz_c_sgz(kc)-1.0_dp)/(ea-1.0_dp)
   else
      eaz_c_sgz(kc)          = 1.0_dp
      eaz_c_quotient_sgz(kc) = zeta_c_sgz(kc)
   end if

end do

!-------- Computation of vz --------

do i=1, IMAX-1
do j=1, JMAX-1

   if (mask(j,i)==3) then   ! floating ice

!  ------ Derivatives of the horizontal velocity

      dvx_dxi =  (insq_g11_g(j,i)*insq_g22_g(j,i)) &
                 *(vx_m(j,i)*sq_g22_sgx(j,i)-vx_m(j,i-1)*sq_g22_sgx(j,i-1)) &
                 *dxi_inv

      dvy_deta = (insq_g11_g(j,i)*insq_g22_g(j,i)) &
                 *(vy_m(j,i)*sq_g11_sgy(j,i)-vy_m(j-1,i)*sq_g11_sgy(j-1,i)) &
                 *deta_inv

!  ------ Basal velocity vz_b

      vz_b(j,i) = ( 0.5_dp*(vx_m(j,i)+vx_m(j,i-1))*dzb_dxi_g(j,i) &
                   +0.5_dp*(vy_m(j,i)+vy_m(j-1,i))*dzb_deta_g(j,i) ) &
                  +dzb_dtau(j,i)-Q_b_tot(j,i)
                  ! kinematic boundary condition at the ice base

!  ------ Velocity at sea level vz_sl

      vz_sl(j,i) = vz_b(j,i) - (z_sl(j,i)-zb(j,i))*(dvx_dxi+dvy_deta)

!  ------ Surface velocity vz_s

      vz_s(j,i) = vz_sl(j,i) - (zs(j,i)-z_sl(j,i))*(dvx_dxi+dvy_deta)

!  ------ Velocity vz_m at the interface between
!                              the upper (kc) and the lower (kt) domain

      if ((n_cts(j,i) == -1).or.(n_cts(j,i) == 0)) then
                        ! cold ice base, temperate ice base

         vz_m(j,i) = vz_b(j,i)

      else   ! n_cts(j,i) == 1, temperate ice layer

         vz_m(j,i) = vz_b(j,i) - H_t(j,i)*(dvx_dxi+dvy_deta)

      end if

!  ------ 3-D velocity vz_c and vz_t

      do kc=0, KCMAX-1
         vz_c(kc,j,i) = vz_sl(j,i) &
                        -(zm(j,i)+eaz_c_quotient_sgz(kc)*H_c(j,i)-z_sl(j,i)) &
                         *(dvx_dxi+dvy_deta)
      end do

      if ((n_cts(j,i) == -1).or.(n_cts(j,i) == 0)) then
                        ! cold ice base, temperate ice base

         do kt=0, KTMAX-1
            vz_t(kt,j,i) = vz_b(j,i)
         end do

      else   ! n_cts(j,i) == 1, temperate ice layer

         do kt=0, KTMAX-1
            vz_t(kt,j,i) = vz_sl(j,i) &
                           -(zb(j,i) &
                             +0.5_dp*(zeta_t(kt)+zeta_t(kt+1))*H_t(j,i) &
                             -z_sl(j,i)) &
                            *(dvx_dxi+dvy_deta)
         end do

      end if

   end if

end do
end do

end subroutine calc_vz_floating

!-------------------------------------------------------------------------------
!> Computation of the vertical velocity vz for static ice.
!-------------------------------------------------------------------------------
subroutine calc_vz_static()

implicit none

integer(i4b) :: i, j, kc, kt

do i=0, IMAX
do j=0, JMAX

   if ((mask(j,i)==0).or.(mask(j,i)==3)) then
                              ! grounded or floating ice

      vz_b(j,i) = dzb_dtau(j,i)-Q_b_tot(j,i)
                  ! kinematic boundary condition at the ice base

      do kt=0, KTMAX-1
         vz_t(kt,j,i) = vz_b(j,i)
      end do

      vz_m(j,i) = vz_b(j,i)

      do kc=0, KCMAX-1
         vz_c(kc,j,i) = vz_b(j,i)
      end do

      vz_s(j,i) = vz_b(j,i)

   else   ! mask(j,i) == (1 or 2)

      vz_b(j,i) = 0.0_dp

      do kt=0, KTMAX-1
         vz_t(kt,j,i) = 0.0_dp
      end do

      vz_m(j,i) = 0.0_dp

      do kc=0, KCMAX-1
         vz_c(kc,j,i) = 0.0_dp
      end do

      vz_s(j,i) = 0.0_dp

   end if

end do
end do

end subroutine calc_vz_static

!-------------------------------------------------------------------------------

end module calc_vz_m
!
