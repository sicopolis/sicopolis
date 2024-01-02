!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ t e m p _ e n t h _ m
!
!> @file
!!
!! Computation of temperature, water content and age with the enthalpy method.
!!
!! @section Copyright
!!
!! Copyright 2013-2024 Ralf Greve, Heinz Blatter
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
!> Computation of temperature, water content and age with the enthalpy method.
!<------------------------------------------------------------------------------
module calc_temp_enth_m

  use sico_types_m
  use sico_variables_m
  use sico_vars_m
  use error_m

  implicit none

  private
  public :: calc_temp_enth

contains

!-------------------------------------------------------------------------------
!> Main subroutine of calc_temp_enth_m:
!! Computation of temperature, water content and age with the enthalpy method.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth(dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
                          dtime_temp)

implicit none

real(dp), intent(in) :: dxi, deta, dzeta_c, dzeta_t, dzeta_r
real(dp), intent(in) :: dtime_temp

integer(i4b) :: i, j, kc, kt, kr, ii, jj
real(dp) :: at1(0:KCMAX), at2_1(0:KCMAX), at2_2(0:KCMAX), &
            at3_1(0:KCMAX), at3_2(0:KCMAX), at4_1(0:KCMAX), &
            at4_2(0:KCMAX), at5(0:KCMAX), at6(0:KCMAX), at7, &
            ai1(0:KCMAX), ai2(0:KCMAX), &
            atr1, acb1, acb2, acb3, acb4, alb1, aqtlde(0:KCMAX), &
            am1, am3(0:KCMAX)
real(dp) :: dtime_temp_inv, dtt_2dxi, dtt_2deta

!-------- Term abbreviations

at7 = 2.0_dp/RHO*dtime_temp

atr1 = KAPPA_R/(RHO_C_R*H_R**2)*dtime_temp/(dzeta_r**2)

if (flag_aa_nonzero) then
   am1 = aa*BETA*dzeta_c/(ea-1.0_dp)
else
   am1 = BETA*dzeta_c
end if

if (flag_aa_nonzero) then
   acb1 = (ea-1.0_dp)/aa/dzeta_c
else
   acb1 = 1.0_dp/dzeta_c
end if

acb2 = KAPPA_R/H_R/dzeta_r
acb3 = RHO*G
acb4 = RHO*G

alb1 = H_R/KAPPA_R*dzeta_r

dtt_2dxi  = 0.5_dp*dtime_temp/dxi
dtt_2deta = 0.5_dp*dtime_temp/deta

dtime_temp_inv = 1.0_dp/dtime_temp

do kc=0, KCMAX

   if (flag_aa_nonzero) then

      at1(kc)   = (ea-1.0_dp)/(aa*eaz_c(kc))*dtime_temp/dzeta_c
      at2_1(kc) = (ea-1.0_dp)/(aa*eaz_c(kc))*dtime_temp/dzeta_c
      at2_2(kc) = (eaz_c(kc)-1.0_dp)/(aa*eaz_c(kc)) &
                  *dtime_temp/dzeta_c
      at3_1(kc) = (ea-1.0_dp)/(aa*eaz_c(kc))*dtime_temp/dzeta_c
      at3_2(kc) = (eaz_c(kc)-1.0_dp)/(aa*eaz_c(kc)) &
                  *dtime_temp/dzeta_c
      at4_1(kc) = (ea-1.0_dp)/(aa*eaz_c(kc))*dtime_temp/dzeta_c
      at4_2(kc) = (eaz_c(kc)-1.0_dp)/(aa*eaz_c(kc)) &
                  *dtime_temp/dzeta_c
      at5(kc)   = (ea-1.0_dp)/(RHO*aa*eaz_c(kc)) &
                  *dtime_temp/dzeta_c
      if (kc /= KCMAX) then
         at6(kc) = (ea-1.0_dp) &
                   /(aa*exp(aa*0.5_dp*(zeta_c(kc)+zeta_c(kc+1)))) &
                   /dzeta_c
      else
         at6(kc) = 0.0_dp
      end if
      ai1(kc) = AGEDIFF*(ea-1.0_dp)/(aa*eaz_c(kc)) &
                *dtime_temp/dzeta_c
      if (kc /= KCMAX) then
         ai2(kc) = (ea-1.0_dp) &
                   /(aa*exp(aa*0.5_dp*(zeta_c(kc)+zeta_c(kc+1)))) &
                   /dzeta_c
      else
         ai2(kc) = 0.0_dp
      end if
      aqtlde(kc) = (aa*eaz_c(kc))/(ea-1.0_dp)*dzeta_c/dtime_temp
      am3(kc)    = (aa*eaz_c(kc))/(ea-1.0_dp)*dzeta_c*BETA

   else

      at1(kc)   = dtime_temp/dzeta_c
      at2_1(kc) = dtime_temp/dzeta_c
      at2_2(kc) = zeta_c(kc) &
                  *dtime_temp/dzeta_c
      at3_1(kc) = dtime_temp/dzeta_c
      at3_2(kc) = zeta_c(kc) &
                  *dtime_temp/dzeta_c
      at4_1(kc) = dtime_temp/dzeta_c
      at4_2(kc) = zeta_c(kc) &
                  *dtime_temp/dzeta_c
      at5(kc)   = 1.0_dp/RHO &
                  *dtime_temp/dzeta_c
      if (kc /= KCMAX) then
         at6(kc) = 1.0_dp &
                   /dzeta_c
      else
         at6(kc) = 0.0_dp
      end if
      ai1(kc) = AGEDIFF &
                *dtime_temp/dzeta_c
      if (kc /= KCMAX) then
         ai2(kc) = 1.0_dp &
                   /dzeta_c
      else
         ai2(kc) = 0.0_dp
      end if
      aqtlde(kc) = dzeta_c/dtime_temp
      am3(kc)    = dzeta_c*BETA

   end if

end do

strain_heating_c = 0.0_dp   ! initialization,
strain_heating_t = 0.0_dp   ! purely diagnostic fields

!-------- Computation loop --------

do i=1, IMAX-1   ! skipping domain margins
do j=1, JMAX-1   ! skipping domain margins

   if (mask(j,i)==0) then   ! glaciated land

!  ------ Old vertical column cold

      if (n_cts(j,i) == -1) then

         n_cts_new(j,i)  = -1
         kc_cts_new(j,i) =  0
         zm_new(j,i)     = zb(j,i)
         H_c_new(j,i)    = H_c(j,i)
         H_t_new(j,i)    = 0.0_dp
         !$AD NOCHECKPOINT
         call calc_temp_enth_1(at1, at2_1, at2_2, at3_1, at3_2, &
                               at4_1, at4_2, at5, at6, at7, &
                               atr1, acb1, acb2, acb3, acb4, alb1, &
                               ai1, ai2, &
                               dtime_temp, dtt_2dxi, dtt_2deta, &
                               dtime_temp_inv, &
                               i, j)

!    ---- Check whether base has become temperate

         if (temp_c_new(0,j,i) > temp_c_m(0,j,i)-eps) then

            n_cts_new(j,i)  = 0
            kc_cts_new(j,i) = 0
            !$AD NOCHECKPOINT
            call calc_temp_enth_2(at1, at2_1, at2_2, at3_1, at3_2, &
                                  at4_1, at4_2, at5, at6, at7, atr1, alb1, &
                                  ai1, ai2, aqtlde, am3, &
                                  dtime_temp, dtt_2dxi, dtt_2deta, &
                                  dtime_temp_inv, &
                                  i, j)

         end if

!  ------ Old vertical column with temperate base

      else if (n_cts(j,i) == 0) then

         n_cts_new(j,i)  = 0
         kc_cts_new(j,i) = kc_cts(j,i)
         zm_new(j,i)     = zb(j,i)
         H_c_new(j,i)    = H_c(j,i)
         H_t_new(j,i)    = H_t(j,i)
         !$AD NOCHECKPOINT
         call calc_temp_enth_2(at1, at2_1, at2_2, at3_1, at3_2, &
                               at4_1, at4_2, at5, at6, at7, atr1, alb1, &
                               ai1, ai2, aqtlde, am3, &
                               dtime_temp, dtt_2dxi, dtt_2deta, &
                               dtime_temp_inv, &
                               i, j)

!    ---- Check whether temperate base becomes cold

         if ( (temp_c_new(1,j,i)-temp_c_new(0,j,i)) < (am1*H_c(j,i)) ) then

            n_cts_new(j,i)  = -1
            kc_cts_new(j,i) =  0
            !$AD NOCHECKPOINT
            call calc_temp_enth_1(at1, at2_1, at2_2, at3_1, at3_2, &
                                  at4_1, at4_2, at5, at6, at7, &
                                  atr1, acb1, acb2, acb3, acb4, alb1, &
                                  ai1, ai2, &
                                  dtime_temp, dtt_2dxi, dtt_2deta, &
                                  dtime_temp_inv, &
                                  i, j)

            if (temp_c_new(0,j,i) > temp_c_m(0,j,i)-eps) then

               n_cts_new(j,i)  = 0
               kc_cts_new(j,i) = 0
               !$AD NOCHECKPOINT
               call calc_temp_enth_2(at1, at2_1, at2_2, at3_1, at3_2, &
                                     at4_1, at4_2, at5, at6, at7, atr1, alb1, &
                                     ai1, ai2, aqtlde, am3, &
                                     dtime_temp, dtt_2dxi, dtt_2deta, &
                                     dtime_temp_inv, &
                                     i, j)

            end if

         end if

      end if

#if (MARGIN==3)

   else if (mask(j,i)==3) then   ! floating ice

      n_cts_new(j,i)  = -1
      kc_cts_new(j,i) =  0
      zm_new(j,i)     = zb(j,i)
      H_c_new(j,i)    = H_c(j,i)
      H_t_new(j,i)    = 0.0_dp
      !$AD NOCHECKPOINT
      call calc_temp_enth_ssa(at1, at2_1, at2_2, at3_1, at3_2, &
                              at4_1, at4_2, at5, at6, at7, atr1, alb1, &
                              ai1, ai2, &
                              dtime_temp, dtt_2dxi, dtt_2deta, &
                              dtime_temp_inv, &
                              i, j)

!  ------ Reset temperatures above melting to the melting point
!         and water contents above zero to zero
!         (should not occur, but just in case)

      do kc=0, KCMAX
         if (temp_c_new(kc,j,i) > temp_c_m(kc,j,i)) &
                    temp_c_new(kc,j,i) = temp_c_m(kc,j,i)
         if (omega_c_new(kc,j,i) > 0.0_dp) &
                    omega_c_new(kc,j,i) = 0.0_dp
      end do

#endif

   else   ! mask(j,i) == 1, 2 (ice-free land or sea point)

      n_cts_new(j,i)  = -1
      kc_cts_new(j,i) =  0
      zm_new(j,i)     = zb(j,i)
      H_c_new(j,i)    = H_c(j,i)
      H_t_new(j,i)    = 0.0_dp
      !$AD NOCHECKPOINT
      call calc_temp_enth_r(atr1, alb1, i, j)

   end if

end do
end do 

!-------- Extrapolate values on margins --------

!  ------ Lower left corner

i=0
j=0

if ( (mask(j,i) == 0).or.(mask(j,i) == 3) ) then
                              ! glaciated land or floating ice
   ii=i+1
   jj=j+1

   do kc=0, KCMAX
      enth_c_new(kc,j,i)  = enth_c_new(kc,jj,ii)
      temp_c_new(kc,j,i)  = temp_c_new(kc,jj,ii)
      omega_c_new(kc,j,i) = omega_c_new(kc,jj,ii)
      age_c_new(kc,j,i)   = age_c_new(kc,jj,ii)
   end do

   do kt=0, KTMAX
            ! redundant, lower (kt) ice layer
      enth_t_new(kt,j,i)  = enth_t_new(kt,jj,ii)
      omega_t_new(kt,j,i) = omega_t_new(kt,jj,ii)
      age_t_new(kt,j,i)   = age_t_new(kt,jj,ii)
   end do

   do kr=0, KRMAX
      temp_r_new(kr,j,i)  = temp_r_new(kr,jj,ii)
   end do

   n_cts_new(j,i)  = n_cts_new(jj,ii)
   kc_cts_new(j,i) = kc_cts_new(jj,ii)
   zm_new(j,i)     = zb(j,i)
   H_c_new(j,i)    = H_c(j,i)
   H_t_new(j,i)    = H_t(j,i)

else   ! mask(j,i) == 1, 2 (ice-free land or sea point)

   n_cts_new(j,i)  = -1
   kc_cts_new(j,i) =  0
   zm_new(j,i)     = zb(j,i)
   H_c_new(j,i)    = H_c(j,i)
   H_t_new(j,i)    = H_t(j,i)

   call calc_temp_enth_r(atr1, alb1, i, j)

end if

!  ------ Lower right corner

i=IMAX
j=0

if ( (mask(j,i) == 0).or.(mask(j,i) == 3) ) then
                              ! glaciated land or floating ice
   ii=i-1
   jj=j+1

   do kc=0, KCMAX
      enth_c_new(kc,j,i)  = enth_c_new(kc,jj,ii)
      temp_c_new(kc,j,i)  = temp_c_new(kc,jj,ii)
      omega_c_new(kc,j,i) = omega_c_new(kc,jj,ii)
      age_c_new(kc,j,i)   = age_c_new(kc,jj,ii)
   end do

   do kt=0, KTMAX
            ! redundant, lower (kt) ice layer
      enth_t_new(kt,j,i)  = enth_t_new(kt,jj,ii)
      omega_t_new(kt,j,i) = omega_t_new(kt,jj,ii)
      age_t_new(kt,j,i)   = age_t_new(kt,jj,ii)
   end do

   do kr=0, KRMAX
      temp_r_new(kr,j,i)  = temp_r_new(kr,jj,ii)
   end do

   n_cts_new(j,i)  = n_cts_new(jj,ii)
   kc_cts_new(j,i) = kc_cts_new(jj,ii)
   zm_new(j,i)     = zb(j,i)
   H_c_new(j,i)    = H_c(j,i)
   H_t_new(j,i)    = H_t(j,i)

else   ! mask(j,i) == 1, 2 (ice-free land or sea point)

   n_cts_new(j,i)  = -1
   kc_cts_new(j,i) =  0
   zm_new(j,i)     = zb(j,i)
   H_c_new(j,i)    = H_c(j,i)
   H_t_new(j,i)    = H_t(j,i)

   call calc_temp_enth_r(atr1, alb1, i, j)

end if

!  ------ Upper left corner

i=0
j=JMAX

if ( (mask(j,i) == 0).or.(mask(j,i) == 3) ) then
                              ! glaciated land or floating ice
   ii=i+1
   jj=j-1

   do kc=0, KCMAX
      enth_c_new(kc,j,i)  = enth_c_new(kc,jj,ii)
      temp_c_new(kc,j,i)  = temp_c_new(kc,jj,ii)
      omega_c_new(kc,j,i) = omega_c_new(kc,jj,ii)
      age_c_new(kc,j,i)   = age_c_new(kc,jj,ii)
   end do

   do kt=0, KTMAX
            ! redundant, lower (kt) ice layer
      enth_t_new(kt,j,i)  = enth_t_new(kt,jj,ii)
      omega_t_new(kt,j,i) = omega_t_new(kt,jj,ii)
      age_t_new(kt,j,i)   = age_t_new(kt,jj,ii)
   end do

   do kr=0, KRMAX
      temp_r_new(kr,j,i)  = temp_r_new(kr,jj,ii)
   end do

   n_cts_new(j,i)  = n_cts_new(jj,ii)
   kc_cts_new(j,i) = kc_cts_new(jj,ii)
   zm_new(j,i)     = zb(j,i)
   H_c_new(j,i)    = H_c(j,i)
   H_t_new(j,i)    = H_t(j,i)

else   ! mask(j,i) == 1, 2 (ice-free land or sea point)

   n_cts_new(j,i)  = -1
   kc_cts_new(j,i) =  0
   zm_new(j,i)     = zb(j,i)
   H_c_new(j,i)    = H_c(j,i)
   H_t_new(j,i)    = H_t(j,i)

   call calc_temp_enth_r(atr1, alb1, i, j)

end if

!  ------ Upper right corner

i=IMAX
j=JMAX

if ( (mask(j,i) == 0).or.(mask(j,i) == 3) ) then
                              ! glaciated land or floating ice
   ii=i-1
   jj=j-1

   do kc=0, KCMAX
      enth_c_new(kc,j,i)  = enth_c_new(kc,jj,ii)
      temp_c_new(kc,j,i)  = temp_c_new(kc,jj,ii)
      omega_c_new(kc,j,i) = omega_c_new(kc,jj,ii)
      age_c_new(kc,j,i)   = age_c_new(kc,jj,ii)
   end do

   do kt=0, KTMAX
            ! redundant, lower (kt) ice layer
      enth_t_new(kt,j,i)  = enth_t_new(kt,jj,ii)
      omega_t_new(kt,j,i) = omega_t_new(kt,jj,ii)
      age_t_new(kt,j,i)   = age_t_new(kt,jj,ii)
   end do

   do kr=0, KRMAX
      temp_r_new(kr,j,i)  = temp_r_new(kr,jj,ii)
   end do

   n_cts_new(j,i)  = n_cts_new(jj,ii)
   kc_cts_new(j,i) = kc_cts_new(jj,ii)
   zm_new(j,i)     = zb(j,i)
   H_c_new(j,i)    = H_c(j,i)
   H_t_new(j,i)    = H_t(j,i)

else   ! mask(j,i) == 1, 2 (ice-free land or sea point)

   n_cts_new(j,i)  = -1
   kc_cts_new(j,i) =  0
   zm_new(j,i)     = zb(j,i)
   H_c_new(j,i)    = H_c(j,i)
   H_t_new(j,i)    = H_t(j,i)

   call calc_temp_enth_r(atr1, alb1, i, j)

end if

!  ------ Lower and upper margins

do i=1, IMAX-1

!    ---- Lower margin

   j=0

   if ( (mask(j,i) == 0).or.(mask(j,i) == 3) ) then
                                 ! glaciated land or floating ice
      ii=i
      jj=j+1

      do kc=0, KCMAX
         enth_c_new(kc,j,i)  = enth_c_new(kc,jj,ii)
         temp_c_new(kc,j,i)  = temp_c_new(kc,jj,ii)
         omega_c_new(kc,j,i) = omega_c_new(kc,jj,ii)
         age_c_new(kc,j,i)   = age_c_new(kc,jj,ii)
      end do

      do kt=0, KTMAX
               ! redundant, lower (kt) ice layer
         enth_t_new(kt,j,i)  = enth_t_new(kt,jj,ii)
         omega_t_new(kt,j,i) = omega_t_new(kt,jj,ii)
         age_t_new(kt,j,i)   = age_t_new(kt,jj,ii)
      end do

      do kr=0, KRMAX
         temp_r_new(kr,j,i)  = temp_r_new(kr,jj,ii)
      end do

      n_cts_new(j,i)  = n_cts_new(jj,ii)
      kc_cts_new(j,i) = kc_cts_new(jj,ii)
      zm_new(j,i)     = zb(j,i)
      H_c_new(j,i)    = H_c(j,i)
      H_t_new(j,i)    = H_t(j,i)

   else   ! mask(j,i) == 1, 2 (ice-free land or sea point)

      n_cts_new(j,i)  = -1
      kc_cts_new(j,i) =  0
      zm_new(j,i)     = zb(j,i)
      H_c_new(j,i)    = H_c(j,i)
      H_t_new(j,i)    = H_t(j,i)

      call calc_temp_enth_r(atr1, alb1, i, j)

   end if

!    ---- Upper margin

   j=JMAX

   if ( (mask(j,i) == 0).or.(mask(j,i) == 3) ) then
                                 ! glaciated land or floating ice
      ii=i
      jj=j-1

      do kc=0, KCMAX
         enth_c_new(kc,j,i)  = enth_c_new(kc,jj,ii)
         temp_c_new(kc,j,i)  = temp_c_new(kc,jj,ii)
         omega_c_new(kc,j,i) = omega_c_new(kc,jj,ii)
         age_c_new(kc,j,i)   = age_c_new(kc,jj,ii)
      end do

      do kt=0, KTMAX
               ! redundant, lower (kt) ice layer
         enth_t_new(kt,j,i)  = enth_t_new(kt,jj,ii)
         omega_t_new(kt,j,i) = omega_t_new(kt,jj,ii)
         age_t_new(kt,j,i)   = age_t_new(kt,jj,ii)
      end do

      do kr=0, KRMAX
         temp_r_new(kr,j,i)  = temp_r_new(kr,jj,ii)
      end do

      n_cts_new(j,i)  = n_cts_new(jj,ii)
      kc_cts_new(j,i) = kc_cts_new(jj,ii)
      zm_new(j,i)     = zb(j,i)
      H_c_new(j,i)    = H_c(j,i)
      H_t_new(j,i)    = H_t(j,i)

   else   ! mask(j,i) == 1, 2 (ice-free land or sea point)

      n_cts_new(j,i)  = -1
      kc_cts_new(j,i) =  0
      zm_new(j,i)     = zb(j,i)
      H_c_new(j,i)    = H_c(j,i)
      H_t_new(j,i)    = H_t(j,i)

      call calc_temp_enth_r(atr1, alb1, i, j)

   end if

end do

!  ------ Left and right margins

do j=1, JMAX-1

!    ---- Left margin

   i=0

   if ( (mask(j,i) == 0).or.(mask(j,i) == 3) ) then
                                 ! glaciated land or floating ice
      ii=i+1
      jj=j

      do kc=0, KCMAX
         enth_c_new(kc,j,i)  = enth_c_new(kc,jj,ii)
         temp_c_new(kc,j,i)  = temp_c_new(kc,jj,ii)
         omega_c_new(kc,j,i) = omega_c_new(kc,jj,ii)
         age_c_new(kc,j,i)   = age_c_new(kc,jj,ii)
      end do

      do kt=0, KTMAX
               ! redundant, lower (kt) ice layer
         enth_t_new(kt,j,i)  = enth_t_new(kt,jj,ii)
         omega_t_new(kt,j,i) = omega_t_new(kt,jj,ii)
         age_t_new(kt,j,i)   = age_t_new(kt,jj,ii)
      end do

      do kr=0, KRMAX
         temp_r_new(kr,j,i)  = temp_r_new(kr,jj,ii)
      end do

      n_cts_new(j,i)  = n_cts_new(jj,ii)
      kc_cts_new(j,i) = kc_cts_new(jj,ii)
      zm_new(j,i)     = zb(j,i)
      H_c_new(j,i)    = H_c(j,i)
      H_t_new(j,i)    = H_t(j,i)

   else   ! mask(j,i) == 1, 2 (ice-free land or sea point)

      n_cts_new(j,i)  = -1
      kc_cts_new(j,i) =  0
      zm_new(j,i)     = zb(j,i)
      H_c_new(j,i)    = H_c(j,i)
      H_t_new(j,i)    = H_t(j,i)

      call calc_temp_enth_r(atr1, alb1, i, j)

   end if

!    ---- Right margin

   i=IMAX

   if ( (mask(j,i) == 0).or.(mask(j,i) == 3) ) then
                                 ! glaciated land or floating ice
      ii=i-1
      jj=j

      do kc=0, KCMAX
         enth_c_new(kc,j,i)  = enth_c_new(kc,jj,ii)
         temp_c_new(kc,j,i)  = temp_c_new(kc,jj,ii)
         omega_c_new(kc,j,i) = omega_c_new(kc,jj,ii)
         age_c_new(kc,j,i)   = age_c_new(kc,jj,ii)
      end do

      do kt=0, KTMAX
               ! redundant, lower (kt) ice layer
         enth_t_new(kt,j,i)  = enth_t_new(kt,jj,ii)
         omega_t_new(kt,j,i) = omega_t_new(kt,jj,ii)
         age_t_new(kt,j,i)   = age_t_new(kt,jj,ii)
      end do

      do kr=0, KRMAX
         temp_r_new(kr,j,i)  = temp_r_new(kr,jj,ii)
      end do

      n_cts_new(j,i)  = n_cts_new(jj,ii)
      kc_cts_new(j,i) = kc_cts_new(jj,ii)
      zm_new(j,i)     = zb(j,i)
      H_c_new(j,i)    = H_c(j,i)
      H_t_new(j,i)    = H_t(j,i)

   else   ! mask(j,i) == 1, 2 (ice-free land or sea point)

      n_cts_new(j,i)  = -1
      kc_cts_new(j,i) =  0
      zm_new(j,i)     = zb(j,i)
      H_c_new(j,i)    = H_c(j,i)
      H_t_new(j,i)    = H_t(j,i)

      call calc_temp_enth_r(atr1, alb1, i, j)

   end if

end do

end subroutine calc_temp_enth

!-------------------------------------------------------------------------------
!> Computation of temperature and age for a cold ice column with the
!! enthalpy method.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_1(at1, at2_1, at2_2, at3_1, at3_2, &
                            at4_1, at4_2, at5, at6, at7, &
                            atr1, acb1, acb2, acb3, acb4, alb1, &
                            ai1, ai2, &
                            dtime_temp, dtt_2dxi, dtt_2deta, &
                            dtime_temp_inv, &
                            i, j)

#if !defined(ALLOW_TAPENADE) /* Normal */
use sico_maths_m, only : tri_sle
#else /* Tapenade */
use sico_maths_m
#endif /* Normal vs. Tapenade */

use enth_temp_omega_m, only : temp_fct_enth

implicit none

integer(i4b), intent(in) :: i, j

real(dp), intent(in) :: at1(0:KCMAX), at2_1(0:KCMAX), at2_2(0:KCMAX), &
                        at3_1(0:KCMAX), at3_2(0:KCMAX), &
                        at4_1(0:KCMAX), at4_2(0:KCMAX), &
                        at5(0:KCMAX), at6(0:KCMAX), at7, &
                        ai1(0:KCMAX), ai2(0:KCMAX), &
                        atr1, acb1, acb2, acb3, acb4, alb1
real(dp), intent(in) :: dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv

integer(i4b) :: kc, kt, kr
real(dp) :: ct1(0:KCMAX), ct2(0:KCMAX), ct3(0:KCMAX), ct4(0:KCMAX), &
            ce5(0:KCMAX), ce6(0:KCMAX), ce7(0:KCMAX), &
            ctr1, ccbe1, ccb2, ccb3, ccb4, clb1
real(dp) :: ct1_sg(0:KCMAX), ct2_sg(0:KCMAX), ct3_sg(0:KCMAX), &
            ct4_sg(0:KCMAX), adv_vert_sg(0:KCMAX), abs_adv_vert_sg(0:KCMAX)
real(dp) :: ci1(0:KCMAX), ci2(0:KCMAX)
real(dp) :: dtt_dxi, dtt_deta
real(dp) :: lgs_a0(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_a1(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_a2(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_x(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_b(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX)

!-------- Check for boundary points --------

if ((i == 0).or.(i == IMAX).or.(j == 0).or.(j == JMAX)) then
   errormsg = ' >>> calc_temp_enth_1: Boundary points not allowed.'
   call error(errormsg)
end if

!-------- Abbreviations --------

call calc_temp_enth_1_a(at1, at2_1, at2_2, at3_1, at3_2, &
                        at4_1, at4_2, at5, at6, at7, &
                        atr1, acb1, acb2, acb3, acb4, alb1, &
                        ai1, ai2, &
                        dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
                        i, j, &
                        ct1, ct2, ct3, ct4, ce5, ce6, ce7, &
                        ctr1, ccbe1, ccb2, ccb3, ccb4, clb1, &
                        ct1_sg, ct2_sg, ct3_sg, ct4_sg, &
                        adv_vert_sg, abs_adv_vert_sg, &
                        ci1, ci2, dtt_dxi, dtt_deta)

!-------- Computation of the bedrock temperature
!         (upper boundary condition: old temperature at the ice base) --------

!  ------ Set-up of the the equations

call calc_temp_enth_1_b(ctr1, clb1, i, j, &
                        lgs_a0, lgs_a1, lgs_a2, lgs_b)

!  ------ Solution of the system of linear equations

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, KRMAX)

!  ------ Assignment of the result (predictor values)

do kr=0, KRMAX
   temp_r_new(kr,j,i) = lgs_x(kr)
end do

!-------- Computation of the ice enthalpy --------

!  ------ Set-up of the the equations

call calc_temp_enth_1_c(ct1, ct2, ct3, ct4, ce5, ce6, ce7, &
                        ct1_sg, ct2_sg, ct3_sg, ct4_sg, &
                        ccbe1, ccb2, ccb3, ccb4, &
                        adv_vert_sg, abs_adv_vert_sg, &
                        dtime_temp, dtt_dxi, dtt_deta, &
                        dtt_2dxi, dtt_2deta, &
                        i, j, &
                        lgs_a0, lgs_a1, lgs_a2, lgs_b)

!  ------ Solution of the system of linear equations

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, KCMAX)

!  ------ Assignment of the result

do kc=0, KCMAX
   enth_c_new(kc,j,i)  = lgs_x(kc)
   temp_c_new(kc,j,i)  = temp_fct_enth(enth_c_new(kc,j,i), temp_c_m(kc,j,i))
   omega_c_new(kc,j,i) = 0.0_dp   ! solution is supposed to be for cold ice
end do

!-------- Water drainage from the non-existing temperate ice --------

Q_tld(j,i) = 0.0_dp

!-------- Set enthalpy and water content in the redundant,
!         lower (kt) ice layer to the value at the ice base --------

do kt=0, KTMAX
   enth_t_new(kt,j,i)  = enth_c_new(0,j,i)
   omega_t_new(kt,j,i) = omega_c_new(0,j,i)
end do

!-------- Computation of the age of ice --------

!  ------ Set-up of the the equations

call calc_temp_enth_1_d(ct1, ct2, ct3, ct4, ci1, ci2, &
                        ct1_sg, ct2_sg, ct3_sg, ct4_sg, &
                        adv_vert_sg, abs_adv_vert_sg, &
                        dtime_temp, dtt_dxi, dtt_deta, &
                        dtt_2dxi, dtt_2deta, &
                        i, j, &
                        lgs_a0, lgs_a1, lgs_a2, lgs_b)

!  ------ Solution of the system of linear equations

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, KCMAX)

!  ------ Assignment of the result,
!         restriction to interval [0, AGE_MAX yr] --------

do kc=0, KCMAX

   age_c_new(kc,j,i) = lgs_x(kc)

   if (age_c_new(kc,j,i) < (AGE_MIN*year2sec)) &
                           age_c_new(kc,j,i) = 0.0_dp
   if (age_c_new(kc,j,i) > (AGE_MAX*year2sec)) &
                           age_c_new(kc,j,i) = AGE_MAX*year2sec

end do

!-------- Age of the ice in the redundant, lower (kt) ice layer --------

do kt=0, KTMAX
   age_t_new(kt,j,i) = age_c_new(0,j,i)
end do

end subroutine calc_temp_enth_1

!-------------------------------------------------------------------------------
!> Computation of temperature and age for a cold ice column with the
!! enthalpy method: Abbreviations.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_1_a(at1, at2_1, at2_2, at3_1, at3_2, &
                              at4_1, at4_2, at5, at6, at7, &
                              atr1, acb1, acb2, acb3, acb4, alb1, &
                              ai1, ai2, &
                              dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
                              i, j, &
                              ct1, ct2, ct3, ct4, ce5, ce6, ce7, &
                              ctr1, ccbe1, ccb2, ccb3, ccb4, clb1, &
                              ct1_sg, ct2_sg, ct3_sg, ct4_sg, &
                              adv_vert_sg, abs_adv_vert_sg, &
                              ci1, ci2, dtt_dxi, dtt_deta)

use ice_material_properties_m, only : ratefac_c_t, kappa_val, c_val, &
                                      creep, viscosity

implicit none

integer(i4b), intent(in) :: i, j

real(dp),     intent(in) :: at1(0:KCMAX), &
                            at2_1(0:KCMAX), at2_2(0:KCMAX), &
                            at3_1(0:KCMAX), at3_2(0:KCMAX), &
                            at4_1(0:KCMAX), at4_2(0:KCMAX), &
                            at5(0:KCMAX), at6(0:KCMAX), at7, &
                            ai1(0:KCMAX), ai2(0:KCMAX), &
                            atr1, acb1, acb2, acb3, acb4, alb1
real(dp),     intent(in) :: dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv

real(dp),    intent(out) :: ct1(0:KCMAX), ct2(0:KCMAX), ct3(0:KCMAX), &
                            ct4(0:KCMAX), ce5(0:KCMAX), ce6(0:KCMAX), &
                            ce7(0:KCMAX), &
                            ctr1, ccbe1, ccb2, ccb3, ccb4, clb1
real(dp),    intent(out) :: ct1_sg(0:KCMAX), ct2_sg(0:KCMAX), &
                            ct3_sg(0:KCMAX), ct4_sg(0:KCMAX), &
                            adv_vert_sg(0:KCMAX), abs_adv_vert_sg(0:KCMAX)
real(dp),    intent(out) :: ci1(0:KCMAX), ci2(0:KCMAX)
real(dp),    intent(out) :: dtt_dxi, dtt_deta

integer(i4b) :: kc
real(dp) :: temp_c_help(0:KCMAX)

!-------- Initialisation --------

ct1             = 0.0_dp
ct2             = 0.0_dp
ct3             = 0.0_dp
ct4             = 0.0_dp
ce5             = 0.0_dp
ce6             = 0.0_dp
ce7             = 0.0_dp
ctr1            = 0.0_dp
clb1            = 0.0_dp
ct1_sg          = 0.0_dp
ct2_sg          = 0.0_dp
ct3_sg          = 0.0_dp
ct4_sg          = 0.0_dp
adv_vert_sg     = 0.0_dp
abs_adv_vert_sg = 0.0_dp
ci1             = 0.0_dp
ci2             = 0.0_dp
dtt_dxi         = 0.0_dp
dtt_deta        = 0.0_dp

!-------- Actual computation --------

ctr1 = atr1

ccbe1 = acb1 &
        *kappa_val(temp_c(0,j,i)) &
        /(c_val(temp_c(0,j,i))*H_c(j,i))
ccb2  = acb2

#if (DYNAMICS==2)
if (.not.flag_shelfy_stream(j,i)) then
#endif

   ccb3 = acb3*0.5_dp*(vx_t(0,j,i)+vx_t(0,j,i-1)) &
              *H_c(j,i)*dzs_dxi_g(j,i)
   ccb4 = acb4*0.5_dp*(vy_t(0,j,i)+vy_t(0,j-1,i)) &
              *H_c(j,i)*dzs_deta_g(j,i)

#if (DYNAMICS==2)
else   ! flag_shelfy_stream(j,i) == .true.

#if !defined(ALLOW_TAPENADE) /* Normal */
   ccb3 = -c_drag(j,i) &
           * sqrt(vx_b_g(j,i)**2  &
                 +vy_b_g(j,i)**2) &
                           **(1.0_dp+p_weert_inv(j,i))
#else /* Tapenade: guarding against undifferentiable sqrt(0) */
   if ((vx_b_g(j,i)**2+vy_b_g(j,i)**2) == 0) then
   ccb3 = 0.0_dp 
   else
   ccb3 = -c_drag(j,i) &
           * sqrt(vx_b_g(j,i)**2  &
                 +vy_b_g(j,i)**2) &
                           **(1.0_dp+p_weert_inv(j,i))
   end if
#endif /* Normal vs. Tapenade */
   ccb4 = 0.0_dp

end if
#endif

clb1 = alb1*q_geo(j,i)

#if (ADV_VERT==1)

do kc=1, KCMAX-1
   ct1(kc) = at1(kc)/H_c(j,i)*0.5_dp*(vz_c(kc,j,i)+vz_c(kc-1,j,i))
end do

kc=0
ct1_sg(kc) = 0.5_dp*(at1(kc)+at1(kc+1))/H_c(j,i)*vz_c(kc,j,i)
             ! only needed for kc=0 ...
kc=KCMAX-1
ct1_sg(kc) = 0.5_dp*(at1(kc)+at1(kc+1))/H_c(j,i)*vz_c(kc,j,i)
             ! ... and kc=KCMAX-1

#elif (ADV_VERT==2 || ADV_VERT==3)

do kc=0, KCMAX-1
   ct1_sg(kc) = 0.5_dp*(at1(kc)+at1(kc+1))/H_c(j,i)*vz_c(kc,j,i)
end do

#endif

do kc=0, KCMAX

   ct2(kc) = ( at2_1(kc)*dzm_dtau(j,i) &
           +at2_2(kc)*dH_c_dtau(j,i) )/H_c(j,i)
   ct3(kc) = ( at3_1(kc)*dzm_dxi_g(j,i) &
           +at3_2(kc)*dH_c_dxi_g(j,i) )/H_c(j,i) &
          *0.5_dp*(vx_c(kc,j,i)+vx_c(kc,j,i-1)) *insq_g11_g(j,i)
   ct4(kc) = ( at4_1(kc)*dzm_deta_g(j,i) &
            +at4_2(kc)*dH_c_deta_g(j,i) )/H_c(j,i) &
          *0.5_dp*(vy_c(kc,j,i)+vy_c(kc,j-1,i)) *insq_g22_g(j,i)
   ce5(kc) = at5(kc)/H_c(j,i)

#if (DYNAMICS==2)
   if (.not.flag_shelfy_stream(j,i)) then
#endif
      ce7(kc) = at7 &
                *enh_c(kc,j,i) &
                *ratefac_c_t(temp_c(kc,j,i), omega_c(kc,j,i), temp_c_m(kc,j,i)) &
                *creep(sigma_c(kc,j,i)) &
                *sigma_c(kc,j,i)*sigma_c(kc,j,i)
#if (DYNAMICS==2)
   else
      ce7(kc) = 2.0_dp*at7 &
                *viscosity(de_c(kc,j,i), &
                           temp_c(kc,j,i), temp_c_m(kc,j,i), 0.0_dp, &
                           enh_c(kc,j,i), 0) &
                *de_c(kc,j,i)**2
   end if
#endif

   strain_heating_c(kc,j,i) = ce7(kc)*dtime_temp_inv

   ci1(kc) = ai1(kc)/H_c(j,i)

end do

#if (ADV_VERT==1)

kc=0
ct2_sg(kc) = 0.5_dp*(ct2(kc)+ct2(kc+1))
ct3_sg(kc) = 0.5_dp*(ct3(kc)+ct3(kc+1))
ct4_sg(kc) = 0.5_dp*(ct4(kc)+ct4(kc+1))
adv_vert_sg(kc) = ct1_sg(kc)-ct2_sg(kc)-ct3_sg(kc)-ct4_sg(kc)
abs_adv_vert_sg(kc) = abs(adv_vert_sg(kc))   ! only needed for kc=0 ...
kc=KCMAX-1
ct2_sg(kc) = 0.5_dp*(ct2(kc)+ct2(kc+1))
ct3_sg(kc) = 0.5_dp*(ct3(kc)+ct3(kc+1))
ct4_sg(kc) = 0.5_dp*(ct4(kc)+ct4(kc+1))
adv_vert_sg(kc) = ct1_sg(kc)-ct2_sg(kc)-ct3_sg(kc)-ct4_sg(kc)
abs_adv_vert_sg(kc) = abs(adv_vert_sg(kc))   ! ... and kc=KCMAX-1

#elif (ADV_VERT==2 || ADV_VERT==3)

do kc=0, KCMAX-1
   ct2_sg(kc) = 0.5_dp*(ct2(kc)+ct2(kc+1))
   ct3_sg(kc) = 0.5_dp*(ct3(kc)+ct3(kc+1))
   ct4_sg(kc) = 0.5_dp*(ct4(kc)+ct4(kc+1))
   adv_vert_sg(kc) = ct1_sg(kc)-ct2_sg(kc)-ct3_sg(kc)-ct4_sg(kc)
   abs_adv_vert_sg(kc) = abs(adv_vert_sg(kc))
end do

#endif

do kc=0, KCMAX-1
   temp_c_help(kc) = 0.5_dp*(temp_c(kc,j,i)+temp_c(kc+1,j,i))
   ce6(kc) = at6(kc) &
             *kappa_val(temp_c_help(kc))/(c_val(temp_c_help(kc))*H_c(j,i))
   ci2(kc) = ai2(kc)/H_c(j,i)
end do

#if (ADV_HOR==3)
dtt_dxi  = 2.0_dp*dtt_2dxi
dtt_deta = 2.0_dp*dtt_2deta
#endif

end subroutine calc_temp_enth_1_a

!-------------------------------------------------------------------------------
!> Computation of temperature and age for a cold ice column with the
!! enthalpy method:
!! Set-up of the equations for the bedrock temperature.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_1_b(ctr1, clb1, i, j, &
                              lgs_a0, lgs_a1, lgs_a2, lgs_b)

implicit none

integer(i4b), intent(in) :: i, j

real(dp),     intent(in) :: ctr1, clb1

real(dp),    intent(out) :: lgs_a0(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
                            lgs_a1(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
                            lgs_a2(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
                            lgs_b(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX)

integer(i4b) :: kr

!-------- Initialisation --------

lgs_a0 = 0.0_dp
lgs_a1 = 0.0_dp
lgs_a2 = 0.0_dp
lgs_b  = 0.0_dp

!-------- Actual computation --------

kr=0
lgs_a1(kr) =  1.0_dp
lgs_a2(kr) = -1.0_dp
lgs_b(kr)  = clb1

#if (Q_LITHO==1)
!   (coupled heat-conducting bedrock)

do kr=1, KRMAX-1
   lgs_a0(kr) = -ctr1
   lgs_a1(kr) = 1.0_dp + 2.0_dp*ctr1
   lgs_a2(kr) = -ctr1
   lgs_b(kr)  = temp_r(kr,j,i)
end do

#elif (Q_LITHO==0)
!   (no coupled heat-conducting bedrock)

do kr=1, KRMAX-1
   lgs_a0(kr) =  1.0_dp
   lgs_a1(kr) =  0.0_dp
   lgs_a2(kr) = -1.0_dp
   lgs_b(kr)  =  2.0_dp*clb1
end do

#endif

kr=KRMAX
lgs_a0(kr) = 0.0_dp
lgs_a1(kr) = 1.0_dp
lgs_b(kr)  = temp_c(0,j,i)

end subroutine calc_temp_enth_1_b

!-------------------------------------------------------------------------------
!> Computation of temperature and age for a cold ice column with the
!! enthalpy method:
!! Set-up of the equations for the ice enthalpy.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_1_c(ct1, ct2, ct3, ct4, ce5, ce6, ce7, &
                              ct1_sg, ct2_sg, ct3_sg, ct4_sg, &
                              ccbe1, ccb2, ccb3, ccb4, &
                              adv_vert_sg, abs_adv_vert_sg, &
                              dtime_temp, dtt_dxi, dtt_deta, &
                              dtt_2dxi, dtt_2deta, &
                              i, j, &
                              lgs_a0, lgs_a1, lgs_a2, lgs_b)

use enth_temp_omega_m, only : enth_fct_temp_omega

implicit none

integer(i4b), intent(in) :: i, j

real(dp),     intent(in) :: ct1(0:KCMAX), ct2(0:KCMAX), ct3(0:KCMAX), &
                            ct4(0:KCMAX), ce5(0:KCMAX), ce6(0:KCMAX), &
                            ce7(0:KCMAX)
real(dp),     intent(in) :: ct1_sg(0:KCMAX), ct2_sg(0:KCMAX), &
                            ct3_sg(0:KCMAX), ct4_sg(0:KCMAX), &
                            ccbe1, ccb2, ccb3, ccb4, &
                            adv_vert_sg(0:KCMAX), abs_adv_vert_sg(0:KCMAX)
real(dp),     intent(in) :: dtime_temp, dtt_dxi, dtt_deta, dtt_2dxi, dtt_2deta

real(dp),    intent(out) :: lgs_a0(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
                            lgs_a1(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
                            lgs_a2(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
                            lgs_b(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX)

integer(i4b) :: kc, kr
real(dp) :: vx_c_help, vy_c_help
real(dp) :: adv_vert_help

!-------- Initialisation --------

lgs_a0 = 0.0_dp
lgs_a1 = 0.0_dp
lgs_a2 = 0.0_dp
lgs_b  = 0.0_dp

!-------- Actual computation --------

kr=KRMAX
kc=0
lgs_a1(kc) = -ccbe1
lgs_a2(kc) =  ccbe1
lgs_b(kc)  =  ccb2*(temp_r_new(kr,j,i)-temp_r_new(kr-1,j,i)) + ccb3 + ccb4

do kc=1, KCMAX-1

#if (ADV_VERT==1)

   lgs_a0(kc) = -0.5_dp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ce5(kc)*ce6(kc-1)
   lgs_a1(kc) = 1.0_dp+ce5(kc)*(ce6(kc)+ce6(kc-1))
   lgs_a2(kc) = 0.5_dp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ce5(kc)*ce6(kc)

#elif (ADV_VERT==2)

   lgs_a0(kc) &
         = -0.5_dp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -ce5(kc)*ce6(kc-1)
   lgs_a1(kc) &
         = 1.0_dp &
           +0.5_dp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -0.5_dp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           +ce5(kc)*(ce6(kc)+ce6(kc-1))
   lgs_a2(kc) &
         =  0.5_dp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           -ce5(kc)*ce6(kc)

#elif (ADV_VERT==3)

   adv_vert_help = 0.5_dp*(adv_vert_sg(kc)+adv_vert_sg(kc-1))

   lgs_a0(kc) &
         = -max(adv_vert_help, 0.0_dp) &
           -ce5(kc)*ce6(kc-1)
   lgs_a1(kc) &
         = 1.0_dp &
           +max(adv_vert_help, 0.0_dp)-min(adv_vert_help, 0.0_dp) &
           +ce5(kc)*(ce6(kc)+ce6(kc-1))
   lgs_a2(kc) &
         =  min(adv_vert_help, 0.0_dp) &
           -ce5(kc)*ce6(kc)

#endif

#if (ADV_HOR==2)

   lgs_b(kc) = enth_c(kc,j,i) + ce7(kc) &
       - ( dtt_2dxi* &
             ( (vx_c(kc,j,i)-abs(vx_c(kc,j,i))) &
               *(enth_c(kc,j,i+1)-enth_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +(vx_c(kc,j,i-1)+abs(vx_c(kc,j,i-1))) &
               *(enth_c(kc,j,i)-enth_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_2deta* &
             ( (vy_c(kc,j,i)-abs(vy_c(kc,j,i))) &
               *(enth_c(kc,j+1,i)-enth_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +(vy_c(kc,j-1,i)+abs(vy_c(kc,j-1,i))) &
               *(enth_c(kc,j,i)-enth_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#elif (ADV_HOR==3)

   vx_c_help = 0.5_dp*(vx_c(kc,j,i)+vx_c(kc,j,i-1))
   vy_c_help = 0.5_dp*(vy_c(kc,j,i)+vy_c(kc,j-1,i))

   lgs_b(kc) = enth_c(kc,j,i) + ce7(kc) &
       - ( dtt_dxi* &
             ( min(vx_c_help, 0.0_dp) &
               *(enth_c(kc,j,i+1)-enth_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +max(vx_c_help, 0.0_dp) &
               *(enth_c(kc,j,i)-enth_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_deta* &
             ( min(vy_c_help, 0.0_dp) &
               *(enth_c(kc,j+1,i)-enth_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +max(vy_c_help, 0.0_dp) &
               *(enth_c(kc,j,i)-enth_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#endif

end do

kc=KCMAX
lgs_a0(kc) = 0.0_dp
lgs_a1(kc) = 1.0_dp
lgs_b(kc)  = enth_fct_temp_omega(temp_s(j,i), 0.0_dp)
                                ! zero water content at the ice surface

end subroutine calc_temp_enth_1_c

!-------------------------------------------------------------------------------
!> Computation of temperature and age for a cold ice column with the
!! enthalpy method:
!! Set-up of the equations for the age of ice.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_1_d(ct1, ct2, ct3, ct4, ci1, ci2, &
                              ct1_sg, ct2_sg, ct3_sg, ct4_sg, &
                              adv_vert_sg, abs_adv_vert_sg, &
                              dtime_temp, dtt_dxi, dtt_deta, &
                              dtt_2dxi, dtt_2deta, &
                              i, j, &
                              lgs_a0, lgs_a1, lgs_a2, lgs_b)

implicit none

integer(i4b), intent(in) :: i, j

real(dp),     intent(in) :: ct1(0:KCMAX), ct2(0:KCMAX), ct3(0:KCMAX), &
                            ct4(0:KCMAX), ci1(0:KCMAX), ci2(0:KCMAX)
real(dp),     intent(in) :: ct1_sg(0:KCMAX), ct2_sg(0:KCMAX), &
                            ct3_sg(0:KCMAX), ct4_sg(0:KCMAX), &
                            adv_vert_sg(0:KCMAX), abs_adv_vert_sg(0:KCMAX)
real(dp),     intent(in) :: dtime_temp, dtt_dxi, dtt_deta, dtt_2dxi, dtt_2deta

real(dp),    intent(out) :: lgs_a0(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
                            lgs_a1(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
                            lgs_a2(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
                            lgs_b(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX)

integer(i4b) :: kc
real(dp) :: vx_c_help, vy_c_help
real(dp) :: adv_vert_help

!-------- Initialisation --------

lgs_a0 = 0.0_dp
lgs_a1 = 0.0_dp
lgs_a2 = 0.0_dp
lgs_b  = 0.0_dp

!-------- Actual computation --------

kc=0                                                 ! adv_vert_sg(0) <= 0
lgs_a1(kc) = 1.0_dp - min(adv_vert_sg(kc), 0.0_dp)   ! (directed downward)
lgs_a2(kc) = min(adv_vert_sg(kc), 0.0_dp)            ! assumed/enforced

#if (ADV_HOR==2)

lgs_b(kc) = age_c(kc,j,i) + dtime_temp &
       - ( dtt_2dxi* &
             ( (vx_c(kc,j,i)-abs(vx_c(kc,j,i))) &
               *(age_c(kc,j,i+1)-age_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +(vx_c(kc,j,i-1)+abs(vx_c(kc,j,i-1))) &
               *(age_c(kc,j,i)-age_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_2deta* &
             ( (vy_c(kc,j,i)-abs(vy_c(kc,j,i))) &
               *(age_c(kc,j+1,i)-age_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +(vy_c(kc,j-1,i)+abs(vy_c(kc,j-1,i))) &
               *(age_c(kc,j,i)-age_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#elif (ADV_HOR==3)

vx_c_help = 0.5_dp*(vx_c(kc,j,i)+vx_c(kc,j,i-1))
vy_c_help = 0.5_dp*(vy_c(kc,j,i)+vy_c(kc,j-1,i))

lgs_b(kc) = age_c(kc,j,i) + dtime_temp &
       - ( dtt_dxi* &
             ( min(vx_c_help, 0.0_dp) &
               *(age_c(kc,j,i+1)-age_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +max(vx_c_help, 0.0_dp) &
               *(age_c(kc,j,i)-age_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_deta* &
             ( min(vy_c_help, 0.0_dp) &
               *(age_c(kc,j+1,i)-age_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +max(vy_c_help, 0.0_dp) &
               *(age_c(kc,j,i)-age_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#endif

do kc=1, KCMAX-1

#if (ADV_VERT==1)

   lgs_a0(kc) = -0.5_dp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ci1(kc)*ci2(kc-1)
   lgs_a1(kc) = 1.0_dp+ci1(kc)*(ci2(kc)+ci2(kc-1))
   lgs_a2(kc) = 0.5_dp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ci1(kc)*ci2(kc)

#elif (ADV_VERT==2)

   lgs_a0(kc) = -0.5_dp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1))
   lgs_a1(kc) =  1.0_dp &
                +0.5_dp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
                -0.5_dp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  )
   lgs_a2(kc) =  0.5_dp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  )

#elif (ADV_VERT==3)

   adv_vert_help = 0.5_dp*(adv_vert_sg(kc)+adv_vert_sg(kc-1))

   lgs_a0(kc) = -max(adv_vert_help, 0.0_dp)
   lgs_a1(kc) =  1.0_dp &
                +max(adv_vert_help, 0.0_dp)-min(adv_vert_help, 0.0_dp)
   lgs_a2(kc) =  min(adv_vert_help, 0.0_dp)

#endif

#if (ADV_HOR==2)

   lgs_b(kc) = age_c(kc,j,i) + dtime_temp &
       - ( dtt_2dxi* &
             ( (vx_c(kc,j,i)-abs(vx_c(kc,j,i))) &
               *(age_c(kc,j,i+1)-age_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +(vx_c(kc,j,i-1)+abs(vx_c(kc,j,i-1))) &
               *(age_c(kc,j,i)-age_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_2deta* &
             ( (vy_c(kc,j,i)-abs(vy_c(kc,j,i))) &
               *(age_c(kc,j+1,i)-age_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +(vy_c(kc,j-1,i)+abs(vy_c(kc,j-1,i))) &
               *(age_c(kc,j,i)-age_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#elif (ADV_HOR==3)

   vx_c_help = 0.5_dp*(vx_c(kc,j,i)+vx_c(kc,j,i-1))
   vy_c_help = 0.5_dp*(vy_c(kc,j,i)+vy_c(kc,j-1,i))

   lgs_b(kc) = age_c(kc,j,i) + dtime_temp &
       - ( dtt_dxi* &
             ( min(vx_c_help, 0.0_dp) &
               *(age_c(kc,j,i+1)-age_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +max(vx_c_help, 0.0_dp) &
               *(age_c(kc,j,i)-age_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_deta* &
             ( min(vy_c_help, 0.0_dp) &
               *(age_c(kc,j+1,i)-age_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +max(vy_c_help, 0.0_dp) &
               *(age_c(kc,j,i)-age_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#endif

end do

kc=KCMAX
if (as_perp(j,i) >= 0.0_dp) then
   lgs_a0(kc) = 0.0_dp
   lgs_a1(kc) = 1.0_dp
   lgs_b(kc)  = 0.0_dp
else
   lgs_a0(kc) = -max(adv_vert_sg(kc-1), 0.0_dp)
   lgs_a1(kc) = 1.0_dp + max(adv_vert_sg(kc-1), 0.0_dp)
                       ! adv_vert_sg(KCMAX-1) >= 0 (directed upward)
                       ! assumed/enforced
#if (ADV_HOR==2)

   lgs_b(kc) = age_c(kc,j,i) + dtime_temp &
       - ( dtt_2dxi* &
             ( (vx_c(kc,j,i)-abs(vx_c(kc,j,i))) &
               *(age_c(kc,j,i+1)-age_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +(vx_c(kc,j,i-1)+abs(vx_c(kc,j,i-1))) &
               *(age_c(kc,j,i)-age_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_2deta* &
             ( (vy_c(kc,j,i)-abs(vy_c(kc,j,i))) &
               *(age_c(kc,j+1,i)-age_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +(vy_c(kc,j-1,i)+abs(vy_c(kc,j-1,i))) &
               *(age_c(kc,j,i)-age_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#elif (ADV_HOR==3)

   vx_c_help = 0.5_dp*(vx_c(kc,j,i)+vx_c(kc,j,i-1))
   vy_c_help = 0.5_dp*(vy_c(kc,j,i)+vy_c(kc,j-1,i))

   lgs_b(kc) = age_c(kc,j,i) + dtime_temp &
       - ( dtt_dxi* &
             ( min(vx_c_help, 0.0_dp) &
               *(age_c(kc,j,i+1)-age_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +max(vx_c_help, 0.0_dp) &
               *(age_c(kc,j,i)-age_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_deta* &
             ( min(vy_c_help, 0.0_dp) &
               *(age_c(kc,j+1,i)-age_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +max(vy_c_help, 0.0_dp) &
               *(age_c(kc,j,i)-age_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#endif

end if

end subroutine calc_temp_enth_1_d

!-------------------------------------------------------------------------------
!> Computation of temperature and age for an ice column with a temperate base
!! with the enthalpy method.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_2(at1, at2_1, at2_2, at3_1, at3_2, &
                            at4_1, at4_2, at5, at6, at7, atr1, alb1, &
                            ai1, ai2, aqtlde, am3, &
                            dtime_temp, dtt_2dxi, dtt_2deta, &
                            dtime_temp_inv, &
                            i, j)

#if !defined(ALLOW_TAPENADE) /* Normal */
use sico_maths_m, only : tri_sle
#else /* Tapenade */
use sico_maths_m
#endif /* Normal vs. Tapenade */

use enth_temp_omega_m, only : enth_fct_temp_omega, &
                              temp_fct_enth, omega_fct_enth

implicit none

integer(i4b), intent(in) :: i, j

real(dp), intent(in) :: at1(0:KCMAX), at2_1(0:KCMAX), at2_2(0:KCMAX), &
                        at3_1(0:KCMAX), at3_2(0:KCMAX), &
                        at4_1(0:KCMAX), at4_2(0:KCMAX), &
                        at5(0:KCMAX), at6(0:KCMAX), at7, &
                        ai1(0:KCMAX), ai2(0:KCMAX), &
                        atr1, alb1, aqtlde(0:KCMAX), am3(0:KCMAX)
real(dp), intent(in) :: dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv

integer(i4b) :: kc, kt, kr
real(dp) :: ct1(0:KCMAX), ct2(0:KCMAX), ct3(0:KCMAX), ct4(0:KCMAX), &
            ce5(0:KCMAX), ce6(0:KCMAX), ce7(0:KCMAX), ctr1, clb1
real(dp) :: ct1_sg(0:KCMAX), ct2_sg(0:KCMAX), ct3_sg(0:KCMAX), &
            ct4_sg(0:KCMAX), adv_vert_sg(0:KCMAX), abs_adv_vert_sg(0:KCMAX)
real(dp) :: ci1(0:KCMAX), ci2(0:KCMAX)
real(dp) :: cqtlde(0:KCMAX), cm3(0:KCMAX)
real(dp) :: dtt_dxi, dtt_deta
real(dp) :: temp_c_val(0:KCMAX), omega_c_val(0:KCMAX)
real(dp) :: lgs_a0(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_a1(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_a2(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_x(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_b(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX)

real(dp), parameter :: eps_omega=1.0e-12_dp

#if defined(ALLOW_TAPENADE) /* Tapenade */
logical :: kcdone
#endif /* Tapenade */

!-------- Check for boundary points --------

if ((i == 0).or.(i == IMAX).or.(j == 0).or.(j == JMAX)) then
   errormsg = ' >>> calc_temp_enth_2: Boundary points not allowed.'
   call error(errormsg)
end if

!-------- Abbreviations --------

call calc_temp_enth_2_a1(at1, at2_1, at2_2, at3_1, at3_2, &
                         at4_1, at4_2, at5, atr1, alb1, &
                         ai1, aqtlde, &
                         dtime_temp, dtt_2dxi, dtt_2deta, i, j, &
                         ct1, ct2, ct3, ct4, ce5, ctr1, clb1, &
                         ct1_sg, ct2_sg, ct3_sg, ct4_sg, &
                         adv_vert_sg, abs_adv_vert_sg, &
                         ci1, cqtlde, dtt_dxi, dtt_deta)

do kc=0, KCMAX
   temp_c_val(kc)  = temp_c(kc,j,i)
   omega_c_val(kc) = omega_c(kc,j,i)
end do

call calc_temp_enth_2_a2(at6, at7, ai2, am3, temp_c_val, omega_c_val, &
                         dtime_temp_inv, &
                         i, j, ce6, ce7, ci2, cm3)

!-------- Computation of the bedrock temperature --------

!  ------ Set-up of the the equations

call calc_temp_enth_2_b(ctr1, clb1, i, j, &
                        lgs_a0, lgs_a1, lgs_a2, lgs_b)

!  ------ Solution of the system of linear equations

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, KRMAX)

!  ------ Assignment of the result

do kr=0, KRMAX
   temp_r_new(kr,j,i) = lgs_x(kr)
end do

!-------- Computation of the ice enthalpy (predictor step) --------

!  ------ Set-up of the equations

call calc_temp_enth_2_c(ct1, ct2, ct3, ct4, ce5, ce6, ce7, &
                        ct1_sg, ct2_sg, ct3_sg, ct4_sg, cm3, &
                        adv_vert_sg, abs_adv_vert_sg, &
                        dtime_temp, dtt_dxi, dtt_deta, &
                        dtt_2dxi, dtt_2deta, &
                        i, j, 0, &
                        lgs_a0, lgs_a1, lgs_a2, lgs_b)

!  ------ Solution of the system of linear equations

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, KCMAX)

!  ------ Assignment of the result

do kc=0, KCMAX
   enth_c_new(kc,j,i)  = lgs_x(kc)
   temp_c_new(kc,j,i)  = temp_fct_enth(enth_c_new(kc,j,i), temp_c_m(kc,j,i))
   omega_c_new(kc,j,i) = omega_fct_enth(enth_c_new(kc,j,i), temp_c_m(kc,j,i))
end do

!  ------ Determination of the CTS

kc_cts_new(j,i) = 0

#if !defined(ALLOW_TAPENADE) /* Normal */

do kc=1, KCMAX-1
   if (omega_c_new(kc,j,i) > eps_omega) then
      kc_cts_new(j,i) = kc
   else
      exit
   end if
end do

#else /* Tapenade */

kcdone = .false.

do kc=1, KCMAX-1
   if (kcdone.eqv..false.) then
      if (omega_c_new(kc,j,i) > eps_omega) then
         kc_cts_new(j,i) = kc
      else
         kcdone = .true.
      end if
   end if 
end do

#endif /* Normal vs. Tapenade */

!-------- Computation of the ice enthalpy
!         (corrector step for the cold-ice domain only
!         in order to fulfull the transition condition at the CTS) --------

#if (CALCMOD==3)   /* ENTM scheme */

if (kc_cts_new(j,i) > 0) then

!  ------ Update of the abbreviations where needed

   do kc=0, KCMAX
      temp_c_val(kc)  = temp_c_new(kc,j,i)
      omega_c_val(kc) = omega_c_new(kc,j,i)
   end do

   call calc_temp_enth_2_a2(at6, at7, ai2, am3, temp_c_val, omega_c_val, &
                            dtime_temp_inv, &
                            i, j, ce6, ce7, ci2, cm3)

!  ------ Set-up of the equations

   call calc_temp_enth_2_c(ct1, ct2, ct3, ct4, ce5, ce6, ce7, &
                           ct1_sg, ct2_sg, ct3_sg, ct4_sg, cm3, &
                           adv_vert_sg, abs_adv_vert_sg, &
                           dtime_temp, dtt_dxi, dtt_deta, &
                           dtt_2dxi, dtt_2deta, &
                           i, j, kc_cts_new(j,i), &
                           lgs_a0, lgs_a1, lgs_a2, lgs_b)

!  ------ Solution of the system of linear equations

   call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, KCMAX)

!  ------ Assignment of the result

   kc=kc_cts_new(j,i)
      enth_c_new(kc,j,i)  = lgs_x(kc)
      temp_c_new(kc,j,i)  = temp_fct_enth(enth_c_new(kc,j,i), temp_c_m(kc,j,i))
      omega_c_new(kc,j,i) = omega_fct_enth(enth_c_new(kc,j,i), temp_c_m(kc,j,i))

   do kc=kc_cts_new(j,i)+1, KCMAX
      enth_c_new(kc,j,i)  = lgs_x(kc)
      temp_c_new(kc,j,i)  = temp_fct_enth(enth_c_new(kc,j,i), temp_c_m(kc,j,i))
      omega_c_new(kc,j,i) = 0.0_dp   ! cold-ice domain
   end do

end if

#elif (CALCMOD==2)   /* ENTC scheme */

!!! continue   ! no corrector step

#else

errormsg = ' >>> calc_temp_enth_2: CALCMOD must be either 2 or 3!'
call error(errormsg)

#endif

!-------- Water drainage from temperate ice (if existing) --------

Q_tld(j,i) = 0.0_dp

do kc=0, kc_cts_new(j,i)

   if (omega_c_new(kc,j,i) > OMEGA_MAX) then

      Q_tld(j,i) = Q_tld(j,i) + cqtlde(kc)*(omega_c_new(kc,j,i)-OMEGA_MAX)

      omega_c_new(kc,j,i) = OMEGA_MAX
      enth_c_new(kc,j,i)  = enth_fct_temp_omega(temp_c_new(kc,j,i), OMEGA_MAX)

   end if

end do

!-------- Set enthalpy and water content in the redundant,
!         lower (kt) ice layer to the value at the ice base --------

do kt=0, KTMAX
   enth_t_new(kt,j,i)  = enth_c_new(0,j,i)
   omega_t_new(kt,j,i) = omega_c_new(0,j,i)
end do

!-------- Computation of the age of ice --------

!  ------ Set-up of the the equations

call calc_temp_enth_2_d(ct1, ct2, ct3, ct4, ci1, ci2, &
                        ct1_sg, ct2_sg, ct3_sg, ct4_sg, &
                        adv_vert_sg, abs_adv_vert_sg, &
                        dtime_temp, dtt_dxi, dtt_deta, &
                        dtt_2dxi, dtt_2deta, &
                        i, j, &
                        lgs_a0, lgs_a1, lgs_a2, lgs_b)

!  ------ Solution of the system of linear equations

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, KCMAX)

!  ------ Assignment of the result
!         restriction to interval [0, AGE_MAX yr]

do kc=0, KCMAX

   age_c_new(kc,j,i) = lgs_x(kc)

   if (age_c_new(kc,j,i) < (AGE_MIN*year2sec)) &
                           age_c_new(kc,j,i) = 0.0_dp
   if (age_c_new(kc,j,i) > (AGE_MAX*year2sec)) &
                           age_c_new(kc,j,i) = AGE_MAX*year2sec

end do

!-------- Age of the ice in the redundant, lower (kt) ice layer --------

do kt=0, KTMAX
   age_t_new(kt,j,i) = age_c_new(0,j,i)
end do

end subroutine calc_temp_enth_2

!-------------------------------------------------------------------------------
!> Computation of temperature and age for an ice column with a temperate base
!! with the enthalpy method: Abbreviations I.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_2_a1(at1, at2_1, at2_2, at3_1, at3_2, &
                               at4_1, at4_2, at5, atr1, alb1, &
                               ai1, aqtlde, &
                               dtime_temp, dtt_2dxi, dtt_2deta, i, j, &
                               ct1, ct2, ct3, ct4, ce5, ctr1, clb1, &
                               ct1_sg, ct2_sg, ct3_sg, ct4_sg, &
                               adv_vert_sg, abs_adv_vert_sg, &
                               ci1, cqtlde, dtt_dxi, dtt_deta)

implicit none

integer(i4b), intent(in) :: i, j

real(dp),     intent(in) :: at1(0:KCMAX), &
                            at2_1(0:KCMAX), at2_2(0:KCMAX), &
                            at3_1(0:KCMAX), at3_2(0:KCMAX), &
                            at4_1(0:KCMAX), at4_2(0:KCMAX), &
                            at5(0:KCMAX), ai1(0:KCMAX), &
                            atr1, alb1, aqtlde(0:KCMAX)
real(dp),     intent(in) :: dtime_temp, dtt_2dxi, dtt_2deta

real(dp),    intent(out) :: ct1(0:KCMAX), ct2(0:KCMAX), ct3(0:KCMAX), &
                            ct4(0:KCMAX), ce5(0:KCMAX), &
                            ctr1, clb1
real(dp),    intent(out) :: ct1_sg(0:KCMAX), ct2_sg(0:KCMAX), &
                            ct3_sg(0:KCMAX), ct4_sg(0:KCMAX), &
                            adv_vert_sg(0:KCMAX), abs_adv_vert_sg(0:KCMAX)
real(dp),    intent(out) :: ci1(0:KCMAX), cqtlde(0:KCMAX)
real(dp),    intent(out) :: dtt_dxi, dtt_deta

integer(i4b) :: kc

!-------- Initialisation --------

ct1             = 0.0_dp
ct2             = 0.0_dp
ct3             = 0.0_dp
ct4             = 0.0_dp
ce5             = 0.0_dp
ctr1            = 0.0_dp
clb1            = 0.0_dp
ct1_sg          = 0.0_dp
ct2_sg          = 0.0_dp
ct3_sg          = 0.0_dp
ct4_sg          = 0.0_dp
adv_vert_sg     = 0.0_dp
abs_adv_vert_sg = 0.0_dp
ci1             = 0.0_dp
cqtlde          = 0.0_dp
dtt_dxi         = 0.0_dp
dtt_deta        = 0.0_dp

!-------- Actual computation --------

ctr1 = atr1
clb1 = alb1*q_geo(j,i)

#if (ADV_VERT==1)

do kc=1, KCMAX-1
   ct1(kc) = at1(kc)/H_c(j,i)*0.5_dp*(vz_c(kc,j,i)+vz_c(kc-1,j,i))
end do

kc=0
ct1_sg(kc) = 0.5_dp*(at1(kc)+at1(kc+1))/H_c(j,i)*vz_c(kc,j,i)
             ! only needed for kc=0 ...
kc=KCMAX-1
ct1_sg(kc) = 0.5_dp*(at1(kc)+at1(kc+1))/H_c(j,i)*vz_c(kc,j,i)
             ! ... and kc=KCMAX-1

#elif (ADV_VERT==2 || ADV_VERT==3)

do kc=0, KCMAX-1
   ct1_sg(kc) = 0.5_dp*(at1(kc)+at1(kc+1))/H_c(j,i)*vz_c(kc,j,i)
end do

#endif

do kc=0, KCMAX

   ct2(kc) = ( at2_1(kc)*dzm_dtau(j,i) &
           +at2_2(kc)*dH_c_dtau(j,i) )/H_c(j,i)
   ct3(kc) = ( at3_1(kc)*dzm_dxi_g(j,i) &
           +at3_2(kc)*dH_c_dxi_g(j,i) )/H_c(j,i) &
          *0.5_dp*(vx_c(kc,j,i)+vx_c(kc,j,i-1)) *insq_g11_g(j,i)
   ct4(kc) = ( at4_1(kc)*dzm_deta_g(j,i) &
            +at4_2(kc)*dH_c_deta_g(j,i) )/H_c(j,i) &
          *0.5_dp*(vy_c(kc,j,i)+vy_c(kc,j-1,i)) *insq_g22_g(j,i)
   ce5(kc) = at5(kc)/H_c(j,i)
   ci1(kc) = ai1(kc)/H_c(j,i)
   cqtlde(kc) = aqtlde(kc)*H_c(j,i)

end do

#if (ADV_VERT==1)

kc=0
ct2_sg(kc) = 0.5_dp*(ct2(kc)+ct2(kc+1))
ct3_sg(kc) = 0.5_dp*(ct3(kc)+ct3(kc+1))
ct4_sg(kc) = 0.5_dp*(ct4(kc)+ct4(kc+1))
adv_vert_sg(kc) = ct1_sg(kc)-ct2_sg(kc)-ct3_sg(kc)-ct4_sg(kc)
abs_adv_vert_sg(kc) = abs(adv_vert_sg(kc))   ! only needed for kc=0 ...
kc=KCMAX-1
ct2_sg(kc) = 0.5_dp*(ct2(kc)+ct2(kc+1))
ct3_sg(kc) = 0.5_dp*(ct3(kc)+ct3(kc+1))
ct4_sg(kc) = 0.5_dp*(ct4(kc)+ct4(kc+1))
adv_vert_sg(kc) = ct1_sg(kc)-ct2_sg(kc)-ct3_sg(kc)-ct4_sg(kc)
abs_adv_vert_sg(kc) = abs(adv_vert_sg(kc))   ! ... and kc=KCMAX-1

#elif (ADV_VERT==2 || ADV_VERT==3)

do kc=0, KCMAX-1
   ct2_sg(kc) = 0.5_dp*(ct2(kc)+ct2(kc+1))
   ct3_sg(kc) = 0.5_dp*(ct3(kc)+ct3(kc+1))
   ct4_sg(kc) = 0.5_dp*(ct4(kc)+ct4(kc+1))
   adv_vert_sg(kc) = ct1_sg(kc)-ct2_sg(kc)-ct3_sg(kc)-ct4_sg(kc)
   abs_adv_vert_sg(kc) = abs(adv_vert_sg(kc))
end do

#endif

#if (ADV_HOR==3)
dtt_dxi  = 2.0_dp*dtt_2dxi
dtt_deta = 2.0_dp*dtt_2deta
#endif

end subroutine calc_temp_enth_2_a1

!-------------------------------------------------------------------------------
!> Computation of temperature and age for an ice column with a temperate base
!! with the enthalpy method: Abbreviations II.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_2_a2(at6, at7, ai2, am3, temp_c_val, omega_c_val, &
                               dtime_temp_inv, &
                               i, j, ce6, ce7, ci2, cm3)

use ice_material_properties_m, only : ratefac_c_t, kappa_val, c_val, &
                                      creep, viscosity

implicit none

integer(i4b), intent(in) :: i, j

real(dp),     intent(in) :: at6(0:KCMAX), at7, ai2(0:KCMAX), am3(0:KCMAX)
real(dp),     intent(in) :: temp_c_val(0:KCMAX), omega_c_val(0:KCMAX)
real(dp),     intent(in) :: dtime_temp_inv

real(dp),    intent(out) :: ce6(0:KCMAX), ce7(0:KCMAX), ci2(0:KCMAX), &
                            cm3(0:KCMAX)

integer(i4b) :: kc
real(dp) :: temp_c_help(0:KCMAX)

!-------- Initialisation --------

ce6 = 0.0_dp
ce7 = 0.0_dp
ci2 = 0.0_dp
cm3 = 0.0_dp

!-------- Actual computation --------

do kc=0, KCMAX

#if (DYNAMICS==2)
   if (.not.flag_shelfy_stream(j,i)) then
#endif
      ce7(kc) = at7 &
                *enh_c(kc,j,i) &
                *ratefac_c_t(temp_c_val(kc), omega_c_val(kc), temp_c_m(kc,j,i)) &
                *creep(sigma_c(kc,j,i)) &
                *sigma_c(kc,j,i)*sigma_c(kc,j,i)
#if (DYNAMICS==2)
   else
      ce7(kc) = 2.0_dp*at7 &
                *viscosity(de_c(kc,j,i), &
                           temp_c_val(kc), temp_c_m(kc,j,i), omega_c_val(kc), &
                           enh_c(kc,j,i), 2) &
                *de_c(kc,j,i)**2
   end if
#endif

   strain_heating_c(kc,j,i) = ce7(kc)*dtime_temp_inv

   cm3(kc) = am3(kc)*H_c(j,i)*c_val(temp_c_val(kc))

end do

do kc=0, kc_cts_new(j,i)-1   ! temperate layer
   ce6(kc) = at6(kc) &
             *NUE/H_c(j,i)
   ci2(kc) = ai2(kc)/H_c(j,i)
end do

do kc=kc_cts_new(j,i), KCMAX-1   ! cold layer
   temp_c_help(kc) = 0.5_dp*(temp_c_val(kc)+temp_c_val(kc+1))
   ce6(kc) = at6(kc) &
             *kappa_val(temp_c_help(kc))/(c_val(temp_c_help(kc))*H_c(j,i))
   ci2(kc) = ai2(kc)/H_c(j,i)
end do

end subroutine calc_temp_enth_2_a2

!-------------------------------------------------------------------------------
!> Computation of temperature and age for an ice column with a temperate base
!! with the enthalpy method:
!! Set-up of the equations for the bedrock temperature.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_2_b(ctr1, clb1, i, j, &
                              lgs_a0, lgs_a1, lgs_a2, lgs_b)

implicit none

integer(i4b), intent(in) :: i, j

real(dp),     intent(in) :: ctr1, clb1

real(dp),    intent(out) :: lgs_a0(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
                            lgs_a1(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
                            lgs_a2(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
                            lgs_b(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX)

integer(i4b) :: kr

!-------- Initialisation --------

lgs_a0 = 0.0_dp
lgs_a1 = 0.0_dp
lgs_a2 = 0.0_dp
lgs_b  = 0.0_dp

!-------- Actual computation --------

kr=0
lgs_a1(kr) = 1.0_dp
lgs_a2(kr) = -1.0_dp
lgs_b(kr)    = clb1

#if (Q_LITHO==1)
!   (coupled heat-conducting bedrock)

do kr=1, KRMAX-1
   lgs_a0(kr) = - ctr1
   lgs_a1(kr) = 1.0_dp + 2.0_dp*ctr1
   lgs_a2(kr) = - ctr1
   lgs_b(kr)    = temp_r(kr,j,i)
end do

#elif (Q_LITHO==0)
!   (no coupled heat-conducting bedrock)

do kr=1, KRMAX-1
   lgs_a0(kr) = 1.0_dp
   lgs_a1(kr) = 0.0_dp
   lgs_a2(kr) = -1.0_dp
   lgs_b(kr)  = 2.0_dp*clb1
end do

#endif

kr=KRMAX
lgs_a0(kr) = 0.0_dp
lgs_a1(kr) = 1.0_dp
lgs_b(kr)  = temp_t_m(0,j,i)

end subroutine calc_temp_enth_2_b

!-------------------------------------------------------------------------------
!> Computation of temperature and age for an ice column with a temperate base
!! with the enthalpy method:
!! Set-up of the equations for the ice enthalpy.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_2_c(ct1, ct2, ct3, ct4, ce5, ce6, ce7, &
                              ct1_sg, ct2_sg, ct3_sg, ct4_sg, cm3, &
                              adv_vert_sg, abs_adv_vert_sg, &
                              dtime_temp, dtt_dxi, dtt_deta, &
                              dtt_2dxi, dtt_2deta, &
                              i, j, kcmin, &
                              lgs_a0, lgs_a1, lgs_a2, lgs_b)

use enth_temp_omega_m, only : enth_fct_temp_omega

implicit none

integer(i4b), intent(in) :: i, j

integer(i4b), intent(in) :: kcmin
real(dp),     intent(in) :: ct1(0:KCMAX), ct2(0:KCMAX), ct3(0:KCMAX), &
                            ct4(0:KCMAX), ce5(0:KCMAX), ce6(0:KCMAX), &
                            ce7(0:KCMAX)
real(dp),     intent(in) :: ct1_sg(0:KCMAX), ct2_sg(0:KCMAX), &
                            ct3_sg(0:KCMAX), ct4_sg(0:KCMAX), &
                            adv_vert_sg(0:KCMAX), abs_adv_vert_sg(0:KCMAX)
real(dp),     intent(in) :: cm3(0:KCMAX)
real(dp),     intent(in) :: dtime_temp, dtt_dxi, dtt_deta, dtt_2dxi, dtt_2deta

real(dp),    intent(out) :: lgs_a0(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
                            lgs_a1(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
                            lgs_a2(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
                            lgs_b(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX)

integer(i4b) :: kc
real(dp) :: vx_c_help, vy_c_help
real(dp) :: adv_vert_help

!-------- Initialisation --------

lgs_a0 = 0.0_dp
lgs_a1 = 0.0_dp
lgs_a2 = 0.0_dp
lgs_b  = 0.0_dp

!-------- Actual computation --------

if (kcmin == 0) then   ! predictor step

   kc=0

   if (kc_cts_new(j,i) == 0) then   ! temperate base without temperate layer

      lgs_a1(kc) = 1.0_dp
      lgs_a2(kc) = 0.0_dp
      lgs_b(kc)  = enth_fct_temp_omega(temp_c_m(kc,j,i), 0.0_dp)

   else   ! kc_cts_new(j,i) > 0, temperate base with temperate layer

      lgs_a1(kc) =  1.0_dp
      lgs_a2(kc) = -1.0_dp
      lgs_b(kc)  =  0.0_dp

   end if

else   ! kcmin > 0, corrector step

   kc=0

   lgs_a1(kc) = 1.0_dp   ! dummy
   lgs_a2(kc) = 0.0_dp   ! setting,
   lgs_b(kc)  = 0.0_dp   ! not needed

   do kc=1, kcmin-1

      lgs_a0(kc) = 0.0_dp   ! dummy
      lgs_a1(kc) = 1.0_dp   ! setting,
      lgs_a2(kc) = 0.0_dp   ! not
      lgs_b(kc)  = 0.0_dp   ! needed

   end do

   kc=kcmin

   lgs_a0(kc) = 0.0_dp
   lgs_a1(kc) = 1.0_dp
   lgs_a2(kc) = -1.0_dp
   lgs_b(kc)  = -cm3(kc)

end if

do kc=kcmin+1, KCMAX-1

#if (ADV_VERT==1)

   lgs_a0(kc) = -0.5_dp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ce5(kc)*ce6(kc-1)
   lgs_a1(kc) = 1.0_dp+ce5(kc)*(ce6(kc)+ce6(kc-1))
   lgs_a2(kc) = 0.5_dp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ce5(kc)*ce6(kc)

#elif (ADV_VERT==2)

   lgs_a0(kc) &
         = -0.5_dp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -ce5(kc)*ce6(kc-1)
   lgs_a1(kc) &
         = 1.0_dp &
           +0.5_dp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -0.5_dp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           +ce5(kc)*(ce6(kc)+ce6(kc-1))
   lgs_a2(kc) &
         =  0.5_dp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           -ce5(kc)*ce6(kc)

#elif (ADV_VERT==3)

   adv_vert_help = 0.5_dp*(adv_vert_sg(kc)+adv_vert_sg(kc-1))

   lgs_a0(kc) &
         = -max(adv_vert_help, 0.0_dp) &
           -ce5(kc)*ce6(kc-1)
   lgs_a1(kc) &
         = 1.0_dp &
           +max(adv_vert_help, 0.0_dp)-min(adv_vert_help, 0.0_dp) &
           +ce5(kc)*(ce6(kc)+ce6(kc-1))
   lgs_a2(kc) &
         =  min(adv_vert_help, 0.0_dp) &
           -ce5(kc)*ce6(kc)

#endif

#if (ADV_HOR==2)

   lgs_b(kc) = enth_c(kc,j,i) + ce7(kc) &
       - ( dtt_2dxi* &
             ( (vx_c(kc,j,i)-abs(vx_c(kc,j,i))) &
               *(enth_c(kc,j,i+1)-enth_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +(vx_c(kc,j,i-1)+abs(vx_c(kc,j,i-1))) &
               *(enth_c(kc,j,i)-enth_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_2deta* &
             ( (vy_c(kc,j,i)-abs(vy_c(kc,j,i))) &
               *(enth_c(kc,j+1,i)-enth_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +(vy_c(kc,j-1,i)+abs(vy_c(kc,j-1,i))) &
               *(enth_c(kc,j,i)-enth_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#elif (ADV_HOR==3)

   vx_c_help = 0.5_dp*(vx_c(kc,j,i)+vx_c(kc,j,i-1))
   vy_c_help = 0.5_dp*(vy_c(kc,j,i)+vy_c(kc,j-1,i))

   lgs_b(kc) = enth_c(kc,j,i) + ce7(kc) &
       - ( dtt_dxi* &
             ( min(vx_c_help, 0.0_dp) &
               *(enth_c(kc,j,i+1)-enth_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +max(vx_c_help, 0.0_dp) &
               *(enth_c(kc,j,i)-enth_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_deta* &
             ( min(vy_c_help, 0.0_dp) &
               *(enth_c(kc,j+1,i)-enth_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +max(vy_c_help, 0.0_dp) &
               *(enth_c(kc,j,i)-enth_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#endif

end do

kc=KCMAX
lgs_a0(kc) = 0.0_dp
lgs_a1(kc) = 1.0_dp
lgs_b(kc)  = enth_fct_temp_omega(temp_s(j,i), 0.0_dp)
                                ! zero water content at the ice surface

end subroutine calc_temp_enth_2_c

!-------------------------------------------------------------------------------
!> Computation of temperature and age for an ice column with a temperate base
!! with the enthalpy method:
!! Set-up of the equations for the age of ice.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_2_d(ct1, ct2, ct3, ct4, ci1, ci2, &
                              ct1_sg, ct2_sg, ct3_sg, ct4_sg, &
                              adv_vert_sg, abs_adv_vert_sg, &
                              dtime_temp, dtt_dxi, dtt_deta, &
                              dtt_2dxi, dtt_2deta, &
                              i, j, &
                              lgs_a0, lgs_a1, lgs_a2, lgs_b)

implicit none

integer(i4b), intent(in) :: i, j

real(dp),     intent(in) :: ct1(0:KCMAX), ct2(0:KCMAX), ct3(0:KCMAX), &
                            ct4(0:KCMAX), ci1(0:KCMAX), ci2(0:KCMAX)
real(dp),     intent(in) :: ct1_sg(0:KCMAX), ct2_sg(0:KCMAX), &
                            ct3_sg(0:KCMAX), ct4_sg(0:KCMAX), &
                            adv_vert_sg(0:KCMAX), abs_adv_vert_sg(0:KCMAX)
real(dp),     intent(in) :: dtime_temp, dtt_dxi, dtt_deta, dtt_2dxi, dtt_2deta

real(dp),    intent(out) :: lgs_a0(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
                            lgs_a1(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
                            lgs_a2(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
                            lgs_b(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX)

integer(i4b) :: kc
real(dp) :: vx_c_help, vy_c_help
real(dp) :: adv_vert_help

!-------- Initialisation --------

lgs_a0 = 0.0_dp
lgs_a1 = 0.0_dp
lgs_a2 = 0.0_dp
lgs_b  = 0.0_dp

!-------- Actual computation --------

kc=0                                                 ! adv_vert_sg(0) <= 0
lgs_a1(kc) = 1.0_dp - min(adv_vert_sg(kc), 0.0_dp)   ! (directed downward)
lgs_a2(kc) = min(adv_vert_sg(kc), 0.0_dp)            ! assumed/enforced

#if (ADV_HOR==2)

lgs_b(kc) = age_c(kc,j,i) + dtime_temp &
       - ( dtt_2dxi* &
             ( (vx_c(kc,j,i)-abs(vx_c(kc,j,i))) &
               *(age_c(kc,j,i+1)-age_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +(vx_c(kc,j,i-1)+abs(vx_c(kc,j,i-1))) &
               *(age_c(kc,j,i)-age_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_2deta* &
             ( (vy_c(kc,j,i)-abs(vy_c(kc,j,i))) &
               *(age_c(kc,j+1,i)-age_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +(vy_c(kc,j-1,i)+abs(vy_c(kc,j-1,i))) &
               *(age_c(kc,j,i)-age_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#elif (ADV_HOR==3)

vx_c_help = 0.5_dp*(vx_c(kc,j,i)+vx_c(kc,j,i-1))
vy_c_help = 0.5_dp*(vy_c(kc,j,i)+vy_c(kc,j-1,i))

lgs_b(kc) = age_c(kc,j,i) + dtime_temp &
       - ( dtt_dxi* &
             ( min(vx_c_help, 0.0_dp) &
               *(age_c(kc,j,i+1)-age_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +max(vx_c_help, 0.0_dp) &
               *(age_c(kc,j,i)-age_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_deta* &
             ( min(vy_c_help, 0.0_dp) &
               *(age_c(kc,j+1,i)-age_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +max(vy_c_help, 0.0_dp) &
               *(age_c(kc,j,i)-age_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#endif

do kc=1, KCMAX-1

#if (ADV_VERT==1)

   lgs_a0(kc) = -0.5_dp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ci1(kc)*ci2(kc-1)
   lgs_a1(kc) = 1.0_dp+ci1(kc)*(ci2(kc)+ci2(kc-1))
   lgs_a2(kc) = 0.5_dp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ci1(kc)*ci2(kc)

#elif (ADV_VERT==2)

   lgs_a0(kc) = -0.5_dp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1))
   lgs_a1(kc) =  1.0_dp &
                +0.5_dp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
                -0.5_dp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  )
   lgs_a2(kc) =  0.5_dp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  )

#elif (ADV_VERT==3)

   adv_vert_help = 0.5_dp*(adv_vert_sg(kc)+adv_vert_sg(kc-1))

   lgs_a0(kc) = -max(adv_vert_help, 0.0_dp)
   lgs_a1(kc) =  1.0_dp &
                +max(adv_vert_help, 0.0_dp)-min(adv_vert_help, 0.0_dp)
   lgs_a2(kc) =  min(adv_vert_help, 0.0_dp)

#endif

#if (ADV_HOR==2)

   lgs_b(kc) = age_c(kc,j,i) + dtime_temp &
       - ( dtt_2dxi* &
             ( (vx_c(kc,j,i)-abs(vx_c(kc,j,i))) &
               *(age_c(kc,j,i+1)-age_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +(vx_c(kc,j,i-1)+abs(vx_c(kc,j,i-1))) &
               *(age_c(kc,j,i)-age_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_2deta* &
             ( (vy_c(kc,j,i)-abs(vy_c(kc,j,i))) &
               *(age_c(kc,j+1,i)-age_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +(vy_c(kc,j-1,i)+abs(vy_c(kc,j-1,i))) &
               *(age_c(kc,j,i)-age_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#elif (ADV_HOR==3)

   vx_c_help = 0.5_dp*(vx_c(kc,j,i)+vx_c(kc,j,i-1))
   vy_c_help = 0.5_dp*(vy_c(kc,j,i)+vy_c(kc,j-1,i))

   lgs_b(kc) = age_c(kc,j,i) + dtime_temp &
       - ( dtt_dxi* &
             ( min(vx_c_help, 0.0_dp) &
               *(age_c(kc,j,i+1)-age_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +max(vx_c_help, 0.0_dp) &
               *(age_c(kc,j,i)-age_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_deta* &
             ( min(vy_c_help, 0.0_dp) &
               *(age_c(kc,j+1,i)-age_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +max(vy_c_help, 0.0_dp) &
               *(age_c(kc,j,i)-age_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#endif

end do

kc=KCMAX
if (as_perp(j,i) >= 0.0_dp) then
   lgs_a0(kc) = 0.0_dp
   lgs_a1(kc) = 1.0_dp
   lgs_b(kc)  = 0.0_dp
else
   lgs_a0(kc) = -max(adv_vert_sg(kc-1), 0.0_dp)
   lgs_a1(kc) = 1.0_dp + max(adv_vert_sg(kc-1), 0.0_dp)
                       ! adv_vert_sg(KCMAX-1) >= 0 (directed upward)
                       ! assumed/enforced
#if (ADV_HOR==2)

   lgs_b(kc) = age_c(kc,j,i) + dtime_temp &
       - ( dtt_2dxi* &
             ( (vx_c(kc,j,i)-abs(vx_c(kc,j,i))) &
               *(age_c(kc,j,i+1)-age_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +(vx_c(kc,j,i-1)+abs(vx_c(kc,j,i-1))) &
               *(age_c(kc,j,i)-age_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_2deta* &
             ( (vy_c(kc,j,i)-abs(vy_c(kc,j,i))) &
               *(age_c(kc,j+1,i)-age_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +(vy_c(kc,j-1,i)+abs(vy_c(kc,j-1,i))) &
               *(age_c(kc,j,i)-age_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#elif (ADV_HOR==3)

   vx_c_help = 0.5_dp*(vx_c(kc,j,i)+vx_c(kc,j,i-1))
   vy_c_help = 0.5_dp*(vy_c(kc,j,i)+vy_c(kc,j-1,i))

   lgs_b(kc) = age_c(kc,j,i) + dtime_temp &
       - ( dtt_dxi* &
             ( min(vx_c_help, 0.0_dp) &
               *(age_c(kc,j,i+1)-age_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +max(vx_c_help, 0.0_dp) &
               *(age_c(kc,j,i)-age_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_deta* &
             ( min(vy_c_help, 0.0_dp) &
               *(age_c(kc,j+1,i)-age_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +max(vy_c_help, 0.0_dp) &
               *(age_c(kc,j,i)-age_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#endif

end if

end subroutine calc_temp_enth_2_d

!-------------------------------------------------------------------------------
!> Computation of temperature, age, water content and enthalpy for an
!! ice-free column.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_r(atr1, alb1, i, j)

#if !defined(ALLOW_TAPENADE) /* Normal */
use sico_maths_m, only : tri_sle
#else /* Tapenade */
use sico_maths_m
#endif /* Normal vs. Tapenade */

use enth_temp_omega_m, only : enth_fct_temp_omega

implicit none

integer(i4b), intent(in) :: i, j

real(dp), intent(in) :: atr1, alb1

integer(i4b) :: kc, kt, kr
real(dp) :: ctr1, clb1
real(dp) :: lgs_a0(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_a1(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_a2(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_x(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_b(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX)
real(dp) :: enth_val

!-------- Abbreviations --------

ctr1 = atr1
clb1 = alb1*q_geo(j,i)

!-------- Set up the equations for the bedrock temperature --------

kr=0
lgs_a1(kr) = 1.0_dp
lgs_a2(kr) = -1.0_dp
lgs_b(kr)    = clb1

#if (Q_LITHO==1)
!   (coupled heat-conducting bedrock)

do kr=1, KRMAX-1
   lgs_a0(kr) = - ctr1
   lgs_a1(kr) = 1.0_dp + 2.0_dp*ctr1
   lgs_a2(kr) = - ctr1
   lgs_b(kr)    = temp_r(kr,j,i)
end do

#elif (Q_LITHO==0)
!   (no coupled heat-conducting bedrock)

do kr=1, KRMAX-1
   lgs_a0(kr) = 1.0_dp
   lgs_a1(kr) = 0.0_dp
   lgs_a2(kr) = -1.0_dp
   lgs_b(kr)  = 2.0_dp*clb1
end do

#endif

kr=KRMAX
lgs_a0(kr) = 0.0_dp
lgs_a1(kr) = 1.0_dp
lgs_b(kr)   = temp_s(j,i)

!-------- Solve system of linear equations --------

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, KRMAX)

!-------- Assign the result --------

do kr=0, KRMAX
   temp_r_new(kr,j,i) = lgs_x(kr)
end do

!-------- Water content, age and enthalpy
!                        in the non-existing lower (kt) ice layer --------

enth_val = enth_fct_temp_omega(temp_s(j,i), 0.0_dp)

do kt=0, KTMAX
   omega_t_new(kt,j,i) = 0.0_dp
   age_t_new(kt,j,i)   = 0.0_dp
   enth_t_new(kt,j,i)  = enth_val
end do

!-------- Temperature, age, water content and enthalpy
!                      in the non-existing upper (kc) ice layer --------

do kc=0, KCMAX
   temp_c_new(kc,j,i)  = temp_s(j,i)
   age_c_new(kc,j,i)   = 0.0_dp
   omega_c_new(kc,j,i) = 0.0_dp
   enth_c_new(kc,j,i)  = enth_val
end do

end subroutine calc_temp_enth_r

!-------------------------------------------------------------------------------
!> Computation of temperature and age for ice shelves (floating ice)
!! with the enthalpy method.
!<------------------------------------------------------------------------------
subroutine calc_temp_enth_ssa(at1, at2_1, at2_2, at3_1, at3_2, &
                              at4_1, at4_2, at5, at6, at7, atr1, alb1, &
                              ai1, ai2, &
                              dtime_temp, dtt_2dxi, dtt_2deta, &
                              dtime_temp_inv, &
                              i, j)

use ice_material_properties_m, only : kappa_val, c_val, viscosity

#if !defined(ALLOW_TAPENADE) /* Normal */
use sico_maths_m, only : tri_sle
#else /* Tapenade */
use sico_maths_m
#endif /* Normal vs. Tapenade */

use enth_temp_omega_m, only : enth_fct_temp_omega, &
                              temp_fct_enth, omega_fct_enth

implicit none

integer(i4b), intent(in) :: i, j

real(dp), intent(in) :: at1(0:KCMAX), at2_1(0:KCMAX), at2_2(0:KCMAX), &
                        at3_1(0:KCMAX), at3_2(0:KCMAX), &
                        at4_1(0:KCMAX), at4_2(0:KCMAX), &
                        at5(0:KCMAX), at6(0:KCMAX), at7, &
                        ai1(0:KCMAX), ai2(0:KCMAX), &
                        atr1, alb1
real(dp), intent(in) :: dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv

integer(i4b) :: kc, kt, kr
real(dp) :: ct1(0:KCMAX), ct2(0:KCMAX), ct3(0:KCMAX), ct4(0:KCMAX), &
            ce5(0:KCMAX), ce6(0:KCMAX), ce7(0:KCMAX), ctr1, clb1
real(dp) :: ct1_sg(0:KCMAX), ct2_sg(0:KCMAX), ct3_sg(0:KCMAX), &
            ct4_sg(0:KCMAX), adv_vert_sg(0:KCMAX), abs_adv_vert_sg(0:KCMAX)
real(dp) :: ci1(0:KCMAX), ci2(0:KCMAX)
real(dp) :: temp_c_help(0:KCMAX)
real(dp) :: vx_c_help, vy_c_help
real(dp) :: adv_vert_help
real(dp) :: dtt_dxi, dtt_deta
real(dp) :: lgs_a0(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_a1(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_a2(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_x(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_b(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX)
real(dp), parameter :: zero=0.0_dp

!-------- Check for boundary points --------

if ((i == 0).or.(i == IMAX).or.(j == 0).or.(j == JMAX)) then
   errormsg = ' >>> calc_temp_enth_ssa: Boundary points not allowed.'
   call error(errormsg)
end if

!-------- Abbreviations --------

ctr1 = atr1
clb1 = alb1*q_geo(j,i)

#if (ADV_VERT==1)

do kc=1, KCMAX-1
   ct1(kc) = at1(kc)/H_c(j,i)*0.5_dp*(vz_c(kc,j,i)+vz_c(kc-1,j,i))
end do

kc=0
ct1_sg(kc) = 0.5_dp*(at1(kc)+at1(kc+1))/H_c(j,i)*vz_c(kc,j,i)
             ! only needed for kc=0 ...
kc=KCMAX-1
ct1_sg(kc) = 0.5_dp*(at1(kc)+at1(kc+1))/H_c(j,i)*vz_c(kc,j,i)
             ! ... and kc=KCMAX-1

#elif (ADV_VERT==2 || ADV_VERT==3)

do kc=0, KCMAX-1
   ct1_sg(kc) = 0.5_dp*(at1(kc)+at1(kc+1))/H_c(j,i)*vz_c(kc,j,i)
end do

#endif

do kc=0, KCMAX

   ct2(kc) = ( at2_1(kc)*dzm_dtau(j,i) &
           +at2_2(kc)*dH_c_dtau(j,i) )/H_c(j,i)
   ct3(kc) = ( at3_1(kc)*dzm_dxi_g(j,i) &
           +at3_2(kc)*dH_c_dxi_g(j,i) )/H_c(j,i) &
          *0.5_dp*(vx_c(kc,j,i)+vx_c(kc,j,i-1)) *insq_g11_g(j,i)
   ct4(kc) = ( at4_1(kc)*dzm_deta_g(j,i) &
            +at4_2(kc)*dH_c_deta_g(j,i) )/H_c(j,i) &
          *0.5_dp*(vy_c(kc,j,i)+vy_c(kc,j-1,i)) *insq_g22_g(j,i)
   ce5(kc) = at5(kc)/H_c(j,i)
   ce7(kc) = 2.0_dp*at7 &
             *viscosity(de_ssa(j,i), &
                        temp_c(kc,j,i), temp_c_m(kc,j,i), 0.0_dp, &
                        enh_c(kc,j,i), 0) &
             *de_ssa(j,i)**2

   strain_heating_c(kc,j,i) = ce7(kc)*dtime_temp_inv

   ci1(kc) = ai1(kc)/H_c(j,i)

end do

#if (ADV_VERT==1)

kc=0
ct2_sg(kc) = 0.5_dp*(ct2(kc)+ct2(kc+1))
ct3_sg(kc) = 0.5_dp*(ct3(kc)+ct3(kc+1))
ct4_sg(kc) = 0.5_dp*(ct4(kc)+ct4(kc+1))
adv_vert_sg(kc) = ct1_sg(kc)-ct2_sg(kc)-ct3_sg(kc)-ct4_sg(kc)
abs_adv_vert_sg(kc) = abs(adv_vert_sg(kc))   ! only needed for kc=0 ...
kc=KCMAX-1
ct2_sg(kc) = 0.5_dp*(ct2(kc)+ct2(kc+1))
ct3_sg(kc) = 0.5_dp*(ct3(kc)+ct3(kc+1))
ct4_sg(kc) = 0.5_dp*(ct4(kc)+ct4(kc+1))
adv_vert_sg(kc) = ct1_sg(kc)-ct2_sg(kc)-ct3_sg(kc)-ct4_sg(kc)
abs_adv_vert_sg(kc) = abs(adv_vert_sg(kc))   ! ... and kc=KCMAX-1

#elif (ADV_VERT==2 || ADV_VERT==3)

do kc=0, KCMAX-1
   ct2_sg(kc) = 0.5_dp*(ct2(kc)+ct2(kc+1))
   ct3_sg(kc) = 0.5_dp*(ct3(kc)+ct3(kc+1))
   ct4_sg(kc) = 0.5_dp*(ct4(kc)+ct4(kc+1))
   adv_vert_sg(kc) = ct1_sg(kc)-ct2_sg(kc)-ct3_sg(kc)-ct4_sg(kc)
   abs_adv_vert_sg(kc) = abs(adv_vert_sg(kc))
end do

#endif

do kc=0, KCMAX-1
   temp_c_help(kc) = 0.5_dp*(temp_c(kc,j,i)+temp_c(kc+1,j,i))
   ce6(kc) = at6(kc) &
             *kappa_val(temp_c_help(kc))/(c_val(temp_c_help(kc))*H_c(j,i))
   ci2(kc) = ai2(kc)/H_c(j,i)
end do

#if (ADV_HOR==3)
dtt_dxi  = 2.0_dp*dtt_2dxi
dtt_deta = 2.0_dp*dtt_2deta
#endif

!-------- Set up the equations for the bedrock temperature --------

kr=0
lgs_a1(kr) = 1.0_dp
lgs_a2(kr) = -1.0_dp
lgs_b(kr)    = clb1

#if (Q_LITHO==1)
!   (coupled heat-conducting bedrock)

do kr=1, KRMAX-1
   lgs_a0(kr) = - ctr1
   lgs_a1(kr) = 1.0_dp + 2.0_dp*ctr1
   lgs_a2(kr) = - ctr1
   lgs_b(kr)    = temp_r(kr,j,i)
end do

#elif (Q_LITHO==0)
!   (no coupled heat-conducting bedrock)

do kr=1, KRMAX-1
   lgs_a0(kr) = 1.0_dp
   lgs_a1(kr) = 0.0_dp
   lgs_a2(kr) = -1.0_dp
   lgs_b(kr)  = 2.0_dp*clb1
end do

#endif

kr=KRMAX
lgs_a0(kr) = 0.0_dp
lgs_a1(kr) = 1.0_dp
lgs_b(kr)  = temp_c_m(0,j,i)-DELTA_TM_SW

!-------- Solve system of linear equations --------

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, KRMAX)

!-------- Assign the result --------

do kr=0, KRMAX
   temp_r_new(kr,j,i) = lgs_x(kr)
end do

!-------- Set up the equations for the ice temperature --------

kc=0
lgs_a1(kc) = 1.0_dp
lgs_a2(kc) = 0.0_dp
lgs_b(kc)  = enth_fct_temp_omega(temp_c_m(kc,j,i)-DELTA_TM_SW, 0.0_dp)
                                                  ! zero water content assumed

do kc=1, KCMAX-1

#if (ADV_VERT==1)

   lgs_a0(kc) = -0.5_dp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ce5(kc)*ce6(kc-1)
   lgs_a1(kc) = 1.0_dp+ce5(kc)*(ce6(kc)+ce6(kc-1))
   lgs_a2(kc) = 0.5_dp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ce5(kc)*ce6(kc)

#elif (ADV_VERT==2)

   lgs_a0(kc) &
         = -0.5_dp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -ce5(kc)*ce6(kc-1)
   lgs_a1(kc) &
         = 1.0_dp &
           +0.5_dp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -0.5_dp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           +ce5(kc)*(ce6(kc)+ce6(kc-1))
   lgs_a2(kc) &
         =  0.5_dp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           -ce5(kc)*ce6(kc)

#elif (ADV_VERT==3)

   adv_vert_help = 0.5_dp*(adv_vert_sg(kc)+adv_vert_sg(kc-1))

   lgs_a0(kc) &
         = -max(adv_vert_help, 0.0_dp) &
           -ce5(kc)*ce6(kc-1)
   lgs_a1(kc) &
         = 1.0_dp &
           +max(adv_vert_help, 0.0_dp)-min(adv_vert_help, 0.0_dp) &
           +ce5(kc)*(ce6(kc)+ce6(kc-1))
   lgs_a2(kc) &
         =  min(adv_vert_help, 0.0_dp) &
           -ce5(kc)*ce6(kc)

#endif

#if (ADV_HOR==2)

   lgs_b(kc) = enth_c(kc,j,i) + ce7(kc) &
       - ( dtt_2dxi* &
             ( (vx_c(kc,j,i)-abs(vx_c(kc,j,i))) &
               *(enth_c(kc,j,i+1)-enth_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +(vx_c(kc,j,i-1)+abs(vx_c(kc,j,i-1))) &
               *(enth_c(kc,j,i)-enth_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_2deta* &
             ( (vy_c(kc,j,i)-abs(vy_c(kc,j,i))) &
               *(enth_c(kc,j+1,i)-enth_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +(vy_c(kc,j-1,i)+abs(vy_c(kc,j-1,i))) &
               *(enth_c(kc,j,i)-enth_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#elif (ADV_HOR==3)

   vx_c_help = 0.5_dp*(vx_c(kc,j,i)+vx_c(kc,j,i-1))
   vy_c_help = 0.5_dp*(vy_c(kc,j,i)+vy_c(kc,j-1,i))

   lgs_b(kc) = enth_c(kc,j,i) + ce7(kc) &
       - ( dtt_dxi* &
             ( min(vx_c_help, 0.0_dp) &
               *(enth_c(kc,j,i+1)-enth_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +max(vx_c_help, 0.0_dp) &
               *(enth_c(kc,j,i)-enth_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_deta* &
             ( min(vy_c_help, 0.0_dp) &
               *(enth_c(kc,j+1,i)-enth_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +max(vy_c_help, 0.0_dp) &
               *(enth_c(kc,j,i)-enth_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#endif

end do

kc=KCMAX
lgs_a0(kc) = 0.0_dp
lgs_a1(kc) = 1.0_dp
lgs_b(kc)  = enth_fct_temp_omega(temp_s(j,i), 0.0_dp)
             ! zero water content assumed

!-------- Solve system of linear equations --------

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, KCMAX)

!-------- Assign the result --------

do kc=0, KCMAX
   enth_c_new(kc,j,i)  = lgs_x(kc)
   temp_c_new(kc,j,i)  = temp_fct_enth(enth_c_new(kc,j,i), temp_c_m(kc,j,i))
   omega_c_new(kc,j,i) = omega_fct_enth(enth_c_new(kc,j,i), temp_c_m(kc,j,i))
end do

!-------- Set enthalpy and water content in the redundant,
!         lower (kt) ice layer to the value at the ice base --------

do kt=0, KTMAX
   enth_t_new(kt,j,i)  = enth_c_new(0,j,i)
   omega_t_new(kt,j,i) = omega_c_new(0,j,i)
end do

!-------- Water drainage from the non-existing temperate ice --------

Q_tld(j,i) = 0.0_dp

!-------- Set up the equations for the age of cold ice --------

kc=0                                                 ! adv_vert_sg(0) <= 0
lgs_a1(kc) = 1.0_dp - min(adv_vert_sg(kc), 0.0_dp)   ! (directed downward)
lgs_a2(kc) = min(adv_vert_sg(kc), 0.0_dp)            ! assumed/enforced

#if (ADV_HOR==2)

lgs_b(kc) = age_c(kc,j,i) + dtime_temp &
       - ( dtt_2dxi* &
             ( (vx_c(kc,j,i)-abs(vx_c(kc,j,i))) &
               *(age_c(kc,j,i+1)-age_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +(vx_c(kc,j,i-1)+abs(vx_c(kc,j,i-1))) &
               *(age_c(kc,j,i)-age_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_2deta* &
             ( (vy_c(kc,j,i)-abs(vy_c(kc,j,i))) &
               *(age_c(kc,j+1,i)-age_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +(vy_c(kc,j-1,i)+abs(vy_c(kc,j-1,i))) &
               *(age_c(kc,j,i)-age_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#elif (ADV_HOR==3)

vx_c_help = 0.5_dp*(vx_c(kc,j,i)+vx_c(kc,j,i-1))
vy_c_help = 0.5_dp*(vy_c(kc,j,i)+vy_c(kc,j-1,i))

lgs_b(kc) = age_c(kc,j,i) + dtime_temp &
       - ( dtt_dxi* &
             ( min(vx_c_help, 0.0_dp) &
               *(age_c(kc,j,i+1)-age_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +max(vx_c_help, 0.0_dp) &
               *(age_c(kc,j,i)-age_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_deta* &
             ( min(vy_c_help, 0.0_dp) &
               *(age_c(kc,j+1,i)-age_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +max(vy_c_help, 0.0_dp) &
               *(age_c(kc,j,i)-age_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#endif

do kc=1, KCMAX-1

#if (ADV_VERT==1)

   lgs_a0(kc) = -0.5_dp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ci1(kc)*ci2(kc-1)
   lgs_a1(kc) = 1.0_dp+ci1(kc)*(ci2(kc)+ci2(kc-1))
   lgs_a2(kc) = 0.5_dp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ci1(kc)*ci2(kc)

#elif (ADV_VERT==2)

   lgs_a0(kc) = -0.5_dp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1))
   lgs_a1(kc) =  1.0_dp &
                +0.5_dp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
                -0.5_dp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  )
   lgs_a2(kc) =  0.5_dp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  )

#elif (ADV_VERT==3)

   adv_vert_help = 0.5_dp*(adv_vert_sg(kc)+adv_vert_sg(kc-1))

   lgs_a0(kc) = -max(adv_vert_help, 0.0_dp)
   lgs_a1(kc) =  1.0_dp &
                +max(adv_vert_help, 0.0_dp)-min(adv_vert_help, 0.0_dp)
   lgs_a2(kc) =  min(adv_vert_help, 0.0_dp)

#endif

#if (ADV_HOR==2)

   lgs_b(kc) = age_c(kc,j,i) + dtime_temp &
       - ( dtt_2dxi* &
             ( (vx_c(kc,j,i)-abs(vx_c(kc,j,i))) &
               *(age_c(kc,j,i+1)-age_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +(vx_c(kc,j,i-1)+abs(vx_c(kc,j,i-1))) &
               *(age_c(kc,j,i)-age_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_2deta* &
             ( (vy_c(kc,j,i)-abs(vy_c(kc,j,i))) &
               *(age_c(kc,j+1,i)-age_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +(vy_c(kc,j-1,i)+abs(vy_c(kc,j-1,i))) &
               *(age_c(kc,j,i)-age_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#elif (ADV_HOR==3)

   vx_c_help = 0.5_dp*(vx_c(kc,j,i)+vx_c(kc,j,i-1))
   vy_c_help = 0.5_dp*(vy_c(kc,j,i)+vy_c(kc,j-1,i))

   lgs_b(kc) = age_c(kc,j,i) + dtime_temp &
       - ( dtt_dxi* &
             ( min(vx_c_help, 0.0_dp) &
               *(age_c(kc,j,i+1)-age_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +max(vx_c_help, 0.0_dp) &
               *(age_c(kc,j,i)-age_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_deta* &
             ( min(vy_c_help, 0.0_dp) &
               *(age_c(kc,j+1,i)-age_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +max(vy_c_help, 0.0_dp) &
               *(age_c(kc,j,i)-age_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#endif

end do

kc=KCMAX
if (as_perp(j,i) >= zero) then
   lgs_a0(kc) = 0.0_dp
   lgs_a1(kc) = 1.0_dp
   lgs_b(kc)  = 0.0_dp
else
   lgs_a0(kc) = -max(adv_vert_sg(kc-1), 0.0_dp)
   lgs_a1(kc) = 1.0_dp + max(adv_vert_sg(kc-1), 0.0_dp)
                       ! adv_vert_sg(KCMAX-1) >= 0 (directed upward)
                       ! assumed/enforced
#if (ADV_HOR==2)

   lgs_b(kc) = age_c(kc,j,i) + dtime_temp &
       - ( dtt_2dxi* &
             ( (vx_c(kc,j,i)-abs(vx_c(kc,j,i))) &
               *(age_c(kc,j,i+1)-age_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +(vx_c(kc,j,i-1)+abs(vx_c(kc,j,i-1))) &
               *(age_c(kc,j,i)-age_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_2deta* &
             ( (vy_c(kc,j,i)-abs(vy_c(kc,j,i))) &
               *(age_c(kc,j+1,i)-age_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +(vy_c(kc,j-1,i)+abs(vy_c(kc,j-1,i))) &
               *(age_c(kc,j,i)-age_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#elif (ADV_HOR==3)

   vx_c_help = 0.5_dp*(vx_c(kc,j,i)+vx_c(kc,j,i-1))
   vy_c_help = 0.5_dp*(vy_c(kc,j,i)+vy_c(kc,j-1,i))

   lgs_b(kc) = age_c(kc,j,i) + dtime_temp &
       - ( dtt_dxi* &
             ( min(vx_c_help, 0.0_dp) &
               *(age_c(kc,j,i+1)-age_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +max(vx_c_help, 0.0_dp) &
               *(age_c(kc,j,i)-age_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_deta* &
             ( min(vy_c_help, 0.0_dp) &
               *(age_c(kc,j+1,i)-age_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +max(vy_c_help, 0.0_dp) &
               *(age_c(kc,j,i)-age_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#endif

end if

!-------- Solve system of linear equations --------

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, KCMAX)

!-------- Assign the result,
!         restriction to interval [0, AGE_MAX yr] --------

do kc=0, KCMAX

   age_c_new(kc,j,i) = lgs_x(kc)

   if (age_c_new(kc,j,i) < (AGE_MIN*year2sec)) &
                           age_c_new(kc,j,i) = 0.0_dp
   if (age_c_new(kc,j,i) > (AGE_MAX*year2sec)) &
                           age_c_new(kc,j,i) = AGE_MAX*year2sec

end do

!-------- Age of the ice in the redundant, lower (kt) ice layer --------

do kt=0, KTMAX
   age_t_new(kt,j,i) = age_c_new(0,j,i)
end do

end subroutine calc_temp_enth_ssa

!-------------------------------------------------------------------------------

end module calc_temp_enth_m
!
