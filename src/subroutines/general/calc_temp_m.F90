!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ t e m p _ m
!
!! Computation of temperature, water content and age.
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
!> Computation of temperature, water content and age.
!-------------------------------------------------------------------------------
module calc_temp_m

  use sico_types_m
  use sico_variables_m

#if (defined(EISMINT) || defined(HEINO) || defined(MOCHO) || defined(NMARS) || defined(SMARS) || defined(XYZ))
  use sico_vars_m
#endif

  use error_m

  implicit none

#if !defined(ALLOW_TAPENADE) /* Normal */
  private
#endif /* Normal */

  public :: calc_temp_poly, calc_temp_cold, calc_temp_const

contains

!-------------------------------------------------------------------------------
!> Computation of temperature, water content and age in polythermal mode.
!-------------------------------------------------------------------------------
subroutine calc_temp_poly(dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
                          dtime_temp)

use ice_material_properties_m, only : kappa_val

implicit none

real(dp), intent(in) :: dxi, deta, dzeta_c, dzeta_t, dzeta_r
real(dp), intent(in) :: dtime_temp

integer(i4b) :: i, j, kc, kt, kr, ii, jj
real(dp) :: at1(0:KCMAX), at2_1(0:KCMAX), at2_2(0:KCMAX), &
            at3_1(0:KCMAX), at3_2(0:KCMAX), at4_1(0:KCMAX), &
            at4_2(0:KCMAX), at5(0:KCMAX), at6(0:KCMAX), at7, &
            acb1, acb2, acb3, acb4, &
            ai1(0:KCMAX), ai2(0:KCMAX), ai3, &
            atr1, am1, am2, alb1
real(dp) :: aw1, aw2, aw3, aw4, aw5, aw7, aw8, aw9, aqtld
real(dp) :: dtime_temp_inv, dtt_2dxi, dtt_2deta
real(dp) :: temp_c_help(0:KCMAX)
real(dp) :: fact_thick
real(dp) :: time_lag_cts
real(dp) :: Vol_t, Vol_t_smooth, korrfakt_t, &
            dH_t_smooth(0:JMAX,0:IMAX)

!-------- Term abbreviations

at7 = 2.0_dp/RHO*dtime_temp

aw1 = dtime_temp/dzeta_t
aw2 = dtime_temp/dzeta_t
aw3 = dtime_temp/dzeta_t
aw4 = dtime_temp/dzeta_t
aw5 = NUE/RHO*dtime_temp/(dzeta_t**2)
aw7 = 2.0_dp/(RHO*L)*dtime_temp
aw8 = BETA**2/(RHO*L) &
      *(kappa_val(0.0_dp)-kappa_val(-1.0_dp))*dtime_temp
aw9 = BETA/L*dtime_temp

ai3 = AGEDIFF*dtime_temp/(dzeta_t**2)

atr1 = KAPPA_R/(RHO_C_R*H_R**2)*dtime_temp/(dzeta_r**2)

if (flag_aa_nonzero) then
   am1 = aa*BETA*dzeta_c/(ea-1.0_dp)
   am2 = aa*L*RHO*dzeta_c/(ea-1.0_dp)
else
   am1 = BETA*dzeta_c
   am2 = L*RHO*dzeta_c
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

aqtld = dzeta_t/dtime_temp

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

   end if

end do

strain_heating_c = 0.0_dp   ! initialization,
strain_heating_t = 0.0_dp   ! purely diagnostic fields

!-------- Computation loop for temperature, water content and age --------

do i=1, IMAX-1   ! skipping domain margins
do j=1, JMAX-1   ! skipping domain margins

   if (mask(j,i)==0) then   ! glaciated land

!  ------ Old vertical column cold

      if (n_cts(j,i) == -1) then

         n_cts_new(j,i) = n_cts(j,i)
         zm_new(j,i)  = zb(j,i)
         H_c_new(j,i) = H_c(j,i)
         H_t_new(j,i) = H_t(j,i)
         !$AD NOCHECKPOINT
         call calc_temp1(at1, at2_1, at2_2, at3_1, at3_2, &
           at4_1, at4_2, at5, at6, at7, atr1, acb1, acb2, &
           acb3, acb4, alb1, ai1, ai2, &
           dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
           i, j)

!    ---- Check whether base has become temperate

         if (temp_c_new(0,j,i) > temp_c_m(0,j,i)) then

            n_cts_new(j,i) = 0
            !$AD NOCHECKPOINT
            call calc_temp2(at1, at2_1, at2_2, at3_1, at3_2, &
                 at4_1, at4_2, at5, at6, at7, atr1, alb1, &
                 ai1, ai2, &
                 dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
                 i, j)

         end if

!    ---- Check whether even temperate layer has formed

         if ( &
              ( n_cts_new(j,i) == 0 ).and. &
              ( (temp_c_new(1,j,i)-temp_c_new(0,j,i)) &
                 > (am1*H_c_new(j,i)) ) &
            ) then

            n_cts_new(j,i) = 1
            zm_new(j,i)  = zb(j,i)+0.001_dp
            H_c_new(j,i) = H_c(j,i)-0.001_dp
            H_t_new(j,i) = H_t(j,i)+0.001_dp
!                 ! CTS preliminarily positioned 1 mm above ice base --------
            !$AD NOCHECKPOINT
            call calc_temp3(at1, at2_1, at2_2, at3_1, at3_2, &
                 at4_1, at4_2, at5, at6, at7, atr1, &
                 am1, am2, alb1, &
                 aw1, aw2, aw3, aw4, aw5, aw7, aw8, aw9, aqtld, &
                 ai1, ai2, ai3, dzeta_t, &
                 dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
                 i, j)
            !$AD NOCHECKPOINT
            call shift_cts_upward(at1, at2_1, at2_2, at3_1, at3_2, &
              at4_1, at4_2, at5, at6, at7, atr1, am1, am2, alb1, &
              aw1, aw2, aw3, aw4, aw5, aw7, aw8, aw9, aqtld, &
              ai1, ai2, ai3, dzeta_t, &
              dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
              i, j)

         end if

!  ------ Old vertical column with temperate base

      else if (n_cts(j,i) == 0) then

         n_cts_new(j,i) = n_cts(j,i)
         zm_new(j,i)  = zb(j,i)
         H_c_new(j,i) = H_c(j,i)
         H_t_new(j,i) = H_t(j,i)
         !$AD NOCHECKPOINT
         call calc_temp2(at1, at2_1, at2_2, at3_1, at3_2, &
              at4_1, at4_2, at5, at6, at7, atr1, alb1, &
              ai1, ai2, &
              dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
              i, j)

!    ---- Check whether temperate base becomes cold

         if ( (temp_c_new(1,j,i)-temp_c_new(0,j,i)) <  (am1*H_c(j,i)) ) then

            n_cts_new(j,i) = -1
            !$AD NOCHECKPOINT
            call calc_temp1(at1, at2_1, at2_2, at3_1, at3_2, &
                 at4_1, at4_2, at5, at6, at7, atr1, acb1, acb2, &
                 acb3, acb4, alb1, ai1, ai2, &
                 dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
                 i, j)

            if (temp_c_new(0,j,i) >= temp_c_m(0,j,i)) then

               n_cts_new(j,i) = 0
               !$AD NOCHECKPOINT
               call calc_temp2(at1, at2_1, at2_2, at3_1, at3_2, &
                    at4_1, at4_2, at5, at6, at7, atr1, alb1, &
                    ai1, ai2, &
                    dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
                    i, j)

            end if

         end if

!    ---- Check whether temperate layer has formed

         if ( &
              ( n_cts_new(j,i) == 0 ).and. &
              ( (temp_c_new(1,j,i)-temp_c_new(0,j,i)) &
                 > (am1*H_c_new(j,i)) ) &
            ) then

            n_cts_new(j,i) = 1
            zm_new(j,i)  = zb(j,i)+0.001_dp
            H_c_new(j,i) = H_c(j,i)-0.001_dp
            H_t_new(j,i) = H_t(j,i)+0.001_dp
!                 ! CTS preliminarily positioned 1 mm above ice base --------
            !$AD NOCHECKPOINT
            call calc_temp3(at1, at2_1, at2_2, at3_1, at3_2, &
                 at4_1, at4_2, at5, at6, at7, atr1, &
                 am1, am2, alb1, &
                 aw1, aw2, aw3, aw4, aw5, aw7, aw8, aw9, aqtld, &
                 ai1, ai2, ai3, dzeta_t, &
                 dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
                 i, j)
            !$AD NOCHECKPOINT
            call shift_cts_upward(at1, at2_1, at2_2, at3_1, at3_2, &
              at4_1, at4_2, at5, at6, at7, atr1, am1, am2, alb1, &
              aw1, aw2, aw3, aw4, aw5, aw7, aw8, aw9, aqtld, &
              ai1, ai2, ai3, dzeta_t, &
              dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
              i, j)

         end if

!  ------ Old vertical column with temperate base and temperate layer

      else   ! n_cts(j,i) == 1

         n_cts_new(j,i) = n_cts(j,i)
         zm_new(j,i)  = zm(j,i)
         H_c_new(j,i) = H_c(j,i)
         H_t_new(j,i) = H_t(j,i)
         !$AD NOCHECKPOINT
         call calc_temp3(at1, at2_1, at2_2, at3_1, at3_2, &
                 at4_1, at4_2, at5, at6, at7, atr1, &
                 am1, am2, alb1, &
                 aw1, aw2, aw3, aw4, aw5, aw7, aw8, aw9, aqtld, &
                 ai1, ai2, ai3, dzeta_t, &
                 dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
                 i, j)

         if ( (temp_c_new(0,j,i)-(-BETA*H_c_new(j,i))) > 0.0_dp ) &
         then
            !$AD NOCHECKPOINT
            call shift_cts_upward(at1, at2_1, at2_2, at3_1, at3_2, &
              at4_1, at4_2, at5, at6, at7, atr1, am1, am2, alb1, &
              aw1, aw2, aw3, aw4, aw5, aw7, aw8, aw9, aqtld, &
              ai1, ai2, ai3, dzeta_t, &
              dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
              i, j)
         else
            !$AD NOCHECKPOINT
            call shift_cts_downward(at1, at2_1, at2_2, at3_1, at3_2, &
              at4_1, at4_2, at5, at6, at7, atr1, am1, am2, alb1, &
              aw1, aw2, aw3, aw4, aw5, aw7, aw8, aw9, aqtld, &
              ai1, ai2, ai3, dzeta_t, &
              dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
              i, j)
         end if

      end if

#if (MARGIN==3)

   else if (mask(j,i)==3) then   ! floating ice

      n_cts_new(j,i) = -1
      zm_new(j,i)  = zb(j,i)
      H_c_new(j,i) = H_c(j,i) + H_t(j,i)
      H_t_new(j,i) = 0.0_dp
      !$AD NOCHECKPOINT
      call calc_temp_ssa(at1, at2_1, at2_2, at3_1, at3_2, &
           at4_1, at4_2, at5, at6, at7, atr1, alb1, &
           ai1, ai2, &
           dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
           i, j)

!  ------ Reset temperatures above melting to the melting point
!         (should not occur, but just in case)

      do kc=0, KCMAX
         if (temp_c_new(kc,j,i) > temp_c_m(kc,j,i)) &
                    temp_c_new(kc,j,i) = temp_c_m(kc,j,i)
      end do

#endif

   else   ! mask(j,i) == 1, 2 (ice-free land or sea point)

      n_cts_new(j,i) = -1
      zm_new(j,i)  = zb(j,i)
      H_c_new(j,i) = H_c(j,i)
      H_t_new(j,i) = H_t(j,i)
      !$AD NOCHECKPOINT
      call calc_temp_r(atr1, alb1, i, j)

endif

end do
end do   ! End of computation loop

!-------- Extrapolate values on margins --------

!  ------ Lower left corner

i=0
j=0

if ( (mask(j,i) == 0).or.(mask(j,i) == 3) ) then
                              ! glaciated land or floating ice
   ii=i+1
   jj=j+1

   do kc=0,KCMAX
      temp_c_new(kc,j,i) = temp_c_new(kc,jj,ii)   ! set cold-ice temperature
      age_c_new(kc,j,i)  = age_c_new(kc,jj,ii)    ! set cold-ice age
   end do

   do kt=0,KTMAX
      omega_t_new(kt,j,i) = omega_t_new(kt,jj,ii) ! set temp.-ice water content
      age_t_new(kt,j,i)   = age_t_new(kt,jj,ii)   ! set temp.-ice age
   end do

   do kr=0,KRMAX
      temp_r_new(kr,j,i) = temp_r_new(kr,jj,ii)   ! set bedrock temperature
   end do

   n_cts_new(j,i) = min(n_cts_new(jj,ii),0)   ! temperate layer excluded
   H_c_new(j,i)   = H_c(j,i)
   H_t_new(j,i)   = H_t(j,i)
   zm_new(j,i)  = zb(j,i)

else   ! mask(j,i) == 1, 2 (ice-free land or sea point)

   n_cts_new(j,i) = -1
   zm_new(j,i)    = zb(j,i)
   H_c_new(j,i)   = H_c(j,i)
   H_t_new(j,i)   = H_t(j,i)

   call calc_temp_r(atr1, alb1, i, j)

end if

!  ------ Lower right corner

i=IMAX
j=0

if ( (mask(j,i) == 0).or.(mask(j,i) == 3) ) then
                              ! glaciated land or floating ice
   ii=i-1
   jj=j+1

   do kc=0,KCMAX
      temp_c_new(kc,j,i) = temp_c_new(kc,jj,ii)   ! set cold-ice temperature
      age_c_new(kc,j,i)  = age_c_new(kc,jj,ii)    ! set cold-ice age
   end do

   do kt=0,KTMAX
      omega_t_new(kt,j,i) = omega_t_new(kt,jj,ii) ! set temp.-ice water content
      age_t_new(kt,j,i)   = age_t_new(kt,jj,ii)   ! set temp.-ice age
   end do

   do kr=0,KRMAX
      temp_r_new(kr,j,i) = temp_r_new(kr,jj,ii)   ! set bedrock temperature
   end do

   n_cts_new(j,i) = min(n_cts_new(jj,ii),0)   ! temperate layer excluded
   H_c_new(j,i)   = H_c(j,i)
   H_t_new(j,i)   = H_t(j,i)
   zm_new(j,i)  = zb(j,i)

else   ! mask(j,i) == 1, 2 (ice-free land or sea point)

   n_cts_new(j,i) = -1
   zm_new(j,i)    = zb(j,i)
   H_c_new(j,i)   = H_c(j,i)
   H_t_new(j,i)   = H_t(j,i)

   call calc_temp_r(atr1, alb1, i, j)

end if

!  ------ Upper left corner

i=0
j=JMAX

if ( (mask(j,i) == 0).or.(mask(j,i) == 3) ) then
                              ! glaciated land or floating ice
   ii=i+1
   jj=j-1

   do kc=0,KCMAX
      temp_c_new(kc,j,i) = temp_c_new(kc,jj,ii)   ! set cold-ice temperature
      age_c_new(kc,j,i)  = age_c_new(kc,jj,ii)    ! set cold-ice age
   end do

   do kt=0,KTMAX
      omega_t_new(kt,j,i) = omega_t_new(kt,jj,ii) ! set temp.-ice water content
      age_t_new(kt,j,i)   = age_t_new(kt,jj,ii)   ! set temp.-ice age
   end do

   do kr=0,KRMAX
      temp_r_new(kr,j,i) = temp_r_new(kr,jj,ii)   ! set bedrock temperature
   end do

   n_cts_new(j,i) = min(n_cts_new(jj,ii),0)   ! temperate layer excluded
   H_c_new(j,i)   = H_c(j,i)
   H_t_new(j,i)   = H_t(j,i)
   zm_new(j,i)  = zb(j,i)

else   ! mask(j,i) == 1, 2 (ice-free land or sea point)

   n_cts_new(j,i) = -1
   zm_new(j,i)    = zb(j,i)
   H_c_new(j,i)   = H_c(j,i)
   H_t_new(j,i)   = H_t(j,i)

   call calc_temp_r(atr1, alb1, i, j)

end if

!  ------ Upper right corner

i=IMAX
j=JMAX

if ( (mask(j,i) == 0).or.(mask(j,i) == 3) ) then
                              ! glaciated land or floating ice
   ii=i-1
   jj=j-1

   do kc=0,KCMAX
      temp_c_new(kc,j,i) = temp_c_new(kc,jj,ii)   ! set cold-ice temperature
      age_c_new(kc,j,i)  = age_c_new(kc,jj,ii)    ! set cold-ice age
   end do

   do kt=0,KTMAX
      omega_t_new(kt,j,i) = omega_t_new(kt,jj,ii) ! set temp.-ice water content
      age_t_new(kt,j,i)   = age_t_new(kt,jj,ii)   ! set temp.-ice age
   end do

   do kr=0,KRMAX
      temp_r_new(kr,j,i) = temp_r_new(kr,jj,ii)   ! set bedrock temperature
   end do

   n_cts_new(j,i) = min(n_cts_new(jj,ii),0)   ! temperate layer excluded
   H_c_new(j,i)   = H_c(j,i)
   H_t_new(j,i)   = H_t(j,i)
   zm_new(j,i)  = zb(j,i)

else   ! mask(j,i) == 1, 2 (ice-free land or sea point)

   n_cts_new(j,i) = -1
   zm_new(j,i)    = zb(j,i)
   H_c_new(j,i)   = H_c(j,i)
   H_t_new(j,i)   = H_t(j,i)

   call calc_temp_r(atr1, alb1, i, j)

end if

!  ------ Lower and upper margins

do i=1, IMAX-1

!    ---- Lower margin

   j=0

   if ( (mask(j,i) == 0).or.(mask(j,i) == 3) ) then
                                 ! glaciated land or floating ice
      ii=i
      jj=j+1

      do kc=0,KCMAX
         temp_c_new(kc,j,i) = temp_c_new(kc,jj,ii)   ! set cold-ice temperature
         age_c_new(kc,j,i)  = age_c_new(kc,jj,ii)    ! set cold-ice age
      end do

      do kt=0,KTMAX
         omega_t_new(kt,j,i) = omega_t_new(kt,jj,ii) ! set temp.-ice water content
         age_t_new(kt,j,i)   = age_t_new(kt,jj,ii)   ! set temp.-ice age
      end do

      do kr=0,KRMAX
         temp_r_new(kr,j,i) = temp_r_new(kr,jj,ii)   ! set bedrock temperature
      end do

      n_cts_new(j,i) = min(n_cts_new(jj,ii),0)   ! temperate layer excluded
      H_c_new(j,i)   = H_c(j,i)
      H_t_new(j,i)   = H_t(j,i)
      zm_new(j,i)  = zb(j,i)

   else   ! mask(j,i) == 1, 2 (ice-free land or sea point)

      n_cts_new(j,i) = -1
      zm_new(j,i)    = zb(j,i)
      H_c_new(j,i)   = H_c(j,i)
      H_t_new(j,i)   = H_t(j,i)

      call calc_temp_r(atr1, alb1, i, j)

   end if

!    ---- Upper margin

   j=JMAX

   if ( (mask(j,i) == 0).or.(mask(j,i) == 3) ) then
                                 ! glaciated land or floating ice
      ii=i
      jj=j-1

      do kc=0,KCMAX
         temp_c_new(kc,j,i) = temp_c_new(kc,jj,ii)   ! set cold-ice temperature
         age_c_new(kc,j,i)  = age_c_new(kc,jj,ii)    ! set cold-ice age
      end do

      do kt=0,KTMAX
         omega_t_new(kt,j,i) = omega_t_new(kt,jj,ii) ! set temp.-ice water content
         age_t_new(kt,j,i)   = age_t_new(kt,jj,ii)   ! set temp.-ice age
      end do

      do kr=0,KRMAX
         temp_r_new(kr,j,i) = temp_r_new(kr,jj,ii)   ! set bedrock temperature
      end do

      n_cts_new(j,i) = min(n_cts_new(jj,ii),0)   ! temperate layer excluded
      H_c_new(j,i)   = H_c(j,i)
      H_t_new(j,i)   = H_t(j,i)
      zm_new(j,i)  = zb(j,i)

   else   ! mask(j,i) == 1, 2 (ice-free land or sea point)

      n_cts_new(j,i) = -1
      zm_new(j,i)    = zb(j,i)
      H_c_new(j,i)   = H_c(j,i)
      H_t_new(j,i)   = H_t(j,i)

      call calc_temp_r(atr1, alb1, i, j)

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

      do kc=0,KCMAX
         temp_c_new(kc,j,i) = temp_c_new(kc,jj,ii)   ! set cold-ice temperature
         age_c_new(kc,j,i)  = age_c_new(kc,jj,ii)    ! set cold-ice age
      end do

      do kt=0,KTMAX
         omega_t_new(kt,j,i) = omega_t_new(kt,jj,ii) ! set temp.-ice water content
         age_t_new(kt,j,i)   = age_t_new(kt,jj,ii)   ! set temp.-ice age
      end do

      do kr=0,KRMAX
         temp_r_new(kr,j,i) = temp_r_new(kr,jj,ii)   ! set bedrock temperature
      end do

      n_cts_new(j,i) = min(n_cts_new(jj,ii),0)   ! temperate layer excluded
      H_c_new(j,i)   = H_c(j,i)
      H_t_new(j,i)   = H_t(j,i)
      zm_new(j,i)  = zb(j,i)

   else   ! mask(j,i) == 1, 2 (ice-free land or sea point)

      n_cts_new(j,i) = -1
      zm_new(j,i)    = zb(j,i)
      H_c_new(j,i)   = H_c(j,i)
      H_t_new(j,i)   = H_t(j,i)

      call calc_temp_r(atr1, alb1, i, j)

   end if

!    ---- Right margin

   i=IMAX

   if ( (mask(j,i) == 0).or.(mask(j,i) == 3) ) then
                                 ! glaciated land or floating ice
      ii=i-1
      jj=j

      do kc=0,KCMAX
         temp_c_new(kc,j,i) = temp_c_new(kc,jj,ii)   ! set cold-ice temperature
         age_c_new(kc,j,i)  = age_c_new(kc,jj,ii)    ! set cold-ice age
      end do

      do kt=0,KTMAX
         omega_t_new(kt,j,i) = omega_t_new(kt,jj,ii) ! set temp.-ice water content
         age_t_new(kt,j,i)   = age_t_new(kt,jj,ii)   ! set temp.-ice age
      end do

      do kr=0,KRMAX
         temp_r_new(kr,j,i) = temp_r_new(kr,jj,ii)   ! set bedrock temperature
      end do

      n_cts_new(j,i) = min(n_cts_new(jj,ii),0)   ! temperate layer excluded
      H_c_new(j,i)   = H_c(j,i)
      H_t_new(j,i)   = H_t(j,i)
      zm_new(j,i)  = zb(j,i)

   else   ! mask(j,i) == 1, 2 (ice-free land or sea point)

      n_cts_new(j,i) = -1
      zm_new(j,i)    = zb(j,i)
      H_c_new(j,i)   = H_c(j,i)
      H_t_new(j,i)   = H_t(j,i)

      call calc_temp_r(atr1, alb1, i, j)

   end if

end do

!-------- Dummy values for omega_c_new and kc_cts_new --------

omega_c_new = 0.0_dp   ! not needed for
kc_cts_new  = 0        ! the polythermal mode

!-------- Smoothing of H_t_new with numerical diffusion --------

!  ------ Volume of temperate ice without smoothing

Vol_t = 0.0_dp
do i=0, IMAX   ! extended to domain margins (22.1.02 -> V1.1)
do j=0, JMAX   ! extended to domain margins (22.1.02 -> V1.1)
   if (n_cts_new(j,i) == 1) then
      Vol_t = Vol_t + H_t_new(j,i)*cell_area(j,i)
   end if
end do
end do

!  ------ Smoothing

do i=1, IMAX-1
do j=1, JMAX-1
   if (n_cts_new(j,i) /= -1) then

      dH_t_smooth(j,i) = NUMDIFF_H_T* ( -4.0_dp*H_t_new(j,i) &
                             +(H_t_new(j,i+1)+H_t_new(j,i-1)) &
                             +(H_t_new(j+1,i)+H_t_new(j-1,i)) )
      if (dH_t_smooth(j,i) > 0.001_dp) n_cts_new(j,i) = 1

   end if
end do
end do

do i=1, IMAX-1
do j=1, JMAX-1
   if (n_cts_new(j,i) == 1) then
      H_t_new(j,i) = H_t_new(j,i) + dH_t_smooth(j,i)
   end if
end do
end do

!  ------ Volume of temperate ice with smoothing

Vol_t_smooth = 0.0_dp
do i=0, IMAX
do j=0, JMAX
   if (n_cts_new(j,i) == 1) then
      Vol_t_smooth = Vol_t_smooth + H_t_new(j,i)*cell_area(j,i)
   end if
end do
end do

!  ------ Correction so that volume is not changed by the smoothing

if (Vol_t_smooth > 0.0_dp) then

   korrfakt_t = Vol_t/Vol_t_smooth
   do i=0, IMAX   ! extended to domain margins (22.1.02 -> V1.1)
   do j=0, JMAX   ! extended to domain margins (22.1.02 -> V1.1)
      if (n_cts_new(j,i) == 1) then
         H_t_new(j,i) = H_t_new(j,i)*korrfakt_t
!               zm_new(j,i)  = zb(j,i) + H_t_new(j,i)
!               H_c_new(j,i) = zs(j,i) - zm_new(j,i)
      end if
   end do
   end do

end if

!-------- Numerical time lag for evolution of H_t_new --------

time_lag_cts = TAU_CTS*year2sec   ! a -> s

do i=0, IMAX   ! extended to domain margins (22.1.02 -> V1.1)
do j=0, JMAX   ! extended to domain margins (22.1.02 -> V1.1)

   if (n_cts_new(j,i) == 1) then

      H_t_new(j,i) = ( time_lag_cts*H_t(j,i) &
                       + dtime_temp*H_t_new(j,i) ) &
                     /(time_lag_cts+dtime_temp)

      zm_new(j,i)  = zb(j,i) + H_t_new(j,i)
      H_c_new(j,i) = zs(j,i) - zm_new(j,i)

   end if

end do
end do

end subroutine calc_temp_poly

!-------------------------------------------------------------------------------
!> Computation of temperature and age in cold-ice mode.
!-------------------------------------------------------------------------------
subroutine calc_temp_cold(dxi, deta, dzeta_c, dzeta_t, dzeta_r, &
                          dtime_temp)

use ice_material_properties_m, only : kappa_val

implicit none

real(dp), intent(in) :: dxi, deta, dzeta_c, dzeta_t, dzeta_r
real(dp), intent(in) :: dtime_temp

integer(i4b) :: i, j, kc, kr, ii, jj
real(dp) :: at1(0:KCMAX), at2_1(0:KCMAX), at2_2(0:KCMAX), &
            at3_1(0:KCMAX), at3_2(0:KCMAX), at4_1(0:KCMAX), &
            at4_2(0:KCMAX), at5(0:KCMAX), at6(0:KCMAX), at7, &
            acb1, acb2, acb3, acb4, &
            ai1(0:KCMAX), ai2(0:KCMAX), ai3, &
            atr1, am1, am2, alb1
real(dp) :: aw1, aw2, aw3, aw4, aw5, aw7, aw8, aw9, aqtld
real(dp) :: dtime_temp_inv, dtt_2dxi, dtt_2deta
real(dp) :: temp_c_help(0:KCMAX)

!-------- Term abbreviations

at7 = 2.0_dp/RHO*dtime_temp

aw1 = dtime_temp/dzeta_t
aw2 = dtime_temp/dzeta_t
aw3 = dtime_temp/dzeta_t
aw4 = dtime_temp/dzeta_t
aw5 = NUE/RHO*dtime_temp/(dzeta_t**2)
aw7 = 2.0_dp/(RHO*L)*dtime_temp
aw8 = BETA**2/(RHO*L) &
      *(kappa_val(0.0_dp)-kappa_val(-1.0_dp))*dtime_temp
aw9 = BETA/L*dtime_temp

ai3 = AGEDIFF*dtime_temp/(dzeta_t**2)

atr1 = KAPPA_R/(RHO_C_R*H_R**2)*dtime_temp/(dzeta_r**2)

if (flag_aa_nonzero) then
   am1 = aa*BETA*dzeta_c/(ea-1.0_dp)
   am2 = aa*L*RHO*dzeta_c/(ea-1.0_dp)
else
   am1 = BETA*dzeta_c
   am2 = L*RHO*dzeta_c
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

aqtld = dzeta_t/dtime_temp

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

   end if

end do

strain_heating_c = 0.0_dp   ! initialization,
strain_heating_t = 0.0_dp   ! purely diagnostic fields

!-------- Computation loop for temperature and age --------

do i=1, IMAX-1   ! skipping domain margins
do j=1, JMAX-1   ! skipping domain margins

   if (mask(j,i)==0) then   ! glaciated land

      n_cts_new(j,i) = -1
      zm_new(j,i)  = zb(j,i)
      H_c_new(j,i) = H_c(j,i)
      H_t_new(j,i) = H_t(j,i)

      call calc_temp1(at1, at2_1, at2_2, at3_1, at3_2, &
           at4_1, at4_2, at5, at6, at7, atr1, acb1, acb2, &
           acb3, acb4, alb1, ai1, ai2, &
           dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
           i, j)

!  ------ Reset temperatures above melting to the melting point,
!         look for the CTS

      kc_cts_new(j,i) = 0

      if (temp_c_new(0,j,i) > temp_c_m(0,j,i)) then
         n_cts_new(j,i)        = 0
         kc_cts_new(j,i)       = 0
         temp_c_new(0,j,i)     = temp_c_m(0,j,i)
         temp_r_new(KRMAX,j,i) = temp_c_m(0,j,i)
      end if

      do kc=1, KCMAX
         if (temp_c_new(kc,j,i) > temp_c_m(kc,j,i)) then
            kc_cts_new(j,i)    = kc
            temp_c_new(kc,j,i) = temp_c_m(kc,j,i)
         end if
      end do

#if (MARGIN==3)

   else if (mask(j,i)==3) then   ! floating ice

      n_cts_new(j,i)  = -1
      kc_cts_new(j,i) =  0
      zm_new(j,i)     = zb(j,i)
      H_c_new(j,i)    = H_c(j,i)
      H_t_new(j,i)    = 0.0_dp

      call calc_temp_ssa(at1, at2_1, at2_2, at3_1, at3_2, &
           at4_1, at4_2, at5, at6, at7, atr1, alb1, &
           ai1, ai2, &
           dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
           i, j)

!  ------ Reset temperatures above melting to the melting point
!         (should not occur, but just in case)

      do kc=0, KCMAX
         if (temp_c_new(kc,j,i) > temp_c_m(kc,j,i)) &
                    temp_c_new(kc,j,i) = temp_c_m(kc,j,i)
      end do

#endif

   else   ! mask(j,i) == 1, 2 (ice-free land or sea point)

      n_cts_new(j,i)  = -1
      kc_cts_new(j,i) =  0
      zm_new(j,i)     = zb(j,i)
      H_c_new(j,i)    = H_c(j,i)
      H_t_new(j,i)    = H_t(j,i)

      call calc_temp_r(atr1, alb1, i, j)

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

   do kc=0,KCMAX
      temp_c_new(kc,j,i) = temp_c_new(kc,jj,ii)   ! set cold-ice temperature
      age_c_new(kc,j,i)  = age_c_new(kc,jj,ii)    ! set cold-ice age
   end do

   do kr=0,KRMAX
      temp_r_new(kr,j,i) = temp_r_new(kr,jj,ii)   ! set bedrock temperature
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

   call calc_temp_r(atr1, alb1, i, j)

end if

!  ------ Lower right corner

i=IMAX
j=0

if ( (mask(j,i) == 0).or.(mask(j,i) == 3) ) then
                              ! glaciated land or floating ice
   ii=i-1
   jj=j+1

   do kc=0,KCMAX
      temp_c_new(kc,j,i) = temp_c_new(kc,jj,ii)   ! set cold-ice temperature
      age_c_new(kc,j,i)  = age_c_new(kc,jj,ii)    ! set cold-ice age
   end do

   do kr=0,KRMAX
      temp_r_new(kr,j,i) = temp_r_new(kr,jj,ii)   ! set bedrock temperature
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

   call calc_temp_r(atr1, alb1, i, j)

end if

!  ------ Upper left corner

i=0
j=JMAX

if ( (mask(j,i) == 0).or.(mask(j,i) == 3) ) then
                              ! glaciated land or floating ice
   ii=i+1
   jj=j-1

   do kc=0,KCMAX
      temp_c_new(kc,j,i) = temp_c_new(kc,jj,ii)   ! set cold-ice temperature
      age_c_new(kc,j,i)  = age_c_new(kc,jj,ii)    ! set cold-ice age
   end do

   do kr=0,KRMAX
      temp_r_new(kr,j,i) = temp_r_new(kr,jj,ii)   ! set bedrock temperature
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

   call calc_temp_r(atr1, alb1, i, j)

end if

!  ------ Upper right corner

i=IMAX
j=JMAX

if ( (mask(j,i) == 0).or.(mask(j,i) == 3) ) then
                              ! glaciated land or floating ice
   ii=i-1
   jj=j-1

   do kc=0,KCMAX
      temp_c_new(kc,j,i) = temp_c_new(kc,jj,ii)   ! set cold-ice temperature
      age_c_new(kc,j,i)  = age_c_new(kc,jj,ii)    ! set cold-ice age
   end do

   do kr=0,KRMAX
      temp_r_new(kr,j,i) = temp_r_new(kr,jj,ii)   ! set bedrock temperature
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

   call calc_temp_r(atr1, alb1, i, j)

end if

!  ------ Lower and upper margins

do i=1, IMAX-1

!    ---- Lower margin

   j=0

   if ( (mask(j,i) == 0).or.(mask(j,i) == 3) ) then
                                 ! glaciated land or floating ice
      ii=i
      jj=j+1

      do kc=0,KCMAX
         temp_c_new(kc,j,i) = temp_c_new(kc,jj,ii)   ! set cold-ice temperature
         age_c_new(kc,j,i)  = age_c_new(kc,jj,ii)    ! set cold-ice age
      end do

      do kr=0,KRMAX
         temp_r_new(kr,j,i) = temp_r_new(kr,jj,ii)   ! set bedrock temperature
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

      call calc_temp_r(atr1, alb1, i, j)

   end if

!    ---- Upper margin

   j=JMAX

   if ( (mask(j,i) == 0).or.(mask(j,i) == 3) ) then
                                 ! glaciated land or floating ice
      ii=i
      jj=j-1

      do kc=0,KCMAX
         temp_c_new(kc,j,i) = temp_c_new(kc,jj,ii)   ! set cold-ice temperature
         age_c_new(kc,j,i)  = age_c_new(kc,jj,ii)    ! set cold-ice age
      end do

      do kr=0,KRMAX
         temp_r_new(kr,j,i) = temp_r_new(kr,jj,ii)   ! set bedrock temperature
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

      call calc_temp_r(atr1, alb1, i, j)

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

      do kc=0,KCMAX
         temp_c_new(kc,j,i) = temp_c_new(kc,jj,ii)   ! set cold-ice temperature
         age_c_new(kc,j,i)  = age_c_new(kc,jj,ii)    ! set cold-ice age
      end do

      do kr=0,KRMAX
         temp_r_new(kr,j,i) = temp_r_new(kr,jj,ii)   ! set bedrock temperature
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

      call calc_temp_r(atr1, alb1, i, j)

   end if

!    ---- Right margin

   i=IMAX

   if ( (mask(j,i) == 0).or.(mask(j,i) == 3) ) then
                                 ! glaciated land or floating ice
      ii=i-1
      jj=j

      do kc=0,KCMAX
         temp_c_new(kc,j,i) = temp_c_new(kc,jj,ii)   ! set cold-ice temperature
         age_c_new(kc,j,i)  = age_c_new(kc,jj,ii)    ! set cold-ice age
      end do

      do kr=0,KRMAX
         temp_r_new(kr,j,i) = temp_r_new(kr,jj,ii)   ! set bedrock temperature
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

      call calc_temp_r(atr1, alb1, i, j)

   end if

end do

!-------- Dummy values for omega_c_new --------

omega_c_new = 0.0_dp   ! not computed in the cold-ice mode

end subroutine calc_temp_cold

!-------------------------------------------------------------------------------
!> Isothermal mode: Setting of the temperature and age to constant values.
!-------------------------------------------------------------------------------
subroutine calc_temp_const()

implicit none

#if (defined(TEMP_CONST))
   if ( TEMP_CONST > -eps ) then
      errormsg = ' >>> calc_temp_const: TEMP_CONST must be negative!'
      call error(errormsg)
   end if
   temp_c_new  = TEMP_CONST
   temp_r_new  = TEMP_CONST
#else
   temp_c_new  = -10.0_dp   ! default value -10 C
   temp_r_new  = -10.0_dp   ! default value -10 C
#endif

temp_c_new = min(temp_c_new, temp_c_m-eps)
             ! keep temperatures below the pressure melting point

omega_t_new = 0.0_dp
omega_c_new = 0.0_dp

Q_tld       = 0.0_dp

#if (defined(AGE_CONST))
   age_c_new   = AGE_CONST *year2sec   ! a -> s
   age_t_new   = AGE_CONST *year2sec   ! a -> s
#else
   age_c_new   = 0.0_dp   ! default value 0
   age_t_new   = 0.0_dp   ! default value 0
#endif

strain_heating_c = 0.0_dp   ! purely diagnostic fields,
strain_heating_t = 0.0_dp   ! not computed in the isothermal mode

n_cts_new   = -1
kc_cts_new  =  0
zm_new      = zb
H_c_new     = H_c
H_t_new     = 0.0_dp

end subroutine calc_temp_const

!-------------------------------------------------------------------------------
!> Computation of temperature and age for a cold ice column.
!-------------------------------------------------------------------------------
subroutine calc_temp1(at1, at2_1, at2_2, at3_1, at3_2, &
   at4_1, at4_2, at5, at6, at7, atr1, acb1, acb2, &
   acb3, acb4, alb1, ai1, ai2, &
   dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
   i, j)

use ice_material_properties_m, only : ratefac_c, kappa_val, c_val, &
                                      creep, viscosity

#if !defined(ALLOW_TAPENADE) /* Normal */
use sico_maths_m, only : tri_sle
#else /* Tapenade */
use sico_maths_m
#endif /* Normal vs. Tapenade */

implicit none

integer(i4b), intent(in)    :: i, j
real(dp), intent(in) :: at1(0:KCMAX), at2_1(0:KCMAX), at2_2(0:KCMAX), &
                        at3_1(0:KCMAX), at3_2(0:KCMAX), at4_1(0:KCMAX), &
                        at4_2(0:KCMAX), at5(0:KCMAX), at6(0:KCMAX), at7, &
                        ai1(0:KCMAX), ai2(0:KCMAX), &
                        atr1, acb1, acb2, acb3, acb4, alb1
real(dp), intent(in) :: dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv

integer(i4b) :: kc, kt, kr
real(dp) :: ct1(0:KCMAX), ct2(0:KCMAX), ct3(0:KCMAX), ct4(0:KCMAX), &
            ct5(0:KCMAX), ct6(0:KCMAX), ct7(0:KCMAX), ctr1, &
            ccb1, ccb2, ccb3, ccb4, clb1
real(dp) :: ct1_sg(0:KCMAX), ct2_sg(0:KCMAX), ct3_sg(0:KCMAX), &
            ct4_sg(0:KCMAX), adv_vert_sg(0:KCMAX), abs_adv_vert_sg(0:KCMAX)
real(dp) :: ci1(0:KCMAX), ci2(0:KCMAX)
real(dp) :: temp_c_help(0:KCMAX)
real(dp) :: vx_c_help, vy_c_help
real(dp) :: adv_vert_help
real(dp) :: dtt_dxi, dtt_deta
real(dp) :: c_val_aux
real(dp) :: lgs_a0(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_a1(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_a2(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_x(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_b(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX)
real(dp), parameter :: zero=0.0_dp

!-------- Check for boundary points --------

if ((i == 0).or.(i == IMAX).or.(j == 0).or.(j == JMAX)) then
   errormsg = ' >>> calc_temp1: Boundary points not allowed!'
   call error(errormsg)
end if

!-------- Abbreviations --------

ctr1 = atr1

ccb1 = acb1 &
   *kappa_val(temp_c(0,j,i)) &
   /H_c(j,i)
ccb2 = acb2

#if (DYNAMICS==2)
if (.not.flag_shelfy_stream(j,i)) then
#endif

   ccb3 = acb3*0.5_dp*(vx_t(0,j,i)+vx_t(0,j,i-1)) &
              *H_c(j,i)*dzs_dxi_g(j,i)
   ccb4 = acb4*0.5_dp*(vy_t(0,j,i)+vy_t(0,j-1,i)) &
              *H_c(j,i)*dzs_deta_g(j,i)

#if (DYNAMICS==2)
else   ! flag_shelfy_stream(j,i) == .true.

   ccb3 = -c_drag(j,i) &
           * sqrt(vx_b_g(j,i)**2  &
                 +vy_b_g(j,i)**2) &
                           **(1.0_dp+p_weert_inv(j,i))
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

   c_val_aux = c_val(temp_c(kc,j,i))

   ct2(kc) = ( at2_1(kc)*dzm_dtau(j,i) &
           +at2_2(kc)*dH_c_dtau(j,i) )/H_c(j,i)
   ct3(kc) = ( at3_1(kc)*dzm_dxi_g(j,i) &
           +at3_2(kc)*dH_c_dxi_g(j,i) )/H_c(j,i) &
          *0.5_dp*(vx_c(kc,j,i)+vx_c(kc,j,i-1)) *insq_g11_g(j,i)
   ct4(kc) = ( at4_1(kc)*dzm_deta_g(j,i) &
            +at4_2(kc)*dH_c_deta_g(j,i) )/H_c(j,i) &
          *0.5_dp*(vy_c(kc,j,i)+vy_c(kc,j-1,i)) *insq_g22_g(j,i)
   ct5(kc) = at5(kc) &
             /c_val_aux &
             /H_c(j,i)

#if (DYNAMICS==2)
   if (.not.flag_shelfy_stream(j,i)) then
#endif

      ct7(kc) = at7 &
                /c_val_aux &
                *enh_c(kc,j,i) &
                *ratefac_c(temp_c(kc,j,i), temp_c_m(kc,j,i)) &
                *creep(sigma_c(kc,j,i)) &
                *sigma_c(kc,j,i)*sigma_c(kc,j,i)

#if (DYNAMICS==2)
   else
      ct7(kc) = 2.0_dp*at7 &
                /c_val_aux &
                *viscosity(de_c(kc,j,i), &
                           temp_c(kc,j,i), temp_c_m(kc,j,i), 0.0_dp, &
                           enh_c(kc,j,i), 0) &
                *de_c(kc,j,i)**2
   end if
#endif

   strain_heating_c(kc,j,i) = c_val_aux*ct7(kc)*dtime_temp_inv

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
   ct6(kc) = at6(kc) &
    *kappa_val(temp_c_help(kc)) &
    /H_c(j,i)
   ci2(kc) = ai2(kc)/H_c(j,i)
end do

#if (ADV_HOR==3)
dtt_dxi  = 2.0_dp*dtt_2dxi
dtt_deta = 2.0_dp*dtt_2deta
#endif

!-------- Set up the temperature equations (ice and bedrock
!         simultaneously) --------

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
kc=0
lgs_a0(kr) = ccb2
lgs_a1(kr) = -(ccb1+ccb2)
lgs_a2(kr) = ccb1
lgs_b(kr)  = ccb3+ccb4

do kc=1, KCMAX-1

#if (ADV_VERT==1)

   lgs_a0(KRMAX+kc) = -0.5_dp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ct5(kc)*ct6(kc-1)
   lgs_a1(KRMAX+kc) = 1.0_dp+ct5(kc)*(ct6(kc)+ct6(kc-1))
   lgs_a2(KRMAX+kc) = 0.5_dp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ct5(kc)*ct6(kc)

#elif (ADV_VERT==2)

   lgs_a0(KRMAX+kc) &
         = -0.5_dp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -ct5(kc)*ct6(kc-1)
   lgs_a1(KRMAX+kc) &
         = 1.0_dp &
           +0.5_dp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -0.5_dp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           +ct5(kc)*(ct6(kc)+ct6(kc-1))
   lgs_a2(KRMAX+kc) &
         =  0.5_dp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           -ct5(kc)*ct6(kc)

#elif (ADV_VERT==3)

   adv_vert_help = 0.5_dp*(adv_vert_sg(kc)+adv_vert_sg(kc-1))

   lgs_a0(KRMAX+kc) &
         = -max(adv_vert_help, 0.0_dp) &
           -ct5(kc)*ct6(kc-1)
   lgs_a1(KRMAX+kc) &
         = 1.0_dp &
           +max(adv_vert_help, 0.0_dp)-min(adv_vert_help, 0.0_dp) &
           +ct5(kc)*(ct6(kc)+ct6(kc-1))
   lgs_a2(KRMAX+kc) &
         =  min(adv_vert_help, 0.0_dp) &
           -ct5(kc)*ct6(kc)

#endif

#if (ADV_HOR==2)

   lgs_b(KRMAX+kc) = temp_c(kc,j,i) + ct7(kc) &
       - ( dtt_2dxi* &
             ( (vx_c(kc,j,i)-abs(vx_c(kc,j,i))) &
               *(temp_c(kc,j,i+1)-temp_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +(vx_c(kc,j,i-1)+abs(vx_c(kc,j,i-1))) &
               *(temp_c(kc,j,i)-temp_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_2deta* &
             ( (vy_c(kc,j,i)-abs(vy_c(kc,j,i))) &
               *(temp_c(kc,j+1,i)-temp_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +(vy_c(kc,j-1,i)+abs(vy_c(kc,j-1,i))) &
               *(temp_c(kc,j,i)-temp_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#elif (ADV_HOR==3)

   vx_c_help = 0.5_dp*(vx_c(kc,j,i)+vx_c(kc,j,i-1))
   vy_c_help = 0.5_dp*(vy_c(kc,j,i)+vy_c(kc,j-1,i))

   lgs_b(KRMAX+kc) = temp_c(kc,j,i) + ct7(kc) &
       - ( dtt_dxi* &
             ( min(vx_c_help, 0.0_dp) &
               *(temp_c(kc,j,i+1)-temp_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +max(vx_c_help, 0.0_dp) &
               *(temp_c(kc,j,i)-temp_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_deta* &
             ( min(vy_c_help, 0.0_dp) &
               *(temp_c(kc,j+1,i)-temp_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +max(vy_c_help, 0.0_dp) &
               *(temp_c(kc,j,i)-temp_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#endif

end do

kc=KCMAX
lgs_a0(KRMAX+kc) = 0.0_dp
lgs_a1(KRMAX+kc) = 1.0_dp
lgs_b(KRMAX+kc)  = temp_s(j,i)

!-------- Solve system of linear equations --------

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, KCMAX+KRMAX)

!-------- Assign the result --------

do kr=0, KRMAX
   temp_r_new(kr,j,i) = lgs_x(kr)
end do

do kc=0, KCMAX
   temp_c_new(kc,j,i) = lgs_x(KRMAX+kc)
end do

!-------- Set water content in the non-existing temperate layer
!         to zero --------

do kt=0, KTMAX
   omega_t_new(kt,j,i) = 0.0_dp
end do

!-------- Water drainage from the non-existing temperate layer --------

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

!-------- Age of the ice in the non-existing temperate layer --------

do kt=0, KTMAX
   age_t_new(kt,j,i) = age_c_new(0,j,i)
end do

end subroutine calc_temp1

!-------------------------------------------------------------------------------
!> Computation of temperature and age for an ice column with a temperate base
!! overlain by cold ice.
!-------------------------------------------------------------------------------
subroutine calc_temp2(at1, at2_1, at2_2, at3_1, at3_2, &
   at4_1, at4_2, at5, at6, at7, atr1, alb1, &
   ai1, ai2, &
   dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
   i, j)

use ice_material_properties_m, only : ratefac_c, kappa_val, c_val, &
                                      creep, viscosity
#if !defined(ALLOW_TAPENADE) /* Normal */
use sico_maths_m, only : tri_sle
#else /* Tapenade */
use sico_maths_m
#endif /* Normal vs. Tapenade */

implicit none

integer(i4b), intent(in)    :: i, j
real(dp), intent(in) :: at1(0:KCMAX), at2_1(0:KCMAX), at2_2(0:KCMAX), &
                        at3_1(0:KCMAX), at3_2(0:KCMAX), at4_1(0:KCMAX), &
                        at4_2(0:KCMAX), at5(0:KCMAX), at6(0:KCMAX), at7, &
                        ai1(0:KCMAX), ai2(0:KCMAX), &
                        atr1, alb1
real(dp), intent(in) :: dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv

integer(i4b) :: kc, kt, kr
real(dp) :: ct1(0:KCMAX), ct2(0:KCMAX), ct3(0:KCMAX), ct4(0:KCMAX), &
            ct5(0:KCMAX), ct6(0:KCMAX), ct7(0:KCMAX), ctr1, clb1
real(dp) :: ct1_sg(0:KCMAX), ct2_sg(0:KCMAX), ct3_sg(0:KCMAX), &
            ct4_sg(0:KCMAX), adv_vert_sg(0:KCMAX), abs_adv_vert_sg(0:KCMAX)
real(dp) :: ci1(0:KCMAX), ci2(0:KCMAX)
real(dp) :: temp_c_help(0:KCMAX)
real(dp) :: vx_c_help, vy_c_help
real(dp) :: adv_vert_help
real(dp) :: dtt_dxi, dtt_deta
real(dp) :: c_val_aux
real(dp) :: lgs_a0(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_a1(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_a2(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_x(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_b(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX)
real(dp), parameter :: zero=0.0_dp

!-------- Check for boundary points --------

if ((i == 0).or.(i == IMAX).or.(j == 0).or.(j == JMAX)) then
   errormsg = ' >>> calc_temp2: Boundary points not allowed!'
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

   c_val_aux = c_val(temp_c(kc,j,i))

   ct2(kc) = ( at2_1(kc)*dzm_dtau(j,i) &
           +at2_2(kc)*dH_c_dtau(j,i) )/H_c(j,i)
   ct3(kc) = ( at3_1(kc)*dzm_dxi_g(j,i) &
           +at3_2(kc)*dH_c_dxi_g(j,i) )/H_c(j,i) &
          *0.5_dp*(vx_c(kc,j,i)+vx_c(kc,j,i-1)) *insq_g11_g(j,i)
   ct4(kc) = ( at4_1(kc)*dzm_deta_g(j,i) &
            +at4_2(kc)*dH_c_deta_g(j,i) )/H_c(j,i) &
          *0.5_dp*(vy_c(kc,j,i)+vy_c(kc,j-1,i)) *insq_g22_g(j,i)
   ct5(kc) = at5(kc) &
             /c_val_aux &
             /H_c(j,i)

#if (DYNAMICS==2)
   if (.not.flag_shelfy_stream(j,i)) then
#endif
      ct7(kc) = at7 &
                /c_val_aux &
                *enh_c(kc,j,i) &
                *ratefac_c(temp_c(kc,j,i), temp_c_m(kc,j,i)) &
                *creep(sigma_c(kc,j,i)) &
                *sigma_c(kc,j,i)*sigma_c(kc,j,i)
#if (DYNAMICS==2)
   else
      ct7(kc) = 2.0_dp*at7 &
                /c_val_aux &
                *viscosity(de_c(kc,j,i), &
                           temp_c(kc,j,i), temp_c_m(kc,j,i), 0.0_dp, &
                           enh_c(kc,j,i), 0) &
                *de_c(kc,j,i)**2
   end if
#endif

   strain_heating_c(kc,j,i) = c_val_aux*ct7(kc)*dtime_temp_inv

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
   ct6(kc) = at6(kc) &
    *kappa_val(temp_c_help(kc)) &
    /H_c(j,i)
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
lgs_b(kr)   = temp_t_m(0,j,i)

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
lgs_b(kc)  = temp_c_m(0,j,i)

do kc=1, KCMAX-1

#if (ADV_VERT==1)

   lgs_a0(kc) = -0.5_dp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ct5(kc)*ct6(kc-1)
   lgs_a1(kc) = 1.0_dp+ct5(kc)*(ct6(kc)+ct6(kc-1))
   lgs_a2(kc) = 0.5_dp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ct5(kc)*ct6(kc)

#elif (ADV_VERT==2)

   lgs_a0(kc) &
         = -0.5_dp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -ct5(kc)*ct6(kc-1)
   lgs_a1(kc) &
         = 1.0_dp &
           +0.5_dp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -0.5_dp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           +ct5(kc)*(ct6(kc)+ct6(kc-1))
   lgs_a2(kc) &
         =  0.5_dp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           -ct5(kc)*ct6(kc)

#elif (ADV_VERT==3)

   adv_vert_help = 0.5_dp*(adv_vert_sg(kc)+adv_vert_sg(kc-1))

   lgs_a0(kc) &
         = -max(adv_vert_help, 0.0_dp) &
           -ct5(kc)*ct6(kc-1)
   lgs_a1(kc) &
         = 1.0_dp &
           +max(adv_vert_help, 0.0_dp)-min(adv_vert_help, 0.0_dp) &
           +ct5(kc)*(ct6(kc)+ct6(kc-1))
   lgs_a2(kc) &
         =  min(adv_vert_help, 0.0_dp) &
           -ct5(kc)*ct6(kc)

#endif

#if (ADV_HOR==2)

   lgs_b(kc) = temp_c(kc,j,i) + ct7(kc) &
       - ( dtt_2dxi* &
             ( (vx_c(kc,j,i)-abs(vx_c(kc,j,i))) &
               *(temp_c(kc,j,i+1)-temp_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +(vx_c(kc,j,i-1)+abs(vx_c(kc,j,i-1))) &
               *(temp_c(kc,j,i)-temp_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_2deta* &
             ( (vy_c(kc,j,i)-abs(vy_c(kc,j,i))) &
               *(temp_c(kc,j+1,i)-temp_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +(vy_c(kc,j-1,i)+abs(vy_c(kc,j-1,i))) &
               *(temp_c(kc,j,i)-temp_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#elif (ADV_HOR==3)

   vx_c_help = 0.5_dp*(vx_c(kc,j,i)+vx_c(kc,j,i-1))
   vy_c_help = 0.5_dp*(vy_c(kc,j,i)+vy_c(kc,j-1,i))

   lgs_b(kc) = temp_c(kc,j,i) + ct7(kc) &
       - ( dtt_dxi* &
             ( min(vx_c_help, 0.0_dp) &
               *(temp_c(kc,j,i+1)-temp_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +max(vx_c_help, 0.0_dp) &
               *(temp_c(kc,j,i)-temp_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_deta* &
             ( min(vy_c_help, 0.0_dp) &
               *(temp_c(kc,j+1,i)-temp_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +max(vy_c_help, 0.0_dp) &
               *(temp_c(kc,j,i)-temp_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#endif

end do

kc=KCMAX
lgs_a0(kc) = 0.0_dp
lgs_a1(kc) = 1.0_dp
lgs_b(kc)  = temp_s(j,i)

!-------- Solve system of linear equations --------

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, KCMAX)

!-------- Assign the result --------

do kc=0, KCMAX
   temp_c_new(kc,j,i) = lgs_x(kc)
end do

!-------- Set water content in the non-existing temperate layer
!         to zero --------

do kt=0, KTMAX
   omega_t_new(kt,j,i) = 0.0_dp
end do

!-------- Water drainage from the non-existing temperate layer --------

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

!-------- Age of the ice in the non-existing temperate layer --------

do kt=0, KTMAX
   age_t_new(kt,j,i) = age_c_new(0,j,i)
end do

end subroutine calc_temp2

!-------------------------------------------------------------------------------
!> Computation of temperature, water content and age for an ice column with a
!! temperate base overlain by a temperate-ice layer.
!-------------------------------------------------------------------------------
subroutine calc_temp3(at1, at2_1, at2_2, at3_1, at3_2, &
   at4_1, at4_2, at5, at6, at7, atr1, am1, am2, alb1, &
   aw1, aw2, aw3, aw4, aw5, aw7, aw8, aw9, aqtld, &
   ai1, ai2, ai3, dzeta_t, &
   dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
   i, j)

use ice_material_properties_m, only : ratefac_c, ratefac_t, kappa_val, c_val, &
                                      creep, viscosity
#if !defined(ALLOW_TAPENADE) /* Normal */
use sico_maths_m, only : tri_sle
#else /* Tapenade */
use sico_maths_m
#endif /* Normal vs. Tapenade */

implicit none

integer(i4b), intent(in)    :: i, j
real(dp), intent(in) :: at1(0:KCMAX), at2_1(0:KCMAX), at2_2(0:KCMAX), &
                        at3_1(0:KCMAX), at3_2(0:KCMAX), at4_1(0:KCMAX), &
                        at4_2(0:KCMAX), at5(0:KCMAX), at6(0:KCMAX), at7, &
                        ai1(0:KCMAX), ai2(0:KCMAX), ai3, &
                        atr1, am1, am2, alb1
real(dp), intent(in) :: aw1, aw2, aw3, aw4, aw5, aw7, aw8, aw9, aqtld
real(dp), intent(in) :: dzeta_t, dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv

integer(i4b) :: kc, kt, kr
real(dp) :: ct1(0:KCMAX), ct2(0:KCMAX), ct3(0:KCMAX), ct4(0:KCMAX), &
            ct5(0:KCMAX), ct6(0:KCMAX), ct7(0:KCMAX), ctr1, cm1, cm2, &
            clb1
real(dp) :: ct1_sg(0:KCMAX), ct2_sg(0:KCMAX), ct3_sg(0:KCMAX), &
            ct4_sg(0:KCMAX), adv_vert_sg(0:KCMAX), abs_adv_vert_sg(0:KCMAX)
real(dp) :: ci1(0:KCMAX), ci2(0:KCMAX), ci3
real(dp) :: cw1(0:KTMAX), cw2(0:KTMAX), cw3(0:KTMAX), cw4(0:KTMAX), &
            cw5, cw7(0:KTMAX), cw8, cw9(0:KTMAX)
real(dp) :: cw1_sg(0:KTMAX), cw2_sg(0:KTMAX), cw3_sg(0:KTMAX), &
            cw4_sg(0:KTMAX), adv_vert_w_sg(0:KTMAX), abs_adv_vert_w_sg(0:KTMAX)
real(dp) :: sigma_c_help(0:KCMAX), sigma_t_help(0:KTMAX), &
            temp_c_help(0:KCMAX)
real(dp) :: vx_c_help, vy_c_help, vx_t_help, vy_t_help
real(dp) :: adv_vert_help, adv_vert_w_help
real(dp) :: dtt_dxi, dtt_deta
real(dp) :: c_val_aux
real(dp) :: lgs_a0(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_a1(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_a2(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_x(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_b(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX)
real(dp), parameter :: zero=0.0_dp

!-------- Check for boundary points --------

if ((i == 0).or.(i == IMAX).or.(j == 0).or.(j == JMAX)) then
   errormsg = ' >>> calc_temp3: Boundary points not allowed!'
   call error(errormsg)
end if

!-------- Abbreviations --------

ctr1 = atr1
cm1  = am1*H_c_new(j,i)
clb1 = alb1*q_geo(j,i)

#if (ADV_VERT==1)

do kc=1, KCMAX-1
   ct1(kc) = at1(kc)/H_c_new(j,i)*0.5_dp*(vz_c(kc,j,i)+vz_c(kc-1,j,i))
end do

kc=0
ct1_sg(kc) = 0.5_dp*(at1(kc)+at1(kc+1))/H_c_new(j,i)*vz_c(kc,j,i)
             ! only needed for kc=0 ...
kc=KCMAX-1
ct1_sg(kc) = 0.5_dp*(at1(kc)+at1(kc+1))/H_c_new(j,i)*vz_c(kc,j,i)
             ! ... and kc=KCMAX-1

#elif (ADV_VERT==2 || ADV_VERT==3)

do kc=0, KCMAX-1
   ct1_sg(kc) = 0.5_dp*(at1(kc)+at1(kc+1))/H_c_new(j,i)*vz_c(kc,j,i)
end do

#endif

do kc=0, KCMAX

   c_val_aux = c_val(temp_c(kc,j,i))

   ct2(kc) = ( at2_1(kc)*dzm_dtau(j,i) &
           +at2_2(kc)*dH_c_dtau(j,i) )/H_c_new(j,i)
   ct3(kc) = ( at3_1(kc)*dzm_dxi_g(j,i) &
           +at3_2(kc)*dH_c_dxi_g(j,i) )/H_c_new(j,i) &
          *0.5_dp*(vx_c(kc,j,i)+vx_c(kc,j,i-1)) *insq_g11_g(j,i)
   ct4(kc) = ( at4_1(kc)*dzm_deta_g(j,i) &
            +at4_2(kc)*dH_c_deta_g(j,i) )/H_c_new(j,i) &
          *0.5_dp*(vy_c(kc,j,i)+vy_c(kc,j-1,i)) *insq_g22_g(j,i)
   ct5(kc) = at5(kc) &
             /c_val_aux &
             /H_c_new(j,i)

   sigma_c_help(kc) &
           = RHO*G*H_c_new(j,i)*(1.0_dp-eaz_c_quotient(kc)) &
             *(dzs_dxi_g(j,i)**2+dzs_deta_g(j,i)**2)**0.5_dp

#if (DYNAMICS==2)
   if (.not.flag_shelfy_stream(j,i)) then
#endif
      ct7(kc) = at7 &
                /c_val_aux &
                *enh_c(kc,j,i) &
                *ratefac_c(temp_c(kc,j,i), temp_c_m(kc,j,i)) &
                *creep(sigma_c_help(kc)) &
                *sigma_c_help(kc)*sigma_c_help(kc)
#if (DYNAMICS==2)
   else
      ct7(kc) = 2.0_dp*at7 &
                /c_val_aux &
                *viscosity(de_c(kc,j,i), &
                           temp_c(kc,j,i), temp_c_m(kc,j,i), 0.0_dp, &
                           enh_c(kc,j,i), 0) &
                *de_c(kc,j,i)**2
   end if
#endif

   strain_heating_c(kc,j,i) = c_val_aux*ct7(kc)*dtime_temp_inv

   ci1(kc) = ai1(kc)/H_c_new(j,i)

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
   ct6(kc) = at6(kc) &
    *kappa_val(temp_c_help(kc)) &
    /H_c_new(j,i)
   ci2(kc) = ai2(kc)/H_c_new(j,i)
end do

cw5 = aw5/(H_t_new(j,i)**2)
cw8 = aw8
ci3 = ai3/(H_t_new(j,i)**2)

#if (ADV_VERT==1)

do kt=1, KTMAX-1
   cw1(kt) = aw1/H_t_new(j,i)*0.5_dp*(vz_t(kt,j,i)+vz_t(kt-1,j,i))
end do

kt=KTMAX
cw1(kt) = aw1/H_t_new(j,i)*0.5_dp*(vz_t(kt-1,j,i)+vz_c(0,j,i))

kt=0
cw1_sg(kt) = aw1/H_t_new(j,i)*vz_t(kt,j,i)
             ! only needed for kt=0 ...
kt=KTMAX-1
cw1_sg(kt) = aw1/H_t_new(j,i)*vz_t(kt,j,i)
             ! ... and kt=KTMAX-1

#elif (ADV_VERT==2 || ADV_VERT==3)

do kt=0, KTMAX-1
   cw1_sg(kt) = aw1/H_t_new(j,i)*vz_t(kt,j,i)
end do

#endif

do kt=1, KTMAX-1
   cw9(kt) = aw9 &
             *c_val(temp_t_m(kt,j,i)) &
    *( dzs_dtau(j,i) &
      + ( 0.5_dp*(vx_t(kt,j,i)+vx_t(kt,j,i-1))*dzs_dxi_g(j,i) &
         +0.5_dp*(vy_t(kt,j,i)+vy_t(kt,j-1,i))*dzs_deta_g(j,i) ) &
      -0.5_dp*(vz_t(kt,j,i)+vz_t(kt-1,j,i)) )
end do

do kt=0, KTMAX

   cw2(kt) = aw2*(dzb_dtau(j,i)+zeta_t(kt)*dH_t_dtau(j,i)) &
             /H_t_new(j,i)
   cw3(kt) = aw3*(dzb_dxi_g(j,i)+zeta_t(kt)*dH_t_dxi_g(j,i)) &
             /H_t_new(j,i) &
             *0.5_dp*(vx_t(kt,j,i)+vx_t(kt,j,i-1)) *insq_g11_g(j,i)
   cw4(kt) = aw4*(dzb_deta_g(j,i)+zeta_t(kt)*dH_t_deta_g(j,i)) &
             /H_t_new(j,i) &
             *0.5_dp*(vy_t(kt,j,i)+vy_t(kt,j-1,i)) *insq_g22_g(j,i)
   sigma_t_help(kt) &
           = sigma_c_help(0) &
             + RHO*G*H_t_new(j,i)*(1.0_dp-zeta_t(kt)) &
               *(dzs_dxi_g(j,i)**2+dzs_deta_g(j,i)**2)**0.5_dp

#if (DYNAMICS==2)
   if (.not.flag_shelfy_stream(j,i)) then
#endif
      cw7(kt) = aw7 &
                *enh_t(kt,j,i) &
                *ratefac_t(omega_t(kt,j,i)) &
                *creep(sigma_t_help(kt)) &
                *sigma_t_help(kt)*sigma_t_help(kt)
#if (DYNAMICS==2)
   else
      cw7(kt) = 2.0_dp*aw7 &
                *viscosity(de_t(kt,j,i), &
                           temp_t_m(kt,j,i), temp_t_m(kt,j,i), &
                           omega_t(kt,j,i), &
                           enh_t(kt,j,i), 0) &
                *de_t(kt,j,i)**2
   end if
#endif

   strain_heating_t(kt,j,i) = L*cw7(kt)*dtime_temp_inv

end do

#if (ADV_VERT==1)

kt=0
cw2_sg(kt) = 0.5_dp*(cw2(kt)+cw2(kt+1))
cw3_sg(kt) = 0.5_dp*(cw3(kt)+cw3(kt+1))
cw4_sg(kt) = 0.5_dp*(cw4(kt)+cw4(kt+1))
adv_vert_w_sg(kt) = cw1_sg(kt)-cw2_sg(kt)-cw3_sg(kt)-cw4_sg(kt)
abs_adv_vert_w_sg(kt) = abs(adv_vert_w_sg(kt))   ! only needed for kt=0 ...
kt=KTMAX-1
cw2_sg(kt) = 0.5_dp*(cw2(kt)+cw2(kt+1))
cw3_sg(kt) = 0.5_dp*(cw3(kt)+cw3(kt+1))
cw4_sg(kt) = 0.5_dp*(cw4(kt)+cw4(kt+1))
adv_vert_w_sg(kt) = cw1_sg(kt)-cw2_sg(kt)-cw3_sg(kt)-cw4_sg(kt)
abs_adv_vert_w_sg(kt) = abs(adv_vert_w_sg(kt))   ! ... and kt=KTMAX-1

#elif (ADV_VERT==2 || ADV_VERT==3)

do kt=0, KTMAX-1
   cw2_sg(kt) = 0.5_dp*(cw2(kt)+cw2(kt+1))
   cw3_sg(kt) = 0.5_dp*(cw3(kt)+cw3(kt+1))
   cw4_sg(kt) = 0.5_dp*(cw4(kt)+cw4(kt+1))
   adv_vert_w_sg(kt) = cw1_sg(kt)-cw2_sg(kt)-cw3_sg(kt)-cw4_sg(kt)
   abs_adv_vert_w_sg(kt) = abs(adv_vert_w_sg(kt))
end do

#endif

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
lgs_b(kr)   = temp_t_m(0,j,i)

!-------- Solve system of linear equations --------

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, KRMAX)

!-------- Assign the result --------

do kr=0, KRMAX
   temp_r_new(kr,j,i) = lgs_x(kr)
end do

!-------- Set up the equations for the water content in
!         temperate ice --------

kt=0
lgs_a1(kt) = 1.0_dp
lgs_a2(kt) = -1.0_dp
lgs_b(kt)  = 0.0_dp

do kt=1, KTMAX-1

#if (ADV_VERT==1)

   lgs_a0(kt) = -0.5_dp*(cw1(kt)-cw2(kt)-cw3(kt)-cw4(kt)) - cw5
   lgs_a1(kt) = 1.0_dp + 2.0_dp*cw5
   lgs_a2(kt) = 0.5_dp*(cw1(kt)-cw2(kt)-cw3(kt)-cw4(kt)) - cw5

#elif (ADV_VERT==2)

   lgs_a0(kt) &
         = -0.5_dp*(adv_vert_w_sg(kt-1)+abs_adv_vert_w_sg(kt-1)) &
           -cw5
   lgs_a1(kt) &
         = 1.0_dp &
           +0.5_dp*(adv_vert_w_sg(kt-1)+abs_adv_vert_w_sg(kt-1)) &
           -0.5_dp*(adv_vert_w_sg(kt)  -abs_adv_vert_w_sg(kt)  ) &
           +2.0_dp*cw5
   lgs_a2(kt) &
         =  0.5_dp*(adv_vert_w_sg(kt)  -abs_adv_vert_w_sg(kt)  ) &
           -cw5

#elif (ADV_VERT==3)

   adv_vert_w_help = 0.5_dp*(adv_vert_w_sg(kt)+adv_vert_w_sg(kt-1))

   lgs_a0(kt) &
         = -max(adv_vert_w_help, 0.0_dp) &
           -cw5
   lgs_a1(kt) &
         = 1.0_dp &
           +max(adv_vert_w_help, 0.0_dp)-min(adv_vert_w_help, 0.0_dp) &
           +2.0_dp*cw5
   lgs_a2(kt) &
         =  min(adv_vert_w_help, 0.0_dp) &
           -cw5

#endif

#if (ADV_HOR==2)

   lgs_b(kt) = omega_t(kt,j,i) + cw7(kt) + cw8 + cw9(kt) &
       - ( dtt_2dxi* &
             ( (vx_t(kt,j,i)-abs(vx_t(kt,j,i))) &
               *(omega_t(kt,j,i+1)-omega_t(kt,j,i)) &
               *insq_g11_sgx(j,i) &
              +(vx_t(kt,j,i-1)+abs(vx_t(kt,j,i-1))) &
               *(omega_t(kt,j,i)-omega_t(kt,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_2deta* &
             ( (vy_t(kt,j,i)-abs(vy_t(kt,j,i))) &
               *(omega_t(kt,j+1,i)-omega_t(kt,j,i)) &
               *insq_g22_sgy(j,i) &
              +(vy_t(kt,j-1,i)+abs(vy_t(kt,j-1,i))) &
               *(omega_t(kt,j,i)-omega_t(kt,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#elif (ADV_HOR==3)

   vx_t_help = 0.5_dp*(vx_t(kt,j,i)+vx_t(kt,j,i-1))
   vy_t_help = 0.5_dp*(vy_t(kt,j,i)+vy_t(kt,j-1,i))

   lgs_b(kt) = omega_t(kt,j,i) + cw7(kt) + cw8 + cw9(kt) &
       - ( dtt_dxi* &
             ( min(vx_t_help, 0.0_dp) &
               *(omega_t(kt,j,i+1)-omega_t(kt,j,i)) &
               *insq_g11_sgx(j,i) &
              +max(vx_t_help, 0.0_dp) &
               *(omega_t(kt,j,i)-omega_t(kt,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_deta* &
             ( min(vy_t_help, 0.0_dp) &
               *(omega_t(kt,j+1,i)-omega_t(kt,j,i)) &
               *insq_g22_sgy(j,i) &
              +max(vy_t_help, 0.0_dp) &
               *(omega_t(kt,j,i)-omega_t(kt,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#endif

end do

kt=KTMAX

#if (!defined(CTS_MELTING_FREEZING) || CTS_MELTING_FREEZING==1)

if (am_perp(j,i) >= zero) then   ! melting condition
   lgs_a0(kt) = 0.0_dp
   lgs_a1(kt) = 1.0_dp
   lgs_b(kt)  = 0.0_dp
else   ! am_perp(j,i) < 0.0, freezing condition
   lgs_a0(kt) = -1.0_dp
   lgs_a1(kt) = 1.0_dp
   lgs_b(kt)  = 0.0_dp
end if

#elif (CTS_MELTING_FREEZING==2)

lgs_a0(kt) = 0.0_dp
lgs_a1(kt) = 1.0_dp   ! melting condition assumed
lgs_b(kt)  = 0.0_dp

#else

errormsg = ' >>> calc_temp3: CTS_MELTING_FREEZING must be either 1 or 2!'
call error(errormsg)

#endif

!-------- Solve system of linear equations --------

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, KTMAX)

!-------- Assign the result, compute the water drainage --------

Q_tld(j,i) = 0.0_dp

do kt=0, KTMAX

   if (lgs_x(kt) < zero) then
      omega_t_new(kt,j,i) = 0.0_dp   ! (as a precaution)
   else if (lgs_x(kt) < OMEGA_MAX) then
      omega_t_new(kt,j,i) = lgs_x(kt)
   else
      omega_t_new(kt,j,i) = OMEGA_MAX
      Q_tld(j,i) = Q_tld(j,i) &
                     +aqtld*H_t_new(j,i)*(lgs_x(kt)-OMEGA_MAX)
   end if

end do

!-------- Set up the equations for the ice temperature --------

!  ------ Abbreviation for the jump of the temperature gradient with
!         the new omega

#if (!defined(CTS_MELTING_FREEZING) || CTS_MELTING_FREEZING==1)

if (am_perp(j,i) >= zero) then   ! melting condition
   cm2 = 0.0_dp
else   ! am_perp(j,i) < 0.0, freezing condition
   cm2  = am2*H_c_new(j,i)*omega_t_new(KTMAX,j,i)*am_perp(j,i) &
          /kappa_val(temp_c(0,j,i))
end if

#elif (CTS_MELTING_FREEZING==2)

cm2 = 0.0_dp   ! melting condition assumed

#else

errormsg = ' >>> calc_temp3: CTS_MELTING_FREEZING must be either 1 or 2!'
call error(errormsg)

#endif

kc=0
lgs_a1(kc) = 1.0_dp
lgs_a2(kc) = -1.0_dp
lgs_b(kc)  = -cm1-cm2

do kc=1, KCMAX-1

#if (ADV_VERT==1)

   lgs_a0(kc) = -0.5_dp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ct5(kc)*ct6(kc-1)
   lgs_a1(kc) = 1.0_dp+ct5(kc)*(ct6(kc)+ct6(kc-1))
   lgs_a2(kc) = 0.5_dp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ct5(kc)*ct6(kc)

#elif (ADV_VERT==2)

   lgs_a0(kc) &
         = -0.5_dp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -ct5(kc)*ct6(kc-1)
   lgs_a1(kc) &
         = 1.0_dp &
           +0.5_dp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -0.5_dp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           +ct5(kc)*(ct6(kc)+ct6(kc-1))
   lgs_a2(kc) &
         =  0.5_dp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           -ct5(kc)*ct6(kc)

#elif (ADV_VERT==3)

   adv_vert_help = 0.5_dp*(adv_vert_sg(kc)+adv_vert_sg(kc-1))

   lgs_a0(kc) &
         = -max(adv_vert_help, 0.0_dp) &
           -ct5(kc)*ct6(kc-1)
   lgs_a1(kc) &
         = 1.0_dp &
           +max(adv_vert_help, 0.0_dp)-min(adv_vert_help, 0.0_dp) &
           +ct5(kc)*(ct6(kc)+ct6(kc-1))
   lgs_a2(kc) &
         =  min(adv_vert_help, 0.0_dp) &
           -ct5(kc)*ct6(kc)

#endif

#if (ADV_HOR==2)

   lgs_b(kc) = temp_c(kc,j,i) + ct7(kc) &
       - ( dtt_2dxi* &
             ( (vx_c(kc,j,i)-abs(vx_c(kc,j,i))) &
               *(temp_c(kc,j,i+1)-temp_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +(vx_c(kc,j,i-1)+abs(vx_c(kc,j,i-1))) &
               *(temp_c(kc,j,i)-temp_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_2deta* &
             ( (vy_c(kc,j,i)-abs(vy_c(kc,j,i))) &
               *(temp_c(kc,j+1,i)-temp_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +(vy_c(kc,j-1,i)+abs(vy_c(kc,j-1,i))) &
               *(temp_c(kc,j,i)-temp_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#elif (ADV_HOR==3)

   vx_c_help = 0.5_dp*(vx_c(kc,j,i)+vx_c(kc,j,i-1))
   vy_c_help = 0.5_dp*(vy_c(kc,j,i)+vy_c(kc,j-1,i))

   lgs_b(kc) = temp_c(kc,j,i) + ct7(kc) &
       - ( dtt_dxi* &
             ( min(vx_c_help, 0.0_dp) &
               *(temp_c(kc,j,i+1)-temp_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +max(vx_c_help, 0.0_dp) &
               *(temp_c(kc,j,i)-temp_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_deta* &
             ( min(vy_c_help, 0.0_dp) &
               *(temp_c(kc,j+1,i)-temp_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +max(vy_c_help, 0.0_dp) &
               *(temp_c(kc,j,i)-temp_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#endif

end do

kc=KCMAX
lgs_a0(kc) = 0.0_dp
lgs_a1(kc) = 1.0_dp
lgs_b(kc)  = temp_s(j,i)

!-------- Solve system of linear equations --------

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, KCMAX)

!-------- Assign the result --------

do kc=0, KCMAX
   temp_c_new(kc,j,i) = lgs_x(kc)
end do

!-------- Set up the equations for the age (cold and temperate ice
!         simultaneously) --------

kt=0                                                  ! adv_vert_w_sg(0) <= 0
lgs_a1(kt) = 1.0_dp - min(adv_vert_w_sg(kt), 0.0_dp)  ! (directed downward)
lgs_a2(kt) = min(adv_vert_w_sg(kt), 0.0_dp)           ! assumed/enforced

#if (ADV_HOR==2)

lgs_b(kt) = age_t(kt,j,i) + dtime_temp &
       - ( dtt_2dxi* &
             ( (vx_t(kt,j,i)-abs(vx_t(kt,j,i))) &
               *(age_t(kt,j,i+1)-age_t(kt,j,i)) &
               *insq_g11_sgx(j,i) &
              +(vx_t(kt,j,i-1)+abs(vx_t(kt,j,i-1))) &
               *(age_t(kt,j,i)-age_t(kt,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_2deta* &
             ( (vy_t(kt,j,i)-abs(vy_t(kt,j,i))) &
               *(age_t(kt,j+1,i)-age_t(kt,j,i)) &
               *insq_g22_sgy(j,i) &
              +(vy_t(kt,j-1,i)+abs(vy_t(kt,j-1,i))) &
               *(age_t(kt,j,i)-age_t(kt,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#elif (ADV_HOR==3)

vx_t_help = 0.5_dp*(vx_t(kt,j,i)+vx_t(kt,j,i-1))
vy_t_help = 0.5_dp*(vy_t(kt,j,i)+vy_t(kt,j-1,i))

lgs_b(kt) = age_t(kt,j,i) + dtime_temp &
       - ( dtt_dxi* &
             ( min(vx_t_help, 0.0_dp) &
               *(age_t(kt,j,i+1)-age_t(kt,j,i)) &
               *insq_g11_sgx(j,i) &
              +max(vx_t_help, 0.0_dp) &
               *(age_t(kt,j,i)-age_t(kt,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_deta* &
             ( min(vy_t_help, 0.0_dp) &
               *(age_t(kt,j+1,i)-age_t(kt,j,i)) &
               *insq_g22_sgy(j,i) &
              +max(vy_t_help, 0.0_dp) &
               *(age_t(kt,j,i)-age_t(kt,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#endif

do kt=1, KTMAX-1

#if (ADV_VERT==1)

   lgs_a0(kt) = -0.5_dp*(cw1(kt)-cw2(kt)-cw3(kt)-cw4(kt)) - ci3
   lgs_a1(kt) = 1.0_dp + 2.0_dp*ci3
   lgs_a2(kt) = 0.5_dp*(cw1(kt)-cw2(kt)-cw3(kt)-cw4(kt)) - ci3

#elif (ADV_VERT==2)

   lgs_a0(kt) = -0.5_dp*(adv_vert_w_sg(kt-1)+abs_adv_vert_w_sg(kt-1))
   lgs_a1(kt) = 1.0_dp &
               +0.5_dp*(adv_vert_w_sg(kt-1)+abs_adv_vert_w_sg(kt-1)) &
               -0.5_dp*(adv_vert_w_sg(kt)  -abs_adv_vert_w_sg(kt)  )
   lgs_a2(kt) =  0.5_dp*(adv_vert_w_sg(kt)  -abs_adv_vert_w_sg(kt)  )

#elif (ADV_VERT==3)

   adv_vert_w_help = 0.5_dp*(adv_vert_w_sg(kt)+adv_vert_w_sg(kt-1))

   lgs_a0(kt) = -max(adv_vert_w_help, 0.0_dp)
   lgs_a1(kt) = 1.0_dp &
               +max(adv_vert_w_help, 0.0_dp)-min(adv_vert_w_help, 0.0_dp)
   lgs_a2(kt) =  min(adv_vert_w_help, 0.0_dp)

#endif

#if (ADV_HOR==2)

   lgs_b(kt) = age_t(kt,j,i) + dtime_temp &
       - ( dtt_2dxi* &
             ( (vx_t(kt,j,i)-abs(vx_t(kt,j,i))) &
               *(age_t(kt,j,i+1)-age_t(kt,j,i)) &
               *insq_g11_sgx(j,i) &
              +(vx_t(kt,j,i-1)+abs(vx_t(kt,j,i-1))) &
               *(age_t(kt,j,i)-age_t(kt,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_2deta* &
             ( (vy_t(kt,j,i)-abs(vy_t(kt,j,i))) &
               *(age_t(kt,j+1,i)-age_t(kt,j,i)) &
               *insq_g22_sgy(j,i) &
              +(vy_t(kt,j-1,i)+abs(vy_t(kt,j-1,i))) &
               *(age_t(kt,j,i)-age_t(kt,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#elif (ADV_HOR==3)

   vx_t_help = 0.5_dp*(vx_t(kt,j,i)+vx_t(kt,j,i-1))
   vy_t_help = 0.5_dp*(vy_t(kt,j,i)+vy_t(kt,j-1,i))

   lgs_b(kt) = age_t(kt,j,i) + dtime_temp &
       - ( dtt_dxi* &
             ( min(vx_t_help, 0.0_dp) &
               *(age_t(kt,j,i+1)-age_t(kt,j,i)) &
               *insq_g11_sgx(j,i) &
              +max(vx_t_help, 0.0_dp) &
               *(age_t(kt,j,i)-age_t(kt,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_deta* &
             ( min(vy_t_help, 0.0_dp) &
               *(age_t(kt,j+1,i)-age_t(kt,j,i)) &
               *insq_g22_sgy(j,i) &
              +max(vy_t_help, 0.0_dp) &
               *(age_t(kt,j,i)-age_t(kt,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#endif

end do

#if (ADV_VERT==1)

kt=KTMAX
kc=0

lgs_a0(kt) = -0.5_dp*(cw1(kt)-cw2(kt)-cw3(kt)-cw4(kt)) - ci3
lgs_a1(kt) = 1.0_dp + 2.0_dp*ci3
lgs_a2(kt) = 0.5_dp*(cw1(kt)-cw2(kt)-cw3(kt)-cw4(kt)) - ci3

#if (ADV_HOR==2)

lgs_b(kt) = age_t(kt,j,i) + dtime_temp &
       - ( dtt_2dxi* &
             ( (vx_t(kt,j,i)-abs(vx_t(kt,j,i))) &
               *(age_t(kt,j,i+1)-age_t(kt,j,i)) &
               *insq_g11_sgx(j,i) &
              +(vx_t(kt,j,i-1)+abs(vx_t(kt,j,i-1))) &
               *(age_t(kt,j,i)-age_t(kt,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_2deta* &
             ( (vy_t(kt,j,i)-abs(vy_t(kt,j,i))) &
               *(age_t(kt,j+1,i)-age_t(kt,j,i)) &
               *insq_g22_sgy(j,i) &
              +(vy_t(kt,j-1,i)+abs(vy_t(kt,j-1,i))) &
               *(age_t(kt,j,i)-age_t(kt,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#elif (ADV_HOR==3)

vx_t_help = 0.5_dp*(vx_t(kt,j,i)+vx_t(kt,j,i-1))
vy_t_help = 0.5_dp*(vy_t(kt,j,i)+vy_t(kt,j-1,i))

lgs_b(kt) = age_t(kt,j,i) + dtime_temp &
       - ( dtt_dxi* &
             ( min(vx_t_help, 0.0_dp) &
               *(age_t(kt,j,i+1)-age_t(kt,j,i)) &
               *insq_g11_sgx(j,i) &
              +max(vx_t_help, 0.0_dp) &
               *(age_t(kt,j,i)-age_t(kt,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_deta* &
             ( min(vy_t_help, 0.0_dp) &
               *(age_t(kt,j+1,i)-age_t(kt,j,i)) &
               *insq_g22_sgy(j,i) &
              +max(vy_t_help, 0.0_dp) &
               *(age_t(kt,j,i)-age_t(kt,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#endif

#elif (ADV_VERT==2 || ADV_VERT==3)

kt=KTMAX
kc=0

if (adv_vert_sg(kc) <= zero) then

   lgs_a0(KTMAX+kc) = 0.0_dp
   lgs_a1(KTMAX+kc) = 1.0_dp - adv_vert_sg(kc)
   lgs_a2(KTMAX+kc) = adv_vert_sg(kc)

#if (ADV_HOR==2)

   lgs_b(KTMAX+kc) = age_c(kc,j,i) + dtime_temp &
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

   lgs_b(KTMAX+kc) = age_c(kc,j,i) + dtime_temp &
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

else if (adv_vert_w_sg(kt-1) >= zero) then

   lgs_a0(kt) = -adv_vert_w_sg(kt-1)
   lgs_a1(kt) = 1.0_dp + adv_vert_w_sg(kt-1)
   lgs_a2(kt) = 0.0_dp

#if (ADV_HOR==2)

   lgs_b(kt) = age_t(kt,j,i) + dtime_temp &
       - ( dtt_2dxi* &
             ( (vx_t(kt,j,i)-abs(vx_t(kt,j,i))) &
               *(age_t(kt,j,i+1)-age_t(kt,j,i)) &
               *insq_g11_sgx(j,i) &
              +(vx_t(kt,j,i-1)+abs(vx_t(kt,j,i-1))) &
               *(age_t(kt,j,i)-age_t(kt,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_2deta* &
             ( (vy_t(kt,j,i)-abs(vy_t(kt,j,i))) &
               *(age_t(kt,j+1,i)-age_t(kt,j,i)) &
               *insq_g22_sgy(j,i) &
              +(vy_t(kt,j-1,i)+abs(vy_t(kt,j-1,i))) &
               *(age_t(kt,j,i)-age_t(kt,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#elif (ADV_HOR==3)

   vx_t_help = 0.5_dp*(vx_t(kt,j,i)+vx_t(kt,j,i-1))
   vy_t_help = 0.5_dp*(vy_t(kt,j,i)+vy_t(kt,j-1,i))

   lgs_b(kt) = age_t(kt,j,i) + dtime_temp &
       - ( dtt_dxi* &
             ( min(vx_t_help, 0.0_dp) &
               *(age_t(kt,j,i+1)-age_t(kt,j,i)) &
               *insq_g11_sgx(j,i) &
              +max(vx_t_help, 0.0_dp) &
               *(age_t(kt,j,i)-age_t(kt,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_deta* &
             ( min(vy_t_help, 0.0_dp) &
               *(age_t(kt,j+1,i)-age_t(kt,j,i)) &
               *insq_g22_sgy(j,i) &
              +max(vy_t_help, 0.0_dp) &
               *(age_t(kt,j,i)-age_t(kt,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#endif

else

   lgs_a0(kt) = -0.5_dp
   lgs_a1(kt) = 1.0_dp
   lgs_a2(kt) = -0.5_dp
   lgs_b(kt)  = 0.0_dp
   ! Makeshift: Average of age_c(kc=1) and age_t(kt=KTMAX-1)

end if

#endif

do kc=1, KCMAX-1

#if (ADV_VERT==1)

   lgs_a0(KTMAX+kc) = -0.5_dp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ci1(kc)*ci2(kc-1)
   lgs_a1(KTMAX+kc) = 1.0_dp+ci1(kc)*(ci2(kc)+ci2(kc-1))
   lgs_a2(KTMAX+kc) = 0.5_dp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ci1(kc)*ci2(kc)

#elif (ADV_VERT==2)

   lgs_a0(KTMAX+kc) = -0.5_dp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1))
   lgs_a1(KTMAX+kc) =  1.0_dp &
                      +0.5_dp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
                      -0.5_dp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  )
   lgs_a2(KTMAX+kc) =  0.5_dp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  )

#elif (ADV_VERT==3)

   adv_vert_help = 0.5_dp*(adv_vert_sg(kc)+adv_vert_sg(kc-1))

   lgs_a0(KTMAX+kc) = -max(adv_vert_help, 0.0_dp)
   lgs_a1(KTMAX+kc) =  1.0_dp &
                      +max(adv_vert_help, 0.0_dp)-min(adv_vert_help, 0.0_dp)
   lgs_a2(KTMAX+kc) =  min(adv_vert_help, 0.0_dp)

#endif

#if (ADV_HOR==2)

   lgs_b(KTMAX+kc) = age_c(kc,j,i) + dtime_temp &
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

   lgs_b(KTMAX+kc) = age_c(kc,j,i) + dtime_temp &
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
   lgs_a0(KTMAX+kc) = 0.0_dp
   lgs_a1(KTMAX+kc) = 1.0_dp
   lgs_b(KTMAX+kc)  = 0.0_dp
else
   lgs_a0(KTMAX+kc) = -max(adv_vert_sg(kc-1), 0.0_dp)
   lgs_a1(KTMAX+kc) = 1.0_dp + max(adv_vert_sg(kc-1), 0.0_dp)
                             ! adv_vert_sg(KCMAX-1) >= 0 (directed upward)
                             ! assumed/enforced
#if (ADV_HOR==2)

   lgs_b(KTMAX+kc) = age_c(kc,j,i) + dtime_temp &
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

   lgs_b(KTMAX+kc) = age_c(kc,j,i) + dtime_temp &
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

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, KCMAX+KTMAX)

!-------- Assign the result,
!         restriction to interval [0, AGE_MAX yr] --------

do kt=0, KTMAX

   age_t_new(kt,j,i) = lgs_x(kt)

   if (age_t_new(kt,j,i) < (AGE_MIN*year2sec)) &
                           age_t_new(kt,j,i) = 0.0_dp
   if (age_t_new(kt,j,i) > (AGE_MAX*year2sec)) &
                           age_t_new(kt,j,i) = AGE_MAX*year2sec

end do

do kc=0, KCMAX

   age_c_new(kc,j,i) = lgs_x(KTMAX+kc)

   if (age_c_new(kc,j,i) < (AGE_MIN*year2sec)) &
                           age_c_new(kc,j,i) = 0.0_dp
   if (age_c_new(kc,j,i) > (AGE_MAX*year2sec)) &
                           age_c_new(kc,j,i) = AGE_MAX*year2sec

end do

end subroutine calc_temp3

!-------------------------------------------------------------------------------
!> Computation of temperature and age for an ice-free column.
!-------------------------------------------------------------------------------
subroutine calc_temp_r(atr1, alb1, i, j)

#if !defined(ALLOW_TAPENADE) /* Normal */
use sico_maths_m, only : tri_sle
#else /* Tapenade */
use sico_maths_m
#endif /* Normal vs. Tapenade */

implicit none

integer(i4b), intent(in)    :: i, j
real(dp), intent(in) :: atr1, alb1

integer(i4b) :: kc, kt, kr
real(dp) :: ctr1, clb1
real(dp) :: lgs_a0(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_a1(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_a2(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_x(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_b(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX)

!-------- Abbreviations --------

ctr1 = atr1
clb1 = alb1*q_geo(j,i)

!-------- Set up the equations for the bedrock temperature --------

kr=0
#if defined(ALLOW_TAPENADE) /* Tapenade */
lgs_a0(kr) = 0.0_dp
#endif /* Tapenade */
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

!-------- Water content and age
!                       in the non-existing lower (kt) ice layer --------

do kt=0, KTMAX
   omega_t_new(kt,j,i) = 0.0_dp
   age_t_new(kt,j,i)   = 0.0_dp
end do

!-------- Temperature and age
!                     in the non-existing upper (kc) ice layer --------

do kc=0, KCMAX
   temp_c_new(kc,j,i)  = temp_s(j,i)
   age_c_new(kc,j,i)   = 0.0_dp
end do

end subroutine calc_temp_r

!-------------------------------------------------------------------------------
!> Upward shifting of the CTS.
!-------------------------------------------------------------------------------
subroutine shift_cts_upward(at1, at2_1, at2_2, at3_1, at3_2, &
              at4_1, at4_2, at5, at6, at7, atr1, am1, am2, alb1, &
              aw1, aw2, aw3, aw4, aw5, aw7, aw8, aw9, aqtld, &
              ai1, ai2, ai3, dzeta_t, &
              dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
              i, j)

implicit none

integer(i4b), intent(in)    :: i, j
real(dp),     intent(in) :: at1(0:KCMAX), at2_1(0:KCMAX), at2_2(0:KCMAX), &
                            at3_1(0:KCMAX), at3_2(0:KCMAX), at4_1(0:KCMAX), &
                            at4_2(0:KCMAX), at5(0:KCMAX), at6(0:KCMAX), at7, &
                            ai1(0:KCMAX), ai2(0:KCMAX), ai3, &
                            atr1, am1, am2, alb1
real(dp),     intent(in) :: aw1, aw2, aw3, aw4, aw5, aw7, aw8, aw9, aqtld
real(dp),     intent(in) :: dzeta_t
real(dp),     intent(in) :: dtime_temp, dtime_temp_inv, dtt_2dxi, dtt_2deta

real(dp) :: zm_shift
real(dp) :: difftemp_a, difftemp_b, interpol

zm_shift = 1.0_dp   ! CTS shift in intervals of 1 m

!-------- Temperature discrepancy from the computation of the main
!         program --------

difftemp_a = temp_c_new(0,j,i)-(-BETA*H_c_new(j,i))
if (difftemp_a <= 0.0_dp) return

!-------- Shift CTS upward until it is too high --------

do while (difftemp_a > 0.0_dp)

   zm_new(j,i)  = zm_new(j,i)  + zm_shift
   if (zm_new(j,i) >= zs(j,i)) then
      zm_new(j,i)  = zm_new(j,i) - zm_shift
      return
   end if
   H_c_new(j,i) = H_c_new(j,i) - zm_shift
   H_t_new(j,i) = H_t_new(j,i) + zm_shift

   dH_t_dtau(j,i) = (zm_new(j,i)-zm(j,i))*dtime_temp_inv
   dzm_dtau(j,i)  = dzb_dtau(j,i)+dH_t_dtau(j,i)
   dH_c_dtau(j,i) = dzs_dtau(j,i)-dzm_dtau(j,i)

   am_perp(j,i) = am_perp_st(j,i) + dzm_dtau(j,i)

   call calc_temp3(at1, at2_1, at2_2, at3_1, at3_2, &
                 at4_1, at4_2, at5, at6, at7, atr1, &
                 am1, am2, alb1, &
                 aw1, aw2, aw3, aw4, aw5, aw7, aw8, aw9, aqtld, &
                 ai1, ai2, ai3, dzeta_t, &
                 dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
                 i, j)

   difftemp_b = difftemp_a
   difftemp_a = temp_c_new(0,j,i)-(-BETA*H_c_new(j,i))

end do

!-------- Interpolate the CTS position from the last (_a) and the
!         last but one (_b) value, weighed with the temperature
!         discrepancies at the CTS --------

interpol = difftemp_a/(difftemp_b-difftemp_a)*zm_shift

zm_new(j,i)  = zm_new(j,i)  + interpol
H_c_new(j,i) = H_c_new(j,i) - interpol
H_t_new(j,i) = H_t_new(j,i) + interpol

dH_t_dtau(j,i) = (zm_new(j,i)-zm(j,i))*dtime_temp_inv
dzm_dtau(j,i)  = dzb_dtau(j,i)+dH_t_dtau(j,i)
dH_c_dtau(j,i) = dzs_dtau(j,i)-dzm_dtau(j,i)

am_perp(j,i) = am_perp_st(j,i) + dzm_dtau(j,i)

call calc_temp3(at1, at2_1, at2_2, at3_1, at3_2, &
                at4_1, at4_2, at5, at6, at7, atr1, &
                am1, am2, alb1, &
                aw1, aw2, aw3, aw4, aw5, aw7, aw8, aw9, aqtld, &
                ai1, ai2, ai3, dzeta_t, &
                dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
                i, j)

end subroutine shift_cts_upward

!-------------------------------------------------------------------------------
!> Downward shifting of the CTS.
!-------------------------------------------------------------------------------
subroutine shift_cts_downward(at1, at2_1, at2_2, at3_1, at3_2, &
              at4_1, at4_2, at5, at6, at7, atr1, am1, am2, alb1, &
              aw1, aw2, aw3, aw4, aw5, aw7, aw8, aw9, aqtld, &
              ai1, ai2, ai3, dzeta_t, &
              dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
              i, j)

implicit none

integer(i4b), intent(in)    :: i, j
real(dp),     intent(in) :: at1(0:KCMAX), at2_1(0:KCMAX), at2_2(0:KCMAX), &
                            at3_1(0:KCMAX), at3_2(0:KCMAX), at4_1(0:KCMAX), &
                            at4_2(0:KCMAX), at5(0:KCMAX), at6(0:KCMAX), at7, &
                            ai1(0:KCMAX), ai2(0:KCMAX), ai3, &
                            atr1, am1, am2, alb1
real(dp),     intent(in) :: aw1, aw2, aw3, aw4, aw5, aw7, aw8, aw9, aqtld
real(dp),     intent(in) :: dzeta_t
real(dp),     intent(in) :: dtt_2dxi, dtt_2deta, dtime_temp, dtime_temp_inv

real(dp) :: zm_shift
real(dp) :: difftemp_a, difftemp_b, interpol

zm_shift = 1.0_dp   ! CTS shift in intervals of 1 m

!-------- Temperature discrepancy from the computation of the main
!         program --------

difftemp_a = temp_c_new(0,j,i)-(-BETA*H_c_new(j,i))
if (difftemp_a >= 0.0_dp) return

!-------- Shift CTS downward until it is too low --------

do while (difftemp_a < 0.0_dp)

   zm_new(j,i)  = zm_new(j,i) - zm_shift

!  ------ Special case: CTS too close to the base

   if (zm_new(j,i) <= zb(j,i)) then

      zm_shift = (zm_new(j,i)+zm_shift)-(zb(j,i)+0.001_dp)
      zm_new(j,i)  = zb(j,i)+0.001_dp
      H_c_new(j,i) = H_c_new(j,i)+H_t_new(j,i)-0.001_dp
      H_t_new(j,i) = 0.001_dp
!                   ! CTS positioned 1 mm above ice base --------

      dH_t_dtau(j,i) = (zm_new(j,i)-zm(j,i))*dtime_temp_inv
      dzm_dtau(j,i)  = dzb_dtau(j,i)+dH_t_dtau(j,i)
      dH_c_dtau(j,i) = dzs_dtau(j,i)-dzm_dtau(j,i)

      am_perp(j,i) = am_perp_st(j,i) + dzm_dtau(j,i)

      call calc_temp3(at1, at2_1, at2_2, at3_1, at3_2, &
                 at4_1, at4_2, at5, at6, at7, atr1, &
                 am1, am2, alb1, &
                 aw1, aw2, aw3, aw4, aw5, aw7, aw8, aw9, aqtld, &
                 ai1, ai2, ai3, dzeta_t, &
                 dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
                 i, j)

      difftemp_b = difftemp_a
      difftemp_a = temp_c_new(0,j,i)-(-BETA*H_c_new(j,i))

      if (difftemp_a >= 0.0_dp) then ! CTS remains above the base

!    ---- Interpolate the CTS position from the last (_a) and the
!         last but one (_b) value, weighed with the temperature
!         discrepancies at the CTS --------

         interpol = difftemp_a/(difftemp_a-difftemp_b)*zm_shift

         zm_new(j,i)  = zm_new(j,i)  + interpol
         H_c_new(j,i) = H_c_new(j,i) - interpol
         H_t_new(j,i) = H_t_new(j,i) + interpol

         dH_t_dtau(j,i) = (zm_new(j,i)-zm(j,i))*dtime_temp_inv
         dzm_dtau(j,i)  = dzb_dtau(j,i)+dH_t_dtau(j,i)
         dH_c_dtau(j,i) = dzs_dtau(j,i)-dzm_dtau(j,i)

         am_perp(j,i) = am_perp_st(j,i) + dzm_dtau(j,i)

         call calc_temp3(at1, at2_1, at2_2, at3_1, at3_2, &
                 at4_1, at4_2, at5, at6, at7, atr1, &
                 am1, am2, alb1, &
                 aw1, aw2, aw3, aw4, aw5, aw7, aw8, aw9, aqtld, &
                 ai1, ai2, ai3, dzeta_t, &
                 dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
                 i, j)

      else   ! CTS disappears

         n_cts_new(j,i) = 0
         zm_new(j,i) = zb(j,i)
         H_c_new(j,i) = H_c_new(j,i)+H_t_new(j,i)
         H_t_new(j,i) = 0.0_dp

         dH_t_dtau(j,i) = (zm_new(j,i)-zm(j,i))*dtime_temp_inv
         dzm_dtau(j,i)  = dzb_dtau(j,i)+dH_t_dtau(j,i)
         dH_c_dtau(j,i) = dzs_dtau(j,i)-dzm_dtau(j,i)

         am_perp(j,i) = am_perp_st(j,i) + dzm_dtau(j,i)

         call calc_temp2(at1, at2_1, at2_2, at3_1, at3_2, &
                 at4_1, at4_2, at5, at6, at7, atr1, alb1, &
                 ai1, ai2, &
                 dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
                 i, j)

      end if

      return
   end if

!  ------ End of treatment of special case

   H_c_new(j,i) = H_c_new(j,i) + zm_shift
   H_t_new(j,i) = H_t_new(j,i) - zm_shift

   dH_t_dtau(j,i) = (zm_new(j,i)-zm(j,i))*dtime_temp_inv
   dzm_dtau(j,i)  = dzb_dtau(j,i)+dH_t_dtau(j,i)
   dH_c_dtau(j,i) = dzs_dtau(j,i)-dzm_dtau(j,i)

   am_perp(j,i) = am_perp_st(j,i) + dzm_dtau(j,i)

   call calc_temp3(at1, at2_1, at2_2, at3_1, at3_2, &
                 at4_1, at4_2, at5, at6, at7, atr1, &
                 am1, am2, alb1, &
                 aw1, aw2, aw3, aw4, aw5, aw7, aw8, aw9, aqtld, &
                 ai1, ai2, ai3, dzeta_t, &
                 dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
                 i, j)

   difftemp_b = difftemp_a
   difftemp_a = temp_c_new(0,j,i)-(-BETA*H_c_new(j,i))

end do

!-------- Interpolate the CTS position from the last (_a) and the
!         last but one (_b) value, weighed with the temperature
!         discrepancies at the CTS --------

interpol = difftemp_a/(difftemp_a-difftemp_b)*zm_shift

zm_new(j,i)  = zm_new(j,i)  + interpol
H_c_new(j,i) = H_c_new(j,i) - interpol
H_t_new(j,i) = H_t_new(j,i) + interpol

dH_t_dtau(j,i) = (zm_new(j,i)-zm(j,i))*dtime_temp_inv
dzm_dtau(j,i)  = dzb_dtau(j,i)+dH_t_dtau(j,i)
dH_c_dtau(j,i) = dzs_dtau(j,i)-dzm_dtau(j,i)

am_perp(j,i) = am_perp_st(j,i) + dzm_dtau(j,i)

call calc_temp3(at1, at2_1, at2_2, at3_1, at3_2, &
                at4_1, at4_2, at5, at6, at7, atr1, &
                am1, am2, alb1, &
                aw1, aw2, aw3, aw4, aw5, aw7, aw8, aw9, aqtld, &
                ai1, ai2, ai3, dzeta_t, &
                dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
                i, j)

end subroutine shift_cts_downward

!-------------------------------------------------------------------------------
!> Computation of temperature and age for ice shelves (floating ice).
!-------------------------------------------------------------------------------
subroutine calc_temp_ssa(at1, at2_1, at2_2, at3_1, at3_2, &
   at4_1, at4_2, at5, at6, at7, atr1, alb1, &
   ai1, ai2, &
   dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv, &
   i, j)

use ice_material_properties_m, only : kappa_val, c_val, viscosity

#if !defined(ALLOW_TAPENADE) /* Normal */
use sico_maths_m, only : tri_sle
#else /* Tapenade */
use sico_maths_m
#endif /* Normal vs. Tapenade */

implicit none

integer(i4b), intent(in)    :: i, j
real(dp), intent(in) :: at1(0:KCMAX), at2_1(0:KCMAX), at2_2(0:KCMAX), &
                        at3_1(0:KCMAX), at3_2(0:KCMAX), at4_1(0:KCMAX), &
                        at4_2(0:KCMAX), at5(0:KCMAX), at6(0:KCMAX), at7, &
                        ai1(0:KCMAX), ai2(0:KCMAX), &
                        atr1, alb1
real(dp), intent(in) :: dtime_temp, dtt_2dxi, dtt_2deta, dtime_temp_inv

integer(i4b) :: kc, kt, kr
real(dp) :: ct1(0:KCMAX), ct2(0:KCMAX), ct3(0:KCMAX), ct4(0:KCMAX), &
            ct5(0:KCMAX), ct6(0:KCMAX), ct7(0:KCMAX), ctr1, clb1
real(dp) :: ct1_sg(0:KCMAX), ct2_sg(0:KCMAX), ct3_sg(0:KCMAX), &
            ct4_sg(0:KCMAX), adv_vert_sg(0:KCMAX), abs_adv_vert_sg(0:KCMAX)
real(dp) :: ci1(0:KCMAX), ci2(0:KCMAX)
real(dp) :: temp_c_help(0:KCMAX)
real(dp) :: vx_c_help, vy_c_help
real(dp) :: adv_vert_help
real(dp) :: dtt_dxi, dtt_deta
real(dp) :: c_val_aux
real(dp) :: lgs_a0(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_a1(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_a2(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_x(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX), &
            lgs_b(0:KCMAX+KTMAX+KRMAX+IMAX+JMAX)
real(dp), parameter :: zero=0.0_dp

!-------- Check for boundary points --------

if ((i == 0).or.(i == IMAX).or.(j == 0).or.(j == JMAX)) then
   errormsg = ' >>> calc_temp_ssa: Boundary points not allowed!'
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

   c_val_aux = c_val(temp_c(kc,j,i))

   ct2(kc) = ( at2_1(kc)*dzm_dtau(j,i) &
           +at2_2(kc)*dH_c_dtau(j,i) )/H_c(j,i)
   ct3(kc) = ( at3_1(kc)*dzm_dxi_g(j,i) &
           +at3_2(kc)*dH_c_dxi_g(j,i) )/H_c(j,i) &
          *0.5_dp*(vx_c(kc,j,i)+vx_c(kc,j,i-1)) *insq_g11_g(j,i)
   ct4(kc) = ( at4_1(kc)*dzm_deta_g(j,i) &
            +at4_2(kc)*dH_c_deta_g(j,i) )/H_c(j,i) &
          *0.5_dp*(vy_c(kc,j,i)+vy_c(kc,j-1,i)) *insq_g22_g(j,i)
   ct5(kc) = at5(kc) &
             /c_val_aux &
             /H_c(j,i)
   ct7(kc) = 2.0_dp*at7 &
             /c_val_aux &
             *viscosity(de_ssa(j,i), &
                        temp_c(kc,j,i), temp_c_m(kc,j,i), 0.0_dp, &
                        enh_c(kc,j,i), 0) &
             *de_ssa(j,i)**2

   strain_heating_c(kc,j,i) = c_val_aux*ct7(kc)*dtime_temp_inv

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
   ct6(kc) = at6(kc) &
    *kappa_val(temp_c_help(kc)) &
    /H_c(j,i)
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
lgs_b(kc)  = temp_c_m(0,j,i)-DELTA_TM_SW

do kc=1, KCMAX-1

#if (ADV_VERT==1)

   lgs_a0(kc) = -0.5_dp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ct5(kc)*ct6(kc-1)
   lgs_a1(kc) = 1.0_dp+ct5(kc)*(ct6(kc)+ct6(kc-1))
   lgs_a2(kc) = 0.5_dp*(ct1(kc)-ct2(kc)-ct3(kc)-ct4(kc)) &
                      -ct5(kc)*ct6(kc)

#elif (ADV_VERT==2)

   lgs_a0(kc) &
         = -0.5_dp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -ct5(kc)*ct6(kc-1)
   lgs_a1(kc) &
         = 1.0_dp &
           +0.5_dp*(adv_vert_sg(kc-1)+abs_adv_vert_sg(kc-1)) &
           -0.5_dp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           +ct5(kc)*(ct6(kc)+ct6(kc-1))
   lgs_a2(kc) &
         =  0.5_dp*(adv_vert_sg(kc)  -abs_adv_vert_sg(kc)  ) &
           -ct5(kc)*ct6(kc)

#elif (ADV_VERT==3)

   adv_vert_help = 0.5_dp*(adv_vert_sg(kc)+adv_vert_sg(kc-1))

   lgs_a0(kc) &
         = -max(adv_vert_help, 0.0_dp) &
           -ct5(kc)*ct6(kc-1)
   lgs_a1(kc) &
         = 1.0_dp &
           +max(adv_vert_help, 0.0_dp)-min(adv_vert_help, 0.0_dp) &
           +ct5(kc)*(ct6(kc)+ct6(kc-1))
   lgs_a2(kc) &
         =  min(adv_vert_help, 0.0_dp) &
           -ct5(kc)*ct6(kc)

#endif

#if (ADV_HOR==2)

   lgs_b(kc) = temp_c(kc,j,i) + ct7(kc) &
       - ( dtt_2dxi* &
             ( (vx_c(kc,j,i)-abs(vx_c(kc,j,i))) &
               *(temp_c(kc,j,i+1)-temp_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +(vx_c(kc,j,i-1)+abs(vx_c(kc,j,i-1))) &
               *(temp_c(kc,j,i)-temp_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_2deta* &
             ( (vy_c(kc,j,i)-abs(vy_c(kc,j,i))) &
               *(temp_c(kc,j+1,i)-temp_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +(vy_c(kc,j-1,i)+abs(vy_c(kc,j-1,i))) &
               *(temp_c(kc,j,i)-temp_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#elif (ADV_HOR==3)

   vx_c_help = 0.5_dp*(vx_c(kc,j,i)+vx_c(kc,j,i-1))
   vy_c_help = 0.5_dp*(vy_c(kc,j,i)+vy_c(kc,j-1,i))

   lgs_b(kc) = temp_c(kc,j,i) + ct7(kc) &
       - ( dtt_dxi* &
             ( min(vx_c_help, 0.0_dp) &
               *(temp_c(kc,j,i+1)-temp_c(kc,j,i)) &
               *insq_g11_sgx(j,i) &
              +max(vx_c_help, 0.0_dp) &
               *(temp_c(kc,j,i)-temp_c(kc,j,i-1)) &
               *insq_g11_sgx(j,i-1) ) &
          +dtt_deta* &
             ( min(vy_c_help, 0.0_dp) &
               *(temp_c(kc,j+1,i)-temp_c(kc,j,i)) &
               *insq_g22_sgy(j,i) &
              +max(vy_c_help, 0.0_dp) &
               *(temp_c(kc,j,i)-temp_c(kc,j-1,i)) &
               *insq_g22_sgy(j-1,i) ) )

#endif

end do

kc=KCMAX
lgs_a0(kc) = 0.0_dp
lgs_a1(kc) = 1.0_dp
lgs_b(kc)  = temp_s(j,i)

!-------- Solve system of linear equations --------

call tri_sle(lgs_a0, lgs_a1, lgs_a2, lgs_x, lgs_b, KCMAX)

!-------- Assign the result --------

do kc=0, KCMAX
   temp_c_new(kc,j,i) = lgs_x(kc)
end do

!-------- Set water content in the non-existing temperate layer
!         to zero --------

do kt=0, KTMAX
   omega_t_new(kt,j,i) = 0.0_dp
end do

!-------- Water drainage from the non-existing temperate layer --------

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

!-------- Age of the ice in the non-existing temperate layer --------

do kt=0, KTMAX
   age_t_new(kt,j,i) = age_c_new(0,j,i)
end do

end subroutine calc_temp_ssa

!-------------------------------------------------------------------------------

end module calc_temp_m
!
