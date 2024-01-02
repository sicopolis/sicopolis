!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  d i s c h a r g e _ w o r k e r s _ m
!
!> @file
!!
!! Ice discharge parameterization for the Greenland ice sheet.
!!
!! References:
!! @li Calov, R., A. Robinson, M. Perrette and A. Ganopolski. 2015.\n
!!     Simulating the Greenland ice sheet under present-day and 
!!     palaeo constraints including a new discharge parameterization.\n
!!     Cryosphere 9, 179-196.
!!
!! @section Copyright
!!
!! Copyright 2011-2024 Reinhard Calov, Andrey Ganopolski
!!                     (Potsdam Institute for Climate Impact Research)
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
!> Ice discharge parameterization for the Greenland ice sheet.
!<------------------------------------------------------------------------------
module discharge_workers_m

  use sico_types_m
  use sico_variables_m
  use sico_vars_m
  use error_m
#if defined(ALLOW_TAPENADE)
  use globals
#endif
 
#if !defined(ALLOW_TAPENADE) /* Normal */
  use compare_float_m
#endif /* Normal */

  implicit none

#if !defined(ALLOW_TAPENADE) /* Normal */
  integer(i4b), private :: disc
  integer(i4b), private :: n_discharge_call, iter_mar_coa
  real(dp),     private :: c_dis_0, s_dis, c_dis_fac
  real(dp),     private :: T_sub_PD, alpha_sub, alpha_o, m_H, m_D, r_mar_eff
  real(dp),     private :: T_sea_freeze

  real(dp),     public  :: dT_glann, dT_sub

  integer(i4b), dimension(0:JMAX,0:IMAX), public  :: mask_mar

  real(dp),     dimension(0:JMAX,0:IMAX), private :: c_dis

  real(dp),     dimension(0:JMAX,0:IMAX), public  :: cst_dist, cos_grad_tc
  real(dp),     dimension(0:JMAX,0:IMAX), public  :: dis_perp

  private
#endif /* Normal */

  public :: disc_param, disc_fields, calc_c_dis_0, discharge

contains


!-------------------------------------------------------------------------------
!> Ice discharge parameters (Greenland).
!! [Assign ice discharge parameters.]
!<------------------------------------------------------------------------------
  subroutine disc_param(dtime)

  ! Author: Reinhard Calov
  ! Institution: Potsdam Institute for Climate Impact Research  
  ! Date: 10.6.16    

  ! Purpose: Calculate discharge parameters.

  implicit none

  real(dp), intent(in) :: dtime

  real(dp) :: dtime_mar_coa

#if !defined(ALLOW_TAPENADE) /* Normal */

  disc = DISC 
  c_dis_0   = C_DIS_0
  c_dis_fac = C_DIS_FAC
  
  m_H = M_H
  m_D = M_D
  
  r_mar_eff = R_MAR_EFF *1000.0_dp   ! km -> m

#else /* Tapenade */

  disc_DW = DISC
  c_dis_0_DW   = C_DIS_0
  c_dis_fac_DW = C_DIS_FAC
  
  m_H_DW = M_H
  m_D_DW = M_D
  
  r_mar_eff_DW = R_MAR_EFF *1000.0_dp   ! km -> m

#endif /* Normal vs. Tapenade */

#if !defined(ALLOW_TAPENADE) /* Normal */

#if (defined(S_DIS))
  s_dis = S_DIS
#else
  s_dis = 1.0_dp   ! default value
#endif

#else /* Tapenade */

#if (defined(S_DIS))
  s_dis_DW = S_DIS
#else
  s_dis_DW = 1.0_dp   ! default value
#endif

#endif /* Normal vs. Tapenade */

#if !defined(ALLOW_TAPENADE) /* Normal */

#if (defined(ALPHA_SUB))
  alpha_sub = ALPHA_SUB
#else
  alpha_sub = 0.5_dp   ! default value
#endif

#else /* Tapenade */

#if (defined(ALPHA_SUB))
  alpha_sub_DW = ALPHA_SUB
#else
  alpha_sub_DW = 0.5_dp   ! default value
#endif

#endif /* Normal vs. Tapenade */

#if !defined(ALLOW_TAPENADE) /* Normal */

#if (defined(ALPHA_O))
  alpha_o = ALPHA_O
#else
  alpha_o = 1.0_dp   ! default value
#endif

  T_sea_freeze = 0.0_dp - DELTA_TM_SW   ! freezing temperature of sea water

#else /* Tapenade */

#if (defined(ALPHA_O))
  alpha_o_DW = ALPHA_O
#else
  alpha_o_DW = 1.0_dp   ! default value
#endif

  T_sea_freeze_DW = 0.0_dp - DELTA_TM_SW   ! freezing temperature of sea water

#endif /* Normal vs. Tapenade */

#if !defined(ALLOW_TAPENADE) /* Normal */
  n_discharge_call = -1
#else /* Tapenade */
  n_discharge_call_DW = -1
#endif /* Normal vs. Tapenade */

#if (defined(DTIME_MAR_COA0))
  dtime_mar_coa = DTIME_MAR_COA0*year2sec   ! a -> s
#else
  errormsg = ' >>> disc_param: DTIME_MAR_COA0 not defined in header file!'
  call error(errormsg)
#endif

#if !defined(ALLOW_TAPENADE) /* Normal */
  if (.not.approx_integer_multiple(dtime_mar_coa, dtime, eps_sp_dp)) then
     errormsg = ' >>> disc_param: dtime_mar_coa must be a multiple of dtime!'
     call error(errormsg)
  end if
#else /* Tapenade */
     print *, ' >>> disc_param: compare_float_m not used in adjoint'
     print *, '                 applications! double check ' 
     print *, '                 dtime_mar_coa and dtime are multiples'
#endif /* Normal vs. Tapenade */

#if !defined(ALLOW_TAPENADE) /* Normal */
  iter_mar_coa = nint(dtime_mar_coa/dtime)
#else /* Tapenade */
  iter_mar_coa_DW = nint(dtime_mar_coa/dtime)
#endif /* Normal vs. Tapenade */

#if (GRID > 1)
  errormsg = ' >>> disc_param: GRID==2 not allowed for this application!'
  call error(errormsg)
#endif

#if (ANF_DAT != 1  && defined(EXEC_MAKE_C_DIS_0))
  errormsg = ' >>> disc_param: EXEC_MAKE_C_DIS_0 defined requires ANF_DAT==1!'
  call error(errormsg)
#endif

  end subroutine disc_param


!-------------------------------------------------------------------------------
!> Constant in ice discharge parameterization (Greenland).
!! [Determine (amount of magnitude of) constant in ice discharge
!! parameterization.]
!<------------------------------------------------------------------------------
  subroutine calc_c_dis_0(dxi, deta)

  ! Author: Reinhard Calov
  ! Institution: Potsdam Institute for Climate Impact Research
  ! Date: 8.6.16    
  ! Purpose: Compute c_dis_0 indirectly from present-day topography 
  !          from the condition that the parameterized total discharge 
  !          equals the discharge give by dis_target. Note: This c_dis_0 could
  !          differ from the optimal one yielded from the siumualted ice thickness.
  !
  !          This is very useful, because c_dis_0 can vary a lot for different powers 
  !          m_D and m_H. This way we easier yield the approximately right value
  !          for c_dis_0 for given powers m_D and m_H.

  ! This routine is not to be used regularly, and it is only executed if the
  ! parameter EXEC_MAKE_C_DIS_0 is defined in the header file.

  implicit none

  real(dp), intent(in) :: dxi, deta

  integer(i4b) :: i, j
  real(dp) :: disc_target, disc_tot, dT_sub

  ! targeted total ice discharge
  ! 350 Gt/yr Calov et al. (2015)   !!! 400 Gt/yr van den Broeke et al (2016)
  disc_target = 350.0_dp ! in Gt/yr  

  H   = zs-zb
  H_c = H
  H_t = 0.0_dp

#if !defined(ALLOW_TAPENADE) /* Normal */
  c_dis_0   =  1.0_dp
  c_dis_fac =  1.0_dp
  s_dis     =  1.0_dp
#else /* Tapenade */
  c_dis_0_DW   =  1.0_dp
  c_dis_fac_DW =  1.0_dp
  s_dis_DW     =  1.0_dp
#endif /* Normal vs. Tapenade */

  call disc_fields()

  ! ensure that we get present-day discharge here (disc=1). 
  dT_glann = 0.0_dp

#if !defined(ALLOW_TAPENADE) /* Normal */
  n_discharge_call = -1
#else /* Tapenade */
  n_discharge_call_DW = -1
#endif /* Normal vs. Tapenade */

  call discharge(dxi, deta)

  !  disc_tot=0.0_dp
  !  do i=0, IMAX
  !  do j=0, JMAX
  !    if(mask(j,i).eq.0.or.mask(j,i).eq.3) then
  !      disc_tot=disc_tot+dis_perp(j,i)*cell_area(j,i)
  !    end if
  !  end do
  !  end do
  !
  !  disc_tot = disc_tot*year2sec         ! m^3/s -> m^3/a

  disc_tot = 0.0_dp

  do i=1, IMAX-1
  do j=1, JMAX-1

     if (mask_mar(j,i) == 1) then
        disc_tot = disc_tot + dis_perp(j,i)*cell_area(j,i)
     end if

  end do
  end do

  disc_tot = disc_tot *year2sec *(RHO/RHO_W)
                  ! m^3/s ice equiv. -> m^3/a water equiv.

  disc_target = disc_target*1.0e+12_dp/RHO_W   ! Gt/a -> m^3/a water equiv.

  write(6, fmt='(a)') ' '
  write(6, fmt='(a)') ' * * * * * * * * * * * * * * * * * * * * * * * * * * * * '

  write(6, fmt='(a)') ' '

#if !defined(ALLOW_TAPENADE) /* Normal */
  write(6, fmt='(3x,a,es12.4)') 'c_dis_0_init = ', c_dis_0
#else /* Tapenade */
  write(6, fmt='(3x,a,es12.4)') 'c_dis_0_init = ', c_dis_0_DW
#endif /* Normal vs. Tapenade */

#if !defined(ALLOW_TAPENADE) /* Normal */
  c_dis_0 = c_dis_0 * disc_target/disc_tot
#else /* Tapenade */
  c_dis_0_DW = c_dis_0_DW * disc_target/disc_tot
#endif /* Normal vs. Tapenade */

  write(6, fmt='(a)') ' '
  write(6, fmt='(3x,a,es12.4)') 'disc_target  = ', disc_target
  write(6, fmt='(3x,a,es12.4)') 'disc_tot     = ', disc_tot

  write(6, fmt='(a)') ' '

#if !defined(ALLOW_TAPENADE) /* Normal */
  write(6, fmt='(3x,a,es12.4)') '--> c_dis_0  = ', c_dis_0
#else /* Tapenade */
  write(6, fmt='(3x,a,es12.4)') '--> c_dis_0  = ', c_dis_0_DW
#endif /* Normal vs. Tapenade */

  end subroutine calc_c_dis_0

!-------------------------------------------------------------------------------
!> Dependence of ice discharge coefficient on latitude (Greenland).
!! [Determine dependence of ice discharge coefficient on latitude.
!! This can be improved. For now I recommend s_dis=1.]
!<------------------------------------------------------------------------------
  subroutine disc_fields()

  ! Author: Reinhard Calov
  ! Institution: Potsdam Institute for Climate Impact Research  
  ! Date: 8.6.16    

  !   Purpose: Assign fields relevant for ice discharge parameterization
  !            1. Discharge coefficient field using linear dependence on 
  !               latitude.
  !            2. Subsurface temperature dependence on longitude and latitude.

  implicit none

  integer(i4b) :: i, j
  real(dp)     :: delta_phi=25.0_dp ! approximately the phi
                                    ! spanning the middle of domain

#if !defined(ALLOW_TAPENADE) /* Normal */

  c_dis = c_dis_0*c_dis_fac &
                 *(1.0_dp-(1.0_dp-s_dis)*(phi*rad2deg-60.0_dp)/delta_phi)

  c_dis = max(c_dis, 0.0_dp)

#else /* Tapenade */

  c_dis_DW = c_dis_0_DW*c_dis_fac_DW &
                 *(1.0_dp-(1.0_dp-s_dis_DW)*(phi*rad2deg-60.0_dp)/delta_phi)

  c_dis_DW = max(c_dis_DW, 0.0_dp)

#endif /* Normal vs. Tapenade */

#if !defined(ALLOW_TAPENADE) /* Normal */
  T_sub_PD = 3.0_dp ! Can depend on lambda, phi later on
#else /* Tapenade */
  T_sub_PD_DW = 3.0_dp ! Can depend on lambda, phi later on
#endif /* Normal vs. Tapenade */

  end subroutine disc_fields

!-------------------------------------------------------------------------------
!> Ice discharge parameterization main formula, controler (general).
!! [Compute ice discharge via a parameterization using distance of ice 
!! margin to coast and ice thickness as parameters.]
!<------------------------------------------------------------------------------
  subroutine discharge(dxi, deta)

#ifdef ALLOW_TAPENADE /* Tapenade */
  use ctrl_m, only: myfloor
#endif /* Tapenade */

  ! Authors: Reinhard Calov, Andrey Ganopolski
  ! Institution: Potsdam Institute for Climate Impact Research  
  ! Date: 11.6.16

  implicit none

  real(dp), intent(in)    :: dxi, deta

  real(dp), parameter :: alpha_tc=60.0_dp ! maximal allowed angle
  real(dp) :: cos_tc

#ifdef ALLOW_TAPENADE /* Tapenade */
  integer(i4b) :: i, j, valmodint
  real(dp) :: tmp, valmod
#endif /* Tapenade */

#if !defined(ALLOW_TAPENADE) /* Normal */
  n_discharge_call = n_discharge_call + 1
#else /* Tapenade */
  n_discharge_call_DW = n_discharge_call_DW + 1
#endif /* Normal vs. Tapenade */

  cos_tc=dcos(alpha_tc*deg2rad)

#if !defined(ALLOW_TAPENADE) /* Normal */

  if ((mod(n_discharge_call, iter_mar_coa)==0).or.(n_discharge_call==1)) then
     write(6, fmt='(10x,a)') 'Computation of mask_mar, cst_dist, cos_grad_tc'
     call marginal_ring(dxi, deta)
     call coastal_distance(dxi, deta)
  end if

#else /* Tapenade */

  ! this is mod broken down algorithmically:
  tmp = real(n_discharge_call_DW)/real(iter_mar_coa_DW)
  call myfloor(tmp, valmodint)
  valmod = real(n_discharge_call_DW) - real(iter_mar_coa_DW) * real(valmodint)
  if ((valmod==0) .or. (n_discharge_call_DW==1)) then
     write(6, fmt='(10x,a)') 'Computation of mask_mar, cst_dist, cos_grad_tc'
     call marginal_ring(dxi, deta)
     call coastal_distance(dxi, deta)
  end if

#endif /* Normal vs. Tapenade */

#if !defined(ALLOW_TAPENADE) /* Normal */

  if(disc.ge.1) then !------------------ disc >= 1
    where(mask_mar.eq.1.and.cos_grad_tc.ge.cos_tc) 
      dis_perp=c_dis*H**m_H/cst_dist**m_D
    elsewhere
      dis_perp=0.0_dp
    end where
    
    if(disc.eq.2) then !----------------- disc = 2
      call calc_sub_oc_dT(T_sub_PD, dT_glann, dT_sub)
      dis_perp=dis_perp*(1.0_dp+alpha_sub*dT_sub) ! actual discharge present-day one corrected 
                                                ! via sub-ocean temperature anomaly
    end if
    dis_perp=max(0.0_dp, dis_perp) ! ensure positive values
  else if(disc.eq.0) then !-------------- disc = 0
    dis_perp=0.0_dp
  end if

#else /* Tapenade */

  if(disc_DW.ge.1) then !------------------ disc >= 1
   do i=0, IMAX
   do j=0, JMAX
      if (mask_mar(j,i).eq.1 .and. cos_grad_tc(j,i).ge.cos_tc) then
         dis_perp(j,i) = c_dis_DW(j,i)*H(j,i)**m_H_DW/cst_dist(j,i)**m_D_DW
      else
         dis_perp(j,i) = 0.0_dp
      end if
   end do
   end do
    if(disc_DW.eq.2) then !----------------- disc = 2
      call calc_sub_oc_dT(T_sub_PD_DW, dT_glann, dT_sub)
      dis_perp=dis_perp*(1.0_dp+alpha_sub_DW*dT_sub) ! actual discharge present-day one corrected 
                                                ! via sub-ocean temperature anomaly
    end if
    dis_perp=max(0.0_dp, dis_perp) ! ensure positive values
  else if(disc_DW.eq.0) then !-------------- disc = 0
    dis_perp=0.0_dp
  end if

#endif /* Normal vs. Tapenade */

  end subroutine discharge

!-------------------------------------------------------------------------------
!> Distance to the coast (general).
!! [Compute distance to the coast for every land point.]
!<------------------------------------------------------------------------------
  subroutine coastal_distance(dxi, deta)

  ! Author: Reinhard Calov
  ! Institution: Potsdam Institute for Climate Impact Research  
  ! Date: 2.11.11
  ! Methode: Options are:
  !          i_search=1: Brute force.
  !          i_search=2: Defining a square with two grid distances side length around the 
  !                      actual grid point and enlanging the square successively by
  !                      one grid. Because this should have been a circle, the
  !                      square is enlarged again if the coast line was found; this time 
  !                      enlargement is by a rectangles a the respective four side of 
  !                      the square. Smaller side lenght equals the rectangle root of 2 of the
  !                      larger side. 

#ifdef ALLOW_TAPENADE /* Tapenade */
  use ctrl_m, only: myceiling
#endif /* Tapenade */

  implicit none

  real(dp), intent(in) :: dxi, deta

  ! --------------------------------------------------
  real(dp), parameter :: grid_dist=5000.0_dp
  real(dp), parameter :: c_smooth=0.2_dp
  integer(i4b), parameter :: i_search=2
  ! --------------------------------------------------
  integer(i4b) :: i, j, i_pos, j_pos, l, l_e, d_l
  real(dp) :: dxi_inv, deta_inv
  real(dp) :: cst_dist_tmp, qrt_1, qrt_2
  logical :: leave_loop

  real(dp), dimension(0:JMAX,0:IMAX) :: zs_tmp, &
                                        grad_zs_x, grad_zs_y, &
                                        grad_dist_x, grad_dist_y

    select case (i_search)

    case(1)
    ! brute force computation of minimal distance to coast for every land point

      do i_pos=0, IMAX
      do j_pos=0, JMAX
        if(mask(j_pos,i_pos).ne.2) then
          cst_dist(j_pos,i_pos)=1.d+20
          do i=0, IMAX
          do j=0, JMAX
            if(mask(j,i).eq.2) then
              cst_dist_tmp=sqrt((xi(i_pos)-xi(i))**2+(eta(j_pos)-eta(j))**2)
              if(cst_dist_tmp.le.cst_dist(j_pos,i_pos)) then
                cst_dist(j_pos,i_pos)=cst_dist_tmp
              end if
            end if
          end do
          end do
        else
          cst_dist(j_pos,i_pos)=grid_dist
        end if
      end do
      end do

    case(2)

      do i_pos=0, IMAX
      do j_pos=0, JMAX
        if(mask(j_pos,i_pos).ne.2) then
          leave_loop=.false.
          cst_dist(j_pos,i_pos)=1.d+20

#if !defined(ALLOW_TAPENADE) /* Normal */
          do l=1, max(IMAX,JMAX)
#else /* Tapenade */
          l = 1
          do while (l<= max(IMAX,JMAX) .and. leave_loop.eqv..false.)
#endif /* Normal vs. Tapenade */

            do i=max(i_pos-l,0), min(i_pos+l,IMAX)
              j=max(j_pos-l,0); j=min(j,JMAX)
              if(mask(j,i).eq.2) then
                leave_loop=.true.
                cst_dist_tmp=sqrt((xi(i_pos)-xi(i))**2+(eta(j_pos)-eta(j))**2)
                if(cst_dist_tmp.le.cst_dist(j_pos,i_pos)) then
                  cst_dist(j_pos,i_pos)=cst_dist_tmp
                end if
              end if
              j=min(j_pos+l,JMAX)
              if(mask(j,i).eq.2) then
                leave_loop=.true.
                cst_dist_tmp=sqrt((xi(i_pos)-xi(i))**2+(eta(j_pos)-eta(j))**2)
                if(cst_dist_tmp.le.cst_dist(j_pos,i_pos)) then
                  cst_dist(j_pos,i_pos)=cst_dist_tmp
                end if
              end if
            end do
            do j=max(j_pos-l+1,0), min(j_pos+l-1,JMAX)
              i=max(i_pos-l,0); i=min(i,IMAX)
              if(mask(j,i).eq.2) then
                leave_loop=.true.
                cst_dist_tmp=sqrt((xi(i_pos)-xi(i))**2+(eta(j_pos)-eta(j))**2)
                if(cst_dist_tmp.le.cst_dist(j_pos,i_pos)) then
                  cst_dist(j_pos,i_pos)=cst_dist_tmp
                end if
              end if
              i=min(i_pos+l,IMAX)
              if(mask(j,i).eq.2) then
                leave_loop=.true.
                cst_dist_tmp=sqrt((xi(i_pos)-xi(i))**2+(eta(j_pos)-eta(j))**2)
                if(cst_dist_tmp.le.cst_dist(j_pos,i_pos)) then
                  cst_dist(j_pos,i_pos)=cst_dist_tmp
                end if
              end if
            end do

#if !defined(ALLOW_TAPENADE) /* Normal */
            if(leave_loop) then
              l_e=l
              d_l=ceiling(l_e*sqrt(2.0_dp)+0.001_dp)
              exit   ! leave loop in l
            end if
#else /* Tapenade */
            if(leave_loop) then
              l_e=l
              call myceiling((l_e*sqrt(2.0_dp)+0.001_dp), d_l)
            end if
          l = l+1
#endif /* Normal vs. Tapenade */

          end do
          ! left
          do i=max(i_pos-(l_e+d_l),0),min(i_pos-(l_e-1),IMAX)
          do j=max(j_pos-l_e,0),min(j_pos+l_e,JMAX)
            if(mask(j,i).eq.2) then
              cst_dist_tmp=sqrt((xi(i_pos)-xi(i))**2+(eta(j_pos)-eta(j))**2)
              if(cst_dist_tmp.le.cst_dist(j_pos,i_pos)) then
                cst_dist(j_pos,i_pos)=cst_dist_tmp
              end if
            end if
          end do
          end do
          ! right
          do i=max(i_pos+l_e+1,0),min(i_pos+l_e+d_l,IMAX)
          do j=max(j_pos-l_e,0),min(j_pos+l_e,JMAX)
            if(mask(j,i).eq.2) then
              cst_dist_tmp=sqrt((xi(i_pos)-xi(i))**2+(eta(j_pos)-eta(j))**2)
              if(cst_dist_tmp.le.cst_dist(j_pos,i_pos)) then
                cst_dist(j_pos,i_pos)=cst_dist_tmp
              end if
            end if
          end do
          end do
          ! lower
          do i=max(i_pos-l_e,0),min(i_pos+l_e,IMAX)
          do j=max(j_pos-(l_e+d_l),0),min(j_pos-(l_e-1),JMAX)
            if(mask(j,i).eq.2) then
              cst_dist_tmp=sqrt((xi(i_pos)-xi(i))**2+(eta(j_pos)-eta(j))**2)
              if(cst_dist_tmp.le.cst_dist(j_pos,i_pos)) then
                cst_dist(j_pos,i_pos)=cst_dist_tmp
              end if
            end if
          end do
          end do
          ! upper
          do i=max(i_pos-l_e,0),min(i_pos+l_e,IMAX)
          do j=max(j_pos+l_e+1,0),min(j_pos+l_e+d_l,JMAX)
            if(mask(j,i).eq.2) then
              cst_dist_tmp=sqrt((xi(i_pos)-xi(i))**2+(eta(j_pos)-eta(j))**2)
              if(cst_dist_tmp.le.cst_dist(j_pos,i_pos)) then
                cst_dist(j_pos,i_pos)=cst_dist_tmp
              end if
            end if
          end do
          end do
        else
          cst_dist(j_pos,i_pos)=grid_dist
        end if
      end do
      end do

    end select

  ! test GIS only

    if(.false.) then
      do i=0, IMAX
      do j=0, JMAX
        if(xi(i).le.-350000.0_dp &
          .and.-2800000.0_dp.le.eta(j) &
          .and.eta(j).le.-2250000.0_dp &
          .or. &
          xi(i).ge.100000.0_dp &
          .and.-2280000.0_dp.le.eta(j) &
          .and.eta(j).le.-1880000.0_dp &
          .or. &
          eta(j).ge.-1250000.0_dp &
          .and.0.0_dp.le.xi(i) &
          .and.xi(i).le.230000.0_dp) then
             cst_dist(j,i)=0.3_dp*cst_dist(j,i)
        end if
      end do
      end do
    end if

  ! forbid zero, something like grid distance is reasonable here.
  ! Dependents of coordinate system, for spherical coordinate system: 
  ! take length on Earth's surface. For now, simply just 5 km are taken.

    cst_dist=max(grid_dist, cst_dist)

  ! forbit to large values
    cst_dist=min(1.0e10_dp, cst_dist)

    zs_tmp=zs

  ! smoothing the surface topography
    do i=1, IMAX-1
    do j=1, JMAX-1
       zs_tmp(j,i) = c_smooth*((zs(j,i-1)+zs(j,i+1))+(zs(j-1,i)+zs(j+1,i))) &
                     +(1.0_dp-4.0_dp*c_smooth)*zs(j,i)
    end do
    end do

  ! gradients (2. order for now ...)

    dxi_inv  = 1.0_dp/dxi
    deta_inv = 1.0_dp/deta

!  ------ x-derivatives

    do i=1, IMAX-1
    do j=0, JMAX
      grad_zs_x(j,i)    = (zs_tmp(j,i+1)-zs_tmp(j,i-1))*0.5_dp*dxi_inv*insq_g11_g(j,i)
      grad_dist_x(j,i)  = (cst_dist(j,i+1)-cst_dist(j,i-1))*0.5_dp*dxi_inv*insq_g11_g(j,i)
    end do
    end do

    do j=0, JMAX
      grad_zs_x(j,0)       = (zs_tmp(j,1)-zs_tmp(j,0))*dxi_inv*insq_g11_g(j,0)
      grad_zs_x(j,IMAX)    = (zs_tmp(j,IMAX)-zs_tmp(j,IMAX-1))*dxi_inv*insq_g11_g(j,IMAX)
      grad_dist_x(j,0)     = (cst_dist(j,1)-cst_dist(j,0))*dxi_inv*insq_g11_g(j,0)
      grad_dist_x(j,IMAX)  = (cst_dist(j,IMAX)-cst_dist(j,IMAX-1))*dxi_inv*insq_g11_g(j,IMAX)
    end do

!  ------ y-derivatives

    do i=0, IMAX
    do j=1, JMAX-1
      grad_zs_y(j,i)    = (zs_tmp(j+1,i)-zs_tmp(j-1,i))*0.5_dp*deta_inv*insq_g22_g(j,i)
      grad_dist_y(j,i)  = (cst_dist(j+1,i)-cst_dist(j-1,i))*0.5_dp*deta_inv*insq_g22_g(j,i)
    end do
    end do

    do i=0, IMAX
      grad_zs_y(0,i)       = (zs_tmp(1,i)-zs_tmp(0,i))*deta_inv*insq_g22_g(0,i)
      grad_zs_y(JMAX,i)    = (zs_tmp(JMAX,i)-zs_tmp(JMAX-1,i))*deta_inv*insq_g22_g(JMAX,i)
      grad_dist_y(0,i)     = (cst_dist(1,i)-cst_dist(0,i))*deta_inv*insq_g22_g(0,i)
      grad_dist_y(JMAX,i)  = (cst_dist(JMAX,i)-cst_dist(JMAX-1,i))*deta_inv*insq_g22_g(JMAX,i)
    end do

  ! angle between topographical gradient and coastal distance gradient

    do i=0, IMAX
    do j=0, JMAX
      qrt_1=grad_zs_x(j,i)*grad_zs_x(j,i)+grad_zs_y(j,i)*grad_zs_y(j,i)
      qrt_2=grad_dist_x(j,i)*grad_dist_x(j,i)+grad_dist_y(j,i)*grad_dist_y(j,i)
      if(qrt_1.ne.0.0_dp.and.qrt_2.ne.0.0_dp) then
  !!! top_cst_alpha(j,i)=dacos( (grad_zs_x(j,i)*grad_dist_x(j,i)+grad_zs_y(j,i)*grad_dist_y(j,i))/ &
  !!! (sqrt(qrt_1)*sqrt(qrt_2)) )

        cos_grad_tc(j,i)= (grad_zs_x(j,i)*grad_dist_x(j,i)+grad_zs_y(j,i)*grad_dist_y(j,i))/ &
        (sqrt(qrt_1)*sqrt(qrt_2))
      else
        cos_grad_tc(j,i)=-1.0_dp
      end if
    end do
    end do

  end subroutine coastal_distance

!-------------------------------------------------------------------------------
!> Ring along an ice sheet margin (general).
!! [Compute marginal ring.]
!<------------------------------------------------------------------------------
  subroutine marginal_ring(dxi, deta)

  ! Author: Reinhard Calov
  ! Institution: Potsdam Institute for Climate Impact Research  
  ! Date: 28.10.11

  ! Purpose: Determine an r_max_eff wide ring arong ice margins 
  ! towards the interior of the ice sheet. Small stripes of land are
  ! not considered as land. For the logics De Morgan is used often.

  ! Methode: Two staggered loops in i,j. The inner i,j loop acts inside
  ! a rectangle defined by r_mar_eff. r_mar_eff should be small compared
  ! to the domain size. 

    implicit none

    real(dp), intent(in) :: dxi, deta

    integer(i4b) :: i, j, i_pos, j_pos, di_eff, dj_eff
    integer(i4b) :: count_tmp
    real(dp) :: r_p

#if !defined(ALLOW_TAPENADE) /* Normal */
    if(r_mar_eff.le.1.0e6_dp) then
#else /* Tapenade */
    if(r_mar_eff_DW.le.1.0e6_dp) then
#endif /* Normal vs. Tapenade */

      mask_mar=0
      do i_pos=1, IMAX-1
      do j_pos=1, JMAX-1
        if((mask(j_pos,i_pos).eq.1.or.mask(j_pos,i_pos).eq.2).and.     &
           .not.(mask(j_pos,i_pos).eq.1.and.mask(j_pos,i_pos-1).eq.0   &
                 .and.mask(j_pos,i_pos+1).eq.0.or. &
                  mask(j_pos,i_pos).eq.1.and.mask(j_pos-1,i_pos).eq.0  &
                  .and.mask(j_pos+1,i_pos).eq.0.or. &
                   mask(j_pos,i_pos).eq.1.and.mask(j_pos,i_pos-1).eq.3 &
                   .and.mask(j_pos,i_pos+1).eq.3.or. &
                   mask(j_pos,i_pos).eq.1.and.mask(j_pos-1,i_pos).eq.3 &
                   .and.mask(j_pos+1,i_pos).eq.3.or. &
                   mask(j_pos,i_pos).eq.1.and.mask(j_pos,i_pos-1).eq.3 &
                   .and.mask(j_pos,i_pos+1).eq.0.or. &
                   mask(j_pos,i_pos).eq.1.and.mask(j_pos-1,i_pos).eq.0 &
                   .and.mask(j_pos+1,i_pos).eq.3)) then ! outside ice sheet, exclude isolated land stripes

#if !defined(ALLOW_TAPENADE) /* Normal */
          di_eff=int(r_mar_eff/dxi)+1; dj_eff=int(r_mar_eff/deta)+1 ! only for grid=0, 1 yet!
#else /* Tapenade */
          di_eff=int(r_mar_eff_DW/dxi)+1; dj_eff=int(r_mar_eff_DW/deta)+1 ! only for grid=0, 1 yet!
#endif /* Normal vs. Tapenade */

          do i=max(i_pos-di_eff,0), min(i_pos+di_eff,IMAX)
          do j=max(j_pos-dj_eff,0), min(j_pos+dj_eff,JMAX)
            r_p=sqrt((xi(i_pos)-xi(i))**2+(eta(j_pos)-eta(j))**2)

#if !defined(ALLOW_TAPENADE) /* Normal */
            if(r_p.le.r_mar_eff.and.(mask(j,i).eq.0.or.mask(j,i).eq.3)) then
#else /* Tapenade */
            if(r_p.le.r_mar_eff_DW.and.(mask(j,i).eq.0.or.mask(j,i).eq.3)) then
#endif /* Normal vs. Tapenade */

              mask_mar(j,i)=1
            end if
          end do
          end do
        end if
      end do
      end do

#if !defined(ALLOW_TAPENADE) /* Normal */
      where(mask.ge.1)
        mask_mar=1
      end where
#else /* Tapenade */
      do i=0,IMAX
      do j=0,JMAX
        if (mask(j,i).ge.1) then
          mask_mar(j,i)=1
        end if
      end do
      end do
#endif /* Normal vs. Tapenade */

    else ! the ring encompassed entire Greenland 
      mask_mar=1
    end if

  end subroutine marginal_ring

!-------------------------------------------------------------------------------
!> Anomaly of subsurface temperature (general).
!! [Compute anomaly of subsurface temperature with a simple parameterization.]
!<------------------------------------------------------------------------------
  subroutine calc_sub_oc_dT(T_sub_PD, dT_glann, dT_sub)

  ! Author: Reinhard Calov
  ! Institution: Potsdam Institute for Climate Impact Research  
  ! Purpose: Compute anomaly of subsurface temperature
  !          with a simple parameterization
  !          using global annual temperature
  !          anomaly (e.g. from CLIMBER)

  implicit none

  real(dp), intent(in)  :: T_sub_PD, dT_glann
  real(dp), intent(out) :: dT_sub

  real(dp) :: T_sub

#if !defined(ALLOW_TAPENADE) /* Normal */
  dT_sub=alpha_o*dT_glann                    ! Sub-ocean temperature anomaly
  T_sub=max(T_sea_freeze, T_sub_PD)+dT_sub   ! Actual sub-ocean temperature
  T_sub=max(T_sea_freeze, T_sub)             ! Ocean temperature should stay
                                             ! above freezing point
#else /* Tapenade */
  dT_sub=alpha_o_DW*dT_glann                    ! Sub-ocean temperature anomaly
  T_sub=max(T_sea_freeze_DW, T_sub_PD)+dT_sub   ! Actual sub-ocean temperature
  T_sub=max(T_sea_freeze_DW, T_sub)             ! Ocean temperature should stay
                                                ! above freezing point
#endif /* Normal vs. Tapenade */

  dT_sub=T_sub-T_sub_PD

  end subroutine calc_sub_oc_dt

!-------------------------------------------------------------------------------

end module discharge_workers_m
!
