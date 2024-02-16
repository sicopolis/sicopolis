!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ t h k _ w a t e r _ b a s _ m
!
!! Computation of the thickness of the water column under the ice base.
!!
!!##### Authors
!!
!! Ralf Greve, Marius Schaefer
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
!> Computation of the thickness of the water column under the ice base.
!-------------------------------------------------------------------------------
module calc_thk_water_bas_m

  use sico_types_m
  use sico_variables_m
  use sico_vars_m
  use error_m

#if (BASAL_HYDROLOGY==1)
  use hydro_m
#endif

  implicit none

  private
  public :: calc_thk_water_bas

contains

!-------------------------------------------------------------------------------
!> Main subroutine of calc_thk_water_bas_m:
!! Computation of the thickness of the water column under the ice base.
!-------------------------------------------------------------------------------
  subroutine calc_thk_water_bas()

  implicit none

  integer(i4b) :: i, j

  logical, save :: firstcall = .true.

#if (BASAL_HYDROLOGY==1)
  real(dp), save                     :: rho_rho_w_ratio
  integer , dimension(0:IMAX,0:JMAX) :: hydro_icemask
  real(dp), dimension(0:IMAX,0:JMAX) :: hydro_topg, hydro_thk, &
                                        hydro_temppabase, hydro_supply, &
                                        hydro_sflux, hydro_vflux, &
                                        hydro_vfluxX, hydro_vfluxY, &
                                        hydro_bwat
  type(hydro_t), save :: hydro
#endif

!-------- Water column --------

#if (BASAL_HYDROLOGY==1)

  if (firstcall) then

     rho_rho_w_ratio = RHO/RHO_W

     call hydro_init(hydro, xi, eta)
     call hydro_gen_conf(hydro, &
          & method='quinn', &
          & avoid_frz=.false., &
          & filter_len=0.0_dp, &
          & rho_seawater=RHO_SW, &
          & rho_freshwater=RHO_W, &
          & rho_ice=RHO)

  end if

  do i=0, IMAX
  do j=0, JMAX

     hydro_topg(i,j)       = zl(j,i)-z_sl(j,i)
     hydro_temppabase(i,j) = temph_b(j,i)

     if (mask(j,i)==0) then   ! grounded ice
        hydro_icemask(i,j) = 1
        hydro_thk(i,j)     = H(j,i)
#if (!defined(MELT_DRAIN) || MELT_DRAIN==0)
        hydro_supply(i,j)  = rho_rho_w_ratio*Q_b_tot(j,i)
#elif (MELT_DRAIN==1)   /* drainage of surface melt water included */
        hydro_supply(i,j)  = rho_rho_w_ratio*(Q_b_tot(j,i) + runoff(j,i))
#else
        errormsg = ' >>> calc_thk_water_bas: MELT_DRAIN must be 0 or 1!'
        call error(errormsg)
#endif
     else
        hydro_icemask(i,j) = 0
        hydro_thk(i,j)     = 0.0_dp
        hydro_supply(i,j)  = 0.0_dp
     end if

  end do
  end do

  call hydro_set_topg(hydro, hydro_topg)
  call hydro_set_thk(hydro, hydro_thk)
  call hydro_set_temppabase(hydro, hydro_temppabase)
  call hydro_set_supply(hydro, hydro_supply)
  call hydro_set_mask(hydro, hydro_icemask)

  call hydro_update(hydro)

  call hydro_get_sflux(hydro, hydro_sflux)
  call hydro_get_vflux(hydro, hydro_vflux, hydro_vfluxX, hydro_vfluxY)
  call hydro_get_bwat(hydro, hydro_bwat)

  do i=0, IMAX
  do j=0, JMAX

     q_w(j,i)   = hydro_vflux(i,j)
     q_w_x(j,i) = hydro_vfluxX(i,j)
     q_w_y(j,i) = hydro_vfluxY(i,j)
     H_w(j,i)   = hydro_bwat(i,j)

  end do
  end do

#else

  do i=0, IMAX
  do j=0, JMAX

     if (mask(j,i)==0) then   ! grounded ice
        q_w(j,i)   = 0.0_dp
        q_w_x(j,i) = 0.0_dp
        q_w_y(j,i) = 0.0_dp
        H_w(j,i)   = 0.0_dp
     end if

  end do
  end do

#endif

  do i=0, IMAX
  do j=0, JMAX

     if (mask(j,i)==2) then   ! ocean
        q_w(j,i)   = 0.0_dp
        q_w_x(j,i) = 0.0_dp
        q_w_y(j,i) = 0.0_dp
        H_w(j,i)   = z_sl(j,i)-zl(j,i)
     else if (mask(j,i)==3) then   ! floating ice
        q_w(j,i)   = 0.0_dp
        q_w_x(j,i) = 0.0_dp
        q_w_y(j,i) = 0.0_dp
        H_w(j,i)   = zb(j,i)-zl(j,i)
     else if (mask(j,i)==1) then   ! ice-free land
        q_w(j,i)   = 0.0_dp
        q_w_x(j,i) = 0.0_dp
        q_w_y(j,i) = 0.0_dp
        H_w(j,i)   = 0.0_dp
     end if

  end do
  end do

  if (firstcall) firstcall = .false.

  end subroutine calc_thk_water_bas

!-------------------------------------------------------------------------------

end module calc_thk_water_bas_m
!
