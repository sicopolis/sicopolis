!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  g l o b a l s
!
!> @file
!!
!! TAPENADE does not like the save statement and instead prefers the global
!! variables to be declared in a separate file.
!!
!! @section Copyright
!!
!! Copyright 2009-2024 Ralf Greve, Shreyas Sunil Gaikwad, Liz Curry-Logan
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
!> TAPENADE does not like the save statement and instead prefers the global
!! variables to be declared in a separate file.
!<++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module globals
        
  use sico_types_m

  implicit none

  real(dp), dimension(0:JMAX,0:IMAX) :: dzs_dx_aux, dzs_dy_aux

  integer(i4b) :: disc_DW
  integer(i4b) :: n_discharge_call_DW, iter_mar_coa_DW
  real(dp)     :: c_dis_0_DW, s_dis_DW, c_dis_fac_DW
  real(dp)     :: T_sub_PD_DW, alpha_sub_DW, alpha_o_DW, m_H_DW, m_D_DW, r_mar_eff_DW
  real(dp)     :: T_sea_freeze_DW
  real(dp)     :: dT_glann, dT_sub
  integer(i4b), dimension(0:JMAX,0:IMAX) :: mask_mar
  real(dp),     dimension(0:JMAX,0:IMAX) :: c_dis_DW
  real(dp),     dimension(0:JMAX,0:IMAX) :: cst_dist, cos_grad_tc
  real(dp),     dimension(0:JMAX,0:IMAX) :: dis_perp 

  real(dp), dimension(-256:255), public       :: c_int_table
  real(dp), dimension(-524288:524287), public :: c_int_inv_table
  integer(i4b), public                        :: n_temp_min
  integer(i4b), public                        :: n_temp_max
  integer(i4b), public                        :: n_enth_min
  integer(i4b), public                        :: n_enth_max
  real(dp), public                            :: L_inv
  real(dp), public                            :: L_eto

  logical      :: firstcall = .true.
  integer(i4b) :: n_year_CE_aux_save = -9999 

  real(dp) :: enh_stream
  logical  :: flag_enh_stream

end module globals
!
