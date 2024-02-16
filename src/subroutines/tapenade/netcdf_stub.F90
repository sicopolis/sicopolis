!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  n e t c d f _ s t u b
!
!! Stub file for NetCDF subroutines for SICOPOLIS-AD v2.
!!
!!##### Authors
!!
!! Shreyas Sunil Gaikwad, Ralf Greve
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
!> Stub file for NetCDF subroutines for SICOPOLIS-AD v2.
!-------------------------------------------------------------------------------
module netcdf

implicit none 

integer, parameter, public :: &
nf90_byte   = 1,            &
nf90_int    = 4,            &
nf90_float  = 5,            &
nf90_double = 6

integer, parameter, public :: nf90_noclobber = 4
integer, parameter, public :: nf90_unlimited = 0
integer, parameter, public :: nf90_global = 0
integer, parameter, public :: nf90_noerr = 0
integer, parameter, public :: nf90_nowrite = 0

interface nf90_put_att
        module procedure nf90_put_att_one_real, nf90_put_att_1D_real, &
                         nf90_put_att_1D_character, nf90_put_att_1D_integer
end interface nf90_put_att

interface nf90_get_att
    module procedure nf90_get_att_text
end interface

interface nf90_put_var
        module procedure nf90_put_var_one_real, nf90_put_var_1D_real, nf90_put_var_2D_real, nf90_put_var_3D_real
        module procedure nf90_put_var_one_integer, nf90_put_var_1D_integer, nf90_put_var_2D_integer, nf90_put_var_3D_integer
        module procedure nf90_put_var_one_logical, nf90_put_var_1D_logical, nf90_put_var_2D_logical, nf90_put_var_3D_logical
end interface nf90_put_var

interface nf90_get_var
        module procedure nf90_get_var_one_real  
        module procedure nf90_get_var_one_real, nf90_get_var_1D_real, nf90_get_var_2D_real, nf90_get_var_3D_real
        module procedure nf90_get_var_one_integer, nf90_get_var_1D_integer, nf90_get_var_2D_integer, nf90_get_var_3D_integer      
end interface nf90_get_var

interface nf90_def_var
        module procedure nf90_def_var_Scalar, nf90_def_var_oneDim, nf90_def_var_ManyDims
end interface nf90_def_var

public :: nf90_inquire_attribute

contains

function nf90_get_att_text(ncid, varid, name, values)
    integer,                          intent( in) :: ncid, varid
    character(len = *),               intent( in) :: name
    character(len = *),               intent(out) :: values
    integer                                       :: nf90_get_att_text
    values = ""
    nf90_get_att_text = 0
end function nf90_get_att_text

function nf90_inquire_attribute(ncid, varid, name, xtype, len, attnum)
    integer,             intent( in)           :: ncid, varid
    character (len = *), intent( in)           :: name
    integer,             intent(out), optional :: xtype, len, attnum
    integer                                    :: nf90_inquire_attribute
    integer                          :: local_xtype, local_len
    nf90_inquire_attribute = 0
end function nf90_inquire_attribute

function nf90_open(path, mode, ncid, bufrsize, cache_size, cache_nelems, &
                   cache_preemption, comm, info)
  implicit none
  character (len = *), intent(in) :: path
  integer, intent(in) :: mode
  integer, intent(out) :: ncid
  integer, optional, intent(inout) :: bufrsize
  integer, optional, intent(in) :: cache_size, cache_nelems
  real, optional, intent(in) :: cache_preemption
  integer, optional, intent(in) :: comm, info
  integer :: nf90_open
end function nf90_open


function nf90_create(path, cmode, ncid, initialsize, chunksize)
    character (len = *), intent(in   ) :: path
    integer,             intent(in   ) :: cmode
    integer,             intent(  out) :: ncid
    integer, optional,   intent(in   ) :: initialsize
    integer, optional,   intent(inout) :: chunksize
    integer                            :: nf90_create
    integer :: fileSize, chunk
    nf90_create = 0
  end function nf90_create

  function nf90_def_dim(ncid, name, len, dimid)
    integer,             intent( in) :: ncid
    character (len = *), intent( in) :: name
    integer,             intent( in) :: len
    integer,             intent(out) :: dimid
    integer                          :: nf90_def_dim
    nf90_def_dim = 0
  end function nf90_def_dim

  function nf90_close(ncid)
    integer, intent( in) :: ncid
    integer              :: nf90_close
    nf90_close = 0
  end function nf90_close

  function nf90_sync(ncid)
    integer, intent( in) :: ncid
    integer              :: nf90_sync
    nf90_sync = 0
  end function nf90_sync

  function nf90_inq_varid(ncid, name, varid)
    integer,             intent( in) :: ncid
    character (len = *), intent( in) :: name
    integer,             intent(out) :: varid
    integer                          :: nf90_inq_varid
    nf90_inq_varid = 0
  end function nf90_inq_varid

  function nf90_enddef(ncid, h_minfree, v_align, v_minfree, r_align)
    integer,           intent( in) :: ncid
    integer, optional, intent( in) :: h_minfree, v_align, v_minfree, r_align
    integer                        :: nf90_enddef
    integer :: hMinFree, VAlign, VMinFree, RAlign
    nf90_enddef = 0
  end function nf90_enddef

  function nf90_inq_dimid(ncid, name, dimid)
    integer,             intent( in) :: ncid
    character (len = *), intent( in) :: name
    integer,             intent(out) :: dimid
    integer                          :: nf90_inq_dimid
    nf90_inq_dimid = 0
  end function nf90_inq_dimid

   function nf90_put_att_one_real(ncid, varid, name, values)
    integer,                                   intent( in) :: ncid, varid
    character(len = *),                        intent( in) :: name
    real (8), intent( in) :: values
    integer                                                :: nf90_put_att_one_real
    nf90_put_att_one_real = 0
   end function nf90_put_att_one_real

   function nf90_put_att_1D_real(ncid, varid, name, values)
    integer,                                   intent( in) :: ncid, varid
    character(len = *),                        intent( in) :: name
    real (8), dimension(:), intent( in) :: values
    integer                                                :: nf90_put_att_1D_real
    nf90_put_att_1D_real = 0
   end function nf90_put_att_1D_real

   function nf90_put_att_1D_integer(ncid, varid, name, values)
    integer,                                   intent( in) :: ncid, varid
    character(len = *),                        intent( in) :: name
    integer(4), dimension(:), intent( in) :: values
    integer                                                :: nf90_put_att_1D_integer
    nf90_put_att_1D_integer = 0
   end function nf90_put_att_1D_integer

   function nf90_put_att_1D_character(ncid, varid, name, values)
    integer,                                   intent( in) :: ncid, varid
    character(len = *),                        intent( in) :: name
    character(len = *), intent( in) :: values
    integer                                                :: nf90_put_att_1D_character
    nf90_put_att_1D_character = 0
   end function nf90_put_att_1D_character

function nf90_put_var_one_real(ncid, varid, values, start, count)
  integer,                         intent( in) :: ncid, varid
           real (8), intent( in) :: values
  integer, dimension(:), optional, intent( in) :: start, count
  integer                                                :: nf90_put_var_one_real
end function nf90_put_var_one_real

function nf90_put_var_1D_real(ncid, varid, values, start, count)
  integer,                         intent( in) :: ncid, varid
           real (8), dimension(:), intent( in) :: values
  integer, dimension(:), optional, intent( in) :: start, count
  integer                                                :: nf90_put_var_1D_real
end function nf90_put_var_1D_real

function nf90_put_var_2D_real(ncid, varid, values, start, count)
  integer,                         intent( in) :: ncid, varid
           real (8), dimension(:,:), intent( in) :: values
  integer, dimension(:), optional, intent( in) :: start, count
  integer                                                :: nf90_put_var_2D_real
end function nf90_put_var_2D_real

function nf90_put_var_3D_real(ncid, varid, values, start, count)
  integer,                         intent( in) :: ncid, varid
           real (8), dimension(:,:,:), intent( in) :: values
  integer, dimension(:), optional, intent( in) :: start, count
  integer                                                :: nf90_put_var_3D_real
end function nf90_put_var_3D_real

function nf90_put_var_one_integer(ncid, varid, values, start, count)
  integer,                         intent( in) :: ncid, varid
           integer (4), intent( in) :: values
  integer, dimension(:), optional, intent( in) :: start, count
  integer                                                :: nf90_put_var_one_integer
end function nf90_put_var_one_integer

function nf90_put_var_1D_integer(ncid, varid, values, start, count)
  integer,                         intent( in) :: ncid, varid
           integer (4), dimension(:), intent( in) :: values
  integer, dimension(:), optional, intent( in) :: start, count
  integer                                                :: nf90_put_var_1D_integer
end function nf90_put_var_1D_integer

function nf90_put_var_2D_integer(ncid, varid, values, start, count)
  integer,                         intent( in) :: ncid, varid
           integer (4), dimension(:,:), intent( in) :: values
  integer, dimension(:), optional, intent( in) :: start, count
  integer                                                :: nf90_put_var_2D_integer
end function nf90_put_var_2D_integer

function nf90_put_var_3D_integer(ncid, varid, values, start, count)
  integer,                         intent( in) :: ncid, varid
           integer (4), dimension(:,:,:), intent( in) :: values
  integer, dimension(:), optional, intent( in) :: start, count
  integer                                                :: nf90_put_var_3D_integer
end function nf90_put_var_3D_integer

function nf90_put_var_one_logical(ncid, varid, values, start, count)
  integer,                         intent( in) :: ncid, varid
           logical, intent( in) :: values
  integer, dimension(:), optional, intent( in) :: start, count
  integer                                                :: nf90_put_var_one_logical
end function nf90_put_var_one_integer

function nf90_put_var_1D_logical(ncid, varid, values, start, count)
  integer,                         intent( in) :: ncid, varid
           logical, dimension(:), intent( in) :: values
  integer, dimension(:), optional, intent( in) :: start, count
  integer                                                :: nf90_put_var_1D_logical
end function nf90_put_var_1D_integer

function nf90_put_var_2D_logical(ncid, varid, values, start, count)
  integer,                         intent( in) :: ncid, varid
           logical, dimension(:,:), intent( in) :: values
  integer, dimension(:), optional, intent( in) :: start, count
  integer                                                :: nf90_put_var_2D_logical
end function nf90_put_var_2D_logical

function nf90_put_var_3D_logical(ncid, varid, values, start, count)
  integer,                         intent( in) :: ncid, varid
           logical, dimension(:,:,:), intent( in) :: values
  integer, dimension(:), optional, intent( in) :: start, count
  integer                                                :: nf90_put_var_3D_logical
end function nf90_put_var_3D_logical

function nf90_get_var_one_real(ncid, varid, values, start, count, stride, map)
  integer,                         intent( in) :: ncid, varid
                                   real (8), intent(out) :: values
  integer, dimension(:), optional, intent( in) :: start, count, stride, map
  integer                                                :: nf90_get_var_one_real
  values = 0.0
  nf90_get_var_one_real = 0
end function nf90_get_var_one_real

function nf90_get_var_1D_real(ncid, varid, values, start, count)
  integer,                         intent( in) :: ncid, varid
           real (8), dimension(:), intent( in) :: values
  integer, dimension(:), optional, intent( in) :: start, count
  integer                                                :: nf90_get_var_1D_real
  values(:) = 0.0
  nf90_get_var_1D_real = 0
end function nf90_get_var_1D_real

function nf90_get_var_2D_real(ncid, varid, values, start, count)
  integer,                         intent( in) :: ncid, varid
           real (8), dimension(:,:), intent( in) :: values
  integer, dimension(:), optional, intent( in) :: start, count
  integer                                                :: nf90_get_var_2D_real
  values(:,:) = 0.0
  nf90_get_var_2D_real = 0
end function nf90_get_var_2D_real

function nf90_get_var_3D_real(ncid, varid, values, start, count)
  integer,                         intent( in) :: ncid, varid
           real (8), dimension(:,:,:), intent( in) :: values
  integer, dimension(:), optional, intent( in) :: start, count
  integer                                                :: nf90_get_var_3D_real
  values(:,:,:) = 0.0
  nf90_get_var_3D_real = 0
end function nf90_get_var_3D_real

function nf90_get_var_one_integer(ncid, varid, values, start, count)
  integer,                         intent( in) :: ncid, varid
           integer (4), intent( in) :: values
  integer, dimension(:), optional, intent( in) :: start, count
  integer                                                :: nf90_get_var_one_integer
  values = 0
  nf90_get_var_one_integer = 0
end function nf90_get_var_one_integer

function nf90_get_var_1D_integer(ncid, varid, values, start, count)
  integer,                         intent( in) :: ncid, varid
           integer (4), dimension(:), intent( in) :: values
  integer, dimension(:), optional, intent( in) :: start, count
  integer                                                :: nf90_get_var_1D_integer
  values(:) = 0
  nf90_get_var_1D_integer = 0
end function nf90_get_var_1D_integer

function nf90_get_var_2D_integer(ncid, varid, values, start, count)
  integer,                         intent( in) :: ncid, varid
           integer (4), dimension(:,:), intent( in) :: values
  integer, dimension(:), optional, intent( in) :: start, count
  integer                                                :: nf90_get_var_2D_integer
  values(:,:) = 0
  nf90_get_var_2D_integer = 0
end function nf90_get_var_2D_integer

function nf90_get_var_3D_integer(ncid, varid, values, start, count)
  integer,                         intent( in) :: ncid, varid
           integer (4), dimension(:,:,:), intent( in) :: values
  integer, dimension(:), optional, intent( in) :: start, count
  integer                                                :: nf90_get_var_3D_integer
  values(:,:,:) = 0
  nf90_get_var_3D_integer = 0
end function nf90_get_var_3D_integer

function nf90_def_var_Scalar(ncid, name, xtype, varid)
  integer, intent(in) :: ncid
  character (len = *), intent(in) :: name
  integer, intent( in) :: xtype
  integer, intent(out) :: varid
  integer                            :: nf90_def_var_Scalar
end function nf90_def_var_Scalar

function nf90_def_var_oneDim(ncid, name, xtype, dimids, varid)
  integer, intent(in) :: ncid
  character (len = *), intent(in) :: name
  integer, intent( in) :: xtype
  integer, intent(in) :: dimids
  integer, intent(out) :: varid
  integer                            :: nf90_def_var_oneDim
end function nf90_def_var_oneDim

function nf90_def_var_ManyDims(ncid, name, xtype, dimids, varid)
  integer, intent(in) :: ncid
  character (len = *), intent(in) :: name
  integer, intent( in) :: xtype
  integer, dimension(:), intent(in) :: dimids
  integer, intent(out) :: varid
  integer                            :: nf90_def_var_ManyDims
end function nf90_def_var_ManyDims

end module netcdf
!
