!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  h y d r o _ m   (including several sub-modules)
!
!! Routing scheme for the basal meltwater.
!!
!! This module uses elements of the code produced by
!! Le Brocq et al. (2006, Comput. Geosci. 32, 1780-1795,
!! doi: 10.1016/j.cageo.2006.05.003), but has been modified
!! and repurposed for this hydrology module of SICOPOLIS.
!!
!!##### Authors
!!
!! Sebastian Beyer, Anne Le Brocq, Thomas Kleiner, Ralf Greve
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
!  In this file, the complete hydrology from libhydro has been merged into
!  one single file, and the precision_m module has been replaced by the
!  sico_types module of SICOPOLIS.
!
#define libhydro_VERSION_MAJOR 2
#define libhydro_VERSION_MINOR 3
#define libhydro_VERSION_PATCH 3
!
! - 2.3.1 by RG: subroutine hydro_get_vflux (return vector volume flux) added.
! - 2.3.2 by RG: Computation of vfluxX and vfluxY added in
!                subroutine flux_correction.
! - 2.3.3 by RG: All sub-modules consistently named hydro_xxx_m,
!                some comments added.
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Module hydro_fortfilt_m of libhydro.
!-------------------------------------------------------------------------------
module hydro_fortfilt_m

  use sico_types_m, only: dp
  use sico_variables_m, only: pi

  implicit none

  private ! everything private by default
  public :: naiveGauss
  public :: fastGauss
  public :: BoxBlur

contains

  subroutine naiveGauss (source, filtered, r)

    implicit none
    integer, intent(in)                  :: r
    real(dp), intent(in)                 :: source(:,:)
    real(dp), intent(out)                :: filtered(:,:)

    integer                              :: i, j, k, l
    integer                              :: ii, jj      ! inside the array
    integer                              :: nx, ny
    integer                              :: dim(2)
    real(dp)                             :: val, wsum
    real(dp)                             :: dsq         ! distance squared
    real(dp)                             :: wght        ! weight

    dim = shape(source)
    nx = dim(1)
    ny = dim(2)

    do i = 1, nx
       do j = 1, ny
          val = 0
          wsum = 0
          do k = i-r, i+r     ! inner loop over kernel
             do l = j-r, j+r  ! inner loop over kernel
                ii = min(nx, max(1,k))   ! make sure i is always inside the grid (this implies values are extendet (stretched at the boundaries))
                jj = min(ny, max(1,l))   ! make sure j is always inside the grid (this implies values are extendet (stretched at the boundaries))
                dsq = (j-l)**2 + (i-k)**2
                wght = exp(-dsq / (2.d0*r**2)) / (2.d0*pi*r**2)
                val = val + source(ii,jj) * wght
                wsum = wsum + wght
                ! print *, i,j, k, l, ii, jj, dsq
             end do
          end do
          filtered(i,j) = val / wsum
       end do
    end do

  end subroutine naiveGauss

  subroutine BoxBlurH (source, filtered, r)
    ! computes horizontal blur
    implicit none
    integer, intent(in)                  :: r
    real(dp), intent(in)                 :: source(:,:)
    real(dp), intent(out)                :: filtered(:,:)

    integer                              :: nx, ny
    integer                              :: dim(2)

    real(dp)                             :: wght  ! weight
    real(dp)                             :: sum
    integer                              :: i,j
    integer                              :: il    ! leftmost  pixel which should be removed from the accumulator
    integer                              :: ir    ! rightmost pixel which should be added   to   the accumulator

    dim = shape(source)
    nx = dim(1)
    ny = dim(2)

    wght = 1.d0 / (2.d0*r+1.d0)

    do j = 1, ny   ! loop over all rows
       ! compute sum at first pixel
       sum = source(1,j)
       do i = 1, r  ! loop over box kernel
          sum = sum + source(1,j) + source(i+1,j) ! always take 1 as left pixel, to not get out of grid
          ! print *, sum
       end do

       ! generate output pixel, then update running sum
       do i = 1, nx
          ! print *, j, i, sum
          filtered(i,j) = sum * wght
          il = max(i-r, 1)     ! make sure we dont get off the grid
          ir = min(i+r+1, nx)  ! make sure we dont get off the grid
          sum = sum + source(ir,j) - source(il,j)
       end do
    end do

  end subroutine BoxBlurH

  subroutine BoxBlur (source, filtered, r)
    ! computes box blur, by calling horizontal boxblur twice, the second time with tansposed matrix
    implicit none
    integer, intent(in)                  :: r
    real(dp), intent(in)                 :: source(:,:)
    real(dp), intent(out)                :: filtered(:,:)

    integer                              :: nx, ny
    integer                              :: dim(2)

    real(dp), allocatable                :: source_transp(:,:)   ! source transposed
    real(dp), allocatable                :: filtered_transp(:,:)   ! filtered transposed

    if (containsNaN (source)) then
      stop "Error: fortfilt: source contains NaNs"
    end if

    dim = shape(source)
    nx = dim(1)
    ny = dim(2)
    allocate(source_transp(ny,nx))
    allocate(filtered_transp(ny,nx))

    ! first horizontal blur
    call BoxBlurH (source, filtered, r)
    ! then transpose result and call horizontal blur again, now really doing vertical blur
    source_transp = transpose (filtered)
    call BoxBlurH (source_transp, filtered_transp, r)
    ! transpose again, to get back initial shape
    filtered = transpose (filtered_transp)

  end subroutine BoxBlur

  subroutine fastGauss (source, filtered, r)
    ! computes a fast approximation to a Gaussian filter, by applying a series of box filters
    implicit none
    integer, intent(in)                  :: r
    real(dp), intent(in)                 :: source(:,:)
    real(dp), intent(out)                :: filtered(:,:)

    integer                              :: nx, ny
    integer                              :: dim(2)

    real(dp), allocatable                :: tmpfilter1(:,:)   ! tmp array to store intermediate results
    real(dp), allocatable                :: tmpfilter2(:,:)   ! tmp array to store intermediate results

    dim = shape(source)
    nx = dim(1)
    ny = dim(2)
    allocate(tmpfilter1(nx,ny))
    allocate(tmpfilter2(nx,ny))

    call BoxBlur (source, tmpfilter1, r)
    call BoxBlur (tmpfilter1, tmpfilter2, r)
    call BoxBlur (tmpfilter2, filtered, r)

  end subroutine fastGauss

  pure function containsNaN (array) result(hasNaN)
    ! checks an array for NaN values
    implicit none
    real(dp), intent(in)     :: array(:,:)
    logical                  :: hasNaN

    if ( any(isnan(array)) ) then
      hasNaN = .true.
    else
      hasNaN = .false.
    end if
  end function containsNaN

end module hydro_fortfilt_m

!-------------------------------------------------------------------------------
!> Module hydro_priority_queue_m of libhydro.
!-------------------------------------------------------------------------------
module hydro_priority_queue_m

  use sico_types_m, only: dp

implicit none

private
public gridCell, queue

type, abstract :: node
end type node

type, extends (node) :: gridCell
  integer            :: i
  integer            :: j
  ! integer            :: z
  real(dp)           :: z
  integer            :: k ! additional integer to guarentee stability of heap
contains
  procedure :: compareTo => compareTo_gridCell
end type

type queue
  type (gridCell), allocatable :: heaplist(:)
  integer                      :: currentSize = 0  ! number of elements in queue
  integer                      :: count = 0        ! counter to increase k
contains
  procedure                    :: isEmpty
  procedure                    :: pop => pop_gridCell
  procedure                    :: push => push_gridCell
  procedure, private           :: siftDown
  procedure                    :: show
end type queue

contains

  pure integer function compareTo_gridCell (self, other) result (res)
    ! compares to grid cells taking into account the height z
    ! and if this leads to a tie then also k
    ! (the order in which the elements where inserted into the heap)
    ! A return value of -1 means the current cell preceds the other cell,
    ! a return value of +1 means the current cell follows the other cell
    class (gridCell), intent(in)     :: self
    class (gridCell), intent(in)     :: other

    if (self%z <= other%z .or. (self%z == other%z .and. self%k < other%k)) then
      res = 1
    else
      res = -1
    end if
  end function compareTo_gridCell

  pure logical function isEmpty (self) result (empty)
    class (queue), intent(in)     :: self

    if ( self%currentSize == 0 ) then
      empty = .true.
    else
      empty = .false.
    end if
  end function isEmpty

  function pop_gridCell (self) result (element)
    class (queue), intent(inout)        :: self
    type (gridCell)                  :: element
    element = self%heaplist(1)
    self%heaplist(1) = self%heaplist(self%currentSize)
    self%currentSize = self%currentSize - 1
    call self%siftDown (1)

  end function pop_gridCell

  integer function push_gridCell (self, element) result (success)
    class (queue), intent(inout)     :: self
    type (gridCell), intent(inout)      :: element
    type (gridCell), allocatable     :: tmp(:)
    integer                          :: i

    ! set k of gridCell to self%count to guarantee stability of heapsort
    ! then increase count
    element%k = self%count
    self%count = self%count + 1

    self%currentSize = self%currentSize + 1

    if (.not.allocated(self%heaplist)) allocate (self%heaplist(1))
    if (size(self%heaplist) < self%currentSize) then
      allocate(tmp(2*size(self%heaplist))) ! allocate whole new row
      tmp(1:self%currentSize-1) = self%heaplist
      call move_alloc(tmp, self%heaplist)
    end if
    self%heaplist(self%currentSize) = element
    i = self%currentSize
    do
      i = i / 2
      if (i==0) exit
      call self%siftDown (i)
      end do

      success = 1

  end function push_gridCell

  subroutine siftDown (self, a)
    class (queue), intent(inout)     :: self
    integer, intent(in)              :: a
    integer                          :: parent, child

     parent = a
     do while (parent*2 <= self%currentSize)
       child = parent*2
       if (child+1 <= self%currentSize) then
         if (self%heaplist(child+1)%compareTo (self%heaplist(child)) > 0) then
           child = child + 1
         end if
       end if
       if (self%heaplist(parent)%compareTo (self%heaplist(child)) < 0) then
         self%heaplist([child, parent]) = self%heaplist([parent, child])
         parent = child
       else
         exit
       end if
     end do
  end subroutine siftDown

  subroutine show (self)
    class (queue), intent(in)     :: self
    integer                       :: i
  do i=1,self%currentSize
    write (*,'(i4, i4, f7.2, i4)') self%heaplist(i)
  end do

  end subroutine

end module hydro_priority_queue_m

!-------------------------------------------------------------------------------
!> Module hydro_priority_flood_m of libhydro.
!-------------------------------------------------------------------------------
module hydro_priority_flood_m

  use sico_types_m, only: dp

  use hydro_priority_queue_m, only: gridCell, queue

implicit none

private
public pf_flood

  ! x offsets of D8 neighbors from a central cell, first 4 are cardinal
  integer, parameter, dimension(8)          :: noff_x = [-1, 0, 1, 0, -1, 1, 1,-1]
  integer, parameter, dimension(8)          :: noff_y = [ 0,-1, 0, 1, -1,-1, 1, 1]

contains

  pure logical function isOnGrid (nx, ny, x, y) result (res)
    integer, intent(in)                       :: nx, ny
    integer, intent(in)                       :: x, y

    if (x <= 0 .or. x > nx .or. y <= 0 .or. y > ny) then
      res = .false.
    else
      res = .true.
    end if
  end function isOnGrid

  subroutine pf_flood (dem, nx,ny, flatmask, tiltflats, nn_o)
    real(dp), dimension(:,:), intent(inout)       :: dem
    integer, intent(in)                           :: nx,ny
    integer, dimension(:,:), intent(out), optional:: flatmask  ! mask where sinks have been filled
    logical, intent(in), optional                 :: tiltflats      ! add a little increment to avoid 'really flat' areas
    integer, intent(in), optional                 :: nn_o          ! number of neighbors considered
    real(dp)                                      :: tiltinc        ! increment for tilting
    type (queue)                                  :: pqopen
    type (gridCell)                               :: cell
    type (gridCell)                               :: c            ! current cell
    logical, allocatable                          :: closed(:,:)
    integer                                       :: nn = 8            ! number of neighbors considered
    integer                                       :: processed_cells=0
    integer                                       :: pitc=0
    integer                                       :: i,j
    integer                                       :: neighbor_i, neighbor_j
    integer                                       :: n
    integer                                       :: success   ! success of push operation

    !------- configuration  ---------
    if (present(tiltflats)) then
      if (tiltflats) then
        ! print *, 'tilting flats is ON'
        tiltinc = 0.01
      else
        ! print *, 'tilting flats is OFF'
        tiltinc = 0
      end if
    end if

    if (present(nn_o)) then
       nn = nn_o
       if ((.not.(nn == 4)) .and. (.not. (nn == 8))) then
          print *, 'pflood: error: number of neighbors has to be 4 or 8 (was ', nn,')'
          stop 1
       end if
    end if
    ! print *, 'priority flood with', nn, ' neighbors'

    ! init flatmask
    if (present(flatmask)) then
      flatmask = 0
    end if

    !------- end configuration  ---------

    allocate(closed(nx,ny))
    closed = .false.

    ! adding edge cells to the priority queue
    do i=1,nx
      cell = gridCell(i,1,dem(i,1),1) ! bottom edge
      success = pqopen%push (cell)
      cell = gridCell(i,ny,dem(i,ny),1) ! top edge
      success = pqopen%push (cell)
      closed(i,1) = .true.
      closed(i,ny) = .true.
    end do
    do j=1,ny
      cell = gridCell(1,j,dem(1,j),1) ! left edge
      success = pqopen%push (cell)
      cell = gridCell(nx,j,dem(nx,j),1) ! right edge
      success = pqopen%push (cell)
      closed(1,j) = .true.
      closed(nx,j) = .true.
    end do

    ! starting the flood algorithm
    do while (.not. pqopen%isEmpty())
      c = pqopen%pop ()
      processed_cells = processed_cells + 1
      ! clever way to cycle through all neighbor cells
      do n=1,nn
        neighbor_i = c%i + noff_x(n)
        neighbor_j = c%j + noff_y(n)
        if (.not.isOnGrid (nx, ny, neighbor_i, neighbor_j)) cycle
        ! check if the cell has already been processed
        if (closed(neighbor_i,neighbor_j)) cycle

        closed(neighbor_i,neighbor_j) = .true.
        if (dem(neighbor_i,neighbor_j) < dem(c%i,c%j)) then
          pitc = pitc + 1
          if (present(flatmask)) flatmask(neighbor_i,neighbor_j) = 1
        end if
        dem(neighbor_i,neighbor_j) = max (dem(neighbor_i,neighbor_j), dem(c%i,c%j)+tiltinc)
        cell = gridCell(neighbor_i,neighbor_j,dem(neighbor_i,neighbor_j),1)
        success = pqopen%push (cell)
      end do
    end do

    ! print *, processed_cells, 'cells processed, ', pitc, 'in pits'
  end subroutine pf_flood

end module hydro_priority_flood_m

!-------------------------------------------------------------------------------
!> Module hydro_balance_flux_m of libhydro.
!-------------------------------------------------------------------------------
module hydro_balance_flux_m

  use sico_types_m, only: dp
  use sico_variables_m, only: pi

  implicit none

  integer, dimension(9), parameter :: dirpp_war = [ 9,8,7,6,5,4,3,2,1 ]
  integer, dimension(9), parameter :: dirpp_tar = [ 4,5,6,3,0,7,2,1,8 ]

  ! HERE BE DRAGONS !

contains

  recursive subroutine dparea_war(i, j, flow_num, nx, ny, x, y, &
       imask, srf, accum, flux)
    implicit none
    integer,intent(in) :: i !! ew_r (x-index,i)
    integer,intent(in) :: j !! ns_r (y-index,j)
    integer,intent(in) :: flow_num !! use 4 or 8 sides
    integer,intent(in) :: nx  !! number of grid points ewn
    integer,intent(in) :: ny  !! number of grid points nsn
    real(dp),dimension(nx),intent(in) :: x
    real(dp),dimension(ny),intent(in) :: y
    integer, dimension(nx,ny),intent(in) :: imask   !!
    real(dp),dimension(nx,ny),intent(in) :: srf     !! surface to compute the flow direction
    real(dp),dimension(nx,ny),intent(in) :: accum   !! accumulation/source
    real(dp),dimension(nx,ny),intent(inout) :: flux !! local computed flux

    real(dp), dimension(9) :: pp,delev
    real(dp), dimension(9) :: dd
    integer :: ii,jj,count,count2,iii,jjj, sumdcount!,check
    real(dp) :: sumd
    real(dp) :: dsq
    real(dp) :: dx,dy

    if(i==1) then
      dx = x(i+1) - x(i)
    else if (i==nx) then
      dx = x(i) - x(i-1)
    else
      dx = 0.5_dp * (x(i+1) - x(i-1))
    end if
    if(j==1) then
      dy = y(j+1) - y(j)
    else if (j==ny) then
      dy = y(j) - y(j-1)
    else
      dy = 0.5_dp * (y(j+1) - y(j-1))
    end if

    !check that it is possible to do the window on this cell...
    !otherwise, just assign it's accumulation and return..
    if (i > 1 .and. i < nx .and. j > 1 .and. j < ny &
        .and. imask(i,j) == 1) then ! org

 ! print *, "Start at ", i,j, flux(i,j)
      dsq = sqrt(dx*dx * dy*dy)
      dd = 1.0_dp / [ dsq, dy, dsq, dx, 1.0_dp, dx, dsq, dy, dsq ]
      dd(5) = 0.0_dp

      !TODO:
      ! change this for non equal spaced grids
      dd([1,3,7,9]) = dd([1,3,7,9]) * sqrt(2.0_dp) / 2.0_dp

      sumd = 0.0_dp

      !Check we haven't been here already... (flux is -1 if we haven't)
      !If we have, then can just return...
      if (flux(i,j) < 0.0_dp) then

        !modified 12/01/05 - flux is now in m3 - so consider it volume rather than flux....
        !        flux(i,j) = accum(i,j)*space*space*0.001    !TK why?
        flux(i,j) = accum(i,j)*dx*dy
        count = 1

        do jj=j-1,j+1
          do ii=i-1,i+1
            pp(:)=0.0_dp
            delev(:) = 0.0_dp
            if(count .ne. 5) then

              !work out proportion
              !first work out p (stored in pp()) for the neighbour
              if (ii > 1 .and. ii < nx .and. jj > 1 .and. jj < ny) then

                !calculate sum of elevation differences
                sumd = 0.0_dp
                count2 = 1

                do jjj=jj-1,jj+1
                  do iii=ii-1,ii+1

                    !4 directions...
                    if (flow_num == 4) then
                      if(count2 .ne. 5 .and. ((count2/2)*2) == count2) then
                        delev(count2) = max(0.0_dp,(srf(ii,jj) - srf(iii,jjj)))
                        sumd = sumd + delev(count2)
                        ! print *, ii,jj,count2,(srf(ii,jj) - srf(iii,jjj)), __LINE__
                      end if

                    else if (flow_num == 8) then
                      if(count2 .ne. 5) then
                        delev(count2) = max(0.0_dp,(srf(ii,jj) - srf(iii,jjj))) * dd(count2)
                        sumd = sumd + delev(count2)
                      end if

                    end if
                    count2 = count2 + 1

                  end do
                end do

                count2 = 1

                if (sumd == 0.0_dp .and. srf(ii,jj) > 0.0_dp) then
                  sumdcount = 0

                  do jjj=jj-1,jj+1
                    do iii=ii-1,ii+1

                      !4 flow directions
                      if (flow_num == 4) then
                        if(count2 .ne. 5 .and. ((count2/2)*2) == count2) then
                          if (srf(iii,jjj) == srf(ii,jj)) then
                            sumdcount = sumdcount + 1
                          end if
                        end if
                        !8 flow directions
                      else if (flow_num == 8) then
                        if(count2 .ne. 5) then
                          if (srf(iii,jjj) == srf(ii,jj)) then
                            sumdcount = sumdcount + 1
                          end if
                        end if
                      end if
                      count2 = count2 + 1

                    end do
                  end do

                  !work out proportions...

                  count2 = 1
                  do jjj=jj-1,jj+1
                    do iii=ii-1,ii+1
                      !4 flow directions
                      if (flow_num == 4) then
                        if(count2 .ne. 5 .and. ((count2/2)*2) == count2) then
                          if (srf(iii,jjj) == srf(ii,jj)) then
                            ! pp(count2) = 1.0_dp/sumdcount !ORG
                            pp(count2) = 1.0_dp/real(sumdcount,KIND=dp)
                          end if
                        end if
                        !8 flow directions
                      else if (flow_num == 8) then
                        if(count2 .ne. 5) then
                          if (srf(iii,jjj) == srf(ii,jj)) then
                            ! pp(count2) = 1.0_dp/sumdcount
                            pp(count2) = 1.0_dp/real(sumdcount,KIND=dp)
                          end if
                        end if
                      end if
                      count2 = count2 + 1
                    end do
                  end do

                  !this bit is the normal case - that there is at least one downslope neighbour
                else
                  do jjj=jj-1,jj+1
                    do iii=ii-1,ii+1

                      if (sumd <= 0.0_dp ) then
                        write (*,'("#WARN: sumd <= 0 at ",2(i4,1x),"neighbor >",2(i4,1x), "sumd =",es16.8)') ii,jj,iii,jjj,sumd
                      end if

                      if (flow_num == 4) then
                        if(count2 .ne. 5 .and. ((count2/2)*2) == count2) then
                          !pp(count2) = delev(count2)/sumd !ORG
                          if (sumd > 0.0_dp) pp(count2) = delev(count2)/sumd
                        end if
                      else if (flow_num == 8) then
                        if(count2 .ne. 5) then
                          !pp(count2) = delev(count2)/sumd !ORG
                          if (sumd > 0.0_dp) pp(count2) = delev(count2)/sumd
                        end if
                      end if

                      count2 = count2 + 1
                    end do
                  end do
                end if

                ! now we have info on the neighbour,
                ! now we need to check if it contributes to our cell
                if (pp(dirpp_war(count)) > 0.0_dp) then
                  call dparea_war(ii, jj, flow_num, nx, ny, x, y, &
                      imask, srf, accum, flux)
                  flux(i,j) = flux(i,j) + pp(dirpp_war(count)) * flux(ii,jj)
                end if

              end if
            end if
            count = count+1
          end do
        end do
      end if
      !for the outside cells?
    else
      flux(i,j) = accum(i,j)*dx*dy
    end if

  end subroutine dparea_war

  recursive subroutine dparea_tar(i, j, flow_num, nx, ny, x, y, &
      imask, srf, accum, flowdir, flux)
    !
    implicit none
    !
    integer,intent(in) :: i !! ew_r (x-index,i)
    integer,intent(in) :: j !! ns_r (y-index,j)
    integer,intent(in) :: flow_num !! use 4 or 8 sides ! ignored in tar
    integer,intent(in) :: nx  !! number of grid points ewn
    integer,intent(in) :: ny  !! number of grid points nsn
    !
    real(dp),dimension(nx),intent(in) :: x !! grid spacing x
    real(dp),dimension(ny),intent(in) :: y !! grid spacing y
    !
    integer, dimension(nx,ny),intent(in) :: imask   !!
    real(dp),dimension(nx,ny),intent(in) :: srf     !! surface to compute the flow direction
    real(dp),dimension(nx,ny),intent(in) :: accum   !! accumulation/source
    real(dp),dimension(nx,ny),intent(in) :: flowdir !! flowdir
    !
    real(dp),dimension(nx,ny),intent(inout) :: flux !! local computed flux
    ! --
    real(dp), dimension(9) :: pp
    integer :: ii,jj,count,sector
    real(dp) :: alpha2
    real(dp) :: alpha
    real(dp) :: dx,dy
    real(dp), parameter :: pi4 = pi/4.0_dp
    real(dp), parameter :: pi2 = pi/2.0_dp
    ! ---

    !check that it is possible to do the window on this cell...
    !otherwise, just assign it's accumulation and return..

    if(i==1) then
      dx = x(i+1) - x(i)
    else if (i==nx) then
      dx = x(i) - x(i-1)
    else
      dx = 0.5_dp * (x(i+1) - x(i-1))
    end if
    if(j==1) then
      dy = y(j+1) - y(j)
    else if (j==ny) then
      dy = y(j) - y(j-1)
    else
      dy = 0.5_dp * (y(j+1) - y(j-1))
    end if

    if (i > 1 .and. i < nx .and. j > 1 .and. j < ny .and. imask(i,j) > 0) then

      !Check we haven't been here already... (flux is -1 if we haven't)
      !If we have, then can just return...
      if (flux(i,j) < 0.0_dp) then

        ! flux(i,j) = accum(i,j)*space*space*0.001
        flux(i,j) = accum(i,j)*dx*dy

        count = 1
        do jj=j-1,j+1
          do ii=i-1,i+1

            pp(:) = 0.0_dp

            ! flowdir < 0 is noflow
            if(count .ne. 5 .and. flowdir(ii,jj) >= 0.0_dp) then

              !work out proportion
              !first work out p (stored in pp()) for the neighbour
              !work out which sector (triangular) the flow falls in

              ! changed tkleiner
              ! changed compared to TARBOTON because of different ordering
              alpha = flowdir(ii,jj) + pi2
              if(alpha > 2.0_dp *pi) alpha = alpha - 2.0_dp * pi

              ! check this TK because we now have 0 2Pi
              ! sector = int((flowdir(ii,jj) + 45.0)/45.0)
              sector = int((alpha + pi4)/pi4)
              if (sector == 9) then
                sector = 1
              end if

              !special case for sector 8, as there is no sector + 1
              !now we know which sector it falls in, we need to work out the proportion
              !that goes in each cell.

              ! check this TK because we now have 0 2Pi
              ! alpha2 = (sector * 45.0) - flowdir(ii,jj)
              alpha2 = (real(sector,dp) * pi4) - alpha
              ! FIX
              if(alpha2 < 0.0_dp) alpha2 = alpha2 + 2.0_dp *pi

              if (sector .ne. 8) then

                ! pp(sector) = alpha2/45.0 ! check this TK because we now have 0 2Pi
                pp(sector) = alpha2/pi4
                pp(sector + 1) = 1.0_dp - pp(sector)

              else

                ! pp(sector) = alpha2/45.0  ! check this TK because we now have 0 2Pi
                pp(sector) = alpha2/pi4  !
                pp(1) = 1.0_dp - pp(sector)

              end if

              ! now we have info on the neighbour,
              ! now we need to check if it contributes to our cell

              if (pp(dirpp_tar(count)) > 0.0_dp ) then

                call dparea_tar(ii, jj, flow_num, nx, ny, x, y, &
                    imask, srf, accum, flowdir, flux)
                flux(i,j) = flux(i,j) + pp(dirpp_tar(count)) * flux(ii,jj)

              end if
            end if
            count = count+1
          end do
        end do

      else
        ! visited earlier
      end if
    else
      flux(i,j) = accum(i,j)*dx*dy
    end if
  end subroutine dparea_tar

  subroutine flow_dir(nx, ny, x, y, imask, srf, srf_x, srf_y, theta)
    implicit none
    integer,intent(in) :: nx  !! number of grid points ewn
    integer,intent(in) :: ny  !! number of grid points nsn
    real(dp),dimension(nx),intent(in) :: x !! grid spacing x
    real(dp),dimension(ny),intent(in) :: y !! grid spacing y
    integer,dimension(nx,ny),intent(in) :: imask     !!
    real(dp),dimension(nx,ny),intent(in) :: srf !! surface to compute the flow direction
    real(dp),dimension(nx,ny),intent(out) :: srf_x     !!
    real(dp),dimension(nx,ny),intent(out) :: srf_y     !!
    real(dp),dimension(nx,ny),intent(out) :: theta     !! flow direction

    integer :: i,j,jstart,jend
    srf_x(:,:) = 0.0_dp
    srf_y(:,:) = 0.0_dp
    theta(:,:) = -9999.0_dp

    jstart = 2
    jend = ny-1

    do j=jstart,jend
      do i = 2,nx-1

        if(imask(i,j) > 0) then

          ! equation (7) in Budd & Warner 1996 (better use CDS from fds.F90)
          srf_x(i,j) = (srf(i+1,j) - srf(i-1,j)) / (x(i+1)-x(i-1))
          srf_y(i,j) = (srf(i,j+1) - srf(i,j-1)) / (y(j+1)-y(j-1))

          ! direction of the slope eq. (9) sin (theta) / cos (theta)
          ! theta(i,j) = atan2(srf_y(i,j),srf_x(i,j)) ! * 180.0_dp/Pi
          theta(i,j) = atan2(-srf_y(i,j),-srf_x(i,j)) ! * 180.0_dp/Pi

          ! atan2 produces results in the range (−Pi, Pi], which can be
          ! mapped to [0, 2Pi) by adding 2Pi to negative values
          if (theta(i,j) < 0.0_dp) then
            theta(i,j) = 2.0_dp * pi + theta(i,j)
          end if
        end if
      end do
    end do
  end subroutine flow_dir

  subroutine tar_flow_dir(nx, ny, x, y, imask, srf, srf_x, srf_y, theta)
    implicit none
    integer,intent(in) :: nx  !! number of grid points ewn
    integer,intent(in) :: ny  !! number of grid points nsn
    real(dp),dimension(nx),intent(in) :: x !! grid spacing x
    real(dp),dimension(ny),intent(in) :: y !! grid spacing y
    integer,dimension(nx,ny),intent(in) :: imask  !!
    real(dp),dimension(nx,ny),intent(in) :: srf    !! surface to compute the flow direction

    real(dp),dimension(nx,ny),intent(out) :: srf_x     !!
    real(dp),dimension(nx,ny),intent(out) :: srf_y     !!
    real(dp),dimension(nx,ny),intent(out) :: theta     !! flow direction

    integer :: i,j
    integer :: max_pos,ct
    real(dp) :: a1, rg, max_val
    real(dp) :: r(8), s(8), dx,dy

    ! af in table 6 Tarboton 1997
    real(dp), parameter :: af(8) = &
        [ 1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp ]

    ! ac in table 6 Tarboton 1997
    real(dp), parameter :: ac(8) = &
        [ 0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 4.0_dp ]

    ! the main body of the calcs. finding the direction of the vector
    ! for each segment and then the angles and then the max vector
    do j = 2,ny-1
      do i = 2,nx-1

        ! equation (7) in Budd & Warner 1996 ! not used here
        srf_x(i,j) = (srf(i+1,j) - srf(i-1,j)) / (x(i+1)-x(i-1))
        srf_y(i,j) = (srf(i,j+1) - srf(i,j-1)) / (y(j+1)-y(j-1))

        if ( imask(i,j) > 0 ) then

          dx = 0.5_dp*(x(i+1)-x(i-1))
          dy = 0.5_dp*(y(j+1)-y(j-1))

          r(:) = 0.0_dp; s(:) = 0.0_dp

          ! table 1 in Tarboton 1997 (index of i,j has different meaning)
          call tar_rs(srf(i,j),srf(i+1,j),srf(i+1,j-1),dx,dy,r(1),s(1))
          call tar_rs(srf(i,j),srf(i,j-1),srf(i+1,j-1),dx,dy,r(2),s(2))
          call tar_rs(srf(i,j),srf(i,j-1),srf(i-1,j-1),dx,dy,r(3),s(3))
          call tar_rs(srf(i,j),srf(i-1,j),srf(i-1,j-1),dx,dy,r(4),s(4))
          call tar_rs(srf(i,j),srf(i-1,j),srf(i-1,j+1),dx,dy,r(5),s(5))
          call tar_rs(srf(i,j),srf(i,j+1),srf(i-1,j+1),dx,dy,r(6),s(6))
          call tar_rs(srf(i,j),srf(i,j+1),srf(i+1,j+1),dx,dy,r(7),s(7))
          call tar_rs(srf(i,j),srf(i+1,j),srf(i+1,j+1),dx,dy,r(8),s(8))

          ! find the maximum downslop vector

          max_val = 0.0_dp

          do ct = 1,8
            if (s(ct) .gt. max_val) then
              max_pos = ct
              max_val = s(ct)
            end if
          end do

          ! in this section the i and j coords of the twp cells
          ! into which the flow goes are recorded so that they tally
          ! with a1 and a2.

          a1 = r(max_pos)
          rg = r(max_pos) * af(max_pos) + ac(max_pos) * pi / 2.0_dp ! eq. 6 in Tarboton 1997

          if (max_val .lt. 0.0_dp) then
            rg = -9999 ! noflow
          else
            ! tkleiner adjust to the different i,j orientation
            if (rg > 0.0_dp)  rg = 2.0_dp*pi - rg
          end if

          ! rg is in radians and anti-clockwise from the East.
          theta(i,j) = rg

        end if
      end do

    end do

  end subroutine tar_flow_dir

  subroutine tar_rs(e0,e1,e2,d1,d2,r,s)
    !
    implicit none

    real(dp), intent(in) :: e0,e1,e2
    real(dp), intent(in) :: d1,d2 !
    real(dp), intent(out) :: r,s

    real(dp) :: s1, s2
    real(dp), save :: r_lim, con
    logical, save :: first = .true.

    if (first) then
      r_lim = atan2(d2,d1)      ! see Tarboton 1997 eq. 5
      con = sqrt(d1**2 + d2**2) ! see Tarboton 1997 eq. 5
      first = .false.
    end if

    s1 = (e0-e1) / d1 ! Tarboton 1997 eq. 1
    s2 = (e1-e2) / d2 ! Tarboton 1997 eq. 2

    r = atan2(s2,s1)        ! r in Tarboton 1997 eq. 3
    s = sqrt(s1**2 + s2**2) ! s in Tarboton 1997 eq. 3

    !
    ! r is not in the range of (0,tan^-1(d2,d1)
    ! Tarboton 1997 eq. 4 and 5
    if (r < 0.0_dp) then
      r = 0.0_dp
      s = s1
    else if (r > r_lim) then
      r = r_lim
      s = (e0-e2) / con
    end if

    return

  end subroutine tar_rs

end module hydro_balance_flux_m

!-------------------------------------------------------------------------------
!> Module hydro_modelstate_flux_m of libhydro.
!-------------------------------------------------------------------------------
module hydro_modelstate_flux_m

  use sico_types_m, only: dp

  implicit none

  type :: modelstate_flux_t
    real(dp), dimension(:,:), allocatable :: vflux  !! vector flux
    real(dp), dimension(:,:), allocatable :: sflux  !! scalar flux
    real(dp), dimension(:,:), allocatable :: vfluxX
    real(dp), dimension(:,:), allocatable :: vfluxY
    real(dp), dimension(:,:), allocatable :: dir    !! direction of flux
    real(dp), dimension(:,:), allocatable :: bwat   !! water layer thickness
  contains
    procedure :: alloc => allocate_flux
    procedure :: dealloc => deallocate_flux
    procedure :: init => init_flux
  end type modelstate_flux_t

contains

  subroutine allocate_flux (self, nx, ny)
    class (modelstate_flux_t), intent(inout)  :: self
    integer, intent(in)                       :: nx, ny
    if (allocated(self%vflux))    deallocate(self%vflux)
    if (allocated(self%sflux))    deallocate(self%sflux)
    if (allocated(self%vfluxX))   deallocate(self%vfluxX)
    if (allocated(self%vfluxY))   deallocate(self%vfluxY)
    if (allocated(self%dir))      deallocate(self%dir)
    if (allocated(self%bwat))     deallocate(self%bwat)

    allocate(self%vflux(nx, ny))
    allocate(self%sflux(nx, ny))
    allocate(self%vfluxX(nx, ny))
    allocate(self%vfluxY(nx, ny))
    allocate(self%dir(nx, ny))
    allocate(self%bwat(nx, ny))
  end subroutine allocate_flux
  
  subroutine deallocate_flux (self)
    class (modelstate_flux_t), intent(inout)  :: self
    
    deallocate(self%vflux)
    deallocate(self%sflux)
    deallocate(self%vfluxX)
    deallocate(self%vfluxY)
    deallocate(self%dir)
    deallocate(self%bwat)
  end subroutine deallocate_flux
  
  subroutine init_flux (self)
    class (modelstate_flux_t), intent(inout)  :: self
    
    self%vflux = -9999
    self%sflux = -9999
    self%vfluxX = -9999
    self%vfluxY = -9999
    self%dir = -9999
    self%bwat = -9999
  end subroutine init_flux

end module hydro_modelstate_flux_m

!-------------------------------------------------------------------------------
!> Module hydro_modelstate_geom_m of libhydro.
!-------------------------------------------------------------------------------
module hydro_modelstate_geom_m

  use sico_types_m, only: dp

  implicit none

  type :: modelstate_geom_t
    integer                               :: nx
    integer                               :: ny
    real(dp)                              :: dx
    real(dp)                              :: dy
    real(dp), dimension(:),   allocatable :: x  ! vector containing x values
    real(dp), dimension(:),   allocatable :: y  ! vector containing y values
    integer,  dimension(:,:), allocatable :: icemask
    real(dp), dimension(:,:), allocatable :: usurf
    real(dp), dimension(:,:), allocatable :: thk
    real(dp), dimension(:,:), allocatable :: topg
    real(dp), dimension(:,:), allocatable :: bmelt
    real(dp), dimension(:,:), allocatable :: temppabase
  contains
    procedure :: alloc => allocate_geom
    procedure :: dealloc => deallocate_geom
    procedure :: init => init_geom
  end type modelstate_geom_t

contains

  subroutine allocate_geom (self, nx, ny)
    class (modelstate_geom_t), intent(inout)  :: self
    integer, intent(in)                       :: nx, ny
    if (allocated(self%x))          deallocate(self%x)
    if (allocated(self%y))          deallocate(self%y)
    if (allocated(self%icemask))    deallocate(self%icemask)
    if (allocated(self%usurf))      deallocate(self%usurf)
    if (allocated(self%thk))        deallocate(self%thk)
    if (allocated(self%topg))       deallocate(self%topg)
    if (allocated(self%bmelt))      deallocate(self%bmelt)
    if (allocated(self%temppabase)) deallocate(self%temppabase)

    allocate(self%x(nx))
    allocate(self%y(ny))
    allocate(self%icemask(nx, ny))
    allocate(self%usurf(nx, ny))
    allocate(self%thk(nx, ny))
    allocate(self%topg(nx, ny))
    allocate(self%bmelt(nx, ny))
    allocate(self%temppabase(nx, ny))
  end subroutine allocate_geom

  subroutine deallocate_geom (self)
    class (modelstate_geom_t), intent(inout)  :: self
    
    deallocate(self%x)
    deallocate(self%y)
    deallocate(self%icemask)
    deallocate(self%usurf)
    deallocate(self%thk)
    deallocate(self%topg)
    deallocate(self%bmelt)
    deallocate(self%temppabase)
  end subroutine deallocate_geom
  
  subroutine init_geom (self)
    !! initialize geometry
    !! set most things to zero or otherwise easy to recognize
    !! 'magic' numbers
    class (modelstate_geom_t), intent(inout)  :: self
    
    self%x = -9999
    self%y = -9999
    self%icemask = 1
    self%usurf = -9999
    self%thk = -9999
    self%topg = -9999
    self%bmelt = -9999
    self%temppabase = -9999
  end subroutine init_geom

end module hydro_modelstate_geom_m

!-------------------------------------------------------------------------------
!> Module hydro_modelstate_pot_m of libhydro.
!-------------------------------------------------------------------------------
module hydro_modelstate_pot_m

  use sico_types_m, only: dp

  implicit none

  type :: modelstate_pot_t
    real(dp), dimension(:,:), allocatable :: pot
    real(dp), dimension(:,:), allocatable :: pot_smooth
    real(dp), dimension(:,:), allocatable :: pot_filled
    real(dp), dimension(:,:), allocatable :: pot_diff
                 !! potential difference (amount filled, e.g. 'lake volume')
  contains
    procedure :: alloc   => allocate_pot
    procedure :: dealloc => deallocate_pot
    procedure :: init    => init_pot
  end type modelstate_pot_t

contains

  subroutine allocate_pot (self, nx, ny)
    class (modelstate_pot_t), intent(inout)  :: self
    integer, intent(in)                      :: nx, ny
    if (allocated(self%pot))         deallocate(self%pot)
    if (allocated(self%pot_smooth))  deallocate(self%pot_smooth)
    if (allocated(self%pot_filled))  deallocate(self%pot_filled)
    if (allocated(self%pot_diff))    deallocate(self%pot_diff)

    allocate(self%pot(nx, ny))
    allocate(self%pot_smooth(nx, ny))
    allocate(self%pot_filled(nx, ny))
    allocate(self%pot_diff(nx, ny))
  end subroutine allocate_pot

  subroutine deallocate_pot (self)
    class (modelstate_pot_t), intent(inout)  :: self

    deallocate(self%pot)
    deallocate(self%pot_smooth)
    deallocate(self%pot_filled)
    deallocate(self%pot_diff)
  end subroutine deallocate_pot

  subroutine init_pot (self)
    class (modelstate_pot_t), intent(inout)  :: self

    self%pot = -9999
    self%pot_smooth = -9999
    self%pot_filled = -9999
    self%pot_diff = -9999
  end subroutine init_pot

end module hydro_modelstate_pot_m

!-------------------------------------------------------------------------------
!> Module hydro_modelstate_m of libhydro.
!-------------------------------------------------------------------------------
module hydro_modelstate_m

  use sico_types_m, only: dp

  use hydro_modelstate_geom_m, only: modelstate_geom_t
  use hydro_modelstate_pot_m, only: modelstate_pot_t
  use hydro_modelstate_flux_m, only: modelstate_flux_t

  implicit none

  type :: hydro_modelstate_t
     type(modelstate_geom_t)        :: geometry
     type(modelstate_pot_t)         :: potential
     type(modelstate_flux_t)        :: flux
  contains
     procedure :: alloc => hydro_modelstate_alloc
     procedure :: dealloc => hydro_modelstate_dealloc
     procedure :: init => hydro_modelstate_init
  end type hydro_modelstate_t

contains

  subroutine hydro_modelstate_alloc (self, nx, ny)
    class (hydro_modelstate_t), intent(inout)    :: self
    integer, intent(in)                          :: nx, ny

    call self%geometry%alloc(nx, ny)
    call self%potential%alloc(nx, ny)
    call self%flux%alloc(nx, ny)
  end subroutine hydro_modelstate_alloc
  
  subroutine hydro_modelstate_dealloc (self)
    class (hydro_modelstate_t), intent(inout)    :: self

    call self%geometry%dealloc ()
    call self%potential%dealloc ()
    call self%flux%dealloc ()
  end subroutine hydro_modelstate_dealloc
  
  subroutine hydro_modelstate_init (self)
    class (hydro_modelstate_t), intent(inout)    :: self

    call self%geometry%init ()
    call self%potential%init ()
    call self%flux%init ()
  end subroutine hydro_modelstate_init
  
end module hydro_modelstate_m

!-------------------------------------------------------------------------------
!> Module hydro_conf_m of libhydro.
!-------------------------------------------------------------------------------
module hydro_conf_m

  use sico_types_m, only: dp

  type hydro_conf_t
     !! Configuration type for hydrive
     !! Contains all necessary values to run the hydro model
     character(256)     :: infile     !! input  NetCDF file
     character(256)     :: outfile    !! output NetCDF file
     character(128)     :: method     !! fluxroutine method (warner, quinn, tarboton)
     character(32)      :: version     !! libhydro version
     character(128)     :: flatresolution !! how flat areas in hyd. pot. should be resolved
     real(dp)           :: global_supply !! water supply globally set from cmd line [m/a]
     logical            :: read_bmelt
     !! if basal melt should be read from netcdf file, this should depend on whether global supply is set
     character(256)     :: ncFormat   !! NetCDF file format of input (test, pism sico, timfd3)
     logical            :: avoid_frz
     !! If frozen areas (temp rel. to pressure melting point below 0) should be avoided
     real(dp)           :: filter_len !! filter length for hyd. potential [m]
     real(dp)           :: rho_freshwater !! used for hydraulic potential
     real(dp)           :: rho_seawater !! used for icemask calculation
     real(dp)           :: rho_ice
     real(dp)           :: gravity_const
     real(dp)           :: min_thickness !! minimum ice thickness to be considered in icemask
     logical            :: comp_mask !! whether the ice mask should be computed in libhydro
   contains
     procedure :: init => init_config
     procedure :: print => print_config
  end type hydro_conf_t

contains

  subroutine init_config(self)
    !! Init config type with default data
    class (hydro_conf_t), intent(out)   :: self

    self%infile = 'unknown'
    self%outfile = 'hydro_out.nc'
    self%method = 'quinn'
    self%version = 'none'
    self%flatresolution = 'fast'
    self%global_supply = 0.0_dp
    self%read_bmelt = .true.
    self%ncFormat = 'test'
    self%avoid_frz = .false.
    self%filter_len = 0.0_dp
    self%rho_freshwater = 1000.0_dp
    self%rho_seawater = 1025.0_dp
    self%rho_ice = 917.0_dp
    self%gravity_const = 9.81_dp
    self%min_thickness = 1.0_dp
    self%comp_mask = .true.
  end subroutine init_config

  subroutine print_config(self)
    class (hydro_conf_t), intent(in)   :: self

    print *, "--- CONFIGURATION:"
    print '(a30, a20)', "Version              ", trim(self%version)
    print '(a30, a20)', "Infile               ", trim(self%infile)
    print '(a30, a20)', "Outfile              ", trim(self%outfile)
    print '(a30, a20)', "NetCDF fmt           ", trim(self%ncFormat)
    print '(a30, a20)', "Fluxrouting method   ", trim(self%method)
    print '(a30, f8.5, a10)',"Global Supply   ", (self%global_supply), "m/year WE"
    print '(a30, l2)',  "read bmelt           ", (self%read_bmelt)
    print '(a30, l2)',  "avoid frozen         ", (self%avoid_frz)
    print '(a30, f8.1, a10)',"Filtering       ", (self%filter_len), " m"
    print '(a30, a20)', "Flat reso:           ", trim(self%flatresolution)
    print '(a30, f8.1, a10)',"Rho freshwater       ", (self%rho_freshwater), " "
    print '(a30, f8.1, a10)',"Rho freshwater       ", (self%rho_seawater), " "
    print '(a30, f8.1, a10)',"Rho ice         ", (self%rho_ice), " "
    print '(a30, f7.2, a10)',"Gravity g       ", (self%gravity_const), " "
    print '(a30, f7.2, a10)',"Min thickness   ", (self%min_thickness), "m"
    print '(a30, l2)',"Compute mask    ", (self%comp_mask), " "

  end subroutine print_config

end module hydro_conf_m

!-------------------------------------------------------------------------------
!> Module hydro_barnes_wrapper_m of libhydro.
!-------------------------------------------------------------------------------
module hydro_barnes_wrapper_m

  use sico_types_m, only: sp, dp

  use hydro_modelstate_pot_m, only: modelstate_pot_t
  use iso_fortran_env, only: error_unit

  implicit none

  private
  public :: barnes_resolve_flats

contains

  subroutine barnes_resolve_flats(nx, ny, potential, nneighbors)
    class(modelstate_pot_t), intent(inout) :: potential  !! hydraulic potential where flats need to be handled
    integer, intent(in)                    :: nneighbors !! number of neighbor cells considered, 4 or 8
    integer, intent(in)                    :: nx, ny
    real(sp), dimension(:,:), allocatable  :: pot_singleprec !! barnes only works in single precision mode, so this holds the values

    write (error_unit, '(a)') "ERROR: libhydro has been compiled without barnes flat resolution!"
    write (error_unit, '(a)') " Don't use that! (Check hydro_gen_conf)"
    stop 1
  end subroutine

end module hydro_barnes_wrapper_m

!-------------------------------------------------------------------------------
!> Main module hydro_m of libhydro: Routing scheme for the basal meltwater.
!-------------------------------------------------------------------------------
module hydro_m

  use sico_types_m, only: dp

  use hydro_modelstate_m, only: hydro_modelstate_t
  use hydro_modelstate_geom_m, only: modelstate_geom_t
  use hydro_modelstate_pot_m, only: modelstate_pot_t
  use hydro_modelstate_flux_m, only: modelstate_flux_t
  use hydro_fortfilt_m, only: BoxBlur
  use hydro_priority_flood_m, only: pf_flood
  use hydro_balance_flux_m, only: dparea_war, dparea_tar, flow_dir, tar_flow_dir
  use hydro_conf_m, only: hydro_conf_t
  use hydro_barnes_wrapper_m, only: barnes_resolve_flats

  implicit none

  private
  ! public api
  public :: hydro_t
  public :: hydro_gen_conf
  public :: hydro_init
  public :: hydro_update
  public :: hydro_get_sflux
  public :: hydro_get_vflux
  public :: hydro_get_bwat
  public :: hydro_set_mask
  public :: hydro_set_supply
  public :: hydro_set_temppabase
  public :: hydro_set_thk
  public :: hydro_set_topg
  ! advanced api
  public :: smooth_potential
  public :: hydro_potential
  public :: fill_sinks
  public :: resolve_flats
  public :: generateIceMask
  public :: balance_flux
  public :: flux_correction
  public :: flux2bwat
  public :: rise_potential_where_freezing
  public :: libhydro_version

  interface smooth_potential
     module procedure smooth_potential_meters, &
                      smooth_potential_pixel
  end interface smooth_potential

  type hydro_t
     !! hydro model type used in the public api
     type (hydro_modelstate_t) :: state
     type (hydro_conf_t)       :: conf
  end type hydro_t

contains

  subroutine hydro_potential (potential, geometry, RHO_FRESHWATER, RHO_ICE, G)
    !! Calculate hydraulic potential $\Phi$ from bedrock elevation
    !! $z_b$ and ice thickness $H$.
    !! $$ \Phi = \rho_{freshwater} g z_b + \rho_{ice} g H $$
    !! Ice surface elevation $z_s$ is not used here.
    !! @todo Some consistency checking here would be nice!
    class(modelstate_pot_t), intent(inout)    :: potential
    class(modelstate_geom_t), intent(in)      :: geometry
    real(dp), intent(in)                      :: RHO_FRESHWATER
    real(dp), intent(in)                      :: RHO_ICE
    real(dp), intent(in)                      :: G

    potential%pot = RHO_FRESHWATER * G * geometry%topg + RHO_ICE &
         & * G * geometry%thk
  end subroutine hydro_potential

  subroutine smooth_potential_meters (potential, fltlen_meters, dx)
    !! Smooth the hydraulic potential
    !! Filter radius is given in meters
    !! This converts the radius to grid pixel and then calls
    !! smooth_potential_pixel
    class(modelstate_pot_t), intent(inout)     :: potential
    real(dp), intent(in)                       :: fltlen_meters
    !! Radius of the filter window in meters
    real(dp), intent(in)                       :: dx
    integer                                    :: fltlen_pixel

    fltlen_pixel = nint (fltlen_meters / dx)
    call smooth_potential_pixel (potential, fltlen_pixel)

  end subroutine smooth_potential_meters

  subroutine smooth_potential_pixel (potential, fltlen_pixel)
    !! Smooth the hydraulic potential
    !! Filter radius is given in pixels
    class(modelstate_pot_t), intent(inout)     :: potential
    integer, intent(in)                        :: fltlen_pixel
    !! Radius of the filter window in pixel

    if (fltlen_pixel > 0) then
       call BoxBlur (potential%pot, potential%pot_smooth, fltlen_pixel)
    else
       ! copy the unchanged data,
       ! because later subroutines expect them to be there
       potential%pot_smooth = potential%pot
    end if

  end subroutine smooth_potential_pixel

  subroutine fill_sinks (potential, routingmethod, tiltflats)
    !! Fill sinks in the hydraulic potential via
    !! priority flooding.
    class(modelstate_pot_t), intent(inout)     :: potential  !! hydraulic potential to be filled
    character(len=*), intent(in)               :: routingmethod !! warner, quinn, tarboton
    logical, intent(in), optional              :: tiltflats  !! add a small increment while filling flats to make sure they drain
    logical                                    :: tiltflats_
    integer                                    :: nneighbors !! number of neighbor cells considered, 4 or 8
    integer                                    :: dim(2)
    integer                                    :: nx, ny

    dim = shape(potential%pot)
    nx = dim(1)
    ny = dim(2)

    nneighbors = method2nneighbors (routingmethod)

    if (present(tiltflats)) then
       tiltflats_ = tiltflats
    end if

    ! priority flood works in-place, so copy the data first
    potential%pot_filled = potential%pot_smooth
    call pf_flood (potential%pot_filled, nx, ny, tiltflats=tiltflats_, nn_o=nneighbors)

    ! compute potential difference (amount filled)
    potential%pot_diff = potential%pot_filled - potential%pot_smooth
  end subroutine fill_sinks

  subroutine resolve_flats_random (potential, nneighbors)
    !! resolve flat areas ($\Phi = 0$) by adding a small
    !! random increment a all points in a flat area.
    !! @warning This does not guarantee drainage and may lead
    !! to loss of water! It is advised to use the Barnes method
    !! or the tilting method instead.
    class(modelstate_pot_t), intent(inout)     :: potential  !! hydraulic potential to be filled
    integer, intent(in)                        :: nneighbors !! number of neighbor cells considered, 4 or 8
    integer                                    :: dim(2)
    integer                                    :: nx, ny
    integer :: cnt,i,j,l
    real(dp) :: grad_x,grad_y,grad,rand01,old

    dim = shape(potential%pot)
    nx = dim(1)
    ny = dim(2)

    outer: do l = 1,10
    cnt = 0
    ! print *, "loop: ",l
    do j = 2,ny-1
       do i = 2,nx-1
          grad_x = ( potential%pot_filled(i+1,j) - potential%pot_filled(i-1,j) )
          grad_y = ( potential%pot_filled(i,j+1) - potential%pot_filled(i,j-1) )
             grad = sqrt(grad_x**2+grad_y**2)

             if (grad <= 0.0_dp) then
                cnt = cnt + 1
                if (grad_x == 0.0_dp) then
                   old = potential%pot_filled(i+1,j)
                   call random_number(rand01)
                   potential%pot_filled(i+1,j) = potential%pot_filled(i+1,j) + (rand01 - 0.5_dp)*1e-3_dp ! +/- 5*10^-7
                end if

                if (grad_y == 0.0_dp) then
                   old = potential%pot_filled(i,j+1)
                   call random_number(rand01)
                   potential%pot_filled(i,j+1) = potential%pot_filled(i,j+1) + (rand01 - 0.5_dp)*1e-3_dp ! +/- 5*10^-7
                end if
             end if
       end do
    end do

    if (cnt > 0 ) then
       write (*,'("# ",i6," points altered (through adding a small random increment) to avoid grad Phi = 0")') cnt
    else
       exit outer
    end if
 end do outer

  end subroutine resolve_flats_random

  subroutine resolve_flats (potential, routingmethod, flatmethod)
    !! Deal with flat areas ($\Phi = 0$) in a sane way.
    !! Options are:
    !! 'random' --- a small random amount is added
    !! to each flat point.
    !!
    !! @Warning
    !! This is not recommended, because flux is not conserved this way!
    !!
    !! 'barnes' --- use the Method of Barnes 2014 to route the
    !! water through a flat area towards the point of lowest elevation
    !! at the boundary
    class(modelstate_pot_t), intent(inout) :: potential  !! hydraulic potential to be filled
    integer                                :: nneighbors !! number of neighbor cells considered, 4 or 8
    character(len=*), intent(in)           :: routingmethod !! warner, quinn, tarboton
    character(len=*), intent(in)           :: flatmethod !! how flat surfaces should be treated to avoid phi=0
    integer                                :: dim(2)
    integer                                :: nx, ny

    dim = shape(potential%pot)
    nx = dim(1)
    ny = dim(2)

    nneighbors = method2nneighbors (routingmethod)

    select case (flatmethod)
    case ('random')
       print *, "WARNING: Random flat resolution is not recommended!"
       print *, "It is only implemented for test/comparison reasons"
       call resolve_flats_random (potential, nneighbors)
    case ('barnes')
        call barnes_resolve_flats (nx, ny, potential, nneighbors)
    case default
       print *, "Method for flat resolution unknown:", flatmethod
       stop 1
    end select
  end subroutine resolve_flats

  subroutine generateIceMask (geometry, RHO_ICE, RHO_SEAWATER, min_thickness)
    !! Generate an ice mask according to ice thickness
    !! and equilibrium condition
    class (modelstate_geom_t), intent(inout)    :: geometry
    real(dp), intent(in)                        :: RHO_ICE
    real(dp), intent(in)                        :: RHO_SEAWATER  !! seawater density
    real(dp), intent(in)                        :: min_thickness !! minimum ice thickness to be considered
    real(dp), dimension(:,:), allocatable       :: H_eq          !! equilibrium height

    allocate(H_eq(geometry%nx, geometry%ny))
    H_eq = ( 1.0_dp - RHO_ICE/RHO_SEAWATER ) * geometry%thk
    geometry%icemask = 0
    where ( geometry%topg + geometry%thk > H_eq ) geometry%icemask = 1
    where ( geometry%thk < min_thickness) geometry%icemask = 0
    deallocate(H_eq)
  end subroutine generateIceMask

  subroutine ensure_positive_supply(geometry)
    !! Set negative supply to zero instead.
    !! For fluxrouting the supply needs to be positive,
    !! so no negative supply is allowed.
    !! @warning Mass is not conserved here!!
    !! @note TODO: do some accounting of this
    !! @note TODO: think about where this should be called
    class (modelstate_geom_t), intent(inout)    :: geometry

    geometry%bmelt = max (geometry%bmelt, 0.0_dp)
  end subroutine ensure_positive_supply

  subroutine balance_flux (flux, potential, geometry, method, mask)
    class (modelstate_flux_t), intent(inout) :: flux
    class (modelstate_pot_t), intent(in)     :: potential
    class (modelstate_geom_t), intent(inout)    :: geometry
    character(len=*), intent(in)             :: method
    integer                                  :: nneighbors
    integer, dimension(:,:), intent(in)      :: mask

    integer                                    :: dim(2)
    integer                                    :: nx, ny
    real(dp), dimension(:,:), allocatable    :: pot_x     !! d pot dx
    real(dp), dimension(:,:), allocatable    :: pot_y     !! d pot dy
    integer                                  :: i,j

    dim = shape(potential%pot)
    nx = dim(1)
    ny = dim(2)

    allocate (pot_x(nx,ny))
    allocate (pot_y(nx,ny))

    ! mark all cells as unvisited
    flux%sflux = -1.0_dp

    call ensure_positive_supply (geometry)

    select case (method)
    case ( 'warner','quinn' )
       if (method == 'warner' ) then
          nneighbors = 4
       elseif (method == 'quinn') then
          nneighbors = 8
       end if

       call flow_dir (nx, ny, geometry%x, geometry%y, mask, potential%pot_filled, pot_x, pot_y, flux%dir)
       ! loop over al the cells
       do j = 2, ny-1
          do i = 2, nx-1
             if (mask(i,j) > 0) then
                call dparea_war (i, j, nneighbors, nx, ny, &
                     geometry%x, geometry%y, mask, potential%pot_filled, &
                     geometry%bmelt, flux%sflux)
             end if
          end do
        end do

    case ( 'tarboton' )
        nneighbors = 8
       call tar_flow_dir (nx, ny, geometry%x, geometry%y, mask, potential%pot_filled, pot_x, pot_y, flux%dir)
       ! loop over al the cells
       do j = 2, ny-1
          do i = 2, nx-1
             if (mask(i,j) > 0) then
                call dparea_tar (i, j, nneighbors, nx, ny, &
                     geometry%x, geometry%y, mask, potential%pot_filled, &
                     geometry%bmelt, flux%dir, flux%sflux)
             end if
          end do
        end do
    case default
       print *, "#ERROR: unknown balance flux method: ", method
       stop 1
    end select
    if(allocated(pot_x)) deallocate (pot_x)
    if(allocated(pot_y)) deallocate (pot_y)

  end subroutine balance_flux

  subroutine flux_correction (flux, potential, geometry, method, mask)
    !! the computed flux field is not a direct measure of the flux
    !! magnitude and need to be scaled due to orientation of the flow
    !! (\see Fricker_et_al_2000, \see Budd_Warner_1996 )
    !!
    !! \see LeBrocq_et_al_2006: The flux routing algorithm calculate the
    !! scalar flux in m3/a.
    !! In order to convert the scalar flux to vector flux (m2/a) the
    !! lenght of the flow cross section needs to be considered.
    !!
    class (modelstate_flux_t), intent(inout) :: flux
    class (modelstate_pot_t), intent(in)     :: potential
    class (modelstate_geom_t), intent(in)    :: geometry
    character(len=*), intent(in)             :: method
    integer, dimension(:,:), intent(in)      :: mask

    integer                                    :: dim(2)
    integer                                    :: nx, ny
    real(dp), dimension(:,:), allocatable    :: pot_x     !! \partial pot / \partial x
    real(dp), dimension(:,:), allocatable    :: pot_y     !! \partial pot / \partial y
    integer                                  :: i,j

    real(dp)                                 :: denom
    real(dp)                                 :: corone
    real(dp)                                 :: corfac
    real(dp)                                 :: dx, dy

    dim = shape(potential%pot)
    nx = dim(1)
    ny = dim(2)

    allocate (pot_x(nx,ny))
    allocate (pot_y(nx,ny))
    select case (method)
    case ( 'warner','quinn' )
       call flow_dir (nx, ny, geometry%x, geometry%y, mask, potential%pot_filled, pot_x, pot_y, flux%dir)
    case ( 'tarboton' )
       call tar_flow_dir (nx, ny, geometry%x, geometry%y, mask, potential%pot_filled, pot_x, pot_y, flux%dir)
    case default
       print *, "#ERROR: unknown balance flux method: ", method
       stop
    end select

    do j = 2,ny-1
       do i = 2,nx-1
        denom = 1.0_dp
        corone = 1.0_dp
        corfac = 1.0_dp
        if (mask(i,j) > 0) then
          ! apply eq. (6) from Fricker_et_al_2000 also in BalanceV2.f
          ! (Budd_Warner_1996)
          corone = abs(pot_x(i,j)) + abs(pot_y(i,j))
          if (corone <= 0.0_dp ) then
            ! we take the flow direction es undifined
            corfac = 1.0_dp
          else
            dx = 0.5_dp*(geometry%x(i+1)-geometry%x(i-1))
            dy = 0.5_dp*(geometry%y(j+1)-geometry%y(j-1))
            !if (dx /= dy ) then
            if ( abs(dx - dy) > 1e-2_dp ) then
              write(*,'("error: dx != dy",f10.2,1x,f10.2)') dx,dy
              stop
            else
              denom = sqrt(pot_x(i,j)**2 + pot_y(i,j)**2)
              if (denom == 0.0_dp) then
                write (*,'(a,2(i4,1x),es13.6)') "denom = 0", i,j,denom
                ! TRY THIS
                ! denom = min(TINY,denom)
                stop
              end if
              ! TODO:
              ! basal water only works for dx=dy
              corfac = dx * corone / denom
            end if
          end if
          if (corfac == 0.0_dp) then
            write (*,'(a,2(i4,1x),es13.6)') "corfac = 0", i,j,corfac
            stop
          end if
          flux%vflux(i,j) = flux%sflux(i,j) / corfac
          flux%vfluxX(i,j) = flux%vflux(i,j) * cos(flux%dir(i,j))
          flux%vfluxY(i,j) = flux%vflux(i,j) * sin(flux%dir(i,j))
        end if
      end do
    end do

    if(allocated(pot_x)) deallocate (pot_x)
    if(allocated(pot_y)) deallocate (pot_y)

  end subroutine flux_correction

  subroutine flux2bwat(flux, potential, geometry, mask)
    !! solve for water layer thickness
    !! Weertman film or Manning-flow eqs. 2, 9 and 10 in
    !! \see LeBrocq_et_al_2009
    !! the viscosity of water can be found in Nye
    implicit none
    class (modelstate_flux_t), intent(inout) :: flux
    class (modelstate_pot_t), intent(in)     :: potential
    class (modelstate_geom_t), intent(in)    :: geometry
    integer, dimension(:,:), intent(in)      :: mask

    integer :: i,j
    real(dp) :: pot_x,pot_y,grad
    ! "Water at the glacier bed", J.F. Nye
    ! Weertman_Birchfield_1983
    real(dp), parameter :: visc_water = 1.8e-3_dp ! Pa*s = 1.8x10^-2 Poise
    integer :: jstart,jend
    integer                                    :: dim(2)
    integer                                    :: nx, ny
    dim = shape(potential%pot)
    nx = dim(1)
    ny = dim(2)
    jstart = 2
    jend = ny-1

    ! init to zero
    flux%bwat = 0.0_dp

    do j=jstart,jend
      do i = 2,nx-1
        grad = 0.0_dp
        pot_x = 0.0_dp
        pot_y = 0.0_dp
        if (mask(i,j) > 0) then
          pot_x = ( potential%pot(i+1,j) - potential%pot(i-1,j) ) / (geometry%x(i+1)-geometry%x(i-1))
          pot_y = ( potential%pot(i,j+1) - potential%pot(i,j-1) ) / (geometry%y(j+1)-geometry%y(j-1))
          grad = sqrt (pot_x**2 + pot_y**2)
          ! todo: find a solution for this
          ! if(cflux(i,j)< 0.0_dp) print *,i,j,cflux(i,j)," error"
          if(grad > 0.0_dp .and. flux%vflux(i,j) > 0.0_dp .and. mask(i,j) > 0) then
            flux%bwat(i,j) = ( 12.0_dp * visc_water * flux%vflux(i,j) / grad )**(1.0_dp/3.0_dp)
          else
            ! ! should be checked
            ! ! flux routing schemes are not possible for grad = 0
            ! ! this is maybe not the correct potential surface
            ! print *, i,j,cflux(i,j), grad
            ! ERRORMSG
            ! stop
            flux%bwat(i,j) = 0.0_dp
          end if
        else
          flux%bwat(i,j) = 0.0_dp
        end if
      end do
    end do

  end subroutine flux2bwat

  subroutine rise_potential_where_freezing (potential, geometry)
    !! The idea is to increase the hydraulic potential to a very high
    !! level (10 times the highest present level) at these points, so
    !! the flux will go around this.  This has to be done after the
    !! hydraulic potential has been computed!

    class (modelstate_pot_t), intent(inout)   :: potential
    class (modelstate_geom_t), intent(in)    :: geometry
    real(dp)                                 :: phi_safe  !! maximum present potential times 10, to be safe

    phi_safe = maxval(potential%pot_filled) * 10

    ! -1 because with 0 there is not much flowing anymore TODO:think about this
    where (geometry%temppabase < -1.0_dp)
      potential%pot_filled = potential%pot_filled + phi_safe
    endwhere
    ! make sure, the potential is low, at places with no ice, so the whole area can drain
    where (geometry%icemask == 0)
      potential%pot_filled = potential%pot_filled - phi_safe
    endwhere

  end subroutine

  pure function method2nneighbors (method) result (nneighbors)
    !! Return the number of considered neighbors for each method
    character(len=*), intent(in)  :: method
    integer                       :: nneighbors

    select case (method)
    case ('warner')
       nneighbors = 4
    case ('quinn')
       nneighbors = 8
    case ('tarboton')
       nneighbors = 8
    case default
       nneighbors = 0
    end select
  end function method2nneighbors

  pure function libhydro_version () result (version_str)
    !! return the version number of libhydro
    integer             :: version_major
    integer             :: version_minor
    integer             :: version_patch
    character(len=10)   :: version_major_str
    character(len=10)   :: version_minor_str
    character(len=10)   :: version_patch_str
    character(len=80)   :: version_str

    version_major = libhydro_VERSION_MAJOR
    version_minor = libhydro_VERSION_MINOR
    version_patch = libhydro_VERSION_PATCH

    write ( version_major_str, '(i0)') version_major
    write ( version_minor_str, '(i0)') version_minor
    write ( version_patch_str, '(i0)') version_patch

    version_str = trim(version_major_str) //'.'// &
         & trim(version_minor_str) //'.'// &
         & trim(version_patch_str)

  end function libhydro_version

  ! -----------
  ! PUBLIC API
  ! -----------

  subroutine hydro_gen_conf (hydro, method, flatresolution, &
       avoid_frz, filter_len, rho_seawater, rho_freshwater, rho_ice, min_thickness)
    type (hydro_t), intent(inout)          :: hydro
    character(len=*), intent(in),optional  :: method
    character(len=*), intent(in),optional  :: flatresolution
    logical, intent(in),optional           :: avoid_frz
    real(dp), intent(in),optional          :: filter_len
    real(dp), intent(in),optional          :: rho_seawater
    real(dp), intent(in),optional          :: rho_freshwater
    real(dp), intent(in),optional          :: rho_ice
    real(dp), intent(in),optional          :: min_thickness

    if (present(method)) hydro%conf%method = method
    if (present(flatresolution)) hydro%conf%flatresolution = flatresolution
    if (present(avoid_frz)) hydro%conf%avoid_frz = avoid_frz
    if (present(filter_len)) hydro%conf%filter_len = filter_len
    if (present(rho_seawater)) hydro%conf%rho_seawater = rho_seawater
    if (present(rho_freshwater)) hydro%conf%rho_freshwater = rho_freshwater
    if (present(rho_ice)) hydro%conf%rho_ice = rho_ice
    if (present(min_thickness)) hydro%conf%min_thickness = min_thickness
  end subroutine hydro_gen_conf

  subroutine hydro_init (hydro, x, y)
    !! initialize hydrology model
    type (hydro_t), intent(inout)    :: hydro
    real(dp), intent(in), dimension(:):: x !! x coordinates
    real(dp), intent(in), dimension(:):: y !! y coordinates

    integer                          :: nx, ny

    nx = size(x)
    ny = size(y)

    call hydro%state%alloc (nx, ny)
    call hydro%state%init ()
    call hydro%conf%init ()
    hydro%state%geometry%x = x
    hydro%state%geometry%y = y

    hydro%state%geometry%nx = nx
    hydro%state%geometry%ny = ny
    hydro%state%geometry%dx = x(2) - x(1)
  end subroutine hydro_init

  subroutine hydro_update (hydro)
    !! update hydrology model: this is the main routine
    !! calculate new potential, route the water and
    !! generate flux and basal water layer thickness
    type (hydro_t), intent(inout)   :: hydro

    if ( hydro%conf%comp_mask ) then
        call generateIceMask  (hydro%state%geometry,  &
            & hydro%conf%rho_ice, &
            & hydro%conf%rho_seawater,  &
            & min_thickness=hydro%conf%min_thickness)
    end if
    call hydro_potential  (hydro%state%potential, hydro%state%geometry, &
        & hydro%conf%rho_freshwater, &
        & hydro%conf%rho_ice, &
        & hydro%conf%gravity_const)
    call smooth_potential (hydro%state%potential, hydro%conf%filter_len, &
         & hydro%state%geometry%dx)
    if (hydro%conf%flatresolution == 'fast') then
       ! if fast is chosen, fill_sinks takes care of resolving flats
       call fill_sinks (hydro%state%potential, hydro%conf%method, &
            & tiltflats=.true.)
    else
       call fill_sinks (hydro%state%potential, hydro%conf%method, &
            & tiltflats=.false.)
       call resolve_flats (hydro%state%potential, hydro%conf%method, &
            & hydro%conf%flatresolution)
    end if
    if (hydro%conf%avoid_frz) then
       call rise_potential_where_freezing (hydro%state%potential, &
            & hydro%state%geometry)
    end if
    call balance_flux (hydro%state%flux, hydro%state%potential, &
         & hydro%state%geometry, hydro%conf%method, hydro%state%geometry%icemask)
    call flux_correction (hydro%state%flux, hydro%state%potential, &
         & hydro%state%geometry, hydro%conf%method, hydro%state%geometry%icemask)
    ! calculate water layer thickness
    call flux2bwat (hydro%state%flux, hydro%state%potential, &
         & hydro%state%geometry, hydro%state%geometry%icemask)

    ! set water layer thickness to certain amount at lake locations
    ! but avoid doing this where the base is frozen
    where ( hydro%state%potential%pot_diff > 1.0_dp .and. &
            hydro%state%geometry%temppabase > -0.001_dp) hydro%state%flux%bwat = 10.0_dp

  end subroutine hydro_update

  ! -----------------------------
  ! setters and getters
  ! -----------------------------
  subroutine hydro_get_bwat (hydro, bwat)
    !! return basal water layer thickness
    type(hydro_t),intent(in)              :: hydro
    real(dp), intent(out), dimension(:,:) :: bwat !! basal water layer thickness (meter)
    bwat = hydro%state%flux%bwat
  end subroutine hydro_get_bwat

  subroutine hydro_get_sflux (hydro, sflux)
    !! return scalar (volume) flux
    type(hydro_t),intent(in)              :: hydro
    real(dp), intent(out), dimension(:,:) :: sflux !! scalar (volume) flux (meter^3/sec)
    sflux = hydro%state%flux%sflux
  end subroutine hydro_get_sflux

  subroutine hydro_get_vflux (hydro, vflux, vfluxX, vfluxY)
    !! return vector (volume) flux
    type(hydro_t),intent(in)              :: hydro
    real(dp), intent(out), dimension(:,:) :: vflux, vfluxX, vfluxY
                                             !! vector (volume) flux (meter^2/sec)
    vflux  = hydro%state%flux%vflux
    vfluxX = hydro%state%flux%vfluxX
    vfluxY = hydro%state%flux%vfluxY
  end subroutine hydro_get_vflux

  subroutine hydro_set_mask (hydro, mask)
    !! set mask in hydro model
    !! This is a integer mask with 1 means routing should
    !! happen and 0 means no routing.
    type(hydro_t),intent(inout)          :: hydro
    integer, intent(in), dimension(:,:) :: mask !! hydro mask, 1 means routing, 0 no routing
    hydro%state%geometry%icemask = mask
    hydro%conf%comp_mask = .false.
    end subroutine hydro_set_mask

  subroutine hydro_set_supply (hydro, supply)
    !! set supply in hydro model
    type(hydro_t),intent(inout)          :: hydro
    real(dp), intent(in), dimension(:,:) :: supply !! water supply (meter per second)
    hydro%state%geometry%bmelt = supply
  end subroutine hydro_set_supply

  subroutine hydro_set_temppabase (hydro, temppabase)
    !! set temppabase in hydro model
    type(hydro_t),intent(inout)          :: hydro
    real(dp), intent(in), dimension(:,:) :: temppabase !! basal temperature rel. to pmp
    hydro%state%geometry%temppabase = temppabase
  end subroutine hydro_set_temppabase

  subroutine hydro_set_thk (hydro, thk)
    !! set ice thickness in hydro model
    type(hydro_t),intent(inout)          :: hydro
    real(dp), intent(in), dimension(:,:) :: thk !! ice thickness (meter)
    hydro%state%geometry%thk = thk
  end subroutine hydro_set_thk

  subroutine hydro_set_topg (hydro, topg)
    !! set bedrock topography in hydro model
    type(hydro_t),intent(inout)          :: hydro
    real(dp), intent(in), dimension(:,:) :: topg !! basal topography (meter)
    hydro%state%geometry%topg = topg
  end subroutine hydro_set_topg

!-------------------------------------------------------------------------------

end module hydro_m
!
