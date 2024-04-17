!-------------------------------------------------------------------------------
!> Check if point is inside a general polygon.
!-------------------------------------------------------------------------------
module point_in_polygon_m

#if (defined(MODEL_SICOPOLIS))
  use sico_types_m, only: dp
#endif

  implicit none

  public

#if (!defined(MODEL_SICOPOLIS))
  integer, parameter :: dp = kind(1.0d0)   ! double precision
#endif

contains

!-------------------------------------------------------------------------------
!> Check if point is inside a general polygon
!!
!! Input parameters:
!!   'n' is the number of sides/vertices defining the polygon
!!   'smalld' is a small double precision number
!!   'larged' is a large double precision number
!!   'xc' is a vector of nodal x-coords (anticlockwise order)
!!   'yc' is a vector of nodal y-coords (anticlockwise order)
!!        Both of these vectors must be of length n
!!   'xq' is the x-coord of the point to be tested (query point)
!!   'yq' is the y-coord of tbe point to be tested (query point)
!!
!! Output parameters:
!!   'mindst' the distance from the point to the nearest point
!!            on the polygon
!!     If 'mindst' is lt zero then point is outside the polygon
!!     If 'mindst' is eq zero then point is on a side of the polygon
!!     If 'mindst' is gt zero then point is inside the polygon
!!
!! Notes:
!!   An improved version of the algorithm of Nordbeck and Rydstedt
!!   by S. W. Sloan (1985, Adv. Eng. Software 7(1), 45-47),
!!   modified and ported to Fortran 90 by R. Greve (December 2018)
!-------------------------------------------------------------------------------
  subroutine point_in_polygon(n, smalld, larged, xc, yc, xq, yq, mindst)

  implicit none

  integer , intent(in)  :: n
  real(dp), intent(in)  :: smalld, larged
  real(dp), intent(in), dimension(:) :: xc, yc
  real(dp), intent(in)  :: xq, yq
  real(dp), intent(out) :: mindst

  integer  :: i, j
  real(dp), dimension(n+2) :: x, y
  real(dp) :: d, area
  real(dp) :: x1, y1, x1p, y1p, x21, y21, t, dx, dy
  character(len=256) :: errmsg
  logical  :: snear

  real(dp), parameter :: c00000 = 0.0_dp, c00001 = 1.0_dp

  ! 'snear' is .true. if distance to nearest side is less than
  !         distance to nearest vertex
  ! 'snear' is .false. if distance to nearest vertex is less than
  !         distance to nearest side
  ! 'mindst' is square of distance to closest point on the polygon

  if ((size(xc) /= n).or.(size(yc) /= n)) then
     errmsg = ' >>> point_in_polygon: length of xc or yc not equal to n!'
     write(6, fmt='(/,a,/)') trim(errmsg)
     stop
  end if

  mindst = larged

  ! Construct vectors x and y of nodal x-/y-coords with length n+2

  do i=1, n
     x(i) = xc(i); y(i) = yc(i)
  end do
  x(n+1) = xc(1); y(n+1) = yc(1)
  x(n+2) = xc(2); y(n+2) = yc(2)

  ! Loop over each side defining polygon

  do i=1, n

     ! Start of side has coords (x1,y1)
     ! End of side has coords (x2,y2)
     ! Point has coords (xq,yq)

     x1  = x(i)
     y1  = y(i)
     x21 = x(i+1)-x1
     y21 = y(i+1)-y1
     x1p = x1-xq
     y1p = y1-yq

     ! Points on infinite line defined by
     !   x=x1+t*(x1-x2)
     !   y=y1+t*(y1-y2)
     ! where
     !   t=0 at (x1,y1)
     !   t=1 at (x2,y2)
     ! Find where normal passing through (xq,yq)
     ! intersects infinite line

     t = -(x1p*x21+y1p*y21)/(x21*x21+y21*y21)

     if (t < c00000) then

        ! Normal does not intersect side
        ! Point is closest to vertex (x1,y1)
        ! Compute square of distance to this vertex

        d = x1p*x1p+y1p*y1p

        if (d < mindst) then
           ! Point is closer to (x1,y1) than any other vertex or side
           snear  = .false.
           mindst = d
           j      = i 
        end if

     else if (t <= c00001) then

        ! Normal intersects side

        dx = x1p+t*x21
        dy = y1p+t*y21
        d  = dx*dx+dy*dy

        if (d < mindst) then
           ! Point is closer to this side than to any other side or vertex
           snear  = .true.
           mindst = d
           j      = i
        end if

     end if

  end do

  mindst = sqrt(mindst)

  if (mindst < smalld) then
     ! Point is on side of polygon
     mindst = c00000 
  else
     if (snear) then
        ! Point is closer to its nearest side than to its nearest
        ! vertex, check if point is to left or right of this side
        ! If point is to left of side it is inside polygon,
        ! else point is outside polygon
        area   = det(x(j), x(j+1), xq, y(j), y(j+1), yq)
        mindst = sign(mindst, area)
     else
        ! Point is closer to its nearest vertex than its nearest side,
        ! check if nearest vertex is concave
        ! If the nearest vertex is concave then point is inside the polygon,
        ! else the point is outside the polygon
        if (j == 1) then
           j=n+1
        end if
        area   = det(x(j+1), x(j), x(j-1), y(j+1), y(j), y(j-1))
        mindst = sign(mindst, area)
     end if
  end if 

  end subroutine point_in_polygon

!-------------------------------------------------------------------------------
!> Compute twice the area of the triangle defined by three points
!! with coords (x1,y1), (x2,y2) and (x3,y3) using determinant formula
!!
!! Input parameters:
!!   'x1,y1' coords of point 1
!!   'x2,y2' coords of point 2
!!   'x3,y3' coords of point 3
!!
!! Output parameters:
!!   'det' twice the area of the triangle defined by the three points
!!
!! Notes:
!!   det is positive if points 1, 2 and 3 define triangle in
!!   anticlockwise order
!!   det is negative if points 1, 2 and 3 define triangle in
!!   clockwise order
!!   det is zero if at least two of the points are coincident
!!   or if all three points are collinear
!-------------------------------------------------------------------------------
  function det(x1, x2, x3, y1, y2, y3)

  implicit none

  real(dp), intent(in) :: x1, x2, x3, y1, y2, y3

  real(dp) :: det

  det = (x1-x3)*(y2-y3)-(x2-x3)*(y1-y3)

  end function det

!-------------------------------------------------------------------------------

end module point_in_polygon_m
!
