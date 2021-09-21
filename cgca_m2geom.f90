!*robodoc*m* CGPACK/cgca_m2geom
!  NAME
!    cgca_m2geom
!  SYNOPSIS

!$Id: cgca_m2geom.f90 379 2017-03-22 09:57:10Z mexas $

module cgca_m2geom

!  DESCRIPTION
!    Module dealing with various 3D geometrical problems
!  AUTHOR 
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    cgca_boxsplit
!  USES
!    cgca_m1co
!  USED BY
!    cgca_m3pfem
!  SOURCE

use cgca_m1co, only : idef
implicit none

private
public :: cgca_boxsplit

contains

!*roboend*


!*robodoc*s* cgca_m2geom/cgca_boxsplit
!  NAME
!    cgca_boxsplit
!  SYNOPSIS

subroutine cgca_boxsplit( lwr, upr, lwr1, upr1, lwr2, upr2 )

!  INPUTS
!    lwr(3) - integer, lower corner of the box
!    upr(3) - integer, upper corner of the box

integer( kind=idef ), intent( in ) :: lwr(3), upr(3)

!  OUTPUTS
!    lwr1(3) - integer, lower corner of new box 1
!    upr1(3) - integer, upper corner of new box 1
!    lwr2(3) - integer, lower corner of new box 2
!    upr2(3) - integer, upper corner of new box 2

integer( kind=idef ), intent( out ) :: lwr1(3), upr1(3), lwr2(3),      &
  upr2(3)

!  DESCRIPTION
!    This routine splits the box, specified by two corner
!    coordinates into two smaller boxes, along the biggest dimension
!    of the original box.
!  SOURCE

integer( kind=idef ) :: boxsize(3), splitdim

! If the box is only a single cell, return immediately
if ( all( lwr .eq. upr ) ) then
  lwr1 = lwr
  lwr2 = lwr
  upr1 = upr
  upr2 = upr
  return
end if

! Find the biggest dimension of the box.
boxsize = upr - lwr + 1
splitdim = maxloc( boxsize, dim=1 ) ! 1, 2 or 3 only

! Set the dimensions of each new box initially equal to
! the old box
  lwr1 = lwr
  upr1 = upr
  lwr2 = lwr
  upr2 = upr

! Change only relevant dimensions
if ( splitdim .eq. 1 ) then
  upr1(1) = ( lwr(1) + upr(1) ) / 2 ! new box 1
  lwr2(1) = upr1(1) + 1             ! new box 2
else if ( splitdim .eq. 2 ) then
  upr1(2) = ( lwr(2) + upr(2) ) / 2 ! new box 1
  lwr2(2) = upr1(2) + 1             ! new box 2
else if ( splitdim .eq. 3 ) then
  upr1(3) = ( lwr(3) + upr(3) ) / 2 ! new box 1
  lwr2(3) = upr1(3) + 1             ! new box 2
end if

end subroutine cgca_boxsplit

!*roboend*

end module cgca_m2geom
