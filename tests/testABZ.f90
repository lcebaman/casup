!*robodoc*u* tests/testABZ
!  NAME
!    testABZ
!  SYNOPSIS

!$Id: testABZ.f90 380 2017-03-22 11:03:09Z mexas $

program testABZ

!  PURPOSE
!    Testing cgca_boxsplit from module cgca_m2geom.
!  NOTE
!    Serial routine, a single image will do.
!  DESCRIPTION
!    cgca_boxsplit splits a box, given by 2 corner coordinates,
!    all integers, into two smaller boxes, along the biggest
!    dimension of the old box.
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  USES
!    cgca testaux
!  USED BY
!    Part of CGPACK test suite
!  SOURCE

use testaux
implicit none

integer( kind=idef) :: lwr(3), upr(3), lwr1(3), upr1(3), lwr2(3),      &
   upr2(3), img, i

! Only image 1 works
img = this_image()
main: if ( img .eq. 1 ) then

! print a banner
call banner("ABZ")

! set some values
lwr = (/ 1, 1, 1 /)
upr = (/ 123, 456, 789 /)

call cgca_boxsplit( lwr, upr, lwr1, upr1, lwr2, upr2 )
write (*,101) "lwr: ", lwr, "upr:", upr
write (*,*) "split into:"
write (*,101) "lwr1: ", lwr1, "upr1:", upr1
write (*,101) "lwr2: ", lwr2, "upr2:", upr2

do i=1,100

  if ( mod(i,2) .eq. 0 ) then
    ! even i
    write (*,*) "choose box 1"
    lwr = lwr1
    upr = upr1
  else
    ! odd i
    write (*,*) "choose box 2"
    lwr = lwr2
    upr = upr2
  end if

  call cgca_boxsplit( lwr, upr, lwr1, upr1, lwr2, upr2 )
  write (*,101) "lwr: ", lwr, "upr:", upr
  write (*,*) "split into:"
  write (*,101) "lwr1: ", lwr1, "upr1:", upr1
  write (*,101) "lwr2: ", lwr2, "upr2:", upr2

  if ( all( lwr1 .eq. upr1 ) .and. all( lwr2 .eq. upr2 ) ) exit

end do

end if main

101 format (2(a5,3(i5),tr1))

end program testABZ

!*roboend*
