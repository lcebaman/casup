!*robodoc*f* cgca_m3clvg/m3clvg_sm2
!  NAME
!    m3clvg_sm2
!  SYNOPSIS

!$Id: m3clvg_sm2.f90 380 2017-03-22 11:03:09Z mexas $

submodule ( cgca_m3clvg ) m3clvg_sm2

!  DESCRIPTION
!    Submodule with aux routine for checking the misorientation
!    thereshold.
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    cgca_tchk
!  USES
!    Variables and routines of module cgca_m3clvg by host
!    association.
!  USED BY
!    Module cgca_m3clvg
!  SOURCE

implicit none

contains

!*roboend*


!*robodoc*s* m3clvg_sm2/cgca_tchk
!  NAME
!    cgca_tchk
!  SYNOPSIS  

module procedure cgca_tchk

!  DESCRIPTION
!    Generate num random unit normal vectors.
!    Calculate the MAX of the MIN for all 26 cell
!    neighbourhood unit vectors.
!    Find the maximum value of all num normal vectors
!    and return it as maxmin.
!    And the opposite - calculate the MIN of the MAX for
!    all 26 cell neighbourhood unit vectors.
!    Find the min value of all num normal vectors
!    and return as minmax.
!  NOTE
!    Any image can call this routine
!  SOURCE

real(    kind=rlrg ) :: n(3), mag, prod, prodmax, prodmin
integer( kind=ilrg ) :: i, x1, x2, x3

maxmin = 0.0_rlrg
minmax = 1.0_rlrg

do i=1,num
  prodmin = 1.0_rlrg
  prodmax = 0.0_rlrg
  call random_number(n)
  mag = sqrt( sum( n**2 ) )
  n = n / mag
  do x3 = -1, 1
  do x2 = -1, 1
  do x1 = -1, 1
    if ( x1 .eq. 0 .and. x2 .eq. 0 .and. x3 .eq. 0) cycle
    prod = abs( dot_product( e(:,x1,x2,x3), n ) )
    if ( prod .lt. prodmin ) prodmin = prod
    if ( prod .gt. prodmax ) prodmax = prod
  end do
  end do
  end do

  if ( prodmin .gt. maxmin) maxmin = prodmin
  if ( prodmax .lt. minmax) minmax = prodmax

end do

end procedure cgca_tchk

!*roboend*

end submodule m3clvg_sm2
