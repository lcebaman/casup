!*robodoc*m* CGPACK/cgca_m4fr
!  NAME
!    cgca_m4fr
!  SYNOPSIS

!$Id: cgca_m4fr.f90 400 2017-05-04 17:47:56Z mexas $

module cgca_m4fr

!  DESCRIPTION
!    Module dealing with fracture
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!  USES
!    cgca_m1co, cgca_m2gb, cgca_
!  USED BY
!  SOURCE

use cgca_m1co
use cgca_m2gb
use cgca_m3clvg
use cgca_m3gbf
implicit none
private

contains

!*roboend*


!*robodoc*s* cgca_m4fr/cgca_fr
!  NAME
!    cgca_fr
!  SYNOPSIS

subroutine cgca_fr( coarray, rt, s1, scrit, periodicbc, iter, &
 heartbeat, debug )

!  INPUTS

integer( kind=iarr ), allocatable, intent( inout ) :: &
  coarray(:,:,:,:)[:,:,:]
real( kind=rdef ), allocatable, intent(inout) :: rt(:,:,:)[:,:,:]
real( kind=rdef ), intent(in) :: s1(3), scrit(3)
logical(kind=ldef),intent(in) :: periodicbc
integer(kind=idef),intent(in) :: iter,heartbeat
logical(kind=ldef),intent(in) :: debug

!  SIDE EFFECTS
!    state of coarray changes
!  DESCRIPTION
!    This routine does one iteration of cleavage propagation, followed
!    by an iteration of grain boundary fracture. It does the halo
!    exchange when required. Then it repeats this cycle for the given
!    number of iterations, "iter".

!integer(kind=iarr),allocatable,save :: array(:,:,:)

integer( kind=idef ) :: &
  lbv(4) ,& ! lower bounds of the complete (plus virtual) coarray
  ubv(4) ,& ! upper bounds of the complete (plus virtual) coarray
  lbr(4) ,& ! lower bounds of the "real" coarray, lower virtual+1
  ubr(4) ,& ! upper bounds of the "real" coarray, upper virtual-1
  x1     ,& ! local coordinates in an array, which are also
  x2     ,& ! do loop counters
  x3,     &
  iteration

integer :: img, nimgs

  img = this_image()
nimgs = num_images()

! Do not check coarray for allocated, as this wastes time.
! Instead let the code fail if coarray is not allocated.

! determine the extents
lbv=lbound(coarray)
ubv=ubound(coarray)
lbr=lbv+1
ubr=ubv-1

! dummy code to suppress the warnings
if ( periodicbc ) then
  write (*,*) s1 * scrit * sum( rt )
end if

itera: do iteration = 1, iter

do x1 = lbr(1), ubr(1)
do x2 = lbr(2), ubr(2)
do x3 = lbr(3), ubr(3)
end do
end do
end do

if ( mod( iteration, heartbeat ) .eq. 0 ) then
  write (*,'(a)') "INFO:"
end if

end do itera

! debug
if ( debug ) then
end if

end subroutine cgca_fr

!*roboend*

end module cgca_m4fr
