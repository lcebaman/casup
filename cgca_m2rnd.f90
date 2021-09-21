!*robodoc*m* CGPACK/cgca_m2rnd
!  NAME
!    cgca_m2rnd
!  SYNOPSIS

!$Id: cgca_m2rnd.f90 380 2017-03-22 11:03:09Z mexas $

module cgca_m2rnd

!  DESCRIPTION
!    Module dealing with random number generation
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    cgca_irs, cgca_ins
!  USES
!    cgca_m1co
!  USED BY
!    cgca
!  SOURCE

use cgca_m1co, only: ldef
implicit none

private
public :: cgca_irs, cgca_ins

contains

!*roboend*

!*robodoc*s* cgca_m2rnd/cgca_irs
!  NAME
!    cgca_irs
!  SYNOPSIS

subroutine cgca_irs(debug)

!  INPUT

logical(kind=ldef),intent(in) :: debug

!  SIDE EFFECTS
!    initialise random seed on all images   
!  DESCRIPTION
!    Initialise random seed based on system_clock intrinsic,
!    adapted from:
!    http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html.
!    Note that the seed is based on THIS_IMAGE intrinsic, hence
!    each image uses a different seed.
!  USES
!    none
!  USED BY
!    none, end user
!  SOURCE

integer :: i, n, clock, errstat=0
integer, allocatable :: seed(:)
         
call random_seed( size = n )

allocate( seed(n), stat=errstat )
if ( errstat .ne. 0 ) then
  write (*,*) "ERROR: cgca_irs/cgca_m2rnd: allocate( seed )," //       &
              " err. stat:", errstat
  error stop
end if

call system_clock(count=clock)
          
seed = int(real(clock)/real(this_image())) +                           &
  999999937*(/ (i - 1, i = 1, n) /)

call random_seed(put = seed)
          
! Debug output
if (debug) write (*,*) "DEBUG: cgca_irs/cgca_m2rnd: this_image():",    &
  this_image(), "; size:", n, "; seed:", seed

deallocate( seed, stat=errstat )
if ( errstat .ne. 0 ) then
  write (*,*) "ERROR: cgca_irs/cgca_m2rnd: deallocate( seed )," //     &
              " err. stat:", errstat
  error stop
end if

end subroutine cgca_irs

!*roboend*


!*robodoc*s* cgca_m2rnd/cgca_ins
!  NAME
!    cgca_irs
!  SYNOPSIS

subroutine cgca_ins(debug)

!  INPUT

logical( kind=ldef ), intent(in) :: debug

!  SIDE EFFECTS
!    Initialise the seed based only on the image number.
!  DESCRIPTION
!    This routine sets reproducible random seeds on each image.
!    If the number of images is kept constant, then the results
!    of the CGPACK simulation should be (need to prove this
!    rigorously!) reproducible. If some changes to the code
!    are not supposed to change the results, use this random
!    seed routine.
!  USES
!    none
!  USED BY
!    none, end user
!  SOURCE

integer :: i, n, errstat=0
integer, allocatable :: seed(:)
         
call random_seed( size = n )

allocate( seed(n), stat=errstat )
if ( errstat .ne. 0 ) then
  write (*,*) "ERROR: cgca_ins/cgca_m2rnd: allocate( seed )," //       &
              " err. stat:", errstat
  error stop
end if

! Set seed to this_image() related values.
seed = this_image() * (/ (i, i=1, n) /)
call random_seed(put = seed)

! Don't need to do this, but just to be double sure, read
! the seed back from the subroutine
call random_seed(get = seed)

! Debug output
if (debug) write (*, "(3(a,i0),9999(tr1,i0))" )                        &
  "DEBUG: cgca_ins/cgca_m2rnd: this_image(): ", this_image(),          &
  " size: ", n, " seed: ", seed

deallocate( seed, stat=errstat )
if ( errstat .ne. 0 ) then
  write (*,*) "ERROR: cgca_ins/cgca_m2rnd: deallocate( seed )," //     &
              " err. stat:", errstat
  error stop
end if

end subroutine cgca_ins

!*roboend*

end module cgca_m2rnd
