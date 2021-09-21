!*robodoc*u* tests/testABH
!  NAME
!    testABH
!  SYNOPSIS

!$Id: testABH.f90 389 2017-03-22 16:31:21Z mexas $

program testABH

!  PURPOSE
!    Checking: cgca_redand, part of cgca_m2red
!  DESCRIPTION
!    Checking collective AND reduction over a logical coarray.
!    Works only when the number of images is 2**p,
!    where p is an integer, so use 2, 4, 8, 16, 32, etc. images. 
!  NOTES
!    The program must be called with 2 command line arguments,
!    both positive integers. These are codimensions along 1 and 2.
!    The number of images must be such that
!    codimension3 = num_images()/( codimension1 * codimension3 )
!    is a positive integer. Example:
!      cafrun -np 16 ./testABH.x 2 2      ! OpenCoarrays
!    or
!      ./testABH.x 2 2                    ! Intel, Cray
!    which will make the third codimension equal to 16/(2*2)=4.
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

real, parameter :: l2 = log(real(2))
logical, parameter :: nodebug = .false.
real :: num
integer(kind=idef) :: p, nimages, img, codim(3)[*]
logical(kind=ldef) :: z[*]

!*********************************************************************72
! first executable statement

nimages=num_images()
img = this_image()

! check than n is a power of 2
p = nint(log(real(nimages))/l2)
if ( 2**p .ne. nimages)                                                &
      error stop "number of images is not a power of 2"

! do a check on image 1
if (img .eq. 1) then
 call getcodim(nimages,codim)
 ! print a banner
 call banner("ABH")
 write (*,'(a,i0,a)') "running on ", nimages, " images in a 3D grid"
end if

! Trying to separate the output
sync all

! initialise random number seed
call cgca_irs(nodebug)

! assign z
call random_number(num)

if (num .gt. 0.5) then
 z = .true.
else
 z = .false.
end if

z = .true.
if (img .eq. nimages) z = .false.

write (*,*) "image", img, "z", z

! Trying to separate the output
sync all

! call collective AND
call cgca_redand(z,p)

write (*,*) "image", img, "answer", z

end program testABH

!*roboend*
