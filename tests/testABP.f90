!*robodoc*u* tests/testABP
!  NAME
!    testABP
!  SYNOPSIS

!$Id: testABP.f90 380 2017-03-22 11:03:09Z mexas $

program testABP
      
!  PURPOSE
!    Checking: cgca_tchk
!  DESCRIPTION
!    Checking the MAXMIN value of the dot product between
!    an arbitrary cleavage plane normal and all 26 unit vectors
!    connecting the central cell with its neighbours.
!    cgca_tchk is a serial routine, so can make only image 1 call it,
!    or even better, make all images execute it, to increase the
!    search space. No sync required.
!  NOTES
!    The program must be called with 2 command line arguments,
!    both positive integers. These are codimensions along 1 and 2.
!    The number of images must be such that
!    codimension3 = num_images()/( codimension1 * codimension3 )
!    is a positive integer. Example:
!      ./testABP.x 2 2
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

integer( kind=idef ) :: nimages, codim(3)[*], image
logical( kind=ldef ) :: image1
real(    kind=rlrg ) :: maxmin, minmax

!*********************************************************************72
! first executable statement

nimages = num_images()
  image = this_image()
 image1 = .false.
if (this_image().eq.1) image1 = .true.

! do a check on image 1
if (image1) then
 call getcodim( nimages, codim )
 ! print a banner
 call banner("ABP")
 ! print the parameter values
 call cgca_pdmp
 write (*,'(a,i0,a)') "running on ", nimages, " images in a 3D grid"
 write (*,*) "codim:", codim
end if

sync all

codim(:) = codim(:)[1]

! initialise random seed
call cgca_irs( debug = .false. )

! check threshold t
call cgca_tchk( 2_ilrg**32, maxmin, minmax ) ! 4,294,967,296
write( *, "(a, i0, 2(a,es20.10))" ) "image: ", image, " maxmin: ",      &
                         maxmin, " minmax: ", minmax

end program testABP

!*roboend*
