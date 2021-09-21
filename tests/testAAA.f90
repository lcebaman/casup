!*robodoc*u* tests/testAAA
!  NAME
!    testAAA
!  SYNOPSIS

!$Id: testAAA.f90 380 2017-03-22 11:03:09Z mexas $

program testAAA

!  PURPOSE
!    Checking: getcodim, cgca_as, cgca_ds 
!  DESCRIPTION
!    Testing allocating and deallocating of a coarray.
!  NOTES
!    The program must be called with 2 command line arguments,
!    both positive integers. These are codimensions along 1 and 2.
!    The number of images must be such that
!    codimension3 = num_images()/( codimension1 * codimension3 )
!    is a positive integer. Example:
!      cafrun -np 16 ./testAAA.x 2 2      ! OpenCoarrays
!    or
!      ./testAAA.x 2 2                    ! Intel, Cray
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

integer( kind=idef ), parameter :: size1=10, size2=10, size3=10
integer( kind=idef ) :: nimgs
integer( kind=idef ) :: codim(3)[*]
integer( kind=iarr ), allocatable :: space1( : , : , : , : ) [:,:,:]

!*********************************************************************72
! first executable statement

nimgs = num_images()

! do a check on image 1
if (this_image() .eq. 1) then
 call getcodim( nimgs, codim )
 ! print a banner
 call banner("AAA")
 ! print the parameter values
 call cgca_pdmp
 write (*, '(a,i0,a)') "running on ", nimgs, " images in a 3D grid"
 write (*, *) "codim:", codim
end if

sync all ! wait for image 1 to set codim

codim(:) = codim(:)[1]

sync all ! wait for each image to read codim from img 1.

! implicit sync all inside
call cgca_as( 1, size1, 1, size2, 1, size3, 1, codim(1), 1, codim(2),  &
              1, 2, space1 )

if ( allocated( space1 ) )                                             &
 write (*, "(2(a,i0), 3(a,4(i0,tr1)), 3(a,3(i0,tr1)) )")               &
  "img: ", this_image(), ". my array, size: ", size( space1 ),         &
  ". shape: " , shape(space1), ". lbound: ", lbound(space1),           &
  ". ubound:", ubound(space1), ". coar index: ", this_image( space1 ), &
 ". lcobound:", lcobound(space1), ". ucobound:", ucobound(space1)

call cgca_ds(space1)
 
if ( .not. allocated(space1) )                                         &
  write (*,'(a,i0,a)')"Image:",this_image(), " space1 not allocated"

end program testAAA 

!*roboend*
