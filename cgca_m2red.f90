!*robodoc*m* CGPACK/cgca_m2red
!  NAME
!    cgca_m2red
!  SYNOPSIS

!$Id: cgca_m2red.f90 380 2017-03-22 11:03:09Z mexas $

module cgca_m2red

!  DESCRIPTION
! Module dealing with collective reduction operations, including
! all required image syncronisation.
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    cgca_redand
!  USES
!    cgca_m1co
!  USED BY
!    cgca
!  SOURCE

use cgca_m1co
implicit none

private
public :: cgca_redand

contains

!*roboend*


!*robodoc*s* cgca_m2red/cgca_redand
!  NAME
!    cgca_redand
!  SYNOPSIS

subroutine cgca_redand(coarray,p)

!  INPUTS

logical(kind=ldef),intent(inout) :: coarray[*]
integer(kind=idef),intent(in) :: p

!  SIDE EFFECTS
!    coarray values change
!  DESCRIPTION
! This routine does collective AND operation over coarray values across
! all images. The result is returned in coarray on every image.
! The result is TRUE
! if coarray values on all images are TRUE, and FALSE otherwise.
! The algorithm implements a divide and conquer scheme that works
! only when the number of images, n, is a power of 2 - n=2**p.
! p is the input to this routine.
!
! If the number of images is 2**p,
! then reduction takes p iterations.
! In this example I have 2**4=16, so it takes 4 iterations.
!
!     img1 img2 img3 img4 img5 img6 img7 img8 img9 img10 img11 img12 img13 img14 img15 img16
!  1. img1 _/   img3 _/   img5 _/   img7 _/   img9 _/    img11 _/    img13 _/    img15 _/
!  2. img1 ______/        img5 ______/        img9 _______/          img13 _______/
!  3. img1 ________________/                  img9 ___________________/
!  4. img1 ____________________________________/
!     img1
!
!  NOTE
!   For efficiency no check is made that n = 2**p. This check must be
!   made in the calling routine or the main program.
!  USES
!    none 
!  USED BY
!    cgca_m2red
!  SOURCE

integer(kind=idef) :: i, img, step, stepold

img = this_image()

step    = 2
stepold = 1

! do the reduction

redu: do i = 1,p

 if (mod(img,step)-1 .eq. 0) then
  sync images (img+stepold)
  coarray = coarray .and. coarray[img+stepold]
 else if (mod(img+stepold,step)-1 .eq. 0) then
  sync images (img-stepold)
 end if

 stepold = step
 step = step * 2

end do redu

! now send the result, which is in z[1] to all images.
! all images wait for image 1, so can use sync images(*),
! but, as the standard suggests, sync images is probably faster.

sync all

coarray = coarray[1]

end subroutine cgca_redand

!*roboend*

end module cgca_m2red
