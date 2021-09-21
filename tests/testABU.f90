!*robodoc*u* tests/testABU
!  NAME
!    testABU
!  SYNOPSIS

!$Id: testABU.f90 526 2018-03-25 23:44:51Z mexas $

program testABU

!  PURPOSE
!    Checking: cgca_imco
!  DESCRIPTION
!    First need to call cgca_gdim, cgca_cadim to
!    calculate all parameters of coarray space.
!  NOTE
!    cgca_gdim and cgca_cadim can be called by any or all images.
!    Their results do not depend on the index of the invoking image.
!    However, cgca_imco must be called by every image,
!    So it makes sense to call all three routines by every image.
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

integer( kind=idef ) :: ir(3), img, nimgs,                             &
  ng,    & ! number of grains in the whole model
  c(3)     ! coarray dimensions
integer( kind=iarr ), allocatable :: space(:,:,:,:) [:,:,:]

real( kind=rdef ), parameter :: zero = 0.0_rdef, one = 1.0_rdef
real( kind=rdef ) ::    &
 lres,                  & ! linear resolution
 qual,                  & ! quality
 bsz0(3),               & ! the given "box" size
 bsz(3),                & ! updated "box" size
 origin(3),             & ! origin of the "box" cs, in FE cs
 rot(3,3),              & ! rotation tensor *from* FE cs *to* CA cs
 dm,                    & ! mean grain size, linear dim, phys units
 res,                   & ! resolutions, cells per grain
 bcol(3), bcou(3)         ! lower and upper phys. coord of the coarray
                          ! on each image

!*********************************************************************72
! first executable statement

! physical dimensions of the box, assume mm
bsz0 = (/ 1.0, 2.0, 3.0 /)

! origin of the box cs, assume mm
origin = (/ 10.0, 11.0, 12.0 /)

! rotation tensor *from* FE cs *to* CA cs.
! The box cs is aligned with the box.
rot = zero
rot(1,1) = one
rot(2,2) = one
rot(3,3) = one

! mean grain size, also mm
dm = 1.0e-1

! resolution
res = 1.0e5 

! In this test set the number of images via the env var
! the code must be able to cope with any value >= 1.
   img = this_image()
 nimgs = num_images()

! do a check on image 1
if ( img .eq. 1 ) then

 ! print a banner
 call banner("ABU")

 ! print the parameter values
 call cgca_pdmp
 write (*,'(a,i0,a)') "running on ", nimgs, " images in a 3D grid"

end if

! want to sync here to make sure the banner is
! printed before the rest.
sync all

! each image calculates the coarray grid dimensions
call cgca_gdim( nimgs, ir, qual )

! calculate the resolution and the actual phys dimensions
! of the box
! subroutine cgca_cadim( bsz, res, dm, ir, c, lres, ng )
bsz = bsz0
call cgca_cadim( bsz, res, dm, ir, c, lres, ng )

write ( *, "(9(a,i0),tr1,g10.3,tr1,i0,3(a,g10.3),a)" )           &
    "img: ", img, " nimgs: ", nimgs,                             &
     " (", c(1), ",", c(2), ",", c(3),                           &
     ")[", ir(1), ",", ir(2), ",", ir(3), "] ", ng,              &
    qual, lres,                                                  &
    " (", bsz(1), ",", bsz(2), ",", bsz(3), ")"

! allocate space coarray with a single layer
call cgca_as(1, c(1), 1, c(2), 1, c(3), 1, ir(1), 1, ir(2), 1, 1,space)

! calculate the extremeties of the box, in the CA cs, on each image
!subroutine cgca_imco( space, lres, bcol, bcou )
call cgca_imco( space, lres, bcol, bcou ) 

write ( *,"(a,i0,2(a,3(g10.3,tr1)),a)" ) "img: ", img,           &
  " CA bcol: (", bcol, ") CA bcou: (", bcou, ")"

! and now in FE cs:
write ( *,"(a,i0,2(a,3(g10.3,tr1)),a)" ) "img: ", img,           &
   " FE bcol: (", matmul( transpose( rot ),bcol ) + origin,      &
  ") FE bcou: (", matmul( transpose( rot ),bcou ) + origin, ")"

! deallocate space
call cgca_ds( space )

end program testABU


!*roboend*
