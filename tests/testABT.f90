!*robodoc*u* tests/testABT
!  NAME
!    testABT
!  SYNOPSIS

!$Id: testABT.f90 526 2018-03-25 23:44:51Z mexas $

program testABT

!  PURPOSE
!    Checking: cgca_gdim, cgca_cadim
!  DESCRIPTION
!    cgca_gdim finds the optimum coarray grid layout for a
!    given total number of images. It also reports the
!    quality of this optimum, from 0 - worst, to 1 - best.
!    cgca_cadim then calculates the coarray dimensions
!    the new updated box size.
!  NOTE
!    Both cgca_gdim and cgca_cadim are serial routines.
!    It makes no sence to run this test at high numbers of images.
!    A single image is enough to test the routines.
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

integer( kind=idef ) :: n, ir(3), nimgs, &
  ng,     & ! number of grains in the whole model
  c(3)      ! coarray dimensions
logical( kind=ldef ) :: image1
real( kind=rdef ) ::    &
 qual,                  & ! quality
 bsz0(3),               & ! the given "box" size
 bsz(3),                & ! updated "box" size
 dm,                    & ! mean grain size, linear dim, phys units
 res,                   & ! resolutions, cells per grain
! tmprnd(3),             & ! array of random numbers
 lres                     ! linear resolution

!*********************************************************************72
! first executable statement

! physical dimensions of the box, assume mm
bsz0 = (/ 10.0, 20.0, 30.0 /)

! mean grain size, also mm
dm = 5.0e-1

! resolution
res = 1.0e5 

 nimgs = num_images()
image1 = .false.
if (this_image() .eq. 1) image1 = .true.

! do a check on image 1
if (image1) then

 ! print a banner
 call banner("ABT")

 ! print the parameter values
 call cgca_pdmp
 write (*,'(a,i0,a)') "running on ", nimgs, " images in a 3D grid"

 ! calculate the coarray grid dimensions
 do n = 1, 2**15
   call cgca_gdim( n, ir, qual )

!   ! choose box sizes at random, max 30 in any dimension
!   call random_number( tmprnd )
!   bsz0 = tmprnd * 30.0

   ! subroutine cgca_cadim( bsz, res, dm, ir, c, lres, ng )
   bsz = bsz0
   call cgca_cadim( bsz, res, dm, ir, c, lres, ng )

   write ( *, "(8(i0,a),es9.2,tr1,es9.2,3(a,es9.2),a)" )               &
    n, "(", c(1), ",", c(2), ",", c(3), ")[" ,                         &
    ir(1), ",", ir(2), ",", ir(3), "] ", ng, " ",                      &
    qual, lres,                                                        &
    " (", bsz(1), ",", bsz(2), ",", bsz(3), ")"
 end do

end if

sync all

end program testABT

!*roboend*
