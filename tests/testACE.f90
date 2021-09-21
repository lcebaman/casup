!*robodoc*u* tests/testACE
!  NAME
!    testACE
!  SYNOPSIS

!$Id: testACE.f90 526 2018-03-25 23:44:51Z mexas $

program testACE

!  PURPOSE
!    Test cgca_fwci/cgca_m2out - a debugging output routine.
!  DESCRIPTION
!    cgca_fwci/cgca_m2out collects coarray data from all images and
!    writes it out from image 1 in a formatted text file. A single
!    line contains the state of a single cells with its full local
!    global coordinates, i.e. its position in the coarray on this
!    image, and the image position in the coarray grid. 
!    The data is dumped out as soon as the solidificaiton is finished.
!
!    Both grain and fracture layers are used and dumped.
!  NOTE
!    The program uses cgca_ins RND seed, to obtain a reproducible
!    model data. If this test is re-run on the same platform with
!    the same number of images, the results must be the same.
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

logical( kind=ldef ), parameter :: yesdebug=.true., nodebug=.false.,   &
 periodicbc=.true.,  noperiodicbc=.false.

real, parameter :: gigabyte=real(2**30)

integer( kind=idef ) :: ir(3), img, nimgs,                             &
  ng,  & ! number of grains in the whole model
  c(3)   ! coarray dimensions

integer( kind=iarr ), allocatable :: space(:,:,:,:) [:,:,:]

integer( kind=ilrg ) :: icells, mcells

real( kind=rdef ) ::    &
 qual,                  & ! quality
 bsz0(3),               & ! the given "box" size
 bsz(3),                & ! updated "box" size
 dm,                    & ! mean grain size, linear dim, phys units
 lres,                  & ! linear resolution, cells per unit of length
 res                      ! resolutions, cells per grain
real( kind=rdef ), allocatable :: grt(:,:,:)[:,:,:]

logical( kind=ldef ) :: solid

!*********************************************************************72
! first executable statement

! physical dimensions of the box, assume mm
bsz0 = (/ 1.0, 1.0, 1.0 /)

! mean grain size, linear dimension, e.g. mean grain diameter, also mm
dm = 5.0e-1

! resolution
res = 1.0e5

! In this test set the number of images via the env var
! the code must be able to cope with any value >= 1.
   img = this_image()
 nimgs = num_images()

! do a check on image 1
if ( img .eq. 1 ) then

 ! print a banner
 call banner("ACE")

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
! c - coarray sizes
! ir - coarray grid sizes
bsz = bsz0
call cgca_cadim( bsz, res, dm, ir, c, lres, ng )

! total number of cells in a coarray
icells = int( c(1), kind=ilrg ) *    &
         int( c(2), kind=ilrg ) *    &
         int( c(3), kind=ilrg )

! total number of cells in the model
mcells = icells * int( nimgs, kind=ilrg )

if ( img .eq. 1 ) then
  write ( *, "(9(a,i0), 2(tr1, es9.2), 3(a, es9.2), a)" )              &
    "img: ", img  , " nimgs: ", nimgs, " (", c(1) ,                    &
    ","    , c(2) , ","       , c(3) , ")[", ir(1),                    &
    ","    , ir(2), ","       , ir(3), "] ", ng   ,                    &
    qual, lres,                                                        &
    " (", bsz(1), ",", bsz(2), ",", bsz(3), ")"
  write (*,*) "dataset sizes for ParaView", c*ir
  write (*,'(a,i0)') "Cells on each image: ", icells
  write (*,"(a, es10.2, a, i0)") "Total cells in the model (real): ",  &
    product( real(c) * real(ir) ), " (int): ", mcells
end if

! allocate space coarray with two layers, implicit SYNC ALL inside
call cgca_as(1, c(1), 1, c(2), 1, c(3), 1, ir(1), 1, ir(2), 1, 2, space)

! initialise the reproducible random number seed
call cgca_ins( yesdebug )

! Set initial values to all layers of space
space( :, :, :, cgca_state_type_grain ) = cgca_liquid_state
space( :, :, :, cgca_state_type_frac  ) = cgca_intact_state

! Allocate rotation tensors, implicit SYNC ALL inside
call cgca_art( 1, ng, 1, ir(1), 1, ir(2), 1, grt )

! Image 1 sets crack nuclei
if ( img .eq. 1 ) then

  ! set a single crack nucleus in the centre of the x3=max(x3) face
  space( c(1)/2, c(2)/2, c(3), cgca_state_type_frac )                  &
        [ ir(1)/2, ir(2)/2, ir(3) ] = cgca_clvg_state_100_edge
end if

! Set grain nuclei, SYNC ALL inside
call cgca_nr( space, ng, yesdebug )

! assign rotation tensors, SYNC ALL inside
call cgca_rt( grt )

! solidify, implicit SYNC ALL inside
!subroutine cgca_sld( coarray, periodicbc, iter, heartbeat, solid )
call cgca_sld( space, noperiodicbc, 0, 10, solid )

! Initiate grain boundaries. cgca_igb has no remote comms. Halo
! exchange is needed to update other images.
call cgca_igb( space )
sync all
call cgca_hxi( space )
sync all

! Smoothen the GB, several iterations.
! cgca_gbs has no remote comms.
call cgca_gbs( space )
sync all
call cgca_hxi( space )
sync all
call cgca_gbs( space )
sync all
call cgca_hxi( space )
sync all

if ( img .eq. 1 ) write (*,*) "dumping model to files"

  ! dump space arrays to files, only image 1 does it, all others
  ! wait at the barrier, hence sync needed
  call cgca_fwci( space, cgca_state_type_grain, "zg0.raw" )
  call cgca_fwci( space, cgca_state_type_grain, "zf0.raw" )

if ( img .eq. 1) write (*,*) "finished dumping model to files"

sync all

! deallocate all coaarrays, implicit sync all.
call cgca_ds( space )
call cgca_drt( grt )

if ( img .eq. 1 ) write (*,*) "Test ACE completed sucessfully"

end program testACE

!*roboend*
