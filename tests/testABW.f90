!*robodoc*u* tests/testABW
!  NAME
!    testABW
!  SYNOPSIS

!$Id: testABW.f90 526 2018-03-25 23:44:51Z mexas $

program testABW

!  PURPOSE
!    Check cgca_clvgp_nocosum (no CO_SUM version).
!    Cleavage propagation in a model with the given "box" size,
!    mean grain size and the microstructure resolution values.
!    Using parallel I/O here!
!  DESCRIPTION
!    First need to call cgca_gdim, cgca_cadim to
!    calculate all parameters of coarray space, including
!    the number of nuclei. Then call solidification and then
!    cleavage.
!  NOTE
!    Using 2 cleavage nuclei to see how the two macro-cracks interact.
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

! specify the total number of cleavage propagation iterations,
! and the number of times the fracture array will be dumped along
! the way
integer( kind=idef ), parameter ::                                     &
 itot = 150_idef,       & ! total cleavage iterations to do
 idmp = 3_idef           ! number of dumps to do along the way 

real, parameter :: gigabyte=real(2**30),                               &
! cleavage stress on 100, 110, 111 planes for BCC,
! see the manual for derivation.
 scrit(3) = (/ 1.05e4, 1.25e4, 4.90e4 /)

!integer :: ierr ! default integer needed for MPI

integer( kind=idef ) :: iter, fiter, ir(3), img, nimgs,                &
  ng,    & ! number of grains in the whole model
  c(3)     ! coarray dimensions
integer( kind=iarr ), allocatable :: space(:,:,:,:) [:,:,:]

real( kind=rdef ) ::    &
 t(3,3),                & ! stress tensor
 qual,                  & ! quality
 bsz0(3),               & ! the given "box" size
 bsz(3),                & ! updated "box" size
 dm,                    & ! mean grain size, linear dim, phys units
 res,                   & ! resolutions, cells per grain
 icells,                & ! total number of cells on each image
 lres                     ! linear resolution, cells per unit of length
real( kind=rdef ), allocatable :: grt(:,:,:)[:,:,:]

logical( kind=ldef ) :: solid
character(6) :: citer

!*********************************************************************72
! first executable statement

! physical dimensions of the box, assume mm
bsz0 = (/ 1.0, 1.0, 1.0 /)

! mean grain size, linear dimension, e.g. mean grain diameter, also mm
dm = 3.0e-1

! resolution
res = 1.0e5

! In this test set the number of images via the env var
! the code must be able to cope with any value >= 1.
   img = this_image()
 nimgs = num_images()

! do a check on image 1
if ( img .eq. 1 ) then

 ! print a banner
 call banner("ABW")

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
icells = product( real( c ) )

! dump some stats from image 1
if ( img .eq. 1 ) then
  write ( *, "(9(a,i0), 2(tr1, es9.2), 3(a, es9.2), a)" )              &
    "img: ", img  , " nimgs: ", nimgs, " (", c(1) ,                    &
    ","    , c(2) , ","       , c(3) , ")[", ir(1),                    &
    ","    , ir(2), ","       , ir(3), "] ", ng   ,                    &
    qual, lres,                                                        &
    " (", bsz(1), ",", bsz(2), ",", bsz(3), ")"
  write (*,'(a,es9.2,a)') "Cells on each image: ", icells
  write (*,*) "dataset sizes for ParaView", c*ir
  write (*,"(a, es10.2, a, i0)") "Total cells in the model (real): ",  &
    product( real(c) * real(ir) ), " (int): ",                         &
    product( int(c, kind=ilrg) * int(ir, kind=ilrg) )
end if

! allocate space coarray with two layers
! implicit sync all
call cgca_as(1, c(1), 1, c(2), 1, c(3), 1, ir(1), 1, ir(2), 1, 2, space)

! initialise random number seed
call cgca_irs( nodebug )

! allocate rotation tensors
call cgca_art( 1, ng, 1, ir(1), 1, ir(2), 1, grt )

! initialise space
space( :, :, :, cgca_state_type_grain ) = cgca_liquid_state
space( :, :, :, cgca_state_type_frac  ) = cgca_intact_state

! assign rotation tensors, sync all inside
call cgca_rt( grt )

! nuclei, sync all inside
call cgca_nr( space, ng, nodebug )

! solidify, implicit sync all inside
call cgca_sld( space, noperiodicbc, 0, 10, solid )

! initiate grain boundaries
call cgca_igb( space )
sync all
call cgca_hxi( space )

! smoothen the GB, several iterations, sync needed
! halo exchange, following smoothing
call cgca_gbs( space )
sync all
call cgca_hxi( space )
sync all
call cgca_gbs( space )
sync all
call cgca_hxi( space )
sync all

! update grain connectivity, local routine, no sync needed
call cgca_gcu( space )

! set a crack nucleus in the centre of the x3=max(x3) face
space( c(1)/2, c(2)/2, c(3), cgca_state_type_frac )  &
        [ ir(1)/2+1, ir(2)/2+1, ir(3) ] = cgca_clvg_state_100_edge

! set a crack nucleus in the centre of the x3=1 face
space( c(1)/2, c(2)/2, 1, cgca_state_type_frac )  &
        [ ir(1)/2+1, ir(2)/2+1, 1 ] = cgca_clvg_state_100_edge

sync all

! all images start MPI, used only for I/O here, so
! no need to start earlier
!call MPI_Init(ierr)

if ( img .eq. 1 ) write (*,*) "dumping model to files"

  ! dump space arrays to files, only image 1 does it, all others
  ! wait at the barrier, hence sync needed
  call cgca_pswci( space, cgca_state_type_grain, "zg0.raw" )
  call cgca_pswci( space, cgca_state_type_frac,  "zf0.raw" )

if ( img .eq. 1) write (*,*) "finished dumping model to files"

sync all

! set the stress tensor
t = 0.0
t(1,1) = 1.0e6
!t(2,2) = -1.0e6

! number of cleavage iterations between file dumps
fiter = itot/idmp

! propagate cleavage, sync inside, dump fracture arrays
! every certain number of increments.
do iter = 1, idmp

 ! Propagate cleavage, sync inside
 ! Run for fiter iterations
 ! subroutine cgca_clvgp_nocosum( coarray, rt, t, scrit, sub,
 !    gcus, periodicbc, iter, heartbeat, debug )
 call cgca_clvgp_nocosum( space, grt, t, scrit, cgca_clvgsd,           &
      cgca_gcupdn, noperiodicbc, fiter, 10, yesdebug )

 ! dump the fracture space array to files, only image 1 does it,
 ! all others wait at the barrier, hence sync needed.
 ! citer is the total number of cleavage fracture iterations,
 ! "c" for character date type.
 write ( citer, "(i0)" ) iter*fiter
 call cgca_pswci( space, cgca_state_type_frac,                         &
                                "zf"//trim(citer)//".raw" )

 sync all

 if ( img .eq. 1 ) write (*,"(a)") &
  "Completed "//trim(citer)//" cleavage iterations"

end do

! deallocate all arrays, implicit sync all.
call cgca_ds( space )
call cgca_dgc
call cgca_drt( grt )

! terminate MPI
!call MPI_Finalize(ierr)

if ( img .eq. 1 ) write (*,*) "Test ABW completed sucessfully"

end program testABW

!*roboend*
