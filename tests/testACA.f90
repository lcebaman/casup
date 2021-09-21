!*robodoc*u* tests/testACA
!  NAME
!    testACA
!  SYNOPSIS

!$Id: testACA.f90 526 2018-03-25 23:44:51Z mexas $

program testACA

!  PURPOSE
!    Testing cgca_m2gl/cgca_ico.
!  DESCRIPTION
!    cgca_ico converts some image index into its cosubscripts.
!    Lots of data has to be set prior to calling cgca_ico. 
!    The test can be run on any number of images.
!    space coarray is established and the cgca_ico is called
!    using its cosubscripts.
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

integer( kind=idef ) :: ir( cgca_scodim ), img, nimgs, flag, ind,      &
 ng,      & ! number of grains in the whole model
 cosub( cgca_scodim ), c( cgca_scodim )  ! coarray dimensions

integer( kind=iarr ), allocatable :: space(:,:,:,:) [:,:,:]

integer( kind=ilrg ) :: icells, mcells

real( kind=rdef ) ::    &
 qual,                  & ! quality
 bsz0(3),               & ! the given "box" size
 bsz(3),                & ! updated "box" size
 dm,                    & ! mean grain size, linear dim, phys units
 lres,                  & ! linear resolution, cells per unit of length
 res                      ! resolutions, cells per grain

!*********************************************************************72
! first executable statement

! physical dimensions of the box, assume mm
bsz0 = (/ 1.0, 2.0, 3.0 /)

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
  call banner("ACA")
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
icells = int( c(1), kind=ilrg ) * int( c(2), kind=ilrg ) *             &
         int( c(3), kind=ilrg )

! total number of cells in the model
mcells = icells * int( nimgs, kind=ilrg )

if ( img .eq. 1 ) then
  write ( *, "(9(a,i0),tr1,g10.3,tr1,g10.3,3(a,g10.3),a)" )            &
    "img: ", img  , " nimgs: ", nimgs, " (", c(1) ,                    &
    ","    , c(2) , ","       , c(3) , ")[", ir(1),                    &
    ","    , ir(2), ","       , ir(3), "] ", ng   ,                    &
    qual, lres,                                                        &
    " (", bsz(1), ",", bsz(2), ",", bsz(3), ")"
  write (*,'(a,i0,a)') "Each image has ",icells, " cells"
  write (*,'(a,i0,a)') "The model has ", mcells, " cells"
end if

! allocate space coarray with a single layer
! implicit sync all
!subroutine cgca_as( l1, u1, l2, u2, l3, u3, col1, cou1, col2, cou2,   &
!                    col3, props, coarray )
call cgca_as(1, c(1), 1, c(2), 1, c(3), 1, ir(1), 1, ir(2), 1, 1, space)

! check cgca_ico on the last image
! make only image 1 do this.
if ( img .eq. 1 ) then

  ! Dump lower and upper cobounds
  write (*,*) "cgca_slcob:", cgca_slcob
  write (*,*) "cgca_sucob:", cgca_sucob

  ! Calculate subscripts for all indices and check
  do ind = 1, num_images()
    call cgca_ico( ind, cosub, flag )
    if ( flag .ne. 0 ) then
      write (*,*) "ERROR: cgca_ico returned flag:", flag
      stop
    end if

    write (*,*) "index=", ind, "cosub:", cosub

    ! check
    if ( image_index( space, cosub ) .ne. ind ) then
      write (*,*) "ERROR: testACA failed: cgca_ico calculated wrong values"
    end if 

  end do

end if

sync all

! deallocate all arrays
call cgca_ds( space )

end program testACA

!*roboend*
