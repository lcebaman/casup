!*robodoc*u* tests/test_hxi
!  NAME
!    test_hxi
!  SYNOPSIS

!$Id: test_hxi.f90 529 2018-03-26 11:25:45Z mexas $

program test_hxi

!  PURPOSE
!    Test cgca_hxi, a halo exchange subroutine.
!  DESCRIPTION
!    Assign coarrays on each image the value of this_image().
!    Then do a hx call and check that the halos are flag.
!  NOTE
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

real,parameter :: gigabyte=real(2**30)
logical( kind=ldef ), parameter :: yesdebug = .true., nodebug = .false.

real( kind=rdef ) ::   &
       qual,           & ! quality
       bsz0(3),        & ! the given "box" size
       bsz(3),         & ! updated "box" size
       dm,             & ! mean grain size, linear dim, phys units
       lres,           & ! linear resolution, cells per unit of length
       res               ! resolutions, cells per grain

integer( kind=idef ) :: ir(3), nimgs, ng, c(3) ! coarray dims

integer( kind=iarr ), allocatable :: space(:,:,:,:)[:,:,:]
integer( kind=iarr ) :: img

integer( kind=ilrg ) :: icells, mcells

integer :: flag

!*********************************************************************72
! first executable statement

! Set to an error value initially

flag = 1 
! physical dimensions of the box, assume mm
  bsz0 = (/ 1.0, 1.1, 0.9 /)

! mean grain size, linear dimension, e.g. mean grain diameter, also mm
  dm = 5.0e-1

! resolution
  res = 1.0e5

  img = int( this_image(), kind=iarr )
nimgs = num_images()

  ! do a check on image 1
  if ( img .eq. 1 ) then

     ! print a banner
     call banner("hxi")

     ! print the parameter values
     call cgca_pdmp

     write (*,'(a,i0,a)') "running on ", nimgs, " images in a 3D grid"

  end if

  ! I want pdmp output appear before the rest.
  ! This might help
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
  icells = int( c(1), kind=ilrg ) * int( c(2), kind=ilrg ) *           &
       int( c(3), kind=ilrg )

  ! total number of cells in the model
  mcells = icells * int( nimgs, kind=ilrg )

  if ( img .eq. 1 ) then
     write ( *, "(9(a,i0),tr1,g10.3,tr1,g10.3,3(a,g10.3),a)" )         &
          "img: ", img  , " nimgs: ", nimgs, " (", c(1) ,              &
          ","    , c(2) , ","       , c(3) , ")[", ir(1),              &
          ","    , ir(2), ","       , ir(3), "] ", ng   ,              &
          qual, lres,                                                  &
          " (", bsz(1), ",", bsz(2), ",", bsz(3), ")"
     write (*,'(a,i0,a)') "Each image has ",icells, " cells"
     write (*,'(a,i0,a)') "The model has ", mcells, " cells"
  end if

! allocate space coarray with a single layer
! implicit sync all
!subroutine cgca_as( l1, u1, l2, u2, l3, u3, col1, cou1, col2, cou2,   &
!                    col3, props, coarray )
call cgca_as(1, c(1), 1, c(2), 1, c(3), 1, ir(1), 1, ir(2), 1, 1, space)

space = img
sync all

call cgca_hxi( space )
sync all

! Test that hx is flag
call cgca_hxic( space, flag )

call co_sum( flag )

if ( img .eq. 1 ) then
  if ( flag .eq. 0 )  write (*,*) "hxi PASS"
end if

sync all

space = img
sync all

call cgca_hxir( space )
sync all

! Test that hx is flag
call cgca_hxic( space, flag )

call co_sum( flag )

if ( img .eq. 1 ) then
  if ( flag .eq. 0 )  write (*,*) "hxir PASS"
end if

! deallocate all arrays
call cgca_ds( space )

end program test_hxi

!*roboend*
