!*robodoc*u* tests/test_hxir
!  NAME
!    test_hxir
!  SYNOPSIS

!$Id: test_hxir.f90 526 2018-03-25 23:44:51Z mexas $

program test_hxir

!  PURPOSE
!    Test cgca_hxir, a halo exchange subroutine with random order of
!    remote calls.
!  DESCRIPTION
!    Two separate solidifications are performed, which should result
!    in bit-by-bit idendical microstructures.
!    Both solidification runs use m3sld_hc/cgca_sld_h, a routine
!    that offers a choice of different halo exchange routines.
!    One solidification is using cgca_hxi, a halo exchange routine with
!    ordered, predictable sequence of remote calls, the same
!    on all images.
!    Another solidification is using cgca_hxir, a halo exchange
!    routine with random order of remote calls.
!    MPI/IO is used for speed. 
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

real( kind=rdef ), allocatable :: grt(:,:,:)[:,:,:]

integer( kind=idef ) :: ir(3), nimgs, img, ng, i, ierr, funit, c(3)

integer( kind=iarr ), allocatable :: space(:,:,:,:)[:,:,:]

integer( kind=ilrg ) :: icells, mcells

logical( kind=ldef ) :: solid

character(len=5) :: fnum

double precision :: t0, t1, tdiff, fsizeb, fsizeg

!*********************************************************************72
! first executable statement

! physical dimensions of the box, assume mm
  bsz0 = (/ 1.0, 1.1, 0.9 /)

! mean grain size, linear dimension, e.g. mean grain diameter, also mm
  dm = 5.0e-1

! resolution
  res = 1.0e5

  img = this_image()
nimgs = num_images()

  ! do a check on image 1
  if ( img .eq. 1 ) then

     ! print a banner
     call banner("hxir")

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
     write (*,*) "dataset sizes for ParaView", c*ir
  end if

  ! Total output file size, in B and in GB.
  fsizeb = real( mcells * storage_size( space, kind=ilrg ) / 8_ilrg )
  fsizeg = fsizeb / gigabyte

! allocate space coarray with a single layer
! implicit sync all
!subroutine cgca_as( l1, u1, l2, u2, l3, u3, col1, cou1, col2, cou2,   &
!                    col3, props, coarray )
call cgca_as(1, c(1), 1, c(2), 1, c(3), 1, ir(1), 1, ir(2), 1, 1, space)

! Start MPI
call MPI_Init( ierr )

! Do 2 solidifications
do i=1,2

  ! initialise the reproducible random number seed
  call cgca_ins( yesdebug )

  ! initialise space
  space = cgca_liquid_state
!  space = img

  ! allocate rotation tensors
  call cgca_art( 1, ng, 1, ir(1), 1, ir(2), 1, grt )

  ! assign rotation tensors, sync all inside
  call cgca_rt( grt )

  ! nuclei, sync all inside
  call cgca_nr( space, ng, yesdebug )

  ! Start timer
  t0 = cgca_benchtime()

  ! Solidify, implicit sync all inside
  !  module subroutine cgca_sld_h( coarray, hx, iter, heartbeat, solid )

  if      ( i .eq. 1 ) then
    call cgca_sld_h( space, cgca_hxi,  2, 1, solid )
  else if ( i .eq. 2 ) then
    call cgca_sld_h( space, cgca_hxir, 2, 1, solid )
  end if

  t1 = cgca_benchtime()

  if (img .eq. 1)  then
    tdiff = t1 - t0
    if      ( i .eq. 1 ) then
      write (*,*) "cgca_hxi  time: ", tdiff, "s"
    else if ( i .eq. 2 ) then
      write (*,*) "cgca_hxir time: ", tdiff, "s"
    end if
  end if

  sync all

  ! MPI/IO
  ! subroutine cgca_pswci2( coarray, stype, fname )
  if      ( i .eq. 1 ) then
    call cgca_pswci2( space, cgca_state_type_grain, "hxi.raw"  )
! Plain text output
write (fnum, '(i0)') img
open( newunit=funit, file="hxi_img" // fnum, status="replace" )
write (funit, *) space
close( funit )

  else if ( i .eq. 2 ) then
    call cgca_pswci2( space, cgca_state_type_grain, "hxir.raw" )
! Plain text output
write (fnum, '(i0)') img
open( newunit=funit, file="hxir_img" // fnum, status="replace" )
write (funit, *) space
close( funit )

  end if

  sync all

end do

! terminate MPI
call MPI_Finalize( ierr )

! deallocate all arrays
call cgca_ds( space )

end program test_hxir

!*roboend*
