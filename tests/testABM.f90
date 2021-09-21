!*robodoc*u* tests/testABM
!  NAME
!    testABM
!  SYNOPSIS

!$Id: testABM.f90 530 2018-03-26 16:10:00Z mexas $

program testABM

!  PURPOSE
!    Testing MPI/IO, cgca_m2mpiio/cgca_pswci
!  DESCRIPTION
!    Timing output of MPI/IO (cgca_pswci) against
!    the serial version (cgca_swci).
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
use mpi

implicit none

real,parameter :: gigabyte=real(2**30), resolution=1.0e-5,             &
 loge2 = log(real(2))
logical(kind=ldef),parameter :: yesdebug = .true., nodebug = .false.

real( kind=rdef ) ::    &
 qual,                  & ! quality
 bsz0(3),               & ! the given "box" size
 bsz(3),                & ! updated "box" size
 dm,                    & ! mean grain size, linear dim, phys units
 lres,                  & ! linear resolution, cells per unit of length
 res                      ! resolutions, cells per grain

integer :: c(3), ir(3), nimgs, img, ng
integer( kind=iarr ), allocatable :: space(:,:,:,:)[:,:,:]

integer( kind=ilrg ) :: icells, mcells

integer :: ierr

real :: time1, time2, fsizeb, fsizeg, tdiff

logical :: iflag

!*********************************************************************72
! first executable statement

! physical dimensions of the box, assume mm
bsz0 = (/ 2.0, 3.0, 1.0 /)

! mean grain size, linear dimension, e.g. mean grain diameter, also mm
dm = 1.0e-1

! resolution
res = 1.0e5

    img = this_image()
nimgs = num_images()

! do a check on image 1
if ( img .eq. 1 ) then
 ! print a banner
 call banner("ABM")
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

! Total output file size, in B and in GB.
fsizeb = real( mcells * storage_size( space, kind=ilrg ) / 8_ilrg )
fsizeg = fsizeb / gigabyte

! allocate space coarray with a single layer
! implicit sync all
!subroutine cgca_as( l1, u1, l2, u2, l3, u3, col1, cou1, col2, cou2,   &
!                    col3, props, coarray )
call cgca_as(1, c(1), 1, c(2), 1, c(3), 1, ir(1), 1, ir(2), 1, 1, space)

! initialise coarray to image number
space = int( img, kind=iarr )

! dump the model, serial
call cpu_time( time1 )
call cgca_swci( space, cgca_state_type_grain, 10, 'serial.raw' )
call cpu_time( time2 )
tdiff = time2-time1
if (img .eq. 1)                                                        &
  write (*,*) "Serial IO: ", tdiff, "s, rate: ", fsizeg/tdiff, "GB/s."

! All images start MPI, used only for I/O here, so no need to start
! earlier. Note that OpenCoarrays doesn't like MPI_Init.
! Probably it's called automatically by the runtime, so
! I put a check for it.
call MPI_Initialized( iflag, ierr )
if ( .not. iflag ) call MPI_Init(ierr)

! dump the model, MPI/IO
call cpu_time( time1 )
call cgca_pswci( space, cgca_state_type_grain, 'mpiio.raw' )
call cpu_time( time2 )
tdiff = time2-time1
if (img .eq. 1)                                                        &
  write (*,*) "MPI/IO: ", tdiff, "s, rate: ", fsizeg/tdiff, "GB/s."

! terminate MPI
! See the note above on MPI_Init
!call MPI_Finalized( iflag, ierr)
!if ( .not. iflag ) call MPI_Finalize(ierr)

! deallocate all arrays
call cgca_ds(space)

end program testABM

!*roboend*
