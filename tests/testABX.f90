!*robodoc*u* tests/testABX
!  NAME
!    testABX
!  SYNOPSIS

!$Id: testABX.f90 526 2018-03-25 23:44:51Z mexas $

program testABX

!  PURPOSE
!    Testing Cray parallel direct access file IO, aka "assign"
!    environment, routine cgca_m2out/cgca_pc.
!  NOTE
!    Compile only on Cray! Don't try to compile with other compilers!
!  DESCRIPTION
!    Verify data integrity and compare timings of
!    cgca_swci and cgca_pc - i.e. a serial writer from image 1
!    and a parallel direct access shared writer. The latter is
!    a non standard Cray extension. So don't try to run this
!    on non-Cray machines.
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

real, parameter :: gigabyte=real(2**30)
logical( kind=ldef ), parameter :: yesdebug = .true., nodebug = .false.

real :: time1, time2, fsizeb, fsizeg, tdiff

integer( kind=idef ) :: ir(3), img, nimgs,                             &
  ng,   & ! number of grains in the whole model
  c(3)    ! coarray dimensions
integer( kind=iarr ), allocatable :: space(:,:,:,:) [:,:,:]

integer( kind=ilrg ) :: icells, mcells

real( kind=rdef ) ::    &
 qual,                  & ! quality
 bsz0(3),               & ! the given "box" size
 bsz(3),                & ! updated "box" size
 dm,                    & ! mean grain size, linear dim, phys units
 lres,                  & ! linear resolution, cells per unit of length
 res                      ! resolutions, cells per grain

logical( kind=ldef ) :: solid

!*********************************************************************72
! first executable statement

! physical dimensions of the box, assume mm
bsz0 = (/ 1.0, 2.0, 3.0 /)

! mean grain size, linear dimension, e.g. mean grain diameter, also mm
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
  call banner("ABX")
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
  write ( *, "(9(a,i0),tr1,g10.3,tr1,i0,3(a,g10.3),a)" )               &
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

! initialise random number seed
call cgca_irs( nodebug )

! initialise space
space( :, :, :, cgca_state_type_grain ) = cgca_liquid_state

! nuclei, sync all inside
call cgca_nr( space, ng, nodebug )

! solidify with co_sum, implicit sync all inside
! subroutine cgca_sld3( coarray, iter, heartbeat, solid )
call cgca_sld3( space, 0, 10, solid )

! dump the model
call cpu_time(time1)
call cgca_pc( space, cgca_state_type_grain, 'craypario.raw' )
call cpu_time(time2)
tdiff = time2-time1
if (img .eq. 1) then
 write (*,*) "File size:", fsizeg, "GB"
 write (*,*) "Cray assign -m on IO: ", tdiff, "s, rate: ",             &
             fsizeg/tdiff, "GB/s."
end if

call cpu_time(time1)
call cgca_swci( space, cgca_state_type_grain, 10, 'serial.raw' )
call cpu_time(time2)
tdiff = time2-time1
if (img .eq. 1) write (*,*) "Serial IO: ", tdiff, "s, rate: ",         &
   fsizeg/tdiff, "GB/s."

! deallocate all arrays
call cgca_ds( space )

end program testABX

!*roboend*
