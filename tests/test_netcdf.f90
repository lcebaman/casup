!*robodoc*u* tests/test_netcdf
!  NAME
!    test_netcdf
!  SYNOPSIS

!$Id: test_netcdf.f90 533 2018-03-30 14:31:26Z mexas $

program test_netcdf

!  PURPOSE
!    Test NetCDF IO timing and integrity against serial IO.
!  DESCRIPTION
!    Serial IO is cgca_swci. NetCDF IO is cgca_pswci3.
!  NOTE
!    Should work on any file system, as long as NetCDF libraries
!    are properly installed.
!  AUTHOR
!    Luis Cebamanos, Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  USES
!    cgca testaux
!  USED BY
!    Part of CGPACK test suite
!  SOURCE

use casup

implicit none

real,parameter :: gigabyte=real(2**30)
logical( kind=ldef ), parameter :: yesdebug = .true., nodebug = .false.

  real( kind=rdef ) ::  &
       qual,            & ! quality
       bsz0(3),         & ! the given "box" size
       bsz(3),          & ! updated "box" size
       dm,              & ! mean grain size, linear dim, phys units
       lres,            & ! linear resolution, cells per unit of length
       res                ! resolutions, cells per grain

  integer( kind=idef ) :: ir(3), nimgs, img, ng, c(3) ! coarray dims

  integer( kind=iarr ), allocatable :: space(:,:,:,:) [:,:,:]

  integer( kind=ilrg ) :: icells, mcells

  double precision :: t0, t1, tdiff, fsizeb, fsizeg

!*********************************************************************72
! first executable statement

! physical dimensions of the box, assume mm
  bsz0 = (/ 1.0, 1.5, 2.0 /)

! mean grain size, linear dimension, e.g. mean grain diameter, also mm
  dm = 3.0e-1

! resolution
  res = 1.0e5

  img = this_image()
nimgs = num_images()

  ! do a check on image 1
  if ( img .eq. 1 ) then

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

sync all

! start MPI
!  call MPI_Init(ierr)

! serial
! subroutine cgca_swci( coarray, stype, iounit, fname )
t0 = cgca_benchtime()
call cgca_swci( space, cgca_state_type_grain, 10, "seri.raw" )
t1 = cgca_benchtime()

sync all

if (img .eq. 1)  then
  tdiff = t1 - t0
  write (*,*) "Serial time: ", tdiff, "s. Rate: ", fsizeg/tdiff, "GB/s."
end if

! netcdf
! subroutine cgca_pswci3( coarray, stype, fname )
t0 = cgca_benchtime()
call cgca_pswci3( space, cgca_state_type_grain, "netcdf.ncdf" )
t1 = cgca_benchtime()

sync all

if (img .eq. 1)  then
  tdiff = t1 - t0
  write (*,*) "NetCDF time:", tdiff, "s, rate: ", fsizeg/tdiff, "GB/s."
end if

! terminate MPI
! call MPI_Finalize(ierr)

  ! deallocate all arrays
  call cgca_ds(space)

end program test_netcdf

!*roboend*
