!*robodoc*u* tests/testACF
!  NAME
!    testACF
!  SYNOPSIS

!$Id: testACF.f90 526 2018-03-25 23:44:51Z mexas $

program testACF

  !  PURPOSE
  !    Testing MPI/IO, hdf5 and netCDF on Lustre FS
  !  DESCRIPTION
  !    Timing output of MPI/IO (cgca_pswci2) against
  !    the netcCDF version (cgca_pswci3) and hdf5
  !    (cgca_hdf5).
  !  NOTE
  !    Works only on lustre file system!
  !  AUTHOR
  !    Luis Cebamanos, Anton Shterenlikht
  !  COPYRIGHT
  !    See LICENSE
  !  USES
  !    cgca testaux
  !  USED BY
  !    Part of CGPACK test suite
  !  SOURCE

  use testaux

  implicit none

  integer, parameter :: maxlen = 64
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

  integer( kind=idef ) :: ir(3), nimgs, img, ng, c(3)  ! coarray dimensions

  integer( kind=iarr ), allocatable :: space(:,:,:,:)[:,:,:]

  integer( kind=ilrg ) :: icells, mcells

  !#################################
  character*(maxlen) :: filename
!  integer, parameter :: numiolayer = 4
!  integer, parameter :: totdim = 4, arrdim = totdim-1, coardim = 3
!  integer, parameter :: numstriping = 3
!  character*(maxlen), dimension(numstriping) :: stripestring
  character*(maxlen), dimension(3)  :: iolayername
!  integer ::  comm, ierr=0, rank=0, mpisize=0, filetype,      &
!       mpi_subarray, fh, funit, i, j, k
  integer ::  ierr=0, i, j, k, errstat

! Anton >>> I think we don't need these vars?
!  integer, dimension(totdim) :: asizehal
!  integer, dimension(arrdim) :: arrsize, arstart, artsize
!  integer, dimension(coardim) :: coarsize, copos

! Add trailing blanks to keep all elements of the array of the same
! length. Max stripe count on ARCHER is 56.
  character( len=maxlen), dimension(9) :: stripe_count = (/            &
    "-c-1 ", "-c0  ", "-c1  ", "-c4  ", "-c8  ",                       &
    "-c16 ", "-c20 ", "-c32 ", "-c40 " /)
  character( len=maxlen), dimension(7) :: stripe_size = (/             &
    "-S1m ", "-S2m ", "-S4m ", "-S8m ", "-S16m",                       &
    "-S32m", "-S64m" /)
  character( len=2*maxlen ) :: dir
  character( len=120 ) :: errmsg
 
  !#################################
  double precision :: t0, t1, tdiff, fsizeb, fsizeg
!*********************************************************************72
  ! first executable statement
  iolayername(1) = 'mpiio.dat'
  iolayername(2) = 'netcdf.dat'
  iolayername(3) = 'hdf5.dat'
! These not used yet
!  iolayername(3) = 'serial.dat'

!  stripestring(1) = 'unstriped'
!  stripestring(2) = 'striped'
!  stripestring(3) = 'defstriped'

  ! physical dimensions of the box, assume mm
  bsz0 = (/ 2.0, 2.0, 3.0 /)

  ! mean grain size, linear dimension, e.g. mean grain diameter, also mm
  dm = 1.0e-1
  !dm = 1.0e0

  ! resolution
  res = 1.0e5

    img = this_image()
  nimgs = num_images()

  ! do a check on image 1
  if ( img .eq. 1 ) then
     ! print a banner
     call banner("ACF")
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

  ! start MPI
  call MPI_Init(ierr)

  ! Loop over lfs stripe counts
  do i = 1, size( stripe_count )

    ! Loop over lfs stripe sizes
    do j = 1, size( stripe_size )

      dir = "lfs" // trim( stripe_count(i) ) // trim( stripe_size(j) )

      ! Image 1 makes a dir with desired lfs settings
      if ( img .eq. 1 ) then

write (*,*) "before mkdir"

        ! Make the dir
        errmsg = ""
        call execute_command_line( command = "mkdir " // trim(dir),    &
          wait = .true. , exitstat = ierr, cmdstat = errstat,          &
          cmdmsg = errmsg ) 

write (*,*) "after mkdir, exitstat:", ierr, "cmdstat:", errstat,       &
            "cmdmsg:", errmsg
if ( ierr .ne. 0 ) error stop

write (*,*) "before lfs setstripe"

        ! Set lfs parameters
        call execute_command_line( command = "lfs setstripe " //       &
          stripe_count(i) // " " // stripe_size(j) // " " //           &
          trim(dir) )

      end if

      ! Loop over IO layers
      do k = 1, size( iolayername )

        filename = trim(dir) // "/" // trim( iolayername(k) )

        sync all

        t0 = cgca_benchtime()

        if ( k .eq. 1 ) then
          ! MPI/IO
          call cgca_pswci2( space, cgca_state_type_grain, filename )
        else if ( k .eq. 2 ) then
          ! NetCDF
          call cgca_pswci3( space, cgca_state_type_grain, filename )
        else if ( k .eq. 3 ) then
          ! HDF5
         call cgca_pswci4( space, cgca_state_type_grain, filename )
        end if

        t1 = cgca_benchtime()

        sync all

        if (img .eq. 1)  then
          tdiff = t1 - t0
          write (*,*) trim( iolayername(k) ), " ",                     &
            trim( stripe_count(i) ), " ", trim( stripe_size(j) ), " ", &
            fsizeg/tdiff
!            tdiff, "s, rate: ", fsizeg/tdiff, "GB/s."
        end if

        sync all

      end do
    end do
  end do

  ! terminate MPI
  call MPI_Finalize(ierr)

  ! deallocate all arrays
  call cgca_ds(space)


contains

subroutine fdelete(filename)

  implicit none

  character *(*) :: filename
  integer :: stat, funit=0

  open(newunit=funit, iostat=stat, file=filename, status='old')
  if (stat.eq.0) close(unit=funit, status='delete')

end subroutine fdelete

end program testACF

!*roboend*
