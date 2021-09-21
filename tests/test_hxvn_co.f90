!*robodoc*u* tests/test_hxvn_co
!  NAME
!    test_hxvn_co
!  SYNOPSIS

!$Id: test_hxvn_co.f90 537 2018-04-03 13:57:55Z mexas $

program test_hxvn_co

!  PURPOSE
!    Test HX routine ca_co_hx_all.
!    Use coarrays for the whole model, not just halos.
!  DESCRIPTION
!    - ca_co_spalloc - user allocates the space coarray, obviously
!      at the start of the simulation.
!    - ca_co_hx_all - a high level routine to do all necessary HX
!      operations, with necessary sync.
!    Must work on any number of images, except when a good
!    decomposition cannot be made.
!    The user needs know nothing about sync.
!  NOTE
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  USES
!    casup
!  USED BY
!    Part of CASUP test suite
!  SOURCE

use casup

implicit none

integer( kind=iarr ), parameter :: huge_iarr = huge(0_iarr)

real( kind=rdef ) :: &
  qual,              & ! quality
  bsz0(3),           & ! the given "box" size
  bsz(3),            & ! updated "box" size
  dm,                & ! mean grain size, linear dim, phys units
  lres,              & ! linear resolution, cells per unit of length
  res                  ! resolutions, cells per grain

integer( kind=iarr ), allocatable :: space(:,:,:) [:,:,:],             &
  space1(:,:,:) [:,:,:]

integer( kind=idef ) :: ir(3), nimgs, img, ng, c(3) ! coarray dimensions

integer( kind=ilrg ) :: icells, mcells

integer :: ierr, d, run

! logical :: flag

!*********************************************************************72
! first executable statement

!bsz0 = (/ 4.0e2, 8.0e2, 6.0e2 /) ! numbers of cells in CA space
bsz0 = (/ 4.0e1, 8.0e1, 6.0e1 /) ! for testing on FreeBSD laptop
  dm = 1.0 ! cell size
 res = 1.0 ! resolution

  img = this_image()
nimgs = num_images()

! do a check on image 1
if ( img .eq. 1 ) then
  write (*,*) "running on", nimgs, "images in a 3D grid"
  write (*,*) "iarr kind:", iarr, "huge(0_iarr):", huge_iarr

  ! In this test space is assigned image numbers - must be big enough
  ! integer kind to avoid inteter overflow.
  if ( nimgs .gt. huge_iarr ) then
    write (*,*) "ERROR: num_images(): ", nimgs,              &
      " is greater than huge(0_iarr)"
    error stop
  end if
end if

! each image calculates the coarray grid dimensions
call cgca_gdim( nimgs, ir, qual )

! calculate the resolution and the actual phys dimensions of the box
! subroutine cgca_cadim( bsz, res, dm, ir, c, lres, ng )
! c - coarray sizes
! ir - coarray grid sizes
bsz = bsz0
call cgca_cadim( bsz, res, dm, ir, c, lres, ng )

! Check that the partition is sane
if ( img .eq. 1 ) then
  if ( any(int(bsz) .ne. int(bsz0) ) ) then
    write (*,*)                                                        &
      "ERROR: bad decomposition - use a 'nicer' number of images"
    write (*,*) "ERROR: wanted         :", int(bsz0)
    write (*,*) "ERROR: but got instead:", int(bsz)
    error stop
  end if
end if

! total number of cells in a coarray
icells = int( c(1), kind=ilrg ) * int( c(2), kind=ilrg ) *             &
         int( c(3), kind=ilrg )

! total number of cells in the model
mcells = icells * int( nimgs, kind=ilrg )

if ( img .eq. 1 ) then
  write ( *, "(8(a,i0),tr1,g10.3,tr1,g10.3,3(a,g10.3),a)" )            &
    "nimgs: ", nimgs, " (", c(1), "," , c(2), "," , c(3), ")[",        &
    ir(1), "," , ir(2), "," , ir(3), "] ", ng, qual, lres,             &
    " (", bsz(1), ",", bsz(2), ",", bsz(3), ")"
     write (*,'(a,i0,a)') "Each image has ",icells, " cells"
     write (*,'(a,i0,a)') "The model has ", mcells, " cells"
end if

! MPI
! init
! lines
! here
! in
! test_mpi_hxvn

! run=1 => ca_iter_tl
! run=2 => ca_iter_dc
! run=3 => ca_iter_omp
outer: do run=1,3

  if ( img .eq. 1 ) then
    select case( run )
    case(1)
      write (*,*) "Checking ca_iter_tl - triple loop"
    case(2)
      write (*,*) "Checking ca_iter_dc - do concurrent"
    case(3)
      write (*,*) "Checking ca_iter_omp - OpenMP"
    end select
  end if

  ! Loop over several halo depths
  ! The max halo depth is 1/4 of the min dimension
  ! of the space CA array
  main: do d=1, int( 0.25 * min( c(1), c(2), c(3) ) )
  
    ! allocate space array coarrays
    !    space - CA array to allocate, with halos!
    !        c - array with space dimensions
    !        d - depth of the halo layer
    !       ir - codimensions
    call ca_co_spalloc( space,  c, d, ir )
    call ca_co_spalloc( space1, c, d, ir )
  
    ! Set space to my image number
    space  = int( img, kind=iarr )
    space1 = space

    ! No need
    ! for separate
    ! HX coarrays!

    ! MPI subarray types in
    ! test_mpi_hxvn

    ! do hx, remote ops
    call ca_co_hx_all( space )
  
    ! halo check, local ops
    ! space - space array, with halos
    ! flag - default integer
    call ca_co_hx_check( space=space, flag=ierr )
    if ( ierr .ne. 0 ) then
      write (*,*) "ERROR: ca_co_hx_check failed: img:", img,           &
                  "flag:", ierr
      error stop
    end if

    ! CA iterations
    ! subroutine ca_co_run( space, hx_sub, iter_sub, kernel, niter )
    select case( run )
    case(1)
      call ca_co_run( space = space, hx_sub = ca_co_hx_all,            &
        iter_sub = ca_iter_tl,  kernel = ca_kernel_copy, niter = 13 )
    case(2)
      call ca_co_run( space = space, hx_sub = ca_co_hx_all,            &
        iter_sub = ca_iter_dc,  kernel = ca_kernel_copy, niter = 13 )
    case(3)
      call ca_co_run( space = space, hx_sub = ca_co_hx_all,            &
        iter_sub = ca_iter_omp, kernel = ca_kernel_copy, niter = 13 )
    end select
  
    ! Must be the same
    if ( any( space(  1:c(1), 1:c(2), 1:c(3) ) .ne.                    &
              space1( 1:c(1), 1:c(2), 1:c(3) ) ) ) then
      write (*,*) "img:", img, "FAIL: space .ne. space1"
      error stop
    end if

    ! No separate halos
    ! to deallocate

    ! Free MPI types in
    ! test_mpi_hxvn

    ! deallocate space
    deallocate( space )
    deallocate( space1 )
  
    if (img .eq. 1 ) write (*,*) "PASS, halo depth:", d

  end do main

end do outer

end program test_hxvn_co

!*roboend*
