!*robodoc*u* tests/test_ising_perf
!  NAME
!    test_ising_perf
!  SYNOPSIS

!$Id: test_ising_perf.f90 536 2018-04-03 12:02:13Z mexas $

program test_ising_perf

!  PURPOSE
!    Test performance of Ising magnetisation.
!    This test uses coarray collectives.
!    Reproducibility on any number of images requires a serial RND
!    for the whole model. This is too slow. Here each image uses
!    its own RND, meaning the results are *not* reproducible when
!    the number of images change. This test is purely for performance
!    analysis.
!  DESCRIPTION
!    See ca_kernel_ising and related routines for details.
!    Note that I use a reproducible RND seed and generate
!    a single sequence of RND values for the whole CA model.
!    Thus the results must be exactly reproducible on any number
!    of images. I include the reference value for the
!    final magnetisation (unscaled, integer). If the test
!    does not produce the same value, it fails.
!    However... the ref magnetisation value is obtained here with
!    gfortran7. It is possible (likely?) that other compliers
!    will produce a different sequence of RND from the same seed.
!    In such cases users need to replace the ref value accordingly.
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  USES
!    cgca
!  USED BY
!    Part of CGPACK test suite
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

integer( kind=iarr ), allocatable :: space(:,:,:), space0(:,:,:)

integer( kind=idef ) :: ir(3), nimgs, img, ng, c(3) ! coarray dimensions

integer( kind=ilrg ) :: icells, mcells

!real, allocatable :: space_ini(:,:,:), rnd_array(:)

integer :: i, iter, seed_size, run
integer( kind=ilrg ) :: energy0, energy1, energy2, magnet0, magnet1,   &
  magnet2
integer, allocatable :: seed_array(:)

real :: time1, time2
real, allocatable :: rnd_arr(:,:,:)

!*********************************************************************72
! first executable statement

 dm = 1.0 ! Linear "size" of one spin cell
res = 1.0 ! resolution, CA cells per spin

! When dm=res=1, then bsz0 is simply CA dimensions in cells! 
bsz0 = (/ 3.0e2, 3.0e2, 3.0e2 /) ! dimensions of the CA model

  img = this_image()
nimgs = num_images()

! do a check on image 1
if ( img .eq. 1 ) then
  write (*,*) "running on", nimgs, "images in a 3D grid"
  write (*,*) "iarr kind:", iarr, "huge(0_iarr):", huge_iarr

  ! In this test sum over all cells in the *model* is done, so the kind
  ! must be big enough to contain the total number of cells in the
  ! model to avoid overflow.
  mcells = nint( product( bsz0 ), kind=ilrg )
  if ( mcells .gt. huge_iarr ) then
    write (*,*) "ERROR: total number of cells in the model:", mcells, &
      "is greater than huge(0_iarr)"
    error stop
  end if
end if

sync all ! Might help to separate the data sent to output_unit

! each image calculates the coarray grid dimensions
call cgca_gdim( nimgs, ir, qual )

! calculate the resolution and the actual phys dimensions of the box
! subroutine cgca_cadim( bsz, res, dm, ir, c, lres, ng )
! c - coarray sizes
! ir - coarray grid sizes
bsz = bsz0
call cgca_cadim( bsz, res, dm, ir, c, lres, ng )

! Check that the partition is sane
!if ( img .eq. 1 ) then
!  if ( any(int(bsz) .ne. int(bsz0) ) ) then
!    write (*,*)                                                        &
!      "ERROR: bad decomposition - use a 'nicer' number of images"
!    write (*,*) "ERROR: wanted         :", int(bsz0)
!    write (*,*) "ERROR: but got instead:", int(bsz)
!    error stop
!  end if
!end if 

! total number of cells in a coarray
icells = int( c(1), kind=ilrg ) * int( c(2), kind=ilrg ) *             &
         int( c(3), kind=ilrg )

! total number of cells in the model
mcells = icells * int( nimgs, kind=ilrg )

if ( img .eq. 1 ) then
  write ( *, "(8(a,i0),2(tr1,g10.3),3(a,g10.3),a)" )                   &
    "nimgs: ", nimgs, " (", c(1), "," , c(2), "," , c(3), ")[",        &
    ir(1), "," , ir(2), "," , ir(3), "] ", ng, qual, lres,             &
    " (", bsz(1), ",", bsz(2), ",", bsz(3), ")"
     write (*,'(a,i0)') "Cells per image: ", icells
     write (*,'(a,i0)') "Cells per model: ", mcells
end if

! allocate space arrays and set all values to zero
!    space - CA array to allocate, with halos!
!        c - array with space dimensions
!        d - depth of the halo layer
call ca_spalloc( space,  c, 1 )
call ca_spalloc( space0, c, 1 )

! allocate hx arrays, implicit sync all
! mask_array is set inside too.
! ir(3) - codimensions
call ca_halloc( ir )

! Init RND
!call cgca_irs( debug = .false. )
! Use a reproducible RND here for verification
call random_seed( size = seed_size )
allocate( seed_array( seed_size ) )
seed_array = (/ (i, i=1,seed_size) /)
call random_seed( put = seed_array )

! Set space arrays
allocate( rnd_arr( lbound(space, dim=1) : ubound(space, dim=1),        &
                   lbound(space, dim=2) : ubound(space, dim=2),        &
                   lbound(space, dim=3) : ubound(space, dim=3) ) )
call random_number( rnd_arr )
space = nint( rnd_arr, kind=iarr )

! Calculate initial energy and magnetisation

! run=1 => ca_iter_tl
! run=2 => ca_iter_dc
! run=3 => ca_iter_omp
do run=1,3
  select case(run)
  case(1)
    call ca_ising_energy_col( space=space, iter_sub = ca_iter_tl,       &
      kernel = ca_kernel_ising_ener, energy=energy0, magnet=magnet0 )
  case(2)
    call ca_ising_energy_col( space=space, iter_sub = ca_iter_dc,       &
      kernel = ca_kernel_ising_ener, energy=energy1, magnet=magnet1 )
  case(3)
    call ca_ising_energy_col( space=space, iter_sub = ca_iter_omp,      &
      kernel = ca_kernel_ising_ener, energy=energy2, magnet=magnet2 )
  end select
end do

if (img .eq. 1 ) then
  write (*,*) "Initial energy and magnetisation"
  write (*,*) "ca_iter_tl :", energy0, magnet0
  write (*,*) "ca_iter_dc :", energy1, magnet1
  write (*,*) "ca_iter_omp:", energy2, magnet2
  if ( energy0 .ne. energy1 .or. magnet0 .ne. magnet1 .or.             &
       energy0 .ne. energy2 .or. magnet0 .ne. magnet2 ) then
    write (*,*) "FAIL: ca_iter_tl, ca_iter_dc, ca_iter_omp differ"
    error stop
  else
    write (*,*) "PASS: ca_iter_tl, ca_iter_dc, ca_iter_omp agree"
  end if
end if

! save old space as space0
space0 = space

! run=1 => ca_iter_tl
! run=2 => ca_iter_dc
! run=3 => ca_iter_omp
main: do run=1,3

  ! Reset space to space0
  space = space0

  ! Start counter
  if ( img .eq. 1 ) call cpu_time( time1 )

  ! CA iterations
  loop: do iter = 1,100
 
    ! Check energy after every iter
    ! subroutine ca_run( space, hx_sub, iter_sub, kernel, niter )
    select case(run)
    case(1)
      call ca_run( space = space, hx_sub = ca_hx_all,                  &
        iter_sub = ca_iter_tl,  kernel = ca_kernel_ising, niter = 1 )
      call ca_ising_energy_col( space=space, iter_sub = ca_iter_tl,     &
        kernel = ca_kernel_ising_ener, energy=energy1, magnet=magnet1 )
    case(2)
      call ca_run( space = space, hx_sub = ca_hx_all,                  &
        iter_sub = ca_iter_dc,  kernel = ca_kernel_ising, niter = 1 )
      call ca_ising_energy_col( space=space, iter_sub = ca_iter_dc,     &
        kernel = ca_kernel_ising_ener, energy=energy1, magnet=magnet1 )
    case(3)
      call ca_run( space = space, hx_sub = ca_hx_all,                  &
        iter_sub = ca_iter_omp, kernel = ca_kernel_ising, niter = 1 )
      call ca_ising_energy_col( space=space, iter_sub = ca_iter_omp,    &
        kernel = ca_kernel_ising_ener, energy=energy1, magnet=magnet1 )
    end select 
  
    if ( img .eq. 1 ) then
      if ( energy1 .ne. energy0 ) then
        write (*,*) "FAIL: energy0:", energy0, "energy1:", energy1
        error stop
      else
  !      write (*,"(a,i0,a,es18.6)") "Magnetisation_after_iter_",      &
  !        iter, ":", real(magnet1) / real(mcells)
      end if
    end if
  
  end do loop

  if ( img .eq. 1 ) then
    ! Stop counter
    call cpu_time( time2 )

    select case(run)
      case(1)
          write (*,*) "TIME ca_iter_tl, s:", time2-time1
      case(2)
          write (*,*) "TIME ca_iter_dc, s:", time2-time1
      case(3)
          write (*,*) "TIME ca_iter_omp,s:", time2-time1
    end select

  end if

end do main

! deallocate halos, implicit sync all
call ca_hdalloc

! deallocate space
deallocate( space )
deallocate( space0 )

end program test_ising_perf

!*roboend*
