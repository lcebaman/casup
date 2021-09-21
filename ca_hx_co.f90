!*robodoc*f* ca_hx/ca_hx_co
!  NAME
!    ca_hx_co
!  SYNOPSIS

!$Id: ca_hx_co.f90 533 2018-03-30 14:31:26Z mexas $

submodule ( ca_hx ) ca_hx_co

!  DESCRIPTION
!    Submodule with routines for whole CA implemented with coarrays,
!    not just halos. 
!  AUTHOR 
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!  USES
!  USED BY
!  SOURCE

implicit none

contains

!*roboend*


!*robodoc*s* ca_hx/ca_co_spalloc
!  NAME
!    ca_co_spalloc
!  SYNOPSIS

module subroutine ca_co_spalloc( space, c, d, ir )

!  INPUT

integer( kind=iarr ), allocatable, intent(inout) :: space(:,:,:) [:,:,:]
integer, intent(in) :: c(3), d, ir(3)

!    space - CA array to allocate, with halos!
!        c - array with space dimensions
!        d - depth of the halo layer
!       ir - codimensions
!  OUTPUT
!    space is allocated and set to zero.
!  SIDE EFFECTS
!    none
!  DESCRIPTION
!    This routine allocates the CA coarray array, with halos of depth d.
!    Also save some vars in this module for future.
!    Also, on first call, allocate a module work local space array, not
!    a coarray (tmp_space) of the same mold as space.
!    This array is used in CA iterations later.
!  USES
!  USED BY
!  SOURCE

integer :: i,j,k

if ( allocated( space ) ) then
  write (*,*) "WARN: ca_hx_co/ca_co_spalloc: image:", this_image(),    &
    "space already allocated, deallocating!"
  deallocate( space, stat=ierr, errmsg=errmsg )
  if ( ierr .ne. 0 ) then
    write (*,*)                                                        &
      "ERROR: ca_hx_co/ca_co_spalloc: deallocate( space ), ierr:",     &
                ierr, "errmsg:", trim(errmsg)
    error stop
  end if
end if

allocate( space( 1-d:c(1)+d, 1-d:c(2)+d, 1-d:c(3)+d )[ir(1), ir(2), *],&
          source=0_iarr, stat=ierr, errmsg=errmsg )
if ( ierr .ne. 0 ) then
  write (*,*)                                                          &
    "ERROR: ca_co_hx/ca_co_spalloc: allocate( space ), ierr:",         &
              ierr, "errmsg:", trim(errmsg)
  error stop
end if

! Calculate once and keep forever
! Need a picture to explain these halo vars.
! The picture is only for dimension 1.
! The same diagram can be drawn for dimensions 2 and 3.
!
! hdepth 
! |    | 
! +----+---------------------+----+--
! |    |                     |    | ^ hdepth
! |    |                     |    | v
! +----+---------------------+----+--
! |<-lhsta  |           |    |    |
! |   0|    |           |    |    |
! |    |1   |           |    |    |
! | hdepth->|           |    |    |
! |    |    |           |<-ihsta  |
! |    |    |           sub->|    |
! |    |    |           |    |<-rhsta
! |    |    |           |  rhend->|
! |    |    |           |    |    |
! +----+---------------------+----+
! |    |                     |    |
! |    |                     |    |
! +----+---------------------+----+

hdepth = d
   sub = c
 ihsta = sub - hdepth + 1
 rhsta = sub + 1
 rhend = sub + hdepth
 lhsta = -hdepth + 1

! Total number of cells in the global CA
total_cells = sub(1)*sub(2)*sub(3)*num_images()

! Tmp space array for CA iterations, allocated on the first call
if ( .not. allocated( tmp_space ) ) then
  allocate( tmp_space( 1-d:c(1)+d, 1-d:c(2)+d, 1-d:c(3)+d ),           &
            source=0_iarr, stat=ierr, errmsg=errmsg )
  if ( ierr .ne. 0 ) then
    write (*,"(a,i0,a,a)") "ERROR: ca_hx_co/ca_co_spalloc: " //        &
      "allocate( tmp_space ), ierr: ", ierr, " errmsg: ", trim(errmsg)
    error stop
  end if
end if

! Mask array, no halos, allocated on the first call
! The values are set in ca_halloc.
if ( .not. allocated( mask_array ) ) then
  allocate( mask_array( c(1), c(2), c(3) ), stat=ierr, errmsg=errmsg )
  if ( ierr .ne. 0 ) then
    write (*,"(a,i0,a,a)") "ERROR: ca_hx_co/ca_co_spalloc: " //        &
      "allocate( mask_array ) ", "ierr:", ierr, "errmsg:", trim(errmsg)
    error stop
  end if
end if

! Calculate once and keep forever
  ci = this_image( space )
ucob = ucobound(   space )

main: associate( d => hdepth, c => sub )

! Now can set mask_array. The mask array must reflect
! the global CA space, i.e. must not be affected by partitioning
! of the model into images. This means the mask array must
! depend on coindex set of this image.
do concurrent( i=1:c(1), j=1:c(2), k=1:c(3) )
  mask_array(i,j,k) = int( mod( (i+j+k + (ci(1)-1)*c(1) +              &
    (ci(2)-1)*c(2) + (ci(3)-1)*c(3) ) , 2 ), kind=iarr )
end do

end associate main

! Calculate coindex sets and image numbers for the 6 neighbours
! Coindex sets
  nei_ci_L1 = (/ ci(1)-1, ci(2), ci(3) /)
  nei_ci_R1 = (/ ci(1)+1, ci(2), ci(3) /)
  nei_ci_L2 = (/ ci(1), ci(2)-1, ci(3) /)
  nei_ci_R2 = (/ ci(1), ci(2)+1, ci(3) /)
  nei_ci_L3 = (/ ci(1), ci(2), ci(3)-1 /)
  nei_ci_R3 = (/ ci(1), ci(2), ci(3)+1 /)

! Image index
  nei_img_L(1) = image_index( space, nei_ci_L1 )
  nei_img_R(1) = image_index( space, nei_ci_R1 )
  nei_img_L(2) = image_index( space, nei_ci_L2 )
  nei_img_R(2) = image_index( space, nei_ci_R2 )
  nei_img_L(3) = image_index( space, nei_ci_L3 )
  nei_img_R(3) = image_index( space, nei_ci_R3 )

end subroutine ca_co_spalloc

!*roboend*


!*robodoc*s* ca_hx/ca_co_hx_all
!  NAME
!    ca_co_hx_all
!  SYNOPSIS

module subroutine ca_co_hx_all( space )

!  INPUT

integer( kind=iarr ), intent(inout), allocatable :: space(:,:,:) [:,:,:]

!    space - coarray array with CA model
!  OUTPUT
!    space is updated, just the halo layers
!  SIDE EFFECTS
!    none
!  DESCRIPTION
!    There is pair-wise handshake sync between images.
!    Each image, except those on the boundary, will send halos to and
!    receive halos from 6 neighbouring images.
!  USES
!  USED BY
!  SOURCE

! An image updates its space coarray halo layer (left side)
! along direction 1 from the last layer in the coarray on an image
! which is 1 lower along codimension 1.
if ( ci(1) .ne. 1 ) then
  sync images( nei_img_L(1) )
  space( lhsta(1):0,           1:sub(2) , 1:sub(3) ) =                 &
    space( ihsta(1) : sub(1) , 1:sub(2) , 1:sub(3) )                   &
      [ nei_ci_L1(1), nei_ci_L1(2), nei_ci_L1(3) ]
end if

! An image updates its space coarray halo layer (right side)
! along direction 1 from the 1st layer in the coarray on an image
! which is 1 higher along codimension 1.
if ( ci(1) .ne. ucob(1) ) then
  sync images( nei_img_R(1) )
  space( rhsta(1) : rhend(1) , 1:sub(2) , 1:sub(3) ) =                 &
    space( 1 : hdepth        , 1:sub(2) , 1:sub(3) )                   &
      [ nei_ci_R1(1), nei_ci_R1(2), nei_ci_R1(3) ]
end if

! An image updates its space coarray halo layer (left side)
! along direction 2 from the last layer in the coarray on an image
! which is 1 lower along codimension 2.
if ( ci(2) .ne. 1 ) then
  sync images( nei_img_L(2) )
  space(   1:sub(1) , lhsta(2) : 0,       1:sub(3) ) =                 &
    space( 1:sub(1) , ihsta(2) : sub(2) , 1:sub(3) )                   &
      [ nei_ci_L2(1), nei_ci_L2(2), nei_ci_L2(3) ]
end if

! An image updates its space coarray halo layer (right side)
! along direction 2 from the 1st layer in the coarray on an image
! which is 1 higher along codimension 1.
if ( ci(2) .ne. ucob(2) ) then
  sync images( nei_img_R(2) )
  space(   1:sub(1) , rhsta(2) : rhend(2) , 1:sub(3) ) =               &
    space( 1:sub(1) , 1 : hdepth          , 1:sub(3) )                 &
      [ nei_ci_R2(1), nei_ci_R2(2), nei_ci_R2(3) ]
end if

! An image updates its space coarray halo layer (left side)
! along direction 3 from the last layer in the coarray on an image
! which is 1 lower along codimension 3.
if ( ci(3) .ne. 1 ) then
  sync images( nei_img_L(3) )
  space(   1:sub(1) , 1:sub(2) , lhsta(3) : 0 ) =                      &
    space( 1:sub(1) , 1:sub(2) , ihsta(3) : sub(3) )                   &
    [ nei_ci_L3(1), nei_ci_L3(2), nei_ci_L3(3) ]
end if

! An image updates its space coarray halo layer (right side)
! along direction 3 from the 1st layer in the coarray on an image
! which is 1 higher along codimension 3.
if ( ci(3) .ne. ucob(3) ) then
  sync images( nei_img_R(3) )
  space(   1:sub(1) , 1:sub(2) , rhsta(3) : rhend(3) ) =               &
    space( 1:sub(1) , 1:sub(2) , 1 : hdepth )                          & 
      [ nei_ci_R3(1), nei_ci_R3(2), nei_ci_R3(3) ]
end if

end subroutine ca_co_hx_all

!*roboend*


!*robodoc*s* ca_hx/ca_co_hx_check
!  NAME
!    ca_co_hx_check
!  SYNOPSIS

module subroutine ca_co_hx_check( space, flag )

!  INPUT

integer( kind=iarr ), intent( in ), allocatable :: space(:,:,:) [:,:,:]

!    space - CA array coarray

!  OUTPUT

integer, intent( out ) :: flag

!    flag .eq. 0 - check passed
!    flag .ne. 0 - check failed
!  SIDE EFFECTS
!    none
!  DESCRIPTION
!    This routine is of very limited use. It checks hx code for only
!    a single case - all images set all their cell values to
!    this_image() and then a single hx step is done.
!    In this case I'm sure what's in my halos. So just check for this
!    and either return flag=0 or flag > 0 indicating which halo failed.
!    Fail values:
!      0  - pass, no failures
!      1  - test 1 failed
!      2  - test 2 failed
!      3  - tests 1,2 failed
!      4  - test 3 failed
!      5  - tests 1,3 failed
!      6  - tests 2,3 failed
!      7  - tests 1,2,3 failed
!      8  - test 4 failed
!      9  - tests 1,4 failed
!      10 - tests 2,4 failed
!      11 - tests 1,2,4 failed
!      12 - tests 3,4 failed
!      13 - tests 1,3,4 failed
!      14 - tests 2,3,4 failed
!      15 - tests 1,2,3,4 failed
!      16 - test 5 failed
!    and so on.
!  USES
!  USED BY
!  SOURCE

! coindex set and the image number of the neighbour
integer :: i(3), n

flag = 0

! Test 1
if ( ci(1) .ne. 1 ) then
  ! This is the neighbour
  i = (/ ci(1)-1, ci(2), ci(3) /) ! neighbour's coindex set
  n = image_index( space, i )    ! neighbour image number
  if ( any( space( lhsta(1):0,        1:sub(2), 1:sub(3) ) .ne. n ) )  &
    flag = flag + 1
end if

! Test 2
if ( ci(1) .ne. ucob(1) ) then
  ! This is the neighbour
  i = (/ ci(1)+1, ci(2), ci(3) /) ! neighbour's coindex set
  n = image_index( space, i )    ! neighbour image number
  if ( any( space( rhsta(1):rhend(1), 1:sub(2), 1:sub(3) ) .ne. n ) )  &
    flag = flag + 2
end if

! Test 3
if ( ci(2) .ne. 1 ) then
  ! This is the neighbour
  i = (/ ci(1), ci(2)-1, ci(3) /) ! neighbour's coindex set
  n = image_index( space, i )    ! neighbour image number
  if ( any( space( 1:sub(1), lhsta(2):0,        1:sub(3) ) .ne. n ) )  &
    flag = flag + 4
end if

! Test 4
if ( ci(2) .ne. ucob(2) ) then
  ! This is the neighbour
  i = (/ ci(1), ci(2)+1, ci(3) /) ! neighbour's coindex set
  n = image_index( space, i )    ! neighbour image number
  if ( any( space( 1:sub(1), rhsta(2):rhend(2), 1:sub(3) ) .ne. n ) )  &
    flag = flag + 8
end if

! Test 5
if ( ci(3) .ne. 1 ) then
  ! This is the neighbour
  i = (/ ci(1), ci(2), ci(3)-1 /) ! neighbour's coindex set
  n = image_index( space, i )    ! neighbour image number
  if ( any( space( 1:sub(1), 1:sub(2), lhsta(3):0 )        .ne. n ) )  &
    flag = flag + 16
end if

! Test 6
if ( ci(3) .ne. ucob(3) ) then
  ! This is the neighbour
  i = (/ ci(1), ci(2), ci(3)+1 /) ! neighbour's coindex set
  n = image_index( space, i )     ! neighbour image number
  if ( any( space( 1:sub(1), 1:sub(2), rhsta(3):rhend(3) ) .ne. n ) )  &
    flag = flag + 32
end if

end subroutine ca_co_hx_check

!*roboend*


!*robodoc*s* ca_hx/ca_co_run
!  NAME
!    ca_co_run
!  SYNOPSIS
  
module subroutine ca_co_run( space, hx_sub, iter_sub, kernel, niter )
    
!  INPUT
    
integer( kind=iarr ), intent(inout), allocatable :: space(:,:,:) [:,:,:]
procedure( hx_co_proto )  :: hx_sub
procedure( iter_proto )   :: iter_sub
procedure( kernel_proto ) :: kernel
integer, intent(in) :: niter
  
!       space - space coarray before iterations start
!      hx_sub - HX routine, e.g.
!             - ca_co_hx_all
!    iter_sub - the subroutine performing a single CA iteration, e.g.
!             - ca_iter_tl - triple nested loop
!             - ca_iter_dc - do concurrent loop
!             - ca_iter_omp - OpenMP loop
!      kernel - a function to be called for every cell inside the loop
!        iter - number of iterations to do
!  OUTPUT
!    space - CA coarray at the end of niter iterations
!  SIDE EFFECTS
!    module array tmp_space is updated
!  DESCRIPTION
!    This is a driver routine for CA iterations. HX is done
!    before each iteration. Then a given number of iterations
!    is performed with a given routine and a given kernel.
!    One iteration is really 2 iterations: odd and even,
!    Hence the upper loop limit is 2*niter.
!  USES
!  USED BY
!  SOURCE

integer :: i

tmp_space = space

do i = 1, 2*niter
  call hx_sub( space )                   ! space updated, with HX
  call iter_sub( space, hdepth, kernel ) ! tmp_space updated, local op
  space = tmp_space                      ! local op
  mask_array = 1_iarr - mask_array       ! Flip the mask array
end do

end subroutine ca_co_run

!*roboend*


!*robodoc*s* ca_hx/ca_co_ising_energy
!  NAME
!    ca_co_ising_energy
!  SYNOPSIS

module subroutine ca_co_ising_energy( space, hx_sub, iter_sub, kernel, &
  energy, magnet )

!  INPUT

integer( kind=iarr ), intent(inout), allocatable :: space(:,:,:)[:,:,:]
procedure( hx_co_proto )  :: hx_sub
procedure( iter_proto )   :: iter_sub
procedure( kernel_proto ) :: kernel

!       space - space coarray before iterations start
!      hx_sub - HX routine
!             - ca_co_hx_all
!    iter_sub - the subroutine performing a single CA iteration, e.g.
!             - ca_co_iter_tl - triple nested loop
!             - ca_co_iter_dc - do concurrent
!             - ca_co_iter_omp - OpenMP
!      kernel - a function to be called for every cell inside the loop
!  OUTPUT

integer( kind=ilrg) , intent(out) :: energy, magnet

!    energy - Total energy of CA system
!    magnet - Total magnetisation of the CA system
!  SIDE EFFECTS
!    module array tmp_space is updated
!  DESCRIPTION
!    Calculate the total energy and the total magnetisation
!    of CA using Ising model. These integer values might be very large
!    so I'm using a large integer kind (ilrg).
!    This routine uses collectives.
!    Magnetisation is defined as the fraction of the 1 spins.
!    The only valid kernel is ca_kernel_ising_ener.
!  USES
!  USED BY
!  SOURCE

call hx_sub( space )                   ! space updated, sync images
call iter_sub( space, hdepth, kernel ) ! tmp_space updated, local op
energy = sum( tmp_space( 1:sub(1), 1:sub(2), 1:sub(3) ) ) 
magnet = sum(     space( 1:sub(1), 1:sub(2), 1:sub(3) ) ) 
call co_sum( energy )
call co_sum( magnet )

end subroutine ca_co_ising_energy

!*roboend*


!*robodoc*s* ca_hx/ca_co_netcdf
!  NAME
!    ca_co_netcdf
!  SYNOPSIS

module subroutine ca_co_netcdf( space, fname )
use netcdf

!  INPUTS

     integer( kind=iarr ), intent( in ), allocatable ::                &
       space(:,:,:) [:,:,:]
     character( len=* ), intent( in ) :: fname

!    coarray - what coarray to dump
!    fname - what file name to use
!  OUTPUTS
!    None
!  SIDE EFFECTS
!    A single binary file is created using netcdf
!    with contents of coarray.
!  AUTHOR
!    Anton Shterenlikht, adapted from the code
!    written by David Henty, Luis Cebamanos, EPCC
!  COPYRIGHT
!   Note that this routine has special Copyright conditions.
!
!    !----------------------------------------------------------------------------!
!    !                                                                            !
!    !  netCDF routine for Fortran Coarrays                                       !
!    !                                                                            !
!    !  David Henty, EPCC; d.henty@epcc.ed.ac.uk                                  !
!    !                                                                            !
!    !  Copyright 2013 the University of Edinburgh                                !
!    !                                                                            !
!    !  Licensed under the Apache License, Version 2.0 (the "License");           !
!    !  you may not use this file except in compliance with the License.          !
!    !  You may obtain a copy of the License at                                   !
!    !                                                                            !
!    !      http://www.apache.org/licenses/LICENSE-2.0                            !
!    !                                                                            !
!    !  Unless required by applicable law or agreed to in writing, software       !
!    !  distributed under the License is distributed on an "AS IS" BASIS,         !
!    !  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  !
!    !  See the License for the specific language governing permissions and       !
!    !  limitations under the License.                                            !
!    !                                                                            !
!    !----------------------------------------------------------------------------!
!
!  DESCRIPTION
!    Parallel Stream Write Coarray of Integers:
!    All images must call this routine!
!
!    MPI must be initialised prior to calling this routine,
!    most probably in the main program.
!    Likewise MPI must be terminated only when no further MPI
!    routines can be called. This will most likely be in the
!    main program. 
!    There are some assumptions about the shape of the passed array.
!
!    The default integer is assumed for the array at present! 
!  USES
!    MPI library, netCDF library
!  USED BY
!    none, end user.
!  SOURCE

integer, parameter :: totdim = 3, arrdim = totdim, coardim = 3

integer :: img, nimgs, comm, ierr=0, rank=0, mpisize=0
          
integer, dimension(totdim) :: asizehal
integer, dimension(arrdim) :: arrsize, arstart, artsize
integer, dimension(coardim) :: coarsize, copos

integer :: ncid, varid, dimids(arrdim)
integer :: x_dimid=0, y_dimid=0, z_dimid=0

  img = this_image()
nimgs = num_images()

asizehal(:) = shape( space )
   copos(:) = this_image( space )

! Subtract halos
 arrsize(:) = asizehal(1:arrdim) - 2*hdepth
coarsize(:) = ucobound(space) - lcobound(space) + 1

! Does the array fit exactly?
if ( product( coarsize ) .ne. nimgs) then
 write(*,*) 'ERROR: ca_hx/ca_co_netcdf: non-conforming coarray'
 error stop
end if

comm = MPI_COMM_WORLD
call MPI_Comm_size( comm, mpisize, ierr )
call MPI_Comm_rank( comm, rank, ierr )

! Sanity check
if ( mpisize .ne. nimgs .or. rank .ne. img-1 ) then
  write(*,*) 'ERROR: ca_hx/ca_co_netcdf: MPI/coarray mismatch'
  error stop
end if
     

! This is the global arrayy
artsize(:) = arrsize(:) * coarsize(:)

! Correspondent portion of this global array 
arstart(:) = arrsize(:) * (copos(:)-1) + 1 ! Use Fortran indexing

! ! debug
! write (*,*) "netCDF-image",img, "asizehal", asizehal, "copos", copos,         &
!  "arrsize", arrsize, "coarsize", coarsize,                            &
!  "artsize", artsize, "arstart", arstart

! Create (i.e. open) the netCDF file. The NF90_NETCDF4 flag causes a 
! HDF5/netCDF-4 type file to be created. The comm and info parameters 
! cause parallel I/O to be enabled. 
ierr = nf90_create( fname, ior(nf90_netcdf4,nf90_mpiio), ncid, &
                    comm=comm, info=MPI_INFO_NULL )
if ( ierr .ne. nf90_noerr) then 
  write (*,*) "ERROR: ca_hx/ca_co_netcdf: nf90_create:",               &
              nf90_strerror( ierr )
  error stop
end if

! Define the dimensions. NetCDF returns an ID for each. Any 
! metadata operations must take place on ALL processors
ierr = nf90_def_dim( ncid, "x", artsize(1), x_dimid )
if ( ierr .ne. nf90_noerr) then 
  write (*,*) "ERROR: ca_hx/ca_co_netcdf: nf90_def_dim X:",            &
              nf90_strerror( ierr )
  error stop
end if

ierr = nf90_def_dim( ncid, "y", artsize(2), y_dimid )
if ( ierr .ne. nf90_noerr) then 
  write (*,*) "ERROR: ca_hx/ca_co_netcdf: nf90_def_dim Y:",            &
              nf90_strerror( ierr )
  error stop
end if

ierr = nf90_def_dim( ncid, "z", artsize(3), z_dimid )
if ( ierr .ne. nf90_noerr) then 
  write (*,*) "ERROR: ca_hx/ca_co_netcdf: nf90_def_dim Z:",            &
              nf90_strerror( ierr )
  error stop
end if

! The dimids array is used to pass the ID's of the dimensions of 
! the variables. 
dimids = (/ x_dimid , y_dimid, z_dimid /)

! Define the variable. The type of the variable in this case is
! NF90_INT (4-byte int).
ierr = nf90_def_var(ncid, "data", NF90_INT, dimids, varid) 
if ( ierr .ne. nf90_noerr) then 
  write (*,*) "ERROR: ca_hx/ca_co_netcdf: nf90_def_var:",              &
              nf90_strerror( ierr )
  error stop
end if

! Make sure file it not filled with default values
! which doubles wrote volume
ierr = nf90_def_var_fill(ncid, varid, 1, 1) 
if ( ierr .ne. nf90_noerr) then 
  write (*,*) "ERROR: ca_hx/ca_co_netcdf: nf90_def_var_fill:",         &
              nf90_strerror( ierr )
  error stop
end if

! End define mode. This tells netCDF we are done defining
! metadata. This operation is collective and all processors will
! write their metadata to disk.
ierr = nf90_enddef(ncid)
if ( ierr .ne. nf90_noerr ) then 
  write (*,*) "ERROR: ca_hx/ca_co_netcdf: nf90_enddef:",               &
              nf90_strerror( ierr )
  error stop
end if

! Parallel access
ierr = nf90_var_par_access( ncid, varid, nf90_collective )
if ( ierr .ne. nf90_noerr ) then 
  write (*,*) "ERROR: ca_hx/ca_co_netcdf: nf90_var_par_access:",       &
              nf90_strerror( ierr )
  error stop
end if

! Write the data to file, start will equal the displacement from the 
! start of the file and count is the number of points each proc writes. 
ierr = nf90_put_var( ncid, varid,                                      &
             space( 1:arrsize(1), 1:arrsize(2), 1:arrsize(3) ),        &
             start = arstart, count = arrsize ) 
if ( ierr .ne. nf90_noerr) then 
  write (*,*) "ERROR: ca_hx/ca_co_netcdf: nf90_put_var:",              &
              nf90_strerror( ierr )
  error stop
end if

! Close the file. This frees up any internal netCDF resources
! associated with the file, and flushes any buffers.
ierr = nf90_close(ncid)
if ( ierr .ne. nf90_noerr) then 
  write (*,*) "ERROR: ca_hx/ca_co_netcdf: nf90_close:",                &
              nf90_strerror( ierr )
  error stop
end if

end subroutine ca_co_netcdf

!*roboend*


!*robodoc*s* ca_hx/ca_co_naive_io
!  NAME
!    ca_co_naive_io
!  SYNOPSIS

module subroutine ca_co_naive_io( coarray, fname )

!  INPUTS
 
integer( kind=iarr ), allocatable, intent( in ) :: coarray(:,:,:)[:,:,:]
character( len=* ),intent( in ) :: fname

!    - coarray - what array to dump
!    - fname - what file name to use
!  SIDE EFFECTS
!    A single binary file is created on image 1 with contents of coarray.
!  DESCRIPTION
!    All images call this routine!
!    However only image 1 does all the work.
!    The other images are waiting. 
!  USES
!    none
!  USED BY
!    none, end user.
!  SOURCE

integer :: ierr=0, coi1, coi2, coi3, i2, i3, funit,                    &
  lb(3),   & ! lower bounds   of the coarray
  ub(3),   & ! upper bounds   of the coarray
  lcob(3), & ! lower cobounds of the coarray
  ucob(3)    ! upper cobounds of the coarray

! Only image1 does this. All other images do nothing.
! So sync all probably should be used after a call to
! this routine in the program.

main: if ( this_image() .eq. 1 ) then

  ! Assume the coarray has halos. Don't write those.
    lb = lbound( coarray ) + hdepth
    ub = ubound( coarray ) - hdepth
  lcob = lcobound( coarray )
  ucob = ucobound( coarray )

  open( newunit=funit, file=fname, form="unformatted",                &
        access="stream", status="replace", iostat=ierr )
  if ( ierr .ne. 0 ) then
    write (*,*) "ERROR: ca_hx/ca_co_naive_io: open:", ierr
    error stop
  end if

 ! nested loops for writing in correct order from all images
 do coi3 = lcob(3), ucob(3)
  do i3 = lb(3), ub(3)
   do coi2 = lcob(2), ucob(2)
    do i2 = lb(2), ub(2)
     do coi1 = lcob(1), ucob(1)
       write( unit=funit, iostat=ierr )                                &
            coarray( lb(1):ub(1), i2, i3 ) [ coi1, coi2, coi3 ]
     end do
    end do
   end do
  end do
 end do

 close( unit=funit, iostat=ierr )
 if ( ierr .ne. 0 ) then
   write (*,*) "ERROR: ca_hx/ca_co_naive_io: close:", ierr
  error stop
 end if

end if main

end subroutine ca_co_naive_io
!*roboend*

end submodule ca_hx_co
