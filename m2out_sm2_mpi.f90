!*robodoc*f* cgca_m2out/m2out_sm2_mpi
!  NAME
!    m2out_sm2_mpi
!  SYNOPSIS

!$Id: m2out_sm2_mpi.f90 380 2017-03-22 11:03:09Z mexas $

submodule ( cgca_m2out ) m2out_sm2_mpi

!  AUTHOR
!    David Henty, modified by Anton Shterenlikht
!  COPYRIGHT
!   Note that this routine has special Copyright conditions.
!
!    !----------------------------------------------------------------------------!
!    !                                                                            !
!    !  MPI-IO routine for Fortran Coarrays                                       !
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
!    Submodule of cgca_m2out dealing with parallel IO using MPI/IO library
!  CONTAINS
!    cgca_pswci
!  USES
!    Variables and parameters of the parent module cgca_m2out
!    via host association.
!  USED BY
!    The parent module cgca_m2out.
!  SOURCE

use mpi
implicit none

contains

!*roboend*


!*robodoc*s* m2out_sm2_mpi/cgca_pswci
!  NAME
!    cgca_pswci
!  SYNOPSIS

!module procedure cgca_pswci
  module subroutine cgca_pswci( coarray, stype, fname )
    ! In submodule m2out_sm2_mpi.f90
    ! Parallel Stream Write Coarray of Integers:
    ! - coarray - what array to dump
    ! - stype - what cell state type to dump
    ! - fname - what file name to use
    integer( kind=iarr ), allocatable, intent( in ) ::                 &
      coarray(:,:,:,:)[:,:,:]
    integer( kind=idef ),intent( in ) :: stype
    character( len=* ), intent( in ) :: fname

!  INPUTS
!    See the interface block in the parent module cgca_m2out.
!  OUTPUTS
!    None
!  AUTHOR
!    Anton Shterenlikht, adapted from the code
!    written by David Henty, EPCC
!  SIDE EFFECTS
!    A single binary file is created using MPI/IO
!    with contents of coarray.
!  DESCRIPTION
!    Parallel Stream Write Coarray of Integers:
!    - coarray - what array to dump
!    - stype - what cell state type to dump
!    - fname - what file name to use
!  NOTES
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
!    cgca_m1co, MPI library
!  USED BY
!    none, end user.
!  SOURCE

integer, parameter :: totdim = 4, arrdim = totdim-1, coardim = 3

integer :: img, nimgs, comm, ierr=0, rank=0, mpisize=0, filetype,      &
           mpi_subarray, fh, funit
integer, dimension(totdim) :: asizehal
integer, dimension(arrdim) :: arrsize, arstart, artsize
integer, dimension(coardim) :: coarsize, copos
integer( kind=MPI_OFFSET_KIND ) :: disp = 0
integer, dimension(MPI_STATUS_SIZE) :: mpistat

character( len=80 ) :: iomsg

  img = this_image()
nimgs = num_images()

asizehal(:) = shape( coarray )
   copos(:) = this_image( coarray )

! Subtract halos
 arrsize(:) = asizehal(1:arrdim) - 2
coarsize(:) = ucobound(coarray) - lcobound(coarray) + 1

! Does the array fit exactly?
if ( product( coarsize ) .ne. nimgs) then
 write(*,*) 'ERROR: m2out_sm2_mpi/cgca_pswci: non-conforming coarray'
 error stop
end if

comm = MPI_COMM_WORLD
call MPI_Comm_size( comm, mpisize, ierr )
call MPI_Comm_rank( comm, rank, ierr )

! Sanity check
if ( mpisize .ne. nimgs .or. rank .ne. img-1 ) then
  write(*,*) 'ERROR: m2out_sm2_mpi/cgca_pswci: MPI/coarray mismatch'
  error stop
end if
     
! Define filetype for this process, ie what portion of the global array
! this process owns. Starting positions use C-indexing
! (ie counting from 0).
artsize(:) = arrsize(:) * coarsize(:)
arstart(:) = arrsize(:) * (copos(:)-1)

! debug
!write (*,*) "image",img, "asizehal", asizehal, "copos", copos,         &
!  "arrsize", arrsize, "coarsize", coarsize,                            &
!  "artsize", artsize, "arstart", arstart, "stype", stype

call MPI_Type_create_subarray( arrdim, artsize, arrsize, arstart,      &
 MPI_ORDER_FORTRAN, MPI_INTEGER, filetype, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
 write (*,*) 'ERROR: m2out_sm2_mpi/cgca_pswci:&
   & MPI_type_create_subarray filetype: rank ', rank
 error stop
end if

call MPI_Type_commit( filetype, ierr )

! Define subarray for this process, ie what portion of the local array
! is to be written (excludes halos); starting positions use C-indexing.

arstart(:) = 1

call MPI_Type_create_subarray( arrdim, asizehal, arrsize, arstart,     &
 MPI_ORDER_FORTRAN, MPI_INTEGER, mpi_subarray, ierr)
if ( ierr .ne. MPI_SUCCESS ) then
 write (*,*) 'ERROR: cgca_pswci: MPI_type_create_subarray&
   & mpi_subarray: rank ', rank
 error stop
end if

call MPI_Type_commit(mpi_subarray, ierr)

!  "Striping information cannot be changed on an existing
! file, so to set the stripe count (and stripe size) for the amount of
! parallelism you want to achieve, the file must be deleted if it exists."
! From: Cray Getting Started on MPI I/O manual, S-2490-40 - Dec 2009:
! http://docs.cray.com/books/S-2490-40/

!if ( rank .eq. 0 ) then
!  call MPI_File_delete( fname, MPI_INFO_NULL, ierr )
!  if ( ierr .ne. MPI_SUCCESS )                                         &
!    error stop "ERROR: cgca_pswci: MPI_file_delete: rank 0"
!end if

! All ranks wait till rank 0 deletes the file
!call MPI_Barrier( comm, ierr )
!if ( ierr .ne. MPI_SUCCESS ) then
!  write (*,*) 'ERROR: cgca_pswci: MPI_file_open: rank ', rank
!  error stop
!end if

! Overwriting MPI/IO files does not involve erasing the file first.
! So if the old file was bigger, the new smaller file will still
! be sized on disk as the old file, with only a part of it overwritten
! with new data. That would be bad. To avoid this
! image 1 removes all previous contents of this file, if it exists.
rm: if ( img .eq. 1 ) then

  ! this should not be necessary, but Cray issues a caution otherwise
  funit = 0

write (*,*) "DEBUG: fname:", fname

  open( newunit=funit, file=fname, status="replace", iostat=ierr,      &
        iomsg=iomsg )
  if ( ierr .ne. 0 ) then
    write (*,*) "ERROR: m2out_sm2_mpi/cgca_pswci: open( fname ),&
      & iostat=", ierr, "iomsg:", iomsg
    error stop
  end if

  write( funit, * , iostat=ierr) ""
  if ( ierr .ne. 0 ) then
    write (*,*) "ERROR: m2out_sm2_mpi/cgca_pswci: write( fname ),&
      & iostat=", ierr
    error stop
  end if

!  flush( funit, iostat=ierr )
!  if ( ierr .ne. 0 ) then
!    write (*,*) "ERROR: cgca_m2mpiio: flush( fname ), iostat=",ierr
!    error stop
!  end if

  close( funit, iostat=ierr )
  if ( ierr .ne. 0 ) then
    write (*,*) "ERROR: m2out_sm2_mpi/cgca_pswci: close( fname ),&
      & iostat=", ierr
    error stop
  end if

end if rm

! all images wait till image 1 erases the previous file
sync all

!  Open the file for writing only and attach to file handle fh
!  No IO hints are passed since MPI_INFO_NULL is specified
call MPI_File_open( comm, fname, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE),&
                    MPI_INFO_NULL, fh, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,*) 'ERROR: m2out_sm2_mpi/cgca_pswci: MPI_file_open: rank ', &
    rank
  error stop
end if

!  Set view for this process using appropriate datatype
call MPI_File_set_view(                                                &
 fh, disp, MPI_INTEGER, filetype, 'native', MPI_INFO_NULL, ierr)
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,*) 'ERROR: m2out_sm2_mpi/cgca_pswci: MPI_file_set_view:&
    & rank ', rank
  error stop
end if

! Write all the data for this process.
! Remove halo data by passing an explicit Fortran subarray
call MPI_File_write_all( fh, coarray(:,:,:,stype), 1, mpi_subarray,    &
                         mpistat, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,*) 'ERROR: m2out_sm2_mpi/cgca_pswci: MPI_file_write_all:&
    & rank ', rank
  error stop
end if

!  Close file
call MPI_File_close( fh, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
 write (*,*) 'ERROR: m2out_sm2_mpi/cgca_pswci: MPI_file_close: rank ', &
   rank
 error stop
end if

call MPI_Type_free( filetype, ierr )
call MPI_Type_free( mpi_subarray, ierr )

!end procedure cgca_pswci
end subroutine cgca_pswci

!*roboend*

end submodule m2out_sm2_mpi
