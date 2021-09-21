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

module cafmpiio
  use mpi
  implicit none

  integer, parameter, private :: ndim = 3

private
public :: cafwrite

! Only meed these six parameters because of Cray compiler bug (see below)

!  integer, parameter, private :: n1 = 128
!  integer, parameter, private :: n2 = 128
!  integer, parameter, private :: n3 =  64
!  integer, parameter, private :: p1 =  8
!  integer, parameter, private :: p2 =  8
!  integer, parameter, private :: p3 = 16

contains

subroutine cafwrite(filename, array)

  character*(*),intent(in) :: filename
  integer, allocatable, intent(in) :: array(:,:,:)[:,:,:]

!  integer, dimension(0:n1+1,0:n2+1,0:n3+1), codimension[p1,p2,*] :: array
! Would like to use the definition below so the subroutine is generic for
! all array sizes, but currently this exposes a Cray compiler bug.
!
!  integer, dimension(:,:,:), codimension[:,:,*] :: array


  integer :: myimage, numimage, comm, ierr, rank, mpisize, &
   filetype, mpi_subarray, fh
  integer, dimension(ndim) :: arraysize, arraystart, arraycosize, copos, &
   arraygsize, arraysubsize
  integer (kind=MPI_OFFSET_KIND) :: disp = 0
  integer, dimension(MPI_STATUS_SIZE) :: status

  myimage  = this_image()
  numimage = num_images()

  arraysize(:)      = shape(array)

! Subtract halos for array subsize

  arraysubsize(:)   = arraysize (:) - 2

! Coshape does not seem to work ...
!  arraycosize(:) = coshape(array)

  arraycosize(:) = ucobound(array) - lcobound(array) + 1

! Does the array fit exactly?

  if (arraycosize(1)*arraycosize(2)*arraycosize(3) /= numimage) then
     write(*,*) 'ERROR: non-conforming coarray!'
     error stop
  end if

  copos(:) = this_image(array)

!  write(*,*) 'coloc, cosize, xsize, ysize = ', &
!              copos(1), copos(2), copos(3), &
!              arraycosize(1), arraycosize(2), arraycosize(3), &
!              arraysubsize(1), arraysubsize(2), arraysubsize(3)

  call MPI_Init(ierr)

  comm = MPI_COMM_WORLD

  call MPI_Comm_size(comm, mpisize, ierr)
  call MPI_Comm_rank(comm, rank, ierr)

! Sanity check

  if (mpisize /= numimage .or. rank /= myimage-1) then
     write(*,*) 'ERROR: MPI / coarray mismatch!'
     stop
  end if
     
!
! Define filetype for this process, ie what portion of the global array
! this process owns; starting positions use C-indexing (ie counting from 0).
!

  arraygsize(:) = arraysubsize(:) * arraycosize(:)
  arraystart(:) = arraysubsize(:) * (copos(:)-1)

  call MPI_Type_create_subarray(ndim, arraygsize, arraysubsize, arraystart, &
                                MPI_ORDER_FORTRAN, MPI_INTEGER, &
                                filetype, ierr)

  call MPI_Type_commit(filetype, ierr)

!
! Define subarray for this process, ie what portion of the local array
! is to be written (excludes halos); starting positions use C-indexing.
!

  arraystart(:) = 1

  call MPI_Type_create_subarray(ndim, arraysize, arraysubsize, arraystart, &
                                MPI_ORDER_FORTRAN, MPI_INTEGER, &
                                mpi_subarray, ierr)

  call MPI_Type_commit(mpi_subarray, ierr)

!
!  Open the file for reading only and attach to file handle fh
!  No IO hints are passed since MPI_INFO_NULL is specified
!

  call MPI_File_open(comm, filename, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
                     MPI_INFO_NULL, fh, ierr)

  if (ierr /= MPI_SUCCESS) write(*,*) 'Open error on rank ', rank

!
!  Set view for this process using appropriate datatype
!

  call MPI_File_set_view(fh, disp, MPI_INTEGER, filetype, 'native', &
                         MPI_INFO_NULL, ierr)

  if (ierr /= MPI_SUCCESS) write(*,*) 'View error on rank ', rank

!
!  Write all the data for this process.
!  Remove halo data by passing an explicit Fortran subarray
!
!
!  call MPI_File_write_all(fh, &
!       array(1:arraysubsize(1), 1:arraysubsize(2), 1:arraysubsize(3)), &
!       arraysubsize(1)*arraysubsize(2)*arraysubsize(3), &
!       MPI_INTEGER, status, ierr)

!
!  Write all the data for this process.
!  Remove halo data by passing an MPI subarray type
!

  call MPI_File_write_all(fh, array, 1, mpi_subarray, status, ierr)

  if (ierr /= MPI_SUCCESS) write(*,*) 'Write error on rank ', rank

!
!  Close file
!

  call MPI_File_close(fh, ierr)

  if (ierr /= MPI_SUCCESS) write(*,*) 'Close error on rank ', rank

  call MPI_Type_free(filetype, ierr)
  call MPI_Type_free(mpi_subarray, ierr)

  call MPI_Finalize(ierr)

end subroutine cafwrite

end module cafmpiio
