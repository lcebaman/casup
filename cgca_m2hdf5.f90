!*robodoc*m* CGPACK/cgca_m2hdf5
!  NAME
!    cgca_m2hdf5
!  SYNOPSIS

!$Id: cgca_m2hdf5.f90 414 2017-05-18 12:55:39Z mexas $

module cgca_m2hdf5

!  AUTHOR
!   Luis Cebamanos
!  COPYRIGHT
!    See LICENSE
!  DESCRIPTION
!    Module with hdf5 related routines
!  CONTAINS
!    Subroutines: cgca_pswci4
!  USES
!    cgca_m1co
!  USED BY
!    cgca
!  SOURCE

use cgca_m1co
use mpi
use hdf5
implicit none
private
public :: cgca_pswci4

contains

!*roboend*


!*robodoc*s* cgca_m2hdf5/cgca_pswci4
!  NAME
!    cgca_pswci4
!  SYNOPSIS

subroutine cgca_pswci4( coarray, stype, fname )

!  INPUTS

    integer( kind=iarr ), allocatable, intent( in ) ::                 &
      coarray(:,:,:,:)[:,:,:]
    integer( kind=idef ),intent( in ) :: stype
    character( len=* ), intent( in ) :: fname

!  OUTPUTS
!    None
!  SIDE EFFECTS
!    A single binary file is created using hdf5
!    with contents of coarray.
!  AUTHOR
!    Luis Cebamanos, adapted from the code
!    written by David Henty, EPCC
!  COPYRIGHT
!   Note that this routine has special Copyright conditions.
!
!    !----------------------------------------------------------------------------!
!    !                                                                            !
!    !  hdf5 routine for Fortran Coarrays                                         !
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
!    Parallel Stream Write Coarray of Integers (PSWC), number 4.
!    Dump the coarray to file in a binary file in HDF5 format:
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
!    cgca_m1co, MPI library, hdf5 library
!  USED BY
!    none, end user.
!  SOURCE

   integer, parameter :: totdim = 4, arrdim = totdim-1, coardim = 3

   integer :: img, nimgs, comm, ierr=0, rank=0, mpisize=0
          
   integer, dimension(totdim) :: asizehal
   integer, dimension(arrdim) :: arrsize, arstart, artsize
   integer, dimension(coardim) :: coarsize, copos
   
   integer :: info = MPI_INFO_NULL
   integer(hsize_t), dimension(coardim) :: dimsf  ! dataset dimensions.
  
   character(len=8), parameter :: dsetname = "IntArray" ! Dataset name
   
   integer(hid_t) :: file_id       ! file identifier 
   integer(hid_t) :: dset_id       ! dataset identifier 
   integer(hid_t) :: filespace     ! dataspace identifier in file 
   integer(hid_t) :: memspace      ! dataspace identifier in memory
   integer(hid_t) :: plist_id      ! property list identifier 
   
   integer(hsize_t), dimension(coardim) :: count  
   integer(hssize_t), dimension(coardim) :: offset 
   
!   integer :: ncid, varid, oldmode, dimids(coardim)
!   integer :: x_dimid, y_dimid, z_dimid

   img = this_image()
 nimgs = num_images()

   comm = MPI_COMM_WORLD
   call MPI_Comm_size( comm, mpisize, ierr )
   call MPI_Comm_rank( comm, rank, ierr )
 
   ! Sanity check
   if ( mpisize .ne. nimgs .or. rank .ne. img-1 ) then
      write(*,*) 'ERROR: cgca_m2hdf5/cgca_pswci4: MPI/coarray mismatch'
      error stop
   end if

   asizehal(:) = shape( coarray )
   copos(:) = this_image( coarray )

   ! Subtract halos
   arrsize(:) = asizehal(1:arrdim) - 2
   coarsize(:) = ucobound(coarray) - lcobound(coarray) + 1
 
   
 
   ! Does the array fit exactly?
   if ( product( coarsize ) .ne. nimgs) then
      write(*,*) 'ERROR: cgca_m2hdf5/cgca_pswci4: non-conforming coarray'
      error stop
   end if
 
   
   ! This is the global arrayy
   artsize(:) = arrsize(:) * coarsize(:)

   ! Correspondent portion of this global array 
   arstart(:) = arrsize(:) * (copos(:)-1) + 1 ! Use Fortran indexing

   dimsf(:) = artsize(:)

   count(1) = arrsize(1)     ! Defines the number of values each proc dumps to 
   count(2) = arrsize(2)
   count(3) = arrsize(3)                                   ! the HDF5 file. 

! ! debug
! write (*,*) "hdf5 - image",img, "asizehal", asizehal, "copos", copos,         &
!  "arrsize", arrsize, "coarsize", coarsize,                            &
!  "artsize", artsize, "arstart", arstart, "stype", stype

  offset(:) = (copos(:)-1) * count(:)  ! Defines the offset used in the HDF5 file

  ! Initialize FORTRAN predefined datatypes
  CALL h5open_f(ierr) 

  ! Setup file access property list with parallel I/O access.
  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
  CALL h5pset_fapl_mpio_f(plist_id, comm, info, ierr)
  ! debug
  if (img .eq. 1) then
     write (*,*) "Setup file access property list"
  end if

  ! Create the file collectively.
  CALL h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, ierr, &
      access_prp = plist_id)
  if (ierr .ne. 0) then
     write(0,*) 'Unable to open: ', trim(fname), ': ', ierr
     call mpi_abort(MPI_COMM_WORLD, 1, ierr)
  endif
  if (img .eq. 1) then
     write (*,*) "Created a file collectively"
  end if
  CALL h5pclose_f(plist_id, ierr)
  if (img .eq. 1) then
     write (*,*) "Close property list"
  end if
  ! Create the data space for the  dataset. 
  CALL h5screate_simple_f(coardim, dimsf, filespace, ierr)
  if (img .eq. 1) then
     write (*,*) "Create data space for the dataset"
  end if
  ! Create the dataset with default properties.
  CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_INTEGER, filespace, &
                      dset_id, ierr)
    if (img .eq. 1) then
       write(*,*) "Create dataset with default properties"
    end if
    CALL h5sclose_f(filespace, ierr)
    if (img .eq. 1) then
       write(*,*) "Close the dataset" 
    end if

  ! Each process defines dataset in memory and writes it to the hyperslab
  ! in the file. 
  CALL h5screate_simple_f(coardim, count, memspace, ierr) 
  if (img .eq. 1) then
     write(*,*) "Each process defines dataset in mem"
  end if
  ! Select hyperslab in the file.
  CALL h5dget_space_f(dset_id, filespace, ierr)
  CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, &
      count, ierr)
  if (img .eq. 1) then
     write(*,*) "Selects a hyperslab in the file"
  end if
     
  ! Create property list for collective dataset write
  CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr) 
  CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
  if (img .eq. 1) then
     write(*,*) "Crete property list for coll dataset"
  end if
  
  ! Write the dataset collectively. 
  CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, coarray(1:arrsize(1),    & 
       1:arrsize(2), 1:arrsize(3), stype), dimsf, ierr,             &
       file_space_id = filespace, mem_space_id = memspace,          &
       xfer_prp = plist_id)
  if (img .eq. 1) then
     write(*,*) "Write dataset colectively"
  end if
  
  ! Close dataspaces.
  CALL h5sclose_f(filespace, ierr)
  CALL h5sclose_f(memspace, ierr)
  if (img .eq. 1) then
     write(*,*) "Close stuff"
  end if
  ! Close the dataset and property list.
  CALL h5dclose_f(dset_id, ierr)
  CALL h5pclose_f(plist_id, ierr)

  ! Close the file.
  CALL h5fclose_f(file_id, ierr)

  ! Close FORTRAN predefined datatypes.
  CALL h5close_f(ierr)

  
end subroutine cgca_pswci4


!*roboend*

end module cgca_m2hdf5

