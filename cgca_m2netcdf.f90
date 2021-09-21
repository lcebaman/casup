!*robodoc*m* CGPACK/cgca_m2netcdf
!  NAME
!    cgca_m2netcdf
!  SYNOPSIS

!$Id: cgca_m2netcdf.f90 380 2017-03-22 11:03:09Z mexas $

module cgca_m2netcdf

!  AUTHOR
!   Luis Cebamanos
!  COPYRIGHT
!    See LICENSE
!  DESCRIPTION
!    Module with netCDF related routines
!  CONTAINS
!    Subroutines: cgca_pswci3
!  USES
!    cgca_m1co
!  USED BY
!    cgca
!  SOURCE

use cgca_m1co
use mpi
use netcdf
implicit none
private
public :: cgca_pswci3

contains

!*roboend*


!*robodoc*s* cgca_m2netcdf/cgca_pswci3
!  NAME
!    cgca_pswci3
!  SYNOPSIS

subroutine cgca_pswci3( coarray, stype, fname )

!  INPUTS

    integer( kind=iarr ), allocatable, intent( in ) ::                 &
      coarray(:,:,:,:)[:,:,:]
    integer( kind=idef ),intent( in ) :: stype
    character( len=* ), intent( in ) :: fname

!  OUTPUTS
!    None
!  SIDE EFFECTS
!    A single binary file is created using netcdf
!    with contents of coarray.
!  AUTHOR
!    Luis Cebamanos, adapted from the code
!    written by David Henty, EPCC
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
!    cgca_m1co, MPI library, netCDF library
!  USED BY
!    none, end user.
!  SOURCE

integer, parameter :: totdim = 4, arrdim = totdim-1, coardim = 3

integer :: img, nimgs, comm, ierr=0, rank=0, mpisize=0
          
integer, dimension(totdim) :: asizehal
integer, dimension(arrdim) :: arrsize, arstart, artsize
integer, dimension(coardim) :: coarsize, copos

integer :: ncid, varid, dimids(arrdim)
integer :: x_dimid=0, y_dimid=0, z_dimid=0

  img = this_image()
nimgs = num_images()

asizehal(:) = shape( coarray )
   copos(:) = this_image( coarray )

! Subtract halos
 arrsize(:) = asizehal(1:arrdim) - 2
coarsize(:) = ucobound(coarray) - lcobound(coarray) + 1

! Does the array fit exactly?
if ( product( coarsize ) .ne. nimgs) then
 write(*,*) 'ERROR: cgca_m2netcdf/cgca_pswci3: non-conforming coarray'
 error stop
end if

comm = MPI_COMM_WORLD
call MPI_Comm_size( comm, mpisize, ierr )
call MPI_Comm_rank( comm, rank, ierr )

! Sanity check
if ( mpisize .ne. nimgs .or. rank .ne. img-1 ) then
  write(*,*) 'ERROR: cgca_m2netcdf/cgca_pswci3: MPI/coarray mismatch'
  error stop
end if
     

! This is the global arrayy
artsize(:) = arrsize(:) * coarsize(:)

! Correspondent portion of this global array 
arstart(:) = arrsize(:) * (copos(:)-1) + 1 ! Use Fortran indexing

! ! debug
! write (*,*) "netCDF-image",img, "asizehal", asizehal, "copos", copos,         &
!  "arrsize", arrsize, "coarsize", coarsize,                            &
!  "artsize", artsize, "arstart", arstart, "stype", stype

! Create (i.e. open) the netCDF file. The NF90_NETCDF4 flag causes a 
! HDF5/netCDF-4 type file to be created. The comm and info parameters 
! cause parallel I/O to be enabled. 
call check( nf90_create(fname, ior(nf90_netcdf4,nf90_mpiio), ncid, &
       comm=comm, info=MPI_INFO_NULL))
  

! Define the dimensions. NetCDF returns an ID for each. Any 
! metadata operations must take place on ALL processors
call check( nf90_def_dim(ncid, "x", artsize(1), x_dimid) )
call check( nf90_def_dim(ncid, "y", artsize(2), y_dimid) )
call check( nf90_def_dim(ncid, "z", artsize(3), z_dimid) )

! The dimids array is used to pass the ID's of the dimensions of 
! the variables. 
dimids = (/ x_dimid , y_dimid, z_dimid /)

! Define the variable. The type of the variable in this case is
! NF90_INT (4-byte int).
call check( nf90_def_var(ncid, "data", NF90_INT, dimids, varid) )

! Make sure file it not filled with default values which doubles wrote volume
call check ( nf90_def_var_fill(ncid, varid, 1, 1) )

! End define mode. This tells netCDF we are done defining
! metadata. This operation is collective and all processors will
! write their metadata to disk.
call check( nf90_enddef(ncid) )

! Parallel access
call check( nf90_var_par_access(ncid, varid, nf90_collective) )

! Write the data to file, start will equal the displacement from the 
! start of the file and count is the number of points each proc writes. 
call check1( nf90_put_var(ncid, varid, coarray(1:arrsize(1), 1:arrsize(2), 1:arrsize(3),stype), &
     start = arstart, count = arrsize) )

! Close the file. This frees up any internal netCDF resources
! associated with the file, and flushes any buffers.
call check( nf90_close(ncid) )

end subroutine cgca_pswci3

subroutine check(status )
  integer, intent ( in) :: status
  if(status /= nf90_noerr) then 
     print *, trim(nf90_strerror(status))
     stop
  end if
end subroutine check

subroutine check1(status)
  integer, intent ( in) :: status
  if(status /= nf90_noerr) then 
     print *, "put_var: ",trim(nf90_strerror(status))
     stop
  end if
end subroutine check1

!*roboend*

end module cgca_m2netcdf

