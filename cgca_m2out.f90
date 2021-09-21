!*robodoc*m* CGPACK/cgca_m2out
!  NAME
!    cgca_m2out
!  SYNOPSIS

!$Id: cgca_m2out.f90 380 2017-03-22 11:03:09Z mexas $

module cgca_m2out

!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  DESCRIPTION
!    Module dealing with output 
!  CONTAINS
!    Subroutines: cgca_swci, cgca_fwci.
!    Submodules: m2out_sm1, m2out_sm2_mpi.
!  USES
!    cgca_m1co
!  USED BY
!    cgca
!  SOURCE

use cgca_m1co
implicit none

private
public :: cgca_swci, cgca_pc, cgca_pswci, cgca_fwci

interface
  module subroutine cgca_pc( coarray, stype, fname )
    ! In submodule m2out_sm1.f90
    ! coarray - what array to dump
    ! stype - what cell state type to dump
    ! fname - what file name to use
    integer( kind=iarr ), allocatable, intent( in ) ::                 &
      coarray(:,:,:,:)[:,:,:]
    integer( kind=idef ),intent( in ) :: stype
    character( len=* ), intent( in ) :: fname
  end subroutine cgca_pc

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
  end subroutine cgca_pswci

end interface

contains

!*roboend*


!*robodoc*s* cgca_m2out/cgca_swci
!  NAME
!    cgca_swci
!  SYNOPSIS

subroutine cgca_swci( coarray, stype, iounit, fname )

!  INPUTS
 
integer( kind=iarr ),allocatable,intent( in ) :: coarray(:,:,:,:)[:,:,:]
integer( kind=idef ),intent( in ) :: stype, iounit
character( len=* ),intent( in ) :: fname

!  SIDE EFFECTS
!    A single binary file is created on image 1 with contents of coarray.
!  DESCRIPTION
!    Stream Write Coarray of Integers:
!    - coarray - what array to dump
!    - stype - what cell state type to dump
!    - iounit - which I/O unit to use
!    - fname - what file name to use
!  NOTES
!    All images call this routine!
!    However only image 1 does all the work.
!    The other images are waiting. 
!  USES
!    none
!  USED BY
!    none, end user.
!  SOURCE

integer :: errstat, coi1, coi2, coi3, i2, i3, &
  lb(4),   & ! lower bounds   of the coarray
  ub(4),   & ! upper bounds   of the coarray
  lcob(3), & ! lower cobounds of the coarray
  ucob(3)    ! upper cobounds of the coarray

! Only image1 does this. All other images do nothing.
! So sync all probably should be used after a call to
! this routine in the program.

main: if ( this_image() .eq. 1 ) then
  errstat = 0

  ! Assume the coarray has halos. Don't write those.
    lb = lbound( coarray ) + 1
    ub = ubound( coarray ) - 1
  lcob = lcobound( coarray )
  ucob = ucobound( coarray )

!write (*,*) "DEBUG: cgca_swci: lb: " , lb
!write (*,*) "DEBUG: cgca_swci: ub: " , ub
!write (*,*) "DEBUG: cgca_swci: lcob: " , lcob
!write (*,*) "DEBUG: cgca_swci: ucob: " , ucob

  open( unit=iounit, file=fname, form="unformatted", access="stream",  &
        status="replace", iostat=errstat )
  if ( errstat .ne. 0 ) then
    write (*,'(a,i0)') "ERROR: cgca_swci/cgca_m2out: " //              &
      "open file for writing: err code: ", errstat
    error stop
  end if

!write (*,*) "DEBUG: cgca_swci: starting data output"

 ! nested loops for writing in correct order from all images
 do coi3 = lcob(3), ucob(3)
   do i3 = lb(3), ub(3)
     do coi2 = lcob(2), ucob(2)
       do i2 = lb(2), ub(2)
         do coi1 = lcob(1), ucob(1)

  write( unit=iounit, iostat=errstat )                                 &
    coarray( lb(1):ub(1), i2, i3, stype ) [ coi1, coi2, coi3 ]
  if (errstat .ne. 0) then
    write (*,'(a,i0)') "ERROR: cgca_swci/cgca_m2out: " //              &
      "write: err code: ", errstat
    error stop
  end if

!write (*,*) "DEBUG: cgca_swci: wrote cells with: i2, i3, coi1, coi2, coi3", i2, i3, coi1, coi2, coi3

         end do
       end do
     end do
   end do
 end do

!write (*,*) "DEBUG: cgca_swci: finished data output"

 close( unit=iounit, iostat=errstat )
 if ( errstat .ne. 0 ) then
   write (*,'(a,i0)') "ERROR: cgca_swci/cgca_m2out: " //               &
     "close file: err code: ", errstat
  error stop
 end if

end if main

end subroutine cgca_swci

!*roboend*


!*robodoc*s* cgca_m2out/cgca_fwci
!  NAME
!    cgca_fwci
!  SYNOPSIS

subroutine cgca_fwci( coarray, stype, fname )

!  INPUTS
 
integer( kind=iarr ),allocatable,intent( in ) :: coarray(:,:,:,:)[:,:,:]
integer( kind=idef ),intent( in ) :: stype
character( len=* ),intent( in ) :: fname

!  SIDE EFFECTS
!    A single formatted file is created on image 1
!    with contents of one layer of the coarray.
!  DESCRIPTION
!    Formatted Write Coarray of Integers:
!    - coarray - what array to dump
!    - stype - what cell state type to dump
!    - fname - what file name to use
!  NOTES
!    The main purpose of this routine is to provide an easy means
!    for checking whether the results are reproducible. Arguably
!    a formatted output gives a better clue than a binary file.
!    This routine is *very* slow - writing a single cell value
!    per line.
!
!    All images call this routine!
!    However only image 1 does all the work.
!    The other images are waiting. 
!  USES
!    none
!  USED BY
!    none, end user.
!  SOURCE

integer :: errstat, coi1, coi2, coi3, i1, i2, i3, funit,               &
  lb(4),   & ! lower bounds   of the coarray
  ub(4),   & ! upper bounds   of the coarray
  lcob(3), & ! lower cobounds of the coarray
  ucob(3)    ! upper cobounds of the coarray

! Only image1 does this. All other images do nothing.
! So sync all probably should be used after a call to
! this routine in the program.

! Give funit any value, just to avoid compiler warnings
funit = 111

main: if ( this_image() .eq. 1 ) then
  errstat = 0

  ! Assume the coarray has halos. Don't write those.
    lb = lbound( coarray ) + 1
    ub = ubound( coarray ) - 1
  lcob = lcobound( coarray )
  ucob = ucobound( coarray )

  open( newunit=funit, file=fname, form="formatted",                   &
        access="sequential", status="replace", iostat=errstat )
  if ( errstat .ne. 0 ) then
    write (*,'(a,i0)') "ERROR: cgca_fwci/cgca_m2out: " //              &
      "open file for writing: err code: ", errstat
    error stop
  end if

  ! nested loops for writing out from all images
  do coi3 = lcob(3), ucob(3)
  do coi2 = lcob(2), ucob(2)
  do coi1 = lcob(1), ucob(1)
    do i3 = lb(3), ub(3)
    do i2 = lb(2), ub(2)
    do i1 = lb(1), ub(1)
      write( unit=funit, iostat=errstat, fmt="(7(i0,tr1))" )           &
        coi1, coi2, coi3, i1, i2, i3, coarray( i1, i2, i3, stype )     &
                                                [ coi1, coi2, coi3 ]
      if (errstat .ne. 0) then
         write (*,'(a,i0)') "ERROR: cgca_fwci/cgca_m2out: " //         &
           "write: err code: ", errstat
         error stop
      end if
    end do
    end do
    end do
  end do
  end do
  end do

  close( unit=funit, iostat=errstat )
  if ( errstat .ne. 0 ) then
    write (*,'(a,i0)') "ERROR: cgca_fwci/cgca_m2out: " //              &
      "close file: err code: ", errstat
    error stop
  end if

end if main

end subroutine cgca_fwci

!*roboend*


end module cgca_m2out
