!*robodoc*f* cgca_m2out/m2out_sm1
!  NAME
!    m2out_sm1
!  SYNOPSIS

!$Id: m2out_sm1.f90 380 2017-03-22 11:03:09Z mexas $

submodule ( cgca_m2out) m2out_sm1

!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  DESCRIPTION
!    Submodule with output routines using Cray parallel
!    coarray IO extensions - basically writing to a shared
!    direct access file.
!  CONTAINS
!    cgca_pc
!  USES
!    Parameters and variables of module cgca_m2out by host association.
!  USED BY
!    Module cgca_m2out
!  SOURCE

implicit none

contains

!*roboend*


!*robodoc*s* m2out_sm1/cgca_pc
!  NAME
!    cgca_pc
!  SYNOPSIS
 
module procedure cgca_pc

!  INPUTS
!    See the interface block in the parent module cgca_m2out.
!  OUTPUTS
!    See the interface block in the parent module cgca_m2out.
!  SIDE EFFECTS
!    creates a binary file from *all images* and writes coarray to it
!  DESCRIPTION
!    Parallel Cray output routine. This routine works only on Cray
!    systems using non-standandard proprietary Cray extension, beware!
!    Refer to Cray "Enhanced I/O: Using the assign Environment",
!    section 13 of Cray Fortran Reference Manual, S-3901-83, June 2014
!    or later version:
!    http://docs.cray.com/cgi-bin/craydoc.cgi?mode=View;id=S-3901-83
!
!    First Cray assign enviroment is established, setting fname
!    for shared IO access. Then image 1 writes the last record in
!    the file. After that all images know the length of the file.
!    Then all images can write their data in the correct place in
!    this file in parallel. The expectation is that this is faster
!    than a serial writer, cgca_swci.
!
!    All images must call this routine!
!    Performance is not guaranteed, use with caution!
!  USES
!    none
!  USED BY
!    parent module cgca_m2out
!  SOURCE

integer :: errstat=0, img, nimgs, i2, i3, fu, reclen, recnum,          &
  help1, help2, help3,                                                 &
  cosub(3),  & ! set of cosubscripts for this image
  lb(4),     & ! lower bounds   of the coarray
  ub(4),     & ! upper bounds   of the coarray
  spsz(4),   & ! size of the space array along each dimension
  lcob(3),   & ! lower cobounds of the coarray
  ucob(3),   & ! upper cobounds of the coarray
  cosz(3)      ! size of the space coarray along each dimension

! Assume the coarray has halos. Don't write those.
   lb = lbound( coarray ) + 1
   ub = ubound( coarray ) - 1
 spsz = ub - lb + 1
 lcob = lcobound( coarray )
 ucob = ucobound( coarray )
 cosz = ucob - lcob + 1
  img = this_image()
cosub = this_image( coarray )
nimgs = num_images()

! Initialise Cray assign environment
! -m on - "Special handling of direct access file that will be accessed
!          concurrently by several processes or tasks"
! -F system - no buffering
! fname - this assign enviroment will apply only to file name "fname".
call asnfile( trim(fname), '-m on -F system', errstat )

! Need to set up record length
inquire( iolength=reclen ) coarray( lb(1):ub(1), 1, 1, stype) [1,1,1]
if ( img .eq. 1 ) then
  write (*,*) "INFO: cgca_pc: asnfile errstat:", errstat
  write (*,*) "INFO: cgca_pc: record length:", reclen
  write (*,*) "INFO: cgca_pc: last record num:", spsz(2)*spsz(3)*nimgs
end if

!give fu any value
fu=-2
! open file on image 1, write the last record to it and close it
if ( img .eq. nimgs ) then
  open( newunit=fu, file=trim(fname), access="direct", recl=reclen,    &
        form="unformatted", status="replace" )
  recnum = spsz(2) * spsz(3) * nimgs 
  write( fu, rec= recnum ) coarray( lb(1):ub(1), ub(2), ub(3), stype )
  close( fu )
end if

! all images wait until the file size is known
sync all

! open file on all images
open( unit=fu, file=fname, access="direct", recl=reclen,               &
      form="unformatted", status="old" )

! Calculate intermediate variables to reduce the FLOPs
! The exact expression for recnum 
!     ( (cosub(3)-1) * spsz(3) + i3 - 1 ) * cosz(2) * spsz(2) * cosz(1) &
!   + ( (cosub(2)-1) * spsz(2) + i2 - 1 ) * cosz(1) + cosub(1)
help3 = (cosub(3)-1) * spsz(3) - 1
help2 = (cosub(2)-1) * spsz(2) - 1
help1 = cosz(2) * spsz(2) * cosz(1)

! write data
do i3 = lb(3), ub(3)
do i2 = lb(2), ub(2)
  recnum = (help3 + i3) * help1 + (help2 + i2) * cosz(1) + cosub(1)
  write( unit=fu, rec= recnum ) coarray( lb(1):ub(1), i2, i3, stype )
end do
end do

! flush data
flush( unit=fu )

! wait till all images wrote data and flushed
sync all

! close the file
close( fu )

end procedure cgca_pc

!*roboend*

end submodule m2out_sm1
