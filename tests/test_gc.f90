!*robodoc*u* tests/testgc
!  NAME
!    testgc
!  SYNOPSIS

!$Id: test_gc.f90 380 2017-03-22 11:03:09Z mexas $

program testgc

!  PURPOSE
!    Checking: grain connectivity, cgca_m2gb
!  DESCRIPTION
!    This program reads the global grain connectivity from stdin
!    and dumps a number of neighbours for each grain to stdout.
!  AUTHOR 
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  USES
!    none
!  USED BY
!    Part of CGPACK test suite
!  SOURCE

implicit none

integer :: errstat=0, data1=0, data2=0, value=0, total=0, grain=0

read (unit=*,fmt=*,iostat=errstat) data1
if (errstat .ne. 0) stop "Cannot read data1"
value=1
total=1
grain=1

! The first value must always be grain 1.
! Stop with an error, if not.
if (data1 .ne. grain) stop "The first grain must be number 1"

do
 read (unit=*,fmt=*,iostat=errstat) data2
 if (errstat .ne. 0) then
  ! this means end of file
  write (*,*) value
  exit
 end if
 if ( data2 .eq. data1 ) then
  ! same grain
  value = value+1 
 else
  write (*,*) value
  value=1
  grain=grain+1
  ! The value must match grain number. Stop if not.
  if (data2 .ne. grain) then
   write (*,*) "This should have been grain", grain
   write (*,*) "Instead I read", data2
   stop "This is an error, and I can do no more for you..."
  end if
 data1 = data2
 end if
 total = total +1
end do

write (*,*) "Total values read: ", total
write (*,*) "Total grains: ", grain

end program testgc

!*roboend*
