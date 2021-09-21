! modified as suggested by Ron:
!   https://software.intel.com/en-us/forums/topic/489052
! to check his theory that remote/read writes
! can be executed out of order, even though they
! are on the same image and within the same segment
 
program z
use iso_fortran_env
implicit none
integer, parameter :: dp = real64
real( kind=dp) :: time1, time2

integer :: x(10)[*], img, nimgs, i
  img = this_image()
nimgs = num_images()
    x = img 

sync all

do i = 2, nimgs
  if ( img .eq. 1 ) then
    time1 = mytime()
    x = x(:)[i]
    sync images ( i )
    time2 = mytime()
    write( *, * ) "Remote read took, s : ", time2-time1
    write( *, "(a,i0,a,10(i0))") "img: ", i, "x:", x(:)[i]
  else if ( img .eq. i ) then
    sync images( 1 )
  end if

  if ( img .eq. 1 ) then
    time1 = mytime()
    x(:)[i] = x
    sync images ( i )
    time2 = mytime()
    write( *, * ) "Remote write took, s : ", time2-time1
    write (*,"(a,i0,a,10(i0))") "img: ", i, "x:", x(:)[i]
  else if ( img .eq. i ) then
    sync images( 1 )
  end if
end do

sync all

write (*,"(a,i0,a,i0,a)") "Image: ", img, " out of ", nimgs, "completed ok"

contains

function mytime() result (tseconds)
        real( dp ) :: tseconds
  integer( INT64 ) :: tsec, trate
  CALL SYSTEM_CLOCK( count=tsec, count_rate=trate )
  tseconds = real(tsec,kind=dp) / real(trate,kind=dp)
end function mytime 

end program z
