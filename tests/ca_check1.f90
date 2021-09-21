program z
use iso_fortran_env
implicit none
integer, parameter :: dp = real64
real( kind=dp)     :: time1, time2

integer :: x(10)[*], img, nimgs, i
  img = this_image()
nimgs = num_images()
    x = img 

sync all

if ( img .eq. 1 ) then
  do i = 2, nimgs
!    time1 = mytime()
!    x = x(:)[i]
!    time2 = mytime()
!    write (*,"(a,g)") "Remote read took, s : ", time2-time1
    time1 = mytime()
    x(:)[i] = x
    time2 = mytime()
    write( *,*) "Remote write took, s : ", time2-time1
    write( *, "(a,i0,a,10(i0))" ) "img: ", i, "x:", x(:)[i]
  end do
end if

sync all

write (*,"(a,i0,a,i0,a)") "Image: ", img, " out of ", nimgs, &
 " completed ok"

contains

function mytime() result (tseconds)
        real( dp ) :: tseconds
  integer( INT64 ) :: tsec, trate
  CALL SYSTEM_CLOCK( count=tsec, count_rate=trate )
  tseconds = real(tsec,kind=dp) / real(trate,kind=dp)
end function mytime 

end program z
