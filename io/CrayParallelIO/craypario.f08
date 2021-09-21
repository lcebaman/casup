program par
implicit none

integer, parameter :: idef = selected_int_kind(10), arrsize=10
integer( kind=idef ) :: img, i( arrsize )[*]
integer :: reclen, fu, errstat
character( len=3 ), parameter :: fname="zzz"

  img = this_image()
    i = img

call asnfile( fname, '-m on -F system', errstat )

if ( img .eq. 1 ) then
  write (*,*) "asnfile errstat=", errstat
end if

inquire( iolength=reclen ) i
write (*,*) "Record length", reclen

if ( img .eq. 1 ) then
  open( newunit=fu, file=fname, access="direct", recl=reclen,          &
        form="unformatted", status="replace" )
  write( fu, rec= num_images() ) i
  close( fu )
end if

sync all

open( fu, file=fname, access="direct", recl=reclen,                    &
      form="unformatted", status="old" )
write( fu, rec=img ) i
flush( fu )

sync all

close( fu )

sync all

open( fu, file=fname, access="direct", recl=reclen,                    &
      form="unformatted", status="old" )
read( fu, rec=img ) i
write (*,*) "img", img, "i", i
flush( fu )

sync all

end program par
