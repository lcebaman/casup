!*robodoc*f* cgca_m2hx/m2hx_hxic
!  NAME
!    m2hx_hxic
!  SYNOPSIS

!$Id: m2hx_hxic.f90 430 2017-06-30 07:39:43Z mexas $

submodule ( cgca_m2hx ) m2hx_hxic

!  DESCRIPTION
!    Submodule of cgca_m2hx with a hx subroutine.
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENCE
!  CONTAINS
!    cgca_hxic
!  USES
!    Variables and parameters from the parent module cgca_m2hx.
!  USED BY
!    The parent module cgca_m2hx.
!  SOURCE

implicit none

contains

!*roboend*


!*robodoc*s* cgca_m2hx/cgca_hxic
!  NAME
!    cgca_hxic
!  SYNOPSIS

module procedure cgca_hxic

!  INPUT
!    See the parent module
!  OUTPUT
!    See INPUT
!  SIDE EFFECTS
!    None
!  DESCRIPTION
!    This routine checks that a prior hx call was done flagly,
!    i.e. that the halo cell states are consistent with the states
!    of the corresponding real boundary cells. This routine can be
!    called for any internal hx algorithm, e.g. cgca_hxi or cgca_hxir.
!  NOTES
!    All images must call this routine! Lots of remote calls.
!  USES
!    All data objects from parent module cgca_m2hx by host association.
!  USED BY
!    module cgca_m2hx
!  SOURCE

integer ::                                                             &
  lbv(4)      , & ! lower bounds of the "virtual" coarray
  ubv(4)      , & ! upper bounds of the "virtual" coarray
  lbr(4)      , & ! lower bounds of the "real" coarray, lbv+1
  ubr(4)      , & ! upper bounds of the "real" coarray, ubv-1
  lcob(3)     , & ! lower cobounds of the coarray
  ucob(3)     , & ! upper cobounds of the coarray
  imgpos(3)       ! position of the image in a coarray grid

! Start with 0. Any error must resuult in a positive value.
flag = 0

! check for allocated 
if ( .not. allocated( coarray ) )                                      &
  error stop "ERROR: cgca_hxic/m2hx_hxic: coarray is not allocated"

   lbv = lbound( coarray )
   ubv = ubound( coarray )
   lbr = lbv + 1
   ubr = ubv - 1
  lcob = lcobound( coarray )
  ucob = ucobound( coarray )
imgpos = this_image( coarray )

! Make sure only the virtual (halo) arrays are assigned to.
! The real array values must never appear on the left
! hand side of the assignment expressions.
! The halo exchange process is copying real array values into halos.
! There must not ever be copying real values
! to real, or halo to halo, or halo to real.
! Also, only local array must appear on the left.
! We are assigning values to the local array's virtual
! (halo) cells using values from real cells from arrays
! in other images.

! Check 2D halos in direction 1

! op 1
if ( imgpos(1) .ne. lcob(1) ) then
  if ( any(                                                            &
     coarray( lbv(1) , lbr(2) : ubr(2) , lbr(3) : ubr(3) , : ) .ne.    &
     coarray( ubr(1) , lbr(2) : ubr(2) , lbr(3) : ubr(3) , : )         &
       [ imgpos(1)-1 , imgpos(2) , imgpos(3) ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

! op 2
if ( imgpos(1) .ne. ucob(1) ) then
  if ( any(                                                            &
    coarray( ubv(1) , lbr(2) : ubr(2) , lbr(3) : ubr(3) , : ) .ne.     &
    coarray( lbr(1) , lbr(2) : ubr(2) , lbr(3) : ubr(3) , : )          &
      [ imgpos(1)+1 , imgpos(2) , imgpos(3) ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

! exchange 2D halos in direction 2

! op 3
if ( imgpos(2) .ne. lcob(2) ) then
  if ( any(                                                            &
    coarray( lbr(1) : ubr(1) , lbv(2) , lbr(3) : ubr(3) , : ) .ne.     &
    coarray( lbr(1) : ubr(1) , ubr(2) , lbr(3) : ubr(3) , : )          &
      [ imgpos(1) , imgpos(2)-1 , imgpos(3) ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

! op 4
if ( imgpos(2) .ne. ucob(2) ) then
  if ( any(                                                            &
    coarray( lbr(1) : ubr(1) , ubv(2) , lbr(3) : ubr(3) , : ) .ne.     &
    coarray( lbr(1) : ubr(1) , lbr(2) , lbr(3) : ubr(3) , : )          &
      [ imgpos(1) , imgpos(2)+1 , imgpos(3) ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

! exchange 2D halos in direction 3

! op 5
if ( imgpos(3) .ne. lcob(3) ) then
  if ( any(                                                            &
    coarray( lbr(1) : ubr(1) , lbr(2) : ubr(2) , lbv(3) , : ) .ne.     &
    coarray( lbr(1) : ubr(1) , lbr(2) : ubr(2) , ubr(3) , : )          &
      [ imgpos(1) , imgpos(2) , imgpos(3)-1 ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

! op 6
if ( imgpos(3) .ne. ucob(3) ) then
  if ( any(                                                            &
  coarray( lbr(1) : ubr(1) , lbr(2) : ubr(2) , ubv(3) , : ) .ne.       &
  coarray( lbr(1) : ubr(1) , lbr(2) : ubr(2) , lbr(3) , : )            &
    [ imgpos(1) , imgpos(2) , imgpos(3)+1 ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

! exchange 1D halos parallel to direction 3

! op 7
if ( imgpos(1) .ne. lcob(1) .and. imgpos(2) .ne. lcob(2) ) then
  if ( any(                                                            &
    coarray( lbv(1) , lbv(2) , lbr(3) : ubr(3) , : ) .ne.              &
    coarray( ubr(1) , ubr(2) , lbr(3) : ubr(3) , : )                   &
      [ imgpos(1)-1 , imgpos(2)-1 , imgpos(3) ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

! op 8
if ( imgpos(1) .ne. ucob(1) .and. imgpos(2) .ne. lcob(2) ) then
  if ( any(                                                            &
    coarray( ubv(1) , lbv(2) , lbr(3) : ubr(3) , :) .ne.               &
    coarray( lbr(1) , ubr(2) , lbr(3) : ubr(3) , :)                    &
      [ imgpos(1)+1 , imgpos(2)-1 , imgpos(3) ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

! op 9
if ( imgpos(1) .ne. ucob(1) .and. imgpos(2) .ne. ucob(2) ) then
  if ( any(                                                            &
    coarray( ubv(1) , ubv(2) , lbr(3) : ubr(3) , : ) .ne.              &
    coarray( lbr(1) , lbr(2) , lbr(3) : ubr(3) , : )                   &
      [ imgpos(1)+1 , imgpos(2)+1 , imgpos(3) ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

! op 10
if ( imgpos(1) .ne. lcob(1) .and. imgpos(2) .ne. ucob(2) ) then
  if ( any(                                                            &
    coarray( lbv(1) , ubv(2) , lbr(3) : ubr(3) , : ) .ne.              &
    coarray( ubr(1) , lbr(2) , lbr(3) : ubr(3) , : )                   &
      [ imgpos(1)-1 , imgpos(2)+1 , imgpos(3) ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

! exchange 1D halos parallel to direction 1

! op 11
if ( imgpos(2) .ne. lcob(2) .and. imgpos(3) .ne. lcob(3) ) then
  if ( any(                                                            &
    coarray( lbr(1) : ubr(1) , lbv(2) , lbv(3) , : ) .ne.              &
    coarray( lbr(1) : ubr(1) , ubr(2) , ubr(3) , : )                   &
      [ imgpos(1) , imgpos(2)-1 , imgpos(3)-1 ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

! op 12
if ( imgpos(2) .ne. lcob(2) .and. imgpos(3) .ne. ucob(3) ) then
  if ( any(                                                            &
    coarray( lbr(1) : ubr(1) , lbv(2) , ubv(3) , : ) .ne.              &
    coarray( lbr(1) : ubr(1) , ubr(2) , lbr(3) , : )                   &
      [ imgpos(1) , imgpos(2)-1 , imgpos(3)+1 ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

! op 13
if ( imgpos(2) .ne. ucob(2) .and. imgpos(3) .ne. ucob(3) ) then
  if ( any(                                                            &
    coarray( lbr(1) : ubr(1) , ubv(2) , ubv(3) , : ) .ne.              &
    coarray( lbr(1) : ubr(1) , lbr(2) , lbr(3) , : )                   &
      [ imgpos(1) , imgpos(2)+1 , imgpos(3)+1 ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

! op 14
if ( imgpos(2) .ne. ucob(2) .and. imgpos(3) .ne. lcob(3) ) then
  if ( any(                                                            &
    coarray( lbr(1) : ubr(1) , ubv(2) , lbv(3) , : ) .ne.              &
    coarray( lbr(1) : ubr(1) , lbr(2) , ubr(3) , : )                   &
      [ imgpos(1) , imgpos(2)+1 , imgpos(3)-1 ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

! exchange 1D halos parallel to direction 2

! op 15
if ( imgpos(1) .ne. lcob(1) .and. imgpos(3) .ne. lcob(3) ) then
  if ( any(                                                            &
    coarray( lbv(1) , lbr(2) : ubr(2) , lbv(3) , : ) .ne.              &
    coarray( ubr(1) , lbr(2) : ubr(2) , ubr(3) , : )                   &
      [ imgpos(1)-1 , imgpos(2) , imgpos(3)-1 ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

! op 16
if ( imgpos(1) .ne. ucob(1) .and. imgpos(3) .ne. lcob(3) ) then
  if ( any(                                                            &
    coarray( ubv(1) , lbr(2) : ubr(2) , lbv(3) , : ) .ne.              &
    coarray( lbr(1) , lbr(2) : ubr(2) , ubr(3) , : )                   &
      [ imgpos(1)+1 , imgpos(2) , imgpos(3)-1 ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

! op 17
if ( imgpos(1) .ne. ucob(1) .and. imgpos(3) .ne. ucob(3) ) then
  if ( any(                                                            &
    coarray( ubv(1) , lbr(2) : ubr(2) , ubv(3) , : ) .ne.              &
    coarray( lbr(1) , lbr(2) : ubr(2) , lbr(3) , : )                   &
      [ imgpos(1)+1 , imgpos(2) , imgpos(3)+1 ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

! op 18
if ( imgpos(1) .ne. lcob(1) .and. imgpos(3) .ne. ucob(3) ) then
  if ( any(                                                            &
    coarray( lbv(1) , lbr(2) : ubr(2) , ubv(3) , : ) .ne.              & 
    coarray( ubr(1) , lbr(2) : ubr(2) , lbr(3) , : )                   &
      [ imgpos(1)-1 , imgpos(2) , imgpos(3)+1 ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

! Exchange 8 scalar halos
! See diagram cgca1 in the manual.

! op 19
if ( (imgpos(1) .ne. lcob(1)) .and. (imgpos(2) .ne. lcob(2)) .and.     &
     (imgpos(3) .ne. lcob(3)) ) then
  if ( any(                                                            &
    coarray( lbv(1) , lbv(2) , lbv(3) , : ) .ne.                       &
    coarray( ubr(1) , ubr(2) , ubr(3) , : )                            &
      [ imgpos(1)-1 , imgpos(2)-1 , imgpos(3)-1 ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

! op 20
if ( (imgpos(1) .ne. ucob(1)) .and. (imgpos(2) .ne. lcob(2)) .and.     &
     (imgpos(3) .ne. lcob(3)) ) then
  if ( any(                                                            &
    coarray( ubv(1) , lbv(2) , lbv(3) , : ) .ne.                       &
    coarray( lbr(1) , ubr(2) , ubr(3) , : )                            &
      [ imgpos(1)+1 , imgpos(2)-1 , imgpos(3)-1 ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

! op 21
if ( (imgpos(1) .ne. ucob(1)) .and. (imgpos(2) .ne. ucob(2)) .and.     &
     (imgpos(3) .ne. lcob(3)) ) then
  if ( any(                                                            &
    coarray( ubv(1) , ubv(2) , lbv(3) , : ) .ne.                       &
    coarray( lbr(1) , lbr(2) , ubr(3) , : )                            &
      [ imgpos(1)+1 , imgpos(2)+1 , imgpos(3)-1 ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

! op 22
if ( (imgpos(1) .ne. lcob(1)) .and. (imgpos(2) .ne. ucob(2)) .and.     &
     (imgpos(3) .ne. lcob(3)) ) then
  if ( any(                                                            &
    coarray( lbv(1) , ubv(2) , lbv(3) , : ) .ne.                       &
    coarray( ubr(1) , lbr(2) , ubr(3) , : )                            &
      [ imgpos(1)-1 , imgpos(2)+1 , imgpos(3)-1 ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

! op 23
if ( (imgpos(1) .ne. lcob(1)) .and. (imgpos(2) .ne. lcob(2)) .and.     &
     (imgpos(3) .ne. ucob(3)) ) then
  if ( any(                                                            &
    coarray( lbv(1) , lbv(2) , ubv(3) , : ) .ne.                       &
    coarray( ubr(1) , ubr(2) , lbr(3) , : )                            &
      [ imgpos(1)-1 , imgpos(2)-1 , imgpos(3)+1 ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

! op 24
if ( (imgpos(1) .ne. ucob(1)) .and. (imgpos(2) .ne. lcob(2)) .and.     &
     (imgpos(3) .ne. ucob(3)) ) then
  if ( any(                                                            &
    coarray( ubv(1) , lbv(2) , ubv(3) , : ) .ne.                       & 
    coarray( lbr(1) , ubr(2) , lbr(3) , : )                            &
      [ imgpos(1)+1 , imgpos(2)-1 , imgpos(3)+1 ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

! op 25
if ( (imgpos(1) .ne. ucob(1)) .and. (imgpos(2) .ne. ucob(2)) .and.     &
     (imgpos(3) .ne. ucob(3)) ) then
  if ( any(                                                            &
    coarray( ubv(1) , ubv(2) , ubv(3) , : ) .ne.                       &
    coarray( lbr(1) , lbr(2) , lbr(3) , : )                            &
      [ imgpos(1)+1 , imgpos(2)+1 , imgpos(3)+1 ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

! op 26
if ( (imgpos(1) .ne. lcob(1)) .and. (imgpos(2) .ne. ucob(2)) .and.     &
     (imgpos(3) .ne. ucob(3)) ) then
  if ( any(                                                            &
    coarray( lbv(1) , ubv(2) , ubv(3) , : ) .ne.                       &
    coarray( ubr(1) , lbr(2) , lbr(3) , : )                            &
      [ imgpos(1)-1 , imgpos(2)+1 , imgpos(3)+1 ] ) ) then
    flag = 1
    ! And return immediately
    return
  end if
end if

end procedure cgca_hxic

!*roboend*

end submodule m2hx_hxic
