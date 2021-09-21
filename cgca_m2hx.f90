!*robodoc*m* CGPACK/cgca_m2hx
!  NAME
!    cgca_m2hx
!  SYNOPSIS

!$Id: cgca_m2hx.f90 431 2017-06-30 13:13:49Z mexas $

module cgca_m2hx

!  DESCRIPTION
!    Module dealing with halo exchange
!  AUTHOR 
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    cgca_hxi, cgca_hxg
!  USES
!    cgca_m1co
!  USED BY
!    cgca_m3clvg cgca_m3sld
!  SOURCE

use cgca_m1co
implicit none

private
public :: cgca_hxi,  &
          cgca_hxir, & ! In submodule m2hx_hxir
          cgca_hxic, & ! In submodule m2hx_hxic
          cgca_hxg

interface
  ! In submodule m2hx_hxir
  module subroutine cgca_hxir( coarray )
    ! State of coarray will change!
    integer( kind=iarr ), allocatable, intent( inout ) ::              &
       coarray(:,:,:,:)[:,:,:]
  end subroutine cgca_hxir

  ! In submodule m2hx_hxic
  module subroutine cgca_hxic( coarray, flag )
    ! State of coarray does not change! 
    integer( kind=iarr ), allocatable, intent( inout ) ::              &
       coarray(:,:,:,:)[:,:,:]
    ! flag is 0 if and only if hx was done correctly, i.e. if
    ! the state of all real cells on the boundaries is consistent with
    ! the states of the corresponding halo cells. 
    integer, intent( out ) :: flag
  end subroutine cgca_hxic
end interface

contains

!*roboend*


!*robodoc*s* cgca_m2hx/cgca_hxi
!  NAME
!    cgca_hxi
!  SYNOPSIS

subroutine cgca_hxi( coarray )

!  INPUT
 
integer( kind=iarr ), allocatable, intent( inout ) ::                  &
 coarray(:,:,:,:)[:,:,:]

!  SIDE EFFECTS
!    coarray is changed
!  DESCRIPTION
!    This routine does internal halo exchange.
!    The routine exchanges halos on *all* cell state types.
!    This is an overkill, as it is likely that only one cell
!    state type needs to be halo exchanged at a time.
!    However, it makes for an easier code, and there is virtually
!    no performance penalty, so we do it this way.
!  NOTES
!    All images must call this routine!
!  USES
!    parameters from cgca_m1co
!  USED BY
!    module cgca_m3clvg: cgca_clgvp
!  SOURCE

integer ::                                                             &
  lbv(4)      , & ! lower bounds of the "virtual" coarray
  ubv(4)      , & ! upper bounds of the "virtual" coarray
  lbr(4)      , & ! lower bounds of the "real" coarray, lbv+1
  ubr(4)      , & ! upper bounds of the "real" coarray, ubv-1
  lcob(3)     , & ! lower cobounds of the coarray
  ucob(3)     , & ! upper cobounds of the coarray
  imgpos(3)       ! position of the image in a coarray grid

! check for allocated 
if ( .not. allocated( coarray ) ) then
  write (*,'(a)') "ERROR: cgca_hxi/cgca_m2hx: coarray is not allocated"
  error stop
end if

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

! exchange 2D halos in direction 1

! op 1
if ( imgpos(1) .ne. lcob(1) )                                          &
  coarray( lbv(1) , lbr(2) : ubr(2) , lbr(3) : ubr(3) , : ) =          &
  coarray( ubr(1) , lbr(2) : ubr(2) , lbr(3) : ubr(3) , : )            &
    [ imgpos(1)-1 , imgpos(2) , imgpos(3) ]

! op 2
if ( imgpos(1) .ne. ucob(1) )                                          &
  coarray( ubv(1) , lbr(2) : ubr(2) , lbr(3) : ubr(3) , : ) =          &
  coarray( lbr(1) , lbr(2) : ubr(2) , lbr(3) : ubr(3) , : )            &
    [ imgpos(1)+1 , imgpos(2) , imgpos(3) ]

! exchange 2D halos in direction 2

! op 3
if ( imgpos(2) .ne. lcob(2) )                                          &
  coarray( lbr(1) : ubr(1) , lbv(2) , lbr(3) : ubr(3) , : ) =          &
  coarray( lbr(1) : ubr(1) , ubr(2) , lbr(3) : ubr(3) , : )            &
    [ imgpos(1) , imgpos(2)-1 , imgpos(3) ]

! op 4
if ( imgpos(2) .ne. ucob(2) )                                          &
  coarray( lbr(1) : ubr(1) , ubv(2) , lbr(3) : ubr(3) , : ) =          &
  coarray( lbr(1) : ubr(1) , lbr(2) , lbr(3) : ubr(3) , : )            &
    [ imgpos(1) , imgpos(2)+1 , imgpos(3) ]

! exchange 2D halos in direction 3

! op 5
if ( imgpos(3) .ne. lcob(3) )                                          &
  coarray( lbr(1) : ubr(1) , lbr(2) : ubr(2) , lbv(3) , : ) =          &
  coarray( lbr(1) : ubr(1) , lbr(2) : ubr(2) , ubr(3) , : )            &
    [ imgpos(1) , imgpos(2) , imgpos(3)-1 ]

! op 6
if ( imgpos(3) .ne. ucob(3) )                                          &
  coarray( lbr(1) : ubr(1) , lbr(2) : ubr(2) , ubv(3) , : ) =          &
  coarray( lbr(1) : ubr(1) , lbr(2) : ubr(2) , lbr(3) , : )            &
    [ imgpos(1) , imgpos(2) , imgpos(3)+1 ]

! exchange 1D halos parallel to direction 3

! op 7
if ( imgpos(1) .ne. lcob(1) .and. imgpos(2) .ne. lcob(2) )             &
  coarray( lbv(1) , lbv(2) , lbr(3) : ubr(3) , : ) =                   &
  coarray( ubr(1) , ubr(2) , lbr(3) : ubr(3) , : )                     &
    [ imgpos(1)-1 , imgpos(2)-1 , imgpos(3) ]

! op 8
if ( imgpos(1) .ne. ucob(1) .and. imgpos(2) .ne. lcob(2) )             &
  coarray( ubv(1) , lbv(2) , lbr(3) : ubr(3) , :) =                    &
  coarray( lbr(1) , ubr(2) , lbr(3) : ubr(3) , :)                      &
    [ imgpos(1)+1 , imgpos(2)-1 , imgpos(3) ]

! op 9
if ( imgpos(1) .ne. ucob(1) .and. imgpos(2) .ne. ucob(2) )             &
  coarray( ubv(1) , ubv(2) , lbr(3) : ubr(3) , : ) =                   &
  coarray( lbr(1) , lbr(2) , lbr(3) : ubr(3) , : )                     &
    [ imgpos(1)+1 , imgpos(2)+1 , imgpos(3) ]

! op 10
if ( imgpos(1) .ne. lcob(1) .and. imgpos(2) .ne. ucob(2) )             &
  coarray( lbv(1) , ubv(2) , lbr(3) : ubr(3) , : ) =                   &
  coarray( ubr(1) , lbr(2) , lbr(3) : ubr(3) , : )                     &
    [ imgpos(1)-1 , imgpos(2)+1 , imgpos(3) ]

! exchange 1D halos parallel to direction 1

! op 11
if ( imgpos(2) .ne. lcob(2) .and. imgpos(3) .ne. lcob(3) )             &
  coarray( lbr(1) : ubr(1) , lbv(2) , lbv(3) , : ) =                   &
  coarray( lbr(1) : ubr(1) , ubr(2) , ubr(3) , : )                     &
    [ imgpos(1) , imgpos(2)-1 , imgpos(3)-1 ]

! op 12
if ( imgpos(2) .ne. lcob(2) .and. imgpos(3) .ne. ucob(3) )             &
  coarray( lbr(1) : ubr(1) , lbv(2) , ubv(3) , : ) =                   &
  coarray( lbr(1) : ubr(1) , ubr(2) , lbr(3) , : )                     &
    [ imgpos(1) , imgpos(2)-1 , imgpos(3)+1 ]

! op 13
if ( imgpos(2) .ne. ucob(2) .and. imgpos(3) .ne. ucob(3) )             &
  coarray( lbr(1) : ubr(1) , ubv(2) , ubv(3) , : ) =                   &
  coarray( lbr(1) : ubr(1) , lbr(2) , lbr(3) , : )                     &
    [ imgpos(1) , imgpos(2)+1 , imgpos(3)+1 ]

! op 14
if ( imgpos(2) .ne. ucob(2) .and. imgpos(3) .ne. lcob(3) )             &
  coarray( lbr(1) : ubr(1) , ubv(2) , lbv(3) , : ) =                   &
  coarray( lbr(1) : ubr(1) , lbr(2) , ubr(3) , : )                     &
    [ imgpos(1) , imgpos(2)+1 , imgpos(3)-1 ]

! exchange 1D halos parallel to direction 2

! op 15
if ( imgpos(1) .ne. lcob(1) .and. imgpos(3) .ne. lcob(3) )             &
  coarray( lbv(1) , lbr(2) : ubr(2) , lbv(3) , : ) =                   &
  coarray( ubr(1) , lbr(2) : ubr(2) , ubr(3) , : )                     &
    [ imgpos(1)-1 , imgpos(2) , imgpos(3)-1 ]

! op 16
if ( imgpos(1) .ne. ucob(1) .and. imgpos(3) .ne. lcob(3) )             &
  coarray( ubv(1) , lbr(2) : ubr(2) , lbv(3) , : ) =                   &
  coarray( lbr(1) , lbr(2) : ubr(2) , ubr(3) , : )                     &
    [ imgpos(1)+1 , imgpos(2) , imgpos(3)-1 ]

! op 17
if ( imgpos(1) .ne. ucob(1) .and. imgpos(3) .ne. ucob(3) )             &
  coarray( ubv(1) , lbr(2) : ubr(2) , ubv(3) , : ) =                   &
  coarray( lbr(1) , lbr(2) : ubr(2) , lbr(3) , : )                     &
    [ imgpos(1)+1 , imgpos(2) , imgpos(3)+1 ]

! op 18
if ( imgpos(1) .ne. lcob(1) .and. imgpos(3) .ne. ucob(3) )             &
  coarray( lbv(1) , lbr(2) : ubr(2) , ubv(3) , : ) =                   &
  coarray( ubr(1) , lbr(2) : ubr(2) , lbr(3) , : )                     &
    [ imgpos(1)-1 , imgpos(2) , imgpos(3)+1 ]

! Exchange 8 scalar halos
! See diagram cgca1 in the manual.

! op 19
if ( (imgpos(1) .ne. lcob(1)) .and. (imgpos(2) .ne. lcob(2)) .and.     &
     (imgpos(3) .ne. lcob(3)) )                                        &
  coarray( lbv(1) , lbv(2) , lbv(3) , : ) =                            &
  coarray( ubr(1) , ubr(2) , ubr(3) , : )                              &
                   [ imgpos(1)-1 , imgpos(2)-1 , imgpos(3)-1 ]

! op 20
if ( (imgpos(1) .ne. ucob(1)) .and. (imgpos(2) .ne. lcob(2)) .and.     &
     (imgpos(3) .ne. lcob(3)) )                                        &
  coarray( ubv(1) , lbv(2) , lbv(3) , : ) =                            &
  coarray( lbr(1) , ubr(2) , ubr(3) , : )                              &
                   [ imgpos(1)+1 , imgpos(2)-1 , imgpos(3)-1 ]

! op 21
if ( (imgpos(1) .ne. ucob(1)) .and. (imgpos(2) .ne. ucob(2)) .and.     &
     (imgpos(3) .ne. lcob(3)) )                                        &
  coarray( ubv(1) , ubv(2) , lbv(3) , : ) =                            &
  coarray( lbr(1) , lbr(2) , ubr(3) , : )                              &
                   [ imgpos(1)+1 , imgpos(2)+1 , imgpos(3)-1 ]

! op 22
if ( (imgpos(1) .ne. lcob(1)) .and. (imgpos(2) .ne. ucob(2)) .and.     &
     (imgpos(3) .ne. lcob(3)) )                                        &
  coarray( lbv(1) , ubv(2) , lbv(3) , : ) =                            &
  coarray( ubr(1) , lbr(2) , ubr(3) , : )                              &
                   [ imgpos(1)-1 , imgpos(2)+1 , imgpos(3)-1 ]

! op 23
if ( (imgpos(1) .ne. lcob(1)) .and. (imgpos(2) .ne. lcob(2)) .and.     &
     (imgpos(3) .ne. ucob(3)) )                                        &
  coarray( lbv(1) , lbv(2) , ubv(3) , : ) =                            &
  coarray( ubr(1) , ubr(2) , lbr(3) , : )                              &
                   [ imgpos(1)-1 , imgpos(2)-1 , imgpos(3)+1 ]

! op 24
if ( (imgpos(1) .ne. ucob(1)) .and. (imgpos(2) .ne. lcob(2)) .and.     &
     (imgpos(3) .ne. ucob(3)) )                                        &
  coarray( ubv(1) , lbv(2) , ubv(3) , : ) =                            &
  coarray( lbr(1) , ubr(2) , lbr(3) , : )                              &
                   [ imgpos(1)+1 , imgpos(2)-1 , imgpos(3)+1 ]

! op 25
if ( (imgpos(1) .ne. ucob(1)) .and. (imgpos(2) .ne. ucob(2)) .and.     &
     (imgpos(3) .ne. ucob(3)) )                                        &
  coarray( ubv(1) , ubv(2) , ubv(3) , : ) =                            &
  coarray( lbr(1) , lbr(2) , lbr(3) , : )                              &
                   [ imgpos(1)+1 , imgpos(2)+1 , imgpos(3)+1 ]

! op 26
if ( (imgpos(1) .ne. lcob(1)) .and. (imgpos(2) .ne. ucob(2)) .and.     &
     (imgpos(3) .ne. ucob(3)) )                                        &
  coarray( lbv(1) , ubv(2) , ubv(3) , : ) =                            &
  coarray( ubr(1) , lbr(2) , lbr(3) , : )                              &
                   [ imgpos(1)-1 , imgpos(2)+1 , imgpos(3)+1 ]

end subroutine cgca_hxi

!*roboend*


!*robodoc*s* cgca_m2hx/cgca_hxg
!  NAME
!    cgca_hxg
!  SYNOPSIS

subroutine cgca_hxg(coarray)

!  INPUT

integer(kind=iarr),allocatable,intent(inout) :: coarray(:,:,:,:)[:,:,:]

!  SIDE EFFECTS
!    coarray is changed
!  DESCRIPTION
!    This routine does the global halo exchange, i.e. on the
!    boundaries of the whole model.
!    In other wordsl it imposes self-similar boundary conditions on
!    the global (super) array.
!    The routine exchanges halos on *all* cell state types.
!    This is an overkill, as it is likely that only one cell
!    state type needs to be halo exchanged at a time.
!    However, it makes for an easier code, and there is
!    virtually no performance penaly, so we do it this way.
!  NOTES
!    All images must call this routine!
!  USES
!    none
!  USED BY
!    module cgca_m3clvg: cgca_clvgp
!  SOURCE

integer(kind=idef) :: &
  lbv(4)      , & ! lower bounds of the "virtual" coarray
  ubv(4)      , & ! upper bounds of the "virtual" coarray
  lbr(4)      , & ! lower bounds of the "real" coarray, lbv+1
  ubr(4)      , & ! upper bounds of the "real" coarray, ubv-1
  lcob(3)     , & ! lower cobounds of the coarray
  ucob(3)     , & ! upper cobounds of the coarray
  imgpos(3)   , & ! position of the image in a coarray grid
  imgpos1mns1 , & ! positions of the neighbouring images
  imgpos1pls1 , & ! along 3 directions
  imgpos2mns1 , & !
  imgpos2pls1 , & !
  imgpos3mns1 , & !
  imgpos3pls1 , & !
  z1,z2,z3        ! temp vars to store coarray grid coordinates 

! check for allocated 
if (.not. allocated(coarray)) &
  error stop "ERROR: cgca_hlg: coarray is not allocated"

lbv=lbound(coarray)
ubv=ubound(coarray)
lbr=lbv+1
ubr=ubv-1
lcob=lcobound(coarray)
ucob=ucobound(coarray)

imgpos=this_image(coarray)
imgpos1mns1=imgpos(1)-1
imgpos1pls1=imgpos(1)+1
imgpos2mns1=imgpos(2)-1
imgpos2pls1=imgpos(2)+1
imgpos3mns1=imgpos(3)-1
imgpos3pls1=imgpos(3)+1

! Make sure only the virtual (halo) arrays are assigned to.
! The real array values must never appear on the left
! hand side of the assignment expressions.
! The halo exchange process is copying real array values
! into halos. There must not ever be copying of real values
! to real, or halo to halo, or halo to real.

!**************************************
! 2D halos
!**************************************

! exchange 2D halos in direction 1

if (imgpos(1) .eq. lcob(1))                          &
  coarray(lbv(1),lbr(2):ubr(2),lbr(3):ubr(3),:) =    &
  coarray(ubr(1),lbr(2):ubr(2),lbr(3):ubr(3),:)      &
    [ucob(1),imgpos(2),imgpos(3)]

if (imgpos(1) .eq. ucob(1))                          &
  coarray(ubv(1),lbr(2):ubr(2),lbr(3):ubr(3),:) =    &
  coarray(lbr(1),lbr(2):ubr(2),lbr(3):ubr(3),:)      &
    [lcob(1),imgpos(2),imgpos(3)]

! exchange 2D halos in direction 2

if (imgpos(2) .eq. lcob(2))                          &
  coarray(lbr(1):ubr(1),lbv(2),lbr(3):ubr(3),:) =    &
  coarray(lbr(1):ubr(1),ubr(2),lbr(3):ubr(3),:)      &
    [imgpos(1),ucob(2),imgpos(3)]

if (imgpos(2) .eq. ucob(2))                          &
  coarray(lbr(1):ubr(1),ubv(2),lbr(3):ubr(3),:) =    &
  coarray(lbr(1):ubr(1),lbr(2),lbr(3):ubr(3),:)      &
    [imgpos(1),lcob(2),imgpos(3)]

! exchange 2D halos in direction 3

if (imgpos(3) .eq. lcob(3))                          &
  coarray(lbr(1):ubr(1),lbr(2):ubr(2),lbv(3),:) =    &
  coarray(lbr(1):ubr(1),lbr(2):ubr(2),ubr(3),:)      &
    [imgpos(1),imgpos(2),ucob(3)]

if (imgpos(3) .eq. ucob(3))                          &
  coarray(lbr(1):ubr(1),lbr(2):ubr(2),ubv(3),:) =    &
  coarray(lbr(1):ubr(1),lbr(2):ubr(2),lbr(3),:)      &
    [imgpos(1),imgpos(2),lcob(3)]

!**************************************
! 1D halos
!**************************************

! exchange 1D halos parallel to direction 3
! the 4 edges of the super array

! operation 1
if (imgpos(1) .eq. lcob(1) .and. imgpos(2) .eq. lcob(2))               &
  coarray(lbv(1),lbv(2),lbr(3):ubr(3),:) =                             &
  coarray(ubr(1),ubr(2),lbr(3):ubr(3),:) [ucob(1),ucob(2),imgpos(3)]

! operation 2
if (imgpos(1) .eq. ucob(1) .and. imgpos(2) .eq. lcob(2))               &
  coarray(ubv(1),lbv(2),lbr(3):ubr(3),:) =                             &
  coarray(lbr(1),ubr(2),lbr(3):ubr(3),:) [lcob(1),ucob(2),imgpos(3)]

! operation 3
if (imgpos(1) .eq. ucob(1) .and. imgpos(2) .eq. ucob(2))               &
  coarray(ubv(1),ubv(2),lbr(3):ubr(3),:) =                             &
  coarray(lbr(1),lbr(2),lbr(3):ubr(3),:) [lcob(1),lcob(2),imgpos(3)]

! operation 4
if (imgpos(1) .eq. lcob(1) .and. imgpos(2) .eq. ucob(2))               &
  coarray(lbv(1),ubv(2),lbr(3):ubr(3),:) =                             &
  coarray(ubr(1),lbr(2),lbr(3):ubr(3),:) [ucob(1),lcob(2),imgpos(3)]

! exchange 1D halos parallel to direction 3
! all intermediate edges, self-similarity along 1

! operation 5
if (imgpos(1) .eq. lcob(1) .and. imgpos(2) .ne. ucob(2))               &
  coarray(lbv(1),ubv(2),lbr(3):ubr(3),:) =                             &
  coarray(ubr(1),lbr(2),lbr(3):ubr(3),:) [ucob(1),imgpos2pls1,imgpos(3)]

! operation 6
if (imgpos(1) .eq. lcob(1) .and. imgpos(2) .ne. lcob(2))               &
  coarray(lbv(1),lbv(2),lbr(3):ubr(3),:) =                             &
  coarray(ubr(1),ubr(2),lbr(3):ubr(3),:) [ucob(1),imgpos2mns1,imgpos(3)]

! operation 7
if (imgpos(1) .eq. ucob(1) .and. imgpos(2) .ne. ucob(2))               &
  coarray(ubv(1),ubv(2),lbr(3):ubr(3),:) =                             &
  coarray(lbr(1),lbr(2),lbr(3):ubr(3),:) [lcob(1),imgpos2pls1,imgpos(3)]

! operation 8
if (imgpos(1) .eq. ucob(1) .and. imgpos(2) .ne. lcob(2))               &
  coarray(ubv(1),lbv(2),lbr(3):ubr(3),:) =                             &
  coarray(lbr(1),ubr(2),lbr(3):ubr(3),:) [lcob(1),imgpos2mns1,imgpos(3)]

! exchange 1D halos parallel to direction 3
! all intermediate edges, self-similarity along 2

! operation 9
if (imgpos(1) .ne. ucob(1) .and. imgpos(2) .eq. lcob(2))               &
  coarray(ubv(1),lbv(2),lbr(3):ubr(3),:) =                             &
  coarray(lbr(1),ubr(2),lbr(3):ubr(3),:) [imgpos1pls1,ucob(2),imgpos(3)]

! operation 10
if (imgpos(1) .ne. lcob(1) .and. imgpos(2) .eq. lcob(2))               &
  coarray(lbv(1),lbv(2),lbr(3):ubr(3),:) =                             &
  coarray(ubr(1),ubr(2),lbr(3):ubr(3),:) [imgpos1mns1,ucob(2),imgpos(3)]

! operation 11
if (imgpos(1) .ne. ucob(1) .and. imgpos(2) .eq. ucob(2))               &
  coarray(ubv(1),ubv(2),lbr(3):ubr(3),:) =                             &
  coarray(lbr(1),lbr(2),lbr(3):ubr(3),:) [imgpos1pls1,lcob(2),imgpos(3)]

! operation 12
if (imgpos(1) .ne. lcob(1) .and. imgpos(2) .eq. ucob(2))               &
  coarray(lbv(1),ubv(2),lbr(3):ubr(3),:) =                             &
  coarray(ubr(1),lbr(2),lbr(3):ubr(3),:) [imgpos1mns1,lcob(2),imgpos(3)]

! exchange 1D halos parallel to direction 1
! the 4 edges of the super array

! operation 1
if (imgpos(2) .eq. lcob(2) .and. imgpos(3) .eq. lcob(3))               &
  coarray(lbr(1):ubr(1),lbv(2),lbv(3),:) =                             &
  coarray(lbr(1):ubr(1),ubr(2),ubr(3),:) [imgpos(1),ucob(2),ucob(3)]

! operation 2
if (imgpos(2) .eq. ucob(2) .and. imgpos(3) .eq. lcob(3))               &
  coarray(lbr(1):ubr(1),ubv(2),lbv(3),:) =                             &
  coarray(lbr(1):ubr(1),lbr(2),ubr(3),:) [imgpos(1),lcob(2),ucob(3)]

! operation 3 
if (imgpos(2) .eq. ucob(2) .and. imgpos(3) .eq. ucob(3))               &
  coarray(lbr(1):ubr(1),ubv(2),ubv(3),:) =                             &
  coarray(lbr(1):ubr(1),lbr(2),lbr(3),:) [imgpos(1),lcob(2),lcob(3)]

! operation 4
if (imgpos(2) .eq. lcob(2) .and. imgpos(3) .eq. ucob(3))               &
  coarray(lbr(1):ubr(1),lbv(2),ubv(3),:) =                             &
  coarray(lbr(1):ubr(1),ubr(2),lbr(3),:) [imgpos(1),ucob(2),lcob(3)]

! exchange 1D halos parallel to direction 1
! intermediate edges, self-similarity along 2

! operation 5
if (imgpos(2) .eq. lcob(2) .and. imgpos(3) .ne. ucob(3))               &
  coarray(lbr(1):ubr(1),lbv(2),ubv(3),:) =                             &
  coarray(lbr(1):ubr(1),ubr(2),lbr(3),:) [imgpos(1),ucob(2),imgpos3pls1]

! operation 6
if (imgpos(2) .eq. lcob(2) .and. imgpos(3) .ne. lcob(3))               &
  coarray(lbr(1):ubr(1),lbv(2),lbv(3),:) =                             &
  coarray(lbr(1):ubr(1),ubr(2),ubr(3),:) [imgpos(1),ucob(2),imgpos3mns1]

! operation 7
if (imgpos(2) .eq. ucob(2) .and. imgpos(3) .ne. ucob(3))               &
  coarray(lbr(1):ubr(1),ubv(2),ubv(3),:) =                             &
  coarray(lbr(1):ubr(1),lbr(2),lbr(3),:) [imgpos(1),lcob(2),imgpos3pls1]

! operation 8
if (imgpos(2) .eq. ucob(2) .and. imgpos(3) .ne. lcob(3))               &
  coarray(lbr(1):ubr(1),ubv(2),lbv(3),:) =                             &
  coarray(lbr(1):ubr(1),lbr(2),ubr(3),:) [imgpos(1),lcob(2),imgpos3mns1]

! exchange 1D halos parallel to direction 1
! intermediate edges, self-similarity along 3

! operation 9
if (imgpos(2) .ne. ucob(2) .and. imgpos(3) .eq. lcob(3))               &
  coarray(lbr(1):ubr(1),ubv(2),lbv(3),:) =                             &
  coarray(lbr(1):ubr(1),lbr(2),ubr(3),:) [imgpos(1),imgpos2pls1,ucob(3)]

! operation 10
if (imgpos(2) .ne. lcob(2) .and. imgpos(3) .eq. lcob(3))               &
  coarray(lbr(1):ubr(1),lbv(2),lbv(3),:) =                             &
  coarray(lbr(1):ubr(1),ubr(2),ubr(3),:) [imgpos(1),imgpos2mns1,ucob(3)]

! operation 11
if (imgpos(2) .ne. ucob(2) .and. imgpos(3) .eq. ucob(3))               &
  coarray(lbr(1):ubr(1),ubv(2),ubv(3),:) =                             &
  coarray(lbr(1):ubr(1),lbr(2),lbr(3),:) [imgpos(1),imgpos2pls1,lcob(3)]

! operation 12
if (imgpos(2) .ne. lcob(2) .and. imgpos(3) .eq. ucob(3))               &
  coarray(lbr(1):ubr(1),lbv(2),ubv(3),:) =                             &
  coarray(lbr(1):ubr(1),ubr(2),lbr(3),:) [imgpos(1),imgpos2mns1,lcob(3)]

! exchange 1D halos parallel to direction 2
! the 4 edges of the super array

! operation 1
if (imgpos(1) .eq. lcob(1) .and. imgpos(3) .eq. lcob(3))               &
  coarray(lbv(1),lbr(2):ubr(2),lbv(3),:) =                             &
  coarray(ubr(1),lbr(2):ubr(2),ubr(3),:) [ucob(1),imgpos(2),ucob(3)]

! operation 2
if (imgpos(1) .eq. lcob(1) .and. imgpos(3) .eq. ucob(3))               &
  coarray(lbv(1),lbr(2):ubr(2),ubv(3),:) =                             &
  coarray(ubr(1),lbr(2):ubr(2),lbr(3),:) [ucob(1),imgpos(2),lcob(3)]

! operation 3
if (imgpos(1) .eq. ucob(1) .and. imgpos(3) .eq. ucob(3))               &
  coarray(ubv(1),lbr(2):ubr(2),ubv(3),:) =                             &
  coarray(lbr(1),lbr(2):ubr(2),lbr(3),:) [lcob(1),imgpos(2),lcob(3)]

! operation 4
if (imgpos(1) .eq. ucob(1) .and. imgpos(3) .eq. lcob(3))               &
  coarray(ubv(1),lbr(2):ubr(2),lbv(3),:) =                             &
  coarray(lbr(1),lbr(2):ubr(2),ubr(3),:) [lcob(1),imgpos(2),ucob(3)]

! exchange 1D halos parallel to direction 2
! intermediate edges, self-similarity along 3

! operation 5
if (imgpos(1) .ne. ucob(1) .and. imgpos(3) .eq. lcob(3))               &
  coarray(ubv(1),lbr(2):ubr(2),lbv(3),:) =                             &
  coarray(lbr(1),lbr(2):ubr(2),ubr(3),:) [imgpos1pls1,imgpos(2),ucob(3)]

! operation 6
if (imgpos(1) .ne. lcob(1) .and. imgpos(3) .eq. lcob(3))               &
  coarray(lbv(1),lbr(2):ubr(2),lbv(3),:) =                             &
  coarray(ubr(1),lbr(2):ubr(2),ubr(3),:) [imgpos1mns1,imgpos(2),ucob(3)]

! operation 7
if (imgpos(1) .ne. ucob(1) .and. imgpos(3) .eq. ucob(3))               &
  coarray(ubv(1),lbr(2):ubr(2),ubv(3),:) =                             &
  coarray(lbr(1),lbr(2):ubr(2),lbr(3),:) [imgpos1pls1,imgpos(2),lcob(3)]

! operation 8
if (imgpos(1) .ne. lcob(1) .and. imgpos(3) .eq. ucob(3))               &
  coarray(lbv(1),lbr(2):ubr(2),ubv(3),:) =                             &
  coarray(ubr(1),lbr(2):ubr(2),lbr(3),:) [imgpos1mns1,imgpos(2),lcob(3)]

! exchange 1D halos parallel to direction 2
! intermediate edges, self-similarity along 1

! operation 9
if (imgpos(1) .eq. lcob(1) .and. imgpos(3) .ne. ucob(3))               &
  coarray(lbv(1),lbr(2):ubr(2),ubv(3),:) =                             &
  coarray(ubr(1),lbr(2):ubr(2),lbr(3),:) [ucob(1),imgpos(2),imgpos3pls1]

! operation 10
if (imgpos(1) .eq. lcob(1) .and. imgpos(3) .ne. lcob(3))               &
  coarray(lbv(1),lbr(2):ubr(2),lbv(3),:) =                             &
  coarray(ubr(1),lbr(2):ubr(2),ubr(3),:) [ucob(1),imgpos(2),imgpos3mns1]

! operation 11
if (imgpos(1) .eq. ucob(1) .and. imgpos(3) .ne. ucob(3))               &
  coarray(ubv(1),lbr(2):ubr(2),ubv(3),:) =                             &
  coarray(lbr(1),lbr(2):ubr(2),lbr(3),:) [lcob(1),imgpos(2),imgpos3pls1]

! operation 12
if (imgpos(1) .eq. ucob(1) .and. imgpos(3) .ne. lcob(3))               &
  coarray(ubv(1),lbr(2):ubr(2),lbv(3),:) =                             &
  coarray(lbr(1),lbr(2):ubr(2),ubr(3),:) [lcob(1),imgpos(2),imgpos3mns1]

!**************************************
! corner halos
!**************************************

! corner 1
if (imgpos(1) .eq. lcob(1) .or.             &
    imgpos(2) .eq. lcob(2) .or.             &
    imgpos(3) .eq. lcob(3))                 &
then
  z1 = imgpos1mns1
  z2 = imgpos2mns1
  z3 = imgpos3mns1
  if (z1 .lt. lcob(1)) z1 = ucob(1)
  if (z2 .lt. lcob(2)) z2 = ucob(2)
  if (z3 .lt. lcob(3)) z3 = ucob(3)
  coarray(lbv(1),lbv(2),lbv(3),:) =         &
  coarray(ubr(1),ubr(2),ubr(3),:) [z1,z2,z3]
end if

! corner 2
if (imgpos(1) .eq. ucob(1) .or.             &
    imgpos(2) .eq. lcob(2) .or.             &
    imgpos(3) .eq. lcob(3))                 &
then
  z1 = imgpos1pls1
  z2 = imgpos2mns1
  z3 = imgpos3mns1
  if (z1 .gt. ucob(1)) z1 = lcob(1)
  if (z2 .lt. lcob(2)) z2 = ucob(2)
  if (z3 .lt. lcob(3)) z3 = ucob(3)
  coarray(ubv(1),lbv(2),lbv(3),:) =         &
  coarray(lbr(1),ubr(2),ubr(3),:) [z1,z2,z3]
end if

! corner 3
if (imgpos(1) .eq. lcob(1) .or.             &
    imgpos(2) .eq. ucob(2) .or.             &
    imgpos(3) .eq. lcob(3))                 &
then
  z1 = imgpos1mns1
  z2 = imgpos2pls1
  z3 = imgpos3mns1
  if (z1 .lt. lcob(1)) z1 = ucob(1)
  if (z2 .gt. ucob(2)) z2 = lcob(2)
  if (z3 .lt. lcob(3)) z3 = ucob(3)
  coarray(lbv(1),ubv(2),lbv(3),:) =         &
  coarray(ubr(1),lbr(2),ubr(3),:) [z1,z2,z3]
end if

! corner 4
if (imgpos(1) .eq. ucob(1) .or.             &
    imgpos(2) .eq. ucob(2) .or.             &
    imgpos(3) .eq. lcob(3))                 &
then
  z1 = imgpos1pls1
  z2 = imgpos2pls1
  z3 = imgpos3mns1
  if (z1 .gt. ucob(1)) z1 = lcob(1)
  if (z2 .gt. ucob(2)) z2 = lcob(2)
  if (z3 .lt. lcob(3)) z3 = ucob(3)
  coarray(ubv(1),ubv(2),lbv(3),:) =         &
  coarray(lbr(1),lbr(2),ubr(3),:) [z1,z2,z3]
end if

! corner 5
if (imgpos(1) .eq. lcob(1) .or.             &
    imgpos(2) .eq. lcob(2) .or.             &
    imgpos(3) .eq. ucob(3))                 &
then
  z1 = imgpos1mns1
  z2 = imgpos2mns1
  z3 = imgpos3pls1
  if (z1 .lt. lcob(1)) z1 = ucob(1)
  if (z2 .lt. lcob(2)) z2 = ucob(2)
  if (z3 .gt. ucob(3)) z3 = lcob(3)
  coarray(lbv(1),lbv(2),ubv(3),:) =         &
  coarray(ubr(1),ubr(2),lbr(3),:) [z1,z2,z3]
end if

! corner 6
if (imgpos(1) .eq. ucob(1) .or.             &
    imgpos(2) .eq. lcob(2) .or.             &
    imgpos(3) .eq. ucob(3))                 &
then
  z1 = imgpos1pls1
  z2 = imgpos2mns1
  z3 = imgpos3pls1
  if (z1 .gt. ucob(1)) z1 = lcob(1)
  if (z2 .lt. lcob(2)) z2 = ucob(2)
  if (z3 .gt. ucob(3)) z3 = lcob(3)
  coarray(ubv(1),lbv(2),ubv(3),:) =         &
  coarray(lbr(1),ubr(2),lbr(3),:) [z1,z2,z3]
end if

! corner 7
if (imgpos(1) .eq. lcob(1) .or.             &
    imgpos(2) .eq. ucob(2) .or.             &
    imgpos(3) .eq. ucob(3))                 &
then
  z1 = imgpos1mns1
  z2 = imgpos2pls1
  z3 = imgpos3pls1
  if (z1 .lt. lcob(1)) z1 = ucob(1)
  if (z2 .gt. ucob(2)) z2 = lcob(2)
  if (z3 .gt. ucob(3)) z3 = lcob(3)
  coarray(lbv(1),ubv(2),ubv(3),:) =         &
  coarray(ubr(1),lbr(2),lbr(3),:) [z1,z2,z3]
end if

! corner 8
if (imgpos(1) .eq. ucob(1) .or.             &
    imgpos(2) .eq. ucob(2) .or.             &
    imgpos(3) .eq. ucob(3))                 &
then
  z1 = imgpos1pls1
  z2 = imgpos2pls1
  z3 = imgpos3pls1
  if (z1 .gt. ucob(1)) z1 = lcob(1)
  if (z2 .gt. ucob(2)) z2 = lcob(2)
  if (z3 .gt. ucob(3)) z3 = lcob(3)
  coarray(ubv(1),ubv(2),ubv(3),:) =         &
  coarray(lbr(1),lbr(2),lbr(3),:) [z1,z2,z3]
end if

end subroutine cgca_hxg

!*roboend*

end module cgca_m2hx
