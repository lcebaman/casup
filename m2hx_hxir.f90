!*robodoc*f* cgca_m2hx/m2hx_hxir
!  NAME
!    m2hx_hxir
!  SYNOPSIS

!$Id: m2hx_hxir.f90 423 2017-06-25 20:54:50Z mexas $

submodule ( cgca_m2hx ) m2hx_hxir

!  DESCRIPTION
!    Submodule of cgca_m2hx with an internal halo exchange routine
!    with random sequence of remote operations.
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENCE
!  CONTAINS
!    cgca_hxir
!  USES
!    Variables and parameters from the parent module cgca_m2hx.
!  USED BY
!    The parent module cgca_m2hx.
!  SOURCE

implicit none

contains

!*roboend*


!*robodoc*s* cgca_m2hx/cgca_hxir
!  NAME
!    cgca_hxir
!  SYNOPSIS

module procedure cgca_hxir

!  INPUT
!    See the parent module
!  OUTPUT
!    See INPUT
!  SIDE EFFECTS
!    coarray is changed
!  DESCRIPTION
!    This routine does internal halo exchange in random order.
!    The routine exchanges halos on *all* cell state types.
!    This is an overkill, as it is likely that only one cell
!    state type needs to be halo exchanged at a time.
!    However, it makes for an easier code, and there is virtually
!    no performance penalty, so we do it this way.
!  NOTES
!    All images must call this routine!
!  USES
!    All data objects from parent module cgca_m2hx by host association.
!  USED BY
!    module cgca_m2hx
!  SOURCE

! Number of groups of remote calls
! This parameter is used to randomoise the order of remote calls
integer, parameter :: ngroups=3

integer ::                                                             &
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
  istart      , & ! starting group of remote calls
  i, idx          ! loop indices

real :: rnd

! check for allocated 
if ( .not. allocated( coarray ) )                                      &
  error stop "ERROR: m2hx_hxir/cgca_hxir: coarray is not allocated"

 lbv = lbound( coarray )
 ubv = ubound( coarray )
 lbr = lbv + 1
 ubr = ubv - 1
lcob = lcobound( coarray )
ucob = ucobound( coarray )

     imgpos = this_image( coarray )
imgpos1mns1 = imgpos(1) - 1
imgpos1pls1 = imgpos(1) + 1
imgpos2mns1 = imgpos(2) - 1
imgpos2pls1 = imgpos(2) + 1
imgpos3mns1 = imgpos(3) - 1
imgpos3pls1 = imgpos(3) + 1

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

! I split all transfers into several groups.
! The starting group is chosen at random.
call random_number( rnd )         ! [ 0 .. 1 )
istart = int( rnd * ngroups ) + 1 ! [ 1 .. ngroups ]

!rcalls: do i = istart , istart + ngroups - 1
rcalls: do i = 1,1

  ! Remote call group index
  idx = i
  if ( idx .gt. ngroups ) idx = idx - ngroups 

  groups: if ( idx .eq. 1 ) then
    ! 1st group of remote calls

    ! exchange 2D halos in direction 1

    if ( imgpos(1) .ne. lcob(1) )                                      &
      coarray( lbv(1) , lbr(2) : ubr(2) , lbr(3) : ubr(3) , : ) =      &
      coarray( ubr(1) , lbr(2) : ubr(2) , lbr(3) : ubr(3) , : )        &
        [ imgpos1mns1 , imgpos(2) , imgpos(3) ]

    if ( imgpos(1) .ne. ucob(1) )                                      &
      coarray( ubv(1) , lbr(2) : ubr(2) , lbr(3) : ubr(3) , : ) =      &
      coarray( lbr(1) , lbr(2) : ubr(2) , lbr(3) : ubr(3) , : )        &
        [ imgpos1pls1 , imgpos(2) , imgpos(3) ]

    ! exchange 1D halos parallel to direction 1
    
    if (imgpos(2) .ne. lcob(2) .and. imgpos(3) .ne. lcob(3))           &
      coarray(lbr(1):ubr(1),lbv(2),lbv(3),:) =                         &
      coarray(lbr(1):ubr(1),ubr(2),ubr(3),:)                           &
        [imgpos(1),imgpos2mns1,imgpos3mns1]
    
    if (imgpos(2) .ne. ucob(2) .and. imgpos(3) .ne. lcob(3))           &
      coarray(lbr(1):ubr(1),ubv(2),lbv(3),:) =                         &
      coarray(lbr(1):ubr(1),lbr(2),ubr(3),:)                           &
        [imgpos(1),imgpos2pls1,imgpos3mns1]
    
    if (imgpos(2) .ne. ucob(2) .and. imgpos(3) .ne. ucob(3))           &
      coarray(lbr(1):ubr(1),ubv(2),ubv(3),:) =                         &
      coarray(lbr(1):ubr(1),lbr(2),lbr(3),:)                           &
        [imgpos(1),imgpos2pls1,imgpos3pls1]
    
    if (imgpos(2) .ne. lcob(2) .and. imgpos(3) .ne. ucob(3))           &
      coarray(lbr(1):ubr(1),lbv(2),ubv(3),:) =                         &
      coarray(lbr(1):ubr(1),ubr(2),lbr(3),:)                           &
        [imgpos(1),imgpos2mns1,imgpos3pls1]

!  else if ( idx .eq. 2 ) then
    ! 2nd group of remote calls

    ! exchange 2D halos in direction 2

    if ( imgpos(2) .ne. lcob(2) )                                      &
      coarray( lbr(1) : ubr(1) , lbv(2) , lbr(3) : ubr(3) , : ) =      &
      coarray( lbr(1) : ubr(1) , ubr(2) , lbr(3) : ubr(3) , : )        &
        [ imgpos(1) , imgpos2mns1 , imgpos(3) ]
    
    if ( imgpos(2) .ne. ucob(2) )                                      &
      coarray( lbr(1) : ubr(1) , ubv(2) , lbr(3) : ubr(3) , : ) =      &
      coarray( lbr(1) : ubr(1) , lbr(2) , lbr(3) : ubr(3) , : )        &
        [ imgpos(1) , imgpos2pls1 , imgpos(3) ]

    ! exchange 1D halos parallel to direction 2
    
    if (imgpos(1) .ne. lcob(1) .and. imgpos(3) .ne. lcob(3))           &
      coarray(lbv(1),lbr(2):ubr(2),lbv(3),:) =                         &
      coarray(ubr(1),lbr(2):ubr(2),ubr(3),:)                           &
        [imgpos1mns1,imgpos(2),imgpos3mns1]
    
    if (imgpos(1) .ne. ucob(1) .and. imgpos(3) .ne. lcob(3))           &
      coarray(ubv(1),lbr(2):ubr(2),lbv(3),:) =                         &
      coarray(lbr(1),lbr(2):ubr(2),ubr(3),:)                           &
        [imgpos1pls1,imgpos(2),imgpos3mns1]
    
    if (imgpos(1) .ne. ucob(1) .and. imgpos(3) .ne. ucob(3))           &
      coarray(ubv(1),lbr(2):ubr(2),ubv(3),:) =                         &
      coarray(lbr(1),lbr(2):ubr(2),lbr(3),:)                           &
        [imgpos1pls1,imgpos(2),imgpos3pls1]
    
    if (imgpos(1) .ne. lcob(1) .and. imgpos(3) .ne. ucob(3))           &
      coarray(lbv(1),lbr(2):ubr(2),ubv(3),:) =                         &
      coarray(ubr(1),lbr(2):ubr(2),lbr(3),:)                           &
        [imgpos1mns1,imgpos(2),imgpos3pls1]

!  else if ( idx .eq. 3 ) then
    ! 3rd group of remote calls

    ! exchange 2D halos in direction 3
    
    if ( imgpos(3) .ne. lcob(3) )                                      &
      coarray( lbr(1) : ubr(1) , lbr(2) : ubr(2) , lbv(3) , : ) =      &
      coarray( lbr(1) : ubr(1) , lbr(2) : ubr(2) , ubr(3) , : )        &
        [ imgpos(1) , imgpos(2) , imgpos3mns1 ]
    
    if ( imgpos(3) .ne. ucob(3) )                                      &
      coarray( lbr(1) : ubr(1) , lbr(2) : ubr(2) , ubv(3) , : ) =      &
      coarray( lbr(1) : ubr(1) , lbr(2) : ubr(2) , lbr(3) , : )        &
        [ imgpos(1) , imgpos(2) , imgpos3pls1 ]

    ! exchange 1D halos parallel to direction 3
    
    if ( imgpos(1) .ne. lcob(1) .and. imgpos(2) .ne. lcob(2) )         &
      coarray( lbv(1) , lbv(2) , lbr(3) : ubr(3) , : ) =               &
      coarray( ubr(1) , ubr(2) , lbr(3) : ubr(3) , : )                 &
        [ imgpos1mns1 , imgpos2mns1 , imgpos(3) ]
    
    if (imgpos(1) .ne. ucob(1) .and. imgpos(2) .ne. lcob(2))           &
      coarray(ubv(1),lbv(2),lbr(3):ubr(3),:) =                         &
      coarray(lbr(1),ubr(2),lbr(3):ubr(3),:)                           &
        [imgpos1pls1,imgpos2mns1,imgpos(3)]
    
    if (imgpos(1) .ne. ucob(1) .and. imgpos(2) .ne. ucob(2))           &
      coarray(ubv(1),ubv(2),lbr(3):ubr(3),:) =                         &
      coarray(lbr(1),lbr(2),lbr(3):ubr(3),:)                           &
        [imgpos1pls1,imgpos2pls1,imgpos(3)]
    
    if (imgpos(1) .ne. lcob(1) .and. imgpos(2) .ne. ucob(2))           &
      coarray(lbv(1),ubv(2),lbr(3):ubr(3),:) =                         &
      coarray(ubr(1),lbr(2),lbr(3):ubr(3),:)                           &
        [imgpos1mns1,imgpos2pls1,imgpos(3)]

    ! exchange 8 scalar halos
    ! first 4 statements are for the lower bound virtual
    ! corners of the coarray along dimension 3.
    
    if ((imgpos(1) .ne. lcob(1)) .and. (imgpos(2) .ne. lcob(2)) .and.  &
        (imgpos(3) .ne. lcob(3)))                                      &
      coarray( lbv(1), lbv(2), lbv(3), : ) =                           &
      coarray( ubr(1), ubr(2), ubr(3), : )                             &
                         [ imgpos1mns1, imgpos2mns1, imgpos3mns1 ]
    
    if ((imgpos(1) .ne. ucob(1)) .and. (imgpos(2) .ne. lcob(2)) .and.  &
        (imgpos(3) .ne. lcob(3)))                                      &
      coarray( ubv(1), lbv(2), lbv(3), : ) =                           &
      coarray( lbr(1), ubr(2), ubr(3), : )                             &
                         [ imgpos1pls1, imgpos2mns1, imgpos3mns1 ]
    
    if ((imgpos(1) .ne. lcob(1)) .and. (imgpos(2) .ne. ucob(2)) .and.  &
        (imgpos(3) .ne. lcob(3)))                                      &
      coarray( lbv(1), ubv(2), lbv(3), : ) =                           &
      coarray( ubr(1), lbr(2), ubr(3), : )                             &
                         [ imgpos1mns1, imgpos2pls1, imgpos3mns1 ]
    
    if ((imgpos(1) .ne. ucob(1)) .and. (imgpos(2) .ne. ucob(2)) .and.  &
        (imgpos(3) .ne. lcob(3)))                                      &
      coarray( ubv(1), ubv(2), lbv(3), : ) =                           &
      coarray( lbr(1), lbr(2), ubr(3), : )                             &
                         [ imgpos1pls1, imgpos2pls1, imgpos3mns1 ]
    
    ! these 4 statements are for the uppper bound virtual
    ! corners of the coarray along dimension 3.
    
    if ((imgpos(1) .ne. lcob(1)) .and. (imgpos(2) .ne. lcob(2)) .and.  &
        (imgpos(3) .ne. ucob(3)))                                      &
      coarray (lbv(1), lbv(2), ubv(3), : ) =                           &
      coarray (ubr(1), ubr(2), lbr(3), : )                             &
                         [ imgpos1mns1, imgpos2mns1, imgpos3pls1 ]
    
    if ((imgpos(1) .ne. ucob(1)) .and. (imgpos(2) .ne. lcob(2)) .and.  &
        (imgpos(3) .ne. ucob(3)))                                      &
      coarray (ubv(1), lbv(2), ubv(3), : ) =                           &
      coarray (lbr(1), ubr(2), lbr(3), : )                             &
                         [ imgpos1pls1, imgpos2mns1, imgpos3pls1 ]
    
    if ((imgpos(1) .ne. lcob(1)) .and. (imgpos(2) .ne. ucob(2)) .and.  &
        (imgpos(3) .ne. ucob(3)))                                      &
      coarray (lbv(1), ubv(2), ubv(3), : ) =                           &
      coarray (ubr(1), lbr(2), lbr(3), : )                             &
                         [ imgpos1mns1, imgpos2pls1, imgpos3pls1 ]
    
    if ((imgpos(1) .ne. ucob(1)) .and. (imgpos(2) .ne. ucob(2)) .and.  &
        (imgpos(3) .ne. ucob(3)))                                      &
      coarray (ubv(1), ubv(2), ubv(3), : ) =                           &
      coarray (lbr(1), lbr(2), lbr(3), : )                             &
                         [ imgpos1pls1, imgpos2pls1, imgpos3pls1 ]

  end if groups

end do rcalls

end procedure cgca_hxir

!*roboend*

end submodule m2hx_hxir
