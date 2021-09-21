!*robodoc*m* CGPACK/cgca_m2gb
!  NAME
!    cgca_m2gb
!  SYNOPSIS

!$Id: cgca_m2gb.f90 529 2018-03-26 11:25:45Z mexas $

module cgca_m2gb

!  DESCRIPTION
!    Module dealing with grain boundaries.
!    Most routines in this module are concerned with creating,
!    updating and printing of the *local* array gc.
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    cgca_gcu, cgca_gcp,
!    cgca_dgc, cgca_gcf, cgca_gcr, cgca_gbs,
!    cgca_agc (private to this module), cgca_gcr_pure, cgca_gcf_pure
!  USES
!    cgca_m1co
!  USED BY
!    cgca_m3clvg
!  SOURCE

use cgca_m1co
implicit none

private
public :: cgca_dgc, cgca_gbs, cgca_gcf, cgca_gcp, cgca_gcr,            &
          cgca_gcu, cgca_igb, &
! Do not use! NOt READY!
 cgca_gcr_pure, cgca_gcf_pure


!*roboend*


!*robodoc*d* cgca_m2gb/gc
!  NAME
!    gc
!  SYNOPSIS

integer( kind=iarr ), allocatable, save :: gc(:,:)

!  DESCRIPTION
!    Local, *not* coarray, grain connectivity (GC) array. This
!    array must be SAVEd.
!  NOTES
!    GC is a private array. It is not accessible from outside of
!    cgca_m2gb module, hence we can use a simple name, with no
!    cgca_ prefix.
!    GC is zero initially!
!  USED BY
!    All routines of module cgca_m2gb:
!    cgca_gcu, cgca_gcp, cgca_agc, cgca_dgc, cgca_gcf, cgca_gcr.
!*roboend*

contains

!*robodoc*s* cgca_m2gb/cgca_gcu
!  NAME
!    cgca_gcu
!  SYNOPSIS

subroutine cgca_gcu( coarray )

!  INPUT

integer( kind=iarr ), intent( in ), allocatable ::                     &
  coarray(:,:,:,:) [:,:,:]

!  SIDE EFFECTS
!    Updated state of gc array
!  DESCRIPTION
!    Update a local grain connectivity (GC) array. This means:
!    scan the whole of the local real (no halos) model coarray.
!    For each cell check all 26 neighbours.
!    When a cell has a neighbour of a different number, this
!    is understood as a grain boundary between the grains
!    denoted by the states of both cells. If this pair is not
!    already in GC, then it is added to it in the right place.
!    GC is sorted, first by the first column, then by the second.
!
!    state(26,2) - 26 possible pairs of different states, and 2 states.
!    In practice, 26 is a stupid value.
!    This can only happen if a grain is only 1 cell, and if each of the
!    26 neighbouring cells has a unique different value.
!    Anyway, defining the array of this size seems fool proof. 
!
!    The 3rd entry in each row of GC is the grain boundary integrity,
!    either cgca_gb_state_intact - intact or cgca_gb_state_fractured -
!    fractured.
!  NOTES
!    Only the (:,:,:,cgca_state_type_grain) values are used in this
!    routine.
!    The cgca_state_type_frac values are *not* used.
!  USES
!    gc, cgca_agc, cgca_dgc
!  USED BY
!    via module cgca_m2gb
!  SOURCE

integer( kind=iarr ), allocatable :: tmp(:,:)
integer( kind=iarr ) :: state(26,2) = 0_iarr
integer :: errstat=0, x1, x2, x3, n1, n2, n3, lbr(4), ubr(4), pair,    &
 con, img, gclen, i, j

img = this_image()

!*********************************************************************72
! checks
!*********************************************************************72

if ( .not. allocated(coarray) ) then
 write (*,'(a,i0)') "ERROR: cgca_gc: coarray not allocated, img: ", img
 error stop
end if

!*********************************************************************72
! end of checks
!*********************************************************************72

! if GC not already allocated, allocate to zero length!
if ( .not.  allocated( gc ) ) call cgca_agc(0)

! Assume the coarray has halos. Ignore those
lbr = lbound( coarray ) + 1
ubr = ubound( coarray ) - 1

       do x3 = lbr(3), ubr(3)
       do x2 = lbr(2), ubr(2)
inner: do x1 = lbr(1), ubr(1)

  ! Loop over 26 neighbours, find all neighbours with different states
  pair = 0
  do n3 = x3-1, x3+1
  do n2 = x2-1, x2+1
  do n1 = x1-1, x1+1

    ! if the states differ but no liquid
    if ( ( coarray(x1,x2,x3,cgca_state_type_grain) .ne.                &
           coarray(n1,n2,n3,cgca_state_type_grain) )                   &
         .and.                                                         &
         (coarray(n1,n2,n3,cgca_state_type_grain) .ne.                 &
          cgca_liquid_state) ) then

      ! this is another pair
      pair = pair + 1

      ! add this pair to state array
      state( pair, 1 ) = coarray( x1, x2, x3, cgca_state_type_grain )
      state( pair, 2 ) = coarray( n1, n2, n3, cgca_state_type_grain )

      ! Make state(pair,1) < state(pair,2)
      if ( coarray(n1,n2,n3,cgca_state_type_grain) .lt.                &
           coarray(x1,x2,x3,cgca_state_type_grain) ) then 

        ! swap them
        state(pair,1) = coarray(n1,n2,n3,cgca_state_type_grain)
        state(pair,2) = coarray(x1,x2,x3,cgca_state_type_grain)
      end if
    end if

  end do
  end do
  end do

  ! At the end of the this loop, "state" array will have
  ! all pairs of different grains containg the central cell.
  ! Now all of these pair, which are not already in GC, need
  ! to be added there.

  ! For all identified pairs of states
  newcon: do con = 1, pair

    ! Check if this connectivity is already in the GC array.
    ! If yes, cycle to the next pair.
    gclen = ubound( gc, dim=1 )

    do i = 1, gclen
      ! each pair repeats twice, so the order doesn't matter
      ! If the pair is already in GC, then cycle
      if ( gc(i,1) .eq. state(con,1) .and. &
           gc(i,2) .eq. state(con,2) ) cycle newcon
    end do

    ! If this connectivity is not yet in the GC array, then extend
    ! the GC array by 2 rows.
    gclen = gclen + 2

! this should be replaced by MOVE_ALLOC!!!

    ! Allocate a temp array with length 2 more than GC, and set
    ! to intact.
    allocate( tmp(gclen, 3), source=0_iarr, stat=errstat )
    if (errstat .ne. 0) then
      write (*,'(2(a,i0))') "ERROR: cgca_gc: img: ", img,              &
                        " cannot allocate tmp, error code: ", errstat
      error stop
    end if

    ! Copy GC to the beginning of the temp array
    tmp(1:gclen-2,:) = gc

    ! Can replace this with move_alloc, but I make use of the
    ! temp array later, so this is better!
    call cgca_dgc
    call cgca_agc( gclen )
    gc = tmp

!********************************************************
! debug output
!
!if (this_image() .eq. 1) then
! do k=1,gclen
!  write (*,*) gc(k,:)
! end do
! write (*,*) "state1,state2",state(con,1),state(con,2)
!end if
!********************************************************

    ! find where to insert the first pair: state(con,1), state(con,2)
    i1: do i=1,gclen

      ! If state(con,1) doesn't exist in the first row at all, but
      ! is smaller that some existing value
      if (state(con,1) .lt. gc(i,1)) then
        tmp(1:gclen-i-1,:) = gc(i:gclen-2,:)
        gc(i,1) = state(con,1)
        gc(i,2) = state(con,2)
        gc(i,3) = cgca_gb_state_intact
        gc(i+1:gclen-1,:) = tmp(1:gclen-i-1,:)
        exit i1
      else if (state(con,1) .eq. gc(i,1)) then

        ! If state(con,1) already exists in the first column, then
        ! need to sort by the second column
        do j=i,gclen
          if (state(con,2) .lt. gc(j,2)) then
            tmp(1:gclen-j-1,:) = gc(j:gclen-2,:)
            gc(j,1) = state(con,1)
            gc(j,2) = state(con,2)
            gc(j,3) = cgca_gb_state_intact
            gc(j+1:gclen-1,:) = tmp(1:gclen-j-1,:)
            exit i1
          else if (state(con,2) .gt. gc(j,2) .and.    &
                   gc(j+1,1) .ne. gc(j,1) ) then

            ! If state(con,2) is greater than all column 2 values,
            ! for the same column 1 value, just add it *after* this
            ! row.
            tmp(1:gclen-j-2,:) = gc(j+1:gclen-2,:)
            gc(j+1,1) = state(con,1)
            gc(j+1,2) = state(con,2)
            gc(j+1,3) = cgca_gb_state_intact
            gc(j+2:gclen-1,:) = tmp(1:gclen-j-2,:)
            exit i1
          end if
        end do
      else if (i .gt. gclen-2) then

        ! If state(con,1) is greater than all column 1 values, just
        ! add it at the end. 
        gc(i,1) = state(con,1)
        gc(i,2) = state(con,2)
        gc(i,3) = cgca_gb_state_intact
        exit i1
      end if
    end do i1 

  ! Find where to insert the second pair:
  ! ( state(con,2), state(con,1) ).
  ! Note that as the first pair has been inserted already,
  ! this loop is not exactly as the previous one. 
  i2: do i=1,gclen
    ! If state(con,2) doesn't exist in the first row at all, but
    ! is smaller that some existing value
    if (state(con,2) .lt. gc(i,1)) then
      tmp(1:gclen-i,:) = gc(i:gclen-1,:)
      gc(i,1) = state(con,2)
      gc(i,2) = state(con,1)
      gc(i,3) = cgca_gb_state_intact
      gc(i+1:gclen,:) = tmp(1:gclen-i,:)
      exit i2
    else if (state(con,2) .eq. gc(i,1)) then
      ! If state(con,2) already exist in the first column, then
      ! need to sort by the second column by state(con,1)

      do j=i,gclen
        if (state(con,1) .lt. gc(j,2)) then
          tmp(1:gclen-j,:) = gc(j:gclen-1,:)
          gc(j,1) = state(con,2)
          gc(j,2) = state(con,1)
          gc(j,3) = cgca_gb_state_intact
          gc(j+1:gclen,:) = tmp(1:gclen-j,:)
          exit i2
        else if (state(con,1) .gt. gc(j,2) .and.   &
                 gc(j+1,1) .ne. gc(j,1) ) then
          ! If state(con,1) is greater than all column 2 values, just
          ! add it *after* this row.
          tmp(1:gclen-j-1,:) = gc(j+1:gclen-1,:)
          gc(j+1,1) = state(con,2)
          gc(j+1,2) = state(con,1)
          gc(j+1,3) = cgca_gb_state_intact
          gc(j+2:gclen,:) = tmp(1:gclen-j-1,:)
          exit i2
        end if
      end do
    else if (i .gt. gclen-1) then
      ! If state(con,2) is greater than all column 1 values, just
      ! add it at the end. 
      gc(i,1) = state(con,2)
      gc(i,2) = state(con,1)
      gc(i,3) = cgca_gb_state_intact
      exit i2
    end if
  end do i2

  deallocate(tmp, stat=errstat )
  if (errstat .ne. 0) then
    write (*,'(2(a,i0))') "ERROR: cgca_gcu: img: ", img,               &
                          " cannot deallocate tmp, err code: ", errstat
    error stop
  end if

end do newcon
end do inner
end do
end do

end subroutine cgca_gcu

!*roboend*


!*robodoc*s* cgca_m2gb/cgca_gcp
!  NAME
!    cgca_gcp
!  SYNOPSIS

subroutine cgca_gcp( ounit, fname )

!  INPUTS

integer( kind=idef ), intent(in) :: ounit
character( len=* ), intent(in) :: fname

!  SIDE EFFECTS
!    create output file on each processor and dump gc array to it
!  DESCRIPTION
!    Print grain connectivity (GC) array to file.
!    Remember that GC is a *local* array, and hence this routine will
!    write only the GC from the processor which called this routine.
!    It is indended to be called from all images with file names supplied
!    linked to the processor/image number via this_image.
!  USES
!    gc
!  USED BY
!    none 
!  SOURCE

integer :: i, errstat=0, img
character( len=10) :: state

img = this_image()

! open file
open( unit=ounit, file=trim(fname), form="formatted",                  &
      status="replace", iostat=errstat )
if ( errstat .ne. 0 ) then
  write (*,'(2(a,i0))') "ERROR: cgca_gcp: img: ", img,                 &
                        " cannot open file, error code: ", errstat
  error stop
end if

! write data, one line at a time
do i = 1, ubound( gc, dim=1 )

  ! convert numerical values of intact and fractured into words
  if ( gc(i,3) .eq. cgca_gb_state_intact) then
    state = "intact"
  else if ( gc(i,3) .eq. cgca_gb_state_fractured) then
    state = "fractured"
  else
    write (*,'(2(a,i0))') "ERROR: cgca_gcp: image: ", img,             &
                   " state is neither intact nor fractured: ", gc(i,3)
    error stop
  end if

  ! actual write
  write ( unit=ounit, fmt="(2(i0,tr1),a10)", iostat=errstat)           &
   gc(i,1:2), state

  if ( errstat .ne. 0 ) then
    write (*,'(2(a,i0))') "ERROR: cgca_gcp: image: ", img,             &
                      " cannot write to file, error code: ", errstat
    error stop
  end if

end do

! close file
close( ounit, iostat=errstat )
if ( errstat .ne. 0 ) then
  write (*,'(2(a,i0))') "ERROR: cgca_gcp: image: ", img,               &
                        " cannot close file, error code: ", errstat
  error stop
end if

end subroutine cgca_gcp

!*roboend*


!*robodoc*s* cgca_m2gb/cgca_agc
!  NAME
!    cgca_agc
!  SYNOPSIS

subroutine cgca_agc( len )

!  INPUT

integer, intent( in ) :: len

!  SIDE EFFECTS
!    Allocate gc array
!  DESCRIPTION
!    Allocate the grain connectivity (GC) array, and set to zero.
!    The first dimension is given by len.
!    The second dimension is 3, because:
!     1 - grain number
!     2 - grain neighbour number
!     3 - grain boundary state - fractured or intact
!  USES
!    gc
!  USED BY
!    module cgca_m2gb: cgca_gcu
!  SOURCE

integer :: errstat=0

allocate( gc(len,3), source=0_iarr, stat=errstat )
if ( errstat .ne. 0 ) then
  write (*,'(2(a,i0))') "ERROR: cgca_agc: img: ", this_image(),        &
                        " cannot allocate gc, error code: ", errstat
  error stop
end if

end subroutine cgca_agc

!*roboend*


!*robodoc*s* cgca_m2gb/cgca_dgc
!  NAME
!    cgca_dgc
!  SYNOPSIS

subroutine cgca_dgc

!  SIDE EFFECTS
!    deallocate gc array
!  USES
!    gc
!  USED BY
!    module cgca_m2gb: cgca_gcu
!  SOURCE

integer :: errstat

deallocate( gc, stat=errstat )
if (errstat .ne. 0) then
  write (*,'(2(a,i0))') "ERROR: cgca_dgc: img: ", this_image(),        &
                        " cannot deallocate gc, err code: ", errstat
  error stop
end if

end subroutine cgca_dgc

!*roboend*


!*robodoc*s* cgca_m2gb/cgca_gcf
!  NAME
!    cgca_gcf
!  SYNOPSIS

subroutine cgca_gcf( g1, g2 )

!  INPUTS

integer( kind=iarr ), intent( in ) :: g1, g2

!  SIDE EFFECTS
!    gc array modified
!  DESCRIPTION
!    This routine changes the state of the grain boundary integrity
!    between the two provided grains to cgca_gb_state_fractured.
!  USES
!    gc
!  USED BY
!    cgca_clvgsd, cgca_clvgsp, cgca_gcupdn
!  SOURCE

integer( kind=iarr ) :: tmp, grain1, grain2
integer :: matches, img, i, gclen, match1

img = this_image()

!*********************************************************************72
! sanity checks
!*********************************************************************72

! grains must be different. However, if there is no sync between images
! it's possible that 2 different images write to gcupd simultaneously,
! thus creating impossible combinations of GB pairs, including
! identical grain numbers for both grains. Let's issue a warning, not
! error for this case, and return, i.e. don't add this impossible pair.
if ( g1 .eq. g2 ) then
 write (*,"(a,i0,tr1,i0,a,i0)")                                        &
   "WARN: cgca_gcf: identical grain numbers: ", g1, g2, " image ", img
 return
 ! error stop
end if

!*********************************************************************72
! end of sanity checks
!*********************************************************************72

! Use the local variables
grain1 = g1
grain2 = g2

! Make grain1 < grain2
if ( grain2 .lt. grain1 ) then
  tmp = grain1
  grain1 = grain2
  grain2 = tmp
end if

! Look for matches
matches = 0
 match1 = 0_iarr
  gclen = ubound(gc,1)

! First match
do i = 1, gclen
  if ( gc(i,1) .eq. grain1 .and. gc(i,2) .eq. grain2 ) then
    gc(i,3) = cgca_gb_state_fractured
    matches = matches+1
    ! line number of the first match
    match1 = i
    exit
  end if
end do

! Second match, start from the line with the first match
do i = match1 + 1, gclen
  if ( gc(i,1) .eq. grain2 .and. gc(i,2) .eq. grain1 ) then
    gc(i,3) = cgca_gb_state_fractured
    matches = matches+1
    exit
  end if
end do

! Sanity check
if ( matches .eq. 0 ) then
 ! don't issue the INFO message, just swamps the stdout
 ! perhaps better logic should be built in future
 ! write (*,'(2(a,i0),tr1,i0,a)') "INFO: cgca_gcf: image ",              &
 ! img, " pair ", g1, g2, " does not exist."
else if ( matches .ne. 2 ) then
 write (*,'(3(a,i0),",",tr1,i0,".")') "ERROR: cgca_gcf: image ",       &
  img, ": found ", matches, " matches for pair: ", g1, g2
 write (*,'(a,i0,a)') "ERROR: cgca_gcf: image ", img,                  &
  ": Should have found exactly two matches or none!"
 write (*,'(a,i0,a)') "ERROR: cgca_gcf: image ", img,                  &
  ": This means the gc array is corrupted. Aborting!"
 error stop
end if

end subroutine cgca_gcf

!*roboend*


!*robodoc*s* cgca_m2gb/cgca_gcf_pure
!  NAME
!    cgca_gcf_pure
!  SYNOPSIS

subroutine cgca_gcf_pure( g1, g2, iflag )

!    Despite the name, this subroutine is actually not PURE!
!    I overlooked a restriction that a pure subroutine cannot
!    modify a variable accissible via host association, which
!    is gc in this case.
!  INPUTS

integer( kind=iarr ), intent(in) :: g1, g2

!  OUTPUTS

integer, intent(out) :: iflag

!  SIDE EFFECTS
!    gc array modified
!  DESCRIPTION
!    This routine changes the state of the grain boundary integrity
!    between the two provided grains to cgca_gb_state_fractured.
!    IFLAG:
!     0 - successful termination
!     1 - the 2 grains are identical - probably fatal error
!     2 - there is no such pair in the array - probably a warning
!     3 - the number of grain pairs is not 2 - probably a fatal error
!  USES
!    gc
!  USED BY
!    cgca_clvgsd, cgca_clvgsp, cgca_gcupdn
!  SOURCE

integer( kind=iarr ) :: tmp, grain1, grain2
integer :: matches, img, i, gclen, match1

img = this_image()

! The two grains must be different. If not set iflag to 1 and return
! immediately.
if ( g1 .eq. g2 ) then
  iflag = 1
  return
end if

! Use the local variables
grain1 = g1
grain2 = g2

! Make grain1 < grain2
if ( grain2 .lt. grain1 ) then
  tmp = grain1
  grain1 = grain2
  grain2 = tmp
end if

! Look for matches
matches = 0
 match1 = 0_iarr
  gclen = ubound(gc,1)

! First match
do i = 1, gclen
  if ( gc(i,1) .eq. grain1 .and. gc(i,2) .eq. grain2 ) then

! CANNOT DO THIS in a PURE routine!!!
! READ the pure rules and fix!!!
!    gc(i,3) = cgca_gb_state_fractured

    matches = matches+1
    ! line number of the first match
    match1 = i
    exit
  end if
end do

! Second match, start from the line with the first match
do i = match1 + 1, gclen
  if ( gc(i,1) .eq. grain2 .and. gc(i,2) .eq. grain1 ) then

! CANNOT DO THIS in a PURE routine!!!
! READ the pure rules and fix!!!
!    gc(i,3) = cgca_gb_state_fractured

    matches = matches+1
    exit
  end if
end do

! Sanity check
if ( matches .eq. 0 ) then
 ! Set iflag to 2 and return immediately
 iflag = 2
 return
else if ( matches .ne. 2 ) then
 ! Set iflag to 3 and return immediately
 iflag = 3
 return
end if

end subroutine cgca_gcf_pure

!*roboend*


!*robodoc*s* cgca_m2gb/cgca_gcr
!  NAME
!    cgca_gcr
!  SYNOPSIS

subroutine cgca_gcr( g1, g2, intact )

!  INPUTS

integer( kind=iarr ), intent( in ) :: g1, g2

!  OUTPUT

logical( kind=ldef ), intent(out) :: intact

!  SIDE EFFECTS
!    none
!  DESCRIPTION
!    Find and return the boundary integrity state between
!    2 grains. Returns .true. if the boundary is intact,
!    .false. if it is fractured. The routine does some
!    simple checks and stops with errors where appropriate.
!  USES
!    gc
!  USED BY
!    module cgca_m3clvg: cgca_clvgsd
!  SOURCE

integer( kind=iarr ) :: tmp, grain1, grain2
integer :: img, i
logical :: match

img = this_image()

!**********************************************************************73
! sanity checks
!**********************************************************************73

! grains must be different
if ( g1 .eq. g2 ) then
  write (*,'(2(a,i0),tr1,i0,".")') "ERROR: cgca_gcr: image ", img, &
    ": The two grain numbers must differ: ", g1, g2
  error stop
end if

!**********************************************************************73
! end of sanity checks
!**********************************************************************73

! Use the local variables
grain1 = g1
grain2 = g2

! Make grain1 < grain2
if ( grain2 .lt. grain1 ) then
     tmp = grain1
  grain1 = grain2
  grain2 = tmp
end if

! Look for matches
match=.false.
do i=1, ubound( gc, 1 )
 mtch: if ( gc(i,1) .eq. grain1 .and. gc(i,2) .eq. grain2 ) then
  match = .true.
  if ( gc(i,3) .eq. cgca_gb_state_intact) then
   intact = .true.
   exit
  else if ( gc(i,3) .eq. cgca_gb_state_fractured) then
   intact = .false.
   exit
  else
   ! Must never end up here. This is an error
   write (*,'(a,i0)')                                                  &
     "ERROR: cgca_gcr: the state is neither intact nor fractured. &
        &Integrity of the gc array is broken. img: ", img
   error stop
  end if
 end if mtch
end do

! Sanity check
if (.not. match) then
  write (*,'(a,i0,a,i0,tr1,i0,".")') "WARN: cgca_gcr: image ", &
    img, ": No match found for given pair: ", g1, g2
  write (*,'(a,i0,a)') "WARN: cgca_gcr: image ", img,      &
    ": Returning .true. in intact!"
  intact = .true.
end if

end subroutine cgca_gcr

!*roboend*


!*robodoc*s* cgca_m2gb/cgca_gcr_pure
!  NAME
!    cgca_gcr_pure
!  SYNOPSIS

pure subroutine cgca_gcr_pure( g1, g2, intact, iflag )

!  INPUTS

integer( kind=iarr ), intent(in) :: g1, g2

!  OUTPUT

logical( kind=ldef ), intent(out) :: intact
integer, intent(out) :: iflag

!  SIDE EFFECTS
!    None, this is a PURE subroutine
!  NOTES
!    Cannot have any external IO in a PURE subroutine. Hence
!    all error conditions are returned back via integer flag,
!    IFLAG. IFLAG meanings:
!      0 - successful completion
!      1 - The 2 grains are identical. This is a fatal condition,
!          the caller routine should probably abort.
!  DESCRIPTION
!    Find and return the boundary integrity state between
!    2 grains. Returns .true. if the boundary is intact,
!    .false. if it is fractured. The routine does some
!    simple checks and stops with errors where appropriate.
!  USES
!    gc
!  USED BY
!    module cgca_m3clvg: cgca_clvgsd
!  SOURCE

integer( kind=iarr ) :: tmp, grain1, grain2
integer :: img, i
logical :: match

img = this_image()

! The 2 grains must be different. If not, set iflag=1, and
! return immediately.
if ( g1 .eq. g2 ) then
  iflag = 1
  return
end if

! Use local variables
grain1 = g1
grain2 = g2

! Make grain1 < grain2
if ( grain2 .lt. grain1 ) then
     tmp = grain1
  grain1 = grain2
  grain2 = tmp
end if

! Look for matches
match=.false.
do i=1, ubound( gc, 1 )
  mtch: if ( gc(i,1) .eq. grain1 .and. gc(i,2) .eq. grain2 ) then
    match = .true.
    if ( gc(i,3) .eq. cgca_gb_state_intact) then
      intact = .true.
      exit
    else if ( gc(i,3) .eq. cgca_gb_state_fractured) then
      intact = .false.
      exit
    else
      ! Must never end up here. This is an error. Set iflag=2 and
      ! return immediately.
      iflag = 2
      return
    end if
  end if mtch
end do

! No match probably means something went wrong. This is probably
! not a fatal condition, but watch out! Set iflag=3 and return
! .true. in inact.
if ( .not. match ) then
  iflag = 3  
  intact = .true.
end if

end subroutine cgca_gcr_pure

!*roboend*


!*robodoc*s* cgca_m2gb/cgca_gbs
!  NAME
!    cgca_gbs
!  SYNOPSIS

subroutine cgca_gbs( coarray )

!  INPUT

integer( kind=iarr ), allocatable, intent(inout) :: &
 coarray(:,:,:,:)[:,:,:]

!  SIDE EFFECTS
!    state of coarray changes
!  DESCRIPTION
!    This routines does Grain Boundary Smoothing (GBS).
!    It works only with cgca_state_type_grain layer.
!    For each cell that has neighbours from other grains the routine
!    substites the cell value (grain number) by that
!    of the grain that has most cells in the (3,3,3) neighbourhood.
!    Example:
!
!         -> 2        5 5 5
!      / |            5 5 5
!     3  V            5 5 5
!        1      5 1 1
!               5 1 1
!               5 5 5
!         1 1 1
!         8 1 1
!         8 1 1
!
!    The central cell is grain 1. In the (3,3,3) neighbourhood,
!    including the central cell, there are 14 cells with grain 5,
!    11 cells with grain 1, and 2 cells with grain 8.
!    The highest number of cells belong to grain 5, hence the
!    central cell state in changed to 5.
!  NOTE
!    There are no remote comms in this routine.
!  SOURCE

integer( kind=iarr ), allocatable :: array(:,:,:)
integer( kind=iarr ) :: neigh(27,2)=0_iarr, newgrain
integer( kind=idef ) :: &
 lbv(4),   & ! lower bounds of the complete (plus virtual) coarray
 ubv(4),   & ! upper bounds of the complete (plus virtual) coarray
 lbr(4),   & ! lower bounds of the "real" coarray, lower virtual+1
 ubr(4),   & ! upper bounds of the "real" coarray, upper virtual-1
 x1,x2,x3, & ! local coordinates in an array, which are also
 n1,n2,n3    ! local coord. of the neighbours [-1,0,1]

integer :: errstat,i

! determine the extents
lbv = lbound(coarray)
ubv = ubound(coarray)
lbr = lbv+1
ubr = ubv-1

! allocate the temp array
allocate( array( lbv(1):ubv(1), lbv(2):ubv(2), lbv(3):ubv(3) ), &
          stat=errstat )
if (errstat .ne. 0) then
  write (*,'(2(a,i0))') "ERROR: cgca_gbs/cgca_m2gb: image: ",          &
    this_image(), " allocate( array ): err. code: ", errstat
  error stop
end if

! copy the grain layer to the temp array
array = coarray(:,:,:,cgca_state_type_grain)

! loop over all cells
do x3 = lbr(3), ubr(3)
do x2 = lbr(2), ubr(2)
do x1 = lbr(1), ubr(1)

  ! Construct the array with all neighbour numbers.
  ! Note, at this stage the array might contain non-grain
  ! phase, e.g. liquid, if the cell is at the model
  ! boundary, with non-grain halo cell neighbours.
  ! It is easier to ignore these cases on calculation.
  ! We'll simply forbid changing into non-grain state later.
  i=0
  do n3 = -1,1
  do n2 = -1,1
  do n1 = -1,1
    i=i+1
    ! The first column of neigh array contains grain numbers
    ! for all 27 cells in the (3,3,3) array
    neigh(i,1) =  array( x1+n1, x2+n2, x3+n3 )
  end do
  end do
  end do
 
  ! Count the number of neighbours of each grain, and
  ! assign to the second column of array neigh.
  do i=1,27
    neigh(i,2) = int( count( neigh(:,1) .eq. neigh(i,1) ), kind=iarr )
  end do
 
  ! Get the grain number that has most cells in the neigh array
  newgrain = neigh( maxloc( neigh(:,2), dim=1 ), 1 )
  if ( newgrain .ne. cgca_liquid_state )               &
    coarray(x1,x2,x3,cgca_state_type_grain) = newgrain

end do
end do
end do

! deallocate the temp array
! The intention is that this routine is called only
! once, or, at best, only a few times, so keeping the
! array allocated for the duration of the execution is a waste.
deallocate( array, stat=errstat )
if (errstat .ne. 0) then
  write (*,'(2(a,i0))') "ERROR: cgca_gbs/cgca_m2gb: image: ",          &
    this_image(), " deallocate( array): err. code: ", errstat
  error stop
end if

end subroutine cgca_gbs

!*roboend*

!*********************************************************************72

!*robodoc*s* cgca_m2gb/cgca_igb
!  NAME
!    cgca_igb
!  SYNOPSIS

subroutine cgca_igb ( coarray )

!  INPUT

integer( kind=iarr ), allocatable, intent(inout) :: &
 coarray(:,:,:,:)[:,:,:]

!  SIDE EFFECTS
!    state of coarray changes
!  DESCRIPTION
!    Initialise Grain Boundary (IGB) cells.
!    Simply scan through the (:,:,:,cgca_state_type_grain) array
!    and mark all cells which have a neighbour of a different state
!    as cgca_gb_state_intact in (:,:,:,cgca_state_type_frac) array.
!    Clearly this routine must be called before any fracture is
!    simulated.
!    Possibly a halo exchange should be called before and/or
!    after calling this routine.
!  NOTES
!    All images must call this routine. No remote comms in this
!    routine.
!  SOURCE

integer( kind=idef ) :: &
 lbv(4),   & ! lower bounds of the complete (plus virtual) coarray
 ubv(4),   & ! upper bounds of the complete (plus virtual) coarray
 lbr(4),   & ! lower bounds of the "real" coarray, lower virtual+1
 ubr(4),   & ! upper bounds of the "real" coarray, upper virtual-1
 x1,x2,x3, & ! local coordinates in an array, which are also
 n1,n2,n3    ! local coord. of the neighbours [-1,0,1]

! determine the extents
lbv = lbound(coarray)
ubv = ubound(coarray)
lbr = lbv+1
ubr = ubv-1

! scan over all cells
outer: do x3 = lbr(3), ubr(3)
       do x2 = lbr(2), ubr(2)
       do x1 = lbr(1), ubr(1)

  ! choose all neighbourhood cells
  inner: do n3 = x3-1, x3+1
         do n2 = x2-1, x2+1
         do n1 = x1-1, x1+1

    ! Ignore the global halo cells
    if ( n1 .eq. lbv(1) .or. n1 .eq. ubv(1) .or.                       &
         n2 .eq. lbv(2) .or. n2 .eq. ubv(2) .or.                       &
         n3 .eq. lbv(3) .or. n3 .eq. ubv(3) ) cycle

    ! If the neighbouring grain .ne. the state of the central
    ! cell, then mark this central cell as GB in fracture
    ! array and exit 
    if ( coarray(n1, n2, n3, cgca_state_type_grain) .ne.               &
         coarray(x1, x2, x3, cgca_state_type_grain) ) then
      coarray(x1, x2, x3, cgca_state_type_frac) = cgca_gb_state_intact
      exit inner
    end if

  end do
  end do
  end do inner

end do
end do
end do outer

end subroutine cgca_igb

!*roboend*

!*********************************************************************72

end module cgca_m2gb
