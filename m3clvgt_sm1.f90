!*robodoc*f* cgca_m3clvg/m3clvgt_sm1
!  NAME
!    m3clvgt_sm1
!  SYNOPSIS

!$Id: m3clvgt_sm1.f90 380 2017-03-22 11:03:09Z mexas $

submodule ( cgca_m3clvgt ) m3clvgt_sm1

!  DESCRIPTION
!    Submodule of cgca_m3clvg with collective routines.
!    This module cannot be used (yet) with ifort 16,
!    so don't build there.
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENCE
!  CONTAINS
!    cgca_clvgp
!  USES
!    Variables and parameters from the parent module cgca_m3clvg.
!  USED BY
!    The parent module cgca_m3clvg.
!  SOURCE

implicit none

contains

!*roboend*


!*robodoc*s* m3clvg_sm3/cgca_clvgpt
!  NAME
!    cgca_clvgpt
!  SYNOPSIS

module procedure cgca_clvgpt

!  INPUTS
!    See the interface in the parent module cgca_m3clvg.
!  OUTPUTS
!    None
!  SIDE EFFECTS
!    Many:
!    - change state of coarray
!    - change state of gc array
!  DESCRIPTION
!    This is a cleavage propagation routine.
!    We copy the model (coarray) into the local array.
!    We then analyse the local array, but update the coarray.
!  NOTES
!    This is a *thread safe* cleavage propagation routine.
!    This means that the order of iterations in the inner loop
!    is not important.
!
!    All images must call this routine
!  USES
!    cgca_clvgs_abstract, cgca_clvgsd, cgca_clvgsp, cgca_clvgn,
!    cgca_hxi, cgca_hxg, cgca_dacf
!  USED BY
!    none, end user
!  SOURCE

real( kind=rdef ) :: n(3)
integer( kind=iarr ), allocatable :: array(:,:,:)
integer( kind=iarr ) :: grain, cstate

integer(kind=idef) ::                                              &
  lbv(4) ,& ! lower bounds of the complete (plus virtual) coarray
  ubv(4) ,& ! upper bounds of the complete (plus virtual) coarray
  lbr(4) ,& ! lower bounds of the "real" coarray, lower virtual+1
  ubr(4) ,& ! upper bounds of the "real" coarray, upper virtual-1
  x1     ,& ! local coordinates in an array, which are also
  x2     ,& ! do loop counters
  x3     ,& !
  iteration ! iteration counter

integer :: thisimage, errstat=0, nimages, ierr
integer, save :: clvgglob[*]

logical(kind=ldef) :: clvgflag

! Make sure to allocate gcupd!
if ( .not. allocated( gcupd ) ) call gcupd_alloc

! Set the global cleavage flag initially to zero on all images,
! i.e. no cleavage
clvgglob = 0

! Set the local cleavage flag to .false.
clvgflag = .false. 

! use local vars to save time
thisimage = this_image()
  nimages = num_images()

! determine the extents
lbv = lbound(coarray)
ubv = ubound(coarray)
lbr = lbv+1
ubr = ubv-1

!*************************************************
! Sanity checks
!*************************************************

! check for allocated
if ( .not. allocated( coarray ) ) then
  write (*,'(a,i0)') "ERROR: cgca_clvgp: image ",thisimage
  write (*,'(a)')    "ERROR: cgca_clvgp: coarray is not allocated"
  error stop
end if
if ( .not. allocated( rt ) ) then
  write (*,'(a,i0)') "ERROR: cgca_clvgp: image ",thisimage
  write (*,'(a)')    "ERROR: cgca_clvgp: rt is not allocated"
  error stop
end if

! check there are no liquid cells in the grain array
if ( any( coarray(lbr(1):ubr(1),lbr(2):ubr(2),lbr(3):ubr(3),           &
  cgca_state_type_grain) .eq. cgca_liquid_state)) then
  write (*,'(a,i0,a)') "ERROR: cgca_clvgp: image ",thisimage,          &
   ": liquid phase in the model"
  error stop
end if

if ( iter .le. 0 ) then
  write (*,'(a,i0,a)') "ERROR: cgca_clvgp: image ",thisimage,          &
   ": negative number of iterations given"
  error stop
end if

!*************************************************
! End of sanity checks
!*************************************************

! allocate the temp array
allocate( array(lbv(1):ubv(1),lbv(2):ubv(2),lbv(3):ubv(3) ),           &
          stat=errstat )
if ( errstat .ne. 0 ) then
  write (*,"(2(a,i0))") "ERROR: cgca_clvgp: image ", thisimage,        &
   " : cannot allocate array, errcode: ", errstat
  error stop
end if

! initialise the iteration counter
iteration = 1

! initial halo exchange, to make sure the coarray is in a correct state
call cgca_hxi( coarray )
if ( periodicbc ) call cgca_hxg( coarray )
sync all

! start the main loop for cleavage iterations
main: do

  ! copy coarray fracture state type into a local array
  array = coarray(:,:,:,cgca_state_type_frac)

  ! propagate cleavage
  do concurrent( x1=lbr(1):ubr(1), x2=lbr(2):ubr(2), x3=lbr(3):ubr(3) )

    ! scan only through undamaged cells
    live: if ( array(x1,x2,x3) .eq. cgca_intact_state .or.             &
               array(x1,x2,x3) .eq. cgca_gb_state_intact) then

      grain = coarray( x1, x2, x3, cgca_state_type_grain )
      call cgca_clvgn_pure( t, grain, rt, scrit, clvgflag, n, cstate,  &
        ierr )

      ! debug
      !   if (debug) write (*,"(a,i0,a,l1)")      &
      !    "DEBUG: cgca_clvgp: img ", thisimage,  &
      !    " clvgflag=", clvgflag

      ! If cleavage conditions are met, propagate cleavage into
      ! this cell. Note that we pass the local array, but return
      ! the new state of the central cell into the coarray.
      ! The sub name is provided as an input to cgca_clvgp.
      ! It can be either the deterministic routine cgca_clvgsd,
      ! or the probabilistic routine cgca_clvgsp.
      if ( clvgflag ) then

        ! Mark that cleavage has occurred. The value is not important,
        ! any non-zero integer will do, but the same on all images.
        clvgglob = 1
        call sub(   array(x1-1:x1+1, x2-1:x2+1, x3-1:x3+1),            &
                  coarray(x1-1:x1+1, x2-1:x2+1, x3-1:x3+1,             &
                                       cgca_state_type_grain),         &
                 n, cstate, debug,                                     &
                    coarray(x1,x2,x3,cgca_state_type_frac) )
      end if
    end if live

  end do

  if ( ierr .ne. 0 ) then
    write (*,'(a,i0)') "ERROR: cgca_clvgp_t/m3clvg_sm3:" //            &
      " cgca_clvgn_pure error, ierr:", ierr
    error stop
  end if

  ! Add together all cleavage identifiers from all images
  ! no sync is required!
  call co_sum( clvgglob )

  ! Check if cleavage happened anywhere in the model.
  if ( clvgglob .eq. 0 ) then
    if ( thisimage .eq. 1 )                                            &
      write (*,*) "INFO: cgca_clvgp: no cleavage anywhere, leaving"
    exit main
  end if

  sync all
   
  ! update all local GC arrays using the given subroutine
  call gcus( periodicbc )

  ! halo exchange after a cleavage propagation step
  call cgca_hxi( coarray )
  if ( periodicbc ) call cgca_hxg( coarray )

  ! deactivate crack flanks, ignore grain boundaries
  call cgca_dacf( coarray, debug=.false. )

  sync all
 
  ! Reset all gcupd
  gcupd = cgca_gb_state_intact

  ! Reset the gcupd counter
  gcucnt = 1

  ! halo exchange after deactivation step
  call cgca_hxi( coarray )
  if ( periodicbc ) call cgca_hxg( coarray )

  sync all

  ! send heartbeat signal to terminal
  if (thisimage .eq. 1 .and. heartbeat .gt. 0) then
    if ( mod( iteration, heartbeat ) .eq. 0) write (*,'(a,i0)')        &
        "INFO: cgca_clvgp: iterations completed: ", iteration
  end if

  if ( iteration .ge. iter ) exit main

  ! increment the iteration counter
  iteration = iteration + 1

end do main

deallocate( array, stat=errstat )
if ( errstat .ne. 0 ) then
  write (*,"(2(a,i0))") "ERROR: cgca_clvgp: image ",thisimage,         & 
   " : cannot deallocate array, errcode: ", errstat
  error stop
end if

! sync before leaving
sync all

end procedure cgca_clvgpt

!*roboend*

end submodule m3clvgt_sm1
