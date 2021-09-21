!*robodoc*f* cgca_m3sld/m3sld_hc
!  NAME
!    m3sld_hc
!  SYNOPSIS

!$Id: m3sld_hc.f90 431 2017-06-30 13:13:49Z mexas $

submodule ( cgca_m3sld ) m3sld_hc

!  DESCRIPTION
!    Submodule with solidification routines which offer a choice
!    of halo exchange routines. Fortran 2015 collectives are used.
!  NOTES
!    Not supported by ifort 16.
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    cgca_sld3
!  USES
!    All date objects of the parent module cgca_m3sld
!    by host association.
!  USED BY
!    The parent module cgca_m3sld.
!  SOURCE

contains

!*roboend*


!*robodoc*s* m3sld_hc/cgca_sld_h
!  NAME
!    cgca_sld_h
!  SYNOPSIS
!    module subroutine cgca_sld_h( coarray, hx, iter, heartbeat, solid )

module procedure cgca_sld_h

!  INPUTS
!    See the parent module cgca_m3sld.
!  OUTPUT
!    See the parent module cgca_m3sld.
!  SIDE EFFECTS
!    State of coarray changed
!  DESCRIPTION
!    This is a simplified version of cgca_sld.
!    Most checks have been removed.
!    We use co_sum reduction.
!    The desired halo exchange routine is passed as input.
!    In addition, it does not support the periodic BC.
!
!    Inputs:
!    - coarray - the model
!    - hx - the halo exchange procedure
!    - iter - number of solidification iterations, if <=0 then do
!      until the coarray has solidified; if >0 then proceed until
!      solid or "iter" iterations have been completed, whichever is sooner
!    - heartbeat - if >0 then dump a simple message every heartbeat
!      iterations
!
!    Outputs:
!    - solid - .true. if the coarray is fully solid, .false. otherwise
!
!    At least one value .ne. cgca_liquid_state must exist on at least
!    one image.
!    At least one cgca_liquid_state value must exist on at least one image.
!  NOTES
!    All images must call this routine!
!  USES
!  USED BY
!    none, end user
!  SOURCE

real :: candidate(3)

integer( kind=iarr ), allocatable :: array(:,:,:)
integer :: errstat
integer( kind=idef ) :: &
  lbv(4)      ,& ! lower bounds of the complete (plus virtual) coarray
  ubv(4)      ,& ! upper bounds of the complete (plus virtual) coarray
  lbr(4)      ,& ! lower bounds of the "real" coarray, lower virtual+1
  ubr(4)      ,& ! upper bounds of the "real" coarray, upper virtual-1
  img,nimages ,& ! to avoid repeated calls to this_image() and num_images()
  x1,x2,x3    ,& ! local coordinates in an array, which are also
  step(3)     ,& ! a random neighbouring cell, out of 26 possibilities
  iteration      ! solidification iteration

! true if the local array has solidified
logical( kind=ldef ) :: finished

character( len=100 ) :: e_message

! finished indicator (findicator) is an integer variable
! to use with CO_SUM.
! It is set to 1 on all images at the beginning of every loop.
! if (finished) findicator = 0.
! Then if co_sum(findicator) is zero, then all images have finished.
integer( kind=idef ), save :: findicator[*]

integer :: flag

! get image number and number of images
    img = this_image()
nimages = num_images()

!*************************************************
! Sanity checks
!*************************************************

! check for allocated
if ( .not. allocated(coarray) ) then
  write (*,*) "ERROR: cgca_sld_h/m3sld_hc: coarray is not allocated"
  error stop
end if

!*************************************************
! End of sanity checks. All seems fine, proceed.
!*************************************************

! determine the extents
lbv = lbound( coarray )
ubv = ubound( coarray )
lbr = lbv + 1
ubr = ubv - 1

! Mark as not solid initially.
solid = .false.

! initialise the iteration counter
iteration = 1

! allocate the temp array
allocate( array( lbv(1):ubv(1) , lbv(2):ubv(2) , lbv(3):ubv(3) ),      &
          source = 0_iarr,                                             &
          stat = errstat, errmsg = e_message )
if ( errstat .ne. 0 ) then
  write (*,'(a,i0,a)')                                                 &
    "ERROR: cgca_sld_h/m3sld_hc: allocate( array ), err stat: ",       &
     errstat, " err message: " // e_message
  error stop
end if

! start the main loop
main: do

 ! copy coarray, grain state type, into a local array
 array = coarray( : , : , : , cgca_state_type_grain )

 ! solidify array
 do x3 = lbr(3),ubr(3)
 do x2 = lbr(2),ubr(2)
 do x1 = lbr(1),ubr(1)
  if (coarray(x1,x2,x3,cgca_state_type_grain) .eq. cgca_liquid_state)  &
    then
      call random_number( candidate ) ! [0..1)
      step = nint( candidate*2-1 )    ! step = [-1 0 1]
      array( x1, x2, x3 ) =                                            &
                     coarray( x1+step(1), x2+step(2), x3+step(3),      &
                                               cgca_state_type_grain )
  end if
 end do
 end do
 end do

 ! update coarray
 coarray( :, :, :, cgca_state_type_grain ) = array

 ! image 1 sends heartbeat signal
 if ( img .eq. 1 ) then
  if ( heartbeat .gt. 0 ) then
   if ( mod(iteration,heartbeat) .eq. 0 ) write (*,'(a,i0)')           &
    "INFO: cgca_sld_h/m3sld_hc: iterations completed: ", iteration
  end if
 end if

 ! increment the iteration counter
 iteration = iteration + 1

 ! exit if done the required number of iterations 
 if (iter .gt. 0 .and. iteration .gt. iter) exit main

 ! set finish to .true. if finished (all solid)
 finished = all( coarray(lbr(1):ubr(1), lbr(2):ubr(2), lbr(3):ubr(3),  &
                  cgca_state_type_grain) .ne. cgca_liquid_state)
 findicator = 1_idef
 if (finished) findicator = 0_idef

 sync all

 ! halo exchange in preparation for the next iteration
 ! no sync inside
 call hx( coarray ) 

 sync all

! Check hx, flag .eq. 0 means ok, flag .ne. 0 is an error. 
call cgca_hxic( coarray, flag )
call co_sum( flag )
if ( img .eq. 1 ) then
  if ( flag .ne. 0 ) then
    write (*,'(a,i0)') "ERROR: cgca_sld_h/m3sld_hc: hx flag: ", flag
    error stop
  end if
end if

 ! do the collective AND on finished
 call co_sum( findicator )

 ! Now all images will have the updated findicator.
 ! Exit if finished
 if ( findicator .eq. 0 ) then
  solid = .true.
  exit main
 end if

end do main

! deallocate all local arrays

deallocate( array, stat=errstat, errmsg = e_message )
if ( errstat .ne. 0 ) then
  write (*,'(a,i0,a)')                                                 &
    "ERROR: cgca_sld_h/m3sld_hc: deallocate( array ), err stat: ",     &
     errstat, " err message: " // e_message
  error stop
end if

end procedure cgca_sld_h

!*roboend*

end submodule m3sld_hc
