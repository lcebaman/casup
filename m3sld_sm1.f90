!*robodoc*f* cgca_m3sld/m3sld_sm1
!  NAME
!    m3sld_sm1
!  SYNOPSIS

!$Id: m3sld_sm1.f90 400 2017-05-04 17:47:56Z mexas $

submodule ( cgca_m3sld ) m3sld_sm1

!  DESCRIPTION
!    Submodule with solidification routines using collectives.
!    These are not supported by ifort 16 yet.
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    cgca_sld3
!  USES
!    All variables and parameters of the parent module cgca_m3sld.
!  USED BY
!    The parent module cgca_m3sld.
!  SOURCE

contains

!*roboend*


!*robodoc*s* m3sld_sm1/cgca_sld3
!  NAME
!    cgca_sld3
!  SYNOPSIS

!module procedure cgca_sld3
  module subroutine cgca_sld3( coarray, iter, heartbeat, solid )
    integer( kind=iarr ), allocatable, intent( inout ) ::              &
      coarray(:,:,:,:)[:,:,:]
    integer( kind=idef ), intent( in ) :: iter,heartbeat
    logical( kind=ldef ), intent( out ) :: solid


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
!    In addition, it does not support the periodic BC.
!
!    Inputs:
!    - coarray - the model
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

! finished indicator (findicator) is an integer variable to use with CO_SUM.
! It is set to 1 on all images at the beginning of every loop.
! if (finished) findicator = 0.
! Then if co_sum(findicator) is zero, then all images have finished.
integer( kind=idef ),save :: findicator[*]

! get image number and number of images
    img = this_image()
nimages = num_images()

!*************************************************
! Sanity checks
!*************************************************

! check for allocated
if ( .not. allocated(coarray) ) &
  error stop "ERROR: cgca_sld3: coarray is not allocated"

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
    "ERROR: cgca_sld3/m3sld_sm1: allocate( array ), err stat: ",       &
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
  if (coarray(x1,x2,x3,cgca_state_type_grain) .eq. cgca_liquid_state) then
   call random_number(candidate)     ! 0 .le. candidate .lt. 1
   step = nint(candidate*2-1)        ! step = [-1 0 1]
   array (x1,x2,x3) = &
    coarray (x1+step(1),x2+step(2),x3+step(3),cgca_state_type_grain)
  end if
 end do
 end do
 end do

 ! update coarray
 coarray(:,:,:,cgca_state_type_grain) = array

 ! image 1 sends heartbeat signal
 if ( img .eq. 1 ) then
  if ( heartbeat .gt. 0 ) then
   if ( mod(iteration,heartbeat) .eq. 0 ) write (*,'(a,i0)') &
    "INFO: cgca_sld3: iterations completed: ", iteration
  end if
 end if

 ! increment the iteration counter
 iteration = iteration + 1

 ! exit if done the required number of iterations 
 if (iter .gt. 0 .and. iteration .ge. iter) exit main

 ! set finish to .true. if finished (all solid)
 finished = all( coarray(lbr(1):ubr(1), lbr(2):ubr(2), lbr(3):ubr(3), &
                  cgca_state_type_grain) .ne. cgca_liquid_state)
 findicator = 1_idef
 if (finished) findicator = 0_idef


 ! halo exchange in preparation for the next iteration
 call cgca_hxi(coarray) 

 ! do the collective AND on finished
 call co_sum(findicator)

 ! Now all images will have the updated findicator.
 ! Exit if finished
 if ( findicator .eq. 0 ) then
  solid = .true.
  exit main
 end if

end do main

! deallocate all local arrays

deallocate(array,stat=errstat)
if (errstat.ne.0) then
  write (*,'(a,i0)') "ERROR: cgca_sld3: image ", img
  write (*,'(a)') "ERROR: cgca_sld3: cannot deallocate array"
  error stop
end if

!end procedure cgca_sld3
end subroutine cgca_sld3

!*roboend*

end submodule m3sld_sm1
