!*robodoc*m* CGPACK/cgca_m3sld
!  NAME
!    cgca_m3sld
!  SYNOPSIS

!$Id: cgca_m3sld.f90 431 2017-06-30 13:13:49Z mexas $

module cgca_m3sld

!  DESCRIPTION
!    Module dealing with solidification
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    Subroutines: cgca_sld, cgca_sld1, cgca_sld2,
!    cgca_sld3 (in submodule m3sld_sm1),
!    cgca_sld_h (in submodule m3sld_hc).
!  USES
!    cgca_m1co, cgca_m2hx, cgca_m2red
!  USED BY
!    cgca
!  SOURCE

use cgca_m1co
use cgca_m2hx
use cgca_m2red

implicit none

private
public :: cgca_sld, cgca_sld1, cgca_sld2, &
          cgca_sld3    , & ! in submodule m3sld_sm1 
          cgca_sld_h       ! in submodule m3sld_hc

abstract interface
  subroutine halo_exchange( array )
    use cgca_m1co
    integer( kind=iarr ), allocatable, intent( inout ) ::              &
                                                array(:,:,:,:)[:,:,:]
  end subroutine halo_exchange
end interface

interface
  ! In submodule m3sld_sm1
  module subroutine cgca_sld3( coarray, iter, heartbeat, solid )
    integer( kind=iarr ), allocatable, intent( inout ) ::              &
      coarray(:,:,:,:)[:,:,:]
    integer( kind=idef ), intent( in ) :: iter,heartbeat
    logical( kind=ldef ), intent( out ) :: solid
  end subroutine cgca_sld3

  ! In submodule m3sld_hc
  module subroutine cgca_sld_h( coarray, hx, iter, heartbeat, solid )
    integer( kind=iarr ), allocatable, intent( inout ) ::              &
      coarray(:,:,:,:)[:,:,:]
    procedure( halo_exchange ) :: hx
    integer( kind=idef ), intent( in ) :: iter, heartbeat
    logical( kind=ldef ), intent( out ) :: solid
  end subroutine cgca_sld_h
end interface

contains

!*roboend*


!*robodoc*s* cgca_m3sld/cgca_sld
!  NAME
!    cgca_sld
!  SYNOPSIS

subroutine cgca_sld( coarray, periodicbc, iter, heartbeat, solid )

!  INPUTS

integer( kind=iarr ), allocatable, intent( inout ) ::                  &
  coarray(:,:,:,:)[:,:,:]
integer( kind=idef ), intent( in ) :: iter, heartbeat
logical( kind=ldef ), intent( in ) :: periodicbc

!  OUTPUT

logical( kind=ldef ), intent( out ) :: solid

!  SIDE EFFECTS
!    None
!  DESCRIPTION
!    This routine scans over all cgca_liquid_state cells and
!    gives them a chance to attach to a grain.
!    Thus the grains grow and the liquid phase decreases.
!    The routine can be run for a number of iterations, or
!    until the whole model has solidified, i.e. no more
!    liquid phase left.
!
!    Inputs:
!    - coarray - the model
!    - periodicbc - if .true. periodic boundary conditions are used,
!      i.e. global halo exchange is called after every iteration
!    - iter - number of solidification iterations, if <=0 then do
!      until the coarray has solidified; if >0 then proceed until
!      solid or "iter" iterations have been completed, whichever is sooner
!    - heartbeat - if >0 then dump a simple message every heartbeat
!      iterations
!
!    Outputs:
!    - solid - .true. if the coarray is fully solid, .false. otherwise
!
!    At least one value .gt. cgca_liquid_state must exist on at least
!    one image.
!    At least one cgca_liquid_state value must exist on at least one image.
!  NOTES
!    All images must call this routine!
!    The grain numbers are *always* greater than the liquid state.
!  USES
!  USED BY
!    none, end user
!  SOURCE

real :: candidate(3)

integer( kind=iarr ), allocatable :: array( : , : , : )
integer :: i, errstat
integer( kind=idef ) :: &
  lbv(4),      & ! lower bounds of the complete (plus virtual) coarray
  ubv(4),      & ! upper bounds of the complete (plus virtual) coarray
  lbr(4),      & ! lower bounds of the "real" coarray, lower virtual+1
  ubr(4),      & ! upper bounds of the "real" coarray, upper virtual-1
  thisimage   ,& ! to avoid repeated calls to this_image()
  numimages   ,& ! to avoid repeated calls to num_images()
  x1          ,& ! local coordinates in an array, which are also
  x2          ,& ! do loop counters
  x3          ,& !
  step(3)     ,& ! a random neighbouring cell, out of 26 possibilities
  iteration      ! solidification iteration

! Must be SAVEd to conform
logical( kind=ldef ), save ::                                          &
  yesnuclei [*],   & ! true if nuclei value exists on any image
  finished [*],    & ! true if solidified on this image
  allfinished [*]    ! true if solidified on all images

! .true. if a grain value in coarray exists at least on one image
logical( kind=ldef ) :: allyesnuc

character( len=100 ) :: e_message

! "finished" is calculated by each image and then the value
! is analysed by image 1 to calculate "allfinished".
! Then all other images read "allfinished" from image 1 for
! synchronisation.
!
! Important: when checking for finished, use only the real
! parts of coarray. Do not analyse the virtual (halo) elements!

! get image number and number of images
thisimage = this_image()
numimages = num_images()

!*************************************************
! Sanity checks
!*************************************************

! check for allocated
if ( .not. allocated( coarray ) ) &
  error stop "ERROR: cgca_sld: coarray is not allocated"

! determine the extents
lbv = lbound( coarray )
ubv = ubound( coarray )
lbr = lbv + 1
ubr = ubv - 1

! Check for at least one nuclei on this image.
yesnuclei = .false.
if ( any( coarray(lbr(1):ubr(1), lbr(2):ubr(2), lbr(3):ubr(3),         &
                  cgca_state_type_grain) .gt. cgca_liquid_state))      &
  yesnuclei=.true.

! Check if all cells on this image have solidified.
finished = .false.
if ( all( coarray(lbr(1):ubr(1), lbr(2):ubr(2), lbr(3):ubr(3),         &
                  cgca_state_type_grain) .gt. cgca_liquid_state))      &
  finished=.true.

         ! yesnuclei and finished defined
sync all ! new execution segment
         ! yesnuclei and finished used

! Image 1 calculates the global values.
if ( thisimage .eq. 1 ) then

  ! Need at least one nuclei on any image.
  allyesnuc = .false.
  do i = 1 , numimages
    if ( yesnuclei[i] ) then
      allyesnuc = .true.
      exit 
    end if
  end do

  ! Need at least one non solidified cell on any image
  allfinished = .true. 
  do i = 1 , numimages
    if ( .not. finished[i] ) then
      allfinished=.false.
      exit 
    end if
  end do

  ! Report an error if no nuclei found
  if ( .not. allyesnuc ) then
    write (*,'(a)') "ERROR: cgca_sld/cgca_m3sld: no nuclei found"
    error stop
  end if

  ! Report an error if no liquid phase found
  if ( allfinished ) then
    write (*,'(a)') "ERROR: cgca_sld/cgca_m3sld: already solid"
    error stop
  end if
    
end if

! All images wait for image 1 here.
sync all

!*************************************************
! End of sanity checks. All seems fine, proceed.
!*************************************************

! Mark as not solid initially
solid = .false.

! allocate the temp array
allocate( array( lbv(1) : ubv(1) , lbv(2) : ubv(2) , lbv(3) : ubv(3) ),&
          stat = errstat,  errmsg = e_message )
if ( errstat .ne. 0 ) then
  write (*,'(a,i0,a)')                                                 &
    "ERROR: cgca_sld/cgca_m3sld: allocate( array ), err stat: ",       &
     errstat, " err message: " // e_message
  error stop
end if

! initialise the iteration counter
iteration = 1

! at this point allfinished=.false. on image 1 and not
! initialised on all other images!

! start the main loop
main: do

  ! halo exchange,        ===>>> remote comms <<<===
  call cgca_hxi( coarray )
  sync all

  ! global halo exchange, ===>>> remote comms <<<===
  if ( periodicbc ) call cgca_hxg( coarray )
  sync all

  ! do the iteration if not finished
  fini: if ( .not. finished ) then

    ! copy coarray, grain state type, into a local array
    array = coarray( : , : , : , cgca_state_type_grain )

    ! solidify array
    do x3 = lbr(3), ubr(3)
    do x2 = lbr(2), ubr(2)
    do x1 = lbr(1), ubr(1)
      if ( coarray( x1, x2, x3, cgca_state_type_grain ) .eq.           &
                                          cgca_liquid_state ) then
        call random_number( candidate )     ! 0 .le. candidate .lt. 1
        step = nint( candidate*2 - 1 )      ! step = [-1 0 1]
        array( x1, x2, x3 ) = coarray( x1 + step(1) ,                  &
                                       x2 + step(2) ,                  &
                                       x3 + step(3) ,                  &
                                            cgca_state_type_grain )
      end if
    end do
    end do
    end do

    ! update coarray on this image
    coarray( : , : , : , cgca_state_type_grain ) = array

    ! see if finished (all solid)
    finished = .true.
    if ( any( coarray( lbr(1):ubr(1) , lbr(2):ubr(2) , lbr(3):ubr(3),  &
          cgca_state_type_grain ) .eq. cgca_liquid_state ) ) then
      finished = .false.
    end if

  end if fini

  ! Global sync here. All images must calculate "finished"
  ! before image 1 reads these values from all images
  ! and calculates the global "allfinished".
  sync all 

  ! image 1 checks if finished and sends heartbeat signal
  img1: if ( thisimage .eq. 1 ) then

    ! assume we are done
    allfinished = .true.
    do i = 1 , numimages

      ! but if any image is not done, then we are not done either
      if ( .not. finished[i] ) then 
        allfinished = .false.
        exit
      end if
    end do

    ! send heartbeat signal to terminal
    if ( heartbeat .gt. 0 ) then
      if ( mod( iteration, heartbeat ) .eq. 0 )                        &
         write (*,'(a,i0)') "INFO: cgca_sld: iterations completed: ",  &
                            iteration
    end if

  end if img1

            ! allfinished calculated on image 1
  sync all  ! new segment
            ! allfinished used on all images

  ! get allfinished from image 1
  allfinished = allfinished [ 1 ]

  ! exit if done
  if ( allfinished ) then
    solid = .true.
    exit main
  end if

  ! Also exit is the max number of iterations has been reached
  if ( iter .gt. 0 .and. iteration .ge. iter ) exit main

  ! increment the iteration counter
  iteration = iteration + 1

end do main

deallocate( array, stat=errstat )
if ( errstat .ne. 0 ) then
  write( *, '(a,i0)' ) "ERROR: cgca_sld/cgca_m3sld:" //                &
    " deallocate( array ), img: ", thisimage
  error stop
end if

end subroutine cgca_sld

!*roboend*


!*robodoc*s* cgca_m3sld/cgca_sld1
!  NAME
!    cgca_sld1
!  SYNOPSIS

subroutine cgca_sld1(coarray,iter,heartbeat,solid)

!  INPUTS

integer(kind=iarr),allocatable,intent(inout) :: coarray(:,:,:,:)[:,:,:]
integer(kind=idef),intent(in) :: iter,heartbeat

!  OUTPUT

logical(kind=ldef),intent(out) :: solid

!  SIDE EFFECTS
!    State of coarray changed
!  DESCRIPTION
!    This is a simplified version of cgca_sld.
!    Most checks have been removed and sync instances
!    reduced. In addition, it does not support the periodic BC.
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

integer(kind=iarr),allocatable :: array(:,:,:)
integer :: errstat
integer(kind=idef) :: &
  lbv(4)      ,& ! lower bounds of the complete (plus virtual) coarray
  ubv(4)      ,& ! upper bounds of the complete (plus virtual) coarray
  lbr(4)      ,& ! lower bounds of the "real" coarray, lower virtual+1
  ubr(4)      ,& ! upper bounds of the "real" coarray, upper virtual-1
  thisimage   ,& ! to avoid repeated calls to this_image() and num_images()
  nimages     ,& !
  x1,x2,x3    ,& ! local coordinates in an array, which are also
  step(3)     ,& ! a random neighbouring cell, out of 26 possibilities
  iteration      ! solidification iteration

! true if the local array has solidified
logical(kind=ldef) :: finished

! number of finished images
! NOTE: if *really* many images are used, then the kind will need to
! be increased!
integer(kind=idef),allocatable :: nfini[:]

! finished is set to .true. on all images at the beginning of each
! iteration. At the end of the iteration finished[1] is set to .false.
! by *all* images if their local finished is .false..
!
! Important: when checking for finished, use only the real
! parts of coarray. Do not analyse the virtual (halo) elements!

! get image number and number of images
thisimage = this_image()
nimages = num_images()

!*************************************************
! Sanity checks
!*************************************************

! check for allocated
if (.not. allocated(coarray)) &
  error stop "ERROR: cgca_sld1: coarray is not allocated"

!*************************************************
! End of sanity checks. All seems fine, proceed.
!*************************************************

! determine the extents
lbv=lbound(coarray)
ubv=ubound(coarray)
lbr=lbv+1
ubr=ubv-1

! allocate the integer array to store the number of finished images
allocate(nfini[*],stat=errstat)
if (errstat.ne.0) then
 write (*,'(a,i0)') "ERROR: cgca_sld1: image ",thisimage
 write (*,'(a)') "ERROR: cgca_sld1: cannot allocate nfini"
 error stop
end if

! allocate the temp array
allocate(array(lbv(1):ubv(1),lbv(2):ubv(2),lbv(3):ubv(3)),stat=errstat)
if (errstat.ne.0) then
  write (*,'(a,i0)') "ERROR: cgca_sld1: image ",thisimage
  write (*,'(a)') "ERROR: cgca_sld1: cannot allocate array"
  error stop
end if

! Mark as not solid initially.
solid = .false.

! set finished to .false.
finished = .false.

! initialise the iteration counter
iteration=1

! start the main loop
main: do

! set the number of finished images to zero.
nfini = 0

! do if not finished
fini: if (.not. finished) then

! copy coarray, grain state type, into a local array
array=coarray(:,:,:,cgca_state_type_grain)

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

! set finish to .true. if finished (all solid)
finished = all( coarray(lbr(1):ubr(1), lbr(2):ubr(2), lbr(3):ubr(3), &
                  cgca_state_type_grain) .ne. cgca_liquid_state)

end if fini

! Each image, but first, waits for the lower image. Then it sets
! finished[1] to .false. if not finished.
if (thisimage .eq. 1) then
 if ( finished ) nfini = nfini + 1
! write (*,*) "image", thisimage, finished, nfini
else
 sync images (thisimage-1)
 if ( finished ) nfini[1] = nfini[1] + 1
! write (*,*) "image", thisimage, finished, nfini[1]
end if

! each image, but last, waits for the next one
if (thisimage .lt. nimages) sync images (thisimage+1)

! halo exchange
call cgca_hxi(coarray) 

! exit if done the required number of iterations 
if (iter .gt. 0 .and. iteration .ge. iter) exit main

! increment the iteration counter
iteration = iteration+1

! image 1 sends heartbeat signal
if (thisimage .eq. 1) then
 if (heartbeat .gt. 0) then
   if (mod(iteration,heartbeat) .eq. 0) write (*,'(a,i0)') &
    "INFO: cgca_sld1: iterations completed: ", &
      iteration
 end if
end if

sync all

! exit if finished
if ( nfini[1] .eq. nimages ) then
 solid = .true.
 exit main
end if

end do main

! deallocate all local arrays

deallocate(array,stat=errstat)
if (errstat.ne.0) then
  write (*,'(a,i0)') "ERROR: cgca_sld1: image ",thisimage
  write (*,'(a)') "ERROR: cgca_sld1: cannot deallocate array"
  error stop
end if

deallocate(nfini,stat=errstat)
if (errstat.ne.0) then
  write (*,'(a,i0)') "ERROR: cgca_sld1: image ",thisimage
  write (*,'(a)') "ERROR: cgca_sld1: cannot deallocate nfini"
  error stop
end if

end subroutine cgca_sld1

!*roboend*


!*robodoc*s* cgca_m3sld/cgca_sld2
!  NAME
!    cgca_sld2
!  SYNOPSIS

subroutine cgca_sld2(coarray,p,iter,heartbeat,solid)

!  INPUTS

integer(kind=iarr),allocatable,intent(inout) :: coarray(:,:,:,:)[:,:,:]
integer(kind=idef),intent(in) :: p,iter,heartbeat

!  OUTPUT

logical(kind=ldef),intent(out) :: solid

!  SIDE EFFECTS
!    State of coarray changed
!  DESCRIPTION
!    This is a simplified version of cgca_sld.
!    Most checks have been removed. cgca_redand is called to check that
!    all images have solidified.
!    In addition, it does not support the periodic BC.
!
!    Inputs:
!    - coarray - the model
!    - p - this routine only works when the number of images is a power
!      of 2. So p is the power: num_images = 2**p. Note that no check
!      for this is made in this routine. This is left up to the calling
!      routine, most probably the main program.
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
!    cgca_redand
!  USED BY
!    none, end user
!  SOURCE

real :: candidate(3)

integer(kind=iarr),allocatable :: array(:,:,:)
integer :: errstat
integer(kind=idef) :: &
  lbv(4)      ,& ! lower bounds of the complete (plus virtual) coarray
  ubv(4)      ,& ! upper bounds of the complete (plus virtual) coarray
  lbr(4)      ,& ! lower bounds of the "real" coarray, lower virtual+1
  ubr(4)      ,& ! upper bounds of the "real" coarray, upper virtual-1
  img,nimages ,& ! to avoid repeated calls to this_image() and num_images()
  x1,x2,x3    ,& ! local coordinates in an array, which are also
  step(3)     ,& ! a random neighbouring cell, out of 26 possibilities
  iteration      ! solidification iteration

! true if the local array has solidified
! It must have SAVE attribute or be ALLOCATABLE. Right now SAVE seems
! faster. The memory freed after this routine, if ALLOCATABLE was used
! instead is insignificant.
logical(kind=ldef),save :: finished[*]

! "finished" is set to .false. on all images at the beginning of
! each iteration. At the end of an iteration, new "finished" value
! is calculated on each image. Then cgca_redall is called. It places
! the result in "finished" on every image. So after that every image
! needs to check just its local "finished". If it is .true., then exit
! the loop. If it is .false., then do another solidification iteration.

! get image number and number of images
img     = this_image()
nimages = num_images()

!*************************************************
! Sanity checks
!*************************************************

! check for allocated
if (.not. allocated(coarray)) &
  error stop "ERROR: cgca_sld2: coarray is not allocated"

!*************************************************
! End of sanity checks. All seems fine, proceed.
!*************************************************

! determine the extents
lbv=lbound(coarray)
ubv=ubound(coarray)
lbr=lbv+1
ubr=ubv-1

! Mark as not solid initially.
solid = .false.

! initialise the iteration counter
iteration=1

! allocate the temp array
! Implicit sync all here
allocate(array(lbv(1):ubv(1),lbv(2):ubv(2),lbv(3):ubv(3)),stat=errstat)
if (errstat.ne.0) then
  write (*,'(a,i0)') "ERROR: cgca_sld2: image ", img
  write (*,'(a)') "ERROR: cgca_sld2: cannot allocate array"
  error stop
end if

! start the main loop
main: do

 ! copy coarray, grain state type, into a local array
 array = coarray(:,:,:,cgca_state_type_grain)

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
    "INFO: cgca_sld2: iterations completed: ", iteration
  end if
 end if

 ! increment the iteration counter
 iteration = iteration+1

 ! exit if done the required number of iterations 
 if (iter .gt. 0 .and. iteration .ge. iter) exit main

 ! set finish to .true. if finished (all solid)
 finished = all( coarray(lbr(1):ubr(1), lbr(2):ubr(2), lbr(3):ubr(3), &
                  cgca_state_type_grain) .ne. cgca_liquid_state)

 ! halo exchange in preparation for the next iteration
 call cgca_hxi(coarray) 

! not sure if I need a global barrier here or not
 sync all

 ! do the collective AND on finished
 call cgca_redand(finished,p)

 ! Now all images will have the updated "finished".
 ! Exit if finished
 if ( finished ) then
  solid = .true.
  exit main
 end if

end do main

! deallocate all local arrays

deallocate(array,stat=errstat)
if (errstat.ne.0) then
  write (*,'(a,i0)') "ERROR: cgca_sld2: image ", img
  write (*,'(a)') "ERROR: cgca_sld2: cannot deallocate array"
  error stop
end if

end subroutine cgca_sld2

!*roboend*

end module cgca_m3sld
