!*robodoc*m* CGPACK/cgca_m3nucl
!  NAME
!    cgca_m3nucl
!  SYNOPSIS

!$Id: cgca_m3nucl.f90 380 2017-03-22 11:03:09Z mexas $

module cgca_m3nucl

!  DESCRIPTION
!    Module dealing with nucleation
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    cgca_nr
!  USES
!    cgca_m2glm
!  USED BY
!    cgca
!  SOURCE

use cgca_m1co
use cgca_m2glm

implicit none

private
public :: cgca_nr

contains

!*roboend*


!*robodoc*s* cgca_m3nucl/cgca_nr
!  NAME
!    cgca_nr
!  SYNOPSIS

subroutine cgca_nr( coarray , number , debug )

!  INPUTS

integer( kind=iarr ), allocatable, intent(inout) ::                    &
 coarray(:,:,:,:)[:,:,:]
integer( kind=idef ), intent( in ) :: number
logical( kind=ldef ), intent( in ) :: debug

!  SIDE EFFECTS
!    State of coarray changed
!  DESCRIPTION
!    This routine randomly scatters the given number of grain
!    nuclei over the model. The grain nuclei are assigned unique
!    numbers starting from zero.
!
!    All elements of the coarray must be in cgca_liquid_state state.
!    If not, the program will stop with a error.
!
!    The number of nuclei must not greater than "critfract" of
!    the model size. I arbitrarily set this to 0.1 for now.
!    However, even this is too high. Although, in principle
!    each cell can be a nuclei, such model would have no physical
!    sense.
!
!    Inputs:
!    - coarray - the model
!    - number - number of nuclei to scatter
!    - debug - if .true. dump some debug output
!  NOTES
!    All images must call this routine!
!    However, the work will be done only by image 1.
!    There's "sync all" at the end of this routine.
!  USES
!    cgca_gl
!  USED BY
!    none, end user
!  SOURCE

real( kind=rdef), parameter :: critfract=0.1

integer( kind=idef ) :: &
  nuc    ,& ! running total of nuclei
  lbr(4) ,& ! lower bounds of the "real" coarray, lower virtual+1
  ubr(4) ,& ! upper bounds of the "real" coarray, upper virtual-1
  szr(3) ,& ! size or the "real" coarray, ubr-lbr+1
  lcob(3),& ! lower cobounds of the coarray
  ucob(3),& ! upper cobounds of the coarray
  supermax(3) ,& ! upper bound of the super array, szr*(ucob-lcob+1)
  supermin(3) ,& ! lower bound of the super array, always 1.
  super(3)    ,& ! coordinates in a super array
  imgpos(3)   ,& ! image position in a grid
  local(3)    ,& ! coordinates within an image
  thisimage   ,& ! this_image()
  nimages        ! num_images()

integer( kind=ilrg ) :: coarsize

real( kind=rdef ) :: candidate(3), frac
logical( kind=ldef ) :: image1

  nimages = num_images()
thisimage = this_image()
   image1 = .false.
if ( thisimage .eq. 1 ) image1 = .true.

!*********************************************************************72
! checks
!*********************************************************************72

if ( .not. allocated( coarray ) ) then
 write( *, '(a,i0)' ) "ERROR: cgca_nr/cgca_m3nucl: coarray is not" //  &
   " allocated, img: ", thisimage
 error stop 
end if

! check that there are only liquid cells in coarray.

lbr = lbound( coarray ) + 1
ubr = ubound( coarray ) - 1

if ( any( coarray(lbr(1):ubr(1),lbr(2):ubr(2),lbr(3):ubr(3),           &
              cgca_state_type_grain) .ne. cgca_liquid_state)) then
 write( *, '(a,i0)' ) "ERROR: cgca_nr/cgca_m3nucl: non-liquid" //      &
   " elements in coarray, img: ", thisimage
 error stop
end if

! check that the number of nuclei is positive

if ( number .lt. 1 ) then
 write (*,'(a,i0)') "ERROR: cgca_nr/cgca_m3nucl: number of nuclei" //  &
   " must be 1 or more, img:", thisimage
 error stop
end if

!*********************************************************************72
! end of checks
!*********************************************************************72

! image 1 must not change values in other images before
! all images pass checks
sync all

img1: if ( image1 ) then

 ! Warn the user if there are too many nuclei
 coarsize = size( coarray( lbr(1) : ubr(1) , lbr(2) : ubr(2) ,         &
                           lbr(3) : ubr(3), cgca_state_type_grain ) ,  &
                  kind=ilrg )

 ! number of grains as a fraction of the model
 frac = number / ( real( nimages ) * real( coarsize ) )

 if ( frac .gt. critfract ) then
   write (*,'(a,g10.3)') "WARN: cgca_nr/cgca_m3nucl: too many " //     &
    "nuclei - no physical sense! nuclei/model size: ", frac
 end if
    
 ! The 4th dimension is the number of cell state types.
 ! It is not relevant here, so don't use it.
      szr = ubr(1:3)-lbr(1:3)+1
     lcob = lcobound(coarray)
     ucob = ucobound(coarray)
 supermax = szr * (ucob-lcob+1)
 supermin = 1
 
 nuc=1
 do
  call random_number(candidate)    ! 0 .le. candidate .lt. 1
  super=int(candidate*supermax)+1  ! 1 .le. super .le. supermax
  
  ! now translate to the image and local coordinates
  call cgca_gl( super, coarray, imgpos, local )
  
  ! If a cell is liquid then assign the running "nuc" number to it.
  ncln: if ( coarray( local(1), local(2), local(3),                    &
             cgca_state_type_grain ) [imgpos(1), imgpos(2), imgpos(3)] &
                .eq. cgca_liquid_state ) then
          coarray( local(1), local(2), local(3),                       &
          cgca_state_type_grain ) [imgpos(1),imgpos(2),imgpos(3)] = nuc
  
    ! If requested, dump some debug output
    if ( debug ) then
      write( *, "(2(a,3(i0,tr1)),a,i0)" ) "DEBUG:" //                  &
        " cgca_nr/cgca_m3nucl: local: ", local, "imgpos: ", imgpos,    &
        " nucleus: ", nuc
    end if
  
    ! Increment the running total of the nuclei generated
    nuc = nuc + 1
  end if ncln
  
  ! If "number" of nuclei have been generated, exit
  if ( nuc .gt. number ) exit
  
 end do

end if img1

! Global sync is required here
sync all

end subroutine cgca_nr

!*roboend*

end module cgca_m3nucl
