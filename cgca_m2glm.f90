!*robodoc*m* CGPACK/cgca_m2glm
!  NAME
!    cgca_m2glm
!  SYNOPSIS

!$Id: cgca_m2glm.f90 379 2017-03-22 09:57:10Z mexas $

module cgca_m2glm

!  DESCRIPTION
!    Module dealing with Global to Local Mapping (glm) and vice versa
!  AUTHOR 
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    cgca_gl, cgca_lg, cgca_ico
!  USES
!    cgca_m1co
!  USED BY
!    cgca_m3nucl
!  SOURCE

use cgca_m1co
implicit none

private
public :: cgca_gl, cgca_lg, cgca_ico, cgca_ico2

contains

!*roboend*


!*robodoc*s* cgca_m2glm/cgca_gl
!  NAME
!    cgca_gl
!  SYNOPSIS

subroutine cgca_gl(super,coarray,imgpos,local)

!  INPUTS

integer(kind=idef),intent(in) :: super(3)
integer(kind=iarr),allocatable,intent(inout) :: coarray(:,:,:,:)[:,:,:]

! OUTPUT

integer(kind=idef),intent(out) :: imgpos(3),local(3)

! DESCRIPTION
!   This routine converts a cell coordinate from a global, super, array
!   to the image coordinates in the coarray grid and the local cell
!   coordinates in this image :
!   - super(3) are cell coordinates in a super array
!   - coarray is the model
!   - imgpos(3) is the image position in the grid
!   - local(3) are cell coordinates in that image's array
! NOTES
!   The global coordinates must start from 1!
!
!   Any image can call this routine
! USES
! USED BY
! SOURCE

integer :: & 
 lbr(4)   ,& ! lower bounds of the "real" coarray, lbv+1
 ubr(4)   ,& ! upper bounds of the "real" coarray, ubv-1
 szr(3)   ,& ! size or the "real" coarray, ubr-lbr+1
 lcob(3)  ,& ! lower cobounds of the coarray
 ucob(3)  ,& ! upper cobounds of the coarray
 usup(3)  ,& ! upper bound of the super array, szr*(ucob-lcob+1)
 thisimage

thisimage = this_image()

! check for allocated

if (.not. allocated(coarray)) then
 write (*,'(a,i0)') "ERROR: cgca_gl: image", thisimage
 write (*,'(a)') "ERROR: cgca_gl: coarray is not allocated"
 error stop
end if

lbr=lbound(coarray)+1
ubr=ubound(coarray)-1

! the 4th dimension is to do with the number of cell state
! types. This is not relevant here.
szr=ubr(1:3)-lbr(1:3)+1

lcob=lcobound(coarray)
ucob=ucobound(coarray)
usup=szr*(ucob-lcob+1)

! check for bounds

if (any(super .gt. usup) .or. any(super .lt. 1)) then
 write (*,'(a,i0)') "ERROR: cgca_gl: image", thisimage
 write (*,'(a)') "ERROR: cgca_gl: one or more super array&
                  & coordinate(s) are ouside the bounds"
 write (*,'(a,3(i0,tr1))') "ERROR: cgca_gl: super array coord: ",super
 write (*,'(a)') "ERROR: cgca_gl: lower bound must be 1"
 write (*,'(a,3(i0,tr1))') "ERROR: cgca_gl: upper bounds: ", usup
 error stop
end if

! actual calculation

imgpos = lcob + (super-1)/szr
local = lbr(1:3) + super-szr*(imgpos-lcob) - 1

! checks after

if (any(imgpos .gt. ucob) .or. any(imgpos .lt. lcob)) then
 write (*,'(a,i0)') "ERROR: cgca_gl: image", thisimage
 write (*,'(a)') "ERROR: cgca_lg: one or more image positions&
                  & are ouside the bounds"
 write (*,'(a,3(i0,tr1))') &
  "ERROR: cgca_gl: image positions: ",imgpos
 write (*,'(a,3(i0,tr1))') &
  "ERROR: cgca_gl: lower image grid bounds: ", lcob
 write (*,'(a,3(i0,tr1))') &
  "ERROR: cgca_gl: upper image grid bounds: ", ucob
 error stop
end if

if (any(local .gt. ubr(1:3)) .or. any(local .lt. lbr(1:3))) then
 write (*,'(a,i0)') "ERROR: cgca_gl: image", thisimage
 write (*,'(a)') "ERROR: cgca_lg: one or more local coordinates &
                  & are ouside the bounds"
 write (*,'(a,3(i0,tr1))') "ERROR: cgca_gl: local coordinates: ",local
 write (*,'(a,3(i0,tr1))') "ERROR: cgca_gl: lower bounds: ", lbr
 write (*,'(a,3(i0,tr1))') "ERROR: cgca_gl: upper bounds: ", ubr
 error stop
end if

end subroutine cgca_gl

!*roboend*


!*robodoc*s* cgca_m2glm/cgca_lg
!  NAME
!    cgca_lg
!  SYNOPSIS

subroutine cgca_lg(imgpos,local,coarray,super)

!  INPUTS

integer(kind=idef),intent(in) :: imgpos(3),local(3)
integer(kind=iarr),allocatable,intent(inout) :: coarray(:,:,:,:)[:,:,:]

! OUTPUT

integer(kind=idef),intent(out) :: super(3)

! DESCRIPTION
!  This routine converts the image coordinates in the grid and the local
!  cell coordinates in this image into the global cell coordinates in
!  the super array:
!   - imgpos(3) is the image position in the grid
!   - local(3) are cell coordinates in that image's array
!   - coarray is the model
!   - super(3) are cell coordinates in a super array
! NOTES
!   The global, super, coordinates must start from 1!
!
!   Any image can call this routine
! USES
! USED BY
!   cgca_gbf1f
! SOURCE

integer :: & 
 lbr(4)   ,& ! lower bounds of the "real" coarray, lbv+1
 ubr(4)   ,& ! upper bounds of the "real" coarray, ubv-1
 szr(3)   ,& ! size or the "real" coarray, ubr-lbr+1
 lcob(3)  ,& ! lower cobounds of the coarray
 ucob(3)  ,& ! upper cobounds of the coarray
 usup(3)  ,& ! upper bound of the super array, szr*(ucob-lcob+1)
 thisimage

thisimage = this_image()

! check for allocated

if (.not. allocated(coarray)) then
 write (*,'(a,i0)') "ERROR: cgca_lg: image", thisimage
 write (*,'(a)') "ERROR: cgca_lg: coarray is not allocated"
 error stop
end if

lbr=lbound(coarray)+1
ubr=ubound(coarray)-1

! the 4th dimension is to do with the number of cell state
! types. This is not relevant here.
szr=ubr(1:3)-lbr(1:3)+1

lcob=lcobound(coarray)
ucob=ucobound(coarray)
usup=szr*(ucob-lcob+1)

! check for bounds

if (any(imgpos .gt. ucob) .or. any(imgpos .lt. lcob)) then
 write (*,'(a,i0)') "ERROR: cgca_lg: image", thisimage
 write (*,'(a)') "ERROR: cgca_lg: one or more image positions&
                  & are ouside the bounds"
 write (*,'(a,3(i0,tr1))') "ERROR: cgca_lg: image positions: ",imgpos
 write (*,'(a,3(i0,tr1))') &
  "ERROR: cgca_lg: lower image grid bounds: ", lcob
 write (*,'(a,3(i0,tr1))') &
  "ERROR: cgca_lg: upper image grid bounds: ", ucob
 error stop
end if

if (any(local .gt. ubr(1:3)) .or. any(local .lt. lbr(1:3))) then
 write (*,'(a,i0)') "ERROR: cgca_lg: image", thisimage
 write (*,'(a)') "ERROR: cgca_lg: one or more local coordinates&
                  & are ouside the bounds"
 write (*,'(a,3(i0,tr1))') &
  "ERROR: cgca_lg: local coordinates: ", local
 write (*,'(a,3(i0,tr1))') &
  "ERROR: cgca_lg: lower bounds: ", lbr
 write (*,'(a,3(i0,tr1))') &
  "ERROR: cgca_lg: upper bounds: ", ubr
 error stop
end if

! actual calculation

super = szr*(imgpos-lcob) + local-lbr(1:3)+1

! check for bounds

if (any(super .gt. usup) .or. any(super .lt. 1)) then
 write (*,'(a,i0)') "ERROR: cgca_lg: image", thisimage
 write (*,'(a)') "ERROR: cgca_lg: one or more super array &
                  & coordinates are ouside the bounds"
 write (*,'(a,3(i0,tr1))') &
  "ERROR: cgca_lg: super array coord: ",super
 write (*,'(a)') "ERROR: cgca_lg: lower bound must be 1"
 write (*,'(a,3(i0,tr1))') "ERROR: cgca_lg: upper bounds: ", usup
 error stop
end if

end subroutine cgca_lg

!*roboend*


!*robodoc*s* cgca_m2glm/cgca_ico
!  NAME
!    cgca_ico
!  SYNOPSIS

subroutine cgca_ico( ind, cosub, flag )

!  INPUTS

integer( kind=idef ), intent(in) :: ind

!  OUTPUT

integer( kind=idef ), intent(out) :: cosub( cgca_scodim ), flag

!  DESCRIPTION
!    This routine converts image index of coarray into its set
!    of cosubscripts.
!  NOTES
!    This is a serial routine, just computation, no inter-image
!    communication is involved.
!    Any and all images can call this routine.
!    flag is set to 0 on normal exit. flag is 1 if the coarray index
!    .lt. 0.
!  SOURCE

integer( kind=idef ) :: codim( cgca_scodim ), step, step2, rem, rem2

! Set as default
flag = 0

! Sanity check
if ( ind .le. 0 ) then

  ! Set the flag and return immediately.
  flag = 1
  return
end if

! codimensions
codim = cgca_sucob - cgca_slcob + 1

! along 1
step = mod( ind, codim(1) )
if ( step .eq. 0 ) step = codim(1)
cosub(1) = cgca_slcob(1) + step - 1

! along 2
! number of full layers
step = ind / ( codim(1) * codim(2) )

! number of images in the last unfilled layer
rem = mod( ind , codim(1) * codim(2) )

! if all layers ar filled, take step2 as codim(2)
if ( rem .eq. 0 ) then
  step2 = codim(2)
else
  ! number of full columns in the last unfilled layer
  step2 = rem / codim(1) 
end if

! number of images in the last unfilled column
rem2 = mod( rem, codim(1) )

! if it's not zero, increment the step
if ( rem2 .ne. 0 ) step2 = step2 + 1

cosub(2) = cgca_slcob(2) + step2 -1 

! along 3
if ( rem .ne. 0 ) step = step + 1
cosub(3) = cgca_slcob(3) + step - 1

end subroutine cgca_ico

!*roboend*


!*robodoc*s* cgca_m2glm/cgca_ico2
!  NAME
!    cgca_ico2
!  SYNOPSIS


subroutine cgca_ico2( lcob, ucob, ind, cosub )

!  INPUTS

integer( kind=idef ), intent(in) :: lcob(:), ucob(:), ind

!  OUTPUT

integer( kind=idef ), intent(out) :: cosub( size(lcob) )

!  DESCRIPTION
!    This routine converts image index of coarray into its set
!    of cosubscripts. It borrows the code from the f2008 standard,
!     http://j3-fortran.org/doc/year/10/10-007r1.pdf
!    Section C.10.1.
!  NOTES
!    This is a serial routine, just computation, no inter-image
!  SOURCE

integer :: n, i, m, ml, extent

n = size( cosub )
m = ind - 1
do i = 1, n-1
  extent = ucob(i) - lcob(i) + 1
  ml = m
  m = m / extent
  cosub( i ) = ml - m * extent + lcob(i)
end do
cosub( n ) = m + lcob( n )

end subroutine cgca_ico2

!*roboend*

end module cgca_m2glm
