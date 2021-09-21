!*robodoc*m* CGPACK/cgca_m2rot
!  NAME
!    cgca_m2rot
!  SYNOPSIS

!$Id: cgca_m2rot.f90 526 2018-03-25 23:44:51Z mexas $

module cgca_m2rot

!  DESCRIPTION
!    Module dealing with tensor rotations, orientations,
!    mis-orientations
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    cgca_rt, cgca_ckrt, rtprint, cgca_csym, cgca_csym_pure,
!    cgca_mis, cgca_miscsym
!  USES
!    cgca_m1co
!  USED BY
!    cgca
!  SOURCE

use cgca_m1co, only : idef, ldef, rdef, cgca_pi
implicit none

private
public :: cgca_ckrt, cgca_csym, cgca_csym_pure, cgca_mis, &
  cgca_miscsym, cgca_rt

contains

!*roboend*


!*robodoc*s* cgca_m2rot/cgca_rt
!  NAME
!    cgca_rt
!  SYNOPSIS

subroutine cgca_rt(r)

!  INPUT
!    Array r to store output
!  OUTPUT

real(kind=rdef),allocatable,intent(inout) :: r(:,:,:)[:,:,:]

!  DESCRIPTION
!    Choose grain rotation tensors at random and return in a
!    coarray r(:,:,:)[:,:,:]. dimension 1 of the coarray is the grain
!    number. dimensions 2 and 3 are for the rotation tensor for the
!    grain, e.g. r(387,2,1) is the rotation tensor component r21 for
!    grain 387.
!
!    Image 1 assigns rotation tensors to all grains. then all other
!    images copy the coarray from image 1.
!
!    The method: choose 3 random angles. Then interpret them as
!    rotations about axes 1, 2 and 3 in turn. The resulting rotation
!    will be quite random.
!
!    Grain axes are the crystallographic directions:
!
!      axis  direction
!        1     [100]
!        2     [010]
!        3     [001]
! 
!    Then the grain orientation tensor for some grain z is defined as:
!
!      r(z,1,1) r(z,1,2) r(z,1,3)
!      r(z,2,1) r(z,2,2) r(z,2,3)
!      r(z,3,1) r(z,3,2) r(z,3,3)
!
!    Index 2  is the grain axis, index 3 is the space axis, i.e.
!    r(z,i,j)= cos(x_i_grain,x_j_space), e.g.
!    r(z,3,2) = cos(3_grain,2_space) is the angle between grain
!    axis 3 and the space axix 2 for grain z.
!
!    This leads to the important convention:
!
!      x_grain =           r  * x_space
!      x_space = transpose(r) * x_grain
!
!    *All* other routines must adhere to this convention.
!  NOTES
!    Call cgca_art before calling this routine.
!  USES
!    none
!  USED BY
!    none, end user
!  SOURCE  

real( kind=rdef ), parameter :: zero = 0.0_rdef, one = 1.0_rdef,       &
  eight=8.0_rdef, twopi = 2 * cgca_pi
real(kind=rdef) :: a(3),q(3,3)
real :: rndnum(3)

integer(kind=idef) :: l(3),u(3),lcob(3),ucob(3),i

! check for allocated
if (.not. allocated(r)) error stop "error in cgca_rt: r is not allocated"

! define cobounds on all images. these are used at the end
! to read the values from image 1.
lcob=lcobound(r)
ucob=ucobound(r)

! Image 1 does all work, otherwise lots of syncs will be needed
image1: if (this_image() .eq. 1) then

! coarray bounds
l=lbound(r)
u=ubound(r)

! check that bounds for dimensions 2 and 3 are 1:3
if (l(2) .ne. 1 .or. l(3) .ne. 1 .or. &
    u(2) .ne. 3 .or. u(3) .ne. 3) then
  write (*,'(a)') &
    "error in cgca_rt: coarray bounds along dimensions 2 or 3 are wrong."
  write (*,'(a)') &
    "error in cgca_rt: these must be [1:3]."
  error stop
end if

! loop over all grains
main: do i=l(1),u(1)

call random_number(rndnum)         ! in [0,1)
a=real(rndnum,kind=rdef)*twopi     ! a in [0,2*pi)

! rotation about 1
 r(i,:,:)=zero
 r(i,1,1)=one
 r(i,2,2)=cos(a(1))
 r(i,2,3)=sin(a(1))
 r(i,3,2)=-r(i,2,3)
 r(i,3,3)=r(i,2,2)

! rotation about 2
 q=zero
 q(1,1)=cos(a(2))
 q(1,3)=sin(a(2))
 q(2,2)=one
 q(3,1)=-q(1,3)
 q(3,3)=q(1,1)

! intermediate compound rotation
 r(i,:,:) = matmul(q,r(i,:,:))

! rotation about 3
 q=zero
 q(1,1)=cos(a(3))
 q(1,2)=sin(a(3))
 q(2,1)=-q(1,2)
 q(2,2)=q(1,1)
 q(3,3)=one

! final compound rotation
 r(i,:,:) = matmul(q,r(i,:,:))

end do main

end if image1 

! global sync here
sync all

! all images read r from image 1
r(:,:,:) = r(:,:,:)[lcob(1),lcob(2),lcob(3)] 

! exit only after all images have the rotation tensors assigned.
sync all

end subroutine cgca_rt

!*roboend*


!*robodoc*s* cgca_m2rot/cgca_ckrt
!  NAME
!    cgca_ckrt
!  SYNOPSIS

subroutine cgca_ckrt(r,debug,flag)

!  INPUTS

real(kind=rdef),intent(in) :: r(3,3)
logical(kind=ldef),intent(in) :: debug

!  OUTPUT

integer(kind=idef),intent(out) :: flag

!  DESCRIPTION
!    Check that the given rotation tensor r(3,3) is
!    orthogonal.
!    If debug .eq. .true., then  verbose diagnostics is
!    printed on error.
!    In this case a non-zero return flag shows the number
!    of the failed test.
!  USES
!    rtprint
!  USED BY
!    none, end user
!  SOURCE

real(kind=rdef),parameter :: one=1.0_rdef,factor=1.0e1_rdef
integer :: i
real :: maxerr, orthogon(3,3), tmp(3,3)

maxerr = factor * epsilon(one)

flag = 0 
tmp = transpose(r)
 do i=1,2
   if ( i .eq. 1 ) orthogon = matmul( r, tmp )
   if ( i .eq. 2 ) orthogon = matmul( tmp, r )

   if (abs(orthogon(1,1)-one) .gt. maxerr) then
     if (debug) call rtprint( r, tmp, orthogon, maxerr )
     flag = 1 * i
     return
   end if

   if (abs(orthogon(1,2))     .gt. maxerr) then
     if (debug) call rtprint( r, tmp, orthogon, maxerr )
     flag = 2 * i
     return
   end if

   if (abs(orthogon(1,3))     .gt. maxerr) then
     if (debug) call rtprint( r, tmp, orthogon, maxerr )
     flag = 3 * i
     return
   end if

   if (abs(orthogon(2,1))     .gt. maxerr) then
     if (debug) call rtprint( r, tmp, orthogon, maxerr )
     flag = 4 * i
     return
   end if

   if (abs(orthogon(2,2)-one) .gt. maxerr) then
     if (debug) call rtprint( r, tmp, orthogon, maxerr )
     flag = 5 * i
     return
   end if

   if (abs(orthogon(2,3))     .gt. maxerr) then
     if (debug) call rtprint( r, tmp, orthogon, maxerr )
     flag = 6 * i
     return
   end if

   if (abs(orthogon(3,1))     .gt. maxerr) then
     if (debug) call rtprint( r, tmp, orthogon, maxerr )
     flag = 7 * i
     return
   end if

   if (abs(orthogon(3,2))     .gt. maxerr) then
     if (debug) call rtprint( r, tmp, orthogon, maxerr )
     flag = 8 * i
     return
   end if

   if (abs(orthogon(3,3)-one) .gt. maxerr) then
     if (debug) call rtprint( r, tmp, orthogon, maxerr )
     flag = 9 * i
     return
   end if

 end do

end subroutine cgca_ckrt

!*roboend*


!*robodoc*s* cgca_m2rot/rtprint
!  NAME
!    rtprint
!  SYNOPSIS

subroutine rtprint(a,b,c,err)

!  INPUT

real(kind=rdef),intent(in) :: a(3,3),b(3,3),c(3,3),err


!  SIDE EFFECTS
!    dumps some text on stdout
!  DESCRIPTION
!    This routine prints on stdout the details of the
!    rotation tensor that failed one of the tests.
!    It prints the tensor itself, as a matrix, the
!    transposed tensor, and then a product.
!  NOTES
!    This routine is not accessible from outside of
!    module cgca_m2rot.
!  USES
!    none
!  USED BY
!    cgca_ckrt
!  SOURCE

integer :: i

write (*,'(a)') "troublesome tensor:"
do i=1,3
  write (*,*) a(i,:)
end do

write (*,*)

write (*,'(a)') "transposed troublesome tensor:"
do i=1,3
  write (*,*) b(i,:)
end do

write (*,*)
write (*,'(a)') "this should have been the unit tensor:"
do i=1,3
  write (*,*) c(i,:)
end do

write (*,*)"max allowed error was:", err

end subroutine rtprint

!*roboend*


!*robodoc*s* cgca_m2rot/cgca_csym
!  NAME
!    cgca_csym
!  SYNOPSIS

subroutine cgca_csym( num, rs )

!  INPUT
 
integer, intent( in ) :: num

!  OUTPUT

real( kind=rdef ), intent( out ) :: rs(3,3)

!  NOTES
!    The symmetry tensors and the misorientation angle
!    equations can be found in "Introduction to texture
!    analysis : macrotexture, microtexture, and
!    orientation mapping", Olaf Engler, Valerie Randle,
!    2nd ed, Boca Raton, CRC Press, 2010, 456 p.
!
!    There is a copy in QBL.
!  DESCRIPTION
!    This routine stores, and outputs on demand, symmetry rotation
!    tensors, for cubic crystals, 24 in total.
!    The first tensor is the unit tensor (trivial case).
!    Then there are 23 non-trivial tensors.
!    There are 24*3*3=216 elements in total:
!    - num - rotation symmetry tensor number
!    - rs - rotation symmetry tensor
!  USES
!    none
!  USED BY
!    cgca_miscsym
!  SOURCE

real( kind=rdef ), parameter :: mnsone = -1.0_rdef,      &
                                  zero =  0.0_rdef,      &
                                   one =  1.0_rdef

! Keep the identity tensor as well, at number zero, if required in future.
! Data are filled columns first!

real( kind=rdef ), parameter :: r(3,3,24) = reshape(     &
(/ one,zero,zero,    zero,one,zero,    zero,zero,one,    &
   one,zero,zero,    zero,zero,mnsone, zero,one,zero,    &
   one,zero,zero,    zero,mnsone,zero, zero,zero,mnsone, &

   one,zero,zero,    zero,zero,one,    zero,mnsone,zero, &
   zero,zero,one,    zero,one,zero,    mnsone,zero,zero, &
   mnsone,zero,zero, zero,one,zero,    zero,zero,mnsone, &

   zero,zero,mnsone, zero,one,zero,    one,zero,zero,    &
   zero,mnsone,zero, one,zero,zero,    zero,zero,one,    &
   mnsone,zero,zero, zero,mnsone,zero, zero,zero,one,    &
   zero,one,zero,    mnsone,zero,zero, zero,zero,one,    &

   zero,zero,one,    one,zero,zero,    zero,one,zero,    &
   mnsone,zero,zero, zero,zero,one,    zero,one,zero,    &
   zero,zero,mnsone, mnsone,zero,zero, zero,one,zero,    &

   zero,mnsone,zero, zero,zero,mnsone, one,zero,zero,    &
   mnsone,zero,zero, zero,zero,mnsone, zero,mnsone,zero, &
   zero,one,zero,    zero,zero,mnsone, mnsone,zero,zero, &

   zero,zero,one,    zero,mnsone,zero, one,zero,zero,    &
   zero,zero,mnsone, zero,mnsone,zero, mnsone,zero,zero, &
   zero,one,zero,    one,zero,zero,    zero,zero,mnsone, &
   zero,mnsone,zero, mnsone,zero,zero, zero,zero,mnsone, &

   zero,zero,one,    mnsone,zero,zero, zero,mnsone,zero, &
   zero,zero,mnsone, one,zero,zero,    zero,mnsone,zero, &
   zero,mnsone,zero, zero,zero,one,    mnsone,zero,zero, &
   zero,one,zero,    zero,zero,one,    one,zero,zero /), &
   (/3,3,24/) )

! sanity check

if (num .lt. 1 .or. num .gt. 24) then 
  write (*,'(a,i0)') "ERROR: cgca_csym: image: ", this_image()
  write (*,'(a)') "ERROR: cgca_csym: num is out of range [1..24]"
  error stop
end if

! simply return the required rotation tensor

rs = r(:,:,num)

end subroutine cgca_csym

!*roboend*


!*robodoc*s* cgca_m2rot/cgca_csym_pure
!  NAME
!    cgca_csym_pure
!  SYNOPSIS

pure subroutine cgca_csym_pure( num, rs, flag )

!  INPUT
!    - num - rotation symmetry tensor number
 
integer( kind=idef ), intent(in) :: num

!  OUTPUT
!    - rs - rotation symmetry tensor

real( kind=rdef ), intent(out) :: rs(3,3)
integer, intent(out) :: flag

!  NOTES
!    This is a PURE subroutine. Hence no external IO is required.
!    Hence I'll report back the errors via a flag:
!      0 - sucessful completion.
!      1 - number is out of range [1..24]. 
!  DESCRIPTION
!    This routine stores, and outputs on demand, symmetry rotation
!    tensors, for cubic crystals, 24 in total.
!    The first tensor is the unit tensor (trivial case).
!    Then there are 23 non-trivial tensors.
!    There are 24*3*3=216 elements in total.

!    The symmetry tensors and the misorientation angle
!    equations can be found in "Introduction to texture
!    analysis : macrotexture, microtexture, and
!    orientation mapping", Olaf Engler, Valerie Randle,
!    2nd ed, Boca Raton, CRC Press, 2010, 456 p.
!  USES
!    none
!  USED BY
!    cgca_miscsym
!  SOURCE

real( kind=rdef ), parameter :: mnsone = -1.0_rdef,      &
                                  zero =  0.0_rdef,      &
                                   one =  1.0_rdef

! Keep the identity tensor as well, at number zero,
! if required in future. Data are filled columns first!

real( kind=rdef ), parameter :: r(3,3,24) = reshape(     &
(/ one,zero,zero,    zero,one,zero,    zero,zero,one,    &
   one,zero,zero,    zero,zero,mnsone, zero,one,zero,    &
   one,zero,zero,    zero,mnsone,zero, zero,zero,mnsone, &

   one,zero,zero,    zero,zero,one,    zero,mnsone,zero, &
   zero,zero,one,    zero,one,zero,    mnsone,zero,zero, &
   mnsone,zero,zero, zero,one,zero,    zero,zero,mnsone, &

   zero,zero,mnsone, zero,one,zero,    one,zero,zero,    &
   zero,mnsone,zero, one,zero,zero,    zero,zero,one,    &
   mnsone,zero,zero, zero,mnsone,zero, zero,zero,one,    &
   zero,one,zero,    mnsone,zero,zero, zero,zero,one,    &

   zero,zero,one,    one,zero,zero,    zero,one,zero,    &
   mnsone,zero,zero, zero,zero,one,    zero,one,zero,    &
   zero,zero,mnsone, mnsone,zero,zero, zero,one,zero,    &

   zero,mnsone,zero, zero,zero,mnsone, one,zero,zero,    &
   mnsone,zero,zero, zero,zero,mnsone, zero,mnsone,zero, &
   zero,one,zero,    zero,zero,mnsone, mnsone,zero,zero, &

   zero,zero,one,    zero,mnsone,zero, one,zero,zero,    &
   zero,zero,mnsone, zero,mnsone,zero, mnsone,zero,zero, &
   zero,one,zero,    one,zero,zero,    zero,zero,mnsone, &
   zero,mnsone,zero, mnsone,zero,zero, zero,zero,mnsone, &

   zero,zero,one,    mnsone,zero,zero, zero,mnsone,zero, &
   zero,zero,mnsone, one,zero,zero,    zero,mnsone,zero, &
   zero,mnsone,zero, zero,zero,one,    mnsone,zero,zero, &
   zero,one,zero,    zero,zero,one,    one,zero,zero /), &
   (/3,3,24/) )

! Set to zero initially
rs = zero

! sanity check
if (num .lt. 1 .or. num .gt. 24) then 
  flag = 1
  return
end if

! If there are no error conditions, simply return the
! required rotation tensor and set the flag to 0.

flag = 0
  rs = r( :, :, num )

end subroutine cgca_csym_pure

!*roboend*


!*robodoc*s* cgca_m2rot/cgca_mis
!  NAME
!    cgca_mis
!  SYNOPSIS

subroutine cgca_mis(r1,r2,angle)

!  INPUTS

!real(kind=rdef),intent(in) :: r1(3,3),r2(3,3)
real(kind=rdef),intent(in) :: r1(:,:),r2(:,:)
!  OUTPUT

real(kind=rdef),intent(out) :: angle

!  DESCRIPTION
!    This routine calculates grain misoreintation.
!    Given 2 orientation tensors, r1 and r2, the
!    misorientation angle (in rad) is:
!    acos((tr(r1*r2^t)-1)/2), where "tr" is tensor trace.
!    The angle is [0,pi].
!  USES
!    none
!  USED BY
!    cgca_miscsym
!  SOURCE

real(kind=rdef),parameter :: one=1.0_rdef
real(kind=rdef) :: misor(3,3),trace,arg

misor=matmul(r1,transpose(r2))
trace=misor(1,1)+misor(2,2)+misor(3,3)
arg = (trace-one) / 2

if (arg .gt.  one) arg= one
if (arg .lt. -one) arg=-one

angle = acos(arg)

if (isnan(angle)) then
  write (*,'(a,i0)') "ERROR: cgca_mis: image: ", this_image()
  write (*,'(a,i0)') "ERROR: cgca_mis: arg: ", arg
  write (*,'(a,i0)') "ERROR: cgca_mis: angle=acos(arg) is NAN"
  error stop
end if

end subroutine cgca_mis

!*roboend*


!*robodoc*s* cgca_m2rot/cgca_miscsym
!  NAME
!    cgca_miscsym
!  SYNOPSIS

subroutine cgca_miscsym(r1,r2,minang)

!  INPUTS

!real(kind=rdef),intent(in) :: r1(3,3),r2(3,3)
real(kind=rdef),intent(in) :: r1(:,:),r2(:,:)
!  OUTPUT

real(kind=rdef),intent(out) :: minang

!  USES
!    cgca_csym, cgca_mis
!  USED BY
!  DESCRIPTION
!    This routine calculates the grain misorientation. angle,
!    taking cubic symmetry into account.
!    Given 2 orientation tensors, r1 and r2, the
!    misorientation angle (in rad) is: acos((tr(r1*r2^t)-1)/2), where "tr"
!    is tensor trace. the angle is [0,pi].
!  SOURCE

real(kind=rdef) :: rot(3,3), angle, tmp(3,3)
integer :: i

minang = 1.0e1 ! any number .gt. PI will do
angle  = 0.0

do i=1,24
  call cgca_csym(i,rot)
  tmp = matmul(rot,r2)
  !call cgca_mis(r1,matmul(rot,r2),angle)
  call cgca_mis(r1,tmp,angle)
  if (angle .lt. minang) minang=angle
end do

end subroutine cgca_miscsym

!*roboend*

end module cgca_m2rot
