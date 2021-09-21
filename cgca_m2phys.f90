!*robodoc*m* CGPACK/cgca_m2phys
!  NAME
!    cgca_m2phys
!  SYNOPSIS

!$Id: cgca_m2phys.f90 530 2018-03-26 16:10:00Z mexas $

module cgca_m2phys

!  DESCRIPTION
!    Module dealing with physical units - length and time.
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    cgca_cadim, cgca_gdim, cgca_imco
!  USES
!    cgca_m1co
!  USED BY
!    perhaps the user directly?
!  SOURCE

use cgca_m1co, only: idef, iarr, rdef
implicit none

private
public :: cgca_gdim, cgca_cadim, cgca_imco

contains

!*roboend*


!*robodoc*s* cgca_m2phys/cgca_gdim
!  NAME
!    cgca_gdim
!  SYNOPSIS

subroutine cgca_gdim( n, ir, qual )

!  INPUT

integer( kind=idef ), intent( in ) :: n

!  OUTPUT

integer( kind=idef ), intent( out ) :: ir(3)
real( kind=rdef ), intent( out ) :: qual

!  SIDE EFFECTS
!    None
!  DESCRIPTION
!    The purpose of this routine is to find 3 coarray
!    grid dimensions, ir(1) >= ir(2) >= ir(3),
!    such that for a given number of images
!    the coarray grid is as "cubic" as possible.
!    In mathematical terms the aim is to find
!    F = min( max( ir(1) - ir(3) ) ).
!    The quality of this minimum is defined as QUAL=1-F/(N-1).
!    Inputs:
!    - N is the total number of images, num_images().
!    Outputs:
!    - ir - the array of 3 coarray grid dimensions.
!    - qual - is the quality of the fitted grid.
!      qual=1 means F=0, i.e. the coarray grid is a cube.
!      qual=0 means F=N-1, i.e. the coarray grid is 1D, i.e. [N,1,1].
!    The outputs of this routine, ir, are used to choose
!    the dimensions of space coarray.
!    The QUAL output is for information only.
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  NOTES
!    If N<1 is given, then the routine returns immediately
!    with QUAL=-1. So the caller can/should check for this condition.
!  USES
!  USED BY
!    Supposed to be called prior to calling cgca_cadim
!  SOURCE

real( kind=rdef ), parameter :: one = 1.0_rdef, zero = 0.0_rdef, &
  third = one/3.0_rdef
integer( kind=idef) :: i, j, k, f, ftrial

! default output
  ir(1) = n
ir(2:3) = 1
   qual = zero

! return immediately if n=1
if ( n .eq. 1 ) return

! return with QUAL=-1 if N<1
if ( n .lt. 1 ) then
  qual = -one
  return
end if
 
! Set the initial value of the objective function
f = n-1

loopi1: do i = n/2, int( real(n)**third ), -1
          if ( mod( n,i ) .ne. 0 ) cycle loopi1
loopi2:   do j = n/i, 2, -1
            if ( j .gt. i ) cycle loopi2 
            if ( mod( n,(i*j) ) .ne. 0 ) cycle loopi2
            k = n/(i*j)
            if ( k .gt. j ) cycle loopi2
            ftrial = i-k
            if ( ftrial .ge. f ) cycle loopi2
            f = ftrial
            ir(1) = i
            ir(2) = j
            ir(3) = k
            if ( f .eq. 0 ) exit loopi1
          end do loopi2
        end do loopi1

qual = one - real(f)/(n-one)

end subroutine cgca_gdim

!*roboend*


!*robodoc*s* cgca_m2phys/cgca_cadim
!  NAME
!    cgca_cadim
!  SYNOPSIS

subroutine cgca_cadim( bsz, res, dm, ir, c, lres, ng )

!  INPUT

real( kind=rdef ), intent( inout ) :: bsz(3)
real( kind=rdef ), intent( in ) :: res, dm
integer, intent( inout ) :: ir(3)

!  OUTPUT

integer, intent( out ) :: c(3), ng
real( kind=rdef ), intent( out ) :: lres

!  SIDE EFFECTS
!    Arrays bsz and ir are input/output. On exit the values of ir are
!    rearranged. The values of bsz are changed.
!  DESCRIPTION
!    Inputs:
!    - bsz - box size, the size of the CA space in physical
!      units of length. The unit itself is not defined.
!      It's use must be consistent across the whole of
!      CGPACK. In particular, speeds will depend on the choice
!      of the length unit.
!    - res - the model resolution, cells per grain. Note that
!      this is *not* spatial resolution. The meaning is that
!      RES cells are required to resolve the shape of a grain.
!      This setting should not depend on the grain size.
!    - dm - the mean grain size, in physical units.
!    - ir - coarray grid dimensions. The intention is that ir
!      is calculate by a call to cgca_gdim. 
!    Outputs:
!    - ir - rearranged coarray grid dimensions 
!    - bsz - new box dimension is calculated, see note 3.
!    - c - numbers of cells in the space coarray
!    - lres - linear resolution, cells per unit of length 
!    - ng - number of grains in the whole model
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  NOTES
!    1. The space coarray should be declared with something like this
!         space( c(1), c(2), c(3), nlayers )[ir(1),ir(2),*]
!    on all images.
!
!    2. An important feature is that the coarray grid
!    dimensions can be rearranged to better suit the physical
!    dimensions of the "box". For example, on input ir(1)
!    is the biggest. However, if the physical dimension of
!    the box is smallest along 1, it makes sense to swap
!    ir(1) with ir(3). This will help achive balanced
!    linear resolution along each dimension.
!
!    3. It not generally possible to have the same
!    linear resolution is all 3 directions and to keep
!    the physical size of the box exactly as given.
!    The decision is made here to give preference to
!    the linear resolution. So the algorithm chooses
!    the same linear resolution in all 3 directions.
!    As a consequence, the physical sizes of the box
!    along all 3 directions can be slightly smaller or
!    bigger than given. The biggest deviations probably
!    arise when the shape is very far from cubic.
!  USES
!  USED BY
!    via module cgca_m2phys
!  SOURCE

real( kind=rdef ), parameter :: one = 1.0_rdef, third = one/3.0_rdef
real( kind=rdef ) :: &
 bvol,               & ! "box" volume in physical units
 dm3                   ! dm**3

! Derived type for easy sorting
type bsz_order
  real( kind=rdef ) :: bsz
  integer :: i
end type bsz_order

type( bsz_order ) :: box(3), tmp 

! Assume that ir is already sorted in decreasing order,
! as returned by cgca_gdim. Check this assertion.
if ( ir(1) .lt. ir(2) .or. ir(2) .lt. ir(3) .or. ir(1) .lt. ir(3) ) then
  write (*,*) "ERROR: cgca_m2phys/cgca_cadim: ir is not sorted on entry"
  error stop
end if

! Set the initial descending order, 1 to 3
box(1) = bsz_order( bsz(1), 1 )
box(2) = bsz_order( bsz(2), 2 )
box(3) = bsz_order( bsz(3), 3 )

!write (*,*) "ir on entry:", ir
!write (*,*) "box:", box

! Now sort box%bsz in decreasing order and the resulting order
! or box%i is what I need.

if ( box(1)%bsz .lt. box(2)%bsz ) then
     tmp = box(1)
  box(1) = box(2)
  box(2) = tmp
end if

if ( box(2)%bsz .lt. box(3)%bsz ) then
     tmp = box(2)
  box(2) = box(3)
  box(3) = tmp
end if

if ( box(1)%bsz .lt. box(2)%bsz ) then
     tmp = box(1)
  box(1) = box(2)
  box(2) = tmp
end if

! Reassign elements in desired order
ir( box%i ) = ir

!write (*,*) "ir on exit:", ir

! box volume
bvol = bsz(1)*bsz(2)*bsz(3)

! grain volume
dm3 = dm**3

! number of grains in the whole model, integer
ng = int( bvol/dm3 )

! linear resolution
! res**third      - cells per mean grain size length, dm
! res**third / dm - cells per unit length
lres = res**third / dm

! numbers of cells
! res**third/dm * bsz(i)       - cells per box length along i
! res**third/dm * bsz(i)/ir(i) - cells per image along i
c = nint( lres * bsz/ir )

! cannot have zero as the array dimension
if ( c(1) .eq. 0 ) c(1) = 1
if ( c(2) .eq. 0 ) c(2) = 1
if ( c(3) .eq. 0 ) c(3) = 1

! warn the user, the box is likely to be very different 
! from the input values.
if ( any( c .eq. 1 ) ) then
  write (*,"(a)")                                                      &
   "WARN: cgca_cadim: the new box sizes are probably wrong, check"
end if

! new box size
bsz = real( ir*c, kind=rdef ) / lres

end subroutine cgca_cadim

!*roboend*


!*robodoc*s* cgca_m2phys/cgca_imco
!  NAME
!    cgca_imco
!  SYNOPSIS

subroutine cgca_imco( space, lres, bcol, bcou )

!  INPUT

integer( kind=iarr ), allocatable, intent( in ) :: &
 space(:,:,:,:)[:,:,:]
real( kind=rdef ), intent( in ) :: lres

!  OUTPUT

real( kind=rdef ), intent( out ) :: bcol(3), bcou(3)

!  DESCRIPTION
!    IMCO stands for IMage COordinates. This routine calculates
!    the lower and the upper physical coordinates of the
!    coarray on this image in CA CS.
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  INPUTS
!    - space - the model coarray. Error will result if it's
!      not allocated. The array is used only to calculate the
!      position of this image within the coarray grid, as
!       this_image( space )
!    - lres - linear resolution of the space cooarray,
!      cells per unit of length.
!  OUTPUTS
!    - bcol - physical coordinates of the lower bounds of the coarray
!      on this image in CA CS.
!    - bcou - physical coordinates of the upper bounds of the coarray
!      on this image in CA CS.
!  SIDE EFFECTS
!      none

real( kind=rdef ) :: one = 1.0_rdef, half = 0.5_rdef

integer( kind=idef ) :: &
 img(3),                & ! the coarray image grid is always 3D
 szsp(3)                  ! size of the space array, on one image!

real( kind=rdef) :: hcsz  ! half of the cell physical size

! this image index
img = this_image( space )

! Size of the space array.
! Remember that space array has 1 halo cell on each
! boundary along each dimension. Don't count those:
szsp(1) = size( space, dim=1 ) - 2_idef
szsp(2) = size( space, dim=2 ) - 2_idef 
szsp(3) = size( space, dim=3 ) - 2_idef

! vectors in CA CS
bcol = real( ((img-1) * szsp(1:3) + 1), kind=rdef) / lres
bcou = real( img * szsp(1:3), kind=rdef ) / lres

! Make sure there are no gaps between the upper
! and the next lower boundary. The gap equals to the
! physical size of a single cell. So add half a cell size
! to the upper boundary and subtract half a cell size from
! the lower boundary.
hcsz = half * ( one / lres )
bcol = bcol - hcsz
bcou = bcou + hcsz

end subroutine cgca_imco

!*roboend*

!*********************************************************************72

end module cgca_m2phys
