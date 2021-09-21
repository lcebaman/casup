!*robodoc*m* CGPACK/cgca_m2alloc
!  NAME
!    cgca_m2alloc
!  SYNOPSIS

!$Id: cgca_m2alloc.f90 381 2017-03-22 11:29:44Z mexas $

module cgca_m2alloc

!  DESCRIPTION
!    Module dealing with the allocation and deallocation of various
!    arrays. Several routines are used because they allocate arrays
!    of different dimensionality, i.e.
!          (:)[:,:,:]
!      (:,:,:)[:,:,:]
!    (:,:,:,:)[:,:,:]
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    cgca_as, cgca_ds, cgca_av, cgca_dv, cgca_art, cgca_drt
!  USES
!    cgca_m1co
!  USED BY
!    cgca
!  SOURCE

use cgca_m1co
implicit none

private
public :: cgca_as, cgca_ds, cgca_av, cgca_dv, cgca_art, cgca_drt

contains

!*roboend*


!*robodoc*s* cgca_m2alloc/cgca_as
!  NAME
!    cgca_as
!  SYNOPSIS

subroutine cgca_as( l1, u1, l2, u2, l3, u3, col1, cou1, col2, cou2,    &
 col3, props, coarray )

!  INPUTS

integer(kind=idef),intent(in) ::                                       &
 l1,   & ! Lower and uppper bounds of the coarray along
 u1,   & ! dimensions 1,2, and 3.
 l2,   &
 u2,   &
 l3,   &
 u3,   &
 col1, & ! Lower and uppper cobounds of the coarray along
 cou1, & ! codimensions 1, 2 and 3. Note the last upper
 col2, & ! cobound is not specified. In allocate it is
 cou2, & ! given as * to allow for arbitrary number of
 col3, & ! images at run time.
 props   ! Number of cell state properties to use.

integer(kind=iarr),allocatable,intent(inout) ::                        &
 coarray(:,:,:,:)[:,:,:]  ! coarray to allocate

!  SIDE EFFECTS
!    coarray becomes allocated with all values assigned
!    cgca_liquid_state.
!
!    Global array variables cgca_slcob, cgca_sucob are assigned values.
!  DESCRIPTION
!    This routine allocates a 4D arrays on each image. The first 3
!    dimensions define a cell, and the last, 4th dimension defines the
!    number or the cell state. The routine increases extent by 2
!    in each direction, thus creating space for storing halos.
!    No check for the validity of the coarray dimension is made here.
!    The user should make sure the coarray dimension values passed
!    to the routine make sense.
!  NOTES
!    We want the coarray to have halos to exchange data between
!    processors. So need to increase the extent by 2 in each dimension.
!  USES
!    none 
!  USED BY
!    cgca_m2alloc
!  SOURCE

integer, parameter :: halo=1
integer :: errstat = 0

if ( .not. allocated(coarray) ) allocate( coarray(                     &
  l1-halo:u1+halo, l2-halo:u2+halo, l3-halo:u3+halo, props)            &
  [col1:cou1, col2:cou2, col3:*], source=cgca_liquid_state,            &
    stat=errstat)
if ( errstat .ne. 0 ) then
  write (*,'(2(a,i0))') 'ERROR: cgca_m2alloc/cgca_as:&
    & allocate( coarray ), image: ', this_image(), " err. code: ",     &
    errstat
  error stop
end if

! Assign the cobounds to the global variables for use by other routines
cgca_slcob = lcobound( coarray )
cgca_sucob = ucobound( coarray )

end subroutine cgca_as

!*roboend*


!*robodoc*s* cgca_m2alloc/cgca_ds
!  NAME
!    cgca_ds
!  SYNOPSIS

subroutine cgca_ds(coarray)

!  INPUT

integer(kind=iarr),allocatable,intent(inout) :: coarray(:,:,:,:)[:,:,:]

!  SIDE EFFECTS
!    coarray becomes deallocated
!  DESCRIPTION
! This routine deallocates a 3D coarray. It first checks whether the
! array is allocated. If the array is not allocated, a warning is
! issued, but the program is *not* terminated.
!  USES
!    none
!  USED BY
!    cgca_m2alloc
!  SOURCE

integer :: errstat

errstat=0

if (allocated(coarray)) then
 deallocate(coarray,stat=errstat)
 if (errstat .ne. 0) then
  write (*,'(a,i0)') "ERROR: cgca_ds: image: ", this_image()
  write (*,'(a)')    "ERROR: cgca_ds: cannot deallocate coarray"
  write (*,'(a,i0)') "ERROR: cgca_ds: error code: ", errstat
  error stop
 end if
else
 write (*,'(a,i0,a)') "WARNING: cgca_ds: image: ", this_image(), &
  ", coarray is not allocated, cannot deallocate"
end if

end subroutine cgca_ds

!*roboend*


!*robodoc*s* cgca_m2alloc/cgca_av
!  NAME
!    cgca_av
!  SYNOPSIS

subroutine cgca_av(l,u,col1,cou1,col2,cou2,col3,coarray)

!  DESCRIPTION
!    This routine allocates a 1D coarray of length l:u.
!    Coarray variable "coarray" becomes allocated, with
!    all values assigned to zero.
!  INPUTS

integer(kind=idef),intent(in) :: l,u,col1,cou1,col2,cou2,col3
integer(kind=ilrg),allocatable,intent(inout) :: coarray(:)[:,:,:]

!  SIDE EFFECTS
!    None
!  USES
!    None
!  USED BY
!    cgca_m2alloc
!  SOURCE

integer :: errstat

errstat=0

if (.not. allocated(coarray))                          &
  allocate(coarray(l:u) [col1:cou1,col2:cou2,col3:*],  &
            source = 0_ilrg, stat=errstat)

if (errstat .ne. 0) then
 write (*,'(a,i0)') "ERROR: cgca_av: image: ", this_image()
 write (*,'(a)') "ERROR: cgca_av: cannot allocate coarray"
 write (*,'(a,i0)') "ERROR: cgca_av: error code: ", errstat
 error stop
end if

end subroutine cgca_av

!*roboend*


!*robodoc*s* cgca_m2alloc/cgca_dv
!  NAME
!    cgca_dv
!  SYNOPSIS

subroutine cgca_dv(coarray)

!  INPUT

integer(kind=ilrg),allocatable,intent(inout) :: coarray(:)[:,:,:]

!  SIDE EFFECTS
!    coarray becomes deallocated
!  DESCRIPTION
!    deallocate volume coarray
!  USES
!    none
!  USED BY
!    cgca_m2alloc
!  SOURCE

integer :: errstat

errstat=0

if (allocated(coarray)) then
 deallocate(coarray,stat=errstat)
 if (errstat .ne. 0) then
  write (*,'(a,i0)') "ERROR: cgca_dv: image: ", this_image()
  write (*,'(a)') "ERROR: cgca_dv: cannot deallocate coarray"
  write (*,'(a,i0)') "ERROR: cgca_dv: error code: ", errstat
  error stop
 end if
end if

! if coarray is not allocated, do nothing

end subroutine cgca_dv

!*roboend*


!*robodoc*s* cgca_m2alloc/cgca_art
!  NAME
!    cgca_art
!  SYNOPSIS

subroutine cgca_art(l,u,col1,cou1,col2,cou2,col3,coarray)

!  INPUTS
 
integer(kind=idef),intent(in) :: l,u,col1,cou1,col2,cou2,col3
real(kind=rdef),allocatable,intent(inout) :: coarray(:,:,:)[:,:,:]

!  SIDE EFFECTS
!    coarray becomes allocated
!  DESCRIPTION
! Allocate rotation tensor array. This is an array (l:u,3,3) defined on
! every image. Note that the first index is the grain number and
! the next two are the rotation tensor. So that the array element
! (153,3,1) is the rotation tensor component R31 for grain number 153.
! (87,:,:) is a 3x3 matrix defining the full (non-symmetric, but
! orthogonal!!!) rotation tensor for grain 87.
!  NOTES
!    This routine must be called prior to calling cgca_rt.
!  USES
!    none
!  USED BY
!    cgca_m2alloc
!  SOURCE

integer :: errstat

errstat=0

if (.not. allocated(coarray))                              &
  allocate(coarray(l:u,3,3) [col1:cou1,col2:cou2,col3:*],  &
            source = 0.0_rdef, stat=errstat)

if (errstat .ne. 0) then
  write (*,'(a,i0)') "ERROR: cgca_art: image: ", this_image()
  write (*,'(a)') "ERROR: cgca_art: cannot allocate coarray"
  error stop
end if

end subroutine cgca_art

!*roboend*


!*robodoc*s* cgca_m2alloc/cgca_drt
!  NAME
!    cgca_drt
!  SYNOPSIS

subroutine cgca_drt(coarray)

!  INPUT

real(kind=rdef),allocatable,intent(inout) :: coarray(:,:,:)[:,:,:]

!  SIDE EFFECTS
!    coarray becomes deallocated
!  DESCRIPTION
!    Deallocate rotation tensor array.
!  USES
!    none
!  USED BY
!    cgca_m2alloc
!  SOURCE

integer :: errstat

errstat=0

if (allocated(coarray)) then
  deallocate(coarray,stat=errstat)
  if (errstat .ne. 0) then
    write (*,'(a,i0)') "ERROR: cgca_drt: image: ", this_image()
    write (*,'(a)') "ERROR: cgca_drt: cannot deallocate coarray"
    error stop
  end if
end if

! if coarray is not allocated, do nothing

end subroutine cgca_drt

!*roboend*

end module cgca_m2alloc
