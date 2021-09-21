!*robodoc*m* CGPACK/cgca_m3gbf
!  NAME
!    cgca_m3gbf
!  SYNOPSIS

!$Id: cgca_m3gbf.f90 529 2018-03-26 11:25:45Z mexas $

module cgca_m3gbf

!  DESCRIPTION
!    Module dealing with grain boundary fractures.
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    cgca_gbf1p, cgca_gbf1f
!  USES
!    cgca_m1co, cgca_m2glm
!  USED BY
!  SOURCE

use cgca_m1co
use cgca_m2glm
implicit none

private
public :: cgca_gbf1p, cgca_gbf1f

contains

!*roboend*


!*robodoc*s* cgca_m3gbf/cgca_gbf1p
!  NAME
!    cgca_gbf1p
!  SYNOPSIS

subroutine cgca_gbf1p(coarray)

!  INPUTS

integer(kind=iarr),allocatable,intent(inout) :: coarray(:,:,:,:)[:,:,:]

!  SIDE EFFECTS
!    state of coarray changes
!  DESCRIPTION
!    This routine does a single iteration of grain boundary fracture
!    propagation assuming periodic boundary conditions.
!  NOTE
!    Use only with periodic BC. For fixed BC use cgca_gbf1f.
!  SOURCE

integer(kind=iarr),allocatable,save :: array(:,:,:)

integer(kind=idef) :: &
  range1(3), range2(3), range3(3), &
  lbv(4) ,& ! lower bounds of the complete (plus virtual) coarray
  ubv(4) ,& ! upper bounds of the complete (plus virtual) coarray
  lbr(4) ,& ! lower bounds of the "real" coarray, lower virtual+1
  ubr(4) ,& ! upper bounds of the "real" coarray, upper virtual-1
  x1     ,& ! local coordinates in an array, which are also
  x2     ,& ! do loop counters
  x3

integer :: thisimage, errstat=0, nimages

! Do not check coarray for allocated, as this wastes time.
! Instead let the code fail if coarray is not allocated.

! use local vars to save time
thisimage = this_image()
nimages = num_images()

! determine the extents
lbv = lbound( coarray )
ubv = ubound( coarray )
lbr = lbv + 1
ubr = ubv - 1

! Allocate the temp array if not already allocated.
! The array has the SAVE attribute, as this routine
! is likely to be called many times.

! Initialise errstat to -1
errstat = -1
 
if (.not. allocated(array)) &
 allocate(array(lbv(1):ubv(1),lbv(2):ubv(2),lbv(3):ubv(3)),stat=errstat)
if (errstat.ne.0) then
 write (*,'(a,i0)') "ERROR: cgca_gbf1p: image ",thisimage
 write (*,'(a)') "ERROR: cgca_gbf1p: cannot allocate array"
 error stop
end if

! Copy coarray fracture state type into a local array
array = coarray(:,:,:,cgca_state_type_frac)

! scan across all cells
do x3=lbr(3),ubr(3)
do x2=lbr(2),ubr(2)
do x1=lbr(1),ubr(1)

 ! Analyse only live cells
 if ( array(x1,x2,x3) .ne. cgca_intact_state ) cycle

 ! set up ranges to save compute time
 range1 = (/ x1-1, x1, x1+1 /)
 range2 = (/ x2-1, x2, x2+1 /)
 range3 = (/ x3-1, x3, x3+1 /)

 ! If the cell
 !  (1) is on the grain boundary
 !  (2) has a fractured neighbour
 ! then it becomes cgca_gb_state_fractured.
 if ( any( coarray(range1,range2,range3,cgca_state_type_grain) .ne. &
           coarray(x1,x2,x3,cgca_state_type_grain) ) .and.          &
      any( array(range1,range2,range3) .ne. cgca_intact_state ) )   &
 then
  coarray(x1,x2,x3,cgca_state_type_frac) = cgca_gb_state_fractured
 end if
    
end do
end do
end do

end subroutine cgca_gbf1p

!*roboend*


!*robodoc*s* cgca_m3gbf/cgca_gbf1f
!  NAME
!    cgca_gbf1f
!  SYNOPSIS

subroutine cgca_gbf1f( coarray )

!  INPUTS

integer( kind=iarr ), allocatable, intent(inout) ::                    &
  coarray(:,:,:,:)[:,:,:]

!  OUTPUTS
!    coarray, as it's intent(INOUT)
!  SIDE EFFECTS
!    None
!  DESCRIPTION
!    This routine does a single iteration of grain boundary fracture
!    propagation assuming fixed boundary conditions.
!  NOTE
!    Use only with fixed BC. For periodic BC use cgca_gbf1p.
!  USES
!    cgca_lg
!  SOURCE

integer( kind=iarr ), allocatable, save :: array( :, :, : )
integer :: range1(3), range2(3), range3(3)

integer(kind=idef) :: &
  lbv(4),    & ! lower bounds of the complete (plus virtual) coarray
  ubv(4),    & ! upper bounds of the complete (plus virtual) coarray
  lbr(4),    & ! lower bounds of the "real" coarray, lower virtual+1
  ubr(4),    & ! upper bounds of the "real" coarray, upper virtual-1
  x1,x2,x3,  & ! local coordinates in an array
  super(3),  & ! global (super) coordinates of a cell
  imgpos(3), & ! image position in the image grid
  local(3),  & ! local coordinates of a cell
  ubsuper(3),& ! upper bounds of the super array
  frnei,     & ! number of fractured neighbours
  i            ! loop counter

integer :: thisimage, errstat=0, nimages

! Do not check coarray for allocated, as this wastes time.
! Instead let the code fail if coarray is not allocated.

! use local vars to save time
thisimage = this_image()
  nimages = num_images()
   imgpos = this_image( coarray )

! determine the extents
lbv = lbound( coarray )
ubv = ubound( coarray )
lbr = lbv + 1
ubr = ubv - 1
ubsuper = (ubr(1:3) - lbr(1:3) + 1) * nimages

! Allocate the temp array if not already allocated.
! The array has the SAVE attribute, as this routine
! is likely to be called many times.

if ( .not. allocated(array) ) then
  allocate( array( lbv(1):ubv(1) , lbv(2):ubv(2) , lbv(3):ubv(3) ),    &
            stat=errstat )
  if ( errstat .ne. 0 ) then
    write (*,'(2(a,i0))') "ERROR: cgca_gbf1f/cgca_m3gbf: image: ",     &
      thisimage, "allocate( array ), stat: ", errstat
    error stop
  end if
end if

! Copy coarray fracture state type into a local array
array = coarray(:,:,:,cgca_state_type_frac)

! scan across all cells
do x3=lbr(3),ubr(3)
do x2=lbr(2),ubr(2)
do x1=lbr(1),ubr(1)

 ! Analyse only live cells
 if ( array(x1,x2,x3) .ne. cgca_intact_state ) cycle

 ! Skip cells adjacent to halo cells
 local = (/ x1, x2, x3 /)
 call cgca_lg(imgpos,local,coarray,super)
 if ( any( super .eq. 1) .or. any( super .eq. ubsuper) ) cycle

 ! set up ranges to save compute time
 range1 = (/ x1-1, x1, x1+1 /)
 range2 = (/ x2-1, x2, x2+1 /)
 range3 = (/ x3-1, x3, x3+1 /)

 ! count fractured neighbours, only cleavage edges and fractured GB
 frnei = 0
 ! first all cleavage edge states
 do i=1,size(cgca_clvg_states_edge)
  frnei = frnei + count(                                            &
   array(range1,range2,range3) .eq. cgca_clvg_states_edge(i) )
 end do
 ! then all GB fractured states
 frnei = frnei + count(                                             &
  array(range1,range2,range3) .eq. cgca_gb_state_fractured )

 ! If the cell
 !  (1) is on the grain boundary
 !  (2) has 5 fractured neighbours
 ! then it becomes cgca_gb_state_fractured.
 if ( any( coarray(range1,range2,range3,cgca_state_type_grain) .ne. &
           coarray(x1,x2,x3,cgca_state_type_grain) ) .and.          &
      (frnei .ge. 5) ) then
  coarray(x1,x2,x3,cgca_state_type_frac) = cgca_gb_state_fractured
 end if
    
end do
end do
end do

end subroutine cgca_gbf1f

!*roboend*

end module cgca_m3gbf
