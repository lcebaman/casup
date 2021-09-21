!*robodoc*f* cgca_m3pfem/m3pfem_sm1
!  NAME
!    m3pfem_sm1
!  SYNOPSIS

!$Id: m3pfem_sm1.f90 380 2017-03-22 11:03:09Z mexas $

submodule ( cgca_m3pfem ) m3pfem_sm1

!  DESCRIPTION
!    Submodule with routines using collectives.
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!   cgca_pfem_map, cgca_pfem_lcentr_dump
!  USES
!    All variables and parameters of module cgca_m3pfem by host
!    association
!  USED BY
!    The host module cgca_m3pfem
!  SOURCE
!*roboend*

contains

!*robodoc*s* m3pfem_sm1/cgca_pfem_map
!  NAME
!    cgca_pfem_map
!  SYNOPSIS

!module procedure cgca_pfem_map
  module subroutine cgca_pfem_map( origin, rot, bcol, bcou )
    real( kind=rdef ), intent( in ) ::                                 &
      origin(3),        & ! origin of the "box" cs, in FE cs
      rot(3,3),         & ! rotation tensor *from* FE cs *to* CA cs
      bcol(3),          & ! lower phys. coords of the coarray on image
      bcou(3)             ! upper phys. coords of the coarray on image


!  INPUTS
!    See interface in the host module cgca_m3pfem.
!  SIDE EFFECTS
!    Array lcentr is changed.
!  DESCRIPTION
!    This routine reads centroids of all elements, in FE coord.
!    system, from all MPI processes and adds those with centroids
!    within its CA "box" to its lcentr array.
!  NOTES
!    This routine is an alternative to cgca_pfem_cenc.
!    While cgca_pfem_cenc implements all-to-all algorithm, this
!    routine uses collectives.
!    This routine must be called only after coarray
!    cgca_pfem_centroid_tmp has been established on all images.
!    This routine *reads* coarrays on other images, hence
!    sync must be used before calling this routine.
!    However, the routine *does not* change coarrays, only reads.
!    So no syncs are required inside this routine, as it constitutes
!    a single execution segment.
!    This routine uses CO_MAX and CO_SUM collective.
!    This routine allocates *large* tmp arrays on every image.
!    The array size is equal or even bigger than the number of FE
!    in the *whole* model, more precisely the array is allocated
!    on every image as
!    ( 5, <max no. of FE on any image>*num_images() ).
!    Hence this routine might give OOM for large models.
!    In that case fall back to cgca_pfem_cenc. 
!  USES
!    lcentr via host association.
!  SOURCE

! Initial length of lcentr array. A good choice will reduce the number
! of deallocate/allocate and will use the memory better.
integer, parameter :: lclenini = 100

real( kind=cgca_pfem_iwp ), allocatable :: tmp(:,:)

! Centroid coords in CA cs
real( kind=cgca_pfem_iwp ) :: cen_ca(3) ! 3D case only

! Temp array to expand/contract lcentr
type( mcen ), allocatable :: lctmp(:)

integer( kind=idef ) :: img, nimgs, maxfe, pos_start, pos_end, lclen,  &
  lcel, j, ctmpsize

integer :: errstat

!*********************************************************************72
! First executable statement
  img = this_image()
nimgs = num_images()

! Calculate the max number of FE stored on this image, i.e. nels_pp.
maxfe = size( cgca_pfem_centroid_tmp%r, dim=2 )

! Save it in a separate var
ctmpsize = maxfe

! Calculate the max number of FE stored on any image.
! Use CO_MAX collective. RESULT_IMAGE is not used.
! Hence the result is assigned to maxfe on all images.
call co_max( maxfe )

!write (*,*) "DEBUG, after co_max, img:", img, " maxfe:", maxfe
 
! Allocate tmp array of length ( maxfe * nimgs )
! 5 real values are stored per each FE:
! 1 - image number, cast to real
! 2 - FE number, cast to real.
! 3-5 - coordinates of the centroid of this FE
! NOTE! Important to set to zero initially, either on allocatation,
! or later, but before use. The following algorithm relies on the
! fact that tmp is zero initially on all images.
allocate( tmp( maxfe * nimgs, 5 ), source = 0.0_cgca_pfem_iwp,         &
          stat=errstat )
if ( errstat .ne. 0 ) then
  write (*,'(2(a,i0))') "ERROR: m3pfem_sm1/cgca_pfem_map: img: ", img, &
    ", allocate( tmp ), stat: ", errstat
  error stop
end if

! Write values in correct places in array tmp on this image.
! Use this_image() as the offset.
pos_start = (img - 1) * maxfe + 1
pos_end = pos_start + ctmpsize - 1

! Write image number
tmp( pos_start : pos_end, 1 ) = real( img, kind=cgca_pfem_iwp ) 

! Write element number
tmp( pos_start : pos_end, 2 ) =                                        &
  real( (/ (j, j = 1, ctmpsize) /), kind=cgca_pfem_iwp ) 

! Write centroid coord
tmp( pos_start : pos_end, 3:5 ) = &
  transpose( cgca_pfem_centroid_tmp%r(:,:) )

! Calculate the sum of tmp arrays over all images.
! Because each image wrote its data in a unique location,
! the sum will just produce the tmp array with data from all images.
! Then this tmp array is delivered back to all images.
! Since RESULT_IMAGE is not used, the result is assigned to tmp
! on all images.
call co_sum( tmp )

! Now each image searches through the whole of tmp array and
! adds all elements with centroids inside its CA box, to its local
! (private) lcentr array.

! Allocate lcentr to the initial guess size.
lclen = lclenini
allocate( lcentr( lclen ), stat=errstat )
if ( errstat .ne. 0 ) then
  write (*,'(2(a,i0))') "ERROR: m3pfem_sm1/cgca_pfem_map: img: ", img, &
    ", allocate( lcentr ), stat: ", errstat
  error stop
end if

! There are no elements yet in lcentr array
lcel = 0

! Loop over all elements in tmp array
elements: do j = 1, size( tmp, dim=1 )

  ! Convert centroid coordinates from FE cs to CA cs.
  ! tmp( 3:5 , j ) - take finite element j, and all centroid
  ! coordinates for it.
  cen_ca = matmul( rot, tmp( j, 3:5 ) - origin )

  ! Check whether CA cs centroid is within the box.
  ! If all CA cs centroid coordinates are greater or equal to
  ! the lower bound of the box, and all of them are also
  ! less of equal to the upper bound of the box, then the centroid
  ! is inside. Then add the new entry.
  inside: if ( all( cen_ca .ge. bcol ) .and.                           &
               all( cen_ca .le. bcou ) ) then

    ! Skip zero elements
    if ( int( tmp(j,1), kind=idef ) .eq. 0_idef .or.                   &
         int( tmp(j,2), kind=idef ) .eq. 0_idef ) cycle elements 

    ! Increment the number of elements
    lcel = lcel + 1

    ! Expand the array if there is no space left to add the new entry.
    expand: if ( lclen .lt. lcel ) then 

      ! Double the length of the array
      lclen = 2 * lclen

      ! Allocate a temp array of this length
      allocate( lctmp( lclen ), stat=errstat )
      if ( errstat .ne. 0 ) then
        write (*,'(2(a,i0))') "ERROR: m3pfem_sm1/cgca_pfem_map: img: ",&
          img, ", allocate( lctmp ) 1, stat: ", errstat
        error stop
      end if

      ! copy lcentr into the beginning of lctmp
      lctmp( 1:size( lcentr ) ) = lcentr

      ! move allocation from the temp array back to lcentr
      call move_alloc( lctmp, lcentr )

    end if expand

    ! Add new entry
    lcentr( lcel ) = mcen( int( tmp(j,1), kind=idef ),                 &
      int( tmp(j,2), kind=idef ), cen_ca )

  end if inside

end do elements

! Can now deallocate tmp
deallocate( tmp, stat=errstat )
if ( errstat .ne. 0 ) then
  write (*,'(2(a,i0))') "ERROR: m3pfem_sm1/cgca_pfem_map: img: ", img, &
    ", deallocate( tmp ), stat: ", errstat
  error stop
end if

! Trim lcentr if it is longer than the number of elements
if ( lclen .gt. lcel ) then

  ! Allocate temp array to the number of elements
  allocate( lctmp( lcel ), stat=errstat )
  if ( errstat .ne. 0 ) then
    write (*,'(2(a,i0))') "ERROR: m3pfem_sm1/cgca_pfem_map: img: ",    &
      img, ", allocate( lctmp ) 2, stat: ", errstat
    error stop
  end if

  ! Copy lcentr elements to the temp array
  lctmp = lcentr( 1 : lcel )

  ! move allocation from lctmp back to lcentr
  call move_alloc( lctmp, lcentr )
end if

!end procedure cgca_pfem_map
end subroutine cgca_pfem_map

!*roboend*


!*robodoc*s* m3pfem_sm1/cgca_pfem_lcentr_dump
!  NAME
!    cgca_pfem_lcentr_dump
!  SYNOPSIS

module procedure cgca_pfem_lcentr_dump

!  SIDE EFFECTS
!    Dump lcentr from this image to OUTPUT_UNIT.
!  DESCRIPTION
!    This routine is used for debugging, i.e. to check that
!    cgca_pfem_cenc and cgca_pfem_map produce identical output,
!    as they should.
!  SOURCE

integer :: i, img

img = this_image()

!write (*,*) "DEBUG: img:", this_image(), "size( lcentr ):", size(lcentr)

!if ( img .le. 50 ) then
  do i = lbound( lcentr, dim=1 ) , ubound( lcentr, dim=1 )
    write (*,"(2(a,i0),tr1,i0,3(tr1,es9.2))") "DEBUG: img: ", img,     &
      " lcentr: ", lcentr(i)
  end do
!end if

end procedure cgca_pfem_lcentr_dump

!*roboend*

end submodule m3pfem_sm1
