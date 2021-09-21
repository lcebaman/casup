!*robodoc*m* CGPACK/cgca_m3pfem
!  NAME
!    cgca_m3pfem
!  SYNOPSIS

!$Id: cgca_m3pfem.f90 380 2017-03-22 11:03:09Z mexas $

module cgca_m3pfem

!  DESCRIPTION
!    Module dealing with interfacing CGPACK with ParaFEM.
!  AUTHOR
!    Anton Shterenlikht, Luis Cebamanos
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    Public coarray variables of derived types:
!    cgca_pfem_centroid_tmp,
!    cgca_pfem_integrity,
!    cgca_pfem_stress.
!
!    Public *local*, non-coarray, variable:
!    cgca_pfem_enew.
!
!    Private *local*, non-coarray, variables:
!    lcentr.
!
!    Public routines:
!    cgca_pfem_boxin, cgca_pfem_cellin,
!    cgca_pfem_cenc, cgca_pfem_cendmp, cgca_pfem_ctalloc,
!    cgca_pfem_ctdalloc, cgca_pfem_ealloc, cgca_pfem_edalloc,
!    cgca_pfem_intcalc1, cgca_pfem_integalloc, cgca_pfem_integdalloc,
!    cgca_pfem_lcentr_dump (in submodule m3pfem_sm1),
!    cgca_pfem_map (in submodule m3pfem_sm1),
!    cgca_pfem_partin,
!    cgca_pfem_salloc, cgca_pfem_sdalloc, cgca_pfem_sdmp,
!    cgca_pfem_simg, cgca_pfem_uym, cgca_pfem_wholein
!  USES
!    Modules cgca_m1co, cgca_m2lnklst, cgca_m2geom
!  USED BY
!    end user?
!  SOURCE

use :: cgca_m1co
use :: cgca_m2lnklst, only : cgca_lnklst_tpayld, cgca_lnklst_node,     &
  cgca_inithead, cgca_addhead, cgca_lstdmp, cgca_rmhead
use :: cgca_m2geom, only : cgca_boxsplit

implicit none

private
public ::                                                              &
! routines
  cgca_pfem_boxin, cgca_pfem_cellin, &
  cgca_pfem_cenc, cgca_pfem_cendmp, cgca_pfem_ctalloc,                 &
  cgca_pfem_ctdalloc, cgca_pfem_ealloc, cgca_pfem_edalloc,             &
  cgca_pfem_integalloc, cgca_pfem_integdalloc, cgca_pfem_intcalc1,     &
  cgca_pfem_lcentr_dump, & ! in submodule m3pfem_sm1
  cgca_pfem_map,         & ! in submodule m3pfem_sm1
  cgca_pfem_partin, &
  cgca_pfem_salloc,  cgca_pfem_sdalloc, cgca_pfem_sdmp, cgca_pfem_simg,&
  cgca_pfem_uym, cgca_pfem_wholein,                         &
! variables 
  cgca_pfem_centroid_tmp, cgca_pfem_enew, cgca_pfem_integrity,         &
  cgca_pfem_stress

! corresponds to typical double precision real.
integer, parameter :: cgca_pfem_iwp = selected_real_kind(15,300)

interface

  module subroutine cgca_pfem_lcentr_dump
  end subroutine cgca_pfem_lcentr_dump

  module subroutine cgca_pfem_map( origin, rot, bcol, bcou )
    real( kind=rdef ), intent( in ) ::                                 &
      origin(3),        & ! origin of the "box" cs, in FE cs
      rot(3,3),         & ! rotation tensor *from* FE cs *to* CA cs
      bcol(3),          & ! lower phys. coords of the coarray on image
      bcou(3)             ! upper phys. coords of the coarray on image
  end subroutine cgca_pfem_map

end interface

!*roboend*


!*robodoc*d* cgca_m3pfem/lcentr
!  NAME
!    lcentr
!  SYNOPSIS

type mcen
  integer( kind=idef ) :: image
  integer( kind=idef ) :: elnum
  real( kind=cgca_pfem_iwp ) :: centr(3)
end type mcen
type( mcen ), allocatable :: lcentr(:)

!  DESCRIPTION
!    A *private* *local* allocatable array of derived type with 3
!    components: (1) image number (2) the local element number on that
!    image and (3) centroid coordinates in CA CS.
!    Each entry in this array
!    corresponds to an FE with centroid coordinates within the coarray
!    "box" on this image.
!
!    Assumption!! This is a 3D problems, so the centroid is
!    defined by 3 coordinates, hence centr(3).
!
!    MCEN stands for Mixed CENtroid data type.
!    LCENTR stands for *Local* array of CENTRoids.
!  NOTE
!    This is *private* array, hence the name does not start
!    with "cgca_pfem".
!  USED BY
!    Many routines of this module.
!*roboend*


!*robodoc*d* cgca_m3pfem/cgca_pfem_centroid_tmp
!  NAME
!    cgca_pfem_centroid_tmp
!  SYNOPSIS
  
type rca
  real( kind=cgca_pfem_iwp ), allocatable :: r(:,:)
end type rca
type( rca ) :: cgca_pfem_centroid_tmp[*]

!  DESCRIPTION
!    RCA stands for Rugged CoArray.
!    cgca_pfem_centroid_tmp is a temporary scalar *coarray* of derived
!    type with allocatable array component, storing centroids of ParaFEM
!    finite elements, in FE coord. system, on this image.
!    The array might be of different length on different images,
!    so have to use an allocatable component of a coarray variable
!    of derived type.
!  USED BY
!    routines of this module + end user
!*roboend*


!*robodoc*d* cgca_m3pfem/cgca_pfem_stress
!  NAME
!    cgca_pfem_stress
!  SYNOPSIS

type type_stress
  real( kind=cgca_pfem_iwp ), allocatable :: stress(:,:,:)
end type type_stress
type( type_stress ) :: cgca_pfem_stress[*]

!  DESCRIPTION
!    This is a coarray with a single allocatable array component,
!    to store all stress components for all integration points
!    for all elements on an image. Have to use a derived type
!    because cgca_pfem_stress%stress can be allocated to different
!    length on different images.
!    This data will be read by all images.
!*roboend*

  
!*robodoc*d* cgca_m3pfem/cgca_pfem_integrity
!  NAME
!    cgca_pfem_integrity
!  SYNOPSIS

type cgca_pfem_integ_type
 real( kind=rdef ), allocatable :: i(:)
end type cgca_pfem_integ_type
type( cgca_pfem_integ_type ) :: cgca_pfem_integrity[*]

!  DESCRIPTION
!    A derived type is needed because the length of the integrity
!    array will differ from image to image. So this is a scalar coarray
!    of derived type with a single component: allocatable array of
!    integrity, i. i=1 means to damage, i=0 means no remaining load
!    bearing capacity.
!    This data will be used to update the Young's modulus
!  NOTE
!    Set i to 1 on allocation to avoid problems later.
!    The reason is that in cases when some FE are not
!    connected to CA, the integrity of these FE will never be
!    set of changed. So setting i to 1 on allocation is fool proof.
!  USED BY
!    cgca_pfem_uym + end user?
!*roboend*
  
  
!*robodoc*d* cgca_m3pfem/cgca_pfem_enew
!  NAME
!    cgca_pfem_enew
!  SYNOPSIS
  
real( kind=cgca_pfem_iwp ), allocatable :: cgca_pfem_enew(:,:)

!  DESCRIPTION
!    Naming: E New as in new Young's modulus. This *local* array
!    stores Young's moduli for each integration point of each
!    FE on this image.
!  USED BY
!    cgca_pfem_uym + end user
!*roboend*


contains


!*robodoc*s* cgca_m3pfem/cgca_pfem_integalloc
!  NAME
!    cgca_pfem_integalloc
!  SYNOPSIS

subroutine cgca_pfem_integalloc( nels_pp )

!  INPUT
!    nels_pp - elements per MPI process (per image).

integer, intent(in) :: nels_pp

!  SIDE EFFECTS
!    Allocatable array component cgca_pfem_integrity%i becomes allocated
!  DESCRIPTION
!    This routine allocates cgca_pfem_integrity%i on this image.
!    This is a *local*, non-coarray, array. Hence this routine can be
!    called by any or all images. It should be called by all images,
!    of course.
!
!    The array is allocated with the length equal to the number FE
!    stored on *that* image.
!
!    Must set i to 1, to take care of cases when some FE are not
!    linked to CA. integrity for such FESuch FE which are 
!  USES
!    cgca_pfem_integrity via host association.
!  USED BY
!    end user?
!  SOURCE

integer :: errstat=0

allocate( cgca_pfem_integrity%i( nels_pp ),  source=1.0_rdef,          &
          stat=errstat )
if ( errstat .ne. 0 ) then
  write (*,"(a,i0)")                                                   &
    "ERROR: cgca_pfem_integalloc: allocate( cgca_pfem_integrity%i )",  &
    errstat
  error stop
end if

end subroutine cgca_pfem_integalloc

!*roboend*
 

!*robodoc*s* cgca_m3pfem/cgca_pfem_integdalloc
!  NAME
!    cgca_pfem_integdalloc
!  SYNOPSIS

subroutine cgca_pfem_integdalloc

!  SIDE EFFECTS
!    Allocatable array component of cgca_pfem_integrity coarray becomes
!    deallocated.
!  DESCRIPTION
!    This routine deallocates allocatable array component of integrity:
!    cgca_pfem_integrity%i.
!  USES
!    cgca_pfem_integrity via host association
!  USED BY
!    end user?
!  SOURCE

integer :: errstat=0

deallocate( cgca_pfem_integrity%i, stat=errstat )
if ( errstat .ne. 0 ) then
  write (*,"(a,i0)")                                                   &
   "ERROR: cgca_pfem_integalloc: deallocate( cgca_pfem_integrity%i )", &
   errstat
  error stop
end if

end subroutine cgca_pfem_integdalloc

!*roboend*


!*robodoc*s* cgca_m3pfem/cgca_pfem_ealloc
!  NAME
!    cgca_pfem_ealloc
!  SYNOPSIS

subroutine cgca_pfem_ealloc( nip, nels_pp )
    
!  INPUTS
!    nip - integer, number of integration points
!    nels_pp - elements per MPI process (per image).

integer, intent( in ) :: nip, nels_pp
    
!  SIDE EFFECTS
!    Allocatable *local* array enew becomes allocated
!  DESCRIPTION
!    This routine allocates an allocatable *local* array.
!    The allocatable array stores the Young's modulus
!    per FE element and integration point, for all FE
!    that are stored on this image.
!  USED BY
!    end user?
!  SOURCE

integer :: errstat=0

allocate( cgca_pfem_enew( nip, nels_pp ), stat=errstat )
if ( errstat .ne. 0 ) then
  write (*,"(a,i0)")                                                   &
   "ERROR: cgca_pfem_ealloc: allocate( cgca_pfem_enew ), err. status", &
    errstat
  error stop
end if

end subroutine cgca_pfem_ealloc

!*roboend*

  
!*robodoc*s* cgca_m3pfem/cgca_pfem_edalloc
!  NAME
!    cgca_pfem_edalloc
!  SYNOPSIS

subroutine cgca_pfem_edalloc
    
!  SIDE EFFECTS
!    Allocatable *local* array cgca_pfem_enew becomes deallocated.
!  DESCRIPTION
!    This routine deallocates an allocatable *local* array used to
!    store the Young's modulus per FE element and integration point
!  USES
!    cgca_pfem_enew via host association.
!  USED BY
!    end user?
!  SOURCE

integer :: errstat=0

deallocate( cgca_pfem_enew, stat=errstat )
if ( errstat .ne. 0 ) error stop                                       &
  "ERROR: cgca_pfem_ealloc: deallocate( cgca_pfem_enew )"
    
end subroutine cgca_pfem_edalloc

!*roboend*


!*robodoc*s* cgca_m3pfem/cgca_pfem_uym
!  NAME
!    cgca_pfem_uym
!  SYNOPSIS

subroutine cgca_pfem_uym( e_orig, nels_pp )

!  INPUTS
!     e_orig - *real* is the original Young's modulus.
!    For now assume a single value, i.e. all int points
!    have identical original value.
!     nels_pp - number of FEs for this image.

real( kind=cgca_pfem_iwp ), intent(in) :: e_orig
integer, intent(in) :: nels_pp

!  SIDE EFFECTS
!    The Young's modulus gets updated with integrity
!  DESCRIPTION
!    UYM stands for Update Young's Modulus
!    This routine updates the value of the Young's modulus, e,
!    e = e_original * integrity.
!    Integrity - integer, cell integrity (from 0.0 to 1.0)
!  NOTES
!    Purely local routine, no coarray operations.
!    It seems the Young's modulus of 0 causes instability.
!    So don't let it get to 0, use a small factor instead, e.g. 1.0e-3.
!  USES
!    cgca_pfem_enew, cgca_pfem_integrity, all via host association.
!  USED BY
!    end user?
!  SOURCE

real( kind=rdef ), parameter :: factor = 1.0e-3_rdef
integer :: fe

do fe = 1, nels_pp

 cgca_pfem_enew( : , fe ) = max( factor*e_orig,                        &
                             e_orig * cgca_pfem_integrity % i( fe ) )
end do

end subroutine cgca_pfem_uym

!*roboend*
  

!*robodoc*s* cgca_m3pfem/cgca_pfem_ctalloc
!  NAME
!    cgca_pfem_ctalloc
!  SYNOPSIS

subroutine cgca_pfem_ctalloc( ndim, nels_pp )

!  INPUTS

integer, intent( in ) :: ndim, nels_pp

!    ndim - integer, number of DOF per node.
!    nels_pp - elements per MPI process (per image).
!  SIDE EFFECTS
!    Allocatable array component cgca_pfem_centroid_tmp%r becomes
!    allocated.
!  DESCRIPTION
!    CTA stands for Centroids Temporary Allocate.
!    This routine allocates an allocatable array component of scalar
!    coarray cgca_pfem_centroid_tmp.
!    The allocatable array stores FE centroid coordinates
!    together with their numbers and MPI ranks where these are stored.
!  NOTES
!    The array component can of different length on different images.
!  USES
!  USED BY
!  SOURCE

integer :: errstat=0

allocate( cgca_pfem_centroid_tmp%r( ndim, nels_pp ),                   &
          source=0.0_cgca_pfem_iwp,                                    &
          stat=errstat )
if ( errstat .ne. 0 ) error stop                                       &
  "ERROR: cgca_pfem_ctalloc: allocate( cgca_pfem_centroid_tmp%r )"

end subroutine cgca_pfem_ctalloc

!*roboend*


!*robodoc*s* cgca_m3pfem/cgca_pfem_ctdalloc
!  NAME
!    cgca_pfem_ctdalloc
!  SYNOPSIS

subroutine cgca_pfem_ctdalloc

!  SIDE EFFECTS
!    Allocatable array component of cgca_pfem_centroid_tmp becomes
!    deallocate
!  DESCRIPTION
!    CTD stands for Centroids Temporary Allocate.
!    This routine deallocates an allocatable array component of coar.
!    This must be done only after all images copied the contents of
!    type( rca ) :: cgca_pfem_centroid_tmp[*] into their local,
!    *not* coarray centroid arrays.
!  USES
!  USED BY
!  SOURCE

integer :: errstat=0

deallocate( cgca_pfem_centroid_tmp%r, stat=errstat )
if ( errstat .ne. 0 ) error stop                                       &
  "ERROR: cgca_pfem_ctd: deallocate( cgca_pfem_centroid_tmp%r )"

end subroutine cgca_pfem_ctdalloc

!*roboend*


!*robodoc*s* cgca_m3pfem/cgca_pfem_cenc
!  NAME
!    cgca_pfem_cenc
!  SYNOPSIS

subroutine cgca_pfem_cenc( origin, rot, bcol, bcou )

!  INPUTS

real( kind=rdef ), intent( in ) ::                                     &
 origin(3),        & ! origin of the "box" cs, in FE cs
 rot(3,3),         & ! rotation tensor *from* FE cs *to* CA cs
 bcol(3),          & ! lower phys. coords of the coarray on image
 bcou(3)             ! upper phys. coords of the coarray on image

!  SIDE EFFECTS
!    Array lcentr is changed.
!  DESCRIPTION
!    CENC stands for CENtroids Collection.
!    This routine reads centroids of all elements, in FE coord.
!    system, from all MPI processes and adds those with centroids
!    within its CA "box" to its lcentr array.
!  NOTES
!    This routine must be called only after coarray
!    cgca_pfem_centroid_tmp has been established on all images.
!    This routine *reads* coarrays on other images, hence
!    sync must be used before calling this routine.
!    However, the routine *does not* change coarrays, only reads.
!    So no syncs are required inside this routine, as it constitutes
!    a single execution segment.
!    This routine uses all-to-all comm pattern. This might be
!    inefficient on large numbers of PEs. In this case on can use
!    cgca_pfem_map (in submodule m3pfem_sm1) instead, which
!    does the same calculation using collectives and large tmp arrays.
!  USES
!    lcentr via host association.
!  SOURCE

! initial length of lcentr array. A good choice will reduce the number
! of deallocate/allocate and will use the memory better.
integer, parameter :: lclenini = 100

integer :: errstat, i, j, nimgs, nelements, img_curr, ndims, rndint, &
 lclen, & ! current length of the lcentr array
 lcel     ! number of elements in lcentr array

! centroid coords in CA cs
real( kind=cgca_pfem_iwp ) :: cen_ca(3) ! 3D case only
real( kind=cgca_pfem_iwp ), allocatable :: tmp(:,:)
real :: rnd

! temp array to expand/contract lcentr
type( mcen ), allocatable :: lctmp(:)

nimgs = num_images()

! Allocate lcentr to the initial guess size.
lclen = lclenini
allocate( lcentr( lclen ), stat=errstat )
if ( errstat .ne. 0 ) then
  write (*,'(a,i0)')                                                   &
   "ERROR: cgca_pfem_cenc: allocate( lcentr ), error code: ", errstat
  error stop
end if

! There are no elements yet in lcentr array
lcel = 0

! Choose the first image at random
! It is assumed that the RND has been initialised by a call
! to cgca_irs earlier on. 
call random_number( rnd )   ! [ 0 .. 1 )
rndint = int( rnd*nimgs )+1 ! [ 1 .. nimgs ]

! loop over all images, starting at a randomly chosen image
images: do i=rndint, rndint+nimgs-1

  ! Get the current image number.
  ! If it's > nimgs, subtract nimgs
  img_curr = i
  if ( img_curr .gt. nimgs ) img_curr = img_curr - nimgs

  ! how many elements
      ndims = size( cgca_pfem_centroid_tmp[ img_curr ] % r, dim=1 )
  nelements = size( cgca_pfem_centroid_tmp[ img_curr ] % r, dim=2 )

  ! use a temp array to pull all centroids data in one call
  allocate( tmp( ndims, nelements ), source=0.0_cgca_pfem_iwp,         &
            stat=errstat )
  if ( errstat .ne. 0 )                                                &
    error stop "ERROR: cgca_m3pfem/cgca_pfem_cenc: allocate( tmp )"
  tmp = cgca_pfem_centroid_tmp[ img_curr ] % r

  ! loop over all elements on that image
  elements: do j = 1, nelements

    ! Convert centroid coordinates from FE cs to CA cs
    ! cgca_pfem_centroid_tmp[i]        - variable on image i
    ! cgca_pfem_centroid_tmp[i]%r      - component that is the centroids
    !                                    real array
    ! cgca_pfem_centroid_tmp[i]%r(:,j) - take finite element j, and all
    !                                    centroid coordinates for it.
    !
    ! old algorithm - lots of small remote calls:
    ! cen_ca = &
    !  matmul( rot, cgca_pfem_centroid_tmp[ img_curr ]%r(:,j) - origin )
    !
    cen_ca = matmul( rot, tmp(:,j) - origin )

    ! Check whether CA cs centroid is within the box.
    ! If all CA cs centroid coordinates are greater or equal to
    ! the lower bound of the box, and all of them are also
    ! less of equal to the upper bound of the box, then the centroid
    ! is inside. Then add the new entry.
    inside: if ( all( cen_ca .ge. bcol ) .and.                         &
                 all( cen_ca .le. bcou ) ) then

      ! Increment the number of elements
      lcel = lcel + 1

      ! Expand the array if there is no space left to add the new entry.
      expand: if ( lclen .lt. lcel ) then 

        ! Double the length of the array
        lclen = 2 * lclen

        ! Allocate a temp array of this length
        allocate( lctmp( lclen ), stat=errstat )
        if ( errstat .ne. 0 ) error stop                               &
          "ERROR: cgca_pfem_cenc: allocate( lctmp ) 1"

        ! copy lcentr into the beginning of lctmp
        lctmp( 1:size( lcentr ) ) = lcentr

        ! move allocation from the temp array back to lcentr
        call move_alloc( lctmp, lcentr )

      end if expand

      ! Add new entry
      lcentr( lcel ) = mcen( img_curr, j, cen_ca )

    end if inside

  end do elements

  deallocate( tmp, stat=errstat )
  if ( errstat .ne. 0 )                                                &
    error stop "ERROR: cgca_pfem_cenc: deallocate( tmp )"

end do images

! Trim lcentr if it is longer than the number of elements
if ( lclen .gt. lcel ) then

  ! Allocate temp array to the number of elements
  allocate( lctmp( lcel ), stat=errstat )
  if ( errstat .ne. 0 ) error stop                                     &
    "ERROR: cgca_pfem_cenc: allocate( lctmp ) 2"

  ! Copy lcentr elements to the temp array
  lctmp = lcentr( 1 : lcel )

  ! move allocation from lctmp back to lcentr
  call move_alloc( lctmp, lcentr )
end if

end subroutine cgca_pfem_cenc

!*roboend*


!*robodoc*s* cgca_m3pfem/cgca_pfem_cendmp
!  NAME
!    cgca_pfem_cendmp
!  SYNOPSIS

subroutine cgca_pfem_cendmp

!  SIDE EFFECTS
!    Dumps some data to stdout
!  DESCRIPTION
!    CENDMP stands for CENtroids array dump.
!    This routine dumps lcentr to stdout.
!  NOTES
!    Must call from all images.
!  SOURCE

integer :: i, img

img = this_image()

do i = 1, size( lcentr )
  write (*,"(3(a,i0),a,3(es10.2,tr1))")                                &
              "CA on img "       , img             ,                   &
              " <-> FE "         , lcentr(i)%elnum ,                   &
              " on img "         , lcentr(i)%image ,                   &
              " centr. in CA cs" , lcentr(i)%centr
end do

end subroutine cgca_pfem_cendmp

!*roboend*


!*robodoc*s* cgca_m3pfem/cgca_pfem_salloc
!  NAME
!    cgca_pfem_salloc
!  SYNOPSIS

subroutine cgca_pfem_salloc( nels_pp, intp, comp )

!  INPUTS
!    nels_pp - number of elements on this image
!    intp - number of integration points per element
!    comp - number of stress tensor components

integer, intent( in ) :: nels_pp, intp, comp

!  SIDE EFFECTS
!    Allocatable component array cgca_pfem_stress%stress becomes
!    allocated
!  DESCRIPTION
!    SALLOC stands for Allocate Stress tensor array.
!    This routine allocates an allocatable array component of coar.
!    The allocatable array stores all stress tensor components,
!    for all integration points on all elements on an image.
!  USES
!    cgca_pfem_iwp, host association
!  USED BY
!    end user
!  SOURCE

integer :: errstat=0

allocate( cgca_pfem_stress%stress( nels_pp, intp, comp ),              &
          source=0.0_cgca_pfem_iwp, stat=errstat )
if ( errstat .ne. 0 ) then
  write (*,"(a,i0)") "ERROR: cgca_pfem_salloc: allocate( &
    &cgca_pfem_stress%stress ), err. status: ", errstat
  error stop
end if

end subroutine cgca_pfem_salloc

!*roboend*


!*robodoc*s* cgca_m3pfem/cgca_pfem_sdalloc
!  NAME
!    cgca_pfem_sdalloc
!  SYNOPSIS

subroutine cgca_pfem_sdalloc

!  SIDE EFFECTS
!    allocatable array cgca_pfem_stress%stress become deallocated
!  DESCRIPTION
!    SDALLOC stands for Deallocate Stress tensor array.
!    This routine deallocates allocatable array component of coar.
!    This routine should be called only when the analysis is complete.
!    Any and every image can call this routine.
!  USES
!    cgca_pfem_stress%stress, host association
!  USED BY
!  SOURCE

integer :: errstat=0

deallocate( cgca_pfem_stress%stress, stat=errstat )
if ( errstat .ne. 0 ) then
  write (*,"(a,i0)") "ERROR: cgca_pfem_sdalloc: deallocate( &
    &cgca_pfem_stress%stress ), err. status: ", errstat
  error stop
end if

end subroutine cgca_pfem_sdalloc

!*roboend*


!*robodoc*s* cgca_m3pfem/cgca_pfem_sdmp
!  NAME
!    cgca_pfem_sdmp
!  SYNOPSIS

subroutine cgca_pfem_sdmp

!  SIDE EFFECTS
!    Dumps some data to stdout
!  DESCRIPTION
!    SDMP stands for Stress tensor dump.
!    This routine dumps stress tensors to stdout.
!  NOTES
!    Must call from all images.
!  SOURCE

integer :: img, nel, nintp, el, intp

  img = this_image()
  nel = size( cgca_pfem_stress%stress, dim=1 )
nintp = size( cgca_pfem_stress%stress, dim=2 )

do el = 1, nel
  do intp = 1, nintp
    write (*,*) "img", img, "FE", el, "int p.", intp, "stress",        &
                cgca_pfem_stress%stress( el, intp, : )
  end do
end do

end subroutine cgca_pfem_sdmp
!*roboend*


!*robodoc*s* cgca_m3pfem/cgca_pfem_simg
!  NAME
!    cgca_pfem_simg
!  SYNOPSIS

subroutine cgca_pfem_simg( simg )

!  OUTPUT
!    simg - mean stress tensor over all integration points on all
!    finite elements linked to CA on this image.
!    Note that I use CGPACK kind, because this var will be input to
!    a CGPACK routine.

real( kind=rdef ), intent(out) :: simg(3,3)

!  DESCRIPTION
!    SIMG stands for mean Stress on an Image.
!    The routine reads all stress tensors from all integration
!    points for all elements which are linked to CA on this image,
!    i.e. from lcentr array, and calculates the mean value.
!    This value is then used to pass to the cleavage routine.
!  NOTE
!    If size( lcentr ) .eq. 0, then there are no FE associated
!    with coarray on this image. The set simg to 0.
!  SOURCE

integer, parameter :: comp=6 ! number of stress components

! Running total stress array
real( kind=rdef ) :: stot( comp )

! Temp stress array, to store remotely read values
real( kind=rdef ), allocatable :: str_tmp( : , : )

integer :: el, nel, rel, nintp, rimg, errstat

! Assertion check
! The number of stress components (6) is the same as dimension 3 of
! cgca_pfem_stress%stress
if ( size( cgca_pfem_stress%stress, dim=3 ) .ne. comp )                &
  error stop "ERROR: cgca_pfem_simg: &
   &size( cgca_pfem_stress%stress, dim=3 ) .ne. comp"

! Total number of elements linked to CA model on this image
! Total number of int. points per element
  nel = size( lcentr )
nintp = size( cgca_pfem_stress%stress, dim=2 )

! If there are no FE linked to this image, set simg to 0
! and return immediately.
if ( nel .eq. 0 ) then
  simg = 0.0_rdef
  return
end if

! Allocate tmp stress array
allocate( str_tmp( nintp, comp ), source=0.0_rdef, stat=errstat )
if ( errstat .ne. 0 )                                                  &
  error stop "ERROR: cgca_pfem_simg: allocate( str_tmp )"

! Add all stress tensors together. Loop over all elements linked
! to CA on this image and over all int. points.
stot = 0.0_rdef
do el=1, nel

  ! Calculate the image and the element numbers to read the stress
  ! data from.
  rimg = lcentr(el) % image
   rel = lcentr(el) % elnum

  ! Remote read of all stress values for this element
  str_tmp =                                                            & 
    real( cgca_pfem_stress[ rimg ] % stress( rel, : , : ), kind=rdef )

  ! Sum over all int. points, i.e. 1st dimension
  stot = stot + sum( str_tmp( : , : ), dim=1 )

end do

! Construct a (3,3) matrix from (6) vector.
! Observe the component order of ParaFEM
!   sx=stress(1)
!   sy=stress(2)
!   sz=stress(3)
!   txy=stress(4)
!   tyz=stress(5)
!   tzx=stress(6)
!   sigm=(sx+sy+sz)/three
!https://code.google.com/p/parafem/source/browse/trunk/parafem/src/modules/shared/new_library.f90
simg(1,1) = stot(1)
simg(2,2) = stot(2)
simg(3,3) = stot(3)
simg(1,2) = stot(4)
simg(2,3) = stot(5)
simg(3,1) = stot(6)
simg(2,1) = simg(1,2)
simg(3,2) = simg(2,3)
simg(1,3) = simg(3,1)

! calculate the mean
simg = simg / real( nel*nintp, kind=rdef )

! deallocate temp stress array
deallocate( str_tmp, stat=errstat )
if ( errstat .ne. 0 ) then
  write (*,'(a,i0)')                                                   &
    "ERROR: cgca_pfem_simg: deallocate( str_tmp ), err. code: ", errstat
  error stop
end if

end subroutine cgca_pfem_simg

!*roboend*


!*robodoc*s* cgca_m3pfem/cgca_pfem_intcalc1
!  NAME
!    cgca_pfem_intcalc1
!  SYNOPSIS

subroutine cgca_pfem_intcalc1( arrsize, fracvol )

!  INPUTS
!     arrsize - contains the 3 sizes of the space coarray.
!    Using the coarray sizes, the characteristic coarray
!    area is calculated.
!     fracvol - *real*, the number of failed (fractured) cells for each
!    image. It is calculated by cgca_fv, which *must* be called prior to
!    calling this routine.

integer( kind=iarr ), intent( in ) :: arrsize(3)
real( kind=rdef), intent( in ) :: fracvol

!  SIDE EFFECTS
!    cgca_pfem_integrity array changes
!  DESCRIPTION
!    All FEs linked to this image get the same value of integrity.
!    These are all FEs in lcentr array. For entry i in this
!    array this is FE cgca_pfem_integrity( lcentr(i)%elnum )
!    on image lcentr(i)%image.
!
!    The integrity is 1 minus the ratio of number of
!    failed cells to the cracteristic coarray area.
!    If integrity < 0, set it to zero.
!  USES
!    lcentr via host association
!  USED BY
!    end user?
!  SOURCE

real :: carea ! characteristic area
integer, parameter :: kind_integ = kind( cgca_pfem_integrity % i )
real( kind=kind_integ ), parameter :: one = 1_kind_integ
integer( kind=idef ) :: i

! Volume, in cells, is the product of 3 coarray sizes.
! Don't forget to remove the halo cells! 
! Characteristic area is volume ** 2/3
carea = product( real( arrsize-2 ) ) ** 0.66666666666666666667

do i = 1, size( lcentr )

 ! integrity is calculate as: i = 1 - min(1,f),
 ! Integrity has the range [1..0], where i=1 for f=0, i=0 for f=1.  
 ! f=fracvol / carea   - Fraction of failed cells, 0 if no fracture,
 !                       1 or above when I consider the CA to have no
 !                       load bearing capacity. 
 ! min( 1, fraction) - To make sure fraction is [0..1].
 cgca_pfem_integrity[ lcentr(i)%image ] % i( lcentr(i)%elnum ) =       &
          one - min( one, fracvol / carea )
end do

end subroutine cgca_pfem_intcalc1

!*roboend*


!*robodoc*s* cgca_m3pfem/cgca_pfem_cellin
!  NAME
!    cgca_pfem_cellin
!  SYNOPSIS

subroutine cgca_pfem_cellin( index, lc, lres, bcol, charlen, debug, flag )
!  INPUTS
!     index - represents each corner cell 
!     lc(3) - integer, local coordinates of a cell in the space
!    coarray on *this* image.
!     lres - real, linear resolution of the model, how many cells per
!    linear physical unit of length. 
!     bcol(3) - real, coordinates of the coarray box on this image
!    in physical units, in CA coord. system.
!     charlen - real, characteristic length of an FE in the model.
!    This parameter is used to determine whether a cell is "close
!    enough" to a centroid of an FE.
!     debug - logical. If .true. will dump some debug info

integer( kind=idef ), intent( in ) :: index
integer( kind=idef ), intent( in ) :: lc(index,3)
real( kind=rdef ), intent( in ) :: lres, bcol(3), charlen
logical( kind=ldef ), intent( in ) :: debug

!  OUTPUTS
!     flag - logical, .true. if the cell in "inside" the FE model,
!    .false. otherwise

logical( kind=ldef ), intent( out ) :: flag

!  SIDE EFFECTS
!    if debug is .true. dumps some output to stdout. Othewise none.
!  DESCRIPTION
!    Scan all FE in lcentr array. If the coordinates of the cell in FE
!    coord. system are close enough to the centroid of at least one
!    element, then the cell is "inside" the FE model. Otherwise it
!    is outside.
!  SOURCE

real( kind=rdef ) :: cacoord(3), cl2
integer( kind=idef ) :: i

! Assume the cell is outside by default. Only if it is proven to be in,
! change the flag to .true.
flag = .false.

! I need characteristic length squared for comparison
cl2 = charlen * charlen

! Calculate the global CA coord. of the given cell from the
! input local coord.
! lc / lres - distance in phys. units from the box lower corner 
! lc / lres + bcol - coord. in phys units in CA CS
cacoord = lc(index,:) / lres + bcol
!cacoord = lc / lres + bcol

if ( debug ) write (*,"(a,i0,a,3(es9.2,a))")                           &
  "DEBUG: cgca_pfem_cellin: img: ", this_image(),                      &
  " cacoord (",  cacoord(1), ",", cacoord(2), ",", cacoord(3), ")"

! Loop over all elements in lcentr
elements: do i = 1, size( lcentr )
  
  ! If the square of the distance between the cell and the centroid
  ! is less than the square of the characteristic length, then it's in
  if ( sum( (cacoord - lcentr(i)%centr(:) )**2 ) .lt. cl2 ) then
    flag = .true.
    exit elements
  end if
end do elements

end subroutine cgca_pfem_cellin

!*roboend*


!*robodoc*s* cgca_m3pfem/cgca_pfem_boxin
!  NAME
!    cgca_pfem_boxin
!  SYNOPSIS

subroutine cgca_pfem_boxin( lowr, uppr, lres, bcol, charlen, debug,    &
  iflag )

!  INPUTS
!      lowr(3) - integer, local coordinates of the lower corner cell
!                in the space coarray on *this* image.
!      uppr(3) - integer, local coordinates of the upper corner cell
!                in the space coarray on *this* image.
!         lres - real, linear resolution of the model, how many cells
!                per linear physical unit of length. 
!      bcol(3) - real, coordinates of the coarray box on this image
!                in physical units, in CA coord. system.
!      charlen - real, characteristic length of an FE in the model.
!                This parameter is used to determine whether a cell
!                is "close enough" to a centroid of an FE.
!        debug - logical. If .true. will dump some dubug output

integer( kind=idef ), intent( in ) :: lowr(3), uppr(3)
real( kind=rdef ), intent( in ) :: lres, bcol(3), charlen
logical( kind=ldef ), intent( in ) :: debug

!  OUTPUT
!    iflag - integer, 1 if all 8 corner cells of the box are inside
!            FE, 2 if all 8 corner cells of the box are outside of FE,
!            3 otherwise.

integer( kind=idef ), intent( out ) :: iflag

!  DESCRIPTION
!    This routine calculates whether 8 corner cells of the given box
!    are inside FE or not. There are 3 possibilities:
!    (1) all 8 cells are inside FE. In this case iflag is set to 1.
!    (2) all 8 cells are outside FE. In this case iflag is set to 2.
!    (3) some cells are inside and others are outside.
!    In this case iflag is set to 3. iflag will be used by a calling
!    routine to decide what to do next.
!
!    The cells in a box are numbered according the Fortran convention,
!    the leftmost index changes first:
!
!                  1                  3
!                  x------------------x   --> 2
!                 /.                 /|
!                / .                / |
!               /  .               /  |
!              /   .              /   |
!          5  /    .           7 /    |
!            x------------------x     | 
!            |     .            |     |
!          / |     x . . . . . .|. . .x 4
!         /  |    . 2           |    /
!        /   |   . .            |   /
!       3    |  .  .            |  /
!            | .   .            | /
!            |.    .            |/
!            x------------------x 
!          6       |              8
!                  |
!                  v
!                  1
!
!    This cell numbering convention is used for very small boxes,
!    i.e. when the box size is 1 along some direction.
!    So if the box size is 1 along 1, then following cells coincide:
!    2 as 1, 4 as 3, 6 as 5 and 8 as 7.
!    If the box size is 1 along 2, then following cells coincide:
!    3 as 1, 4 as 2, 7 as 5 and 8 as 6.
!    If the box size is 1 along 3, then following cells coincide:
!    5 as 1, 6 as 2, 7 as 3 and 8 as 4.
!    This logic is used in the code. Logical array same(3) is used
!    to show which box size is 1. same(i) is .true. if the box
!    is of size 1 cell along dimension i. The box size is 1 when
!    the lower cood. matches the upper coord. So same is calculated
!    simly as: same = lowr .eq. uppr.
!  SIDE EFFECTS
!    If debug is .true. dumps some output to stdout. Othewise none.
!  USES
!    cgca_pfem_cellin
!  SOURCE

integer( kind=idef ) :: i, local(8,3), img
logical( kind=ldef ) :: same(3), flagarr(8)

img = this_image()

! From given corner cells, calculate 8 corner local coordinates
local( 1, : ) = (/ lowr(1), lowr(2), lowr(3) /)
local( 2, : ) = (/ uppr(1), lowr(2), lowr(3) /)
local( 3, : ) = (/ lowr(1), uppr(2), lowr(3) /)
local( 4, : ) = (/ uppr(1), uppr(2), lowr(3) /)
local( 5, : ) = (/ lowr(1), lowr(2), uppr(3) /)
local( 6, : ) = (/ uppr(1), lowr(2), uppr(3) /)
local( 7, : ) = (/ lowr(1), uppr(2), uppr(3) /)
local( 8, : ) = (/ uppr(1), uppr(2), uppr(3) /)

! Take care of repeated cells, i.e. when one or more box dimensions
! is 1. same(i) is .true. if the box lower and upper coord. are
! the same along dimension i.
same = lowr .eq. uppr

! Call cgca_pfem_cellin for each corner cell.
main: do i = 1, 8

  ! Cell 1 will always be evaluated, i.e. will call cgca_pfem_cellin.
  ! All other cells will be evaluated only if they are unique.

  ! Check direction 1
  dir1: if ( same(1) ) then
    ! Box is single cell long along 1. This means the following
    ! cells have the same flag: 2 as 1, 4 as 3, 6 as 5, 8 as 7.
    if ( ( i .eq. 2 ) .or. ( i .eq. 4 ) .or.                           &
         ( i .eq. 6 ) .or. ( i .eq. 8 ) ) then
      flagarr(i) = flagarr(i-1)
      cycle main
    end if 
  end if dir1
  
  ! Check direction 2
  dir2: if ( same(2) ) then
    ! Box is single cell long along 2. This means the following
    ! cells have the same flag: 3 as 1, 4 as 2, 7 as 5, 8 as 6.
    if ( ( i .eq. 3 ) .or. ( i .eq. 4 ) .or.                           &
         ( i .eq. 7 ) .or. ( i .eq. 8 ) ) then
      flagarr(i) = flagarr(i-2)
      cycle main
    end if 
  end if dir2

  ! Check direction 3
  dir3: if ( same(3) ) then
    ! Box is single cell long along 3. This means the following
    ! cells have the same flag: 5 as 1, 6 as 2, 7 as 3, 8 as 4.
    if ( ( i .eq. 5 ) .or. ( i .eq. 6 ) .or.                           &
         ( i .eq. 7 ) .or. ( i .eq. 8 ) ) then
      flagarr(i) = flagarr(i-4)
      cycle main
    end if 
  end if dir3

  ! subroutine cgca_pfem_cellin( lc, lres, bcol, charlen, debug, flag )
  ! flag .eq. .true. if inside FE
  ! flag .eq. .false. if outside FE
  call cgca_pfem_cellin( i, local, lres, bcol, charlen, debug,    &
    flagarr(i) )

  if ( debug ) write (*,"(4(a,i0),a,l1)")                              &
    "DEBUG: cgca_pfem_boxin: img: ", img, " local cell coord (",       &
     local( i, 1 ), ",", local( i, 2 ), ",", local( i, 3 ),            &
     ") flag: ", flagarr(i)

end do main

! If all flags are .true. then the box is inside, set iflag to 1
! If all flags are .false. then the box is outside, set iflag to 2
! Otherwise, part of the box is in and part is out, set iflag to 3
if ( all( flagarr ) ) then
  iflag = 1
elseif ( all( .not. flagarr ) ) then
  iflag = 2
else
  iflag = 3
end if

end subroutine cgca_pfem_boxin

!*roboend*


!*robodoc*s* cgca_m3pfem/cgca_pfem_wholein
!  NAME
!    cgca_pfem_wholein
!  SYNOPSIS

subroutine cgca_pfem_wholein( coarray )

!  INPUT
!    coarray - main model coarray

integer( kind=iarr ), allocatable, intent( inout ) ::                  &
  coarray( : , : , : , : ) [ : , : , : ]

!  SIDE EFFECTS
!    state of coarray changed
!  DESCRIPTION
!    This is a very primitive routine to decide if cells are "inside"
!    the FE model or not. It work on the whole coarray on an image.
!    If there are no FE linked to coarray on this image, i.e. if
!    lcentr array is of zero length, then all cells
!    in the fracture layer of the coarray on this image are turned to
!    cgca_state_null.
!  SOURCE     

if ( size( lcentr ) .eq. 0 )                                           &
  coarray( : , : , : , cgca_state_type_frac ) = cgca_state_null
end subroutine cgca_pfem_wholein

!*roboend*


!*robodoc*s* cgca_m3pfem/cgca_pfem_partin
!  NAME
!    cgca_pfem_partin
!  SYNOPSIS

subroutine cgca_pfem_partin( coarray, cadim, lres, bcol, charlen, debug)

!  INPUT
!      coarray - main model coarray
!        cadim - coarray dimensions. Can calculate from coarray, but
!                these will be known already anyway, so makes sense
!                to pass as input
!         lres - real, linear resolution of the model, how many cells
!                per linear physical unit of length. 
!      bcol(3) - real, coordinates of the coarray box on this image
!                in physical units, in CA coord. system.
!      charlen - real, characteristic length of an FE in the model.
!                This parameter is used to determine whether a cell
!                is "close enough" to a centroid of an FE.
!        debug - logical. If .true. will dump some dubug output

integer( kind=iarr ), allocatable, intent( inout ) ::                  &
  coarray( : , : , : , : ) [ : , : , : ]
integer( kind=iarr ), intent( in ) :: cadim(3)
real( kind=rdef ), intent( in ) :: lres, bcol(3), charlen
logical( kind=ldef ), intent( in ) :: debug

!  SIDE EFFECTS
!    state of coarray changed
!  DESCRIPTION
!    This is the most thorough routine to decide if cells are "inside"
!    the FE model or not. It starts by checking boxes the size of the
!    whole coarray on this image. If a box in partially in and
!    partially out, it is split into two smaller boxes and the process
!    continues until each box is either fully in, or fully out.
!    If it is fully out, all fracture level cells of this box are
!    set to state cgca_state_null.
!  NOTES
!    This routine calls cgca_pfem_boxin, which in turn calls
!    cgca_pfem_cellin, which accesses lcentr array on its own image.
!    This routine updates coarray *locally*, on its own image only.
!    So there are no remote read or write operations in this routine.
!    No sync is needed inside this routine. Sync is likely needed
!    before and after. In particular, the coarray and
!    lcentr array must be defined on all images
!    prior to calling this routine on any image.
!  USES
!    cgca_pfem_boxin, cgca_boxsplit, cgca_inithead, cgca_addhead,
!    cgca_rmhead, cgca_lstdmp
!  USED BY
!    end user?
!  SOURCE     

integer( kind=idef ) :: lwr(3), upr(3), iflag, stat, lwr1(3),          &
  upr1(3), lwr2(3), upr2(3)
type( cgca_lnklst_tpayld ) :: payload
type( cgca_lnklst_node ), pointer :: head

integer :: iter, img

img = this_image()

! Start with a box the size of the whole coarray on this image
! Note a conversion from iarr to idef
lwr = 1
upr = int( cadim, kind=idef )

! Initialise the linked list with this box as the head node
! Returns the pointer to the head node
payload%lwr = lwr
payload%upr = upr
call cgca_inithead( head, payload )

! Start iteration counter
iter = 1

! Check all nodes on the list
list: do

  ! debug output
  if ( debug ) then
    write (*,'(2(a,i0),a)') "DEBUG: cgca_pfem_partin: img: ", img,     &
                         " iter: ", iter, " linked list dump:"
    call cgca_lstdmp( head )
  end if

  ! Get the payload from the head node on the list
  payload = head%value
  lwr = payload%lwr
  upr = payload%upr 

  ! Initialise iflag before it's called for the first time
  iflag=-1

  ! Check if this box is in/out
  !subroutine cgca_pfem_boxin( lowr, uppr, lres, bcol, charlen, debug,
  !   iflag )
  call cgca_pfem_boxin( lwr, upr, lres, bcol, charlen, debug, iflag )

  ! Remove this box from the linked list in any case
  ! Check stat value later.
  call cgca_rmhead( head, stat )

  ! The whole box in, iflag=1
  chkiflag: if ( iflag .eq. 1 ) then

    ! Exit if the list has zero nodes
    if ( stat .eq. 1 ) exit list

  ! The whole box out, iflag=2
  else if ( iflag .eq. 2 ) then
 
    ! Mark all cells of this box, on *this* image, in the fracture
    ! layer as cgca_state_null.
    coarray( lwr(1):upr(1) , lwr(2):upr(2), lwr(3):upr(3) ,            &
             cgca_state_type_frac ) = cgca_state_null

    ! Exit if the list has zero nodes
    if ( stat .eq. 1 ) exit list

  ! Part in/part out
  else if ( iflag .eq. 3 ) then

    ! Split the box in two along the biggest dimension of the box
    call cgca_boxsplit( lwr, upr, lwr1, upr1, lwr2, upr2 )

    ! Add two new boxes on top of head.
    ! The first box.
    payload%lwr = lwr1
    payload%upr = upr1

    ! For the first node, if the head is not associated, i.e.
    ! stat is 1, then use cgca_inithead instead of cgca_addhead.
    ! The head will not be associated if the initial box, i.e.
    ! the whole CA on this image has to be split, which must always
    ! happen in the beginning of the process.
    if ( stat .eq. 0 ) then
      call cgca_addhead( head, payload )
    else
      call cgca_inithead( head, payload )
    end if

    ! The second box.
    payload%lwr = lwr2
    payload%upr = upr2
    call cgca_addhead( head, payload )

  end if chkiflag

  ! increase the iteration counter
  iter = iter + 1

end do list

end subroutine cgca_pfem_partin

!*roboend*

end module cgca_m3pfem
