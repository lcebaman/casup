!*robodoc*f* cgca_m3clvg/m3clvg_sm1
!  NAME
!    m3clvg_sm1
!  SYNOPSIS

!$Id: m3clvg_sm1.f90 491 2018-02-20 22:22:58Z mexas $

submodule ( cgca_m3clvg ) m3clvg_sm1

!  DESCRIPTION
!    Submodule of module cgca_m3clvg. It contains subroutines dealing
!    with updating the grain connectivity array.
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    cgca_gcupda, cgca_gcupdn
!  USES
!    All variables and parameters of module cgca_m3clvg by host
!    association
!  USED BY
!    The parent module cgca_m3clvg.
!  SOURCE

implicit none

contains

!*roboend*


!*robodoc*s* m3clvg_sm1/cgca_gcupda
!  NAME
!    cgca_gcupda
!  SYNOPSIS  

  module subroutine cgca_gcupda( periodicbc )
    logical( kind=ldef ), intent( in ) :: periodicbc

!  SIDE EFFECTS
!    State of GB array in module cgca_m2gb is updated.
!  DESCRIPTION
!    This routine reads gcupd from *all* images and
!    adds the pairs to the local GB array on this image.
!    In other words it implements an all-to-all communication pattern.
!    If you want to use just the nearest neighbouring images,
!    use cgca_gcupdn instead.
!    Synchronisation must be used before and after
!    calling this routine, to comply with the standard.
!  NOTES
!    logical var periodicbc that is passed as input is ignored.
!    This input var is provided only to make the interface this routine
!    identical to cgca_gcupdn. Identical interfaces allow creating an
!    abstract interface for both gcupd and gcupdn and thus pass
!    the name of the routine as input to cgca_clvgp* cleavage
!    propagation routines.
!  USES
!    cgca_gcf, cgca_ico
!  SOURCE

integer( kind=idef ) :: i, j, img, nimgs, img_curr, rndint,            &
  cosub( cgca_scodim ), flag
integer( kind=kind(gcupd) ) ::                                         &
  gcupd_local( cgca_gcupd_size1, cgca_gcupd_size2 )
real :: rnd

logical :: l_tmp

! Just to suppress the unused compiler warning
l_tmp = periodicbc

  img = this_image()
nimgs = num_images()

! choose the first image at random
call random_number( rnd )   ! [ 0 .. 1 )
rndint = int( rnd*nimgs )+1 ! [ 1 .. nimgs ]

! Initialise flag with any value
flag = 0

! loop over all images, starting at a randomly chosen image
images: do j = rndint, rndint+nimgs-1

  ! Get the current image number.
  ! If it's > nimgs, subtract nimgs
  img_curr = j
  if ( img_curr .gt. nimgs ) img_curr = img_curr - nimgs

  ! Skip this image, because the GC array has already been updated
  ! in cgca_clvgsd.
  if ( img_curr .eq. img ) cycle images

  ! Calculate cosubscripts:
  call cgca_ico( img_curr, cosub, flag )
  if ( flag .ne. 0 ) then
    write (*,*) &
      "ERROR: m3clvg/cgca_gcupda: cgca_ico exit with error: flag:", flag
    error stop
  end if

  ! copy gcupd from image j into a local var
  gcupd_local( : , : ) = gcupd( : , : ) [ cosub(1), cosub(2), cosub(3) ]

  gcarray: do i = 1, cgca_gcupd_size1

    ! The gcupd array is filled with fractured pairs from the beginning
    ! so exit as soon as the GB state is intact.
    if ( gcupd_local( i , 3 ) .eq. cgca_gb_state_intact ) exit gcarray

!write (*,*) "DEBUG: cgca_gcupd: img:", img, &
!            "gcupd_local(i,:):", gcupd_local( i , : )

    ! add the pair to the GC array on this image
    call cgca_gcf( gcupd_local( i , 1 ), gcupd_local( i , 2 ) )

 end do gcarray
end do images

end subroutine cgca_gcupda

!*roboend*


!*robodoc*s* m3clvg_sm1/cgca_gcupdn
!  NAME
!    cgca_gcupdn
!  SYNOPSIS  

  module subroutine cgca_gcupdn( periodicbc )
    logical( kind=ldef ), intent( in ) :: periodicbc

!  INPUT
!    periodicbc - logical, .true. if the CA space has periodic BC,
!    and .false. otherwise.
!  SIDE EFFECTS
!    State of GB array in module cgca_m2gb is updated
!  DESCRIPTION
!    This routine reads gcupd from the *nearest neighbouring*
!    images only, and adds the pairs to the local GB array on this
!    image. If you want to read from all images, use
!    cgca_gcupd.
!    Synchronisation must be used before and after
!    calling this routine, to comply with the standard.
!  NOTES
!    This routine must be used only after gcupd has been
!    allocated. A runtime error will result if gcupd has not
!    been allocated yet. 
!  USES
!    cgca_gcf
!  SOURCE

integer( kind=idef ) :: i, j, k, s, mycod( cgca_scodim ),              &
 neicod( cgca_scodim )
integer( kind=kind(gcupd) ) ::                                     &
 gcupd_local( cgca_gcupd_size1, cgca_gcupd_size2 )

! Get my coindex set
if ( .not. allocated( gcupd ) ) then
  write (*,'(a)')                                                      &
     "ERROR: cgca_m3clvg/cgca_gcupdn: gcupd not allocated"
  error stop
end if
mycod = this_image( gcupd )

! Loop over all nearest neighbours, taking special attention of
! the images at the edges of the model
do i = -1 , 1
do j = -1 , 1
inner: do k = -1 , 1

  ! Get the coindex set of the neighbour
  neicod = mycod + (/ i, j, k /)

  ! Skip this image
  if ( all( neicod .eq. mycod ) ) cycle inner

  ! Dealing with edges
  ! Loop over all codimensions
  edges: do s = 1 , cgca_scodim

    ! If the neighbour is below the lower edge
    if ( neicod( s ) .lt. cgca_slcob( s ) ) then
      if ( periodicbc ) then
        ! If periodic BC are in use, take the data from the opposite
        ! edge.
        neicod( s ) = cgca_sucob( s )       
      else
        ! Otherwise, do not pull data from this neighbour, move to
        ! the next one.
        cycle inner
      end if
    end if     

    ! If the neighbour is above the upper edge
    if ( neicod( s ) .gt. cgca_sucob( s ) ) then
      if ( periodicbc ) then
        ! If periodic BC are in use, take the data from the opposite
        ! edge.
        neicod( s ) = cgca_slcob( s )       
      else
        ! Otherwise, do not pull data from this neighbour, move to
        ! the next one.
        cycle inner
      end if
    end if     

  end do edges

  ! Now the coindex set of the neighbour has been obtained.
  ! Pull its data

  ! Copy gcupd from the neighbouring image into a local var.
  ! Remote read.
  gcupd_local( : , : ) =                                               &
    gcupd( : , : ) [ neicod(1), neicod(2), neicod(3) ]

  ! Scan all values in gcupd. Can reuse loop index s.
  gcarray: do s = 1, cgca_gcupd_size1

    ! The gcupd array is filled with fractured pairs from the beginning
    ! so exit as soon as the GB state is intact.
    if ( gcupd_local( s , 3 ) .eq. cgca_gb_state_intact ) exit gcarray

    !write (*,*) "DEBUG: cgca_gcupd: img:", img, &
    !            "gcupd_local(i,:):", gcupd_local( i , : )

    ! add the pair to the GC array on this image
    call cgca_gcf( gcupd_local( s , 1 ), gcupd_local( s , 2 ) )

  end do gcarray

end do inner
end do
end do

end subroutine cgca_gcupdn

!*roboend*

end submodule m3clvg_sm1
