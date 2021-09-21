!*robodoc*f* ca_hx/ca_hx_mpi
!  NAME
!    ca_hx_mpi
!  SYNOPSIS

!$Id: ca_hx_mpi.f90 528 2018-03-26 09:02:14Z mexas $

submodule ( ca_hx ) ca_hx_mpi

!  DESCRIPTION
!    Submodule of module ca_hx with MPI related routines.
!    To aid portability, the module works only with default integer
!    kind, i.e. MPI_integer. Other MPI integer kinds might not be
!    widely available, meaning that other Fortran integer kinds might
!    be less portable. So make sure that space array kind is the same
!    as default integer. This is likely to be the case with 
!      integer, parameter :: iarr = selected_int_kind( 8 )
!    iarr is set in cgca_m1co.
!
!    Creation/release (free) of MPI types is left as the user's
!    responsibility. This is because the user might want to change
!    halo depth in the same program. This is hard/impossible to keep
!    completely invisible to the user.
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    ca_mpi_halo_type_create, ca_mpi_halo_type_free, ca_mpi_hxvn1m, &
!    ca_mpi_hxvn1p, ca_mpi_hxvn2m, ca_mpi_hxvn2p, ca_mpi_hxvn3m, &
!    ca_mpi_hxvn3p, ca_mpi_hx_all
!  USES
!  USED BY
!  SOURCE
!
! For reference
!
! MPI_SEND(  BUF, COUNT, DATATYPE, DEST,   TAG, COMM,          IERROR )
! MPI_RECV(  BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS,  IERROR )
!
! MPI_ISEND( BUF, COUNT, DATATYPE, DEST,   TAG, COMM, REQUEST, IERROR )
! MPI_IRECV( BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, REQUEST, IERROR )

implicit none

! Tags for sending messages in 6 directions.
integer, parameter :: TAG1L = 1, TAG1R = 2, TAG2L = 3, TAG2R = 4,      &
  TAG3L = 5, TAG3R = 6

integer, save :: &
  rank,          & ! MPI rank
  !status( MPI_STATUS_SIZE ), & ! used in MPI_RECV, etc.
  mpi_h1_LV,     & ! MPI halo, dim 1, left virtual
  mpi_h1_LR,     & ! MPI halo, dim 1, left real
  mpi_h1_RR,     & ! MPI halo, dim 1, right real
  mpi_h1_RV,     & ! MPI halo, dim 1, right virtual
  mpi_h2_LV,     & ! MPI halo, dim 2, left virtual
  mpi_h2_LR,     & ! MPI halo, dim 2, left real
  mpi_h2_RR,     & ! MPI halo, dim 2, right real
  mpi_h2_RV,     & ! MPI halo, dim 2, right virtual
  mpi_h3_LV,     & ! MPI halo, dim 3, left virtual
  mpi_h3_LR,     & ! MPI halo, dim 3, left real
  mpi_h3_RR,     & ! MPI halo, dim 3, right real
  mpi_h3_RV,     & ! MPI halo, dim 3, right virtual
  mpi_ca_integer,& ! MPI matching type for iarr
  errcode,       & ! Need to preserve ierr
  errlen           ! The length of the output error message

! A flag to track the state of MPI types for halos.
! Set initially to .false.
! Calling ca_mpi_halo_type_create sets it to .true.
! Calling ca_mpi_halo_type_free sets it to .false. again.

logical, save :: halo_type_created = .false.

contains

!*roboend*


!*robodoc*s* ca_hx/ca_mpi_halo_type_create
!  NAME
!    ca_mpi_halo_type_create
!  SYNOPSIS

module subroutine ca_mpi_halo_type_create( space )

!  INPUT

integer( kind=iarr), intent( inout ), allocatable :: space(:,:,:)

!    space - the CA array
!  OUTPUT
!    none
!  SIDE EFFECTS
!    12 MPI halo types, module variables, are created.
!  DESCRIPTION
!    For each direction there are 4 MPI halo data types:
!    - array elements in the halo part of the array to the
!      left of the real data,
!    - array elements of halo thickness inside the real part of the
!      array on its left side,
!    - array elements of halo thickness inside the real part of the
!      array on its right side,
!    - array elements in the hallp part of the array to the right
!      of the real data.
!    Refer to the diagram in ca_hx/ca_spalloc.
!  NOTES
!    Call this routine after ca_halloc.
!    All images must call this routine!
!    Pay particular attention to the starts of all arrays.
!    Refer to the details in e.g:
!      https://www.open-mpi.org/doc/v3.0/man3/MPI_Type_create_subarray.3.php
!    In particular:
!      In a Fortran program with arrays indexed starting from 1,
!      if the starting coordinate of a particular dimension
!      of the subarray is n, then the entry in array of starts
!      for that dimension is n-1. 
!    A diagram is probably needed for starts, because it's different
!    from that in ca_hx/ca_spalloc.
!    Using only a single dimension, e.g. 1.
!
!      +------+------+----------------+------+------+
!      |  LV  |  LR  |                |  RR  |  RV  |
!      +------+------+----------------+------+------+
!      ^      ^                       ^      ^
!      |      |                       |      |
!      0      hdepth                  |      sizes(1)-hdepth
!                                     sizes(1)-2*hdepth
!
!      starts for 4 halo arrays along dim 1
! 
!  USES
!  USED BY
!  SOURCE

!MPI_TYPE_CREATE_SUBARRAY(NDIMS, ARRAY_OF_SIZES, ARRAY_OF_SUBSIZES,
!    ARRAY_OF_STARTS, ORDER, OLDTYPE, NEWTYPE, IERROR)
!
!    INTEGER    NDIMS, ARRAY_OF_SIZES(*), ARRAY_OF_SUBSIZES(*),
!    ARRAY_OF_STARTS(*), ORDER, OLDTYPE, NEWTYPE, IERROR

integer :: sizes(3), subsizes(3), starts(3)

! Set MPI rank, keep forever
call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr)
!write (*,*) "my rank:", rank, "my img:", this_image()

! Set MPI matching type for iarr: mpi_ca_integer.
! Set once, keep forever.
call MPI_TYPE_CREATE_F90_INTEGER( ca_range, mpi_ca_integer, ierr )

! The sizes is just the shape of the space array, for all cases
sizes = shape( space )

! Dimension 1

subsizes = (/ hdepth, sub(2), sub(3) /)

! 1. dimension 1, left virtual (LV) type

starts = (/ 0, hdepth, hdepth /)

!write (*,"(3(a,3(i0,tr1)))") "sizes: ", sizes, " subsizes: ", subsizes, " starts: ", starts

call MPI_TYPE_CREATE_SUBARRAY( 3, sizes, subsizes, starts,             &
  MPI_ORDER_FORTRAN, mpi_ca_integer, mpi_h1_LV, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  errcode = ierr
  call MPI_ERROR_STRING( errcode, errmsg, errlen, ierr )
  write (*,"(a,i0,a)") "ERROR ca_hx_mpi/ca_mpi_halo_type_create: " //  &
    "MPI_TYPE_CREATE_SUBARRAY: dim 1: left virtual (LV): error: ",     &
    errcode, " error message: " // trim(errmsg)
  error stop
end if

call MPI_TYPE_COMMIT( mpi_h1_LV, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type_create: " //    &
    "MPI_TYPE_COMMIT: dim 1: left virtual (LV): ierr: ", ierr
  error stop
end if

! 2. dimension 1, left real (LR) type

starts = (/ hdepth, hdepth, hdepth /)

call MPI_TYPE_CREATE_SUBARRAY( 3, sizes, subsizes, starts,             &
  MPI_ORDER_FORTRAN, mpi_ca_integer, mpi_h1_LR, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type_create: " //    &
    "MPI_TYPE_CREATE_SUBARRAY: dim 1: left real (LR): ierr: ", ierr
  error stop
end if

call MPI_TYPE_COMMIT( mpi_h1_LR, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type_create: " //    &
    "MPI_TYPE_COMMIT: dim 1: left real (LR): ierr: ", ierr
  error stop
end if

! 3. dimension 1, right real (RR) type

starts = (/ sizes(1) - 2*hdepth, hdepth, hdepth /)

call MPI_TYPE_CREATE_SUBARRAY( 3, sizes, subsizes, starts,             &
  MPI_ORDER_FORTRAN, mpi_ca_integer, mpi_h1_RR, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type_create: " //    &
    "MPI_TYPE_CREATE_SUBARRAY: dim 1: right real (RR): ierr: ", ierr
  error stop
end if

call MPI_TYPE_COMMIT( mpi_h1_RR, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type_create: " //    &
    "MPI_TYPE_COMMIT: dim 1: right real (RR): ierr: ", ierr
  error stop
end if

! 4. dimension 1, right virtual (RV) type

starts = (/ sizes(1) - hdepth, hdepth, hdepth /)

call MPI_TYPE_CREATE_SUBARRAY( 3, sizes, subsizes, starts,             &
  MPI_ORDER_FORTRAN, mpi_ca_integer, mpi_h1_RV, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type_create: " //    &
    "MPI_TYPE_CREATE_SUBARRAY: dim 1: right virtual (RV): ierr: ", ierr
  error stop
end if

call MPI_TYPE_COMMIT( mpi_h1_RV, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type_create: " //    &
    "MPI_TYPE_COMMIT: dim 1: right virtual (RV): ierr: ", ierr
  error stop
end if

! Dimension 2

subsizes = (/ sub(1), hdepth, sub(3) /)

! 5. dimension 2, left virtual (LV) type

starts = (/ hdepth, 0, hdepth /)

call MPI_TYPE_CREATE_SUBARRAY( 3, sizes, subsizes, starts,             &
  MPI_ORDER_FORTRAN, mpi_ca_integer, mpi_h2_LV, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type: " //           &
    "MPI_TYPE_CREATE_SUBARRAY: dim 2: left virtual (LV): ierr: ", ierr
  error stop
end if

call MPI_TYPE_COMMIT( mpi_h2_LV, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type: " //           &
    "MPI_TYPE_COMMIT: dim 2: left virtual (LV): ierr: ", ierr
  error stop
end if

! 6. dimension 2, left real (LR) type

starts = (/ hdepth, hdepth, hdepth /)

call MPI_TYPE_CREATE_SUBARRAY( 3, sizes, subsizes, starts,             &
  MPI_ORDER_FORTRAN, mpi_ca_integer, mpi_h2_LR, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type: " //           &
    "MPI_TYPE_CREATE_SUBARRAY: dim 2: left real (LR): ierr: ", ierr
  error stop
end if

call MPI_TYPE_COMMIT( mpi_h2_LR, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type: " //           &
    "MPI_TYPE_COMMIT: dim 2: left real (LR): ierr: ", ierr
  error stop
end if

! 7. dimension 2, right real (RR) type

starts = (/ hdepth, sizes(2) - 2*hdepth, hdepth /)

call MPI_TYPE_CREATE_SUBARRAY( 3, sizes, subsizes, starts,             &
  MPI_ORDER_FORTRAN, mpi_ca_integer, mpi_h2_RR, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type: " //           &
    "MPI_TYPE_CREATE_SUBARRAY: dim 2: right real (RR): ierr: ", ierr
  error stop
end if

call MPI_TYPE_COMMIT( mpi_h2_RR, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type: " //           &
    "MPI_TYPE_COMMIT: dim 2: right real (RR): ierr: ", ierr
  error stop
end if

! 8. dimension 2, right virtual (RV) type

starts = (/ hdepth, sizes(2) - hdepth, hdepth /)

call MPI_TYPE_CREATE_SUBARRAY( 3, sizes, subsizes, starts,             &
  MPI_ORDER_FORTRAN, mpi_ca_integer, mpi_h2_RV, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type: " //           &
    "MPI_TYPE_CREATE_SUBARRAY: dim 2: right virtual (RV): ierr: ", ierr
  error stop
end if

call MPI_TYPE_COMMIT( mpi_h2_RV, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type: " //           &
    "MPI_TYPE_COMMIT: dim 2: right virtual (RV): ierr: ", ierr
  error stop
end if

! Dimension 3

subsizes = (/ sub(1), sub(2), hdepth /)

! 9. dimension 3, left virtual (LV) type

starts = (/ hdepth, hdepth, 0 /)

call MPI_TYPE_CREATE_SUBARRAY( 3, sizes, subsizes, starts,             &
  MPI_ORDER_FORTRAN, mpi_ca_integer, mpi_h3_LV, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type: " //           &
    "MPI_TYPE_CREATE_SUBARRAY: dim 3: left virtual (LV): ierr: ", ierr
  error stop
end if

call MPI_TYPE_COMMIT( mpi_h3_LV, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type: " //           &
    "MPI_TYPE_COMMIT: dim 3: left virtual (LV): ierr: ", ierr
  error stop
end if

! 10. dimension 3, left real (LR) type

starts = (/ hdepth, hdepth, hdepth /)

call MPI_TYPE_CREATE_SUBARRAY( 3, sizes, subsizes, starts,             &
  MPI_ORDER_FORTRAN, mpi_ca_integer, mpi_h3_LR, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type: " //           &
    "MPI_TYPE_CREATE_SUBARRAY: dim 3: left real (LR): ierr: ", ierr
  error stop
end if

call MPI_TYPE_COMMIT( mpi_h3_LR, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type: " //           &
    "MPI_TYPE_COMMIT: dim 3: left real (LR): ierr: ", ierr
  error stop
end if

! 11. dimension 3, right real (RR) type

starts = (/ hdepth, hdepth, sizes(3) - 2*hdepth /)

call MPI_TYPE_CREATE_SUBARRAY( 3, sizes, subsizes, starts,             &
  MPI_ORDER_FORTRAN, mpi_ca_integer, mpi_h3_RR, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type: " //           &
    "MPI_TYPE_CREATE_SUBARRAY: dim 3: right real (RR): ierr: ", ierr
  error stop
end if

call MPI_TYPE_COMMIT( mpi_h3_RR, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type: " //           &
    "MPI_TYPE_COMMIT: dim 3: right real (RR): ierr: ", ierr
  error stop
end if

! 12. dimension 3, right virtual (RV) type

starts = (/ hdepth, hdepth, sizes(3) - hdepth /)

call MPI_TYPE_CREATE_SUBARRAY( 3, sizes, subsizes, starts,             &
  MPI_ORDER_FORTRAN, mpi_ca_integer, mpi_h3_RV, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type: " //           &
    "MPI_TYPE_CREATE_SUBARRAY: dim 3: right virtual (RV): ierr: ", ierr
  error stop
end if

call MPI_TYPE_COMMIT( mpi_h3_RV, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type: " //           &
    "MPI_TYPE_COMMIT: dim 3: right virtual (RV): ierr: ", ierr
  error stop
end if

! MPI types for halos have been created.
! Set the corresponding flag to .true.
halo_type_created = .true.

end subroutine ca_mpi_halo_type_create

!*roboend*


!*robodoc*s* ca_hx/ca_mpi_halo_type_free
!  NAME
!    ca_mpi_halo_type_free
!  SYNOPSIS

module subroutine ca_mpi_halo_type_free

!  INPUT
!    none
!  OUTPUT
!    none
!  SIDE EFFECTS
!    12 MPI halo types, module variables, are freed.
!  DESCRIPTION
!    Refer to ca_mpi_halo_type_create for details of these 12 types.
!    Need to call this routine if want to re-create the halo types,
!    perhaps with different halo depth, or for a different space
!    array.
!  NOTES
!    Will give an error if data types are not committed.
!    All images must call this routine!
!  USES
!  USED BY
!  SOURCE

! Dimension 1

call MPI_TYPE_FREE( mpi_h1_LV, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type_free: " //      &
    "MPI_TYPE_FREE: dim 1: left virtual (LV): ierr: ", ierr
  error stop
end if

call MPI_TYPE_FREE( mpi_h1_LR, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type_free: " //      &
    "MPI_TYPE_FREE: dim 1: left real (LR): ierr: ", ierr
  error stop
end if

call MPI_TYPE_FREE( mpi_h1_RR, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type_free: " //      &
    "MPI_TYPE_FREE: dim 1: right real (RR): ierr: ", ierr
  error stop
end if

call MPI_TYPE_FREE( mpi_h1_RV, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type_free: " //      &
    "MPI_TYPE_FREE: dim 1: right virtual (RV): ierr: ", ierr
  error stop
end if

! Dimension 2

call MPI_TYPE_FREE( mpi_h2_LV, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type_free: " //      &
    "MPI_TYPE_FREE: dim 2: left virtual (LV): ierr: ", ierr
  error stop
end if

call MPI_TYPE_FREE( mpi_h2_LR, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type_free: " //      &
    "MPI_TYPE_FREE: dim 2: left real (LR): ierr: ", ierr
  error stop
end if

call MPI_TYPE_FREE( mpi_h2_RR, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type_free: " //      &
    "MPI_TYPE_FREE: dim 2: right real (RR): ierr: ", ierr
  error stop
end if

call MPI_TYPE_FREE( mpi_h2_RV, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type_free: " //      &
    "MPI_TYPE_FREE: dim 2: right virtual (RV): ierr: ", ierr
  error stop
end if

! Dimension 3

call MPI_TYPE_FREE( mpi_h3_LV, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type_free: " //      &
    "MPI_TYPE_FREE: dim 3: left virtual (LV): ierr: ", ierr
  error stop
end if

call MPI_TYPE_FREE( mpi_h3_LR, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type_free: " //      &
    "MPI_TYPE_FREE: dim 3: left real (LR): ierr: ", ierr
  error stop
end if

call MPI_TYPE_FREE( mpi_h3_RR, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type_free: " //      &
    "MPI_TYPE_FREE: dim 3: right real (RR): ierr: ", ierr
  error stop
end if

call MPI_TYPE_FREE( mpi_h3_RV, ierr )
if ( ierr .ne. MPI_SUCCESS ) then
  write (*,"(a,i0)") "ERROR ca_hx_mpi/ca_mpi_halo_type_free: " //      &
    "MPI_TYPE_FREE: dim 3: right virtual (RV): ierr: ", ierr
  error stop
end if

! MPI types have been freed.
! Reset the flag back to .false.
! Will need to re-create MPI types for halos *before* any HX.

halo_type_created = .false.

end subroutine ca_mpi_halo_type_free

!*roboend*


!*robodoc*s* ca_hx/ca_mpi_hxvn1m
!  NAME
!    ca_mpi_hxvn1m
!  SYNOPSIS

module subroutine ca_mpi_hxvn1m( space )

!  INPUT

integer( kind=iarr), intent(inout), allocatable :: space(:,:,:)

!    space - the CA array
!  OUTPUT
!    space is updated
!  SIDE EFFECTS
!    none
!  DESCRIPTION
!    Use non-blocking send/receive.
!    An image does 2 remote ops:
!    - Send its space array real halo layer, left side (mpi_h1_LR)
!      along dimension 1 into a virtual halo layer, right side
!      (mpi_h1_RV) on an image which is 1 lower along codimension 1.
!      Tag this message with TAG1L.
!    - Receive its space array virtual halo layer, left side
!      (mpi_h1_LV) along dimension 1 from a real halo layer,
!      right side (mpi_h1_RR) on an image which is 1 lower along
!      codimension 1. Tag this message with TAG1R.
!
!    Schematic diagram, only showing what is relevant for HX along
!    dimension 1:
!
!        ----------> dimension 1
!  
!       image P / rank P+1       |       image Q / rank Q+1
!                                |
!                                |
!                                | image Q, TAG1L, send data type mpi_h1_LR
!                       +--------|------------+
!                       |        |            |
!                       |        |            |
!      +----------------|-+      |      +-----|------------+ 
!      |                | |      |      |     |            |
!      |   +----------+-V-+      |      +---+-|------------+
!      |   |          |   |      |      |   | ^            |
!      |   |          | h |      |      | h |              |
!      |   |   real   | a |      |      | a |   real       |
!      |   |          | l |      |      | l |              |
!      |   |          | o |      |      | o |              |
!      |   |         V|   |      |      |   |              |
!      |   +---------|+---+      |      +-^-+--------------+
!      |             |    |      |      | |                |
!      +-------------|----+      |      +-+----------------+ 
!                    |           |        |
!                    |           |        |
!                    +-----------|--------+
!                                | image Q, TAG1R, receive data type mpi_h1_LV
!                                |
!                                |
!
!  USES
!  USED BY
!    ca_mpi_hx_all
!  SOURCE

integer :: reqs1m(2), stats(MPI_STATUS_SIZE, 2)

if ( ci(1) .ne. 1 ) then

  ! Rank is image number -1.

  ! Receive from the left neighbour, tag = TAG1R
  call MPI_IRECV( space, 1, mpi_h1_LV, nei_img_L(1)-1, TAG1R,          &
    MPI_COMM_WORLD, reqs1m(1), ierr )

  ! Send to the left neighbour, tag = TAG1L
  call MPI_ISEND( space, 1, mpi_h1_LR, nei_img_L(1)-1, TAG1L,          &
    MPI_COMM_WORLD, reqs1m(2), ierr )

  call MPI_WAITALL( 2, reqs1m, stats, ierr )

end if

end subroutine ca_mpi_hxvn1m

!*roboend*


!*robodoc*s* ca_hx/ca_mpi_hxvn1p
!  NAME
!    ca_mpi_hxvn1p
!  SYNOPSIS

module subroutine ca_mpi_hxvn1p( space )

!  INPUT

integer( kind=iarr), intent(inout), allocatable :: space(:,:,:)

!    space - the CA array
!  OUTPUT
!    space is updated
!  SIDE EFFECTS
!    none
!  DESCRIPTION
!    Use non-blocking send/receive.
!    An image does 2 remote ops:
!    - Send its space array real halo layer, right side (mpi_h1_RR)
!      along dimension 1 into a virtual halo layer, left side
!      (mpi_h1_LV) on an image which is 1 higher along codimension 1.
!      Tag this message with TAG1R.
!    - Receive its space array virtual halo layer, right side
!      (mpi_h1_RV) along dimension 1 from a real halo layer,
!      right side (mpi_h1_LR) on an image which is 1 higher along
!      codimension 1. Tag this message with TAG1L.
!
!    Schematic diagram, only showing what is relevant for HX along
!    dimension 1:
!
!        ----------> dimension 1
!  
!       image P / rank P+1       |       image Q / rank Q+1
!                                |
!                                |
!       image P, TAG1L, receive data type mpi_h1_RV
!                       +--------|------------+
!                       |        |            |
!                       |        |            |
!      +----------------|-+      |      +-----|------------+ 
!      |                | |      |      |     |            |
!      |   +----------+-V-+      |      +---+-|------------+
!      |   |          |   |      |      |   | ^            |
!      |   |          | h |      |      | h |              |
!      |   |   real   | a |      |      | a |   real       |
!      |   |          | l |      |      | l |              |
!      |   |          | o |      |      | o |              |
!      |   |         V|   |      |      |   |              |
!      |   +---------|+---+      |      +-^-+--------------+
!      |             |    |      |      | |                |
!      +-------------|----+      |      +-+----------------+ 
!                    |           |        |
!                    |           |        |
!                    +-----------|--------+
!       image P, TAG1R, send data type mpi_h1_RR
!
!  USES
!  USED BY
!    ca_mpi_hx_all
!  SOURCE

integer :: reqs1p(2), stats(MPI_STATUS_SIZE, 2)

if ( ci(1) .ne. ucob(1) ) then

  ! Rank is image number -1.

  ! Receive from the right neighbour, tag = TAG1L
  call MPI_IRECV( space, 1, mpi_h1_RV, nei_img_R(1)-1, TAG1L,          &
    MPI_COMM_WORLD, reqs1p(1), ierr )

  ! Send to the right neighbour, tag = TAG1R
  call MPI_ISEND( space, 1, mpi_h1_RR, nei_img_R(1)-1, TAG1R,          &
    MPI_COMM_WORLD, reqs1p(2), ierr )

  call MPI_WAITALL( 2, reqs1p, stats, ierr )

end if

end subroutine ca_mpi_hxvn1p

!*roboend*


!*robodoc*s* ca_hx/ca_mpi_hxvn2m
!  NAME
!    ca_mpi_hxvn2m
!  SYNOPSIS

module subroutine ca_mpi_hxvn2m( space )

!  INPUT

integer( kind=iarr ), intent(inout), allocatable :: space(:,:,:)

!    space - the CA array
!  OUTPUT
!    space is updated
!  SIDE EFFECTS
!    none
!  DESCRIPTION
!    HX along dimension 2. See ca_mpi_hxvn1m.
!  USES
!  USED BY
!    ca_mpi_hx_all
!  SOURCE

integer :: reqs2m(2), stats(MPI_STATUS_SIZE, 2)

if ( ci(2) .ne. 1 ) then

  ! Rank is image number -1.

  ! Receive from the left neighbour, tag = TAG2R
  call MPI_IRECV( space, 1, mpi_h2_LV, nei_img_L(2)-1, TAG2R,          &
    MPI_COMM_WORLD, reqs2m(1), ierr )

  ! Send to the left neighbour, tag = TAG2L
  call MPI_ISEND( space, 1, mpi_h2_LR, nei_img_L(2)-1, TAG2L,          &
    MPI_COMM_WORLD, reqs2m(2), ierr )

  call MPI_WAITALL( 2, reqs2m, stats, ierr )

end if

end subroutine ca_mpi_hxvn2m

!*roboend*


!*robodoc*s* ca_hx/ca_mpi_hxvn2p
!  NAME
!    ca_mpi_hxvn2p
!  SYNOPSIS

module subroutine ca_mpi_hxvn2p( space )

!  INPUT

integer( kind=iarr ), intent(inout), allocatable :: space(:,:,:)

!    space - the CA array
!  OUTPUT
!    space is updated
!  SIDE EFFECTS
!    none
!  DESCRIPTION
!    HX along dimension 2. See ca_mpi_hxvn1p.
!  USES
!  USED BY
!    ca_mpi_hx_all
!  SOURCE

integer :: reqs2p(2), stats(MPI_STATUS_SIZE, 2)

if ( ci(2) .ne. ucob(2) ) then

  ! Rank is image number -1.

  ! Receive from the right neighbour, tag = TAG2L
  call MPI_IRECV( space, 1, mpi_h2_RV, nei_img_R(2)-1, TAG2L,          &
    MPI_COMM_WORLD, reqs2p(1), ierr )

  ! Send to the right neighbour, tag = TAG2R
  call MPI_ISEND( space, 1, mpi_h2_RR, nei_img_R(2)-1, TAG2R,          &
    MPI_COMM_WORLD, reqs2p(2), ierr )

  call MPI_WAITALL( 2, reqs2p, stats, ierr )

end if

end subroutine ca_mpi_hxvn2p

!*roboend*


!*robodoc*s* ca_hx/ca_mpi_hxvn3m
!  NAME
!    ca_mpi_hxvn3m
!  SYNOPSIS

module subroutine ca_mpi_hxvn3m( space )

!  INPUT

integer( kind=iarr ), intent(inout), allocatable :: space(:,:,:)

!    space - the CA array
!  OUTPUT
!    space is updated
!  SIDE EFFECTS
!    none
!  DESCRIPTION
!    HX along dimension 3. See ca_mpi_hxvn1m.
!  USES
!  USED BY
!    ca_mpi_hx_all
!  SOURCE

integer :: reqs3m(2), stats(MPI_STATUS_SIZE, 2)

if ( ci(3) .ne. 1 ) then

  ! Rank is image number -1.

  ! Receive from the left neighbour, tag = TAG3R
  call MPI_IRECV( space, 1, mpi_h3_LV, nei_img_L(3)-1, TAG3R,          &
    MPI_COMM_WORLD, reqs3m(1), ierr )

  ! Send to the left neighbour, tag = TAG3L
  call MPI_ISEND( space, 1, mpi_h3_LR, nei_img_L(3)-1, TAG3L,          &
    MPI_COMM_WORLD, reqs3m(2), ierr )

  call MPI_WAITALL( 2, reqs3m, stats, ierr )

end if

end subroutine ca_mpi_hxvn3m

!*roboend*


!*robodoc*s* ca_hx/ca_mpi_hxvn3p
!  NAME
!    ca_mpi_hxvn3p
!  SYNOPSIS

module subroutine ca_mpi_hxvn3p( space )

!  INPUT

integer( kind = iarr ), intent(inout), allocatable :: space(:,:,:)

!    space - the CA array
!  OUTPUT
!    space is updated
!  SIDE EFFECTS
!    none
!  DESCRIPTION
!    HX along dimension 3. See ca_mpi_hxvn1p.
!  USES
!  USED BY
!    ca_mpi_hx_all
!  SOURCE

integer :: reqs3p(2), stats(MPI_STATUS_SIZE, 2)

if ( ci(3) .ne. ucob(3) ) then

  ! Rank is image number -1.

  ! Receive from the right neighbour, tag = TAG3L
  call MPI_IRECV( space, 1, mpi_h3_RV, nei_img_R(3)-1, TAG3L,          &
    MPI_COMM_WORLD, reqs3p(1), ierr )

  ! Send to the right neighbour, tag = TAG3R
  call MPI_ISEND( space, 1, mpi_h3_RR, nei_img_R(3)-1, TAG3R,          &
    MPI_COMM_WORLD, reqs3p(2), ierr )

  call MPI_WAITALL( 2, reqs3p, stats, ierr )

end if

end subroutine ca_mpi_hxvn3p

!*roboend*


!*robodoc*s* ca_hx/ca_mpi_hx_all
!  NAME
!    ca_mpi_hx_all
!  SYNOPSIS

module subroutine ca_mpi_hx_all( space )

!  INPUT

integer( kind=iarr ), intent(inout), allocatable :: space(:,:,:)

!    space - non coarray array with CA model
!  OUTPUT
!    space is changed
!  SIDE EFFECTS
!    none
!  DESCRIPTION
!    Do all MPI HX.
!    To avoid problems I don't allow the user to call individual
!    hx routines. These are private to this module.
!    The user only calls this routine.
!  NOTE
!    ca_mpi_halo_type_create must be called prior to calling this
!    routine.
!    Note! This routine will only work if iarr is the *default* integer.
!    This is because MPI_INTEGER is used for space, as other MPI
!    integer kinds might not be implemented.
!  USES
!    ca_mpi_hxvn1m, ca_mpi_hxvn1p, ca_mpi_hxvn2m, ca_mpi_hxvn2p,
!    ca_mpi_hxvn3m, ca_mpi_hxvn3p
!  USED BY
!  SOURCE

! Make sure (some) MPI halo types have been created.
if ( .not. halo_type_created ) then
  write (*,"(a)") "ERROR ca_hx_mpi/ca_mpi_hx_all: Need to create " //  &
    "MPI types. Call ca_mpi_halo_type_create first!"
  error stop
end if

call ca_mpi_hxvn1m( space )
call ca_mpi_hxvn1p( space )
call ca_mpi_hxvn2m( space )
call ca_mpi_hxvn2p( space )
call ca_mpi_hxvn3m( space )
call ca_mpi_hxvn3p( space )

end subroutine ca_mpi_hx_all

!*roboend*


!*robodoc*s* ca_hx/ca_mpi_ising_energy
!  NAME
!    ca_mpi_ising_energy
!  SYNOPSIS

module subroutine ca_mpi_ising_energy( space, iter_sub, kernel,        &
  energy, magnet )

!  INPUT

integer( kind=iarr ), intent(inout), allocatable :: space(:,:,:)
procedure( iter_proto ) :: iter_sub
procedure( kernel_proto ) :: kernel

!       space - space array before iterations start
!    iter_sub - the subroutine performing a single CA iteration, e.g.
!             - ca_iter_tl - triple nested loop
!             - ca_iter_dc - do concurrent
!             - ca_iter_omp - OpenMP
!      kernel - a function to be called for every cell inside the loop
!  OUTPUT

integer( kind=ilrg ), intent( out ) :: energy, magnet

!    energy - Total energy of CA system
!    magnet - Total magnetisation of the CA system
!  SIDE EFFECTS
!    module array tmp_space is updated
!  DESCRIPTION
!    Calculate the total energy and the total magnetisation
!    of CA using Ising model. Note that I'm passing integers of kind
!    ilrg to MPI_INTEGER8. This should work as long as ilrg is 8 bytes
!    long. So set ilrg to selected_int_kind( 10 ).
!    This routine uses MPI_ALLREDUCE with MPI_SUM.
!    Magnetisation is defined as the fraction of the 1 spins.
!    The only valid kernel is ca_kernel_ising_ener.
!  USES
!  USED BY
!  SOURCE

integer( kind=ilrg ) :: img_energy, img_magnet

call ca_mpi_hx_all( space )    ! space updated, sync images
! tmp_space updated, local op
call iter_sub( space=space, halo=hdepth, kernel=kernel )
img_energy = sum( tmp_space( 1:sub(1), 1:sub(2), 1:sub(3) ) )
img_magnet = sum(     space( 1:sub(1), 1:sub(2), 1:sub(3) ) )

! write (*,*) "img:", this_image(), "img_energy:", img_energy, "img_magnet:", img_magnet

call MPI_ALLREDUCE( img_energy, energy, 1, MPI_INTEGER8, MPI_SUM,      &
  MPI_COMM_WORLD, ierr) 
call MPI_ALLREDUCE( img_magnet, magnet, 1, MPI_INTEGER8, MPI_SUM,      &
  MPI_COMM_WORLD, ierr) 

end subroutine ca_mpi_ising_energy

!*roboend*

end submodule ca_hx_mpi
