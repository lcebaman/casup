!*robodoc*m* CASUP/ca_hx
!  NAME
!    ca_hx
!  SYNOPSIS

!$Id: ca_hx.f90 536 2018-04-03 12:02:13Z mexas $

module ca_hx

!  DESCRIPTION
!    Module with halo exchange for CASUP, the 2nd gen
!    CA library for SUPercomputers. In this module halos are
!    separate arrays from the central part of CA, which is
!    not a coarray. 
!  AUTHOR 
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    Public subroutines: ca_spalloc, ca_halloc, ca_hdalloc, ca_hx_all,
!    ca_hx_check, ca_iter_dc, ca_iter_omp, ca_iter_tl, ca_run,
!    ca_ising_energy, ca_ising_energy_col, ca_set_space_rnd,
!    ca_mpi_halo_type_create, ca_mpi_halo_type_free,
!    ca_mpi_hx_all, ca_mpi_ising_energy,
!    ca_co_spalloc, ca_co_hx_all, ca_co_hx_check, ca_co_run,
!    ca_co_ising_energy, ca_co_netcdf, ca_co_naive_io.
!    Public pure functions: ca_kernel_copy, ca_kernel_ising,
!    ca_kernel_ising_ener.
!    Private subroutines: ca_hxvn1m, ca_hxvn1p, ca_hxvn2m, ca_hxvn2p,
!    ca_hxvn3m, ca_hxvn3p, ca_hx_ready; then in submodule ca_hx_mpi:
!    ca_mpi_hxvn1m, ca_mpi_hxvn1p, ca_mpi_hxvn2m, ca_mpi_hxvn2p,
!    ca_mpi_hxvn3m, ca_mpi_hxvn3p.
!  USES
!    cgca_m1co
!  USED BY
!  SOURCE

use cgca_m1co, only : ca_range, iarr, idef, ilrg
use mpi
implicit none

private
public :: ca_spalloc, ca_halloc, ca_hdalloc, ca_hx_all, ca_hx_check, &
  ca_iter_dc, ca_iter_omp, ca_iter_tl, ca_run, ca_kernel_copy,       &
  ca_kernel_ising, ca_kernel_ising_ener,  &
  ca_ising_energy, ca_ising_energy_col, &
  ca_set_space_rnd, &
  ca_mpi_halo_type_create, ca_mpi_halo_type_free,  &
  ca_mpi_hx_all, ca_mpi_ising_energy, &
  ca_co_spalloc, ca_co_hx_all, ca_co_hx_check, ca_co_run, &
  ca_co_ising_energy, ca_co_netcdf, ca_co_naive_io

! These are halo coarrays. Parts of the CA array are designated
! "halo", but those are not coarrays.
!
! h1 - halo along dimension 1.
! h2 - halo along dimension 1.
! h3 - halo along dimension 1.
! "minu" - minus, "plus" - plus. 
! "t" - tmp arrays.
integer( kind=iarr ), allocatable ::          &
  h1minu(:,:,:)[:,:,:], h1plus(:,:,:)[:,:,:], &
  h2minu(:,:,:)[:,:,:], h2plus(:,:,:)[:,:,:], &
  h3minu(:,:,:)[:,:,:], h3plus(:,:,:)[:,:,:], &
  tmp_space(:,:,:), mask_array(:,:,:)

integer :: hdepth,     & ! halo depth
  ci(3),         & ! coindex set of my image
  ucob(3),       & ! upper cobounds of halo coarrays
  ihsta(3),      & ! inner right halo start
  rhsta(3),      & ! outer right halo start
  rhend(3),      & ! outer right halo end
  lhsta(3),      & ! outer left halo start
  sub(3),        & ! upper bounds of space array (without halo)
  ierr,          & ! error variable
  total_cells,   & ! total number of cells in the global CA.
  nei_ci_L1(3),  & ! coindex set of left neighbour along 1
  nei_ci_R1(3),  & ! coindex set of right neighbour along 1
  nei_ci_L2(3),  & ! coindex set of left neighbour along 2
  nei_ci_R2(3),  & ! coindex set of right neighbour along 2
  nei_ci_L3(3),  & ! coindex set of left neighbour along 3
  nei_ci_R3(3),  & ! coindex set of right neighbour along 3
  nei_img_L(3),  & ! image indices for 3 left neighbours
  nei_img_R(3)     ! image indices for 3 right neighbours

! A scalar to calculate the total energy of CA
! These will not be needed when collectives can be used.
integer( kind=iarr ) :: co_energy[*], co_magnet[*]

character( len=500 ) :: errmsg

abstract interface

!****************************************************************
! Cray bug in 8.5.8!
! Should be able to use hdepth from the module via IMPORT!
! When fixed, remove "halo" and change back to use "hdepth"
! in this module, ca_hx_co and ca_hx_mpi!
!****************************************************************

  !*******************************************************************72
  ! For non-coarray CA model
  !*******************************************************************72

  ! For a kernel function
  pure function kernel_proto( space, halo, coord )
    use cgca_m1co
!    import hdepth
    implicit none
    integer, intent(in) :: halo
    integer( kind=iarr ), intent(in), contiguous ::                    &
      space( 1-halo: , 1-halo: , 1-halo: )
    integer, intent(in) :: coord(3)
    integer( kind=iarr ) kernel_proto
  end function kernel_proto

  ! For a subroutine doing a single iteration
  subroutine iter_proto( space, halo, kernel )
    use cgca_m1co
    import :: kernel_proto
    implicit none
    integer, intent( in ) :: halo
    integer( kind=iarr ), intent(in), contiguous ::                    &
      space( 1-halo: , 1-halo: , 1-halo: )
    procedure( kernel_proto ) :: kernel
  end subroutine iter_proto

  ! For a HX routine
  subroutine hx_proto( space )
    use cgca_m1co
    implicit none
    integer( kind=iarr ), intent( inout ), allocatable :: space(:,:,:)
  end subroutine hx_proto

  !*******************************************************************72
  ! For coarray CA model
  !*******************************************************************72

  ! For a HX routine
  subroutine hx_co_proto( space )
    use cgca_m1co
    implicit none
    integer( kind=iarr ), intent( inout ), allocatable ::              &
      space(:,:,:)[:,:,:]
  end subroutine hx_co_proto

end interface

! Interfaces for submodule procedures.
 
interface

  !*******************************************************************72
  ! In submodule ca_hx_mpi
  !*******************************************************************72

  ! NOTE! space is *default* integer!

  module subroutine ca_mpi_halo_type_create( space )
    integer( kind=iarr ), intent( inout ), allocatable :: space(:,:,:)
  end    subroutine ca_mpi_halo_type_create

  module subroutine ca_mpi_halo_type_free
  end    subroutine ca_mpi_halo_type_free

  module subroutine ca_mpi_hxvn1m( space )
    integer( kind=iarr ), intent(inout), allocatable :: space(:,:,:)
  end    subroutine ca_mpi_hxvn1m

  module subroutine ca_mpi_hxvn1p( space )
    integer( kind=iarr ), intent(inout), allocatable :: space(:,:,:)
  end    subroutine ca_mpi_hxvn1p

  module subroutine ca_mpi_hxvn2m( space )
    integer( kind=iarr ), intent(inout), allocatable :: space(:,:,:)
  end    subroutine ca_mpi_hxvn2m

  module subroutine ca_mpi_hxvn2p( space )
    integer( kind=iarr ), intent(inout), allocatable :: space(:,:,:)
  end    subroutine ca_mpi_hxvn2p

  module subroutine ca_mpi_hxvn3m( space )
    integer( kind=iarr ), intent(inout), allocatable :: space(:,:,:)
  end    subroutine ca_mpi_hxvn3m

  module subroutine ca_mpi_hxvn3p( space )
    integer( kind=iarr ), intent(inout), allocatable :: space(:,:,:)
  end    subroutine ca_mpi_hxvn3p

! Note! this will only work if iarr is the *default* integer.
! This is because MPI_INTEGER is used for space, as other MPI
! integer kinds might not be implemented.
  module subroutine ca_mpi_hx_all( space )
    integer( kind=iarr ), intent(inout), allocatable :: space(:,:,:)
  end    subroutine ca_mpi_hx_all

  module subroutine ca_mpi_ising_energy( space, iter_sub, kernel,      &
                                         energy, magnet )
    integer( kind=iarr ), intent(inout), allocatable :: space(:,:,:)
    procedure( iter_proto ) :: iter_sub
    procedure( kernel_proto ) :: kernel
    integer( kind=ilrg ), intent(out) :: energy, magnet
  end    subroutine ca_mpi_ising_energy

  !*******************************************************************72
  ! In submodule ca_hx_co
  !*******************************************************************72

  module subroutine ca_co_spalloc( space, c, d, ir )
    integer( kind=iarr ), intent( inout ), allocatable ::              &
      space(:,:,:) [:,:,:]
    integer, intent(in) :: c(3), d, ir(3)
  end subroutine ca_co_spalloc

  module subroutine ca_co_hx_all( space )
    integer( kind=iarr ), intent( inout ), allocatable ::              &
      space(:,:,:) [:,:,:]
  end subroutine ca_co_hx_all

  module subroutine ca_co_hx_check( space, flag )
    integer( kind=iarr ), intent( in ), allocatable ::                 &
      space(:,:,:) [:,:,:]
    integer, intent( out ) :: flag
  end subroutine ca_co_hx_check

  module subroutine ca_co_run( space, hx_sub, iter_sub, kernel, niter )
    integer( kind=iarr ), intent( inout ), allocatable ::              &
      space(:,:,:) [:,:,:]
    procedure( hx_co_proto )  :: hx_sub
    procedure( iter_proto )   :: iter_sub
    procedure( kernel_proto ) :: kernel
    integer, intent(in) :: niter
  end subroutine ca_co_run

  module subroutine ca_co_ising_energy( space, hx_sub, iter_sub,       &
       kernel, energy, magnet )
    integer( kind=iarr ), intent( inout ), allocatable ::              &
      space(:,:,:)[:,:,:]
    procedure( hx_co_proto ) :: hx_sub
    procedure( iter_proto ) :: iter_sub
    procedure( kernel_proto ) :: kernel
    integer( kind=ilrg) , intent(out) :: energy, magnet
  end subroutine ca_co_ising_energy

  module subroutine ca_co_netcdf( space, fname )
    integer( kind=iarr ), intent( in ), allocatable ::                &
      space(:,:,:) [:,:,:]
    character( len=* ), intent( in ) :: fname
  end subroutine ca_co_netcdf

  module subroutine ca_co_naive_io( coarray, fname )
    integer( kind=iarr ), intent( in ), allocatable ::                &
      coarray(:,:,:) [:,:,:]
    character( len=* ),intent( in ) :: fname
  end subroutine ca_co_naive_io

end interface

contains

!*roboend*


!*robodoc*s* ca_hx/ca_kernel_copy
!  NAME
!    ca_kernel_copy
!  SYNOPSIS

pure function ca_kernel_copy( space, halo, coord )

!  INPUTS

integer, intent( in ) :: halo
integer( kind=iarr ), intent( in ), contiguous ::                      &
  space( 1-halo: , 1-halo: , 1-halo: )
integer, intent(in) :: coord(3)

!  OUTPUT

integer( kind=iarr ) ca_kernel_copy

!  SIDE EFFECTS
!    None
!  DESCIPTION
!    This is a simplest CA kernel function. Simply copy the previous
!    state of each cell. This function is used only to test the
!    library.
!  USES
!  USED BY
!    ca_iter_tl, ca_run
!  SOURCE

ca_kernel_copy = space( coord(1), coord(2), coord(3) )

end function ca_kernel_copy

!*roboend*


!*robodoc*s* ca_hx/ca_kernel_ising
!  NAME
!    ca_kernel_ising
!  SYNOPSIS

pure function ca_kernel_ising( space, halo, coord )

!  INPUTS

integer, intent( in ) :: halo
integer( kind=iarr ), intent( in ), contiguous ::                      &
  space( 1-halo: , 1-halo: , 1-halo: )
integer, intent( in ) :: coord(3)

!  OUTPUT

integer( kind=iarr ) ca_kernel_ising

!  SIDE EFFECTS
!    None
!  DESCIPTION
!    Ising magnetisation CA kernel function.
!    The CA state is magnetic spin, either 0 (down) or 1 (up).
!    The sign of the spin (CA cell state) is changed if and only if
!    this spin has the same number of parallel and anti-parallel
!    neighbours, i.e. the sum of 6 neighbours is exactly 3 -
!    three neighbours of spin up (3*1=3) and three neighbours of
!    spin down (0). The multiplier alternates between 0 and 1 from one
!    cell to the next. This preserves the energy.
!    For more details see Sec. 2.2.3 "The Q2R rule" in:
!    B. Chopard, M. Droz "Cellular Automata Modeling of Physical
!    Systems", Cambridge, 1998.
!  USES
!  USED BY
!    ca_iter_tl, ca_iter_dc, ca_iter_omp, ca_run, ca_co_run
!  SOURCE

integer( kind=iarr ) :: n

associate( s => space, i => coord(1), j => coord(2), k => coord(3) )

n = s(i-1,j,k) + s(i+1,j,k) + s(i,j-1,k) + s(i,j+1,k) + s(i,j,k-1) +   &
    s(i,j,k+1)

if ( n .eq. 3 .and. mask_array(i,j,k) .eq. 1 ) then
  ! If the sum of 6 neighbours is exactly 3 and the mask value is 1
  ! then flip the state.
  ca_kernel_ising = 1_iarr - s(i,j,k)
else
  ! Otherwise no change
  ca_kernel_ising = s(i,j,k)
end if

end associate

end function ca_kernel_ising

!*roboend*


!*robodoc*s* ca_hx/ca_kernel_ising_ener
!  NAME
!    ca_kernel_ising_ener
!  SYNOPSIS

pure function ca_kernel_ising_ener( space, halo, coord )

!  INPUTS

integer, intent( in ) :: halo
integer( kind=iarr ), intent( in ), contiguous ::                      &
  space( 1-halo: , 1-halo: , 1-halo: )
integer, intent( in ) :: coord(3)

!  OUTPUT

integer( kind=iarr ) ca_kernel_ising_ener

!  SIDE EFFECTS
!    None
!  DESCIPTION
!    Ising magnetisation CA kernel function for energy calculation.
!    Each neighbour of the same spin adds -1 to the energy.
!    Each neighbour of the opposite spin adds +1 to the energy.
!    For more details see Sec. 2.2.3 "The Q2R rule" in:
!    B. Chopard, M. Droz "Cellular Automata Modeling of Physical
!    Systems", Cambridge, 1998.
!
!    To avoid double counting, each cell only checks links with cells
!    to its *left* in all 3 dimensions. This leaves the last cells
!    in the global CA space along each dimension,
!    which must also look to the right, i.e. to the global right halo.
!    This is not satisfactory because not all cells are processed
!    in the same way. However, I cannot think of a better algorithm.
!  USES
!  USED BY
!    ca_iter_tl, ca_iter_dc, ca_iter_omp, ca_ising_energy
!  SOURCE

integer( kind=iarr ) :: count(6)

count=0

associate( s => space, i => coord(1), j => coord(2), k => coord(3) )

! s(i,j,k) - s(i-1,j,k)                 => 0, -1, 1
! abs( s(i,j,k) - s(i-1,j,k) )          => 0, 1
! 2 * abs( s(i,j,k) - s(i-1,j,k) )      => 0, 2
! 2 * abs( s(i,j,k) - s(i-1,j,k) ) - 1  => -1, 1

! Neighbours to the left along 3 directions
count(1) = 2_iarr * abs( s(i,j,k) - s(i-1,j,k) ) - 1_iarr
count(2) = 2_iarr * abs( s(i,j,k) - s(i,j-1,k) ) - 1_iarr
count(3) = 2_iarr * abs( s(i,j,k) - s(i,j,k-1) ) - 1_iarr

! Neighbours to the right, only for the globally last cells
if ( ci(1) .eq. ucob(1) .and. i .eq. sub(1) ) then
  count(4) = 2_iarr * abs( s(i,j,k) - s(i+1,j,k) ) - 1_iarr
end if
if ( ci(2) .eq. ucob(2) .and. j .eq. sub(2) ) then
  count(5) = 2_iarr * abs( s(i,j,k) - s(i,j+1,k) ) - 1_iarr
end if
if ( ci(3) .eq. ucob(3) .and. k .eq. sub(3) ) then
  count(6) = 2_iarr * abs( s(i,j,k) - s(i,j,k+1) ) - 1_iarr
end if

ca_kernel_ising_ener = sum( count )

! write (*,"(a,4(i0,tr1))") "i,j,k,energy:", i,j,k,ca_kernel_ising_ener

end associate

end function ca_kernel_ising_ener

!*roboend*


!*robodoc*s* ca_hx/ca_spalloc
!  NAME
!    ca_spalloc
!  SYNOPSIS

subroutine ca_spalloc( space, c, d )

!  INPUT

integer( kind=iarr ), allocatable, intent(inout) :: space(:,:,:) 
integer, intent(in) :: c(3), d

!    space - CA array to allocate, with halos!
!        c - array with space dimensions
!        d - depth of the halo layer
!  OUTPUT
!    space is allocated and set to zero.
!  SIDE EFFECTS
!    none
!  DESCRIPTION
!    This routine allocates the CA array space, with halos of depth d.
!    Also save some vars in this module for future.
!    Also, on first call, allocate a module work space array
!    (tmp_space) of the same mold as space.
!    This array is used in CA iterations later.
!  USES
!  USED BY
!  SOURCE

if ( allocated( space ) ) then
  write (*,*) "WARN: ca_hx/ca_spalloc: image:", this_image(), "space", &
              "already allocated, deallocating!"
  deallocate( space, stat=ierr, errmsg=errmsg )
  if ( ierr .ne. 0 ) then
    write (*,*) "ERROR: ca_hx/ca_spalloc: deallocate( space ), ierr:", &
                ierr, "errmsg:", trim(errmsg)
    error stop
  end if
end if

allocate( space( 1-d:c(1)+d, 1-d:c(2)+d, 1-d:c(3)+d ), source=0_iarr,  &
          stat=ierr, errmsg=errmsg )
if ( ierr .ne. 0 ) then
  write (*,*) "ERROR: ca_hx/ca_spalloc: allocate( space ), ierr:",     &
              ierr, "errmsg:", trim(errmsg)
  error stop
end if

! Calculate once and keep forever
! Need a picture to explain these halo vars.
! The picture is only for dimension 1.
! The same diagram can be drawn for dimensions 2 and 3.
!
! hdepth 
! |    | 
! +----+---------------------+----+--
! |    |                     |    | ^ hdepth
! |    |                     |    | v
! +----+---------------------+----+--
! |<-lhsta  |           |    |    |
! |   0|    |           |    |    |
! |    |1   |           |    |    |
! | hdepth->|           |    |    |
! |    |    |           |<-ihsta  |
! |    |    |           sub->|    |
! |    |    |           |    |<-rhsta
! |    |    |           |  rhend->|
! |    |    |           |    |    |
! +----+---------------------+----+
! |    |                     |    |
! |    |                     |    |
! +----+---------------------+----+

hdepth = d
   sub = c
 ihsta = sub - hdepth + 1
 rhsta = sub + 1
 rhend = sub + hdepth
 lhsta = -hdepth + 1

! Total number of cells in the global CA
total_cells = sub(1)*sub(2)*sub(3)*num_images()

! Tmp space array for CA iterations, allocated on the first call
if ( .not. allocated( tmp_space ) ) then
  allocate( tmp_space, source=space, stat=ierr, errmsg=errmsg )
  if ( ierr .ne. 0 ) then
    write (*,"(a,i0,a,a)") "ERROR: ca_hx/ca_spalloc: " //              &
      "allocate( tmp_space ), ierr: ", ierr, " errmsg: ", trim(errmsg)
    error stop
  end if
end if

! Mask array, no halos, allocated on the first call
! The values are set in ca_halloc.
if ( .not. allocated( mask_array ) ) then
  allocate( mask_array( c(1), c(2), c(3) ), stat=ierr, errmsg=errmsg )
  if ( ierr .ne. 0 ) then
    write (*,"(a,i0,a,a)") "ERROR: ca_hx/ca_spalloc: " //              &
      "allocate( mask_array ) ", "ierr:", ierr, "errmsg:", trim(errmsg)
    error stop
  end if
end if

end subroutine ca_spalloc

!*roboend*


!*robodoc*s* ca_hx/ca_halloc
!  NAME
!    ca_halloc
!  SYNOPSIS

subroutine ca_halloc( ir )

!  INPUT

integer :: ir(3)

!    ir - codimensions
!  OUTPUT
!    none
!  SIDE EFFECTS
!    halo coarrays (module variables) are allocated
!  DESCRIPTION
!    This routine allocates 6 halo coarrays for von Neumann
!    6-neighbourhood. The coarrays are module variables.
!    Coarray allocation is an implicit sync.
!    Halos have depth hdepth.
!    Halo coarrays have the same cobounds on all images.
!  NOTES
!    All images must call this routine!
!  USES
!  USED BY
!  SOURCE

integer :: i,j,k

main: associate( d => hdepth, c => sub )

allocate( h1minu( d, c(2), c(3) ) [ ir(1), ir(2), * ], stat=ierr,      &
          errmsg=errmsg )
if ( ierr .ne. 0 ) then
  write (*,*) "ERROR: ca_hx/ca_halloc: allocate( h1minu ). ierr: ",    &
              ierr, "errmsg:", trim(errmsg)
  error stop
end if

allocate( h1plus( d, c(2), c(3) ) [ ir(1), ir(2), * ], stat=ierr,      &
          errmsg=errmsg )
if ( ierr .ne. 0 ) then
  write (*,*) "ERROR: ca_hx/ca_halloc: allocate( h1plus ). ierr: ",    &
              ierr, "errmsg:", trim(errmsg)
  error stop
end if

allocate( h2minu( c(1), d, c(3) ) [ ir(1), ir(2), * ], stat=ierr,      &
          errmsg=errmsg )
if ( ierr .ne. 0 ) then
  write (*,*) "ERROR: ca_hx/ca_halloc: allocate( h2minu ). ierr: ",    &
              ierr, "errmsg:", trim(errmsg)
  error stop
end if

allocate( h2plus( c(1), d, c(3) ) [ ir(1), ir(2), * ], stat=ierr,      &
          errmsg=errmsg )
if ( ierr .ne. 0 ) then
  write (*,*) "ERROR: ca_hx/ca_halloc: allocate( h2plus ). ierr: ",    &
              ierr, "errmsg:", trim(errmsg)
  error stop
end if

allocate( h3minu( c(1), c(2), d ) [ ir(1), ir(2), * ], stat=ierr,      &
          errmsg=errmsg )
if ( ierr .ne. 0 ) then
  write (*,*) "ERROR: ca_hx/ca_halloc: allocate( h3minu ). ierr: ",    &
              ierr, "errmsg:", trim(errmsg)
  error stop
end if

allocate( h3plus( c(1), c(2), d ) [ ir(1), ir(2), * ], stat=ierr,      &
          errmsg=errmsg)
if ( ierr .ne. 0 ) then
  write (*,*) "ERROR: ca_hx/ca_halloc: allocate( h3plus ). ierr: ",    &
              ierr, "errmsg:", trim(errmsg)
  error stop
end if

! Calculate once and keep forever
  ci = this_image( h1minu )
ucob = ucobound(   h1minu )

! Now can set mask_array. The mask array must reflect
! the global CA space, i.e. must not be affected by partitioning
! of the model into images. This means the mask array must
! depend on coindex set of this image.
do concurrent( i=1:c(1), j=1:c(2), k=1:c(3) )
  mask_array(i,j,k) = int( mod( (i+j+k + ( ci(1)-1)*c(1) +             &
    (ci(2)-1)*c(2) + (ci(3)-1)*c(3) ) , 2 ), kind=iarr ) 
end do

end associate main

! Calculate coindex sets and image numbers for the 6 neighbours
! Coindex sets
  nei_ci_L1 = (/ ci(1)-1, ci(2), ci(3) /)
  nei_ci_R1 = (/ ci(1)+1, ci(2), ci(3) /)
  nei_ci_L2 = (/ ci(1), ci(2)-1, ci(3) /)
  nei_ci_R2 = (/ ci(1), ci(2)+1, ci(3) /)
  nei_ci_L3 = (/ ci(1), ci(2), ci(3)-1 /)
  nei_ci_R3 = (/ ci(1), ci(2), ci(3)+1 /)

! Image index
  nei_img_L(1) = image_index( h1plus, nei_ci_L1 )
  nei_img_R(1) = image_index( h1plus, nei_ci_R1 )
  nei_img_L(2) = image_index( h1plus, nei_ci_L2 )
  nei_img_R(2) = image_index( h1plus, nei_ci_R2 )
  nei_img_L(3) = image_index( h1plus, nei_ci_L3 )
  nei_img_R(3) = image_index( h1plus, nei_ci_R3 )

end subroutine ca_halloc

!*roboend*


!*robodoc*s* ca_hx/ca_hdalloc
!  NAME
!    ca_hdalloc
!  SYNOPSIS

subroutine ca_hdalloc

!  INPUT
!    none
!  OUTPUT
!    none
!  SIDE EFFECTS
!    halo coarrays (module variables) are deallocated
!  DESCRIPTION
!    This routine deallocates 6 halo coarrays for von Neumann
!    6-neighbourhood. The coarrays are module variables.
!    Coarray deallocation is an implicit sync.
!    Halo coarrays have the same cobounds on all images.
!  NOTES
!    All images must call this routine!
!  USES
!  USED BY
!  SOURCE

deallocate( h1minu, stat=ierr, errmsg=errmsg )
if ( ierr .ne. 0 ) then
  write (*,*) "ERROR: ca_hx/ca_hdalloc: deallocate( h1minu ). ierr: ", &
              ierr, "errmsg:", trim(errmsg)
  error stop
end if

deallocate( h1plus, stat=ierr, errmsg=errmsg )
if ( ierr .ne. 0 ) then
  write (*,*) "ERROR: ca_hx/ca_hdalloc: deallocate( h1plus ). ierr: ", &
              ierr, "errmsg:", trim(errmsg)
  error stop
end if

deallocate( h2minu, stat=ierr, errmsg=errmsg )
if ( ierr .ne. 0 ) then
  write (*,*) "ERROR: ca_hx/ca_hdalloc: deallocate( h2minu ). ierr: ", &
              ierr, "errmsg:", trim(errmsg)
  error stop
end if

deallocate( h2plus, stat=ierr, errmsg=errmsg )
if ( ierr .ne. 0 ) then
  write (*,*) "ERROR: ca_hx/ca_hdalloc: deallocate( h2plus ). ierr: ", &
              ierr, "errmsg:", trim(errmsg)
  error stop
end if

deallocate( h3minu, stat=ierr, errmsg=errmsg )
if ( ierr .ne. 0 ) then
  write (*,*) "ERROR: ca_hx/ca_hdalloc: deallocate( h3minu ). ierr: ", &
              ierr, "errmsg:", trim(errmsg)
  error stop
end if

deallocate( h3plus, stat=ierr, errmsg=errmsg)
if ( ierr .ne. 0 ) then
  write (*,*) "ERROR: ca_hx/ca_hdalloc: deallocate( h3plus ). ierr: ", &
              ierr, "errmsg:", trim(errmsg)
  error stop
end if

end subroutine ca_hdalloc

!*roboend*


!*robodoc*s* ca_hx/ca_hx_ready
!  NAME
!    ca_hx_ready
!  SYNOPSIS

subroutine ca_hx_ready( space )

!  INPUT

integer( kind=iarr ), intent(in), allocatable :: space(:,:,:)

!    space - CA model array
!  OUTPUT
!    none
!  SIDE EFFECTS
!    halo coarrays are updated
!  DESCRIPTION
!    This is the first step of the 2 step hx process.
!    I prepare my 6 coarray halos for use by my neighbours.
!    These are copies of the slabs of space array, i.e. top/bottom,
!    left/right and front/back.
!    Note that these are not space array halos!
!    Refer to the schematic in ca_spalloc for details.
!    When done, I declare that my 6 coarray halos "ready" to be read -
!    ready for remote calls.
!    This is a local routine, which must be run before any
!    of hx routines are called.
!    This routine is called when all cells on my image have been
!    processed, so doesn't make sense to separate separate assignments
!    into different routines - no performance gain at all. 
!  USES
!  USED BY
!    ca_hx_all
!  SOURCE

if ( ci(1) .ne. 1       ) then
  h1minu(:,:,:) = space( 1 : hdepth        , 1:sub(2) , 1:sub(3) )
end if
if ( ci(1) .ne. ucob(1) ) then
  h1plus(:,:,:) = space( ihsta(1) : sub(1) , 1:sub(2) , 1:sub(3) )
end if
if ( ci(2) .ne. 1       ) then
  h2minu(:,:,:) = space( 1:sub(1) , 1 : hdepth        , 1:sub(3) )
end if
if ( ci(2) .ne. ucob(2) ) then
  h2plus(:,:,:) = space( 1:sub(1) , ihsta(2) : sub(2) , 1:sub(3) )
end if
if ( ci(3) .ne. 1       ) then
  h3minu(:,:,:) = space( 1:sub(1) , 1:sub(2) , 1 : hdepth        )
end if
if ( ci(3) .ne. ucob(3) ) then
  h3plus(:,:,:) = space( 1:sub(1) , 1:sub(2) , ihsta(3) : sub(3) )
end if

end subroutine ca_hx_ready

!*roboend*


!*robodoc*s* ca_hx/ca_hxvn1m
!  NAME
!    ca_hxvn1m
!  SYNOPSIS

subroutine ca_hxvn1m( space )

!  INPUT

integer( kind=iarr ), intent(inout), allocatable :: space(:,:,:)

!    space - the CA array
!  OUTPUT
!    space is updated
!  SIDE EFFECTS
!    none
!  DESCRIPTION
!    An image updates its space array halo layer (left side)
!    along direction 1 from a coarray halo (h1plus) on an image
!    which is 1 lower along codimension 1.
!  USES
!  USED BY
!    ca_hx_all
!  SOURCE

! coindex set and the image number of the neighbour
!integer :: i(3), n

if ( ci(1) .ne. 1 ) then
  ! This is the neighbour
!  i = (/ ci(1)-1, ci(2), ci(3) /) ! neighbour's coindex set
!  n = image_index( h1plus, i )    ! neighbour image number
!  sync images( n )
!  space( lhsta(1):0, 1:sub(2), 1:sub(3) ) =                            &
!    h1plus(:,:,:) [ i(1), i(2), i(3) ]
  sync images( nei_img_L(1) )
  space( lhsta(1):0, 1:sub(2), 1:sub(3) ) =                            &
    h1plus(:,:,:) [ nei_ci_L1(1), nei_ci_L1(2), nei_ci_L1(3) ]
end if

end subroutine ca_hxvn1m

!*roboend*


!*robodoc*s* ca_hx/ca_hxvn1p
!  NAME
!    ca_hxvn1p
!  SYNOPSIS

subroutine ca_hxvn1p( space )

!  INPUT

integer( kind=iarr ), intent(inout), allocatable :: space(:,:,:)

!    space - the CA array
!  OUTPUT
!    space is updated
!  SIDE EFFECTS
!    none
!  DESCRIPTION
!    An image updates its space array halo layer (right side)
!    along direction 1 from a coarray halo (h1minu) on an image
!    which is 1 higher along codimension 1.
!  USES
!  USED BY
!    ca_hx_all
!  SOURCE

! coindex set and the image number of the neighbour
!integer :: i(3), n

if ( ci(1) .ne. ucob(1) ) then
!  ! This is the neighbour
!  i = (/ ci(1)+1, ci(2), ci(3) /) ! neighbour's coindex set
!  n = image_index( h1plus, i )    ! neighbour image number
!  sync images( n )
!  space( rhsta(1) : rhend(1), 1:sub(2) , 1:sub(3) ) =                  &
!    h1minu(:,:,:) [ i(1), i(2), i(3) ]
  sync images( nei_img_R(1) )
  space( rhsta(1) : rhend(1), 1:sub(2) , 1:sub(3) ) =                  &
    h1minu(:,:,:) [ nei_ci_R1(1), nei_ci_R1(2), nei_ci_R1(3) ]
end if

end subroutine ca_hxvn1p

!*roboend*


!*robodoc*s* ca_hx/ca_hxvn2m
!  NAME
!    ca_hxvn2m
!  SYNOPSIS

subroutine ca_hxvn2m( space )

!  INPUT

integer( kind=iarr ), intent(inout), allocatable :: space(:,:,:)

!    space - the CA array
!  OUTPUT
!    space is updated
!  SIDE EFFECTS
!    none
!  DESCRIPTION
!    An image updates its space array halo layer (left side)
!    along direction 2 from a coarray halo (h1plus) on an image
!    which is 1 lower along codimension 2.
!  USES
!  USED BY
!    ca_hx_all
!  SOURCE

! coindex set and the image number of the neighbour
!integer :: i(3), n

if ( ci(2) .ne. 1 ) then
!  ! This is the neighbour
!  i = (/ ci(1), ci(2)-1, ci(3) /) ! neighbour's coindex set
!  n = image_index( h2plus, i )    ! neighbour image number
!  sync images( n )
!  space( 1:sub(1) , lhsta(2) : 0, 1:sub(3) ) =                         &
!    h2plus(:,:,:) [ i(1), i(2), i(3) ]
  sync images( nei_img_L(2) )
  space( 1:sub(1) , lhsta(2) : 0, 1:sub(3) ) =                         &
    h2plus(:,:,:) [ nei_ci_L2(1), nei_ci_L2(2), nei_ci_L2(3) ]
end if

end subroutine ca_hxvn2m

!*roboend*


!*robodoc*s* ca_hx/ca_hxvn2p
!  NAME
!    ca_hxvn2p
!  SYNOPSIS

subroutine ca_hxvn2p( space )

!  INPUT

integer( kind=iarr ), intent(inout), allocatable :: space(:,:,:)

!    space - the CA array
!  OUTPUT
!    space is updated
!  SIDE EFFECTS
!    none
!  DESCRIPTION
!    An image updates its space array halo layer (right side)
!    along direction 2 from a coarray halo (h2minu) on an image
!    which is 1 higher along codimension 1.
!  USES
!  USED BY
!    ca_hx_all
!  SOURCE

! coindex set and the image number of the neighbour
!integer :: i(3), n 

if ( ci(2) .ne. ucob(2) ) then
!  ! This is the neighbour
!  i = (/ ci(1), ci(2)+1, ci(3) /) ! neighbour's coindex set
!  n = image_index( h2plus, i )    ! neighbour image number
!  sync images( n )
!  space( 1:sub(1) , rhsta(2) : rhend(2) , 1:sub(3) ) =                 &
!    h2minu(:,:,:) [ i(1), i(2), i(3) ]
  sync images( nei_img_R(2) )
  space( 1:sub(1) , rhsta(2) : rhend(2) , 1:sub(3) ) =                 &
    h2minu(:,:,:) [ nei_ci_R2(1), nei_ci_R2(2), nei_ci_R2(3) ]
end if

end subroutine ca_hxvn2p

!*roboend*


!*robodoc*s* ca_hx/ca_hxvn3m
!  NAME
!    ca_hxvn3m
!  SYNOPSIS

subroutine ca_hxvn3m( space )

!  INPUT

integer( kind=iarr ), intent(inout), allocatable :: space(:,:,:)

!    space - the CA array
!  OUTPUT
!    space is updated
!  SIDE EFFECTS
!    none
!  DESCRIPTION
!    An image updates its space array halo layer (left side)
!    along direction 3 from a coarray halo (h3plus) on an image
!    which is 1 lower along codimension 3.
!  USES
!  USED BY
!    ca_hx_all
!  SOURCE

! coindex set and the image number of the neighbour
!integer :: i(3), n

if ( ci(3) .ne. 1 ) then
!  ! This is the neighbour
!  i = (/ ci(1), ci(2), ci(3)-1 /) ! neighbour's coindex set
!  n = image_index( h3plus, i )    ! neighbour image number
!  sync images( n )
!  space( 1:sub(1) , 1:sub(2) , lhsta(3) : 0 ) =                        &
!    h3plus(:,:,:) [ i(1), i(2), i(3) ]
  sync images( nei_img_L(3) )
  space( 1:sub(1) , 1:sub(2) , lhsta(3) : 0 ) =                        &
    h3plus(:,:,:) [ nei_ci_L3(1), nei_ci_L3(2), nei_ci_L3(3) ]
end if

end subroutine ca_hxvn3m

!*roboend*


!*robodoc*s* ca_hx/ca_hxvn3p
!  NAME
!    ca_hxvn3p
!  SYNOPSIS

subroutine ca_hxvn3p( space )

!  INPUT

integer( kind=iarr ), intent(inout), allocatable :: space(:,:,:)

!    space - the CA array
!  OUTPUT
!    space is updated
!  SIDE EFFECTS
!    none
!  DESCRIPTION
!    An image updates its space array halo layer (right side)
!    along direction 3 from a coarray halo (h3minu) on an image
!    which is 1 higher along codimension 3.
!  USES
!  USED BY
!    ca_hx_all
!  SOURCE

! coindex set and the image number of the neighbour
!integer :: i(3), n

if ( ci(3) .ne. ucob(3) ) then
!  ! This is the neighbour
!  i = (/ ci(1), ci(2), ci(3)+1 /) ! neighbour's coindex set
!  n = image_index( h3plus, i )    ! neighbour image number
!  sync images( n )
!  space( 1:sub(1) , 1:sub(2) , rhsta(3) : rhend(3) ) =                 &
!    h3minu(:,:,:) [ i(1), i(2), i(3) ]
  sync images( nei_img_R(3) )
  space( 1:sub(1) , 1:sub(2) , rhsta(3) : rhend(3) ) =                 &
    h3minu(:,:,:) [ nei_ci_R3(1), nei_ci_R3(2), nei_ci_R3(3) ]
end if

end subroutine ca_hxvn3p

!*roboend*


!*robodoc*s* ca_hx/ca_hx_all
!  NAME
!    ca_hx_all
!  SYNOPSIS

subroutine ca_hx_all( space )

!  INPUT

integer( kind=iarr ), intent(inout), allocatable :: space(:,:,:)

!    space - non coarray array with CA model
!  OUTPUT
!    none
!  SIDE EFFECTS
!    halo coarrays are updated
!  DESCRIPTION
!    HX is a 2-step process. In step 1 I copy parts of my space array,
!    not space halos!, into my halo coarrays. This is done by
!    ca_hx_ready. Step 2 is remote comms - I pull neighbour halo
!    coarrays into my local space array halo sections.
!
!    Important! the 6 routines
!    must be called in the same order on all images to avoid deadlocks.
!    So to make it fool proof I don't allow the user to call individual
!    hx routines. These are private to this module.
!    The user only calls this routine.
!    The same applies to ca_hx_ready, because it prepares coarray halos.
!  USES
!    ca_hx_ready, ca_hxvn1m, ca_hxvn1p, ca_hxvn2m, ca_hxvn2p,
!    ca_hxvn3m, ca_hxvn3p
!  USED BY
!  SOURCE

! HX is a 2-step process
! Step 1
call ca_hx_ready( space ) ! halos are prepared and ready for use

! Step 2
call ca_hxvn1m( space )   ! halos are used - controled by sync images
call ca_hxvn1p( space )
call ca_hxvn2m( space )
call ca_hxvn2p( space )
call ca_hxvn3m( space )
call ca_hxvn3p( space )

end subroutine ca_hx_all

!*roboend*


!*robodoc*s* ca_hx/ca_hx_check
!  NAME
!    ca_hx_check
!  SYNOPSIS

subroutine ca_hx_check( space, flag )

!  INPUT

integer( kind=iarr ), intent(in), allocatable :: space(:,:,:)

!    space - CA array
!  OUTPUT

integer, intent( out ) :: flag

!    flag .eq. 0 - check passed
!    flag .ne. 0 - check failed
!  SIDE EFFECTS
!    none
!  DESCRIPTION
!    This routine is of very limited use. It checks hx code for only
!    a single case - all images set all their cell values to
!    this_image() and then a single hx step is done.
!    In this case I'm sure what's in my halos. So just check for this
!    and either return flag=0 or flag > 0 indicating which halo failed.
!    Fail values:
!      0  - pass, no failures
!      1  - test 1 failed
!      2  - test 2 failed
!      3  - tests 1,2 failed
!      4  - test 3 failed
!      5  - tests 1,3 failed
!      6  - tests 2,3 failed
!      7  - tests 1,2,3 failed
!      8  - test 4 failed
!      9  - tests 1,4 failed
!      10 - tests 2,4 failed
!      11 - tests 1,2,4 failed
!      12 - tests 3,4 failed
!      13 - tests 1,3,4 failed
!      14 - tests 2,3,4 failed
!      15 - tests 1,2,3,4 failed
!      16 - test 5 failed
!    and so on.
!  USES
!  USED BY
!  SOURCE

! coindex set and the image number of the neighbour
integer :: i(3), n

flag = 0

! Test 1
if ( ci(1) .ne. 1 ) then
  ! This is the neighbour
  i = (/ ci(1)-1, ci(2), ci(3) /) ! neighbour's coindex set
  n = image_index( h1plus, i )    ! neighbour image number
  if ( any( space( lhsta(1):0,        1:sub(2), 1:sub(3) ) .ne. n ) )  &
    flag = flag + 1
end if

! Test 2
if ( ci(1) .ne. ucob(1) ) then
  ! This is the neighbour
  i = (/ ci(1)+1, ci(2), ci(3) /) ! neighbour's coindex set
  n = image_index( h1plus, i )    ! neighbour image number
  if ( any( space( rhsta(1):rhend(1), 1:sub(2), 1:sub(3) ) .ne. n ) )  &
    flag = flag + 2
end if

! Test 3
if ( ci(2) .ne. 1 ) then
  ! This is the neighbour
  i = (/ ci(1), ci(2)-1, ci(3) /) ! neighbour's coindex set
  n = image_index( h1plus, i )    ! neighbour image number
  if ( any( space( 1:sub(1), lhsta(2):0,        1:sub(3) ) .ne. n ) )  &
    flag = flag + 4
end if

! Test 4
if ( ci(2) .ne. ucob(2) ) then
  ! This is the neighbour
  i = (/ ci(1), ci(2)+1, ci(3) /) ! neighbour's coindex set
  n = image_index( h1plus, i )    ! neighbour image number
  if ( any( space( 1:sub(1), rhsta(2):rhend(2), 1:sub(3) ) .ne. n ) )  &
    flag = flag + 8
end if

! Test 5
if ( ci(3) .ne. 1 ) then
  ! This is the neighbour
  i = (/ ci(1), ci(2), ci(3)-1 /) ! neighbour's coindex set
  n = image_index( h1plus, i )    ! neighbour image number
  if ( any( space( 1:sub(1), 1:sub(2), lhsta(3):0 )        .ne. n ) )  &
    flag = flag + 16
end if

! Test 6
if ( ci(3) .ne. ucob(3) ) then
  ! This is the neighbour
  i = (/ ci(1), ci(2), ci(3)+1 /) ! neighbour's coindex set
  n = image_index( h1plus, i )     ! neighbour image number
  if ( any( space( 1:sub(1), 1:sub(2), rhsta(3):rhend(3) ) .ne. n ) )  &
    flag = flag + 32
end if

end subroutine ca_hx_check

!*roboend*


!*robodoc*s* ca_hx/ca_iter_tl
!  NAME
!    ca_iter_tl
!  SYNOPSIS

subroutine ca_iter_tl( space, halo, kernel )

!  INPUT

integer, intent( in ) :: halo
integer( kind=iarr ), intent( in ), contiguous ::                      &
  space( 1-halo: , 1-halo: , 1-halo: )
procedure( kernel_proto ) :: kernel

!    space - CA array at the end of the previous iteration
!  OUTPUT
!    none
!  SIDE EFFECTS
!    module array tmp_space, CA state at the end of this iteration,
!    is updated
!  DESCRIPTION
!    This routine does a single CA iteration with triple nested loops.
!    The space array from the previous iteration is read only.
!    The space array at the end of this iteration, tmp_space, from
!    this module, is written only.
!  USES
!  USED BY
!    ca_run, ca_co_run, ca_ising_energy
!  SOURCE

integer :: i, j, k

! Do not include halo! So start at 1 and end at sub.
do k = 1, sub(3)
do j = 1, sub(2)
do i = 1, sub(1)

  tmp_space( i,j,k ) = kernel( space = space, halo = hdepth,           &
    coord = (/ i , j , k /) )

!  write (*, "(a,i0,a,4(i0,tr1))") &
!    "img: ", this_image(), " i,j,k,energy: ", i, j, k, tmp_space(i,j,k)

end do
end do
end do

end subroutine ca_iter_tl
!*roboend*


!*robodoc*s* ca_hx/ca_iter_dc
!  NAME
!    ca_iter_dc
!  SYNOPSIS

subroutine ca_iter_dc( space, halo, kernel )

!  INPUT

integer, intent( in ) :: halo
integer( kind=iarr ), intent( in ), contiguous ::                      &
  space( 1-halo: , 1-halo: , 1-halo: )
procedure( kernel_proto ) :: kernel

!    space - CA array at the end of the previous iteration
!  OUTPUT
!    none
!  SIDE EFFECTS
!    module array tmp_space, CA state at the end of this iteration,
!    is updated
!  DESCRIPTION
!    This routine does a single CA iteration with DO CONCURRENT.
!    The space array from the previous iteration is read only.
!    The space array at the end of this iteration, tmp_space, from
!    this module, is written only.
!  USES
!  USED BY
!    ca_run, ca_ising_energy
!  SOURCE

integer :: i, j, k

! Do not include halo! So start at 1 and end at sub.
do concurrent( k = 1:sub(3), j = 1:sub(2), i = 1:sub(1) )
  tmp_space( i,j,k ) = kernel( space = space, halo = hdepth,           &
    coord = (/ i , j , k /) )
end do

end subroutine ca_iter_dc
!*roboend*


!*robodoc*s* ca_hx/ca_iter_omp
!  NAME
!    ca_iter_omp
!  SYNOPSIS

subroutine ca_iter_omp( space, halo, kernel )

!  INPUT

integer, intent( in ) :: halo
integer( kind=iarr ), intent( in ), contiguous ::                      &
  space( 1-halo: , 1-halo: , 1-halo: )
procedure( kernel_proto ) :: kernel

!    space - CA array at the end of the previous iteration
!  OUTPUT
!    none
!  SIDE EFFECTS
!    module array tmp_space, CA state at the end of this iteration,
!    is updated
!  DESCRIPTION
!    This routine does a single CA iteration with triple nested loops.
!    The space array from the previous iteration is read only.
!    The space array at the end of this iteration, tmp_space, from
!    this module, is written only.
!  USES
!  USED BY
!    ca_run, ca_ising_energy
!  SOURCE

integer :: i, j, k

! Do not include halo! So start at 1 and end at sub.

!$omp parallel do default(none)                                        &
!$omp private(i,j,k) shared(sub,space,hdepth,tmp_space)
! !$omp collapse(3)
do k = 1, sub(3)
do j = 1, sub(2)
do i = 1, sub(1)
  tmp_space( i,j,k ) = kernel( space = space, halo = hdepth,           &
    coord = (/ i , j , k /) )
end do
end do
end do
!$omp end parallel do

end subroutine ca_iter_omp
!*roboend*


!*robodoc*s* ca_hx/ca_run
!  NAME
!    ca_run
!  SYNOPSIS

subroutine ca_run( space, hx_sub, iter_sub, kernel, niter )

!  INPUT

integer( kind=iarr ), intent(inout), allocatable :: space(:,:,:)
procedure( hx_proto ) :: hx_sub
procedure( iter_proto ) :: iter_sub
procedure( kernel_proto ) :: kernel
integer, intent(in) :: niter

!       space - space array before iterations start
!      hx_sub - HX routine, e.g.
!             - ca_hx_all
!             _ ca_mpi_hx_all
!    iter_sub - the subroutine performing a single CA iteration, e.g.
!             - ca_iter_tl - triple nested loop
!             - ca_iter_dc - do concurrent
!             - ca_iter_omp - OpenMP
!      kernel - a function to be called for every cell inside the loop
!        iter - number of iterations to do
!  OUTPUT
!    space - CA array at the end of niter iterations
!  SIDE EFFECTS
!    module array tmp_space is updated
!  DESCRIPTION
!    This is a driver routine for CA iterations. HX is done
!    before each iteration. Then a given number of iterations
!    is performed with a given routine and a given kernel.
!    One iteration is really 2 iterations: odd and even,
!    Hence the upper loop limit is 2*niter.
!  USES
!  USED BY
!  SOURCE

integer :: i

tmp_space = space

do i = 1, 2*niter
  call hx_sub( space )                ! space updated, with HX
  ! tmp_space updated, local op
  call iter_sub( space=space, halo=hdepth, kernel=kernel )
  space = tmp_space                   ! local op
  mask_array = 1_iarr - mask_array    ! Flip the mask array
end do

end subroutine ca_run

!*roboend*


!*robodoc*s* ca_hx/ca_ising_energy
!  NAME
!    ca_ising_energy
!  SYNOPSIS

subroutine ca_ising_energy( space, iter_sub, kernel, energy, magnet )

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

integer, intent(out) :: energy, magnet

!    energy - Total energy of CA system
!    magnet - Total magnetisation of the CA system
!  SIDE EFFECTS
!    module array tmp_space is updated
!  DESCRIPTION
!    Calculate the total energy and the total magnetisation
!    of CA using Ising model.
!    This routine does not use collectives.
!    Magnetisation is defined as the fraction of the 1 spins.
!    The only valid kernel is ca_kernel_ising_ener.
!  USES
!  USED BY
!  SOURCE

integer( kind=iarr) :: my_energy, my_magnet

co_energy = 0
co_magnet = 0

call ca_hx_all( space )        ! space updated, sync images
call iter_sub( space, hdepth, kernel ) ! tmp_space updated, local op
my_energy = sum( tmp_space( 1:sub(1), 1:sub(2), 1:sub(3) ) ) 
my_magnet = sum(     space( 1:sub(1), 1:sub(2), 1:sub(3) ) ) 

!write (*,"(a,i0,a,999i3)") "img: ", this_image(), " tmp_space: ", tmp_space
!write (*,*) "img:", this_image(), "my energy:", my_energy 

! This is a tmp version on systems with no CO_SUM!
! Change on Cray!
!
! Image 1 calculates the total values in a *coarray* variable!
critical
  co_energy[1] = co_energy[1] + my_energy
  co_magnet[1] = co_magnet[1] + my_magnet
end critical

! I read the total values from image 1
! Magnetisation is real value, scaled by the total number of
! cells in the global model.
sync all
energy = co_energy[1]
magnet = co_magnet[1]

! Better do the scaling in the end user program.
! I want the routine to return an integer magnetisation
! so that the results can be exactly reproducible.
!magnet = real( co_magnet[1] ) / real( total_cells )

end subroutine ca_ising_energy

!*roboend*


!*robodoc*s* ca_hx/ca_ising_energy_col
!  NAME
!    ca_ising_energy_col
!  SYNOPSIS

subroutine ca_ising_energy_col( space, iter_sub, kernel, energy, magnet)

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

integer( kind=ilrg) , intent(out) :: energy, magnet

!    energy - Total energy of CA system
!    magnet - Total magnetisation of the CA system
!  SIDE EFFECTS
!    module array tmp_space is updated
!  DESCRIPTION
!    Calculate the total energy and the total magnetisation
!    of CA using Ising model. These integer values might be very large
!    so I'm using a large integer kind (ilrg).
!    This routine uses collectives.
!    Magnetisation is defined as the fraction of the 1 spins.
!    The only valid kernel is ca_kernel_ising_ener.
!  USES
!  USED BY
!  SOURCE

call ca_hx_all( space )        ! space updated, sync images
call iter_sub( space, hdepth, kernel ) ! tmp_space updated, local op
energy = sum( tmp_space( 1:sub(1), 1:sub(2), 1:sub(3) ) ) 
magnet = sum(     space( 1:sub(1), 1:sub(2), 1:sub(3) ) ) 
call co_sum( energy )
call co_sum( magnet )

end subroutine ca_ising_energy_col

!*roboend*


!*robodoc*s* ca_hx/ca_set_space_rnd
!  NAME
!    ca_set_space_rnd
!  SYNOPSIS

subroutine ca_set_space_rnd( seed, frac1, space )

!  INPUT

integer :: seed(:)
real :: frac1

!    seed - RND seed array
!   frac1 - fraction of "1" cells, spin up cells, [0..1].
!  OUTPUT

integer( kind=iarr ), intent( inout ), contiguous ::                   &
  space( 1-hdepth: , 1-hdepth:, 1-hdepth: )

!   space - CA array is set to reproducible RND values, and
!           in a way that the global CA array is not dependent
!           on the partition or the number of images 
!  DESCRIPTION
!    Must call this routine after ca_halloc.
!    Generate a 1D array or RND from a given seed.
!    The length of the array is the size of the global CA model,
!    i.e. don't generate values for the halos.
!    Image 1 writes this array to file. Then, inside a critical
!    region each image reads its own data.
!    This is a slow routine.
!  SIDE EFFECTS
!    A tmp external file is created from image 1 to store RND values.
!    It is deleted at the exit from this routine, also by image 1.
!  SOURCE

integer :: co_k, co_j, co_i, k, j, i, funit

integer( kind=iarr ) :: value

integer( kind=iarr ), allocatable :: rnd_arr_int(:)

real, allocatable :: rnd_arr(:)

character(:), allocatable :: fname

fname = "tmp_rnd_file"

! Image 1 sets the file
make_file: if ( this_image() .eq. 1 ) then

  ! Sanity check
  if ( frac1 .lt. 0.0 .or. frac1 .gt. 1.0 ) then
    write (*,*) "ERROR: ca_hx/ca_set_space_rnd: frac1 outside of" //   &
      " admissable range [0..1]: frac1:", frac1
    error stop
  end if

    write (*,*) 
  ! total_cells would have been set already in ca_spalloc
  ! No halo cells are included in this number!
  allocate( rnd_arr(     total_cells ) )
  allocate( rnd_arr_int( total_cells ) )

  call random_seed( put = seed )
  call random_number( rnd_arr )
!  rnd_arr_int = nint( rnd_arr )

  rnd_arr_int = 0
  where ( rnd_arr .lt. frac1 )
    rnd_arr_int = 1
  end where

  open( newunit=funit, file=fname, status="replace", access="stream",  &
        form="unformatted", iostat=ierr )
  if ( ierr .ne. 0 ) then
    write (*,*) "ERROR: ca_hx/ca_set_space_rnd: " //                   &
      "open( fname ), ierr:", ierr
    error stop
  end if

  write (funit) rnd_arr_int

  close(funit, status="keep", iostat=ierr )
  if ( ierr .ne. 0 ) then
    write (*,*) "ERROR: ca_hx/ca_set_space_rnd: " //                   &
      "close( fname ), ierr:", ierr
    error stop
  end if

  deallocate( rnd_arr )
  deallocate( rnd_arr_int )

end if make_file

! All images wait for img1 to write data to file and close it.
sync all

! Each image in turn reads the data from the file
critical

  open( newunit=funit, file=fname, status="old", access="stream",      &
        form="unformatted", iostat=ierr )
  if ( ierr .ne. 0 ) then
    write (*,*) "ERROR: ca_hx/ca_set_space_rnd: " //                   &
      "open( fname ), ierr:", ierr
  end if

  do co_k = 1, ucob(3)
  do    k = 1,  sub(3)
  do co_j = 1, ucob(2)
  do    j = 1,  sub(2)
  do co_i = 1, ucob(1)
  do    i = 1,  sub(1)
    read( funit ) value
    if ( co_k .eq. ci(3) .and. co_j .eq. ci(2) .and.                   &
         co_i .eq. ci(1) ) then
       space(i,j,k) = value
    end if
  end do
  end do
  end do
  end do
  end do
  end do

  close(funit, status="keep", iostat=ierr )
  if ( ierr .ne. 0 ) then
    write (*,*) "ERROR: ca_hx/ca_set_space_rnd: " //                   &
      "close( fname ), ierr:", ierr
  end if

end critical

! All images must read their data before image 1 deletes the file
sync all

! Image 1 deletes the file
file_delete: if ( this_image() .eq. 1 ) then

  open( newunit=funit, file=fname, status="old", iostat=ierr )
  if ( ierr .ne. 0 ) then
    write (*,*) "ERROR: ca_hx/ca_set_space_rnd: " //                   &
      "open( fname ), ierr:", ierr
    error stop
  end if

  close(funit, status="delete", iostat=ierr )
  if ( ierr .ne. 0 ) then
    write (*,*) "ERROR: ca_hx/ca_set_space_rnd: " //                   &
      "close( fname ), ierr:", ierr
    error stop
  end if

end if file_delete

end subroutine ca_set_space_rnd

!*roboend*

end module ca_hx
