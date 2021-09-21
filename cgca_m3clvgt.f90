!*robodoc*m* CGPACK/cgca_m3clvgt
!  NAME
!    cgca_m3clvgt
!  SYNOPSIS

!$Id: cgca_m3clvgt.f90 380 2017-03-22 11:03:09Z mexas $

module cgca_m3clvgt

!  DESCRIPTION
!    Module containing only *thread safe* cleavage routines.
!    In particular, all routines are PURE, designed to be usable from
!    DO CONCURRENT. The module contains variables, submodules and
!    subroutines.
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    cgca_clvgpt (in submodule m3clvgt_sm1).
!  USES
!    cgca_m1co, cgca_m2gb, cgca_m2hx, cgca_m2rot, cgca_m2rnd,
!    cgca_m2glm
!  USED BY
!    cgca
!  SOURCE

use cgca_m1co
use cgca_m2gb,  only: cgca_gcr_pure, cgca_gcf_pure
use cgca_m2glm, only: cgca_ico
use cgca_m2hx,  only: cgca_hxi, cgca_hxg
use cgca_m2rot, only: cgca_csym_pure
use cgca_m2rnd, only: cgca_irs

implicit none

private
public :: cgca_clvgpt

! This array is used to update local GC arrays.
! The components are as follows:
!
!  gcupd(:,1) - grain
!  gcupd(:,2) - neighbour
!  gcupd(:,3) - state, either cgca_gb_state_intact    or
!                             cgca_gb_state_fractured
!
! The idea is that this array is updated every time a grain
! boundary is crossed. Then all local arrays are updated using
! cgca_gcf and this coarray.

integer( kind=iarr ), allocatable :: gcupd(:,:) [:,:,:]

! Counter of a gcupd pair. Set to 1 initially.
integer( kind=idef ) :: gcucnt=1

real( kind=rdef ), parameter ::    &
     zero = 0.0_rdef,              &
      one = 1.0_rdef,              &
    sqrt2 = sqrt(2.0_rdef),        &
 onesqrt2 = one/sqrt2,             &
    sqrt3 = sqrt(3.0_rdef),        &
 onesqrt3 = one/sqrt3,             &

 ! 27 unit vectors connecting the central cell with all its
 ! 26 neighbours + itself, which is a zero vector.
 e(3, -1:1, -1:1, -1:1) =          &
   reshape( (/                     &
 -onesqrt3,-onesqrt3,-onesqrt3,    & ! -1 -1 -1
  zero,    -onesqrt2,-onesqrt2,    & !  0 -1 -1
  onesqrt3,-onesqrt3,-onesqrt3,    & !  1 -1 -1
 
 -onesqrt2, zero,    -onesqrt2,    & ! -1  0 -1 
  zero,     zero,    -one,         & !  0  0 -1
  onesqrt2, zero,    -onesqrt2,    & !  1  0 -1
 
 -onesqrt3, onesqrt3,-onesqrt3,    & ! -1  1 -1
  zero,     onesqrt2,-onesqrt2,    & !  0  1 -1
  onesqrt3, onesqrt3,-onesqrt3,    & !  1  1 -1
 
 -onesqrt2,-onesqrt2, zero,        & ! -1 -1  0
  zero,    -one,      zero,        & !  0 -1  0
  onesqrt2,-onesqrt2, zero,        & !  1 -1  0
 
 -one,      zero,     zero,        & ! -1  0  0
  zero,     zero,     zero,        & !  0  0  0
  one,      zero,     zero,        & !  1  0  0
 
 -onesqrt2, onesqrt2, zero,        & ! -1  1  0
  zero,     one,      zero,        & !  0  1  0
  onesqrt2, onesqrt2, zero,        & !  1  1  0
 
 -onesqrt3,-onesqrt3, onesqrt3,    & ! -1 -1  1
  zero,    -onesqrt2, onesqrt2,    & !  0 -1  1
  onesqrt3,-onesqrt3, onesqrt3,    & !  1 -1  1
 
 -onesqrt2, zero,     onesqrt2,    & ! -1  0  1
  zero,     zero,     one,         & !  0  0  1
  onesqrt2, zero,     onesqrt2,    & !  1  0  1
 
 -onesqrt3, onesqrt3, onesqrt3,    & ! -1  1  1
  zero,     onesqrt2, onesqrt2,    & !  0  1  1
  onesqrt3, onesqrt3, onesqrt3     & !  1  1  1
  /), (/ 3,3,3,3 /) )


! Abstract interface is a fortran 2003 feature.
! This interface is for *thread safe* cleavage "change state" routines.
! The interface for such routines - cgca_clvgsdt and
! cgca_clvgspt, must match it.
! Using this interface, the cleavage "change state"
! routines can be passed as actual arguments to the cleavage
! propagation routines, e.g. cgca_clvgpt, where
! the dummy arguments for these routines are defined by
!  procedure(cgca_clvgst_abstract)

abstract interface
 subroutine cgca_clvgst_abstract( farr, marr, n, cstate, debug,        &
   newstate )
  use cgca_m1co
  integer, parameter :: l=-1, centre=l+1, u=centre+1
  integer( kind=iarr ), intent(in) :: farr(l:u,l:u,l:u),               &
    marr(l:u,l:u,l:u), cstate
  real( kind=rdef ), intent(in) :: n(3)
  logical( kind=ldef ), intent(in) :: debug
  integer( kind=iarr ), intent(out) :: newstate
 end subroutine cgca_clvgst_abstract

 subroutine gcupd_abstract( periodicbc )
  import :: ldef
  logical( kind=ldef ), intent( in ) :: periodicbc
 end subroutine
end interface

! Interfaces for submodule procedures.

interface

  ! In submodule m3clvgt_sm1
  module subroutine cgca_clvgpt( coarray, rt, t, scrit, sub, gcus,    &
      periodicbc, iter, heartbeat, debug )
    ! Inputs:
    ! coarray - cellular array
      integer( kind=iarr ), allocatable, intent(inout) ::              &
        coarray(:,:,:,:)[:,:,:]
    ! rt - rotation tensor coarray
      real( kind=rdef ), allocatable, intent(inout) :: rt(:,:,:)[:,:,:]
    ! t - stress tensor in spatial CS
    ! - scrit - critical values of cleavage stress on 100,
    !   110 and 111 planes
    ! Clearly t and scrit must be in the same units!
      real( kind=rdef ), intent(in) :: t(3,3), scrit(3)
    ! sub - name of the cleavage state calculation routine,
    !   either cgca_clvgsd, or cgca_clvgsp.
      procedure( cgca_clvgst_abstract ) :: sub
    ! gcus - name of the grain connectivity update subroutine, either
    !        cgca_gcupda - all-to-all, or cgca_gcupdn - nearest
    ! neighbour.
      procedure( gcupd_abstract ) :: gcus
    ! periodicbc - if .true. periodic boundary conditions are used,
    !   i.e. global halo exchange is called before every iteration
      logical( kind=ldef ), intent(in) :: periodicbc
    !      iter - number of cleavage iterations, if <=0 then error
    ! heartbeat - if >0 then dump a simple message every
    !             heartbeat iterations
      integer( kind=idef ), intent(in) :: iter, heartbeat
    ! debug - if .true. then will call cgca_dacf with debug
      logical( kind=ldef ), intent(in) :: debug
  end subroutine cgca_clvgpt
end interface

contains

!*roboend*


!*robodoc*s* cgca_m3clvg/gcupd_alloct
!  NAME
!    gcupd_alloct
!  SYNOPSIS

subroutine gcupd_alloct

!  DESCRIPTION
!    This is a private routine, hence the name does not start with
!    cgca_.
!    This routine allocates gcupd array coarray. Therefore it
!    involves implicit sync all.
!    If gcupd is allocated already,
!    it is first deallocated and then reallocated with new
!    codimensions.
!  SIDE EFFECTS
!    gcupd is (re)allocated.
!  USES
!  USED BY
!  SOURCE

integer :: errstat = 0

! Deallocate if already allocated
if ( allocated( gcupd ) ) then
  deallocate( gcupd, stat=errstat )
  if ( errstat .ne. 0 ) then
    write (*,'(a,i0)') "ERROR: cgca_m3clvg/gcupd_alloc:&
     & deallocate( gcupd ), err. code:", errstat
    error stop
  end if
end if

! Allocate and set all values to cgca_gb_state_intact.
allocate( gcupd( cgca_gcupd_size1, cgca_gcupd_size2 )              &
  [ cgca_slcob(1):cgca_sucob(1), cgca_slcob(2):cgca_sucob(2),          &
    cgca_slcob(3):* ], source = cgca_gb_state_intact, stat=errstat )
if ( errstat .ne. 0 ) then
  write (*,'(a,i0)') "ERROR: cgca_m3clvg/gcupd_alloc:&
     & allocate( gcupd ), err. code:", errstat
  error stop
end if

end subroutine gcupd_alloct

!*roboend*


!*robodoc*s* cgca_m3clvg/cgca_clvgn_pure
!  NAME
!    cgca_clvgn_pure
!  SYNOPSIS  

pure subroutine cgca_clvgn_pure( t, grain, r, tcrit, flag, n, cstate,  &
  ierr )

!  INPUTS
!    r - rotation array with extents (number of grains, 3, 3)

real( kind=rdef ), intent(in) :: t(3,3), r(:,:,:), tcrit(3)
integer( kind=iarr), intent(in) :: grain

!  OUTPUTS

logical( kind=ldef ), intent(out) :: flag
real( kind=rdef ), intent(out) :: n(3)
integer( kind=iarr ), intent(out) :: cstate
integer, intent(out) :: ierr

!  DESCRIPTION
!    Given the stress tensor in the spatial CS, (t)
!    and the crystal rotation tensor (r), this routine first
!    calculates whether cleavage happens or not (flag). If cleavage
!    conditions are met (flag=.true.), then the routine calculates the
!    acting cleavage plane normal in spatial coord. system (n) and the
!    type of the cleavage plane (cstate).
!
!    Inputs:
!    - t - stress tensor in spatial CS
!    - r - crystal rotation tensor
!    - tcrit - critical values of cleavage stress on 100, 110
!              and 111 planes
!
!    Outputs:
!    - flag - .true. if cleavage conditions are met, .false. otherwise
!    - n - unit normal vector defining the acting cleavage plane
!    - cstate - cell state, one of the cleavage type.
!    - ierr - error flag. ierr=0 means a successful execution.
!      Any positive value means an error. The value is taken from
!      the output flag of cgca_csym_pure. 
!
!    On output, if flag=.false. then n=0, cstate=cgca_instact_state
!  SIDE EFFECTS
!    None, this is a PURE subroutine.
!  NOTES
!    This is a PURE subroutine, no side effects.
!    It should be used instead of
!    cgca_clvgn where a pure subroutine is required, e.g. from
!    a do concurrent.
!  USES
!    cgca_csym
!  USED BY
!    cgca_clvgsd, cgca_clvgsp
!  NOTES
!    Not accessible from outside of the cgca_clvg module
!  SOURCE

! unit vectors defining 3 cleavage plane families
real( kind=rdef ), parameter ::                & 
 n100(3) = (/ one,      zero,     zero     /), &
 n110(3) = (/ onesqrt2, onesqrt2, zero     /), &
 n111(3) = (/ onesqrt3, onesqrt3, onesqrt3 /)

integer( kind=iarr ) :: i
integer :: iflag

real( kind=rdef ) :: &
 rlocal(3,3), & ! array to store r(grain,:,:), just a convenience
 tc(3,3),     & ! stress tensor in crystal CS
 rsym(3,3),   & ! rotation symmetry tensor
 n100rot(3),  & ! normal to a 100 plane
 n110rot(3),  & ! normal to a 110 plane
 n111rot(3),  & ! normal to a 111 plane
 ttmp(3),     & ! normal stress to a crystallographic plane
 tmax(3),     & ! max normal stress to a crystallographic plane
 n100max(3),  & ! normal to the 100 plane that has max normal stress
 n110max(3),  & ! normal to the 110 plane that has max normal stress
 n111max(3),  & ! normal to the 111 plane that has max normal stress
 p(3),        & ! cleavage stress ratios
 pmax           ! max cleavage stress ratio

! Set rlocal
rlocal = r( grain, :, : )

! Set important values to zero, and initialise the default
! values for outputs
n100max = zero
n110max = zero
n111max = zero
   flag = .false.
      n = zero
 cstate = cgca_intact_state
   tmax = zero
   ierr = 0

! stress tensor in crystal CS
! By our convention, r rotates a vector from the cryst.
! coord. system to spatial. Here we rotate the other
! way round, from spatial to crystal system, hence the
! transpose: ! tc = r^T . t . r
tc = matmul( matmul( transpose( rlocal ), t ), rlocal )

! Find the max normal stresses to three cleavage plane families.
do i = 1, 24

  ! Pick a rotation tensor. If iflag .ne. 0, then an error
  ! condition has occurred. Then abandon the computation, set
  ! flag accordingly and return immediately.
  call cgca_csym_pure( i, rsym, iflag )
  if ( iflag .ne. 0 ) then
    ierr = iflag
    return
  end if

  ! normals to partucular {100}, {110}, {111} planes
  n100rot = matmul( rsym, n100 )
  n110rot = matmul( rsym, n110 )
  n111rot = matmul( rsym, n111 )

  ! Projections of the stress tensor to normals, i.e.
  ! the normal stresses to crystallographic planes. 
  ! (tc . n) . n
  ttmp(1) = dot_product( n100rot, matmul( tc, n100rot) )
  ttmp(2) = dot_product( n110rot, matmul( tc, n110rot) )
  ttmp(3) = dot_product( n111rot, matmul( tc, n111rot) )

  ! Choose the plane of each type, that has the max normal stress
  ! Signs are taken into account. Since tmax is zero initially,
  ! negative stresses are less than tmax, and thus do not cause
  ! cleavage.

  if ( ttmp(1) .gt. tmax(1) ) then   ! 100 plane
    tmax(1) = ttmp(1)                ! stress
    n100max = n100rot                ! normal to a particular 100 plane
  end if
  if ( ttmp(2) .gt. tmax(2) ) then   ! 110 plane
    tmax(2) = ttmp(2)                ! stress
    n110max = n110rot                ! normal to a particular 110 plane
  end if
  if ( ttmp(3) .gt. tmax(3) ) then   ! 111 plane
    tmax(3) = ttmp(3)                ! stress
    n111max = n111rot                ! normal to a particular 111 plane
  end if

end do

! calculate the cleavage factors (ratios)
p = tmax / tcrit
pmax = maxval(p)

! check for cleavage
if ( pmax .ge. one ) then

  ! cleavage is happening, on 100 by default
    flag = .true.
       n = n100max
  cstate = cgca_clvg_state_100_edge 

  if ( p(2) .gt. p(1) ) then
    ! cleavage on 110
         n = n110max
    cstate = cgca_clvg_state_110_edge 
  end if
  
  if ( p(3) .gt. p(1) .and. p(3) .gt. p(2) ) then
    ! cleavage on 111
         n = n111max
    cstate = cgca_clvg_state_111_edge 
  end if

  ! rotate n back into the spatial coord. system
  n = matmul( rlocal, n ) 

end if

end subroutine cgca_clvgn_pure

!*roboend*


!*robodoc*s* cgca_m3clvg/cgca_clvgsdt
!  NAME
!    cgca_clvgsdt
!  SYNOPSIS

pure subroutine cgca_clvgsdt( farr, marr, n, cstate, debug, newstate )

!  PARAMETERS

integer, parameter :: l=-1, centre=l+1, u=centre+1
real( kind=rdef ), parameter :: trshld = 0.17325932611400130485_rdef

!  INPUTS

integer( kind=iarr ), intent(in) :: farr(l:u,l:u,l:u), &
 marr(l:u,l:u,l:u), cstate
real(    kind=rdef ), intent(in) :: n(3)
logical( kind=ldef ), intent(in) :: debug

!  OUTPUT

integer( kind=iarr ), intent(out) :: newstate

!  USES
!    cgca_gcr, cgca_gcf
!  USED BY
!    cgca_clvgp
!  DESCRIPTION
!    This routine determines the cleavage (CLVG) state (S) of the
!    central cell. This is a deterministic (D) routine.
!    Hence the name is CLVGSD.
!    If there is a cleaved neighbour, such that the vector connecting
!    it to the central cell is on or close to the cleavage plane,
!    then the central cell state is changed to the given cleavage state.
!
!    If the cleaved neighbour belongs to another grain, the analysis
!    takes into account the grain boundary state. If the grain
!    boundary is marked as intact, then the crack can cross it,
!    and the central cell takes the values of the cleaved neighbour.
!    If the grain boundary is marked as fractured, then the crack
!    cannot cross it. This means that even if there is a cleaved
!    neighbour on the cleavage plane, the central cell will still not
!    change state.
!
!    If the grain boudary is intact, and crack crosses it,
!    mark this GB as fractured immediately in the local CG array
!    and add this GB fracture to gcupd on this image.
!
!    Inputs:
!    - farr - (3,3,3) array of failed states, cell state type
!             cgca_state_type_frac
!    - marr - (3,3,3) array of material states, cell state type
!             cgca_state_type_grain
!    - n - vector defining the cleavage plane
!    - cstate - cleavage state 
!    - debug - if .true. will dump some debug info
!
!    Outputs:
!    - newstate - updated central cell state
!
!  NOTES
!    The threshold analysis is explained in
!    A. Shterenlikht, L. Margetts, Three-dimensional cellular automata
!    modelling of cleavage propagation across crystal boundaries in
!    polycrystalline microstructures, Proc. Roy. Soc. A 471:20150039,
!    http://dx.doi.org/10.1098/rspa.2015.0039 .
!    There are always at least 2 neighbouring cells for which
!    dot_product(e,n) is less than ~0.1732. However, 2 cells
!    is not enough! If threshold is left at that, 2 cells
!    lead to linear 1D cracks, which are not reasonable.
!  SOURCE

integer( kind=idef ) :: &
 x1, x2, x3,            & ! coordinates within (3,3,3) neighbourhood
 img,                   & ! this_image()
 nimgs                    ! num_images()

integer :: iflag

logical :: intact

  img = this_image()
nimgs = num_images()

! initially the new central cell state is the same as old
newstate = farr( centre, centre, centre )

! scan all neighbourhood cells
outer: do x3 = l, u
       do x2 = l, u
       do x1 = l, u

  ! Check whether the neighbour has cleaved. The central cell is never
  ! cleaved, so the check will always fail for the central cell.
  clv: if ( any( cgca_clvg_states .eq. farr(x1,x2,x3)) ) then

    ! If the neighbour is from the same grain, then just check the
    ! orientation of the cleavage plane
    same: if ( marr(centre,centre,centre) .eq. marr(x1,x2,x3) ) then

      ! debug
      !  if (debug) then
      !   call cgca_gcr(1,2, intact)
      !   if (.not. intact) write (*,"(2(a,i0),a,3(tr1,f7.4),a,f7.4)") &
      !    "zzz: grain=",marr(centre,centre,centre), &
      !    ", cstate=", cstate, &
      !    ", n=(", n, ") ", &
      !    abs(dot_product(e(:,x1,x2,x3),n))
      !  end if

      ! If the dot product is smaller than the threshold, i.e. if
      ! the vector connecting the central cell with the neighbour
      ! is in the cleavage plane, then the central cell has cleaved
      if ( abs( dot_product( e(:,x1,x2,x3), n ) ) .lt. trshld ) then
        newstate = cstate  ! change its state
        exit outer         ! and exit the main loop
      end if

    else

      ! The neighbour is from a different grain. In this case
      ! the boundary must be intact in addition to the cleavage
      ! plane orientation criterion.
      intact = .false.

      ! Get the recorded state of the boundary between the two grains.
      ! cgca_gcr_pure returns .TRUE. if the boundary is intact
      ! and .FALSE. if the boundary is fractured. If iflag .ne. 0,
      ! then some error condition occurred, potentially fatal.
      ! However, to keep this routine PURE, we cannot abort here,
      ! so not sure what I can do with iflag, perhaps combine it
      ! somehow with the success flag of this subroutine?
      call cgca_gcr_pure( marr(centre,centre,centre), marr(x1,x2,x3),  &
                          intact, iflag )

      ! same as in the previous case, but with the added constraint
      ! that the grain boundary must be intact
      crossgb: if ( intact .and.                                       &
           abs( dot_product( e(:,x1,x2,x3), n ) ) .lt. trshld ) then

        ! change its state
        newstate = cstate

        ! Mark GB as fractured on this image
        ! The return diagnostic iflag is ignored for now.
        call cgca_gcf_pure( marr(centre, centre, centre),              &
          marr(x1, x2, x3), iflag )

        ! Add the (grain, neighbour) pair to gcupd coarray
        ! on this image
        gcupd( gcucnt, : ) = (/ marr( centre , centre , centre ) , &
                                marr(   x1   ,   x2   ,   x3   ) ,     &
                                cgca_gb_state_fractured /)

        ! debug output, not allowed in a pure procedure
!        if (debug)                                                     &
!          write (*,"(4(a,i0),2(tr1,i0),')',a,27(i0,tr1),').')")        &
!            "DEBUG: cgca_clvgsd: img: ", img, ": newstate=", newstate, & 
!            ", gcucnt=", gcucnt,                                       &
!            ", calling cgca_gcf, gcupd=(", gcupd( gcucnt , : ),&
!            ", marr=(", marr

        ! increment the gcupd pair counter
        gcucnt = gcucnt + 1

        ! IO not allowed in a pure procedure, so just terminate without
        ! issuing the error.
        ! Issue fatal error if the length of the gcupd has been
        ! exceeded.
        if ( gcucnt .gt. cgca_gcupd_size1 ) then
!          write( *, '(a,i0)' ) "ERROR: cgca_m3clvg/cgca_clvgsd:&
!            & gcucnt .gt. cgca_gcupd_size1, image: ", img
          error stop
        end if

        ! now exit the main loop
        exit outer

      end if crossgb

    end if same

 end if clv

 ! now check another neighbouring cell, i.e. increment the loop counter

end do
end do
end do outer

end subroutine cgca_clvgsdt

!*roboend*


!*robodoc*s* cgca_m3clvg/cgca_clvgspt
!  NAME
!    cgca_clvgspt
!  SYNOPSIS

pure subroutine cgca_clvgspt( farr, marr, n, cstate, debug, newstate )

!  PARAMETERS

integer,parameter :: l=-1, centre=l+1, u=centre+1
real(kind=rdef),parameter ::            &
 ltrshld = 0.17325932611400130485_rdef, & ! see cgca_clvgpd
 utrshld = 0.27_rdef, interval = utrshld-ltrshld

!  INPUTS

integer(kind=iarr),intent(in) :: farr(l:u,l:u,l:u), &
 marr(l:u,l:u,l:u), cstate
real(kind=rdef),intent(in) :: n(3)
logical(kind=ldef),intent(in) :: debug

!  OUTPUT

integer(kind=iarr),intent(out) :: newstate

!  DESCRIPTION
!    This routine determines the cleavage (CLVG) state (S) of the central
!    cell. This is a probabilistic (P) routine. Hence the name is CLVGSP.
!    If there is a cleaved neighbour, such that the vector connecting it to
!    the central cell is on or close to the cleavage plane, then the
!    central cell has a probability of changing state to the given cleavage
!    state.
!
!    If the cleaved neighbour belongs to another grain, the analysis takes
!    into account the grain boundary state. If the grain boundary is
!    marked as intact, then the crack can cross it, and the central cell
!    takes the values of the cleaved neighbour. If the grain boundary
!    is marked as fractured, then the crack cannot cross it. This means
!    that even if there is a cleaved neighbour on the cleavage plane,
!    the central cell will still not change state.
!
!    Inputs:
!    - farr - array of failed states, cell state type cgca_state_type_frac
!    - marr - array of material states, cel state type cgca_state_type_grain
!    - n - vector defining the cleavage plane
!    - cstate - cleavage state 
!    - debug - if .true. will dump some debug info
!
!    Outputs:
!    - newstate - updated central cell state
!  NOTES
!    This routine has two thresholds, the upper and the lower.
!  USES
!    cgca_gcr, cgca_gcf
!  USED BY
!  SOURCE

integer(kind=idef) :: x1,x2,x3
real(kind=rdef) :: rnd, proj, prob

logical :: intact

! initially the new central cell state is the same as old
newstate = farr(centre,centre,centre)

! scan all neighbourhood cells
outer: do x3=l,u
do x2=l,u
do x1=l,u

 ! Check whether the neighbour has cleaved. The central cell is never
 ! cleaved, so the check will always fail for the central cell.
 clv: if ( any( cgca_clvg_states .eq. farr(x1,x2,x3)) ) then

  proj = abs(dot_product(e(:,x1,x2,x3),n))

  same: if ( marr(centre,centre,centre) .eq. marr(x1,x2,x3) ) then
   ! The neighbour is from the same grain.
   ! If the central cell is close to the cleavage plane,
   ! change its state. This is a deterministic check for the
   ! lower threshold. If the central cell is further from the
   ! cleavage plane, but not too far, it has some probability
   ! to cleave. This is a probabilistic check for upper threshold.
   z1: if ( proj .lt. ltrshld ) then
    newstate = cstate
    exit outer
   else if ( proj .lt. utrshld ) then
    ! The power must be .ge. 1. If the power is 1,
    ! then the probability is a linear function of proj.
    ! If the power is > 1 then the probability is a power function.
    ! The higher the power, the steeper the descent.
    ! In other words, the higher the power, the lower the chances
    ! of cleavage for proj values greater than the lower threshold. 
    prob = ((utrshld-proj) / interval)**1
    call random_number(rnd)
    if ( prob .gt. rnd) newstate = cstate
    exit outer
   end if z1

  else

   ! The neighbour is from another grain. As above, but the
   ! additional condition is that the grain boundary must be intact.
   intact = .false.
   call cgca_gcr(marr(centre,centre,centre), marr(x1,x2,x3), intact)
   inta: if ( intact ) then
    z2: if ( proj .lt. ltrshld ) then
     newstate = cstate

     ! Mark GB as fractured straight away
     call cgca_gcf( marr(centre,centre,centre), marr(x1,x2,x3) )

     ! debug output
     if (debug) write (*,"(2(a,i0),2(tr1,i0),').')")                   &
       "DEBUG: cgca_clvgsd: image ", this_image(),                     &
       ": called _gcf, set gcupd=(", gcupd

     exit outer
    else if ( proj .lt. utrshld ) then
     prob = ((utrshld-proj) / interval)**1
     call random_number(rnd)
     if ( prob .gt. rnd) then
      newstate = cstate
      ! Mark GB as fractured straight away
      call cgca_gcf( marr(centre,centre,centre), marr(x1,x2,x3) )
      ! debug output
      if (debug) write (*,"(2(a,i0),2(tr1,i0),').')")                  &
       "DEBUG: cgca_clvgsd: image ", this_image(),                     &
       ": called _gcf, set gcupd=(", gcupd
      exit outer
     end if
    end if z2
   end if inta

  end if same
 end if clv
end do
end do
end do outer

end subroutine cgca_clvgspt

!*roboend*


!*robodoc*s* cgca_m3clvg/cgca_dacft
!  NAME
!    cgca_dacft
!  SYNOPSIS

subroutine cgca_dacft( coarray , debug )

!  INPUTS

integer( kind=iarr ), allocatable,intent(inout) ::                     &
 coarray(:,:,:,:)[:,:,:]
logical( kind=ldef ), intent(in) :: debug

!  SIDE EFFECTS
!    Changed state of coarray
!  DESCRIPTION
!    This routine DeActivates Crack Flanks, hence the name DACF.
!    The idea is that we must distinguish crack edge cells, which
!    can attract new cracked (cleaved) cells, and crack flanks, which
!    are inactive, i.e. cells representing crack flanks have very low
!    SIF (stress intensity factor), and maybe low stresses too.
!    Crack flank cells cannot attract new cleaved cells.
!    The cell states are defined in module cgca_m1co.
!
!    The distinction is made based on the number of cleaved
!    neighbours. If there are too many cleaved neighbours, then
!    the central cell is a crack flank.
!
!    The intention is that this routine is called after every cleavage
!    propagation increment, to prevent cracks becoming large 3D bodies.
!
!    This routines runs only once over the coarray. So we don't put the
!    sync here. But a sync all probably should be used before and after
!    a call to this routine. Also, the halo exchange probably should
!    be done after running this routine.
!
!    If debug=.true. then will dump *lots* of debug output
!  USES
!  USED BY
!    cgca_clvgp
!  SOURCE

integer, parameter ::                                                  &
 lclvg_states = lbound( cgca_clvg_states , dim=1 ),                    &
 uclvg_states = ubound( cgca_clvg_states , dim=1 ),                    &
 l=-1, centre=l+1, u=centre+1

integer( kind=idef ) :: &
 lbv(4) ,& ! lower bounds of the complete (plus virtual) coarray
 ubv(4) ,& ! upper bounds of the complete (plus virtual) coarray
 lbr(4) ,& ! lower bounds of the "real" coarray, lower virtual+1
 ubr(4) ,& ! upper bounds of the "real" coarray, upper virtual-1
 x1     ,& ! local coordinates in an array, which are also
 x2     ,& ! do loop counters
 x3
integer( kind=iarr ), allocatable, save :: array(:,:,:)
integer( kind=iarr ) :: neiarr(l:u,l:u,l:u)
integer :: thisimage, errstat=0, ncount, i
logical( kind=ldef ) :: clvnei(l:u,l:u,l:u), &
 samegrain(l:u,l:u,l:u), csg(l:u,l:u,l:u)

! get image number
thisimage = this_image()

! no sanity checks for speed

! determine the extents
lbv = lbound(coarray)
ubv = ubound(coarray)
lbr = lbv+1
ubr = ubv-1

! allocate the temp array
if (.not. allocated(array)) allocate( &
  array(lbv(1):ubv(1), lbv(2):ubv(2), lbv(3):ubv(3)), stat=errstat)
if (errstat.ne.0) then
  write (*,'(a,i0)') "ERROR: cgca_dacf: image ",thisimage
  write (*,'(a)') "ERROR: cgca_dacf: cannot allocate array"
  error stop
end if

! copy coarray state type 2, fracture states, to temp array
array = coarray(:,:,:,cgca_state_type_frac)

! loop over all cells
do x3 = lbr(3), ubr(3)
do x2 = lbr(2), ubr(2)
do x1 = lbr(1), ubr(1)

 ! analyse only crack edge cells
 cleav: if ( any (coarray(x1,x2,x3,cgca_state_type_frac) .eq.        &
          cgca_clvg_states_edge )) then

   ! count the number of cleaved neighbours of the same grain
   neiarr = coarray( x1-1:x1+1, x2-1:x2+1, x3-1:x3+1,                &
                     cgca_state_type_frac )

   ! Do not count the central cell, so set its state to some
   ! value that is not in the cgca_clvg_states. I use the intact
   ! state here.
   neiarr(centre,centre,centre) = cgca_intact_state

   ! logical array: .true. if the values are identical,
   ! .false. otherwise. samegrain is a (3,3,3) neighbourhood
   ! array around the central cell in question.
   samegrain = coarray( x1, x2, x3, cgca_state_type_grain ) .eq.     &
               coarray( x1-1:x1+1, x2-1:x2+1, x3-1:x3+1,             &
                          cgca_state_type_grain )

   ! debug
   if (debug) write (*,'(a,i0,a,27(tr1,i0),a,27(tr1,l1))')           &
    "DEBUG: cgca_dacf: image ", thisimage, ": neiarr=", neiarr,      &
    ", samegrain=", samegrain

   ncount = 0

   ! for all cleavage states
   do i = lclvg_states, uclvg_states
     clvnei = cgca_clvg_states(i) .eq. neiarr
        csg = clvnei .and. samegrain
     ncount = ncount + count(csg)

     ! debug
     if (debug) write (*,'(4(a,i0),2(a,27(tr1,l1)),a,i0)')           &
       "DEBUG: cgca_dacf: image ", thisimage, ": i=", i,             &
       ", cgca_clvg_states(",i,")=", cgca_clvg_states(i),            &
       ", clvnei=", clvnei, ", csg=", csg, ", ncount=", ncount

   end do

   ! if a cell has 6 cleaved neighbours or more,
   ! then mark it as crack flank. 
   nei: if ( ncount .ge. 6) then
   
     ! debug
     if (debug) write (*,'(3(a,i0),2(tr1,i0),a)')                    &
       "DEBUG: cgca_dacf: image ", thisimage,                        &
       ": grain=", coarray(x1,x2,x3,cgca_state_type_grain),          &
       ", crack front cell x1,x2,x3=", x1, x2, x3, " deactivated."

     ! change the state preserving the cleavage plane family
     ! Note: we are reading from "coarray" but writing into
     ! the temp "array"
     if ( coarray( x1, x2 ,x3, cgca_state_type_frac ) .eq.           &
                      cgca_clvg_state_100_edge ) then
       array(x1,x2,x3) = cgca_clvg_state_100_flank
     else if ( coarray( x1, x2, x3, cgca_state_type_frac ) .eq.      &
                      cgca_clvg_state_110_edge ) then
       array(x1,x2,x3) = cgca_clvg_state_110_flank
     else 
       array(x1,x2,x3) = cgca_clvg_state_111_flank
     end if

   end if nei
 end if cleav
end do
end do
end do

! write array to coarray
coarray(:,:,:,cgca_state_type_frac) = array

! do not deallocate array. Let it exist until the program
! terminates.

end subroutine cgca_dacft

!*roboend*


!*robodoc*s* cgca_m3clvg/cgca_dacf1t
!  NAME
!    cgca_dacf1t
!  SYNOPSIS

subroutine cgca_dacf1t(coarray,debug)

!  INPUTS

integer( kind=iarr ), allocatable, intent(inout) ::                    &
 coarray(:,:,:,:)[:,:,:]
logical( kind=ldef ), intent(in) :: debug

!  SIDE EFFECTS
!    Changed state of coarray
!  DESCRIPTION
!    Same as cgca_dacf, but no attention is paid to which
!    grain a cell belongs to. So when we count the number
!    of fractured neighbours, these can be from any grain.
!    What this means is that, when a crack propagates
!    from one grain to another, and the crack planes
!    IN THE MODEL coinside, this routine will deactivate
!    the fractured cells on the interface. This helps the
!    grain boundary fracture analysis. Since there are no
!    crack fronts at the boundary, there is no bondary to
!    fracture! In other words, if a crack really slices
!    from one grain to another with no deviation, then
!    the grain boundary fracture does not happen, the
!    two grain system in question is already fully separated
!    into two different bodies.
!  USES
!  USED BY
!    cgca_clvgp
!  SOURCE

integer,parameter :: lclvg_states = lbound(cgca_clvg_states,dim=1), &
 uclvg_states = ubound(cgca_clvg_states,dim=1),                     &
 l=-1, centre=l+1, u=centre+1

integer(kind=idef) :: &
  lbv(4) ,& ! lower bounds of the complete (plus virtual) coarray
  ubv(4) ,& ! upper bounds of the complete (plus virtual) coarray
  lbr(4) ,& ! lower bounds of the "real" coarray, lower virtual+1
  ubr(4) ,& ! upper bounds of the "real" coarray, upper virtual-1
  x1     ,& ! local coordinates in an array, which are also
  x2     ,& ! do loop counters
  x3        !
integer(kind=iarr),allocatable,save :: array(:,:,:)
integer(kind=iarr) :: neiarr(l:u,l:u,l:u)
integer :: thisimage, errstat=0, ncount, i
logical(kind=ldef) :: clvnei(l:u,l:u,l:u)

! get image number
thisimage=this_image()

! no sanity checks for speed

! determine the extents
lbv = lbound(coarray)
ubv = ubound(coarray)
lbr = lbv + 1
ubr = ubv - 1

! allocate the temp array

if ( .not. allocated( array ) )                                                &
 allocate( array( lbv(1):ubv(1), lbv(2):ubv(2), lbv(3):ubv(3) ), stat=errstat )
if (errstat.ne.0) then
  write (*,'(2(a,i0))') "****** ERROR: cgca_dacf1: image ",thisimage,          &
     ": cannot allocate array, error code ", errstat
  error stop
end if

! copy coarray state type 2, fracture states, to temp array
array = coarray(:,:,:,cgca_state_type_frac)

! loop over all cells
do x3 = lbr(3),ubr(3)
do x2 = lbr(2),ubr(2)
do x1 = lbr(1),ubr(1)

 ! analyse only crack edge cells
 cleav: if ( any (coarray(x1,x2,x3,cgca_state_type_frac) .eq.      &
          cgca_clvg_states_edge )) then

  ! Count the number of cleaved neighbours. Do not count the central
  ! cell, so set its state to some value that is not in the
  ! cgca_clvg_states. I use the intact state here.
  neiarr = coarray(x1-1:x1+1, x2-1:x2+1, x3-1:x3+1, cgca_state_type_frac)
  neiarr(centre,centre,centre) = cgca_intact_state

  ncount = 0
  do i=lclvg_states,uclvg_states
   clvnei = cgca_clvg_states(i) .eq. neiarr
   ncount = ncount + count(clvnei)

   ! debug
   if (debug) write (*,'(4(a,i0),a,27(tr1,l1),a,i0)')              &
    "DEBUG: cgca_dacf: image ", thisimage, ": i=", i,              &
    ", cgca_clvg_states(",i,")=", cgca_clvg_states(i),             &
    ", clvnei=", clvnei, ", ncount=", ncount

  end do

  ! if a cell has 5 cleaved neighbours or more,
  ! then mark it as crack flank. 
  nei: if ( ncount .ge. 5) then
   
   ! debug
   if (debug) write (*,'(3(a,i0),2(tr1,i0),a)')                    &
    "DEBUG: cgca_dacf1: image ", thisimage,                        &
    ": grain=", coarray(x1,x2,x3,cgca_state_type_grain),           &
    ", crack front cell x1,x2,x3=", x1, x2, x3, " deactivated."

   ! change the state preserving the cleavage plane family
   ! Note: we are reading from "coarray" but writing into the temp
   ! "array"
   if ( coarray(x1,x2,x3,cgca_state_type_frac) .eq.                &
         cgca_clvg_state_100_edge ) then
     array(x1,x2,x3) = cgca_clvg_state_100_flank
   else if ( coarray(x1,x2,x3,cgca_state_type_frac) .eq.           &
          cgca_clvg_state_110_edge ) then
     array(x1,x2,x3) = cgca_clvg_state_110_flank
   else 
     array(x1,x2,x3) = cgca_clvg_state_111_flank
   end if

  end if nei
 end if cleav
end do
end do
end do

! write array to coarray
coarray(:,:,:,cgca_state_type_frac) = array

! do not deallocate array. Let it exist until the program
! terminates.

end subroutine cgca_dacf1t

!*roboend*

end module cgca_m3clvgt
