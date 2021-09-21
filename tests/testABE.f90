!*robodoc*u* tests/testABE
!  NAME
!    testABE
!  SYNOPSIS

!$Id: testABE.f90 380 2017-03-22 11:03:09Z mexas $

program testABE

!  PURPOSE
!    Checking: cgca_gbs
!  DESCRIPTION
!    Checking grain boundary smoothing routine
!  NOTES
!    The program must be called with 2 command line arguments,
!    both positive integers. These are codimensions along 1 and 2.
!    The number of images must be such that
!    codimension3 = num_images()/( codimension1 * codimension3 )
!    is a positive integer. Example:
!      cafrun -np 16 ./testABE.x 2 2      ! OpenCoarrays
!    or
!      ./testABE.x 2 2                    ! Intel, Cray
!    which will make the third codimension equal to 16/(2*2)=4.
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  USES
!    cgca testaux
!  USED BY
!    Part of CGPACK test suite
!  SOURCE

use testaux

implicit none

logical(kind=ldef),parameter :: yesdebug=.true., nodebug=.false., &
  periodicbc=.true.,  noperiodicbc=.false.

real,parameter :: gigabyte=real(2**30), resolution=1.0e-5, &
! cleavage stress on 100, 110, 111 planes for BCC,
! see the manual for derivation.
  scrit(3) = (/ 1.05e4, 1.25e4, 4.90e4 /)

 real(kind=rdef),allocatable :: grt(:,:,:)[:,:,:]

integer(kind=idef) :: l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3, &
 cou3, nuc,nimages,codim(3)[*], i
integer(kind=iarr),allocatable :: space(:,:,:,:)[:,:,:]
integer(kind=ilrg) :: icells,mcells

real :: image_storage

logical(kind=ldef) :: solid, image1

!*********************************************************************72
! first executable statement

nimages = num_images()
 image1 = .false.
if (this_image() .eq. 1) image1 = .true.

! do a check on image 1
if (image1) then
 call getcodim(nimages,codim)
 ! print a banner
 call banner("ABE")
 ! print the parameter values
 call cgca_pdmp
 write (*,'(a,i0,a)') "Running on ", nimages, " images in a 3D grid."
end if

sync all

codim(:) = codim(:)[1]

l1 = 1
l2 = l1
l3 = l1

u1 = 10
u2 = u1
u3 = u1

col1 = 1
cou1 = codim(1)-col1+1
col2 = 1
cou2 = codim(2)-col2+1
col3 = 1
cou3 = codim(3)-col3+1

! total number of cells in a coarray
icells = int(u1-l1+1,kind=ilrg) * int(u2-l2+1,kind=ilrg) * &
  int(u3-l3+1,kind=ilrg)

! total number of cells in the model
mcells = icells * int(codim(1),kind=ilrg) * int(codim(2),kind=ilrg) * &
  int(codim(3),kind=ilrg)

! total number of nuclei
!nuc = resolution*mcells
nuc = 2

if (image1) then
  write (*,'(a,2(i0,":",i0,","),i0,":",i0,")")') &
    "bounds: (",l1,u1,l2,u2,l3,u3
  write (*,'(a,2(i0,":",i0,","),i0,":",i0,")")') &
    "cobounds: (",col1,cou1,col2,cou2,col3,cou3

  ! An absolute minimum of storage, in GB, per image.
  ! A factor of 2 is used because will call _sld, which
  ! allocates another array of the same size and kind as
  ! coarray.
  image_storage = real(2 * icells*storage_size(space)/8)/gigabyte

  write (*,'(a,i0,a)') "Each image has ",icells, " cells"
  write (*,'(a,i0,a)') "The model has ", mcells, " cells"
  write (*,'(a,i0,a)') "The model has ", nuc, " nuclei"
  write (*,'(a,es9.2,a)') "Each image will use at least ", &
    image_storage, " GB memory"
end if

! initialise random number seed
call cgca_irs(nodebug)

! allocate space with two layers
call cgca_as(l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3,2,space)

! allocate rotation tensors
call cgca_art(1,nuc,col1,cou1,col2,cou2,col3,grt)

! initialise space
space(:,:,:,cgca_state_type_grain) = cgca_liquid_state
space(:,:,:,cgca_state_type_frac ) = cgca_intact_state

! nuclei, sync all inside
call cgca_nr(space,nuc,yesdebug)

! solidify, implicit sync all inside
call cgca_sld(space,noperiodicbc,0,10,solid)

! dump space grain layer to file
call cgca_swci(space,cgca_state_type_grain,10,"zg0.raw")
sync all

! smooth the grain boundary, many iterations
do i=1,16
call cgca_gbs(space)
sync all

if (image1) write (*,*) "done GB smoothing, iter=",i

if (i .eq. 1) call cgca_swci(space,cgca_state_type_grain,10,"zg1.raw")
if (i .eq. 2) call cgca_swci(space,cgca_state_type_grain,10,"zg2.raw")
if (i .eq. 4) call cgca_swci(space,cgca_state_type_grain,10,"zg4.raw")
if (i .eq. 8) call cgca_swci(space,cgca_state_type_grain,10,"zg8.raw")
! no sync here, because the halo exchange changes
! only *local* coarray values.

! internal halo exchange
call cgca_hxi(space)
sync all
end do

! dump space grain layer to file
! all others wait at the barrier, hence sync needed
call cgca_swci(space,cgca_state_type_grain,10,"zg16.raw")

! However, since there's nothing more to do, no sync is needed.

! deallocate all arrays, implicit sync all.
call cgca_ds(space)

end program testABE

!*roboend*
