!*robodoc*u* tests/testAAV
!  NAME
!    testAAV
!  SYNOPSIS

!$Id: testAAV.f90 380 2017-03-22 11:03:09Z mexas $

program testAAV

!  PURPOSE
!    Checking: cgca_clvgp_nocosum (no CO_SUM version), cgca_clvgsd
!  DESCRIPTION
!*********************************************************************72
!    Checking deterministic cleavage propagation with 1 cleavage
!    nucleus in the middle of image 1.  Single grain with random
!    orientation.  Several cleavage increments.
!    Check how crack propagates across image boundary.
!  NOTES
!    The program must be called with 2 command line arguments,
!    both positive integers. These are codimensions along 1 and 2.
!    The number of images must be such that
!    codimension3 = num_images()/( codimension1 * codimension3 )
!    is a positive integer. Example:
!      cafrun -np 16 ./testAAV.x 2 2      ! OpenCoarrays
!    or
!      ./testAAV.x 2 2                    ! Intel, Cray
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
   periodicbc=.true.
 real,parameter :: gigabyte=real(2**30), resolution=1.0e-5, &

! cleavage stress on 100, 110, 111 planes for BCC,
! see the manual for derivation.
  scrit(3) = (/ 1.05e4, 1.25e4, 4.90e4 /)

 real(kind=rdef),allocatable :: grt(:,:,:)[:,:,:]
 real(kind=rdef) :: t(3,3)  ! stress tensor

 integer(kind=idef) :: l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3, &
   cou3,   &
   nuc,    & ! number of nuclei in the model
   nimages,codim(3)[*]
 integer(kind=iarr),allocatable :: space(:,:,:,:)[:,:,:]
 integer(kind=ilrg),allocatable :: grainvol(:)[:,:,:],fracvol(:)[:,:,:]
 integer(kind=ilrg) :: icells,mcells

 logical(kind=ldef) :: image1

 real :: image_storage

!*********************************************************************72
! first executable statement

nimages=num_images()
image1=.false.
if (this_image().eq.1) image1=.true.

! do a check on image 1
if (image1) then
 call getcodim(nimages,codim)
 ! print a banner
 call banner("AAV")
 ! print the parameter values
 call cgca_pdmp
 write (*,'(a,i0,a)') "running on ", nimages, " images in a 3D grid"
 write (*,*) "codim:", codim
end if

sync all

codim(:) = codim(:)[1]

l1=1
l2=l1
l3=l1

! The array size is only controlled by this value
! in this program.
u1=50
u2=u1
u3=u1

col1=1
cou1=codim(1)-col1+1
col2=1
cou2=codim(2)-col2+1
col3=1
cou3=codim(3)-col3+1

! total number of cells in a coarray
icells = int(u1-l1+1,kind=ilrg) * int(u2-l2+1,kind=ilrg) * &
  int(u3-l3+1,kind=ilrg)

! total number of cells in the model
mcells = icells * int(codim(1),kind=ilrg) * int(codim(2),kind=ilrg) * &
  int(codim(3),kind=ilrg)

! total number of nuclei
!nuc = resolution*mcells
! just 1 nuclei in this test
nuc=1

if (image1) then
  write (*,'(a,2(i0,":",i0,","),i0,":",i0,")")') &
    "bounds: (",l1,u1,l2,u2,l3,u3
  write (*,'(a,2(i0,":",i0,","),i0,":",i0,")")') &
    "cobounds: (",col1,cou1,col2,cou2,col3,cou3

  ! An absolute minimum of storage, in GB, per image.
  image_storage = real( &
    ! 2nd space array is allocated in _sld
    2 * icells*storage_size(space) &
    ! grain volumes
    + nuc*storage_size(grainvol) &
    ! rotation tensors
    + nuc*9*storage_size(grt))/8/gigabyte 

  write (*,'(a,i0,a)') "Each image has ",icells, " cells"
  write (*,'(a,i0,a)') "The model has ", mcells, " cells"
  write (*,'(a,i0,a)') "The model has ", nuc, " nuclei"
  write (*,'(a,es9.2,a)') "Each image will use at least ", &
    image_storage, " GB memory"
end if

! initialise random number seed
call cgca_irs(nodebug)

! allocate coarray
call cgca_as(l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3,2,space)

! initialise coarrays: a single grain and no damage
space(:,:,:,cgca_state_type_grain) = 1
space(:,:,:,cgca_state_type_frac) = cgca_intact_state

! set a single crack nucleus at the upper end of image1 array
space(u1,u2,u3,cgca_state_type_frac)[col1,col2,col3] = &
 cgca_clvg_state_100_edge

! Allocate grain and fracture volumes
call cgca_av(1,nuc,col1,cou1,col2,cou2,col3,grainvol)
call cgca_av(1,nuc,col1,cou1,col2,cou2,col3,fracvol)

! allocate rotation tensors
call cgca_art(1,nuc,col1,cou1,col2,cou2,col3,grt)

! set the stress tensor
t = 0.0
t(1,1) = 1.0e6
t(2,2) = -1.0e6

! assign rotation tensors, sync all inside
call cgca_rt(grt)


! propagate cleavage, sync inside
! subroutine cgca_clvgp_nocosum( coarray, rt, t, scrit, sub,
!    gcus, periodicbc, iter, heartbeat, debug )
call cgca_clvgp_nocosum( space, grt, t, scrit, cgca_clvgsd,            &
     cgca_gcupdn, .false., 50, 10, nodebug )

! dump the array
call cgca_swci(space,cgca_state_type_frac,10,"zf.raw")
sync all

! calculate volumes, sync all inside
call cgca_gv(space,grainvol)

! dump grain volumes
if (image1) write (*,"(i0)") grainvol

! deallocate all arrays
call cgca_ds(space)
call cgca_dv(grainvol)
call cgca_dv(fracvol)
call cgca_drt(grt)

end program testAAV

!*roboend*
