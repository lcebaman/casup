!*robodoc*u* tests/testAAS
!  NAME
!    testAAS
!  SYNOPSIS

!$Id: testAAS.f90 380 2017-03-22 11:03:09Z mexas $

program testAAS
      
!  PURPOSE
!    Checking: cgca_dacf
!  DESCRIPTION
!    Checking deactivation of crack flanks.
!    Note that in this test we need to use at least 2 state types
!    in the space coarray.
!
!    In image1 coarray only, cells on one plane are given "crack edge"
!    values. Then cgca_dacf is called. As all interior cells on this planes
!    are switched to "crack flanks".
!  NOTES
!    The program must be called with 2 command line arguments,
!    both positive integers. These are codimensions along 1 and 2.
!    The number of images must be such that
!    codimension3 = num_images()/( codimension1 * codimension3 )
!    is a positive integer. Example:
!      cafrun -np 16 ./testAAS.x 2 2      ! OpenCoarrays
!    or
!      ./testAAS.x 2 2                    ! Intel, Cray
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
 real(kind=rdef) :: s1(3)

 integer(kind=idef) :: l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3, &
   cou3,   &
   nuc,    & ! number of nuclei in the model
   nimages,codim(3)[*]
 integer(kind=iarr),allocatable :: space(:,:,:,:)[:,:,:]
 integer(kind=ilrg),allocatable :: grainvol(:)[:,:,:], fracvol(:)[:,:,:]
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
 call banner("AAS")
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
u1=10
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

! initialise coarray: a single grain and no damage
space(:,:,:,cgca_state_type_grain) = 1
space(:,:,:,cgca_state_type_frac) = cgca_intact_state

! modify the fracture state array:
! set image 1 cells on some plane to an edge state
space(l1+1:u1-1,l2+1:u2-1,u3/2,cgca_state_type_frac)[col1,col2,col3] = &
 cgca_clvg_state_100_edge

! allocate grain and fracture volumes
call cgca_av(1,nuc,col1,cou1,col2,cou2,col3,grainvol)
call cgca_av(1,nuc,col1,cou1,col2,cou2,col3,fracvol)

! allocate rotation tensors
call cgca_art(1,nuc,col1,cou1,col2,cou2,col3,grt)

! assign rotation tensors, sync all inside
call cgca_rt(grt)

! set s1, although not used in this test
s1 = (/ 1.0, 1.0, 1.0 /)

! dump both the grain and the fracture arrays
call cgca_swci(space,cgca_state_type_grain,10,"zg1.raw")
call cgca_swci(space,cgca_state_type_frac,10,"zf1.raw")
sync all

! deactivate flanks
call cgca_dacf(space,yesdebug)
sync all

! dump both the grain and the fracture arrays
call cgca_swci(space,cgca_state_type_grain,10,"zg2.raw")
call cgca_swci(space,cgca_state_type_frac,10,"zf2.raw")

! deallocate all arrays
call cgca_ds(space)
call cgca_dv(grainvol)
call cgca_dv(fracvol)
call cgca_drt(grt)

end program testAAS

!*roboend*
