!*robodoc*u* tests/testABK
!  NAME
!    testABK
!  SYNOPSIS

!$Id: testABK.f90 380 2017-03-22 11:03:09Z mexas $

program testABK

!  PURPOSE
!    Checking: cgca_av, cgca_dv, cgca_gv
!  DESCRIPTION
!    Checking and timing grain volume array calculation on Cray.
!    Use with CRAY only!
!  NOTES
!    The program must be called with 2 command line arguments,
!    both positive integers. These are codimensions along 1 and 2.
!    The number of images must be such that
!    codimension3 = num_images()/( codimension1 * codimension3 )
!    is a positive integer. Example:
!      ./testABK.x 2 2
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
 real,parameter :: gigabyte=real(2**30), resolution=1.0e-5

 integer(kind=idef) :: l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3,cou3, &
   nuc,    & ! number of nuclei in the model
   nimages, img, codim(3)[*], i
 integer(kind=iarr),allocatable :: space(:,:,:,:)[:,:,:]
 integer(kind=ilrg),allocatable :: grainvol(:)[:,:,:]
 integer(kind=ilrg) :: icells, mcells

 logical(kind=ldef) :: solid

 real :: image_storage

 integer :: pat_status

!*********************************************************************72
! first executable statement

nimages = num_images()
img = this_image()

! do a check on image 1
if ( img .eq. 1 ) then
 call getcodim(nimages,codim)
 ! print a banner
 call banner("ABK")
 ! print the parameter values
 call cgca_pdmp
 write (*,'(a,i0,a)') "running on ", nimages, " images in a 3D grid"
end if

! all images read codim from image 1
if ( img .eq. 1 ) then
 sync images(*)
else
 sync images(1)
 codim(:) = codim(:)[1]
end if

l1=1
l2=l1
l3=l1

! The array size is only controlled by this value
! in this program.
u1=128
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
nuc = int( resolution * real(mcells) )

if ( img .eq. 1 ) then
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

! allocate space and grain volume coarrays
call cgca_as(l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3,1,space)
call cgca_av(1,nuc,col1,cou1,col2,cou2,col3,grainvol)

! initialise coarray
space = cgca_liquid_state

! populate nuclei
call cgca_nr(space,nuc,nodebug)

! solidify
call cgca_sld(space,periodicbc,0,1,solid)

! calculate volumes with my routine
call pat_region_begin(1,"_gv",pat_status)
call cgca_gv(space,grainvol)
call pat_region_end(1,pat_status)

! dump grain volumes
if ( img .eq. 1 ) then
 write (*,*) "results from cgca_gv"
 do i=lbound(grainvol,dim=1),ubound(grainvol,dim=1)
  write (*,"(i0,tr1,i0)") i, grainvol(i)
 end do
end if

! now calculate volumes using CO_SUM intrinsic
call pat_region_begin(2,"_gvl+co_sum",pat_status)
call cgca_gvl(space,grainvol)
call co_sum(grainvol)
call pat_region_end(2,pat_status)

! dump grain volumes
if ( img .eq. 1 ) then
 write (*,*) "results from cgca_gvl + co_sum"
 do i=lbound(grainvol,dim=1),ubound(grainvol,dim=1)
  write (*,"(i0,tr1,i0)") i, grainvol(i)
 end do
end if

! deallocate all arrays
call cgca_ds(space)
call cgca_dv(grainvol)

end program testABK

!*roboend*
