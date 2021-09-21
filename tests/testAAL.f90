!*robodoc*u* tests/testAAL
!  NAME
!    testAAL
!  SYNOPSIS

!$Id: testAAL.f90 380 2017-03-22 11:03:09Z mexas $

program testAAL

!  PURPOSE
!    Checking: cgca_sld, cgca_nr, cgca_irs, cgca_hxg
!  DESCRIPTION
!    Still solidification, but this program is designed for
!    timing analysis.
!  NOTES
!    The program must be called with 2 command line arguments,
!    both positive integers. These are codimensions along 1 and 2.
!    The number of images must be such that
!    codimension3 = num_images()/( codimension1 * codimension3 )
!    is a positive integer. Example:
!      cafrun -np 16 ./testAAL.x 2 2      ! OpenCoarrays
!    or
!      ./testAAL.x 2 2                    ! Intel, Cray
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
real, parameter :: gigabyte=real(2**30), resolution=1.0e-5

 integer(kind=idef) :: l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3, &
   cou3, &
   nuc,    & ! number of nuclei in the model
   nimages,myimage,codim(3)[*]
 integer(kind=iarr),allocatable :: space1(:,:,:,:)[:,:,:]
 integer(kind=ilrg) :: icells,mcells

 logical(kind=ldef) :: solid

 real :: image_storage

!*********************************************************************72
! first executable statement

nimages=num_images()
myimage=this_image()

! do a check on image 1
if (myimage .eq. 1) then
 call getcodim(nimages,codim)
 ! print a banner
 call banner("AAL")
 ! print the parameter values
 call cgca_pdmp
 write (*,'(a,i0,a)') "running on ", nimages, " images in a 3D grid"
 write (*,*) "codim:", codim
end if

sync all

codim(:) = codim(:)[1]

if (myimage .eq. 2) call system("sleep 1")

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

! total number of cells in an image
icells = int(u1-l1+1,kind=ilrg) * int(u2-l2+1,kind=ilrg) * &
  int(u3-l3+1,kind=ilrg)

! total number of cells in the model
mcells = icells * int(codim(1),kind=ilrg) * int(codim(2),kind=ilrg) * &
  int(codim(3),kind=ilrg)

! total number of nuclei
nuc = int( resolution * real( mcells ) )

if (myimage .eq. 1) then
  write (*,'(a,2(i0,":",i0,","),i0,":",i0,")")') &
    "bounds: (",l1,u1,l2,u2,l3,u3
  write (*,'(a,2(i0,":",i0,","),i0,":",i0,")")') &
    "cobounds: (",col1,cou1,col2,cou2,col3,cou3

  ! An absolute minimum of storage, in GB, per image.
  ! A factor of 2 is used because will call _sld, which
  ! allocates another array of the same size and kind as
  ! coarray.
  image_storage = real(2 * icells*storage_size(space1)/8)/gigabyte 

  write (*,'(a,i0,a)') "Each image has ",icells, " cells"
  write (*,'(a,i0,a)') "The model has ", mcells, " cells"
  write (*,'(a,i0,a)') "The model has ", nuc, " nuclei"
  write (*,'(a,es9.2,a)') "Each image will use at least ", &
    image_storage, " GB memory"
end if

! initialise random number seed
call cgca_irs(nodebug)

! allocate coarray
call cgca_as(l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3,1,space1)

! initialise coarray to liquid
space1(:,:,:,cgca_state_type_grain) = cgca_liquid_state

sync all

! populate nuclei
call cgca_nr(space1,nuc,nodebug)

!call cgca_swci(space1,10,"z0.raw")
!sync all

! solidify 1
call cgca_sld(space1,periodicbc,0,0,solid)

call cgca_swci(space1,cgca_state_type_grain,10,"z9end.raw")

call cgca_ds(space1)

end program testAAL

!*roboend*
