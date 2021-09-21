!*robodoc*u* tests/testABG
!  NAME
!    testABG
!  SYNOPSIS

!$Id: testABG.f90 380 2017-03-22 11:03:09Z mexas $

program testABG

!  PURPOSE
!    Timing cgca_sld1
!  DESCRIPTION
!    cgca_sld1 is the second attempt of solidification. It uses
!    SYNC IMAGES to enforce the order of a collective operation.
!    However, like cgca_sld, it still does it one image at a time.
!    It is more complex to program than cgca_sld, but hopefully
!    slightly more efficient.
!
!    This is timing test only. Use other tests to check the correctness
!    of cgca_sld1.
!  NOTES
!    The program must be called with 2 command line arguments,
!    both positive integers. These are codimensions along 1 and 2.
!    The number of images must be such that
!    codimension3 = num_images()/( codimension1 * codimension3 )
!    is a positive integer. Example:
!      cafrun -np 16 ./testABG.x 2 2      ! OpenCoarrays
!    or
!      ./testABG.x 2 2                    ! Intel, Cray
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

real,parameter :: gigabyte=real(2**30), resolution=1.0e-5
logical(kind=ldef),parameter :: yesdebug = .true., nodebug = .false.

integer(kind=idef) :: l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3,cou3, &
 nuc,    & ! number of nuclei in the model
 nimages, codim(3)[*]
integer(kind=iarr),allocatable :: space(:,:,:,:)[:,:,:]
integer(kind=ilrg) :: icells, mcells, img

logical(kind=ldef) :: solid

real :: image_storage, time1, time2

!*********************************************************************72
! first executable statement

img = this_image()
nimages = num_images()

! do a check on image 1
if ( img .eq. 1 ) then
 call getcodim(nimages,codim)
 ! print a banner
 call banner("ABG")
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
u1=512
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
nuc = int( resolution * real( mcells ) )

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

! allocate coarray
call cgca_as(l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3,1,space)

! initialise coarray to liquid
space = cgca_liquid_state 

! populate nuclei
call cgca_nr(space,nuc,nodebug)

! solidify
call cpu_time(time1)
call cgca_sld1(space,0,1,solid)
call cpu_time(time2)
write (*,*) "img", img, "time, s", time2-time1

! dump the model
!call cgca_swci(space,cgca_state_type_grain,10,'z.raw')

! deallocate all arrays
call cgca_ds(space)

end program testABG

!*roboend*
