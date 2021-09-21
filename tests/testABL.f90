!*robodoc*u* tests/testABL
!  NAME
!    testABL
!  SYNOPSIS

!$Id: testABL.f90 380 2017-03-22 11:03:09Z mexas $

program testABL

!  PURPOSE
!    Timing and checking cgca_sld3
!  DESCRIPTION
!    Timing solidification where the check for the complete
!    solidification of the whole model is done with CO_SUM.
!  NOTES
!    The program must be called with 2 command line arguments,
!    both positive integers. These are codimensions along 1 and 2.
!    The number of images must be such that
!    codimension3 = num_images()/( codimension1 * codimension3 )
!    is a positive integer. Example:
!      cafrun -np 16 ./testABL.x 2 2      ! OpenCoarrays
!    or
!      ./testABL.x 2 2                    ! Intel, Cray
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

real,parameter :: gigabyte=real(2**30), resolution=1.0e-5, &
 loge2 = log(real(2))
logical(kind=ldef),parameter :: yesdebug = .true., nodebug = .false.

integer(kind=idef) :: l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3,cou3, &
 nuc,    & ! number of nuclei in the model
 nimages, codim(3)[*], p
integer(kind=iarr),allocatable :: space(:,:,:,:)[:,:,:]
integer(kind=ilrg) :: icells, mcells, img

logical(kind=ldef) :: solid

real :: image_storage, time1, time2

!*********************************************************************72
! first executable statement

img = this_image()
nimages = num_images()

! check than nimages is a power of 2
p = nint(log(real(nimages))/loge2)
if ( 2**p .ne. nimages) error stop "number of images is not a power of 2"

! do a check on image 1
if ( img .eq. 1 ) then
 call getcodim(nimages,codim)
 ! print a banner
 call banner("ABL")
 ! print the parameter values
 call cgca_pdmp
 write (*,'(a,i0,a)') "running on ", nimages, " images in a 3D grid"
 ! dump the value of p
 write (*,"(a,i0)") "p=",p
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
u1=16
u2=32
u3=32

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
call cgca_sld3(space,0,1,solid)
call cpu_time(time2)
write (*,*) "img", img, "time, s", time2-time1

! dump the model
!call cgca_swci(space,cgca_state_type_grain,10,'z.raw')

! deallocate all arrays
call cgca_ds(space)

end program testABL

!*roboend*
