!*robodoc*u* tests/testAAK
!  NAME
!    testAAK
!  SYNOPSIS

!$Id: testAAK.f90 529 2018-03-26 11:25:45Z mexas $

program testAAK

!  PURPOSE
!    Checking: cgca_sld, cgca_nr, cgca_irs, cgca_hxg
!  DESCRIPTION
!    This program is designed to assess the max coarray dimension,
!    assuming cubic coarray, and cubic coarray grid.
!    Run it until it fails to fit in memory.
!    The array is not written to disk.
!
!    Note that on HPC systems it makes sense to maximise the
!    number of processors, and run as quickly as possible, rather
!    then maximise the memory used by each processor. So the purpose
!    of this test is limited to computers with limited number of
!    processors, where achieving a big model requires using all
!    memory.
!
!    On the other hand, CGPACK, at present, does not scale well,
!    so using lots of memory per node and fewer nodes is best
!    for performance for now.
!  NOTES
!    The program must be called with 2 command line arguments,
!    both positive integers. These are codimensions along 1 and 2.
!    The number of images must be such that
!    codimension3 = num_images()/( codimension1 * codimension3 )
!    is a positive integer. Example:
!      cafrun -np 16 ./testAAK.x 2 2      ! OpenCoarrays
!    or
!      ./testAAK.x 2 2                    ! Intel, Cray
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
 real,parameter :: gigabyte=real(2**30)

 integer(kind=idef) :: l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3, &
   cou3,nimages,myimage,codim(3)[*]
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
 call banner("AAK")
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

u1=10
u2=u1
u3=u1

col1=1
cou1=codim(1)-col1+1
col2=1
cou2=codim(2)-col2+1
col3=1
cou3=codim(3)-col3+1

! initialise random number seed
call cgca_irs(nodebug)

main: do

if (myimage .eq. 1) then
  write (*,'(a,2(i0,":",i0,","),i0,":",i0,")")') &
    "bounds: (",l1,u1,l2,u2,l3,u3
  write (*,'(a,2(i0,":",i0,","),i0,":",i0,")")') &
    "cobounds: (",col1,cou1,col2,cou2,col3,cou3

  ! total number of cells in a coarray
  icells = int(u1-l1+1,kind=ilrg) * int(u2-l2+1,kind=ilrg) * &
    int(u3-l3+1,kind=ilrg)

  ! total number of cells in the model
  mcells = icells * int(codim(1),kind=ilrg) * int(codim(2),kind=ilrg) * &
    int(codim(3),kind=ilrg)

  ! An absolute minimum of storage, in GB, per image.
  ! A factor of 2 is used because will call _sld, which
  ! allocates another array of the same size and kind as
  ! coarray.
  image_storage = real(2*icells) * real(storage_size(space1)/8)/gigabyte 

  write (*,'(a,i0,a)') "Each image has ",icells, " cells"
  write (*,'(a,i0,a)') "The model has ", mcells, " cells"
  write (*,'(a,es9.2,a)') "Each image will use at least ", &
    image_storage, " GB memory"
end if

! allocate coarray
call cgca_as(l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3,1,space1)

! initialise coarray to all solid
space1 = int( this_image(), kind=iarr )

sync all

! make only one liquid cell
space1(l1,l2,l3,cgca_state_type_grain)[col1,col2,col3] = &
 cgca_liquid_state

!call cgca_swci(space1,10,"z0.raw")
sync all

! solidify 1
call cgca_sld(space1,periodicbc,0,0,solid)

!call cgca_swci(space1,10,"z9end.raw")

call cgca_ds(space1)

! if all is fine, continue
if (myimage .eq. 1) write (*,*) "ok"

! increase the size of coarray
u1 = u1+1
u2 = u1
u3 = u1

end do main

end program testAAK

!*roboend*
