!*robodoc*u* tests/testAAZ
!  NAME
!    testAAZ
!  SYNOPSIS

!$Id: testAAZ.f90 526 2018-03-25 23:44:51Z mexas $

program testAAZ

!  PURPOSE
!    Checking: cgca_gcf, cgca_gcr
!  DESCRIPTION
!    Checking the grain connectivity boundary status.
!  NOTES
!    The program must be called with 2 command line arguments,
!    both positive integers. These are codimensions along 1 and 2.
!    The number of images must be such that
!    codimension3 = num_images()/( codimension1 * codimension3 )
!    is a positive integer. Example:
!      cafrun -np 16 ./testAAZ.x 2 2      ! OpenCoarrays
!    or
!      ./testAAZ.x 2 2                    ! Intel, Cray
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
logical(kind=ldef),parameter :: yesdebug=.true., nodebug=.false., &
  periodicbc=.true.,  noperiodicbc=.false.

integer(kind=idef) :: l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3, &
 cou3, nuc,nimages,codim(3)[*]
integer(kind=iarr),allocatable :: space(:,:,:,:)[:,:,:]
integer(kind=ilrg) :: icells,mcells

real :: image_storage

logical(kind=ldef) :: solid,image1, intact
character(6) :: image

!*********************************************************************72
! first executable statement

nimages = num_images()
 image1 = .false.
if (this_image() .eq. 1) image1 = .true.

! do a check on image 1
if (image1) then
 call getcodim(nimages,codim)
 ! print a banner
 call banner("AAZ")
 ! print the parameter values
 call cgca_pdmp
 write (*,'(a,i0,a)') "running on ", nimages, " images in a 3D grid"
end if

sync all

codim(:) = codim(:)[1]

l1 = 1
l2 = l1
l3 = l1

u1 = 50
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
nuc = int( resolution * real( mcells ) )

100 format( a,3(i0,":",i0,a) )
 
if (image1) then
  write (*,100) "bounds: (", l1,u1,",", l2,u2,",", l3,u3, ")"
  write (*,100) "cobounds: [", &
    col1,cou1,",", col2,cou2,",", col3,cou3, "]"

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

! allocate space 
call cgca_as(l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3,1,space)

! initialise random number seed
call cgca_irs(nodebug)

! initialise the grain coarray to liquid
space(:,:,:,cgca_state_type_grain) = cgca_liquid_state

sync all

! nuclei, sync all inside
call cgca_nr(space,nuc,yesdebug)

! solidify, implicit sync all inside
call cgca_sld(space,noperiodicbc,0,1,solid)

! update grain connectivity
call cgca_gcu(space)
sync all

! mark gb between two grains as intact
! read the value
call cgca_gcr( 3_iarr, 2_iarr, intact )
if (.not. intact) write (*,*) "image", this_image(), intact 
! set the value
call cgca_gcf( 3_iarr, 2_iarr )
! read the value
call cgca_gcr( 3_iarr, 2_iarr, intact)
if (.not. intact) write (*,*) "image", this_image(), intact 
sync all

! dump grain connectivity to files
write (image,"(i0)") this_image()
call cgca_gcp(ounit=10,fname="cgca_gcp_out"//image)

! dump space array to file
call cgca_swci(space,cgca_state_type_grain,10,"z9end.raw")
sync all

! deallocate all arrays
call cgca_ds(space)
call cgca_dgc

end program testAAZ

!*roboend*
