!*robodoc*u* tests/testAAC
!  NAME
!    testAAC
!  SYNOPSIS

!$Id: testAAC.f90 536 2018-04-03 12:02:13Z mexas $

program testAAC

!  PURPOSE
!    Checking: cgca_hxi
!  DESCRIPTION
!    Checking internal halo calculations.
!  NOTES
!    The program must be called with 2 command line arguments,
!    both positive integers. These are codimensions along 1 and 2.
!    The number of images must be such that
!    codimension3 = num_images()/( codimension1 * codimension3 )
!    is a positive integer. Example:
!      cafrun -np 16 ./testAAC.x 2 2      ! OpenCoarrays
!    or
!      ./testAAC.x 2 2                    ! Intel, Cray
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

 integer( kind=idef ),parameter :: size1=10, size2=10, size3=10
 integer( kind=idef ) :: nimages,nbhd,i,j,k
 integer( kind=idef ) :: codim(3)[*]
 real :: scaling
 integer(kind=iarr),allocatable :: space1(:,:,:,:)[:,:,:], &
   space2(:,:,:,:)[:,:,:]
logical(kind=ldef) :: image1

!*********************************************************************72
! first executable statement

nimages=num_images()
image1 = .false.
if (this_image() .eq. 1) image1 = .true.

! do a check on image 1
if (image1) then
 call getcodim(nimages,codim)
 ! print a banner
 call banner("AAC")
 ! print the parameter values
 call cgca_pdmp
 write (*,'(a,i0,a)') "running on ", nimages, " images in a 3D grid"
 write (*,*) "codim:", codim
end if

sync all

nbhd=26
scaling=1.0/real(nbhd+1)

codim(:) = codim(:)[1]

if (this_image() .eq. 2) call system("sleep 1")

call cgca_as(1,size1,1,size2,1,size3,1,codim(1),1,codim(2),1,2,space1)
call cgca_as(1,size1,1,size2,1,size3,1,codim(1),1,codim(2),1,2,space2)

if (allocated(space1)) then
  write (*,'(a,i0,a)')"Image:",this_image(), " space1 allocated"
  write (*,'(a,i0,a,3(i0,tr1),a)') &
    "Image: ",this_image()," is ",this_image(space1)," in the grid"
end if

if (allocated(space2)) then
  write (*,'(a,i0,a)')"Image:",this_image(), " space2 allocated"
end if

space1(:,:,:,cgca_state_type_grain) = int( this_image(), kind=iarr )
space1(:,:,:,cgca_state_type_frac) = cgca_intact_state
space2=0

sync all

call cgca_swci(space1,cgca_state_type_grain,10,'z1.raw')
sync all

call cgca_hxi(space1)
sync all

do k=1,size3
do j=1,size2
do i=1,size1
  space2( i, j, k, cgca_state_type_grain ) =                           &
  space1( i, j, k, cgca_state_type_grain ) -                           &
    nint( scaling * sum( space1( i-1:i+1 , j-1:j+1 , k-1:k+1 ,         &
                                 cgca_state_type_grain ) ), kind=iarr )   
end do
end do
end do

sync all

if (this_image() .eq. 3) call system("sleep 2")

call cgca_swci(space2,cgca_state_type_grain,10,'z2.raw')

call cgca_ds(space1)
if (.not. allocated(space1)) &
  write (*,'(a,i0,a)')"Image:",this_image(), " space1 not allocated"

call cgca_ds(space2)
if (.not. allocated(space2)) &
  write (*,'(a,i0,a)')"Image:",this_image(), " space2 not allocated"

end program testAAC

!*roboend*
