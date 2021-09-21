!*robodoc*u* tests/testAAD
!  NAME
!    testAAD
!  SYNOPSIS

!$Id: testAAD.f90 380 2017-03-22 11:03:09Z mexas $

program testAAD

!  PURPOSE
!    Checking: cgca_gl, cgca_lg.
!  DESCRIPTION
!    Checking the mapping routines.
!  NOTES
!    The program must be called with 2 command line arguments,
!    both positive integers. These are codimensions along 1 and 2.
!    The number of images must be such that
!    codimension3 = num_images()/( codimension1 * codimension3 )
!    is a positive integer. Example:
!      cafrun -np 16 ./testAAD.x 2 2      ! OpenCoarrays
!    or
!      ./testAAD.x 2 2                    ! Intel, Cray
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

 integer(kind=idef) :: l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3, &
   nimages,super_in(3),super_out(3),imgpos(3),local(3),test,  &
   codim(3)[*]
 integer(kind=iarr),allocatable :: space1(:,:,:,:)[:,:,:]

!*********************************************************************72
! first executable statement

nimages=num_images()

! do a check on image 1
if (this_image() .eq. 1) then
 call getcodim(nimages,codim)
 ! print a banner
 call banner("AAD")
 ! print the parameter values
 call cgca_pdmp
 write (*,'(a,i0,a)') "running on ", nimages, " images in a 3D grid"
 write (*,*) "codim:", codim
end if

sync all

codim(:) = codim(:)[1]

if (this_image() .eq. 2) call system("sleep 1")

!***********************************************************
! test 1

test = 1

l1=1
u1=10
l2=1
u2=10
l3=1
u3=10
col1=1
cou1=2
col2=1
cou2=2
col3=1

if (this_image() .eq. 1) then
  write (*,*) "test 1"
  write (*,'(a,2(i0,":",i0,","),i0,":",i0,")")') &
    "bounds: (",l1,u1,l2,u2,l3,u3
  write (*,'(a,2(i0,":",i0,","),i0,":",i0,")")') &
    "cobounds: (",col1,cou1,col2,cou2,col3, &
     nimages/((cou1-col1+1)*(cou2-col2+1))+col3-1
end if

call cgca_as(l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3,1,space1)

super_in=(/20,10,1/)

if (this_image() .eq. 1) then
  call cgca_gl(super_in,space1,imgpos,local)
  write (*,'(a,3(i0,tr1))') "super in ", super_in
  write (*,'(a,3(i0,tr1))') "imgpos ", imgpos
  write (*,'(a,3(i0,tr1))') "local  ", local
  call cgca_lg(imgpos,local,space1,super_out)
  write (*,'(a,3(i0,tr1))') "super out ", super_out

  if (all(super_out .eq. super_in)) then
    write (*,'(a,i0,a)') "test ", test, " passed ok"
  else
    write (*,'(a,i0,a)') "****** ERROR: TEST ", test, " FAILED"
    error stop
  end if

end if
sync all

call cgca_ds(space1)

!***********************************************************
! test 2

test = 2

l1=-7
u1=-3
l2=-10
u2=10
l3=-5
u3=0
col1=-2
cou1=-1
col2=-1
cou2=0
col3=8

if (this_image() .eq. 1) then
  write (*,*) "test  2"
  write (*,'(a,2(i0,":",i0,","),i0,":",i0,")")') &
    "bounds: (",l1,u1,l2,u2,l3,u3
  write (*,'(a,2(i0,":",i0,","),i0,":",i0,")")') &
    "cobounds: (",col1,cou1,col2,cou2,col3, &
     nimages/((cou1-col1+1)*(cou2-col2+1))+col3-1
end if

call cgca_as(l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3,2,space1)

super_in=(/5,4,3/)

if (this_image() .eq. 1) then
  call cgca_gl(super_in,space1,imgpos,local)
  write (*,'(a,3(i0,tr1))') "super in ", super_in
  write (*,'(a,3(i0,tr1))') "imgpos ", imgpos
  write (*,'(a,3(i0,tr1))') "local  ", local
  call cgca_lg(imgpos,local,space1,super_out)
  write (*,'(a,3(i0,tr1))') "super out ", super_out

  if (all(super_out .eq. super_in)) then
    write (*,'(a,i0,a)') "test ", test, " passed ok"
  else
    write (*,'(a,i0,a)') "****** ERROR: TEST ", test, " FAILED"
    error stop
  end if

end if
sync all

call cgca_ds(space1)

!***********************************************************
! test 3 

test = 3

l1=73
u1=100
l2=15
u2=20
l3=-99
u3=-90
col1=15
cou1=16
col2=0
cou2=1
col3=-10

if (this_image() .eq. 1) then
  write (*,*) "test 3"
  write (*,'(a,2(i0,":",i0,","),i0,":",i0,")")') &
    "bounds: (",l1,u1,l2,u2,l3,u3
  write (*,'(a,2(i0,":",i0,","),i0,":",i0,")")') &
    "cobounds: (",col1,cou1,col2,cou2,col3, &
     nimages/((cou1-col1+1)*(cou2-col2+1))+col3-1
end if

call cgca_as(l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3,3,space1)

super_in=(/5,4,3/)

if (this_image() .eq. 1) then
  call cgca_gl(super_in,space1,imgpos,local)
  write (*,'(a,3(i0,tr1))') "super in ", super_in
  write (*,'(a,3(i0,tr1))') "imgpos ", imgpos
  write (*,'(a,3(i0,tr1))') "local  ", local
  call cgca_lg(imgpos,local,space1,super_out)
  write (*,'(a,3(i0,tr1))') "super out ", super_out

  if (all(super_out .eq. super_in)) then
    write (*,'(a,i0,a)') "test ", test, " passed ok"
  else
    write (*,'(a,i0,a)') "****** ERROR: TEST ", test, " FAILED"
    error stop
  end if

end if
sync all

call cgca_ds(space1)

end program testAAD

!*roboend*
