!*robodoc*u* tests/testABY
!  NAME
!    testABY
!  SYNOPSIS

!$Id: testABY.f90 380 2017-03-22 11:03:09Z mexas $

program testABY

!  PURPOSE
!    Testing routines from the linked list module, cgca_m2lnklst.
!  NOTE
!    All routines are serial, so no need to use multiple images.
!    A single image will do.
!  DESCRIPTION
!    Verify data integrity and compare timings of
!    cgca_swci and cgca_pc - i.e. a serial writer from image 1
!    and a parallel direct access shared writer. The latter is
!    a non standard Cray extension. So don't try to run this
!    on non-Cray machines.
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

type( cgca_lnklst_node ), pointer :: head, tmp
type( cgca_lnklst_tpayld ) :: box
integer :: img, i, stat

! Only image 1 works
img = this_image()
main: if ( img .eq. 1 ) then

! print a banner
call banner("ABY")

! initialise head
box%lwr = (/ 1, 2, 3 /)
box%upr = (/ 4, 5, 6 /)
call cgca_inithead( head, box )

! dump the list
write (*,*) "The list after initialising head"
call cgca_lstdmp( head )

! add several nodes on top of head
box%lwr = (/ 7, 8, 9 /)
box%upr = (/ 10, 11, 12 /)
call cgca_addhead( head, box )
box%lwr = (/ 13, 14, 15 /)
box%upr = (/ 16, 17, 18 /)
call cgca_addhead( head, box )
box%lwr = (/ 19, 20, 21 /)
box%upr = (/ 22, 23, 24 /)
call cgca_addhead( head, box )
box%lwr = (/ 25, 26, 27 /)
box%upr = (/ 28, 29, 30 /)
call cgca_addhead( head, box )
box%lwr = (/ 31, 32, 33 /)
box%upr = (/ 34, 35, 36 /)
call cgca_addhead( head, box )

! dump the list
write (*,*) "The list after adding 5 nodes on top of head"
call cgca_lstdmp( head )

! add few nodes 3 levels lower than head
tmp => head 
do i=1,2
  tmp => tmp%next 
end do
box%lwr = (/ 101, 102, 103 /)
box%upr = (/ 104, 105, 106 /)
call cgca_addmiddle( tmp, box ) 
box%lwr = (/ 107, 108, 109 /)
box%upr = (/ 110, 111, 112 /)
call cgca_addmiddle( tmp, box ) 
box%lwr = (/ 113, 114, 115 /)
box%upr = (/ 116, 117, 118 /)
call cgca_addmiddle( tmp, box ) 
box%lwr = (/ 119, 120, 121 /)
box%upr = (/ 122, 123, 124 /)
call cgca_addmiddle( tmp, box ) 

! dump the list
write (*,*) "The list after adding 4 nodes 3 levels lower than head"
call cgca_lstdmp( head )

! remove the head node several times
do i=1,2
  call cgca_rmhead( head, stat )
  if ( stat .eq. 1 ) then
    write (*,*) "Reached NULL"
    exit
  end if
end do

! dump the list
write (*,*) "Removed two head nodes, now the list is:"
call cgca_lstdmp( head )

! remove few middle nodes 3 levels below head
tmp => head 
do i=1,2
  tmp => tmp%next 
end do
do i = 1,3
  call cgca_rmmiddle( tmp, stat )
  if ( stat .eq. 1 ) write (*,*) "WARN: cgca_rmmiddle: Reached NULL"
end do

! dump the list
write (*,*) "Removed 3 middle nodes 3 levels below head, the new list:"
call cgca_lstdmp( head )

! continue removing till NULL has been reached
write (*,*) "Continue removing from that point, till NULL has been reached."
remove: do
  call cgca_rmmiddle( tmp, stat )
!write (*,*) "stat: ", stat
  if ( stat .eq. 1 ) then
    write (*,*) "Reached NULL"
    exit remove
  end if
end do remove

! dump the list
write (*,*) "The list after removing all nodes till NULL"
call cgca_lstdmp( head )

! add 3 more nodes on top of head
box%lwr = (/ 2201, 2201, 2201 /)
box%upr = (/ 2333, 2333, 2333 /)
call cgca_addhead( head, box )
box%lwr = (/ 3201, 3201, 3201 /)
box%upr = (/ 3333, 3333, 3333 /)
call cgca_addhead( head, box )
box%lwr = (/ 4201, 4201, 4201 /)
box%upr = (/ 4333, 4333, 4333 /)
call cgca_addhead( head, box )

! dump the list
write (*,*) "The list after adding 3 more nodes on top of head"
call cgca_lstdmp( head )

! remove all head nodes till NULL has been reached
write (*,*) "Removing all head nodes"
do
  call cgca_rmhead( head, stat )
  if ( stat .eq. 1 ) then
    write (*,*) "Reached NULL. associated( head ) = ", associated( head )
    exit
  end if
  if ( .not. associated( head%next ) ) then
    write (*,*) "The list when head%next is not associated, just a single node left"
    call cgca_lstdmp( head )
  end if 
end do

! dump the list
write (*,*) "The list after removing all head nodes"
call cgca_lstdmp( head )

end if main

end program testABY

!*roboend*
