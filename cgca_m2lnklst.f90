!*robodoc*m* CGPACK/cgca_m2lnklst
!  NAME
!    cgca_m2lnklst
!  SYNOPSIS

!$Id: cgca_m2lnklst.f90 380 2017-03-22 11:03:09Z mexas $

module cgca_m2lnklst

!  DESCRIPTION
!    Module with link list types and routines.
!    The module is mainly useful for linking CGPAK to FE.
!    In case the CA box is sticking outside of the FE model.
!    Routines of this module help effectively find all cells
!    which are inside and outside of the FE model.
!  AUTHOR 
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    All public. Derived types cgca_lnklst_tpayld, cgca_lnklst_node.
!    Routines
!    cgca_inithead, cgca_addhead, cgca_addmiddle, cgca_rmhead,
!    cgca_rmmiddle, cgca_lstdmp.
!  USES
!    cgca_m1co
!  USED BY
!    cgca_m3pfem
!  SOURCE

use cgca_m1co, only : idef
implicit none

private
public ::                                                              &
! derived types
          cgca_lnklst_tpayld, cgca_lnklst_node, &
! routines
          cgca_inithead, &
          cgca_addhead, &
          cgca_addmiddle, cgca_rmhead, cgca_rmmiddle, cgca_lstdmp

!*roboend*


!*robodoc*d* cgca_m2lnklst/cgca_lnklst_tpayld
!  NAME
!    cgca_lnklst_tpayld
!  SYNOPSIS

type cgca_lnklst_tpayld
  integer :: lwr(3), upr(3)
end type cgca_lnklst_tpayld

! DESCRIPTION
!   Payload type for all link list routines.
!   The payload consists of two integer arrays of length 3.
!   The arrays contain lower and upper corner coordinates of a CA box
!   in local CA coord. system. 
! USED BY
!   All routines of module cgca_m2lnklst.
!*roboend*


!*robodoc*d* cgca_m2lnklst/cgca_lnklst_node
!  NAME
!    cgca_lnklst_node
!  SYNOPSIS

type cgca_lnklst_node
  type( cgca_lnklst_tpayld ) :: value
  type( cgca_lnklst_node ), pointer :: next
end type cgca_lnklst_node

! DESCRIPTION
!   A derived type for a node in the linked list. A very traditional
!   type. The payload is of derived type cgca_lnklst_tpayld.
! USED BY
!   All routines of module cgca_m2lnklst.
!*roboend*


contains


!*robodoc*s* cgca_m2lnklst/cgca_inithead
!  NAME
!    cgca_inithead
!  SYNOPSIS

subroutine cgca_inithead( head, payload )

!  OUTPUT

type( cgca_lnklst_node ), pointer, intent( out ) :: head

!  INPUT

type( cgca_lnklst_tpayld ), intent( in ) :: payload

!  SIDE EFFECTS
!    Memory for one entity of type "node" is allocated.
!    The pointer to this memory is returned as "head".
!    The value of that memory is set to "payload".
!  DESCRIPTION
!    This routine initialises the head node of the linked list.
!    The head node is the node at the very top of the list.
!    The head node has nothing higher. This is the only node that
!    can be accessed directly. Access to all other nodes is via
!    the head node and pointers therein.
!  SOURCE

  allocate( head )
  head%value = payload
  head%next => null()
end subroutine cgca_inithead

!*roboend*


!*robodoc*s* cgca_m2lnklst/cgca_addhead
!  NAME
!    cgca_addhead
!  SYNOPSIS

subroutine cgca_addhead( head, payload )

!  INPUTS

type( cgca_lnklst_node ), pointer, intent( inout ) :: head
type( cgca_lnklst_tpayld ), intent( in ) :: payload

!  SIDE EFFECTS
!    Memory for one entity of type "node" is allocated.
!    The value of that memory is set to "payload".
!    The pointer to this memory is returned as "head".
!    The new head points to the old head.
!  DESCRIPTION
!    This routine adds another node on top of the head,
!    i.e. puts another node higher than current head.
!    The new node becomes the new head. It points to the
!    old head.
!  SOURCE

  type( cgca_lnklst_node ), pointer :: tmp

  allocate( tmp )
  tmp%value = head%value
  tmp%next => head%next
  allocate( head )
  head%value = payload
  head%next => tmp
end subroutine cgca_addhead

!*roboend*


!*robodoc*s* cgca_m2lnklst/cgca_addmiddle
!  NAME
!    cgca_addmiddle
!  SYNOPSIS

subroutine cgca_addmiddle( node, payload )

!  INPUTS

type( cgca_lnklst_node ), pointer, intent( in ) :: node
type( cgca_lnklst_tpayld ), intent( in ) :: payload

!  SIDE EFFECTS
!    Memory for one entity of type cgca_lnklst_node is allocated.
!    The value of that memory is set to "payload".
!    The new entity points to where "node" pointed before.
!    Node now points to the new entity.
!    A schematic diagram:
!
!      Before                         After
!
!      node   (next,value)     node     (next,value)
!               \                         |_____
!                 \                             v 
!                   \          new node (next,value)
!                    \                    |_____ 
!                     v                         v
!      node X (next,value)     node X   (next,value)
!
!  DESCRIPTION
!    This routine adds another node *lower* than the given node.
!    Lower here means further from the head and closer to NULL.
!    The new node points to where the old node pointed.
!    The old node points to the new node.
!    So the list length is +1.
!  SOURCE

  type( cgca_lnklst_node ), pointer :: tmp

  allocate( tmp )
  tmp%value = payload
   tmp%next => node%next
  node%next => tmp
end subroutine cgca_addmiddle

!*roboend*


!*robodoc*s* cgca_m2lnklst/cgca_rmhead
!  NAME
!    cgca_rmhead
!  SYNOPSIS

subroutine cgca_rmhead( head, stat )

!  INPUT

type( cgca_lnklst_node ), pointer, intent( inout ) :: head

!  OUTPUT
!    stat - integer, 0 if no problem, 1 if the head node is NULL.

integer( kind=idef ), intent( out ) :: stat

!  SIDE EFFECTS
!    Memory for one entity of type "node" is freed.
!    The pointer to the old head now points to where the old head was
!    pointing. This pointer is returned as "head".
!  DESCRIPTION
!    This routine removes the head node. The list length decreases by 1.
!    The pointer to the old head is given on entry. On exit this pointer
!    points to where the old head was pointing, i.e. one node closer to
!    NULL. If there was only a single node on top of head, then the
!    head will return null (unassociated) and stat will be 1.
!    If there is no head node already, i.g. head is not associated
!    already, head will not be changed and stat of 1 will be returned.
!  SOURCE

type( cgca_lnklst_node ), pointer :: tmp
stat = 0

if ( associated( head ) ) then
  tmp => head
  head => head%next
  deallocate( tmp )
end if

! This pointer is not associated only if NULL has been reached. 
! Do nothing and set the output flag accordingly.
if ( .not. associated( head ) ) stat = 1

end subroutine cgca_rmhead

!*roboend*


!*robodoc*s* cgca_m2lnklst/cgca_rmmiddle
!  NAME
!    cgca_rmmiddle
!  SYNOPSIS

subroutine cgca_rmmiddle( node, stat )

!  INPUT

type( cgca_lnklst_node ), pointer, intent( in ) :: node

!  OUTPUT

integer( kind=idef), intent( out ) :: stat

!  SIDE EFFECTS
!    Memory for one entity of type cgca_lnklst_node is freed.
!
!      Before                            After
!
!      node           (next,value)       node   (next,value)    
!                       |_____                    \       
!                             v                     \       
!      node to remove (next,value)                    \    
!                       |_____                         \           
!                             v                         v
!      node X         (next,value)       node X (next,value)
!
!  DESCRIPTION
!    The node below the given node is removed. Below here means
!    further from the head and closer to NULL. The node that
!    pointed to the node to remove before, now points to where
!    the node to remove was pointing.
!  NOTES
!    On output stat=0 means no problem.
!    If stat=1, then an attempt has been made to remove NULL.
!  SOURCE

  type( cgca_lnklst_node ), pointer :: tmp

  stat = 0
  tmp => node%next
  if ( associated( tmp ) ) then
    node%next => tmp%next
    deallocate( tmp )
  else
    ! This pointer is not associated only if NULL has been reached. 
    ! Do nothing but set the output flag accordingly.
    stat = 1
  end if
end subroutine cgca_rmmiddle

!*roboend*


!*robodoc*s* cgca_m2lnklst/cgca_lstdmp
!  NAME
!    cgca_lstdmp
!  SYNOPSIS

subroutine cgca_lstdmp( head )

!  INPUT

type( cgca_lnklst_node ), pointer, intent( in ) :: head

!  SIDE EFFECTS
!    Values of all nodes are dumped to stdout.
!  DESCRIPTION
!    This routine dumps all nodes, one per line, starting from HEAD,
!    till it reaches NULL.
!  SOURCE

  type( cgca_lnklst_node ), pointer :: tmp

  if ( .not. associated( head ) ) return
  tmp => head
  do
    write (*,*) tmp%value
    tmp => tmp%next
    if ( .not. associated( tmp ) ) exit
  end do
end subroutine cgca_lstdmp

!*roboend*

end module cgca_m2lnklst
