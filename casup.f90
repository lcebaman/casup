!*robodoc*m* CASUP/casup
!  NAME
!    casup
!  SYNOPSIS

!$Id: casup.f90 533 2018-03-30 14:31:26Z mexas $

module casup

!  DESCRIPTION
!    The top level module for CASUP.
!  NOTES
!    The lowest level modules are level 1, e.g. cgca_m1co.
!    Level 1 modules use no other modules.
!    Level 2 modules use only level 1 modules.
!    Level 3 modules use some level 2 modules and possibly
!    also level 1 modules. And so on. All modules except the
!    top level, this one, are named
!     cgca_mX<name>
!    where X is the level number, starting from 1, and
!    name is the module name.
!
!    The modules group routines dealing with a particular
!    data structure or a problem, e.g. cgca_m3clvg deals
!    with cleavage propagation, cgca_m2gb deals with
!    the grain connectivity array. 
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    All modules in the library, which are not used
!    by other modules directly
!  USED BY
!    none, top level module 
!  SOURCE

use cgca_m1co
use cgca_m1clock

use ca_hx
use cgca_m2alloc
use cgca_m2gb
use cgca_m2geom
use cgca_m2glm
use cgca_m2hdf5
use cgca_m2hx
use cgca_m2lnklst
use cgca_m2mpiio
use cgca_m2netcdf
use cgca_m2out
use cgca_m2pck
use cgca_m2phys
use cgca_m2red
use cgca_m2rnd
use cgca_m2rot
use cgca_m2stat

use cgca_m3clvg
! use cgca_m3clvgt - does not build
use cgca_m3gbf
use cgca_m3nucl
!use cgca_m3pfem
use cgca_m3sld

!use cgca_m4fr

end module casup

!*roboend*
