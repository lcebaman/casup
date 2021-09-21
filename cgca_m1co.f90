!*robodoc*m* CGPACK/cgca_m1co
!  NAME
!    cgca_m1co
!  SYNOPSIS

!$Id: cgca_m1co.f90 537 2018-04-03 13:57:55Z mexas $

module cgca_m1co

!  DESCRIPTION
!    Lowest level module, contains named global parameters
!    and variables.
!    Contains routines which do not use modules (level 1 routines).
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    Various global parameters and variables
!  USED BY
!    Most/all higher level modules.
!  SOURCE

use iso_fortran_env
implicit none

!*roboend*

!*********************************************************************72
! parameters
!*********************************************************************72

!*robodoc*p* cgca_m1co/rdef
!  PARAMETER
integer, parameter :: rdef = selected_real_kind( 6, 30 )
!  DESCRIPTION
!    Default real kind
!*roboend*


!*robodoc*p* cgca_m1co/rlrg
!  PARAMETER

integer, parameter :: rlrg = selected_real_kind( 15, 300 )

!  DESCRIPTION
!    High precision real kind, most likely will be double
!    precision
!*roboend*

!*robodoc*p* cgca_m1co/idef
!  PARAMETER

integer, parameter :: idef = selected_int_kind( 8 )

!  DESCRIPTION
!    Default integer kind
!*roboend*

!*robodoc*p* cgca_m1co/ca_range
!  PARAMETER

integer, parameter :: ca_range = 8

!  DESCRIPTION
!    Range argument for CA data, used in SELECTED_INT_KIND
!    and in MPI MPI_Type_create_f90_integer.
!    The use of the same range parameter will ensure that the
!    integer data kinds (MPI types) are matching.
!    Typical values (no guarantee) are:
!      ca_range     integer size, byte
!       10            8 
!       8             4, default integer
!       4             2
!       2             1
!*roboend*

!*robodoc*p* cgca_m1co/iarr
!  PARAMETER

integer, parameter :: iarr = selected_int_kind( ca_range )

!  DESCRIPTION
!    Integer kind for CA arrays
!*roboend*

!*robodoc*p* cgca_m1co/ilrg
!  PARAMETER

integer, parameter :: ilrg = selected_int_kind( 10 )

! DESCRIPTION
!   Integer kind for large numbers, e.g. volumes, total
!   number of cells, etc.
!*roboend*

!*robodoc*p* cgca_m1co/ldef
!  PARAMETER

integer, parameter :: ldef = kind( .true. )

!  DESCRIPTION
!    Default logical kind
!*roboend*

!*robodoc*p* cgca_m1co/pi
!  PARAMETER

real( kind=rdef ), parameter :: cgca_pi = 3.14159265358979323846264338_rdef

!  DESCRIPTION
!    pi
!*roboend*

!*robodoc*p* cgca_m1co/cgca_state_type_grain
!  PARAMETER

integer( kind=idef ), parameter :: cgca_state_type_grain = 1_idef

! DESCRIPTION
!   Cell state type for grains
!*roboend*

!*robodoc*p* cgca_m1co/cgca_state_type_frac
!  PARAMETER

integer( kind=idef ), parameter :: cgca_state_type_frac = 2_idef

! DESCRIPTION
!   Cell state type for fractures
!*roboend*

!*robodoc*p* cgca_m1co/cgca_liquid_state
!  PARAMETER

integer( kind=iarr ), parameter :: cgca_liquid_state = 0_iarr

!  DESCRIPTION
!    Liquid phase, cell state of type cgca_state_type_grain.
!    All states of the same type must be unique.
!  NOTE
!    All grains must be greater than the liquid state, i.e. all
!    grains must be positive.
!*roboend*

!*********************************************************************72
! Fracture layer states in decreasing value.
! All states in a layer must be unique.
!*********************************************************************72

!*robodoc*p* cgca_m1co/cgca_state_null
!  PARAMETER
integer( kind=iarr ), parameter :: cgca_state_null = huge( 0_iarr )
!  DESCRIPTION
!    An inactive (null, void, nonexistent) state of a cell in the
!    fracture layer, cell state of type cgca_state_type_frac.
!    This state is given to cells which are outside of the FE model.
!    These cells are not analysed at all in any fracture routines.
!*roboend*


!*robodoc*p* cgca_m1co/cgca_gb_state_intact
!  PARAMETER
integer( kind=iarr ), parameter :: cgca_gb_state_intact = 2_iarr
!  DESCRIPTION
!    Intact grain boundary, cell state of type cgca_state_type_frac.
!*roboend*


!*robodoc*p* cgca_m1co/cgca_gb_state_fractured
!  PARAMETER
integer( kind=iarr ), parameter :: cgca_gb_state_fractured = 1_iarr
!  DESCRIPTION
!    Fractured grain boundary, cell state of type cgca_state_type_frac.
!*roboend*


!*robodoc*p* cgca_m1co/cgca_intact_state
!  PARAMETER
integer( kind=iarr ), parameter :: cgca_intact_state = 0_iarr
!  DESCRIPTION
!    Intact state for fracture array, cell state of type
!    cgca_state_type_frac.
!*roboend*


!*robodoc*p* cgca_m1co/cgca_clvg_state_100_flank
!  PARAMETER
integer( kind=iarr ), parameter :: cgca_clvg_state_100_flank = -1_iarr
!  DESCRIPTION
!    Flanks of a cleavage crack of {100} family.
!    Cell state of type cgca_state_type_frac.
!*roboend*


!*robodoc*p* cgca_m1co/cgca_clvg_state_100_edge
!  PARAMETER
integer( kind=iarr ), parameter :: cgca_clvg_state_100_edge = -2_iarr
!  DESCRIPTION
!    Edges of a cleavage crack of {100} family.
!    Cell state of type cgca_state_type_frac.
!*roboend*


!*robodoc*p* cgca_m1co/cgca_clvg_state_110_flank
!  PARAMETER
integer( kind=iarr ), parameter :: cgca_clvg_state_110_flank = -3_iarr
!  DESCRIPTION
!    Flanks of a cleavage crack of {110} family.
!    Cell state of type cgca_state_type_frac.
!*roboend*


!*robodoc*p* cgca_m1co/cgca_clvg_state_110_edge
!  PARAMETER
integer( kind=iarr ), parameter :: cgca_clvg_state_110_edge = -4_iarr
!  DESCRIPTION
!    Edges of a cleavage crack of {110} family.
!    Cell state of type cgca_state_type_frac.
!*roboend*


!*robodoc*p* cgca_m1co/cgca_clvg_state_111_flank
!  PARAMETER
integer( kind=iarr ), parameter :: cgca_clvg_state_111_flank = -5_iarr
!  DESCRIPTION
!    Flanks of a cleavage crack of {111} family.
!    Cell state of type cgca_state_type_frac.
!*roboend*


!*robodoc*p* cgca_m1co/cgca_clvg_state_111_edge
!  PARAMETER
integer( kind=iarr ), parameter :: cgca_clvg_state_111_edge = -6_iarr
!  DESCRIPTION
!    Edges of a cleavage crack of {111} family.
!    Cell state of type cgca_state_type_frac.
!*roboend*


!*robodoc*p* cgca_m1co/cgca_clvg_lowest_state
!  PARAMETER
integer( kind=iarr ), parameter ::                                     &
 cgca_clvg_lowest_state = cgca_clvg_state_111_edge
!  DESCRIPTION
!    The the lowest cleavage state, used for sizing the lower
!    bound of e.g. the grain volume array.
!    Cell state of type cgca_state_type_frac.
!*roboend*


!*robodoc*p* cgca_m1co/cgca_clvg_states_flank
!  PARAMETER
integer(kind=iarr), parameter :: cgca_clvg_states_flank(3) =           &
  (/ cgca_clvg_state_100_flank,                                        & 
     cgca_clvg_state_110_flank,                                        &
     cgca_clvg_state_111_flank /)
!  DESCRIPTION
!    Array to store all flank cleavage states.
!    Cell state of type cgca_state_type_frac.
!*roboend*


!*robodoc*p* cgca_m1co/cgca_clvg_states_edge
!  PARAMETER
integer(kind=iarr), parameter :: cgca_clvg_states_edge(3)  =           &
  (/ cgca_clvg_state_100_edge,                                         &
     cgca_clvg_state_110_edge,                                         &
     cgca_clvg_state_111_edge /)
!  DESCRIPTION
!    Array to store all edge cleavage states.
!    Cell state of type cgca_state_type_frac.
!*roboend*


!*robodoc*p* cgca_m1co/cgca_clvg_states
!  PARAMETER
integer(kind=iarr), parameter ::                                       &
 cgca_clvg_states( size(cgca_clvg_states_flank) +                      &
                   size(cgca_clvg_states_edge) ) =                     &
  (/ cgca_clvg_states_flank, cgca_clvg_states_edge /)
!  DESCRIPTION
!    Array to store all cleavage states, flanks and edges.
!    Cell state of type cgca_state_type_frac.
!*roboend*


!*robodoc*p* cgca_m1co/cgca_frac_states
!  PARAMETER
integer( kind=iarr ), parameter ::                                     &
 cgca_frac_states( size(cgca_clvg_states) + 1 ) =                      &
  (/ cgca_gb_state_fractured, cgca_clvg_states /)
!  DESCRIPTION
!    Array to store all fracture states: cleavage, fractured GB, etc.
!    Cell state of type cgca_state_type_frac.
!*roboend*


!*robodoc*p* cgca_m1co/cgca_lowest_state
!  PARAMETER
integer(kind=iarr), parameter :: cgca_lowest_state = -huge(0_iarr)
!  DESCRIPTION
!    Lowest possible state in the model
!*roboend*


!*robodoc*p* cgca_m1co/cgca_gcupd_size1
!  PARAMETER
integer( kind=idef ), parameter :: cgca_gcupd_size1 = 100_idef
!  DESCRIPTION
!    This is size1 for the allocatable array coarray gcupd_arr
!    in module cgca_m3clvg. The value is the maximum number
!    of fractured grain boundaries on an image in a single
!    CA iteration. 100 is probably an overkill. Perhaps 2-3
!    will be enough, but it's not too big a waste.
!*roboend*


!*robodoc*p* cgca_m1co/cgca_gcupd_size2
!  PARAMETER
integer( kind=idef ), parameter :: cgca_gcupd_size2 = 3_idef
!  DESCRIPTION
!    This is size2 for the allocatable array coarray gcupd_arr
!    in module cgca_m3clvg. The value is the number of useful
!    fields for each entry. At present only 3 entries planned:
!    (1) grain 1, (2) grain 2, (3) state of grain boundary
!    between grains 1 and 2.
!*roboend*


!*robodoc*p* cgca_m1co/cgca_scodim
!  PARAMETER
integer( kind=idef ), parameter :: cgca_scodim = 3
!  DESCRIPTION
!    Number of spatial codimensions for the main space coarray
!    and other coarrays with more that 1 codimensions.
!*roboend*

!*********************************************************************72
! variables
!*********************************************************************72

!*robodoc*v* cgca_m1co/cgca_slcob
!  PARAMETER
integer( kind=idef ) :: cgca_slcob( cgca_scodim )
!  DESCRIPTION
!    Lower cobounds of the space coarray and other coarrays with
!    more that 1 codimensions.
!*roboend*

!*robodoc*v* cgca_m1co/cgca_sucob
!  PARAMETER
integer( kind=idef ) :: cgca_sucob( cgca_scodim )
!  DESCRIPTION
!    Upper cobounds of the space coarray and other coarrays with
!    more that 1 codimensions.
!*roboend*

end module cgca_m1co
