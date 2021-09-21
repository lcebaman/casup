!*robodoc*m* CGPACK/cgca_m2pck
!  NAME
!    cgca_m2pck
!  SYNOPSIS

!$Id: cgca_m2pck.f90 380 2017-03-22 11:03:09Z mexas $

module cgca_m2pck

!  DESCRIPTION
!    Module dealing with checking consistency of the many various
!    global CGPACK parameters set in cgca_m1co.
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    cgca_as, cgca_ds, cgca_av, cgca_dv, cgca_art, cgca_drt
!  USES
!    cgca_m1co
!  USED BY
!    cgca
!  SOURCE

use cgca_m1co
implicit none

private
public :: cgca_pdmp

contains

!*roboend*


!*robodoc*s* cgca_m2pck/cgca_pdmp
!  NAME
!    cgca_pdmp
!  SYNOPSIS

subroutine cgca_pdmp

!  DESCRIPTION
!    Dump global CGPACK parameters from cgca_m1co to stdout.
!    The user might want to see all values in one place.
!  USES
!    cgca_m1co
!  USED BY
!    cgca_m2alloc
!  SOURCE

! ifort 16 still does not support this
! Cray apparently has these 2 functions, but *not* in iso_fortran_env!
! So just comment out for now. If really needed, these can be
! added to the main program.
!write (* , "(a,a)" ) "CGPACK compiled with: ", compiler_version()
!write (* , "(a,a)" ) "CGPACK compiler options: ", compiler_options()
write (*,"(a)") "CGPACK cell state types:"
write (*,"(a40,i0)") "cgca_state_type_grain: ", cgca_state_type_grain
write (*,"(a40,i0)") "cgca_state_type_frac: ", cgca_state_type_frac
write (*,*)
write (*,"(a)") "CGPACK grain layer states:"
write (*,"(a40,i0)") "cgca_liquid_state: ", cgca_liquid_state
write (*,*)
write (*,"(a)") "CGPACK fracture layer states:"
write (*,"(a40,i0)") "cgca_state_null: ", cgca_state_null
write (*,"(a40,i0)") "cgca_gb_state_intact: ", cgca_gb_state_intact
write (*,"(a40,i0)") "cgca_gb_state_fractured: ",                      &
                      cgca_gb_state_fractured
write (*,"(a40,i0)") "cgca_intact_state: ", cgca_intact_state
write (*,"(a40,i0)") "cgca_clvg_state_100_flank: ",                    &
                      cgca_clvg_state_100_flank
write (*,"(a40,i0)") "cgca_clvg_state_100_edge: ",                     &
                      cgca_clvg_state_100_edge
write (*,"(a40,i0)") "cgca_clvg_state_110_flank: ",                    &
                      cgca_clvg_state_110_flank
write (*,"(a40,i0)") "cgca_clvg_state_110_edge: ",                     &
                      cgca_clvg_state_110_edge
write (*,"(a40,i0)") "cgca_clvg_state_111_flank: ",                    &
                      cgca_clvg_state_111_flank
write (*,"(a40,i0)") "cgca_clvg_state_111_edge: ",                     &
                      cgca_clvg_state_111_edge
write (*,"(a40,999(i0,tr1))") "cgca_clvg_states_flank: ",              &
                               cgca_clvg_states_flank
write (*,"(a40,999(i0,tr1))") "cgca_clvg_states_edge: ",               &
                               cgca_clvg_states_edge
write (*,"(a40,999(i0,tr1))") "cgca_clvg_states: ", cgca_clvg_states
write (*,"(a40,999(i0,tr1))") "cgca_frac_states: ", cgca_frac_states
write (*,"(a40,i0)") "cgca_clvg_lowest_state: ", cgca_clvg_lowest_state
write (*,*)
write (*,"(a)") "CGPACK lowest state for both types:"
write (*,"(a40,i0)") "cgca_lowest_state: ", cgca_lowest_state
write (*,*)
write (*,"(a)") "CGPACK kinds:"
write (*,"(a40,i0)") "iarr: ", iarr
write (*,"(a40,i0)") "idef: ", idef
write (*,"(a40,i0)") "ilrg: ", ilrg
write (*,"(a40,i0)") "ldef: ", ldef
write (*,"(a40,i0)") "rdef: ", rdef
write (*,"(a40,i0)") "rlrg: ", rlrg
write (*,*)
write (*,"(a)") "CGPACK other parameters:"
write (*,"(a40,f20.10)") "pi: ", cgca_pi

end subroutine cgca_pdmp

!*roboend*

end module cgca_m2pck
