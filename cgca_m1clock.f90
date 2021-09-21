!*robodoc*m* CGPACK/cgca_m1clock
!  NAME
!    cgca_m1clock
!  SYNOPSIS

!$Id: cgca_m1clock.f90 379 2017-03-22 09:57:10Z mexas $

module cgca_m1clock

!  DESCRIPTION
!    Module with timing clocks
!  AUTHOR
!    Luis Cebamanos, modified by Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    Public functions:
!    cgca_benchtime
!  SOURCE

  implicit none
  private
  public :: cgca_benchtime

  integer, parameter :: dp = kind(1.0d0),                              &
    int32kind = selected_int_kind( 9),                                 &
    int64kind = selected_int_kind(18),                                 &
    intkind = int64kind

  logical,       save :: firstcall = .true.
  real(kind=dp), save :: ticktime = 0.0_dp

!
!  Select high resolution clock
!

  integer(kind = intkind) :: count, rate

contains

real(kind=dp) function benchtick()

  if (firstcall) then

     firstcall = .false.
     call system_clock(count, rate)
     ticktime = 1.0_dp / real(rate, kind=dp)

  end if

  benchtick = ticktime

end function benchtick

!*roboend*


!*robodoc*s* cgca_m1clock/cgca_benchtime
!  NAME
!    cgca_benchtime
!  SYNOPSIS

real(kind=dp) function cgca_benchtime()

!  SOURCE

  real(kind=dp) :: dummy

! Ensure clock is initialised  

  if (firstcall) dummy = benchtick()

  call system_clock(count)

  cgca_benchtime  = real(count, kind=dp)*ticktime

end function cgca_benchtime

!*roboend*

end module cgca_m1clock
