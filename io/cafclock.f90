!  COPYRIGHT
!   Note that this routine has special Copyright conditions.
!
!    !----------------------------------------------------------------------------!
!    !                                                                            !
!    !  MPI-IO routine for Fortran Coarrays                                       !
!    !                                                                            !
!    !  David Henty, EPCC; d.henty@epcc.ed.ac.uk                                  !
!    !                                                                            !
!    !  Copyright 2013 the University of Edinburgh                                !
!    !                                                                            !
!    !  Licensed under the Apache License, Version 2.0 (the "License");           !
!    !  you may not use this file except in compliance with the License.          !
!    !  You may obtain a copy of the License at                                   !
!    !                                                                            !
!    !      http://www.apache.org/licenses/LICENSE-2.0                            !
!    !                                                                            !
!    !  Unless required by applicable law or agreed to in writing, software       !
!    !  distributed under the License is distributed on an "AS IS" BASIS,         !
!    !  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  !
!    !  See the License for the specific language governing permissions and       !
!    !  limitations under the License.                                            !
!    !                                                                            !
!    !----------------------------------------------------------------------------!

module cafclock
  implicit none

  logical,          save, private :: firstcall = .true.
  double precision, save, private :: ticktime = 0.0
  integer, parameter :: int32kind = selected_int_kind( 9)
  integer, parameter :: int64kind = selected_int_kind(18)

private
public :: caftime, caftick

!
!  Select high resolution clock
!

  integer, parameter :: intkind = int64kind
  integer(kind = intkind) :: count,rate

contains

double precision function caftime()
  double precision :: dummy
  ! Ensure clock is initialised  
  if (firstcall) dummy = caftick()
  call system_clock(count)
  caftime  = dble(count)*ticktime
end function caftime


double precision function caftick()
  if (firstcall) then
     firstcall = .false.
     call system_clock(count, rate)
     ticktime = 1.0d0/dble(rate)
  end if
  caftick = ticktime
end function caftick

end module cafclock
