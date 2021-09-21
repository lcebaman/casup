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

module cafserialio
  use mpi
  implicit none

!  integer, parameter, private :: n1 = 128
!  integer, parameter, private :: n2 = 128
!  integer, parameter, private :: n3 =  64
!  integer, parameter, private :: p1 =  8
!  integer, parameter, private :: p2 =  8
!  integer, parameter, private :: p3 = 16

  integer, parameter, private :: iounit = 12
  integer, parameter, private :: ndim = 3
  integer, private :: i1, i2, i3, coi1, coi2, coi3, myimage
  integer, private, dimension(ndim) :: lb, ub, lcob, ucob

private
public :: cafwritestream, cafwritedirect

contains

!*************************************************************************
subroutine cafwritestream(filename, array)

  character*(*), intent(in) :: filename
  integer, allocatable, intent(in) :: array(:,:,:)[:,:,:]

! Would like to use the definition below so the subroutine is generic for
! all array sizes, but currently this causes the Cray compiler to crash
!
! integer :: array(:,:,:)[:,:,*]

  myimage  = this_image()

  if (myimage == 1) then

     open(unit=iounit, file=filename, form='unformatted', access="stream")

     lb = lbound(array)+1
     ub = ubound(array)-1
     lcob = lcobound(array)
     ucob = ucobound(array)

     do coi3=lcob(3),ucob(3)
        do i3=lb(3),ub(3)
           do coi2=lcob(2),ucob(2)
              do i2=lb(2),ub(2)
                 do coi1=lcob(1),ucob(1)

                    write(iounit) array(lb(1):ub(1),i2,i3)[coi1,coi2,coi3]

                 end do
              end do
           end do
        end do
     end do

     close(unit=iounit)

  end if

end subroutine cafwritestream


!*************************************************************************
subroutine cafwritedirect(filename, array)

  character*(*), intent(in) :: filename
  integer, allocatable, intent(in) :: array(:,:,:)[:,:,:]

! Would like to use the definition below so the subroutine is generic for
! all array sizes, but currently this causes the Cray compiler to crash
!
! integer :: array(:,:,:)[:,:,*]

  integer, parameter :: intsize=4

  integer :: irec, reclength

  myimage  = this_image()

  if (myimage == 1) then

     lb = lbound(array)+1
     ub = ubound(array)-1
     lcob = lcobound(array)
     ucob = ucobound(array)

     reclength = (ub(1)-lb(1)+1)*intsize

     open(unit=iounit, file=filename, form='unformatted', &
          access='direct', recl=reclength)

     irec = 1

     do coi3=lcob(3),ucob(3)
        do i3=lb(3),ub(3)
           do coi2=lcob(2),ucob(2)
              do i2=lb(2),ub(2)
                 do coi1=lcob(1),ucob(1)

                    write(iounit,rec=irec) &
                         array(lb(1):ub(1),i2,i3)[coi1,coi2,coi3]

                    irec = irec + 1
                    
                 end do
              end do
           end do
        end do
     end do

     close(unit=iounit)

  end if

end subroutine cafwritedirect

end module cafserialio
