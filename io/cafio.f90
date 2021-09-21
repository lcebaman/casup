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

program cafio
  use cafclock
  use cafmpiio
  use cafserialio
  implicit none

  integer, parameter :: ndim = 3
  integer, parameter :: n1 = 128, n2=n1, n3=n1
  integer, parameter :: p1 = 2,  p2=p1, p3=p1
  integer, parameter :: l1 = n1*p1
  integer, parameter :: l2 = n2*p2
  integer, parameter :: l3 = n3*p3
  integer, parameter :: numimage = p1*p2*p3
  integer, parameter :: iounit = 12
  integer, parameter :: mb = 1024*1024

  integer :: i1, i2, i3, j1, j2, j3, myimage, mypos(ndim)
  integer, allocatable :: ca(:,:,:)[:,:,:]

  double precision :: t0, t1, time, iorate, mbdata

!********************************************************************

allocate ( ca(0:n1+1, 0:n2+1, 0:n3+1) [p1, p2, *] )

  myimage = this_image()
  if (numimage /= num_images()) then
    if (myimage == 1) then
      write(*,*) 'Error: compiled for ', numimage, & 
                 ' image(s), but running on ', num_images()
    end if
    error stop
  end if

  mbdata = float(4*n1*n2*n3)*float(p1*p2*p3)/float(mb)

  if (myimage == 1) then
     write(*,*) 'Running on ', numimage, ' image(s)'
     write(*,*) 'Image grid is [', p1, ', ', p2, ', ', p3, ']'
     write(*,*) 'Array size is (', n1, ', ', n2, ', ', n3, ')'
     write(*,*)
     write(*,*) 'Total amount of data = ', mbdata, ' MB'
     write(*,*)
     write(*,*) 'Clock resolution is ', caftick()*1.0e6, ', usecs'
     write(*,*)
  end if
  
  mypos(:) = this_image(ca)

! Set halos to illegal values

  ca(:,:,:) = -1
  
! Set ca core to have unique values
 
!allocate ( ca(0:n1+1, 0:n2+1, 0:n3+1) [p1, p2, *] )
ca(1:n1,1:n2,1:n3) = this_image() 

!sync all
!  t0 = caftime()
!  call cafwritestream('stream.dat', ca)
!  t1 = caftime()
!  time = t1 - t0
!  iorate = mbdata/time
!  if (myimage == 1) &
!       write(*,*) 'time = ', time, ', stream IO rate = ', iorate, ' MB/s'

!sync all
!  t0 = caftime()
!  call cafwritedirect('direct.dat', ca)
!  t1 = caftime()
!  time = t1 - t0
!  iorate = mbdata/time
!  if (myimage == 1) &
!       write(*,*) 'time = ', time, ', direct IO rate = ', iorate, ' MB/s'

sync all

  t0 = caftime()
  call cafwrite('native.dat', ca)
  t1 = caftime()
  time = t1 - t0
  iorate = mbdata/time
  if (myimage == 1) &
       write(*,*) 'time = ', time, ', MPI-IO IO rate = ', iorate, ' MB/s'
  
end program cafio
