!*robodoc*u* tests/future_ca_omp1
!  NAME
!    future_ca_omp1
!  SYNOPSIS

!$Id: future_ca_omp1.f90 520 2018-03-13 18:02:06Z mexas $

program future_ca_omp1

!  PURPOSE
!    The future* tests are not part of CASU. These are to test
!    emerging capabilities. This test checks coarrays inside
!    an OpenMP parallel region.
!  DESCRIPTION
!    Run on 2 images only! This is just for demo purposes.
!    A 1D integer array coarray is set to 0 on both images.
!    The last element on image 2 is set to 1.
!    The kernel copies the value to the right to itself.
!    So gradually all values change to 1. The HX is implemented
!    using sync images.
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  USES
!  USED BY
!  SOURCE

use :: omp_lib
implicit none
integer, parameter :: n=20
!integer, external :: omp_get_num_threads, omp_get_thread_num
integer :: a(0:n+1)[*], b(0:n+1), i, img, iter, tmp, tid, nthr, nimgs

img = this_image()
nimgs = num_images()

if ( nimgs .ne. 2 ) then
  write (*,*) "ERROR: this demo program runs only on 2 images"
  error stop
end if

! Set b=0 on both images, except b(n+1)=1 on image 2
if (img .eq. 1 ) b = 0
if (img .eq. 2 ) then
  b = 0
  b(n+1) = 1
end if

! 2*n iterations are required to propagate 1 across both
! images.
main: do iter = 1, 2*n
  a = b
  !$omp parallel do default(none) private(i,tmp,tid) &
  !$omp shared(img,a,b,nthr)
  loop: do i=1, n
    nthr = omp_get_num_threads()
    if (img .eq. 1 .and. i .eq. n ) then
      tid = omp_get_thread_num()
      write (*,"(a,3(i0,tr1))") "img, nthr, tid: ", img, nthr, tid
      sync images (2)
      tmp = a(1)  [2]
      sync images (2)
      a(n+1) = tmp
    end if
    if (img .eq. 2 .and. i .eq. 1 ) then
      tid = omp_get_thread_num()
      write (*,"(a,3(i0,tr1))") "img, nthr, tid: ", img, nthr, tid
      sync images  (1)
      tmp = a(n) [1]
      sync images  (1)
      a(0) = tmp
    end if
    b(i) = max( a(i+1), a(i-1) )
  end do loop
  !$omp end parallel do
  write (*,"(a,i0,tr1,i0,tr1,999i1)") "iter, img, b: ", iter, img, b(1:n)
end do main

end program future_ca_omp1

!*roboend*
