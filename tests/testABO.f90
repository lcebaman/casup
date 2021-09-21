!*robodoc*u* tests/testABO
!  NAME
!    testABO
!  SYNOPSIS

!$Id: testABO.f90 418 2017-06-07 14:02:10Z mexas $

program testABO

!  PURPOSE
!    Checking: cgca_clvgp, cgca_clvgsd, cgca_pswci, cgca_sld3
!  DESCRIPTION
!    Checking cleavage propagation across grain boundary
!    with many grains. Still a single cleavage nucleus.
!    Put it at the centre of one of the faces of the model.
!    Solidification is done with Cray reduction, CO_SUM (cgca_sld3).
!    Output is done with MPI/IO (cgca_pswci) and with
!    the serial version (cgca_swci). The timings of both
!    IO routines are done for comparison.
!  NOTES
!    The program must be called with 2 command line arguments,
!    both positive integers. These are codimensions along 1 and 2.
!    The number of images must be such that
!    codimension3 = num_images()/( codimension1 * codimension3 )
!    is a positive integer. Example:
!      ./testABO.x 2 2
!    which will make the third codimension equal to 16/(2*2)=4.
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

logical( kind=ldef ), parameter :: yesdebug=.true., nodebug=.false., &
  periodicbc=.true.,  noperiodicbc=.false.

real, parameter :: gigabyte=real(2**30), resolution=1.0e-5,          &
! cleavage stress on 100, 110, 111 planes for BCC,
! see the manual for derivation.
  scrit(3) = (/ 1.05e4, 1.25e4, 4.90e4 /)

real( kind=rdef ), allocatable :: grt(:,:,:)[:,:,:]
real( kind=rdef ) :: t(3,3)   ! stress tensor

integer( kind=idef ) :: l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3,  &
 cou3, nuc,nimages,codim(3)[*]
integer( kind=iarr ), allocatable :: space(:,:,:,:)[:,:,:]
integer( kind=ilrg ) :: icells, mcells

real :: time1, time2, fsizeb, fsizeg, tdiff

logical( kind=ldef ) :: solid, image1

!*********************************************************************72
! first executable statement

nimages=num_images()
image1=.false.
if (this_image() .eq. 1) image1 = .true.

! do a check on image 1
if ( image1 ) then
 call getcodim(nimages,codim)
 ! print a banner
 call banner("ABO")
 ! print the parameter values
 call cgca_pdmp
 write (*,'(a,i0,a)') "running on ", nimages, " images in a 3D grid"
end if

sync all

codim(:) = codim(:)[1]

l1 = 1
l2 = l1
l3 = l1

! The array size is only controlled by this value in this program.
u1 = 2**7 ! 128
u2 = u1
u3 = u1

col1 = 1
cou1 = codim(1) - col1 + 1
col2 = 1
cou2 = codim(2) - col2 + 1
col3 = 1
cou3 = codim(3) - col3 + 1

! total number of cells in a coarray
icells = int( u1-l1+1, kind=ilrg ) *    &
         int( u2-l2+1, kind=ilrg ) *    &
         int( u3-l3+1, kind=ilrg )

! total number of cells in the model
mcells = icells * int( nimages, kind=ilrg )

! total number of nuclei
! nuc should not exceed resolution*mcells
nuc = 100

! total number of nuclei
nuc = int( resolution * real( mcells ) )

if ( image1 ) then
  write (*,'(a,2(i0,":",i0,","),i0,":",i0,")")') &
    "bounds: (",l1,u1,l2,u2,l3,u3
  write (*,'(a,2(i0,":",i0,","),i0,":",i0,")")') &
    "cobounds: (",col1,cou1,col2,cou2,col3,cou3

! Total output file size, in B and in GB.
fsizeb = real( mcells*storage_size(space,kind=ilrg)/8_ilrg )
fsizeg = fsizeb / gigabyte

  write (*,'(a,i0,a)') "Each image has ",icells, " cells"
  write (*,'(a,i0,a)') "The model has ", mcells, " cells"
  write (*,'(a,i0,a)') "The model has ", nuc, " nuclei"
  write (*,'(2(a,es9.2),a)') "The output file size is ", fsizeb, &
   " B, or ", fsizeg, "GB."
end if

! initialise random number seed
call cgca_irs( nodebug )

! allocate space with two layers
call cgca_as( l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3,2,space )

! allocate rotation tensors
call cgca_art( 1, nuc, col1, cou1, col2, cou2, col3, grt )

! initialise space
space( :, :, :, cgca_state_type_grain ) = cgca_liquid_state
space( :, :, :, cgca_state_type_frac ) = cgca_intact_state

! nuclei, sync all inside
call cgca_nr( space, nuc, nodebug )

! assign rotation tensors, sync all inside
call cgca_rt( grt )

! set a single crack nucleus in the centre of the front face
space( u1/2, u2/2, u3, cgca_state_type_frac )[ cou1/2, cou2/2, cou3] = &
 cgca_clvg_state_100_edge

! solidify
! subroutine cgca_sld3(coarray,iter,heartbeat,solid)
call cgca_sld3( space, 0, 1, solid )

! update grain connectivity, local routine, no sync needed
call cgca_gcu( space )

! global sync needed to wait for cgca_gcu to complete on all images
sync all

! set the stress tensor
t = 0.0
t(1,1) = 1.0e6
t(2,2) = -1.0e6

! propagate cleavage, sync inside
! subroutine cgca_clvgp( coarray, rt, t, scrit, sub, gcus,
!    periodicbc, iter, heartbeat, debug )
call cgca_clvgp( space, grt, t, scrit, cgca_clvgsd, cgca_gcupdn,       &
     noperiodicbc, 200, 10, nodebug )

! dump the model
call cpu_time(time1)
call cgca_pswci( space, cgca_state_type_frac, 'frac_mpi.raw' )
call cpu_time(time2)
tdiff = time2-time1
if ( image1 ) write (*,*) "MPI/IO: ", tdiff, "s, rate: ", fsizeg/tdiff, "GB/s."

call cpu_time(time1)
call cgca_swci( space, cgca_state_type_frac, 10, 'frac_ser.raw' )
call cpu_time(time2)
tdiff = time2-time1
if ( image1 ) write (*,*) "Serial IO: ", tdiff, "s, rate: ", fsizeg/tdiff, "GB/s."

! However, since there's nothing more to do, no sync is needed.

! deallocate all arrays, implicit sync all.
call cgca_ds( space )
call cgca_dgc
call cgca_drt( grt )

end program testABO

!*roboend*
