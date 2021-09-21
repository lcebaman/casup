PROGRAM cleavage1

implicit none

INTEGER PGOPEN

external :: pgopen,pgclose

 IF (PGOPEN('cleavage1.ps/ps') .LE. 0) STOP

 call pgsch(2.0)
 call pgslw(3)
! CALL PGEX0
 CALL PGEX1
 CALL PGCLOS

END

!****************************************************************

SUBROUTINE PGEX0

implicit none

external :: pgqinf,pgqvsz
CHARACTER*64 VALUE
INTEGER LENGTH
REAL X, Y, X1, X2, Y1, Y2

!
! Information available from PGQINF:
!
      CALL PGQINF('version',  VALUE, LENGTH)
      WRITE (*,*) 'version=', VALUE(:LENGTH)
      CALL PGQINF('state',    VALUE, LENGTH)
      WRITE (*,*) 'state=',   VALUE(:LENGTH)
      CALL PGQINF('user',     VALUE, LENGTH)
      WRITE (*,*) 'user=',    VALUE(:LENGTH)
      CALL PGQINF('now',      VALUE, LENGTH)
      WRITE (*,*) 'now=',     VALUE(:LENGTH)
      CALL PGQINF('device',   VALUE, LENGTH)
      WRITE (*,*) 'device=',  VALUE(:LENGTH)
      CALL PGQINF('file',     VALUE, LENGTH)
      WRITE (*,*) 'file=',    VALUE(:LENGTH)
      CALL PGQINF('type',     VALUE, LENGTH)
      WRITE (*,*) 'type=',    VALUE(:LENGTH)
      CALL PGQINF('dev/type', VALUE, LENGTH)
      WRITE (*,*) 'dev/type=',VALUE(:LENGTH)
      CALL PGQINF('hardcopy', VALUE, LENGTH)
      WRITE (*,*) 'hardcopy=',VALUE(:LENGTH)
      CALL PGQINF('terminal', VALUE, LENGTH)
      WRITE (*,*) 'terminal=',VALUE(:LENGTH)
      CALL PGQINF('cursor',   VALUE, LENGTH)
      WRITE (*,*) 'cursor=',  VALUE(:LENGTH)
!
! Get view surface dimensions:
!
      CALL PGQVSZ(1, X1, X2, Y1, Y2)
      X = X2-X1
      Y = Y2-Y1
      WRITE (*,100) X, Y, X*25.4, Y*25.4
  100 FORMAT (' Plot dimensions (x,y; inches): ',F9.2,', ',F9.2, &
       '                          (mm): ',F9.2,', ',F9.2)

END

!****************************************************************

SUBROUTINE PGEX1

implicit none

external :: pgenv,pglab,pgline

INTEGER I,j
REAL XS(5),YS(5), XR(100), YR(100),ang,dang,pi,pi2,a(2),b(2),&
  xnew(100),ynew(100)

      pi=4.0*atan(1.0)
      pi2=pi/2.0
      dang=pi/180.0
      a=(/1.0,0.0/)
      b=0.0

!
! Call PGENV to specify the range of the axes and to draw a box, and
! PGLAB to label it. The x-axis runs from 0 to 10, and y from 0 to 20.
!
      CALL PGENV(0.0, 1.0, 0.0, 1.0, 0, 1)
      CALL PGLAB('abs(dot_product), (z)','(-atan((z-0.3)*10*N)+pi/2)/pi', &
                 & 'Prob. of cleavage')

 do j=1,10
      ang=0.0

      DO i=1,91
          b(1)= cos(ang)
          b(2)= sin(ang)
          XR(I)= dot_product(a,b)
          YR(I)=(-atan((abs(xr(i))-0.3)*j*10)+pi2)/pi
          ang=ang+dang
      end do

      CALL PGLINE(91,XR,YR)

 end do

END
