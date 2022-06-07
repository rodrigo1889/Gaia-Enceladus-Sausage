  module constants
  !N-Body simulation's code
    implicit none

    real(kind=8), parameter :: yr = 3.15d7, msun=2.d33
    real(kind=8), parameter :: au = 1.49d13,G=6.67d-8
    real(kind=8), parameter :: pc = 3.08d18
    real(kind=8), parameter :: r0=5.d-2*pc, dtmax=100.   !cgs
    integer, parameter      :: N=250

  end module constants
  !
  Program MainPrueba
    use constants
    Implicit None
 !   
    !   Final time is tm and the output time is calling dtp
    !
  real(kind=8), parameter         :: tm=1000000000.001, dtp=50000000. ! in yr
  real(kind=8)                    :: dt,t,tmax
  real(kind=8), dimension(N)      :: x,y,z,vx,vy,vz,m,aa
  real(kind=8), dimension(N)      :: x0,y0,z0,vx0,vy0,vz0
  real(kind=8), dimension(3,N)    :: F
  real(kind=8), dimension(3,N,N)  :: FF
  real(kind=8)     :: c,rx,ry,r3,rz,tprint,dtprint,tprintp,dtprintp
  real(kind=8)     :: um,ul,ut,e0,utt,uttt,uv
  integer  :: i,IT,itprint,ifll
  !
  itprint=0
  ifll=0
  tprint=0.
  tprintp=0.
 !
  t=0.0
  !
  ! The bodies initial conditions.
  !
  call OPENFILE(X,Y,Z,VX,VY,VZ,M)
  dt=10.       ! [yr]
  tmax=tm*yr
  dtprint=dtp*yr
  dt=dt*yr
  X=X
  Y=Y
  Z=Z
  VX=VX
  VY=VY
  VZ=VZ
  M=M
  dtprintp=dtp*yr
  call output(x,y,z,vx,vy,vz,m,t,itprint)
!
  do while (t <= tmax)
     !
     !
     call NBODY(X,Y,Z,VX,VY,VZ,M,aa,dt)
     !
     if (t .gt. 2.*dt) call TimeStep(x,y,z,vx,vy,vz,aa,dt)
      if (t+dt .gt. tmax .and. ifll == 0) then
        dt=tmax-t
        ifll=1
     endif
     write(*,*) '----> Time , dt, tmax', t/yr, dt/yr, tmax/yr
     if(t.ge.tprint) then
        write(*,*) 'output------->',itprint,t/yr

        call output(x,y,z,vx,vy,vz,m,t,itprint)
        tprint=tprint+dtprint
        itprint=itprint+1
     end if
     if(t.ge.tprintp) then
!
        tprintp=tprintp+dtprintp
     endif
  t=t+dt
  enddo

     !   
1000 format(701(' ',e12.4))
1001 format('---->',3ES12.4)
!
END PROGRAM MAINPRUEBA

!//////////////////////////////////////////////////////////////////////////////////////////
Subroutine NBODY(X,Y,Z,VX,VY,VZ,M,aa,dt)
  !
  ! The position and velocity evolution
  !
  use constants, only:N
  implicit none
!
  real(kind=8)                 :: m1,m2,m3,t
  real(kind=8), dimension(N)   :: x,v,xp,vp,y,vy,vyp,rp,vx,vxp,yp,z,zp,vzp,vz,r,m,mm,aa
  real(kind=8), dimension(N)   :: x0,y0,z0,vx0,vy0,vz0
  real(kind=8), dimension(3,N) :: F

  real(kind=8)     ::  dt, S
  real(kind=8)     :: ul,um,ut,e0,utt,uttt,uv
  integer  :: Nstep,step
  integer  :: i,j,IT
! 

     !Predictor
     call force(X,Y,Z,M,F,aa)
     do i=1,N                  
        vxp(i)=vx(i)+(f(1,i)/m(i))*dt*0.5d0
 !   where f(1,i) is the force over x axis of the  particule 1 over all others i's particles
        vyp(i)=vy(i)+(f(2,i)/m(i))*dt*0.5d0
!   where f(2,i) is the force over y axis of the  particule 1 over all others i's particles
        vzp(i)=vz(i)+(f(3,i)/m(i))*dt*0.5d0
!   where f(3,i) is the force over z axis of the  particule 1 over all others i's particles
        xp(i)=x(i)+vxp(i)*dt
        yp(i)=y(i)+vyp(i)*dt
        zp(i)=z(i)+vzp(i)*dt
     end do
     !Corrector
     call force(xp,yp,zp,M,F,aa)
     do i=1,n
        vx0(i)=vxp(i)+(f(1,i)/m(i))*dt*0.5d0
        vy0(i)=vyp(i)+(f(2,i)/m(i))*dt*0.5d0
        vz0(i)=vzp(i)+(f(3,i)/m(i))*dt*0.5d0
     end do

     x=xp
     y=yp
     z=zp
     vx=vx0
     vy=vy0
     vz=vz0

end SUBROUTINE NBODY
!////////////////////////////////////////////////////////////////////////////////////////////////// 
 subroutine force (X,Y,Z,M,F,aa)
   ! 
   use constants
  implicit none
  Integer                       :: i,j
  Real(Kind=8)                          :: rx,ry,rz,r3,mm
  Real(Kind=8), dimension(N)            :: x,y,z,m,aa
  Real(Kind=8), dimension (3,N,N)       :: FF
  Real(Kind=8), dimension (3,N)         :: F
    ! The particules are not 
  DO I=1,N
     FF(1,I,I)=0.d0
     FF(2,I,I)=0.d0
     FF(3,I,I)=0.d0 
  END DO
  !
  DO I=1,N
     DO J=1,I-1
        RX=X(I)-X(J)
        RY=Y(I)-Y(J)
        RZ=Z(I)-Z(J)
        R3=(RX**2+RY**2+RZ**2)**0.5d0   
        MM=m(I)*m(J)
        FF(1,I,J)=(-G*MM*RX)/(R3+R0)**3d0
        FF(2,I,J)=(-G*MM*RY)/(R3+R0)**3d0!  print *,x
        FF(3,I,J)=(-G*MM*RZ)/(R3+R0)**3d0
        FF(1,J,I)=-FF(1,I,J)
        FF(2,J,I)=-FF(2,I,J)
        FF(3,J,I)=-FF(3,I,J)
        !
     END DO
  END DO
  !Como condicion inicial sobre la fuerza
  DO I=1,N  
     F(1,I)=0.d0
     F(2,I)=0.d0
     F(3,I)=0.d0
     !Aqui se calcula la fuerza total
     DO J=1,N
        F(1,I)=F(1,I)+FF(1,I,J)
        F(2,I)=F(2,I)+FF(2,I,J)
        F(3,I)=F(3,I)+FF(3,I,J)
     END DO
  ENDDO
  !
  Do i=1,N
     aa(i)=sqrt(F(1,I)**2d0+F(2,I)**2d0+F(3,I)**2d0)/m(i)
  Enddo
  !
end subroutine FORCE
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////
Subroutine TimeStep (x,y,z,vx,vy,vz,aa,dt)
  use constants
  Implicit None
!  Integer, parameter     :: N=500
  Integer                :: i,j,p
  Real(Kind=8)                   :: rx,ry,rz,r3,rmin,vmax,dt,amax
  Real(Kind=8), dimension(N)     :: x,y,z,vv,vx,vy,vz,rr,aa
  Real, parameter :: A0=0.15d0
  !
  Rmin=1d50
  DO I=1,N
     DO J=1,I-1
        RX=X(I)-X(J)
        RY=Y(I)-Y(J)
        RZ=Z(I)-Z(J)
        R3=sqrt(RX**2+RY**2+RZ**2) 
        if (r3.lt.rmin) Rmin=R3
     enddo
     VV(i)=(VX(i)**2d0+VY(i)**2d0+VZ(i)**2d0)**0.5d0
 
  enddo
  !
  !
  amax=maxval(aa(2:N))
  !
!  
    rmin=max(rmin,r0)
    !
     dt=A0*sqrt(Rmin/amax)
!  
End Subroutine TimeStep
!
!////////////////////////////////////////////////////////////////////////////
SUBROUTINE OPENFILE(X,Y,Z,VX,VY,VZ,M)
  use constants, only:N
  IMPLICIT NONE
!  INTEGER,parameter           :: N=200
  REAL(KIND=8), DIMENSION(N)         :: X,Y,Z,VX,VY,VZ,M,VW,MDOT,TS
  CHARACTER (40)              :: File1
  Integer                     :: i,NN
!
  FILE1='posest.dat'
1000 FORMAT(A)
  
  OPEN(UNIT=10,FILE=FILE1,STATUS='OLD')
   read(10,*)NN
  Do i=1,N
     READ(10,*,end=2) X(i),Y(i),Z(i),VX(i),VY(i),VZ(i),M(i)
  EndDo
!
2 continue


  !
  CLOSE(UNIT=10)
END SUBROUTINE OPENFILE
!////////////////////////////////////////////////////////////////////////////
SUBROUTINE output(X,Y,Z,VX,VY,VZ,M,t,itprint)
  use constants, only : N
  IMPLICIT NONE
!
  REAL(KIND=8), DIMENSION(N)         :: X,Y,Z,VX,VY,VZ,M
  real(kind=8)                :: t,vt,sig
  CHARACTER (40)              :: File1, File2
  Integer                     :: i,NN,itprint
  !
  !
  IF(itprint.LT.1000) THEN 
     WRITE(FILE1,1001)itprint
     WRITE(FILE2,1003)itprint
  endif
  
!
  OPEN(UNIT=10,FILE=FILE1,STATUS='UNKNOWN')    
  OPEN(UNIT=11,FILE=FILE2,STATUS='UNKNOWN')
  write(10,1002)
  write(11,1002)
  Do i=1,200
     sig=-vx(i)/abs(vx(i))
     vt=sig*sqrt(vx(i)**2.+vy(i)**2.+vz(i)**2.)
     write(10,1000) X(i),Y(i),Z(i),vt   !
  EndDo
  Do i=200,250
     sig=-vx(i)/abs(vx(i))
     vt=sig*sqrt(vx(i)**2.+vy(i)**2.+vz(i)**2.)
     write(11,1000) X(i),Y(i),Z(i),vt   !
  EndDo
  
2 continue
  !
  CLOSE(UNIT=10)
1000 FORMAT(4ES14.3)
1001 FORMAT('./DAT3/outputa',I4.4,'.3D')
1003 FORMAT('./DAT3/outputb',I4.4,'.3D')
  1002 FORMAT('    X    ','    Y    ','    Z    ','    v    ')   
!
  END SUBROUTINE OUTPUT


  



  
