*********************************************
*  SUBROUTINES                              *
*     cmtran                                *
*     dircos                                *
*     efromv                                *
*     exterp                                *
*     exterp2 (30/07/2001)                  *
*     labtran                               *
*     plngth                                *
*     rcke                                  *
*     rvect                                 *
*     transvec                              *
*     velocity                              *
*     (last modified by satoh on 17/08/2000)*
*                                           *
*********************************************
c==========================================================
C   This is file CMTRAN.FOR
c
      Subroutine CMTRAN(N1,N2,V, Vcom,E1,E2)
c
      Real*8 B,Bcom,E,Em1,Em2,D1,Etot,G,Gem1,P1
      Common /MASSES/ Emass(21)
      Data D1/1.0D+0/, C/2.997925E+10/
c
      Em1=Emass(N1)
      Em2=Emass(N2)
c
c -----------------------------------------------------------
c Initialize. (added by d.satoh on 2005.01.07 at JAERI)
      Vcom = 0.
      E1   = 0.
      E2   = 0.
c -----------------------------------------------------------
c check by d.satoh (2005.01.06@jaeri)
c      write(17,*) 'N1= ',N1
c      write(17,*) 'N2= ',N2
c      write(17,*) 'V = ',V
c      write(17,*) 'Vc= ',Vcom
c      write(17,*) 'E1= ',E1
c      write(17,*) 'E2= ',E2
c      write(17,*) 'Em1=',Em1
c      write(17,*) 'Em2=',Em2
c -----------------------------------------------------------
c       Let particle of mass Em1 approach particle of mass Em2 with
c        velocity V.  Particle of mass Em2 is at rest.  Then compute
c        velocity of the center of mass of the Em1-Em2 system with respect
c        to the mass-Em2-at-rest system, and compute the relativistic
c        kinetic energies, E1 and E2, of the two particles (respectively)
c        in the center-of-mass system.
c
      Beta=V/C
      B=Beta
      G=D1/DSQRT(D1-B*B)
      Gem1=G*Em1
      Etot=Gem1+Em2
      P1=Gem1*B
      Bcom=P1/Etot
      Beta=Bcom
      Vcom=Beta*C
      G=D1/DSQRT(D1-Bcom*Bcom)
      E=G*Em2 - Em2
      E2=E
      E=G*(Gem1 - Bcom*P1) - Em1
      E1=E
c -----------------------------------------------------------
c check by d.satoh (2005.01.06@jaeri)
c      write(17,*) 'N1= ',N1
c      write(17,*) 'N2= ',N2
c      write(17,*) 'V = ',V
c      write(17,*) 'Vc= ',Vcom
c      write(17,*) 'E1= ',E1
c      write(17,*) 'E2= ',E2
c      write(17,*) 'Em1=',Em1
c      write(17,*) 'Em2=',Em2
c -----------------------------------------------------------
      Return
      END
c==========================================================
c    This is file DIRCOS.FOR
c
c      Purpose is to return with direction cosines, C1,C2,C3, for input
c       velocities V1,V2,V3.


c
      SUBROUTINE DIRCOS(V1,V2,V3, C1,C2,C3)
c
      Real*8 C,Csq,X,Y,Z
      Vsq=V1*V1 + V2*V2 + V3*V3
      V=SQRT(Vsq)
    1 C1=V1/V
      C2=V2/V
      C3=V3/V
      X=C1
      Y=C2
      Z=C3
      Csq=X*X + Y*Y + Z*Z
      IF (Csq .LE. 1.0D+0) Return
      C=DSQRT(Csq)
      CC=C
      IF (CC .LT. 1.00001) CC=1.00001
      C1=C1/CC
      C2=C2/CC
      C3=C3/CC
      Return
      END
c==========================================================
c   This is file EFROMV.FOR
c
      Function EFROMV(N,V)
c       Purpose is to compute particle energy (in MeV) from input velocity V
c        (in cm/sec)


c       N = same as in function VELOCITY


c
      Real*8 B,G,D1
cs      Common /MASSES/ Emass(19)
      Common /MASSES/ Emass(21)
      Data D1/1.0D+0/,C/2.997925E+10/
c
      E=0.0
cs      IF (N.GE.1 .AND. N.LE.19) goto 1
      IF (N.GE.1 .AND. N.LE.20) goto 1
csatoh
c     Type 4, N,V
      write(*,4)N,V
    4 Format(/'  *** Error in Function EFROMV; N,V = 'I12,1PE10.3/)
      Goto 7
    1 M=N
      Beta=V/C
      B=Beta
      G=D1/DSQRT(D1 - B*B)
      G=G-D1
      Gc=G
      E=Emass(M)*Gc
c_satoh------------------------------
c      write(18,*)'mass=',Emass(M),M
c------------------------------------
    7 Efromv=E
      Return
      End


c==========================================================
c    This is file EXTERP.FOR
c
c      Purpose is to provide interpolation between two points in a table
c      of independent data, A, to get a value from a table of dependent
c      data, B, for an input value C(A).
c          N=dimension of tables A and B
c          Ntype variable:                 FOR:     A-array   B-array
c                                       Ntype =1    Linear    Linear
c                                             =2    Linear    Log
c                                             =3    Log       Linear
c                                             =4    Log       Log
c
      FUNCTION EXTERP(A,B,C,N,Ntype)
      Real*8 Ac,Ae,Af,Amin,R,Xb,Xc,Y,Z
      Dimension A(1),B(1)
      Data Amin/3.72008D-24/
c
      Exterp=B(1)
      If (C.LE.A(1)) Return
      Exterp=B(N)
      If (C.GE.A(N)) Return
c
c       First 4 statements take care of case where input value of C was
c         outside of the limits of the tabular values of the A array.
c
      Do 2 I=2,N
      IF (C - A(I)) 3,9,2
    9 Exterp=B(I)
      Return
    2 Continue
    3 J=I
      I=J-1
c
      L=Ntype
      Ac=C
      Ae=A(I)
      Af=A(J)
      Xb=B(I)
      Xc=B(J)
      If (L.LE.2) goto 4
      IF (Ac .LT. Amin) Ac=Amin
      Ac=Dlog(Ac)
      IF (Ae .LT. Amin) Ae=Amin
      Ae=Dlog(Ae)
      IF (Af .LT. Amin) Af=Amin
      Af=Dlog(Af)
    4 IF (L.LE.1 .OR. L.EQ.3) goto 7
      IF (Xb .LT. Amin) Xb=Amin
      Xb=Dlog(Xb)
      IF (Xc .LT. Amin) Xc=Amin
      Xc=Dlog(Xc)
    7 R=(Ac-Ae)/(Af-Ae)
      Y=Xb+R*(Xc-Xb)
      Z=Y
      IF (L.EQ.2 .OR. L.EQ.4) Z=Dexp(Y)
      Exterp=Z
      Return
      END
c==========================================================


      FUNCTION EXTERP2(A,B,C,N,Ntype)
      Real*8 Ac,Ae,Af,Amin,R,Xb,Xc,Y,Z
c_modified by d.satoh
      Real*8 dc,x1,x2,y1,y2,S,dy,By
c_satoh_end
      Dimension A(1),B(1)
      Data Amin/3.72008D-24/
c
      Exterp2=B(1)
      If (C.LE.A(1)) Return
      Exterp2=B(N)
c_modified by d.satoh
c_y-axis is log scale!!
      If (C.GE.A(N)) then


         dc = C


         x1 = A(N)
         y1 = B(N)
         y1 = dlog(y1)


         x2 = A(N-1)
         y2 = B(N-1)
         y2 = dlog(y2)


         S = ( y1 - y2 )/( x1 - x2 )
         dy = S*( dc-x1 )
         By = y1 + dy
         Exterp2 = exp(By)
         return
      end if
c_satoh_end
c
c       First 4 statements take care of case where input value of C was
c         outside of the limits of the tabular values of the A array.
c
      Do 2 I=2,N
      IF (C - A(I)) 3,9,2
    9 Exterp2=B(I)
      Return
    2 Continue
    3 J=I
      I=J-1
c
      L=Ntype
      Ac=C
      Ae=A(I)
      Af=A(J)
      Xb=B(I)
      Xc=B(J)
      If (L.LE.2) goto 4
      IF (Ac .LT. Amin) Ac=Amin
      Ac=Dlog(Ac)
      IF (Ae .LT. Amin) Ae=Amin
      Ae=Dlog(Ae)
      IF (Af .LT. Amin) Af=Amin
      Af=Dlog(Af)
    4 IF (L.LE.1 .OR. L.EQ.3) goto 7
      IF (Xb .LT. Amin) Xb=Amin
      Xb=Dlog(Xb)
      IF (Xc .LT. Amin) Xc=Amin
      Xc=Dlog(Xc)
    7 R=(Ac-Ae)/(Af-Ae)
      Y=Xb+R*(Xc-Xb)
      Z=Y
      IF (L.EQ.2 .OR. L.EQ.4) Z=Dexp(Y)
      Exterp2=Z
      Return
      END
c==========================================================


C   This is file LABTRAN.FOR
c
      Subroutine LABTRAN(Vxp,Vyp,Vzp, V, Vx,Vy,Vz)
c
c       Given velocity components Vxp, Vyp, and Vzp in the "primed"


c        coordinate system which is moving with velocity V along the


c        Z-axis of the laboratory system, compute the relativistically


c        correct velocity components Vx, Vy, and Vz of the velocity


c        in the laboratory system.


c
      Real*8 Vxx,Vyy,Vzz,Vv,Beta,C,Denom,D1,Sqr
      Data D1/1.0D+0/, C/2.99792458D+10/
c
      Vxx=Vxp
      Vyy=Vyp
      Vzz=Vzp
      Vv=V
      Beta=Vv/C
      Denom=D1+Beta*Vzz/C
      Sqr=DSQRT(D1-Beta*Beta)
      Vxx=Vxx*Sqr/Denom
      Vyy=Vyy*Sqr/Denom
      Vzz=(Vzz+Vv)/Denom
      Vx=Vxx
      Vy=Vyy
      Vz=Vzz
      Return
      END
c==========================================================
c   This is file PLNGTH.FOR
c
c   April 1986 ---
c       Purpose is to compute length of a vector with direction cosines
c       AX,AY,AZ from origin X,Y,Z to a surface of a cylinder of outer radius
c       R and length H.  Cartesian coordinates are used wherein the X-Y plane
c       is the 'front' surface of the cylinder, and the Z-axis is the
c       cylinder axis.
c
      Function PLNGTH(Ax,Ay,Az,  X,Y,Z,  R,H)
c
        IF (Az) 1,2,3
   1    Akzet=-Z/Az
        Goto 4


   2    Akzet=10000.
C           (It is assumed that 10000 is >>> cylinder dimensions.)


        IF (Z .LT. 0.0) Akzet=Z


        Goto 4


   3    Akzet=(H-Z)/Az
   4    Pl=Akzet
        Arsq= Ax*Ax+Ay*Ay


        IF (Arsq .LE. 1.0E-10) goto 5


        C=Ay*X-Ax*Y


        Argu=Arsq*R*R-C*C


        IF (Argu .GE. 0.0) goto 6


csatoh
c       Type 10, Argu,Arsq,C
        write(*,10)Argu,Arsq,C
  10    Format(//20x,'Problem in -PLNGTH-  Arguments='1PE10.3/
     x   25x,'ARSQ='E10.3,', C='E10.3,'.  Set Argument = 0.0')
csatoh
c       Type 11, X,Y,Z, Ax,Ay,Az


        write(*,11)X,Y,Z,Ax,Ay,Az
  11    Format(10x,'    Input X, Y, Z ='1P3E10.3/
     x   8x,'      and Dir. Cos. ='1P3E10.3/)
        Argu=0.0


   6    St=SQRT(Argu)
        Akrad=(St-Ax*X-Ay*Y)/Arsq


        IF (Akrad .LT. 0.0) Akrad=0.0


c         Last is to prevent a (small) negative value -- due to single-
c           precision round-off errors -- when values of X,Y place the
c           point essentially on the (curved) surface of the cylinder.


        Pl=Akrad


        IF (Akrad .GT. Akzet) Pl=Akzet


   5    PLNGTH=Pl
        Return
        END
c==========================================================
C   This is file RCKE.FOR
c
      Function RCKE(N1,N2,TKE)
c
c   RCKE = Relativistically Correct Kinetic Energy
c       TKE is the total kinetic energy available in the center of mass


c       of the system of particles having masses EM1 and EM2.  RCKE
c       computes the kinetic energy for the particle having mass EM1.


c
      Real*8 Em1,Em2,E,Etot,Dhalf
      Common /MASSES/ Emass(21)
      Data Dhalf/0.5E+0/
c
      Em1=Emass(N1)
      Em2=Emass(N2)
      E=TKE
      Etot=E+Em1+Em2
      E=Dhalf*(Etot*Etot + Em1*Em1 - Em2*Em2)/Etot
      E=E-Em1
      RCKE=E
      Return
      END
c==========================================================
C   This is file RVECT.FOR
c
c       Purpose is to obtain a random unit vector with components


c         (i.e. direction cosines) Xn, Yn, Zn


      Subroutine RVECT(Xn,Yn,Zn)
c
        Real*8 Zet,S,D1
        Common /RANDM/ Ixx
        Data D1/1.0D+0/
c
        Z=1.0-2.0*RAN(IXX)


        ZN=Z


        Zet=Z


        THETA=6.283185*RAN(IXX)


        S=DSQRT(1.0-Zet*Zet)


        SINZ=S


        XN=COS(THETA)*SINZ


        YN=SIN(THETA)*SINZ


      Return
      END
c==========================================================
C   This is file TRANSVEC.FOR
c
      Subroutine TRANSVEC(Zx,Zy,Zz)
c
c        Given a vector with direction cosines Vxp,Vyp, and Vzp in a
c        coordinate system having its Z-axis direction cosines as
c        Zx, Zy, and Zz in a second coordinate system, determine the
c        directions cosines of the vector in the second coordinate
c        system.  Since the X and Y axes corresponding to the given
c        Z axis are not given, it is assumed that the Y axis lies
c        in the x-y plane of the second coordinate system.  The
c        resulting direction cosines of the vector in the second
c        coordinate system are Vx, Vy, and Vz.


c
      Real*8 D1,Dum1,Dum2,dcut
      Common /VECTOR/ Vxp,Vyp,Vzp, Vx,Vy,Vz
      Data D1/1.0D+0/
      data dcut/-1.0d-2/
c
      IF (Zz.GE.1.0) goto 2
      Dum1=Zz
      Dum2=DSQRT(D1-Dum1*Dum1)
      Xz=-Dum2
cs      IF (Xz .GT. -0.0004) goto 2
      IF (Xz .GT. dcut ) goto 2
      Yx=Zy/Xz
      Yy=-Zx/Xz
      Xx=Yy*Zz
      Xy=-Yx*Zz
c
      Vx=Xx*Vxp + Yx*Vyp + Zx*Vzp
      Vy=Xy*Vxp + Yy*Vyp + Zy*Vzp
      Vz=Xz*Vxp          + Zz*Vzp
      Return
c
    2 Vx=Vxp
      Vy=Vyp
      Vz=Vzp
      Return
      END
c==========================================================
      Function VELOCITY(N,E)
c
c       N=1  Particle=neutron
c         2     "    =proton
c         3     "    =alpha
c         4     "    =9-Be
c         5     "    =11-B
c         6     "    =12-B
c         7     "    =12-C
c         8     "    =8-Be
c         9     "    =11-C
c        10     "    =deuteron
c        11     "    =5-He
c        12     "    =10-B
c        13     "    =triton
c        14     "    =11-Be
c        15     "    =10-Be
c        16     "    =8-Li
c        17     "    =7-Li
c        18     "    =6-Li
c        19     "    =3-He
c        20     "    =pion+ (Added by satoh)
c        21     "    =pion- (Added by satoh)
c
      Real*8 G,Gsq,B,D1
      Common /MASSES/ Emass(21)
      Data D1/1.0D+0/,C/2.997925E+10/
      Data Emass /939.583,938.272,3728.43,8394.86,10255.2,11191.4,
     & 11178.0,7456.95,10257.17,1876.14,4668.9,9327.05,2809.45,
     & 10266.69,9327.62,7472.96,6535.42,5603.097,2809.44,139.570,
     & 139.570/


c
      V=0.0
      IF (E .LE. 0.0) goto 2
      IF (N.GE.1 .AND. N.LE.20) goto 1
    2 write(*,4)N,E
    4 Format(/'  *** Error in Function VELOCITY; N, E ='I12,1PE10.3/)
      Goto 7
    1 M=N
      Gamma=E/Emass(M)
      G=Gamma
      G=G+D1
      Gsq=G*G
      B=DSQRT(D1-D1/Gsq)
      Beta=B
      V=C*Beta
    7 Velocity=V
      Return
      END
c==========================================================


