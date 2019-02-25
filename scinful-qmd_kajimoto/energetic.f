***************************************
*  SUBROUTINES                        *
*    alfapd                           *
*    alfapt                           *
*    dpli7                            *
*    eli8dk                           *
*                                     *
***************************************
c===========================================
c   This is file ALFAPD.FOR.
c
      Subroutine ALFAPD(Excita,V6Li)
c
c          6-Li --> d + alpha.
      Common /COLLIS/ Nelm,Echrg(6)
c
      Data Kd/10/, Ka/3/, He4pd/1.474/
c
      Ta=Excita-He4pd
      Ed=RCKE(Kd,Ka,Ta)
      Vd=VELOCITY(Kd,Ed)
      CALL RVECT(Xd,Yd,Zd)
      Vdx=Xd*Vd
      Vdy=Yd*Vd
      Vdz=Zd*Vd
      CALL LABTRAN(Vdx,Vdy,Vdz, V6Li, Vx,Vy,Vz)
      Vdl=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Echrg(3)=EFROMV(Kd,Vdl)
c       That's it for the deuteron.  Now get alpha energy.


      Eac=Ta-Ed
      Vac=VELOCITY(Ka,Eac)
      Vax=-Xd*Vac
      Vay=-Yd*Vac
      Vaz=-Zd*Vac
      CALL LABTRAN(Vax,Vay,Vaz, V6Li, Vx,Vy,Vz)
      Va=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Echrg(5)=EFROMV(Ka,Va)
      Return
c       At exit have E(deuteron) in Echrg(3) and E(alpha) in Echrg(5).


      END
c===========================================
C   This is file ALFAPT.FOR
c
      Subroutine ALFAPT(Excitn,V7Li)
c
c       Does energetics for  7-Li --> alpha + triton.


c
      Common /COLLIS/ Nelm,Echrg(6)
c
      Data Ka/3/, Kt/13/, He4pt/2.467/
c
      Ta=Excitn-He4pt
      Et=RCKE(Kt,Ka,Ta)
      CALL RVECT(Xt,Yt,Zt)
      Vtc=VELOCITY(Kt,Et)
      Vtx=Xt*Vtc
      Vty=Yt*Vtc
      Vtz=Zt*Vtc
      CALL LABTRAN(Vtx,Vty,Vtz, V7Li, Vx,Vy,Vz)
      Vt=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Echrg(6)=EFROMV(Kt,Vt)
c       Triton done.  Get alpha energy.


      Ealf=Ta-Et
      Valf=VELOCITY(Ka,Ealf)
      Vax=-Xt*Valf
      Vay=-Yt*Valf
      Vaz=-Zt*Valf
      CALL LABTRAN(Vax,Vay,Vaz, V7Li, Vx,Vy,Vz)
      Va=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Echrg(5)=EFROMV(Ka,Va)
c         and that gets the alpha energy.


      Return
      END
c===========================================
c   This is file DPLI7.FOR.
c
      Subroutine DPLI7(Excit,V9Be)
c
c            Does energetics for  9-Be --> d + 7-Li
c
      Common /COLLIS/  Nelm, Echrg(6)
      Data  Kd/10/, K7Li/17/, Dp7Li/16.6965/
c
      Ta=Excit- Dp7Li
      Ed=RCKE(Kd,K7Li,Ta)
      CALL RVECT(Xd,Yd,Zd)
      Vdc=VELOCITY(Kd,Ed)
      Vdx=Xd*Vdc
      Vdy=Yd*Vdc
      Vdz=Zd*Vdc
      CALL LABTRAN(Vdx,Vdy,Vdz, V9Be, Vx,Vy,Vz)
      Vd=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Echrg(4) = -EFROMV(Kd,Vd)
c       Negative value of E(deuteron) to alert Subroutine BANKER.
      E7Li=Ta-Ed
      V7c=VELOCITY(K7Li,E7Li)
      V7x=-Xd*V7c
      V7y=-Yd*V7c
      V7z=-Zd*V7c
      CALL LABTRAN(V7x,V7y,V7z, V9Be, Vx,Vy,Vz)
      V7=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Echrg(3)=EFROMV(K7Li,V7)
c       That gets the 7-Li ion energy.
      Return
      END
c===========================================
C   This is file ELI8DK.FOR
c
      Subroutine ELI8DK(Ex,V8Li,Eg)
c
c       Particle decay of highly excited 8-Li ion.


c
c       Decay modes considered in this subroutine are:


c         (a)   8-Li --> n + 7-Li; then,


c         (b)   7-Li --> gamma + 7-Li (ground state);


c      or (b')  7-Li --> triton + alpha;
c      or (b")  7-Li --> n + 6-Li; then,
c         (c")  6-Li --> gamma + 6-Li (ground state);


c      or (c'") 6-Li --> d + alpha.
c
      Real*4 Li7pn, Li6pn
      Common /COLLIS/ Nelm,Echrg(6)
      Common /MASSES/ Emass(21)
      Common /RANDM/  Irx
      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
c               At entry, Xn,Yn,Zn are direction cosines of the moving


c                 8-Li ion in "detector" laboratory coordinates.


      Common /NEUTRN/ Eneut, Dc1,Dc2,Dc3, Xneut,Yneut,Zneut
      Common /NEUTR2/ Eneut2, U2,V2,W2, X2,Y2,Z2
      Common /EXC7LI/ E7Li(6),Li6pn
      Common /EXC6LI/ ExLi6(6),He4pd
c          Array ExLi6 and variable He4pd given in Subroutine TENBDK


c
      Data Kn/1/, K7Li/17/, K6Li/18/
      Data Li7pn/2.033/, Li6pn/7.251/
c         Next array represents level structure of 7-Li:


      Data E7Li/0.0, 0.478, 4.63, 6.68, 7.456, 9.7/
c         For 7-Li ground state is stable, 0.478-MeV level decays by


c         gamma emission, 4.63- and 6.68-MeV levels decay by triton +


c         alpha, and higher-lying levels decay (primarily) by neutron


c         emission to levels in 6-Li.


c
      Eg=0.0
      Ta=Ex
c       First decay (a) 8-Li --> n + 7-Li:


      Taex=Ta-Li7pn
      IF (Taex .LE. 0.0) Return
      Enmax=RCKE(Kn,K7Li,Taex)
      Ehat=Enmax+5.0
c         (No magic to 5.0; it's just a value chosen to somewhat lessen


c          the "chosen" 8-Li excitation energies, on the average.)


      Fn=0.065 + 0.001*Ehat
      Temp=Fn*Ehat
      Eloww=0.0
c       Now check for existing 2d neutron waiting to be processed.  If


c         there is one, restrict present decay to neutron-bound levels


c         of 7-Li.


      IF (Eneut2 .LE. 0.0) goto 5
      Tryy=Enmax-Li6pn
      IF (Tryy .GT. 0.0) Eloww=Tryy
    5 Enn=CHOOSN(Eloww,Enmax,Temp)
      Try=Enn*(Emass(Kn) + Emass(K7Li))/Emass(K7Li)
      Entry=RCKE(Kn,K7Li,Try)
      Try=Try*Enn/Entry
      Excitn=Taex-Try
      Level=7
c       Possibly pair -EXCITN- with level in 7-Li.


      IF (Excitn .GT. 11.0) goto 10
      K=1
      Do 8 I=2,6
      IF (Excitn .LT. E7Li(I)) goto 9
    8 K=K+1
    9 Excitn=E7Li(K)
      Level=K
   10 Taex=Taex-Excitn
      Enn=RCKE(Kn,K7Li,Taex)
      CALL RVECT(Xp,Yp,Zp)
      Vnn=VELOCITY(Kn,Enn)
      Vnx=Xp*Vnn
      Vny=Yp*Vnn
      Vnz=Zp*Vnn
      CALL LABTRAN(Vnx,Vny,Vnz, V8Li, Vx,Vy,Vz)
      Vneut=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Eneut=EFROMV(Kn,Vneut)
c       Rotate neutron coordinates to "detector" system.


      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
c       Done.  Store information in NEUTRN common area.


      Dc1=Xn
      Dc2=Yn
      Dc3=Zn
c       Now get 7-Li ion information.


      E7Lic=Taex-Enn
      V7Lic=VELOCITY(K7Li,E7Lic)
      V7x=-Xp*V7Lic
      V7y=-Yp*V7Lic
      V7z=-Zp*V7Lic
      CALL LABTRAN(V7x,V7y,V7z, V8Li, Vx,Vy,Vz)
      V7Li=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
c       Test for 7-Li decay mode (b') or (b"):


      IF (Level .GE. 3) goto 15
c       Program counter to here means ground state or 0.478-MeV level


c         so tidy up 7-Li.


      Eg=Excitn
      Echrg(3)=EFROMV(K7Li,V7Li)
      Return
c
c       Next, decide between decay modes (b') and (b").


   15 IF (Level .GE. 5) goto 20
c
c       Next section for  7-Li --> alpha + triton.


      CALL ALFAPT(Excitn,V7Li)
      Return
c       Ends section on  7-Li --> alpha + triton.


c
c       Next is for second neutron emission to levels in 6-Li.


c       First transform motion of 7-Li ion to detector coordinates.


   20 CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
c       Done -- now have 7-Li ion motion in "detector" coordinates.
      Ta=Excitn-Li6pn
      Enmax=RCKE(Kn,K6Li,Ta)
      Ehat=Enmax+5.0
      Fn=0.065 + 0.001*Ehat
      Temp=Fn*Ehat
      Eloww=0.0
      Enn=CHOOSN(Eloww,Enmax,Temp)
      Try=Enn*(Emass(Kn) + Emass(K6Li))/Emass(K6Li)
      Entry=RCKE(Kn,K6Li,Try)
      Try=Try*Enn/Entry
      Excitn2=Ta-Try
c       Maybe pair -EXCITN2- with energy of a level in 6-Li.


      Level=7
      IF (Excitn2 .GT. 20.0) goto 25
      K=1
      Do 22 I=2,6
      IF (Excitn2 .LT. ExLi6(I)) goto 23
   22 K=K+1
   23 Excitn2=ExLi6(K)
      Level=K
      IF (K. EQ. 3) Eg=Excitn2
c
   25 Taex=Ta-Excitn2
      En2c=RCKE(Kn,K6Li,Taex)
      Vn2c=VELOCITY(Kn,En2c)
      CALL RVECT(Xp,Yp,Zp)
      Vn2x=Xp*Vn2c
      Vn2y=Yp*Vn2c
      Vn2z=Zp*Vn2c
      CALL LABTRAN(Vn2x,Vn2y,Vn2z, V7Li, Vx,Vy,Vz)
      Vn2=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      En2=EFROMV(Kn,Vn2)
      Eneut2=-En2
c       Negative value to alert Subroutine BANKER; will be corrected


c         at exit of calling routine.


c       Now rotate neutron coordinates into "detector" system.


      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
c       Done.  Save information in NEUTR2 common area.


      U2=Xn
      V2=Yn
      W2=Zn
      X2=Xneut
      Y2=Yneut
      Z2=Zneut
c       Okay.  Now get information on the 6-Li ion.


      E6Li=Taex-En2c
      V6Lic=VELOCITY(K6Li,E6Li)
      V6x=-Xp*V6Lic
      V6y=-Yp*V6Lic
      V6z=-Zp*V6Lic
      CALL LABTRAN(V6x,V6y,V6z, V7Li, Vx,Vy,Vz)
      V6Li=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
c       See subroutine TENBDK for decay characteristics of levels of 6-Li.


      IF (Level .EQ. 2) goto 30
      IF (Level .GE. 4) goto 30
c       Program counter to here for Level=1 or 3.  Tidy up 6-Li.


      Echrg(3)=EFROMV(K6Li,V6Li)
      Return
c
c       Last section: 6-Li --> d + alpha.


   30 CALL ALFAPD(Excitn2,V6Li)
      Return
c
c   At exit from this routine the -Echrg- array will have information on
c       decay of 8-Li as follows:


c
c   Final Decay Mode    |   Echrg(3)    |   Echrg(5)    |   Echrg(6)


c   ------------------- | ------------- | ------------- | -------------
c       n + 7-Li        |    E(7-Li)    |     0.0       |     0.0


c      2n + 6-Li        |    E(6-Li)    |     0.0       |     0.0


c     n + t + alpha     |      0.0      |   E(alpha)    |    E(t)


c    2n + d + alpha     |     E(d)      |   E(alpha)    |     0.0


c
c      Echrg(4) has energy of alpha from 12-B --> alpha + 8-Li, or
c                   energy of proton from 9-Be --> p + 8-Li.


c
      END
c===========================================
