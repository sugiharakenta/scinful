***************************************
*  SUBROUTINES                        *
*     ap7li                           *
*     ex8li                           *
*     exli8                           *
*     pp8li                           *
*     tenbdk                          *
*                                     *
***************************************
c==========================================
C  This is file AP7LI.FOR
c
      Subroutine AP7LI(Ta,V11B,Eg)
C  Purpose is to follow  11-B --> 7-Li + alpha  decay, including
c       any subsequent 7-Li excited ion decay.


c
      Real*4 Li6pn
      Common /RANDM/  Irx
      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
c       Xn,Yn,Zn at entry are direction cosines of 11-B ion in


c           "detector" coordinates.


      Common /MASSES/ Emass(21)
      Common /NEUTRN/ Eneut, U,V,W,X,Y,Z
      Common /NEUTR2/ Eneut2, U2,V2,W2,X2,Y2,Z2
      Common /COLLIS/ Nelm,Echrg(6)
      Common /EXC7LI/ E7Li(6),Li6pn
c       See Subroutine ELI8DK for values in EXC7LI common.


      Common /EXC6LI/ E6Li(6),He4pd
c       See Subroutine TENBDK for values in EXC6LI common.


c
      Data Ka/3/, K7Li/17/, Kn/1/, K6Li/18/
c
      Eax=RCKE(Ka,K7Li,Ta)
      Eamin=0.0
      Ndx=0
      IF (Eneut .LE. 0.0) goto 5
c       Program counter to here means call to AP7LI came from NPX subroutine.


c         Otherwise call came from subroutine ND.


      Ndx=1
c       Check to see if there is a 2d neutron waiting to be processed.  If so


c         then fix 11-B decay only to n-stable states of 7-Li.


      IF (Eneut2 .LE. 0.0) goto 5
c         Fix 11-B ion decay by limiting minimum value of energy of out-


c         going alpha particle.


      Eamin=Ta-Li6pn
    5 Temp=4.0
c       See comments just prior to statement 38 in subroutine NN3ALF for


c          estimating "alpha" continuum.


    6 Eatry=CHOOSP(Eax,Temp)
      IF (Eatry .LT. Eamin) goto 6
      Try=Eatry*(Emass(Ka) + Emass(K7Li))/Emass(K7Li)
      Etry=RCKE(Ka,K7Li,Try)
      Try=Try*Eatry/Etry
      Excit=Ta-Try
c       Now, maybe, set  -EXCIT-  to specific energy level in 7-Li.


      Level=7
      IF (Excit .GT. 11.0) goto 10
      K=1
      Do 8 I=2,6
      IF (Excit .LT. E7Li(I)) goto 9
    8 K=K+1
    9 Excit=E7Li(K)
      Level=K
   10 Taex=Ta-Excit
      Eaa=RCKE(Ka,K7Li,Taex)
      CALL RVECT(Xp,Yp,Zp)
      Vaa=VELOCITY(Ka,Eaa)
      Vax=Xp*Vaa
      Vay=Yp*Vaa
      Vaz=Zp*Vaa
      CALL LABTRAN(Vax,Vay,Vaz, V11B, Vx,Vy,Vz)
      Valf=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Ealf=EFROMV(Ka,Valf)
      Echrg(4)=Ealf
c       That gets the alpha information.


c       Now get 7-Li ion information


      E7Lic=Taex-Eaa
      V7Lic=VELOCITY(K7Li,E7Lic)
      V7x=-Xp*V7Lic
      V7y=-Yp*V7Lic
      V7z=-Zp*V7Lic
      CALL LABTRAN(V7x,V7y,V7z, V11B, Vx,Vy,Vz)
      V7Li=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      IF (Level .GE. 3) goto 15
c       Program to here have bound state of 7-Li to tidy up.


      Eg=Excit
      Echrg(3)=EFROMV(K7Li,V7Li)
      Return
c
c       Next: decide 7-Li particle decay mode.


   15 IF (Level .GE. 5) goto 20
c       If program to here have  7-Li --> alpha + triton  decay mode.


      CALL ALFAPT(Excit,V7Li)
      Return
c
c       Next is for  7-Li --> n + 6-Li  decay mode.  First get 7-Li motion


c         into detector coordinates.


   20 CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
C       Done.  Now get neutron information.


      Ta=Excit-Li6pn
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
c       Maybe pair -EXCITN2-  with energy of a level in 6-Li


      Level=7
      IF (Excitn2 .GT. 20.0) goto 25
      K=1
      Do 22 I=2,6
      IF (Excitn2 .LT. E6Li(I)) goto 23
   22 K=K+1
   23 Excitn2=E6Li(K)
      Level=K
      IF (K .EQ. 3) Eg=Excitn2
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
c       Transform neutron dir. cos. into "detector" coordinates.


      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
c       Now recall whether call to this routine was from Subroutine NPX


c          or subroutine ND.


      IF (Ndx .EQ. 1) goto 27
c       Program counter to here, Subroutine ND was calling routine, so


c         store information into /NEUTRN/ common variables.


      Eneut=En2
      U=Xn
      V=Yn
      W=Zn
      Goto 28
c       Next is when calling subroutine was NPX; store information into


c         /NEUTR2/ common variables.


   27 Eneut2=-En2
c          Negative value to alert Subroutine BANKER; it will be corrected


c            at the exit of the calling routine.


      U2=Xn
      V2=Yn
      W2=Zn
      X2=X
      Y2=Y
      Z2=Z
c       Done.  Now get 6-Li ion information.


   28 E6Lix=Taex-En2c
      V6Lic=VELOCITY(K6Li,E6Lix)
      V6x=-Xp*V6Lic
      V6y=-Yp*V6Lic
      V6z=-Zp*V6Lic
      CALL LABTRAN(V6x,V6y,V6z, V7Li, Vx,Vy,Vz)
      V6Li=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
c       See Subroutine TENBDK for decay characteristics of levels of 6-Li


      IF (Level .EQ. 2) goto 30
      IF (Level .GE. 4) goto 30
c       Program to here, bound state of 6-Li.  Tidy up.


      Echrg(3)=EFROMV(K6Li,V6Li)
      Return
c
c      Last section -- 6-Li --> d + alpha
   30 CALL ALFAPD(Excitn2,V6Li)
      Return
c       See end of subroutine ELI8DK for -Echrg- array at end of this


c         routine.


      END
c==========================================
c   This is file EX8LI.FOR
c
c     Contains routines involving "decay" of Li ions.
c
      Function EX8LI(Exc12B)
cc
c     Purpose: for alpha decay of highly-excited 12-B, determine 8-Li
c       excitation energy.  Present programming mocks up TNG results.
c
      Real*4 Li8pa
      Common /RANDM/  Irx
      Common /EXC8LI/ E8Li(10)
c
      Dimension Ex(28),P(10),P1(28),P2(28),P3(28),P4(28),P5(28),
     x    P6(28),P7(28),P8(28),P9(28),P10(28)
c
c     Next array represents level structure of 8-Li:
      Data E8Li/0.0,0.981,2.255, 3.21, 5.4, 6.1, 6.53, 7.1, 9.0, 10.82/
      Data N/28/, Nterp/1/, Temp/4.0/, Li8pa/10.002/
c
c     Next arrays relate probabilities for population of individual levels
c       vs excitation energy (Ex) in 12-B.
      Data Ex/ 11.0,  11.4,  12.3,  13.3,  14.2,  15.0,  15.5,  16.5,
     a  17.5,  18.5,  19.5,  20.0,  21.0,  22.0,  23.0,  24.0,  25.0,
     b  26.0,  27.0,  28.0,  30.0,  32.0,  34.0,  36.0,  38.0,  40.0,
     c  44.0,  48.0/
      Data P1/  1.0,  .9815, .985,  .8606, .6478, .574,  .513,  .41,
     a  .341,  .288,  .3448, .338,  .315,  .28,   .25,   .226,  .215,
     b  .2066, .2004, .186,  .1582, .141,  .1235, .1087, .0933, .089,
     c  .079,  .0665/
      Data P2/  1.0,   1.0,   1.0,  .9996, .9808, .899,  .813,  .66,
     a  .553,  .483,  .5268, .518,  .49,   .45,   .415,  .388,  .373,
     b  .3616, .3514, .333,  .3002, .275,  .2505, .2257, .1983, .179,
     c  .149,  .1115/
      Data P3/  1.0,   1.0,   1.0,   1.0,  .9998, .989,  .983,  .95,
     a  .893,  .783,  .7798, .758,  .71,   .66,   .616,  .583,  .56,
     b  .5416, .5264, .502,  .4572, .42,   .3845, .3507, .3133, .284,
     c  .235,  .1815/
      Data P4/  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  .984,
     a  .968,  .903,  .9048, .888,  .838,  .7843, .733,  .6943, .665,
     b  .6406, .6204, .591,  .5382, .494,  .4535, .4147, .3733, .338,
     c  .276,  .2115/
      Data P5/  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     a  .999,  .988,  .9798, .958,  .918,  .8773, .838,  .8093, .786,
     b  .7656, .7454, .715,  .6532, .595,  .5405, .4977, .4433, .4,
     c  .322,  .2445/
      Data P6/  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     a   1.0,  .994,  .9898, .973,  .947,  .9233, .903,  .8833, .866,
     b  .8486, .8284, .796,  .7294, .666,  .6075, .5557, .5033, .455,
     c  .366,  .2785/
      Data P7/  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     a   1.0,   1.0,  .9998, .988,  .97,   .9553, .943,  .9313, .921,
     b  .9086, .8924, .863,  .7954, .73,   .6675, .6127, .5583, .507,
     c  .414,  .3235/
      Data P8/  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     a   1.0,   1.0,   1.0,  .998,  .99,   .9843, .978,  .9723, .9665,
     b  .9576, .9424, .913,  .8444, .777,  .7115, .6527, .5953, .542,
     c  .4433, .3475/
      Data P9/  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     a   1.0,   1.0,   1.0,   1.0,   1.0,  .9998, .9985, .9973, .9947,
     b  .9886, .9754, .9472, .8804, .8152, .7495, .6897, .6303, .5762,
     c  .4733, .3755/
      Data P10/ 1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     a   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  .9999, .9995, .998,
     b  .993,  .981,  .954,  .889,  .825,  .76,   .7,    .64,   .585,
     c  .48,   .38/
c
      E=Exc12B
      P(1) =EXTERP(Ex, P1,E,N,Nterp)
      P(2) =EXTERP(Ex, P2,E,N,Nterp)
      P(3) =EXTERP(Ex, P3,E,N,Nterp)
      P(4) =EXTERP(Ex, P4,E,N,Nterp)
      P(5) =EXTERP(Ex, P5,E,N,Nterp)
      P(6) =EXTERP(Ex, P6,E,N,Nterp)
      P(7) =EXTERP(Ex, P7,E,N,Nterp)
      P(8) =EXTERP(Ex, P8,E,N,Nterp)
      P(9) =EXTERP(Ex, P9,E,N,Nterp)
      P(10)=EXTERP(Ex,P10,E,N,Nterp)
c
      Pran=RAN(Irx)
      Nb=1
      Do 1 J=1,10
      IF (Pran .LE. P(J)) goto 2
    1 Nb=Nb+1
    2 IF (Nb .EQ. 11) goto 3
      Ex8Li=E8Li(Nb)
      Return
c
c     Program to here--highly excited 8-Li, E(level) from "continuum"
    3 Ex8Li=E8Li(10)
      Emax=E-22.5
      IF (Emax .LE. 0.0) Return
      Eatry=CHOOSP(Emax,Temp)
      Ex8Li=E-Eatry-Li8pa
      Return
      END
c==========================================
C   This is file EXLI8.FOR
c
      Function EXLI8(Exc9Be)
cc
c       Purpose: for proton decay of highly-excited 9-Be, determine 8-Li


c         excitation energy...program mocks up TNG results.


c
      Real*4 Li8pp
      Common /RANDM/  Irx
      Common /EXC8LI/ E8Li(10)
c           E8Li array in Function EX8LI


c
      Dimension Ex(28),P(10),P1(28),P2(28),P3(28),P4(28),P5(28),
     x     P6(28),P7(28),P8(28),P9(28),P10(28)
      Data N/28/, Nterp/1/, Li8pp/16.888/
c
c       Next arrays relate probabilities for population of individual


c         levels vs excitation energy (Ex) in 9-Be.


      Data Ex/ 17.9,  18.25, 19.2,  19.6,  20.15, 20.6,  22.35, 22.7,
     a  23.05, 23.45, 24.15, 25.0,  25.9,  26.8,  27.7,  28.7,  30.0,
     b  32.0,  34.0,  36.0,  28.0,  40.0,  42.0,  44.0,  46.0,  48.0,
     c  50.0,  55.0/
      Data P1/  1.0,  .99,   .67,   .629,  .585,  .551,  .514,  .497,
     a  .417,  .388,  .3367, .323,  .31,   .3042, .3042, .2822, .28,
     b  .236,  .195,  .151,  .113,  .08,   .056,  .038,  .025,  .0165,
     c  .0106, .003/
      Data P2/  1.0,   1.0,   1.0,  .989,  .97,   .936,  .834,  .795,
     a  .697,  .655,  .5667, .525,  .483,  .4472, .4282, .3842, .358,
     b  .288,  .2295, .273,  .127,  .089,  .0616, .0415, .0271, .0177,
     c  .0113, .003/
      Data P3/  1.0,   1.0,   1.0,   1.0,   1.0,  .991,  .93,   .908,
     a  .885,  .857,  .7847, .745,  .697,  .6472, .6102, .5442, .493,
     b  .386,  .2985, .22,   .158,  .109,  .0743, .0492, .0318, .0204,
     c  .0128, .003/
      Data P4/  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  .998,
     b  .995,  .989,  .9687, .942,  .894,  .8352, .7702, .6962, .623,
     c  .484,  .3695, .27,   .193,  .1325, .0909, .0597, .0388, .0248,
     d  .0156, .0035/
      Data P5/  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     a   1.0,  .9976, .9917, .98,   .952,  .9142, .8672, .8062, .741,
     b  .594,  .4564, .335,  .241,  .1665, .1149, .0762, .0503, .0328,
     c  .021,  .0054/
      Data P6/  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     a   1.0,   1.0,  .998,  .992,  .977,  .9562, .9292, .8922, .841,
     b  .697,  .5485, .407,  .294,  .2045, .1419, .0947, .0628, .0416,
     c  .0268, .0074/
      Data P7/  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     a   1.0,   1.0,   1.0,  .998,  .992,  .9822, .9687, .9462, .908,
     b  .769,  .6135, .46,   .335,  .2363, .1659, .1127, .0758, .0511,
     c  .0336, .0104/
      Data P8/  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     a   1.0,   1.0,   1.0,   1.0,   1.0,  .9997, .9967, .9882, .966,
     b  .8375, .6815, .518,  .381,  .2713, .1916, .1307, .0878, .0591,
     c  .0386, .012/
      Data P9/  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     a   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  .9997, .991,
     b  .8855, .7395, .576,  .431,  .3108, .2216, .1532, .1044, .0711,
     c  .0471, .0153/
      Data P10/ 1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     a   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  .992,
     b  .8917, .7532, .5925, .447,  .3252, .2336, .1627, .1114, .0759,
     c  .0501, .0163/
c
      E=Exc9Be
      P(1) =EXTERP(Ex, P1,E,N,Nterp)
      P(2) =EXTERP(Ex, P2,E,N,Nterp)
      P(3) =EXTERP(Ex, P3,E,N,Nterp)
      P(4) =EXTERP(Ex, P4,E,N,Nterp)
      P(5) =EXTERP(Ex, P5,E,N,Nterp)
      P(6) =EXTERP(Ex, P6,E,N,Nterp)
      P(7) =EXTERP(Ex, P7,E,N,Nterp)
      P(8) =EXTERP(Ex, P8,E,N,Nterp)
      P(9) =EXTERP(Ex, P9,E,N,Nterp)
      P(10)=EXTERP(Ex,P10,E,N,Nterp)
c
      Pran=RAN(Irx)
      M=1
      Do 1 J=1,10
      IF (Pran .LE. P(J)) goto 2
    1 M=M+1
    2 IF (M .EQ. 11) goto 3
      ExLi8=E8Li(M)
      Return
c
c   Or get highly-excited 8-Li, E(Level) from "continuum"
    3 ExLi8=E8Li(10)
      Emax=E-29.4
      IF (Emax .LE. 0.0) Return
      Ehat=Emax+5.0
      F=0.1245 + 0.001*ABS(Ehat-45.0)
      Temp=F*Ehat
      Ep=CHOOSP(Emax,Temp)
      ExLi8=E-Ep-Li8pp
      Return
      END
c==========================================
C   This is file PP8LI.FOR
c
      Subroutine PP8LI(Excit,V9Be,Eg)
c
c       Next portion of programming for 9-Be --> p + 8-Li decay.


c
      Real*4 Li8pp,Li7pn
      Common /MASSES/ Emass(21)
      Common /PROTON/NOPROT,RNGEN(144),PRAN(144),DRAN(144),ARAN(144),
     a              HRAN(144),TRAN(144)
c      Common /PROTON/ Noprot,Range(95),Proten(95)
      Common /NEUTRN/ Eneut, U,V,W, X,Y,Z
      Common /RANDM/  Irx
      Common /NAID/   Rad,Ht,R1,R2
      Common /COLLIS/ Nelm,Echrg(6)
      Common /VECTOR/ Xb,Yb,Zb, Xn,Yn,Zn
c       At entry, Xn,Yn,Zn are dir. cos. of 9-Be ion motion in "detector"


c          coordinates.


c
      Data Li8pp/16.888/, Li7pn/2.033/, Kp/2/, K8Li/16/
c
      Ta=Excit-Li8pp
c       Test for possible alert to this routine to do only proton decay,


c         i.e., no further decay of 8-Li to  n + 7-Li.


      IF (V9Be .GT. 0.0) goto 4
      V9Be=-V9Be
      Excitp=0.0
c       Choose between ground state and first-excited state of 8-Li.


      IF (RAN(Irx) .LT. 0.3) Excitp=0.9808
      Goto 6
c           Next is for no limits on 8-Li excitation energy.


    4 Excitp=EXLI8(Excit)
c       That gets 8-Li excitation energy.


      Ix=1
      IF (Excitp .GT. Li7pn) goto 7
    6 Eg=Excitp
      Ix=0
    7 Taex=Ta-Excitp
      Epp=RCKE(Kp,K8Li,Taex)
      Vpp=VELOCITY(Kp,Epp)
      CALL RVECT(Xp,Yp,Zp)
      Vpx=Xp*Vpp
      Vpy=Yp*Vpp
      Vpz=Zp*Vpp
      CALL LABTRAN(Vpx,Vpy,Vpz, V9Be, Vx,Vy,Vz)
      Vprot=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Eprot=EFROMV(Kp,Vprot)
      Echrg(4)=Eprot
c          Check for proton escape:


      CALL DIRCOS(Vx,Vy,Vz, Xb,Yb,Zb)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
      Protpl=PLNGTH(Xn,Yn,Zn, X,Y,Z, Rad,Ht)
      Rprot=EXTERP(RNGEN,PRAN, Eprot, Noprot,4)
      RangeQ=Rprot-Protpl
      IF (RangeQ .LE. 1.0200E-04) goto 8
c       Proton escapes.

      Echrg(2)=EXTERP(PRAN,RNGEN,RangeQ,Noprot,4) ! kajimoto 2011/03/18
c
c       Now get 8-Li ion:

    8 E8Li=Taex-Epp
      V8Lic=VELOCITY(K8Li,E8Li)
      V8x=-Xp*V8Lic
      V8y=-Yp*V8Lic
      V8z=-Zp*V8Lic
      CALL LABTRAN(V8x,V8y,V8z, V9Be, Vx,Vy,Vz)
      V8Li=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      IF (Ix .EQ. 1) goto 2
c        If Ix = 0, then no additional neutron decay. So tidy up and exit.


      Echrg(3)=EFROMV(K8Li,V8Li)
      Return
c
c       Set up for 8-Li decay.


    2 CALL DIRCOS(Vx,Vy,Vz, Xb,Yb,Zb)
      CALL TRANSVEC(Zx,Zy,Zz)
c       Set up done.  Get 8-Li decay.


      CALL ELI8DK(Excitp,V8Li,Egamma)
      Eg=Egamma
c       That's that!


      Return
      END
c==========================================
c   This is file TENBDK.FOR
c
c       Purpose is to follow the decay of particle-unstable excited


c         states of 10-B.  This routine is used in Subroutines ND, NT,


c         and NPX.


c
c       Choices of decay modes are:


c         (a)   10-B --> p + 9-Be (ground state)


c       or,


c         (a')  10-B --> d + 8-Be


c         (b')   8-Be--> alpha + alpha


c       or,


c         (a")  10-B --> alpha + 6-Li


c         (b")   6-Li--> gamma + 6-Li (ground state)


c       or,


c         (b")   6-Li--> d + alpha


c
c       or, under certain circumstances,


c         (a3)  10-B --> p + 9-Be (excited)


c         (b3)   9-Be--> n + 8-Be


c         (c3)   8-Be--> alpha + alpha


c       or, fractionally,


c         (b3')  9-Be--> d + 7-Li


c
      Subroutine TENBDK(Excit,V10B,Eg)
c
      Real*4 Li6pa
      Common /COLLIS/ Nelm,Echrg(6)
      Common /MASSES/ Emass(21)
      Common /RANDM/  Irx
      Common /NAID/   Rad,Ht,R1,R2
      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
c          At entry to this routine Xn,Yn,Zn are direction cosines of


c          moving 10-B ion in "detector" lab. coordinates.


      Common /NEUTRN/ Eneut, Dc1,Dc2,Dc3, Xneut,Yneut,Zneut
csatoh20041214@jaeri
      COMMON/PROTON/NOPROT,RNGEN(144),PRAN(144),DRAN(144),ARAN(144),
     a              HRAN(144),TRAN(144)
c      Common /PROTON/ Noprot,Range(95),Proten(95)
      Common /EXC6LI/ ExLi6(6),He4pd
      Common /BENINE/ Ex9Be(6),Be8pn
c          See Subroutine N3He for values in BENINE common.


c
c         Next array represents level structure of 6-Li.
      Data ExLi6/0.0,2.185, 3.563, 4.31, 5.37, 5.65/
c
      Data Kp/2/, K6Li/18/, Kd/10/, Ka/3/, K9Be/4/, Kn/1/, K8Be/8/
      Data Qnt/18.93/, Be8pd/6.026/, Be9pp/6.585/, Be8aa/0.092/
      Data He4pd/1.474/, Li6pa/4.46/, Ntrp/1/
c
      Eg=0.0
      Ta=Excit
c       Choose which decay mode to follow.


      IF (Excit .LE. Be8pd) goto 10
      IF (Excit .LE. Be9pp) goto 5
c       Such experimental information as exists suggests that about half


c         of the highly-excited levels in 10-B decay predominantly by


c         proton emission.


      IF (RAN(Irx) .GE. 0.5) goto 5
c
c   Study p + 9-Be mode first.
      Taex=Ta-Be9pp
      Try=Taex
      Ep=RCKE(Kp,K9Be,Taex)
      Npgo=-1
      Excit=0.0
      IF (RAN(Irx) .LE. 0.04) goto 4
c       Last is to slightly enhance 9-Be ground state -- 4% is ad hoc.


c
c       Next test variable -ENEUT-; if it's zero then can set up for


c         possible 9-Be --> n + 8-Be calculation.


      IF (Eneut .GT. 0.0) goto 4
c       -ENEUT- is zero, so look for decay mode (a3), above.


      Ehat=Ep+5.0
      F=0.1245 + 0.001*ABS(Ehat-45.0)
      T=F*Ehat
      Eprot=CHOOSP(Ep,T)
      Try=Eprot*(Emass(Kp) + Emass(K9Be))/Emass(K9Be)
      Eptry=RCKE(Kp,K9Be,Try)
      Try=Try*Eprot/Eptry
      Excit=Taex-Try
c          (No test for 9-Be --> p + 8-Li as already have one proton's info.)


c       Next check possible  9-Be --> d + 7-Li decay mode.


      Rani=RAN(Irx)
      IF (Excit .LT. 17.4 .OR. Rani .GT. 0.09) goto 2
      Npgo=1
      Goto 9
    2 Npgo=0
      IF (Excit .GE. 16.0) goto 9
      K=1
      Do 6 I=2,6
      IF (Excit .LT. Ex9Be(I)) goto 7
    6 K=K+1
    7 Excit=Ex9Be(K)
c       Recheck for 9-Be ground state at this point.


      IF (K .EQ. 1) Npgo=-1
      Try=Taex-Excit
    9 Ep=RCKE(Kp,K9Be,Try)
    4 Vp=VELOCITY(Kp,Ep)
c         "Fission" of 10-B is isotropic in center of mass coords.


      CALL RVECT(Xp,Yp,Zp)
      Vpx=Xp*Vp
      Vpy=Yp*Vp
      Vpz=Zp*Vp
      CALL LABTRAN(Vpx,Vpy,Vpz, V10B, Vx,Vy,Vz)
      Vpl=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Epl=EFROMV(Kp,Vpl)
      Echrg(5)=Epl
      Echrg(6)=0.0
c       Check for possible proton escape from detector:


      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Za=Xn
      Zb=Yn
      Zc=Zn
      CALL TRANSVEC(Za,Zb,Zc)
      Protpl=PLNGTH(Xn,Yn,Zn, Xneut,Yneut,Zneut, Rad,Ht)
      Rprot=EXTERP(RNGEN,PRAN,Epl,Noprot,4)
      RangeQ=Rprot-Protpl
      IF (RangeQ .LE. 1.0200E-04) goto 3
      Echrg(6)=EXTERP(PRAN,RNGEN,RangeQ,Noprot,4) ! kajimoto 2011/03/18
c       Next: Get information for 9-Be ion.


    3 E9Bec=Try-Ep
      V9Bec=VELOCITY(K9Be,E9Bec)
      V9x=-Xp*V9Bec
      V9y=-Yp*V9Bec
      V9z=-Zp*V9Bec
      CALL LABTRAN(V9x,V9y,V9z, V10B, Vx,Vy,Vz)
      V9Be=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
C
      IF (Npgo .EQ. -1) goto 13
      IF (Npgo .EQ. 1) goto 25
c       Program counter to here have 9-Be --> n + 8-Be decay mode


c         First rotate 9-Be ion motion into "detector" coordinates:


      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Za,Zb,Zc)
c       Done.


      Tax=Excit-Be8aa-Be8pn
      Enc=RCKE(Kn,K8Be,Tax)
      Vnc=VELOCITY(Kn,Enc)
      CALL RVECT(Xp,Yp,Zp)
      Vnx=Xp*Vnc
      Vny=Yp*Vnc
      Vnz=Zp*Vnc
      CALL LABTRAN(Vnx,Vny,Vnz, V9Be, Vx,Vy,Vz)
      Vnlab=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Enn=EFROMV(Kn,Vnlab)
c       That's the neutron's energy.  Store it in /NEUTRN/ common area.


      Eneut=Enn
c       Rotate neutron's motion into "detector" coordinates and store them.


      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Za=Xn
      Zb=Yn
      Zc=Zn
      CALL TRANSVEC(Za,Zb,Zc)
      Dc1=Xn
      Dc2=Yn
      Dc3=Zn
c       Now get 8-Be ion motion.


      E8Bec=Tax-Enc
      V8Bec=VELOCITY(K8Be,E8Bec)
      V8x=-Xp*V8Bec
      V8y=-Yp*V8Bec
      V8z=-Zp*V8Bec
      CALL LABTRAN(V8x,V8y,V8z, V9Be, Vx,Vy,Vz)
      V8Be=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
c       Done.  Now get  8-Be --> 2 alphas.


      Iechrg=3
      Goto 16
c
c       Next it to tidy up  10-B --> p + 9-Be (ground state) decay mode.


   13 Echrg(3)=EFROMV(K9Be,V9Be)
      Return
c       END SECTION on p + 9-Be


c
c  Choosing between the other two modes of 10-B decay appears a toss up.
    5 IF (RAN(Irx) .GE. 0.5) goto 10
c
c       Do 10-B --> d + 8-Be decay mode next.  Ignore 8-Be excited states


c          as we are already high in (effective) 12-C excitation, and these


c          branches of decay are relatively minor.


c
      Taex=Ta-Be8pd
      Ed=RCKE(Kd,K8Be,Taex)
      Vd=VELOCITY(Kd,Ed)
      CALL RVECT(Xd,Yd,Zd)
      Vdx=Xd*Vd
      Vdy=Yd*Vd
      Vdz=Zd*Vd
      CALL LABTRAN(Vdx,Vdy,Vdz, V10B, Vx,Vy,Vz)
      Vdl=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Echrg(3)=EFROMV(Kd,Vdl)
c       and as above for the proton, ignore possible deuteron escape from


c       the detector.


      E8Bec=Taex-Ed
      V8Bec=VELOCITY(K8Be,E8Bec)
      V8x=-Xd*V8Bec
      V8y=-Yd*V8Bec
      V8z=-Zd*V8Bec
      CALL LABTRAN(V8x,V8y,V8z, V10B, Vx,Vy,Vz)
      V8Be=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
c       Now get 8-Be --> alpha + alpha(2)


      Iechrg=4
   16 Ta=Be8aa
      Ea=0.5*Ta
c           Alphas share the available energy.


      Va=VELOCITY(Ka,Ea)
      CALL RVECT(Xc,Yc,Zc)
      Vax=Xc*Va
      Vay=Yc*Va
      Vaz=Zc*Va
      CALL LABTRAN(Vax,Vay,Vaz, V8Be, Vx,Vy,Vz)
      Valph=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Echrg(Iechrg)=EFROMV(Ka,Valph)
c         Recall symmetry so that Va2x=-Vax, etc, for the other


c         alpha.


      Va2x=-Vax
      Va2y=-Vay
      Va2z=-Vaz
      CALL LABTRAN(Va2x,Va2y,Va2z, V8Be, Vx,Vy,Vz)
      Valph2=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Echrg(Iechrg+1)=EFROMV(Ka,Valph2)
      Return
c       END SECTION for 10-B --> d + 8-Be


c               and for 10-B --> p + n + 8-Be.


c
c  That leaves 10-B --> alpha + 6-Li decay modes to do.
   10 Taex=Ta-Li6pa
      Eax=RCKE(Ka,K6Li,Taex)
      Eg=0.0
      Excita=0.0
      K=1
      IF (RAN(Irx) .LE. 0.03) goto 18
c       As done above, a slight (3%, ad hoc) enhancement of ground-state.


      Temp=4.0
c   See comments for alpha "continuum" computation just prior to statement
c       no. 38 in subroutine NN3ALF.


      Eatry=CHOOSP(Eax,Temp)
      Try=Eatry*(Emass(Ka) + Emass(K6Li))/Emass(K6Li)
      Etry=RCKE(Ka,K6Li,Try)
      Try=Try*Eatry/Etry
      Excita=Taex-Try
c          Now pair -Excita- with a level in 6-Li
      K=7
      IF (Excita .GT. 20.0) goto 18
      K=1
      Do 14 I=2,6
      IF (Excita .LT. ExLi6(I)) goto 15
   14 K=K+1
   15 Excita=ExLi6(K)
c       6-Li level decay scheme is as follows:


c               K=1     Ex=0.0          Stable


c               K=2     Ex=2.185        d + alpha (very weak gamma)


c               K=3     Ex=3.563        ground-state gamma ray


c               K=4     Ex=4.31         d + alpha


c               K=5     Ex=5.366        weak g.s. gamma reported, but


c                                        p or n decay most likely.


c               K=6     Ex=5.65         d + alpha


c               K>6     Ex>20           triton + 3-He


c   For the present programming, take K=3 gamma decay and the remaining
c       K>1 excited state decay to be d + alpha.


c
      IF (K .EQ. 3) Eg=Excita
c       Now get information on alpha.


   18 Taex=Taex-Excita
      Eaa=RCKE(Ka,K6Li,Taex)
      Vaa=VELOCITY(Ka,Eaa)
      CALL RVECT(Xa,Ya,Za)
      Vax=Xa*Vaa
      Vay=Ya*Vaa
      Vaz=Za*Vaa
      CALL LABTRAN(Vax,Vay,Vaz, V10B, Vx,Vy,Vz)
      Va=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Echrg(4)=EFROMV(Ka,Va)
c       Now set up 6-Li ion.


      E6Lic=RCKE(K6Li,Ka,Taex)
      V6Lic=VELOCITY(K6Li,E6Lic)
      V6x=-Xa*V6Lic
      V6y=-Ya*V6Lic
      V6z=-Za*V6Lic
      CALL LABTRAN(V6x,V6y,V6z, V10B, Vx,Vy,Vz)
      V6Li=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      IF(K .EQ. 2) goto 20
      IF (K. GE. 4) goto 20
c       For K=1 or 3 tidy up 6-Li energy.


      Echrg(3)=EFROMV(K6Li,V6Li)
      Return
c       Penultimate section: 6-Li --> d + alpha
   20 CALL ALFAPD(Excita,V6Li)
      Return
c
c       Last section:  9-Be --> d + 7-Li  (no deuteron escape considered)


   25 CALL DPLI7(Excit,V9Be)
      Return
c
c   At exit from this subroutine the -Echrg- array will have information on
c       decay of 10-B as follows:


c
c   Final decay mode    |   Echrg(3)    |   Echrg(4)    |   Echrg(5)


c    ------------------ | ------------- | ------------- | -------------
c       p + 9-Be        |   E(9-Be)     |     0.0       |    E(p)


c    p + n + 2 alpha    |   E(alpha)    |   E(alpha)    |    E(p)
c     alpha + 6-Li      |   E(6-Li)     |   E(alpha)    |    0.0


c      d + 2 alpha      |    E(d)       |   E(alpha)    |   E(alpha)


c     p + d + 7-Li      |   E(7-Li)     |    -E(d)      |    0.0
c
c   For reactions  p + 9-Be, p + d + 7-Li, and p + n + 2alpha, Echrg(6)
c        will have proton escape information, if any.


c       For the other two reactions, Echrg(6) = 0.


      END
c==========================================
