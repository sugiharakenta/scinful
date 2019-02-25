*****************************************************
*  SUBROUTINES                                      *
*     hydrog                                        *
*     pscat                                         *
*     inelas                                        *
*     nalpha                                        *
*     nn3alf                                        *
*     n3he                                          *
*     npx                                           *
*     nd                                            *
*     nt                                            *
*     n2n                                           *
*     np                                            *
*     cscat                                         *
*     proton (add at 04.12.14,modified at 04.12.15) *
*                                                   *
*****************************************************
c This is file HYDROG.FOR
c     Purpose is to do the calculations for the neutron-proton
c       collision.  Some of the programming is taken from the O5S
c       coding.
c
cs_4v
      SUBROUTINE HYDROG
C
C    HYDROG CONTAINS THE KINEMATICS FOR THE N+H = N+H COLLISION,
C       THE POSSIBILITY OF THE RECOIL PROTON LEAKING OUT OF
C       THE SCINTILLATOR IS CONSIDERED.
C
C
      Common /RANDM/ IRX
      Common /NEUTRN/ SPDSQ, U, V, W, X, Y, Z
      COMMON/PROTON/NOPROT,RNGEN(144),PRAN(144),DRAN(144),ARAN(144),
     a              HRAN(144),TRAN(144)
c          /PROTON/ common data in File PROTN.FOR


      Common /MASSES/ Emass(21)
      Common /NAID/ Rdet,Ht,RC,Rz
c
c by d.satoh for Linux (2005.01.11@jaeri)
c      Common /COLLIS/ Nelm, Echrg(3)
      Common /COLLIS/ Nelm, Echrg(6)
      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
c
ckaji 2011/03/18
      common /RES/ iswt, bias(5), icont, ictrl, iesc, ilight
c
C      Dimension Energy(16),Alpha(16),Beta(16),Gamma(16)
      Dimension Energy(13),Alpha(13),Beta(13),Gamma(13)
      Dimension Delta(13),EPSIR(13)
      Dimension XMU(21),XMUSQ(21),XMUCUB(21),XMUTRI(21),XMUFOR(21)


      DIMENSION PMU(21)
c
      DATA  ENTOLD/0.0/, CKVALU/13.7/, Kn/1/, Kp/2/
C      DATA NUMBER/16/
      DATA NUMBER/13/
C      DATA (ENERGY(I),I=1,16)/13.7,14.1,17.9,22.5,27.5,32.5,37.5,42.5,
C     &47.5,52.5,57.5,62.5,70.0,80.0,90.0,100./
      DATA ENERGY /13.7   ,14.1   ,17.9  ,20.0   ,25.0   ,30.0   ,
     &     40.0   ,50.0   ,60.0   ,70.0  ,80.0   ,90.0   ,100.0/


C      DATA (ALPHA(I),I=1,16)/56.4,52.2,40.0,32.9,27.0,21.8,17.8,13.8,
C     &13.3,11.7,9.4, 7.6, 5.6, 5.4, 4.4, 3.1/
      Data Alpha /56.4    ,52.2   ,40.0  ,37.3878,29.1367,23.3544,
     &     15.7512,11.3439,8.63440,6.75356,5.45053,4.5356,3.86055/
C      DATA (BETA(I),I=1,16)/ 0.0,-3.4,-3.4,0.15,0.05,0.85,.3,0.30,.4,1.,
C     10.80,0.80,0.05,-.65,-.15,-.40/
       Data Beta /0.0     ,-3.4   ,-3.4   ,-1.18 ,-.75780,-.412017,
     &     .021592,.132066,-4.63e-3,-.24424,-.50353,-.74459,-.94257/
C      DATA (GAMMA(I),I=1,16)/0.0,2.0,1.4,3.25,3.05,5.55,5.1,7.0,6.3,6.8,
C     1 7.4, 7.2,8.45,7.05,8.05, 7.4/
       Data Gamma/0.0     ,2.0    ,1.4    ,1.86507,2.90311,3.74780,
     &     4.87744,5.49871,5.78747,5.78322,5.61946,5.34592,4.93463/
       Data Delta/0.0     ,0.0    ,0.0   ,-1.06337,-.92175,-.731767,
     &   -.3693726,-.13831,-.0372414,6.1575e-3,4.10505e-2,9.0466e-2,
     &   1.682531e-1/
       Data Epsir/0.0     ,0.0    ,0.0   ,1.028816,.906036,.772477,
     &   .5947638,.6040201,.7999593,1.103231,1.485119,1.933848,
     &   2.434463/


       DATA XMU /-1.00,-0.90,-0.80,-0.70,
     &-0.60,-0.50,-0.40,-0.30,-0.20,-0.10,-0.00, 0.10,
     & 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00/
        DATA XMUSQ /1.000, .810, .640, .490,
     & .360, .250, .160, .090, .040, .010, .000, .010,
     & .040, .090, .160, .250, .360, .490, .640, .810,1.000/
       DATA XMUCUB /-1.0000,-0.7290,-0.5120,-0.3430,
     &-0.2160,-0.1250,-0.0640,-0.0270,-0.0080,-0.0010,-0.0000, 0.0010,
     & .0080, 0.0270, 0.0640, 0.1250, 0.2160, 0.3430, 0.5120, 0.7290,1./
       DATA XMUTRI/1.0,0.6561,0.4096,0.2401,0.1296,
     & 0.0625,0.0256,0.0081,0.0016,0.0001,0.,0.0001,0.0016,0.0081,
     & 0.0256,0.0625,0.1296,0.2401,0.4096,0.6561,1.0/
       DATA XMUFOR/-1.0,-0.59049,-0.32768,-0.16807,-0.07776,
     &-0.03125,-0.01024,-0.00243,-0.00032,-0.00001,0.,0.00001,0.00032,
     & 0.00243,0.01024,0.03125,0.07776,0.16807,0.32768,0.59049,1.0/
c
      En=SPDSQ
      Vn=VELOCITY(Kn,En)
      CALL CMTRAN(Kn,Kp,Vn, Vcom,Enc,Epc)
      Vnc=VELOCITY(Kn,Enc)
      Vpc=VELOCITY(Kp,Epc)
C
C * * * IF SPDSQ IS LESS THAN CKVALU USE ISOTROPIC (CM) SCHEME FOR
C * * * DETERMINING THE DIRECTION COSINES OF THE OUTGOING NEUTRON.
C
      IF ( SPDSQ .LT. CKVALU )   GO TO 140
C
C * * * DETERMINE THE VALUES OF A, B AND C TO BE USED IN CALCULATING
C * * * THE ANGULAR DISTRIBUTION FUNCTION.  ( THE ANGULAR DISTRIBUTION
C * * * IS APPROXIMATED BY P(E,MU) = A(E) + B(E) * MU + C(E)* MU**2 ).
C
C       CALCULATION OF XMU AND ITS POWERS HAS BEEN REMOVED, AND


C       A LOOK-UP TABLE IS USED INSTEAD.  (ZWBell addition, I think.)
C
      IF ( SPDSQ .EQ. ENTOLD )   GO TO 110
      ENTOLD = SPDSQ
C
      Nterp=1
      A=EXTERP(Energy,Alpha,En,Number,Nterp)
      B=EXTERP(Energy,Beta,En,Number,Nterp)
      C=EXTERP(Energy,Gamma,En,Number,Nterp)
      D=EXTERP(Energy,DELTA,En,Number,Nterp)
      E=EXTERP(Energy,EPSIR,En,Number,Nterp)
c
C * * * INTEGRATE THE ANGULAR DISTRIBUTION FUNCTION.
C
   80 B2 = B / 2.
      C3 = C / 3.
      D4 = D/  4.
      E5 = E/  5.
      DO  90  K = 1, 21
      PMU(K) = A * ( XMU(K) + 1. ) + B2 * ( XMUSQ(K) - 1. ) +
     A         C3 * ( XMUCUB(K) + 1. )
     A       + D4 * ( XMUTRI(K) - 1. )
     A       + E5 * ( XMUFOR(K) + 1. )
   90 CONTINUE
C
C * * * PICK MU FROM THE DISTRIBUTION ON INTERVAL [0..PMU(21)]
C
  110 RAND = RAN(IRX)*PMU(21)
      Nterp=1
      Xxmu=EXTERP(Pmu,Xmu,Rand,21,Nterp)
C
C * * * CHOOSE AZIMUTHAL DIRECTION.
C
        RANN = 6.28318530*RAN(IRX)
        SINPHI = SIN(RANN)
        COSPHI = COS(RANN)
C
      Sintheta=SQRT(1.0-Xxmu*Xxmu)
      Xp=Sintheta*Cosphi
      Yp=Sintheta*Sinphi
      Zp=Xxmu
      Goto 15
c
  140 CALL RVECT(Xp,Yp,Zp)
c       Direction cosines of outgoing neutron are Xp, Yp, Zp.
   15 Vnx=Xp*Vnc
      Vny=Yp*Vnc
      Vnz=Zp*Vnc
c       Vnx,Vny,Vnz are components of neutron velocity in the "neutron"
c         center-of-mass coordinate system.  Get the same for the proton.
      Vpx=-Xp*Vpc
      Vpy=-Yp*Vpc
      Vpz=-Zp*Vpc
c        For the neutron -- get velocity components in "neutron" lab.
c          coordinate system.
      CALL LABTRAN(Vnx,Vny,Vnz, Vcom, Vx,Vy,Vz)
      Vneut=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
c        Now get neutron's energy, then get its direction cosines in the
c          "neutron" lab. system, finally rotate those into "detector"
c          coordinate system using TRANSVEC.
      Eneut=EFROMV(Kn,Vneut)
      SPDSQ=Eneut
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Zx=U
      Zy=V
      Zz=W
      CALL TRANSVEC(Zx,Zy,Zz)
      U=Xn
      V=Yn
      W=Zn
c
c       Now determine proton's energy and see if the proton stays in the
c          detector or gets out.
      CALL LABTRAN(Vpx,Vpy,Vpz, Vcom, Vx,Vy,Vz)
      Vprot=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Eprot=EFROMV(Kp,Vprot)
      Echrg(1)=Eprot
      Echrg(2)=0.0
c      write(26,*) Echrg(1)
c
ckaji      IF (Eprot .LE. 0.1) goto 20
      IF (Eprot .LE. 0.07) goto 20
c          Eprot < 100 keV, range a few micrometers.
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
c       Xn,Yn,Zn are the proton's direction cosines in lab coordinates.
      Protpl=PLNGTH(Xn,Yn,Zn, X,Y,Z, Rdet,Ht)
      Rprot=EXTERP(RNGEN,PRAN,Eprot,Noprot,4)


c       Now compare proton's range in detector (RPROT) with the path length
c         of the vector having direction cosines Xn, Yn, Zn from position
c         X, Y, Z (PROTPL).
      RangeQ=Rprot-Protpl
ckaji      IF (RangeQ .LE. 0.1) goto 20
      IF (RangeQ .LE. 1.0200E-04) goto 20
c         If it goes to next step, the proton got out of the detector.
      Eprotf=EXTERP(PRAN,RNGEN,RangeQ,Noprot,4) ! kajimoto 2011/03/18
      Echrg(2)=Eprotf
   20 Nelm=1
      CALL BANKR2
      Return
      END
c------------------------------------------------------------------
      SUBROUTINE PSCAT
C
      Common /RANDM/ IRX
      Common /NEUTRN/ Eneut, U, V, W, X1, Y, Z
      Common /LGNDRE/ IC, E(220), F(6,220)
c       See file LEGENDRE.FOR for the data in common /LGNDRE/ arrays.


      Common /COLLIS/ Nelm, Echrg(6)
      Common /CSCATT/ Xp,Yp,Zp, Vnca,Vcca,Vcm
ckaji 2011/03/18
      common /RES/ iswt, bias(5), icont, ictrl, iesc, ilight
c
      Dimension Fi(6),A1(220),B1(220),B2(220),B3(220),B4(220)
      Dimension B5(220),B6(220)
      Data Nthru/0/, Nterp/3/, Kn/1/, KC/7/, Nfi/6/
c
C     For the present programming to obtain the neutron's polar scattering
C     angle the Legendre polynomial "technique" is used.  The BLOCK DATA have
C       coefficients for E(neut) up to 20 MeV from the ENDF/B evaluation;


c       whether the values are correct or not I can't say.  For larger


c       E(neut) values are included which have been deduced from comparisons


c       with experimental angular distributions as 20.8, 26 and 40 MeV.


c
c     It appears, from the O5S coding, that the values of the Legendre
c       coefficients in the Block Data are tabulated for specific incident


c       neutron energies in MeV in the laboratory frame of reference, and the


c       values themselves are for angular distributions in the center of mass.


c
      If (Nthru .GT. 0) goto 10
c               Set-up needed only once -- the first time through.
c==CHECK!===============
c          write(17,*)'PSCAT_IN: Eneut=',Eneut
c==CHECK!===============


      Do 2 I=1,Ic
      A1(I)=E(I)
      B1(I)=F(1,I)
      B2(I)=F(2,I)
      B3(I)=F(3,I)
      B4(I)=F(4,I)
      B5(I)=F(5,I)
    2 B6(I)=F(6,I)
      Nthru=1
c
C * * * CHOOSE THE COSINE OF THE POLAR ANGLE OF THE OUTGOING NEUTRON.
c
   10 En=Eneut
      Vn=VELOCITY(Kn,En)
      CALL CMTRAN(Kn,KC,Vn, Vcom,Enc,Ecc)
      Vcm=Vcom
      Vnca=VELOCITY(Kn,Enc)
      Vcca=VELOCITY(KC,Ecc)
c         Vnca, Vcca = velocities of neutron and Carbon ion in center of mass


c
      Fi(1)=EXTERP(A1,B1,En,Ic,Nterp)
      Fi(2)=EXTERP(A1,B2,En,Ic,Nterp)
      Fi(3)=EXTERP(A1,B3,En,Ic,Nterp)
      Fi(4)=EXTERP(A1,B4,En,Ic,Nterp)
      Fi(5)=EXTERP(A1,B5,En,Ic,Nterp)
      Fi(6)=EXTERP(A1,B6,En,Ic,Nterp)
c
c       Get polar scattering angle from function CHOOSL


      FMU = CHOOSL(Fi,NFi)
c
   24 SINPSI = SQRT ( 1. - FMU * FMU )
        RANN = 6.28318530*RAN(IRX)


        SINETA = SIN(RANN)


        COSETA = COS(RANN)


c
      Xp=Sinpsi*Coseta
      Yp=Sinpsi*Sineta
      Zp=Fmu
c
      CALL CSCAT
c               CSCAT finishes up -- gets new neutron energy, dir. cosines.
      Nelm=2
      CALL BANKR2
      Return
      END
c------------------------------------------------------------------
C   This is file INELAS.FOR
c
      SUBROUTINE INELAS
C  INELAS DOES THE CALCULATIONS FOR THE 12C INELASTIC REACTION.
C
C
      Common /RANDM/ IRX
      Common /NEUTRN/ Eneut, Ax, Ay, Az, X, Y, Z
      Common /COLLIS/ Nelm,Echrg(6)
      Common /CSCATT/ Xpn,Ypn,Zpn,Vnca,Vcca,Vcm
c
      Dimension A1(10),B1(10),B2(10),B3(10),B4(10),Fi(4)
c
      Data Ipoly/10/, Nterp/1/, Nfi/4/
      Data Q/4.433/, Rmass/0.922461/, Kn/1/, K12C/7/, TenMeV/10.0/
c
c         Next arrays to get anisotropic scattering distribution.
c               (See Glasgow, et al, Nuclear Science and Engineering, 61
c                       (1976) 521 for data between 9.19 and 13 MeV.)
c
      Data A1/ 6.0, 8.56, 9.19, 10.69, 10.96, 11.73,
     m  12.95, 14.6, 20.8, 26.0/
      Data B1/0.07784, 0.03034, -0.0396, 0.166, 0.2027, 0.258,
     m 0.1733, 0.21885, 0.39757, 0.55621/
      Data B2/0.04819, 0.17527, 0.1836, 0.203, 0.2504, 0.2406,
     m 0.2843, 0.209, 0.20374, 0.24515/
      Data B3/0.004335, -0.03078, -0.00117, 0.0587, 0.0376, 0.0448,
     m 0.0489, 0.03678, 0.07678, 0.09518/
      Data B4/0.0, 0.01355, 0.0182, 0.0186, 0.0276, 0.0128,
     m 0.0233, -0.001315, 0.01595, 0.01973/
c
C       Enter with incident neutron energy, direction cosines.
C       Exit with final neutron energy, carbon-ion energy, new direction
C        cosines of scattered neutron.  Spot of interaction is unchanged.
c
      En=Eneut
      Vn=VELOCITY(Kn,En)
      CALL CMTRAN(Kn,K12C,Vn, Vcom,Enc,Ecc)
      Vcm=Vcom
      Tec=Enc+Ecc
      Ta=Tec-Q
      IF (Ta .GT. 0.0) goto 10
csatoh
c     Type 5, En
      write(*,5) en
    5 Format(/'  *** Error Subroutine INELAS; E(neut) at entry = '
     V    1PE11.4/)
      Eneut=0.0
      Echrg(1)=0.0
      Return
c
   10 Ena=Rmass*Ta
      IF (Ta .GT. TenMeV) Ena=RCKE(Kn,K12C,Ta)
      Eca=Ta-Ena
      Vnca=VELOCITY(Kn,Ena)
      Vcca=VELOCITY(K12C,Eca)
c
c       O5S treats inelastic scattering as isotropic for all incident


c        neutron energies.  We will not do the same but will instead


c        use a relatively crude grid of Legendre coefficients to get


c        a handle on the angular distribution of the scattered neutrons.


c
      Fi(1)=EXTERP(A1,B1,En,Ipoly,Nterp)
      Fi(2)=EXTERP(A1,B2,En,Ipoly,Nterp)
      Fi(3)=EXTERP(A1,B3,En,Ipoly,Nterp)
      Fi(4)=EXTERP(A1,B4,En,Ipoly,Nterp)
      Fmu=CHOOSL(Fi,Nfi)
      Sinpsi=SQRT(1.0 - Fmu*Fmu)
      Phi=6.283185*RAN(Irx)
      Xpn=Sinpsi*COS(Phi)
      Ypn=Sinpsi*SIN(Phi)
      Zpn=Fmu
      CALL CSCAT
c
      Nelm=3
      CALL BANKR2
c       Now to check for possible Compton scattering of photon in detector


      Egamma=Q
      CALL PHOTON(Egamma)
      Return
      END
c------------------------------------------------------------------
c   This is file NALPHA.FOR
C
       SUBROUTINE NALPHA
c
C  NALPHA DOES THE CALCULATIONS FOR THE REACTION N + 12C -> ALPHA + 9BE.
C   12-C(N,ALPHA)9-BE COLLISION
C
      Common /RANDM/ IRX
      Common /NEUTRN/ Eneut,U,V,W,X,Y,Z
      Common /COLLIS/ Nelm,Echrg(6)
ckaji 2011/03/18
      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
      Common /NAID/ Rdet,Ht,RC,Rz
      common /RES/ iswt, bias(5), icont, ictrl, iesc, ilight
      Common /PROTON/NOPROT,RNGEN(144),PRAN(144),DRAN(144),ARAN(144),
     a              HRAN(144),TRAN(144)
      Common /MASSES/ Emass(21) ! added by kajimoto 06/13
c
      Dimension A1(10),B1(10),B2(10),B3(10),B4(10),B5(10),B6(10),Fi(6)
c
      Data Q/5.71/, Kn/1/, K12C/7/, K9Be/4/, Ka/3/
      Data Ipoly/10/, Nfi/6/, Nterp/1/
C
c       Next are normalized Legendre polynomial coefficients, B1 thru B6,


c        for En = A1.  Data for En .LE. 9.83 from G. Dietze et al. in Nuclear


c        Data for Science and Technology, 6-10 Sept 1982, Antwerp.


c        Data for En = 13.9 and 15.6 deduced from angular distributions


c        measured by


c        Data for En = 14.1 MeV deduced from angular distribution measured


c        by Haight et al, Nuclear Science & Eng. 87 (1984) 41.


c        Data for En = 11.5 MeV included to bridge the gap and has no


c        experimental basis.


      Data A1/8.0, 8.64, 8.99, 9.22, 9.41, 9.83, 11.5, 13.9,14.1,15.6/
      Data B1/0.2447, 0.226, -0.1538, -0.1006, -0.1708, -0.04401,
     X 0.06217, 0.1525, 0.24192, 0.06217/
      Data B2/0.07407, -0.1602, 0.06007, 0.1265, 0.1833, 0.1678,
     X -0.058, 0.01694, 0.10631, -0.058/
      Data B3/0.04101, -0.0905, 0.00231, 0.0048, 0.02966, 0.01752,
     X -0.0324, -0.00457, 0.05494, -0.03239/
      Data B4/0.0236, 0.02001, 0.08777, 0.138, 0.1263, 0.04874,
     X 0.0497, 0.04144, 0.09454, 0.04966/
      Data B5/ 0.0, 0.02022, 0.0, 0.00375, -0.01591, -0.105,
     X 0.0259, 0.05025, 0.07588, 0.02585/
      Data B6/0.0, 0.0, 0.0, 0.01643, 0.03383, -0.0293, -0.02105,
     X -0.0102, 0.015443, -0.02105/
c
      En=Eneut
      Eneut=0.0
      Vn=VELOCITY(Kn,En)
      CALL CMTRAN(Kn,K12C,Vn, Vcom,Enc,Ecc)
      Tec=Enc+Ecc
      Ta=Tec-Q
      IF (Ta .GT. 0.0) goto 10
csatoh
c     Type 5, En
      write(*,5)En
    5 Format(/'   *** Error in Subroutine NALPHA; E(neut) at entry ='
     Y 1PE11.3/10x,'Set E(Alpha) = 0 and exit'/)
      Echrg(1)=0.0
      Echrg(2)=0.0
      Return
c
   10 Ealpha=RCKE(Ka,K9Be,Ta)
      Vac=VELOCITY(Ka,Ealpha)
c added by kajimoto 06/13
      IF (RAN(Irx) .LT. 0.04) goto 11 ! added by kajimoto 06/13
      Temp=0.06137*En+4.1 !(same 3He emission)
      Ealpha=CHOOSP(Ealpha,Temp)
      Try=Ealpha*(Emass(Ka) + Emass(K9Be))/Emass(K9Be)
      Ealpha=RCKE(Ka,K9Be,Try)
      Vac=VELOCITY(Ka,Ealpha)
c added by kajimoto 06/13
c  Now get cosine of polar scattering angle.
   11 Fi(1)=EXTERP(A1,B1,En,Ipoly,Nterp)
      Fi(2)=EXTERP(A1,B2,En,Ipoly,Nterp)
      Fi(3)=EXTERP(A1,B3,En,Ipoly,Nterp)
      Fi(4)=EXTERP(A1,B4,En,Ipoly,Nterp)
      Fi(5)=EXTERP(A1,B5,En,Ipoly,Nterp)
      Fi(6)=EXTERP(A1,B6,En,Ipoly,Nterp)
      Fmu=CHOOSL(Fi,Nfi)
c added by kajimoto 06/13
      if(En .GE. 40.0) then
        Fmu=akalbach(En,Ealpha,13.0,6.0,0.0,1.0,2.0,2.0)
      endif
c added by kajimoto 06/13
      Sinpsi=SQRT(1.0 - Fmu*Fmu)
      Phi=6.283185*RAN(Irx)
      Xa=Sinpsi*COS(Phi)
      Ya=Sinpsi*SIN(Phi)
      Za=Fmu
      Vax=Xa*Vac
      Vay=Ya*Vac
      Vaz=Za*Vac
      CALL LABTRAN(Vax,Vay,Vaz, Vcom, Vx,Vy,Vz)
      Va=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Ealfa=EFROMV(Ka,Va)
      Echrg(1)=Ealfa
ckaji
c      IF (Ealfa .LE. 0.1) goto 20
      IF (Ealfa .LE. 0.55) goto 20
c          Eprot < 100 keV, range a few micrometers.
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(U,V,W)
      Alphpl=PLNGTH(Xn,Yn,Zn,X,Y,Z, Rdet,Ht)
      Ralpha=EXTERP(RNGEN,ARAN,Ealfa,Noprot,4)
c       Now compare alpha's range in detector with the path length
c         of the vector having direction cosines Xn, Yn, Zn from position
c         X, Y, Z.
      RangeQ=Ralpha-Alphpl
      IF (RangeQ .LE. 3.5000E-04) goto 20
c         If it goes to next step, the proton got out of the detector.
      Echrg(1)=EXTERP(ARAN,RNGEN,RangeQ,Noprot,4) ! kajimoto 2011/03/18
      if(iesc.eq.2) Echrg(1)=0.0 ! kajimoto 2011/03/18
   20 continue
ckaji end
      E9Bec=Ta-Ealpha
      V9Bec=VELOCITY(K9Be,E9Bec)
      V9x=-Xa*V9Bec
      V9y=-Ya*V9Bec
      V9z=-Za*V9Bec
      CALL LABTRAN(V9x,V9y,V9z, Vcom, Vx,Vy,Vz)
      V9Be=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      E9Be=EFROMV(K9Be,V9Be)
      Echrg(2)=E9Be
ckaji 04/14/2011
      if(iesc.eq.2.and.Echrg(1).eq.0.0) Echrg(2)=0.0
ckaji end
      Nelm=4
      CALL BANKR2
      Return
      END
c------------------------------------------------------------------
C   This is file NN3ALF.FOR
c
c     Purpose is to compute the reaction  n + 12-C --> n' + 3 alphas
c       At exit should have the n' energy and its direction cosines in


c       the appropriate variables in the NEUTRN Common area, and the


c       energies of the 3 alphas in the array Echrg in the


c       COLLIS Common area.


c
c       The calculation follows one of several reaction schemes:


c
c    1. (a)        n + 12-C  -->  n' + 12-C (excited)
c       (b)  12-C (excited)  -->  alpha + 8-Be (ground state)
c       (c)  8-Be (grnd st)  -->  2 alphas.
c
c    2. (a)  same as 1. (a)
c       (b)  12-C (excited)  -->  alpha + 8-Be (excited state at 3.0 MeV)
c       (c)  8-Be (excited)  -->  2 alphas.
c
c    3. (a)  same as 1. (a)
c       (b)  12-C (excited)  -->  3 alphas via 3-body breakup
c
c    4. (a)        n + 12-C  -->  alpha + 9-Be (excited)
c       (b)  9-Be (excited)  -->  n + 8-Be (ground state)
c       (c)  same as 1. (c)
c
c    5. (a)  same as 4. (a)
c       (b)  9-Be (excited)  -->  n + 8-Be (excited state at 3.0 MeV)
c       (c)  same as 2. (c)
c
c    6. (a)  same as 4. (a)
c       (b)  9-Be (excited)  -->  alpha + 5-He
c       (c)  5-He            -->  n + alpha
c
c         Added 5/87:


c    7. (a)  same as 4. (a)
c       (b)  9-Be (excited!) --> p + 8-Li
c       (c)  8-Li --> several modes.
c
      Subroutine NN3ALF
c
      Real*4 Li8pp, Li7pn
      Common /RANDM/ Irx
      Common /P3ALF/ Pbar(13)
      Common /NEUTRN/ Eneut, Vnprx, Vnpry, Vnprz, X, Y, Z
      Common /NEUTR2/ Eneut2, U2,V2,W2, X2,Y2,Z2
      Common /COLLIS/ Nelm, Echrg(6)
      Common /VECTOR/ Xb,Yb,Zb, Xn,Yn,Zn
      Common /MASSES/ Emass(21)
      Common /PPLI8/  Nx9,Ex9(14),Pprot(14)
ckaji 2011/03/18
      common /RES/ iswt, bias(5), icont, ictrl, iesc, ilight
c
      Dimension Q(14)
c
      Data Rmass/0.922461/, Rmass2/0.666670/
      Data K5He/11/, K9Be/4/, Kn/1/, KC/7/, K8Be/8/, Ka/3/
      Data TenMeV/10.0/, Q5He/0.88/, Qcont/17.9/
      Data Q/7.656,9.641,8.13,8.54,10.84,10.4,11.84,12.46,14.08,
     X   16.1,0.0,16.99,19.7,0.0/
      Data Qnn/7.3665/, E8Bex/3.0/, QBa/2.46/, QBn/1.666/
      Data E8Begs/0.092/, Qaa/22.5/, Qna/5.71/
C
c       Next arrays for determining proton branching following


c         n + 12-C --> alpha + 9-Be.


      Data Ex9/ 16.89, 17.11, 17.33, 18.1, 18.7, 18.94, 19.7,
     a  19.94, 21.1, 26.0, 32.0, 38.6, 44.0, 51.3/
      Data Pprot/ 0.0, .000036, .0134, .0484, .1226, .0965, .156,
     a  .1237, .2,   .22,  .175, .18,   .2,  .25/
      Data Nx9/14/, Nterp/1/, K8Li/16/, Li8pp/16.888/, Li7pn/2.033/
      Data Kp/2/
c
c
c       Initial set up.


      Do 1 I=1,6
    1 Echrg(I)=0.0
      Egamma=0.0
c
c       First step: transform to center-of-mass coordinates


      En=Eneut
      Vn=VELOCITY(Kn,En)
      CALL CMTRAN(Kn,KC,Vn, Vcom,Enc,Ecc)
      Tec=Enc+Ecc
      Ta=Tec-Q(1)
      IF (Ta .GT. 0.0) goto 10
csatoh
c     Type 5, En
      write(*,5)En
    5 Format(/'   *** Error in Subroutine NN3ALF; E(neut) at entry ='
     & 1PE11.3/10x,'Set E(neut) = E(all alphas) = 0 and exit'/)
      Eneut=0.0
      Return
c               Last is an "error" return.


c
c         Next is to determine which of 14 reactions to use in rest of


c           computation.


   10 CALL P3ALPH(En)
c          That sets up the -Pbar- array of the P3ALF Common area


c         Choose branching by random number from this array.


      Pran=RAN(Irx)
      Nbrnch=1
      Do 12 J=1,13
      IF (Pran .LE. Pbar(J)) goto 14
   12 Nbrnch=Nbrnch+1
   14 Continue
c
c      Now have Nbrnch between 1 and 14.  The chosen reaction is
c               as follows:


c
c       Nbrnch    From (n,n')   or   From (n,alpha)


c       ------    -----------        --------------


c          1      Ex = 7.65 MeV
c          2          9.64
c          3                         Ex = 2.43 MeV


c          4                           2.8 + 3.05


c          5         10.84


c          6                              4.70


c          7      11.8+12.7+13.3


c          8                              6.76


c          9         14.08


c         10      16.1 - 18.0


c         11       Continuum
c         12                           11.5 group


c         13                           14.0 group


c         14                           Continuum


c
c
      Goto (20,20,40,40,20,40,20,40,20,20,15,40,40,38), Nbrnch
c
c       Next step is for the 12-C continuum, i.e. excitation of a 12-C


c         "excited state" having Ex > 18 MeV.  We get an approximate


c         n' energy from Function CHOOSN, and then get from that an


c         effective excitation energy for a "level" in the continuum


c         presumably excited in the reaction.


c
   15 Ta=Tec-Qcont
      Encom=RCKE(Kn,KC,Ta)
c               That's the maximum energy for the "continuum" neutron.


      F=0.065+0.001*En
      Temp=F*En
c               -Temp- determined empirically so that the CHOOSN function


c                 deals with a neutron continuum spectrum similar to that


c                 deduced by the nuclear model code TNG.


      Elow=0.0
      Entry=CHOOSN(Elow,Encom,Temp)
      Try=Entry*(Emass(Kn) + Emass(KC))/Emass(KC)
      Q(11)=Tec-Try
c
c         Next step is to get the outgoing neutron energy and direction


c          cosines in detector (laboratory) coordinates.


   20 Ta=Tec-Q(Nbrnch)
      Ena=Rmass*Ta
      IF (Ta .GT. TenMeV) Ena=RCKE(Kn,KC,Ta)
      Eca=Ta-Ena
      Vnca=VELOCITY(Kn,Ena)
      Vcca=VELOCITY(KC,Eca)
c               Choose neutron velocity direction in the center-of-mass


c                coordinates by random number.  An isotropic


c                distribution is assumed.


      CALL RVECT(Xb,Yb,Zb)
      Zx=Vnprx
      Zy=Vnpry
      Zz=Vnprz
c          Get neutron velocity components in "neutron" center-of-mass


c               coordinates:


      Vnxp=Xb*Vnca
      Vnyp=Yb*Vnca
      Vnzp=Zb*Vnca
c         Get carbon velocity components in "neutron" c.o.m. coordinates.


      Vcxp=-Xb*Vcca
      Vcyp=-Yb*Vcca
      Vczp=-Zb*Vcca
c
c          Transform neutron velocity components into "neutron" laboratory


c               coordinates:


      CALL LABTRAN(Vnxp,Vnyp,Vnzp, Vcom, Vnwx,Vnwy,Vnwz)
      Vneut=SQRT(Vnwx*Vnwx + Vnwy*Vnwy +Vnwz*Vnwz)
c               Now can get E(Neut) and dir. cosines in lab coordinates


c                 for the outgoing neutron.


      Eneut=EFROMV(Kn,Vneut)
      CALL DIRCOS(Vnwx,Vnwy,Vnwz, Xb,Yb,Zb)
c           Xb,Yb,Zb are dir. cosines in "neutron" laboratory coordinates.


c               Now get dir. cosines in "detector" laboratory coordinates:


      CALL TRANSVEC(Zx,Zy,Zz)
      Vnprx=Xn
      Vnpry=Yn
      Vnprz=Zn
c       That takes care of the outgoing neutron, n', of reaction 1. (a)


c
c         Next step: the "fission" of the excited carbon ion.


c           The center-of-mass for this step is


c           the moving excited carbon ion.


      CALL LABTRAN(Vcxp,Vcyp,Vczp, Vcom, Vcx,Vcy,Vcz)
      Vcarb=SQRT(Vcx*Vcx + Vcy*Vcy + Vcz*Vcz)
c
c         Now test to see if ground state or excited state of 8-Be is


c               involved.


      Ta=Q(Nbrnch)-Qnn-E8Begs
      N8Bex=0
      IF (Nbrnch-5) 25,17,18
   17 IF (RAN(Irx) - 0.6) 25,25,19
   18 IF (Nbrnch .EQ. 11) goto 25
   19 Ta=Ta-E8Bex
      N8Bex=1
c         N8Bex=0 means ground state involved; N8Bex=1 means excited state.


   25 Eaa=Rmass2*Ta
      IF (Ta .GT. TenMeV) Eaa=RCKE(Ka,K8Be,Ta)
c          Check to see if a contiuum reaction -- if so test for the


c               3-body breakup mode -- if so get an alpha energy (by


c               random number in the Function CHOOSA).


      L3body=0
      IF (Nbrnch .LT. 11) goto 22
      Pcomp=(En - 20.0)/50.
c        The assumption of variable Pcomp: a simple increase in the
c         probability of 3-body breakup with increasing E(neutron)


      IF (RAN(Irx) .GE. Pcomp) goto 22
c           If it gets to here it's a 3-body breakup reaction.


      L3body=1
      Emax=Eaa
      Eaa=CHOOSA(Emax)
   22 Eba=Ta-Eaa
      Vaa=VELOCITY(Ka,Eaa)
c         The next portion gets the energy of the alpha from the fission


c       reaction 1. (b) or 2. (b) or of the the first alpha from the 3-body


c       breakup reaction 3. (b).


c
c          The angular distribution of the "fission" products is, perforce,


c               randomly chosen.


      CALL RVECT(Xc,Yc,Zc)
c          Velocities are slow enough for these heavier ions to do the


c               transformations non-relativistically.


      Va1x=Xc*Vaa
      Va1y=Yc*Vaa
      Va1z=Zc*Vaa + Vcarb
c       Va1x, Va1y, Va1z are components of alpha velocity in the "carbon-ion"


c           laboratory coordinate system.  Nonrelativistic transformation.


      Valph1=SQRT(Va1x*Va1x + Va1y*Va1y + Va1z*Va1z)
c       Valph1 is the same in the "detector" laboratory coordinate system as


c               it is in the "carbon-ion" laboratory coordinate system.


      Ealph1=EFROMV(Ka,Valph1)
      Echrg(1)=Ealph1
c         and that's it for the "first" alpha particle.


c
c    Now recheck for the 3-body breakup reaction again.
      IF (L3body .EQ. 0) goto 26
c
c   Next is for 3-body breakup reaction 3. (b).  Get 2d alpha information.
      Ea2=0.5*Eba
      Va2=VELOCITY(Ka,Ea2)
      Vopp=0.5*Vaa
c               (That's the component of velocity of the 2d and 3rd


c                alphas along a z-axis defined by the velocity of the


c                first alpha in the 3-body breakup.)


      Dirzet=Vopp/Va2
      Sintheta=SQRT(1.0-Dirzet*Dirzet)
      Phi=3.1415926*Ran(Irx)
c               For the 2d alpha; the 3rd will scatter at Phi + Pi


      Sinphi=SIN(Phi)
      Cosphi=COS(Phi)
c          Need to transform coordinates from those defined by the z-axis


c           determined by the velocity vector of the first alpha to the


c           "carbon-ion" coordinates.


      Xb=Sintheta*Cosphi
      Yb=Sintheta*Sinphi
      Zb=Dirzet
      Zx=-Xc
      Zy=-Yc
      Zz=-Zc
      CALL TRANSVEC(Zx,Zy,Zz)
c        Xn,Yn,Zn are now direction cosines of 2d alpha in "carbon-ion"


c          center-of-mass coordinates.


      Va2x=Va2*Xn
      Va2y=Va2*Yn
      Va2z=Va2*Zn + Vcarb
c               (Non-relativistic coord. transformation, again.)


      Valph2=SQRT(Va2x*Va2x + Va2y*Va2y + Va2z*Va2z)
      Ealph2=EFROMV(Ka,Valph2)
      Echrg(2)=Ealph2
c         And that's it for the 2d alpha.  Now for the third alpha --


c               its velocity, Va3, equals Va2 -- but in a


c               different direction


      Phi=Phi+3.1415926
      Sinphi=SIN(Phi)
      Cosphi=COS(Phi)
c         Get dir. cosines of 3d alpha in "carbon-ion" c.o.m. coordinates.


      Xb=Sintheta*Cosphi
      Yb=Sintheta*Sinphi
c         Zb is the same as for 2d alpha; so are Zx, Zy, and Zz.


      CALL TRANSVEC(Zx,Zy,Zz)
c         Now get velocity components of 3d alpha in "carbon-ion" lab. coords.


      Va3x=Va2*Xn
      Va3y=Va2*Yn
      Va3z=Va2*Zn + Vcarb
      Valph3=SQRT(Va3x*Va3x + Va3y*Va3y + Va3z*Va3z)
      Goto 60
c  End of 3-body breakup energetics and kinematics.
c
c      Next branch is "fission" of the moving 8-Be ion, 1. (b) or 2. (b).
c         The center of mass is the moving 8-Be ion.


   26 Vba=VELOCITY(K8Be,Eba)
      Vbwx=-Xc*Vba
      Vbwy=-Yc*Vba
      Vbwz=-Zc*Vba + Vcarb
      VBe = SQRT(Vbwx*Vbwx + Vbwy*Vbwy + Vbwz*Vbwz)
c         Now recall whether Ground state or excited state of 8-Be


c          is involved.


   21 Ta=E8Begs
      IF (N8Bex .EQ. 1) Ta=Ta+E8Bex
   28 Ea2=0.5*Ta
   29 Va2=VELOCITY(Ka,Ea2)
c          Again, perforce, isotropy in angular distribution of the two


c            alphas in the center-of-mass coordinates.


      CALL RVECT(Xc,Yc,Zc)
      Va2x=Xc*Va2
      Va2y=Yc*Va2
      Va2z=Zc*Va2+Vbe
      Valph2=SQRT(Va2x*Va2x + Va2y*Va2y + Va2z*Va2z)
      Ealph2=EFROMV(Ka,Valph2)
      Echrg(2)=Ealph2
c         and that's it for the second alpha.


c
c         Va3x=-Va2x and Va3y=-Va2y


      Va3z=-Zc*Va2+Vbe
      Valph3=SQRT(Va2x*Va2x + Va2y*Va2y + Va3z*Va3z)
      Goto 60
c  All done there, except to tidy up at very end.
c
c
c   Next sections are for the reactions  n + 12-C --> alpha + 9-Be.
c
c       Next 5 statements obtain an effective 9-Be "continuum" energy


c       level very approximately.  Two rationalizations: (1) probability


c       for this branch is small, a few percent of the  n + 3 alpha


c       reaction; and (2) the "continuum" is not known experimentally,


c       and the calculated "continuum" using TNG is very similar to that


c       calculated for the proton "continuum".


   38 Ta=Tec-Qaa
      Eamax=RCKE(Ka,K9Be,Ta)
      Temp=4.0
      Eatry=CHOOSP(Eamax,Temp)
      Try=Eatry*(Emass(Ka) + Emass(K9Be))/Emass(K9Be)
      Etry=RCKE(Ka,K9Be,Try)
      Try=Try*Eatry/Etry
      Q(14)=Tec-Try
c
c       Before continuing on alpha + 9-Be, determine if subsequent 9-Be


c         decay will be by proton emission.


      Npgo=0
      Excit9=Q(14)-Qna
      IF (Excit9 .LE. Ex9(1)) goto 40
      Prob=EXTERP(Ex9,Pprot,Excit9,Nx9,Nterp)
      IF (RAN(Irx) .LE. Prob) Npgo=1
c       OK.  Npgo=1 means proton decay of 9-Be; Npgo = 0 means alpha decay


c         of 9-Be.  In either case first work on  alpha + 9-Be energetics.


c       Get information on the alpha.


   40 Ta=Tec-Q(Nbrnch)
      Ealpha=RCKE(Ka,K9Be,Ta)
      Va1=VELOCITY(Ka,Ealpha)
      CALL RVECT(Xa,Ya,Za)
      Va1x=Xa*Va1
      Va1y=Ya*Va1
      Va1z=Za*Va1
c       Now transform to neutron laboratory coordinates.


      CALL LABTRAN(Va1x,Va1y,Va1z, Vcom, Vx,Vy,Vz)
c       ... done.  Vx, Vy, Vz are velocity components of alpha in the


c          neutron lab. coordinates.


      Va=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Ea1=EFROMV(Ka,Va)
      Echrg(1)=Ea1
c          That's it for the (first) alpha.


c
c       Now get 9-Be motion into detector coordinates.  We have to get back


c         to detector coordinates because we still have to get information


c         on the neutron in those coordinates.


      E9Be=Ta-Ealpha
      V9Bec=VELOCITY(K9Be,E9Be)
      V9xc=-Xa*V9Bec
      V9yc=-Ya*V9Bec
      V9zc=-Za*V9Bec
      CALL LABTRAN(V9xc,V9yc,V9zc, Vcom, V9x,V9y,V9z)
      V9Be=SQRT(V9x*V9x + V9y*V9y + V9z*V9z)
      CALL DIRCOS(V9x,V9y,V9z, Xb,Yb,Zb)
      Zx=Vnprx
      Zy=Vnpry
      Zz=Vnprz
      CALL TRANSVEC(Zx,Zy,Zz)
c       Now have 9-Be in detector laboratory coordinates; its direction


c        cosines are  Xn, Yn, and Zn.


c
c       First decide if next reaction is reaction 7. (b).


      IF (Npgo .EQ. 1) goto 75
c       Program counter to here, then not 7. (b), so try for 6. (b).


      IF (Nbrnch .EQ. 3) goto 48
      IF (Nbrnch .EQ. 6  .AND.  Ran(Irx) .LT. 0.87) goto 50
      IF (Nbrnch .EQ. 8  .AND.  Ran(Irx) .LT. 0.45) goto 50
      IF (Nbrnch .EQ. 12 .AND.  Ran(Irx) .LT. 0.88) goto 50
      IF (Nbrnch .EQ. 13 .AND.  Ran(Irx) .LT. 0.60) goto 50
c       If program gets to here it wasn't  6. (b).


c       Choice of  4. (b)  or  5. (b)  depends only on the 9-Be level.


      Qq=Qnn+E8Begs
c       So, next is "fission" of 9-Be --> n + 8-Be.  Determine which


c                8-Be state is involved.


      N8Bex=0
      IF (Nbrnch .LT. 8) goto42
      N8Bex=1
      Qq=Qq+E8Bex
   42 Ta=Q(Nbrnch)-Qq
c        Get information on the neutron first.


      Enn=RCKE(Kn,K8Be,Ta)
      Vnn=VELOCITY(Kn,Enn)
      CALL RVECT(Xa,Ya,Za)
      Vnx=Xa*Vnn
      Vny=Ya*Vnn
      Vnz=Za*Vnn
c       Transform velocity components to  9-Be  laboratory coordinates.


      CALL LABTRAN(Vnx,Vny,Vnz, V9Be, Vx,Vy,Vz)
c         ... done.  Vx, Vy, Vz are the desired velocity components.


      Vneutr=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Eneut=EFROMV(Kn,Vneutr)
      CALL DIRCOS(Vx,Vy,Vz, Xb,Yb,Zb)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
c        and now have neutron's direction cosines in detector coordinates


      Vnprx=Xn
      Vnpry=Yn
      Vnprz=Zn
c         so that takes care of the neutron for reactions  4. (b) and 5. (b)


c       Now need to get motion of 8-Be in the 9-Be lab. coordinates.


      E8Be=Ta-Enn
      V8Be=VELOCITY(K8Be,E8Be)
      V8x=-Xa*V8Be
      V8y=-Ya*V8Be
      V8z=-Za*V8Be + V9Be
      VBe=SQRT(V8x*V8x + V8y*V8y + V8z*V8z)
c       Rest of computation same as for 1. (c) or 2. (c)


      Goto 21
c
c     Coming down the home stretch!!  Next is for  9-Be --> alpha + 5-He
c       reaction  6. (b).


c        First look at case where the Ex=2.43 MeV state of 9-Be decays into


c        alpha + 5-He.  This is clearly an unusual decay mode, since the Q-


c        value for decay of 9-Be into alpha + 5-He is 2.46 MeV, or some 30


c        keV larger than Ex.  However both the level in 9-Be and the ground


c        state in 5-He are quite broad, and the decay mechanism must take


c        these widths into account.  Programming in such widths seems a bit


c        much for the present purpose.  Instead we will simply note that on


c        the average there is no residual kinetic energy in the 9-Be center


c        of mass to be shared by the alpha and the 5-He ion after the decay.


c        Hence V(alpha2) = V(5-He) = V(9-Be) in the lab. coordinates


   48 Ea2=EFROMV(Ka,V9Be)
      Echrg(2)=Ea2
      V5He=V9Be
c         Also the dir. cosines of the 5-He ion in the detector lab. coords.


c          are the same as the dir. cosines of the 9-Be (i.e., Xn,Yn,Zn).


      Goto 55
c
c       Decay of other 9-Be levels will involve some c.o.m. kinetic energy.


   50 Qq=Qnn+QBa-QBn
      Ta=Q(Nbrnch)-Qq
      Ea2c=RCKE(Ka,K5He,Ta)
c       That's the energy of the 2d alpha in the 9-Be center of mass.


      Va2c=VELOCITY(Ka,Ea2c)
      CALL RVECT(Xa,Ya,Za)
c         Again, isotropy in the decay process in the center of mass.


      Va2x=Xa*Va2c
      Va2y=Ya*Va2c
      Va2z=Za*Va2c + V9Be
c         and again, a non-relativistic transformation to lab. coordinates.


      Va2=SQRT(Va2x*Va2x + Va2y*Va2y + Va2z*Va2z)
      Ea2=EFROMV(Ka,Va2)
      Echrg(2)=Ea2
c         So -- finished with alpha from 9-Be decay.  Now set up the motion


c          of the 5-He ion in detector lab. coordinates (remember, we still


c          have the neutron's parameters to get into detector lab coords.).


      E5Hec=Ta-Ea2c
      V5Hec=VELOCITY(K5He,E5Hec)
      V5Hex=-Xa*V5Hec
      V5Hey=-Ya*V5Hec
      V5Hez=-Za*V5Hec + V9Be
c         and once more a non-relativistic transformation, this time to get


c           the motion of the 5-He ion in the 9-Be lab. coordinates.


      V5He=SQRT(V5Hex*V5Hex + V5Hey*V5Hey + V5Hez*V5Hez)
      CALL DIRCOS(V5Hex,V5Hey,V5Hez, Xb,Yb,Zb)
c        Recall 9-Be ion has dir. cosines Xn,Yn,Zn in detector lab. coords.


      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
C       Now we have the 5-He ion dir. cosines in detector lab. coordinates;


c         the direction cosines are Xn, Yn, Zn.  So now to the decay of the


c         5-He into a neutron and an alpha.  Get neutron information first.


   55 Ta=Q5He
      Ena=RCKE(Kn,Ka,Ta)
      Vnc=VELOCITY(Kn,Ena)
      CALL RVECT(Xa,Ya,Za)
c          once more, isotropy in the decay.


      Vnxc=Xa*Vnc
      Vnyc=Ya*Vnc
      Vnzc=Za*Vnc
c        Transform neutron's velocity components from 5-He c.o.m. coords.


c         to 5-He lab. coords.


      CALL LABTRAN(Vnxc,Vnyc,Vnzc, V5He, Vx,Vy,Vz)
c        Got them.  Vx,Vy,Vz are neutron's velocity components in 5-He c.o.m.


      Vneut=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Eneut=EFROMV(Kn,Vneut)
c        Eneut=E(neutron) in lab. coords.  Now get dir. cosines in detector


c          lab. coordinates.


      CALL DIRCOS(Vx,Vy,Vz, Xb,Yb,Zb)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
      Vnprx=Xn
      Vnpry=Yn
      Vnprz=Zn
c         Whew!!  That's it for the neutron.  Now for the final alpha.


      Ea3c=Ta-Ena
      Va3c=VELOCITY(Ka,Ea3c)
      Va3x=-Xa*Va3c
      Va3y=-Ya*Va3c
      Va3z=-Za*Va3c
      CALL LABTRAN(Va3x,Va3y,Va3z, V5He, Vx,Vy,Vz)
      Valph3=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
c
c   All that's left is to finish up on the third alpha.
   60 Ealph3=EFROMV(Ka,Valph3)
      Echrg(3)=Ealph3
      Goto 90
c       END of (n,n 3-alpha) reaction computation.


c
c       Next portion of programming for 9-Be --> p + 8-Li decay.


c         Set outgoing neutron energy to zero to begin with.


   75 Eneut=0.0
      Egamma=0.0
      Eleftov=0.0
      CALL PP8LI(Excit9,V9Be,Egamma)
c
c       That's it.  At the end of this section the Echrg array should be:


c               Echrg(1) = First alpha from n + 12-C --> alpha + 9-Be


c               Echrg(2) = Escaped proton information, if it escaped.


c               Echrg(4) = Proton from 9-Be --> p + 8-Li


c               Echrg(3), (5), (6) from routine ELI8DK.


c
c       Final bit of business.


   90 Nelm=5
      CALL BANKR2
      Eneut2=ABS(Eneut2)
      IF (Egamma .GT. 0.0) CALL PHOTON(Egamma)
      Return
      END
c------------------------------------------------------------------
C   This is file N3HE.FOR
c
c       Purpose is to compute energetics for reaction


c         n + 12-C --> 3-He + 10-Be, with possible subsequent decay


c             10-Be --> n + 9-Be.            As of 6/87.


c       Added in 9-Be decay modes -- 8/87.


c
      Subroutine N3HE
c
      Common /COLLIS/ Nelm,Echrg(6)
      Common /NEUTRN/ Eneut, U,V,W, X,Y,Z
      Common /MASSES/ Emass(21)
      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
      Common /RANDM/  Irx
      Common /NAID/   Rad, Ht, Rr1, Rr2
      Common /EXCTBE/ Ex10Be(4),Be9pn
c       Array Ex10Be is in subroutine NPX.


      Common /NEUTR2/ Eneut2, U2,V2,W2, X2,Y2,Z2
      Common /PPLI8/  Nx9,Ex9(14),Pp(14)
      Common /BENINE/ Exc9Be(6),Be8pn
ckaji 2011/03/18
      Common /PROTON/NOPROT,RNGEN(144),PRAN(144),DRAN(144),ARAN(144),
     a              HRAN(144),TRAN(144)
      common /RES/ iswt, bias(5), icont, ictrl, iesc, ilight
c       See subroutine NN3ALF for values of variables in -PPLI8- common.


c
      Dimension Epoly(3),B1(3),B2(3),B3(3),B4(3),B5(3),B6(3),Fi(6)
c
      Data Kn/1/, K12C/7/, K3He/19/, K10Be/15/, K9Be/4/, K8Be/8/
      Data Q/19.47/, Nfi/6/, Be8pn/1.6652/, Be8gs/0.092/, Ka/3/
c   Next array represents 9-Be level structure well enough for present use.
      Data Exc9Be/0.0, 2.9, 4.7, 6.8, 11.5, 14.0/
c
c       Next arrays will be used to determine angular distributions of


c         outgoing 3-He ions using Legendre coefficients determined from


c         study of the data provided by Subramanian, et al, Physical


c         Review C28, 521 (1983).


      Data Ipoly/3/, Epoly/27.4, 39.7, 60.7/
      Data B1/0.4204, 0.4207, 0.3993/,  B2/0.209, 0.1761, 0.1569/
      Data B3/0.0925, 0.046,  0.09874/, B4/0.033, 0.0234, 0.04945/
      Data B5/0.001,  0.00243,0.01223/, B6/0.004, 0.01271,0.00558/
c
c   Start computation.  Initialize variables.
      Do 1 I=1,6
    1 Echrg(I)=0.0
      En2=0.0
c       Go to center-of-mass and check E(neutron)


      En=Eneut
      Eneut=0.0
      Vn=VELOCITY(Kn,En)
      CALL CMTRAN(Kn,K12C,Vn, Vcom,Enc,Ecc)
      Tec=Enc+Ecc
      Ta=Tec-Q
      IF (Ta .GT. 0.0) goto 2
csatoh
c     Type 100, En
      write(*,100)E
  100 Format(/'   ***  Error in Subroutine N3HE, En = '1PE11.4/
     x  10x,'Set E(neutron) = 0.0 and Exit... ')
      Return
c       Last is an "error" return.


c
c       Set up Legendre coeff. computation.


    2 Fi(1)=EXTERP(Epoly,B1,En,Ipoly,1)
      Fi(2)=EXTERP(Epoly,B2,En,Ipoly,1)
      Fi(3)=EXTERP(Epoly,B3,En,Ipoly,1)
      Fi(4)=EXTERP(Epoly,B4,En,Ipoly,1)
      Fi(5)=EXTERP(Epoly,B5,En,Ipoly,1)
      Fi(6)=EXTERP(Epoly,B6,En,Ipoly,1)
c
c       Get center-of-mass energies:


      Ehcom=RCKE(K3He,K10Be,Ta)
      Ehe=Ehcom
      Excit=0.0
      Egamma=0.0
      Eleftov=0.0
      Level=1
      IF (Ehcom .LT. Ex10Be(2)) goto 8
      IF (RAN(Irx) .LT. 0.04) goto 8
c         Last is to slightly enhance ground-state population of 10-Be.


      Level=5
ckaji      Temp=3.0
      Temp=0.06137*En+4.1 ! modified by kajimoto 06/08
c         Last is an assumption, and probably not very critical since the


c           3-He production is small.


      Ehe=CHOOSP(Ehcom,Temp)
      Try=Ehe*(Emass(K3He) + Emass(K10Be))/Emass(K10Be)
      Etry=RCKE(K3He,K10Be,Try)
      Try=Try*Ehe/Etry
c         As usual, one iteration on -Try- for improved precision.


      Excit=Ta-Try
      IF (Excit .GT. Be9pn) goto 6
c       If -Excit- is < 9-Be + n, then pair -Excit- with "level" in 10-Be.


      K=1
      Do 3 I=2,4
      IF (Excit .LT. Ex10Be(I)) goto 4
    3 K=K+1
    4 Excit=Ex10Be(K)
      Level=K
c       See subroutine NPX, comments following statement no. 42, for decay


c        characteristics of 10-Be levels.


      IF (K.GE.2) Egamma=Ex10Be(2)
      IF (K.EQ.3 .AND. RAN(Irx).GE.0.5) Egamma=Excit
      Eleftov=Excit-Egamma
c
    6 Taex=Ta-Excit
      Ehe=RCKE(K3He,K10Be,Taex)
c
    8 Qex=Q+Excit
c       Now get 3-He scattering angle by random choice from Legendre


c         coefficients tabulated above.

      Ctheta=CHOOSL(Fi,Nfi)
c added by kajimoto 06/13
      if(En .GE. 40.0) then
        Ctheta=akalbach(En,Ehe,13.0,6.0,0.0,1.0,2.0,1.0)
      endif
c added by kajimoto 06/13
c       Done.  Now choose azimuthal angle and get direction cosines of 3-He.


      Rann=6.283185*RAN(Irx)
      Sinphi=SIN(Rann)
      Cosphi=COS(Rann)
      Sintheta=SQRT(1.0-Ctheta*Ctheta)
      Cx=Sintheta*Cosphi
      Cy=Sintheta*Sinphi
      Cz=Ctheta
c       Done.  Get velocity components of 3-He ion in center-of-mass coords.


      Vhe=VELOCITY(K3He,Ehe)
      Vhx=Cx*Vhe
      Vhy=Cy*Vhe
      Vhz=Cz*Vhe
c       Done.  Transform to velocity components in laboratory coords.


      CALL LABTRAN(Vhx,Vhy,Vhz, Vcom, Vx,Vy,Vz)
c       Done.  Get Velocity and Energy of 3-He ion in laboratory coords.


      V3He=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      E3He=EFROMV(K3He,V3He)
      Echrg(1)=E3He
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(U,V,W)
c       Xn,Yn,Zn are the proton's direction cosines in lab coordinates.
      Protpl=PLNGTH(Xn,Yn,Zn, X,Y,Z, Rad,Ht)
      Rprot=EXTERP(RNGEN,HRAN,E3He,Noprot,4)
c       Now compare proton's range in detector (RPROT) with the path length
c         of the vector having direction cosines Xn, Yn, Zn from position
c         X, Y, Z (PROTPL).
      RangeQ=Rprot-Protpl
ckaji      IF (RangeQ .LE. 0.1) goto 99
      IF (RangeQ .LE. 0.2) goto 99
c         If it goes to next step, the proton got out of the detector.
      E3He=EXTERP(HRAN,RNGEN,RangeQ,Noprot,4) ! kajimoto 2011/03/18
      Echrg(1)=E3He
      if(iesc .eq. 2) Echrg(1)=0.0
c       That's the energy of the 3-He ion.  Now study 10-Be ion motion.

   99 continue
      Taex=Tec-Qex
      EBe=RCKE(K10Be,K3He,Taex)
      VBe=VELOCITY(K10Be,EBe)
      Vbx=-Cx*VBe
      Vby=-Cy*VBe
      Vbz=-Cz*Vbe
      CALL LABTRAN(Vbx,Vby,Vbz, Vcom, Vx,Vy,Vz)
      V10Be=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
c       Check for 10-Be decay.  If not, save energy of 10-Be ion and exit.


      IF (Level .LE. 4) goto 20
c       Program counter to here, study 10-Be --> n + 9-Be.


c
c       First rotate 10-Be motion into "detector" lab. coordinates.


      CALL DIRCOS(Vbx,Vby,Vbz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(U,V,W)
c       Done.  Dir. cosines of 3-He motion are Xn,Yn,Zn in "detector"


c         laboratory coordinate system.


      Ta=Excit-Be9pn
      Enmax=RCKE(Kn,K9Be,Ta)
      Enn=Enmax
      Ehat=En-Q
      Fn=0.065+0.001*Ehat
      Temp=Fn*Ehat
      Elow=0.0
      Entry=CHOOSN(Elow,Enmax,Temp)
      Try=Entry*(Emass(Kn) + Emass(K9Be))/Emass(K9Be)
      Etry=RCKE(Kn,K9Be,Try)
      Try=Try*Entry/Etry
      Excit=Ta-Try
c       At this point check for possible proton decay.


      Npgo=0
      IF (Excit .LE. Ex9(1)) goto 10
      Prob=EXTERP(Ex9,Pp,Excit,Nx9,1)
      IF (RAN(Irx) .GT. Prob) goto 10
      Npgo=1
c       Program to here, have  9-Be --> p + 8-Li


      Goto 14
c
c       Next check for 2d neutron already waiting to be processed.  If


c         so, then cannot have further neutron emission from 9-Be, so


c         go directly to ground state of 9-Be.


   10 IF (Eneut2 .GT. 0.0) goto 15
c       Program counter to here, then have  9-Be --> n + 8-Be


c               followed by                 8-Be --> 2 alphas.


c       Restrict 8-Be to ground state for simplicity -- justification is


c         small cross sections just to get this far, so won't be very


c         many of these events.


      Npgo=-1
      IF (Excit .GE. 16.0) goto 14
      K=1
      Do 12 I=2,6
      IF (Excit .LT. Exc9Be(I)) goto 13
   12 K=K+1
   13 Excit=Exc9Be(K)
c       Check for  -EXCIT-  =0, i.e., ground state of 9-Be.


      IF (K .EQ. 1) Npgo=0
      Try=Ta-Excit
   14 Enn=RCKE(Kn,K9Be,Try)
   15 Vn=VELOCITY(Kn,Enn)
C          Isotropy of  n + 9-Be  in 10-Be center-of-mass coordinates.


      CALL RVECT(Xp,Yp,Zp)
      Vnx=Xp*Vn
      Vny=Yp*Vn
      Vnz=Zp*Vn
c       Transform velocity components to 10-Be lab. coords.


      CALL LABTRAN(Vnx,Vny,Vnz, V10Be, Vx,Vy,Vz)
c       Done.  Get velocity, energy of neutron in 10-Be lab. coords.


      Vnlab=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Enlab=EFROMV(Kn,Vnlab)
      Eneut=Enlab
c       Done.  Rotate neutron velocity into "detector" lab. coords.


      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Z1=Xn
      ZB=Yn
      Z3=Zn
      CALL TRANSVEC(Z1,ZB,Z3)
c         Done.  Put new dir. cosines into /NEUTRN/ Common variables.


      U=Xn
      V=Yn
      W=Zn
c       (Same position coordinates as for  n + 12-C assumed.)


c       Now get information on energy of 9-Be to save.


      IF (Npgo .EQ. 0) Try=Ta
      E9Bec=Try-Enn
      V9Bec=VELOCITY(K9Be,E9Bec)
      V9x=-Xp*V9Bec
      V9y=-Yp*V9Bec
      V9z=-Zp*V9Bec
c       Transform those c.o.m. components to 10-Be laboratory coordinates.


      CALL LABTRAN(V9x,V9y,V9z, V10Be, Vx,Vy,Vz)
c       Done.  Get energy of 9-Be.


      V9Be=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
c       Check for further neutron or proton decay of 9-Be


      IF (Npgo .EQ. 0) goto 19
c       Program to here, further particle decay.  First rotate 9-Be


c         velocity vector to "detector" coordinates.


      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Z1,ZB,Z3)
c       Done.  Set up and test for proton decay of 9-Be.


      Egamma=0.0
      Eleftov=0.0
      IF (Npgo .EQ. 1) goto 16
c       Program to here, it's  9-Be --> n + 8-Be.


      Taex=Excit-Be8gs-Be8pn
      Enc=RCKE(Kn,K8Be,Taex)
      Vnc=VELOCITY(Kn,Enc)
      CALL RVECT(Xp,Yp,Zp)
      Vnx=Xp*Vnc
      Vny=Yp*Vnc
      Vnz=Zp*Vnc
      CALL LABTRAN(Vnx,Vny,Vnz, V9Be, Vx,Vy,Vz)
      Vnlab=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Eneut2=EFROMV(Kn,Vnlab)
c       Rotate neutron velocity into "detector" lab. coordinates.


      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Z1=Xn
      ZB=Yn
      Z3=Zn
      CALL TRANSVEC(Z1,ZB,Z3)
c       Done.  Put values into  /NEUTR2/  common variables.


      U2=Xn
      V2=Yn
      W2=Zn
      X2=X
      Y2=Y
      Z2=Z
C       Now get 8-Be ion motion.


      E8Bec=Taex-Enc
      V8Bec=VELOCITY(K8Be,E8Bec)
      V8x=-Xp*V8Bec
      V8y=-Yp*V8Bec
      V8z=-Zp*V8Bec
      CALL LABTRAN(V8x,V8y,V8z, V9Be, Vx,Vy,Vz)
      V8Be=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
c       Done.  Now  8-Be --> 2 alphas.


      Ta=Be8gs
      Ea=0.5*Ta
      Vac=VELOCITY(Ka,Ea)
      CALL RVECT(Xp,Yp,Zp)
      Vax=Xp*Vac
      Vay=Yp*Vac
      Vaz=Zp*Vac
      CALL LABTRAN(Vax,Vay,Vaz, V8Be, Vx,Vy,Vz)
      Valph=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Echrg(2)=EFROMV(Ka,Valph)
c       That's one alpha.  Now get the other.


      Vax2=-Vax
c          (By symmetry.)


      Vay2=-Vay
      Vaz2=-Vaz
      CALL LABTRAN(Vax2,Vay2,Vaz2, V8Be, Vx,Vy,Vz)
      Valph2=SQRT(Vx*Vx + Vy*Vy +Vz*Vz)
      Echrg(3)=EFROMV(Ka,Valph2)
c           Done!!


      Goto 25
c
c       Next is proton decay, i.e.,  9-Be --> p + 8-Li


c         A complexity at this point.  If there is a second neutron


c          waiting to be processed, then Subroutine PP8LI must be


c          alerted to inhibit further neutron decay of 8-Li.  In


c          addition, Subroutine ELI8DK (which is called by PP8LI)


c          puts the information on the neutron from  8-Li --> n + 7-Li


c          into the /NEUTRN/ common; however, this common area already


c          has information for the neutron from  10-Be --> n + 9-Be


c          reaction just completed.  So need to do some preparation.


   16 V9=V9Be
      IF (Eneut2 .GT. 0.0) goto 17
c       If program counter to here, no second neutron waiting to be


c         processed, so shift first neutron information into /NEUTR2/ common.


      Eneut2 = Eneut
      U2 = U
      V2 = V
      W2 = W
      X2 = X
      Y2 = Y
      Z2 = Z
      Eneut = 0.0
      I2=0
      Goto 18
c
c    But if there is a second neutron, alert -PP8LI- through variable V9Be.
   17 V9=-V9Be
      I2=1
   18 CALL PP8LI(Excit,V9,Egamma)
c       If both neutrons determined by present routine, set Eneut2 to
c        negative to alert -BANKR2- ; will return to >0 before exit from
c        routine.


      IF (I2 .EQ. 0) Eneut2=-Eneut2
      Goto 25
c
   19 Echrg(3)=EFROMV(K9Be,V9Be)
c       That's it for 9-Be ground-state results!


      Goto 25
c          END  n + 9-Be SECTION.


c
c       Tidy up 10-Be ion.


   20 Echrg(3)=EFROMV(K10Be,V10Be)
c
   25 Nelm=12
      CALL BANKR2
      IF (Eneut2 .LT. 0.0) Eneut2=-Eneut2
c       Check for gamma ray interaction.


      IF (Egamma .GT. 0.0) CALL PHOTON(Egamma)
c       (See P.S.D. comment following statement 26 in subroutine NT.)


      IF (Egamma .LE. 0.0) Return
      IF (Eleftov.GT. 0.0) CALL PHOTON(Eleftov)
      Return
c
c
c          At exit from this routine particle energies are stored:


c-----------------------------------------------------------------------------
c   Decay Mode                  Echrg array                     Eneut   Eneut2


c                (1)     (2)     (3)     (4)     (5)     (6)


c-----------------------------------------------------------------------------
c   3-He+10-Be   3-He    0.0    10-Be    0.0     0.0     0.0     0.0      --
c  3-He+n+9-Be   3-He    0.0     9-Be    0.0     0.0     0.0      En      --
c   3-He+2n+2a   3-He   alpha   alpha    0.0     0.0     0.0     En1     En2
c 3-He+n+p+8-Li  3-He  p(esc)    8-Li  proton    0.0     0.0      En      --
c 3-He+2n+p+7-Li 3-He  p(esc)    7-Li  proton    0.0     0.0     En1     En2
c 3-He+2n+p+t+a  3-He  p(esc)    0.0   proton   alpha   triton   En1     En2
c-----------------------------------------------------------------------------
c  N.B.--cannot get to 6-Li in this routine as that would entail a
c       third neutron


c
      END
c------------------------------------------------------------------
C   This is file NPX.FOR
c
      Subroutine NPX
c       Purpose to compute energetics of the reactions starting with
c       n + 12-C --> p + 12-B.  Program first checks for populating the
c       residual 12-B ion in a low-lying (particle-bound) energy level,
c       and if so program reverts to Subroutine NP.  If not, then
c       test for possible second proton emission, or possible alpha
c       emission, and if neither of those then consider neutron emission.
c       The several possible "next" reactions are:
c          12-B --> p + 11-Be,
c          12-B --> alpha + 8-Li, or
c          12-B --> n + 11-B.
c       Following each of these reactions the program will consider
c       possible further particle emissions (if the excitation energy
c       is large enough) of the residual nuclei 11-Be, 8-Li or 11-B.
C
      Real*4 Li6pa, Li7pa, Li7pn, Li8pa
      Common /RANDM/  Irx
      Common /COLLIS/ Nelm,Echrg(6)
      Common /NAID/   Rad,Ht,Rc,Rz
      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
      Common /NEUTRN/ SPDSQ, U, V, W, Xneut, Yneut, Zneut
      Common /MASSES/ Emass(21)
      COMMON/PROTON/NOPROT,RNGEN(144),PRAN(144),DRAN(144),ARAN(144),
     a              HRAN(144),TRAN(144)
c      Common /PROTON/ Noprot,Range(95),Proten(95)
      Common /PPOLAR/ X(37),Ptheta(37),N37
c       The arrays X and Ptheta are given in Subroutine NP
      Common /EXC11B/ Exc11(14),TenBpn,Li7pa
      Common /NEUTR2/ Eneut2, U2,V2,W2, X2,Y2,Z2
      Common /EXC10B/ Ex10B(9)
c       Array Ex10B is in subroutine NT.
      Common /EXC8LI/ E8Li(10)
c       Array E8Li is in Function EX8LI.
      Common /EXCTBE/ Ex10Be(4),Q9Bepn
ckaji 2011/03/18
      common /RES/ iswt, bias(5), icont, ictrl, iesc, ilight
c
      Dimension Ex12B(18), PP2P(18)
      Dimension Ex12B1(12),Ppalph(12)
c
c       Next array has a representation of excited states of 11-B.
      Data Exc11/0.0,2.125,4.445,5.02,6.743,6.792,7.286,7.978,
     a 8.56,8.92,9.2,9.88,10.33,11.25/
c          The last 5 decay by both alpha and gamma.
c
      Data Qnp/12.613/, Be11pp/14.096/, Be10pn/0.503/, Ka/3/
      Data Q9Bepn/6.812/, Nterp/1/, Npp2p/18/, B11pn/3.37/
      Data Kn/1/, K12C/7/, Kp/2/, K12B/6/, K11B/5/, K11Be/14/, K10Be/15/
      Data K10B/12/, TenBpn/11.46/, Li7pa/8.665/, Li6pa/4.46/, K9Be/4/
      Data Li7pn/2.033/, Li8pa/10.002/, K8Li/16/, Nppalf/12/
c
c       Next data array represent particle-bound states of 10-Be.
      Data Ex10Be/0.0, 3.368, 5.959, 6.263/
c
c       Next arrays to determine probability of second proton emission.
c         Calculations using TNG show that  (n,2p)/(n,p)  is a function
c         of the 12-B excitation, but essentially independent of incident
c         neutron energy.
      Data Ex12B/14.1, 14.45, 15.06,  16.0,  16.8,  18.0,  19.2, 20.0,
     a  22.0, 24.0, 25.6, 27.6, 29.4, 32.4, 34.6, 36.4, 39.0, 45.0/
      Data PP2P/0.0,0.0012, 0.0039, 0.011, 0.019, 0.032, 0.042,0.05,
     a 0.068, 0.085, 0.1, 0.13, 0.17, 0.24, 0.30, 0.35, 0.40, 0.42/
c
c       Next arrays to determine probability of alpha emission following
c         first proton emission.  Again, TNG calculations show dependence
c         on 12-B excitation, but essentially not on incident neutron energy.
      Data Ex12B1/10.003, 10.14, 11.4, 12.5, 13.6, 19.6, 22.0, 26.0,
     a   30.0, 38.0, 50.0, 100.0/
      Data Ppalph/0.0, 0.00036, 0.047, 0.1, 0.15, 0.15, 0.1, 0.075,
     a   0.048, 0.03, 0.024, 0.018/
c
C   The reaction scheme is:
c
c       (a)  n + 12-C --> p + 12-B;      then:
c               Decide if residual 12-B is in a bound state, and if so
c                 go to Subroutine NP.  If not, then:
c       (b)      12-B --> n + 11-B;
c               Check for possible further neutron emission; if so, then:
c       (c)      11-B --> n + 10-B;
c
c       New additions (5/87):
c       (b')     12-B --> p + 11-Be, with possible further emission:
c       (c')     11-Be --> n + 10-Be, with possible further emission:
c       (d')     10-Be --> n + 9-Be (ground state, only)
c
c       Further additions (5/87):
c       (b")     12-B --> alpha + 8-Li, with possible further emission:
c       (c")     8-Li --> gamma or neutron decay.
c
c       Added 9/87:
c       (c3)    11-B --> alpha + 7-Li, with possible further 7-Li breakup.
c
c       Initial set-up:
      En2=0.0
      Do 1 I=1,6
   1  Echrg(I)=0.0
c               Go to Center-of-mass coordinates
      En=SPDSQ
      Vn = VELOCITY(Kn,En)
      CALL CMTRAN(Kn,K12C,Vn, Vcom,Eneut,Ecc)
      Tec=Ecc+Eneut
      Ta=Tec-Qnp
      IF (Ta .GT. 0.0) goto 2
csatoh
c     Type 100, En
      write(*,222)En
  222 Format(/'   ***  Error in Subroutine NPX; E(neut) at entry ='
     b 1PE11.3/10x,'Set E(Neut) = E(Prot) = 0 and exit.'/)
      SPDSQ=0.0
      Return
c
c
c    First step -- decide if particle-stable state of 12-B is involved.
    2 Ratio=SIGCNP(En)/(SIGCNP(En) + SIGCNPN(En))
      IF (RAN(Irx) .GT. Ratio) goto 7
    3 CALL NP
      Return
c          (The point to this last bit is that originally the (n,p) and
c            (n,pn) reactions were considered entirely separately.  The
c            present amalgamation is, thus, a hybrid...)
c
c       Set up for first proton emission:
C              CHOOSE C.M. PROTON ENERGY using Function CHOOSP
    7 Tapn=Ta-B11pn
      IF (Tapn .LE. 0.0) goto 3
c       (Last statement just a check; probably shouldn't be needed.)
C
c   At this point we consider the proton spectrum for En = 60 MeV neutron
c       interactions with 12-C as reported by Brady et al, Journal of Physics
c       G 10 (1984) 363, which shows sufficient enhancement of the population
c       of (unbound) levels at about 4.5 and 8.3 MeV in 12-B to be worth
c       programming.  Thus ...
      IF (En .LE. 34.0) goto 8
      Enhanc=Ratio*(En-34.0)/50.0
c      Enhanc=0.053 ! kajimoto 06/10
      if(En .GE. 40.0) Enhanc=0.025 ! kajimoto 06/20
      if(En .GE. 50.0) Enhanc=0.015
      if(En .GE. 80.0) Enhanc=0.015
      Excit=4.5
      Try=Ta-Excit
      Rani=RAN(Irx)
      IF (Rani .LE. Enhanc) goto 12
      IF (En .LE. 39.0) goto 8
      Excit=8.3
      Try=Ta-Excit
      Enhanc2=Enhanc+Ratio*(En-39.0)/30.0
      IF (Rani .LE. Enhanc2) goto 12
c
c       Program counter to here -- treat proton output as "continuum"
    8 Epcom=RCKE(Kp,K12B,Tapn)
c        Epcom is the maximum energy the proton can have in the center of mass
      F=0.1245+0.001*ABS(En-45.0)
      Temp=F*En
      if(En .GE. 37.0) Temp=0.29672*En-6.1121 ! kajimoto 06/10
c       Temp value chosen empirically so that CHOOSP gives a reasonable
c        representation of the proton "continuum" as calculated by the
c        nuclear model code TNG of C. Y. Fu.
      Eprot = CHOOSP(Epcom,Temp)
c        The difference in energy between Epcom and Eprot could be
c         interpreted as the energy of a "virtual" energy level in the 12-B
c         nucleus.  The effect of this "level" should be to increase
c         the value of Q by the "virtual" excitation energy.
      Try=Eprot*(Emass(Kp)+Emass(K12B))/Emass(K12B)
c               That's approximately the "total available energy" that
c                would yield Eprot from the function RCKE.
      Eptry=RCKE(Kp,K12B,Try)
      Try=Try*Eprot/Eptry
c               (Iterate once on TRY for slightly more accuracy.)
      Excit=Ta-Try
c               EXCIT  should be close to the "virtual" excitation energy
c                 of the residual 12-B ion.
      IF (Excit .LE. B11pn) goto 3
c       Last is a check, mainly on possible single-precision round-off error.
c
   12 Qex=Qnp + Excit
      Eprot=RCKE(Kp,K12B,Try)
c
C         FIND THE COSINE  OF THE POLAR DIRECTION OF THE PROTON(C.M.)
   15 Rand=RAN(IRX)
   18 Ctheta = EXTERP(Ptheta,X,Rand,N37,Nterp)
c added by kajimoto 06/13
      if(En .GE. 40.0) then
        Ctheta=akalbach(En,Eprot,13.0,6.0,0.0,1.0,1.0,0.0)
        goto 20
      endif
c added by kajimoto 06/13
      IF (Tec .GE. 51.7) goto 20
      Ciso=1.0-2.0*Rand
      Slope=(Ctheta-Ciso)/(51.7-Qex)
      Ctheta=Slope*(Tec-Qex)
C          FIND THE SINE AND COSINE OF THE AZIMUTHAL DIRECTION OF THE
C           PROTON (C.M.) -- By Random-number Choice.
   20 Rann = 6.2831853*RAN(Irx)
      Sinphi = SIN(Rann)
      Cosphi = COS(Rann)
c         FIND THE C. M. DIRECTION COSINES OF THE PROTON.
      Sintheta=SQRT(1.-Ctheta*Ctheta)
      Xp=Sintheta*Cosphi
      Yp=Sintheta*Sinphi
      Zp=Ctheta
      Vpr=VELOCITY(Kp,Eprot)
C        NOW GET C. M. VELOCITY COMPONENTS OF THE PROTON.
      Vpx=Xp*Vpr
      Vpy=Yp*Vpr
      Vpz=Zp*Vpr
      CALL LABTRAN(Vpx,Vpy,Vpz, Vcom, Vx,Vy,Vz)
C        THAT HAS GOT THE VELOCITY COMPONENTS OF THE PROTON IN THE
C               "NEUTRON" LABORATORY COORDINATE SYSTEM.
      Vprot=SQRT(Vx*Vx +Vy*Vy +Vz*Vz)
C         FIND THE LAB ENERGY OF THE PROTON.
      Epl=EFROMV(Kp,Vprot)
c      Store information; then get distance from interaction point to
c       detector surface in the direction of the proton's velocity.
      Echrg(1)=Epl
      Echrg(2)=0.0
c     Get proton dir. cosines in "detector" laboratory coord. system.
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(U,V,W)
c       Done.  Now test for proton escape from detector.
ckaji      IF (Epl .LE. 0.1) goto 70
      IF (Epl .LE. 0.07) goto 70
      Protpl=PLNGTH(Xn,Yn,Zn, Xneut,Yneut,Zneut, Rad,Ht)
C       FIND RANGE OF THE PROTON.
      RPROT = EXTERP(RNGEN,PRAN,Epl,Noprot,4)
C  DETERMINE THE AMOUNT OF THIS RANGE THAT LIES WITHIN THE SCINTILLATOR.
      RangeQ=Rprot-Protpl
      IF (RangeQ .LE. 1.02e-04) goto 70
c           If it gets to here, proton "escapes" ...
      Eprotf=EXTERP(PRAN,RNGEN,RangeQ,Noprot,4) ! kajimoto 2011/03/18
      Echrg(2)=Eprotf
C
c
C       FIND THE COMPONENTS OF THE 12B VELOCITY IN THE LAB.
c       First adjust total available energy for the energy of the
c             "virtual" excited state.
   70 Taex=Tec-Qex
      E12B=RCKE(K12B,Kp,Taex)
      V12B=VELOCITY(K12B,E12B)
      Vbcx=-Xp*V12B
      Vbcy=-Yp*V12B
      Vbcz=-Zp*V12B
c        (Recall that the direction cosines of the 12-B ion in the center
c         of mass are negatives of the direction cosines of the proton.)
      CALL LABTRAN(Vbcx,Vbcy,Vbcz, Vcom, Vbx,Vby,Vbz)
      V12B=SQRT(Vbx*Vbx + Vby*Vby + Vbz*Vbz)
c
      CALL DIRCOS(Vbx,Vby,Vbz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(U,V,W)
c     Now variables Xn,Yn,Zn = direction cosines of the motion of the
c        12-B in "detector" laboratory coordinate system.
c
c     New addition 5/87: check for second proton or alpha emission.
      IF (Excit .LE. Ex12B(1)) goto 72
c
c           SECOND PROTON OR ALPHA EMISSION
c       If program gets to here, enough excitation in 12-B for second
c        proton or possibly alpha emission.
      Fr=EXTERP(Ex12B,PP2P,Excit,Npp2p,Nterp)
      Rani=RAN(Irx)
      IF (Rani .LE. Fr) goto 71
      Fr=Fr+EXTERP(Ex12B1,Ppalph,Excit,Nppalf,Nterp)
      IF (Rani .GT. Fr) goto 72
c
C
C
C
C   ***********************************************************************
C   *                                                                     *
C   *                 ALPHA PLUS 8-LI                                     *
C   *                                                                     *
c   * Program counter to here means alpha emission:  12-B --> alpha + 8-Li
C   *                                                             *
      Ta=Excit-Li8pa
      Excita=EX8LI(Excit)
      Egamma=0.0
      Eleftov=0.0
      IF (Excita .GT. Li7pn) goto 50
      Egamma=Excita
c       Either ground state or 1st excited state of 8-Li, and no neutron.
      SPDSQ=0.0
   50 Taex=Ta-Excita
      Epa=RCKE(Ka,K8Li,Taex)
      CALL RVECT(Xpa,Ypa,Zpa)
      Vpa=VELOCITY(Ka,Epa)
      Vacx=Xpa*Vpa
      Vacy=Ypa*Vpa
      Vacz=Zpa*Vpa
      CALL LABTRAN(Vacx,Vacy,Vacz, V12B, Vx,Vy,Vz)
      Valf=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Ealf=EFROMV(Ka,Valf)
      Echrg(4)=Ealf
c     That's it for the alpha; now get information on 8-Li ion.
      E8Lic=Taex-Epa
      V8Lic=VELOCITY(K8Li,E8Lic)
      Vax=-Xpa*V8Lic
      Vay=-Ypa*V8Lic
      Vaz=-Zpa*V8Lic + V12B
      V8Li=SQRT(Vax*Vax + Vay*Vay + Vaz*Vaz)
c     Check for further neutron emission by 8-Li ion.
      IF (SPDSQ .GT. 0.0) goto 52
c     Program counter to here = no neutron emission, so tidy up and exit.
      Echrg(3)=EFROMV(K8Li,V8Li)
      Goto 90
c
c      Next, set up for decay of 8-Li.
   52 CALL DIRCOS(Vax,Vay,Vaz, Xpn,Ypn,Zpn)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
c     Ready to go to 8-Li decay subroutine.
      CALL ELI8DK(Excita,V8Li,Egamma)
c     Finished breaking up 8-Li, so finish off.
      Goto 90
c
c
C
C   ************************************************************************
C   *                                                                      *
C   *                 PROTON PLUS 11-BE                                    *
C   *                                                                      *
c   *    Got to here: 2d PROTON emission to consider.                      *
C   *                                                                      *
   71 Ta=Excit-Be11pp
      Ep2com=RCKE(Kp,K11Be,Ta)
      Ehat=En-Qnp
      F=0.1245+0.001*ABS(Ehat-45.0)
      Temp=F*Ehat
      Eprot2=CHOOSP(Ep2com,Temp)
      Try2=Eprot2*(Emass(Kp) + Emass(K11Be))/Emass(K11Be)
      Ep2try=RCKE(Kp,K11Be,Try2)
c       and as done before, iterate once on Try2 for better precision
      Try2=Try2*Eprot2/Ep2try
      Excitp=Ta-Try2
      Egamma=0.0
      Eleftov=0.0
      IF (Excitp .GT. 1.6) goto 30
c     Low-energy excitation of 11-Be.  Check for E(gamma).
      IF (Excitp .LT. 0.3918) goto 29
      Excitp = 0.3198
      Egamma = 0.3198
c     Only one excited bound state of 11-Be.  However, no subsequent
c       neutron emission -- just (n,2p).
   29 SPDSQ=0.0
   30 Taex=Ta-Excitp
c     Get information on outgoing 2d proton.
      Epp=RCKE(Kp,K11Be,Taex)
c     Isotropy for 2d proton emission in center-of-mass coordinates (of
c        moving 12-B ion)
      CALL RVECT(Xp,Yp,Zp)
      Vpp=VELOCITY(Kp,Epp)
      Vpcx=Xp*Vpp
      Vpcy=Yp*Vpp
      Vpcz=Zp*Vpp
      CALL LABTRAN(Vpcx,Vpcy,Vpcz, V12B, Vpx,Vpy,Vpz)
      Vprot=SQRT(Vpx*Vpx + Vpy*Vpy + Vpz*Vpz)
      Eprot=EFROMV(Kp,Vprot)
      Echrg(5)=Eprot
      Echrg(6)=0.0
c     Now check for 2d proton escape.
      CALL DIRCOS(Vpx,Vpy,Vpz, Xpn,Ypn,Zpn)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
c     Now have 2d proton direction cosines in "detector" lab. coord. system.
c       Its direction cosines are Xn, Yn, Zn.
      Protpl=PLNGTH(Xn,Yn,Zn, Xneut,Yneut,Zneut, Rad,Ht)
      Rprot=EXTERP(RNGEN,PRAN, Eprot, Noprot,Nterp)
      RangeQ=Rprot-Protpl
      IF (RangeQ .LE. 0.0) goto 32
c     If program gets to here, 2d proton escapes detector.
      Echrg(6)=EXTERP(PRAN,RNGEN,RangeQ,Noprot,Nterp) ! kajimoto 2011/03/18
c     Now get 11-Be velocity and energy in lab. coord. system
   32 E11Be=Taex-Epp
      V11Be=VELOCITY(K11Be,E11Be)
      Vbx=-Xp*V11Be
      Vby=-Yp*V11Be
      Vbz=-Zp*V11Be + V12B
c     Again, non-relativistic transformation from c.o.m. to lab. coordinates
c       of the moving 12-B ion.
      V11Be=SQRT(Vbx*Vbx + Vby*Vby + Vbz*Vbz)
      IF (SPDSQ .GT. 0.0) goto 35
c     Program counter to here means no further neutron emission.  We
c      excited either the ground state or the 320-keV state in 11-Be.
   33 E11Be=EFROMV(K11Be,V11Be)
      Echrg(3)=E11Be
      Goto 90
c
C   *************************************************************************
C   *                                                                       *
C   *          NEUTRON PLUS 10-BE                                           *
C   *                                                                       *
c   * If program goes to next step, then have (n,2pn) reaction to finish up.
c   *  Schematic will look something like the drawing in the N2N subroutine.
C   *                                                                       *
   35 CALL DIRCOS(Vbx,Vby,Vbz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
c     Now have 11-Be ion dir. cosines (i.e. Xn,Yn,Zn) in "detector" lab.
c      coordinate system.
      Ta=Excitp-Be10pn
      Encom=RCKE(Kn,K10Be,Ta)
      Ehat=Ehat-Be11pp
      Fn=0.065+0.001*Ehat
      Temp=Fn*Ehat
      Eloww=0.0
c     Test if Eneut2 has a value left over from a previous "2n" reaction --
c      if so then force present neutron emission to populate a state in
c      10-Be having Ex < binding energy of 9-Be + n.
      IF (Eneut2 .LE. 0.0) goto 36
      Tryy=Encom-Q9Bepn
      IF (Tryy .GT. 0.0) Eloww=Tryy
   36 Enn=CHOOSN(Eloww,Encom,Temp)
      IF (Enn.GT.0.0) goto 37
      SPDSQ=0.0
      Goto 33
   37 Try=Enn*(Emass(Kn) + Emass(K10Be))/Emass(K10Be)
c     Once again, an iteration for precision.
      Entry=RCKE(Kn,K10Be,Try)
      Try=Try*Enn/Entry
      Excitn=Ta-Try
      Egamma=0.0
      Eleftov=0.0
      IF (Excitn .GT. Q9Bepn) goto 44
c
      K=1
c     Now pair -Excitn- with a level energy for 10-Be.
      Do 40 I=2,4
      IF (Excitn .LT. Ex10Be(I)) goto 42
   40 K=K+1
   42 Excitn=Ex10Be(K)
c       Set up gamma decay characteristics:
      IF (K .GE. 2) Egamma=Ex10Be(2)
      IF (K.EQ.3 .AND. RAN(Irx).GE.0.5) Egamma=Excitn
      Eleftov=Excitn-Egamma
c     Now get lab. coord. energy of emitted neutron.
c       Again, isotropy of emission in center of mass (of the 11-Be ion).
   44 Taex=Ta-Excitn
      Enn=RCKE(Kn,K10Be,Taex)
      CALL RVECT(Xp,Yp,Zp)
      Vnn=VELOCITY(Kn,Enn)
      Vncx=Xp*Vnn
      Vncy=Yp*Vnn
      Vncz=Zp*Vnn
      CALL LABTRAN(Vncx,Vncy,Vncz, V11Be, Vx,Vy,Vz)
      Vneut=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Eneut=EFROMV(Kn,Vneut)
      SPDSQ=Eneut
c     Done.  Now get neutron dir. cosines into "Detector" lab. coordinates.
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
c     Done.  Xn,Yn,Zn are dir. cosines of neutron in "detector" lab. coords.
      U=Xn
      V=Yn
      W=Zn
c     Now get 10-Be energy in laboratory coordinates.
      E10Be=Taex-Enn
      V10Be=VELOCITY(K10Be,E10Be)
      Vbx=-Xp*V10Be
      Vby=-Yp*V10Be
      Vbz=-Zp*V10Be + V11Be
c     Non-relativistic transformation from c.o.m. to lab. coords.
      V10Bel=SQRT(Vbx*Vbx + Vby*Vby + Vbz*Vbz)
c
C
C   **************************************************************************
C   *                                                                        *
C   *                NEUTRON PLUS 9-BE                                       *
C   *                                                                        *
c   *   Test for possible SECOND NEUTRON emission.                           *
C   *                                                                        *
      IF (Excitn .LE. Q9Bepn) goto 60
c
c           SECOND NEUTRON:  10-Be --> n + 9-Be
      CALL DIRCOS(Vbx,Vby,Vbz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
c     Now Xn,Yn,Zn are dir. cosines of 10-Be ion in "detector" lab. coords.
      Ta=Excitn-Q9Bepn
      En2c=RCKE(Kn,K9Be,Ta)
      CALL RVECT(Xn2,Yn2,Zn2)
      Vn2=VELOCITY(Kn,En2c)
      Vn2x=Xn2*Vn2
      Vn2y=Yn2*Vn2
      Vn2z=Zn2*Vn2
      CALL LABTRAN(Vn2x,Vn2y,Vn2z, V10Bel, Vx,Vy,Vz)
      Vneut=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      En2=EFROMV(Kn,Vneut)
      Eneut2=-En2
c     This negative is temporary, just to alert BANKER routine.
c     Now get direction cosines of 2d neutron into "detector" coords.
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
c     Done.  Transfer to 2d neutron variables.
      U2=Xn
      V2=Yn
      W2=Zn
      X2=Xneut
      Y2=Yneut
      Z2=Zneut
c     and that's done.  Get energy of 9-Be.
      E9be=Ta-En2c
      V9Be=VELOCITY(K9Be,E9Be)
      Vbx=-Xn2*V9Be
      Vby=-Yn2*V9Be
      Vbz=-Zn2*V9Be
      CALL LABTRAN(Vbx,Vby,Vbz, V10Bel, Vx,Vy,Vz)
      V9Belab=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Echrg(3)=EFROMV(K9Be,V9Belab)
      Goto 90
c        END SECTION ON 10-Be --> n + 9-Be
c
c     Finish case where 10-Be is final heavy ion.
   60 Echrg(3)=EFROMV(K10Be,V10Bel)
      Goto 90
c        END OF THIS SECTION
c
c
c
C   ************************************************************************
C   *                                                                      *
C   *                NEUTRON PLUS 11-B                                     *
C   *                                                                      *
c   *         Next section for 12-B --> n + 11-B                     *
c   *                                                                      *
   72 Ta=Excit-B11pn
      Encom=RCKE(Kn,K11B,Ta)
      Ehat=En-Qnp
      F=0.065 + 0.001*Ehat
      Temp=F*Ehat
c       Temp value sort of empirically determined so that CHOOSN gives a
c        reasonable representation of the neutron "continuum" as calculated
c        by the nuclear model code TNG.
      Eloww=0.0
c     Test if Eneut2 has a value left over from a previous "2n" reaction --
c      if so then force the present neutron emission to populate a state
c      in 11-B having Ex < binding energy of 10-B + n.
      IF (Eneut2 .LE. 0.0) goto 73
      Tryy = Encom-TenBpn
      IF (Tryy .GT. 0.0) Eloww=Tryy
   73 Enn=CHOOSN(Eloww,Encom,Temp)
c
c   AND As observed when choosing E(proton) from a distribution (as done above)
c     the difference in energy between Encom and Enn could be interpreted as
c     a "Virtual" energy level in 11-B.
c
      Try=Enn*(Emass(Kn)+Emass(K11B))/Emass(K11B)
      Entry=RCKE(Kn,K11B,Try)
      Try=Try*Enn/Entry
      Excitn=Ta-Try
c       If EXCITN is relatively small, pair it up with a known excited
c           state in 11-B.
      Egamma=0.0
      Eleftov=0.0
      IF (Excitn .LT. TenBpn) goto 75
      IF (Excitn .GT. TenBpn+0.5) goto 80
      Excitn=TenBpn+0.25
      K=15
      Goto 79
c     (That covers alpha emission in competition with neutron emission...)
   75 K=1
      Do 74 I=2,14
      IF (Excitn .LT. Exc11(I)) goto 76
   74 K=K+1
   76 Excitn=Exc11(K)
      Egamma=Exc11(K)
c     Most 11-B levels decay predominantly by ground-state transitions;
c      however, for three levels the dominant decay transition is through
c      an excited state.  In such cases, there are at least two gamma rays
c      to be followed, including the ground-state transition of
c      the intermediate excited state.
      IF (K.EQ.8 .OR. K.EQ.9) Egamma=Exc11(2)
      IF (K .EQ. 11) Egamma=Exc11(3)
      Eleftov=Exc11(K)-Egamma
c     That's the "other" gamma radiation.
c
   79 Taex=Ta-Excitn
      Enn=RCKE(Kn,K11B,Taex)
c
   80 CALL RVECT(Xp,Yp,Zp)
      Vnn=VELOCITY(Kn,Enn)
      Vncx=Xp*Vnn
      Vncy=Yp*Vnn
      Vncz=Zp*Vnn
      CALL LABTRAN(Vncx,Vncy,Vncz, V12B, Vnx,Vny,Vnz)
      Vneut=SQRT(Vnx*Vnx + Vny*Vny + Vnz*Vnz)
      Eneut=EFROMV(Kn,Vneut)
      SPDSQ=Eneut
c        Now get "new" neutron direction cosines in the
c         "detector" laboratory coordinates.
      CALL DIRCOS(Vnx,Vny,Vnz, Xpn,Ypn,Zpn)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
      U=Xn
      V=Yn
      W=Zn
c     and, finally, get the 11-B energy.
      Taex=Ta-Excitn
      E11B=RCKE(K11B,Kn,Taex)
      V11B=VELOCITY(K11B,E11B)
      Vbx=-Xp*V11B
      Vby=-Yp*V11B
      Vbz=-Zp*V11B
      CALL LABTRAN(Vbx,Vby,Vbz, V12B, Vx,Vy,Vz)
      V11B=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
c
c     At this point want to look for second neutron emission,
c        i.e. (n,p2n) reaction.
c     First check -EXCITN- of 11-B for possible alpha emission.
      IF (Excitn .LE. Li7pa) goto 88
      IF (K .EQ. 10) goto 88
      Tryn=RAN(Irx)
      IF (K.EQ.11 .AND. Tryn.LT.0.9) goto 88
      IF (K.EQ.13 .AND. Tryn.LT.0.0001) goto 88
c     If program counter gets to here, could be alpha or second
c       neutron emission.
c     First get moving 11-B ion into "detector" coordinates.
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
c     Done.  Xn,Yn,Zn are dir. cosines of 11-B motion in "detector" coords.
      IF (K.EQ.15 .AND. Tryn.LT.0.56) goto 81
      IF (Excitn .GT. TenBpn) goto 82
C
C
C
C   ************************************************************************
C   *                                                                      *
C   *             ALPHA PLUS 7-LI                        *
C   *                                                    *
c   *    Program counter to here, it's  11-B --> 7-Li + alpha              *
C   *                                                    *
   81 Egamma=0.0
      Eleftov=0.0
      Ta=Excitn-Li7pa
      CALL AP7LI(Ta,V11B,Egamma)
      Goto 90
c     END SECTION ON ALPHA PRODUCTION
C
C
C
C   ************************************************************************
C   *                                                                      *
C   *                NEUTRON PLUS 10-B                                     *
C   *                                                                      *
C   *   Next is case of second neutron emission, 11-B --> n + 10-B.        *
C   *                                                                      *
   82 Ta=Excitn-TenBpn
      En2com=RCKE(Kn,K10B,Ta)
      Ehat=Ehat-B11pn
      Fn=0.065 + 0.001*Ehat
      Temp=Fn*Ehat
      Eloww=0.0
      En2=CHOOSN(Eloww,En2com,Temp)
      Try=En2*(Emass(Kn) + Emass(K10B))/Emass(K10B)
c     Iterate once for greater precision (since it costs very little)
      Entry=RCKE(Kn,K10B,Try)
      Try=Try*En2/Entry
      Excitn2=Ta-Try
      Egamma=0.0
      Eleftov=0.0
      Level=10
      IF (Excitn2 .GT. 6.11) goto 86
      K=1
c     Pair -EXCITN2- with 10-B energy level.
      Do 84 I=2,9
      IF (Excitn2 .LT. Ex10B(I)) goto 85
   84 K=K+1
   85 Excitn2=Ex10B(K)
      Level=K
c     See discussion following statement No. 18 of Subroutine NT for
c       gamma decay of levels in 10-B.
      IF (K .GT. 1) Egamma=Ex10B(2)
      Eleftov=Ex10B(K)-Egamma
c     Now get information on 2d neutron into "detector" lab. coords.
   86 Taex=Ta-Excitn2
      Enn=RCKE(Kn,K10B,Taex)
      CALL RVECT(Xn2,Yn2,Zn2)
      Vn2=VELOCITY(Kn,Enn)
      V2x=Xn2*Vn2
      V2y=Yn2*Vn2
      V2z=Zn2*Vn2
      CALL LABTRAN(V2x,V2y,V2z, V11B, Vx,Vy,Vz)
      Vneut2=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      En2=EFROMV(Kn,Vneut2)
      Eneut2=-En2
c     The negative is temporary, just to alert BANKER routine.
c     Now get direction cosines of 2d neutron
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
c     Done.  Xn,Yn,Zn are the needed direction cosines
c     Fill in rest of data for 2d neutron.
      U2=Xn
      V2=Yn
      W2=Zn
      X2=Xneut
c     Assumption: various heavy ions (12-C, 12-B, ... ) travel a very
c      short path compared to detector dimensions.
      Y2=Yneut
      Z2=Zneut
c     Done with 2d neutron; now get 10-B ion information.
      E10B=Taex-Enn
      V10B=VELOCITY(K10B,E10B)
      Vbx=-Xn2*V10B
      Vby=-Yn2*V10B
      Vbz=-Zn2*V10B + V11B
c     Again, non-relativistic transformation of c.o.m. --> lab. coords.
      V10Blab=SQRT(Vbx*Vbx + Vby*Vby + Vbz*Vbz)
c     Test for 10-B particle-emission decay:
      IF (Level .LE. 6) goto 87
      IF (Level.EQ.8 .AND. RAN(Irx).LT. 0.5) goto 87
c     Program counter to here, prepare for 10-B particle decay.
      Eleftov=0.0
      CALL DIRCOS(Vbx,Vby,Vbz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
c     All prepared.  Go get 'em, tiger.
      CALL TENBDK(Excitn2,V10Blab,Egamma)
      Goto 90
c     Tidy up (n,p2n) results:
   87 Echrg(3)=EFROMV(K10B,V10Blab)
      Goto 90
c        END THIS SECTION OF SECOND NEUTRON EMISSION
c
c     Tidy up (n,pn) results:
   88 E11B=EFROMV(K11B,V11B)
      Echrg(3)=E11B
c
   90 Nelm=7
      CALL BANKR2
c     Now correct the Eneut2 variable.
      Eneut2=ABS(Eneut2)
c     Now check for possible gamma-ray contributions
      IF (Egamma .GT. 0.0) CALL PHOTON(Egamma)
      IF (Egamma .LE. 0.0) Return
      IF (Eleftov .GT. 0.0) CALL PHOTON(Eleftov)
      Return
c
c
c        At end of this routine, Echrg(1) has E(proton) and
c        Echrg (2) has energy lost if proton escapes.
c      Rest of Echrg array will look like:
c
c    Reaction       (3)      (4)       (5)       (6)     SPDSQ   Eneut2
c  --------------------------------------------------------------------
c      p + n      E(11-B)    0.0       0.0       0.0     E(n1)    ---
c     p + 2n      E(10-B)    0.0       0.0       0.0     E(n1)   -E(n2)
c    p + alpha    E(8-Li)  E(alpha)    0.0       0.0      0.0     ---
c  p + n + alpha  E(7-Li)  E(alpha)    0.0       0.0     E(n1)    ---
c  p + 2n+ alpha  E(6-Li)  E(alpha)    0.0       0.0     E(n1)   -E(n2)
c  p+d+2n+2alpha    E(d)   E(alpha)  E(alpha)    0.0     E(n1)   -E(n2)
c  p+t +n+2alpha    0.0    E(alpha)  E(alpha)   E(t)     E(n1)    ---
c      2p         E(11-Be)   0.0      E(p2)    p2-esc     0.0     ---
c    2p + n       E(10-Be)   0.0      E(p2)    p2-esc    E(n1)    ---
c    2p + 2n      E(9-Be)    0.0      E(p2)    p2-esc    E(n1)   -E(n2)
c
      END
c
c------------------------------------------------------------------
C   This is file ND.FOR
c
c       Purpose to compute energies of the deuteron and 11-B coming from
c        the collision  n + 12-C --> d + 11-B.
c         One assumption is that the deuteron "continuum" has the
c         same relative energy distribution as computed for protons by
c         the nuclear model code TNG.  In addition, the programming
c         specifies some amount of ground-state (n,d) collisions for
c         all En, more or less consistent with measured spectra.
c
c       Added in 5/87: Consider possibility of 11-B --> n + 10-B
c         following the (n,d) reaction to a highly-excited "state" in 10-B.
c
c
      Subroutine ND
c
      Real*4 Li7pa
c
c by d.satoh for Linux (2005.01.11 @ jaeri)
c      Common /COLLIS/ Nelm,Echrg(5)
      Common /COLLIS/ Nelm,Echrg(6)
      Common /NEUTRN/ Eneut, U,V,W,X,Y,Z
      Common /MASSES/ Emass(21)
      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
      Common /PPOLAR/ Xmu(37), Pth(37), Npolar
      Common /NAID/   Rad, Ht, Rc, Rz
      Common /RANDM/  Irx
      COMMON/PROTON/NOPROT,RNGEN(144),PRAN(144),DRAN(144),ARAN(144),
     a              HRAN(144),TRAN(144)
c      Common /PROTON/ Noprot,Range(95),Proten(95)
      Common /EXC11B/ Exc11(14),TenBpn,Li7pa
c          Exc11 array tabulated in subroutine NPX
      Common /EXC10B/ Ex10B(9)
c          Ex10B array tabulated in Subroutine NT
ckaji 2011/03/18
      common /RES/ iswt, bias(5), icont, ictrl, iesc, ilight
c
      Dimension Ptheta(37)
c       The next array is used to determine the polar scattering
c        angle of the deuteron and is taken from the measurements
c        in Nuclear Instruments & Methods 129 (1975) 241 for
c        En(lab)=56 MeV.
      Data Ptheta/0.0, 0.013, 0.0596, 0.1318, 0.2128, 0.2959, 0.3768,
     Q  0.4527, 0.522, 0.584, 0.6387, 0.6861, 0.7267, 0.7611, 0.7901,
     P  0.8141, 0.8341, 0.8504, 0.8636, 0.8755, 0.8873, 0.8989,
     N  0.9103, 0.9212, 0.9318, 0.9418, 0.9513, 0.960, 0.9681,
     M  0.9753, 0.9817, 0.9872, 0.9918, 0.9954, 0.9979,0.9995, 1.0/
      Data Q/13.732/, Kn/1/, K12C/7/, Kd/10/, K11B/5/, Nterp/1/
      Data K10B/12/
c
c       Start computation.  Zero charged-particle array.
      Do 1 I=1,6
    1 Echrg(I)=0.0
c       Go to center-of-mass coordinates.  Check to be sure incident
c        neutron energy is large enough.
      En=Eneut
      Eneut=0.0
      Vn=VELOCITY(Kn,En)
      CALL CMTRAN(Kn,K12C,Vn, Vcom,Enc,Ecc)
      Tec=Enc+Ecc
      Ta=Tec-Q
      IF (Ta .GT. 0.0) goto 2
      write(*,333)En
  333 Format(/'   ***  Error in Subroutine ND; E(neut) at entry = '
     a  1PE11.3/10x,'Set E(neut) = E(deuteron) = 0 and exit.'/)
      Return
c
c       Get center-of-mass energies, then transform to lab coord. energies
    2 Edcom=RCKE(Kd,K11B,Ta)
      Edeut=Edcom
      Excit=0.0
      Egamma=0.0
      Eleftov=0.0
      IF (Edcom .LT. 2.22) goto 9
      GSenhanc=0.00543*En
      IF (En .GT. 30.0) GSenhanc=0.08*SQRT(En-25.85)
      if(En .GT. 52.0) GSenHanc=4.2578*En**(-0.59252) ! added by kajimoto 06/08
      IF (RAN(Irx) .LT. GSenhanc) goto 9
c       Last is to enhance ground-state transition, and is ad hoc.
      Temp=0.13*En
      IF (En.LT.30.0) Temp=12.6-1.07*En+0.026*En*En
      if(En.GT.37.0) Temp=3.723e-01*En-8.964 ! added by kajimoto 06/08
c       Ad hoc formula for "nuclear temperature" to be used in CHOOSP
      Edeut=CHOOSP(Edcom,Temp)
c           Assume same "continuum" distribution of deuterons as we
c            have for protons.  It's an assumption.
      Try=Edeut*(Emass(Kd) + Emass(K11B))/Emass(K11B)
      Edtry=RCKE(Kd,K11B,Try)
c       Same "continuum" level excitation energy determination procedure
c         as used for protons in subroutine NPX.
      Try=Try*Edeut/Edtry
      Excit=Ta-Try
c        If variable Excit is relatively small, pair it with a known
c        excited state in 11-B.
      IF (Excit .LT. TenBpn)  goto 4
      IF (Excit .GT. TenBpn+0.5) goto 9
      Excit=TenBpn+0.25
      K=15
      Goto 8
    4 K=1
      Do 17 I=2,14
      IF (Excit .LT. Exc11(I)) goto 18
   17 K=K+1
   18 Excit=Exc11(K)
c       For gamma-ray energy assignments see comment following statement
c         No. 76 in Subroutine NPN
      Egamma=Exc11(K)
      IF (K.EQ.8 .OR. K.EQ.9) Egamma=Exc11(2)
      IF (K.EQ. 11) Egamma=Exc11(3)
      Eleftov=Exc11(K)-Egamma
c
    8 Taex=Ta-Excit
      Edeut=RCKE(Kd,K11B,Taex)
c
    9 Qex=Q + Excit
c
c        Next get polar scattering angle for the deuteron from array
c         Ptheta, above, and array Xmu (=array X in subroutine NP).
      Rand=RAN(Irx)
      Ctheta=EXTERP(Ptheta,Xmu,Rand,Npolar,Nterp)
      IF (Tec .GE. 51.7) goto 20
      Ciso=1.0-2.0*Rand
      Slope=(Ctheta-Ciso)/(51.7-Qex)
c       For En(c.o.m.) < 51.7 MeV interpolate between Ptheta results
c         and isotropy to get Ctheta.
      Ctheta=Slope*(Tec-Qex)
      Ctheta= EXTERP(Ptheta,Xmu,Rand,Npolar,Nterp)
c added by kajimoto 06/13
      if(En .GE. 40.0) then
        Ctheta=akalbach(En,Edeut,13.0,6.0,0.0,1.0,1.0,1.0)
      endif
c added by kajimoto 06/13
c
c       Now get azimuthal scattering and direction cosines
   20 Rann = 6.283185*RAN(Irx)
      Sinphi=SIN(Rann)
      Cosphi=COS(Rann)
      Sintheta=SQRT(1.0 - Ctheta*Ctheta)
      Cx=Sintheta*Cosphi
      Cy=Sintheta*Sinphi
      Cz=Ctheta
      Vdeut=VELOCITY(Kd,Edeut)
      Vdx=Cx*Vdeut
      Vdy=Cy*Vdeut
      Vdz=Cz*Vdeut
      CALL LABTRAN(Vdx,Vdy,Vdz, Vcom, Vx,Vy,Vz)
      Vd=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Edt=EFROMV(Kd,Vd)
c
      Echrg(1)=Edt
c               That's the energy for the deuteron.
c          Now check for deuteron escape in same manner as done for
c                protons.
      Echrg(2)=0.0
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(U,V,W)
c       Xn,Yn,Zn are deuteron's direction cosines in detector coords.
ckaji      IF (Edt .LE. 0.2) goto 10
      IF (Edt .LE. 0.13) goto 10
      Deutpl=PLNGTH(Xn,Yn,Zn, X,Y,Z, Rad,Ht)
      Epofd=Edt
      Rdeut=EXTERP(RNGEN,DRAN,Epofd,Noprot,4)
c        That gets the range of the deuteron using the proton range data.
      RangeQ=Rdeut-Deutpl
ckaji      IF (RangeQ .LE. 0.0) goto 10
      IF (RangeQ .LE. 1.96e-4) goto 10
c        If it gets to next step, deuteron escaped.
      Echrg(2)=EXTERP(DRAN,RNGEN,RangeQ,Noprot,4) ! kajimoto 2011/03/18
c      RangeQ=0.5*RangeQ
c      Eout=2.0*EXTERP(Range,Proten,RangeQ,Noprot,Nterp)
c      Echrg(2)=Eout
c  Okay -- settled deuteron escape.
c       Next: get energy of recoil 11-B ion.
   10 Taex=Tec-Qex
      E11B=RCKE(K11B,Kd,Taex)
      V11B=VELOCITY(K11B,E11B)
      Vbx=-Cx*V11B
      Vby=-Cy*V11B
      Vbz=-Cz*V11B
      CALL LABTRAN(Vbx,Vby,Vbz, Vcom, Vx,Vy,Vz)
      V11B=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
c  Okey to here for straight (n,d).  Now test for subsequent neutron decay,
c       i.e., (n,dn) reaction leading to bound states of 7-Li
      IF (Excit .LE. Li7pa) goto 40
      IF (K .EQ. 10) goto 40
      Tryn=RAN(Irx)
      IF (K.EQ.11 .AND. Tryn.LT.0.9) goto 40
      IF (K.EQ.13 .AND. Tryn.LT.0.0001) goto 40
c       Program counter to here, have alpha or neutron emission.
c       First get 11-B ion into "detector" coordinates.
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(U,V,W)
c       Done.  Xn,Yn,Zn = dir. cosines of moving 11-B ion in "detector"
c         laboratory coordinates.
      IF (K.EQ.15 .AND. Tryn.LT.0.56) goto 14
      IF (Excit .GT. TenBpn) goto 15
c       Program to here, have 11-B --> 7-Li + alpha decay
   14 Egamma=0.0
      Eleftov=0.0
      Ta=Excit-Li7pa
      CALL AP7LI(Ta,V11B,Egamma)
      Goto 45
C       END ALPHA DECAY SECTION
c
c      Program to here have 11-B --> n + 10-B decay to consider.
   15 Ta=Excit-TenBpn
      Encom=RCKE(Kn,K10B,Ta)
      Ehat=En-Q
      Fn=0.065 + 0.001*Ehat
      Temp=Fn*Ehat
c       Assumption of last three statements is that the "continuum"
c         neutron distribution from decay of the excited 11-B ion is the
c         same as used for decay of an excited 12-C ion after adjusting
c         the "incident" neutron energy.
      Eloww=0.0
      Enn=CHOOSN(Eloww,Encom,Temp)
      Try=Enn*(Emass(Kn) + Emass(K10B))/Emass(K10B)
c               Iterate once for improved accuracy on -Try-.
      Entry=RCKE(Kn,K10B,Try)
      Excitn=Ta-Try
      Egamma=0.0
      Eleftov=0.0
      Level=10
      IF (Excitn .GE. 6.11) goto 26
      K=1
      Do 24 I=2,9
      IF (Excitn .LT. Ex10B(I)) goto 25
   24 K=K+1
   25 Excitn=Ex10B(K)
      Level=K
c       For gamma decay of excited states of 10-B see comments following
c         Statement No. 18 in Subroutine NT
      IF (K.GT.1) Egamma=Ex10B(2)
      Eleftov=Ex10B(K)-Egamma
c
c       Now get information on neutron into "detector" laboratory coords.
   26 Taex=Ta-Excitn
      Enn=RCKE(Kn,K10B,Taex)
      CALL RVECT(Xnn,Ynn,Znn)
      Vnn=VELOCITY(Kn,Enn)
      Vxc=Xnn*Vnn
      Vyc=Ynn*Vnn
      Vzc=Znn*Vnn
      CALL LABTRAN(Vxc,Vyc,Vzc, V11B, Vx,Vy,Vz)
      Vneutl=Sqrt(Vx*Vx + Vy*Vy + Vz*Vz)
      Eneut=EFROMV(Kn,Vneutl)
c       That gets outgoing neutron energy.  Now get dir. cosines.
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
c       Done.  Save in /NEUTRN/ Common variable locations.
      U=Xn
      V=Yn
      W=Zn
c       Now get information on 10-B ion.
      E10B=Taex-Enn
      V10B=VELOCITY(K10B,E10B)
      Vbx=-Xnn*V10B
      Vby=-Ynn*V10B
      Vbz=-Znn*V10B + V11B
c       Non-relativistic c.o.m. --> lab. coord. transformation for
c          heavy ion motion.
      V10Blab=SQRT(Vbx*Vbx + Vby*Vby + Vbz*Vbz)
c       Now test for possible (highly-excited) 10-B ion particle decay.
      IF (Level .LE. 6) goto 37
      IF (Level.EQ.8 .AND. RAN(Irx).LT.0.5) goto 37
      Eleftov=0.0
c       Get 10-B ion dir. cosines in "detector" lab. coordinates
      CALL DIRCOS(Vbx,Vby,Vbz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
c       Done.  Now do 10-B decay.
      CALL TENBDK(Excitn,V10Blab,Egamma)
      Goto 45
c
   37 Echrg(3)=EFROMV(K10B,V10Blab)
      Goto 45
c       END SECTION for neutron emission.
c
c       Next statement is the tag end of the (n,d) computation.
   40 EB=EFROMV(K11B,V11B)
      Echrg(3)=EB
c               and that's it for the the 11-B ion
   45 Nelm=9
c
c satoh-test-2002.12.15
c      write(26,*) echrg(2)
c
      CALL BANKR2
c       Now check for possible gamma-ray contributions:
      IF (Egamma .GT. 0.0) CALL PHOTON(Egamma)
c       See Comment end of Subroutine NT for next two statements.
      IF (Egamma .LE. 0.0) Return
      IF (Eleftov .GT. 0.0) CALL PHOTON(Eleftov)
      Return
      END
c------------------------------------------------------------------
C   This is file NT.FOR
c
      Subroutine NT
C  Purpose to compute energies of the triton and 10-B ion coming from
c       the collision reaction  n + 12-C --> t + 10-B which (reaction)


c       has cross sections in the 30 mb region around En = 40 MeV.


c       If the c.o.m. total kinetic energy is not too large, the


c       program will fix the 10-B excitation at one of the known low-


c       lying energy levels.  The triton scattering angle (or more


c       precisely, the cosine of the scattering angle) is chosen by


c       random number from Legendre polynomial fits to data for En


c       between 27 and 61 MeV from Subramanian, et al, Physical


c       Review C28, 521.


c
      Common /COLLIS/ Nelm, Echrg(6)
      Common /NEUTRN/ Eneut, U,V,W,X,Y,Z
      Common /MASSES/ Emass(21)
      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
      Common /RANDM/  Irx
      Common /EXC10B/ Exc10(9)
ckaji 2011/03/18
      COMMON/PROTON/NOPROT,RNGEN(144),PRAN(144),DRAN(144),ARAN(144),
     a              HRAN(144),TRAN(144)
      common /RES/ iswt, bias(5), icont, ictrl, iesc, ilight
      Common /NAID/ Rdet,Ht,RC,Rz
c
      Dimension Epoly(3),B1(3),B2(3),B3(3),B4(3),B5(3),B6(3),Fi(6)
c
c       Next array represents low-lying excited states of 10-B


      Data Exc10/0.0, 0.7183, 1.7402, 2.154, 3.587, 4.774,
     a   5.11, 5.17, 5.92/
c
c       Next arrays are for Legendre polynomials representing angular


c         distributions of tritons.  These distributions are forward


c         peaked at all En.


      Data Ipoly/3/, Epoly/27.4, 39.7, 60.7/, B1/0.42035, 0.3074,
     x  0.3681/, B2/0.209, 0.1459, 0.1619/, B3/0.0925, 0.1106, 0.1018/,
     y  B4/0.03296, 0.0228, 0.03733/, B5/0.01292, -0.01512, -0.00044/,
     z  B6/0.00382, -0.00692, -0.00318/
      Data Nfi/6/, Q/18.93/, Kn/1/, K12C/7/, Kt/13/, K10B/12/,Nterp/1/
c
c   START COMPUTATION.
      Do 1 I=1,6
    1 Echrg(I)=0.0
c       Go to center of mass coordinates.


      En=Eneut
      Eneut=0.0
      Vn=VELOCITY(Kn,En)
      CALL CMTRAN(Kn,K12C,Vn, Vcom,Enc,Ecc)
      Tec=Enc+Ecc
      Ta=Tec-Q
      IF (Ta .GT. 0.0) goto 2
csatoh
c     Type 100, En
      write(*,444)En
  444 Format(/'   ***  Error in Subroutine NT; E(neut) at entry = '
     a  1PE11.3/10x,'Set E(neut) = E(Triton) = 0 and exit.'/)
      Return
c
c   Last is an error return.
c
c      Get center-of-mass energies, then transform to lab coordinate
c       energies.  First set up Legendre coeff. for En.


    2 Fi(1)=EXTERP(Epoly,B1,En,Ipoly,Nterp)
      Fi(2)=EXTERP(Epoly,B2,En,Ipoly,Nterp)
      Fi(3)=EXTERP(Epoly,B3,En,Ipoly,Nterp)
      Fi(4)=EXTERP(Epoly,B4,En,Ipoly,Nterp)
      Fi(5)=EXTERP(Epoly,B5,En,Ipoly,Nterp)
      Fi(6)=EXTERP(Epoly,B6,En,Ipoly,Nterp)
      Etcom=RCKE(Kt,K10B,Ta)
      Etrit=Etcom
      Excit=0.0
      Egamma=0.0
      Eleftov=0.0
      Level=1
      IF (Etcom .LT. 0.72) goto 9
ckaji      IF (RAN(Irx) .LT. 0.05) goto 9
      IF (RAN(Irx) .LT. 0.06) goto 9
c   Last is to enhance slightly the ground-state transition,
c               and is ad hoc.


      Level=10
ckaji      Temp=0.5+0.04*En
      Temp=6.1370e-02*En+1.3271 ! modified by kajimoto 06/08
c         Assume same "continuum" distribution of tritons as is used for


c          protons.  It's an assumption!


      Etrit=CHOOSP(Etcom,Temp)
      Try=Etrit*(Emass(Kt) + Emass(K10B))/Emass(K10B)
      Etry=RCKE(Kt,K10B,Try)
c       Same level excitation energy determination procedure as used


c         for deuterons in subroutine ND.


      Try=Try*Etrit/Etry
      Excit=Ta-Try
c         If variable Excit is relatively small, pair it with a known


c          excited state in 10-B


      IF (Excit .GT. 6.11) goto 9
      K=1
      Do 17 I=2,9
      IF (Excit .LT. Exc10(I)) goto 18
   17 K=K+1
   18 Excit=Exc10(K)
      Level=K
c       Gamma decay of levels in 10-B strongly proceed via the first-excited


c        state of 10-B.  This means at least two photons per decay of an


c        excited state, one of which corresponds to the ground-state


c        transition of the first excited state.


      IF (K .GT. 1) Egamma=Exc10(2)
      Eleftov=Exc10(K)-Egamma
      Taex=Ta-Excit
      Etrit=RCKE(Kt,K10B,Taex)
c
    9 Qex=Q + Excit
c       Next get polar scattering angle for the triton from array Fi


      Ctheta=CHOOSL(Fi,Nfi)
c added by kajimoto 06/13
      if(En .GE. 40.0) then
        Ctheta=akalbach(En,Etrit,13.0,6.0,0.0,1.0,1.0,2.0)
      endif
c added by kajimoto 06/13
c       Now get azimuthal scattering angle.


   20 Rann = 6.283185*RAN(Irx)
      Sinphi=SIN(Rann)
      Cosphi=COS(Rann)
      Sintheta=SQRT(1.0 - Ctheta*Ctheta)
      Cx=Sintheta*Cosphi
      Cy=Sintheta*Sinphi
      Cz=Ctheta
      Vtri=VELOCITY(Kt,Etrit)
      Vtx=Cx*Vtri
      Vty=Cy*Vtri
      Vtz=Cz*Vtri
      CALL LABTRAN(Vtx,Vty,Vtz, Vcom, Vx,Vy,Vz)
      Vt=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Ett=EFROMV(Kt,Vt)
      Echrg(1)=Ett
ckaji      IF (Ett .LE. 1.0) goto 10
      IF (Ett .LE. 0.2) goto 10
c          Eprot < 100 keV, range a few micrometers.
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(U,V,W)
c       Xn,Yn,Zn are the proton's direction cosines in lab coordinates.
      Protpl=PLNGTH(Xn,Yn,Zn, X,Y,Z, Rdet,Ht)
      Rprot=EXTERP(RNGEN,TRAN,Ett,Noprot,4)


c       Now compare proton's range in detector (RPROT) with the path length
c         of the vector having direction cosines Xn, Yn, Zn from position
c         X, Y, Z (PROTPL).
      RangeQ=Rprot-Protpl
      IF (RangeQ .LE. 3.02E-04) goto 10
c         If it goes to next step, the proton got out of the detector.
      if(iesc.ne.2) then ! kajimoto 2011/03/18
        Echrg(1)=EXTERP(TRAN,RNGEN,RangeQ,Noprot,4)
      else
        Echrg(1)=0.0
      endif ! kajimoto 2011/03/18
c
c       That's the energy for the triton.


c         Triton escape is assumed to be negligable.


c       So next get the energy of the recoil 10-B ion.


   10 Taex=Tec-Qex
      E10B=RCKE(K10B,Kt,Taex)
      V10B=VELOCITY(K10B,E10B)
      Vbx=-Cx*V10B
      Vby=-Cy*V10B
      Vbz=-Cz*V10B + Vcom
c         Last transforms to "neutron" Lab. coords. non-relativistically.


      V10B=SQRT(Vbx*Vbx + Vby*Vby + Vbz*Vbz)
c
c       Test for possible (highly-excited) 10-B ion particle decay.


      IF (Level .LE. 6) goto 12
      IF (Level.EQ.8 .AND. RAN(Irx).LT.0.5) goto 12
      Eleftov=0.0
c       Get 10-B ion direction cosines in "detector" lab. coordinates.


      CALL DIRCOS(Vbx,Vby,Vbz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(U,V,W)
c       Done.  Now do decay of highly excited 10-B ion.


      CALL TENBDK(Excit,V10B,Egamma)
      Goto 16
c
   12 Echrg(3)=EFROMV(K10B,V10B)
c         and that's it for the 10-B ion.


   16 Nelm=11
      CALL BANKR2
c       Now check for possible gamma-ray contribution.


      IF (Egamma .GT. 0.0) CALL PHOTON(Egamma)
c       If "P.S.D." is "on" (see Subroutine Input), and if the last call


c         to PHOTON triggered an event subtraction, then on return


c         Egamma < 0.0; so Return.  Otherwise check for possible second


c         gamma ray.


      IF (Egamma .LE. 0.0) Return
      IF (Eleftov .GT. 0.0) CALL PHOTON(Eleftov)
      Return
      End
c------------------------------------------------------------------
c   This is file N2N.FOR
c
c       Purpose is to compute the energies and directions of scatter
c        of the two neutrons and the energy of the 11-C ion from the
c        reaction  n + 12-C --> n' + n' + 11-C.
c
c       New features added 4/87: (1) fix excitation of residual 11-C to
c        be less than (a third) neutron separation energy; (2) allow
c        residual 11-C to decay by proton emission, if there is a
c        sufficient amount of excitation energy; and (3) in any case,
c        check for possible gamma-ray decay of particle-bound excited
c        states in either 11-C or 10-B, the latter if proton emission
c        is allowed.
c
      Subroutine N2N
C
      Common /RANDM/  IRX
      Common /COLLIS/ Nelm,Echrg(6)
      Common /NEUTR2/ Eneut2, U2,V2,W2, X2,Y2,Z2
      Common /NAID/   Rad,Ht,Rc,Rz
      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
      Common /NEUTRN/ SPDSQ, U, V, W, Xneut, Yneut, Zneut
      Common /MASSES/ Emass(21)
      Common /EXC10B/ Exc10(9)
c          Exc10 array in Subroutine NT
      COMMON /PROTON/NOPROT,RNGEN(144),PRAN(144),DRAN(144),ARAN(144),
     a              HRAN(144),TRAN(144)
c      Common /PROTON/ Noprot,Range(95),Proten(95)
      Common /PPOLAR/ X(37),Pt(37),N37
ckaji 2011/03/18
      common /RES/ iswt, bias(5), icont, ictrl, iesc, ilight
c
      Dimension Ptheta(37), Exc11C(11)
c       Next array is used to determine the polar scattering angle for
c        the first neutron.  It was taken from the formalism of the O5S
c        code, and represents a mildly forward-peaked scattering.  The
c        X-array which is its companion is given in Subroutine NP.
      Data Ptheta/0.0, 0.0332, 0.06605, 0.0986, 0.1309, 0.1628,
     A 0.1944, 0.2258, 0.2568, 0.2875, 0.3179, 0.348, 0.3778, 0.4073,
     B 0.4364, 0.4653, 0.4938, 0.5221, 0.550, 0.5776, 0.6049, 0.6319,
     C 0.6586, 0.685, 0.7111, 0.7369, 0.7624, 0.7875, 0.8124, 0.8369,
     D 0.8611, 0.885, 0.9086, 0.9319, 0.9549, 0.9776, 1.0/
c       Next array gives energies of particle-stable excited states in 11-C
      Data Exc11C/0.0, 2.0, 4.32, 4.8, 6.34, 6.48, 6.9, 7.5, 8.1, 8.42,
     a   8.66/
      Data Q/18.72/, Nterp/1/, TenCpn/13.124/, TenBpp/8.691/
      Data Kn/1/, Kp/2/, K12C/7/, K11C/9/, K10B/12/
c
C   The reaction scheme for purposes of computation is:
c
c       n + 12-C --> n + 12-C;      then:
c           12-C --> n + 11-C;
c               We take E(neutron) from CHOOSN.
c
c       Addition 4/87.  Consider possibility of 11-C --> p + 10-B
c
c               Start computation: go to Center-of-mass coordinates
      En=SPDSQ
      Vn = VELOCITY(Kn,En)
c -------------------------------------------------
c check by d.satoh (2005.01.06@jaeri)
c      write(17,*) 'En = ', En
c      write(17,*) 'Kn = ', Kn
c      write(17,*) 'Vn = ', Vn
c -------------------------------------------------
      CALL CMTRAN(Kn,K12C,Vn, Vcom,Eneut,Ecc)
      Tec=Ecc+Eneut
      Ta=Tec-Q
      Do 80 I=1,6
   80 Echrg(I)=0.0
      Epl=0.0
c
      IF (Ta .GT. 0.0) goto 2
c -------------------------------------------------
c check by d.satoh (2005.01.06@jaeri)
c      write(17,*) 'N1 = ', Kn
c      write(17,*) 'N2 = ', K12C
c      write(17,*) 'V  = ', Vn
c      write(17,*) 'Vc = ', Vcom
c      write(17,*) 'E1 = ', Eneut
c      write(17,*) 'E2 = ', Ecc
c      write(17,*) 'Tec= ', Tec
c      write(17,*) 'Ta = ', Ta
c      write(17,*) '  '
c      write(17,*) '----------------'
c -------------------------------------------------
csatoh
c     Type 100, En
      write(*,555)En
  555 Format(/'   ***  Error in Subroutine N2N; E(neut) at entry ='
     b 1PE11.3/10x,'Set E(neut1) = E(neut2) = 0 and exit.'/)
      SPDSQ=0.0
      Eneut2=0.0
      Return
c
C              CHOOSE C.M. NEUTRON ENERGY


    2 Encom=RCKE(Kn,K12C,Ta)
c    Encom is the maximum energy the neutron can have in the center of mass
c         for the (n,2n) reaction.
      Eloww=0.0
      F=0.065+0.001*En
      Temp=F*En
c          Temp value chosen so that the "continuum" neutron spectrum calcu-
c            lated in CHOOSN is reasonably representative of that computed by
c            the nuclear model code TNG.
      Enc1 = CHOOSN(Eloww,Encom,Temp)
c        The difference in energy between Encom and Enc1 could be
c         interpreted as the added energy (larger than Q) of a "virtual"
c         energy level in the 12-C nucleus decaying by neutron emission.
      Try=Enc1*(Emass(Kn)+Emass(K12C))/Emass(K12C)
c               (That's approximately the "total available energy" that
c                would yield Enc1 from the function RCKE.)
      Entry=RCKE(Kn,K12C,Try)
      Try=Try*Enc1/Entry
c               (Iterate once on TRY for slightly more accuracy.)
      Excit=Ta-Try
c         (EXCIT  should be close to the added "virtual" excitation energy.)
      Qex=Q + Excit
c
C         FIND THE COSINE OF THE POLAR DIRECTION OF THE First Neutron (C.M.)
   15 RAND=RAN(IRX)
   18 CTHETA = EXTERP(Ptheta,X,Rand,N37,Nterp)
C          FIND THE SINE AND COSINE OF THE AZIMUTHAL DIRECTION OF THE
C           Neutron (C.M.) -- By Random-number Choice.
      RANN = 6.2831853*RAN(Irx)
      SINPHI = SIN(RANN)
      COSPHI = COS(RANN)
C         FIND THE C. M. VELOCITY COMPONENTS OF THE Neutron.
      Sintheta=SQRT(1.-Ctheta*Ctheta)
      Xp=Sintheta*Cosphi
      Yp=Sintheta*Sinphi
      Zp=Ctheta
      Vnr=VELOCITY(Kn,Enc1)
      Vnx=Xp*Vnr
      Vny=Yp*Vnr
      Vnz=Zp*Vnr
      CALL LABTRAN(Vnx,Vny,Vnz, Vcom, Vx,Vy,Vz)
      Vneut=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
C         FIND THE LAB ENERGY OF THE Neutron.
      Enl=EFROMV(Kn,Vneut)
c        Store information
      SPDSQ=Enl
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Zx=U
      Zy=V
      Zz=W
      CALL TRANSVEC(Zx,Zy,Zz)
      U=Xn
      V=Yn
      W=Zn
c       That takes care of neutron number one.
c
C         FIND THE COMPONENTS OF THE 12-C VELOCITY IN THE LAB.
c           First adjust total available energy for the energy of the
c                 "virtual" excited state.
   70 Taex=Tec-Qex
      E12C=RCKE(K12C,Kn,Taex)
      V12C=VELOCITY(K12C,E12C)
      Vbcx=-Xp*V12C
      Vbcy=-Yp*V12C
      Vbcz=-Zp*V12C
c          (Recall that the direction cosines of the 12-C ion in the center
c           of mass are negatives of the direction cosines of the neutron.)
      CALL LABTRAN(Vbcx,Vbcy,Vbcz, Vcom, Vbx,Vby,Vbz)
      V12C=SQRT(Vbx*Vbx + Vby*Vby + Vbz*Vbz)
c
c    Now worry about the second neutron energy.  The maximum center-of-mass
c       kinetic energy available is given by the variable EXCIT where now the
c       center of mass is the 12-C ion, and its motion is the center-of-mass
c       motion.  Need to get motion of 12-C ion in detector coordinates.
c
      CALL DIRCOS(Vbx,Vby,Vbz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
c       Now variables Xn,Yn,Zn = direction cosines of the motion of the
c          12-C in "detector" laboratory coordinate system.
c
      Ta=Excit
      Encom=RCKE(Kn,K11C,Ta)
      Ehat=En-Q
      F=0.065+0.001*Ehat
      Temp=F*Ehat
c   Now the added assumption -- following the neutron decay of the highly
c       excited 12-C nucleus, the residual 11-C ion is stable against further
c       neutron emission.  However, the possibility of proton decay of the
c       11-C to (as it turns out) a particle-bound state in 10-B will
c       be included.
c
      Tryy = Encom-TenCpn
      IF (Tryy .GT. 0.0) Eloww=Tryy
      Enn=CHOOSN(Eloww,Encom,Temp)
c
c      AND, once again, when choosing from a distribution the difference
c       in energy between Encom and Enn could be interpreted as a "Virtual"
c       energy level in 11-C.
c
      Try=Enn*(Emass(Kn)+Emass(K11C))/Emass(K11C)
      Entry=RCKE(Kn,K11C,Try)
      Try=Try*Enn/Entry
      Excitn=Ta-Try
c
c       If variable -Excitn- is small enough, pair it with a known particle-
c         stable state in 11-C
      Egamma=0.0
      IF (Excitn .GT. TenBpp) goto 9
      K=1
      Do 7 I=2,11
      IF (Excitn .LT. Exc11C(I)) goto 8
    7 K=K+1
    8 Excitn=Exc11C(K)
      Egamma=Excitn
      Taex=Ta-Excitn
      Enn=RCKE(Kn,K11C,Taex)
    9 Continue
c
c               SECOND NEUTRON EMISSION
c       For second neutron emission, the center-of-mass scattering is
c          isotropic.
      CALL RVECT(Xp,Yp,Zp)
      Vnn=VELOCITY(Kn,Enn)
      Vncx=Xp*Vnn
      Vncy=Yp*Vnn
      Vncz=Zp*Vnn
      CALL LABTRAN(Vncx,Vncy,Vncz, V12C, Vnx,Vny,Vnz)
      Vneut=SQRT(Vnx*Vnx + Vny*Vny + Vnz*Vnz)
      Eneut=EFROMV(Kn,Vneut)
      Eneut2=Eneut
c          Now get "new" neutron direction cosines in the
c           "detector" laboratory coordinates.
      CALL DIRCOS(Vnx,Vny,Vnz, Xpn,Ypn,Zpn)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
      U2=Xn
      V2=Yn
      W2=Zn
c       now get the 11-C energy.
      Taex=Ta-Excitn
      E11C=RCKE(K11C,Kn,Taex)
      V11C=VELOCITY(K11C,E11C)
      Vbx=-Xp*V11C
      Vby=-Yp*V11C
      Vbz=-Zp*V11C+V12C
c         For Vbz including "LABTRAN" nonrelativistically.
c          Now have components of 11-C velocity in 12-C ion laboratory
c          coordinate system.
      V11C=SQRT(Vbx*Vbx + Vby*Vby + Vbz*Vbz)
      E11C=EFROMV(K11C,V11C)
      IF (Excitn .LE. TenBpp) goto 30
c
c               PROTON EMISSION
c         If program gets to here, 11-C is energetic enough to decay by
c        proton emission (actually, also by alpha emission, but we don't
c        attempt that case) to particle-stable levels in 10-B.
c       First get 11-C motion into "detector" laboratory coordinate system.
c        This step is needed to determine possible proton leakage from the
c        detector later in the program.
      CALL DIRCOS(Vbx,Vby,Vbz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
c         Done!  Xn,Yn,Zn = direction cosines of 11-C in "detector"
c           laboratory coordinate system.
      Ta=Excitn-TenBpp
      Epcom=RCKE(Kp,K10B,Ta)
      F=0.1245+0.001*ABS(En-45.0)
      Temp=F*En
      Eprot=CHOOSP(Epcom,Temp)
c         and as done before, one iteration for slightly improved accuracy..
      Try=Eprot*(Emass(Kp) + Emass(K10B))/Emass(K10B)
      Eptry=RCKE(Kp,K10B,Try)
      Try=Try*Eprot/Eptry
      Excit=Ta-Try
c       Now pair with an excited state in 10-B (see Function NT).
      Eleftov=0.0
      Egamma=0.0
      K=1
c       By present construction only 5 levels in 10-B are energetically
c          available, and they are all particle stable.
      Do 22 I=2,5
      IF (Excit .LT. Exc10(I)) goto 23
   22 K=K+1
   23 Excit=Exc10(K)
c       See comment following statement no. 18 in Subroutine NT for


c         gamma decay of levels in 10-B.


      IF (K .GT. 1) Egamma=Exc10(2)
      Eleftov=Exc10(K)-Egamma
      Taex=Ta-Excit
      Eprot=RCKE(Kp,K10B,Taex)
c
c       Get proton velocity in 11-C center-of-mass coordinates.


   25 Vprot=VELOCITY(Kp,Eprot)
c       Done.  "Fission" of 11-C --> p + 10-B is isotropic, perforce.


      CALL RVECT(Xp, Yp, Zp)
      Vpx=Xp*Vprot
      Vpy=Yp*Vprot
      Vpz=Zp*Vprot
      CALL LABTRAN(Vpx,Vpy,Vpz, V11C, Vx,Vy,Vz)
      Vprotl=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
c       Vprotl=velocity of proton in 11-C laboratory coordinates.


      Epl=EFROMV(Kp,Vprotl)
c       Store energy information; check for proton escape from detector.


      Echrg(2)=Epl
      Echrg(3)=0.0
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
c       Now Xn,Yn,Zn = direction cosines of proton in detector coordinates.

      IF (Epl .LE. 0.07) goto 34
      Protpl=PLNGTH(Xn,Yn,Zn, Xneut,Yneut,Zneut, Rad,Ht)
c       Get proton range to compare with length of path to detector surface.


      Rprot=EXTERP(RNGEN,PRAN,Epl,Noprot,4)
      RangeQ=Rprot-Protpl
      IF (RangeQ .LE. 1.02e-4) goto 34
c       If program gets to here proton escapes.

      Eprotf=EXTERP(PRAN,RNGEN,RangeQ,Noprot,4) ! kajimoto 2011/03/18
      Echrg(3)=Eprotf
   34 E10B=Taex-Eprot
      V10B=VELOCITY(K10B,E10B)
      Vbcx=-Xp*V10B
      Vbcy=-Yp*V10B
      Vbcz=-Zp*V10B +V11C
c       again, non-relativistic transformation, this time to 11-C laboratory


c        coordinate system.


      V10Blab=SQRT(Vbcx*Vbcx + Vbcy*Vbcy + Vbcz*Vbcz)
      E10B=EFROMV(K10B,V10Blab)
      Echrg(1)=E10B
      Goto 31
   30 Echrg(1)=E11C
   31 Nelm=8
      CALL BANKR2
      IF (Egamma .GT. 0.0) CALL PHOTON(Egamma)
c       See comment regarding "P.S.D." following statement no 16 in


c         Subroutine NT.


      IF (Egamma .LE. 0.0) Return
      IF (Eleftov .GT. 0.0) CALL PHOTON(Eleftov)
      Return
      END
c------------------------------------------------------------------
C   This is file NP.FOR
c
c   Purpose is to compute the energy of the outgoing proton from the
c       reaction  n + 12-C --> p + 12B and to check for its possible


c       escape from the detector and account for the energy loss if it does.


c       The ground state and first four excited states of 12-B are stable


c       against particle emission; the programming chooses which of these


c       five states is the state excited by the incident neutron from a


c       table of probabilities generated from cross sections derived from


c       the Hauser-Feshback code TNG developed by C. Y. Fu.


c
      Subroutine NP
C
      Common /RANDM/  Irx
      Common /COLLIS/ Nelm,Echrg(6)
      Common /NAID/   Rad,Ht,Rc,Rz
      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
      Common /NEUTRN/ SPDSQ, U, V, W, Xneut, Yneut, Zneut
      Common /MASSES/ Emass(21)
      Common /PROTON/NOPROT,RNGEN(144),PRAN(144),DRAN(144),ARAN(144),
     a              HRAN(144),TRAN(144)
c      Common /PROTON/ Noprot,Range(95),Proten(95)
      Common /PPOLAR/ X(37),Ptheta(37),Npolar
ckaji 2011/03/18
      common /RES/ iswt, bias(5), icont, ictrl, iesc, ilight
C

      Dimension Ei(12),Pgs(12),P1ex(12),P2ex(12),P3ex(12),Qex(5)
      Dimension Px(4)
      Equivalence (Qex(1),Q)
c
c       Next are Q values for the 12-B ground and first 4 excited states.


      Data Qex/12.613, 13.566, 14.287, 15.233, 15.333/
c
c       Next arrays are used to determine which bound state of 12-B was


c         populated.  Pgs, P1ex, P2ex and P3ex are cumulative probabilities.


      Data Ei/14.75, 15.5, 16.5, 17.0, 18.0, 18.25, 19.0, 20.0, 21.0,
     Z  23.0, 25.0, 40.0/
      Data Pgs/1.0, 0.691, 0.468, 0.368, 0.287, 0.279, 0.246, 0.224,
     Z  0.215, 0.203, 0.198, 0.180/
      Data P1ex/1.0, 1.0, 0.85, 0.779, 0.658, 0.646, 0.592, 0.552,
     Z  0.53, 0.511, 0.502, 0.475/
      Data P2ex/1.0, 1.0, 1.0, 0.985, 0.904, 0.894, 0.857, 0.827,
     Z  0.808, 0.79, 0.782, 0.762/
      Data P3ex/1.0, 1.0, 1.0, 0.998, 0.976, 0.974, 0.966, 0.959,
     Z  0.955, 0.95, 0.947, 0.939/
c
      Data Npolar/37/
c       Next array is used also in Subroutines NPN and N2N.


      Data X /1.0, .99619, .98481, .965926, .939693,
c               This table is COS(Theta) for theta every 5 degrees.


     &  0.90631, .866025, .81915, .76604, .707107, .64279, .57358,
     &  0.5, .42262, .34202, .25882, .17305, .081557, 0.0, -.081557,
     &  -.17305, -.25882, -.34202, -.42262, -0.5, -.57358, -.64279,
     &  -.707107, -.76604, -.81915, -.866025, -.90631, -.939693,
     &  -.965926, -.98481, -.99619, -1.0/
c
c       The next array is used to determine the polar scattering angle


c        of the proton, and is taken from measurements reported in


c        Nuclear Instruments and Methods 129 (1975) 241.  These data are


c        for En(lab)=56 MeV corresponding to En(c.o.m.)=51.7 MeV; the


c        measured angular distribution is strongly forward peaked for


c        all measured proton energies.


      Data Ptheta/0.0, 0.0091, 0.0385, 0.091, 0.1601, 0.2333, 0.3069,
     S 0.3782, 0.4456, 0.5078, 0.5644, 0.6151, 0.660, 0.6994, 0.7335,
     T 0.7629, 0.7879, 0.8091, 0.8269, 0.842, 0.857, 0.8717, 0.8861,
     U 0.9001, 0.9134, 0.9262, 0.9382, 0.9495, 0.9595, 0.9687, 0.9768,
     V 0.9838, 0.9896, 0.9941, 0.9974, 0.9993, 1.0/
C
      Data Npx/12/, Nterp/1/
      Data Kn/1/, K12C/7/, Kp/2/, K12B/6/
c
c               Go to Center-of-mass coordinates


      En=SPDSQ
      Vn = VELOCITY(Kn,En)
      CALL CMTRAN(Kn,K12C,Vn, Vcom,Eneut,Ecc)
      Tec=Ecc+Eneut
      Ta=Tec-Q
      IF (Ta .GT. 0.0) goto 2
csatoh
c     Type 100, En
      write(*,100)En
  100 Format(/'   ***  Error in Subroutine NP; E(neut) at entry ='
     b 1PE11.3/10x,'Set E(Neut) = E(Proton) = 0 and exit'/)
      SPDSQ=0.0
      Echrg(1)=0.0
      Return
c
C              CHOOSE C.M. PROTON ENERGY


    2 Px(1)=EXTERP(Ei,Pgs,En,Npx,Nterp)
      Px(2)=EXTERP(Ei,P1ex,En,Npx,Nterp)
      Px(3)=EXTERP(Ei,P2ex,En,Npx,Nterp)
      Px(4)=EXTERP(Ei,P3ex,En,Npx,Nterp)
      Eran=RAN(Irx)
      Eleftov=0.0
      K=1
      Do 4 J=1,4
      IF (Eran .LE. Px(J)) goto 6
    4 K=K+1
    6 Ta=Tec-Qex(K)
      Egamma=Qex(K)-Qex(1)
      IF (K .NE. 4) goto 7
c       Third excited state of 12-B decays primarily to the first excited


c         state.  Other excited states decay primarily to the ground state.


      Egamma=Qex(2)-Qex(1)
      Eleftov=Qex(4)-Qex(2)
    7 Eprot=RCKE(Kp,K12B,Ta)
c
C         FIND THE COSINE  OF THE POLAR DIRECTION OF THE PROTON(C.M.)
   15 Rand=RAN(Irx)
   18 Ctheta = EXTERP(Ptheta,X,Rand,Npolar,Nterp)
c added by kajimoto 06/13
      if(En .GE. 40.0) then
        Ctheta=akalbach(En,Eprot,13.0,6.0,0.0,1.0,1.0,0.0)
        goto 20
      endif
c added by kajimoto 06/13
c       Assumed angular distribution at En(c.o.m.) .GE. 51.7 MeV, but
c        modified toward isotropy for En(c.o.m.) < 51.7 MeV.
      IF (Tec .GE. 51.7) goto 20
      Ciso=1.0-2.0*Rand
      Slope=(Ctheta-Ciso)/(51.7-Qex(K))
      Ctheta=Slope*Ta
C          FIND THE SINE AND COSINE OF THE AZIMUTHAL DIRECTION OF THE


C           PROTON (C.M.) -- By Random-number Choice.


   20 Rann = 6.2831853*RAN(Irx)
      Sinphi = SIN(Rann)
      Cosphi = COS(Rann)
C         FIND THE C. M. DIRECTION COSINES OF THE PROTON.


      Sintheta=SQRT(1.-Ctheta*Ctheta)
      Xp=Sintheta*Cosphi
      Yp=Sintheta*Sinphi
      Zp=Ctheta
      Vpr=VELOCITY(Kp,Eprot)
C        NOW GET C. M. VELOCITY COMPONENTS OF THE PROTON.


      Vpx=Xp*Vpr
      Vpy=Yp*Vpr
      Vpz=Zp*Vpr
      CALL LABTRAN(Vpx,Vpy,Vpz, Vcom, Vx,Vy,Vz)
C        THAT HAS GOT THE VELOCITY COMPONENTS OF THE PROTON IN THE
C               "NEUTRON" LABORATORY COORDINATE SYSTEM.


      Vprot=SQRT(Vx*Vx +Vy*Vy +Vz*Vz)
C         FIND THE LAB ENERGY OF THE PROTON.


      Epl=EFROMV(Kp,Vprot)
c        Store information; then get distance from interaction point to


c         detector surface in the direction of the proton's velocity.


      Echrg(1)=Epl
      Echrg(2)=0.0
c       First check for very low energy proton; skip if so.


      IF (Epl .LE. 0.07) goto 70
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(U,V,W)
      Protpl=PLNGTH(Xn,Yn,Zn, Xneut,Yneut,Zneut, Rad,Ht)
C         FIND RANGE OF THE PROTON.

      Rprot = EXTERP(RNGEN,PRAN,Epl,Noprot,4)
C  DETERMINE THE AMOUNT OF THIS RANGE THAT LIES WITHIN THE SCINTILLATOR.
      RangeQ=Rprot-Protpl
      IF (RangeQ .LE. 1.02e-4) goto 70
c               If it gets to here, proton "escapes" ...


      Eprotf=EXTERP(PRAN,RNGEN,RangeQ,Noprot,4) ! kajimoto 2011/03/18
      Echrg(2)=Eprotf
c
C         FIND THE COMPONENTS OF THE 12B VELOCITY IN THE LAB.


   70 E12B=RCKE(K12B,Kp,Ta)
      V12B=VELOCITY(K12B,E12B)
      Vbcx=-Xp*V12B
      Vbcy=-Yp*V12B
      Vbcz=-Zp*V12B
c          (Recall that the direction cosines of the 12-B ion in the center


c           of mass are negatives of the direction cosines of the proton.)


      CALL LABTRAN(Vbcx,Vbcy,Vbcz, Vcom, Vbx,Vby,Vbz)
      V12B=SQRT(Vbx*Vbx + Vby*Vby + Vbz*Vbz)
C         FIND THE LAB ENERGY OF THE 12B.


      E12B = EFROMV(K12B,V12B)
      Echrg(3)=E12B
      Nelm=6
      SPDSQ=0.0
      CALL BANKR2
c       Check for possible gamma-ray interactions:


      IF (Egamma .GT. 0.0) CALL PHOTON(Egamma)
      IF (Egamma .LE. 0.0) Return
      IF (Eleftov .GT. 0.0) CALL PHOTON(Eleftov)
      Return
      END
c------------------------------------------------------------------
c   This is file CSCAT.FOR
c
      Subroutine CSCAT
c
c         Finish up  n + 12-C  elastic and inelastic scattering.
c
      Common /CSCATT/ Xp,Yp,Zp, Vnca,Vcca,Vcom
      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
      Common /COLLIS/ Nelm,Echrg(6)
      Common /NEUTRN/ Eneut,U,V,W,X,Y,Z
c
      Data Kn/1/, K12C/7/
c
c       Get neutron velocity components in "neutron" center-of-mass


c         coordinate system:


      Vnx=Xp*Vnca
      Vny=Yp*Vnca
      Vnz=Zp*Vnca
c       Get, similarly, carbon-ion velocity components:


      Vcx=-Xp*Vcca
      Vcy=-Yp*Vcca
      Vcz=-Zp*Vcca
c
c       Transform neutron velocity components from "neutron" center-of-mass


c         coordinate system to "neutron" laboratory coordinate system:


      CALL LABTRAN(Vnx,Vny,Vnz, Vcom, Vx,Vy,Vz)
c       Now get V(neutron) and E(neutron).  V(neutron) has the same value


c         in "neutron" lab. coordinates as in "detector" lab. coordinates.


      Vneut=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Eneut=EFROMV(Kn,Vneut)
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
c       Xpn,Ypn and Zpn are neutron's dir. cosines in "neutron" lab. coords.


c       Rotate neutron's dir. cosines into "detector" lab. coordinates


      CALL TRANSVEC(U,V,W)
      U=Xn
      V=Yn
      W=Zn
c          and the information is saved in common block labelled NEUTRN.


c       Now get information on the Carbon ion:


      CALL LABTRAN(Vcx,Vcy,Vcz, Vcom, Vx,Vy,Vz)
      Vcarb=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Ecarb=EFROMV(K12C,Vcarb)
      Echrg(1)=Ecarb
c               All Done!


      Return
      END
c------------------------------------------------------------------
C  This is file PROTN.FOR
cc
c       Contains the block data for computing proton ranges.
c               (Additions 6/87.)
c
c       Modified by d.satoh to expand data tables. Extra data were
c       calculated by Bethe formula. (2004.12.15@JAERI)
c
c      BLOCK DATA PROTN
c      COMMON/PROTON/NOPROT,RANGE(95),PROTEN(95)
c      DATA NOPROT/95/
c       Ranges are in cm for NE-213.  These values are changed in the INPUT
c         subroutine for NE-110.
c      DATA (RANGE(I),I=1,95)/ 2.500E-5, 1.0400E-4,
c     a 1.5500E-4, 1.9500E-4, 2.3500E-4, 2.9500E-4, 4.5000E-4, 6.6000E-4,
c     1 9.0000E-4, 1.1500E-3, 1.4400E-3, 1.8000E-3, 2.1500E-3, 2.5000E-3,
c     2 3.5000E-3, 4.5000E-3, 5.6000E-3, 6.9000E-3, 8.3000E-3, 9.7000E-3,
c     4 0.0113,    0.013,     0.0148,    0.0167,    0.0187,    0.0208,
c     5 0.0230,    0.0252,    0.0276,    0.03,      0.0325,    0.0352,
c     6 0.0379,    0.0407,    0.0461,    0.0466,    0.0496,    0.0528,
c     7 0.0559,    0.0627,    0.0698,    0.0772,    0.0849,    0.0929,
c     7 0.1012,    0.1098,    0.1142,    0.1255,    0.1374,    0.1906,
c     8 0.2517,    0.3202,    0.3958,    0.4785,    0.568,     0.6642,
c     9 0.7672,    0.8768,    0.993,     1.1062,    1.2349,    1.3699,
c     A 1.5110,    1.6583,    1.8117,    1.9711,    2.2000,    2.3085,
c     B 2.4856,    2.6683,    2.8568,    3.0509,    3.2506,    3.4557,
c     C 3.6664,    3.8824,    4.1039,    4.3306,    4.5626,    8.6700,
c     + 10.330,   12.0400,   13.8700,   15.8000,   17.8300,   19.9500,
c     + 22.170,   24.4800,   26.8800,   29.3600,   31.9200,   34.5600,
c     + 37.270,   40.0600,   42.9300/
c
c       Next are proton energies, in MeV.
c      DATA (PROTEN(I),I=1,95)/0.01,     0.06,
c     a 0.1,   0.13,  0.16,  0.2,   0.3,   0.4,   0.5,   0.6,   0.7,
c     1 0.8,   0.9,   1.0,   1.2,   1.4,   1.6,   1.8,   2.0,   2.2,
c     2 2.4,   2.6,   2.8,   3.0,   3.2,   3.4,   3.6,   3.8,   4.0,
c     3 4.2,   4.4,   4.6,   4.8,   5.0,   5.2,   5.4,   5.6,   5.8,
c     4 6.0,   6.4,   6.8,   7.2,   7.6,   8.0,   8.4,   8.8,   9.0,
c     5 9.5,  10.0,  12.0,  14.0,  16.0,  18.0,  20.0,  22.0,  24.0,
c     6 26.0, 28.0,  30.0,  32.0,  34.0,  36.0,  38.0,  40.0,  42.0,
c     7 44.0, 46.0,  48.0,  50.0,  52.0,  54.0,  56.0,  58.0,  60.0,
c     8 62.0, 64.0,  66.0,  68.0,  70.0, 100.0,
c     + 110., 120., 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0,
c     + 200., 210., 220.0, 230.0, 240.0, 250.0/
c      END
c
      BLOCK DATA PROTN
      COMMON/PROTON/NOPROT,RNGEN(144),PRAN(144),DRAN(144),ARAN(144),
     a              HRAN(144),TRAN(144)
      DATA NOPROT/144/
c       Ranges are in cm for NE-213.  These values are changed in the INPUT
c         subroutine for NE-110.
      DATA (PRAN(I),I=1,144)/
     A 2.4680e-05,2.6490e-05,2.8240e-05,2.9930e-05,3.1590e-05,
     A 3.3200e-05,3.4770e-05,3.6310e-05,3.7820e-05,4.0760e-05,
     A 4.4310e-05,4.7740e-05,5.1070e-05,5.4320e-05,5.7510e-05,
     A 6.0630e-05,6.3710e-05,6.6740e-05,7.2710e-05,7.8570e-05,
     A 8.4360e-05,9.0110e-05,9.5830e-05,1.0200e-04,1.1300e-04,
     A 1.2500e-04,1.3600e-04,1.4800e-04,1.6100e-04,1.7300e-04,
     A 1.8600e-04,2.0000e-04,2.1300e-04,2.2700e-04,2.4200e-04,
     A 2.7200e-04,3.1200e-04,3.5500e-04,4.0000e-04,4.4800e-04,
     A 4.9900e-04,5.5100e-04,6.0600e-04,6.6400e-04,7.8600e-04,
     A 9.1700e-04,1.0560e-03,1.2050e-03,1.3610e-03,1.5260e-03,
     A 1.8780e-03,2.2610e-03,2.6720e-03,3.1100e-03,3.5700e-03,
     A 4.0570e-03,4.5720e-03,5.1130e-03,5.6820e-03,6.2770e-03,
     A 6.8990e-03,8.2200e-03,1.0017e-02,1.1973e-02,1.4085e-02,
     A 1.6351e-02,1.8769e-02,2.1336e-02,2.4050e-02,2.6911e-02,
     A 3.3058e-02,3.9772e-02,4.7043e-02,5.4862e-02,6.3222e-02,
     A 7.2115e-02,9.1455e-02,1.1300e-01,1.3600e-01,1.6200e-01,
     A 1.8900e-01,2.1800e-01,2.5000e-01,2.8300e-01,3.1800e-01,
     A 3.5400e-01,3.9300e-01,4.7600e-01,5.8900e-01,7.1300e-01,
     A 8.4700e-01,9.9200e-01,1.1470e+00,1.3120e+00,1.4870e+00,
     A 1.6720e+00,2.0700e+00,2.5050e+00,2.9760e+00,3.4830e+00,
     A 4.0250e+00,4.6000e+00,5.8470e+00,7.2200e+00,8.7150e+00,
     A 1.0327e+01,1.2051e+01,1.3882e+01,1.5818e+01,1.7855e+01,
     A 1.9988e+01,2.2216e+01,2.4534e+01,2.9427e+01,3.6003e+01,
     A 4.3056e+01,5.0550e+01,5.8454e+01,6.6739e+01,7.5380e+01,
     A 8.4351e+01,9.3632e+01,1.1300e+02,1.3300e+02,1.5500e+02,
     A 1.7700e+02,2.0000e+02,2.2300e+02,2.7200e+02,3.2200e+02,
     A 3.7400e+02,4.2700e+02,4.8000e+02,5.3500e+02,5.9000e+02,
     A 6.4600e+02,7.0200e+02,7.5900e+02,8.1600e+02,9.3100e+02,
     A 1.0750e+03,1.2190e+03,1.3640e+03,1.5090e+03/
C
C      deuteron
      DATA (DRAN(I),I=1,144)/
     A 3.0780e-05,3.3170e-05,3.5480e-05,3.7720e-05,3.9900e-05,
     A 4.2030e-05,4.4100e-05,4.6130e-05,4.8110e-05,5.1950e-05,
     A 5.6560e-05,6.0980e-05,6.5240e-05,6.9350e-05,7.3350e-05,
     A 7.7230e-05,8.1020e-05,8.4730e-05,9.1920e-05,9.8860e-05,
     A 1.0600e-04,1.1200e-04,1.1900e-04,1.2500e-04,1.3700e-04,
     A 1.4900e-04,1.6100e-04,1.7300e-04,1.8400e-04,1.9600e-04,
     A 2.0700e-04,2.1900e-04,2.3000e-04,2.4200e-04,2.5300e-04,
     A 2.7700e-04,3.0700e-04,3.3800e-04,3.7000e-04,4.0400e-04,
     A 4.3800e-04,4.7400e-04,5.1100e-04,5.4900e-04,6.3000e-04,
     A 7.1600e-04,8.0600e-04,9.0200e-04,1.0030e-03,1.1090e-03,
     A 1.3340e-03,1.5780e-03,1.8410e-03,2.1210e-03,2.4170e-03,
     A 2.7310e-03,3.0600e-03,3.4060e-03,3.7670e-03,4.1430e-03,
     A 4.5340e-03,5.3580e-03,6.4600e-03,7.6390e-03,8.9010e-03,
     A 1.0248e-02,1.1680e-02,1.3195e-02,1.4793e-02,1.6472e-02,
     A 2.0071e-02,2.3988e-02,2.8219e-02,3.2758e-02,3.7600e-02,
     A 4.2741e-02,5.3895e-02,6.6208e-02,7.9655e-02,9.4218e-02,
     A 1.1000e-01,1.2700e-01,1.4400e-01,1.6300e-01,1.8300e-01,
     A 2.0400e-01,2.2600e-01,2.7300e-01,3.3700e-01,4.0800e-01,
     A 4.8400e-01,5.6600e-01,6.5400e-01,7.4800e-01,8.4800e-01,
     A 9.5300e-01,1.1790e+00,1.4270e+00,1.6970e+00,1.9870e+00,
     A 2.2980e+00,2.6280e+00,3.3480e+00,4.1450e+00,5.0160e+00,
     A 5.9600e+00,6.9750e+00,8.0590e+00,9.2110e+00,1.0428e+01,
     A 1.1711e+01,1.3056e+01,1.4463e+01,1.7457e+01,2.1525e+01,
     A 2.5940e+01,3.0688e+01,3.5753e+01,4.1121e+01,4.6781e+01,
     A 5.2719e+01,5.8926e+01,7.2092e+01,8.6211e+01,1.0100e+02,
     A 1.1700e+02,1.3400e+02,1.5100e+02,1.8700e+02,2.2600e+02,
     A 2.6700e+02,3.1000e+02,3.5400e+02,4.0000e+02,4.4700e+02,
     A 4.9500e+02,5.4400e+02,5.9400e+02,6.4400e+02,7.4800e+02,
     A 8.8000e+02,1.0160e+03,1.1540e+03,1.2930e+03/
C
C      triton
      DATA (TRAN(I),I=1,144)/
     A 3.4360e-05,3.7150e-05,3.9860e-05,4.2490e-05,4.5060e-05,
     A 4.7560e-05,5.0010e-05,5.2390e-05,5.4730e-05,5.9270e-05,
     A 6.4700e-05,6.9910e-05,7.4930e-05,7.9770e-05,8.4450e-05,
     A 8.9000e-05,9.3430e-05,9.7750e-05,1.0600e-04,1.1400e-04,
     A 1.2200e-04,1.2900e-04,1.3700e-04,1.4400e-04,1.5700e-04,
     A 1.7100e-04,1.8300e-04,1.9600e-04,2.0800e-04,2.2000e-04,
     A 2.3200e-04,2.4400e-04,2.5600e-04,2.6700e-04,2.7900e-04,
     A 3.0200e-04,3.3100e-04,3.5900e-04,3.8900e-04,4.1800e-04,
     A 4.4800e-04,4.7900e-04,5.1000e-04,5.4200e-04,6.0900e-04,
     A 6.7800e-04,7.5100e-04,8.2800e-04,9.0700e-04,9.9100e-04,
     A 1.1670e-03,1.3570e-03,1.5610e-03,1.7770e-03,2.0060e-03,
     A 2.2480e-03,2.5020e-03,2.7680e-03,3.0450e-03,3.3340e-03,
     A 3.6340e-03,4.2670e-03,5.1180e-03,6.0340e-03,7.0120e-03,
     A 8.0510e-03,9.1420e-03,1.0282e-02,1.1475e-02,1.2725e-02,
     A 1.5394e-02,1.8289e-02,2.1406e-02,2.4742e-02,2.8295e-02,
     A 3.2061e-02,4.0218e-02,4.9202e-02,5.8998e-02,6.9590e-02,
     A 8.0965e-02,9.3113e-02,1.0600e-01,1.2000e-01,1.3400e-01,
     A 1.4900e-01,1.6500e-01,1.9900e-01,2.4500e-01,2.9600e-01,
     A 3.5100e-01,4.1000e-01,4.7400e-01,5.4100e-01,6.1200e-01,
     A 6.8800e-01,8.5100e-01,1.0290e+00,1.2230e+00,1.4310e+00,
     A 1.6540e+00,1.8920e+00,2.4100e+00,2.9850e+00,3.6130e+00,
     A 4.2950e+00,5.0300e+00,5.8160e+00,6.6520e+00,7.5370e+00,
     A 8.4710e+00,9.4520e+00,1.0481e+01,1.2674e+01,1.5665e+01,
     A 1.8927e+01,2.2449e+01,2.6222e+01,3.0238e+01,3.4491e+01,
     A 3.8971e+01,4.3672e+01,5.3706e+01,6.4551e+01,7.6163e+01,
     A 8.8500e+01,1.0200e+02,1.1500e+02,1.4400e+02,1.7600e+02,
     A 2.0900e+02,2.4400e+02,2.8100e+02,3.2000e+02,3.6000e+02,
     A 4.0100e+02,4.4400e+02,4.8700e+02,5.3200e+02,6.2300e+02,
     A 7.4300e+02,8.6600e+02,9.9300e+02,1.1230e+03/
c
C      alpha
      DATA (ARAN(I),I=1,144)/
     A 1.8460e-05,2.0110e-05,2.1720e-05,2.3300e-05,2.4840e-05,
     A 2.6360e-05,2.7840e-05,2.9300e-05,3.0730e-05,3.3510e-05,
     A 3.6860e-05,4.0080e-05,4.3190e-05,4.6190e-05,4.9110e-05,
     A 5.1950e-05,5.4710e-05,5.7410e-05,6.2640e-05,6.7660e-05,
     A 7.2520e-05,7.7230e-05,8.1800e-05,8.6250e-05,9.4820e-05,
     A 1.0300e-04,1.1100e-04,1.1800e-04,1.2600e-04,1.3300e-04,
     A 1.3900e-04,1.4600e-04,1.5200e-04,1.5900e-04,1.6500e-04,
     A 1.7700e-04,1.9100e-04,2.0400e-04,2.1800e-04,2.3000e-04,
     A 2.4300e-04,2.5500e-04,2.6800e-04,2.8000e-04,3.0300e-04,
     A 3.2700e-04,3.5000e-04,3.7300e-04,3.9600e-04,4.2000e-04,
     A 4.6700e-04,5.1500e-04,5.6500e-04,6.1700e-04,6.7000e-04,
     A 7.2400e-04,7.8000e-04,8.3800e-04,8.9800e-04,9.6000e-04,
     A 1.0230e-03,1.1550e-03,1.3300e-03,1.5160e-03,1.7140e-03,
     A 1.9220e-03,2.1420e-03,2.3740e-03,2.6160e-03,2.8700e-03,
     A 3.4100e-03,3.9940e-03,4.6220e-03,5.2930e-03,6.0080e-03,
     A 6.7650e-03,8.4050e-03,1.0198e-02,1.2134e-02,1.4226e-02,
     A 1.6474e-02,1.8875e-02,2.1426e-02,2.4126e-02,2.6973e-02,
     A 2.9966e-02,3.3102e-02,3.9797e-02,4.8950e-02,5.8961e-02,
     A 6.9817e-02,8.1504e-02,9.4009e-02,1.0700e-01,1.2100e-01,
     A 1.3600e-01,1.6800e-01,2.0400e-01,2.4200e-01,2.8300e-01,
     A 3.2700e-01,3.7400e-01,4.7600e-01,5.9000e-01,7.1400e-01,
     A 8.4900e-01,9.9400e-01,1.1500e+00,1.3160e+00,1.4920e+00,
     A 1.6770e+00,1.8720e+00,2.0770e+00,2.5140e+00,3.1110e+00,
     A 3.7630e+00,4.4690e+00,5.2270e+00,6.0360e+00,6.8940e+00,
     A 7.8000e+00,8.7530e+00,1.0793e+01,1.3008e+01,1.5390e+01,
     A 1.7930e+01,2.0623e+01,2.3462e+01,2.9547e+01,3.6149e+01,
     A 4.3230e+01,5.0753e+01,5.8686e+01,6.7001e+01,7.5671e+01,
     A 8.4673e+01,9.3984e+01,1.0400e+02,1.1300e+02,1.3400e+02,
     A 1.6100e+02,1.8900e+02,2.1800e+02,2.4800e+02/
C
C     3Helium
      DATA (HRAN(I),I=1,144)/
     A 1.7550e-05,1.9050e-05,2.0520e-05,2.1940e-05,2.3340e-05,
     A 2.4700e-05,2.6030e-05,2.7340e-05,2.8620e-05,3.1110e-05,
     A 3.4100e-05,3.6970e-05,3.9750e-05,4.2440e-05,4.5050e-05,
     A 4.7590e-05,5.0070e-05,5.2500e-05,5.7200e-05,6.1730e-05,
     A 6.6110e-05,7.0340e-05,7.4450e-05,7.8430e-05,8.6080e-05,
     A 9.3350e-05,1.0000e-04,1.0700e-04,1.1300e-04,1.2000e-04,
     A 1.2600e-04,1.3200e-04,1.3700e-04,1.4300e-04,1.4800e-04,
     A 1.5900e-04,1.7200e-04,1.8500e-04,1.9700e-04,2.0900e-04,
     A 2.2100e-04,2.3200e-04,2.4400e-04,2.5600e-04,2.7900e-04,
     A 3.0200e-04,3.2500e-04,3.4900e-04,3.7300e-04,3.9800e-04,
     A 4.4800e-04,5.0100e-04,5.5600e-04,6.1300e-04,6.7200e-04,
     A 7.3400e-04,7.9800e-04,8.6400e-04,9.3300e-04,1.0040e-03,
     A 1.0780e-03,1.2320e-03,1.4380e-03,1.6590e-03,1.8950e-03,
     A 2.1460e-03,2.4120e-03,2.6920e-03,2.9870e-03,3.2960e-03,
     A 3.9580e-03,4.6770e-03,5.4520e-03,6.2840e-03,7.1660e-03,
     A 8.0930e-03,1.0098e-02,1.2311e-02,1.4727e-02,1.7343e-02,
     A 2.0154e-02,2.3159e-02,2.6354e-02,2.9738e-02,3.3307e-02,
     A 3.7060e-02,4.0996e-02,4.9401e-02,6.0900e-02,7.3488e-02,
     A 8.7146e-02,1.0200e-01,1.1800e-01,1.3400e-01,1.5200e-01,
     A 1.7100e-01,2.1100e-01,2.5600e-01,3.0400e-01,3.5600e-01,
     A 4.1200e-01,4.7100e-01,6.0000e-01,7.4300e-01,9.0000e-01,
     A 1.0700e+00,1.2530e+00,1.4490e+00,1.6570e+00,1.8780e+00,
     A 2.1110e+00,2.3560e+00,2.6120e+00,3.1590e+00,3.9060e+00,
     A 4.7200e+00,5.5980e+00,6.5400e+00,7.5430e+00,8.6040e+00,
     A 9.7220e+00,1.0896e+01,1.3401e+01,1.6108e+01,1.9007e+01,
     A 2.2088e+01,2.5340e+01,2.8755e+01,3.6035e+01,4.3877e+01,
     A 5.2228e+01,6.1042e+01,7.0279e+01,7.9902e+01,8.9879e+01,
     A 1.0000e+02,1.1100e+02,1.2200e+02,1.3300e+02,1.5600e+02,
     A 1.8500e+02,2.1600e+02,2.4800e+02,2.8000e+02/
C
c       Next are proton energies, in MeV.
      DATA (RNGEN(I),I=1,144)/
     A 1.0000e-02,1.1000e-02,1.2000e-02,1.3000e-02,1.4000e-02,
     A 1.5000e-02,1.6000e-02,1.7000e-02,1.8000e-02,2.0000e-02,
     A 2.2500e-02,2.5000e-02,2.7500e-02,3.0000e-02,3.2500e-02,
     A 3.5000e-02,3.7500e-02,4.0000e-02,4.5000e-02,5.0000e-02,
     A 5.5000e-02,6.0000e-02,6.5000e-02,7.0000e-02,8.0000e-02,
     A 9.0000e-02,1.0000e-01,1.1000e-01,1.2000e-01,1.3000e-01,
     A 1.4000e-01,1.5000e-01,1.6000e-01,1.7000e-01,1.8000e-01,
     A 2.0000e-01,2.2500e-01,2.5000e-01,2.7500e-01,3.0000e-01,
     A 3.2500e-01,3.5000e-01,3.7500e-01,4.0000e-01,4.5000e-01,
     A 5.0000e-01,5.5000e-01,6.0000e-01,6.5000e-01,7.0000e-01,
     A 8.0000e-01,9.0000e-01,1.0000e+00,1.1000e+00,1.2000e+00,
     A 1.3000e+00,1.4000e+00,1.5000e+00,1.6000e+00,1.7000e+00,
     A 1.8000e+00,2.0000e+00,2.2500e+00,2.5000e+00,2.7500e+00,
     A 3.0000e+00,3.2500e+00,3.5000e+00,3.7500e+00,4.0000e+00,
     A 4.5000e+00,5.0000e+00,5.5000e+00,6.0000e+00,6.5000e+00,
     A 7.0000e+00,8.0000e+00,9.0000e+00,1.0000e+01,1.1000e+01,
     A 1.2000e+01,1.3000e+01,1.4000e+01,1.5000e+01,1.6000e+01,
     A 1.7000e+01,1.8000e+01,2.0000e+01,2.2500e+01,2.5000e+01,
     A 2.7500e+01,3.0000e+01,3.2500e+01,3.5000e+01,3.7500e+01,
     A 4.0000e+01,4.5000e+01,5.0000e+01,5.5000e+01,6.0000e+01,
     A 6.5000e+01,7.0000e+01,8.0000e+01,9.0000e+01,1.0000e+02,
     A 1.1000e+02,1.2000e+02,1.3000e+02,1.4000e+02,1.5000e+02,
     A 1.6000e+02,1.7000e+02,1.8000e+02,2.0000e+02,2.2500e+02,
     A 2.5000e+02,2.7500e+02,3.0000e+02,3.2500e+02,3.5000e+02,
     A 3.7500e+02,4.0000e+02,4.5000e+02,5.0000e+02,5.5000e+02,
     A 6.0000e+02,6.5000e+02,7.0000e+02,8.0000e+02,9.0000e+02,
     A 1.0000e+03,1.1000e+03,1.2000e+03,1.3000e+03,1.4000e+03,
     A 1.5000e+03,1.6000e+03,1.7000e+03,1.8000e+03,2.0000e+03,
     A 2.2500e+03,2.5000e+03,2.7500e+03,3.0000e+03/
      END
c
