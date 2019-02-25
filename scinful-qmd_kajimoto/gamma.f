*************************************
*  SUBROUTINES                      *
*     betarg                        *
*     comdcs                        *
*     compsig                       *
*     comptn                        *
*     photon                        *
*                                   *
*************************************
c========================================
c    This is the BETARG.FOR file.
c
c       This file was also taken from my Monte Carlo code to compute
c         responses in sodium iodide.  However, see comment below
c         for present use.
c
      Function BETARG(N,Z)
c
C       N=1 FIND RANGE  FOR ENERGY Z (in MeV)
C        =2 FIND ENERGY FOR RANGE Z (in cm)
c
      Common /DENS/ Apb(2)
      Data Nthru/1/, Avoinv/1660.3/
      Data A/-7.139648/,B/2.534449/,C/-.0907533/
c       The above constants do well for electron ranges in sodium iodide,
c         for example.  It is known that electron ranges, expressed in
c         mg/cm**2, vary only slightly with the atomic number, Z, of the
c         stopping medium, and so these constants are approximately
c         correct for the C:H material we are dealing with, but not likely
c         correct to the number of decimal places given.  Calculated
c         ranges don't have to be too accurate, however, since the number
c         of electrons per neutron is likely quite small, on the order of
c         a few percent.
c
      IF (Nthru .NE. 1) goto 4
      R=Apb(1)/Apb(2)
      W=12.0+R
      Spgr=W*Avoinv*Apb(2)
c         Spgr should be the material density in mg/cm**3


      Nthru=2
    4 IF (Z.LE.0.) GO TO 1
      GO TO (2,3,1),N
c_satoh@jaeri
c    1 Type 100,N,Z
    1 write(*,100) N,Z
  100 Format (19H0**Error BETARG, N=I4,6H, ARG=E12.3/)
      BETARG=0.
      Return
c
    2 Zet=1000.0*Z
      ELN=ALOG(Zet)
      RLN=A+ELN*(B+C*ELN)
      BETARG=EXP(RLN)/SPGR
      Return
c
    3 SQ=B*B+4.*C*(ALOG(SPGR*Z)-A)
      IF (Sq .GE. 0.0)  goto 14
c_satoh@jaeri
c      Type 10, Z,Sq
      write(*,10) Z, Sq
  10  Format(//20X,'Problem in -BETARG-   input range='1pe10.3,
     x  ' cm, leading to a negative arg ='e10.3,'.  Set it = 0.0'/)
      Sq=0.0
   14 S=0.5*(SQRT(SQ)-B)/C
      BETARG=0.001*EXP(S)
      Return
      END
c========================================
C   This is file COMDCS.FOR
c
c   The following routine was taken from my Monte Carlo code to compute
c       responses of photons by NaI (or BaF2) scintillation detectors.


c
      Subroutine COMDCS(Egam,Ebeta)
c
C  GIVEN VECTOR R IN DETECTOR COORDINATES FROM X0,Y0,Z0                 COMDC 02
C     TO X,Y,Z, AND ENERGY OF COMPTON SCATTERED BETA,                   COMDC 03
C     COMPUTE POLAR SCATTERING ANGLES FOR BETA AND                      COMDC 04
C     SCATTERED GAMMA, CHOOSE AZIMUTHAL ANGLE BY                        COMDC 05
C     MONTE CARLO, THEN COMPUTE DIRECTION COSINES FOR                   COMDC 06
C     BETA AND GAMMA IN DETECTOR COORDINATES.                           COMDC 07
c
      Real*8 Dr,Dx,Dy,Dz
      COMMON /SPOT/ X,Y,Z,X0,Y0,Z0                                      
      COMMON /CBET/ XB,YB,ZB,XG,YG,ZG                                  
      Common /RANDM/ Ixx
c
      DATA EMCSQ/0.511/, PI/3.14159/
      EG=EGAM                                                           
      EB=EBETA                                                          
      A=EG/EMCSQ                                                        
      A2=(1.+A)**2
      A3=2.*A*EG/EB-1.-2.*A                                             
      TPB=ATAN (SQRT(A3/A2))                                            
      IF (TPB.LT.0.) TPB=TPB+PI                                         
      TPG=ACOS(1.-2./(1.+A3))
      IF (TPG.LT.0.) TPG=TPG+PI                                         
      PPB=2.*PI*RAN(Ixx)                                                
C         TPB and PPB are the polar and azimuthal angles of the scattered
C           electron with respect to the z-axis defined by the direction of
C           incident photon.
      DX=X-X0                                                           
      DY=Y-Y0                                                           
      DZ=Z-Z0                                                           
      DR=DSQRT(DX*DX+DY*DY+DZ*DZ)                                       
      R=Dr
      TR=ACOS(DZ/R)                                                     
      IF (TR.LT.0.) TR=TR+PI                                            
      IF (ABS(DX)  .GT. .001*ABS(DY)) GO TO 1                           
      PR=0.5*PI                                                         
      GO TO 2                                                           
    1 PR=ATAN(DY/DX)                                                    
      IF (PR.LT.0.) PR=PR+PI                                            
    2 IF (DY.LT.0.) PR=PR+PI                                            
c
C      TR and PR are the polar and azimuthal angles of the vector
C        representing the direction of the incident photon with respect
C        to the detector coordinate axis.
c
C      What follows is the transformation to get the vector representing
C        the direction of the scattered electron also into the detector
C        coordinate systems.  Its direction cosines are then determined.
c
      ERI=SIN (TR)*COS (PR)                                             
      ERJ=SIN (TR)*SIN (PR)                                             
      ERK=COS (TR)                                                      
      EPI=-SIN (PR)                                                     
      EPJ=COS (PR)                                                      
      ETI=COS (PR)*COS (TR)                                             
      ETJ=SIN (PR)*COS (TR)                                             
      ETK=-SIN (TR)                                                     
      XB=ERI*COS (TPB)+SIN (TPB)*(EPI*SIN (PPB)+ETI*COS (PPB))          
      YB=ERJ*COS (TPB)+SIN (TPB)*(EPJ*SIN (PPB)+ETJ*COS (PPB))          
      ZB=ERK*COS (TPB)+SIN (TPB)*ETK*COS (PPB)                          
C       And that is it for the electron.  Now do the scattered photon vector.
c
      PPG=PPB-PI                                                        
      XG=ERI*COS (TPG)+SIN (TPG)*(EPI*SIN (PPG)+ETI*COS (PPG))          
      YG=ERJ*COS (TPG)+SIN (TPG)*(EPJ*SIN (PPG)+ETJ*COS (PPG))          
      ZG=ERK*COS (TPG)+SIN (TPG)*ETK*COS (PPG)                          


      RETURN
      END
c========================================
C  This is file COMPSIG.FOR
c
      Function COMPSIG(Egamma)
c
c       Purpose is to return a photon scattering cross section in cm**-1


c
      Common /DENS/ Apb(2)
      Dimension Eg(25),SigH(25),SigC(25)
c
      Data Eg/0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.15,
     x 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0,
     y 6.0, 8.0, 10.0, 15.0/
c       Next arrays are photon scattering cross sections in barns/atom


      Data SigH/0.609, 0.606, 0.592, 0.576, 0.559, 0.544, 0.516, 0.492,
     w 0.443, 0.406, 0.353, 0.317, 0.289, 0.267, 0.235, 0.211, 0.172,
     x 0.146, 0.115, 0.0962, 0.0831, 0.0734, 0.0601, 0.0511, 0.0379/
      Data SigC/13.6, 7.33, 4.4, 3.72, 3.45, 3.3, 3.09, 2.94, 2.66,
     y 2.43, 2.12, 1.9, 1.73, 1.6, 1.41, 1.27, 1.03, 0.88, 0.71, 0.61,
     z 0.54, 0.49, 0.42, 0.38, 0.33/
      Data N/25/, Nterp/4/
c
      E=Egamma
      Hsig=EXTERP(Eg,SigH,E,N,Nterp)
      Csig=EXTERP(Eg,SigC,E,N,Nterp)
      C=Hsig*Apb(1) + Csig*Apb(2)
      Compsig=C
      Return
      END
c========================================
C   This is file COMPTN.FOR
c
      Subroutine COMPTN(Egamma,Ebeta)
C      Routine to obtain electron energy, Ebeta, for an input photon
c       energy, Egamma.  Energies are in MeV.


      Common /RANDM/ Ixx
      Dimension S(100), E(100)
      Data EMCSQ/0.511/, Eold/0.0/, Ndat/100/, Nterp/1/
c
      IF (Egamma .EQ. Eold) goto 2
      Eold=Egamma
      Eg=Egamma
      Alpha=Eg/EMCSQ
      Alphsq=Alpha*Alpha
      Emax=Eg*2.*Alpha/(1.0+2.0*Alpha)
      Edel=Emax/FLOAT(Ndat)
      Halfd=0.5*Edel
      E(1)=0.0
      S(1)=0.0
      DO 1 I=2,Ndat
      Ek=E(I-1)+Halfd
      DE=Eg-Ek
      TODT=Ek/DE
      SIG=(2.+TODT*(TODT/ALPHSQ+Ek/EG-2./ALPHA))/ALPHSQ
      S(I)=S(I-1)+Sig
      E(I)=Ek+Halfd
    1 Continue
      Sum=S(Ndat)
c
    2 W=Sum*RAN(Ixx)
      Ebeta=EXTERP(S,E,W,Ndat,Nterp)
      Return
      END
c========================================
C  This is file PHOTON.FOR
c
c   Purpose is to following scattering of a photon in the detector until
c     the photon escapes or is totally absorbed.
c
      Subroutine PHOTON(Egamma)
c
c   Primary subroutine -- follows photon scattering, and at each interaction
c       determines energy of scattered electron and converts to light units.


c       Also determines if the electron escapes the detector and therefore


c       only part of its energy is deposited in the detector.  Compton


c       scattering only is considered in the scattering process until the


c       photon energy becomes < 30 keV in which case total absorbtion of


c       the remaining energy is assumed.  Pair production is not considered.


c
      Common /NEUTRN/ En, U,V,W,Xneut,Yneut,Zneut
      Common /NAID/ Rdet, Ht, Rc, Rz
      Common /SPOT/ X,Y,Z, Xnew,Ynew,Znew
      Common /CBET/ Bx,By,Bz, Xg,Yg,Zg
      Common /RANDM/ Ixx
      Common /COLLIS/ Nelm,Echrg(6)
      Common /GFLAG/  Igflag
c
      E=Egamma
      X=Xneut
      Y=Yneut
      Z=Zneut
      Rdetsq=Rdet*Rdet
      Amtlt=0.0
c       Start off with assumed isotropic distribution of gamma ray at the


c        neutron interaction spot.  Gamma-ray angular distributions are not,


c        in general, isotropic, but they are symmetric about 90 deg, and


c        the non-isotropy is generally smaller than a factor of two.


      CALL RVECT(Cx,Cy,Cz)
c                    Start loop on photon scattering:


    1 Sig=COMPSIG(E)
      Path=PLNGTH(Cx,Cy,Cz, X,Y,Z, Rdet,Ht)
      D=EXP(-Path*Sig)
      T=RAN(Ixx)
      IF (T .LT. D) goto 5
c          If program counter gets to here, a photon scattering took place.


      IF (E .LT. 0.03) goto 4
      Dzet=-ALOG(T)/Sig
      Xnew=X+Cx*Dzet
      Ynew=Y+Cy*Dzet
      Dsq=Xnew*Xnew + Ynew*Ynew
      IF (Dsq .GT. Rdetsq) goto 5
c       Last 2 steps may be unnecessary check on single-precision arithmetic.


      Znew=Z+Cz*Dzet
      CALL COMPTN(E,Ebeta)
      CALL COMDCS(E,Ebeta)
      Amtlt=0.842*Ebeta + Amtlt
      Bpath=PLNGTH(Bx,By,Bz, Xnew,Ynew,Znew, Rdet,Ht)
      N=1
      Brange=BETARG(N,Ebeta)
      RangeQ=Brange-Bpath
      IF (RangeQ .LE. 0.0) goto 3
c       If program counter gets to here, the electron escaped before losing


c         all of its energy.  (Actually not very accurate, for even


c         fast electrons get banged about and hardly ever travel in a


c         straight line.  But for the present program, this inaccuracy


c         isn't very important.)


      N=2
      Eblost=BETARG(N,RangeQ)
      Amtlt=Amtlt-0.842*Eblost
c       Get scattered photon's energy, reset position coordinates and


c         direction cosines, and go back through the loop.
    3 E=E-Ebeta
      X=Xnew
      Y=Ynew
      Z=Znew
      Cx=Xg
      Cy=Yg
      Cz=Zg
      Goto 1
c       Next is terminal step for photon with energy < 30 KeV.


    4 Amtlt=Amtlt+0.842*E
c       End of routine -- deposit accumulated light in BANKER routine


    5 IF (Amtlt .LE. 0.0) Return
      Nelm=10
      Esave=Echrg(1)
      Echrg(1)=Amtlt
      CALL BANKR2
      Echrg(1)=Esave
c   If "P.S.D." is "on" (see BANKER subroutines) return with a negative
c       value for Egamma as a signature to calling subroutine.


      IF (Igflag .EQ. 1) Egamma=-Egamma
      Return
      END
c========================================
