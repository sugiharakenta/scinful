c   This is File DOMEGA.FOR
c       Provide Initial Information for Monte-Carlo
c       Response Routine for Cylindrical Detector.
c
c     Right-Circular Cylinder of radius RAD and length HT.
c
c     Purposes: (a) Deduce subtended solid angle; and
c               (b) Determine and record coordinates and direction
c                   cosines of initial interactions in the detector
c                   for NHIST events.
c
      FUNCTION DOMEGA(NHIST)
      Double Precision Rdp,Dsdp,Dsdpp,D1,Xsdp,Fdp
      Common /TRY/  TRIES(8,1000000)
      Common /NAID/  RAD,HT,RCOLIM,R
      Common /INIT/ NH,Nhitot
      Common /RANDM/  IXX
cs
      Common /MAXWEL/ Esourc,E1,T
cs
      Data Pi/3.1415927/
      Data Nmax/1000000/, D1/1.0D+0/
c
c  Steps in calculation are as follows:
c
c  1.  Determine a comparatively small solid angle which completely encloses
c      the detector as observed by the source.  This is the TRIAL solid angle.
c  2.  By random choice, obtain some fraction of the overall TRIAL solid angle.
c      (The source is assumed to radiate isotropically into this solid angle.)
c  3.  Pick a direction (vector) by random number.
c  4.  Test to see if this vector from the source intercepts some place
c      on the detector surface that is "seen" by the source.
c      a.  If it does, record the coordinates of this spot and deduce and
c          record the direction cosines of the vector from the source spot.
c          Keep this information in the  TRIES(I,J)  array.
c      b.  If it doesnt, try again.  Keep track of total number of attempts.
c      c.  If it does, then see if there was an interaction.
c      d.  If no interaction, try again.
c  5.  Keep choosing until the total number of interactions equals the
c      initially desired number of histories.
c  6.  The solid angle subtended by the detector is the number of hits
c      divided by the total number of attempts times the  TRIAL  solid angle.
c  7.  At exit the  TRIES(I,J)  array has information on the position in the
c      detector for the initial NHIST interactions in the detector.
c
c       TRIES(1,J) = X Coordinate
c       TRIES(2,J) = Y Coordinate
c       TRIES(3,J) = Z Coordinate
c       TRIES(4,J) = X Direction Cosine
c       TRIES(5,J) = Y Direction Cosine
c       TRIES(6,J) = Z Direction Cosine
c       TRIES(7,J) = Energy of the interacting neutron
c       TRIES(8,J) = Index of the type of interaction
c
c  The Cartesian coordinate system has its origin at the center of the front
c      face of the detector such that the positive Z axis is into the
c      detector along the cylindrical symmetry axis.  The back surface of the
c      detector is at Z = +HT.  The Front Face of the detector is in the X-Y
c      plane.
c
c      If the position of the source is offset from the Z axis, the
c      offset is taken as along the positive X axis, i.e., Y=0 by symmetry.
c


cs -----------------------------------------------------------------------------
c Initialize (add by d.satoh '05.02.01)
      Nint=0
      Nhitot=0
c
      R = 0.
      NINDX = 0
      SINZ = 0.
      AZ = 0.
      H = 0.
      RD = 0.
      XD = 0.
      XP = 0.
      YP = 0.
      PHI = 0.
      RPRIME = 0.
      DWP = 0.
      NCOL = 0
      NHIT = 0
      NTOT = 0
      CHIMAX = 0.
      DOMSS = 0.
      DOMRR = 0.
      DOM = 0.
      F = 0.
      DS = 0.
      FDP = 0.
      XSDP = 0.
      DSDPP = 0.
      DSDP = 0.
      RDP = 0.
c
cs -----------------------------------------------------------------------------
      Zstart=Tries(3,1)
      Xstart= SQRT(Tries(1,1)**2 + Tries(2,1)**2)
c
      Ntries=Nhist
      IF (Ntries. GT. Nmax)  Ntries=Nmax
      Ncontl=1
   3  IF (Zstart) 4,10,20
c
c   If ZSTART < 0.0 then the front face of the detector is illuminated.
c
   4  IF (Xstart .GT. Rad)  go to 40
c
c      If the Offset of the Source is > detector radius, then the cylindrical
c         surface of the detector is illuminated.  However, if a collimator
c       is designated, the collimator is assumed to cover all of the
c       cylindrical surface and perhaps also some of the front face of the
c       detector.  For ease in computation the collimator is taken to be
c       infinitely thin and totally absorbing.
c
c  Geometry of calculation in the next portions of the program:
c
CCCCCCCCC  FRONT FACE    FRONT FACE    FRONT FACE    FRONT FACE  CCCCCCCCC
C                                                                        C
C                         DETECTOR.DETECTOR.DETECTOR.DETECTOR.D          C
C                         F                                   D          C
C   SOURCE                R                                   D          C
C     :                   O                                   D          C
C     :                   N                                   D          C
C     --------------------T - - - - - - - - - - - - - - - - - D---> +Z   C
C                         F                                   D          C
C                         A                                   D          C
C                         C                                   D          C
C                         E                                   D          C
C                         DETECTOR.DETECTOR.DETECTOR.DETECTOR.D          C
C                                                                        C
C     I<------ DS ------->I<----------------- HT ------------>I          C
C                                                                        C
CCCCCCCCC  FRONT FACE    FRONT FACE    FRONT FACE    FRONT FACE  CCCCCCCCC
c
   2  Ds=-Zstart
      R=Rad
      IF (Rcolim.LT.Rad .AND. Rcolim.GT.0.0)  R=Rcolim
   1  Rdp=R
c           Selective double-precision arithmetic here
      Dsdp=Ds
        Dsdpp=Ds
      Xsdp=Xstart
      Fdp=D1-Dsdp/DSQRT(Dsdp*Dsdp + (Rdp+Xsdp)**2)
      F=Fdp
      Dom=2.0*Pi*F
c           Solid angle = 2*Pi*(1-COS(ATAN(Rad/Distance)))
      Domrr=Dom
      Domss=0.0
      Chimax=Pi
      IF (Xstart .LE. R) go to 9
      Fdp=D1-Dsdp/DSQRT(Dsdp*Dsdp + (Xsdp-Rdp)**2)
      F=Fdp
      Domss=2.*Pi*F
      Domrr=Dom-Domss
      Chimax=ASIN(R/Xstart)
        Dom=Domrr*Chimax/Pi
    9 IF (Ncontl .NE. 1)  go to 41
      Ntot=1
      Nhit=1
      Ncol=1
   5  Dwp=0.5*Domrr*RAN(IXX)/Pi
cs
c        write(80,*)'RAN--------->',RAN(IXX)
      Dwp=Dwp+0.5*Domss/Pi
      Fdp=Dwp
c           Selective double precisions here, too
      Rdp=Dsdpp*DSQRT(D1/(D1-Fdp)**2 - D1)
      Rprime=Rdp
      Phi=2.*Chimax*(0.5-RAN(Ixx)) + Pi
      Yp=Rprime*SIN(Phi)
      IF (ABS(Yp) .GT. R)  go to 6
      Xp=Rprime*COS(Phi)
      Xd=Xp+Xstart
      Rd=SQRT(Yp*Yp + Xd*Xd)
      IF (Rd .GT. R)  go to 6
      H=SQRT(Ds*Ds + Rprime*Rprime)
      Az=Ds/H
      Sinz=SQRT(1.-Az*Az)
      Tries(1,Ncol)=Xd
      Tries(2,Ncol)=Yp
      Tries(3,Ncol)=0.0
      Tries(4,Ncol)=Sinz*COS(Phi)
      Tries(5,Ncol)=Sinz*SIN(Phi)
      Tries(6,Ncol)=Az
cs
c        write(80,*) 'X= ',Tries(1,Ncol)
c        write(80,*) 'Y= ',Tries(2,Ncol)
c        write(80,*) 'Z= ',Tries(3,Ncol)


   8  IF (Ncontl .NE. 1) go to 49
      Nhit=Nhit+1
cs
c_expansion
c        if (Esourc.gt.80) then
        if (Esourc.gt.150.0) then
        call INTERACT2(Ncol,Nindx)
        Nint=Nint+1
        else
      Call INTERACT(Ncol,Nindx)
        end if
        IF (Nindx .LE. 0) goto 6
        Ncol=Ncol+1
      If (Ncol.GT.Ntries) go to 37
   6  Ntot=Ntot+1
      Go to 5
c
c      If the program goes to the next statement it's because the initially-
c       chosen Z=0 puts the source on the plane of the face of the detector.
c
  10  IF (Xstart .GT. Rad)  go to 30
c
c      If XSTART is smaller than RAD then the source is ON the front face of
c         the detector.  Need only to choose initial direction cosines.  In
c       this case if the Z-Direction Cosine is chosen negative, can simply
c       reverse direction by 180 degrees with no loss of randomness.
c
      Ntot=1
      Ncol=1
      Nhit=1
  11  CALL RVECT(Ax,Ay,Az)
      IF (Az)  12,15,13
  12  Tries(4,Ncol)=-Ax
      Tries(5,Ncol)=-Ay
      Tries(6,Ncol)=-Az
      Go to 14
  13  Tries(4,Ncol)=Ax
      Tries(5,Ncol)=Ay
      Tries(6,Ncol)=Az
  14  Tries(1,Ncol)=Xstart
      Tries(2,Ncol)=0.0
      Tries(3,Ncol)=0.0
      Nhit=Nhit+1
cs
c_expansion
c        if (Esourc.gt.80) then
        if (Esourc.gt.150.0) then
      Call INTERACT2(Ncol,Nindx)
        Nint=Nint+1
        else
      Call INTERACT(Ncol,Nindx)
        end if


      IF (Nindx .LE. 0) goto 15
      Ncol=Ncol+1
      IF (Ncol .GT. Ntries) go to 16
  15  Ntot=Ntot+1
      Go to 11
c        write(80,*)'Nint16------->',Nint
  16  Domega=2.*Pi*FLOAT(Ntries)/FLOAT(Ntot)
  18  Nhitot=Nhit
      RETURN
c
c  Next is for ZSTART > 0.0, i.e. curved side of the detector is "seen."
c
  20  IF (Zstart .LT. Ht)  go to 21
c
c      If Program goes to here then the source is "behind" the detector.
c       In this case, reverse the coordinates to put the source in front
c       of the detector -- i.e. ignore P.M. tube and environs.
c
      Zstart=Ht-Zstart
      Go to 3
c
  21  IF (Xstart-Rad)  22,23,30
c
c     >>>>>>  IF  XSTART  is smaller than RAD ----
c             then the source is "inside" the detector, and that
c             situation is typed out before continuing.
csatoh
c  22 Type 99, Xstart, Zstart, Rad
22       write(*,99) Xstart, Zstart, Rad


  99  Format(/'  ****  Source Inside Detector'/10x,
     x            'X-START, Z-START, RADIUS ='1P3E10.3)
csatoh
c     Type 98
        write(*,98)
  98  Format(4x,'If Okay Type 1, Else 0'/3H ? )
csatoh
c     Accept 97, Nokay
        read(*,97)Nokay
  97  Format(I1)
      IF (Nokay .NE. 1) STOP 'Fix INPUT.DATA File.'
      Write (21,99) Xstart,Zstart,Rad
      Domega=4.*Pi
      Nhit=1
  96  Tries(1,Ncol)=Xstart
      Tries(2,Ncol)=0.0
      Tries(3,Ncol)=Zstart
      Call RVECT(Ax,Ay,Az)
      Tries(4,Ncol)=Ax
      Tries(5,Ncol)=Ay
      Tries(6,Ncol)=Az
      Nhit=Nhit+1
cs
c_expansion
c        if (Esourc.gt.80) then
        if (Esourc.gt.150.0) then
      Call INTERACT2(Ncol,Nindx)
        Nint=Nint+1
        else
      Call INTERACT(Ncol,Nindx)
        end if


      IF (Nindx .LE. 0) goto 28
      Ncol=Ncol+1
      IF (Ncol .GT. Ntries) goto 18
      Go to 96
cc
cc
cc
c      If XSTART = RADius  then the source is on the curved (cylindrical)
c       surface of the detector.  See Comment following STATEMENT
c       labelled 10, above, only now  AX  must be negative.
c
  23  Ntot=1
      Ncol=1
      Nhit=1
  25  Call RVECT(Ax,Ay,Az)
      IF (Ax) 26,28,24
  24  Tries(4,Ncol)=-Ax
      Tries(5,Ncol)=-Ay
      Tries(6,Ncol)=-Az
      Go to 27
  26  Tries(4,Ncol)=Ax
      Tries(5,Ncol)=Ay
      Tries(6,Ncol)=Az
  27  Tries(1,Ncol)=Xstart
      Tries(2,Ncol)=0.0
      Tries(3,Ncol)=Zstart
      Nhit=Nhit+1
cs
c_expansion
c        if (Esourc.gt.80) then
        if (Esourc.gt.150.0) then
      Call INTERACT2(Ncol,Nindx)
        Nint=Nint+1
        else
      Call INTERACT(Ncol,Nindx)
        end if
      IF (Nindx .LE. 0) goto 28
      Ncol=Ncol+1
      IF (Ncol .GT. Ntries) go to 16
  28  Ntot=Ntot+1
      Go to 25
cc
cc
cc
c  Next is source "seeing" curved (cylindrical) surface of the detector.
c    Example of geometry for side entry follows:
c
CCCCCCCCC  SIDE ENTRY    SIDE ENTRY    SIDE ENTRY    SIDE ENTRY  CCCCCCCCC
C                                                                        C
C            SOURCE                                                      C
C              :                                                         C
C              :                                                         C
C              :                                                         C
C              :                                                         C
C              :                                                         C
C              :                                                         C
C              :                                                         C
C              :                                                         C
C         DETECTOR.DETECTOR.DETECTOR.DETECTOR.DETECTOR.D                 C
C         D                                            D                 C
C         D                                            D                 C
C         D                                            D                 C
C         D                                            D                 C
C     --- * ----------------------------------------   D ---> +Z AXIS    C
C         D                                            D                 C
C         D                                            D                 C
C         D                                            D                 C
C         D                                            D                 C
C         DETECTOR.DETECTOR.DETECTOR.DETECTOR.DETECTOR.D                 C
C                                                                        C
CCCCCCCCC  SIDE ENTRY    SIDE ENTRY    SIDE ENTRY    SIDE ENTRY  CCCCCCCCC
c
c
  30  R1=Ht-Zstart
      IF (R1 .LT. Zstart)  R1=Zstart
      Dss=Xstart-Rad
      Thet=ASIN(Rad/Xstart)
      Ymax=Dss*TAN(Thet)
      Rss=SQRT(R1*R1 + Ymax*Ymax)
c       Selective double-precision arithmetic, again.
        Rdp=Rss
        Dsdp=Dss
        Fdp=D1-Dsdp/DSQRT(Dsdp*Dsdp + Rdp*Rdp)
        F=Fdp
      Dom =2.*Pi*F
        Domr=Dom
        Dom2=0.0
        Phimax=Pi
        IF (Zstart .GE. 0.0) goto 34
c        Again double-precision arithmetic
        Rdp=Zstart
        Dsdp=Xstart
        Fdp=D1-Dsdp/DSQRT(Dsdp*Dsdp + Rdp*Rdp)
        F=Fdp
        Doms=2.*Pi*F
        Domr=Dom-Doms
        Phimax=ATAN(-Rad/Zstart)
        Dom=Domr*Phimax/Pi
  34    IF (Ncontl .NE. 1) go to 42
      Ntot=1
      Ncol=1
      Nhit=1
  31  Dwp=0.5*Domr*RAN(IXX)/Pi
cs
c        write(80,*)'RAN31------>',RAN(IXX)
        Dwp=Dwp+0.5*Doms/Pi
      Rprime=Dss*SQRT(1./(1.-Dwp)**2 -1.)
      Phi=2.*Phimax*(0.5-RAN(Ixx))
      Zbar=Rprime*COS(Phi)
      Zp=Zbar+Zstart
      IF (Zp.LT.0.0) go to 33
      IF (Zp.GT. Ht) go to 33
      Yp=Rprime*SIN(Phi)
      IF (ABS(Yp) .GT. Ymax)  go to 33
      IF (ABS(Yp) .GT. 0.02*Rad)  go to 35
c
c  If program goes to here, entry is close enough to the X axis to
c      take Y=0 for the starting point of entry.
c
      Tries(1,Ncol)=Rad
      Tries(2,Ncol)=0.0
      Tries(3,Ncol)=Zp
      Tries(5,Ncol)=0.0
      V=SQRT(Dss*Dss + Rprime*Rprime)
  32  Tries(4,Ncol)=-Dss/V
      Tries(6,Ncol)=Zbar/V
      IF (Ncontl .NE. 1) go to 46
      Nhit=Nhit+1
cs
c_expansion
c        if (Esourc.gt.80) then
        if (Esourc.gt.150.0) then
      Call INTERACT2(Ncol,Nindx)
        Nint=Nint+1
        else
      Call INTERACT(Ncol,Nindx)
        end if
      IF (Nindx .LE. 0) goto 33
      Ncol=Ncol+1
      IF (Ncol .GT. Ntries) go to 37
  33  Ntot=Ntot+1
      Go to 31
  37  Domega=Dom*FLOAT(Nhit-1)/FLOAT(Ntot)
c
      Nhitot=Nhit
      RETURN
c
  35  Thetpr=ATAN(Yp/Dss)
      Alph=ASIN((Xstart*SIN(Thetpr))/Rad)
      Beta=Alph-Thetpr
      Xs=Rad*COS(Beta)
      Xd=Xstart-Xs
      Ys=Yp*Xd/Dss
      Zs=Zstart+Zbar*Xd/Dss
      IF (Zs .GE. Ht)  go to 33
      IF (Zs .LT.0.0)  go to 33
      Tries(1,Ncol)=Xs
      Tries(2,Ncol)=Ys
      Tries(3,Ncol)=Zs
      V=SQRT(Dss*Dss + Yp*Yp + Zbar*Zbar)
      Tries(5,Ncol)=Yp/V
      Go to 32
c
c
c   For what follows the source is at a position to illuminate
c     the front face of the detector as well as the cylindrical
c     side.  First check to see if there is a collimator.  If
c     there is, assume that the collimator is centered on the
c     Z axis shielding the cylindrical side and perhaps some of
c     the front face.
c
  40  IF (Rcolim .GT. 0.0  .AND.  Rcolim .LE. Rad)  go to 2
c
c   Now to the most general case -- can illuminate the front face
c      AND the cylindrical side.  We will do a different trick here.
c      Instead of trying to deduce an all-encompassing  TRIAL  solid
c      angle, we will get One  TRIAL  solid angle for the detector face
c      and Another for the cylindrical surface using the algebra
c      already in place, above.  Then, for the first "Hit" we choose
c      between the two cases by random number on the basis of the
c      relative sizes of the two  TRIAL  solid angles.  We then strobe
c      whichever is chosen until a "Hit" is registered, keeping track
c      of the total number of tries taken.  Then we do an immediate
c      adjustment of the  TRIAL  solid angle for that case to get
c      a "Running" Value of the solid angle for that case.  This
c      "Running" solid angle will be the one compared with the  TRIAL
c      solid angle of the case that wasnt run to determine which case
c      will get the second "Hit" case.  The idea is that after a few
c      "Hits" the "Running" solid angles will approximate the actual
c      solid angles for each case well enough such that choosing which
c      case to study by random number will represent the actual
c      situation.
c
      Ncontl=2
      Go to 2
c
  41  Domft=Dom
c               That's the  TRIAL  solid angle for the face,
      go to 30
c
  42  Domst=Dom
c               and that's it for the cylindrical side.
c
      Ncol=1
      Nhitf=0
      Nhits=0
      Ntotf=0
      Ntots=0
      Domfr=Domft
      Domsr=Domst
c               DOMFR and DOMSR are the "Running" solid angles
  44  Nhit=Nhitf+Nhits+1
      IF (Nhitf .LE. 0)  go to 45
      Domfr=Domft*FLOAT(Nhitf)/FLOAT(Ntotf)
  45  IF (Nhits .LE. 0)  go to 47
      Domsr=Domst*FLOAT(Nhits)/FLOAT(Ntots)
  47  Ratio =Domfr/(Domfr+Domsr)
      U=RAN(IXX)
      IF (U .LE. Ratio) go to 48
c
c        "Go to 48" means Go To "Face" of the Detector
c        Else GO TO "Side" of the detector, which follows immediately
c
      Ntot=Ntots
      Dom=Domst
      Go to 31
  46  Ntots=Ntot+1
      Nhits=Nhits+1
      Go to 50
c
  48  Dom=Domft
      Ntot=Ntotf
      Go to 5
  49  Ntotf=Ntot+1
      Nhitf=Nhitf+1
cs
c_expansion
c  50    if (Esourc.gt.80) then
  50    if (Esourc.gt.150.0) then
      Call INTERACT2(Ncol,Nindx)
c        write(80,*) 'INTERACT2.1'
        Nint=Nint+1
        else
      Call INTERACT(Ncol,Nindx)
        end if


      IF (Nindx .LE. 0) goto 44
      Ncol=Ncol+1
      IF (Ncol .LE. Ntries) goto 44
      Domega=Domfr+Domsr
      Nhitot=Nhit
c
      Return
      END
