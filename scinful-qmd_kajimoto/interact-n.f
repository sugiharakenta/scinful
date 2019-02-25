C   This is file TOTALX.FOR
c
c     Purpose is to determine total attenuation cross section (in units of
c      cm**-1) and also the fractional parts of this total due to the separate
c      partial cross sections for incident neutron of energy Eneut.
c
      Function TOTALX(Eneut)
c
      real*8 emev, sigct, sigcr,sigcs
      integer incp,ia,iz
c
      Common /DENS/ Apb(2)
      Common /PROB/ P(12), indx
      Dimension S(10)
c
      En=Eneut
      S(1)=SIGHYD(En)*Apb(1)
c -----------------------------------------------------------------------
c modified by d.satoh (06.01.25)
c
      if( En .le. 110.0 ) then
        S(2)=Apb(2)*SIGCELAS(En)
      else
        incp = 2              ! neutron
        emev = dble(En) ! incident energy (MeV)
        ia = 12               ! mass number
        iz =6                       ! charge number
        call sigrc(incp,emev,ia,iz,sigct,sigcr,sigcs)    ! Niita
        sige = sngl(sigcs)
        S(2)=Apb(2)*sige
      end if
c -----------------------------------------------------------------------
      S(3)=Apb(2)*SIGCINEL(En)
      S(4)=Apb(2)*SIGCNAL(En)
      S(5)=Apb(2)*SIGCNN3A(En)
c     include (n,3-He); 6/87.
      S(6)=Apb(2)*SIGCN3HE(En)
c     change 5/87: amalgamate (n,p) and (n,pn)
      S(7)=Apb(2)*(SIGCNPN(En)+SIGCNP(En))
      S(8)=Apb(2)*SIGCN2N(En)
      S(9)=Apb(2)*SIGCND(En)
      S(10)=Apb(2)*SIGCNT(En)
c
      Sigt=S(1)
      Do 1 J=2,10
      P(J+2)=1.0
    1 Sigt=Sigt+S(J)
c     Sigt = total "attenuation" cross section desired
      P(1)=2.0001
      P(2)=Sigt
      P(3)=S(1)/Sigt
      Do 4 J=2,10
      P(J+2)=P(J+1)+S(J)/Sigt
      IF (P(J+2) .LE. 0.999999) P(1)=P(1)+1.0
    4 Continue
c
c added by d.satoh (05.08.12) -------------------------------------------
cc      if( En .gt. 80.0 .and. En .lt. 150.0 ) then
c bug-fix by d.satoh (06.01.24)
c      if( En .gt. 80.0 .and. En .le. 150.0 ) then
c           incp = 2
c           emev = dble(En)
c           ia = 12
c           iz = 6
c           call sigrc(incp,emev,ia,iz,sigct,sigcr,sigcs)    ! Niita parametrization.
c          sigct2 = sngl(sigct) *0.95    ! correction factor for Carbon nucleus by daiki.
c           P(2) = S(1) + Apb(2)*sigct2
c           Totalx = P(2)
c           Return
c     end if
c -----------------------------------------------------------------------
      Totalx=Sigt
      Return
      END
c------------------------------------------------------------------------
c   This is file IBOX.FOR
c
      INTEGER FUNCTION IBOX(Path,Xsect,D4)
c     Purpose is to determine by random number if an interaction took place
c        along the PATH given an "attenuation" cross section XSECT.  If not,
c        set IBOX=0.  If so, return with IBOX = index (1-10) of the event
c        type and determine the distance along PATH to the interaction point.
c        Return with this distance as variable D4.
c
      Real*8 D,DP,U,X
      Common /RANDM/ Irx
      Common /PROB/ P(12), indx
c
      D4=0.0
      DP=dble(Path)
      X=dble(Xsect)
      Itype=0
      U=RAN(Irx)
      D=DEXP(-DP*X)
      IF (U .LT. D) goto 15
c
      Ip1=IFIX(P(1)) -1
      V=RAN(Irx)
      K=2
      Do 8 I=1,Ip1
      IF (V .LE. P(I+2)) goto 10
    8 K=K+1
   10 Itype=K-1
      D=-DLOG(U)/X
      D4=sngl(D)
   15 Ibox=Itype
      Return
      END
c------------------------------------------------------------------------
c    This is file TOTALX2.FOR
c
c     Purpose is to determine total attenuation cross section (in units of
c      cm**-1) and also the fractional parts of this total due to the separate
c      partial cross sections for incident neutron of energy Eneut.
c
      FUNCTION TOTALX2(Ene)
      Common /DENS/ Apb(2)
cs
      common /PROB2/p
      Dimension S(2)
csA      Dimension Apb(2)
*
c      apb(1)=0.04833
c      apb(2)=0.03984
*
      En=Ene
      S(1)=Sighyd2(En)*apb(1)
      S(2)=Sigctot(En)*apb(2)
c
      Sigt=S(1)+S(2)
      P=S(1)/Sigt
      Totalx2=Sigt
csd
c              write(13,*)'totalx2,p,apb(1),apb(2),s(1),s(2),en'
c     +                 ,totalx2,p,apb,s,en
      Return
      end
c------------------------------------------------------------------------
c   This is file IBOX2.FOR
c
      INTEGER FUNCTION IBOX2(Path,Xsect,D4)
c     Purpose is to determine by random number if an interaction took place
c        along the PATH given an "attenuation" cross section XSECT.  If not,
c        set IBOX=0.  If so, return with IBOX2 = index (1-2) of the event
c        type and determine the distance along PATH to the interaction point.
c        Return with this distance as variable D4.


c
      Real*8 D,DP,U,X
c
      common /RANDM/ Irx
cs
      common /PROB2/ P
c
      D4=0.0
cA00      DP=Path
      Dp=dble(Path)
cA00      X=Xsect
      X=dble(Xsect)
      Itype=0
      U=RAN(Irx)
ctest_plag-in
c      U=0.99
ctest
      D=DEXP(-DP*X)
cs               write(13,*)'D,U',D,U
      IF (U .LT. D) goto 15
c
      V=RAN(Irx)


      k=1


      IF (V .LE. P) goto 10
      k=2
   10 Itype=k
c                write(13,*)'x= ',x
      D=-DLOG(U)/X
cA00      D4=D
      D4=sngl(D)
   15 Ibox2=Itype
csd
c                write(13,*)'ibox2,d4,u,x,p,v',ibox2,d4,u,x,p,v
      Return
      END
c------------------------------------------------------------------------
c    This is file INTERACT.FOR
c
c        Purpose is to determine if an interaction took place along the
c         initial flight path of the neutron (from the source) in the
c         detector.  If not, NINDX = 0; if so, NINDX = 1, and the position
c         of the interaction and the neutron energy and type of interaction
c         will be stored in the TRIES array.
c
      Subroutine INTERACT(N,Nindx)
c               (The variable N is the running index of interactions)
      Common /RANDM/ Irx
      Common /TRY/ Tries (8,1000000)
      Common /MAXWEL/ Eu,Elow,T
      Common /NAID/ Rdet,Ht,Rcolim,R
      Common /PROB/ P(12), indx
      Dimension Enn(41),Sn(41)
      Data Indx/0/, Pi/3.1415927/
c
      Nindx=0
      M=N
      X=Tries(1,M)
      Y=Tries(2,M)
      Z=Tries(3,M)
      U=Tries(4,M)
      V=Tries(5,M)
      W=Tries(6,M)
      Flgt=PLNGTH(U,V,W, X,Y,Z, Rdet,Ht)
      Totx=P(2)
      IF (Indx .GT. 0) goto 6
c               If INDX is > 0 then using the same neutron energy as the last
c               time through this routine, so skip over intermediate steps.
      Eneut=Eu
      IF (T .EQ. 0.0) goto 5
c               Next steps are to choose the initial neutron energy from a
c               preselected distribution.  For T > 0.0 the energy is chosen
c               from a Maxwellian distribution.  For T < 0.0 the energy is
c               chosen from a uniform distribution.  However, there is a
c               slight variation in the latter in this particular coding;
c               a small skew toward smaller energies is included so as to
c               better represent the source flux of neutrons from the ORELA
c               for comparisons with response measurements taken with this
c               flux.
c
      IF (T .GT. 0.0) goto 2
      R1=RAN(Irx)
      R1=R1**1.2
c               The power 1.2 gives the skewing close enough.
      En=(Eu-Elow)*R1 + Elow
      Goto 4
c
c               Next is Maxwellian choice:
    2 En=CHOOSM(Elow,Eu,T)
c
    4 Eneut=En
    5 Totx=TOTALX(Eneut)
c
c       Now check for interactions along path -- use function IBOX
    6 Indx=1
      Itype=IBOX(Flgt,Totx,Dist)
      IF (Itype .LE. 0) Return
cs121201
c      if (Eneut .le. 0.0) return
cs121201
      Tries(1,M)=X+U*Dist
      Tries(2,M)=Y+V*Dist
      Tries(3,M)=Z+W*Dist
      Tries(7,M)=Eneut
      Tries(8,M)=FLOAT(Itype)+0.00001
      Nindx=1
      IF (T .NE. 0.0) Indx=0
      Return
      END
c------------------------------------------------------------------------
      Subroutine INTERACT2(N,Nindx)
c               (The variable N is the running index of interactions)


      Common /RANDM/ Irx
      Common /TRY/ Tries (8,1000000)
      Common /MAXWEL/ Eu,Elow,T
      Common /NAID/ Rdet,Ht,Rcolim,R
      Common /PROB/ PX(12), indx
      Dimension Enn(41),Sn(41)
      Data Indx/0/, Pi/3.1415927/
c---------------------------------------------
c
      Nindx=0
      M=N
      X=Tries(1,M)
      Y=Tries(2,M)
      Z=Tries(3,M)
      U=Tries(4,M)
      V=Tries(5,M)
      W=Tries(6,M)
      Flgt=PLNGTH(U,V,W, X,Y,Z, Rdet,Ht)
cs      Totx=P(2)
c_00/8/22      Totx=TOTALX2(Eu)
      Totx=TOTALX2(Eu)
c---------------------------------------------
c3_modified by satoh on 22.08.00
c3      ip=1
c3      totx = totalx3(ip,eu)
c---------------------------------------------
      IF (Indx .GT. 0) goto 6
c               If INDX is > 0 then using the same neutron energy as the last
c               time through this routine, so skip over intermediate steps.


      Eneut=Eu
      IF (T .EQ. 0.0) goto 5
c               Next steps are to choose the initial neutron energy from a
c               preselected distribution.  For T > 0.0 the energy is chosen
c               from a Maxwellian distribution.  For T < 0.0 the energy is
c               chosen from a uniform distribution.  However, there is a
c               slight variation in the latter in this particular coding;
c               a small skew toward smaller energies is included so as to
c               better represent the source flux of neutrons from the ORELA
c               for comparisons with response measurements taken with this
c               flux.
c
      IF (T .GT. 0.0) goto 2
      R1=RAN(Irx)
      R1=R1**1.2
c               The power 1.2 gives the skewing close enough.
      En=(Eu-Elow)*R1 + Elow
      Goto 4
c
c               Next is Maxwellian choice:


    2 En=CHOOSM(Elow,Eu,T)
c
    4 Eneut=En
cs    5 Totx=TOTALX(Eneut)
    5 Totx=totalx2(Eneut)
c3      ip=1
c3    5 Totx=totalx3(ip,Eneut)
c
c       Now check for interactions along path -- use function IBOX
    6 Indx=1
cs      Itype=IBOX(Flgt,Totx,Dist)
      Itype=IBOX2(Flgt,Totx,Dist)
      IF (Itype .LE. 0) Return
      Tries(1,M)=X+U*Dist
      Tries(2,M)=Y+V*Dist
      Tries(3,M)=Z+W*Dist
      Tries(7,M)=Eneut
      Tries(8,M)=FLOAT(Itype)+0.00001
      Nindx=1
      IF (T .NE. 0.0) Indx=0
      Return
      END
c------------------------------------------------------------------------
