c   This is file BANKER.FOR
c     Purpose is to do the bookkeeping, updating collision counters,
c       incrementing light-unit-bin counters, etcetera
c
c     Programming follows organization of the O5S Code.
c
c       Added (n,d) reaction channel 2/87.
c       Added (n,t) reaction channel 3/87.
c       Added photon reactions tally 3/87.
c       Added (n,2np) reaction tally 4/87.
c       Added (n,2p) and (n,2pn) reactions -- 5/87.
c       Added (n,dn) reaction -- 5/87.
c       Added decay of highly excited 10-B to give
c         multi-body breakup reaction -- 5/87.
c       Added decay of highly excited 8-Li to give
c         additional multi-body breakup -- 5/87.
c       Added (n,3-He) and (n,3-He n) reactions -- 6/87.
c       Attenuation of light by the scintillator added 6/87.
c       Z-interaction distribution included 7/87.
c
c       Added QMD reaction -- 8/05
c---------------------------------------------------------------------
c     Last modification at August/2005 by d.satoh
c---------------------------------------------------------------------
      Subroutine BANKR
c
      Integer SUMA,SUMB,SUMC,SUMCP,SUMD,SUME,SUMF,SUMG,SUMH,SUMJ,SUMK
      Integer SBOXT(999),SBNNPR(999),SB1H(999),SBCNLY(999),SB1A0H(999),
     A  SBND(999), SB3A0H(999), SBNP(999), SBNPN(999), SBN2N(999),
     B  Sphot(999),S3He(999),Alspec(1000),Dspec(1000),Pspec(1000),Sboxes
     C  ,S3HeX,Sumt,SBNT(999),Tspec(1000),Hspec(1000),Zsave(10)
      Integer CountLi,CountBe,CountB,maxl
c add by d.satoh at 05/08/09
c light-output bins for QMD
      integer SBPQ(999),SBDQ(999),SBTQ(999),SBHQ(999),SBAQ(999),SPIQ(999)
      real BINS, Bno
c
      Common /LTABLE/ ITAB,ENE(127),HYDL(127),CARL(127),ALPL(127),
     x Dlight(127)
      Common /NBOXES/ Wbox
      Common /NEUTRN/ Eneut, U, V, W, X, Y, Zint
      Common /NAID/   Rdet,Ht,Rjk,Rkn
      Common /COLLIS/ Nelm,Echrg(6)
      Common /NCOLLS/ Nclsns(21),Ndie,Nfstcl,Ntcoll
      Common /MAXWEL/ Esourc,Elow,Temp
      Common /GFLAG/  Igflag
      Common /NEUTR2/ En2, U2,V2,W2, X2,Y2,Z2
      Common /GAMATN/ Grho
c
csatoh 01/7/22
      common /effnaka/  bins, sboxt, maxl
csatoh 04/12/14
      common /INIT/ Nhist, Nhits
csatoh 05/02/01
      common /try/ tries(8,1000000)
      common /prob/ px(12), indx
      common /RES/ iswt, bias(5), icont, ictrl, iesc, ilight
csatoh 05/02/21
      common /Flightcom/ Flight
c
      Dimension BINS(1000),Bno(1000)
cc      Dimension LSAVE(30)
      Dimension LSAVE(50)
C
      Data Lboxes/999/, Sboxes/1000/
C     (Variable -LBOXES- is the number of light-unit bins.)
C     (Variable -Sboxes- is the number of energy bins.)


C  BANKR0 INITIALIZES CASE PARAMETERS.
      ENTRY BANKR0


C  Compute values of Light Bins


      call effn


      Bno(1)=1.0


      DO I = 1,Lboxes
         Bno(I+1)=Bno(I)+1.0
      end do


      Do 1020 I=1,10
 1020 Zsave(I)=0
C
C  INITIALIZE COUNTERS
C     NTCOL COUNTS THE TOTAL COLLISIONS.
C     NFSTCL COUNTS THE NEUTRONS HAVING AT LEAST ONE COLLISION.
C     NPRTES COUNTS THE PROTONS WHICH LEAK OUT OF THE SCINTILLATOR.
C     NDLOST counts leaking deuterons.
C     NDIE counts number of neutrons which have E < Ecutoff after one
C     or more scatterings (but are still in the detector).
C     NHYD COUNTS THE (N,H) REACTIONS.
C     NELAS COUNTS THE INELASTIC CARBON COLLISIONS.
C     NALP COUNTS THE (N, ALPHA) COLLISIONS.
C     N3ALP COUNTS THE (N, NPRIME, 3 ALPHA) COLLISIONS.
c     NONA1 counts (n,alpha p) Reactions -- see subroutine NN3ALF
c     NONA2 counts (n,alpha pn) Reactions -- also NN3ALF routine.
c     NONA3 counts (n,alpha p2n) Reactions -- also NN3ALF routine.
c     NONA4 counts (n,2alpha pnt) Reactions -- also NN3ALF routine.
c     NONA5 counts (n,2alpha p2nd) Reactions -- also NN3ALF routine.
C     NONP COUNTS THE (N,P) REACTIONS.
C     NOND Counts the (N,D) Reactions (added in after O5S).
C     NONDA count (n,d alpha 7-Li) events, Function ND
c     NONDN Counts the (N,DN) Reactions -- see Function ND
c     NONDN1 counts (n,dnp 9-Be) multibody reactions.
c     NONDN2 counts (n,dnd 2alpha) multibody reactions.
c     NONDN3 counts (n,dn alpha 6-Li) multibody reactions.
C     NONT Counts the (N,T) Reactions (added in after N,D)
c     NONT1 counts (n,tp 9-Be) multibody reactions.
c     NONT2 counts (n,td 2alpha) multibody reactions.
c     NONT3 counts (n,t alpha 6-Li) multibody reactions.
c     NONT4 counts (n,t pn 2alpha) reactions.
c     NONT5 counts (n,t pd 7-Li) reactions.
C     NONPN COUNTS THE (N,PN) REACTIONS.
c     NONP2N counts (n,p2n) reactions -- in Funtion NPX
c     NON2P counts the (n,2p) reactions -- see Function NPX
c     NON2PN counts the (n,2pn) reactions -- also Function NPX
c     NONP0 counts (n,p alpha) reactions -- also Function NPX
c     NONP1 counts (n,pn alpha) reactions -- also Function NPX
c     NONP2 counts (n,p2nd 2alpha) reactions -- also Function NPX
c     NONP3 counts (n,p2n alpha) reactions -- also Function NPX
c     NONP4 counts (n,2p2n) reactions -- also Function NPX
c     NONP5 counts (n,pnt 2alpha) reactions -- also Function NPX
C     NON2N COUNTS THE (N,2N) COLLISIONS.
c     NON2NP counts the (N,2NP) Collisions -- in Function N2N;
c     see also (n,p2n) reaction counter above.
c     NCLSNS counts the collision distribution.
c     NHONCE counts the number of neutrons having at least one
c     hydrogen collision.
c     NCONCE counts the number of neutrons having at least one carbon
c     collision.
c     NPHOT counts gamma rays detected.  Note: gamma rays MUST have
c       an accompanying (N,X) reaction, and so will not be included
c       in NTCOL.
c     NON3HE counts the (n,3-He) Collisions -- 6/87 addition.
c     NON3HEN counts the (n,3-He n) collisions -- also 6/87.
c     NONH1 counts (n,3-He 2n 2alpha) collisions;
c     NONH2 counts (n,3-He np 8-Li) collisions;
c     NONH3 counts (n,3-He 2np 7-Li) collisions; and
c     NONH4 counts (n,3-He 2npt alpha) collisions -- 7/87.
C
C
 1060 NTCOL  = 0
      NON3HE = 0
      NON3HEN= 0
      NONH1  = 0
      NONH2  = 0
      NONH3  = 0
      NONH4  = 0
      NFSTCL = 0
      NPRTES = 0
      Ndlost = 0
      Ndie   = 0
      NHYD   = 0
      NCARB  = 0
      NELAS  = 0
      NALP   = 0
      N3ALP = 0
      NONA1 = 0
      NONA2 = 0
      NONA3 = 0
      NONA4 = 0
      NONA5 = 0
      NONP  = 0
      NOND  = 0
      Nonda = 0
      NONDN = 0
      Nondn1 = 0
      Nondn2 = 0
      Nondn3 = 0
      NONT  = 0
      Nont1 = 0
      Nont2 = 0
      Nont3 = 0
      Nont4 = 0
      Nont5 = 0
      NONPN = 0
      NONP2N = 0
      NON2P  = 0
      NON2PN = 0
      NONP0 = 0
      NONP1 = 0
      NONP2 = 0
      NONP3 = 0
      NONP4 = 0
      NONP5 = 0
      NON2N = 0
      NON2NP = 0
      NHONCE = 0
      NCONCE = 0
      Nphot = 0
      DO 1061 L=1,21
1061  NCLSNS(L)=0
C INITIALIZE REACTION LIGHT BINS.
      DO  1070  L = 1, LBOXES
      SBOXT(L) = 0
      SB1H(L) = 0
      SBCNLY(L) = 0
      SBNNPR(L) = 0
      SBND(L) = 0
      SBNT(L)  = 0
      SB1A0H(L) = 0
      SB3A0H(L) = 0
      SBNP(L) = 0
      SBNPN(L) = 0
      SBN2N(L) = 0
      Sphot(L) = 0
      S3He(L)=0
c add 05/Aug
      SBPQ(L) = 0
      SBDQ(L) = 0
      SBTQ(L) = 0
      SBHQ(L) = 0
      SBAQ(L) = 0
      SPIQ(L) = 0
 1070 CONTINUE
c     Added in determinations of alpha, deuteron, and proton spectra.
c       First initialize spectra arrays -- bins will be 1 MeV in
c       response energy.
c     Added in triton spectrum, too.
c     And, finally, 3-He spectrum.
c     Keep count of 'heavy' ions -- Li, Be, and B
      Do 1071 L=1,Sboxes
      Dspec(L)=0
      Tspec(L)=0
      Pspec(L)=0
      Hspec(L)=0
 1071 Alspec(L)=0
      CountLi=0
      CountBe=0
      CountB=0
csatoh 01/07/22
      maxl = 0
csatoh 05/02/01
      do m = 1, 6
            Echrg(m) = 0.0
      end do
c
      do mm = 1,12
            px(mm) = 0.0
      end do
c
      indx = 0
c
      RETURN
C
C
C
C
C                    INDIVIDUAL HISTORY INITIALIZATION
      ENTRY BANKR1
c
C       INITIALIZE THE LIGHT SUMMATION PARAMETER FOR THIS NEUTRON.
      FLIGHT=0.
C       INITIALIZE THE COLLISION INDICATION PARAMETERS.
      NTEST=0
      NCTST=0
C       INITIALIZE COLLISION COUNTERS
c          These will be used in the "P.S.D." analysis for NE-213, if
c          requested.
c Tatsu 04/6/28
      FLtotal=0.0 ! Total Light output before considering quenching effect
      FLstop=0.0  ! Total Light Output by Stopped Particle & Residual Nuclide
      Ethre=0.0   ! Escaped proton above this energy identified as gamma
c
      DO  2010  I = 1, 50
 2010 Lsave(I) = 0
      Npsd=1
      Nprej=Nprtes
      Ndrej=Ndlost
      RETURN
C
C
C
C
C                    COLLISION ANALYSIS
C  BANKR2 ANALYZES THE LIGHT GIVEN IN THE FOLLOWING REACTIONS
C    1. NELM=1, N+H-->N+H, PROTON LEAKAGE CONSIDERED.
C    2. NELM=2, N+12C-->N+12C, ELASTIC.
C    3. NELM=3, N+12C-->N+12C, INELASTIC.
C    4. NELM=4, N+12C-->ALPHA+9BE(ground state only)
C    5. NELM=5, N+12C-->N+3*ALPHA
C    6. NELM=6, N+12C-->P+12B, PROTON LEAKAGE CONSIDERED.
C    7. NELM=7, N+12C-->P+N+11B, PROTON LEAKAGE CONSIDERED; also
c                    possible N+12C-->P+P+11Be and
c                    possible N+12C-->P+P+N+10Be reactions.
C    8. NELM=8, N+12C-->N+N+11C, with possible 11-C --> p + 10-B.
C    9. NELM=9, N+12C-->D+11B  (Added in after O5S)
C                    with subsequent 11B-->N+10B possibility.
C   10. NELM=10.  Photon interaction in the detector (added after O5S).
C   11. NELM=11.  N+12C-->T+10B  (Added in after O5S)
C   12. NELM=12.  N+12C --> 3He+10-Be  (Added in 6/87; more in 7/87).
C   13. NELM=-1.  N+12C --> pions      (Added in 11/99; by d.satoh).
C   14. NELM=-2.  N+12C --> alphas     (Added in 11/99; by d.satoh).
C   15. NELM=-3.  N+12C --> deuterons  (Added in 11/99; by d.satoh).
C   16. NELM=-4.  N+12C --> tritons    (Added in 01/00; by d.satoh).
C   17. NELM=-5.  N+12C --> 3He        (Added in 01/00; by d.satoh).
C   18. NELM=-6.  proton production via QMD  (Added in 02/05; by d.satoh).
C
      ENTRY BANKR2
C  INITIALIZE LIGHT SUMMER.
      AMTLT = 0.
c
C  BRANCH ON COLLISION TYPE.
      GO TO (10, 20, 20, 40, 50, 10, 10, 20, 60, 70, 75, 80), NELM
      if (NELM.eq.-1) goto 666
      if (NELM.eq.-2) goto 777
      if (NELM.eq.-3) goto 888
      if (NELM.eq.-4) goto 999
      if (NELM.eq.-5) goto 555
      if (NELM.eq.-6) goto 444
      if (NELM.eq.-10) return ! neutron
c
      write(*,*) 'Error: NELM=',NELM
      return
c
c ----------------------------------------------
c     PROTON
c ----------------------------------------------
c
  444 En=Echrg(1)
c     Get proton light, then check for possible proton escape.
      N=2
      Amtlt=ETOL(N,En)
      if(iesc.eq.0.and.Echrg(2).le.Ethre) FLstop=FLstop+ETOL(N,En) ! Tatsu 04/6/28
      Ep=Echrg(2)
      IF (Ep .LE. 0.0 ) goto 121
c if program encounts this section, proton leack from scintillator.
      Amtlt=Amtlt-ETOL(N,Ep)
      if(iesc.eq.0.and.Echrg(2).le.Ethre) FLstop=FLstop-ETOL(N,Ep) ! Tatsu 04/6/28
  121 continue
      L=1+IFIX(Echrg(1))
      IF (L.GT.Sboxes) L=Sboxes
      Pspec(L)=Pspec(L)+1
      if(iesc.eq.2 .and. Echrg(2).gt.0.0) then ! kajimoto 2011/03/18
         Amtlt=0.0
         goto 3400
      end if ! kajimoto 2011/03/18
c
      if( Amtlt .lt. 0.0 ) then
            write(*,*) 'Error: Amtlt is negative!!', NELM, Echrg(1), Echrg(2), Amtlt
            Amtlt=0.0
      end if
      Goto 3400
c
c ----------------------------------------------
c     PION
c ----------------------------------------------
c
  666 Epi=Echrg(1)
c       Get pion light, then check for possible pion escape.
      if(iesc.eq.2 .and. Echrg(2).gt.0.0) then ! kajimoto 2011/03/18
         Amtlt=0.0
         goto 3400
      end if ! kajimoto 2011/03/18
      Amtlt=-4.9008+0.82915*Epi
      if(iesc.eq.0.and.Echrg(2).le.Ethre) FLstop=FLstop+Amtlt  ! Tatsu 04/6/28
      Epif=Echrg(2)
      IF (Epif .LE. 0.0) goto 111
      Amtlt=Amtlt+4.9008-0.82915*Epp
      if(iesc.eq.0.and.Echrg(2).le.Ethre) FLstop=FLstop
     &+4.9008-0.82915*Epp  ! Tatsu 04/6/28
  111 continue
c_modified by satoh on 23/08/2000
      if( Amtlt .lt. 0.0 ) then
            write(*,*) 'Error: Amtlt is negative!!', NELM, Echrg(1), Amtlt
            Amtlt=0.0
      end if
      Goto 3400
c
c ----------------------------------------------
c     ALPHA   (made by d.satoh '99.11.11)
c ----------------------------------------------
c
  777 Eal=Echrg(1)
c       Get alpha light, then check for possible alpha escape.
c               Increment alpha spectrum array.
      L=1+IFIX(Eal)
      IF (L .GT. Sboxes) L=Sboxes
      Alspec(L)=Alspec(L)+1
c               Get alpha light --
      N=6
      if(iesc.eq.2 .and. Echrg(2).gt.0.0) then ! kajimoto 2011/03/18
         Amtlt=0.0
         goto 3400
      end if ! kajimoto 2011/03/18
      Amtlt=ETOL(N,Eal)
      if(iesc.eq.0.and.Echrg(2).le.Ethre) then  ! Tatsu 04/6/28
      FLstop=FLstop+ETOL(N,Eal)
      endif
c ----
      Ealf=Echrg(2)
      IF (Ealf .LE. 0.0) goto 3400
      Amtlt=Amtlt-ETOL(N,Ealf)
      if(iesc.eq.0.and.Echrg(2).le.Ethre) FLstop=FLstop-ETOL(N,Ealf) ! Tatsu 04/6/28
c ----
c_modified by satoh on 23/08/2000
      if( Amtlt .lt. 0.0 ) then
            write(*,*) 'Error: Amtlt is negative!!', NELM, Echrg(1), Amtlt
            Amtlt=0.0
      end if
      Goto 3400
c
c ----------------------------------------------
c     DEUTERON   (made by d.satoh '99.11.11)
c ----------------------------------------------
c
  888 Ed=Echrg(1)
c               First increment deuteron spectrum array --
      L=1+IFIX(Ed)
      IF (L .GT. Sboxes) L=Sboxes
      Dspec(L)=Dspec(L)+1
c               Next get deuteron light --
      N=8
      if(iesc.eq.2 .and. Echrg(2).gt.0.0) then ! kajimoto 2011/03/18
         Amtlt=0.0
         goto 3400
      end if ! kajimoto 2011/03/18
      Amtlt=ETOL(N,Ed)
      if(iesc.eq.0.and.Echrg(2).le.Ethre) then ! Tatsu 04/6/28
       FLstop=FLstop+ETOL(N,Ed)
      endif
c ----
      Edf=Echrg(2)
      if (Edf .le. 0.0) goto 3400
      Amtlt=Amtlt-ETOL(N,Edf)
      if(iesc.eq.0.and.Echrg(2).le.Ethre) FLstop=FLstop-ETOL(N,Edf) ! Tatsu 04/6/28
c ----
c_modified by satoh on 23/08/2000
      if( Amtlt .lt. 0.0 ) then
            write(*,*) 'Error: Amtlt is negative!!', NELM, Echrg(1), Amtlt
            Amtlt=0.0
      end if
      goto 3400
c
c ----------------------------------------------
c     TRITON   (made by d.satoh '99.11.11)
c ----------------------------------------------
c
  999 Etr=Echrg(1)
c               First increment triton spectrum array --
      L=1+IFIX(Etr)
      IF (L .GT. Sboxes) L=Sboxes
      Tspec(L)=Tspec(L)+1
c               Next get triton light --
      N=8
      if(iesc.eq.2 .and. Echrg(2).gt.0.0) then ! kajimoto 2011/03/18
         Amtlt=0.0
         goto 3400
      end if ! kajimoto 2011/03/18
      Amtlt=0.8*ETOL(N,Etr)
      if(iesc.eq.0) FLstop=FLstop+0.8*ETOL(N,Etr)  ! Tatsu 04/6/28
c
      Etf=Echrg(2)
      if (Etf .le. 0.0) goto 3400
c_modified by satoh on 23/08/2000
      if( Amtlt .lt. 0.0 ) then
            write(*,*) 'Error: Amtlt is negative!!', NELM, Echrg(1), Amtlt
            Amtlt=0.0
      end if
      goto 3400
c
c ----------------------------------------------
c     3He   (made by d.satoh '99.11.11)
c ----------------------------------------------
c
  555 EHe=Echrg(1)
c               First increment 3He spectrum array --
      L=1+IFIX(EHe)
      IF (L .GT. Sboxes) L=Sboxes
      Hspec(L)=Hspec(L)+1
c               Next get 3He light --
      N=6
      if(iesc.eq.2 .and. Echrg(2).gt.0.0) then ! kajimoto 2011/03/18
         Amtlt=0.0
         goto 3400
      end if ! kajimoto 2011/03/18
      Amtlt=1.25*ETOL(N,EHe)
      if(iesc.eq.0) FLstop=FLstop+1.25*ETOL(N,EHe)  ! Tatsu 04/6/28
c
      EHf=Echrg(2)
      if (EHf .le. 0.0) goto 3400
c_modified by satoh on 23/08/2000
      if( Amtlt .lt. 0.0 ) then
            write(*,*) 'Error: Amtlt is negative!!', NELM, Echrg(1), Amtlt
            Amtlt=0.0
      end if
      goto 3400
C
c======================================================================c
C  NELM=1.  IN ORDER TO APPROPRIATELY CALCULATE THE AMOUNT OF LIGHT
C       GIVEN INSIDE THE SCINTILLATOR, FIRST FIND THE AMOUNT
C       CORRESPONDING TO THE ENERGY ORIGINALLY GIVEN THE PROTON, THEN
C       FIND THE AMOUNT CORRESPONDING TO THE ENERGY WITH WHICH
C       THE PROTON ESCAPED THE SCINTILLATOR.  THE ACTUAL AMOUNT IS THE
C       DIFFERENCE OF THESE QUANTITIES.    (NOTE---  THIS SAME
C       PROCEDURE IS USED FOR NELM = 6 AND NELM = 7.)
C
   10 En=Echrg(1)
c     Get proton light, then check for possible proton escape.
      N=2
      Amtlt=ETOL(N,En)
      if(iesc.eq.0.and.Echrg(2).le.Ethre) FLstop=FLstop+ETOL(N,En) ! Tatsu 04/6/28
      Ep=Echrg(2)
      IF (Ep .LE. 0.0 ) goto 12
c if program encounts this section, proton leack from scintillator.
      Nprtes=Nprtes+1
      Amtlt=Amtlt-ETOL(N,Ep)
      if(iesc.eq.0.and.Echrg(2).le.Ethre) FLstop=FLstop-ETOL(N,Ep) ! Tatsu 04/6/28
   12 if(iesc.eq.2.and.Echrg(2).gt.0.0) Amtlt=0.0 ! kajimoto 2011/03/18
      IF (Nelm .GE. 6) goto 15
      NHYD = NHYD + 1
c     Keep tabs on interaction spot.
check_A  (Added in 01/00 by d.satoh)
c      Increament proton spectrum -- n+p->n+p reactions (SCINFUL)
c      Increament proton spectrum -- every p emission reactions.(SCINFUL-QMD)
      L=1+IFIX(Echrg(1))
      IF (L.GT.Sboxes) L=Sboxes
      Pspec(L)=Pspec(L)+1
      Depth = 1.0 + 10.0*Zint/Ht
      Idepth=IFIX(Depth)
      IF (Idepth .GT. 10) Idepth=10
      Zsave(Idepth)=Zsave(Idepth) + 1
      IF (Ntest .GT. 0) goto 3400
      Ntest=1
      Nhonce=Nhonce+1
      Goto 3400
c
c        Increment proton spectrum -- only protons from 12-C reactions.
   15 L=1+IFIX(Echrg(1))
      IF (L.GT.Sboxes) L=Sboxes
      Pspec(L)=Pspec(L)+1
c         Now finish up with carbon light ---
c     Test for tertiary reaction channel active
      IF (Nelm .EQ. 7) goto 16
c     It's not, so finish up 12-B particle-stable system
      N=4
      En=Echrg(3)
      IF (En .GT. 0.0) Amtlt=Amtlt + 1.2*ETOL(N,En)
      if(iesc.eq.0) FLstop=FLstop+1.2*ETOL(N,En)  ! Tatsu 04/6/28
      if(iesc.eq.2.and.Echrg(2).gt.0.07) Amtlt=0.0 !kajimoto 04/14/2011
      CountB = CountB + 1
c     Boron-ion light somewhat larger than Carbon-ion light (see Nucl.
c      Instrum. Methods 138 (1976) 93).
      Nonp=Nonp+1
      Goto 3400
c
c     NELM=7, tertiary reaction: first check for triton reaction
   16 IF (Echrg(3) .LE. 0.0) goto 63
c     Not that; check for alpha emission
      IF (Echrg(4) .GT. 0.0) goto 74
c     Not that; try for second proton emission
      IF (Echrg(5) .GT. 0.0) goto 67
c     No second proton; take care of residual boron-ion energy.
      En=Echrg(3)
      N=4
      Amtlt=Amtlt+1.2*ETOL(N,En)
      if(iesc.eq.0) FLstop=FLstop+1.2*ETOL(N,En)  ! Tatsu 04/6/28
      CountB = CountB + 1
      IF (En2 .GE. 0.0) Nonpn=Nonpn+1
      IF (En2 .LT. 0.0) Nonp2n=Nonp2n+1
      Goto 3400
   11 IF (Eneut .GT. 0.0) goto 13
      NON2P=NON2P + 1
      Goto 3400
   13 IF (En2 .GE. 0.0) Non2pn=Non2pn+1
      IF (En2 .LT. 0.0) Nonp4 = Nonp4 + 1
      Goto 3400
c       Next is for n + 12-C --> p + alpha + x-Li + 0,1,2 neutrons
   14 IF (Eneut .GT. 0.0) goto 19
      Nonp0=Nonp0+1
      Goto 3400
   19 IF (En2 .GE. 0.0) Nonp1=Nonp1+1
      IF (En2 .LT. 0.0) Nonp3=Nonp3+1
      Goto 3400
c       Next is for n + 12-C --> p + alpha + n + triton + alpha
   18 Nonp5=Nonp5+1
      En=Echrg(6)
      L=1+IFIX(En)
      IF (L .GT. Sboxes) L=Sboxes
      Tspec(L)=Tspec(L)+1
      N=8
      Amtlt=Amtlt+0.8*ETOL(N,En)
      if(iesc.eq.0) FLstop=FLstop+0.8*ETOL(N,En)  ! Tatsu 04/6/28
      Goto 3400
c
C  NELM=2, 3, or 8.  GET LIGHT GIVEN BY CARBON.
   20 En=Echrg(1)
      N=4
      Amtlt=ETOL(N,En)
      if(iesc.eq.0) FLstop=FLstop+ETOL(N,En)  ! Tatsu 04/6/28
      IF (Nelm-3) 25,30,35
   25 NCARB=NCARB+1
      IF (Nctst .GT. 0) goto 3400
      Nctst=1
      Nconce=Nconce+1
      Goto 3400
   30 NELAS=NELAS+1
      Goto 3400
c   Next is for (N,2N) with added feature of subsequent 11-C decay
   35 Enn=Echrg(2)
      IF (Enn .GT. 0.0) goto 36
      Non2n=Non2n+1
      Goto 3400
   36 Amtlt=Amtlt +0.2*ETOL(N,En)
      if(iesc.eq.0) FLstop=FLstop+0.2*ETOL(N,En)  ! Tatsu 04/6/28
c     That increments from 'Carbon' to 'Boron' light
      CountB=CountB+1
      N=2
      if(iesc.ne.2.or.Echrg(3).eq.0.0)  Amtlt=Amtlt+ETOL(N,Enn) ! kajimoto 2011/03/18
      if(iesc.eq.0.and.Echrg(3).le.Ethre) FLstop=FLstop+ETOL(N,Enn)  ! Add by Tatsu 04/06/28
c     Increment proton spectrum
      L=1+IFIX(En)
      IF (L .GT. Sboxes) L=Sboxes
      Pspec(L)=Pspec(L)+1
      Non2np=Non2np+1
c     Now check for proton escape.
      En=Echrg(3)
      IF (En .LE. 0.0) goto 3400
      Nprtes=Nprtes+1
      Amtlt=Amtlt-ETOL(N,En)
      if(iesc.eq.0.and.Echrg(3).le.Ethre) FLstop=FLstop-ETOL(N,En)  ! Add by Tatsu 04/06/28
      Goto 3400
c
C  NELM=4.  INCREMENT (N,ALPHA) COUNTER.
   40 NALP = NALP + 1
      En=Echrg(1)
c           Increment alpha spectrum array.
      L=1+IFIX(En)
      IF (L .GT. Sboxes) L=Sboxes
      Alspec(L)=Alspec(L)+1
c           Get alpha light --
      if(iesc.eq.2 .and. En.eq.0.0) goto 3400 ! kajimoto 2011/04/14
      N=6
      Amtlt=ETOL(N,En)
      Adaiki=Amtlt
      if(iesc.eq.0) FLstop=FLstop+ETOL(N,En)  ! Tatsu 04/6/28
c     AMTLT is the alpha light; now get the Be light.
   41 En=Echrg(2)
      N=4
      Amtlt=Amtlt+1.33*ETOL(N,En)
      Asatoh=Amtlt-Adaiki
      if(iesc.eq.0) FLstop=FLstop+1.33*ETOL(N,En)  ! Tatsu 04/6/28
      CountBe=CountBe+1
c     Be light is about 1.33 x C light (G. Dietze and H. Klein,
c       PTB-Bericht ND-22 (1982)).
      Goto 3400
c
c  NELM=5.  INCREMENT (N, 3 ALPHA) COUNTER.
   50 IF (Echrg(4) .GT. 0.0) goto 53
      N3ALP = N3ALP + 1
      Amtlt=0.0
      I1=1
      I2=3
c       The alpha spectrum will have 3 alphas for each event.
   51 N=6
      Do 52 I=I1,I2
      En=Echrg(I)
      L=1+IFIX(En)
      IF (L .GT. Sboxes) L=Sboxes
      Alspec(L)=Alspec(L)+1
c       However, the alpha light is summed -- 1 count for each event.
      Amtlt=Amtlt+ETOL(N,En)
   52 if(iesc.eq.0) FLstop=FLstop+ETOL(N,En)  ! Tatsu 04/6/28
      Goto 3400
c
c     Account for possible 8-Li breakup reactions.
   53 En=Echrg(4)
      N=2
      L=1+IFIX(En)
      IF (L .GT. Sboxes) L=Sboxes
      Pspec(L)=Pspec(L)+1
      Amtlt=ETOL(N,En)
      if(iesc.eq.0.and.Echrg(2).le.Ethre) FLstop=FLstop+ETOL(N,En)  ! Add by Tatsu 04/06/28
      En=Echrg(2)
c     Test for proton escape.
      IF (En .GT. 0.0) Amtlt=Amtlt-ETOL(N,En)
      if(iesc.eq.0.and.Echrg(2).le.Ethre) FLstop=FLstop-ETOL(N,En) ! Add by Tatsu 04/06/28
      IF (En .GT. 0.0) Nprtes=Nprtes+1
c     Check for triton decay mode.
      En = Echrg(6)
      IF (En .LE. 0.0) goto 55
c     Program counter to here, had triton + (second) alpha.
      N=8
      Amtlt=Amtlt+0.8*ETOL(N,En)
      if(iesc.eq.0) FLstop=FLstop+0.8*ETOL(N,En)  ! Tatsu 04/6/28
      Nona4=Nona4+1
      L=1+IFIX(En)
      IF (L .GT. Sboxes) L=Sboxes
      Tspec(L)=Tspec(L)+1
c     (Shift Echrg array for convenience for final two alphas.)
   54 Exx=Echrg(1)
      Echrg(1)=Echrg(4)
      Echrg(4)=Exx
      I1=4
      I2=5
      Goto 51
c     Next programming--not triton, so test for deuteron mode of
c       decay of 8-Li.
   55 IF (Echrg(5) .LE. 0.0) goto 56
c     Program counter to here, deuteron mode encountered.  Get deuteron
c       information counted.
      Nona5=Nona5+1
      En=Echrg(3)
      L=1+IFIX(En)
      IF (L .GT. Sboxes) L=Sboxes
      Dspec(L)=Dspec(L)+1
      N=8
      Amtlt=Amtlt+ETOL(N,En)
      if(iesc.eq.0) FLstop=FLstop+ETOL(N,En)  ! Tatsu 04/6/28
c     Now finish up the two alphas:
      Goto 54
c     Program counter to here means decay to 8-Li, 7-Li or 6-Li.
   56 N=4
      En=Echrg(3)
      C2Li=2.6-0.006*En
      Amtlt=Amtlt+C2Li*ETOL(N,En)
      if(iesc.eq.0) FLstop=FLstop+C2Li*ETOL(N,En)  ! Tatsu 04/6/28
      CountLi=CountLi+1
c     Variable C2Li represents very approximate relationship between
c       Li-ion light and carbon-ion light; see Nucl. Instrum. Methods
c       138 (1976) 93 for measurements on NE-102.
      I1=1
      I2=1
      IF (Eneut .GT. 0.0) goto 58
      Nona1=Nona1+1
      Goto 51
   58 IF (En2 .LT. 0.0) Nona3=Nona3+1
      IF (En2 .GE. 0.0) Nona2=Nona2+1
      Goto 51
c     END 8-Li decay.
c
c  NELM=9.  Increment (N,D) reaction
   60 En=Echrg(1)
c           First increment deuteron spectrum array --
      L=1+IFIX(En)
      IF (L .GT. Sboxes) L=Sboxes
      Dspec(L)=Dspec(L)+1
c           Next get deuteron light --
      N=8
      Amtlt=ETOL(N,En)
      if(iesc.eq.0.and.Echrg(2).le.Ethre) FLstop=FLstop+Amtlt  ! Tatsu 04/6/28
 9888 format(/'  Amount of D-light='1pe10.3/'    Echrg-array ='
     a        1p6e10.3/'   Eneut ='1pe11.4/)
c       See Nuclear Instrum. & Methods 138 (1976) 93 for
c       measurements of deuteron light in NE-102.
      En = Echrg(2)
c       Test for escaped deuteron just like proton case, above.
      IF (En .LE. 0.0) goto 62
      Ndlost=Ndlost+1
      Amtlt=Amtlt-ETOL(N,En)
      if(iesc.eq.0.and.Echrg(2).le.Ethre) FLstop=FLstop-ETOL(N,En)  ! Tatsu 04/6/28
   62 if(iesc.eq.2.and.Echrg(2).gt.0.0) Amtlt=0.0 ! kajimoto 2011/03/18
      IF (Eneut .GT. 0.0) goto 61
      IF (Echrg(4) .LE. 0.0) goto 59
c     Program counter to here had (n,d alpha 7-Li) event.
      Nonda=Nonda+1
      En = Echrg(3)
      N=4
      C2Li=2.6-0.006*En
      Amtlt=Amtlt+C2Li*ETOL(N,En)
      if(iesc.eq.0) FLstop=FLstop+C2Li*ETOL(N,En)  ! Tatsu 04/6/28
      CountLi=CountLi+1
c     See comment after statement 56 for Li light.
      I1=3
      I2=3
      Goto 51
   59 Nond=Nond+1
      Goto 68
c
c     Check for multibody breakup reactions following (n,dn)
c     See Subroutine TENBDK for "signatures" of different reactions.
   61 IF (Echrg(4) .LE. 0.0) goto 66
      IF (Echrg(5) .LE. 0.0) goto 65
c       Next portion for 10-B --> d + 2 alpha breakup.
   63 N=6
      Do 64 I=4,5
      En=Echrg(I)
      L=1+IFIX(En)
      IF (L .GT. Sboxes) L=Sboxes
      Alspec(L)=Alspec(L)+1
   64 Amtlt=Amtlt+ETOL(N,En)
      if(iesc.eq.0) FLstop=FLstop+ETOL(N,En)  ! Tatsu 04/6/28
      IF (Nelm .EQ. 7  .AND.  Echrg(3) .LE. 0.0) goto 18
c     That gets the contributions from the 2 alphas; now get deuteron
c       (or triton if Nelm=7 and Echrg(3) = 0) contribution.
      N=8
      En=Echrg(3)
      L=1+IFIX(En)
      IF (L .GT. Sboxes) L=Sboxes
      Dspec(L)=Dspec(L)+1
      Amtlt=Amtlt+ETOL(N,En)
      if(iesc.eq.0) FLstop=FLstop+ETOL(N,En)  ! Tatsu 04/6/28
      IF (Nelm .EQ. 7) Nonp2=Nonp2+1
      IF (Nelm .EQ. 9) Nondn2=Nondn2+1
      IF (Nelm .EQ. 11) Nont2=Nont2+1
      Goto 3400
c
c     10-B --> alpha + 6-Li is next.
   65 N=6
      En=Echrg(4)
      L=1+IFIX(En)
      IF (L .GT. Sboxes) L=Sboxes
      Alspec(L)=Alspec(L)+1
      Amtlt=Amtlt+ETOL(N,En)
      if(iesc.eq.0) FLstop=FLstop+ETOL(N,En)  ! Tatsu 04/6/28
      En=Echrg(3)
      N=4
      C2Li=2.6-0.006*En
      Amtlt=Amtlt+C2Li*ETOL(N,En)
      if(iesc.eq.0) FLstop=FLstop+C2Li*ETOL(N,En)  ! Tatsu 04/6/28
c     See comment above after statement 56 about Li-ion light.
      CountLi=CountLi+1
      IF (Nelm .EQ. 7) goto 14
      IF (Nelm .EQ. 9) Nondn3=Nondn3+1
      IF (Nelm .EQ. 11) Nont3=Nont3+1
      Goto 3400
c
   66 IF (Echrg(5) .LE. 0.0) goto 69
c     10-B --> p + 9-Be is next bit of programming.
   67 N=2
      En=Echrg(5)
      L=1+IFIX(En)
      IF (L .GT. Sboxes) L=Sboxes
      Pspec(L)=Pspec(L)+1
      if(iesc.eq.2.and.Echrg(6).gt.0.0) goto 71 ! kajimoto 2011/03/18
      Amtlt=Amtlt+ETOL(N,En)
      if(iesc.eq.0.and.Echrg(6).le.Ethre) FLstop=FLstop+ETOL(N,En)  ! Tatsu 04/6/28
      En=Echrg(6)
      IF (En .LE. 0.0) goto 71
      Nprtes=Nprtes+1
      Amtlt=Amtlt-ETOL(N,En)
      if(iesc.eq.0.and.Echrg(6).le.Ethre) FLstop=FLstop-ETOL(N,En) ! Tatsu 04/6/28
c     Check for  d + 7-Li  reaction following (n,2p2n)
      IF (Nelm.EQ.7 .AND. Echrg(4).LT. 0.0) goto 79
c     Not that, so get Be-ion energy and light.
   71 N=4
      En=Echrg(3)
      Amtlt=Amtlt+1.33*ETOL(N,En)
      if(iesc.eq.0) FLstop=FLstop+1.33*ETOL(N,En)  ! Tatsu 04/6/28
      CountBe=CountBe+1
c     See comment for NELM=8 regarding Be light.
      IF (Nelm .EQ. 7) goto 11
      IF (Nelm .EQ. 9) Nondn1=Nondn1+1
      IF (Nelm .EQ. 11) Nont1=Nont1+1
      Goto 3400
c
c     10-B didn't further breakup for next programming.
   69 Nondn=Nondn+1
   68 En=Echrg(3)
      N=4
      Amtlt=Amtlt+1.2*ETOL(N,En)
      if(iesc.eq.0) FLstop=FLstop+1.2*ETOL(N,En)  ! Tatsu 04/6/28
      CountB=CountB+1
c     Assume Boron-ion light = 1.2 times Carbon-ion light.  It's an
c        assumption which won't matter much unless vastly wrong.
      Goto 3400
c
c
c  NELM=10.  Photon interaction recorded in the detector.  The amount
c     of light was determined in the analysis in Subroutine PHOTON
c     and saved in variable Echrg(1).
   70 Nphot=Nphot+1
c      Now check to see if rejecting results when a gamma ray is detected.
      IF (Igflag .EQ. 1) goto 4400
      if(iesc.ne.0.and.igflag.eq.1) goto 4400 ! Tatsu 04/6/28
      Amtlt=Echrg(1)
      Goto 3400
c
c  NELM=11.  Triton reaction in detector.  Increment (N,T) counter, etc.
   75 En=Echrg(1)
      L=1+IFIX(En)
      IF (L .GT. Sboxes) L=Sboxes
      Tspec(L)=Tspec(L)+1
      N=8
      Amtlt=0.8*ETOL(N,En)
      if(iesc.eq.2 .and. En.eq.0.0)  Amtlt=0.0 ! kajimoto 2011/03/18
      if(iesc.eq.0) FLstop=FLstop+0.8*ETOL(N,En)  ! Tatsu 04/6/28
c     Triton light assumed less than deuteron light
c  Check for further 10-B ion breakup.
      IF (Echrg(4) .EQ. 0.0) goto 76
      IF (Echrg(4) .LT. 0.0) goto 79
      IF (Eneut .GT. 0.0) goto 77
   74 IF (Echrg(5) .LE. 0.0) goto 65
      Goto 63
   76 IF (Echrg(5) .GT. 0.0) goto 67
c     If program counter gets to here, 10-B didn't break up.
      Nont=Nont+1
      En=Echrg(3)
      N=4
      Amtlt=Amtlt+1.2*ETOL(N,En)
      if(iesc.eq.0) FLstop=FLstop+1.2*ETOL(N,En)  ! Tatsu 04/6/28
      CountB=CountB+1
      Goto 3400
c     Program counter to here, have 10-B --> n + p + 2 alpha breakup.
   77 Nont4 = Nont4 + 1
      N=2
      En=Echrg(5)
      L=1 + IFIX(En)
      IF (L .GT. Sboxes) L=Sboxes
      Pspec(L)=Pspec(L)+1
      if(iesc.eq.2 .and. Echrg(6).gt.0.0) goto 78 ! kajimoto 2011/03.18
      Amtlt=Amtlt+ETOL(N,En)
      if(iesc.eq.0.and.Echrg(6).le.Ethre) FLstop=FLstop+ETOL(N,En)  ! Tatsu 04/6/28
      En=Echrg(6)
c     check for proton escape.
      IF (En .LE. 0.0) goto 78
      Nprtes=Nprtes+1
      Amtlt=Amtlt-ETOL(N,En)
      if(iesc.eq.0.and.Echrg(6).le.Ethre) FLstop=FLstop-ETOL(N,En)  ! Tatsu 04/6/28
c     Now finish up the 2 alphas:
   78 I1=3
      I2=4
      Goto 51
c
c   Next is for  d + 7-Li  breakup of 9-Be.
   79 N=8
      En=-Echrg(4)
c     (Recall that in subroutine TENBDK, E(deuteron) was set negative.)
      L=1+IFIX(En)
      IF (L .GT. Sboxes) L=Sboxes
      Dspec(L) = Dspec(l) + 1
      Amtlt=Amtlt+ETOL(N,En)
      if(iesc.eq.0) FLstop=FLstop+ETOL(N,En)  ! Tatsu 04/6/28
      En=Echrg(3)
      C2Li=2.6-0.006*En
c     (For Li-ion light output see comment after statement 56, above.)
      N=4
      Amtlt=Amtlt + C2Li*ETOL(N,En)
      if(iesc.eq.0) FLstop=FLstop+C2Li*ETOL(N,En)  ! Tatsu 04/6/28
      CountLi=CountLi+1
      IF (Nelm .EQ. 11) Nont5=Nont5+1
      Goto 3400
c
c  NELM=12.  3-He reaction in detector with possible added neutron.
   80 En=Echrg(1)
      N=6
      Amtlt=1.25*ETOL(N,En)
      if(iesc.eq.0) FLstop=FLstop+1.25*ETOL(N,En)  ! Tatsu 04/6/28
c     3-He light is somewhat greater than alpha light for same E;
c       factor of 1.25 is approximately the difference in 3-He and 4-He
c       light in NE-102 reported by Becchetti, Thorn and Levine in
c       Nuclear Instrum. Methods 138 (1976) 93-104 (to within my ability
c       to read their fig. 5).
      L=1+IFIX(En)
      IF (L .GT. Sboxes) L=Sboxes
      Hspec(L)=Hspec(L)+1
c     That increments 3-He spectrum.  Test for multibody breakup reactions.
      IF (Echrg(4) .GT. 0.0) goto 84
      IF (Echrg(2) .GT. 0.0) goto 82
c     Program counter to here, not multibody breakup.
c     Therefore, get Be-ion light -- 10-Be light same as 9-Be light ...
      En=Echrg(3)
      N=4
      Amtlt=Amtlt+1.33*ETOL(N,En)
      if(iesc.eq.0) FLstop=FLstop+1.33*ETOL(N,En)  ! Tatsu 04/6/28
      CountBe=CountBe+1
c       Now test for which counter to increment.
      IF (Eneut .GT. 0.0) goto 88
      Non3He=Non3He+1
      Goto 3400
c
c   Next is multibody breakup reactions following (n,3-He) two-body collision.
c
c     Program counter to here, we have 10-Be --> 9-Be --> 8-Be --> 2 alphas.
   82 Nonh1=Nonh1+1
      I1=2
      I2=3
      Goto 51
c     Next portion the 9-Be decayed by proton emission to 8-Li
   84 N=2
      En = Echrg(4)
      L=1+IFIX(En)
      IF (L .GT. Sboxes) L=Sboxes
      Pspec(L)=Pspec(L)+1
c        That updates proton spectrum.
      if(iesc.eq.2 .and. Echrg(2).gt.0.0) goto 85 ! kajimoto 2011/03/18
      Amtlt=Amtlt+ETOL(N,En)
      if(iesc.eq.0.and.Echrg(2).le.Ethre) FLstop=FLstop+ETOL(N,En)  ! Tatsu 04/6/28
c       Proton light added -- check for proton escape.
      En=Echrg(2)
      IF (En .LE. 0.0) goto 85
      Nprtes=Nprtes+1
      Amtlt=Amtlt-ETOL(N,En)
      if(iesc.eq.0.and.Echrg(2).le.Ethre) FLstop=FLstop-ETOL(N,En)  ! Tatsu 04/6/28
   85 En=Echrg(6)
c     Test for 8-Li --> 7-Li --> triton + alpha.  If so, Echrg(6) > 0.0
      IF (En .GT. 0.0) goto 87
c     Program counter to here, have either 8-Li or 7-Li ion light to get.
c      To determine which test both neutron energies.
      N=4
      En=Echrg(3)
      C2Li=2.6-0.006*En
      Amtlt=Amtlt+C2Li*ETOL(N,En)
      if(iesc.eq.0) FLstop=FLstop+C2Li*ETOL(N,En)  ! Tatsu 04/6/28
      CountLi=CountLi+1
c        (See statement 56 for approximate Li-ion light.)
      IF (En2.LT.0.0 .AND. Eneut.GT.0.0) goto 86
      Nonh2=Nonh2+1
      Goto 3400
   86 Nonh3=Nonh3+1
      Goto 3400
c     Next is 7-Li --> triton + alpha situation.
   87 N=8
      Amtlt=Amtlt+0.8*ETOL(N,En)
      if(iesc.eq.0) FLstop=FLstop+0.8*ETOL(N,En)  ! Tatsu 04/6/28
c        (See statement 75 for approximate triton light.)
      Nonh4=Nonh4+1
      L=1+IFIX(En)
      IF (L .GT. Sboxes) L=Sboxes
      Tspec(L)=Tspec(L)+1
c     That takes care of the triton.  Now finish off with alpha.
      I1=5
      I2=5
      Goto 51
c
c     Program counter to statement 88 to get index for (n,3-He n) reaction.
   88 Non3Hen=Non3Hen+1
c
c
c   Now have light output, unattentuated by scintillator, for the specific
c     reaction studied.  Determine attenuation factor for light loss in the
c     scintillator, then correct light output for this loss.
 3400 Attenu=-Grho*(Ht-Zint)
      FLtotal=FLtotal+amtlt  ! Light Output before considering quenching by Tatsu 04/6/28
      Amtlt=Amtlt*EXP(Attenu)
      FLIGHT = FLIGHT + AMTLT
      IF (Nelm .NE. 10) NTCOL = NTCOL + 1
c     That's the total number of collisions (excluding photon interactions).
      ELL=EXTERP(Bins,Bno,Amtlt,1000,1)
      L=IFIX(ELL)
c check type out-----light unit bins expansion problem
      maxB=999
      if (L .le. maxB) goto 321
      write(*,*) 'ERROR=====>light-unit bins problem'
      write(*,*) 'L= ',L,' light output= ',Amtlt
c daiki (06.01.23)
      L = L -1
  321 continue
c        Save info on this scattering in case a deletion is needed due to
c         photon detection and P.S.D. rejection.  (If 30 such saved
c         scatterings take place, terminate event!)
c add by d.satoh ('05.02.09)
      if( Nelm .le. 0 ) goto 5445
      Lsave(Npsd)=1000*Nelm + L
c check type out  --- problem with calculations for En 6 MeV on up with
c     too many counts in the 475th light-unit box.  Where do they
c      come from???
c     IF (L.LT.999) goto 3399
c      Do 3398 Jpsd=1,Npsd
c 3398            write(6,3401)jpsd,Lsave(jpsd)
c 3401      Format('   *** For L=1000 -->  LSAVE('I2,')= 'I8)
c        write(6,3402)Echrg(1),Amtlt,Eneut
c 3402      format('    Echrg(1), Amtlt, Eneut ='1p3e12.4/)
c 3399      Continue
      Npsd=Npsd+1
      IF (Npsd .GT. 50) Eneut=0.0
c
C  BRANCH ON THE INDICATOR WHICH CONTAINS THE COLLISION-TYPE THAT
C    THE NEUTRON IS TO BE RECORDED AS.  THEN INCREMENT THE
C    APPROPRIATE COUNTERS.  Count Photon interactions, too.
C
 5445 continue
      Goto(230,240,245,250,260,270,280,290,275,300,277,249),NELM
      if (NELM.eq.-1) goto 301
      if (NELM.eq.-2) goto 302
      if (NELM.eq.-3) goto 303
      if (NELM.eq.-4) goto 304
      if (NELM.eq.-5) goto 305
      if (NELM.eq.-6) goto 306
  230 SB1H(L) = SB1H(L) + 1
      Return
  240 SBCNLY(L) = SBCNLY(L) + 1
      Return
  245 SBNNPR(L) = SBNNPR(L) + 1
      Return
  249 S3He(L) = S3He(L) + 1
      Return
  250 SB1A0H(L) = SB1A0H(L) + 1
      Return
  260 SB3A0H(L) = SB3A0H(L) + 1
      Return
  270 SBNP(L) = SBNP(L) + 1
      Return
  275 Sbnd(L) = Sbnd(L) + 1
      Return
  277 Sbnt(L) = Sbnt(L) + 1
      Return
  280 SBNPN(L) = SBNPN(L) + 1
      Return
  290 SBN2N(L) = SBN2N(L) + 1
      RETURN
  300 Sphot(L) = Sphot(L)+1
      Return
c add 05/Aug.
  301 SPIQ(L) = SPIQ(L) + 1
      return
  302 SBAQ(L) = SBAQ(L) + 1
      return
  303 SBDQ(L) = SBDQ(L) + 1
      return
  304 SBTQ(L) = SBTQ(L) + 1
      return
  305 SBHQ(L) = SBHQ(L) + 1
      return
  306 SBPQ(L) = SBPQ(L) + 1
      return
C
c   Next is the "P.S.D." event rejection -- deletion of information already
c     stored for this history.  Note that this deletion means that the
c     number of originally desired histories will be reduced by one on each
c     pass through this portion of the programming.
 4400 Npsd=Npsd-1
      Do 4490 L=1,Npsd
      I=Lsave(L)
      K=I/1000
      M=I-1000*K
c daiki(06.01.24)
      if( M .eq. 0 ) then
      write(*,*) 'Error in PSD!!'
      end if
c
      Goto (410,415,420,4490,425,430,4435,4440,4445,4490,4450,460), K
  410 SB1H(M)=SB1H(M)-1
      Goto 4490
  415 SBCNLY(M)=SBCNLY(M)-1
      Goto 4490
  420 SBNNPR(M)=SBNNPR(M)-1
      Goto 4490
  425 SB3A0H(M)=SB3A0H(M)-1
      Goto 4490
  430 SBNP(M)=SBNP(M)-1
      Nprtes=Nprej
      Goto 4490
 4435 SBNPN(M)=SBNPN(M)-1
      Nprtes=Nprej
      Goto 4490
 4440 SBN2N(M)=SBN2N(M)-1
      Goto 4490
 4445 SBND(M)=SBND(M)-1
      Ndlost=Ndrej
      Goto 4490
 4450 Sbnt(M)=Sbnt(M)-1
      Goto 4490
  460 S3He(M)=S3He(M)-1
 4490 Continue
c
      Flight=0.0
      Eneut=0.0
      En2 = 0.0
      Return
C
C
c           NEUTRON HISTORY ANALYSIS
C    BANKR3 RECORDS THE LIGHT PRODUCED DURING EACH NEUTRON HISTORY
c
      ENTRY BANKR3
c
      IF (Flight.LE.0.0) Return


c  Tatsu 7/4/04 for P.S.D.
      FLmax=30.0   ! Event with delayed light output above this value (MeVee) is identified as neutr
      Thratio=0.95 ! Event with FLratio (see below) above this value is identified as neutron event
      FLratio=FLstop/FLtotal  ! Delayed light output ratio
      if(iesc.eq.0.and.FLratio.lt.Thratio.and.FLstop.lt.FLmax) then
       Flight=0.0
       return
      endif


c  Tatsu 7/4/04 for output PSD result
c      if(flight*1.25.ge.0.2576) then ! above neutron threshold
c       Eve0=Eve0+1     ! Total Event
c       if(FLratio.ge.1.00.or.FLStop.ge.FLmax) Eve1=Eve1+1  ! Most Severe Case
c       if(FLratio.ge.0.95.or.FLStop.ge.FLmax) Eve2=Eve2+1  ! Not so severe Case
c       if(FLratio.ge.0.80.or.FLStop.ge.FLmax) Eve3=Eve3+1  ! Normal Case
c       if(FLratio.ge.0.50.or.FLStop.ge.FLmax) Eve4=Eve4+1  ! Normal Case
c       if(FLratio.ge.0.25.or.FLStop.ge.FLmax) Eve5=Eve5+1  ! Normal Case
c       if(FLratio.ge.0.01.or.FLStop.ge.FLmax) Eve6=Eve6+1  ! Generous Case
c      endif
cs
      Ell=EXTERP(Bins,Bno,Flight,1000,1)
      L=IFIX(Ell)
c daiki (05.11.25)
      if( L .eq. 1000 ) then
      L = L - 1
      end if
cs_01.07.22--------------------------------
**  Determin the max bin No.
      if( L .ge. maxl ) then
        maxl = L
      end if
cs_01.07.22--------------------------------


C     Increment the "Totals" counter.
      Sboxt(L)=Sboxt(L)+1
c
      Return
C
C
C
C                    OUTPUT CASE RESULTS
      ENTRY BANKR4
C       CALCULATE THE NUMBER OF COLLISIONS PER COLLIDING NEUTRON.
      FC = NTCOL
      FA=FLOAT(Nfstcl)
      IF (FA .GT. 0.0) COLPN = FC / FA
C       OUTPUT NEUTRON, COLLISION-TYPE, AND LIGHT INFORMATION.
c      Test to see if there were "P.S.D." deletions, and if so,
c        note that information first.
      IF (Nphot .LE. 0) goto 5005
      IF (Igflag .EQ. 1) Write(21,5000) Nphot
 5000 Format(/ 39H  Ne-213 PSD On; No. of photon rejects=,I10)
 5005 WRITE(21,5010) NFSTCL,Ht,NHONCE,NHYD,
     z  (Zsave(I),I=1,10),NCONCE,NCARB,NELAS
      Nna=Nalp+N3alp+Nona1+Nona2+Nona3+Nona4+Nona5
      IF (Nna .LE. 0) goto 5007
      Write(21,5001) NALP,
     a  N3ALP, Nona1, Nona2, Nona3, Nona4, Nona5
c
 5007 Nnp=Nonp+Nonpn+Nonp0+Nonp1+Nonp2n+Nonp2+Nonp5+Non2p
      Nnt=Nont+Nont1+Nont2+Nont3+Nont4
      Nnd=Nond+Nonda+Nondn+Nondn1+Nondn2+Nondn3
      Nn3h=Non3He+Non3Hen+Nonh1
      IF (Nnp+Nnt+Nnd+Nn3h .GE. 1)  goto 5009
      Write (21,5008) Nphot, Nprtes, Ndie, Colpn, Nclsns
 5008 Format(1H0,6x,30HNumber of PHOTON Interactions=, I12/
     a  7x,27HNumber of Escaped PROTONS =, I15/
     c  7x,32HNumber of NEUTRONS below Cutoff=, I10//
     d  41H No. of Collisions per Colliding NEUTRON=,0PF9.4/
     e  24H COLLISION DISTRIBUTION:/7x,2H 0,5x,1H1,5x,1H2,5x,1H3,5x,1H4,
     f  60H     5     6     7     8     9    10    11    12    13    14,
     g  36H    15    16    17    18    19   >19/I9,18I6,2I5/)
      Goto 5019
c
 5009 Write(21,5011) NONP, Nonp0, Nonpn, NONT, Nonp1, Nont1, Nonp2n,
     W Nont2, Nonp3, Nont3, Nonp2, Nont4, Nonp5, Nont5, Non2p,
     X Non2pn, Non3He, Nonp4, Non3Hen, Nonh1
c
      Write(21,5014) NOND, Nonh2, Nonda, Nonh3, Nondn, Nonh4,
     Z Nondn1, Nondn2, Nphot, Nondn3, Nprtes
c
      Write(21,5012)  NON2N, Ndlost, Non2np, Ndie,
     a   COLPN, NCLSNS
c
 5010 FORMAT (/ 38H NEUTRONS HAVING AT LEAST 1 COLLISION=,I11,6x,
     d  54HHydrogen-Collision Depth Distribution in Steps of 0.1*,
     r  F6.1,3H cm/
     e  38H NEUTRONS HAVING A HYDROGEN COLLISION=,I11,5x,
     E  10HFront Face, 3(6x,3H-->,6x), 9HP.M. Tube /
     f  39H NUMBER OF ELASTIC HYDROGEN COLLISIONS=,I10,I14,9I6/
     1  38H NEUTRONS HAVING ELASTIC C COLLISIONS=,I11/
     2  37H NUMBER OF ELASTIC CARBON COLLISIONS=,I12/
     D  39H NUMBER OF INELASTIC CARBON COLLISIONS=,I10)
 5001 Format(40H NEUTRONS HAVING AN (N,ALPHA) COLLISION=,I9/
     F  43H NEUTRONS HAVING AN (N,N-3ALPHA) COLLISION=,I6/
     f  41H   Neutrons having (N,Alpha P) Reactions=,I8/
     G  42H   Neutrons having (N,Alpha PN) Reactions=,I7/
     h  43H   Neutrons having (N,Alpha P2N) Reactions=,I6/
     I  44H   Neutrons having (N,2Alpha PNT) Reactions=,I5/
     J  45H   Neutrons having (N,2Alpha P2ND) Reactions=,I4)
 5011 FORMAT( 41H NEUTRONS HAVING AN (N,P 12-B) COLLISION=,I8/
     x  41H   Neutrons having (N,P ALPHA) Reactions=,I8/
     a  42H NEUTRONS HAVING AN (N,PN 11-B) COLLISION=,I7, 10x,
     y  40HNEUTRONS HAVING AN (N,T 10-B) COLLISION=,I13/
     A  42H   Neutrons having (N,PN ALPHA) Reactions=,I7,12x,
     z  40HNeutrons Having an (N,TP 9-Be) Reaction=,I11/
     B  39H   Neutrons having an (N,P2N) Reaction=,I10,12x,
     w  42HNeutrons Having an (N,TD 2ALPHA) Reaction=,I9/
     b  43H   Neutrons having (N,P2N ALPHA) Reactions=,I6,12x,
     v  45HNeutrons Having an (N,T ALPHA 6-Li) Reaction=,I6/
     C  43H   Neutrons with (N,PD2N 2ALPHA) Reactions=,I6,12x,
     u  43HNeutrons Having an (N,TPN 2ALPHA) Reaction=,I8/
     y  42H   Neutrons with (N,PNT 2ALPHA) Reactions=,I7,12x,
     Z  41HNeutrons having an (N,TPD 7-Li) Reaction=,I10/
     c  42H   Neutrons having (N,2P 11-Be) Reactions=,I7/
     D  43H   Neutrons having (N,2PN 10-Be) Reactions=,I6,10x,
     e  43HNEUTRONS HAVING AN (N,3-He 10-Be) REACTION=,I10/
     E  43H   Neutrons having (N,2P2N 9-Be) Reactions=,I6,12x,
     f  42HNeutrons having (N,3-He N 9-Be) Reactions=,I9/
     g  61x,45HNeutrons having (N,3-He 2N 2Alpha) Reactions=,I6)
 5014  Format(41H NEUTRONS HAVING AN (N,D 11-B) COLLISION=,I8,12x,
     G  43HNeutrons having (N,3-He NP 8-Li) Reactions=,I8/
     H  41H   Neutrons having (N,D Alpha) Reactions=,I8,12x,
     A  44HNeutrons having (N,3-He 2NP 7-Li) Reactions=,I7/
     f  43H   Neutrons having an (N,DN 10-B) Reaction=,I6,12x,
     B  46HNeutrons having (N,3-He 2NPT Alpha) Reactions=,I5/
     v  42H   Neutrons having (N,DNP 9-Be) Reactions=, I7/
     b  44H   Neutrons having (N,DND 2ALPHA) Reactions=, I5,10x,
     E  30HNumber of PHOTON Interactions=,I12/
     c  42H   Neutrons having (N,DN ALPHA) Reactions=, I7,10x,
     D  27HNumber of Escaped PROTONS =,I15)
 5012 Format(37H NEUTRONS HAVING AN (N,2N) COLLISION=,I12,
     A  10x, 29HNumber of Escaped DEUTERONS =I13/
     a  40H   Neutrons Having an (N,2NP) Collision=,I9,
     J  10x, 32HNumber of NEUTRONS below Cutoff=,I10//
     K  41H NO. OF COLLISIONS PER COLLIDING NEUTRON=0PF9.4/
     N  24H COLLISION DISTRIBUTION:/7X,2H 0,
     O  4X,2H 1,4X,2H 2,4X,2H 3,4X,2H 4,4X,2H 5,4X,2H 6,4X,2H 7,
     P  4X,2H 8,4X,2H 9,4X,2H10,4X,2H11,4X,2H12,4X,2H13,4X,2H14,
     Q  4X,2H15,4X,2H16,4X,2H17,4X,2H18,3X,2H19,2X,3H>19/
     R  1X,I8,18I6,2I5/)
c
 5019 WRITE(21,5020)
 5020 FORMAT (1H1,5X,11HLIGHT RANGE,3X,5HTOTAL,2X,5HHYDGN,1X,6HC-ELAS,
     A 14H  INEL  n,alfa,5x,10Hn,n3a  n,p,3x,8Hn,d  n,t
     B 3x,4Hn,pX,5x,5HGamma,4x,23HEquivalent Energy (MeV)/52x,
     C 6Hn,3-He,28x,10Hn,2n  Rays,3x,6HProtons,5x,5hAlpha4x,6HCarbon)
c
c     Suppress print-out of zeroes for largest "light-unit" bins.
      Kwit = Lboxes
      DO  5030  L = 1, LBOXES
      IF ( SBOXT(Kwit)  .NE.  0 )   GO TO 5040
 5030 Kwit=Kwit-1
 5040 Kwit=Kwit+2
      IF (Kwit .GT. Lboxes) Kwit=Lboxes
      Write (22,5045) Kwit,Ntcoll,Esourc,Elow,Temp
 5045 Format (2I10,1P3E10.3)
c
      Linect=1
      Do 5080 L=1,Kwit
c           (Added in new code -- equivalent particle energies output)
      Bavg=0.5*(Bins(L) + Bins(L+1))
      M=1
      Plt=ETOL(M,Bavg)
      M=3
      Clt=ETOL(M,Bavg)
      M=5
      Alt=ETOL(M,Bavg)
c           (Another addition -- "E" format for small light units)
      IF (Bins(L+1) .GE. 0.0011) goto 5050
      Bsl=Bins(L)*10000.
      Bsm=Bins(L+1)*10000.
      Write(21,5055) Bsl,Bsm,Sboxt(L),Sb1h(L),Sbcnly(L),Sbnnpr(L),
     A  Sb1a0h(L), S3He(L), Sb3a0h(L), Sbnp(L), Sbnd(L), Sbnt(L),
     B  Sbnpn(L), Sbn2n(L), Sphot(L), Plt,Alt,Clt
 5055 Format(F6.1,'E-4 -',F4.1,3HE-4,3I7,2I6,I5,3I6,
     A  I5,I6,2I5,2x,0p2F10.6,F9.5)
      Goto 5062
 5050 WRITE(21,5060) BINS(L),BINS(L+1),SBOXT(L),SB1H(L),
     1 SBCNLY(L),SBNNPR(L), SB1A0H(L), S3He(L), SB3A0H(L),
     2 SBNP(L), Sbnd(L), Sbnt(L), SBNPN(L), SBN2N(L), Sphot(L),
     3 Plt, Alt, Clt
 5060 FORMAT (F9.4,2H -,F7.4,3I7,2I6,I5,3I6,I5,I6,2I5,0P3F10.4)
 5062 Linect=Linect+1
      IF (Linect .LE. 56) goto 5065
      Write (21,5020)
      Linect=1
c comment out by daiki.
cc 5065 WRITE(22,5070) L,BINS(L),BINS(L+1),SBOXT(L),SB1H(L),
cc     1 SBCNLY(L),SBNNPR(L), SB1A0H(L), S3He(L),SB3A0H(L), SBNP(L),
cc     2 Sbnd(L), SBNPN(L), SBN2N(L), Sphot(L)
 5065 continue
c daiki modified...-----------------------------------------------------
      sarea = 3.141592e+0*Rdet**2
      denom = (BINS(L+1)-BINS(L))*1.25*REAL(Nhits) / sarea
      write(25, 5555) BINS(L)*1.25,BINS(L+1)*1.25,Sboxt(L)/denom,
     + SB1H(L)/denom, SBCNLY(L)/denom, SBNNPR(L)/denom, S3He(L)/denom,
     + SB1A0H(L)/denom,
     + SB3A0H(L)/denom, SBNP(L)/denom, SBND(L)/denom, SBNT(L)/denom,
     + SBNPN(L)/denom, SBN2N(L)/denom, SPHOT(L)/denom,
     + SBPQ(L)/denom, SBDQ(L)/denom, SBTQ(L)/denom, SBHQ(L)/denom,
     + SBAQ(L)/denom, SPIQ(L)/denom
 5555 format (1pe11.3,1x,1pe11.3,1x,1pe11.3,1x,1pe11.3,1x,1pe11.3,1x,
     &        1pe11.3,1x,1pe11.3,1x,1pe11.3,1x,1pe11.3,1x,1pe11.3,1x,
     &        1pe11.3,1x,1pe11.3,1x,1pe11.3,1x,1pe11.3,1x,1pe11.3,1x,
     &        1pe11.3,1x,1pe11.3,1x,1pe11.3,1x,1pe11.3,1x,1pe11.3,1x,
     &        1pe11.3)
c ----------------------------------------------------------------------
 5070 FORMAT (I4,1PE11.3,1X,1H-,1X,1PE10.3,5(1X,I5),7I5)
 5080 CONTINUE
      WRITE(21,5095)KWIT
 5095 FORMAT(/'  N(Rows)=',I4)
c
C       INITIALIZE SUMMATION VARIABLES, SUM EACH OF THE COLUMNS
C          FROM THE TOTAL ON,PRINT THE TOTALS AND RETURN.
      SUMA = 0
      SUMB = 0
      SUMC = 0
      SUMCP= 0
      SUMD = 0
      SUME = 0
      SUMF = 0
      SUMG = 0
      SUMH = 0
      SUMJ = 0
      Sumk = 0
      Sumt = 0
      S3HeX =0
c
      XSUM = 0.
      DO  5100  L = 1, KWIT
      SUMA = SUMA+SBOXT(L)
      SUMB = SUMB+SB1H(L)
      SUMC = SUMC + SBCNLY(L)
      SUMCP= SUMCP + SBNNPR(L)
      SUMD = SUMD+SB1A0H(L)
      SUME = SUME+SB3A0H(L)
      SUMF = SUMF + SBNP(L)
      SUMG = SUMG + SBNPN(L)
      SUMH = SUMH + SBN2N(L)
      SUMJ = SUMJ + Sbnd(L)
      Sumt = Sumt + Sbnt(L)
      Sumk = Sumk + Sphot(L)
      S3HeX= S3HeX+ S3He(L)
      XSUM=XSUM+(BINS(L+1)+BINS(L))*0.5*FLOAT(SBOXT(L))
c           (Last corrects a small error in original O5S programming)
 5100 CONTINUE
      WRITE(21,5110)SUMA,SUMB,SUMC,SUMCP,SUMD,S3HeX,SUME,SUMF,SUMJ,
     A Sumt,SUMG, SUMH,Sumk
 5110 FORMAT (15H0SUM OF COLUMNS,3X,3I7,2I6,I5,3I6,I5,I6,2I5)
      IF (Suma.EQ.0) STOP '    Stopped because SUMA = 0 in BANKR4.'
c           (That last is just in case --- 12/86)
      XLITE=XSUM/SUMA
      WRITE(21,5115)XLITE
 5115 FORMAT(' AVG LIGHT/INTERACTING NEUTRON=',F10.6)
c
c       Added feature -- output of charged-particle spectra in 1 MeV bins.
 6000 Format(29H1    Charged-particle spectra / 10x, 6HE-Resp, 4x,
     x  6HAlphas, 2x, 9HDeuterons, 3x, 7HProtons, 3x, 7HTritons,
     g  4x,4H3-He)
      Lmax=Sboxes
 6005 Continue
      J=Hspec(Lmax)+Alspec(Lmax)+Dspec(Lmax)+Pspec(Lmax)+Tspec(Lmax)
      IF (J.GT.0) goto 6010
      Lmax=Lmax-1
      IF (Lmax .LE. 0) Return
      Goto 6005
 6010 IF (Lmax .LT. Sboxes) Lmax=Lmax+1
      Suma=0
      Sumb=0
      Sumc=0
      Sumd=0
      Sume=0
      Write (21,6000)
      Do 6020 L=1,Lmax
      K=L-1
cs      IF (L .EQ. 57) Write(21,6000)
      Write (21,6015) K,L,Alspec(L),Dspec(L),Pspec(L),
     a Tspec(L),Hspec(L)
 6015 Format(I12,2H -,I2,I9,4I10)
      Suma=Suma+Alspec(L)
      Sumb=Sumb+Dspec(L)
      Sumc=Sumc+Pspec(L)
      Sumd=Sumd+Tspec(L)
      Sume=Sume+Hspec(L)
 6020 Continue
      Write(21,6030) Suma,Sumb,Sumc,Sumd,Sume
 6030 Format(13x,5HSums=,I7,4I10)
      Write(21,6040) CountLi,CountBe,CountB
 6040 Format(/10x,'No(Li ions)='I6/10x,'No(Be ions)='I6/
     x  10x,'No(B ions)='I7)


      SUMA=0
      SUMB=0
      SUMC=0
      SUMCP=0
      SUMD=0
      SUME=0
      SUMF=0
      SUMG=0
      SUMH=0
      SUMJ=0
      Sumk=0
cc    Sboxes=0
      S3HeX=0
      Sumt=0


cs
      do 6098 ij=1,999
      SBOXT(ij)=0
      SBNNPR(ij)=0
      SB1H(ij)=0
      SBCNLY(ij)=0
      SB1A0H(ij)=0
      SBND(ij)=0
      SB3A0H(ij)=0
      SBNP(ij)=0
      SBNPN(ij)=0
      SBN2N(ij)=0
      Sphot(ij)=0
      S3He(ij)=0
      Sbnt(ij)=0
c Add 05/Aug.
      SBPQ(ij)=0
      SBDQ(ij)=0
      SBTQ(ij)=0
      SBHQ(ij)=0
      SBAQ(ij)=0
      SPIQ(ij)=0
 6098 continue
      do 6099 ij=1,1000
      Alspec(ij)=0
      Dspec(ij)=0
      Pspec(ij)=0
      Tspec(ij)=0
      Hspec(ij)=0
 6099 continue
      do 6097 ij=1,10
      Zsave(ij)=0
 6097 continue


      CountLi=0
      CountBe=0
      CountB=0


c add by d.satoh ('05.02.01)
c Purpose is to initialize the TRIES array for continuous mode.


      do ijk = 1, 8
            do kji = 1, Nhist
                  tries(ijk,kji) = 0.0
            end do
      end do


      Return
      END
