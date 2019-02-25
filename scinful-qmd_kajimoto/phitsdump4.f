c  This subroutine is for using PHITS dump file as source data
      subroutine phitsdump(ireact)
c     ireact = 0: no neutron reaction occurred, 1: neutron reaction occurred
      parameter(npara=10) ! number of parameter (kf,x,y,z,u,v,w,e,hist,weight)
      common /dumpara/ dump(npara)
      real*8 dump
      common /incidnt/ Enp,u,v,w,x,y,z
      Common /TRY/ Tries (8,1000000)
      Common /NAID/ Rdet,Ht,RC,Rz
      Common /PROTON/NOPROT,RNGEN(144),PRAN(144),DRAN(144),ARAN(144),
     a              HRAN(144),TRAN(144)
c      Common /PROTON/ NOPROT, RANGE(95), PROTEN(95)
      Common /COLLIS/ Nelm, Echrg(6)
      Common /Flightcom/ Flight
      Common /RANDM/ Ixx
      Common /russian/ ikill
      Common /NEUTRN/ Entmp,Utmp,Vtmp,Wtmp,Xtmp,Ytmp,Zint


      ireact=0
      N=0           ! Charged Particle Index
      Nterp=4       ! Logarithemic Interporation
      amtlt=0.0     ! Amount light of this event (only for charged particle)
      Echrg(1)=0.0      ! Incident Charged Particle Energy (MeV)
      Echrg(2)=0.0      ! Escaped Charged Particle Energy (MeV)


      kf=int(dump(1))  ! kf code of particle
      do ip=1,7
       tries(ip,1)=dump(ip+1)
      enddo


c      write(59,*) kf,dump(9),tries(3,1),tries(7,1),dump(10)


      if(ikill.eq.0) then  ! new history event
       if(dump(9).lt.ran(ixx)) then  ! Russian Roulette is performed by every history (NOT event!)
        ikill=1  ! killed history
        return
       endif
      endif
      ikill=2  ! survived history


c Neutron
      if(kf.eq.2112) then ! neutron
       call interact3(1,nindx)
       if(nindx.ne.0) ireact=1 ! neutron reaction occured
       return


c Charged Particle
      else  ! for charged particle


c Proton
       if(kf.eq.2212) then ! proton, this part is taken from kinema3.f
        Nelm=1  ! Index for Proton
        Echrg(1)=Tries(7,1) ! Proton Energy (MeV)
        Protpl=PLNGTH(Tries(4,1),Tries(5,1),Tries(6,1),
     1  Tries(1,1),Tries(2,1),Tries(3,1), Rdet,Ht)    ! Get Possible Flight Distance
        Rprot=EXTERP(RNGEN,PRAN,Echrg(1),Noprot,Nterp) ! Get Range of Proton
        RangeQ=Rprot-Protpl
        Zint=Tries(3,1)+Tries(6,1)/2.0*min(Rprot,Protpl) ! effective interaction point
        IF (RangeQ .gt. 0.0) then  ! Proton Escape !!
         Echrg(2)=EXTERP(RNGEN,PRAN,RangeQ,Noprot,Nterp)
c         write(59,*) 'Espcape',echrg(1),echrg(2),rprot,protpl
        endif


c Deuteron
       elseif(kf.eq.1000002) then  ! deuteron part taken from kinema3.f, line 2124
        Nelm=-3  ! Index for Deuteron
        Echrg(1)=Tries(7,1)*2.0  ! PHITS output is in (MeV/n), so convert to (MeV)
        IF (Echrg(1) .gt. 0.2) then
         Deutpl=PLNGTH(Tries(4,1),Tries(5,1),Tries(6,1),
     1   Tries(1,1),Tries(2,1),Tries(3,1), Rdet,Ht)    ! Get Possible Flight Distance
         Epofd=0.5*Echrg(1)
         Rdeut=2.0*EXTERP(RNGEN,DRAN,Epofd,Noprot,Nterp)
c        That gets the range of the deuteron using the proton range data.
         Zint=Tries(3,1)+Tries(6,1)/2.0*min(Rdeut,Deutpl) ! effective interaction point
         RangeQ=Rdeut-Deutpl
         IF (RangeQ .gt. 0.0) then  ! Deuteron Escape !!
          RangeQ=RangeQ
          Echrg(2)=EXTERP(RNGEN,DRAN,RangeQ,Noprot,Nterp)
         endif
        endif


c Triton
       elseif(kf.eq.1000003) then  ! Triton
        Nelm=-4  ! Index for Triton (same as deuteron)
        Echrg(1)=Tries(7,1)*3.0  ! PHITS output is in (MeV/n), so convert to (MeV)
        Zint=Tries(3,1) ! effective interaction point


c He-3
       elseif(kf.eq.2000003) then  ! He-3
        Nelm=-5  ! Index for He-3 (same as alpha)
        Echrg(1)=Tries(7,1)*3.0  ! PHITS output is in (MeV/n), so convert to (MeV)
        Zint=Tries(3,1) ! effective interaction point


c Alpha
       elseif(kf.eq.2000004) then
        Nelm=-2  ! Index for Alpha
        Echrg(1)=Tries(7,1)*4.0  ! PHITS output is in (MeV/n), so convert to (MeV)
        Zint=Tries(3,1) ! effective interaction point


c Pion+ & Pion-
       elseif(abs(kf).eq.211) then
        Nelm=-1   ! Index for Pion
        Echrg(1)=Tries(7,1)
        if(Echrg(1).gt.10.0) then
         path=PLNGTH(Tries(4,1),Tries(5,1),Tries(6,1),
     1   Tries(1,1),Tries(2,1),Tries(3,1), Rdet,Ht)    ! Get Possible Flight Distance
         isw = 1  ! Range Calculation
         if(kf.gt.0) then  ! Pion+
          id=20
          ic=1
         else              ! Pion-
          id=21
          ic=-1
         endif
         call dedxl(isw,id,ic,Echrg(1),rg,tls)
         Zint=Tries(3,1)+Tries(6,1)/2.0*min(rg,path) ! effective interaction point
         rangeq = rg - path
         if(rangeq.gt.0.0) then ! Pion Escape !!
          isw = 2  ! Energy Loss Calculation
          tls = 0.0 ! Total Energy Loss
          call dedxl(isw,id,ic,Echrg(1),path,tls)
          Echrg(2) = Echrg(1) - tls
         endif
        endif


c electron (11), positron (-11), photon(22), neutrino(+-12,14), muon(+-13)
       elseif(abs(kf).eq.11.or.kf.eq.22.or.abs(kf).eq.13.or.
     1 abs(kf).eq.12.or.abs(kf).eq.14) then
        return  ! currently ignore those event
       else
        write(*,*) 'Unknown Particle with kf=',kf
        return
       endif


c Calculate Charged Particle Light Output
       call bankr2
c       write(59,*) kf,Echrg(1),Echrg(2),dump(10),flight,Zint  ! for debug
      endif


      end


      subroutine readdump(ihan)
c     ihan = 0:normal, 1:with same history number, 2:End of file
      parameter(npara=10) ! number of parameter (kf,x,y,z,u,v,w,e,hist,weight)
      common /dumpara/ dump(npara)
      Common /Flightcom/ Flight
      Common /russian/ ikill
      real*8 dump


      ihan=0
      iold=int(dump(10)+0.01) ! old history number
c      write(59,*) iold,flight
   10 read(235,end=999) (dump(ip),ip=1,npara)
      inew=int(dump(10)+0.01) ! new history number
      if(iold.eq.inew) then  ! same history
       if(ikill.eq.1) then   ! ikill=1 means already killed history, so read new event until next hi
        goto 10
       elseif(ikill.eq.2) then ! ikill=2 means already survived history
        ihan=1
       else
        write(6,*) 'Russian Roulette Error, ikill=',ikill
        stop
       endif
      else
       ikill=0   ! for new history, reset ikill
      endif
      return
  999 ihan=2
      return
      end


      Subroutine INTERACT3(N,Nindx)
c               (The variable N is the running index of interactions)
      Common /RANDM/ Irx
      Common /TRY/ Tries (8,1000000)
      Common /MAXWEL/ Eu,Elow,T
      Common /NAID/ Rdet,Ht,Rcolim,R
      Common /PROB/ P(12), indx
      Dimension Enn(41),Sn(41)
      Data Pi/3.1415927/
c
      Nindx=0
      M=N
      X=Tries(1,M)
      Y=Tries(2,M)
      Z=Tries(3,M)
      U=Tries(4,M)
      V=Tries(5,M)
      W=Tries(6,M)
      Eneut=Tries(7,M)
      if(Eneut.lt.0.1.or.Eneut.gt.3000.0) return
      Flgt=PLNGTH(U,V,W, X,Y,Z, Rdet,Ht)
      if(Eneut.le.150.0) then       ! Original
       Totx=TOTALX(Eneut)
       Itype=IBOX(Flgt,Totx,Dist)
      else                          ! JQMD mode
       Totx=TOTALX2(Eneut)
       Itype=IBOX2(Flgt,Totx,Dist)
      endif
      IF (Itype .LE. 0) Return
      Tries(1,M)=X+U*Dist
      Tries(2,M)=Y+V*Dist
      Tries(3,M)=Z+W*Dist
      Tries(7,M)=Eneut
      Tries(8,M)=FLOAT(Itype)+0.00001
      Nindx=1
      Return
      END


cccccccccccccccccccccccccccccccc
      subroutine efft2   ! output only response function
cccccccccccccccccccccccccccccccc
      integer sboxt(999),total, pboxt(999)
      real effic(5),count(5), bias2(5)
      common /effnaka/  bins(1000), sboxt, maxl
      common /INIT/ Nhist,Nhits
      common /MAXWEL/ Esourc,E1,T
      common /RES/ iswt, bias(5), icont, ictrl, iesc, ilight


      data count / 0.0, 0.0, 0.0, 0.0, 0.0/


      Nhits=1  ! NOT normalized by source number because PHITS source is not equal to input data
      Nhist=1  ! NOT normalized by source number because PHITS source is not equal to input data
      do 5000 ijk=1,5
      bias2(ijk)=0.0
      effic(ijk)=0.0
      count(ijk)=0.0
5000  continue


c+++++ response function +++++
c_2003/06/16@jaeri


      if( iswt .eq. 2 ) then
          write(33,*) 'These data NOT normalized by source and fluence'
          write(33,*) 'LO_low(MeVee), LO_up(MeVee), Response(/MeVee),
     +Error(/MeVee)'


          do kk = 1,maxl
              serr = 0.
              if (sboxt(kk) .gt. 0) serr = sqrt(float(sboxt(kk)))
              binw = ( bins(kk+1) - bins(kk) ) * 1.25


              dlight = (bins(kk) + bins(kk+1) )/2. * 1.25        ! MeVee.
              cnt = float(sboxt(kk))
              resp = cnt / binw / float(nhits)
              err = serr / binw / float(nhits)


c daiki modified at jaeri 2004.12.13 -------------------------------------------
              write(33,40) bins(kk)*1.25,bins(KK+1)*1.25,resp,err
   40         format(1pe11.3,1x,1pe11.3,1x,1pe11.3,1x,1pe11.3)
c
          end do


      end if
c+++++ response function +++++


      return
      end


