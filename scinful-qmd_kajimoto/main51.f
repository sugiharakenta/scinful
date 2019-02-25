*********************************************************************************
*  SSSSSSSS   CCCCCCCC  IIIIIIIII NNN       NNN FFFFFFFFFF UUU     UUU LLL      *
* SSSSSSSSSS CCCCCCCCCC IIIIIIIII NNNN      NNN FFFFFFFFFF UUU     UUU LLL      *
* SSS    SSS CCCC   CCC    III    NNNNN     NNN FFF        UUU     UUU LLL      *
* SSSS       CCC           III    NNN NN    NNN FFF        UUU     UUU LLL      *
*  SSSSSS    CCC           III    NNN  NN   NNN FFFFFFFFFF UUU     UUU LLL      *
*    SSSSSS  CCC           III    NNN   NN  NNN FFF        UUU     UUU LLL      *
*       SSSS CCC           III    NNN    NN NNN FFF        UUU     UUU LLL      *
* SSS    SSS CCCC   CCC    III    NNN     NNNNN FFF        UUUU   UUUU LLL      *
* SSSSSSSSSS CCCCCCCCCC IIIIIIIII NNN      NNNN FFF         UUUUUUUUU  LLLLLLLL *
*  SSSSSSSS   CCCCCCCC  IIIIIIIII NNN       NNN FFF          UUUUUUU   LLLLLLLL *
*                                                                               *
*                                               QQQQQ    MMM     MMM DDDDDDD    *
*                                             QQQQQQQQQ  MMMM   MMMM DDDDDDDD   *
*                                            QQQ     QQQ MMMMM MMMMM DDD   DDD  *
*                               ==========   QQQ     QQQ MMM MMM MMM DDD    DDD *
*                               ==========   QQQ     QQQ MMM  M  MMM DDD    DDD *
*                                            QQQ     QQQ MMM     MMM DDD    DDD *
*                                            QQQ QQQ QQQ MMM     MMM DDD    DDD *
*                                             QQQQQQQQQ  MMM     MMM DDD   DDD  *
*                                               QQQQQQ   MMM     MMM DDDDDDDD   *
*                                                   QQQ  MMM     MMM DDDDDDD    *
*********************************************************************************
c                This is main program for SCINFUL-QMD code                      *
c                         programed by SATOH Daiki                              *
c                     Japan Atomic Energy Agency (JAEA)                         *
c                      mailto: satoh.daiki@jaea.go.jp                           *
c                                                                               *
c                           Last reviced in 2006.                               *
*********************************************************************************
      common /incidnt/ Enp,u,v,w,x,y,z
      common /clst/ no
      common /maxwel/ Esourc,e1,t
      common /ncll/ ncll
      common /init/ nhist,nhits
      common /randm/ irx
      common /try/ tries(8,1000000)
      common /naid/ R,H,rcollim,rz
      common /cutoff/ ecutoff
      common /cascade/ ids(20),xs(20),ys(20),zs(20),ks(20),
     +           us(20),vs(20),ws(20),ens(20)
      common /RES/ iswt, bias(5), icont, ictrl, iesc, ilight
*--------------------------------------------
      dimension ids2(20),xs2(20),ys2(20),zs2(20),ks2(20),
     +           us2(20),vs2(20),ws2(20),ens2(20)
*--------------------------------------------
      open(unit=20,file='Input.data',status='old')
      open(unit=21,file='NEW.lpt',status='old')
      open(unit=22,file='NEW.pch',status='old')
      open(unit=30,file='scin_bias.data',status='old')
      open(unit=31,file='scinful-qmd.txt')
      open(unit=23,file='energy.inp',status='old')
      open(unit=25,file='spectra.out',status='unknown')
*--------------------------------------------
      do ji = 1, 5
            bias(ji) = 0.0
      end do
c
cc    read(30,9) (bias(in),in=1,5)
cc    9     format(5(1x,f8.3))
      read(30,*) (bias(in),in=1,5)
*--------------------------------------------
c_Read the Input.data!!!
      ictrl = 1
  500 continue
c
      call input
      rewind(20)
c
      if ( icont .eq. 2 ) then
            read(23,*,end=501) Ene
            write(*,480) ene
  480       format(/,'c ++++++++++++++++++++++++++++++++++++++++++++++++'/
     +      ' Incident neutron energy =',f7.4)
            Esourc = ene
      end if
c
      Enp=Esourc
      n=nhist
      d=domega(n)
      e=float(nhist)/float(nhits)   ! zero-bias efficiency
*--------------------------------------------
      write(*,99) d
   99 Format(1h ,'Finished first part.  Domega = '1PE12.4,' sr.')
c
      write(*,11) Nhits,Nhist,esourc
   11 format(' Nhits = ',i10,' events,  Nhist = ',i7,
     +    ' events,  Energy = ',f9.4,' MeV')
c
      write(*,10) (bias(i),i=1,5)
c
      if( icont .eq. 1 ) then
            write(31,121) Nhist,Nhits,e,d
            write(21,121) Nhist,Nhits,e,d
  121       Format(1H ,'No. of Neutrons traced =',I10/
     X             1H ,'No. of Neutrons used   =',I10/
     X           1H ,'Zero-bias efficiency   =',0PF13.6/
     X             1H ,'Subtended Solid angle  =',1PE13.5,3H sr)
            write(31,*) 'c ------------------------------------------------------'
      end if
c
      if( ictrl .eq. 2 ) then
            write(31,10) (bias(i),i=1,5)
   10       format(1h ,'Bias ========> ',5(2x,f8.3))
c
            if( icont .eq. 2 ) then
            write(31,*)'c ------------------------------------------------------'
            write(31,*) 'Energy(MeV)  Eff(1)   Eff(2)   Eff(3)   Eff(4)   Eff(5)'
            write(31,*)'c ------------------------------------------------------'
            end if
c
      end if
c
*--------------------------------------------
c_Initiate all counting registers!!
      call bankr0
*--------------------------------------------
c         MULTI RUN CONTROL
*--------------------------------------------
***** Add by Tatsu 2005/2/20  ***************
      if(icont.eq.3) then
       n=100000000   ! for PHITSdump mode, n should be infinite
       open(unit=235,file='phits.dmp',form='unformatted',status='old')
       call readdump(ihan)
       if(ihan.ne.0) then
        write(*,*) 'There is no dump data from PHITS!'
        stop
       endif
      endif
***** End of Modification *******************


c_Repeat following routine until j=n
c in order to reduce statistical error!!
************* start *************************
      do 2000 j=1,n
      ip=1    !incident particle type (ex.,1:neutron)
      no=1
***** Changed by Tatsu *****
      if(ihan.ne.1) call bankr1  ! ihan=1 means the same history event, so need not initilization
***** End of Modification **
      ncll=0
*--------------------------------------------
      do ij1 =1,20
      ids(ij1)=0
      ks(ij1) =0.
      ens(ij1)=0.
      xs(ij1) =0.
      ys(ij1) =0.
      zs(ij1) =0.
      us(ij1) =0.
      vs(ij1) =0.
      ws(ij1) =0.
      end do
*--------------------------------------------
      ids(1)=ip
c ***** Changed by Tatsu for using PHITS dump file 2005/2/20 **
      if(icont.ne.3) then  ! normal mode
       x=tries(1,j)
       y=tries(2,j)
       z=tries(3,j)
       u=tries(4,j)
       v=tries(5,j)
       w=tries(6,j)
       Enp=tries(7,j)
       itype=tries(8,j)
      else                ! use PHITS dump file
       call phitsdump(ireact)
       if(ireact.eq.0) goto 300 ! no reaction occurred
       x=tries(1,1)
       y=tries(2,1)
       z=tries(3,1)
       u=tries(4,1)
       v=tries(5,1)
       w=tries(6,1)
       Enp=tries(7,1)
       itype=tries(8,1)
      endif
      xs(1)=x
      ys(1)=y
      zs(1)=z
      us(1)=u
      vs(1)=v
      ws(1)=w
      ens(1)=Enp
      ks(1)=itype
c ****** End of Modification by Tatsu ************************
*--------------------------------------------
      if(enp.le.150) then


cs121201
      if(enp .le. 0.0) then
          write(*,*) ' *** Error @ main.f; enp = ',enp
          goto 2000
      end if
cs121201


      call scin (itype,Enp,u,v,w,x,y,z)


      else
*--------------------------------------------
  222 continue
*--------------------------------------------
      do 123 ij2 =1,20
      ids2(ij2)=0
      ks2(ij2) =0.
      ens2(ij2)=0.
      xs2(ij2) =0.
      ys2(ij2) =0.
      zs2(ij2) =0.
      us2(ij2) =0.
      vs2(ij2) =0.
      ws2(ij2) =0.
  123 continue
*--------------------------------------------
      do 100 i=1,no
      ids2(i)=ids(i)
      ks2(i)=ks(i)
      ens2(i)=ens(i)
      xs2(i)=xs(i)
      ys2(i)=ys(i)
      zs2(i)=zs(i)
      us2(i)=us(i)
      vs2(i)=vs(i)
      ws2(i)=ws(i)
*--------------------------------------------
  100 continue
*--------------------------------------------
      do 124 ijk =1,20
      ids(ijk)=0
      ks(ijk) =0.
      ens(ijk)=0.
      xs(ijk) =0.
      ys(ijk) =0.
      zs(ijk) =0.
      us(ijk) =0.
      vs(ijk) =0.
      ws(ijk) =0.
  124 continue
*--------------------------------------------
c_initialize!
      nn=no
      no=0
*--------------------------------------------
      do 111 ii=1,nn
      ip=ids2(ii)
      itype=ks2(ii)
      enp=ens2(ii)
      x=xs2(ii)
      y=ys2(ii)
      z=zs2(ii)
      u=us2(ii)
      v=vs2(ii)
      w=ws2(ii)
*--------------------------------------------
      call jqmd (ip,itype)
      irx = irx * 663608941
  111 continue
*--------------------------------------------
      if( no .ge. 1 ) goto 222
*--------------------------------------------
      endif
*--------------------------------------------


c *** Changed by Tatsu for using PHITS dump file ******
  300 if(icont.eq.3) then
       call readdump(ihan)  ! ihan=0:history over, =1:next event has the same history, =2:End of Fil
       if(ihan.eq.2) then  ! End of dump file
        call bankr3
        call efft2  ! need not to call bankr4 because that causes an error
        stop
       endif
      endif
c_Increment counting registers, if one history is over
      if(ihan.ne.1) call bankr3
c *** End of Modification ******************************


*--------------------------------------------
 2000 continue
************* finish ************************


*--------------------------------------------
c_Calculation of neutron detection efficiency
c for 5-bias.
      call efft
      call bankr4
*--------------------------------------------
      if( icont .eq. 2 ) then
            goto 500
      end if


  501 continue
*--------------------------------------------
      close(20)
      close(21)
      close(22)
      close(23)
      close(30)
      close(31)
      close(33)
      close(25)
      stop
      end
