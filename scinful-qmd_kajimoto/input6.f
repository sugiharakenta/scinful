C   This is file INPUT.FOR
C     Purpose is to read input data file to start the computation!!
c     ( modified by D.Satoh @ 2001.07.19 )
c     ( modified by D.Satoh @ 2004.12.14 )
c     ( modified by D.Satoh @ 2005.02.02 )
c     ( modified by D.Satoh @ 2006.07.07 )
c--------------------------------------------------------------------
      Subroutine INPUT


      common /init/ Nhist,Nhits
      Common /MAXWEL/ Esourc,E1,T
      Common /RANDM/ Irx
      Common /LTABLE/ Itab,Ene(127),Hydl(127),Carl(127),Alpl(127),
     x Dlight(127)
      Common /NBOXES/ Wbox
csatoh20041214@jaeri
      COMMON /PROTON/NOPROT,RNGEN(144),PRAN(144),DRAN(144),ARAN(144),
     a              HRAN(144),TRAN(144)
c      Common /PROTON/ Nprt,Range(95),Epr(95)
      Common /TRY/ Tries(8,1000000)
      Common /NAID/ Radius,Height,Rcolim,Rdummy
      Common /DENS/ Apb(2)
      Common /CUTOFF/ Ecutoff
      Common /GFLAG/ Igflag
      Common /GAMATN/ Grho
      common /RES/ iswt, bias(5), icont, ictrl, iesc, ilight
      common /SCINT/ Scntyp


      dimension Title(20)
      Dimension Apb110(2),Apb213(2)
      Data Index/2/, Im1/-1/
      Data Apb110/0.05237,0.04743/


c       "Attenuation" densities for H and C in (atoms/barn)/cm for NE-110
c        based on a density of 1.032 gm/cm**3 and H:C ratio of 1.104


      Data Apb213/0.04833,0.03984/


c     Similarly for NE-213 based on a density of 0.874 gm/cm**3 and
c       H:C ratio of 1.213


CARD READ NO. 1   ---  a title, up to 80 characters
      Read (20,100,End=99) (Title(I),I=1,20)
  100 Format(20A4)
c
CARD READ NO. 2   --- (1)  No. of histories;
c           (2) Scintillator type - either 110. or 213. (or -213.; see
c               comment below) for NE-110 or NE-213, respectively;
c           (3) Light attenuation factor (in cm**-1) for light attenuation
c               by the scintillator, based on a report by Kuijper et al.
c               Nucl. Instrum. Methods 42 (l966) 56-60.
      Igflag=0
c
      Read (20,*,End=99) Nhist,Scntyp,Grho
c
      IF (Scntyp .EQ. 213.) goto 3
      IF (Scntyp .EQ. 110.) goto 92
c
c     If program counter gets to here, input SCNTYP was -213., indicating
c       a desire to use "pulse-shape-discrimination" computer methods to
c       exclude events which result in gamma rays being detected.
      Igflag=1
      Scntyp=213.
      Goto 3
c
   92 Apb(1)=Apb110(1)
      Apb(2)=Apb110(2)
c ----
c     d.satoh('06.07.03atJAEA)
c     Change the values of data-tables for NE-110.
c
      if( icntl .ne. 1 ) goto 4    ! The tables have already updated.
c
c           Change -RANGE- values for NE-110
c      Do 1 I=1,Nprt
ckajimoto    1 Range(I)=0.8546*Range(I)
c           Change -GRHO- for denser NE-110
      Grho=Grho*1.19
c
c      Do 33 I=2,124
      Do 33 I=2,127
      Aux=3.0
      IF (Ene(I) .GT. 0.081) goto 2
      Aux=1.0+2.0*Ene(I)/0.081
c        GLMorgan had added the factor of 3 increase in Carbon light;
c         spectra obtained with a 12-mm thick NE-110 detector support
c         some increase in carbon light, but how much is uncertain.
c         They also support some increase in alpha light.
    2 Carl(I)=Carl(I)*Aux
      Alpl(I)=1.75*Alpl(I)
      Dlight(I)=1.52*Dlight(I)
   33 Continue
c ----
c
      Goto 4
    3 Apb(1)=Apb213(1)
      Apb(2)=Apb213(2)
c
CARD READ NO. 3   ---  Neutron source energy information
c
    4 Read (20,*,End=99) Esourc,E1,T,Ecutoff
c    4 Read (20,104,End=99) Esourc,E1,T,Ecutoff
c  104      format(4F)
c
c     If input neutron energy < threshold for photon production
c       set flag IGFLAG.
      IF (Esourc .LE. 4.805) Igflag=1
      Isc=IFIX(Scntyp+0.0001)
      IF (T .EQ. 0.0) goto 6
      IF (Esourc .GE. E1) goto 7
      Zz=E1
      E1=Esourc
      Esourc=Zz
      IF (Esourc .LE. 4.805) Igflag=1
    7 IF (Esourc .LT. 0.15) Write (21,105) Esourc, E1
  105 Format(23H Upper Energy of Group ,1PE13.5,4H MeV/
     X       23H Lower Energy of Group ,1PE13.5,4H MeV)
      IF (Esourc .GE. 0.15) Write (21,205) Esourc, E1
  205 Format(23H Upper Energy of Group ,0PF9.4,4H MeV/
     x       23H Lower Energy of Group ,0PF9.4,4H MeV)
      IF (T.LT.0.0) goto 8
      Write (21,106) T
  106 Format(23H Maxwellian Avg. Temp.  ,0PF9.4,4H MeV)
      Go to 8
    6 IF (Esourc .LT. 0.1) Write (21,107) Esourc
  107 Format(32H Monoenergetic Source Energy is  1PE13.5,4H MeV)
      IF (Esourc .GE. 0.1) Write (21,207) Esourc
  207 Format(32H Monoenergetic Source Energy is  0PF9.4,4H MeV)
    8 IF (Ecutoff .LT. 10000.0) Write (21,108) Isc,Grho,Ecutoff
  108 Format(22H Scintillator Type NE-,I3,12x,
     W            'Light-attenuation factor is'F7.4,'/cm'/
     X       22H Low-energy Cut-off = ,0PF9.4,4H MeV)
      IF (Ecutoff .GE. 10000.0) Write (21,208) Isc, Grho, Ecutoff
  208 Format (22H Scintillator Type NE-,I3,12x,
     W            27HLight-attenuation factor is,F7.4,3H/cm /
     X        22H Low-energy Cut-off = ,1PE12.5,4H MeV/
     Y     6x,25HWhich Seems AWFULLY High  ,14x,9H<---***--  )
c
CARD READ NO. 4   ---  Starting random number (if < 0 program chooses one)
csatoh
      Read (20,*,End=99) Irx
c_satoh@jaeri
   10 Write (21,110) Irx
  110 Format(23H Initial Random Seed =  ,I11,10H (Base 10) )
csatoh
CARD READ NO. 5
      Read (20,*) Xsourc,Ysourc,Zsourc
c
      Tries(1,1)=Xsourc
      Tries(2,1)=Ysourc
      Tries(3,1)=Zsourc
c
CARD READ NO. 6   ---   detector information
csatoh
      Read (20,*,End=99) Radius, Height, Rcolim
c
CARD READ NO. 7   ---  response function   ! (add@2001.07.19)
      read(20,*) iswt, icont, iesc, ilight
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     iswt;
c      1 -> efficiency mode
c      2 -> response mode
ccc
c      icont;
c      1 -> single mode
c      2 -> multiple mode
c      3 -> tatsuhiko-san special (for PHITS)
ccc
c      iesc;
c      1 -> normal
c      2 -> eliminate escaping proton event
c      0 -> tatsuhiko-san special
ccc
c      ilight;
c      1 -> light output database in original SCINFUL
c      2 -> light output function proposed by Nakao et al.
c      3 -> light output function proposed by Satoh et al.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if ( iswt .lt. 1 .or. iswt .gt. 2 ) then
        write(*,*) 'ERROR @ CARD No. 7: iswt= ',iswt
        stop
      end if
c
      if ( icont .lt. 1 .or. icont .gt. 3 ) then
        write(*,*) 'ERROR @ CARD No. 7: icont= ',icont
        stop
      end if
c
      if ( iesc .lt. 0 .or. iesc .gt. 2 ) then
        write(*,*) 'ERROR @ CARD No. 7: iesc= ',iesc
        stop
      end if
c
      if ( ilight .lt. 1 .or. ilight .gt. 3 ) then
        write(*,*) 'ERROR @ CARD No. 7: ilight= ',ilight
        stop
      end if
c--------------------------------------------------------------------
c WRITING SECTION
      if ( ictrl .eq. 1 ) then
c ++ FOR SCREEN +++
            write(*,*)'  '
            write(*,*)'*** Input configuration!! ************'
            write(*,1000) Nhist,Scntyp,Grho
 1000       format(i7,f6.1,1pe13.5)
            write(*,1001) Esourc,E1,T,Ecutoff
 1001       format(3(f7.2),1pe13.5)
            write(*,1002) Irx
 1002       format(i10)
            Write(*,1003) Xsourc,Ysourc,Zsourc
 1003       format(3(f8.2))
            Write(*,1003) Radius,Height,Rcolim
            write(*,1004) iswt, icont, iesc, ilight
 1004       format( 4(i2))
ccc
          if ( iswt .eq. 1 ) then
              write(*,*) '* EFFICIENCY'
            else if ( iswt .eq. 2 ) then
                  write(*,*) '* RESPONSE'
                  open(unit=33,file='res.out',status='unknown')
                  write(33,*) '--- RESPONSE FUNCTION --------------------------'
c added by kajimoto 10/08/05
                   open(unit=34,file='res1.out',status='unknown')
            end if
c
            if ( icont .eq. 1 ) then
                  write(*,*) '* SINGLE calculation'
            else if( icont .eq. 2) then
                  write(*,*) '* MULTIPLE calculation'
            else
                  write(*,*) '* PHITS DUMP MODE'
            end if
c
            if ( iesc .eq. 2 ) then
              write(*,*) '* ELIMINATION of escaping proton events'
            end if
c
          if ( ilight .eq. 1 ) then
              write(*,*) '* LIGHT-OUTPUT DATABASE in SCINFUL'
          else if ( ilight .eq. 2 ) then
              write(*,*) '* LIGHT-OUTPUT FUNCTION proposed by Nakao et al.'
          else
              write(*,*) '* LIGHT-OUTPUT FUNCTION proposed by Satoh et al.'
          end if
c
            write(*,*)'**************************************'
c
c ++ FOR FILES +++
            write(31,*) 'c ------------------------------------------------------'
            Write(21,101) (Title(I),I=1,20)
            Write(31,101) (Title(I),I=1,20)
  101       format (1H ,20A4)
            write(31,34) Esourc
   34       format (1h ,'Energy     =', 1x,f7.2)
            write(31,35) tries(1,1),tries(2,1),tries(3,1)
   35       format (1h ,'Source Pos =', f8.2,1x,f8.2,1x,f8.2)
c
            if ( iesc .eq. 2 ) then
                write(31,*) '# ELIMINATION of escaping proton events'
          end if
c
          if ( ilight .eq. 1 ) then
              write(31,*) '# LIGHT-OUTPUT DATABASE in SCINFUL'
          else if ( ilight .eq. 2 ) then
              write(31,*) '# LIGHT-OUTPUT FUNCTION proposed by Nakao et al.'
          else
              write(31,*) '# LIGHT-OUTPUT FUNCTION proposed by Satoh et al.'
          end if
            write(31,*) 'c ------------------------------------------------------'
c
      end if
c
      ictrl = ictrl + 1
c
c--------------------------------------------------------------------
      Return
   99 Stop ' No more input data '
      END
