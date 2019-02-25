************************************************************************
*                                                                      *
*        PART 2: Initialization of QMD part                            *
*                                                                      *
*                                                                      *
*   List of subprograms in rough order of relevance with main purpose  *
*      ( s = subroutine, f = function, b = block data, e = entry )     *
*                                                                      *
*                                                                      *
*  s  readcfg   to read the input parameters from configuration file   *
*  s  qmdint    to initialize the QMD part                             *
*  s  kcomp     to determine the parameters of the interactions        *
*  s  idname    to identifies the particle or nucleus                  *
*  s  chname    to give the particle and Nucleus name as character     *
*  f  ulmass    to give the mass of the partilce from kf code          *
*  f  luchge    to give the charge of the partilce from kf code        *
*                                                                      *
*                                                                      *
************************************************************************


************************************************************************
*                                                                      *
      subroutine readcfg
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to read the following parameters from                   *
*              input configuration file                                *
*                                                                      *
*----------------------------------------------------------------------*
*        icfg should be defined in user program.                       *
*                                                                      *
*           icfg = 0;  no input file                                   *
*           icfg = 1;  read input file of default name                 *
*           icfg = 2;  read input file name from unit 5                *
*           icfg = 3;  no input file and                               *
*                      take input data from main by characters         *
*           icfg = 4;  no input file and                               *
*                      take input data from main by variables          *
*                                                                      *
*----------------------------------------------------------------------*
*        icfg = 1, 2, 4                                                *
*                                                                      *
*              input       : content                     ; variables   *
*                                                                      *
*              'proj'      : projectile                  ; mstq1(1)    *
*                                                        ; mstq1(2)    *
*                                                        ; mstq1(3)    *
*              'targ'      : target                      ; mstq1(4)    *
*                                                        ; mstq1(5)    *
*                                                        ; mstq1(6)    *
*              'event'     : number of events            ; mstq1(7)    *
*              'tstep'     : total number of time step   ; mstq1(8)    *
*              'frame'     : reference frame             ; mstq1(9)    *
*                                                                      *
*                                                                      *
*              'win'       : incident energy or momentum ; parq1(1)    *
*                                                        ; parq1(2)    *
*              'bmin'      : minimum impact parameter    ; parq1(3)    *
*              'bmax'      : maximum impact parameter    ; parq1(4)    *
*              'dt'        : time step                   ; parq1(5)    *
*                                                                      *
*              'fname(i)'  : file name                   ; fname(i)    *
*              'mstq1(i)'  : integer parameters          ; mstq1(i)    *
*              'parq1(i)'  : real parameters             ; parq1(i)    *
*                                                                      *
*----------------------------------------------------------------------*
*        icfg = 3                                                      *
*                                                                      *
*              projc       : projectile                  (a80)         *
*              targc       : target                      (a80)         *
*              framec      : reference frame             (a80)         *
*              winc        : incident energy or momentum (a80)         *
*              bminc       : minimum impact parameter    (a80)         *
*              bmaxc       : maximum impact parameter    (a80)         *
*              dtc         : time step                   (a80)         *
*              tstepc      : total number of time step   (a80)         *
*              eventc      : number of events            (a80)         *
*                                                                      *
*----------------------------------------------------------------------*
*        store input parameters                                        *
*                                                                      *
*              icnum       : number of input data                      *
*              indata      : character input data                      *
*              indlng      : length of character input data            *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


      common /input1/ mstq1(mxpa1), parq1(mxpa1)
      common /input3/ projc, targc, framec, winc,
     &                bminc, bmaxc, dtc, tstepc,
     &                eventc
      common /swich2/ icfg, imany, icpus, idatm
      common /inecho/ icnum, indata(200), indlng(200)
      common /fnamid/ fname(mxfnm), idf(mxfnm), ifnl(mxfnm)


*-----------------------------------------------------------------------


      character       fname*80
      character       indata*80


      character*80    projc, targc, framec, winc,
     &                bminc, bmaxc, dtc, tstepc,
     &                eventc


*-----------------------------------------------------------------------


      character       chinp0*80, chinp*80, pname*80
      character       chalc*26, chauc*26


      logical         exex


      character       chind0(10)*80, chind(10)*80, chindw(10)*80
      dimension       lchind(10), mchind(10)


      data chindw / 'proj = ','targ = ','frame = ','win = ',
     &              'bmin = ','bmax = ','dt = ','tstep = ',
     &              'event = ','  '/


      data mchind / 7, 7, 8, 6, 7, 7, 5, 11, 8, 0 /


*-----------------------------------------------------------------------


      data chalc / 'abcdefghijklmnopqrstuvwxyz' /
      data chauc / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' /


*-----------------------------------------------------------------------
*        output unit for error message
*-----------------------------------------------------------------------


            ieo   = 6


*-----------------------------------------------------------------------


            ierr  = 0
            icnum = 0


*-----------------------------------------------------------------------
*     icfg = 2;  read input file name from unit 5
*-----------------------------------------------------------------------


         if( icfg .eq. 2 ) then


            write(ieo,*) ' *** Input Configuration File ***'


            read(5,'(a80)') fname(1)


         end if


*-----------------------------------------------------------------------
*     icfg = 3; input data in character base
*-----------------------------------------------------------------------


         if( icfg .eq. 3 ) then


               chind0(1) = projc
               chind0(2) = targc
               chind0(3) = framec
               chind0(4) = winc
               chind0(5) = bminc
               chind0(6) = bmaxc
               chind0(7) = dtc
               chind0(8) = tstepc
               chind0(9) = eventc


            do j = 1, 9


               do i = 1, 80


                  if( chind0(j)(i:i) .ne. ' ' ) goto 70


               end do


   70             i1 = i


               do i = 80, 1, -1


                  if( chind0(j)(i:i) .ne. ' ' ) goto 80


               end do


   80             i2 = i


                  lchind(j) = i2 - i1 + 1


               if( lchind(j) .gt. 0 ) then


                  chind0(j)(1:lchind(j)) = chind0(j)(i1:i2)


               else


                  lchind(j) = 1
                  chind0(j)(1:lchind(j)) = ' '


               end if


                  chind(j)(1:lchind(j)+mchind(j))
     &                    = chindw(j)(1:mchind(j))//
     &                      chind0(j)(1:lchind(j))
                  lchind(j) = lchind(j) + mchind(j)


            end do


         end if


*-----------------------------------------------------------------------
*        character length of file name
*-----------------------------------------------------------------------


            do j = 1, mxfnm


               do i = 1, 80


                  if( fname(j)(i:i) .ne. ' ' ) goto 50


               end do


   50             i1 = i


               do i = 80, 1, -1


                  if( fname(j)(i:i) .ne. ' ' ) goto 60


               end do


   60             i2 = i


                  ifnl(j) = i2 - i1 + 1


               if( ifnl(j) .gt. 0 ) then


                  fname(j)(1:ifnl(j)) = fname(j)(i1:i2)


               else


                  ifnl(j) = 0


               end if


               do i = ifnl(j) + 1, 80


                  fname(j)(i:i) = ' '


               end do


            end do




*-----------------------------------------------------------------------
*        parameters are given in main by mstq1 and parq1
*-----------------------------------------------------------------------


         if( icfg .eq. 4 ) return


*-----------------------------------------------------------------------
*        Open Configuration file
*-----------------------------------------------------------------------


         if( icfg .le. 2 ) then


            inquire( file = fname(1), exist = exex )


            if( exex .eqv. .false. ) then


               write(ieo,'('' Error: Input CFG File Does Not Exist'')')
               write(ieo,'('' ====='')')


               ierr = 1
               goto 200


            end if


            idi = idf(1)


            open( idi, file = fname(1) , status = 'old' )


         end if


*-----------------------------------------------------------------------
*     Read One Line from the file fname(1)
*-----------------------------------------------------------------------


               ntry  = 0


   11 continue


*-----------------------------------------------------------------------


         if( icfg .le. 2 ) then


               read( idi, '(a80)', end = 200 ) chinp


         end if


*-----------------------------------------------------------------------


         if( icfg .eq. 3 ) then


               incha = ntry + 1


               if( incha .gt. 9 ) goto 200


               chinp = chind(incha)(1:lchind(incha))


         end if


*-----------------------------------------------------------------------


               ntry = ntry + 1


*-----------------------------------------------------------------------
*        Convert character from upper to lower.
*-----------------------------------------------------------------------


               chinp0 = chinp


            do il = 1, 80


               do m = 1, 26


               if( chauc(m:m) .eq. chinp(il:il) )
     &             chinp(il:il) = chalc(m:m)


               end do


            end do


*-----------------------------------------------------------------------
*        Check chinp
*-----------------------------------------------------------------------


                  l1  = 0
                  l2  = 0
                  l3  = 0
                  l4  = 0


                  lrc = 0
                  llc = 0
                  l   = 1


                  ich = 0


            do while( l .le. 73 )


               if(      chinp(l:l+4) .eq. 'mstq1' .or.
     &                  chinp(l:l+4) .eq. 'parq1' .or.
     &                  chinp(l:l+4) .eq. 'fname' .or.
     &                  chinp(l:l+4) .eq. 'event' .or.
     &                  chinp(l:l+4) .eq. 'tstep' .or.
     &                  chinp(l:l+4) .eq. 'frame' ) then


                  l1 = l
                  l2 = l + 4
                  l  = l + 5


                  ich = l2


               else if( chinp(l:l+3) .eq. 'bmin' .or.
     &                  chinp(l:l+3) .eq. 'bmax' .or.
     &                  chinp(l:l+3) .eq. 'proj' .or.
     &                  chinp(l:l+3) .eq. 'targ' ) then


                  l1  = l
                  l2  = l + 3
                  l   = l + 4
                  ich = l2


               else if( chinp(l:l+1) .eq. 'dt' ) then


                  l1  = l
                  l2  = l + 1
                  l   = l + 2
                  ich = l2


               else if( chinp(l:l+2) .eq. 'win' ) then


                  l1  = l
                  l2  = l + 2
                  l   = l + 3
                  ich = l2


               else if( chinp(l:l) .eq. '#' .or.
     &                  chinp(l:l) .eq. '!' ) then


                  if( ich .gt. 0 ) ich = l


                  goto 130


               else if( chinp(l:l) .eq. '=' ) then


                  leq = l
                  ich = l


                  goto 110


               else if( chinp(l:l) .eq. '(' ) then


                  ich = l
                  llc = l + 1
                  l   = l + 1


               else if( chinp(l:l) .eq. ')' ) then


                  ich = l
                  lrc = l - 1
                  l   = l + 1


               else if( chinp(l:l+2) .eq. 'end' ) then


                  goto 200


               else if( chinp(l:l) .ne. ' ' .and.
     &                  chinp(l:l) .ne. ' ' ) then


                  ich = l
                  l   = l + 1


               else


                  l = l + 1


               end if


            end do


                  goto 130




*-----------------------------------------------------------------------
*        Get column number l3 and l4 of input value after '='.
*-----------------------------------------------------------------------


  110    continue


               if( l1 .eq. 0 .or. l2 .eq. 0 ) goto 130


               ibrank = 1


         do l = leq + 1, 80




            if( chinp(l:l) .ne. ' ' .and.
     &          chinp(l:l) .ne. '   ' .and.
     &          ibrank .eq. 1 ) then


               l3     = l
               ibrank = ibrank + 1
               ich    = l


            else
     &      if( chinp(l:l) .ne. ' ' .and.
     &          chinp(l:l) .ne. '   ' .and.
     &          ibrank .ne. 1 ) then


               ich = l


            else
     &      if( ( chinp(l:l) .eq. ' ' .or.
     &            chinp(l:l) .eq. ' ' ) .and.
     &            ibrank .ne. 1 ) then


               l4  = l - 1
               ich = l - 1


               goto 120


            end if


         end do


*-----------------------------------------------------------------------


            if( l3 .eq. 0 .or. l4 .eq. 0 ) goto 130


*-----------------------------------------------------------------------
*        Get information from the line
*-----------------------------------------------------------------------


  120    continue


*-----------------------------------------------------------------------
*        reference frame : mstq1(9)
*-----------------------------------------------------------------------


         if( chinp(l1:l1+4) .eq. 'frame' ) then


               if( chinp(l3:l3+1) .eq. 'cm' ) then


                  mstq1(9) = 1


               else if( chinp(l3:l3+2) .eq. 'lab' ) then


                  mstq1(9) = 0


               else if( chinp(l3:l3+1) .eq. 'nn' ) then


                  mstq1(9) = 2


               else


                  goto 999


               end if


*-----------------------------------------------------------------------
*        identification of projectile and targe
*
*        projectile : mstq1(1), mstq1(2), mstq1(3)
*        target     : mstq1(4), mstq1(5), mstq1(6)
*-----------------------------------------------------------------------


         else if( chinp(l1:l1+3) .eq. 'proj' ) then


                  pname = chinp(l3:l4)


               call idname(ieo,ierr,pname,l4-l3+1,idnum,idmas,idchg)


                  if( ierr .eq. 1 ) goto 998


                  mstq1(1) = idnum
                  mstq1(2) = idmas
                  mstq1(3) = idchg


         else if( chinp(l1:l1+3) .eq. 'targ' ) then


                  pname = chinp(l3:l4)


               call idname(ieo,ierr,pname,l4-l3+1,idnum,idmas,idchg)


                  if( ierr .eq. 1 ) goto 998


                  mstq1(4) = idnum
                  mstq1(5) = idmas
                  mstq1(6) = idchg


*-----------------------------------------------------------------------
*        incident energy or momentum
*
*        parq1(1), parq1(2)
*-----------------------------------------------------------------------


         else if( chinp(l1:l1+2) .eq. 'win' ) then


                  elab = -1.0
                  plab = -1.0


               if( chinp(l4-3:l4) .eq. 'gevc' ) then


                  read( chinp(l3:l4-4), '(g25.0)', err = 999 ) plab


               else if( chinp(l4-2:l4) .eq. 'gev' ) then


                  read( chinp(l3:l4-3), '(g25.0)', err = 999 ) elab


               else if( chinp(l4-3:l4) .eq. 'mevc' ) then


                  read( chinp(l3:l4-4), '(g25.0)', err = 999 ) plab


                  plab = plab / 1000.


               else if( chinp(l4-2:l4) .eq. 'mev' ) then


                  read( chinp(l3:l4-3), '(g25.0)', err = 999 ) elab


                  elab  = elab / 1000.


               else


                  goto 999


               end if


                  parq1(1) = elab
                  parq1(2) = plab


*-----------------------------------------------------------------------
*        bmin, bmax, dt, tstep, event
*
*        parq1(3), parq1(4), parq1(5)
*        mstq1(7), mstq1(8)
*-----------------------------------------------------------------------


         else if( chinp(l1:l1+3) .eq. 'bmin' ) then


               read( chinp(l3:l4), '(g20.0)', err = 999 ) parq1(3)


         else if( chinp(l1:l1+3) .eq. 'bmax' ) then


               read( chinp(l3:l4), '(g20.0)', err = 999 ) parq1(4)


         else if( chinp(l1:l1+1) .eq. 'dt' ) then


               read( chinp(l3:l4), '(g20.0)', err = 999 ) parq1(5)


         else if( chinp(l1:l1+4) .eq. 'event' ) then


               read( chinp(l3:l4), '(i12)', err = 999 ) mstq1(7)


         else if( chinp(l1:l1+4).eq.'tstep') then


               read( chinp(l3:l4), '(i12)', err = 999 ) mstq1(8)


*-----------------------------------------------------------------------
*        file name
*        fname(i)
*-----------------------------------------------------------------------


         else if( chinp(l1:l1+4) .eq. 'fname' ) then


               read( chinp(llc:lrc), '(i12)', err = 999 ) i


               if( i .gt. mxfnm ) goto 999


               fname(i) = chinp0(l3:l4)


               ifnl(i) = l4 - l3 + 1


*-----------------------------------------------------------------------
*        mstq1(i), parq1(i)
*-----------------------------------------------------------------------


         else if( chinp(l1:l1+4) .eq. 'mstq1' ) then


               read( chinp(llc:lrc), '(i12)', err = 999 ) i


               if( i .gt. mxpa1 ) goto 999


               read( chinp(l3:l4),   '(i12)', err = 999 ) ntmp


               mstq1(i) = ntmp


         else if( chinp(l1:l1+4) .eq. 'parq1' ) then


               read( chinp(llc:lrc), '(i12)',   err = 999 ) i


               if( i .gt. mxpa1 ) goto 999


               read( chinp(l3:l4),   '(g20.0)', err = 999 ) tmp


               parq1(i) = tmp


*-----------------------------------------------------------------------
*        error
*-----------------------------------------------------------------------


         else


               goto 999


         end if


*-----------------------------------------------------------------------
*        store input lines
*-----------------------------------------------------------------------


               icnum = icnum + 1


               l2m = max(l2,lrc+1)


               indata(icnum)(1:l2m-l1+1)  = chinp0(l1:l2m)
               indata(icnum)(l2m-l1+2:11) = ' '
               indata(icnum)(12:17)       = '& = & '
               indata(icnum)(18:l4-l3+18) = chinp0(l3:l4)


               indlng(icnum) = l4-l3+18


               goto 11


*-----------------------------------------------------------------------
*        Warning for unrecognized parameter name
*-----------------------------------------------------------------------


  130    continue


            if( ich .gt. 1 ) then


                  write(ieo,'('' Warning: Unrecognized input'',
     &                        '' parameter name'')')
                  write(ieo,'('' ========'')')


               if( icfg .le. 2 ) then


                  write(ieo,'('' '',a15,'' '',i3,'': '',80a1)')
     &            fname(1)(1:15), ntry, ( chinp0(j:j), j = 1, ich )
                  write(ieo,*)


               else if( icfg .eq. 3 ) then


                  write(ieo,'('' in main : '',80a1)')
     &            ( chind(ntry)(j:j), j = 1, lchind(ntry) )
                  write(ieo,*)


               end if


            end if


               goto 11


*-----------------------------------------------------------------------
*     Error occured in input parameter description
*-----------------------------------------------------------------------


  999 continue


                  write(ieo,'('' Error: Unrecognized input'',
     &                        '' parameter description'')')
                  write(ieo,'('' ======'')')


  998 continue


               if( icfg .le. 2 ) then


                  write(ieo,'('' '',a15,'' '',i3,'': '',80a1)')
     &            fname(1)(1:15), ntry, ( chinp0(j:j), j = 1, ich )
                  write(ieo,*)


               else if( icfg .eq. 3 ) then


                  write(ieo,'('' in main : '',80a1)')
     &            ( chind(ntry)(j:j), j = 1, lchind(ntry) )
                  write(ieo,*)


               end if


                  ierr = 1


*-----------------------------------------------------------------------
*     End of input data
*-----------------------------------------------------------------------


  200 continue


               if( icfg .le. 2 ) then


                  close(idi)


               end if


*-----------------------------------------------------------------------
*     Stop by error code
*-----------------------------------------------------------------------


         if( ierr .eq. 1 ) then


            write(ieo,'(/
     &          '' Program is stopped at subroutine READCFG'')')


            stop


         end if


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine qmdint
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to initialize the input parameters for QMD              *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


      common /input1/ mstq1(mxpa1), parq1(mxpa1)


      common /const0/ idnta, idnpr, massta, masspr, mstapr, msprpr
      common /const1/ elab, rdist, bmin, bmax, ibch, ibin
      common /const2/ dt, ntmax, iprun, iprun0
      common /const3/ nfreq, nfrec, nfred


      common /const4/ plab, srtcm, ylabb, pincm
      common /const5/ pzpr, pxpr, rzpr, rxpr
      common /const6/ pzta, pxta, rzta, rxta
      common /const7/ betpr, gampr, prmas, radpr
      common /const8/ betta, gamta, tamas, radta


      common /vriab0/ massal, massba, nmeson
      common /vriab1/ b, llnow, ntnow
      common /rannum/ iseed, iseed0, iseed1


      common /impact/ bval(1000), ibnum(1000), bweight(1000), bdef


      common /poten0/ aaa, bbb, rpot, esymm
      common /poten1/ gamm, c0, c3, cs, cl, wl
      common /poten2/ t0, t3, rkk


      common /swich0/ icoul, irelcr, iavoid
      common /swich1/ ipot, insys, irkg, icolt
      common /swich2/ icfg, imany, icpus, idatm
      common /swich3/ ielst, jelst, kelst
      common /swich4/ ifin, ifout


      common /grndc0/ dsam, ddif, dsam2, ddif2
      common /grndc1/ cdp, c0p, c3p, csp, clp
      common /grndc2/ r00, r01, saa, rada, radb
      common /grndc3/ ipchs, mntry, dtg, fric


      common /gradu0/ c0g, c3g, csg, pag
      common /caldis/ c0w, c3w, clw, c0sw
      common /pauli0/ cpw, cph, cpc


      common /startf/ iday0,imon0,iyer0,ihor0,imin0,isec0


      common /coultr/ eccm, pzcc, rmax0, zeroz
      common /framtr/ betafr(0:2), gammfr(0:2)


      common /cood2b/ rha(nnn,nnn), rhe(nnn,nnn), rhc(nnn,nnn)
      common /cood2r/ rr2(nnn,nnn), rbij(nnn,nnn), pp2(nnn,nnn)


      common /fnamid/ fname(mxfnm), idf(mxfnm), ifnl(mxfnm)
      character       fname*80


      common /summ04/ ianal


      logical         exex


*-----------------------------------------------------------------------
*           output unit for error message
*-----------------------------------------------------------------------


                  ieo  = 6


*-----------------------------------------------------------------------
*           Basic Input [1-9]
*-----------------------------------------------------------------------


                  idnpr  = mstq1(1)
                  masspr = mstq1(2)
                  msprpr = mstq1(3)
                  idnta  = mstq1(4)
                  massta = mstq1(5)
                  mstapr = mstq1(6)
                  iprun  = mstq1(7)
                  ntmax  = mstq1(8)
                  insys  = mstq1(9)


                  elab   = parq1(1)
                  plab   = parq1(2)
                  bmin   = parq1(3)
                  bmax   = parq1(4)
                  dt     = parq1(5)


*-----------------------------------------------------------------------
*           Numerics and Control [10-29]
*-----------------------------------------------------------------------


                  iseed  = mstq1(10)
                  ibch   = mstq1(11)
                  ibin   = mstq1(12)
                  irkg   = mstq1(13)
                  imany  = mstq1(14)
                  ifout  = mstq1(15)
                  ifin   = mstq1(16)
                  ielst  = mstq1(17)


*-----------------------------------------------------------------------
*           Interaction [30-59]
*-----------------------------------------------------------------------


                  ipot   = mstq1(30)
                  icoul  = mstq1(31)
                  irelcr = mstq1(32)


*-----------------------------------------------------------------------


                  wl     = parq1(30)
                  rpot0  = parq1(31)
                  esymm  = parq1(32) * 0.001


*-----------------------------------------------------------------------
*           Collision [60-89]
*-----------------------------------------------------------------------


                  icolt  = mstq1(60)
                  iavoid = mstq1(61)


*-----------------------------------------------------------------------
*           Ground state [90-119]
*-----------------------------------------------------------------------


                  ipchs  = mstq1(90)
                  mntry  = mstq1(91)


*-----------------------------------------------------------------------


                  saa    = parq1(90)
                  r00    = parq1(91)
                  r01    = parq1(92)


                  rada   = parq1(93)
                  radb   = parq1(94)


                  dtg    = parq1(95)
                  fric   = parq1(96)


                  rdist  = parq1(97)


*-----------------------------------------------------------------------
*           Analysis [150-200]      detect cpu time, date and time.
*                                   the others are in sm_init
*-----------------------------------------------------------------------


                  icpus  = mstq1(169)
                  idatm  = mstq1(170)


                  ianal  = mstq1(175)


*-----------------------------------------------------------------------
*        Initial vlues of llnow, ntnow and iprun0
*-----------------------------------------------------------------------


                  llnow  = 0
                  ntnow  = 0
                  iprun0 = iprun


*-----------------------------------------------------------------------
*        Start Total CPU time, Date and Time
*-----------------------------------------------------------------------


                  call cputime(1)
                  call datetime(iyer0,imon0,iday0,ihor0,imin0,isec0)


*-----------------------------------------------------------------------
*        Initialization of random number
*-----------------------------------------------------------------------


                  call ranint


*-----------------------------------------------------------------------
*        Intialization Multi Run Control
*-----------------------------------------------------------------------


                  call howmany(0)


*-----------------------------------------------------------------------
*        Set the mass of projectile and target particle
*-----------------------------------------------------------------------


                  prmas  = ulmass( idnpr )
                  tamas  = ulmass( idnta )


*-----------------------------------------------------------------------
*        Set beam energy and momentum.
*-----------------------------------------------------------------------


               if( elab .gt. 0.0 .and. plab .le. 0.0 ) then


                     plab = sqrt( elab * ( 2.0 * prmas + elab ) )


               else if( elab .le. 0.0 .and. plab .gt. 0.0 ) then


                     elab = sqrt( prmas**2 + plab**2 ) - prmas


               else if( elab .gt. 0.0 .and. plab .gt. 0.0 ) then


                     plabp = sqrt( elab * ( 2.0 * prmas + elab ) )


                  if( abs( plab - plabp ) .gt. 0.000001 ) then


                     write(ieo,'('' Error: elab and plab'',
     &                           '' is mismatched'')')
                     write(ieo,'('' ====='')')
                     write(ieo,'(/
     &                '' Program is stopped at subroutine qmdint'')')


                     stop


                  else


                     plab = plabp


                  end if


               else


                     elab = 0.0
                     plab = 0.0


               end if


*-----------------------------------------------------------------------
*        Set parameters according to the choice of frame.
*-----------------------------------------------------------------------


                  n1 = masspr
                  n2 = massta


*-----------------------------------------------------------------------
*         no target is error
*-----------------------------------------------------------------------


         if( n2 .eq. 0 ) then


                  write(ieo,'('' Error: Target is not specified'')')
                  write(ieo,'('' ====='')')
                  write(ieo,'(/
     &                '' Program is stopped at subroutine qmdint'')')


                  stop


*-----------------------------------------------------------------------


         else if( n1 .gt. 0 ) then


*-----------------------------------------------------------------------


                  ee   = elab + prmas + tamas
                  srt  = sqrt( ee**2 - plab**2 )
                  pstn = pcmsr( srt, prmas, tamas )


                  ptot = plab * n1
                  etot = elab * n1 + prmas * n1 + tamas * n2
                  stot = sqrt( etot**2 - ptot**2 )
                  pstt =  pcmsr( stot, prmas*n1, tamas*n2 )


                  pzcc = pstt
                  eccm = stot - ( prmas * n1 + tamas * n2 )


*-----------------------------------------------------------------------
*           Events defined in the CM frame.
*-----------------------------------------------------------------------


            if( insys .eq. 1 ) then


                  p1   =  pstt / n1
                  p2   = -pstt / n2
                  e1   =  sqrt( prmas**2 + p1**2 )
                  e2   =  sqrt( tamas**2 + p2**2 )


                  ylab = 0.5 * log( ( etot + ptot ) / ( etot - ptot ) )


                  betacm = 0.0
                  gammcm = 1.0


                  betalb = - ptot / etot
                  gammlb =   etot / stot


                  betann = ( p1 + p2 ) / ( e1 + e2 )
                  gammnn = ( e1 + e2 )
     &                   / sqrt( ( e1 + e2 )**2 - ( p1 + p2 )**2 )


*-----------------------------------------------------------------------
*           Events defined in Lab ( fixed target ) frame.
*-----------------------------------------------------------------------


            else if( insys .eq. 0 ) then


                  p1   = plab
                  p2   = 0.0
                  e1   = sqrt( prmas**2 + p1**2 )
                  e2   = sqrt( tamas**2 + p2**2 )


                  ylab = 0.0


                  betacm = ( p1 * n1 + p2 * n2 )
     &                   / ( e1 * n1 + e2 * n2 )
                  gammcm = ( e1 * n1 + e2 * n2 )
     &                   / sqrt( ( e1 * n1 + e2 * n2 )**2
     &                         - ( p1 * n1 + p2 * n2 )**2 )


                  betalb = 0.0
                  gammlb = 1.0


                  betann = ( p1 + p2 ) / ( e1 + e2 )
                  gammnn = ( e1 + e2 )
     &                   / sqrt( ( e1 + e2 )**2 - ( p1 + p2 )**2 )


*-----------------------------------------------------------------------
*           Frame defined in NN CM.
*-----------------------------------------------------------------------


            else if( insys .eq. 2 ) then


                  p1   =  pstn
                  p2   = -pstn


                  e1   = sqrt( prmas**2 + p1**2 )
                  e2   = sqrt( tamas**2 + p2**2 )


                  ylab = 0.5 * log( ( ee + plab ) / ( ee - plab ) )


                  betacm = ( p1 * n1 + p2 * n2 )
     &                   / ( e1 * n1 + e2 * n2 )
                  gammcm = ( e1 * n1 + e2 * n2 )
     &                   / sqrt( ( e1 * n1 + e2 * n2 )**2
     &                         - ( p1 * n1 + p2 * n2 )**2 )


                  betalb = - plab / ee
                  gammlb =   ee / srt


                  betann = 0.0
                  gammnn = 1.0


*-----------------------------------------------------------------------
*           Unrecognize frame : Error
*-----------------------------------------------------------------------


            else


                  write(ieo,'('' Error: Unrecognized input'',
     &                        '' coordinate frame'')')
                  write(ieo,'('' ====='')')


                  write(ieo,'(/
     &                '' Program is stopped at subroutine qmdint'')')


                  stop


            end if


*-----------------------------------------------------------------------
*           Set values
*-----------------------------------------------------------------------


                  srtcm  = srt
                  ylabb  = ylab
                  pincm  = pstn


                  pzpr   = p1
                  pxpr   = 0.0
                  betpr  = p1 / e1
                  gampr  = e1 / prmas


                  pzta   = p2
                  pxta   = 0.0
                  betta  = p2 / e2
                  gamta  = e2 / tamas


                  betafr(0) = betalb
                  gammfr(0) = gammlb
                  betafr(1) = betacm
                  gammfr(1) = gammcm
                  betafr(2) = betann
                  gammfr(2) = gammnn


*-----------------------------------------------------------------------
*        Target only
*-----------------------------------------------------------------------


         else if( n1 .eq. 0 ) then


                  elab   = 0.0
                  plab   = 0.0


                  srtcm  = 0.0
                  ylabb  = 0.0
                  pincm  = 0.0


                  pzpr   = 0.0
                  pxpr   = 0.0
                  betpr  = 0.0
                  gampr  = 1.0


                  pzta   = 0.0
                  pxta   = 0.0
                  betta  = 0.0
                  gamta  = 1.0


                  eccm   = 0.0
                  pzcc   = 0.0


                  betafr(0) = 0.0
                  gammfr(0) = 1.0
                  betafr(1) = 0.0
                  gammfr(1) = 1.0
                  betafr(2) = 0.0
                  gammfr(2) = 1.0


         end if


*-----------------------------------------------------------------------
*        massal, massba
*-----------------------------------------------------------------------


                  massal = masspr + massta


                  massba = 0.0
                  if( idnpr .eq. 0 ) massba = massba + masspr
                  if( idnta .eq. 0 ) massba = massba + massta




*-----------------------------------------------------------------------
*        Radius of Nucleus
*-----------------------------------------------------------------------


                  radpr = 0.0
                  radta = 0.0


                  if( n1 .gt. 2 ) radpr = r00 * float( n1 )**(1./3.)
                  if( n2 .gt. 2 ) radta = r00 * float( n2 )**(1./3.)


*-----------------------------------------------------------------------
*        Initial Distance : rmax0
*-----------------------------------------------------------------------
*
*        rdist ( D = -1.0 fm ) : initial distance
*
*            < 0.0:
*                 elab < 200 MeV
*                       rmax0 = radta/gamta + radpr/gampr + 6.0
*                 200 MeV < elab < 1 GeV
*                       rmax0 = radta/gamta + radpr/gampr + 4.0
*                 elab > 1 GeV
*                       rmax0 = radta/gamta + radpr/gampr + 4.0
*                 elab > 10 GeV
*                       rmax0 = radta/gamta + radpr/gampr + 2.0
*            >=0.0: rmax0 = radta/gamma + radpr/gamma + rdist
*
*-----------------------------------------------------------------------


               if( rdist .lt. 0.1 ) then


                  if( elab .ge. 10.0 ) then


                     rminm = 2.0


                  else if( elab .ge. 1.0 ) then


                      rminm = 4.0


                  else if( elab .ge. 0.2 ) then


                     rminm = 4.0


                  else


                     rminm = 6.0


                  end if


               else


                     rminm = rdist


               end if


                  rmax0 = radpr / gampr + radta / gamta + rminm


*-----------------------------------------------------------------------
*        Shift distance of z and x.
*-----------------------------------------------------------------------


               rzpr  = -rmax0 * float( n2 ) / float( n2 + n1 )
               rxpr  =  0.0


               rzta  =  rmax0 * float( n1 ) / float( n2 + n1 )
               rxta  =  0.0


*-----------------------------------------------------------------------
*        shift distance for Lab. system
*-----------------------------------------------------------------------


            if( insys .eq. 0 ) then


               zeroz = rmax0 * float( n1 ) / float( n2 + n1 )
     &               + 16.0 - radpr - rmax0


            else


               zeroz = 0.0


            end if


*-----------------------------------------------------------------------
*        check of impact parameter
*-----------------------------------------------------------------------


            if( bmin .gt. bmax .or. bmin .lt. 0.0 ) then


               write(ieo,'('' Error: invalid input bmin or bmax'')')
               write(ieo,'('' ====='')')


               write(ieo,'(/
     &             '' Program is stopped at subroutine qmdint'')')


               stop


            end if


*-----------------------------------------------------------------------
*        check of irkg
*-----------------------------------------------------------------------


            if( irkg .ne. 2 .and. irkg .ne. 4 ) irkg = 2




*-----------------------------------------------------------------------
*        Choice of the impact parameter bin
*-----------------------------------------------------------------------


                  if( ibch .lt. 0  .or.  ibch .gt. 2 ) ibch = 1


                  if( ifout .gt. 0 .and. ibch .eq. 2 ) ibch = 1
                  if( ifin  .gt. 0 .and. ibch .eq. 1 ) ibch = 2


*-----------------------------------------------------------------------
*           initialization of impact parameter and weight
*-----------------------------------------------------------------------


               if( ibch .eq. 0 ) then


                     ibin = 1
                     bdef = ( bmax - bmin ) / float(ibin)


                     ibnum(1)   = 0
                     bval(1)    = bmin + bdef / 2.0
                     bweight(1) = 2.0 * pi * bval(1) * bdef * 10.0


               else if( ibch .eq. 1 ) then


                     if( iprun .lt. ibin ) ibin = iprun


                     jprun  = iprun / ibin
                     iprun  = jprun * ibin
                     iprun0 = iprun


                     bdef = ( bmax - bmin ) / float(ibin)


                  do i = 1, ibin


                     ibnum(i)   = jprun
                     bval(i)    = bmin + bdef / 2.0 + float(i-1) * bdef
                     bweight(i) = 2.0 * pi * bval(i) * bdef * 10.0
     &                          * float(ibin)


                  end do


               else if( ibch .eq. 2 ) then


                     bdef = ( bmax - bmin ) / float(ibin)


                  do i = 1, ibin


                     ibnum(i)   = 0
                     bval(i)    = bmin + bdef / 2.0 + float(i-1) * bdef
                     bweight(i) = 2.0 * pi * bval(i) * bdef * 10.0


                  end do


               end if


*-----------------------------------------------------------------------
*        Choice of the judgement of elastic collision
*-----------------------------------------------------------------------


                  if( ielst .gt. 3 .or. ielst .lt. 0 ) ielst = 3


*-----------------------------------------------------------------------
*           Reading data from files
*-----------------------------------------------------------------------


               if( ifin .gt. 0 ) then


                  do i = 1, ifin


                     if( ifnl(i+10) .le. 0 ) then


                        write(ieo,'('' Error: output data file is'',
     &                              '' not defined'')')
                        write(ieo,'('' ====='')')


                        write(ieo,'(/
     &                  '' Program is stopped at subroutine qmdint'')')


                        stop


                     end if


                        inquire( file = fname(i+10), exist = exex )


                     if( exex .eqv. .false. ) then


                        write(ieo,'('' Error: input data file'',
     &                              '' does not exist:'',
     &                              '' filename = '',80a1)')
     &                              (fname(i+10)(k:k), k=1,ifnl(i+10))
                        write(ieo,'('' ====='')')


                        stop


                     end if


                        open( idf(i+10), file = fname(i+10),
     &                        status = 'old' )


                  end do




                        call qmdreadc




                  do i = 1, ifin


                        rewind( idf(i+10) )


                  end do


               end if


*-----------------------------------------------------------------------
*           Open file for QMD output
*-----------------------------------------------------------------------


               if( ifout .gt. 0 ) then


                  open( idf(10), file = fname(10), status='unknown' )


               end if


*-----------------------------------------------------------------------
*        QMD with only skyrme : Set potential parameters
*-----------------------------------------------------------------------


               if( ipot .eq. 1 ) then


                  rpot = 0.333333


               else if( ipot .eq. 2 ) then


                  rpot = 0.666667


               else if( ipot .eq. 3 ) then


                  rpot = 1.0


               else


                  write(ieo,'('' Error: invalid input mstq1(30)'')')
                  write(ieo,'('' ====='')')


                  write(ieo,'(/
     &                '' Program is stopped at subroutine qmdint'')')


                  stop


               end if


               if( rpot0 .gt. 0.0 ) rpot = rpot0


*-----------------------------------------------------------------------


               gamm =  rpot + 1.0


               call kcomp(rpot,t3,t0,aaa,bbb,rkk)




*-----------------------------------------------------------------------
*        Local potentials.
*-----------------------------------------------------------------------


               c0 = aaa/(rho0*(4*pi*wl)**1.5*2.0)
               c3 = bbb/(rho0**gamm*(4.0*pi*wl)**(1.5*gamm)*(gamm+1.0))
               cs = esymm/(rho0*(4.0*pi*wl)**1.5*2.0)
               cl = ccoul/2.0 * icoul


*-----------------------------------------------------------------------
*        Pauli
*-----------------------------------------------------------------------


               cpw = 1.0 / 2.0 / wl
               cph = 2.0 * wl / hbc**2
               cpc = 4.0


*-----------------------------------------------------------------------
*        Ground State
*-----------------------------------------------------------------------


               dsam  = 1.5
               ddif  = 1.0


               dsam2 = dsam * dsam
               ddif2 = ddif * ddif


               cdp = 1.0 / ( 4.0 * pi * wl )**1.5


               c0p = c0 * 2.0
               c3p = c3 * ( gamm + 1.0 )
               csp = cs * 2.0
               clp = cl * 2.0


            if( masspr .eq. 0 ) then


               rzpr = 0.0


               bmin = 0.0
               bmax = 0.0


            end if


*-----------------------------------------------------------------------
*        gradu
*-----------------------------------------------------------------------


               c0g = - c0 / ( 2.0 * wl )
               c3g = - c3 / ( 4.0 * wl ) * gamm
               csg = - cs / ( 2.0 * wl )
               pag =   gamm - 1.0


*-----------------------------------------------------------------------
*        caldis
*-----------------------------------------------------------------------


               c0w  = 1.0 /   4.0 / wl
               c3w  = 1.0 /   4.0 / wl
               clw  = 2.0 / ( 4.0 * pi * wl )**0.5
               c0sw = sqrt(c0w)


*-----------------------------------------------------------------------
*        set diagonal elements of two-body quantities to be zero
*-----------------------------------------------------------------------


            do i = 1, nnn
*
               rha(i,i)  =   0.0
               rhe(i,i)  =   0.0
               rhc(i,i)  =   0.0
               rr2(i,i)  =   0.0
               pp2(i,i)  =   0.0
              rbij(i,i)  =   0.0
*
            end do


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine kcomp(rpot,t3,t0,aaa,bbb,rkk)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 26                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to set the parameter of Skyrme force and                *
*              calculate incompressibility in the case of simple       *
*              Skyrme interaction.                                     *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              rpot        pawer of density dependent                  *
*              t0, t3,     coefficients of Skyrme                      *
*              aaa, bbb    coefficients of Skyrme                      *
*              rkk         incompressibility                           *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------
*        Calculate nuclear matter properties.
*-----------------------------------------------------------------------


               ebin  = ebinm * 0.001
               pfer  = hbc * ( 3./2. * pi**2 * rho0 )**(1./3.)
               efer  = pfer**2 / 2. / rmass


*-----------------------------------------------------------------------
*        t0, t3 and aaa, bbb
*-----------------------------------------------------------------------


               t3   =  8./3./rpot/rho0**(1.+rpot)*(efer/5.-ebin)
               t0   = -16./15.*efer/rho0 - (1.+rpot)*t3*rho0**rpot


               aaa  =  3./4.*t0*rho0
               bbb  =  3./8.*t3*(2.+rpot)*rho0**(1.+rpot)


*-----------------------------------------------------------------------
*        Incompressibility
*-----------------------------------------------------------------------


               rkk  =  6./5.*efer+9./4.*t0*rho0
     &              + 9./8.*(1.+rpot)*(2.+3.*rpot)*t3*rho0**(1.+rpot)




*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine idname(ieo,ierr,pname,len,idnum,idmas,idchg)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 26                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to identifies the particle or nucleus                   *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              ieo         : file unit for error message               *
*              ierr        : error flag                                *
*                                                                      *
*              pname       : input pname in character*8                *
*              len         : length of the character name              *
*                                                                      *
*              idnum       : output id number                          *
*                            0  : non, nucleon, nucleus                *
*                            kf : others                               *
*              idmas       : output mass number                        *
*              idchg       : output charge number                      *
*                                                                      *
*        Called subroutine and function :                              *
*                                                                      *
*              function   luchge(kf)                                   *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


      parameter ( mxbeam =   6 )
      parameter ( mxnucl = 104 )


*-----------------------------------------------------------------------


      character       pname*80


      character       num(0:9)*1, cnuc*3


      character       elmnt(mxnucl)*3


      character       chcde(mxbeam)*8
      dimension       kcode(mxbeam), lchcd(mxbeam)


*-----------------------------------------------------------------------


      data chcde / 'n0      ','p+      ','pi+     ','pi-     ',
     &             'pi0     ','non     ' /


      data lchcd /          2,         2,         3,         3,
     &                      3,         3 /


      data kcode /       2112,      2212,       211,      -211,
     &                    111,     -9999 /


      data elmnt /
     &    'h  ','he ','li ','be ','b  ','c  ','n  ','o  ','f  ','ne ',
     &    'na ','mg ','al ','si ','p  ','s  ','cl ','ar ','k  ','ca ',
     &    'sc ','ti ','v  ','cr ','mn ','fe ','co ','ni ','cu ','zn ',
     &    'ga ','ge ','as ','se ','br ','kr ','rb ','sr ','y  ','zr ',
     &    'nb ','mo ','te ','ru ','rh ','pd ','ag ','cd ','in ','sn ',
     &    'sb ','te ','j  ','xe ','cs ','ba ','la ','ce ','pr ','nd ',
     &    'pm ','sm ','eu ','gd ','tb ','dy ','ho ','er ','tm ','yb ',
     &    'lu ','hf ','ta ','w  ','re ','os ','ir ','pt ','au ','hg ',
     &    'tl ','pb ','bi ','po ','at ','rn ','fr ','ra ','ac ','th ',
     &    'pa ','u  ','np ','pu ','am ','cm ','bk ','cf ','es ','fm ',
     &    'md ','no ','lw ','ku ' /


      data num / '0','1','2','3','4','5','6','7','8','9' /


*-----------------------------------------------------------------------


         ierr = 0


*-----------------------------------------------------------------------
*     particle or nucleus name,
*
*     Pre-Procedure : Modify the particle names if needed.
*-----------------------------------------------------------------------


                  inuc = 0
                  jnuc = 0


                  idnum = 0
                  idmas = 0
                  idchg = 0


*-----------------------------------------------------------------------
*           'n ' -> 'n0'
*-----------------------------------------------------------------------


            if( pname(1:2) .eq. 'n ' ) then


                  pname(1:8) = 'n0      '


*-----------------------------------------------------------------------
*           'p ' -> 'p+'
*-----------------------------------------------------------------------


            else if( pname(1:2) .eq. 'p ' ) then


                  pname(1:8) = 'p+      '


*-----------------------------------------------------------------------
*           'non' -> 'non     '
*-----------------------------------------------------------------------


            else if( pname(1:3) .eq. 'non' ) then


                  pname(1:8) = 'non     '


            end if


*-----------------------------------------------------------------------
*        Check the length of number in the front of the parameter
*-----------------------------------------------------------------------


               do ln = 0, 9


                  if( pname(1:1) .eq. num(ln) ) jnuc = 1


               end do


            if( jnuc .eq. 1 ) then


               do l  = 1, 3
               do ln = 0, 9


                  if( pname(l:l) .eq. num(ln) )
     &                inuc = inuc + 1


               end do
               end do


            end if


*-----------------------------------------------------------------------
*        Nucleus case, i.e.  start from number
*-----------------------------------------------------------------------


            if( inuc .ge. 1 ) then


                  read( pname(1:inuc), '(i12)' ) nmass


            end if


*-----------------------------------------------------------------------
*        hadron beam
*-----------------------------------------------------------------------


            if( inuc .eq. 0 ) then


               do j = 1, mxbeam


                  ll = lchcd(j)


                  if( pname(1:ll) .eq. chcde(j)(1:ll) ) then


                     idnum = kcode(j)


                  end if


               end do


               if( idnum .eq. 0 ) then


                  write(ieo,'('' Error: Unrecognized input'',
     &                        '' projectile or target particle'')')
                  write(ieo,'('' ====='')')


                  ierr = 1
                  return


               end if


               if( idnum .ne. -9999 ) then


                  idmas = 1
                  idchg = luchge(idnum) / 3


               end if


*-----------------------------------------------------------------------
*        Nucleon ; idnum = 0
*-----------------------------------------------------------------------


               if( idnum .eq. -9999 .or.
     &             idnum .eq.  2112 .or.
     &             idnum .eq.  2212 ) then


                  idnum = 0


               end if


*-----------------------------------------------------------------------
*        Nucleus ; idnum = 0
*-----------------------------------------------------------------------


            else if( inuc .gt. 0 ) then


                  cnuc = pname(inuc+1:len)


                  if( len - inuc .eq. 1 ) cnuc = cnuc(1:1)//'  '
                  if( len - inuc .eq. 2 ) cnuc = cnuc(1:2)//' '


               do j = 1, mxnucl


                  if( cnuc(1:3) .eq. elmnt(j)(1:3) ) then


                     nchag = j


                     goto 155


                  end if


               end do


                  write(ieo,'('' Error: Unrecognized input'',
     &                        '' projectile or target Nucleus'')')
                  write(ieo,'('' ====='')')


                  ierr = 1
                  return


  155          continue


                  idnum  = 0
                  idmas  = nmass
                  idchg  = nchag


            end if


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine chname(mnds0,mmas,mchg,chau,nchau)
*                                                                      *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to give the particle and Nucleus name                   *
*              as character string                                     *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              mnds0       : inds of the particle                      *
*              mmas        : mass of Nucleus                           *
*              mchg        : charge of Nucleus                         *
*              chau        : output of the name as character           *
*              nchau       : length of the name as character           *
*                                                                      *
*        Called subroutine and function :                              *
*                                                                      *
*              non                                                     *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      parameter ( mxpart =  11 )
      parameter ( mxnucl = 104 )


*-----------------------------------------------------------------------


      character       chau*16
      character       nucl*3
      character       nmas*4
      character       num(0:9)*1


      character       elmnt(mxnucl)*3
      dimension       nlmnt(mxnucl)


      character       chcde(0:mxpart)*12
      dimension       lchcd(0:mxpart)


*-----------------------------------------------------------------------


      data chcde / 'Gamma       ',
     &             'p           ','n           ',
     &             '\Delta^{++} ','\Delta^{+}  ',
     &             '\Delta^{0}  ','\Delta^{-}  ',
     &             'N^{*+}      ','N^{*0}      ',
     &             '\pi^{+}     ','\pi^{0}     ',
     &             '\pi^{-}     ' /


      data lchcd /              5,
     &                          1,            1,
     &                         10,            9,
     &                          9,            9,
     &                          6,            6,
     &                          7,            7,
     &                          7 /




      data elmnt /
     &    'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ',
     &    'Na ','Mg ','Al ','Si ','P  ','S  ','Cl ','Ar ','K  ','Ca ',
     &    'Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',
     &    'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ',
     &    'Nb ','Mo ','Te ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ',
     &    'Sb ','Te ','J  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',
     &    'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ',
     &    'Lu ','Hf ','Ta ','W  ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     &    'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',
     &    'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ',
     &    'Md ','No ','Lw ','Ku ' /


      data nlmnt /
     &    1,2,2,2,1,1,1,1,1,2,
     &    2,2,2,2,1,1,2,2,1,2,
     &    2,2,1,2,2,2,2,2,2,2,
     &    2,2,2,2,2,2,2,2,1,2,
     &    2,2,2,2,2,2,2,2,2,2,
     &    2,2,1,2,2,2,2,2,2,2,
     &    2,2,2,2,2,2,2,2,2,2,
     &    2,2,2,1,2,2,2,2,2,2,
     &    2,2,2,2,2,2,2,2,2,2,
     &    2,1,2,2,2,2,2,2,2,2,
     &    2,2,2,2/


      data num / '0','1','2','3','4','5','6','7','8','9' /


*-----------------------------------------------------------------------
*        id of the particle
*-----------------------------------------------------------------------


            if( mnds0 .eq. 0 ) then


               if( mmas .eq. 1 ) then


                  mnds = 1


               else


                  mnds = 0


               end if


            else


               if( mnds0 .eq.  211 .or.
     &             mnds0 .eq. -211 .or.
     &             mnds0 .eq.  111 ) then


                  mnds = 4


               else


                  mnds = mnds0


               end if


            end if


*-----------------------------------------------------------------------
*        Nucleus
*-----------------------------------------------------------------------


            if( mnds .eq. 0 ) then


               if( mchg .le. 104 ) then


                  nucl(1:3) = elmnt( mchg )
                  namel = nlmnt( mchg )


               else


                  i1   = mod( mchg,       10 )
                  i10  = mod( mchg / 10,  10 )
                  i100 = mod( mchg / 100, 10 )


                  nucl(1:3) = num(i100)//num(i10)//num(i1)
                  namel = 3


               end if


                  i1   = mod( mmas,       10 )
                  i10  = mod( mmas / 10,  10 )
                  i100 = mod( mmas / 100, 10 )


               if( mmas .lt. 10 ) then


                  nmas(1:1) = num(i1)
                  lmas      = 1


               else if( mmas .lt. 100 ) then


                  nmas(1:2) = num(i10)//num(i1)
                  lmas      = 2


               else


                  nmas(1:3) = num(i100)//num(i10)//num(i1)
                  lmas      = 3


               end if


                  chau  = '^{'//nmas(1:lmas)//'}'//nucl
                  nchau = lmas + 3 + namel




*-----------------------------------------------------------------------
*        Gamma
*-----------------------------------------------------------------------


            else if( mnds .eq. -1 ) then


                  chau  = chcde(0)(1:lchcd(0))
                  nchau = lchcd(0)


*-----------------------------------------------------------------------
*        Nucleon
*-----------------------------------------------------------------------


            else if( mnds .eq. 1 ) then


               if( mchg .eq. 1 ) then


                  chau  = chcde(1)(1:lchcd(1))
                  nchau = lchcd(1)


               else if( mchg .eq. 0 ) then


                  chau  = chcde(2)(1:lchcd(2))
                  nchau = lchcd(2)


               end if


*-----------------------------------------------------------------------
*        Delta
*-----------------------------------------------------------------------


            else if( mnds .eq. 2 ) then


               if( mchg .eq. 2 ) then


                  chau  = chcde(3)(1:lchcd(3))
                  nchau = lchcd(3)


               else if( mchg .eq. 1 ) then


                  chau  = chcde(4)(1:lchcd(4))
                  nchau = lchcd(4)


               else if( mchg .eq. 0 ) then


                  chau  = chcde(5)(1:lchcd(5))
                  nchau = lchcd(5)


               else if( mchg .eq. -1 ) then


                  chau  = chcde(6)(1:lchcd(6))
                  nchau = lchcd(6)


               end if






*-----------------------------------------------------------------------
*        N-star
*-----------------------------------------------------------------------


            else if( mnds .eq. 3 ) then


               if( mchg .eq. 1 ) then


                  chau  = chcde(7)(1:lchcd(7))
                  nchau = lchcd(7)


               else if( mchg .eq. 0 ) then


                  chau  = chcde(8)(1:lchcd(8))
                  nchau = lchcd(8)


               end if


*-----------------------------------------------------------------------
*        Pion
*-----------------------------------------------------------------------


            else if( mnds .eq. 4 ) then


               if( mchg .eq. 1 ) then


                  chau  = chcde(9)(1:lchcd(9))
                  nchau = lchcd(9)


               else if( mchg .eq. 0 ) then


                  chau  = chcde(10)(1:lchcd(10))
                  nchau = lchcd(10)


               else if( mchg .eq. -1 ) then


                  chau  = chcde(11)(1:lchcd(11))
                  nchau = lchcd(11)


               end if




            end if


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      function ulmass(kf)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 26                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to give the mass of the partilce with id = kf           *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              kf  : particle kf id                                    *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


            if( kf .eq. 0 ) then


                  ulmass = rmass


            else if( kf .eq. 2112 .or. kf .eq. 2212 ) then


                  ulmass = rmass


            else if( kf .eq.  211 .or.
     &               kf .eq. -211 .or.
     &               kf .eq.  111 ) then


                  ulmass = pmass


            else


                  ulmass = 0.0


            end if


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      function luchge(kf)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 26                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to give the charge of the partilce with id = kf         *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              kf  : particle kf id                                    *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


            if( kf .eq. 2212 .or. kf .eq. 211 ) then


                  luchge =  3


            else if( kf .eq. -211 ) then


                  luchge = -3


            else


                  luchge =  0


            end if


*-----------------------------------------------------------------------


      return
      end




