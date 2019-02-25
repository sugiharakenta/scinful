************************************************************************
*                                                                      *
*        PART 7: Output on file
*                                                                      *
*                                                                      *
*   List of subprograms in rough order of relevance with main purpose  *
*      ( s = subroutine, f = function, b = block data, e = entry )     *
*                                                                      *
*                                                                      *
*  s  sm_sumo   to write final output on file                          *
*  s  sm_hedo   to write header, input echo and summary                *
*  s  sm_coll   to sum up the collision number and                     *
*               to write collision history                             *
*  s  sm_grdo   to write ground state properties for masspr = 0        *
*  s  sm_cpuo   to write cputime on screen and/or file                 *
*  s  sm_init   to initiallize the summary and input echo              *
*  s  sm_evnt   to initialize the summary at each event                *
*  s  sm_timd   to output the phase space information for QMDDISP      *
*  s  sm_timc   to output the phase space information for summary      *
*  s  sm_timq   to store some values at the time steps for summary     *
*  s  sm_lstp   to count the last pion number for summary              *
*  s  sm_coln   to count the number of collisions for summary          *
*  s  sm_mdis   to write mass and charge distribution of QMD / SDM     *
*  s  sm_mass   to store the mass distribution of the clusters         *
*  s  sm_norm   to normarize the cross section by total event number   *
*  s  sm_clse   to close files for summary                             *
*                                                                      *
************************************************************************
*                                                                      *
************************************************************************
*                                                                      *
      subroutine sm_sumo
*                                                                      *
*                                                                      *
*        Last Revised:     1998 12 09                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to write final output on idsp(20)                       *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


      common /input1/ mstq1(mxpa1), parq1(mxpa1)


      common /verqmd/ versn, lastr, iyeav, imonv, idayv


      common /const0/ idnta, idnpr, massta, masspr, mstapr, msprpr
      common /const1/ elab, rdist, bmin, bmax, ibch, ibin
      common /const2/ dt, ntmax, iprun, iprun0
      common /const3/ nfreq, nfrec, nfred
      common /const4/ plab, srtcm, ylabb, pincm
      common /const5/ pzpr, pxpr, rzpr, rxpr
      common /const6/ pzta, pxta, rzta, rxta
      common /const7/ betpr, gampr, prmas, radpr
      common /const8/ betta, gamta, tamas, radta


      common /rannum/ iseed, iseed0, iseed1


      common /swich1/ ipot, insys, irkg, icolt
      common /swich2/ icfg, imany, icpus, idatm
      common /swich3/ ielst, jelst, kelst
      common /swich4/ ifin, ifout


      common /vriab0/ massal, massba, nmeson
      common /vriab1/ b, llnow, ntnow
      common /vriab2/ icevnt, idevnt


      common /poten0/ aaa, bbb, rpot, esymm
      common /poten1/ gamm, c0, c3, cs, cl, wl
      common /poten2/ t0, t3, rkk


      common /coultr/ eccm, pzcc, rmax0, zeroz


      common /coodrp/ r(3,nnn),  p(5,nnn)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)


      common /coln00/ lcoll(30)


      common /fnamid/ fname(mxfnm), idf(mxfnm), ifnl(mxfnm)
      character       fname*80
      common /inecho/ icnum, indata(200), indlng(200)
      character       indata*80


      common /startf/ iday0,imon0,iyer0,ihor0,imin0,isec0
      common /startt/ iday1,imon1,iyer1,ihor1,imin1,isec1
      common /cputim/ stime(30), cputm(30)


*-----------------------------------------------------------------------


      common /framtr/ betafr(0:2), gammfr(0:2)


      common /summas/ sumas(3,0:maxpt,0:maxnt)
      common /summch/ tmas(3,0:maxpt+maxnt), cdis(3,0:maxpt)


      common /sdmcut/ sdmemin


      common /sdmsw0/ issdm
      common /sdmsw1/ iswids, isevap, isfiss, isgrnd
      common /sdmsw2/ imengb, imbarr, imangm, imlevd, imgamm


      common /clustx/ tdecay(4)


*-----------------------------------------------------------------------


      common /summ00/ jdsp(50), idsp(50), fdsp(50), ifds(50)
      character    fdsp*80
      common /summ01/ iefw, ised, mulg, icdp(4), scald, isef
      common /summ02/ mdimg, pdimg, mminx, mmaxx, mminy, mmaxy,
     &                pminx, pmaxx, pminy, pmaxy, scalo, itdg
      common /summ03/ ireac(0:10), rcross(4)


      common /summ05/ npin(0:mmm), ndel(0:mmm), nres(0:mmm)
      common /summ06/ ictime, iltime
      common /summ07/ s_time(0:mmm)
      common /summ08/ s_ebin(0:mmm),s_epot(0:mmm),s_ekin(0:mmm),
     &                s_emas(0:mmm),s_epin(0:mmm),s_flow(0:mmm),
     &                s_rmsr(0:mmm),s_outp(0:mmm)
      common /summ09/ epott0, epotp0, ekint0, ekinp0, ebint0, ebinp0
      common /summ10/ adpir, anpir, anpid, ardpi, arnpi, adnpi
      common /summ11/ acoll, abloc, acnne, aeblc, aelpn, aelnd, aeldd,
     &                acnnd, acnnr, acndd, acndn, acnrn, acddn, acddd,
     &                acrdr, acdrd, acrrr, acdnd, acdnr, acrnd, acrnr
      common /summ12/ ncoll(30), nall(30,0:mmm)
      common /summ13/ denmom(0:mmm), denden(0:mmm), xyz(2,0:50),
     &                denpat(0:mmm)
      common /summ14/ dentim(0:mmm,0:50), denerr(2,0:50), dentm(0:50),
     &                denmtm(0:mmm,0:50), denmer(2,0:50)
      common /summ15/ radtac, jmax, facp
      common /summ16/ csw, ccw, kmax, imax


*-----------------------------------------------------------------------


      if( jdsp(20) .eq. 0 ) return


*-----------------------------------------------------------------------


               ida = idsp(20)


               inew = 0


*-----------------------------------------------------------------------
*     Header, Input Echo and Summary
*-----------------------------------------------------------------------


      if( jdsp(21) .ne. 0 ) then


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Header, Input Echo and Summary''/
     &          ''*'',71(''-''))')


         if( inew .gt. 0 ) then


            write(ida,'(/''newpage:'')')


         end if


            write(ida,'( /''infl: {'',80a1)')
     &           ( fdsp(21)(k:k), k = 1, ifds(21) ), '}'


            inew = inew + 1


      end if


*-----------------------------------------------------------------------
*     Collision History
*-----------------------------------------------------------------------


      if( jdsp(22) .ne. 0 ) then


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Collision History''/
     &          ''*'',71(''-''))')


         if( inew .gt. 0 ) then


            write(ida,'(/''newpage:'')')


         end if


            write(ida,'( /''infl: {'',80a1)')
     &           ( fdsp(22)(k:k), k = 1, ifds(22) ), '}'


            inew = inew + 1


      end if


*-----------------------------------------------------------------------
*     Ground State Properties
*-----------------------------------------------------------------------


      if( jdsp(23) .ne. 0 ) then


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Ground State Properties''/
     &          ''*'',71(''-''))')


         if( inew .gt. 0 ) then


            write(ida,'(/''newpage:'')')


         end if


            write(ida,'( /''infl: {'',80a1)')
     &           ( fdsp(23)(k:k), k = 1, ifds(23) ), '}'


            inew = inew + 1


      end if


*-----------------------------------------------------------------------
*     Plot of ANGEL display
*-----------------------------------------------------------------------


      if( jdsp(1) .ne. 0 ) then


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Snapshot of R-Space in Three directions''/
     &          ''*'',71(''-''))')


         if( inew .gt. 0 ) then


            write(ida,'(/''newpage:'')')


         end if


            write(ida,'( /''infl: {'',80a1)')
     &           ( fdsp(1)(k:k), k = 1, ifds(1) ), '}'


            inew = inew + 1


      end if


*-----------------------------------------------------------------------


      if( jdsp(2) .ne. 0 ) then


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Snapshot of P-Space in Three directions''/
     &          ''*'',71(''-''))')


         if( inew .gt. 0 ) then


            write(ida,'(/''newpage:'')')


         end if


            write(ida,'( /''infl: {'',80a1)')
     &           ( fdsp(2)(k:k), k = 1, ifds(2) ), '}'


            inew = inew + 1


      end if


*-----------------------------------------------------------------------


      if( jdsp(3) .ne. 0 ) then


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Snapshot of R and P-Space in Three directions''/
     &          ''*'',71(''-''))')


         if( inew .gt. 0 ) then


            write(ida,'(/''newpage:'')')


         end if


            write(ida,'( /''infl: {'',80a1)')
     &           ( fdsp(3)(k:k), k = 1, ifds(3) ), '}'


            inew = inew + 1


      end if


*-----------------------------------------------------------------------


      if( jdsp(4) .ne. 0 ) then


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Snapshot of R-Space''/
     &          ''*'',71(''-''))')


         if( inew .gt. 0 ) then


            write(ida,'(/''newpage:'')')


         end if


            write(ida,'( /''infl: {'',80a1)')
     &           ( fdsp(4)(k:k), k = 1, ifds(4) ), '}'


            inew = inew + 1


      end if


*-----------------------------------------------------------------------


      if( jdsp(5) .ne. 0 ) then


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Time Evolution of R-Space''/
     &          ''*'',71(''-''))')


         if( inew .gt. 0 ) then


            write(ida,'(/''newpage:'')')


         end if


            write(ida,'( /''infl: {'',80a1)')
     &           ( fdsp(5)(k:k), k = 1, ifds(5) ), '}'


            inew = inew + 1


      end if


*-----------------------------------------------------------------------


      if( jdsp(6) .ne. 0 ) then


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Snapshot of P-Space''/
     &          ''*'',71(''-''))')


         if( inew .gt. 0 ) then


            write(ida,'(/''newpage:'')')


         end if


            write(ida,'( /''infl: {'',80a1)')
     &           ( fdsp(6)(k:k), k = 1, ifds(6) ), '}'


            inew = inew + 1


      end if


*-----------------------------------------------------------------------


      if( jdsp(7) .ne. 0 ) then


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Time Evolution of P-Space''/
     &          ''*'',71(''-''))')


         if( inew .gt. 0 ) then


            write(ida,'(/''newpage:'')')


         end if


            write(ida,'( /''infl: {'',80a1)')
     &           ( fdsp(7)(k:k), k = 1, ifds(7) ), '}'


            inew = inew + 1


      end if


*-----------------------------------------------------------------------


      if( jdsp(8) .ne. 0 ) then


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Snapshot of R and P-Space''/
     &          ''*'',71(''-''))')


         if( inew .gt. 0 ) then


            write(ida,'(/''newpage:'')')


         end if


            write(ida,'( /''infl: {'',80a1)')
     &           ( fdsp(8)(k:k), k = 1, ifds(8) ), '}'


            inew = inew + 1


      end if


*-----------------------------------------------------------------------


      if( jdsp(9) .ne. 0 ) then


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Time Evolution of R and P-Space''/
     &          ''*'',71(''-''))')


         if( inew .gt. 0 ) then


            write(ida,'(/''newpage:'')')


         end if


            write(ida,'( /''infl: {'',80a1)')
     &           ( fdsp(9)(k:k), k = 1, ifds(9) ), '}'


            inew = inew + 1


      end if


*-----------------------------------------------------------------------


      if( jdsp(10) .ne. 0 ) then


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Snapshot R-Space with Color Plot''/
     &          ''*'',71(''-''))')


         if( inew .gt. 0 ) then


            write(ida,'(/''newpage:'')')


         end if


            write(ida,'( /''infl: {'',80a1)')
     &           ( fdsp(10)(k:k), k = 1, ifds(10) ), '}'


            inew = inew + 1


      end if


*-----------------------------------------------------------------------


      if( jdsp(11) .ne. 0 ) then


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Time Evolution of R-Space with Color Plot''/
     &          ''*'',71(''-''))')


         if( inew .gt. 0 ) then


            write(ida,'(/''newpage:'')')


         end if


            write(ida,'( /''infl: {'',80a1)')
     &           ( fdsp(11)(k:k), k = 1, ifds(11) ), '}'


            inew = inew + 1


      end if


*-----------------------------------------------------------------------


      if( jdsp(12) .ne. 0 ) then


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Snapshot R-Space with Contour and Color Plot''/
     &          ''*'',71(''-''))')


         if( inew .gt. 0 ) then


            write(ida,'(/''newpage:'')')


         end if


            write(ida,'( /''infl: {'',80a1)')
     &           ( fdsp(12)(k:k), k = 1, ifds(12) ), '}'


            inew = inew + 1


      end if


*-----------------------------------------------------------------------


      if( jdsp(13) .ne. 0 ) then


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Time Evolution of R-Space'',
     &               '' with Contour and Color Plot''/
     &          ''*'',71(''-''))')


         if( inew .gt. 0 ) then


            write(ida,'(/''newpage:'')')


         end if


            write(ida,'( /''infl: {'',80a1)')
     &           ( fdsp(13)(k:k), k = 1, ifds(13) ), '}'


            inew = inew + 1


      end if


*-----------------------------------------------------------------------
*     Mass Distribution of QMD
*-----------------------------------------------------------------------


      if( jdsp(30) .ne. 0 ) then


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Mass Distribution of QMD''/
     &          ''*'',71(''-''))')


         if( inew .gt. 0 ) then


            write(ida,'(/''newpage:'')')


         end if


            write(ida,'( /''infl: {'',80a1)')
     &           ( fdsp(30)(k:k), k = 1, ifds(30) ), '}'


            inew = inew + 1


      end if


*-----------------------------------------------------------------------
*     Mass Distribution of SDM
*-----------------------------------------------------------------------


      if( jdsp(31) .ne. 0 ) then


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Mass Distribution of SDM''/
     &          ''*'',71(''-''))')


         if( inew .gt. 0 ) then


            write(ida,'(/''newpage:'')')


         end if


            write(ida,'( /''infl: {'',80a1)')
     &           ( fdsp(31)(k:k), k = 1, ifds(31) ), '}'


            inew = inew + 1


      end if


*-----------------------------------------------------------------------
*     Mass Distribution of fission
*-----------------------------------------------------------------------


      if( jdsp(32) .ne. 0 ) then


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Mass Distribution of fission''/
     &          ''*'',71(''-''))')


         if( inew .gt. 0 ) then


            write(ida,'(/''newpage:'')')


         end if


            write(ida,'( /''infl: {'',80a1)')
     &           ( fdsp(32)(k:k), k = 1, ifds(32) ), '}'


            inew = inew + 1


      end if


*-----------------------------------------------------------------------




      return
      end




************************************************************************
*                                                                      *
      subroutine sm_hedo
*                                                                      *
*                                                                      *
*        Last Revised:     1998 12 09                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*             to write header, input echo and summary                  *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


      common /verqmd/ versn, lastr, iyeav, imonv, idayv


      common /const0/ idnta, idnpr, massta, masspr, mstapr, msprpr
      common /const1/ elab, rdist, bmin, bmax, ibch, ibin
      common /const2/ dt, ntmax, iprun, iprun0
      common /const4/ plab, srtcm, ylabb, pincm
      common /const5/ pzpr, pxpr, rzpr, rxpr
      common /const6/ pzta, pxta, rzta, rxta
      common /const7/ betpr, gampr, prmas, radpr
      common /const8/ betta, gamta, tamas, radta


      common /rannum/ iseed, iseed0, iseed1


      common /swich1/ ipot, insys, irkg, icolt
      common /swich2/ icfg, imany, icpus, idatm
      common /swich3/ ielst, jelst, kelst
      common /swich4/ ifin, ifout


      common /vriab2/ icevnt, idevnt


      common /poten0/ aaa, bbb, rpot, esymm
      common /poten1/ gamm, c0, c3, cs, cl, wl
      common /poten2/ t0, t3, rkk


      common /inecho/ icnum, indata(200), indlng(200)
      character       indata*80


      common /startf/ iday0,imon0,iyer0,ihor0,imin0,isec0
      common /startt/ iday1,imon1,iyer1,ihor1,imin1,isec1
      common /cputim/ stime(30), cputm(30)


*-----------------------------------------------------------------------


      common /sdmcut/ sdmemin


      common /sdmsw0/ issdm
      common /sdmsw1/ iswids, isevap, isfiss, isgrnd
      common /sdmsw2/ imengb, imbarr, imangm, imlevd, imgamm


      common /clustx/ tdecay(4)


*-----------------------------------------------------------------------


      common /summ00/ jdsp(50), idsp(50), fdsp(50), ifds(50)
      character    fdsp*80


      common /summ03/ ireac(0:10), rcross(4)


      common /summ05/ npin(0:mmm), ndel(0:mmm), nres(0:mmm)
      common /summ06/ ictime, iltime
      common /summ07/ s_time(0:mmm)
      common /summ08/ s_ebin(0:mmm),s_epot(0:mmm),s_ekin(0:mmm),
     &                s_emas(0:mmm),s_epin(0:mmm),s_flow(0:mmm),
     &                s_rmsr(0:mmm),s_outp(0:mmm)
      common /summ09/ epott0, epotp0, ekint0, ekinp0, ebint0, ebinp0
      common /summ10/ adpir, anpir, anpid, ardpi, arnpi, adnpi
      common /summ11/ acoll, abloc, acnne, aeblc, aelpn, aelnd, aeldd,
     &                acnnd, acnnr, acndd, acndn, acnrn, acddn, acddd,
     &                acrdr, acdrd, acrrr, acdnd, acdnr, acrnd, acrnr


*-----------------------------------------------------------------------


      character       chta*16, chpr*16


*-----------------------------------------------------------------------
*     Header, Input Echo and Summary
*-----------------------------------------------------------------------


      if( jdsp(21) .eq. 0 ) return


*-----------------------------------------------------------------------


            ida = idsp(21)


*-----------------------------------------------------------------------


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Header, Input Echo and Summary'' /
     &          ''*'',71(''-'') //
     &          ''msul:{\rm\it JQMD  plotted by  \ANGEL}''/
     &          ''msdc:{\rm\Large \- \page  \-}''/
     &          ''msur:{\rm\it Header}'')')


            write(ida,'(/
     &          ''p: secp mssg notf'')')


*-----------------------------------------------------------------------
*     Title of JQMD
*-----------------------------------------------------------------------


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Title of JQMD'' /
     &          ''*'',71(''-''))')


            write(ida,'(/
     &          ''wt: x(8) y(24.5) ix(2) iy(3) f(5) s(0.6)'',
     &          '' box(OvalshadowBox)''/
     &          ''\vspace{1.5}{\Huge\Huge\Huge\hv\bf   JQMD  }''/
     &          ''\vspace{1.5}'',
     &          ''{\color{r}\LARGE\bf J}aeri '',
     &          ''{\color{r}\LARGE\bf Q}uantum \ ''/
     &          ''{\color{r}\LARGE\bf M}olecular '',
     &          ''{\color{r}\LARGE\bf D}ynamics''/
     &          ''for Nuclear Reactions''/
     &          ''Version ='',f5.2,/
     &          ''made by''/
     &          ''Japan Atomic Energy Research Institute''/
     &          ''Last Revised  '',i4,1x,i2.2,1x,i2.2/
     &          ''e:'')')
     &          versn, iyeav,imonv, idayv


*-----------------------------------------------------------------------
*     Box for parameters
*-----------------------------------------------------------------------


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Box for Parameters'' /
     &          ''*'',71(''-''))')


            write(ida,'(/
     &          ''box: x(9.2) y(7.0) s(0.6) sx(27) sy(40)'',
     &          '' box(doublebox)'')')


*-----------------------------------------------------------------------
*     System and Energy
*-----------------------------------------------------------------------


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     System and Energy'' /
     &          ''*'',71(''-''))')


            write(ida,'(/
     &          ''wt: x(1) y(17.5) ix(1) iy(3) f(5) s(0.7)'',
     &          '' box(shadowBox)'')')


*-----------------------------------------------------------------------
*        normal reactions
*-----------------------------------------------------------------------


         if( masspr .gt. 0 ) then


                  call chname(idnta,massta,mstapr,chta,nchta)
                  call chname(idnpr,masspr,msprpr,chpr,nchpr)


               write(ida,'(
     &               ''\vspace{0.5} \ ''/
     &               '' & {\hv\bf\Huge{\Large '',20a1)')
     &               (chpr(k:k),k=1,nchpr),'}','\ '
               write(ida,'(
     &               ''   on  {\Large '',20a1)')
     &               (chta(k:k),k=1,nchta),'}',' ','}'


               write(ida,'(''\vspace{1.0}\ '')')


            if( elab .gt. 1.0 ) then


               write(ida,'(
     &               '' E_{Lab} & = '',f8.2,'' A GeV,    '',
     &               '' P_{Lab}   = '',f8.2,'' A GeV/c''/
     &               '' E^{nn}_{cm}  & = '',f8.2,'' GeV'')')
     &               elab, plab, srtcm - prmas - tamas


            else


               write(ida,'(
     &               '' E_{Lab} & = '',f8.2,'' A MeV,    '',
     &               '' P_{Lab}   = '',f8.2,'' A MeV/c''/
     &               '' E^{nn}_{cm}  & = '',f8.2,'' MeV'')')
     &               elab*1000., plab*1000.,
     &               ( srtcm - prmas - tamas ) * 1000.


            end if


               write(ida,'(''\vspace{-0.8}'')')


*-----------------------------------------------------------------------
*        ground state properties
*-----------------------------------------------------------------------


         else


                  call chname(idnta,massta,mstapr,chta,nchta)


               write(ida,'(
     &               ''\vspace{1.0}\hspace{2.0}\ ''/
     &               ''{\hv\bf\Huge\Large '',20a1)')
     &               (chta(k:k),k=1,nchta),'}'


               write(ida,'(''\vspace{1.0}\ '')')


               write(ida,'(
     &               ''Ground State Properties''/)')


         end if


               write(ida,'(''e:'')')


*-----------------------------------------------------------------------
*     Write parameters and CPU times
*-----------------------------------------------------------------------


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Parameters and CPU time'' /
     &          ''*'',71(''-''))')


*-----------------------------------------------------------------------
*        Reference Frame
*-----------------------------------------------------------------------


         if( insys .eq. 1 ) then


            write(ida,'(/''w:Calculated in C.M. frame :'',
     &                 ''/x(5) y(13.2) s(0.6) f(5)'')')


         else if( insys .eq. 0 ) then


            write(ida,'(/''w:Calculated in Lab. frame :'',
     &                 ''/x(5) y(13.2) s(0.6) f(5)'')')


         else if( insys .eq. 2 ) then


            write(ida,'(/''w:Calculated in N-N C.M. frame :'',
     &                 ''/x(5) y(13.2) s(0.6) f(5)'')')


         end if


*-----------------------------------------------------------------------
*        Incident velocity, gamma factor, momentum and position
*-----------------------------------------------------------------------


            write(ida,'(''wtab: x(6) y(13.0) ix(1) iy(3) '',
     &                 ''f(5) s(0.6) tab{lcrr}'')')


            write(ida,'(
     &            ''              &   & '',
     &            ''     Beam  &     Target'')')


            write(ida,'(
     &            '' Velocity / c & : & '',f9.7,''  &  '',f9.7)')
     &            betpr, betta


         if( elab .gt. 1.0 ) then


            write(ida,'(
     &            '' Gamma factor & : & '',f9.3,''  &  '',f9.3)')
     &            gampr, gamta


            write(ida,'(
     &            '' p_z  (GeV/c) & : & ''f9.3,''  &  '',f9.3)')
     &            pzpr, pzta


         else


            write(ida,'(
     &            '' Gamma factor & : & '',f9.5,''  &  '',f9.5)')
     &            gampr, gamta


            write(ida,'(
     &            '' p_z  (MeV/c) & : & '',f9.3,''  &  '',f9.3)')
     &            pzpr*1000., pzta*1000.


         end if


            write(ida,'(
     &            '' r_z  (fm)    & : & '',f9.3,''  &  '',f9.3)')
     &            rzpr, rzta


            write(ida,'(''e:'')')


*-----------------------------------------------------------------------
*        Impact Parameter
*-----------------------------------------------------------------------


            write(ida,'(/''w:Impact Parameter Range :'',
     &                 ''/x(5) y(10.0) s(0.6) f(5)'')')


            write(ida,'(''wtab: x(6) y(9.8) ix(1) iy(3) '',
     &                 ''f(5) s(0.6) tab{lcccl}'')')


            write(ida,'(
     &            f7.3,'' & < & b & < & ''f7.3,'' (fm)'')')
     &            bmin, bmax


            write(ida,'(''e:'')')


*-----------------------------------------------------------------------
*        Event number, time step and seed of rumdom number
*-----------------------------------------------------------------------


         if( iprun .eq. iprun0 ) then


            write(ida,'(/
     &            ''w:Time Evolution :  {\color{b}'',
     &            ''[ Calculation is Finished Normally ]}\ ''/
     &            ''/x(5) y(8.8) s(0.6) f(5)'')')


         else


            write(ida,'(/
     &            ''w:Time Evolution :  {\color{r}'',
     &            ''[ Run is Terminated at Event = '',i6,'' ]}\ ''/
     &            ''/x(5) y(8.8) s(0.6) f(5)'')') iprun


         end if


            write(ida,'(
     &            ''wtab: x(6) y(8.6) ix(1) iy(3) f(5)'',
     &            '' s(0.6) tab{lcl}'')')


            write(ida,'(
     &            '' Input Number of Events & = & '',i9)') iprun0


            write(ida,'(
     &            '' Number of Time Step    & = & '',i9)') ntmax


            write(ida,'(
     &            '' Time step (fm/c)       & = & '',f9.3)') dt


            write(ida,'(
     &            '' Seed of random number  & = & '',i12)') iseed0


            write(ida,'(''e:'')')


*-----------------------------------------------------------------------
*        Date and elaspe time
*-----------------------------------------------------------------------


            write(ida,'(/
     &            ''w:Date and Elaspe time:'',
     &            ''/x(5) y(6.1) s(0.6) f(5)'')')


            write(ida,'(
     &      ''wtab: x(6) y(5.9) ix(1) iy(3) f(5) s(0.6) tab{rcrr}'')')


         if( idatm .eq. 1 ) then


            write(ida,'('' Starting time & = & '',
     &                 i4.4,''-'',i2.2,''-'',i2.2,''    & '',
     &                 i2.2,'':'',i2.2,'':'',i2.2)')
     &                 iyer0,imon0,iday0,ihor0,imin0,isec0


            write(ida,'(''   Ending time & = & '',
     &                 i4.4,''-'',i2.2,''-'',i2.2,''    & '',
     &                 i2.2,'':'',i2.2,'':'',i2.2)')
     &                 iyer1,imon1,iday1,ihor1,imin1,isec1


         else


            write(ida,'('' Starting time & = & not detected &'')')


            write(ida,'(''   Ending time & = & not detected &'')')


         end if


         if( icpus .eq. 1 ) then


               ilpsh = int( stime(1) / 3600.0 )
               ilpsm = int( ( stime(1) - ilpsh * 3600.0 ) / 60.0 )
               elpss = stime(1) - ilpsh * 3600.0 - ilpsm * 60.0


            write(ida,'(''   Elapse time & = & '',i4,'' h '',
     &                  i4,'' m &  '',f6.2,'' sec'')')
     &                  ilpsh, ilpsm, elpss


         else


            write(ida,'(''   Elapse time & = & not detected &'')')


         end if


            write(ida,'(''e:'')')


*-----------------------------------------------------------------------
*        CPU time
*-----------------------------------------------------------------------


         if( icpus .eq. 1 ) then


            write(ida,'(/
     &            ''w:CPU time:/x(5) y(3.9) s(0.6) f(5)'')')


            write(ida,'(
     &      ''wtab: x(6) y(3.7) ix(1) iy(3) f(5) s(0.6) tab{rcrcrl}'')')


                  ipcp2 = int(cputm(2)/cputm(1)*100.0)
                  ipcp3 = int(cputm(3)/cputm(1)*100.0)
                  ipcp4 = int(cputm(4)/cputm(1)*100.0)
                  ipcp5 = int(cputm(5)/cputm(1)*100.0)
                  ipcp6 = int(cputm(6)/cputm(1)*100.0)


                  ipcpr = 100 - ipcp2 - ipcp3
     &                        - ipcp4 - ipcp5 - ipcp6
                  cputr = cputm(1) - cputm(2) - cputm(3)
     &                  - cputm(4) - cputm(5) - cputm(6)


            write(ida,'(''Total CPU time & = &        &     &'',
     &                  f12.2,'' & sec'')') cputm(1)


            write(ida,'(''    Mean Field & = &    '',
     &                  i3,'' & % : &'',f12.2,'' & sec'')')
     &                  ipcp2, cputm(2)
            write(ida,'(''     Collision & = &    '',
     &                  i3,'' & % : &'',f12.2,'' & sec'')')
     &                  ipcp3, cputm(3)
            write(ida,'(''  Ground state & = &    '',
     &                  i3,'' & % : &'',f12.2,'' & sec'')')
     &                  ipcp4, cputm(4)
            write(ida,'(''       Summary & = &    '',
     &                  i3,'' & % : &'',f12.2,'' & sec'')')
     &                  ipcp5, cputm(5)
            write(ida,'(''  Final Decays & = &    '',
     &                  i3,'' & % : &'',f12.2,'' & sec'')')
     &                  ipcp6, cputm(6)
            write(ida,'(''        Others & = &    '',
     &                  i3,'' & % : &'',f12.2,'' & sec'')')
     &                  ipcpr, cputr


         else


            write(ida,'(/
     &            ''w:CPU time:  not detected/'',
     &            ''x(5) y(3.9) s(0.6) f(5)'')')


         end if


            write(ida,'(''e:'')')


*-----------------------------------------------------------------------
*     Second Page
*-----------------------------------------------------------------------


            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Second Page''/
     &            ''*'',71(''-''))')


            write(ida,'(/
     &            ''newpage:'')')


            write(ida,'(/
     &            ''msul:{\rm\it JQMD  plotted by  \ANGEL}''/
     &            ''msdc:{\rm\Large \- \page  \-}''/
     &            ''msur:{\rm\it Input Echo and Summary}'')')


            write(ida,'(/
     &            ''p: secp mssg notf'')')


*-----------------------------------------------------------------------
*     Input Echo
*-----------------------------------------------------------------------


            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Input Echo''/
     &            ''*'',71(''-''))')


            write(ida,'(/
     &            ''w:Input Echo/x(1.5) y(23.5) ix(2) s(0.8) f(5)'')')


            write(ida,'(
     &            ''wtab: x(1.5) y(23.2) ix(2) iy(3) f(5) s(0.6)'',
     &            '' tab{|rcl|}''/''\hline''/)')


      if( icnum .gt. 0 ) then


         if( icnum .gt. 46 ) then


            icnm = 46


         else


            icnm = icnum


         end if


         do i = 1, icnm


            write(ida,'('' '',80a1)') ( indata(i)(j:j), j=1,indlng(i) )


         end do


         if( icnum .gt. 46 ) then


            write(ida,'('' .......... & & '')')
            write(ida,'('' .......... & & more input data'')')


         end if


      else


            write(ida,'('' Input Parameters & are given in & Main'')')


      end if


            write(ida,'(/''\hline''/''e:'')')


*-----------------------------------------------------------------------
*        Mean Field Parameters
*-----------------------------------------------------------------------


            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Mean Field Parameters''/
     &            ''*'',71(''-''))')


            write(ida,'(/
     &       ''w:Mean Field'',
     &       ''/x(6.5) y(23.5) ix(2) s(0.8) f(5)'')')


            write(ida,'(
     &       ''wtab: x(6.5) y(23.2) ix(2) iy(3) f(5) s(0.6)'',
     &       '' tab{|lcrl|}''/''\hline''/)')


            write(ida,'('' A     & = & '',f8.2,'' & (MeV)'')') aaa*1000.
            write(ida,'('' B     & = & '',f8.2,'' & (MeV)'')') bbb*1000.
            write(ida,'('' \tau  & = & '',f8.4,'' & ''     )') gamm
            write(ida,'('' C_s   & = & '',f8.2,'' & (MeV)'')')
     &                                                       esymm*1000.
            write(ida,'('' K     & = & '',f8.2,'' & (MeV)'')') rkk*1000.
            write(ida,'('' L     & = & '',f8.2,'' & (fm^2)'')') wl


            write(ida,'(/''\hline''/''e:'')')


*-----------------------------------------------------------------------
*     Reading Data from Files
*-----------------------------------------------------------------------


      if( ifin .gt. 0 ) then


            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Reading Data from Files''/
     &            ''*'',71(''-''))')


            write(ida,'(/
     &         ''w:Data from Files\ ''/
     &         ''/x(13.0) y(23.5) ix(2) s(0.8) f(5)'')')


            write(ida,'(
     &         ''wtab: x(13.0) y(23.2) ix(2) iy(3) f(5) s(0.6) '',
     &         ''tab{|l|c|}''/''\hline'')')


            write(ida,'(
     &         ''                         & # of Events''/
     &         ''\hline''/
     &         '' Input Number            & '',i9/
     &         '' All in Files            & '',i9/
     &         '' All within b-range      & '',i9/
     &         '' Selected within b-range & '',i9)')
     &         iprun0, icevnt, idevnt, iprun


            write(ida,'(/''\hline''/''e:'')')


*-----------------------------------------------------------------------
*     QMD calculation
*-----------------------------------------------------------------------


      else if( ifin .eq. 0 ) then


*-----------------------------------------------------------------------
*        Average ground state energy
*-----------------------------------------------------------------------


            rmn  = float(iprun)


            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Average ground state energy''/
     &            ''*'',71(''-''))')


            write(ida,'(/
     &         ''w:Ground State Energy\ ''/
     &         ''/x(13.0) y(23.5) ix(2) s(0.8) f(5)'')')


            write(ida,'(
     &         ''wtab: x(13.0) y(23.2) ix(2) iy(3) f(5) s(0.6) '',
     &         ''tab{|c|r|r|}''/''\hline'')')


              write(ida,'(
     &         '' Energy (MeV/A) &     Target  &    Projectile''/
     &         ''\hline\hline''/
     &         '' Potential      & '',f10.3,'' & '',f10.3/
     &         '' Kinetic        & '',f10.3,'' & '',f10.3/
     &         '' Binding        & '',f10.3,'' & '',f10.3)')
     &           epott0*1000./rmn, epotp0*1000./rmn,
     &           ekint0*1000./rmn, ekinp0*1000./rmn,
     &           ebint0*1000./rmn, ebinp0*1000./rmn


            write(ida,'(/''\hline''/''e:'')')


*-----------------------------------------------------------------------
*        Check of total energy consevation
*-----------------------------------------------------------------------


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Check of total energy consevation''/
     &          ''*'',71(''-''))')


            write(ida,'(/''c: plot of total time evolution'')')
            write(ida,'( ''c: Time (fm/c)'')')
            write(ida,'( ''c: Energy/u (MeV)'')')
            write(ida,'(
     &          ''c:   t       edif        ekin        epot    '',
     &          ''    emas        epion'')')


            write(ida,'(
     &          ''n:   x     y(edif),l0  y(kin),l0r  y(pot),l0b'',
     &          ''  y(mas),l0c  y(pin),l0r'')')


            write(ida,'(f7.2,1x,5e12.4)') 0.0, 0.0,
     &                                    s_ekin(0)*1000./rmn,
     &                                    s_epot(0)*1000./rmn,
     &                                    s_emas(0)*1000./rmn,
     &                                    s_epin(0)*1000./rmn


            write(ida,'(f7.2,1x,5e12.4)')
     &      ( s_time(i),
     &       (s_ebin(i)+s_emas(i)+s_epin(i)
     &       -s_ebin(0)-s_emas(0)-s_epin(0))*1000./rmn,
     &        s_ekin(i)*1000./rmn,
     &        s_epot(i)*1000./rmn,
     &        s_emas(i)*1000./rmn,
     &        s_epin(i)*1000./rmn,  i = 1, ictime )


            write(ida,'(/
     &         ''w:Energy Conservation'',
     &         ''/x(13.0) y(20.3) ix(2) s(0.8) f(5)'')')


            write(ida,'(
     &         ''w:{\singleBox E_{diff}   = '',1p1e12.4,
     &         ''   (MeV/A)\ ''/
     &         ''/x(13.0) y(19.5) ix(2) f(5) s(0.6)'')')
     &         (s_ebin(ictime)+s_emas(ictime)+s_epin(ictime)
     &         -s_ebin(0)-s_emas(0)-s_epin(0))*1000./rmn




*-----------------------------------------------------------------------
*        Number of Collisions
*-----------------------------------------------------------------------


         if( masspr .gt. 0 ) then


*-----------------------------------------------------------------------


            rmn  = float(iprun)


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Number of Collisions''/
     &          ''*'',71(''-''))')


            write(ida,'(/
     &         ''w:Number of Collisions\ ''/
     &         ''/x(13.0) y(18.3) ix(2) s(0.8) f(5)'')')


            write(ida,'(
     &         ''wtab: x(13.0) y(18.0) ix(2) iy(3) f(5) s(0.6) '',
     &         ''tab{|c|r|c|r|}''/''\hline'')')


            write(ida,'('' All Collisions   & '',f9.3,'' & '',
     &                  '' Blocked          & '',f9.3)')
     &                  acoll/rmn, abloc/rmn + aeblc/rmn


            write(ida,'('' Inelastic        & '',f9.3,'' & '',
     &                  '' Elastic          & '',f9.3)')
     &                   acoll/rmn - abloc/rmn - aeblc/rmn - acnne/rmn,
     &                   acnne/rmn


            write(ida,'(''\hline\hline'')')


            write(ida,'('' NN \to  ND       & '',f9.3,'' & '',
     &                  '' ND \to  NN       & '',f9.3)')
     &                  acnnd/rmn, acndn/rmn


            write(ida,'('' NN \to  NR       & '',f9.3,'' & '',
     &                  '' NR \to  NN       & '',f9.3)')
     &                  acnnr/rmn, acnrn/rmn


            write(ida,'('' NN \to  DD       & '',f9.3,'' & '',
     &                  '' DD \to  NN       & '',f9.3)')
     &                  acndd/rmn, acddn/rmn


            write(ida,'(''\hline\hline'')')


            write(ida,'('' ND \to  DD       & '',f9.3,'' & '',
     &                  '' DD \to  ND       & '',f9.3)')
     &                  acddd/rmn, acdnd/rmn


            write(ida,'('' NR \to  DR       & '',f9.3,'' & '',
     &                  '' DR \to  NR       & '',f9.3)')
     &                  acrdr/rmn, acdnr/rmn


            write(ida,'('' ND \to  RD       & '',f9.3,'' & '',
     &                  '' RD \to  ND       & '',f9.3)')
     &                  acdrd/rmn, acrnd/rmn


            write(ida,'('' NR \to  RR       & '',f9.3,'' & '',
     &                  '' RR \to  NR       & '',f9.3)')
     &                  acrrr/rmn, acrnr/rmn


            write(ida,'(''\hline\hline'')')


            write(ida,'('' D \to  N \+ \pi  & '',f9.3,'' & '',
     &                  '' N \+ \pi  \to  D & '',f9.3)')
     &                  adnpi/rmn, anpid/rmn


            write(ida,'('' R \to  N \+ \pi  & '',f9.3,'' & '',
     &                  '' N \+ \pi  \to  R & '',f9.3)')
     &                  arnpi/rmn, anpir/rmn


            write(ida,'('' R \to  D \+ \pi  & '',f9.3,'' & '',
     &                  '' D \+ \pi  \to  R & '',f9.3)')
     &                  ardpi/rmn, adpir/rmn


            write(ida,'(''\hline\hline'')')


            write(ida,'('' Final \pi        & '',f9.3,'' & & '')')
     &                  float(npin(ictime+1))/rmn


            write(ida,'(''\hline''/''e:'')')


*-----------------------------------------------------------------------


         end if


*-----------------------------------------------------------------------


      end if


*-----------------------------------------------------------------------
*     Reaction Cross Section
*-----------------------------------------------------------------------


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Reaction Cross Section''/
     &          ''*'',71(''-''))')


            write(ida,'(/
     &         ''w:Reaction Cross Section\ ''/
     &         ''/x(13.0) y(11.0) ix(2) s(0.8) f(5)'')')


            write(ida,'(
     &         ''wtab: x(13.0) y(10.7) ix(2) iy(3) f(5) s(0.6) '',
     &         ''tab{|l|r|}''/''\hline'')')


            write(ida,'('' Reaction Type'',11x''& # of Events'')')
            write(ida,'(''\hline'')')
            write(ida,'('' Elastic. without Coll.  & '',i9)') ireac(0)


            if( ielst .eq. 1 ) write(ida,'(''\hline\hline'')')


            write(ida,'('' Elastic. with Collision & '',i9)') ireac(1)


            if( ielst .eq. 2 ) write(ida,'(''\hline\hline'')')


            write(ida,'('' Inelast. without Coll.  & '',i9)') ireac(2)


            if( ielst .eq. 3 ) write(ida,'(''\hline\hline'')')


            write(ida,'('' Inelast. with Collision & '',i9)') ireac(3)




            write(ida,'(''\hline\hline'')')


            write(ida,'('' Input total b-range (mb) &   '',f8.2)')
     &                     rcross(1)


            write(ida,'('' Beam X-section (mb)      &   '',f8.2)')
     &                     rcross(2)


            write(ida,'('' Reaction X-section (mb)  &   '',f8.2)')
     &                     rcross(3)


            write(ida,'(''\hline''/''e:'')')




               nelst = 0


            do issc = 0, ielst - 1


               nelst = nelst + ireac(issc)


            end do


               rinel = max( 1.0, float( iprun - nelst ) )


*-----------------------------------------------------------------------
*     Parameters for SDM
*-----------------------------------------------------------------------


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Parameters for SDM''/
     &          ''*'',71(''-''))')


*-----------------------------------------------------------------------
      if( issdm .gt. 0 ) then
*-----------------------------------------------------------------------


            write(ida,'(/
     &         ''w:Statistical Decay is included\ ''/
     &         ''/x(13.0) y(6.0) ix(2) s(0.8) f(5)'')')


            write(ida,'(
     &         ''wtab: x(13.0) y(5.7) ix(2) iy(3) f(5) s(0.6) '',
     &         ''tab{|r|c|}''/''\hline'')')


            write(ida,'('' SDM / QMD        & '',i5)') issdm




         if( iswids .eq. 0 ) then


            write(ida,'('' decay width      & without angular mom.'')')


         else


            write(ida,'('' decay width      & with angular mom.'')')


         end if


         if( isevap .eq. 1 ) then


            write(ida,'('' particle decay   &'',
     &                  '' n, p, d, t, ^3He, \alpha'')')


         else


            write(ida,'('' particle decay   &'',
     &                  '' not included'')')


         end if


         if( iswids .eq. 1 .and. imgamm .eq. 1 ) then


            write(ida,'('' gamma decay      & included'')')


         else


            write(ida,'('' gamma decay      & not included'')')


         end if


         if( isgrnd .eq. 1 ) then


            write(ida,'('' ground decay     & included'')')


         else


            write(ida,'('' ground decay     & not included'')')


         end if


         if( isfiss .eq. 1 ) then


            write(ida,'('' fission          & included'')')


         else


            write(ida,'('' fission          & not included'')')


         end if


            write(ida,'('' cutoff energy    & '',f6.2,'' MeV'')')
     &                  sdmemin


*-----------------------------------------------------------------------


         if( iswids .eq. 1 ) then


            write(ida,'('' # of energy bin  & '',i3)') imengb


            if( imbarr .eq. 1 ) then


               write(ida,'('' barrier       & modified'')')


            else


               write(ida,'('' barrier       & simple'')')


            end if


            if( imangm .eq. 1 ) then


               write(ida,'('' angular mom.  & included'')')


            else


               write(ida,'('' angular mom.  & not included'')')


            end if


            if( imlevd .eq. 1 ) then


               write(ida,'('' level density & with angular mom.'')')


            else


               write(ida,'('' level density & without angular mom.'')')


            end if


         end if


*-----------------------------------------------------------------------


            write(ida,'(''\hline\hline'')')


            write(ida,'('' Decay Mode       & # / inel. events'')')


            write(ida,'(''\hline'')')


            write(ida,'('' particle decay   & '',f9.3)')
     &                  tdecay(1) / rinel


            write(ida,'('' gamma decay      & '',f9.3)')
     &                  tdecay(2) / rinel


            write(ida,'('' ground decay     & '',f9.3)')
     &                  tdecay(3) / rinel


            write(ida,'('' fission          & '',f9.3)')
     &                  tdecay(4) / rinel


*-----------------------------------------------------------------------


            write(ida,'(''\hline''/''e:'')')


*-----------------------------------------------------------------------
       else
*-----------------------------------------------------------------------


            write(ida,'(/
     &         ''w:Statistical Decay is not included\ ''/
     &         ''/x(13.0) y(6.0) ix(2) s(0.8) f(5)'')')




*-----------------------------------------------------------------------
       end if
*-----------------------------------------------------------------------




      return
      end




************************************************************************
*                                                                      *
      subroutine sm_coll
*                                                                      *
*                                                                      *
*        Last Revised:     1998 12 09                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*             to sum up the collision number                           *
*             and write collision history on file                      *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


      common /const0/ idnta, idnpr, massta, masspr, mstapr, msprpr
      common /const2/ dt, ntmax, iprun, iprun0


      common /summ00/ jdsp(50), idsp(50), fdsp(50), ifds(50)
      character    fdsp*80


      common /summ05/ npin(0:mmm), ndel(0:mmm), nres(0:mmm)
      common /summ06/ ictime, iltime
      common /summ07/ s_time(0:mmm)
      common /summ08/ s_ebin(0:mmm),s_epot(0:mmm),s_ekin(0:mmm),
     &                s_emas(0:mmm),s_epin(0:mmm),s_flow(0:mmm),
     &                s_rmsr(0:mmm),s_outp(0:mmm)
      common /summ10/ adpir, anpir, anpid, ardpi, arnpi, adnpi
      common /summ11/ acoll, abloc, acnne, aeblc, aelpn, aelnd, aeldd,
     &                acnnd, acnnr, acndd, acndn, acnrn, acddn, acddd,
     &                acrdr, acdrd, acrrr, acdnd, acdnr, acrnd, acrnr
      common /summ12/ ncoll(30), nall(30,0:mmm)


*-----------------------------------------------------------------------


      character       squ*1
      data squ       /"'"/


*-----------------------------------------------------------------------


      if( jdsp(20) .eq. 0 ) return


*-----------------------------------------------------------------------


            if( masspr .gt. 0 ) then


                     acoll = 0.0
                     abloc = 0.0
                     acnne = 0.0
                     aeblc = 0.0


                     aelpn = 0.0
                     aelnd = 0.0
                     aeldd = 0.0


                     acnnd = 0.0
                     acnnr = 0.0
                     acndd = 0.0


                     acndn = 0.0
                     acnrn = 0.0
                     acddn = 0.0


                     acddd = 0.0
                     acrdr = 0.0
                     acdrd = 0.0
                     acrrr = 0.0


                     acdnd = 0.0
                     acdnr = 0.0
                     acrnd = 0.0
                     acrnr = 0.0


               do i = 0, ictime


                     acoll = acoll + float(nall( 1,i))
                     aelnd = aelnd + float(nall( 2,i))
                     aelpn = aelpn + float(nall( 3,i))
                     aeldd = aeldd + float(nall( 4,i))
                     abloc = abloc + float(nall( 5,i))
                     aeblc = aeblc + float(nall( 6,i))
                     acnne = acnne + float(nall( 7,i))
                     acnnd = acnnd + float(nall( 8,i))
                     acndn = acndn + float(nall( 9,i))
                     acnnr = acnnr + float(nall(10,i))
                     acnrn = acnrn + float(nall(11,i))
                     acndd = acndd + float(nall(12,i))
                     acddn = acddn + float(nall(13,i))
                     acddd = acddd + float(nall(14,i))
                     acrdr = acrdr + float(nall(15,i))
                     acdrd = acdrd + float(nall(16,i))
                     acrrr = acrrr + float(nall(17,i))
                     acdnd = acdnd + float(nall(18,i))
                     acdnr = acdnr + float(nall(19,i))
                     acrnd = acrnd + float(nall(20,i))
                     acrnr = acrnr + float(nall(21,i))


               end do


                     adnpi = 0.0
                     arnpi = 0.0
                     ardpi = 0.0


                     anpid = 0.0
                     anpir = 0.0
                     adpir = 0.0


                  do i = 0, ictime + 1


                     adnpi = adnpi + float(nall(22,i))
                     arnpi = arnpi + float(nall(23,i))
                     ardpi = ardpi + float(nall(24,i))


                     anpid = anpid + float(nall(25,i))
                     anpir = anpir + float(nall(26,i))
                     adpir = adpir + float(nall(27,i))


                  end do


            end if




*-----------------------------------------------------------------------
*     Collision History 01
*-----------------------------------------------------------------------


      if( jdsp(22) .eq. 0 ) return


*-----------------------------------------------------------------------


            ida = idsp(22)


*-----------------------------------------------------------------------


            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Collision History 01''/
     &            ''*'',71(''-''))')


            write(ida,'(/
     &            ''msul:{\rm\it JQMD  plotted by  \ANGEL}''/
     &            ''msdc:{\rm\Large \- \page  \-}''/
     &            ''msur:{\rm\it Collision History}'')')


            write(ida,'(/
     &            ''p: port nofr'')')


*-----------------------------------------------------------------------
*     Plot of Total Time Evolution
*-----------------------------------------------------------------------


            rmn  = float(iprun)


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Plot of Total Time Evolution''/
     &          ''*'',71(''-''))')


            write(ida,'( /''p: scal(0.4) xorg(0.0) yorg(4.1) nosx'')')


            write(ida,'( /''p: xmin(0) xmax('',f8.2,'')'')')
     &            s_time(ictime+1)


            write(ida,'( /a1,''Total Time Evolution''a1)') squ, squ
            write(ida,'( /''x: Time (fm/c)'')')
            write(ida,'(  ''y: E/A (MeV)'')')
            write(ida,'( /
     &          ''c:   t       epot        emas        epin    '',
     &                    ''    ekin'')')
            write(ida,'(  ''h:   x   y(E_{pot}),l0r y(mass),l0g'',
     &            '' y(pion),l0c y(E_{total}),l0b'')')


            write(ida,'(f7.2,1x,4e12.4)')
     &         ( s_time(i),
     &           s_epot(i)*1000./rmn,
     &          (s_epot(i)+s_emas(i)-s_emas(0))*1000./rmn,
     &          (s_epot(i)+s_emas(i)-s_emas(0)+s_epin(i))*1000./rmn,
     &          (s_epot(i)+s_emas(i)-s_emas(0)+s_epin(i)
     &           +s_ekin(i))*1000./rmn,
     &           i = 0, ictime )


*-----------------------------------------------------------------------
*     Plot of Collision History
*-----------------------------------------------------------------------


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Plot of Collision History''/
     &          ''*'',71(''-''))')


            write(ida,'( /''z: xorg(1.7)'')')
            write(ida,'( /a1,''Collision History'',a1)') squ, squ


            write(ida,'( /''p: xmin(0) xmax('',f8.2,'')'')')
     &            s_time(ictime+1)


            write(ida,'( /''x: Time (fm/c)'')')
            write(ida,'(  ''y: Collisions / Time'')')


            write(ida,'( /
     &          ''c:   t   all coll     blocked     inelastic'',
     &          ''   elastic     e-block'')')




            write(ida,'(''c:total'',5g12.4)')
     &                 acoll/rmn, abloc/rmn + aeblc/rmn,
     &                 acoll/rmn - abloc/rmn - aeblc/rmn - acnne/rmn,
     &                 acnne/rmn, aeblc/rmn


            write(ida,'(''h: x y(all),h0 y(blocked),h0b '',
     &          ''y(inelastic),h0r y(elastic),h0dc n(e-block)'')')


            write(ida,'(f7.2,5e12.4)')
     &            ( s_time(i-1),
     &              float(nall(1,i))/rmn,
     &              float(nall(5,i))/rmn
     &            + float(nall(6,i))/rmn,
     &              max( 0.0,
     &              float(nall(1,i))/rmn
     &            - float(nall(5,i))/rmn
     &            - float(nall(6,i))/rmn
     &            - float(nall(7,i))/rmn ),
     &              float(nall(7,i))/rmn,
     &              float(nall(6,i))/rmn,  i = 1, ictime )


            write(ida,'(f7.2,5e12.4)')
     &              s_time(ictime),
     &              float(nall(1,ictime))/rmn,
     &              float(nall(5,ictime))/rmn
     &            + float(nall(6,ictime))/rmn,
     &              max( 0.0,
     &              float(nall(1,ictime))/rmn
     &            - float(nall(5,ictime))/rmn
     &            - float(nall(6,ictime))/rmn
     &            - float(nall(7,ictime))/rmn ),
     &              float(nall(7,ictime))/rmn,
     &              float(nall(6,ictime))/rmn


*-----------------------------------------------------------------------
*     Plot of Elastic Collision History
*-----------------------------------------------------------------------


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Plot of Elastic Collision History''/
     &          ''*'',71(''-''))')


            write(ida,'( /''c: '',a1,''Elastic Collision History'',a1)')
     &                      squ, squ


            write(ida,'( /''p: xmin(0) xmax('',f8.2,'')'')')
     &            s_time(ictime+1)


            write(ida,'( /''c: Time (fm/c)'')')
            write(ida,'(  ''c: Collisions / Time'')')


            write(ida,'( /
     &          ''c:   t     pn          N-DR        DR-DR'')')




            write(ida,'(''c:total'',3g12.4)')
     &                aelpn/rmn, aelnd/rmn, aeldd/rmn


            write(ida,'(
     &          ''n:   x   y(pn),h0l   y(ND),h0m   y(DD),h0d'')')


            write(ida,'(F7.2,3E12.4)')
     &            ( s_time(i-1),
     &              float(nall( 3,i))/rmn,
     &              float(nall( 2,i))/rmn,
     &              float(nall( 4,i))/rmn, i = 1, ictime )


            write(ida,'(f7.2,3e12.4)')
     &              s_time(ictime),
     &              float(nall( 3,ictime))/rmn,
     &              float(nall( 2,ictime))/rmn,
     &              float(nall( 4,ictime))/rmn


*-----------------------------------------------------------------------
*     Plot of NN Inelastic Collision History
*-----------------------------------------------------------------------


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Plot of NN Inelastic Collision History''/
     &          ''*'',71(''-''))')


            write(ida,'( /''z: xorg(-1.7) yorg(-1.7)'')')
            write(ida,'( /a1,''NN Inelastic Collision and Inverse'',
     &                    a1)') squ, squ


            write(ida,'( /''p: xmin(0) xmax('',f8.2,'')'')')
     &            s_time(ictime+1)


            write(ida,'( /''x: Time (fm/c)'')')
            write(ida,'(  ''y: Collisions / Time'')')


            write(ida,'( /
     &          ''c:   t     NN->ND      NN->NR      NN->DD'')')


            write(ida,'(''c:total'',3g12.4)')
     &              acnnd/rmn, acnnr/rmn, acndd/rmn


            write(ida,'(''h: x y(NN\toND),h0l y(NN\toNR),h0m '',
     &                  ''y(NN\toDD),h0d'')')


            write(ida,'(f7.2,3e12.4)')
     &            ( s_time(i-1),
     &              float(nall( 8,i))/rmn,
     &              float(nall(10,i))/rmn,
     &              float(nall(12,i))/rmn, i = 1, ictime )


            write(ida,'(f7.2,3e12.4)')
     &              s_time(ictime),
     &              float(nall( 8,ictime))/rmn,
     &              float(nall(10,ictime))/rmn,
     &              float(nall(12,ictime))/rmn


*-----------------------------------------------------------------------
*     Plot of Inverse NN Inelastic Collision History
*-----------------------------------------------------------------------


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Plot of Inverse NN Inelastic Collision'',
     &               '' History''/
     &          ''*'',71(''-''))')


            write(ida,'( /''p: xmin(0) xmax('',f8.2,'')'')')
     &            s_time(ictime+1)


            write(ida,'( /''c: Time (fm/c)'')')
            write(ida,'(  ''c: Collisions / Time'')')


            write(ida,'( /
     &          ''c:   t     ND->NN      NR->NN      DD->NN'')')




            write(ida,'(''c:total'',3g12.4)')
     &              acndn/rmn, acnrn/rmn, acddn/rmn


            write(ida,'(''h: x y(ND\toNN),h0lr y(NR\toNN),h0mb '',
     &                  ''y(DD\toNN),h0dc'')')


            write(ida,'(f7.2,3e12.4)')
     &            ( s_time(i-1),
     &              float(nall( 9,i))/rmn,
     &              float(nall(11,i))/rmn,
     &              float(nall(13,i))/rmn, i = 1, ictime )


            write(ida,'(f7.2,3e12.4)')
     &              s_time(ictime),
     &              float(nall( 9,ictime))/rmn,
     &              float(nall(11,ictime))/rmn,
     &              float(nall(13,ictime))/rmn


*-----------------------------------------------------------------------
*     Plot of ND and NR Inelastic Collision History
*-----------------------------------------------------------------------


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Plot of ND and NR Inelastic Collision History''/
     &          ''*'',71(''-''))')


            write(ida,'( /''z: yorg(-1.7)'')')
            write(ida,'( /a1,
     &          ''ND and NR Inelastic Collision and Inverse'',
     &                    a1)') squ, squ


            write(ida,'( /''p: xmin(0) xmax('',f8.2,'')'')')
     &            s_time(ictime+1)


            write(ida,'( /''x: Time (fm/c)'')')
            write(ida,'(  ''y: Collisions / Time'')')


            write(ida,'( /
     &          ''c:   t     ND->DD      NR->DR      ND->RD'',
     &          ''      NR->RR'')')


            write(ida,'(''c:total'',4g12.4)')
     &              acddd/rmn, acrdr/rmn, acdrd/rmn, acrrr/rmn


            write(ida,'(''h: x y(ND\toDD),h0 Y(NR\toDR),h0m '',
     &                  ''Y(ND\toRD),h0d  y(NR\toRR),h0p'')')


            write(ida,'(f7.2,4e12.4)')
     &            ( s_time(i-1),
     &              float(nall(14,i))/rmn,
     &              float(nall(15,i))/rmn,
     &              float(nall(16,i))/rmn,
     &              float(nall(17,i))/rmn, i = 1, ictime )


            write(ida,'(f7.2,4e12.4)')
     &              s_time(ictime),
     &              float(nall(14,ictime))/rmn,
     &              float(nall(15,ictime))/rmn,
     &              float(nall(16,ictime))/rmn,
     &              float(nall(17,ictime))/rmn


*-----------------------------------------------------------------------
*     Plot of Inverse ND and NR Inelastic Collision History
*-----------------------------------------------------------------------


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Plot of Inverse ND and NR Inelastic '',
     &          ''Collision History''/
     &          ''*'',71(''-''))')




            write(ida,'( /''p: xmin(0) xmax('',f8.2,'')'')')
     &            s_time(ictime+1)


            write(ida,'( /''c: Time (fm/c)'')')
            write(ida,'(  ''c: Collisions / Time'')')


            write(ida,'( /''c:   t     DD->ND      DR->NR      RD->ND'',
     &                    ''      RR->NR'')')


            write(ida,'(''c:total'',4g12.4)')
     &              acdnd/rmn, acdnr/rmn, acrnd/rmn, acrnr/rmn


            write(ida,'(''h: x y(DD\toND),h0r y(DR\toNR),h0mb '',
     &                  ''y(RD\toND),h0dc y(RR\toNR),h0pr'')')


            write(ida,'(f7.2,4e12.4)')
     &            ( s_time(i-1),
     &              float(nall(18,i))/rmn,
     &              float(nall(19,i))/rmn,
     &              float(nall(20,i))/rmn,
     &              float(nall(21,i))/rmn, i = 1, ictime )


            write(ida,'(f7.2,4e12.4)')
     &              s_time(ictime),
     &              float(nall(18,ictime))/rmn,
     &              float(nall(19,ictime))/rmn,
     &              float(nall(20,ictime))/rmn,
     &              float(nall(21,ictime))/rmn


*-----------------------------------------------------------------------
*     Plot of Pion Production History
*-----------------------------------------------------------------------


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Plot of Pion Production History''/
     &          ''*'',71(''-''))')


            write(ida,'( /''z: xorg(1.7) yorg(1.7)'')')
            write(ida,'( /a1,''\pi  Production and Absorption'',
     &                    a1)') squ, squ


            write(ida,'( /''p: xmin(0) xmax('',f8.2,'')'')')
     &            s_time(ictime+1)


            write(ida,'( /''x: Time (fm/c)'')')
            write(ida,'(  ''y: Collisions / Time'')')


            write(ida,'( /
     &          ''c:   t    D->N+pi     R->N+pi     R->D+pi'')')


            write(ida,'(''c:total'',3g12.4)')
     &              adnpi/rmn, arnpi/rmn, ardpi/rmn


            write(ida,'(''h: x y(D\toN+\pi),h0 y(R\toN+\pi),h0m '',
     &                  ''y(R\toD+\pi),h0d'')')


            write(ida,'(f7.2,3e12.4)')
     &            ( s_time(i-1),
     &              float(nall(22,i))/rmn,
     &              float(nall(23,i))/rmn,
     &              float(nall(24,i))/rmn, i = 1, ictime + 1 )


            write(ida,'(f7.2,3e12.4)')
     &              s_time(ictime+1),
     &              float(nall(22,ictime+1))/rmn,
     &              float(nall(23,ictime+1))/rmn,
     &              float(nall(24,ictime+1))/rmn


*-----------------------------------------------------------------------
*     Plot of Pion Absorption History
*-----------------------------------------------------------------------


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Plot of Pion Absorption History''/
     &          ''*'',71(''-''))')


            write(ida,'( /''p: xmin(0) xmax('',f8.2,'')'')')
     &            s_time(ictime+1)


            write(ida,'( /''c: Time (fm/c)'')')
            write(ida,'(  ''c: Collisions / Time'')')


            write(ida,'( /
     &          ''c:   t    N+pi->D     N+pi->R     D+pi->R'')')


            write(ida,'(''c:total'',3g12.4)')
     &              anpid/rmn, anpir/rmn, adpir/rmn


            write(ida,'(''h: x y(N+\pi\toD),h0r '',
     &                 ''y(N+\pi\toR),h0mb '',
     &                 ''y(D+\pi\toR),h0dc'')')


            write(ida,'(F7.2,3E12.4)')
     &            ( s_time(i-1),
     &              float(nall(25,i))/rmn,
     &              float(nall(26,i))/rmn,
     &              float(nall(27,i))/rmn,  i = 1, ictime )


            write(ida,'(f7.2,3e12.4)')
     &              s_time(ictime),
     &              float(nall(25,ictime))/rmn,
     &              float(nall(26,ictime))/rmn,
     &              float(nall(27,ictime))/rmn


*-----------------------------------------------------------------------
*     Plot of pion, Delta and N* History
*-----------------------------------------------------------------------


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Plot of Pion, Delta and N* History''/
     &          ''*'',71(''-''))')


            write(ida,'( /''z: yorg(-1.7)'')')
            write(ida,'( /a1,''\pi, \Delta  and N^*'',
     &                    a1)') squ, squ


            write(ida,'( /''p: xmin(0) xmax('',f8.2,'')'')')
     &            s_time(ictime+1)


            write(ida,'( /''x: Time (fm/c)'')')
            write(ida,'(  ''y: Number of \pi, \Delta  and N^*'')')


            write(ida,'( /
     &          ''c:   t    all pion    all delt    all resn'')')


            write(ida,'(''c:total'',3g12.4)')
     &              float(npin(ictime+1))/rmn,float(ndel(ictime+1))/rmn,
     &              float(nres(ictime+1))/rmn


            write(ida,'(
     &          ''h:  x  y(\pi),h0r y(\Delta),h0b y(N^*),h0c'')')


            write(ida,'(f7.2,3e12.4)')
     &            ( s_time(i-1),
     &              float(npin(i))/rmn,
     &              float(ndel(i))/rmn,
     &              float(nres(i))/rmn,  i = 1, ictime + 1 )


            write(ida,'(f7.2,3e12.4)')
     &              s_time(ictime+1),
     &              float(npin(ictime+1))/rmn,
     &              float(ndel(ictime+1))/rmn,
     &              float(nres(ictime+1))/rmn


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine sm_grdo
*                                                                      *
*                                                                      *
*        Last Revised:     1998 12 09                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*             to write ground state properties for masspr = 0          *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


      common /const0/ idnta, idnpr, massta, masspr, mstapr, msprpr
      common /const2/ dt, ntmax, iprun, iprun0


      common /summ00/ jdsp(50), idsp(50), fdsp(50), ifds(50)
      character    fdsp*80


      common /summ06/ ictime, iltime
      common /summ07/ s_time(0:mmm)
      common /summ08/ s_ebin(0:mmm),s_epot(0:mmm),s_ekin(0:mmm),
     &                s_emas(0:mmm),s_epin(0:mmm),s_flow(0:mmm),
     &                s_rmsr(0:mmm),s_outp(0:mmm)


      common /summ13/ denmom(0:mmm), denden(0:mmm), xyz(2,0:50),
     &                denpat(0:mmm)
      common /summ14/ dentim(0:mmm,0:50), denerr(2,0:50), dentm(0:50),
     &                denmtm(0:mmm,0:50), denmer(2,0:50)
      common /summ15/ radtac, jmax, facp
      common /summ16/ csw, ccw, kmax, imax


*-----------------------------------------------------------------------


      character       squ*1
      data squ       /"'"/


*-----------------------------------------------------------------------
*     Ground State Properties for masspr = 0
*-----------------------------------------------------------------------


      if( jdsp(23) .eq. 0 ) return


*-----------------------------------------------------------------------


            ida = idsp(23)


*-----------------------------------------------------------------------
*     Ground State Properties (1)
*-----------------------------------------------------------------------


            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Ground State Properties (1)''/
     &            ''*'',71(''-''))')


            write(ida,'(/
     &            ''msul:{\rm\it JQMD  plotted by  \ANGEL}''/
     &            ''msdc:{\rm\Large \- \page  \-}''/
     &            ''msur:{\rm\it Ground State Properties (1)}'')')


            write(ida,'(/
     &            ''p: port nofr'')')


*-----------------------------------------------------------------------
*     Plot of Root Mean Squar Radius
*-----------------------------------------------------------------------


            rmn  = float(iprun)


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Plot of Root Mean Squar Radius''/
     &          ''*'',71(''-''))')


            write(ida,'( /''p: scal(0.7) xorg(0.1) yorg(1.9) nosx'')')
            write(ida,'(  ''p: ymin(0)'')')


            sek = 0.0


         do i = 0, ictime


            sek = sek + s_rmsr(i)


         end do


            av_rmsr = sek / float(ictime+1) / rmn


            write(ida,'(/a1,''\left<r^2\right>^{1/2} ='',g12.4,
     &                      '' fm'',a1)')
     &                  squ, av_rmsr, squ


            write(ida,'( /''x: Time (fm/c)'')')
            write(ida,'(  ''y: \left<r^2\right>^{1/2}  (fm)'')')


            write(ida,'( /''c:   t       rmsr'')')
            write(ida,'(  ''h:   x   '',
     &                    ''y(\left<r^2\right>^{1/2}),l0r'')')


            write(ida,'(F7.2,1X,E12.4)')
     &      ( s_time(i), s_rmsr(i)/rmn, i = 0, ictime )


*-----------------------------------------------------------------------
*     Plot of Exit Particle
*-----------------------------------------------------------------------


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Plot of Out Going Particle''/
     &          ''*'',71(''-''))')


            write(ida,'( /''c:   t       outp'')')
            write(ida,'(  ''h:   x     y(Out),hl0b'')')


            write(ida,'(f7.2,1x,e12.4)')
     &      ( s_time(i), s_outp(i)/rmn, i = 0, ictime )


*-----------------------------------------------------------------------
*     Plot of Binding Energy and Potential Energy
*-----------------------------------------------------------------------


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Plot of Binding Energy and Potential Energy''/
     &          ''*'',71(''-''))')


            write(ida,'( /''z: yorg(-1.7)'')')


            eliq0  = - bndeng(mstapr,massta-mstapr) / float(massta)


         if( eliq0 .lt. 0.0 ) then


            write(ida,'( /''p: ymax(0)'')')


         end if


            write(ida,'(/a1,''E_{bin} =  '',g10.4,''(  '',g10.4,
     &          '') MeV'',a1)') squ,
     &              (s_epot(0)+s_emas(0)-s_emas(0)+s_epin(0)
     &               +s_ekin(0))*1000./rmn, eliq0, squ




            write(ida,'( /''x: Time (fm/c)'')')
            write(ida,'(  ''y: E/A (MeV)'')')


            write(ida,'( /
     &          ''c:   t      epot        emas        epin     '',
     &                             ''   ekin'')')


            write(ida,'(  ''h:   x    y(E_{pot}),l0r  n'',
     &            ''          n       y(E_{total}),l0b'')')


            write(ida,'(f7.2,1x,4e12.4)')
     &      ( s_time(i),
     &        s_epot(i)*1000./rmn,
     &       (s_epot(i)+s_emas(i)-s_emas(0))*1000./rmn,
     &       (s_epot(i)+s_emas(i)-s_emas(0)+s_epin(i))*1000./rmn,
     &       (s_epot(i)+s_emas(i)-s_emas(0)+s_epin(i)
     &        +s_ekin(i))*1000./rmn,
     &        i = 0, ictime )


*-----------------------------------------------------------------------
*     Ground State Properties (2)
*-----------------------------------------------------------------------


            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Ground State Properties (2)''/
     &            ''*'',71(''-''))')


            write(ida,'(/
     &            ''newpage:'')')


            write(ida,'(/
     &            ''msul:{\rm\it JQMD  plotted by  \ANGEL}''/
     &            ''msdc:{\rm\Large \- \page  \-}''/
     &            ''msur:{\rm\it Ground State Properties (2)}'')')


            write(ida,'(/
     &            ''p: port nofr'')')


*-----------------------------------------------------------------------
*     Plot of Density Distribution
*-----------------------------------------------------------------------


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Plot of Density Distribution''/
     &          ''*'',71(''-''))')


            write(ida,'( /''p: scal(0.7) xorg(0.1) yorg(1.9) nosx'')')


            rrmax = sqrt(radtac)


            write(ida,'( ''p: xmax('',g12.4,'')'')') rrmax


            write(ida,'(/''x: r (fm)'')')


            write(ida,'( ''y: \rho(r)'')')
            write(ida,'(/''h:    x  y(QMD smeared),lst0b'',
     &                   '' d+           d-'')')


            do i = 0, 50


               denerr(1,i) = 0.0
               denerr(2,i) = 0.0


            end do


         do j = 0, ictime
         do i = 0, kmax - 1


            if( ( dentim(j,i) / rmn .gt.
     &            denden(i) / rmn / float(ictime+1) ) .and.
     &          ( dentim(j,i) / rmn - denden(i) / rmn / float(ictime+1)
     &       .gt. denerr(1,i) ) )
     &            denerr(1,i) = dentim(j,i) / rmn
     &                        - denden(i) / rmn / float(ictime+1)


            if( ( dentim(j,i) / rmn .lt.
     &            denden(i) / rmn / float(ictime+1) ) .and.
     &          ( dentim(j,i) / rmn - denden(i) / rmn / float(ictime+1)
     &       .lt. denerr(2,i) ) )
     &            denerr(2,i) = dentim(j,i) / rmn
     &                        - denden(i) / rmn / float(ictime+1)


         end do
         end do


         do i = 0, kmax - 1


            rabs = ( xyz(2,i) + xyz(1,i) ) / 2.0


            if( i .eq. 0 ) rabs = 0.0
*
            write(ida,'(4e13.5)')
     &            rabs, denden(i) / rmn / float(ictime+1),
     &            denerr(1,i), -denerr(2,i)


         end do




         do i = 0, imax


            if( i .eq. 0 ) then


               volum = 4.0 * pi/3.0 * 0.3535533


            else


               rmin  = sqrt( float(i) - 0.5 )
               rmax  = sqrt( float(i) + 0.5 )
               rdif  = rmax - rmin
               volum = 4.0 * pi * ( rmin**2 * rdif
     &                            + rmin * rdif**2
     &                            + rdif**3 / 3.0 )


            end if


               denpat(i) = denpat(i) / volum


         end do


            write(ida,'(/''h:    x          y(QMD Particle),ht0r'')')


         do i = 0, imax


            if( i .eq. 0 ) then


               rabs = 0.0


            else


               rabs = sqrt( ( float(i) + float(i-1) ) / 2.0 )


            end if


               write(ida,'(2e13.5)')
     &               rabs, denpat(i) / rmn / float(ictime+1)


         end do


         if( massta .eq. 40 ) then


            write(ida,'(/
     &       ''c: Plot of Density Distribution of HF for 40 Ca'',//
     &       ''h: x   y(HF \( ^{40}Ca\)),mst0b'',/
     &       ''   .000      .16867    ''/
     &       ''   .200      .16810    ''/
     &       ''   .400      .16653    ''/
     &       ''   .600      .16412    ''/
     &       ''   .800      .16123    ''/
     &       ''  1.000      .15826    ''/
     &       ''  1.200      .15560    ''/
     &       ''  1.400      .15352    ''/
     &       ''  1.600      .15211    ''/
     &       ''  1.800      .15120    ''/
     &       ''  2.000      .15041    ''/
     &       ''  2.200      .14919    ''/
     &       ''  2.400      .14693    ''/
     &       ''  2.600      .14303    ''/
     &       ''  2.800      .13704    ''/
     &       ''  3.000      .12865    ''/
     &       ''  3.200      .11782    ''/
     &       ''  3.400      .10476    ''/
     &       ''  3.600      .90021E-01''/
     &       ''  3.800      .74423E-01''/
     &       ''  4.000      .58983E-01''/
     &       ''  4.200      .44730E-01''/
     &       ''  4.400      .32488E-01''/
     &       ''  4.600      .22692E-01''/
     &       ''  4.800      .15351E-01''/
     &       ''  5.000      .10148E-01''/
     &       ''  5.200      .66156E-02''/
     &       ''  5.400      .42882E-02''/
     &       ''  5.600      .27836E-02''/
     &       ''  5.800      .18209E-02''/
     &       ''  6.000      .12078E-02'')')


         end if


            rrad = 1.124 * float(massta)**(1./3.)


            write(ida,'(/,
     &          ''c: Plot of Density Distributon for Matter'')')


            write(ida,'(''h: x        y(Matter),d0'')')


            write(ida,'(''   0.0      0.168''/
     &                       g12.4, ''0.168''/
     &                       g12.4, ''0.0'')') rrad, rrad


*-----------------------------------------------------------------------
*     Plot of Momentum Distribution
*-----------------------------------------------------------------------


            write(ida,'(/
     &          ''*'',71(''-'') /
     &          ''*     Plot of Momentum Distribution''/
     &          ''*'',71(''-''))')


            write(ida,'( /''z: yorg(-1.7)'')')
            write(ida,'( /''p: xmax(1.6)'')')


         do i = 0, jmax


            if( i .eq. 0 ) then


               volum = 4.0 * pi / 3.0 * (0.7071 / facp )**3 / hbc**3


            else


               pmin  = sqrt( float(i) - 0.5 ) / facp
               pmax  = sqrt( float(i) + 0.5 ) / facp
               pdif  = pmax - pmin
               volum = 4.0 * pi
     &               * ( pmin**2 * pdif
     &                 + pmin * pdif**2
     &                 + pdif**3 / 3.0 ) / hbc**3


            end if


               denmom(i) = denmom(i) / volum


            do j = 0, ictime


               denmtm(j,i) = denmtm(j,i) / volum


            end do


         end do


            do i = 0, 50


               denmer(1,i) = 0.0
               denmer(2,i) = 0.0


            end do


         do j = 0, ictime
         do i = 0, kmax - 1


            if( ( denmtm(j,i) / rmn .gt.
     &            denmom(i) / rmn / float(ictime+1) ) .and.
     &          ( denmtm(j,i) / rmn - denmom(i) / rmn / float(ictime+1)
     &       .gt. denmer(1,i) ) )
     &            denmer(1,i) = denmtm(j,i) / rmn
     &                        - denmom(i) / rmn / float(ictime+1)


            if( ( denmtm(j,i) / rmn .lt.
     &            denmom(i) / rmn / float(ictime+1) ) .and.
     &          ( denmtm(j,i) / rmn - denmom(i) / rmn / float(ictime+1)
     &       .lt. denmer(2,i) ) )
     &            denmer(2,i) = denmtm(j,i) / rmn
     &                        - denmom(i) / rmn / float(ictime+1)


         end do
         end do




            sek = 0.0


         do i = 0, ictime


            sek = sek + s_ekin(i)


         end do


            av_ekin = sek * 1000. / float(ictime+1) / rmn


            write(ida,'(/a1,''E_{kin} ='',g12.4,'' MeV'',a1)')
     &                  squ, av_ekin, squ


            write(ida,'(/''x: p (1/fm)'')')
            write(ida,'( ''y: \rho(p)'')')


            write(ida,'(/''h:    x          y(QMD),ht0r'',
     &                   ''  d+           d-'')')


         do i = 0, jmax


            if( i .eq. 0 ) then


               pabs = 0.0


            else


               pabs = sqrt( (float(i) + float(i-1) ) / 2.0 )
     &              / facp / hbc


            end if


               write( ida, '(4E13.5)')
     &                pabs, denmom(i) / rmn / float(ictime+1),
     &                denmer(1,i), -denmer(2,i)


         end do


      if( massta .eq. 40 ) then


         write(ida,'(/
     &       ''c: Plot of Momentum Distribution of HF for 40 Ca'',//
     &       ''h: x   y(HF \( ^{40}Ca\)),mst0b'',/
     &       ''  0.0    12.37     ''/
     &       ''  0.1    11.96     ''/
     &       ''  0.2    10.98     ''/
     &       ''  0.3     9.98     ''/
     &       ''  0.4     9.32     ''/
     &       ''  0.5     8.938    ''/
     &       ''  0.6     8.498    ''/
     &       ''  0.7     7.704    ''/
     &       ''  0.8     6.523    ''/
     &       ''  0.9     5.118    ''/
     &       ''  1.0     3.710    ''/
     &       ''  1.1     2.476    ''/
     &       ''  1.2     1.504    ''/
     &       ''  1.3     .8144    ''/
     &       ''  1.4     .3829    ''/
     &       ''  1.5     .1506    ''/
     &       ''  1.6     .04705   ''/
     &       ''  1.7     .01064   ''/
     &       ''  1.8     .001389  ''/
     &       ''  1.9     .0002762 ''/
     &       ''  2.0     .00003717'')')


      end if


         rhop = 3. * float(massta) / ( 4. * pi * 2.30 )


            write(ida,'(/,
     &          ''c: Plot of Momentum Distribution of Matter''//
     &          ''h: x      y(Matter),d0'')')


            write(ida,'('' 0.0    '',g12.4/
     &                  '' 1.32   '',g12.4/
     &                  '' 1.32      0.0'')') rhop, rhop


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine sm_cpuo
*                                                                      *
*                                                                      *
*        Last Revised:     1998 12 09                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*             to write cputime on screen and/or file                   *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      common /const2/ dt, ntmax, iprun, iprun0
      common /swich2/ icfg, imany, icpus, idatm


      common /startf/ iday0,imon0,iyer0,ihor0,imin0,isec0
      common /startt/ iday1,imon1,iyer1,ihor1,imin1,isec1
      common /cputim/ stime(30), cputm(30)


      common /summ00/ jdsp(50), idsp(50), fdsp(50), ifds(50)
      character    fdsp*80


*-----------------------------------------------------------------------


      if( jdsp(20) .eq. 0 ) return


*-----------------------------------------------------------------------
*     write cpu time on 6 and idf(4)
*-----------------------------------------------------------------------


      do 4000 jj = 1, 2


*-----------------------------------------------------------------------
*     File unit of output
*-----------------------------------------------------------------------


         if( jj .eq. 1 ) then


            if( jdsp(24) .eq. 0 ) goto 4000


            ida = 6


         else


            ida = idsp(20)


         end if


            write(ida,'(//16x,
     &          ''****************************************'')')


         if( iprun .eq. iprun0 ) then


            write(ida,'( 16x,
     &          ''* Qmd Calculation is Finished Normally *'')')


         else


            write(ida,'( 16x,''* Run is Terminated at Event = '',
     &                   i6,''  *'')') iprun


         end if


         if( idatm .eq. 1 ) then


            write(ida,'( 16x,
     &          ''*--------------------------------------*'')')
            write(ida,'( 16x,''* Starting time = '',
     &                       i4.4,''-'',i2.2,''-'',i2.2,''  '',
     &                       i2.2,'':'',i2.2,'':'',i2.2,'' *'')')
     &                       iyer0,imon0,iday0,ihor0,imin0,isec0


            write(ida,'( 16x,''*   Ending time = '',
     &                       i4.4,''-'',i2.2,''-'',i2.2,''  '',
     &                       i2.2,'':'',i2.2,'':'',i2.2,'' *'')')
     &                       iyer1,imon1,iday1,ihor1,imin1,isec1


         end if


         if( icpus .eq. 1 ) then


               ilpsh = int( stime(1) / 3600.0 )
               ilpsm = int( ( stime(1) - ilpsh * 3600.0 ) / 60.0 )
               elpss = stime(1) - ilpsh * 3600.0 - ilpsm * 60.0


            write(ida,'( 16x,
     &          ''*--------------------------------------*'')')
            write(ida,'( 16x,''*   Elapse time ='',i4,'' h'',
     &                                             i3,'' m'',
     &                                           f6.2,'' sec *'')')
     &                   ilpsh, ilpsm, elpss


            write(ida,'( 16x,''*     TOTAL CPU = '',f16.2,'' sec *'')')
     &                   cputm(1)
            write(ida,'( 16x,
     &          ''*--------------------------------------*'')')


                  ipcp2 = int(cputm(2)/cputm(1)*100.0)
                  ipcp3 = int(cputm(3)/cputm(1)*100.0)
                  ipcp4 = int(cputm(4)/cputm(1)*100.0)
                  ipcp5 = int(cputm(5)/cputm(1)*100.0)
                  ipcp6 = int(cputm(6)/cputm(1)*100.0)


                  ipcpr = 100 - ipcp2 - ipcp3
     &                        - ipcp4 - ipcp5 - ipcp6
                  cputr = cputm(1) - cputm(2) - cputm(3)
     &                  - cputm(4) - cputm(5) - cputm(6)


            write(ida,'( 16x,''*    Mean Field = '',
     &                   i3,'' %'',f11.2,'' sec *'')') ipcp2, cputm(2)
            write(ida,'( 16x,''*     Collision = '',
     &                   i3,'' %'',f11.2,'' sec *'')') ipcp3, cputm(3)
            write(ida,'( 16x,''*  Ground State = '',
     &                   i3,'' %'',f11.2,'' sec *'')') ipcp4, cputm(4)
            write(ida,'( 16x,''*       Summary = '',
     &                   i3,'' %'',f11.2,'' sec *'')') ipcp5, cputm(5)
            write(ida,'( 16x,''*  Final Decays = '',
     &                   i3,'' %'',f11.2,'' sec *'')') ipcp6, cputm(6)


            write(ida,'( 16x,''*        Others = '',
     &                   i3,'' %'',f11.2,'' sec *'')') ipcpr, cputr


         end if


            write(ida,'( 16x,
     &          ''****************************************'')')






*-----------------------------------------------------------------------


 4000 continue


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine sm_init
*                                                                      *
*                                                                      *
*        Last Revised:     1998 12 11                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to initiallize the summary and input echo               *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


      common /input1/ mstq1(mxpa1), parq1(mxpa1)


      common /verqmd/ versn, lastr, iyeav, imonv, idayv


      common /const0/ idnta, idnpr, massta, masspr, mstapr, msprpr
      common /const1/ elab, rdist, bmin, bmax, ibch, ibin
      common /const2/ dt, ntmax, iprun, iprun0
      common /const3/ nfreq, nfrec, nfred
      common /const4/ plab, srtcm, ylabb, pincm
      common /const5/ pzpr, pxpr, rzpr, rxpr
      common /const6/ pzta, pxta, rzta, rxta
      common /const7/ betpr, gampr, prmas, radpr
      common /const8/ betta, gamta, tamas, radta


      common /poten1/ gamm, c0, c3, cs, cl, wl


      common /rannum/ iseed, iseed0, iseed1


      common /swich1/ ipot, insys, irkg, icolt
      common /swich4/ ifin, ifout


      common /fnamid/ fname(mxfnm), idf(mxfnm), ifnl(mxfnm)
      character       fname*80
      common /inecho/ icnum, indata(200), indlng(200)
      character       indata*80


      common /startf/ iday0,imon0,iyer0,ihor0,imin0,isec0
      common /startt/ iday1,imon1,iyer1,ihor1,imin1,isec1


      common /summas/ sumas(3,0:maxpt,0:maxnt)


      common /summ00/ jdsp(50), idsp(50), fdsp(50), ifds(50)
      character    fdsp*80
      common /summ01/ iefw, ised, mulg, icdp(4), scald, isef
      common /summ02/ mdimg, pdimg, mminx, mmaxx, mminy, mmaxy,
     &                pminx, pmaxx, pminy, pmaxy, scalo, itdg
      common /summ03/ ireac(0:10), rcross(4)


      common /summ05/ npin(0:mmm), ndel(0:mmm), nres(0:mmm)
      common /summ08/ s_ebin(0:mmm),s_epot(0:mmm),s_ekin(0:mmm),
     &                s_emas(0:mmm),s_epin(0:mmm),s_flow(0:mmm),
     &                s_rmsr(0:mmm),s_outp(0:mmm)
      common /summ09/ epott0, epotp0, ekint0, ekinp0, ebint0, ebinp0
      common /summ12/ ncoll(30), nall(30,0:mmm)
      common /summ13/ denmom(0:mmm), denden(0:mmm), xyz(2,0:50),
     &                denpat(0:mmm)
      common /summ14/ dentim(0:mmm,0:50), denerr(2,0:50), dentm(0:50),
     &                denmtm(0:mmm,0:50), denmer(2,0:50)
      common /summ15/ radtac, jmax, facp
      common /summ16/ csw, ccw, kmax, imax


*-----------------------------------------------------------------------


      character       chta*16, chpr*16


*-----------------------------------------------------------------------


                  nfreq  = mstq1(171)
                  nfrec  = mstq1(172)
                  nfred  = max( 1, mstq1(173) )


*-----------------------------------------------------------------------
*        Parameters for Display
*-----------------------------------------------------------------------


                  jdsp(20) = mstq1(150)


                  jdsp(21) = mstq1(151)
                  jdsp(22) = mstq1(152)
                  jdsp(23) = mstq1(153)


                  jdsp(26) = mstq1(157)
                  jdsp(27) = mstq1(158)


                  jdsp(30) = mstq1(159)
                  jdsp(31) = mstq1(160)
                  jdsp(32) = mstq1(161)


                  jdsp(33) = mstq1(174)


                  jchois   = mstq1(154)


               do i = 1, 13


                  idsp(i) = 0
                  jdsp(i) = 0


                  if( jchois .eq. i ) jdsp(i) = 1


               end do


                  jdsp(24) = mstq1(155)
                  jdsp(25) = mstq1(156)


               if( masspr .eq. 0 ) then


                  jdsp(22) = 0


               else


                  jdsp(23) = 0


               end if


*-----------------------------------------------------------------------
*        for ifin > 0 : without QMD calculation
*-----------------------------------------------------------------------


               if( ifin .gt. 0 ) then


                     jdsp(22) = 0
                     jdsp(23) = 0


                     jdsp(25) = 0
                     jdsp(26) = 0


                  do i = 1, 13


                     jdsp(i) = 0


                  end do


               end if


*-----------------------------------------------------------------------
*        No output of summary
*-----------------------------------------------------------------------


               if( jdsp(20) .eq. 0 ) then


                     jdsp(21) = 0
                     jdsp(22) = 0
                     jdsp(23) = 0
                     jdsp(24) = 0
                     jdsp(25) = 0
                     jdsp(26) = 0
                     jdsp(30) = 0
                     jdsp(31) = 0
                     jdsp(32) = 0
                     jdsp(33) = 0


                  do i = 1, 13


                     jdsp(i) = 0


                  end do


                     goto 6000


               end if


*-----------------------------------------------------------------------
*        Some Parameters from Input
*-----------------------------------------------------------------------


                  iefw   = mstq1(166)
                  ised   = mstq1(167)
                  mulg   = mstq1(168)


                     if( ised .lt. 1 .or. ised .gt. 3 ) ised = 1


                  icdp(1) = mstq1(162)
                  icdp(2) = mstq1(163)
                  icdp(3) = mstq1(164)
                  icdp(4) = mstq1(165)


*-----------------------------------------------------------------------


                  rdimg   = parq1(151)
                  fdimg   = parq1(152)


                  xrmin   = parq1(153)
                  xrmax   = parq1(154)
                  yrmin   = parq1(155)
                  yrmax   = parq1(156)


                  xpmin   = parq1(157)
                  xpmax   = parq1(158)
                  ypmin   = parq1(159)
                  ypmax   = parq1(160)


                  scalo   = parq1(161)
                  scald   = parq1(162)


*-----------------------------------------------------------------------
*           size of display
*-----------------------------------------------------------------------


*              size of display for ground state


                  mdimg = nint( rdimg )
                  pdimg = fdimg


*              size of display for reactions


                  mminx = nint( xrmin )
                  mmaxx = nint( xrmax )
                  mminy = nint( yrmin )
                  mmaxy = nint( yrmax )


                  pminx = float(nint( xpmin ))
                  pmaxx = float(nint( xpmax ))
                  pminy = float(nint( ypmin ))
                  pmaxy = float(nint( ypmax ))


*-----------------------------------------------------------------------
*           Time digit selection
*-----------------------------------------------------------------------


               if( nint( nfrec * dt ) * 10 .eq.
     &             nint( 10.0 * nfrec * dt ) ) then


                  itdg = 0


               else


                  itdg = 1


               end if


*-----------------------------------------------------------------------
*        summary for cluster
*-----------------------------------------------------------------------


               do k = 1, 3
               do i = 0, maxpt
               do j = 0, maxnt


                  sumas(k,i,j) = 0.0


               end do
               end do
               end do


*-----------------------------------------------------------------------
*        summary of reaction type
*-----------------------------------------------------------------------


               do i = 0, 10


                  ireac(i) = 0


               end do


                  rcross(1) = ( bmax**2 - bmin**2 ) * pi * 10.0
                  rcross(2) = 0.0
                  rcross(3) = 0.0
                  rcross(4) = 0.0


*-----------------------------------------------------------------------
*        Some Constants
*-----------------------------------------------------------------------


                  radtac = ( 1.124 * float(massta)**(1./3.) * 1.7 ) ** 2
                  rrad  =  1.124 * float(massta)**(1./3.)


                  csw = 1.0 / 2.0 / wl
                  ccw = 1.0 / ( 2.0 * pi * wl ) ** ( 3. / 2. )


                  jmax  = 40
                  facp  = 14.0
                  facp  = 10.0


                  imax = ( nint(rrad) + 4 )**2


                  kmax = 24
                  rr0  = 1.87
                  vv0  = 4./3. * pi * rr0**3


               do i = 0, kmax


                  rr = float(i) * 0.5 + 1.0


                  if( i .eq. 0 ) then


                     xyz(1,i) = 0.0
                     xyz(2,i) = rr0


                  else


                     xyz(1,i) = ( rr**3 - 0.5 * rr0**3 )**(1./3.)
                     xyz(2,i) = ( rr**3 + 0.5 * rr0**3 )**(1./3.)


                  end IF


               end do


               do i = 0, mmm


                  do j = 1, 27


                     nall(j,i) = 0


                  end do


                  npin(i) = 0
                  ndel(i) = 0
                  nres(i) = 0


                  s_epot(i) = 0.0
                  s_ekin(i) = 0.0
                  s_emas(i) = 0.0
                  s_epin(i) = 0.0
                  s_ebin(i) = 0.0


                  s_outp(i) = 0.0
                  s_rmsr(i) = 0.0
                  s_flow(i) = 0.0


               end do


                  epott0 = 0.0
                  epotp0 = 0.0
                  ekint0 = 0.0
                  ekinp0 = 0.0
                  ebint0 = 0.0
                  ebinp0 = 0.0


               do i = 0, 50


                     denmom(i) = 0.0
                     denden(i) = 0.0
                     denpat(i) = 0.0


                  do j = 0, mmm


                     dentim(j,i) = 0.0
                     denmtm(j,i) = 0.0


                  end do


               end do




*-----------------------------------------------------------------------
*     Open files
*-----------------------------------------------------------------------
*        choice of display
*-----------------------------------------------------------------------
*
*           jdsp(21)   ; file for header, input echo and summary
*           jdsp(22)   ; file for the collision history
*           jdsp(23)   ; file for the ground state properties
*
*           jdsp(30)   ; file for the mass distribution of QMD
*           jdsp(31)   ; file for the mass distribution of SDM
*           jdsp(32)   ; file for the mass distribution of fission
*
*           jdsp(33)   ; file for QMDDISP
*
*           jdsp(1)    ; Snapshot of R-space in Threee Directions
*                        for Ground state.
*           jdsp(2)    ; Snapshot of P-space in Threee Directions
*                        for Ground state.
*           jdsp(3)    ; Snapshot of R and P-space in Threee Directions
*                        for Ground state.
*           jdsp(4)    ; Snapshot of R-space
*           jdsp(5)    ; Time Evolution of R-Space
*           jdsp(6)    ; Snapshot of P-space
*           jdsp(7)    ; Time Evolution of P-Space
*           jdsp(8)    ; Snapshot of R and P-space
*           jdsp(9)    ; Time Evolution of R and P-Space
*           jdsp(10)   ; Snapshot R-Space with Color Plot
*           jdsp(11)   ; Time Evolution of R-Space with Color Plot
*           jdsp(12)   ; Snapshot R-Space with Contour and Color Plot
*           jdsp(13)   ; Time Evolution of R-Space
*                        with Contour and Color Plot
*
*-----------------------------------------------------------------------
*
*           iefw    ; write one data of graph in one file
*
*           icdp(i) ; display particle on color plot
*
*           icdp(1) ; Display Nucleon on the Color Culuster Plot
*           icdp(2) ; Display Delta   on the Color Culuster Plot
*           icdp(3) ; Display N*      on the Color Culuster Plot
*           icdp(4) ; Display Pion    on the Color Culuster Plot
*
*           ised    ; three direction of the display
*                   1  ; z-x
*                   2  ; z-y
*                   3  ; x-y
*
*           mulg    ; number of graphs in a page
*
*           scalo   ; global scale
*
*-----------------------------------------------------------------------
*     Open file for the summary
*-----------------------------------------------------------------------


            fdsp(20) = fname(4)(1:ifnl(4))//'-all0.ang'


            idsp(20) = idf(4)


            ifds(20) = ifnl(4) + 9


            open( idsp(20), file = fdsp(20), status='unknown' )


*-----------------------------------------------------------------------
*     Open file for Header, Input Echo and Summary
*-----------------------------------------------------------------------


         if( jdsp(21) .ne. 0 ) then


            fdsp(21) = fname(4)(1:ifnl(4))//'-head.ang'


            idsp(21) = 60


            ifds(21) = ifnl(4) + 9


            open( idsp(21), file = fdsp(21), status='unknown' )


         end if


*-----------------------------------------------------------------------
*     Open file for Collision History
*-----------------------------------------------------------------------


         if( jdsp(22) .ne. 0 ) then


            fdsp(22) = fname(4)(1:ifnl(4))//'-col0.ang'


            idsp(22) = 61


            ifds(22) = ifnl(4) + 9


            open( idsp(22), file = fdsp(22), status='unknown' )


         end if


*-----------------------------------------------------------------------
*     Open file for Ground State Properties
*-----------------------------------------------------------------------


         if( jdsp(23) .ne. 0 ) then


            fdsp(23) = fname(4)(1:ifnl(4))//'-grd0.ang'


            idsp(23) = 62


            ifds(23) = ifnl(4) + 9


            open( idsp(23), file = fdsp(23), status='unknown' )


         end if


*-----------------------------------------------------------------------
*     Open file for Mass Distribution of QMD
*-----------------------------------------------------------------------


         if( jdsp(30) .ne. 0 ) then


            fdsp(30) = fname(4)(1:ifnl(4))//'-mqmd.ang'


            idsp(30) = 63


            ifds(30) = ifnl(4) + 9


            open( idsp(30), file = fdsp(30), status='unknown' )


         end if


*-----------------------------------------------------------------------
*     Open file for Mass Distribution of SDM
*-----------------------------------------------------------------------


         if( jdsp(31) .ne. 0 ) then


            fdsp(31) = fname(4)(1:ifnl(4))//'-msdm.ang'


            idsp(31) = 64


            ifds(31) = ifnl(4) + 9


            open( idsp(31), file = fdsp(31), status='unknown' )


         end if


*-----------------------------------------------------------------------
*     Open file for Mass Distribution of fission
*-----------------------------------------------------------------------


         if( jdsp(32) .ne. 0 ) then


            fdsp(32) = fname(4)(1:ifnl(4))//'-mfis.ang'


            idsp(32) = 65


            ifds(32) = ifnl(4) + 9


            open( idsp(32), file = fdsp(32), status='unknown' )


         end if


*-----------------------------------------------------------------------
*     Open file for QMDDISP
*-----------------------------------------------------------------------


         if( jdsp(33) .ne. 0 ) then


            fdsp(33) = fname(4)(1:ifnl(4))//'-disp.shn'


            idsp(33) = 66


            ifds(33) = ifnl(4) + 9


            open( idsp(33), file = fdsp(33), status='unknown' )


*-----------------------------------------------------------------------


            io = idsp(33)


            ievnt = min( iprun, jdsp(33) )


            write(io,'(7i4,3e13.5)')
     &                  ievnt,
     &                  massta, masspr, mstapr, msprpr,
     &                  0, 0, 0.0, 0.0, elab*1000.0


            write(io,'(2i6)')
     &                  int( float( ntmax / nfred * nfred ) * dt ),
     &                  int( dt * nfred )


         end if


*-----------------------------------------------------------------------
*     Open file for grand coordinate
*-----------------------------------------------------------------------


         if( jdsp(1) .ne. 0 ) then


            isef = 0


            idsp(1) = 70


            fdsp(1) = fname(4)(1:ifnl(4))//'-gr00.ang'


            ifds(1) = ifnl(4) + 9


            open( idsp(1), file = fdsp(1), status='unknown' )


         end if


*-----------------------------------------------------------------------
*     Open file for grand momentum
*-----------------------------------------------------------------------


         if( jdsp(2) .ne. 0 ) then


            isef = 0


            idsp(2) = 72


            fdsp(2) = fname(4)(1:ifnl(4))//'-gp00.ang'


            ifds(2) = ifnl(4) + 9


            open( idsp(2), file = fdsp(2), status='unknown' )


         end if


*-----------------------------------------------------------------------
*     Open file for grand coordinate and momentum
*-----------------------------------------------------------------------


         if( jdsp(3) .ne. 0 ) then


            isef = 1


            idsp(1) = 74
            idsp(2) = 74
            idsp(3) = 74


            fdsp(3) = fname(4)(1:ifnl(4))//'-grp0.ang'


            ifds(3) = ifnl(4) + 9


            open( idsp(3), file = fdsp(3), status='unknown' )


         end if


*-----------------------------------------------------------------------
*     Open file for coordinate
*-----------------------------------------------------------------------


         if( jdsp(4) .ne. 0 ) then


            isef = 0


            idsp(4) = 76


            fdsp(4) = fname(4)(1:ifnl(4))//'-r000.ang'


            ifds(4) = ifnl(4) + 9


            open( idsp(4), file = fdsp(4), status='unknown' )


         end if


*-----------------------------------------------------------------------
*     Open file for coordinate time evolution
*-----------------------------------------------------------------------


         if( jdsp(5) .ne. 0 ) then


            isef = 2


            idsp(4) = 76
            idsp(5) = 76


            fdsp(5) = fname(4)(1:ifnl(4))//'-rt00.ang'


            ifds(5) = ifnl(4) + 9


            open( idsp(5), file = fdsp(5), status='unknown' )


         end if


*-----------------------------------------------------------------------
*     Open file for momentum
*-----------------------------------------------------------------------


         if( jdsp(6) .ne. 0 ) then


            isef = 0


            idsp(6) = 78


            fdsp(6) = fname(4)(1:ifnl(4))//'-p000.ang'


            ifds(6) = ifnl(4) + 9


            open( idsp(6), file = fdsp(6), status='unknown' )


         end if


*-----------------------------------------------------------------------
*     Open file for momentum time evolution
*-----------------------------------------------------------------------


         if( jdsp(7) .ne. 0 ) then


            isef = 2


            idsp(6) = 78
            idsp(7) = 78


            fdsp(7) = fname(4)(1:ifnl(4))//'-pt00.ang'


            ifds(7) = ifnl(4) + 9


            open( idsp(7), file = fdsp(7), status='unknown' )


         end if


*-----------------------------------------------------------------------
*     Open file for coordinate and momentum
*-----------------------------------------------------------------------


         if( jdsp(8) .ne. 0 ) then


            isef = 1


            idsp(4) = 78
            idsp(6) = 78
            idsp(8) = 78


            fdsp(8) = fname(4)(1:ifnl(4))//'-pr00.ang'


            ifds(8) = ifnl(4) + 9


            open( idsp(8), file = fdsp(8), status='unknown' )


         end if


*-----------------------------------------------------------------------
*     Open file for coordinate and momentum time evolution
*-----------------------------------------------------------------------


         if( jdsp(9) .ne. 0 ) then


            isef = 3


            idsp(4) = 78
            idsp(6) = 78
            idsp(9) = 78


            fdsp(9) = fname(4)(1:ifnl(4))//'-prt0.ang'


            ifds(9) = ifnl(4) + 9


            open( idsp(9), file = fdsp(9), status='unknown' )


         end if


*-----------------------------------------------------------------------
*     Open file for color coordinate
*-----------------------------------------------------------------------


         if( jdsp(10) .ne. 0 ) then


            isef = 0


            idsp(10) = 80


            fdsp(10) = fname(4)(1:ifnl(4))//'-rc00.ang'


            ifds(10) = ifnl(4) + 9


            open( idsp(10), file = fdsp(10), status='unknown' )


         end if


*-----------------------------------------------------------------------
*     Open file for color coordinate time evolution
*-----------------------------------------------------------------------


         if( jdsp(11) .ne. 0 ) then


            isef = 2


            idsp(10) = 80
            idsp(11) = 80


            fdsp(11) = fname(4)(1:ifnl(4))//'-rtc0.ang'


            ifds(11) = ifnl(4) + 9


            open( idsp(11), file = fdsp(11), status='unknown' )


         end if


*-----------------------------------------------------------------------
*     Open file for coordinate and color coordinate
*-----------------------------------------------------------------------


         if( jdsp(12) .ne. 0 ) then


            isef = 4


            idsp( 4) = 80
            idsp(10) = 80
            idsp(12) = 80


            fdsp(12) = fname(4)(1:ifnl(4))//'-rrc0.ang'


            ifds(12) = ifnl(4) + 9


            open( idsp(12), file = fdsp(12), status='unknown' )


         end if


*-----------------------------------------------------------------------
*     Open file for coordinate and color coordinate time evolution
*-----------------------------------------------------------------------


         if( jdsp(13) .ne. 0 ) then


            isef = 5


            idsp( 4) = 80
            idsp(10) = 80
            idsp(13) = 80


            fdsp(13) = fname(4)(1:ifnl(4))//'-rrtc.ang'


            ifds(13) = ifnl(4) + 9


            open( idsp(13), file = fdsp(13), status='unknown' )


         end if




*-----------------------------------------------------------------------
*     write JQMD logo and reaction parameters on the front page
*     of the output file and/or on the screan
*-----------------------------------------------------------------------
*        write title on 6, idf(4), idf(10)
*-----------------------------------------------------------------------


 6000 continue


      do 5000 jj = 1, 3


*-----------------------------------------------------------------------
*     File unit of output
*-----------------------------------------------------------------------


         if( jj .eq. 1 ) then


            if( jdsp(24) .eq. 0 ) goto 5000


               ida = 6


         else if( jj .eq. 2 ) then


            if( jdsp(20) .eq. 0 ) goto 5000


               ida = idsp(20)


         else if( jj .eq. 3 ) then


            if( ifout .eq. 0 ) goto 5000


               ida = idf(10)


         end if


*-----------------------------------------------------------------------
*     Logo of JQMD
*-----------------------------------------------------------------------


      write(ida,'(/
     &   12X,   ''       JJJ     QQQQQ     MMM     MMM  DDDDDDD''/
     &   12X,   ''       JJJ   QQQQQQQQQ   MMMM   MMMM  DDDDDDDD''/
     &   12X,   ''       JJJ  QQQ     QQQ  MMMMM MMMMM  DDD   DDD''/
     &   12X,   ''       JJJ  QQQ     QQQ  MMM MMM MMM  DDD    DDD''/
     &   12X,   ''       JJJ  QQQ     QQQ  MMM  M  MMM  DDD    DDD''/
     &   12X,   ''       JJJ  QQQ     QQQ  MMM     MMM  DDD    DDD''/
     &   12X,   ''       JJJ  QQQ QQQ QQQ  MMM     MMM  DDD    DDD''/
     &   12X,   ''JJJ    JJJ   QQQQQQQQQ   MMM     MMM  DDD   DDD''/
     &   12X,   '' JJJJJJJJ      QQQQQQ    MMM     MMM  DDDDDDDD''/
     &   12X,   ''   JJJJ            QQQ   MMM     MMM  DDDDDDD'')')


      write(ida,'(//
     &   17x,      ''    Jaeri Quantum Molecular Dynamics''/
     &   17x,      ''         for Nuclear Reactions''//
     &
     &   17x,      ''           Version ='',f5.2,//
     &
     &   17x,      ''              made by''//
     &
     &   17x,      ''Japan Atomic Energy Research Institute''//
     &
     &   17x,      ''      Last Revised  '',i4,1x,i2.2,1x,i2.2)')
     &
     &                        versn, iyeav,imonv, idayv


*-----------------------------------------------------------------------
*     Input Echo
*-----------------------------------------------------------------------


      if( icnum .gt. 0 ) then


            write(ida,'(//16x,
     &         ''*------------- Input Echo -------------*''/)')


         do i = 1, icnum


            write(ida,'(24x,80a1)')
     &                 ( indata(i)(j:j), j =  1, 11 ),
     &                 ( indata(i)(j:j), j = 13, 15 ),
     &                 ( indata(i)(j:j), j = 17, indlng(i) )


         end do


            write(ida,'(/16x,
     &         ''*------------- Input Echo -------------*''/)')


      end if


*-----------------------------------------------------------------------
*     Reaction parameters
*-----------------------------------------------------------------------


         write(ida,'(/3x,67(''*'')/3x,''*'',65X,''*'')')


      if( masspr .gt. 0 ) then


*-----------------------------------------------------------------------


            write(ida,'(
     &            3x,''*'',10x,''Reaction :'',45x,''*''/
     &            3x,''*'',65x,''*'')')


               call chname(idnta,massta,mstapr,chta,nchta)
               call chname(idnpr,masspr,msprpr,chpr,nchpr)


            write(ida,'(
     &            3x,''*'',14x,a16,''  on  '',a16,13x,''*''/
     &            3x,''*'',65x,''*'')') chpr, chta


                  msprne = masspr-msprpr
                  mstane = massta-mstapr


                  if( idnta .ne. 0 ) msprne = 0
                  if( idnpr .ne. 0 ) mstane = 0


            write(ida,'(
     &            3x,''*'',14x,''mass'',i4,
     &            ''('',i3,'','',i3,'')'',
     &            '' ==> mass'',i4,''('',i3,'','',i3,'')'',12x,''*'')')
     &            masspr,msprpr,msprne,
     &            massta,mstapr,mstane


         if( elab .gt. 1.0 ) then


            write(ida,'(
     &            3x,''*'',65x,''*''/
     &            3x,''*'',16x,
     &            ''Beam energy   = '',f7.2,'' A GeV'',20x,''*''/
     &            3x,''*'',16x,
     &            ''Beam momentum = '',f7.2,'' A GeV/c'',18x,''*''/
     &            3x,''*'',16x,
     &            ''NN CM energy  = '',f7.2,''   GeV'',20x,''*'')')
     &            elab, plab, srtcm - prmas - tamas


         else


            write(ida,'(
     &            3x,''*'',65x,''*''/
     &            3x,''*'',16x,
     &            ''Beam energy   = '',f7.2,'' A MeV'',20x,''*''/
     &            3x,''*'',16x,
     &            ''Beam momentum = '',f7.2,'' A MeV/c'',18x,''*''/
     &            3x,''*'',16x,
     &            ''NN CM energy  = '',f7.2,''   MeV'',20x,''*'')')
     &            elab*1000., plab*1000.,
     &            ( srtcm - prmas - tamas ) * 1000.


         end if


*-----------------------------------------------------------------------


         if( insys .eq. 1 ) then


            write(ida,'(
     &            3x,''*'',65x,''*''/3x,''*'',10x,
     &            ''Calculated in C.M. frame :'',29x,''*'')')


         else if( insys .eq. 0 ) then


            write(ida,'(
     &            3x,''*'',65x,''*''/3x,''*'',10x,
     &            ''Calculated in Lab. frame :'',29x,''*'')')


         else if( insys .eq. 2 ) then


            write(ida,'(
     &            3x,''*'',65x,''*''/3x,''*'',10x,
     &            ''Calculated in N-N C.M. frame :'',25x,''*'')')


         end if


*-----------------------------------------------------------------------


            write(ida,'(
     &            3x,''*'',65x,''*'')')


            write(ida,'(
     &            3x,''*'',38x,
     &            ''Beam'',8x,''Target'',9x,''*'')')


            write(ida,'(
     &            3x,''*'',16x,
     &            ''Velocity / c :'',3x,f9.7,5x,f9.7,9x,''*'')')
     &            betpr, betta


         if( elab .gt. 1.0 ) then


            write(ida,'(
     &            3x,''*'',16x,
     &            ''Gamma factor :'',3x,f9.3,5x,f9.3,9x,''*'')')
     &            gampr, gamta


            write(ida,'(
     &            3x,''*'',16x,
     &            ''p_z  (GeV/c) :'',3x,f9.3,5x,f9.3,9x,''*'')')
     &            pzpr, pzta


         else


            write(ida,'(
     &            3x,''*'',16x,
     &            ''Gamma factor :'',3x,f9.5,5x,f9.5,9x,''*'')')
     &            gampr, gamta


            write(ida,'(
     &            3x,''*'',16x,
     &            ''p_z  (MeV/c) :'',3x,f9.3,5x,f9.3,9x,''*'')')
     &            pzpr*1000., pzta*1000.


         end if


            write(ida,'(
     &            3x,''*'',16x,
     &            ''r_z  (fm)    :'',3x,f9.3,5x,f9.3,9x,''*'')')
     &            rzpr, rzta


*-----------------------------------------------------------------------


            write(ida,'(
     &            3x,''*'',65x,''*''/3x,''*'',10x,
     &            ''Impact Parameter Range :'',31x,''*'')')


            write(ida,'(
     &            3x,''*'',65x,''*''/3x,''*'',16x,
     &            f7.3,''   < b < ''f7.3,'' (fm)'',21x,''*'')')
     &            bmin, bmax


*-----------------------------------------------------------------------


      else if( masspr .eq. 0 ) then


*-----------------------------------------------------------------------


            write(ida,'(
     &             3x,''*'',21x,''Nucleus'',i5,
     &             ''('',i3,'','',i3,'')'',23x,''*'')')
     &             massta,mstapr,massta-mstapr
            write(ida,'(
     &             3x,''*'',65x,''*''/3x,''*'',20x,
     &             ''Ground State Properties :'',20x,''*'')')


*-----------------------------------------------------------------------


      end if


*-----------------------------------------------------------------------


            write(ida,'(
     &            3x,''*'',65x,''*''/3x,''*'',10x,
     &            ''Time Evolution :'',39x,''*'')')


            write(ida,'(
     &            3x,''*'',65x,''*''/3x,''*'',16x,
     &            ''Number of Events      = '',i12,13x,''*'')') iprun


            write(ida,'(
     &            3x,''*'',16x,
     &            ''Number of Time Step   = '',i12,13x,''*'')') ntmax


            write(ida,'(
     &            3x,''*'',16x,
     &            ''Time step (fm/c)      = '',f12.3,13x,''*'')') dt


            write(ida,'(
     &            3x,''*'',16x,
     &            ''Seed of random number = '',i12,13x,''*'')') iseed0




*-----------------------------------------------------------------------


            write(ida,'(
     &            3x,''*'',65x,''*''/3x,''*'',10x,
     &            ''Date :'',49x,''*'')')


            write(ida,'(
     &            3x,''*'',65x,''*''/3x,''*'',16x,
     &            ''Calculated at   '',
     &            i4.4,''-'',i2.2,''-'',i2.2,''   '',
     &            i2.2,'':'',i2.2,'':'',i2.2,12x''*'')')
     &            iyer0,imon0,iday0,ihor0,imin0,isec0


*-----------------------------------------------------------------------


            write(ida,'(3x,''*'',65x,''*''/3x,67(''*''))')


*-----------------------------------------------------------------------
*        Header for output file of QMD results
*-----------------------------------------------------------------------


         if( jj .eq. 3 ) then


            write(ida,'()')


            write(ida,'(''  ik  iz  in  ic  jj  id  is  iq  im'',
     &                   6x,''px'',12x,''py'',12x,''pz'',
     &                  12x,''et'',12x,''rm'',12x,''ex'')')


         end if


*-----------------------------------------------------------------------


 5000 continue


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine sm_evnt
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 26                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to initialize the summary at each event                 *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'


*-----------------------------------------------------------------------


      common /const0/ idnta, idnpr, massta, masspr, mstapr, msprpr


      common /summ00/ jdsp(50), idsp(50), fdsp(50), ifds(50)
      character    fdsp*80


      common /summ06/ ictime, iltime
      common /summ09/ epott0, epotp0, ekint0, ekinp0, ebint0, ebinp0
      common /summ12/ ncoll(30), nall(30,0:mmm)


*-----------------------------------------------------------------------


      dimension       it(0:nnn)
      dimension       rs(3,nnn), ps(3,nnn)


*-----------------------------------------------------------------------


      if( jdsp(20) .eq. 0 ) return


*-----------------------------------------------------------------------
*        time index
*-----------------------------------------------------------------------


               ictime = -1


*-----------------------------------------------------------------------
*        for collision summary
*-----------------------------------------------------------------------


            do i = 1, 27


               ncoll(i)  = 0


            end do


*-----------------------------------------------------------------------
*        ground state energy after boosted
*-----------------------------------------------------------------------


            do i = 1, massta


               it(i) = i


            end do


               it(0) = massta


               call etotal(1,it,rs,ps,ekint,epott,ebint,emast,epin,jj)


            do i = 1, masspr


               it(i) = i + massta


            end do


               it(0) = masspr


               call etotal(1,it,rs,ps,ekinp,epotp,ebinp,emasp,epin,jj)


               epott0 = epott0 + epott
               epotp0 = epotp0 + epotp
               ekint0 = ekint0 + ekint
               ekinp0 = ekinp0 + ekinp
               ebint0 = ebint0 + ebint
               ebinp0 = ebinp0 + ebinp


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine sm_timd
*                                                                      *
*                                                                      *
*        Last Revised:     1998 12 09                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to output the phase space information for QMDDISP       *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'


*-----------------------------------------------------------------------


      common /const2/ dt, ntmax, iprun, iprun0


      common /vriab0/ massal, massba, nmeson
      common /vriab1/ b, llnow, ntnow


      common /coodrp/ r(3,nnn),  p(5,nnn)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)


      common /summ00/ jdsp(50), idsp(50), fdsp(50), ifds(50)
      character    fdsp*80


*-----------------------------------------------------------------------


         if( jdsp(33) .eq. 0 ) return


         if( llnow .gt. jdsp(33) ) return


            io = idsp(33)


*-----------------------------------------------------------------------
*     output QMDDISP
*-----------------------------------------------------------------------


         if( ntnow .eq. 0 ) then


            write(io,'(e13.5)') b


         end if


            write(io,'(2i6)') int( ntnow * dt ), massal


         do i = 1, massal


               idpar = inds(i)
               if( idpar .eq. 4 ) idpar = 50


            write(io,'('' '',1p3e13.5,3i4)')
     &          r(1,i), r(2,i), r(3,i), idpar, ichg(i), 0


         end do


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine sm_timc
*                                                                      *
*                                                                      *
*        Last Revised:     1998 12 09                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to output the phase space information for summary       *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      common /const2/ dt, ntmax, iprun, iprun0
      common /const3/ nfreq, nfrec, nfred


      common /vriab1/ b, llnow, ntnow


      common /summ00/ jdsp(50), idsp(50), fdsp(50), ifds(50)
      character    fdsp*80
      common /summ01/ iefw, ised, mulg, icdp(4), scald, isef
      common /summ02/ mdimg, pdimg, mminx, mmaxx, mminy, mmaxy,
     &                pminx, pmaxx, pminy, pmaxy, scalo, itdg


*-----------------------------------------------------------------------


      save            scalp, formp


*-----------------------------------------------------------------------


         if( jdsp(20) .eq. 0 ) return


*-----------------------------------------------------------------------


         if( idsp(1)  .eq. 0 .and.
     &       idsp(2)  .eq. 0 .and.
     &       idsp(4)  .eq. 0 .and.
     &       idsp(6)  .eq. 0 .and.
     &       idsp(10) .eq. 0 .and.
     &       jdsp(25) .eq. 0 .and.
     &       jdsp(26) .eq. 0 ) return


         if( llnow .gt. jdsp(27) ) return


*-----------------------------------------------------------------------
*     cluster analysis
*-----------------------------------------------------------------------


               ptime = dt * float(ntnow)


               call cldist


*-----------------------------------------------------------------------
*     display particle position on terminal and/or file
*-----------------------------------------------------------------------


            if( jdsp(25) .ne. 0 .or. jdsp(26) .ne. 0 ) then


               call disp06(ptime)


            end if


*-----------------------------------------------------------------------
*     write the summary on the file
*-----------------------------------------------------------------------


            if(  ntnow + nfrec .gt. ntmax ) then


               itim = 1


            else


               itim = 0


            end if


            if( ptime .lt. 0.0001 ) then


               scalp = -1.0
               formp = -1.0


            end if


*-----------------------------------------------------------------------


         if( idsp(1) .ne. 0 ) then


            call disp01(mdimg,scalo,scalp,isef,
     &                  ptime,itim,itdg,iefw)


         end if


*-----------------------------------------------------------------------


         if( idsp(2) .ne. 0 ) then


            call disp07(pdimg,scalo,scalp,isef,
     &                  ptime,itim,itdg,iefw)


         end if


*-----------------------------------------------------------------------


         if( idsp(4) .ne. 0 ) then


            call disp02(mminx,mmaxx,mminy,mmaxy,
     &                  scalo,scalp,formp,ised,isef,mulg,
     &                  ptime,itim,itdg,iefw)


         end if


*-----------------------------------------------------------------------


         if( idsp(6) .ne. 0 ) then


            call disp08(pminx,pmaxx,pminy,pmaxy,
     &                  scalo,scalp,formp,ised,isef,mulg,
     &                  ptime,itim,itdg,iefw)


         end if


*-----------------------------------------------------------------------


         if( idsp(10) .ne. 0 ) then


            call disp04(mminx,mmaxx,mminy,mmaxy,
     &                  scalo,scalp,formp,ised,isef,mulg,
     &                  ptime,itim,itdg,icdp,iefw)


         end if


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine sm_timq
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 26                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to store some values at the time steps for summary      *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


      common /const0/ idnta, idnpr, massta, masspr, mstapr, msprpr
      common /const2/ dt, ntmax, iprun, iprun0
      common /const3/ nfreq, nfrec, nfred


      common /vriab0/ massal, massba, nmeson
      common /vriab1/ b, llnow, ntnow
      common /swich1/ ipot, insys, irkg, icolt


      common /poten1/ gamm, c0, c3, cs, cl, wl


      common /coodrp/ r(3,nnn),  p(5,nnn)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)


      common /summ00/ jdsp(50), idsp(50), fdsp(50), ifds(50)
      character    fdsp*80


      common /summ05/ npin(0:mmm), ndel(0:mmm), nres(0:mmm)
      common /summ06/ ictime, iltime
      common /summ07/ s_time(0:mmm)
      common /summ08/ s_ebin(0:mmm),s_epot(0:mmm),s_ekin(0:mmm),
     &                s_emas(0:mmm),s_epin(0:mmm),s_flow(0:mmm),
     &                s_rmsr(0:mmm),s_outp(0:mmm)
      common /summ12/ ncoll(30), nall(30,0:mmm)
      common /summ13/ denmom(0:mmm), denden(0:mmm), xyz(2,0:50),
     &                denpat(0:mmm)
      common /summ14/ dentim(0:mmm,0:50), denerr(2,0:50), dentm(0:50),
     &                denmtm(0:mmm,0:50), denmer(2,0:50)
      common /summ15/ radtac, jmax, facp
      common /summ16/ csw, ccw, kmax, imax


*-----------------------------------------------------------------------


         if( jdsp(20) .eq. 0 ) return


*-----------------------------------------------------------------------
*        time variables
*-----------------------------------------------------------------------


                  ictime = ictime + 1


            if( llnow .eq. 1 ) then


                  s_time(ictime) = dt * float(ntnow)


               if( ntnow + nfreq .gt. ntmax ) then


                  iltime = ictime


                  s_time(iltime+1) = s_time(ictime)
     &                             + s_time(1) - s_time(0)


               end if


            end if


*-----------------------------------------------------------------------
*        potential, kinetic energy
*-----------------------------------------------------------------------


                  call epotall(epot)


                  ekin  = 0.0
                  emas  = 0.0
                  epin  = 0.0


            do i = 1, massal


               if( inds(i) .eq. 4 ) then


                  epin = epin + p(4,i)


               else


                  emas = emas + p(5,i)
                  ekin = ekin + p(4,i) - p(5,i)


               end if


            end do


               s_epot(ictime) = s_epot(ictime) + epot / float( massba )
               s_ekin(ictime) = s_ekin(ictime) + ekin / float( massba )
               s_emas(ictime) = s_emas(ictime) + emas / float( massba )
               s_epin(ictime) = s_epin(ictime) + epin / float( massba )
               s_ebin(ictime) = s_ebin(ictime)
     &                        +  ( ekin + epot ) / float( massba )


*-----------------------------------------------------------------------
*        collision summary
*-----------------------------------------------------------------------


         if( icolt .eq. 1 ) then


            if( ictime .eq. 0 ) then


               if( idnta .ne. 0 ) npin(0) = npin(0) + 1
               if( idnpr .ne. 0 ) npin(0) = npin(0) + 1


            else if( ictime .gt. 0 ) then


               do j = 1, 27


                  nall(j,ictime) = nall(j,ictime) + ncoll(j)


               end do


                  isekp = 0
                  isekd = 0
                  isekr = 0


               do i = 1, massal


                  if( inds(i) .eq. 2 ) isekd = isekd + 1
                  if( inds(i) .eq. 3 ) isekr = isekr + 1
                  if( inds(i) .eq. 4 ) isekp = isekp + 1


               end do


                  ndel(ictime) = ndel(ictime) + isekd
                  nres(ictime) = nres(ictime) + isekr
                  npin(ictime) = npin(ictime) + isekp


               do j = 1, 27


                  ncoll(j) = 0


               end do


            end if


         end if


*-----------------------------------------------------------------------
*     ground state properties
*-----------------------------------------------------------------------


      if( masspr .eq. 0 ) then


*-----------------------------------------------------------------------


            sek = 0.0


            sekx = 0.0
            seky = 0.0
            sekz = 0.0


         do i = 1, massta


            sekx = sekx + r(1,i)
            seky = seky + r(2,i)
            sekz = sekz + r(3,i)


         end do


            sekx = sekx / float(massta)
            seky = seky / float(massta)
            sekz = sekz / float(massta)


         do i = 1, massta


            rsr = ( r(1,i) - sekx )**2
     &          + ( r(2,i) - seky )**2
     &          + ( r(3,i) - sekz )**2


            if( rsr .gt. radtac )
     &           s_outp(ictime) = s_outp(ictime) + 1.0


            sek = sek + rsr




            ii = nint( (p(1,i)**2 + p(2,i)**2 + p(3,i)**2 ) * facp**2 )


            if( ii .le. jmax ) then


               denmom(ii) = denmom(ii) + 1.0
               denmtm(ictime,ii) = denmtm(ictime,ii) + 1.0


            end if


            ii = nint( r(1,i)**2 + r(2,i)**2 + r(3,i)**2 )
            if( ii .le. imax )  denpat(ii) = denpat(ii) + 1.0


         end do


            s_rmsr(ictime) = s_rmsr(ictime)
     &                     + sqrt( sek / float(massta) + 3.0 * wl )




            lmax = 20


         do k = 0, kmax - 1


            sek = 0.0


            do l = 1, lmax


                  rr0 = xyz(1,k) + ( xyz(2,k) - xyz(1,k) ) * rn()


                  the = pi * rn()
                  phi = 2. * pi * rn()


                  rrx = rr0 * sin(the) * cos(phi)
                  rry = rr0 * sin(the) * sin(phi)
                  rrz = rr0 * cos(the)


               do i = 1, massta


                  rdis2 = ( rrx - r(1,i) + sekx )**2 +
     &                    ( rry - r(2,i) + seky )**2 +
     &                    ( rrz - r(3,i) + sekz )**2


                  expa1 = - rdis2 * csw


                  if( expa1 .gt. epsx ) then


                     rh1 = exp( expa1 )


                  else


                     rh1 = 0.0


                  end if


                     sek = sek + rh1


               end do


            end do


               dentm(k) = sek / float(lmax) * ccw


         end do




               sek = 0.0


         do i = 0, kmax - 1


            if( i .eq. 0 ) then


               dmax = ( xyz(2,i+1) + xyz(1,i+1) ) / 4.0
               vv0  = 4. / 3. * pi * dmax**3


            else


               dmin =  dmax
               dmax = ( xyz(2,i)   + xyz(1,i)   ) / 4.0
     &              + ( xyz(2,i+1) + xyz(1,i+1) ) / 4.0
               vv0 = 4. / 3. * pi * ( dmax**3 - dmin**3 )


            end if


               sek = sek + dentm(i) * vv0


         end do


         do i = 0, kmax - 1


            dentm(i) = dentm(i) * float(massta) / sek


            dentim(ictime,i) = dentim(ictime,i) + dentm(i)


            denden(i) = denden(i) + dentm(i)


         end do


*-----------------------------------------------------------------------


      end if


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine sm_lstp
*                                                                      *
*                                                                      *
*        Last Revised:     1998 12 09                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to count the last pion number for summary               *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'


*-----------------------------------------------------------------------


      common /swich1/ ipot, insys, irkg, icolt
      common /vriab0/ massal, massba, nmeson


      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)


      common /coln00/ lcoll(30)


      common /summ00/ jdsp(50), idsp(50), fdsp(50), ifds(50)
      character    fdsp*80
      common /summ05/ npin(0:mmm), ndel(0:mmm), nres(0:mmm)
      common /summ06/ ictime, iltime
      common /summ12/ ncoll(30), nall(30,0:mmm)


*-----------------------------------------------------------------------


      if( jdsp(20) .eq. 0 ) return


*-----------------------------------------------------------------------


            if( icolt .ne. 1 )  return


                  isekp = 0
                  isekd = 0
                  isekr = 0


               do i = 1, massal


                  if( inds(i) .eq. 2 ) isekd = isekd + 1
                  if( inds(i) .eq. 3 ) isekr = isekr + 1
                  if( inds(i) .eq. 4 ) isekp = isekp + 1


               end do


                  ndel(ictime+1) = ndel(ictime+1) + isekd
                  nres(ictime+1) = nres(ictime+1) + isekr
                  npin(ictime+1) = npin(ictime+1) + isekp


                  nall(22,ictime+1) = nall(22,ictime+1) + lcoll(22)
                  nall(23,ictime+1) = nall(23,ictime+1) + lcoll(23)
                  nall(24,ictime+1) = nall(24,ictime+1) + lcoll(24)


*-----------------------------------------------------------------------




      return
      end




************************************************************************
*                                                                      *
      subroutine sm_coln
*                                                                      *
*                                                                      *
*        Last Revised:     1998 12 09                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to count the number of collisions for summary           *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'


*-----------------------------------------------------------------------


      common /coln00/ lcoll(30)


      common /summ00/ jdsp(50), idsp(50), fdsp(50), ifds(50)
      character    fdsp*80


      common /summ12/ ncoll(30), nall(30,0:mmm)


*-----------------------------------------------------------------------


      if( jdsp(20) .eq. 0 ) return


*-----------------------------------------------------------------------


         do i = 1, 27


            ncoll(i) = ncoll(i) + lcoll(i)


         end do


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine sm_mdis(icd)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 12 09                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to write mass and charge distribution of QMD / SDM      *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              icd         : = 30 ; QMD                                *
*                              31 ; SDM                                *
*                              32 ; fission mass distribution          *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'


*-----------------------------------------------------------------------


      common /summas/ sumas(3,0:maxpt,0:maxnt)
      common /summch/ tmas(3,0:maxpt+maxnt), cdis(3,0:maxpt)


      common /summ00/ jdsp(50), idsp(50), fdsp(50), ifds(50)
      character    fdsp*80


*-----------------------------------------------------------------------


      character       squ*1
      data squ       /"'"/


*-----------------------------------------------------------------------
*     Mass Distribution
*-----------------------------------------------------------------------


      if( jdsp(icd) .eq. 0 ) return


*-----------------------------------------------------------------------


            ida = idsp(icd)


*-----------------------------------------------------------------------
*     Mass Distribuion of QMD
*-----------------------------------------------------------------------


         if( icd .eq. 30 ) then


            kk = 1


            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Mass and Charge Distribution of QMD''/
     &            ''*'',71(''-''))')


            write(ida,'(/
     &            ''msul:{\rm\it JQMD  plotted by  \ANGEL}''/
     &            ''msdc:{\rm\Large \- \page  \-}''/
     &            ''msur:{\rm\it Mass Distribution of QMD}'')')


         else if( icd .eq. 31 ) then


            kk= 2


            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Mass and Charge Distribution of SDM''/
     &            ''*'',71(''-''))')


            write(ida,'(/
     &            ''msul:{\rm\it JQMD  plotted by  \ANGEL}''/
     &            ''msdc:{\rm\Large \- \page  \-}''/
     &            ''msur:{\rm\it Mass Distribution of SDM}'')')


         else if( icd .eq. 32 ) then


            kk= 3


            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Mass and Charge Distribution of fission''/
     &            ''*'',71(''-''))')


            write(ida,'(/
     &            ''msul:{\rm\it JQMD  plotted by  \ANGEL}''/
     &            ''msdc:{\rm\Large \- \page  \-}''/
     &            ''msur:{\rm\it Mass Distribution of fission}'')')


         end if




            write(ida,'(/
     &            ''p: port nofr'')')


*-----------------------------------------------------------------------


            do i = 0, maxpt + maxnt


               tmas(kk,i) = 0.0


            end do


            do i = 0, maxpt


               cdis(kk,i) = 0.0


            end do


*-----------------------------------------------------------------------


               itmax  = 0
               icmax  = 0
               inmax  = 0


            do iz = 0, maxpt
            do in = 0, maxnt


               if( sumas(kk,iz,in) .gt. 0.0 ) then


                  tmas(kk,iz+in) = tmas(kk,iz+in) + sumas(kk,iz,in)
                  cdis(kk,iz)    = cdis(kk,iz) + sumas(kk,iz,in)


                  if( iz + in .gt. itmax ) itmax = iz + in
                  if( iz      .gt. icmax ) icmax = iz
                  if(      in .gt. inmax ) inmax =      in


               end if


            end do
            end do


               cdis(kk,0) = cdis(kk,0) - tmas(kk,0)


*-----------------------------------------------------------------------
*        Write Mass Distribution
*-----------------------------------------------------------------------


         if( icd .eq. 30 ) then


            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Plot of Mass Distribution of QMD''/
     &            ''*'',71(''-''))')


            write(ida,'(/''p: scal(0.7) xorg(0.1) yorg(1.9)'')')
            write(ida,'( ''p: ylog'')')


            write(ida,'( /a1,''Mass Distribution of QMD''a1)') squ, squ


            write(ida,'(/''x: A_{f}'')')
            write(ida,'( ''y: Yield (mb)'')')


            write(ida,'(/''c:  mass       qmd'')')
            write(ida,'( ''h:  x-0.5    y(QMD),lh0'')')


         else if( icd .eq. 31 ) then


            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Plot of Mass Distribution of SDM''/
     &            ''*'',71(''-''))')


            write(ida,'(/''p: scal(0.7) xorg(0.1) yorg(1.9)'')')
            write(ida,'( ''p: ylog'')')


            write(ida,'( /a1,''Mass Distribution of SDM''a1)') squ, squ


            write(ida,'(/''x: A_{f}'')')
            write(ida,'( ''y: Yield (mb)'')')


            write(ida,'(/''c:  mass       sdm'')')
            write(ida,'( ''h:  x-0.5    y(SDM),lh0'')')


         else if( icd .eq. 32 ) then


            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Plot of Mass Distribution of fission''/
     &            ''*'',71(''-''))')


            write(ida,'(/''p: scal(0.7) xorg(0.1) yorg(1.9)'')')
            write(ida,'( ''p: ylog'')')


            write(ida,'( /a1,''Mass Distribution of fission''a1)')
     &                    squ, squ


            write(ida,'(/''x: A_{f}'')')
            write(ida,'( ''y: Yield (mb)'')')


            write(ida,'(/''c:  mass       fission'')')
            write(ida,'( ''h:  x-0.5    y(fission),lh0'')')


         end if




            write(ida,'(f7.1,e16.5)') 1.0, 0.0


         do i = 1, itmax + 1


            write(ida,'(f7.1,e16.5)') float(i), tmas(kk,i)


         end do




*-----------------------------------------------------------------------
*        Write Charge Distribution of QMD
*-----------------------------------------------------------------------


         if( icd .eq. 30 ) then


            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Plot of Charge Distribution of QMD''/
     &            ''*'',71(''-''))')


            write(ida,'(/''z: yorg(-1.7)'')')
            write(ida,'( ''p: ylog'')')


            write(ida,'(/a1,''Charge Distribution of QMD''a1)')
     &                    squ, squ


            write(ida,'(/''x: Z_{f}'')')
            write(ida,'( ''y: Yield (mb)'')')


            write(ida,'(/''c:  charge     qmd'')')
            write(ida,'( ''h:  x-0.5    y(QMD),lh0'')')


         else if( icd .eq. 31 ) then


            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Plot of Charge Distribution of SDM''/
     &            ''*'',71(''-''))')


            write(ida,'(/''z: yorg(-1.7)'')')
            write(ida,'( ''p: ylog'')')


            write(ida,'(/a1,''Charge Distribution of SDM''a1)')
     &                    squ, squ


            write(ida,'(/''x: Z_{f}'')')
            write(ida,'( ''y: Yield (mb)'')')


            write(ida,'(/''c:  charge     sdm'')')
            write(ida,'( ''h:  x-0.5    y(SDM),lh0'')')


         else if( icd .eq. 32 ) then


            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Plot of Charge Distribution of fission''/
     &            ''*'',71(''-''))')


            write(ida,'(/''z: yorg(-1.7)'')')
            write(ida,'( ''p: ylog'')')


            write(ida,'(/a1,''Charge Distribution of fission''a1)')
     &                    squ, squ


            write(ida,'(/''x: Z_{f}'')')
            write(ida,'( ''y: Yield (mb)'')')


            write(ida,'(/''c:  charge     fission'')')
            write(ida,'( ''h:  x-0.5    y(fission),lh0'')')


         end if




            write(ida,'(f7.1,e16.5)') 1.0, 0.0


         do i = 1, icmax + 1


            write(ida,'(f7.1,e16.5)') float(i), cdis(kk,i)


         end do


*-----------------------------------------------------------------------
*        write mass distribution in 2-Dim. of QMD
*-----------------------------------------------------------------------


         if( icd .eq. 30 ) then


            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Plot of Mass Distribution in 2-Dim. of QMD''/
     &            ''*'',71(''-''))')


            write(ida,'(/''newpage:'')')


            write(ida,'(/
     &            ''msul:{\rm\it JQMD  plotted by  \ANGEL}''/
     &            ''msdc:{\rm\Large \- \page  \-}''/
     &            ''msur:{\rm\it 2-Dim. Mass Distribution of QMD}'')')


         else if( icd .eq. 31 ) then


            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Plot of Mass Distribution in 2-Dim. of SDM''/
     &            ''*'',71(''-''))')


            write(ida,'(/''newpage:'')')


            write(ida,'(/
     &            ''msul:{\rm\it JQMD  plotted by  \ANGEL}''/
     &            ''msdc:{\rm\Large \- \page  \-}''/
     &            ''msur:{\rm\it 2-Dim. Mass Distribution of SDM}'')')


         else if( icd .eq. 32 ) then


            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Plot of Mass Distribution in 2-Dim. '',
     &            ''of fission''/
     &            ''*'',71(''-''))')


            write(ida,'(/''newpage:'')')


            write(ida,'(/
     &            ''msul:{\rm\it JQMD  plotted by  \ANGEL}''/
     &            ''msdc:{\rm\Large \- \page  \-}''/
     &            ''msur:{\rm\it 2-Dim. Mass Distribution '',
     &            ''of fission}'')')


         end if




            write(ida,'(/
     &            ''p: port nofr'')')


*-----------------------------------------------------------------------


            dxmax = float(inmax+2)
            dymax = float(icmax+2)
            dform = dymax / dxmax


*-----------------------------------------------------------------------


            write(ida,'( /''p: nosp zlog'')')
            write(ida,'(  ''p: xorg(-0.05) yorg(0.2) afac(0.6)'')')
            write(ida,'(  ''p: xmin(0) xmax('',f5.1,'')'')') dxmax
            write(ida,'(  ''p: ymin(0) ymax('',f5.1,'')'')') dymax
            write(ida,'(  ''p: form('',f7.3,'')'')') dform


            write(ida,'(/''x: Neutron Number {\it N}'')')
            write(ida,'( ''y: Proton Number {\it Z}'')')


            write(ida,'( ''hc: y = '',i3,'' to 1 by -1 ;'',
     &                      '' x = 1 to '',i3,'' by 1 ;'')')
     &                     icmax+2, inmax+2


         do i = icmax+2, 1, -1


            write(ida,'(10e11.3)')
     &           ( sumas(kk,i,l), l = 1, inmax+2 )


         end do


*-----------------------------------------------------------------------
*        write magic numbers
*-----------------------------------------------------------------------


            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Plot of Magic Numbers in 2-Dim.''/
     &            ''*'',71(''-''))')


         if( 11 .lt. inmax+2 .and. 20 .lt. icmax+2 ) then
            write(ida,'(/''w:20/y(20) x(11) iy(2) ix(3) s(0.7)'')')
         else
            write(ida,'(/''n:20/y(20) x(11) iy(2) ix(3) s(0.7)'')')
         end if
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  12  20.5'')')
            write(ida,'( ''  52  20.5'')')
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  12  19.5'')')
            write(ida,'( ''  52  19.5'')')


         if( 17 .lt. inmax+2 .and. 28 .lt. icmax+2 ) then
            write(ida,'(/''w:28/y(28) x(17) iy(2) ix(3) s(0.7)'')')
         else
            write(ida,'(/''n:28/y(28) x(17) iy(2) ix(3) s(0.7)'')')
         end if
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  18  28.5'')')
            write(ida,'( ''  84  28.5'')')
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  18  27.5'')')
            write(ida,'( ''  84  27.5'')')


         if( 47 .lt. inmax+2 .and. 50 .lt. icmax+2 ) then
            write(ida,'(/''w:50/y(50) x(47) iy(2) ix(3) s(0.7)'')')
         else
            write(ida,'(/''n:50/y(50) x(47) iy(2) ix(3) s(0.7)'')')
         end if
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  48   50.5'')')
            write(ida,'( ''  128  50.5'')')
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  48   49.5'')')
            write(ida,'( ''  128  49.5'')')


         if( 79 .lt. inmax+2 .and. 82 .lt. icmax+2 ) then
            write(ida,'(/''w:82/y(82) x(79) iy(2) ix(3) s(0.7)'')')
         else
            write(ida,'(/''n:82/y(82) x(79) iy(2) ix(3) s(0.7)'')')
         end if
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  80   82.5'')')
            write(ida,'( ''  150  82.5'')')
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  80   81.5'')')
            write(ida,'( ''  150  81.5'')')


         if( 20 .lt. inmax+2 .and.  7 .lt. icmax+2 ) then
            write(ida,'(/''w:20/x(20) y(7) ix(2) iy(3) s(0.7)'')')
         else
            write(ida,'(/''n:20/x(20) y(7) ix(2) iy(3) s(0.7)'')')
         end if
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  20.5   8'')')
            write(ida,'( ''  20.5   30'')')
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  19.5   8'')')
            write(ida,'( ''  19.5   30'')')


         if( 28 .lt. inmax+2 .and.  7 .lt. icmax+2 ) then
            write(ida,'(/''w:28/x(28) y(7) ix(2) iy(3) s(0.7)'')')
         else
            write(ida,'(/''n:28/x(28) y(7) ix(2) iy(3) s(0.7)'')')
         end if
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  28.5   8'')')
            write(ida,'( ''  28.5   36'')')
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  27.5   8'')')
            write(ida,'( ''  27.5   36'')')


         if( 50 .lt. inmax+2 .and. 17 .lt. icmax+2 ) then
            write(ida,'(/''w:50/x(50) y(17) ix(2) iy(3) s(0.7)'')')
         else
            write(ida,'(/''n:50/x(50) y(17) ix(2) iy(3) s(0.7)'')')
         end if
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  50.5   18'')')
            write(ida,'( ''  50.5   52'')')
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  49.5   18'')')
            write(ida,'( ''  49.5   52'')')


         if( 82 .lt. inmax+2 .and. 25 .lt. icmax+2 ) then
            write(ida,'(/''w:82/x(82) y(25) ix(2) iy(3) s(0.7)'')')
         else
            write(ida,'(/''n:82/x(82) y(25) ix(2) iy(3) s(0.7)'')')
         end if
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  82.5   26'')')
            write(ida,'( ''  82.5   84'')')
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  81.5   26'')')
            write(ida,'( ''  81.5   84'')')


         if( 126 .lt. inmax+2 .and. 47 .lt. icmax+2 ) then
            write(ida,'(/''w:126/x(126) y(47) ix(2) iy(3) s(0.7)'')')
         else
            write(ida,'(/''n:126/x(126) y(47) ix(2) iy(3) s(0.7)'')')
         end if
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  126.5   48'')')
            write(ida,'( ''  126.5   100'')')
            write(ida,'( ''h: x y,l0zzz'')')
            write(ida,'( ''  125.5   48'')')
            write(ida,'( ''  125.5   100'')')


*-----------------------------------------------------------------------
*        write stable nuclei
*-----------------------------------------------------------------------


            write(ida,'(/
     &            ''*'',71(''-'') /
     &            ''*     Plot of Stable Nuclei in 2-Dim.''/
     &            ''*'',71(''-''))')


         if( inmax .gt. 100 ) then


            write(ida,'( ''h: x y,n6xxxxx'')')


         else if( inmax .gt. 75 ) then


            write(ida,'( ''h: x y,n6xxxx'')')


         else if( inmax .gt. 50 ) then


            write(ida,'( ''h: x y,n6xxx'')')


         else if( inmax .gt. 30 ) then


            write(ida,'( ''h: x y,n6xx'')')


         else


            write(ida,'( ''h: x y,n6'')')


         end if


            write(ida,'( ''  0   1'')')
            write(ida,'( ''  1   1'')')
            write(ida,'( ''  1   2'')')
            write(ida,'( ''  2   2'')')
            write(ida,'( ''  3   3'')')
            write(ida,'( ''  4   3'')')
            write(ida,'( ''  5   4'')')
            write(ida,'( ''  5   5'')')
            write(ida,'( ''  6   5'')')
            write(ida,'( ''  6   6'')')
            write(ida,'( ''  7   6'')')
            write(ida,'( ''  7   7'')')
            write(ida,'( ''  8   7'')')
            write(ida,'( ''  8   8'')')
            write(ida,'( ''  9   8'')')
            write(ida,'( ''  10  8'')')
            write(ida,'( ''  10  9'')')
            write(ida,'( ''  10  10'')')
            write(ida,'( ''  11  10'')')
            write(ida,'( ''  12  10'')')
            write(ida,'( ''  12  11'')')
            write(ida,'( ''  12  12'')')
            write(ida,'( ''  13  12'')')
            write(ida,'( ''  14  12'')')
            write(ida,'( ''  14  13'')')
            write(ida,'( ''  14  14'')')
            write(ida,'( ''  15  14'')')
            write(ida,'( ''  16  14'')')
            write(ida,'( ''  16  15'')')
            write(ida,'( ''  16  16'')')
            write(ida,'( ''  17  16'')')
            write(ida,'( ''  18  16'')')
            write(ida,'( ''  18  17'')')
            write(ida,'( ''  18  18'')')
            write(ida,'( ''  20  16'')')
            write(ida,'( ''  20  17'')')
            write(ida,'( ''  20  18'')')
            write(ida,'( ''  20  19'')')
            write(ida,'( ''  20  20'')')
            write(ida,'( ''  21  19'')')
            write(ida,'( ''  22  18'')')
            write(ida,'( ''  22  19'')')
            write(ida,'( ''  22  20'')')
            write(ida,'( ''  23  20'')')
            write(ida,'( ''  24  20'')')
            write(ida,'( ''  24  21'')')
            write(ida,'( ''  24  22'')')
            write(ida,'( ''  25  22'')')
            write(ida,'( ''  26  20'')')
            write(ida,'( ''  26  22'')')
            write(ida,'( ''  26  24'')')
            write(ida,'( ''  27  22'')')
            write(ida,'( ''  27  23'')')
            write(ida,'( ''  28  20'')')
            write(ida,'( ''  28  22'')')
            write(ida,'( ''  28  23'')')
            write(ida,'( ''  28  24'')')
            write(ida,'( ''  28  26'')')
            write(ida,'( ''  29  24'')')
            write(ida,'( ''  30  24'')')
            write(ida,'( ''  30  25'')')
            write(ida,'( ''  30  26'')')
            write(ida,'( ''  30  28'')')
            write(ida,'( ''  31  26'')')
            write(ida,'( ''  32  26'')')
            write(ida,'( ''  32  27'')')
            write(ida,'( ''  32  28'')')
            write(ida,'( ''  33  28'')')
            write(ida,'( ''  34  28'')')
            write(ida,'( ''  34  29'')')
            write(ida,'( ''  34  30'')')
            write(ida,'( ''  36  28'')')
            write(ida,'( ''  36  29'')')
            write(ida,'( ''  36  30'')')
            write(ida,'( ''  37  30'')')
            write(ida,'( ''  38  30'')')
            write(ida,'( ''  38  31'')')
            write(ida,'( ''  38  32'')')
            write(ida,'( ''  40  30'')')
            write(ida,'( ''  40  31'')')
            write(ida,'( ''  40  32'')')
            write(ida,'( ''  40  34'')')
            write(ida,'( ''  41  32'')')
            write(ida,'( ''  42  32'')')
            write(ida,'( ''  42  33'')')
            write(ida,'( ''  42  34'')')
            write(ida,'( ''  42  36'')')
            write(ida,'( ''  43  34'')')
            write(ida,'( ''  44  32'')')
            write(ida,'( ''  44  34'')')
            write(ida,'( ''  44  35'')')
            write(ida,'( ''  44  36'')')
            write(ida,'( ''  46  34'')')
            write(ida,'( ''  46  35'')')
            write(ida,'( ''  46  36'')')
            write(ida,'( ''  46  38'')')
            write(ida,'( ''  47  36'')')
            write(ida,'( ''  48  34'')')
            write(ida,'( ''  48  36'')')
            write(ida,'( ''  48  37'')')
            write(ida,'( ''  48  38'')')
            write(ida,'( ''  49  38'')')
            write(ida,'( ''  50  36'')')
            write(ida,'( ''  50  37'')')
            write(ida,'( ''  50  38'')')
            write(ida,'( ''  50  39'')')
            write(ida,'( ''  50  40'')')
            write(ida,'( ''  50  42'')')
            write(ida,'( ''  51  40'')')
            write(ida,'( ''  52  40'')')
            write(ida,'( ''  52  41'')')
            write(ida,'( ''  52  42'')')
            write(ida,'( ''  52  44'')')
            write(ida,'( ''  53  42'')')
            write(ida,'( ''  54  40'')')
            write(ida,'( ''  54  42'')')
            write(ida,'( ''  54  44'')')
            write(ida,'( ''  55  42'')')
            write(ida,'( ''  55  44'')')
            write(ida,'( ''  56  40'')')
            write(ida,'( ''  56  42'')')
            write(ida,'( ''  56  44'')')
            write(ida,'( ''  56  46'')')
            write(ida,'( ''  57  44'')')
            write(ida,'( ''  58  42'')')
            write(ida,'( ''  58  44'')')
            write(ida,'( ''  58  45'')')
            write(ida,'( ''  58  46'')')
            write(ida,'( ''  58  48'')')
            write(ida,'( ''  59  46'')')
            write(ida,'( ''  60  44'')')
            write(ida,'( ''  60  46'')')
            write(ida,'( ''  60  47'')')
            write(ida,'( ''  60  48'')')
            write(ida,'( ''  62  46'')')
            write(ida,'( ''  62  47'')')
            write(ida,'( ''  62  48'')')
            write(ida,'( ''  62  50'')')
            write(ida,'( ''  63  47'')')
            write(ida,'( ''  64  46'')')
            write(ida,'( ''  64  48'')')
            write(ida,'( ''  64  49'')')
            write(ida,'( ''  64  50'')')
            write(ida,'( ''  65  48'')')
            write(ida,'( ''  65  50'')')
            write(ida,'( ''  66  48'')')
            write(ida,'( ''  66  49'')')
            write(ida,'( ''  66  50'')')
            write(ida,'( ''  67  50'')')
            write(ida,'( ''  68  48'')')
            write(ida,'( ''  68  50'')')
            write(ida,'( ''  68  52'')')
            write(ida,'( ''  69  50'')')
            write(ida,'( ''  70  50'')')
            write(ida,'( ''  70  51'')')
            write(ida,'( ''  70  52'')')
            write(ida,'( ''  70  54'')')
            write(ida,'( ''  71  52'')')
            write(ida,'( ''  72  51'')')
            write(ida,'( ''  72  52'')')
            write(ida,'( ''  72  53'')')
            write(ida,'( ''  72  54'')')
            write(ida,'( ''  73  52'')')
            write(ida,'( ''  74  50'')')
            write(ida,'( ''  74  52'')')
            write(ida,'( ''  74  53'')')
            write(ida,'( ''  74  54'')')
            write(ida,'( ''  74  56'')')
            write(ida,'( ''  75  54'')')
            write(ida,'( ''  76  52'')')
            write(ida,'( ''  76  54'')')
            write(ida,'( ''  76  56'')')
            write(ida,'( ''  77  54'')')
            write(ida,'( ''  78  52'')')
            write(ida,'( ''  78  54'')')
            write(ida,'( ''  78  55'')')
            write(ida,'( ''  78  56'')')
            write(ida,'( ''  78  58'')')
            write(ida,'( ''  79  56'')')
            write(ida,'( ''  80  54'')')
            write(ida,'( ''  80  56'')')
            write(ida,'( ''  80  58'')')
            write(ida,'( ''  81  56'')')
            write(ida,'( ''  81  57'')')
            write(ida,'( ''  82  54'')')
            write(ida,'( ''  82  56'')')
            write(ida,'( ''  82  57'')')
            write(ida,'( ''  82  58'')')
            write(ida,'( ''  82  59'')')
            write(ida,'( ''  82  60'')')
            write(ida,'( ''  82  62'')')
            write(ida,'( ''  83  60'')')
            write(ida,'( ''  84  58'')')
            write(ida,'( ''  84  60'')')
            write(ida,'( ''  85  60'')')
            write(ida,'( ''  85  62'')')
            write(ida,'( ''  86  60'')')
            write(ida,'( ''  86  62'')')
            write(ida,'( ''  87  62'')')
            write(ida,'( ''  88  60'')')
            write(ida,'( ''  88  62'')')
            write(ida,'( ''  88  63'')')
            write(ida,'( ''  88  64'')')
            write(ida,'( ''  90  60'')')
            write(ida,'( ''  90  62'')')
            write(ida,'( ''  90  63'')')
            write(ida,'( ''  90  64'')')
            write(ida,'( ''  90  66'')')
            write(ida,'( ''  91  64'')')
            write(ida,'( ''  92  62'')')
            write(ida,'( ''  92  64'')')
            write(ida,'( ''  92  66'')')
            write(ida,'( ''  93  64'')')
            write(ida,'( ''  94  64'')')
            write(ida,'( ''  94  65'')')
            write(ida,'( ''  94  66'')')
            write(ida,'( ''  94  68'')')
            write(ida,'( ''  95  66'')')
            write(ida,'( ''  96  64'')')
            write(ida,'( ''  96  66'')')
            write(ida,'( ''  96  68'')')
            write(ida,'( ''  97  66'')')
            write(ida,'( ''  98  66'')')
            write(ida,'( ''  98  67'')')
            write(ida,'( ''  98  68'')')
            write(ida,'( ''  98  70'')')
            write(ida,'( ''  99  68'')')
            write(ida,'( ''  100 68'')')
            write(ida,'( ''  100 69'')')
            write(ida,'( ''  100 70'')')
            write(ida,'( ''  101 70'')')
            write(ida,'( ''  102 68'')')
            write(ida,'( ''  102 70'')')
            write(ida,'( ''  102 72'')')
            write(ida,'( ''  103 70'')')
            write(ida,'( ''  104 70'')')
            write(ida,'( ''  104 71'')')
            write(ida,'( ''  105 71'')')
            write(ida,'( ''  105 72'')')
            write(ida,'( ''  106 70'')')
            write(ida,'( ''  106 72'')')
            write(ida,'( ''  106 74'')')
            write(ida,'( ''  107 72'')')
            write(ida,'( ''  107 73'')')
            write(ida,'( ''  108 72'')')
            write(ida,'( ''  108 73'')')
            write(ida,'( ''  108 74'')')
            write(ida,'( ''  108 76'')')
            write(ida,'( ''  109 74'')')
            write(ida,'( ''  110 74'')')
            write(ida,'( ''  110 75'')')
            write(ida,'( ''  110 76'')')
            write(ida,'( ''  111 76'')')
            write(ida,'( ''  112 74'')')
            write(ida,'( ''  112 75'')')
            write(ida,'( ''  112 76'')')
            write(ida,'( ''  112 78'')')
            write(ida,'( ''  113 76'')')
            write(ida,'( ''  114 76'')')
            write(ida,'( ''  114 77'')')
            write(ida,'( ''  114 78'')')
            write(ida,'( ''  116 76'')')
            write(ida,'( ''  116 77'')')
            write(ida,'( ''  116 78'')')
            write(ida,'( ''  116 80'')')
            write(ida,'( ''  117 78'')')
            write(ida,'( ''  118 78'')')
            write(ida,'( ''  118 79'')')
            write(ida,'( ''  118 80'')')
            write(ida,'( ''  119 80'')')
            write(ida,'( ''  120 78'')')
            write(ida,'( ''  120 80'')')
            write(ida,'( ''  121 80'')')
            write(ida,'( ''  122 80'')')
            write(ida,'( ''  122 81'')')
            write(ida,'( ''  122 82'')')
            write(ida,'( ''  124 80'')')
            write(ida,'( ''  124 81'')')
            write(ida,'( ''  124 82'')')
            write(ida,'( ''  125 82'')')
            write(ida,'( ''  126 82'')')
            write(ida,'( ''  126 83'')')


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine sm_mass(iss,iz1,in1,sfact)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 12 09                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to store the mass distribution of the clusters          *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              iss            =1 ; after QMD                           *
*                             =2 ; after SDM                           *
*                             =3 ; before fission                      *
*                                                                      *
*              iz1            : proton number                          *
*              in1            : neutron number                         *
*                                                                      *
*              sfact          : weight facter                          *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'


      common /summas/ sumas(3,0:maxpt,0:maxnt)


*-----------------------------------------------------------------------


               ipro = iz1
               ineu = in1


            if( ipro .le. maxpt .and. ineu .le. maxnt ) then


               sumas(iss,ipro,ineu) = sumas(iss,ipro,ineu) + sfact


            end if




*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine sm_norm
*                                                                      *
*                                                                      *
*        Last Revised:     1998 12 09                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*             to normarize the cross section by the total event number *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'


*-----------------------------------------------------------------------


      common /const2/ dt, ntmax, iprun, iprun0
      common /summ03/ ireac(0:10), rcross(4)
      common /summas/ sumas(3,0:maxpt,0:maxnt)


*-----------------------------------------------------------------------
*        normalization of total event number : iprun
*-----------------------------------------------------------------------


                  rmn  = float(iprun)


                  rcross(2) = rcross(2) / rmn
                  rcross(3) = rcross(3) / rmn
                  rcross(4) = rcross(4) / rmn


               do k = 1, 3
               do i = 0, maxpt
               do j = 0, maxnt


                  sumas(k,i,j) = sumas(k,i,j) / rmn


               end do
               end do
               end do


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine sm_clse
*                                                                      *
*                                                                      *
*        Last Revised:     1998 12 09                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*             to close files for summary                               *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'


*-----------------------------------------------------------------------


      common /fnamid/ fname(mxfnm), idf(mxfnm), ifnl(mxfnm)
      character       fname*80


      common /summ00/ jdsp(50), idsp(50), fdsp(50), ifds(50)
      character    fdsp*80


      common /swich4/ ifin, ifout


*-----------------------------------------------------------------------
*        close the files
*-----------------------------------------------------------------------


         do i = 1, 13


            if( jdsp(i) .ne. 0 ) then


               close(idsp(i))


            end if


         end do


         do i = 20, 23


            if( jdsp(i) .ne. 0 ) then


               close(idsp(i))


            end if


         end do


         do i = 30, 33


            if( jdsp(i) .ne. 0 ) then


               close(idsp(i))


            end if


         end do


         if( ifout .gt. 0 ) then


               close(idf(10))


         else if( ifin .gt. 0 ) then


            do i = 1, ifin


               close(idf(i+10))


            end do


         end if


*-----------------------------------------------------------------------


      return
      end




