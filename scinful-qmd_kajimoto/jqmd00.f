************************************************************************
*                                                                      *
*        PART 1: Main steering routines                                *
*                                                                      *
*                                                                      *
*   List of subprograms in rough order of relevance with main purpose  *
*      ( s = subroutine, f = function, b = block data, e = entry )     *
*                                                                      *
*                                                                      *
*  s  jqmd00    to control total JQMD code                             *
*  s  qmdevent  to simulate one event by QMD                           *
*  s  finsumry  to summarize the overall events                        *
*  s  qmdjudge  to determine the weight of one event of QMD and        *
*               to judge the elastic or inelastic reaction type        *
*  s  qmdsum    to summarize one event of QMD                          *
*  s  sdmsum    to summarize one event of SDM                          *
*  s  qmdreadd  to read one event result of QMD from the file          *
*  s  qmdwrite  to write one event result of QMD on the file           *
*  s  qmdreadc  to read and check QMD results from the files           *
*                                                                      *
*                                                                      *
************************************************************************


************************************************************************
*                                                                      *
                          subroutine jqmd00
*                                                                      *
*                                                                      *
*                 JJJ     QQQQQ     MMM     MMM  DDDDDDD               *
*                 JJJ   QQQQQQQQQ   MMMM   MMMM  DDDDDDDD              *
*                 JJJ  QQQ     QQQ  MMMMM MMMMM  DDD   DDD             *
*                 JJJ  QQQ     QQQ  MMM MMM MMM  DDD    DDD            *
*                 JJJ  QQQ     QQQ  MMM  M  MMM  DDD    DDD            *
*                 JJJ  QQQ     QQQ  MMM     MMM  DDD    DDD            *
*                 JJJ  QQQ QQQ QQQ  MMM     MMM  DDD    DDD            *
*          JJJ    JJJ   QQQQQQQQQ   MMM     MMM  DDD   DDD             *
*           JJJJJJJJ      QQQQQQ    MMM     MMM  DDDDDDDD              *
*             JJJJ            QQQ   MMM     MMM  DDDDDDD               *
*                                                                      *
*                                                                      *
*                 Jaeri Quantum Molecular Dynamics                     *
*                      for Nuclear Reactions                           *
*                                                                      *
                    parameter ( Version = 1.00 )
*                                                                      *
*                            made by                                   *
*                                                                      *
*                Japan Atomic Energy Research Institute                *
*                                                                      *
                 parameter ( Last Revised = 1998 12 03 )
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'


*-----------------------------------------------------------------------


      common /const2/ dt, ntmax, iprun, iprun0
      common /vriab1/ b, llnow, ntnow
      common /swich4/ ifin, ifout
      common /sdmsw0/ issdm


*-----------------------------------------------------------------------
*     Initialization
*-----------------------------------------------------------------------


                  call jqmdver( Version, Last Revised )
                  call readcfg
                  call qmdint
                  call sdmint
                  call anal_int
                  call sm_init


*-----------------------------------------------------------------------
*     Many Events
*-----------------------------------------------------------------------


         do 5000 ll = 1, iprun


                  llnow = ll


                  if(mod(ll,1000).eq.0) write(*,*) 'event=', ll


*-----------------------------------------------------------------------
*           one event by QMD or from Files
*-----------------------------------------------------------------------


               if( ifin .eq. 0 ) then


                  call qmdevent


               else if( ifin .gt. 0 ) then


                  call qmdreadd


               end if


*-----------------------------------------------------------------------
*           write one event on file
*-----------------------------------------------------------------------


               if( ifout .gt. 0 ) then


                  call qmdwrite


               end if


*-----------------------------------------------------------------------
*           one event weight and judgement of inelastic reaction
*           and summary of QMD
*-----------------------------------------------------------------------


                  call qmdjudge
                  call qmdsum


*-----------------------------------------------------------------------
*           SDM (statistical decay) of clusters and summary of SDM
*-----------------------------------------------------------------------


            if( issdm .gt. 0 ) then


               do iss = 1, issdm


                  call sdment
                  call sdmsum


               end do


            end if


*-----------------------------------------------------------------------


 5000    continue


*-----------------------------------------------------------------------
*        Final Summary
*-----------------------------------------------------------------------


                  call finsumry


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine qmdevent
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 26                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to simulate one event by QMD                            *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      common /const1/ elab, rdist, bmin, bmax, ibch, ibin
      common /const2/ dt, ntmax, iprun, iprun0
      common /const3/ nfreq, nfrec, nfred
      common /vriab1/ b, llnow, ntnow
      common /swich1/ ipot, insys, irkg, icolt


      common /impact/ bval(1000), ibnum(1000), bweight(1000), bdef


      common /rannum/ iseed, iseed0, iseed1
      common /coln01/ iccoll


      data ibsf /0/
      save ibsf


*-----------------------------------------------------------------------
*        Initialization of one QMD event
*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*           event number, initial time, and initial randum seed
*           total collision flag
*-----------------------------------------------------------------------


                  ntnow  = 0
                  iseed1 = iseed
                  iccoll = 0


*-----------------------------------------------------------------------
*           impact parameter and multi run control
*-----------------------------------------------------------------------


               if( ibch .eq. 0 ) then


                     call howmany(1)


                     b = sqrt( max( 0.0,
     &                   bmin**2 + ( bmax**2 - bmin**2 ) * rn() ) )




               else if( ibch .eq. 1 ) then


                     ibsf = ibsf + 1


                  if( ibsf .gt. ibin ) then


                     call howmany(1)


                     ibsf = 1


                  end if


                     b = bval(ibsf) - bdef / 2.0 + bdef * rn()


               end if


*-----------------------------------------------------------------------
*           make ground state and boost
*-----------------------------------------------------------------------


                  call cputime(4)


                  call ground
                  call rboost


                  call cputime(4)


*-----------------------------------------------------------------------
*           Summary for the initial time
*-----------------------------------------------------------------------


                  call cputime(5)


                  call sm_evnt
                  call sm_timq
                  call sm_timc
                  call sm_timd


                  call cputime(5)




*-----------------------------------------------------------------------
*        Time Evolution
*-----------------------------------------------------------------------


         do 100 nt = 1, ntmax


               ntnow = nt


*-----------------------------------------------------------------------
*           time integration
*-----------------------------------------------------------------------


                  call cputime(2)


               if( irkg .eq. 2 ) then


                  call rk12(dt)


               else if( irkg .eq. 4 ) then


                  call rkg4(dt)


               end if


                  call cputime(2)


*-----------------------------------------------------------------------
*           collision term
*-----------------------------------------------------------------------


               if( icolt .eq. 1 ) then


                  call cputime(3)


                  call pionem(dt)
                  call relcol
                  call pionab
                  call sm_coln


                  call cputime(3)


               end if


*-----------------------------------------------------------------------
*           summary at the certain time interval
*-----------------------------------------------------------------------


                  call cputime(5)


               if( ( ntnow / nfreq ) * nfreq .eq. ntnow ) then


                  call sm_timq


               end if


               if( ( ntnow / nfrec ) * nfrec .eq. ntnow ) then


                  call sm_timc


               end if


               if( ( ntnow / nfred ) * nfred .eq. ntnow ) then


                  call sm_timd


               end if


                  call cputime(5)


*-----------------------------------------------------------------------


  100    continue


*-----------------------------------------------------------------------
*        Final pion decay and Final analysis of clusters
*-----------------------------------------------------------------------


                  call cputime(6)


                  call fpidecay
                  call cldist


                  call cputime(6)


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine finsumry
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 26                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to summarize the overall runs                           *
*                                                                      *
************************************************************************


      common /startt/ iday1,imon1,iyer1,ihor1,imin1,isec1


*-----------------------------------------------------------------------
*        normal summary
*-----------------------------------------------------------------------


                  call cputime(5)


                  call sm_norm


                  call sm_grdo
                  call sm_coll
                  call sm_mdis(30)
                  call sm_mdis(31)
                  call sm_mdis(32)


                  call cputime(5)


*-----------------------------------------------------------------------
*        Final Total CPU time, Date and Time
*-----------------------------------------------------------------------


                  call datetime(iyer1,imon1,iday1,ihor1,imin1,isec1)
                  call cputime(1)


*-----------------------------------------------------------------------
*              final summary
*-----------------------------------------------------------------------


                  call howmany(1)


                  call sm_hedo
                  call sm_cpuo
                  call sm_sumo


                  call sm_clse


*-----------------------------------------------------------------------
*        call user subroutine for summary
*-----------------------------------------------------------------------


                  call anal_fin


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine qmdjudge
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 26                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to determin weight of the event, [qmdfac] and           *
*              to judge the elastic or inelastic reaction type and     *
*              to sum up the reaction cross section                    *
*                                                                      *
*                                                                      *
*        Variables: in common block /swich3/                           *
*                                                                      *
*              ielst       : input flag, jelst < ielst: elastic        *
*              jelst       : elastic or inelastic flag                 *
*                     = 0  : elastic without collision                 *
*                     = 1  : elastic with collision                    *
*                     = 2  : inelastic without collision               *
*                     = 3  : inelastic with collision                  *
*              kelst       : 0 -> elastic, 1-> inelastic               *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'


*-----------------------------------------------------------------------


      common /const0/ idnta, idnpr, massta, masspr, mstapr, msprpr
      common /const1/ elab, rdist, bmin, bmax, ibch, ibin


      common /vriab3/ qmdfac, sdmfac


      common /clustf/ nclst, iclust(nnn)
      common /clustg/ jclust(0:7,nnn),  qclust(0:7,nnn)


      common /coln01/ iccoll
      common /swich3/ ielst, jelst, kelst
      common /sdmcut/ sdmemin


      common /impact/ bval(1000), ibnum(1000), bweight(1000), bdef


      common /summ03/ ireac(0:10), rcross(4)


*-----------------------------------------------------------------------


               call cputime(6)


*-----------------------------------------------------------------------
*        weight for QMD
*-----------------------------------------------------------------------


               bimp = qclust(0,1)


               if( bdef .ne. 0.0 ) then
                  ib = int( ( bimp - bmin ) / bdef ) + 1
               else
                  ib = 1
               end if


               ib = max( ib, 1 )
               ib = min( ib, ibin )


               qmdfac = bweight(ib)


*-----------------------------------------------------------------------
*        judgement of elastic or inelastic reaction
*-----------------------------------------------------------------------


                     iels = 1


            if( nclst .eq. 2 ) then


               if( ( jclust(1,1) .eq. mstapr .and.
     &               jclust(2,1) .eq. massta - mstapr .and.
     &               jclust(1,2) .eq. msprpr .and.
     &               jclust(2,2) .eq. masspr - msprpr )     .or.
     &             ( jclust(1,2) .eq. mstapr .and.
     &               jclust(2,2) .eq. massta - mstapr .and.
     &               jclust(1,1) .eq. msprpr .and.
     &               jclust(2,1) .eq. masspr - msprpr ) )   then


                  if( qclust(6,1) .lt. sdmemin .and.
     &                qclust(6,2) .lt. sdmemin  ) then


                     iels = 0


                  end if


               end if


            end if


            if( iccoll .eq. 0 ) then


               if( iels .eq. 0 ) then


                  jelst = 0


               else


                  jelst = 2


               end if


            else


               if( iels .eq. 0 ) then


                  jelst = 1


               else


                  jelst = 3


               end if


            end if


            ireac(jelst) = ireac(jelst) + 1


*-----------------------------------------------------------------------
*        elastic flag
*-----------------------------------------------------------------------


            if( jelst .ge. ielst ) then


               kelst = 1


            else


               kelst = 0


            end if


*-----------------------------------------------------------------------
*        reaction cross section
*-----------------------------------------------------------------------


               rcross(2) = rcross(2) + qmdfac


               if( kelst .eq. 1 )
     &         rcross(3) = rcross(3) + qmdfac


*-----------------------------------------------------------------------


               call cputime(6)


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine qmdsum
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 26                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to summarize one QMD event:                             *
*              to store the mass distribution of the clusters          *
*                 after QMD calculation (inelastic only)               *
*              to call user program for summary of QMD event           *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'


      common /vriab3/ qmdfac, sdmfac
      common /swich3/ ielst, jelst, kelst


      common /clustf/ nclst, iclust(nnn)
      common /clustg/ jclust(0:7,nnn),  qclust(0:7,nnn)


*-----------------------------------------------------------------------


               call cputime(5)


*-----------------------------------------------------------------------
*        do loop for the clusters and particles of QMD
*-----------------------------------------------------------------------


         do i = 1, nclst


                  ik = iclust(i)


                  jj = jclust(0,i)
                  iz = jclust(1,i)
                  in = jclust(2,i)
                  id = jclust(3,i)
                  is = jclust(4,i)
                  ic = jclust(5,i)
                  iq = jclust(6,i)
                  im = jclust(7,i)


                  bi = qclust(0,i)
                  px = qclust(1,i)
                  py = qclust(2,i)
                  pz = qclust(3,i)
                  et = qclust(4,i)
                  rm = qclust(5,i)
                  ex = qclust(6,i)


*-----------------------------------------------------------------------
*        store the mass distribution of the clusters
*        after QMD calculation only for inelastic reaction
*-----------------------------------------------------------------------


               if( kelst .eq. 1 .and. ik .le. 2 ) then


                  call sm_mass(1,iz,in,qmdfac)


               end if


*-----------------------------------------------------------------------
*        call user subroutine for summary of QMD
*-----------------------------------------------------------------------


               call anal_qmd(ik,jj,iz,in,id,is,ic,iq,im,
     &                       bi,px,py,pz,et,rm,ex)


*-----------------------------------------------------------------------


         end do


*-----------------------------------------------------------------------


               call cputime(5)


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine sdmsum
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 26                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to summarize one SDM event:                             *
*              to store the mass distribution of the clusters          *
*                 after SDM calculation (inelastic only)               *
*              to call user program for summary of SDM event           *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'


      common /vriab3/ qmdfac, sdmfac
      common /swich3/ ielst, jelst, kelst


      common /clustt/ nclsts, iclusts(nnn)
      common /clustw/ jclusts(0:7,nnn),  qclusts(0:7,nnn)


*-----------------------------------------------------------------------


               call cputime(5)


*-----------------------------------------------------------------------
*        do loop for the clusters and particles of SDM
*-----------------------------------------------------------------------


         do i = 1, nclsts


                  ik = iclusts(i)


                  jj = jclusts(0,i)
                  iz = jclusts(1,i)
                  in = jclusts(2,i)
                  id = jclusts(3,i)
                  is = jclusts(4,i)
                  ic = jclusts(5,i)
                  iq = jclusts(6,i)
                  im = jclusts(7,i)


                  bi = qclusts(0,i)
                  px = qclusts(1,i)
                  py = qclusts(2,i)
                  pz = qclusts(3,i)
                  et = qclusts(4,i)
                  rm = qclusts(5,i)
                  ex = qclusts(6,i)


*-----------------------------------------------------------------------
*        store the mass distribution of the clusters
*        after SDM calculation
*-----------------------------------------------------------------------


               if( kelst .eq. 1 .and. ik .le. 2 ) then


                  call sm_mass(2,iz,in, qmdfac * sdmfac )


               end if


*-----------------------------------------------------------------------
*        call user subroutine for summary of SDM
*-----------------------------------------------------------------------


               call anal_sdm(ik,jj,iz,in,id,is,ic,iq,im,
     &                       bi,px,py,pz,et,rm,ex)


*-----------------------------------------------------------------------


         end do


*-----------------------------------------------------------------------


               call cputime(5)


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine qmdreadd
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 26                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to read one event result of QMD from files              *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'


      common /const1/ elab, rdist, bmin, bmax, ibch, ibin
      common /swich4/ ifin, ifout
      common /coln01/ iccoll


      common /clustf/ nclst, iclust(nnn)
      common /clustg/ jclust(0:7,nnn),  qclust(0:7,nnn)


      common /fnamid/ fname(mxfnm), idf(mxfnm), ifnl(mxfnm)
      character       fname*80
      common /rannum/ iseed, iseed0, iseed1


      data jfile /10/
      data kfile /0/


      save jfile, kfile


      character dhead*20


*-----------------------------------------------------------------------
*        output unit for error message
*-----------------------------------------------------------------------


               ieo = 6


*-----------------------------------------------------------------------
*        save initial random number
*-----------------------------------------------------------------------


               iseed1 = iseed


*-----------------------------------------------------------------------


               call cputime(5)


*-----------------------------------------------------------------------
*        initial time of the file
*-----------------------------------------------------------------------


  600    continue


         if( kfile .eq. 0 ) then


               jfile = jfile + 1


               io = idf(jfile)


  300          read(io,'(a20)') dhead


               if( dhead .ne. '  ik  iz  in  ic  jj' ) goto 300


               kfile = 1


         end if


*-----------------------------------------------------------------------
*        read data
*-----------------------------------------------------------------------


         if( kfile .eq. 1 ) then


               read(io,'(9i4,1p6e14.6)', end = 500 )
     &            ik, iz, in, ic, jj, id, is, iq, im,
     &            px, py, pz, et, rm, ex


               goto 700


  500       continue


               kfile = 0


               goto 600


  700       continue


*-----------------------------------------------------------------------


                  nclst  = jj
                  iccoll = id


                  bimp        = px
                  qclust(0,1) = px


            do i = 1, nclst


               read(io,'(9i4,1p6e14.6)')
     &            ik, iz, in, ic, jj, id, is, iq, im,
     &            px, py, pz, et, rm, ex


                  iclust(i)   = ik


                  jclust(0,i) = jj
                  jclust(1,i) = iz
                  jclust(2,i) = in
                  jclust(3,i) = id
                  jclust(4,i) = is
                  jclust(5,i) = ic
                  jclust(6,i) = iq
                  jclust(7,i) = im


                  qclust(0,i) = qclust(0,1)
                  qclust(1,i) = px
                  qclust(2,i) = py
                  qclust(3,i) = pz
                  qclust(4,i) = et
                  qclust(5,i) = rm
                  qclust(6,i) = ex


            end do


            if( bimp .gt. bmax .or. bimp .lt. bmin ) goto 600


         end if


*-----------------------------------------------------------------------


               call cputime(5)


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine qmdwrite
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 26                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to write one event result of QMD on file                *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'


      common /coln01/ iccoll


      common /clustf/ nclst, iclust(nnn)
      common /clustg/ jclust(0:7,nnn),  qclust(0:7,nnn)


      common /fnamid/ fname(mxfnm), idf(mxfnm), ifnl(mxfnm)
      character       fname*80


      data  zero /0.0/
      data izero /0/


*-----------------------------------------------------------------------
*        event header
*-----------------------------------------------------------------------


               call cputime(5)


               io = idf(10)


               bimp = qclust(0,1)


               icolnm = iccoll
               if( icolnm .gt. 9999 ) icolnm = 9999


               write(io,'(9i4,1p6e14.6)')
     &            izero, izero, izero, izero,
     &            nclst, icolnm, izero, izero, izero,
     &            bimp, zero, zero, zero, zero, zero


*-----------------------------------------------------------------------
*        do loop for the clusters and particles of QMD
*-----------------------------------------------------------------------


         do i = 1, nclst


                  ik = iclust(i)


                  jj = jclust(0,i)


                  iz = jclust(1,i)
                  in = jclust(2,i)


                  id = jclust(3,i)
                  is = jclust(4,i)
                  ic = jclust(5,i)
                  iq = jclust(6,i)


                  im = 0


                  px = qclust(1,i)
                  py = qclust(2,i)
                  pz = qclust(3,i)


                  et = qclust(4,i)
                  rm = qclust(5,i)


                  ex = qclust(6,i)


               write(io,'(9i4,1p6e14.6)')
     &            ik, iz, in, ic, jj, id, is, iq, im,
     &            px, py, pz, et, rm, ex




         end do


*-----------------------------------------------------------------------


               call cputime(5)


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine qmdreadc
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 26                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to read and check QMD results from files                *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'
      include 'param02.inc'


      common /const1/ elab, rdist, bmin, bmax, ibch, ibin
      common /const2/ dt, ntmax, iprun, iprun0


      common /vriab2/ icevnt, idevnt


      common /swich4/ ifin, ifout


      common /impact/ bval(1000), ibnum(1000), bweight(1000), bdef


      common /fnamid/ fname(mxfnm), idf(mxfnm), ifnl(mxfnm)
      character       fname*80


      data jfile /10/
      data kfile /0/


      save jfile, kfile


      character dhead*20


*-----------------------------------------------------------------------
*        output unit for error message
*-----------------------------------------------------------------------


               ieo = 6


*-----------------------------------------------------------------------


               call cputime(5)


*-----------------------------------------------------------------------
*     initial time of the file
*-----------------------------------------------------------------------


               idevnt = 0


               ll = 0


*-----------------------------------------------------------------------


 5000 continue


               ll = ll + 1


*-----------------------------------------------------------------------


  600    continue


         if( kfile .eq. 0 ) then


               jfile = jfile + 1


            if( jfile - 10 .gt. ifin ) then


               ll = ll - 1


               goto 1000


            end if


               io = idf(jfile)


  300          read(io,'(a20)', end = 200 ) dhead


               if( dhead .ne. '  ik  iz  in  ic  jj' ) goto 300


               goto 100


  200       continue


                  write(ieo,'('' Error: there is no data'',
     &                        '' or no header in input data file:'',
     &                        '' filename = '',80a1)')
     &                        (fname(jfile)(k:k), k=1,ifnl(jfile))


                  write(ieo,'('' ====='')')


                  stop


  100          kfile = 1


         end if


*-----------------------------------------------------------------------
*        read data
*-----------------------------------------------------------------------


         if( kfile .eq. 1 ) then


               read(io,'(9i4,1p6e14.6)', end = 500 )
     &            ik, iz, in, ic, jj, id, is, iq, im,
     &            px, py, pz, et, rm, ex


            if( ik .ne. 0 .or.
     &          iz .ne. 0 .or.
     &          in .ne. 0 .or.
     &          ic .ne. 0 ) then


               write(ieo,'(/'' **** Reading data is terminated,''/
     &                      '' event header is missing at event = '',
     &                   i6/'' filename = '',80a1)') ll,
     &                        (fname(jfile)(k:k), k=1,ifnl(jfile))


               stop


            end if


               goto 700


  500       continue


               kfile = 0


               goto 600


  700       continue


*-----------------------------------------------------------------------


                  bimp = px


            if( bimp .le. bmax .and. bimp .ge. bmin ) then


                  idevnt = idevnt + 1


               if( ibch .eq. 2 .and. idevnt .le. iprun ) then


                  if( bdef .ne. 0.0 ) then
                     ib = int( ( bimp - bmin ) / bdef ) + 1
                  else
                     ib = 1
                  end if


                  ib = max( ib, 1 )
                  ib = min( ib, ibin )


                     ibnum(ib) = ibnum(ib) + 1


               end if


            end if


*-----------------------------------------------------------------------


                  nclst = jj


            do i = 1, nclst


               read(io,'(9i4,1p6e14.6)', end = 800 )
     &            ik, iz, in, ic, jj, id, is, iq, im,
     &            px, py, pz, et, rm, ex


               if( ik .eq. 0 .and.
     &             iz .eq. 0 .and.
     &             in .eq. 0 .and.
     &             ic .eq. 0 ) goto 800


            end do


                  goto 900


  800       continue


               write(ieo,'(/'' **** Reading data is terminated,''/
     &                      '' number of clusters is wrong'',
     &                      '' at event = '',
     &                   i6/'' filename = '',80a1)') ll,
     &                        (fname(jfile)(k:k), k=1,ifnl(jfile))


               stop


  900       continue


         end if


*-----------------------------------------------------------------------


         goto 5000


*-----------------------------------------------------------------------


 1000    continue


               icevnt = ll


               if( idevnt .lt. iprun ) iprun = idevnt


            if( ibch .eq. 2 ) then


               do i = 1, ibin


                  if( ibnum(i) .gt. 0 ) then


                     bweight(i) = bweight(i) / float(ibnum(i))
     &                          * iprun


                  else


                     bweight(i) = 0.0


                  end if


               end do


            end if




*-----------------------------------------------------------------------


               call cputime(5)


*-----------------------------------------------------------------------


      return
      end




