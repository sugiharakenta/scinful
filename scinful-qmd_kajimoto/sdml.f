************************************************************************
*                                                                      *
*        PART 6: Statistical Decay Model                               *
*                                                                      *
*                                                                      *
*   List of subprograms in rough order of relevance with main purpose  *
*      ( s = subroutine, f = function, b = block data, e = entry )     *
*                                                                      *
*                                                                      *
*  s  sdmtest   a short main for test of SDM                           *
*  s  sdment    to entry of SDM from QMD results                       *
*  s  sdmint    to initialize SDM                                      *
*  s  sdmexec   to execute statistical particle decay / fission.       *
*  s  sdmwid0   to calculate decay width and decay products            *
*               without angular momentum.                              *
*  s  sdmwid1   to calculate decay width and decay products            *
*               with angular momentum.                                 *
*  f  sdmlev    to calculate level density                             *
*  s  sdjsum    to sum up angular momentum.                            *
*  s  sdmfisw   to calculate fission width                             *
*  s  sdmfiss   to determine the mass and charge of fission fragment   *
*  f  sdmfisb   to determine the fission barrier.                      *
*  s  sdmlist   to list event record and particle data                 *
*  f  bndeng    to give binding energy per baryon (MeV)                *
*  f  eliq      to calculate liquid drop binding energy (MeV)          *
*  s  sdmtable  to read mass table from block data                     *
*                                                                      *
************************************************************************
*                                                                      *
************************************************************************
*                                                                      *
      subroutine sdmtest
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              A short main for test.                                  *
*                                                                      *
*                                                                      *
************************************************************************


      include 'param01.inc'


*-----------------------------------------------------------------------


      common /input1/ mstq1(mxpa1), parq1(mxpa1)
      common /rannum/ iseed, iseed0, iseed1


*-----------------------------------------------------------------------


            iseed  = 1234567


            mstq1(120) = 1
            mstq1(121) = 0
            mstq1(122) = 1
            mstq1(123) = 1
            mstq1(124) = 1


            parq1(120) = 1.0


            mstq1(125) = 30
            mstq1(126) = 1
            mstq1(127) = 1
            mstq1(128) = 1
            mstq1(129) = 1




            read(*,*) iseed


*-----------------------------------------------------------------------


         call sdmint
         call sdmexec(82,126,10,1.5*208., 0.,0.,0.,ierr)     ! 208Pb
         call sdmlist(6,1)


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine sdment
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to entry of SDM from QMD results                        *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'


*-----------------------------------------------------------------------


      common /sdmsw0/ issdm
      common /vriab3/ qmdfac, sdmfac


      common /clusti/ itc(0:nnn,0:nnn)
      common /clustf/ nclst, iclust(nnn)
      common /clustg/ jclust(0:7,nnn),  qclust(0:7,nnn)


*-----------------------------------------------------------------------


      common /clusts/ nclust, kclust(3,nnn)
      common /clustu/ lclust(0:7,nnn), sclust(0:7,nnn)


*-----------------------------------------------------------------------


      common /clustt/ nclsts, iclusts(nnn)
      common /clustw/ jclusts(0:7,nnn),  qclusts(0:7,nnn)


      common /clustv/ kdecay(4)
      common /clustx/ tdecay(4)


*-----------------------------------------------------------------------


         if( issdm .eq. 0 ) return


            call cputime(6)


*-----------------------------------------------------------------------
*        do loop for the clusters of QMD
*-----------------------------------------------------------------------


                  nclsts = 0


         do i = 1, nclst


*-----------------------------------------------------------------------
*        SDM ; statistical decay of cluster
*-----------------------------------------------------------------------


         if( iclust(i) .eq. 0 ) then


*-----------------------------------------------------------------------


                  jj = jclust(0,i)
                  iz = jclust(1,i)
                  in = jclust(2,i)
                  id = jclust(3,i)
                  is = jclust(4,i)
                  ic = jclust(5,i)


                  bi = qclust(0,i)
                  px = qclust(1,i)
                  py = qclust(2,i)
                  pz = qclust(3,i)


                  ex = qclust(6,i)


*-----------------------------------------------------------------------
*              for the case the cluster includes the Delta or N*
*-----------------------------------------------------------------------


               if( id .gt. 0 .or. is .gt. 0 ) then


                  jc = min( 0, ic )
                  in = max( 0, iz + in + id + is - jc )
                  iz = jc


               end if


*-----------------------------------------------------------------------
*              nucleus : statistical decay, call sdmexec
*-----------------------------------------------------------------------


                  ierr = 0


               call sdmexec(iz,in,jj,ex,px,py,pz,ierr)


                  if( ierr .ne. 0 ) goto 100


ccc            call sdmlist(6,0)


*-----------------------------------------------------------------------
*              total decay mode counter
*-----------------------------------------------------------------------


                  do l = 1, 4


                     tdecay(l) = tdecay(l)
     &                         + float(kdecay(l)) * sdmfac


                  end do


*-----------------------------------------------------------------------
*              new booking
*-----------------------------------------------------------------------


                  do j = 1, nclust


                     if( kclust(1,j) .ge. 100 ) then


                           nclsts = nclsts + 1


                        if( lclust(1,j) .eq. 0 .and.
     &                      lclust(2,j) .eq. 0 ) then


                           iclusts(nclsts) = 6


                        else if( lclust(1,j) .eq. 1 .and.
     &                           lclust(2,j) .eq. 0 ) then


                           iclusts(nclsts) = 1


                        else if( lclust(1,j) .eq. 0 .and.
     &                           lclust(2,j) .eq. 1 ) then


                           iclusts(nclsts) = 2


                        else


                           iclusts(nclsts) = 0


                        end if


                           jclusts(0,nclsts) = lclust(0,j)
                           jclusts(1,nclsts) = lclust(1,j)
                           jclusts(2,nclsts) = lclust(2,j)
                           jclusts(3,nclsts) = 0
                           jclusts(4,nclsts) = 0
                           jclusts(5,nclsts) = lclust(1,j)
                           jclusts(6,nclsts) = 0
                           jclusts(7,nclsts) = kclust(3,j)


                           qclusts(0,nclsts) = bi


                        do l = 1, 7


                           qclusts(l,nclsts) = sclust(l,j)


                        end do


                     end if


                  end do


*-----------------------------------------------------------------------


         end if


*-----------------------------------------------------------------------
*              particle : pass through
*-----------------------------------------------------------------------


  100       continue


            if( iclust(i) .gt. 0 .or.
     &          ierr .ne. 0 ) then


                          nclsts = nclsts + 1


                          iclusts(nclsts)   = iclust(i)


                       do j = 0, 7


                          jclusts(j,nclsts) = jclust(j,i)
                          qclusts(j,nclsts) = qclust(j,i)


                       end do


            end if


*-----------------------------------------------------------------------


         end do


*-----------------------------------------------------------------------


            call cputime(6)


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine sdmint
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to initialize SDM                                       *
*                                                                      *
*----------------------------------------------------------------------*
*     switch for SDM                                                   *
*----------------------------------------------------------------------*
*                                                                      *
*        issdm      : mstq1(120)                                       *
*                                                                      *
*----------------------------------------------------------------------*
*     minimum energy cut off energy for evapolation and fission        *
*----------------------------------------------------------------------*
*                                                                      *
*        sdmemin    : parq1(161)                                       *
*                     ( D = 1.0 MeV )                                  *
*                                                                      *
*----------------------------------------------------------------------*
*     switches for decay mode                                          *
*----------------------------------------------------------------------*
*                                                                      *
*        iswids     : mstq1(121)                                       *
*               = 0 : simple decay width without gamma nor angular mom.*
*                 1 : decay width with gammma and angular momentum.    *
*                                                                      *
*        isevap     : mstq1(122)                                       *
*               = 0 : without particle evapolation                     *
*                 1 : with particle evapolation                        *
*                                                                      *
*        isfiss     : mstq1(123)                                       *
*               = 0 : without fission decay                            *
*                 1 : with fission decay                               *
*                                                                      *
*        isgrnd     : mstq1(124)                                       *
*               = 0 : without ground state decay                       *
*                 1 : with ground state decay                          *
*                                                                      *
*----------------------------------------------------------------------*
*     parameters for sdmwid1                                           *
*----------------------------------------------------------------------*
*                                                                      *
*        imengb     : mstq1(125)                                       *
*                   : (D=30) number of energy bin                      *
*                                                                      *
*        imbarr     : mstq1(126)                                       *
*               = 0 : simple barrier                                   *
*                 1 : modified barrier of sepc                         *
*                                                                      *
*        imangm     : mstq1(127)                                       *
*               = 0 : without angular momentum                         *
*                 1 : with angular momentum                            *
*                                                                      *
*        imlevd     : mstq1(128)                                       *
*               = 0 : level density without angular momentum           *
*                 1 : level density with angular momentum              *
*                                                                      *
*        imgamm     : mstq1(129)                                       *
*               = 0 : without gamma decay                              *
*                 1 : with gamma decay                                 *
*                                                                      *
************************************************************************


      include 'param01.inc'


*-----------------------------------------------------------------------


      common /input1/ mstq1(mxpa1), parq1(mxpa1)


      common /vriab3/ qmdfac, sdmfac


      common /sdmcut/ sdmemin


      common /sdmsw0/ issdm
      common /sdmsw1/ iswids, isevap, isfiss, isgrnd
      common /sdmsw2/ imengb, imbarr, imangm, imlevd, imgamm


      common /clustx/ tdecay(4)


*-----------------------------------------------------------------------
*        switches for decay modes
*-----------------------------------------------------------------------


            issdm  = mstq1(120)


               sdmfac = 1.0


               if( issdm .gt. 0 ) sdmfac = 1.0 / float(issdm)


            iswids = mstq1(121)
            isevap = mstq1(122)
            isfiss = mstq1(123)
            isgrnd = mstq1(124)


*-----------------------------------------------------------------------
*        minimum energy cut off energy for evapolation and fission
*-----------------------------------------------------------------------


            sdmemin = parq1(120)


*-----------------------------------------------------------------------
*        parameters for sdmwid1
*-----------------------------------------------------------------------


            imengb = mstq1(125)
            imbarr = mstq1(126)
            imangm = mstq1(127)
            imlevd = mstq1(128)
            imgamm = mstq1(129)


            if( imlevd .eq. 0 ) imangm = 0




*-----------------------------------------------------------------------
*        initialize mass table
*-----------------------------------------------------------------------


            call sdmtable


*-----------------------------------------------------------------------
*        total decay mode counter
*-----------------------------------------------------------------------


            do i = 1, 4


               tdecay(i) = 0


            end do


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine sdmexec(iz,in,jang,ex,px,py,pz,ierr)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to execute statistical particle decay / fission.        *
*                                                                      *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              iz          : proton number                             *
*              in          : neutron number                            *
*              jang        : angular momentum                          *
*              px,py,pz    : momentum in GeV                           *
*              ierr        : error flag                                *
*                                                                      *
*                                                                      *
************************************************************************


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


      common /vriab3/ qmdfac, sdmfac
      common /sdmsw0/ issdm


      common /clusts/ nclust, kclust(3,nnn)
      common /clustu/ lclust(0:7,nnn), sclust(0:7,nnn)
      common /clustv/ kdecay(4)


      common /evapol/ ievmax, izip(0:100), inip(0:100), jip(0:100)


      common /summ03/ ireac(0:10), rcross(4)


*-----------------------------------------------------------------------


      common /masmx0/ imax0
      common /masid0/ inuc0(0:maxpt,0:maxnt)
      common /maspn0/ nzz0(0:nnuc), nnn0(0:nnuc), jgs0(0:nnuc)


      common /grnddc/ idgrnd(0:nnuc)


      common /sdmcut/ sdmemin


      common /sdmsw1/ iswids, isevap, isfiss, isgrnd
      common /sdmsw2/ imengb, imbarr, imangm, imlevd, imgamm


*-----------------------------------------------------------------------


      dimension pwid(0:100)
      dimension pcm(5)


*-----------------------------------------------------------------------
*        Error flag
*-----------------------------------------------------------------------


            ierr = 0


*-----------------------------------------------------------------------
*        Decay mode counter
*-----------------------------------------------------------------------


            do i = 1, 4


               kdecay(i) = 0


            end do


               jfiss = 0


*-----------------------------------------------------------------------
*        check of initial nucleus
*-----------------------------------------------------------------------


            if( iz .gt. maxpt .or. in .gt. maxnt ) then


               write(*,*) ' **** Error at sdmexec, too large iz, in'
               write(*,*) ' ======================================='
               write(*,*) ' iz    = ', iz,    ' in    = ', in
               write(*,*) ' maxpt = ', maxpt, ' maxnt = ', maxnt


               ierr = 1
               return


            end if


*-----------------------------------------------------------------------
*        Booking.  Store initial nucleus in K, L and S vectors.
*-----------------------------------------------------------------------


            nclust = 1


            kclust(1,1) = 0
            kclust(2,1) = 0
            kclust(3,1) = 0


            lclust(0,1) = jang


            lclust(1,1) = iz
            lclust(2,1) = in


            sclust(1,1) = px
            sclust(2,1) = py
            sclust(3,1) = pz


            be = bndeng(iz,in)


            pm = ( iz * rpmass + in * rnmass - be + ex ) / 1000.0


            sclust(4,1) = sqrt( pm**2 + px**2 + py**2 + pz**2 )
            sclust(5,1) = pm


            sclust(6,1) = ex


*-----------------------------------------------------------------------
*        minimum energy fot particle decay / fission
*-----------------------------------------------------------------------


               emin = sdmemin


*-----------------------------------------------------------------------
*     Start decay.
*-----------------------------------------------------------------------


               num    = 1
               inum   = 0


*-----------------------------------------------------------------------
 2000    continue
*-----------------------------------------------------------------------


               inum = inum + 1
               i1   = inum


*-----------------------------------------------------------------------
*           check dimension
*-----------------------------------------------------------------------


            if( num + 2 .gt. nnn ) then


               write(*,*) ' **** Error at sdmexec, too many products'
               write(*,*) ' ========================================='
               write(*,*) ' num + 2 = ', num + 2, ' nnn  = ', nnn


               ierr = 1
               return


            end if


*-----------------------------------------------------------------------
*        mother id
*-----------------------------------------------------------------------


               iz1 = lclust(1,i1)
               in1 = lclust(2,i1)
                m1 = iz1 + in1


               e1     = sclust(6,i1)


               em1 = iz1 * rpmass
     &             + in1 * rnmass
     &             - bndeng(iz1,in1)
     &             + e1


               j1     = lclust(0,i1)
               idata  = inuc0(iz1,in1)


               pcm(1) = sclust(1,i1)
               pcm(2) = sclust(2,i1)
               pcm(3) = sclust(3,i1)
               pcm(4) = sclust(4,i1)
               pcm(5) = sclust(5,i1)


*-----------------------------------------------------------------------
*           Gamma case
*-----------------------------------------------------------------------


            if( m1 .eq. 0 ) then


                  kclust(1,i1) = kclust(1,i1) + 100


                  goto 3000


            end if


*-----------------------------------------------------------------------
*        Reset fission and evaporation width.
*-----------------------------------------------------------------------


               tcs  = 0.0
               tcsf = 0.0


*-----------------------------------------------------------------------
*     Decay only if minimam excitation energy has.
*-----------------------------------------------------------------------


      if( e1 .gt. emin ) then


*-----------------------------------------------------------------------
*        particle evapolation width
*-----------------------------------------------------------------------


          if( isevap .eq. 1 ) then


               if( iswids .eq. 0 ) then


                  call sdmwid0(iz1,in1,e1,j1,
     &                         iz2,in2,e2,j2,
     &                         iz3,in3,e3,j3,
     &                         tcs,pwid,0,ipdec)


               else if( iswids .eq. 1 ) then


                  call sdmwid1(iz1,in1,e1,j1,
     &                         iz2,in2,e2,j2,
     &                         iz3,in3,e3,j3,
     &                         tcs,pwid,0,ipdec)


               end if


                  tcn = pwid(1)


          else


               do i = 0, ievmax


                  pwid(i) = 0.0


               end do


                  tcn = 1.0


          end if


*-----------------------------------------------------------------------
*        fission width,
*-----------------------------------------------------------------------


         if( isfiss .eq. 1 ) then


                  call sdmfisw(iz1,in1,e1,j1,tcn,tcsf,barf,saf)


          end if


*-----------------------------------------------------------------------


      end if


*-----------------------------------------------------------------------
*        Total decay width and branches
*-----------------------------------------------------------------------


               totwid = tcs + tcsf


               if( totwid .gt. 0.0 ) then


                  jdec = 1


               else if( idgrnd(idata) .ne. 0 .and.
     &                  isgrnd .eq. 1 ) then


                  jdec = 2


               else


                  jdec = 0


               end if


*-----------------------------------------------------------------------
*     Decay channel
*-----------------------------------------------------------------------




      if( jdec .eq. 1 .or. jdec .eq. 2) then


               ifiss = 0


*-----------------------------------------------------------------------
*        evaporation or fission
*-----------------------------------------------------------------------


         if( jdec .eq. 1 ) then


*-----------------------------------------------------------------------
*           decide evaporation or fission
*-----------------------------------------------------------------------


               xran  = rn() * totwid


*-----------------------------------------------------------------------
*           fission channnel
*-----------------------------------------------------------------------


            if( tcsf .ge. xran ) then


               call sdmfiss(iz1,in1,e1,j1,
     &                      iz2,in2,e2,j2,
     &                      iz3,in3,e3,j3,
     &                      barf,saf,ierr)


               if( ierr .ne. 0 ) return


               ifiss = 1
               jfiss = jfiss + 1


               fisfac = qmdfac * sdmfac


               call sm_mass(3,iz1,in1,fisfac)


            else


*-----------------------------------------------------------------------
*           Evaporation channels, find decay branch.
*-----------------------------------------------------------------------


                  tcip0 = tcsf


               do ipdec = 0, ievmax


                  tcip0 = tcip0 + pwid(ipdec)


                  if( tcip0 .ge. xran ) goto 101


               end do


*-----------------------------------------------------------------------
*                 Decay branch could not be found.
*                 Choice maxmum width mode
*-----------------------------------------------------------------------


                        ipdmx = 0
                        tcmx  = 0


                  do ipdec = 0, ievmax


                     if( tcmx .lt. pwid(ipdec) ) then


                        ipdmx = ipdec
                        tcmx  = pwid(ipdec)


                     end if


                  end do


                        ipdec = ipdmx


*-----------------------------------------------------------------------


  101          continue


*-----------------------------------------------------------------------
*              Determine final nuclei.
*-----------------------------------------------------------------------


               if( iswids .eq. 0 ) then


                  call sdmwid0(iz1,in1,e1,j1,
     &                         iz2,in2,e2,j2,
     &                         iz3,in3,e3,j3,
     &                         tcs,pwid,1,ipdec)


               else if( iswids .eq. 1 ) then


                  call sdmwid1(iz1,in1,e1,j1,
     &                         iz2,in2,e2,j2,
     &                         iz3,in3,e3,j3,
     &                         tcs,pwid,1,ipdec)




               end if


*-----------------------------------------------------------------------


            end if


*-----------------------------------------------------------------------
*           Mass of produced particle in MeV.
*-----------------------------------------------------------------------


            if( iz3 + in3 .gt. 0 ) then


               m2  = iz2 + in2
               m3  = iz3 + in3


               em2 = iz2 * rpmass
     &             + in2 * rnmass
     &             - bndeng(iz2,in2)
     &             + e2


               em3 = iz3 * rpmass
     &             + in3 * rnmass
     &             - bndeng(iz3,in3)
     &             + e3


            else


               m2  = m1
               m3  = 0


               em2 = iz2 * rpmass
     &             + in2 * rnmass
     &             - bndeng(iz2,in2)
     &             + e2


               em3 = 0.0


            end if


*-----------------------------------------------------------------------
*        ground state decay
*-----------------------------------------------------------------------


         else if( jdec .eq. 2 ) then


*-----------------------------------------------------------------------


               iz2 = nzz0( idgrnd(idata) )
               in2 = nnn0( idgrnd(idata) )
                m2 = iz2 + in2
                j2 = 0
                e2 = 0


               iz3 = iz1 - iz2
               in3 = in1 - in2
                m3 = iz3 + in3
                j3 = 0
                e3 = 0


               em2 = iz2 * rpmass
     &             + in2 * rnmass
     &             - bndeng(iz2,in2)
     &             + e2


               em3 = iz3 * rpmass
     &             + in3 * rnmass
     &             - bndeng(iz3,in3)
     &             + e3


*-----------------------------------------------------------------------


         end if


*-----------------------------------------------------------------------
*        line number for two daughters
*-----------------------------------------------------------------------


               i2 = num + 1
               i3 = num + 2


               num    = num + 2
               nclust = nclust + 2


*-----------------------------------------------------------------------
*        Determine momenta of decay products.
*        Lorentz-transformation into reference frame.
*-----------------------------------------------------------------------


               em1  = em1 / 1000.
               em2  = em2 / 1000.
               em3  = em3 / 1000.


               pr   = pcmsr( pcm(5), em2, em3 )


               cos1 = 1.0 - 2.0 * rn()
               sin1 = sqrt( 1.0 - cos1**2 )
               phi1 = 2.0 * pi * rn()


               pxr  = pr * sin1 * cos(phi1)
               pyr  = pr * sin1 * sin(phi1)
               pzr  = pr * cos1


               pcs    = pcm(1) * pxr + pcm(2) * pyr + pcm(3) * pzr


               ecm1   = sqrt( em2**2 + pxr**2 + pyr**2+ pzr**2 )
               trans1 = ( pcs / ( pcm(4) + pcm(5) ) + ecm1 ) / pcm(5)


               ecm2   = sqrt( em3**2 + pxr**2 + pyr**2 + pzr**2 )
               trans2 = ( -pcs / ( pcm(4) + pcm(5) ) + ecm2 ) / pcm(5)


*-----------------------------------------------------------------------
*           counting decay mode and booking
*-----------------------------------------------------------------------


            if( ifiss .eq. 0 ) then


               if( jdec .eq. 1 ) then


                  if( m3 .ne. 0 ) then


                     kdecay(1) = kdecay(1) + 1


                     kclust(1,i3) = kclust(1,i1) / 10 * 10 + 1
                     kclust(1,i2) = kclust(1,i1) / 10 * 10 + 2


                  else


                     kdecay(2) = kdecay(2) + 1


                     kclust(1,i3) = kclust(1,i1) / 10 * 10 + 3
                     kclust(1,i2) = kclust(1,i1) / 10 * 10 + 4


                  end if


               else


                     kdecay(3) = kdecay(3) + 1


                     kclust(1,i3) = kclust(1,i1) / 10 * 10 + 5
                     kclust(1,i2) = kclust(1,i1) / 10 * 10 + 6


               end if




            else


                     kdecay(4) = kdecay(4) + 1


                     kclust(1,i2) = kclust(1,i1) / 10 * 10 + 10
                     kclust(1,i3) = kclust(1,i1) / 10 * 10 + 10


                     if( kclust(1,i2) .ge. 100 ) kclust(1,i2) = 90
                     if( kclust(1,i3) .ge. 100 ) kclust(1,i3) = 90


            end if


*-----------------------------------------------------------------------
*           determine the quantities of two daughters
*-----------------------------------------------------------------------


               lclust(0,i2) = j2
               lclust(1,i2) = iz2
               lclust(2,i2) = in2


               kclust(2,i2) = i1
               kclust(3,i2) = kclust(3,i1) + 1


               sclust(1,i2) = pxr + pcm(1) * trans1
               sclust(2,i2) = pyr + pcm(2) * trans1
               sclust(3,i2) = pzr + pcm(3) * trans1
               sclust(4,i2) = sqrt( em2**2
     &                      + sclust(1,i2)**2
     &                      + sclust(2,i2)**2
     &                      + sclust(3,i2)**2 )
               sclust(5,i2) = em2
               sclust(6,i2) = e2


*-----------------------------------------------------------------------


               lclust(0,i3) = j3
               lclust(1,i3) = iz3
               lclust(2,i3) = in3


               kclust(2,i3) = i1
               kclust(3,i3) = kclust(3,i1) + 1


               sclust(1,i3) = -pxr + pcm(1) * trans2
               sclust(2,i3) = -pyr + pcm(2) * trans2
               sclust(3,i3) = -pzr + pcm(3) * trans2
               sclust(4,i3) = sqrt( em3**2
     &                      + sclust(1,i3)**2
     &                      + sclust(2,i3)**2
     &                      + sclust(3,i3)**2 )
               sclust(5,i3) = em3
               sclust(6,i3) = e3


*-----------------------------------------------------------------------
*     final state
*-----------------------------------------------------------------------


      else if( jdec .eq. 0 ) then




                  kclust(1,i1) = kclust(1,i1) + 100




      end if


*-----------------------------------------------------------------------


 3000    if( num .gt. inum ) goto 2000


*-----------------------------------------------------------------------
*        fission cross section
*-----------------------------------------------------------------------


            if( jfiss .gt. 0 ) then


               rcross(4) = rcross(4) + qmdfac * sdmfac


            end if


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine sdmwid0(iz1,in1,ex1,j1,
     &                   iz2,in2,ex2,j2,
     &                   iz3,in3,ex3,j3,
     &                   totwid,pwid,ifdec,ipdec)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate decay width and decay products             *
*              without angular momentum.                               *
*                                                                      *
*              statistical decay of light particles                    *
*              according to the temperature.                           *
*              neutron proton deutron triton 3he 4he                   *
*              are considered                                          *
*                                                                      *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              izi         : charge of fragment                        *
*              ini         : nuetron number of fragment                *
*              exi         : excitation energy                         *
*              ji          : angular momentum                          *
*                                                                      *
*              totwid      : total decay width                         *
*              pwid        : partial decay width                       *
*                                                                      *
*              ifdec       : switch for this subroutine                *
*                      = 0 : calcuate decay width                      *
*                      = 1 : calculate decay product                   *
*                                                                      *
*              ipdec       : decay channel for the decay products      *
*                                                                      *
*                                                                      *
************************************************************************


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


      parameter      ( ddct = - 30.0 )
      parameter      ( mcrit = 1000 )


      parameter      ( denpa = 8.0 )


*-----------------------------------------------------------------------


      common /evapol/ ievmax, izip(0:100), inip(0:100), jip(0:100)
      common /evaval/ sepc(2,0:100)


      dimension   pwid(0:100)
      dimension   qvalp(0:100), eclbp(0:100)
      dimension   smi(0:100), ey1(0:100), s(0:100)


      save        qvalp, eclbp, smi, ey1, s


*-----------------------------------------------------------------------


               e1 = ex1
               m1 = iz1 + in1


*-----------------------------------------------------------------------


         if( ifdec .eq. 1 ) goto 3000


*-----------------------------------------------------------------------


                pwid(0) = 0.0


            do ipid = 1, ievmax


               eclbp(ipid) = 0.0
               qvalp(ipid) = 0.0
                pwid(ipid) = 0.0


            end do


*-----------------------------------------------------------------------
*     Calculate total and partial decay width.
*-----------------------------------------------------------------------


               s0 = 2.0 * sqrt( float(m1) / denpa * e1 )


               totwid = 0.0


*-----------------------------------------------------------------------


      do 100 ipid = 1, ievmax


*-----------------------------------------------------------------------


               s(ipid) = 0.0


*-----------------------------------------------------------------------
*         decay of 1 => 2 + 3
*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*           Set emission particle
*-----------------------------------------------------------------------


               iz3 = izip(ipid)
               in3 = inip(ipid)
                m3 = iz3 + in3


*-----------------------------------------------------------------------
*           Set residual particle
*-----------------------------------------------------------------------


               iz2 = iz1 - iz3
               in2 = in1 - in3
                m2 = iz2 + in2


*-----------------------------------------------------------------------
*           check
*-----------------------------------------------------------------------


               if( iz3 .gt. iz2 .or.
     &              m3 .gt.  m2 .or.
     &             iz2 .ge.  m2 .or.
     &             iz2 .lt.   0 .or.
     &             in2 .lt.   0       ) goto 100


*-----------------------------------------------------------------------
*           Q-value and Coulomb barrier
*-----------------------------------------------------------------------


               qvalp(ipid) = bndeng(iz1,in1)
     &                     - bndeng(iz2,in2)
     &                     - bndeng(iz3,in3)


               a23 = float( m2 )**(1./3.)


               eclbp(ipid) =  ccoul * 1000.0
     &                     * float( iz2 ) * float( iz3 )
     &                     / ( sepc(1,ipid) + sepc(2,ipid) * a23 )


*-----------------------------------------------------------------------
*           Level density
*-----------------------------------------------------------------------


               ss = e1 - eclbp(ipid) - qvalp(ipid)


                  if( ss .le. 0.0 ) goto 100


               sa = float(m2) / denpa


               s(ipid)  = 2.0 * sqrt( sa * ss )


            if( s(ipid) - s0 .lt. ddct ) then


                spl = 0.0


            else


                spl = exp( s(ipid) - s0 )


            end if


            if( s(ipid) .gt. -ddct ) then


                smi(ipid) = 0.0


            else


                smi(ipid) = exp( - s(ipid) )


            end if




               ey1(ipid)  = ( ( 2.0 * s(ipid)**2 - 6.0 * s(ipid) + 6.0 )
     &                      + ( s(ipid)**2 - 6.0 ) * smi(ipid) )


*-----------------------------------------------------------------------
*           Partial decay width and total decay width
*-----------------------------------------------------------------------


               pwid(ipid) = float(m3) * a23**2
     &                    / sa**2 * spl * ey1(ipid)


               if( pwid(ipid) .lt. 0.0 ) pwid(ipid) = 0.0




               totwid = totwid + pwid(ipid)


*-----------------------------------------------------------------------


  100 continue


*-----------------------------------------------------------------------


         return






*-----------------------------------------------------------------------
*     Set decay particles:  ifdec = 1
*-----------------------------------------------------------------------


 3000 continue


*-----------------------------------------------------------------------


               iz3 = izip(ipdec)
               in3 = inip(ipdec)
               ex3 = 0.0


               iz2 = iz1 - iz3
               in2 = in1 - in3


               m2 = iz2 + in2
               m3 = iz3 + in3


               j3 = jip(ipdec)
               j2 = max( 0, j1 - j3 )


*-----------------------------------------------------------------------
*        Determine the energy and momentum.
*-----------------------------------------------------------------------


            if( e1 - qvalp(ipdec) - eclbp(ipdec) .lt. 0.1 .or.
     &          m2 .eq. 1 ) then


               erel = max( 0.0, e1 - qvalp(ipdec) )


            else


               sa = float(m2) / denpa


               t  = ( ( s(ipdec)**3 - 6.0 * s(ipdec)**2
     &                + 15.0 * s(ipdec) -15.0 )
     &              + smi(ipdec) / 8.0
     &              * ( s(ipdec)**4 - 12.0 * s(ipdec)**2
     &                + 15.0 * s(ipdec) ) )
     &              / sa / ey1(ipdec)


               if( t .le. 0.0 ) then


                  erel = max( 0.0, e1 - qvalp(ipdec) )
                  goto 330


               end if


               ermax = e1 - qvalp(ipdec)
               smax  = t / 2.718282




               itry = 0
   20          itry = itry + 1


               if( itry .ge. mcrit ) goto 330


               erel = eclbp(ipdec) + ( ermax - eclbp(ipdec) ) * rn()
               erc  = erel - eclbp(ipdec)


               if( erc * exp( -erc / t ) .lt. smax * rn() ) goto 20


            end if


  330          continue


               ex2  = max( 0.0, e1 - erel - qvalp(ipdec) )




*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine sdmwid1(iz1,in1,e14,j1,
     &                   iz2,in2,e24,j2,
     &                   iz3,in3,e34,j3,
     &                   tcs04,tcip04,ifdec,ipdec)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate decay width and decay products             *
*              with angular momentum.                                  *
*                                                                      *
*              statistical decay of light particles                    *
*              according to the temperature.                           *
*              neutron proton deutron triton 3he 4he                   *
*              are considered                                          *
*                                                                      *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              izi         : charge of fragment                        *
*              ini         : nuetron number of fragment                *
*              exi         : excitation energy                         *
*              ji          : angular momentum                          *
*                                                                      *
*              tcs04       : total decay width                         *
*              tcip04      : partial decay width                       *
*                                                                      *
*              ifdec       : switch for this subroutine                *
*                      = 0 : calcuate decay width                      *
*                      = 1 : calculate decay product                   *
*                                                                      *
*              ipdec       : decay channel for the decay products      *
*                                                                      *
*                                                                      *
************************************************************************


      include 'param01.inc'


*-----------------------------------------------------------------------


      implicit real*8(a-h,o-z)


*-----------------------------------------------------------------------


      real*4  e14, e24, e34, tcip04(0:100), tcs04, rn, bndeng
      real*4  sepc, sdmemin


*-----------------------------------------------------------------------


      common /evapol/ ievmax, izip(0:100), inip(0:100), jip(0:100)
      common /evaval/ sepc(2,0:100)


      common /sdmcut/ sdmemin


      common /sdmsw2/ imengb, imbarr, imangm, imlevd, imgamm


*-----------------------------------------------------------------------
*     Some constants
*-----------------------------------------------------------------------


      parameter ( rmass = 938.3d0, hc  = 197.3d0 )
      parameter ( ccoul = 1.439767d0 )
      parameter ( pi = 3.141592653589793d0 )
      parameter ( xi1 = 0.5d-8, xi2 = 0.3d-9 )


*-----------------------------------------------------------------------


      parameter( maxtc = 5200 )


      dimension tcsave(maxtc), tcip(0:100), ntcoffset(0:100)
      dimension j2idx(maxtc),ieidx(maxtc),lmaxsv(maxtc)


      dimension de(1000), ee(1000), we(1000)
      dimension itmp1(1000),itmp2(1000),ltmp(1000),iesv(1000)


*-----------------------------------------------------------------------


      save     tcsave, ntcoffset
      save     j2idx,ieidx,lmaxsv,ng


*-----------------------------------------------------------------------
*           iemax : number of energy bin
*           ifj   : switch of angular momentum
*           emin  : minimum energy of decay particles
*-----------------------------------------------------------------------


            iemax = imengb
            ifj   = imangm
            emin  = dble(sdmemin)


*-----------------------------------------------------------------------


            e1    = dble(e14)
            j1    = j1 * ifj


*-----------------------------------------------------------------------


            if( ifdec .eq. 1 ) then


                  tcs    = tcip(ipdec)
                  ifsave = 0
                  ifcalc = 0


            else


                  ifsave = 1
                  ifcalc = 1


               do ip = 0, ievmax


                  tcip(ip) = 0


               end do


            end if


*-----------------------------------------------------------------------


               ntc    = 0
               random = dble(rn())


*-----------------------------------------------------------------------


               m1     = iz1 + in1
               r1     = 1.2 * float( m1 )**(1./3.)
               de1    = emin
               dl1    = sdmlev(iz1,in1,e1,de1,j1)
               denomi = dl1 * 2 * pi * de1


*-----------------------------------------------------------------------
*        dl1 = 0.0 => parent is in the wrong state. but continue.
*-----------------------------------------------------------------------


            if( denomi .lt. 0.0001 ) then


               goto 888


            end if


*-----------------------------------------------------------------------
*        gamma deacy
*-----------------------------------------------------------------------


               ip    = 0
               tcip0 = 0


               ntcoffset(ip) = ntc


         if( imgamm .eq. 1 .and.
     &     ( ifdec .eq. 0 .or. ip .eq. ipdec ) ) then


               iz2 = iz1
               in2 = in1
                m2 = m1


               iz3 = 0
               in3 = 0
                m3 = 0
                e3 = 0
                j3 = 0


*-----------------------------------------------------------------------
*           make energy bin
*-----------------------------------------------------------------------


                  emax = e1
                  pow  = 1.0 / ( 0.4 * sqrt( emax * m1 / 8.0 ) )


               do ie2 = 1, iemax


                  yy = ( ie2 - 0.5 ) / ( iemax - 1 )
                  if( ie2 .eq. iemax ) yy = 1.0


                  we(ie2) = emax * yy**pow


               end do


                  ee(1) = we(1) / 2.0
                  de(1) = we(1)


               do ie2 = 2, iemax


                  ee(ie2) = ( we(ie2) + we(ie2-1) ) / 2.0
                  de(ie2) =   we(ie2) - we(ie2-1)


               end do


*-----------------------------------------------------------------------
*           gamma width
*-----------------------------------------------------------------------


                  itmp  = ifj * ( j1 + 2 ) - ifj * max( 0, j1 - 2 ) + 1
                  max10 = iemax * itmp


*-----------------------------------------------------------------------
            if( ifdec .eq. 0 ) then
*-----------------------------------------------------------------------


               do 10 i = 1, max10


                  ie2 = ( i - 1 ) / itmp + 1
                  j2  = ifj * max( 0, j1 - 2 ) + mod( i - 1 , itmp)


                  ieidx(i) = ie2
                  j2idx(i) = j2


                  de2  = de(ie2)
                  e2   = ee(ie2)
                  erel = e1 - e2


                  dl = sdmlev(iz1,in1,e2,de2,j2)


                  if( ( abs( j2 - j1 ) .le. 2 ) .and.
     &                ( j2 + j1 + mod(m1,2) .ge. 2 ) ) then


                        tcoeff = xi2 * erel**5 * dl * de2 / denomi


                     if( ( abs( j2 - j1 ) .le. 1 ) .and.
     &                   ( j2 + j1 + mod(m1,2) .ge. 1 ) ) then


                        tcoeff = xi2 * erel**5 * dl * de2 / denomi
     &                         + xi1 * erel**3 * dl * de2 / denomi


                     end if


                  else


                        tcoeff = 0.0


                  end if


                        tcsave(i) = tcoeff
                        tcip0     = tcip0 + tcoeff


   10          continue


                        ntc = max10


*-----------------------------------------------------------------------
            else
*-----------------------------------------------------------------------


              do 101 i = 1, max10


                    ie2   = ieidx(i)
                    j2    = j2idx(i)
                    de2   = de(ie2)
                    e2    = ee(ie2)
                    erel  = e1 - e2
                    tcip0 = tcip0 + tcsave(i)


                    if( tcip0 .ge. tcs * random ) goto 999


  101         continue


            end if


*-----------------------------------------------------------------------
         end if
*-----------------------------------------------------------------------


                     tcip(ip) = tcip0




*-----------------------------------------------------------------------
*     particle decay
*-----------------------------------------------------------------------


      do 30 ip = 1, ievmax


               if( ifdec .eq. 0 ) then


                  ntcoffset(ip) = ntc


               else if( ifdec .eq. 1 ) then


                  ntc = ntcoffset(ip)


               end if


*-----------------------------------------------------------------------


         if( ifdec .eq. 0 .or. ip .eq. ipdec ) then


*-----------------------------------------------------------------------


                  iz3 = izip(ip)
                  in3 = inip(ip)
                   j3 = jip(ip) * ifj


                   m3 = iz3 + in3
                   e3 = 0


                  iz2 = iz1 - iz3
                  in2 = in1 - in3
                   m2 = iz2 + in2


                  tcip0 = 0


*-----------------------------------------------------------------------


            if(  m2 .gt. 0 .and.
     &          iz2 .ge. 0 .and.
     &          in2 .ge. 0 ) then


                  r2 = 1.2 * float( m2 )**(1./3.)
                  r3 = 1.2 * float( m3 )**(1./3.)


                  qval12 = dble( bndeng(iz1,in1)
     &                         - bndeng(iz2,in2)
     &                         - bndeng(iz3,in3) )


               if( imbarr .eq. 1 ) then


                  a23 = float( m2 )**(1./3.)
                  b23 = iz2 * iz3 * ccoul
     &                / ( dble(sepc(1,ip)) + dble(sepc(2,ip)) * a23 )


               else


                  b23 = iz2 * iz3 * ccoul / ( r2 + r3 )


               end if


                  e2max = e1 - qval12 - b23


*-----------------------------------------------------------------------


            if( e2max .gt. 0 ) then


*-----------------------------------------------------------------------
*              make energy bin
*-----------------------------------------------------------------------


                  emax = e2max
                  pow  = 1.0 / ( 0.4 * sqrt( emax * m1 / 8.0 ) )


               do ie2 = 1, iemax


                  yy = ( ie2 - 0.5 ) / ( iemax - 1 )
                  if( ie2 .eq. iemax ) yy = 1.0


                  we(ie2) = emax * yy**pow


               end do


                  ee(1) = we(1) / 2.0
                  de(1) = we(1)


               do ie2 = 2, iemax


                  ee(ie2) = ( we(ie2) + we(ie2-1) ) / 2.0
                  de(ie2) =   we(ie2) - we(ie2-1)


               end do




*-----------------------------------------------------------------------
            if( ifdec .eq. 0 ) then
*-----------------------------------------------------------------------


                           ntcold = ntc
                           ntctmp = 0


                  do 40 ie2 = 1, iemax


                           de2  = de(ie2)
                           e2   = ee(ie2)
                           erel = e1 - e2 - qval12


                           angmx2 = 2.0 * rmass * m3 * m2 / m1
     &                            * ( r2 + r3 )**2
     &                            * ( erel - b23 )


                        if( angmx2 .ge. 0 ) then


                           iesv(ie2)  = 1
                           ltmp(ie2)  = int( sqrt(angmx2) / hc )
                           itmp1(ie2) = ifj
     &                                * max( 0, j1 - ltmp(ie2) - 1 )
                           itmp2(ie2) = ifj * ( J1 + ltmp(ie2) + 1 )


                           ntctmp = ntctmp + itmp2(ie2) - itmp1(ie2) + 1


                        else


                           iesv(ie2) = 0


                        end if


  40              continue


                     if( ntctmp + ntc .gt. maxtc ) then


                           ng = ip


                           write(*,*)
     &                     ' **** Warning at sdmwid1, over maxtc'
                           write(*,*)
     &                     ' ==================================='
                           write(*,*) ' over maxtc : ntc = ', ntc


                     else


                           ng = -1


                     end if




               if( ng .lt. 0 ) then


                  do 41 ie2 = 1, iemax


                     if( iesv(ie2) .eq. 1 ) then


                        do 50 j2 = itmp1(ie2), itmp2(ie2)


                           ntc = ntc + 1


                           j2idx(ntc)  = j2
                           ieidx(ntc)  = ie2
                           lmaxsv(ntc) = ltmp(ie2)


  50                    continue


                     end if


  41              continue


                  do 444 i = ntcold + 1, ntc


                           ie2  = ieidx(i)
                           j2   = j2idx(i)
                           de2  = de(ie2)
                           e2   = ee(ie2)
                           lmax = lmaxsv(i)


                        if( ifj .eq. 1 ) then


                           call sdjsum( 2 * J1 + mod(m1,2),
     &                                  2 * J2 + mod(m2,2),
     &                                  2 * J3 + mod(m3,2),
     &                                  lmax, w )


                        else


                           w = 1


                        end if


                           dl = sdmlev(iz2,in2,e2,de2,j2)


                           tcoeff    = w * dl * de2 / denomi
                           tcip0     = tcip0 + tcoeff
                           tcsave(i) = tcoeff


 444              continue




               else




                  do 411 ie2 = 1, iemax


                     if( iesv(ie2) .eq. 1 ) then


                        do 501 j2 = itmp1(ie2), itmp2(ie2)


                           ntc  = ntc + 1
                           lmax = ltmp(ie2)
                           de2  = de(ie2)
                           e2   = ee(ie2)


                           if( ifj .eq. 1 ) then


                              call sdjsum( 2 * j1 + mod(m1,2),
     &                                     2 * J2 + mod(m2,2),
     &                                     2 * J3 + mod(m3,2),
     &                                     lmax, w )


                           else


                              W = 1


                           end if


                           dl = sdmlev(iz2,in2,e2,de2,j2)


                           tcoeff = w * dl * de2 / denomi
                           tcip0  = tcip0 + tcoeff


 501                    continue


                     end if


 411              continue


               end if


                           ntcoffset(ip+1) = ntc


*-----------------------------------------------------------------------
            else
*-----------------------------------------------------------------------


               if( ng .lt. ip ) then


*vocl loop,scalar


                  do 441 i = ntcoffset(ip) + 1, ntcoffset(ip+1)


                           tcip0 = tcip0 + tcsave(i)


                     if( tcip0 .ge. tcs * random ) then


                           ie2  = ieidx(i)
                           de2  = de(ie2)
                           e2   = ee(ie2)
                           erel = e1 - e2 - qval12
                           lmax = lmaxsv(i)
                           j2   = j2idx(i)


                           goto 999


                     end if


  441             continue


               else


*vocl loop,scalar


                  do 44 ie2 = 1, iemax


                           de2  = de(ie2)
                           e2   = ee(ie2)
                           erel = e1 - e2 - qval12


                           angmx2 = 2.0 * rmass * m3 * m2 / m1
     &                            * ( r2 + r3 )**2
     &                            * ( erel - b23 )


                     if( angmx2 .ge. 0 ) then


                           lmax = int( sqrt(angmx2) / hc )
*vocl loop,scalar


                        do 55 J2 = ifj * max( 0, j1 - lmax - 1 ),
     &                             ifj * ( j1 + lmax + 1 )


                           ntc = ntc + 1


                           if( ifj .eq. 1 ) then


                              call sdjsum( 2 * j1 + mod(m1,2),
     &                                     2 * j2 + mod(m2,2),
     &                                     2 * j3 + mod(m3,2),
     &                                     lmax, w )


                           else


                              w = 1


                           end if


                           dl = sdmlev(iz2,in2,e2,de2,j2)


                           tcoeff = w * dl * de2 / denomi
                           tcip0  = tcip0 + tcoeff


                           if( tcip0 .ge. tcs * random ) goto 999


  55                    continue


                     end if


  44              continue


               end if


            end if




*-----------------------------------------------------------------------


            end if


            end if


                  tcip(ip) = tcip0


*-----------------------------------------------------------------------


         end if


*-----------------------------------------------------------------------
   30 continue
*-----------------------------------------------------------------------




            if( ifdec .eq. 1 ) then


               write(*,*) ' **** Warning at sdmwid1, no decay channel'
               write(*,*) ' ========================================='
               write(*,*) '  iz1,in1,e1,j1,tcs = ', iz1,in1,e1,j1,tcs
               write(*,*) '  ipdec,tcs,tcip(ipdec) = ',
     &                       ipdec,tcs,tcip(ipdec)


            end if


*-----------------------------------------------------------------------


  888    continue


                  tcs   = 0.0
                  tcs04 = 0.0


               do 100 ip = 0, ievmax


                  tcs = tcs + tcip(ip)
                  tcip04(ip) = real(tcip(ip))
                  tcs04 = tcs04 + tcip04(ip)


  100          continue


               return


*-----------------------------------------------------------------------


  999 continue


               e24   = real(e2)
               e34   = real(e3)




*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      function sdmlev(iz,in,e,de,j)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate level density according to                 *
*              Ref. F.Puhlhofer, Nucl. Phys. A280(1977)267.            *
*                                                                      *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              iz          : proton number                             *
*              in          : neutron number                            *
*              e           : excitation energy (MeV)                   *
*              de          : energy bin (MeV)                          *
*              j           : angluar momentum                          *
*                                                                      *
*                                                                      *
************************************************************************


      include 'param01.inc'


*-----------------------------------------------------------------------


      implicit real*8(a-h,o-z)


      real*4 sdmemin


*-----------------------------------------------------------------------


      common /sdmsw2/ imengb, imbarr, imangm, imlevd, imgamm


      common /sdmcut/ sdmemin


*-----------------------------------------------------------------------


      parameter( rmass = 938.3d0, hc  = 197.3d0 )
      parameter( denpa = 8.0d0 )


*-----------------------------------------------------------------------
*              initial value
*-----------------------------------------------------------------------


               denl = 0


               if( e .lt. dble(sdmemin) .and.
     &             j .eq. 0 ) denl = 1.d0 / de


               m    = iz + in
               amas = m * 1.0d0
               a    = amas / denpa


*-----------------------------------------------------------------------
*     a =< 4
*-----------------------------------------------------------------------


            if( m .le. 4 ) then


                  sdmlev = denl


                  return


            end if


*-----------------------------------------------------------------------
*     a > 4
*-----------------------------------------------------------------------
*        simple level density
*-----------------------------------------------------------------------


         if( imlevd .eq. 0 ) then


               denl = exp( 2 * sqrt( a * e ) ) / 1000


*-----------------------------------------------------------------------
*        Level density from F. Puhlhofer.
*-----------------------------------------------------------------------


         else


            if( ( m / 2 ) * 2 .ne. m ) then


               delta = -0.7d0


            else if( ( iz / 2 ) * 2. ne. iz ) then


               delta = -2.0d0


            else


               delta = 0.7d0


            end if


               r = 4 * ( amas**(5.d0/3.d0) )
     &           * rmass * ( 1.2d0**2 ) / 5.0d0 / hc / hc / a


               u1 = e - ( j + mod(m,2) / 2.0d0 )**2 / a / r - delta
               u2 = e - ( j + mod(m,2) / 2.0d0 + 1 )**2 / a / r - delta


            if( u2 .gt. 0 ) then


               t1 = ( 1.5d0 + sqrt( 1.5d0 * 1.5d0 + 4 * a * u1 ) )
     &            / 2.0d0 / a


               t2 = ( 1.5d0 + sqrt( 1.5d0 * 1.5d0 + 4 * a * u2 ) )
     &            / 2.0d0 / a


               denl = 1.0d0 / 12.0d0 / sqrt(r) / a / a
     &              * ( exp( 2.0d0 * sqrt( a * u1 ) )
     &              / t1 / t1 / t1
     &              - exp( 2.0d0 * sqrt( a * u2 ) )
     &              / t2 / t2 / t2 )


            end if


         end if




               sdmlev = denl


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine sdjsum(j1,j2,j3,lmax,w)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to sum up angular momentum.                             *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              j1, j2, j3  : angular momentum                          *
*              lmax        : maximum relative angular momentum         *
*              w           : number of states                          *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      implicit real*8(a-h,o-z)


*-----------------------------------------------------------------------


            w = 0.0


            if( mod( j1 + j2 + j3, 2 ) .eq. 1 ) return


            jj1 = abs( j2 - j3 )
            jj2 = j2 + j3


*-----------------------------------------------------------------------


            if( jj2 - jj1 .eq. 0 ) then


               l1 = abs( j1 - jj1 ) / 2
               l2 = min( ( j1 + jj1 ) / 2, lmax )


               w  = dim( l2 + 1 - l1, 0 )


            else if( jj2 - jj1 .eq. 2 ) then


               l1 = abs( j1 - jj1 ) / 2
               l2 = min( ( j1 + jj1 ) / 2, lmax )
               l3 = abs( j1 - ( jj1 + 2 ) ) / 2
               l4 = min( ( j1 + jj1 + 2 ) / 2, lmax )


               w  = dim( l2 + 1 - l1, 0 ) + dim( l4 + 1 - l3, 0 )


            else


               l1 = abs( j1 - jj1 ) / 2
               l2 = min( ( j1 + jj1 ) / 2, lmax )
               l3 = abs( j1 - ( jj1 + 2 ) ) / 2
               l4 = min( ( j1 + jj1 + 2 ) / 2, lmax )
               l5 = abs( j1 - ( jj1 + 4 ) ) / 2
               l6 = min(( j1 + jj1 + 4 ) / 2, lmax )


               w  = dim( l2 + 1 - l1, 0 )
     &            + dim( l4 + 1 - l3, 0 )
     &            + dim( l6 + 1 - l5, 0 )


            end if


*-----------------------------------------------------------------------


      return
      end


************************************************************************
*                                                                      *
      subroutine sdmfisw(iz1,in1,e1,j1,tcn,tcsf,barf,saf)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate fission width                              *
*                                                                      *
*                                                                      *
************************************************************************


      include 'param01.inc'


*-----------------------------------------------------------------------


      parameter      ( ddct = - 30.0 )


      parameter      ( sanf  =  10.0,
     &                 saa   =  7.9169e-08,
     &                 sbb   = -1.6143e-04,
     &                 scc   =  1.0922,
     &                 rk0   =  13.0 )


*-----------------------------------------------------------------------


               tcsf  = 0.0


               m1 = iz1 + in1


               if( m1 .le. 20 ) return


*-----------------------------------------------------------------------


                     san   =  float(m1) / sanf
                     saf   = ( saa * e1**2 + sbb * e1 + scc ) * san


                     qvaln = bndeng(iz1,in1) - bndeng(iz1,in1-1)
                     ss    = e1 - qvaln


                     barf  = sdmfisb(m1,iz1)


               if( e1 - qvaln .gt. 0.0 .and.
     &             e1 - barf  .gt. 0.0 ) then


                        ssn  = 2.0 * sqrt( san * ( e1 - qvaln ) )
                        ssf  = 2.0 * sqrt( saf * ( e1 - barf  ) )


                        ssfn = ssf - ssn


                     if( ssfn .gt. -ddct ) ssfn = -ddct


                        ratfn = 0.0


                     if( ssfn .ge. ddct ) then


                        ratfn = rk0 * san * ( ssf - 1.0 )
     &                        / ( 4.0 * float(m1)**(2./3.)
     &                          * saf * ( e1 - qvaln )  )
     &                        * exp( ssfn )
                     end if


                        if( ratfn .gt. 0.0 ) tcsf = tcn * ratfn


               end if


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine sdmfiss(iz1,in1,e1,j1,
     &                   iz2,in2,e2,j2,
     &                   iz3,in3,e3,j3,
     &                   barf,saf,ierr)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to determine the mass and charge of fission fragment    *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              iz1         : charge of fragment                        *
*              in1         : nuetron number of fragment                *
*              e1          : excitation energy                         *
*              j1          : angular momentum                          *
*                                                                      *
*              (iz2,in2,e2,i2)                                         *
*              (iz3,in3,e3,j3)  : fragments' quantities                *
*                                                                      *
*              barf        : fission barrier                           *
*              saf         : level density parameter of fission        *
*              ierr        : error flag                                *
*                                                                      *
*                                                                      *
************************************************************************


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


      parameter      ( ddct = - 30.0 )
      parameter      ( bnorm = 1.2011 )
      parameter      ( mcrit = 1000 )


*-----------------------------------------------------------------------
*           Error Frag
*-----------------------------------------------------------------------


               ierr = 0


*-----------------------------------------------------------------------
*           temperature and width of the mass distribution
*-----------------------------------------------------------------------


            m1   = iz1 + in1
            rm1  = float(m1)
            rz1  = float(iz1)


c           temp = sqrt( ( e1 - barf ) / saf )
            temp = sqrt( abs( ( e1 - barf ) / saf ) )


            width = ( e1 - barf + 7.0 ) * bnorm


*     -------------------------------------------------------
*     this width is questionable for high excitation energy !
*     -------------------------------------------------------


*-----------------------------------------------------------------------
*     mass of fission fragment for subactinides
*-----------------------------------------------------------------------


      if( iz1 .lt. 89 ) then


            imsa = 0
   50       imsa = imsa + 1


            if( imsa .lt. mcrit ) then


               rm2  = rn() * ( rm1 - 1.0 )
               rarg = - ( rm2 - rm1/2.0 )**2 / width**2
               m2   = nint( rm2 )


               if( rarg .lt. ddct )            goto 50
               if( rn() .gt. exp( rarg ) )    goto 50
               if( m2 .eq. 0 .or. m2 .eq. m1 ) goto 50


            else


               m2 = m1 / 2


            end if




*-----------------------------------------------------------------------
*     mass of fission fragment for actinides
*-----------------------------------------------------------------------


      else


            bm1 = 0.4 * rm1
            bm2 = 0.5 * rm1
            bm3 = 0.6 * rm1


            exm = e1 + 6.0


            alpa = 19.98160
            beta = 78.61184


         if( exm .le. 25.0 ) then


            alpa = exp( 0.5991 * exm - 13.1869 )
            beta = exp( 0.7013 * exm - 17.5325 )


         else if( exm .le. 40.0 ) then


            alpa = exp( 0.2008 * exm**0.8 - 0.8451 )
            beta = exp( 2.2672 * sqrt(exm)- 11.3431 )


         else if( exm .le. 48.0 ) then


            beta = exp( 2.2672 * sqrt(exm)- 11.3431 )


         end if


            rdm1 = alpa + beta * exp( - ( bm1 - bm2 )**2 / width**2 )
     &           + alpa * exp( - ( bm1 - bm3 )**2 / width**2 )


            rdm2 = beta
     &           + 2.0 * alpa * exp( - ( bm1 - bm2 )**2 / width**2 )


            rdmax = max( rdm1, rdm2 ) * 1.1




            imsa = 0
   60       imsa = imsa + 1


            if( imsa .lt. mcrit ) then


               rm2  = rn() * ( rm1 - 1.0 )


                  ear1 = 0.0
                  ear2 = 0.0
                  ear3 = 0.0


                  rar1 = - ( rm2 - bm1 )**2 / width**2
                  if( rar1 .gt. ddct ) ear1 = exp( rar1 )


                  rar2 = - ( rm2 - bm2 )**2 / width**2
                  if( rar2 .gt. ddct ) ear2 = exp( rar2 )


                  rar3 = - ( rm2 - bm3 )**2 / width**2
                  if( rar3 .gt. ddct ) ear3 = exp( rar3 )


               m2 = nint( rm2 )


               if( rn() * rdmax .gt. ear1 + ear2 + ear3 ) goto 60
               if( m2 .eq. 0 .or. m2 .eq. m1 )             goto 60


          else


               m2 = m1 / 2


          end if


      end if


*-----------------------------------------------------------------------
*     deetermine the charge of fission fragment
*-----------------------------------------------------------------------


         rho  = 1.1
         r0   = 1.2
         bets = -31.4506
         phi  = 44.2355


         s    = 0.1 * ccoul * 1000.0 / r0 / bets
     &          * ( rm1 / 2. )**(2./3.) * ( 1.0 - 5. / 8. / rho )


         zbar = - s * rz1 / 2. + ( 1. + s ) * float( m2 ) * rz1 / rm1




         widi = - 16.0 * bets / rm1 / temp
     &        * ( 1. + phi / bets * ( 2. / rm1 )**(1./3.)
     &            - 0.055 * ccoul * 1000.0 / r0 / bets * rm1**(2./3.) )


            imsa = 0
   70       imsa = imsa + 1


            if( imsa .lt. mcrit ) then


               rz2  = rn() * ( rz1 - 1.0 )
               rarg = - ( rz2 - zbar )**2 * widi


               iz2  = nint( rz2 )


               if( rarg .lt. ddct )         goto 70
               if( rn() .gt. exp( rarg ) ) goto 70
               if( iz2 .eq.  0 .or.
     &             iz2 .ge. m2 .or.
     &             iz2 .eq. iz1 )           goto 70


          else


               iz2 = min( m2 / 2 , iz1 / 2 )


          end if




*-----------------------------------------------------------------------
*     now fission fragments are determined
*-----------------------------------------------------------------------


            m3  = m1 - m2


            iz2 = iz2 - 1
   80       iz2 = iz2 + 1


            if( iz2 .lt. m2 .and. iz2 .le. iz1 ) then


               iz3 = iz1 - iz2
               if( iz3 .ge. m3 ) goto 80


            else


               write(*,*) ' **** Error at sdmfiss,',
     &                    ' iz3 can not be determined'
               write(*,*) ' ======================',
     &                    '=========================='
               write(*,*) ' iz1 = ', iz1, ' in1 = ', in1


               ierr = 1
               return


            end if


            in2 = m2 - iz2
            in3 = m3 - iz3


            rm = float( m2 * m3 ) / float( m1 )




*-----------------------------------------------------------------------
*     determine the kinetic energy and excitation energy
*-----------------------------------------------------------------------


            rad2 = 1.2 * float(m2)**(1./3.)
            rad3 = 1.2 * float(m3)**(1./3.)


            rad  = rad2 + rad3


            amoi2  = 2./5. * m2 * rad2**2
            amoi3  = 2./5. * m3 * rad3**2
            amoi23 = rm * rad**2 + amoi2 + amoi3


            xj2 = j1 * amoi2 / amoi23
            xj3 = j1 * amoi3 / amoi23


            j2  = xj2
            j3  = xj3


*-----------------------------------------------------------------------
*           rotation energy and relative motion
*-----------------------------------------------------------------------


            er = ( j1 - xj2 - xj3)**2 / 2.0
     &         / ( amoi23 - amoi2 - amoi3 )


            er = er + 0.1071 * rz1**2 / rm1**(1./3.) + 22.2


            qvalf = bndeng(iz1,in1)
     &            - bndeng(iz2,in2)
     &            - bndeng(iz3,in3)


            etf =  max( 0.0, e1 - qvalf )
            efx = etf - er


            if( efx .le. 0.0 ) then


                er  = etf
                efx = 0.0


            end if


            e2 = efx / rm1 * m2
            e3 = efx / rm1 * m3


*-----------------------------------------------------------------------


       return
       end




************************************************************************
*                                                                      *
      function sdmfisb(m1,iz1)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to determine the fission barrier.                       *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              m1          : mass of parent nucleus                    *
*              iz1         : charge of parent nucleus                  *
*              sdmfisb     : fission barrier                           *
*                                                                      *
*                                                                      *
*        Comments:                                                     *
*                                                                      *
*              These double-humped fission barrier data are            *
*              taken from  Kupriyanof et al.                           *
*              Sov. J. Nucl. Phys. 32 (1980) 184                       *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      dimension  fba(64), fbb(64)


*-----------------------------------------------------------------------


      data fba / 5.69, 5.50, 5.08, 5.07, 5.69, 5.42, 5.26, 5.45, 5.59,
     &           5.42, 5.48, 5.45, 5.41, 5.29, 5.24, 5.57, 5.45, 5.65,
     &           5.48, 5.67, 5.79, 5.97, 5.92, 5.91, 5.83, 5.89, 5.65,
     &           5.63, 5.42, 5.69, 6.10, 6.02, 6.17, 5.96, 6.17, 5.82,
     &           6.25, 6.22, 6.40, 6.16, 6.17, 5.94, 5.92, 5.71, 5.67,
     &           6.40, 6.59, 6.34, 6.44, 6.09, 6.26, 5.82, 5.92, 5.37,
     &           6.56, 6.45, 6.53, 6.41, 6.54, 6.32, 6.32, 6.10, 5.89,
     &           5.48 /


      data fbb / 7.89, 7.68, 7.30, 7.34, 7.35, 7.14, 7.04, 6.58, 6.79,
     &           6.68, 6.80, 6.84, 6.86, 6.79, 6.80, 6.26, 6.21, 6.48,
     &           6.38, 5.75, 5.95, 6.20, 6.23, 6.29, 6.28, 6.40, 6.23,
     &           6.26, 6.12, 5.21, 5.69, 5.68, 5.93, 5.79, 6.08, 5.79,
     &           5.34, 5.39, 5.65, 5.48, 5.46, 5.41, 5.52, 5.32, 5.34,
     &           4.87, 5.15, 4.98, 5.16, 4.89, 5.13, 4.77, 4.94, 4.45,
     &           4.50, 4.38, 4.54, 4.50, 4.72, 4.57, 4.65, 4.50, 4.36,
     &           4.02 /


*-----------------------------------------------------------------------


      if( m1 .le. 90 ) then


            sdmfisb = 52.0 * exp( - ( ( m1 - 90.0 ) / 84.7  )**2 )


      else if( m1 .le. 200 ) then


            sdmfisb = 52.0 * exp( - ( ( m1 - 90.0 ) / 110.0 )**2 )


      else if( m1 .le. 224 ) then


            sdmfisb = 23.0 * exp( - ( ( m1 - 210.0 ) / 13.2 )**2 )


      else


            sdmfisb = 6.0


            ifis = 0


         if( iz1 .eq. 88 )      then


            if( m1 .le. 228 )                   ifis = 1  + m1 - 225


         else if( iz1 .eq. 89 ) then


            if( m1 .ge. 226 .and. m1 .le. 228 ) ifis = 5  + m1 - 226


         else if( iz1 .eq. 90 ) then


            if( m1 .ge. 227 .and. m1 .le. 234 ) ifis = 8  + m1 - 227


         else if( iz1 .eq. 91 ) then


            if( m1 .ge. 230 .and. m1 .le. 233 ) ifis = 16 + m1 - 230


         else if( iz1 .eq. 92 ) then


            if( m1 .ge. 231 .and. m1 .le. 240 ) ifis = 20 + m1 - 231


         else if( iz1 .eq. 93 ) then


            if( m1 .ge. 233 .and. m1 .le. 239 ) ifis = 30 + m1 - 233


         else if( iz1 .eq. 94 ) then


            if( m1 .ge. 237 .and. m1 .le. 245 ) ifis = 37 + m1 - 237


         else if( iz1 .eq. 95 ) then


            if( m1 .ge. 239 .and. m1 .le. 247 ) ifis = 46 + m1 - 239


         else if( iz1 .eq. 96 ) then


            if( m1 .ge. 241 .and. m1 .le. 250 ) ifis = 55 + m1 - 241


         end if




         if( ifis .ne. 0 ) then


            sdmfisb = max( fba(ifis), fbb(ifis) )


         end if




      end if


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine sdmlist(io,icd)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to list event record and particle data                  *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              io          : output unit                               *
*              icd         : = 1 ; write information kclust            *
*                            = 0 ; no information kclust               *
*                                                                      *
*                                                                      *
************************************************************************


      include 'param01.inc'


*-----------------------------------------------------------------------


      common /clusts/ nclust, kclust(3,nnn)
      common /clustu/ lclust(0:7,nnn), sclust(0:7,nnn)
      common /clustv/ kdecay(4)


*-----------------------------------------------------------------------


      character chau*16


*-----------------------------------------------------------------------
*     Output nuclear cluster
*-----------------------------------------------------------------------


                  ntt = 0
                  nzt = 0
                  nnt = 0


            do i = 1, nclust


               if( kclust(1,i) .ge. 100 ) then


                  ntt = ntt + 1
                  nzt = nzt + lclust(1,i)
                  nnt = nnt + lclust(2,i)


               end if


            end do


*-----------------------------------------------------------------------


            write(io,'(/72(''=''))')


            write(io,'(/
     &         ''   Event Summary of SDM '')')


            write(io,'(/
     &         ''   Final total products = '',i3)') ntt


            write(io,'(
     &         ''   Evapolation events   = '',i3)') kdecay(1)
            write(io,'(
     &         ''   Gamma decay          = '',i3)') kdecay(2)
            write(io,'(
     &         ''   Ground state decay   = '',i3)') kdecay(3)
            write(io,'(
     &         ''   Fission events       = '',i3)') kdecay(4)


*-----------------------------------------------------------------------


         if( icd .eq. 1 ) then


            write(io,'(/72(''=''))')


            write(io,'(/
     &         ''   kc1 =     0 : initial mother''/
     &         ''       = + 100 : final products''/
     &         ''       = +  10 : fission product''/
     &         ''       =     1 : evapolation particle''/
     &         ''       =     2 : evapolation residue''/
     &         ''       =     3 : gamma ray''/
     &         ''       =     4 : gamma decay residue''/
     &         ''       =     5 : ground state decay particle''/
     &         ''       =     6 : ground state decay residue'')')


            write(io,'(/
     &         ''   kc2 : line number of parent nucleus'')')


            write(io,'(/
     &         ''   kc3 : rank'')')


         end if


*-----------------------------------------------------------------------


            write(io,'(/72(''=''))')


            write(io,'(/
     &         ''   i  kc1 kc2 kc3     Name       z    n   jj'',
     &         ''      Ex(MeV)       Ek(MeV)''/)')


*-----------------------------------------------------------------------


         do ii = 1, 3


                  nprint = 0


            do i = 1, nclust


               if     ( ii .eq. 1 .and.
     &                ( kclust(1,i) .eq. 0 .or.
     &                  kclust(1,i) .eq. 100 ) ) then


                  iprint = 1


               else if( ii .eq. 2 .and.
     &                  kclust(1,i) .gt.   0 .and.
     &                  kclust(1,i) .lt. 100 ) then


                  iprint = 1


               else if( ii .eq. 3 .and.
     &                  kclust(1,i) .gt. 100 ) then


                  iprint = 1


               else


                  iprint = 0


               end if


*-----------------------------------------------------------------------


               if( iprint .eq. 1 ) then


                  nprint = nprint + 1


                     jj = lclust(0,i)
                     np = lclust(1,i)
                     nn = lclust(2,i)


                     mas = nn + np


                  if( mas .eq. 0 ) then


                     mnds = -1


                  else if( mas .eq. 1 ) then


                     mnds = 1


                  else


                     mnds = 0


                  end if


                     call chname(mnds,mas,np,chau,len)




                  if( mas .gt. 1 .and.
     &                sclust(6,i) .gt. 0.001 ) then


                     chau = chau(1:len)//'*'


                  end if


                     ekin = ( sclust(4,i) - sclust(5,i) ) * 1000.0
                     ex   = sclust(6,i)


                  write(io,'(i4,1x,i4,1x,i3,1x,i3,3x,
     &                       a9,2x,i3,2x,i3,2x,i3,1x,
     &                      2(1x,3g14.5))')
     &                  i, kclust(1,i), kclust(2,i), kclust(3,i),
     &                  chau(1:10), np, nn, jj,
     &                  ex, ekin




               end if


*-----------------------------------------------------------------------


            end do


               if( nprint .ge. 1 ) write(io,'(72(''=''))')


         end do


               write(io,'(3x,''Total'',23x,i3,2x,i3)') nzt, nnt




*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      function bndeng(nz,nn)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to give binding energy per baryon (MeV)                 *
*              Experimental Data and                                   *
*              Liquid Drop Mass Formula (BW) are used.                 *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              nz          : proton  number                            *
*              nn          : neutron number                            *
*                                                                      *
************************************************************************


      include 'param01.inc'


*-----------------------------------------------------------------------


      common /sdmsw1/ iswids, isevap, isfiss, isgrnd


      common /masmx0/ imax0


      common /maspn0/ nzz0(0:nnuc), nnn0(0:nnuc), jgs0(0:nnuc)
      common /masbe0/ be0(0:nnuc)
      common /masid0/ inuc0(0:maxpt,0:maxnt)


*-----------------------------------------------------------------------


               bndeng = 0.0
               mass   = nz + nn


*-----------------------------------------------------------------------
*        check
*-----------------------------------------------------------------------


            if( nz  .lt. 0 .or.
     &          nn  .lt. 0 .or.
     &         mass .le. 1 ) then


               return


            end if


*-----------------------------------------------------------------------
*     sdmwid0 case
*-----------------------------------------------------------------------


         if( iswids .eq. 0 ) then




                  bndeng = eliq(nz,nn)




*-----------------------------------------------------------------------
*     sdmwid1 case
*-----------------------------------------------------------------------


         else if( iswids .eq. 1 ) then


*-----------------------------------------------------------------------
*        Super Heavy Nuclei or out of table
*        ---> Use Liquid Drop Mass Formula
*-----------------------------------------------------------------------


            if( nz .gt. maxpt .or.
     &          nn .gt. maxnt ) then


                  bndeng = eliq(nz,nn)


            else if( inuc0(nz,nn) .eq. 0 ) then


                  bndeng = eliq(nz,nn)


*-----------------------------------------------------------------------
*        Normal Nucleus
*-----------------------------------------------------------------------


            else


                  bndeng = be0( inuc0(nz,nn) )


            end if


         end if


*-----------------------------------------------------------------------


      return
      end






************************************************************************
*                                                                      *
      function eliq(nz,nn)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate liquid drop binding energy ( MeV )         *
*                                                                      *
*        Variables:                                                    *
*              nz    :   proton number                                 *
*              nn    :   neutron number                                *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      parameter( bvol=15.56, bsur=17.23, bsym=46.57, coul=1.2 )


*-----------------------------------------------------------------------


         eliq = 0.0


         if( nn + nz .le. 1 ) return


         mm = nn + nz
         a  = 1.0 * nn + nz
         rc = 1.24 * a**(1./3.)


         eliq = ( bvol * a
     &          - bsur * a**(2./3.)
     &          - bsym / 2. * ( nn - nz ) * ( nn - nz ) / a
     &          - 3./5. * nz * coul * nz * coul / rc )


*-----------------------------------------------------------------------
*        special nucleus
*-----------------------------------------------------------------------


c        if( mm.eq. 2 .and. nz.eq. 1 ) eliq = 2.224
c        if( mm.eq. 3 .and. nz.eq. 1 ) eliq = 8.482
c        if( mm.eq. 3 .and. nz.eq. 2 ) eliq = 7.718
c        if( mm.eq. 4 .and. nz.eq. 1 ) eliq = 3.280
c        if( mm.eq. 4 .and. nz.eq. 2 ) eliq = 28.296
c        if( mm.eq. 4 .and. nz.eq. 3 ) eliq = 4.8
c        if( mm.eq. 5 .and. nz.eq. 2 ) eliq = 27.338
c        if( mm.eq. 5 .and. nz.eq. 3 ) eliq = 26.331
c        if( mm.eq. 6 .and. nz.eq. 4 ) eliq = 26.932
c        if( mm.eq. 6 .and. nz.eq. 2 ) eliq = 29.265
c        if( mm.eq. 6 .and. nz.eq. 3 ) eliq = 31.992
c        if( mm.eq. 8 .and. nz.eq. 4 ) eliq = 56.497
c        if( mm.eq.11 .and. nz.eq. 3 ) eliq = 45.54
c        if( mm.eq.12 .and. nz.eq. 6 ) eliq = 92.162
c        if( mm.eq.16 .and. nz.eq. 8 ) eliq = 127.62


         if( mm.eq. 2 .and. nz.eq. 1 ) eliq = 1.11 * a
         if( mm.eq. 3 .and. nz.eq. 1 ) eliq = 2.83 * a
         if( mm.eq. 3 .and. nz.eq. 2 ) eliq = 2.57 * a
         if( mm.eq. 4 .and. nz.eq. 1 ) eliq = 1.40 * a
         if( mm.eq. 4 .and. nz.eq. 2 ) eliq = 7.07 * a
         if( mm.eq. 4 .and. nz.eq. 3 ) eliq = 1.20 * a
         if( mm.eq. 5 .and. nz.eq. 2 ) eliq = 5.48 * a
         if( mm.eq. 5 .and. nz.eq. 3 ) eliq = 5.27 * a
         if( mm.eq. 6 .and. nz.eq. 4 ) eliq = 4.49 * a
         if( mm.eq. 6 .and. nz.eq. 2 ) eliq = 4.88 * a
         if( mm.eq. 6 .and. nz.eq. 3 ) eliq = 5.33 * a
         if( mm.eq. 8 .and. nz.eq. 4 ) eliq = 7.06 * a
         if( mm.eq.11 .and. nz.eq. 3 ) eliq = 4.14 * a
         if( mm.eq.12 .and. nz.eq. 6 ) eliq = 7.68 * a
         if( mm.eq.16 .and. nz.eq. 8 ) eliq = 7.98 * a




*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      block data evtable
*                                                                      *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              data table of evapolation particle                      *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      common /evapol/ ievmax, izip(0:100), inip(0:100), jip(0:100)
      common /evaval/ sepc(2,0:100)


*-----------------------------------------------------------------------
*     data for evapolation particle
*-----------------------------------------------------------------------


      data ievmax / 6 /


      data ( izip(i), i = 1, 6 ) / 0,1,1,1,2,2 /
      data ( inip(i), i = 1, 6 ) / 1,0,1,2,1,2 /
      data (  jip(i), i = 1, 6 ) / 0,0,1,0,0,0 /


      data ( sepc(1,i), sepc(2,i), i = 1, 6 ) /
     &            1.0000, 1.00000,  !neutron DUMMY
     &            8.0606, 0.50836,  !proton
     &            7.9869, 0.62142,  !deutron
     &            7.9869, 0.62142,  !trition
     &            6.4355, 0.78601,  !3-helium
     &            6.1333, 0.88761/  !alpha particle


*-----------------------------------------------------------------------


      end




************************************************************************
*                                                                      *
      subroutine sdmtable
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to read mass table from block data                      *
*                                                                      *
*                                                                      *
************************************************************************


      include 'param01.inc'


*-----------------------------------------------------------------------


      common /masmx0/ imax0


      common /maspn0/ nzz0(0:nnuc), nnn0(0:nnuc), jgs0(0:nnuc)
      common /masbe0/ be0(0:nnuc)
      common /masid0/ inuc0(0:maxpt,0:maxnt)


      common /grnddc/ idgrnd(0:nnuc)


      common /evapol/ ievmax, izip(0:100), inip(0:100), jip(0:100)


*-----------------------------------------------------------------------


                  idgrnd(0) = 0


               do iz = 0, maxpt
               do in = 0, maxnt


                  inuc0(iz,in) = 0


               end do
               end do




*-----------------------------------------------------------------------
*           Normal Nucleus data from block data
*-----------------------------------------------------------------------


               do i = 1, imax0


                  iz = nzz0(i)
                  in = nnn0(i)


                  inuc0(iz,in) = i


                  idgrnd(i) = 0


               end do


*-----------------------------------------------------------------------
*     table of grand state decay nucleus
*
*               4H  ->  3H  + n
*               4Li ->  3He + p
*               5H  ->  4H  + n
*               5He ->  4He + n
*               5Li ->  4He + p
*               7B  ->  6Be + p
*               7He ->  6He + n
*              10Li ->  9Li + n
*              13Be -> 12Be + n
*              16B  -> 15B  + n
*
*               6H  ->  5H  + n
*               7H  ->  6H  + n
*               8H  ->  7H  + n
*
*               5Be ->  4Li + p
*               6Be ->  5Li + p
*               6B  ->  5Be + p
*               7C  ->  6B  + p
*               9He ->  8He + p
*
*-----------------------------------------------------------------------


               idgrnd( 6) =  4
               idgrnd( 8) =  5
               idgrnd( 9) =  6
               idgrnd(10) =  7
               idgrnd(11) =  7
               idgrnd(15) = 12
               idgrnd(17) = 13
               idgrnd(31) = 27
               idgrnd(44) = 39
               idgrnd(59) = 54


               idgrnd(1993) = 9
               idgrnd(1994) = 1993
               idgrnd(1995) = 1994


               idgrnd(1996) = 8
               idgrnd(1997) = 11
               idgrnd(1998) = 1996
               idgrnd(1999) = 1998
               idgrnd(2000) = 22




*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      block data mstable
*                                                                      *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              mass table for the normal                               *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'


*-----------------------------------------------------------------------


      common /masmx0/ imax0


*-----------------------------------------------------------------------


      common /maspn0/ nzz0(0:nnuc), nnn0(0:nnuc), jgs0(0:nnuc)
      common /masbe0/ be0(0:nnuc)


*-----------------------------------------------------------------------


      data imax0 / 2000 /


*-----------------------------------------------------------------------
*     Normal Nucleus
*-----------------------------------------------------------------------


      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=   1,  10)/
     &   0,  1,   0.000,  2 ! n
     &,  1,  0,   0.000,  2 ! p
     &,  1,  1,   2.225,  3 ! 2H
     &,  1,  2,   8.483,  2 ! 3H
     &,  2,  1,   7.719,  2 ! 3He
     &,  1,  3,   5.584,  5 ! 4H
     &,  2,  2,  28.297,  1 ! 4He
     &,  3,  1,   4.809,  5 ! 4Li
     &,  1,  4,   5.786,  2 ! 5H
     &,  2,  3,  27.404,  4/! 5He
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=  11,  20)/
     &   3,  2,  26.331,  4 ! 5Li
     &,  4,  2,  26.925,  1 ! 6Be
     &,  2,  4,  29.268,  1 ! 6He
     &,  3,  3,  31.995,  3 ! 6Li
     &,  5,  2,  24.649,  4 ! 7B
     &,  4,  3,  37.602,  4 ! 7Be
     &,  2,  5,  28.826,  4 ! 7He
     &,  3,  4,  39.246,  4 ! 7Li
     &,  5,  3,  37.739,  5 ! 8B
     &,  4,  4,  56.502,  1/! 8Be
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=  21,  30)/
     &   6,  2,  24.794,  1 ! 8C
     &,  2,  6,  31.400,  1 ! 8He
     &,  3,  5,  41.279,  5 ! 8Li
     &,  5,  4,  56.317,  4 ! 9B
     &,  4,  5,  58.167,  4 ! 9Be
     &,  6,  3,  39.038,  2 ! 9C
     &,  3,  6,  45.343,  2 ! 9Li
     &,  5,  5,  64.753,  7 ! 10B
     &,  4,  6,  64.979,  1 ! 10Be
     &,  6,  4,  60.319,  1/! 10C
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=  31,  40)/
     &   3,  7,  44.539,  1 ! 10Li
     &,  7,  3,  35.740,  1 ! 10N
     &,  5,  6,  76.208,  4 ! 11B
     &,  4,  7,  65.483,  2 ! 11Be
     &,  6,  5,  73.444,  4 ! 11C
     &,  3,  8,  45.501,  2 ! 11Li
     &,  7,  4,  58.081,  2 ! 11N
     &,  5,  7,  79.578,  3 ! 12B
     &,  4,  8,  68.700,  1 ! 12Be
     &,  6,  6,  92.165,  1/! 12C
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=  41,  50)/
     &   7,  5,  74.045,  3 ! 12N
     &,  8,  4,  58.531,  1 ! 12O
     &,  5,  8,  84.458,  2 ! 13B
     &,  4,  9,  66.902,  2 ! 13Be
     &,  6,  7,  97.112,  2 ! 13C
     &,  7,  6,  94.109,  2 ! 13N
     &,  8,  5,  75.567,  2 ! 13O
     &,  5,  9,  85.434,  1 ! 14B
     &,  4, 10,  68.904,  1 ! 14Be
     &,  6,  8, 105.289,  1/! 14C
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=  51,  60)/
     &   9,  5,  72.351,  1 ! 14F
     &,  7,  7, 104.663,  3 ! 14N
     &,  8,  6,  98.736,  1 ! 14O
     &,  5, 10,  87.633,  2 ! 15B
     &,  6,  9, 106.507,  2 ! 15C
     &,  9,  6,  96.373,  2 ! 15F
     &,  7,  8, 115.497,  2 ! 15N
     &,  8,  7, 111.960,  2 ! 15O
     &,  5, 11,  87.235,  1 ! 16B
     &,  6, 10, 110.759,  1/! 16C
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=  61,  70)/
     &   9,  7, 111.413,  5 ! 16F
     &,  7,  9, 117.988,  5 ! 16N
     &, 10,  6,  97.212,  1 ! 16Ne
     &,  8,  8, 127.624,  1 ! 16O
     &,  5, 12,  88.036,  2 ! 17B
     &,  6, 11, 111.464,  2 ! 17C
     &,  9,  8, 128.225,  6 ! 17F
     &,  7, 10, 123.871,  2 ! 17N
     &, 10,  7, 112.916,  2 ! 17Ne
     &,  8,  9, 131.769,  6/! 17O
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=  71,  80)/
     &   6, 12, 115.226,  1 ! 18C
     &,  9,  9, 137.376,  3 ! 18F
     &,  7, 11, 126.539,  1 ! 18N
     &, 11,  7, 111.363,  1 ! 18Na
     &, 10,  8, 132.147,  1 ! 18Ne
     &,  8, 10, 139.814,  1 ! 18O
     &,  6, 13, 114.237,  2 ! 19C
     &,  9, 10, 147.807,  2 ! 19F
     &,  7, 12, 132.285,  2 ! 19N
     &, 11,  8, 131.825,  2/! 19Na
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=  81,  90)/
     &  10,  9, 143.787,  2 ! 19Ne
     &,  8, 11, 143.771,  2 ! 19O
     &,  9, 11, 154.409,  1 ! 20F
     &, 12,  8, 134.476,  1 ! 20Mg
     &,  7, 13, 133.757,  1 ! 20N
     &, 11,  9, 145.983,  1 ! 20Na
     &, 10, 10, 160.652,  1 ! 20Ne
     &,  8, 12, 151.375,  1 ! 20O
     &,  9, 12, 162.510,  2 ! 21F
     &, 12,  9, 149.204,  2/! 21Mg
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=  91, 100)/
     &   7, 14, 137.078,  2 ! 21N
     &, 11, 10, 163.084,  4 ! 21Na
     &, 10, 11, 167.414,  4 ! 21Ne
     &,  8, 13, 155.126,  2 ! 21O
     &, 13,  9, 149.195,  1 ! 22Al
     &,  9, 13, 167.709,  1 ! 22F
     &, 12, 10, 168.582,  1 ! 22Mg
     &, 11, 11, 174.154,  7 ! 22Na
     &, 10, 12, 177.779,  1 ! 22Ne
     &,  8, 14, 161.827,  1/! 22O
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 101, 110)/
     &  13, 10, 168.709,  2 ! 23Al
     &,  9, 14, 175.257,  2 ! 23F
     &, 12, 11, 181.730,  4 ! 23Mg
     &, 11, 12, 186.571,  4 ! 23Na
     &, 10, 13, 182.979,  2 ! 23Ne
     &,  8, 15, 161.439,  2 ! 23O
     &, 13, 11, 183.600,  9 ! 24Al
     &,  9, 15, 178.028,  1 ! 24F
     &, 12, 12, 198.262,  1 ! 24Mg
     &, 11, 13, 193.531,  9/! 24Na
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 111, 120)/
     &  10, 14, 191.845,  1 ! 24Ne
     &, 14, 10, 172.026,  1 ! 24Si
     &, 13, 12, 200.533,  6 ! 25Al
     &,  9, 16, 181.910,  2 ! 25F
     &, 12, 13, 205.593,  6 ! 25Mg
     &, 11, 14, 202.542,  6 ! 25Na
     &, 10, 15, 196.118,  2 ! 25Ne
     &, 14, 11, 187.014,  6 ! 25Si
     &, 13, 13, 211.899, 11 ! 26Al
     &, 12, 14, 216.687,  1/! 26Mg
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 121, 130)/
     &  11, 15, 208.145,  1 ! 26Na
     &, 10, 16, 202.229,  1 ! 26Ne
     &, 15, 11, 186.867,  1 ! 26P
     &, 14, 12, 206.052,  1 ! 26Si
     &, 13, 14, 224.958,  6 ! 27Al
     &, 12, 15, 223.131,  2 ! 27Mg
     &, 11, 16, 214.958,  2 ! 27Na
     &, 10, 17, 203.441,  2 ! 27Ne
     &, 15, 12, 206.789,  2 ! 27P
     &, 14, 13, 219.366,  6/! 27Si
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 131, 140)/
     &  13, 15, 232.683,  1 ! 28Al
     &, 12, 16, 231.634,  1 ! 28Mg
     &, 11, 17, 218.530,  1 ! 28Na
     &, 15, 13, 221.430,  1 ! 28P
     &, 16, 12, 209.298,  1 ! 28S
     &, 14, 14, 236.544,  1 ! 28Si
     &, 13, 16, 242.119,  2 ! 29Al
     &, 12, 17, 235.439,  2 ! 29Mg
     &, 11, 18, 222.812,  2 ! 29Na
     &, 15, 14, 239.291,  2/! 29P
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 141, 150)/
     &  16, 13, 224.719,  2 ! 29S
     &, 14, 15, 245.018,  2 ! 29Si
     &, 13, 17, 247.869,  1 ! 30Al
     &, 17, 13, 224.009,  1 ! 30Cl
     &, 12, 18, 242.551,  1 ! 30Mg
     &, 11, 19, 225.164,  1 ! 30Na
     &, 15, 15, 250.618,  1 ! 30P
     &, 16, 14, 243.693,  1 ! 30S
     &, 14, 16, 255.628,  1 ! 30Si
     &, 13, 18, 255.150,  2/! 31Al
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 151, 160)/
     &  17, 14, 243.990,  2 ! 31Cl
     &, 12, 19, 244.733,  2 ! 31Mg
     &, 11, 20, 231.005,  2 ! 31Na
     &, 15, 16, 262.925,  2 ! 31P
     &, 16, 15, 256.747,  2 ! 31S
     &, 14, 17, 262.217,  2 ! 31Si
     &, 13, 19, 259.412,  1 ! 32Al
     &, 18, 14, 246.420,  1 ! 32Ar
     &, 17, 15, 258.321,  1 ! 32Cl
     &, 12, 20, 251.794,  1/! 32Mg
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 161, 170)/
     &  11, 21, 233.277,  1 ! 32Na
     &, 15, 17, 270.862,  1 ! 32P
     &, 16, 16, 271.790,  1 ! 32S
     &, 14, 18, 271.432,  1 ! 32Si
     &, 13, 20, 265.564,  2 ! 33Al
     &, 18, 15, 261.666,  2 ! 33Ar
     &, 17, 16, 274.067,  2 ! 33Cl
     &, 12, 21, 252.846,  2 ! 33Mg
     &, 15, 18, 280.966,  2 ! 33P
     &, 16, 17, 280.432,  2/! 33S
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 171, 180)/
     &  14, 19, 275.981,  2 ! 33Si
     &, 13, 21, 268.415,  1 ! 34Al
     &, 18, 16, 278.732,  1 ! 34Ar
     &, 17, 17, 285.574,  1 ! 34Cl
     &, 19, 15, 261.051,  1 ! 34K
     &, 15, 19, 287.250,  1 ! 34P
     &, 16, 18, 291.849,  1 ! 34S
     &, 14, 20, 283.733,  1 ! 34Si
     &, 13, 22, 273.177,  2 ! 35Al
     &, 18, 17, 291.474,  2/! 35Ar
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 181, 190)/
     &  17, 18, 298.221,  2 ! 35Cl
     &, 19, 16, 278.811,  2 ! 35K
     &, 15, 20, 295.712,  2 ! 35P
     &, 16, 19, 298.836,  2 ! 35S
     &, 14, 21, 286.595,  2 ! 35Si
     &, 18, 18, 306.728,  1 ! 36Ar
     &, 20, 16, 281.581,  1 ! 36Ca
     &, 17, 19, 306.801,  1 ! 36Cl
     &, 19, 17, 293.140,  1 ! 36K
     &, 15, 21, 299.614,  1/! 36P
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 191, 200)/
     &  16, 20, 308.727,  1 ! 36S
     &, 14, 22, 292.296,  1 ! 36Si
     &, 18, 19, 315.516,  2 ! 37Ar
     &, 20, 17, 296.167,  2 ! 37Ca
     &, 17, 20, 317.112,  2 ! 37Cl
     &, 19, 18, 308.585,  2 ! 37K
     &, 15, 22, 305.925,  2 ! 37P
     &, 16, 21, 313.041,  2 ! 37S
     &, 14, 23, 294.708,  2 ! 37Si
     &, 18, 20, 327.355,  1/! 38Ar
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 201, 210)/
     &  20, 18, 313.135,  1 ! 38Ca
     &, 17, 21, 323.220,  1 ! 38Cl
     &, 19, 19, 320.659,  1 ! 38K
     &, 15, 23, 309.547,  1 ! 38P
     &, 16, 22, 321.067,  1 ! 38S
     &, 21, 17, 294.752,  1 ! 38Sc
     &, 18, 21, 333.952,  2 ! 39Ar
     &, 20, 19, 326.429,  2 ! 39Ca
     &, 17, 22, 331.297,  2 ! 39Cl
     &, 19, 20, 333.735,  2/! 39K
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 211, 220)/
     &  15, 24, 315.359,  2 ! 39P
     &, 16, 23, 325.276,  2 ! 39S
     &, 21, 18, 312.444,  2 ! 39Sc
     &, 18, 22, 343.823,  1 ! 40Ar
     &, 20, 20, 342.065,  1 ! 40Ca
     &, 17, 23, 337.106,  1 ! 40Cl
     &, 19, 21, 341.536,  1 ! 40K
     &, 16, 24, 332.588,  1 ! 40S
     &, 21, 19, 326.963,  1 ! 40Sc
     &, 22, 18, 314.693,  1/! 40Ti
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 221, 230)/
     &  18, 23, 349.923,  2 ! 41Ar
     &, 20, 21, 350.428,  2 ! 41Ca
     &, 17, 24, 345.037,  2 ! 41Cl
     &, 19, 22, 351.632,  2 ! 41K
     &, 16, 25, 336.520,  2 ! 41S
     &, 21, 20, 343.151,  2 ! 41Sc
     &, 22, 19, 329.505,  2 ! 41Ti
     &, 18, 24, 359.347,  1 ! 42Ar
     &, 20, 22, 361.905,  1 ! 42Ca
     &, 17, 25, 350.129,  1/! 42Cl
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 231, 240)/
     &  19, 23, 359.167,  1 ! 42K
     &, 21, 21, 354.700,  1 ! 42Sc
     &, 22, 20, 346.919,  1 ! 42Ti
     &, 23, 19, 329.034,  1 ! 42V
     &, 18, 25, 364.978,  2 ! 43Ar
     &, 20, 23, 369.839,  2 ! 43Ca
     &, 17, 26, 356.921,  2 ! 43Cl
     &, 19, 24, 368.804,  2 ! 43K
     &, 21, 22, 366.836,  2 ! 43Sc
     &, 22, 21, 359.192,  2/! 43Ti
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 241, 250)/
     &  23, 20, 347.106,  2 ! 43V
     &, 18, 26, 373.341,  1 ! 44Ar
     &, 20, 24, 380.971,  1 ! 44Ca
     &, 24, 20, 349.875,  1 ! 44Cr
     &, 19, 25, 376.094,  1 ! 44K
     &, 21, 23, 376.533,  1 ! 44Sc
     &, 22, 22, 375.486,  1 ! 44Ti
     &, 23, 21, 361.008,  1 ! 44V
     &, 18, 27, 378.872,  2 ! 45Ar
     &, 20, 25, 388.386,  2/! 45Ca
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 251, 260)/
     &  24, 21, 363.907,  2 ! 45Cr
     &, 19, 26, 384.970,  2 ! 45K
     &, 21, 24, 387.861,  2 ! 45Sc
     &, 22, 23, 385.016,  2 ! 45Ti
     &, 23, 22, 377.108,  2 ! 45V
     &, 18, 28, 386.943,  1 ! 46Ar
     &, 20, 26, 398.787,  1 ! 46Ca
     &, 24, 22, 381.979,  1 ! 46Cr
     &, 19, 27, 391.851,  1 ! 46K
     &, 25, 21, 364.206,  1/! 46Mn
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 261, 270)/
     &  21, 25, 396.621,  1 ! 46Sc
     &, 22, 24, 398.206,  1 ! 46Ti
     &, 23, 23, 390.372,  1 ! 46V
     &, 20, 27, 406.063,  2 ! 47Ca
     &, 24, 23, 395.208,  2 ! 47Cr
     &, 19, 28, 400.200,  2 ! 47K
     &, 25, 22, 382.728,  2 ! 47Mn
     &, 21, 26, 407.268,  2 ! 47Sc
     &, 22, 25, 407.086,  2 ! 47Ti
     &, 23, 24, 403.374,  2/! 47V
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 271, 280)/
     &  20, 28, 416.008,  1 ! 48Ca
     &, 24, 24, 411.480,  1 ! 48Cr
     &, 19, 29, 404.794,  1 ! 48K
     &, 25, 23, 397.049,  1 ! 48Mn
     &, 21, 27, 415.507,  1 ! 48Sc
     &, 22, 26, 418.714,  1 ! 48Ti
     &, 23, 25, 413.917,  1 ! 48V
     &, 20, 29, 421.149,  2 ! 49Ca
     &, 24, 25, 422.063,  2 ! 49Cr
     &, 26, 23, 399.639,  2/! 49Fe
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 281, 290)/
     &  25, 24, 413.564,  2 ! 49Mn
     &, 21, 28, 425.636,  2 ! 49Sc
     &, 22, 27, 426.857,  2 ! 49Ti
     &, 23, 26, 425.473,  2 ! 49V
     &, 20, 30, 427.507,  1 ! 50Ca
     &, 24, 26, 435.063,  1 ! 50Cr
     &, 26, 24, 417.670,  1 ! 50Fe
     &, 25, 25, 426.648,  1 ! 50Mn
     &, 21, 29, 431.692,  1 ! 50Sc
     &, 22, 28, 437.802,  1/! 50Ti
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 291, 300)/
     &  23, 27, 434.807,  1 ! 50V
     &, 27, 24, 417.760,  2 ! 51Co
     &, 24, 27, 444.325,  2 ! 51Cr
     &, 26, 25, 431.540,  2 ! 51Fe
     &, 25, 26, 440.334,  2 ! 51Mn
     &, 21, 30, 438.444,  2 ! 51Sc
     &, 22, 29, 444.175,  2 ! 51Ti
     &, 23, 28, 445.858,  2 ! 51V
     &, 27, 25, 432.831,  1 ! 52Co
     &, 24, 28, 456.364,  1/! 52Cr
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 301, 310)/
     &  26, 26, 447.716,  1 ! 52Fe
     &, 25, 27, 450.870,  1 ! 52Mn
     &, 21, 31, 443.436,  1 ! 52Sc
     &, 22, 30, 451.983,  1 ! 52Ti
     &, 23, 29, 453.170,  1 ! 52V
     &, 27, 26, 449.313,  2 ! 53Co
     &, 24, 29, 464.304,  2 ! 53Cr
     &, 26, 27, 458.400,  2 ! 53Fe
     &, 25, 28, 462.925,  2 ! 53Mn
     &, 28, 25, 435.300,  2/! 53Ni
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 311, 320)/
     &  22, 31, 457.475,  2 ! 53Ti
     &, 23, 30, 461.666,  2 ! 53V
     &, 27, 27, 462.754,  1 ! 54Co
     &, 24, 30, 474.023,  1 ! 54Cr
     &, 26, 28, 471.778,  1 ! 54Fe
     &, 25, 29, 471.864,  1 ! 54Mn
     &, 28, 26, 453.172,  1 ! 54Ni
     &, 22, 32, 463.987,  1 ! 54Ti
     &, 23, 31, 467.804,  1 ! 54V
     &, 27, 28, 476.840,  2/! 55Co
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 321, 330)/
     &  24, 31, 480.270,  2 ! 55Cr
     &, 29, 26, 452.781,  2 ! 55Cu
     &, 26, 29, 481.077,  2 ! 55Fe
     &, 25, 30, 482.091,  2 ! 55Mn
     &, 28, 27, 467.368,  2 ! 55Ni
     &, 23, 32, 474.956,  2 ! 55V
     &, 27, 29, 486.925,  1 ! 56Co
     &, 24, 32, 488.500,  1 ! 56Cr
     &, 29, 27, 467.823,  1 ! 56Cu
     &, 26, 30, 492.275,  1/! 56Fe
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 331, 340)/
     &  25, 31, 489.362,  1 ! 56Mn
     &, 28, 28, 484.007,  1 ! 56Ni
     &, 23, 33, 480.228,  1 ! 56V
     &, 27, 30, 498.302,  2 ! 57Co
     &, 24, 33, 494.097,  2 ! 57Cr
     &, 29, 28, 485.015,  2 ! 57Cu
     &, 26, 31, 499.921,  2 ! 57Fe
     &, 25, 32, 498.012,  2 ! 57Mn
     &, 28, 29, 494.276,  2 ! 57Ni
     &, 30, 27, 469.242,  2/! 57Zn
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 341, 350)/
     &  27, 31, 506.875,  1 ! 58Co
     &, 24, 34, 501.429,  1 ! 58Cr
     &, 29, 29, 497.128,  1 ! 58Cu
     &, 26, 32, 509.966,  1 ! 58Fe
     &, 25, 33, 504.806,  1 ! 58Mn
     &, 28, 30, 506.473,  1 ! 58Ni
     &, 30, 28, 486.944,  1 ! 58Zn
     &, 27, 32, 517.329,  2 ! 59Co
     &, 29, 30, 509.890,  2 ! 59Cu
     &, 26, 33, 516.547,  2/! 59Fe
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 351, 360)/
     &  25, 34, 512.146,  2 ! 59Mn
     &, 28, 31, 515.473,  2 ! 59Ni
     &, 30, 29, 500.346,  2 ! 59Zn
     &, 27, 33, 524.821,  1 ! 60Co
     &, 29, 31, 519.953,  1 ! 60Cu
     &, 26, 34, 525.394,  1 ! 60Fe
     &, 25, 35, 517.690,  1 ! 60Mn
     &, 28, 32, 526.862,  1 ! 60Ni
     &, 30, 30, 515.011,  1 ! 60Zn
     &, 27, 34, 534.143,  2/! 61Co
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 361, 370)/
     &  29, 32, 531.662,  2 ! 61Cu
     &, 26, 35, 531.039,  2 ! 61Fe
     &, 31, 30, 515.867,  2 ! 61Ga
     &, 28, 33, 534.683,  2 ! 61Ni
     &, 30, 31, 525.479,  2 ! 61Zn
     &, 27, 35, 540.822,  1 ! 62Co
     &, 29, 33, 540.549,  1 ! 62Cu
     &, 26, 36, 539.031,  1 ! 62Fe
     &, 31, 31, 527.958,  1 ! 62Ga
     &, 28, 34, 545.281,  1/! 62Ni
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 371, 380)/
     &  30, 32, 538.140,  1 ! 62Zn
     &, 27, 36, 549.240,  2 ! 63Co
     &, 29, 34, 551.403,  2 ! 63Cu
     &, 31, 32, 540.950,  2 ! 63Ga
     &, 32, 31, 531.027,  2 ! 63Ge
     &, 28, 35, 552.120,  2 ! 63Ni
     &, 30, 33, 547.253,  2 ! 63Zn
     &, 27, 37, 555.252,  1 ! 64Co
     &, 29, 35, 559.320,  1 ! 64Cu
     &, 31, 33, 551.168,  1/! 64Ga
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 381, 390)/
     &  32, 32, 545.979,  1 ! 64Ge
     &, 28, 36, 561.777,  1 ! 64Ni
     &, 30, 34, 559.115,  1 ! 64Zn
     &, 33, 32, 545.988,  2 ! 65As
     &, 29, 36, 569.230,  2 ! 65Cu
     &, 31, 34, 563.057,  2 ! 65Ga
     &, 32, 33, 556.031,  2 ! 65Ge
     &, 28, 37, 567.875,  2 ! 65Ni
     &, 30, 35, 567.095,  2 ! 65Zn
     &, 33, 33, 558.430,  1/! 66As
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 391, 400)/
     &  29, 37, 576.297,  1 ! 66Cu
     &, 31, 35, 572.198,  1 ! 66Ga
     &, 32, 34, 569.314,  1 ! 66Ge
     &, 28, 38, 576.843,  1 ! 66Ni
     &, 30, 36, 578.156,  1 ! 66Zn
     &, 33, 34, 571.632,  2 ! 67As
     &, 29, 38, 585.417,  2 ! 67Cu
     &, 31, 36, 583.425,  2 ! 67Ga
     &, 32, 35, 578.214,  2 ! 67Ge
     &, 28, 39, 582.364,  2/! 67Ni
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 401, 410)/
     &  34, 33, 561.279,  2 ! 67Se
     &, 30, 37, 585.209,  2 ! 67Zn
     &, 33, 35, 581.823,  1 ! 68As
     &, 29, 39, 591.573,  1 ! 68Cu
     &, 31, 37, 591.704,  1 ! 68Ga
     &, 32, 36, 590.808,  1 ! 68Ge
     &, 34, 34, 576.441,  1 ! 68Se
     &, 30, 38, 595.407,  1 ! 68Zn
     &, 33, 36, 594.245,  2 ! 69As
     &, 35, 34, 576.180,  2/! 69Br
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 411, 420)/
     &  29, 40, 600.195,  2 ! 69Cu
     &, 31, 38, 602.012,  2 ! 69Ga
     &, 32, 37, 599.004,  2 ! 69Ge
     &, 34, 35, 586.643,  2 ! 69Se
     &, 30, 39, 601.890,  2 ! 69Zn
     &, 33, 37, 603.536,  1 ! 70As
     &, 35, 35, 588.772,  1 ! 70Br
     &, 29, 41, 605.717,  1 ! 70Cu
     &, 31, 39, 609.667,  1 ! 70Ga
     &, 32, 38, 610.541,  1/! 70Ge
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 421, 430)/
     &  34, 36, 600.004,  1 ! 70Se
     &, 30, 40, 611.104,  1 ! 70Zn
     &, 33, 38, 615.161,  2 ! 71As
     &, 35, 36, 602.194,  2 ! 71Br
     &, 31, 40, 618.975,  2 ! 71Ga
     &, 32, 39, 617.957,  2 ! 71Ge
     &, 36, 35, 591.421,  2 ! 71Kr
     &, 34, 37, 609.576,  2 ! 71Se
     &, 30, 41, 616.940,  2 ! 71Zn
     &, 33, 39, 623.572,  1/! 72As
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 431, 440)/
     &  35, 37, 612.705,  1 ! 72Br
     &, 31, 41, 625.496,  1 ! 72Ga
     &, 32, 40, 628.705,  1 ! 72Ge
     &, 36, 36, 606.863,  1 ! 72Kr
     &, 34, 38, 622.452,  1 ! 72Se
     &, 30, 42, 625.822,  1 ! 72Zn
     &, 33, 40, 634.361,  2 ! 73As
     &, 35, 38, 625.517,  2 ! 73Br
     &, 31, 42, 634.707,  2 ! 73Ga
     &, 32, 41, 635.488,  2/! 73Ge
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 441, 450)/
     &  36, 37, 618.044,  2 ! 73Kr
     &, 34, 39, 630.838,  2 ! 73Se
     &, 30, 43, 630.789,  2 ! 73Zn
     &, 33, 41, 642.343,  1 ! 74As
     &, 35, 39, 635.214,  1 ! 74Br
     &, 31, 43, 641.068,  1 ! 74Ga
     &, 32, 42, 645.688,  1 ! 74Ge
     &, 36, 38, 631.156,  1 ! 74Kr
     &, 34, 40, 642.914,  1 ! 74Se
     &, 30, 44, 639.501,  1/! 74Zn
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 451, 460)/
     &  33, 42, 652.589,  2 ! 75As
     &, 35, 40, 647.149,  2 ! 75Br
     &, 31, 44, 649.680,  2 ! 75Ga
     &, 32, 43, 652.194,  2 ! 75Ge
     &, 36, 39, 641.368,  2 ! 75Kr
     &, 37, 38, 633.935,  2 ! 75Rb
     &, 34, 41, 650.942,  2 ! 75Se
     &, 30, 45, 644.363,  2 ! 75Zn
     &, 33, 43, 659.918,  1 ! 76As
     &, 35, 41, 656.365,  1/! 76Br
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 461, 470)/
     &  31, 45, 655.632,  1 ! 76Ga
     &, 32, 44, 661.623,  1 ! 76Ge
     &, 36, 40, 654.380,  1 ! 76Kr
     &, 37, 39, 645.107,  1 ! 76Rb
     &, 34, 42, 662.104,  1 ! 76Se
     &, 30, 46, 652.524,  1 ! 76Zn
     &, 33, 44, 669.614,  2 ! 77As
     &, 35, 42, 667.375,  2 ! 77Br
     &, 31, 46, 663.674,  2 ! 77Ga
     &, 32, 45, 667.695,  2/! 77Ge
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 471, 480)/
     &  36, 41, 663.587,  2 ! 77Kr
     &, 37, 40, 657.679,  2 ! 77Rb
     &, 34, 43, 669.522,  2 ! 77Se
     &, 38, 39, 649.746,  2 ! 77Sr
     &, 30, 47, 656.956,  2 ! 77Zn
     &, 33, 45, 676.510,  1 ! 78As
     &, 35, 43, 675.663,  1 ! 78Br
     &, 31, 47, 669.015,  1 ! 78Ga
     &, 32, 46, 676.313,  1 ! 78Ge
     &, 36, 42, 675.573,  1/! 78Kr
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 481, 490)/
     &  37, 41, 667.730,  1 ! 78Rb
     &, 34, 44, 680.019,  1 ! 78Se
     &, 38, 40, 663.708,  1 ! 78Sr
     &, 30, 48, 664.198,  1 ! 78Zn
     &, 33, 46, 685.562,  2 ! 79As
     &, 35, 44, 686.347,  2 ! 79Br
     &, 31, 48, 676.217,  2 ! 79Ga
     &, 32, 47, 682.194,  2 ! 79Ge
     &, 36, 43, 683.934,  2 ! 79Kr
     &, 37, 42, 679.572,  2/! 79Rb
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 491, 500)/
     &  34, 45, 686.980,  2 ! 79Se
     &, 38, 41, 673.390,  2 ! 79Sr
     &, 33, 47, 691.974,  1 ! 80As
     &, 35, 45, 694.240,  1 ! 80Br
     &, 31, 49, 681.009,  1 ! 80Ga
     &, 32, 48, 690.126,  1 ! 80Ge
     &, 36, 44, 695.463,  1 ! 80Kr
     &, 37, 43, 688.974,  1 ! 80Rb
     &, 34, 46, 696.893,  1 ! 80Se
     &, 38, 42, 686.391,  1/! 80Sr
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 501, 510)/
     &  33, 48, 700.625,  2 ! 81As
     &, 35, 46, 704.396,  2 ! 81Br
     &, 32, 49, 695.108,  2 ! 81Ge
     &, 36, 45, 703.345,  2 ! 81Kr
     &, 37, 44, 700.301,  2 ! 81Rb
     &, 34, 47, 703.594,  2 ! 81Se
     &, 38, 43, 695.533,  2 ! 81Sr
     &, 33, 49, 706.247,  1 ! 82As
     &, 35, 47, 711.990,  1 ! 82Br
     &, 32, 50, 702.630,  1/! 82Ge
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 511, 520)/
     &  36, 46, 714.301,  1 ! 82Kr
     &, 37, 45, 709.140,  1 ! 82Rb
     &, 34, 48, 712.861,  1 ! 82Se
     &, 38, 44, 708.144,  1 ! 82Sr
     &, 39, 43, 699.272,  1 ! 82Y
     &, 33, 50, 714.079,  2 ! 83As
     &, 35, 48, 721.589,  2 ! 83Br
     &, 36, 47, 721.766,  2 ! 83Kr
     &, 37, 46, 719.986,  2 ! 83Rb
     &, 34, 49, 718.756,  2/! 83Se
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 521, 530)/
     &  38, 45, 716.953,  2 ! 83Sr
     &, 39, 44, 711.874,  2 ! 83Y
     &, 33, 51, 718.360,  1 ! 84As
     &, 35, 49, 728.395,  1 ! 84Br
     &, 36, 48, 732.285,  1 ! 84Kr
     &, 37, 47, 728.823,  1 ! 84Rb
     &, 34, 50, 727.360,  1 ! 84Se
     &, 38, 46, 728.929,  1 ! 84Sr
     &, 39, 45, 721.198,  1 ! 84Y
     &, 40, 44, 718.163,  1/! 84Zr
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 531, 540)/
     &  35, 50, 737.377,  2 ! 85Br
     &, 36, 49, 739.397,  2 ! 85Kr
     &, 37, 48, 739.301,  2 ! 85Rb
     &, 34, 51, 732.060,  2 ! 85Se
     &, 38, 47, 737.455,  2 ! 85Sr
     &, 39, 46, 733.412,  2 ! 85Y
     &, 40, 45, 727.925,  2 ! 85Zr
     &, 35, 51, 742.739,  1 ! 86Br
     &, 36, 50, 749.259,  1 ! 86Kr
     &, 41, 45, 731.424,  1/! 86Nb
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 541, 550)/
     &  37, 49, 747.952,  1 ! 86Rb
     &, 34, 52, 738.421,  1 ! 86Se
     &, 38, 48, 748.944,  1 ! 86Sr
     &, 39, 47, 742.888,  1 ! 86Y
     &, 40, 46, 740.807,  1 ! 86Zr
     &, 35, 52, 749.061,  2 ! 87Br
     &, 36, 51, 754.775,  2 ! 87Kr
     &, 41, 46, 744.586,  2 ! 87Nb
     &, 37, 50, 757.881,  2 ! 87Rb
     &, 38, 49, 757.372,  2/! 87Sr
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 551, 560)/
     &  39, 48, 754.728,  2 ! 87Y
     &, 40, 47, 750.368,  2 ! 87Zr
     &, 35, 53, 754.012,  1 ! 88Br
     &, 36, 52, 761.829,  1 ! 88Kr
     &, 42, 46, 750.365,  1 ! 88Mo
     &, 41, 47, 754.647,  1 ! 88Nb
     &, 37, 51, 763.959,  1 ! 88Rb
     &, 38, 50, 768.485,  1 ! 88Sr
     &, 39, 49, 764.090,  1 ! 88Y
     &, 40, 48, 762.631,  1/! 88Zr
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 561, 570)/
     &  36, 53, 767.001,  2 ! 89Kr
     &, 42, 47, 760.737,  2 ! 89Mo
     &, 41, 48, 766.920,  2 ! 89Nb
     &, 37, 52, 771.146,  2 ! 89Rb
     &, 38, 51, 774.850,  2 ! 89Sr
     &, 39, 50, 775.559,  2 ! 89Y
     &, 40, 49, 771.941,  2 ! 89Zr
     &, 36, 54, 773.463,  1 ! 90Kr
     &, 42, 48, 773.755,  1 ! 90Mo
     &, 41, 49, 777.025,  1/! 90Nb
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 571, 580)/
     &  37, 53, 777.071,  1 ! 90Rb
     &, 38, 52, 782.653,  1 ! 90Sr
     &, 43, 47, 763.776,  1 ! 90Tc
     &, 39, 51, 782.416,  1 ! 90Y
     &, 40, 50, 783.918,  1 ! 90Zr
     &, 36, 55, 778.125,  2 ! 91Kr
     &, 42, 49, 783.859,  2 ! 91Mo
     &, 41, 50, 789.079,  2 ! 91Nb
     &, 37, 54, 783.542,  2 ! 91Rb
     &, 38, 53, 788.456,  2/! 91Sr
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 581, 590)/
     &  43, 48, 776.858,  2 ! 91Tc
     &, 39, 52, 790.357,  2 ! 91Y
     &, 40, 51, 791.117,  2 ! 91Zr
     &, 36, 56, 783.577,  1 ! 92Kr
     &, 42, 50, 796.539,  1 ! 92Mo
     &, 41, 51, 796.962,  1 ! 92Nb
     &, 37, 55, 788.764,  1 ! 92Rb
     &, 38, 54, 795.754,  1 ! 92Sr
     &, 43, 49, 787.885,  1 ! 92Tc
     &, 39, 53, 796.901,  1/! 92Y
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 591, 600)/
     &  40, 52, 799.753,  1 ! 92Zr
     &, 36, 57, 787.418,  2 ! 93Kr
     &, 42, 51, 804.606,  2 ! 93Mo
     &, 41, 52, 805.795,  2 ! 93Nb
     &, 37, 56, 794.636,  2 ! 93Rb
     &, 44, 49, 793.548,  2 ! 93Ru
     &, 38, 55, 801.213,  2 ! 93Sr
     &, 43, 50, 800.631,  2 ! 93Tc
     &, 39, 54, 804.378,  2 ! 93Y
     &, 40, 53, 806.485,  2/! 93Zr
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 601, 610)/
     &  42, 52, 814.287,  1 ! 94Mo
     &, 41, 53, 813.025,  1 ! 94Nb
     &, 37, 57, 799.247,  1 ! 94Rb
     &, 44, 50, 806.881,  1 ! 94Ru
     &, 38, 56, 807.965,  1 ! 94Sr
     &, 43, 51, 809.249,  1 ! 94Tc
     &, 39, 55, 810.605,  1 ! 94Y
     &, 40, 54, 814.704,  1 ! 94Zr
     &, 42, 53, 821.659,  2 ! 95Mo
     &, 41, 54, 821.516,  2/! 95Nb
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 611, 620)/
     &  37, 58, 804.409,  2 ! 95Rb
     &, 45, 50, 809.939,  2 ! 95Rh
     &, 44, 51, 815.834,  2 ! 95Ru
     &, 38, 57, 812.217,  2 ! 95Sr
     &, 43, 52, 819.177,  2 ! 95Tc
     &, 39, 56, 817.527,  2 ! 95Y
     &, 40, 55, 821.175,  2 ! 95Zr
     &, 42, 54, 830.813,  1 ! 96Mo
     &, 41, 55, 828.409,  1 ! 96Nb
     &, 37, 59, 808.701,  1/! 96Rb
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 621, 630)/
     &  45, 51, 819.304,  1 ! 96Rh
     &, 44, 52, 826.529,  1 ! 96Ru
     &, 38, 58, 818.218,  1 ! 96Sr
     &, 43, 53, 827.057,  1 ! 96Tc
     &, 39, 57, 822.796,  1 ! 96Y
     &, 40, 56, 829.028,  1 ! 96Zr
     &, 42, 55, 837.635,  2 ! 97Mo
     &, 41, 56, 836.484,  2 ! 97Nb
     &, 46, 51, 824.720,  2 ! 97Pd
     &, 45, 52, 830.303,  2/! 97Rh
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 631, 640)/
     &  44, 53, 834.595,  2 ! 97Ru
     &, 38, 59, 822.300,  2 ! 97Sr
     &, 43, 54, 836.532,  2 ! 97Tc
     &, 39, 58, 828.718,  2 ! 97Y
     &, 40, 57, 834.609,  2 ! 97Zr
     &, 42, 56, 846.277,  1 ! 98Mo
     &, 41, 57, 842.474,  1 ! 98Nb
     &, 46, 52, 836.302,  1 ! 98Pd
     &, 45, 53, 838.982,  1 ! 98Rh
     &, 44, 54, 844.823,  1/! 98Ru
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 641, 650)/
     &  38, 60, 828.672,  1 ! 98Sr
     &, 43, 55, 843.813,  1 ! 98Tc
     &, 39, 59, 833.699,  1 ! 98Y
     &, 40, 58, 841.019,  1 ! 98Zr
     &, 47, 52, 838.831,  2 ! 99Ag
     &, 42, 57, 852.203,  2 ! 99Mo
     &, 41, 58, 849.362,  2 ! 99Nb
     &, 46, 53, 845.216,  2 ! 99Pd
     &, 45, 54, 849.403,  2 ! 99Rh
     &, 44, 55, 852.288,  2/! 99Ru
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 651, 660)/
     &  43, 56, 852.777,  2 ! 99Tc
     &, 39, 60, 840.081,  2 ! 99Y
     &, 40, 59, 845.689,  2 ! 99Zr
     &, 47, 53, 848.323,  1 ! 100Ag
     &, 42, 58, 860.494,  1 ! 100Mo
     &, 41, 59, 855.048,  1 ! 100Nb
     &, 46, 54, 856.405,  1 ! 100Pd
     &, 45, 55, 857.550,  1 ! 100Rh
     &, 44, 56, 861.962,  1 ! 100Ru
     &, 43, 57, 859.542,  1/! 100Tc
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 661, 670)/
     &  40, 60, 852.470,  1 ! 100Zr
     &, 47, 54, 859.795,  2 ! 101Ag
     &, 48, 53, 853.212,  2 ! 101Cd
     &, 42, 59, 865.893,  2 ! 101Mo
     &, 41, 60, 862.109,  2 ! 101Nb
     &, 46, 55, 864.675,  2 ! 101Pd
     &, 45, 56, 867.440,  2 ! 101Rh
     &, 44, 57, 868.764,  2 ! 101Ru
     &, 43, 58, 867.921,  2 ! 101Tc
     &, 40, 61, 856.992,  2/! 101Zr
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 671, 680)/
     &  47, 55, 868.866,  1 ! 102Ag
     &, 48, 54, 865.184,  1 ! 102Cd
     &, 49, 53, 855.101,  1 ! 102In
     &, 42, 60, 874.011,  1 ! 102Mo
     &, 41, 61, 867.591,  1 ! 102Nb
     &, 46, 56, 875.244,  1 ! 102Pd
     &, 45, 57, 874.878,  1 ! 102Rh
     &, 44, 58, 877.984,  1 ! 102Ru
     &, 43, 59, 874.266,  1 ! 102Tc
     &, 47, 56, 879.408,  2/! 103Ag
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 681, 690)/
     &  48, 55, 874.426,  2 ! 103Cd
     &, 49, 54, 867.143,  2 ! 103In
     &, 42, 61, 879.130,  2 ! 103Mo
     &, 41, 62, 874.713,  2 ! 103Nb
     &, 46, 57, 882.868,  2 ! 103Pd
     &, 45, 58, 884.197,  2 ! 103Rh
     &, 44, 59, 884.217,  2 ! 103Ru
     &, 43, 60, 882.648,  2 ! 103Tc
     &, 47, 57, 887.830,  1 ! 104Ag
     &, 48, 56, 885.747,  1/! 104Cd
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 691, 700)/
     &  49, 55, 876.965,  1 ! 104In
     &, 42, 62, 887.092,  1 ! 104Mo
     &, 46, 58, 892.862,  1 ! 104Pd
     &, 45, 59, 891.197,  1 ! 104Rh
     &, 44, 60, 893.126,  1 ! 104Ru
     &, 50, 54, 871.482,  1 ! 104Sn
     &, 43, 61, 888.510,  1 ! 104Tc
     &, 47, 58, 897.826,  2 ! 105Ag
     &, 48, 57, 894.305,  2 ! 105Cd
     &, 49, 56, 888.526,  2/! 105In
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 701, 710)/
     &  42, 63, 891.804,  2 ! 105Mo
     &, 46, 59, 899.956,  2 ! 105Pd
     &, 45, 60, 900.171,  2 ! 105Rh
     &, 44, 61, 899.037,  2 ! 105Ru
     &, 50, 55, 881.494,  2 ! 105Sn
     &, 43, 62, 896.421,  2 ! 105Tc
     &, 47, 59, 905.752,  1 ! 106Ag
     &, 48, 58, 905.172,  1 ! 106Cd
     &, 49, 57, 897.844,  1 ! 106In
     &, 46, 60, 909.519,  1/! 106Pd
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 711, 720)/
     &  45, 61, 906.760,  1 ! 106Rh
     &, 44, 62, 907.503,  1 ! 106Ru
     &, 51, 55, 881.883,  1 ! 106Sb
     &, 50, 56, 893.466,  1 ! 106Sn
     &, 43, 63, 901.983,  1 ! 106Tc
     &, 47, 60, 915.299,  2 ! 107Ag
     &, 48, 59, 913.099,  2 ! 107Cd
     &, 49, 58, 908.830,  2 ! 107In
     &, 46, 61, 916.048,  2 ! 107Pd
     &, 45, 62, 915.320,  2/! 107Rh
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 721, 730)/
     &  44, 63, 912.952,  2 ! 107Ru
     &, 51, 56, 894.165,  2 ! 107Sb
     &, 50, 57, 902.947,  2 ! 107Sn
     &, 43, 64, 909.535,  2 ! 107Tc
     &, 47, 61, 922.568,  1 ! 108Ag
     &, 48, 60, 923.435,  1 ! 108Cd
     &, 49, 59, 917.501,  1 ! 108In
     &, 46, 62, 925.272,  1 ! 108Pd
     &, 45, 63, 921.551,  1 ! 108Rh
     &, 44, 64, 921.134,  1/! 108Ru
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 731, 740)/
     &  51, 57, 904.237,  1 ! 108Sb
     &, 50, 58, 914.519,  1 ! 108Sn
     &, 52, 56, 896.374,  1 ! 108Te
     &, 47, 62, 931.760,  2 ! 109Ag
     &, 48, 61, 930.796,  2 ! 109Cd
     &, 49, 60, 927.997,  2 ! 109In
     &, 46, 63, 931.427,  2 ! 109Pd
     &, 45, 64, 929.713,  2 ! 109Rh
     &, 44, 65, 926.196,  2 ! 109Ru
     &, 51, 58, 916.028,  2/! 109Sb
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 741, 750)/
     &  50, 59, 923.311,  2 ! 109Sn
     &, 52, 57, 906.596,  2 ! 109Te
     &, 47, 63, 938.566,  1 ! 110Ag
     &, 48, 62, 940.676,  1 ! 110Cd
     &, 49, 61, 935.954,  1 ! 110In
     &, 46, 64, 940.227,  1 ! 110Pd
     &, 45, 65, 935.605,  1 ! 110Rh
     &, 51, 59, 925.410,  1 ! 110Sb
     &, 50, 60, 934.596,  1 ! 110Sn
     &, 52, 58, 918.957,  1/! 110Te
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 751, 760)/
     &  47, 64, 947.408,  2 ! 111Ag
     &, 48, 63, 947.653,  2 ! 111Cd
     &, 49, 62, 946.022,  2 ! 111In
     &, 46, 65, 945.994,  2 ! 111Pd
     &, 45, 66, 943.276,  2 ! 111Rh
     &, 51, 60, 936.892,  2 ! 111Sb
     &, 50, 61, 942.775,  2 ! 111Sn
     &, 52, 59, 928.739,  2 ! 111Te
     &, 47, 65, 953.873,  1 ! 112Ag
     &, 48, 64, 957.049,  1/! 112Cd
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 761, 770)/
     &  53, 59, 929.998,  1 ! 112I
     &, 49, 63, 953.688,  1 ! 112In
     &, 46, 66, 954.362,  1 ! 112Pd
     &, 51, 61, 945.863,  1 ! 112Sb
     &, 50, 62, 953.564,  1 ! 112Sn
     &, 52, 60, 940.891,  1 ! 112Te
     &, 47, 66, 962.365,  2 ! 113Ag
     &, 48, 65, 963.593,  2 ! 113Cd
     &, 53, 60, 942.070,  2 ! 113I
     &, 49, 64, 963.132,  2/! 113In
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 771, 780)/
     &  46, 67, 959.747,  2 ! 113Pd
     &, 51, 62, 956.638,  2 ! 113Sb
     &, 50, 63, 961.309,  2 ! 113Sn
     &, 52, 61, 949.953,  2 ! 113Te
     &, 47, 67, 968.557,  1 ! 114Ag
     &, 48, 66, 972.634,  1 ! 114Cd
     &, 53, 61, 951.772,  1 ! 114I
     &, 49, 65, 970.408,  1 ! 114In
     &, 51, 63, 965.137,  1 ! 114Sb
     &, 50, 64, 971.609,  1/! 114Sn
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 781, 790)/
     &  52, 62, 961.674,  1 ! 114Te
     &, 54, 60, 945.009,  1 ! 114Xe
     &, 47, 68, 976.378,  2 ! 115Ag
     &, 48, 67, 978.779,  2 ! 115Cd
     &, 53, 62, 963.393,  2 ! 115I
     &, 49, 66, 979.444,  2 ! 115In
     &, 51, 64, 975.343,  2 ! 115Sb
     &, 50, 65, 979.156,  2 ! 115Sn
     &, 52, 63, 969.976,  2 ! 115Te
     &, 54, 61, 954.691,  2/! 115Xe
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 791, 800)/
     &  47, 69, 982.160,  1 ! 116Ag
     &, 48, 68, 987.475,  1 ! 116Cd
     &, 55, 61, 955.910,  1 ! 116Cs
     &, 53, 63, 972.455,  1 ! 116I
     &, 49, 67, 986.228,  1 ! 116In
     &, 51, 65, 983.340,  1 ! 116Sb
     &, 50, 66, 988.719,  1 ! 116Sn
     &, 52, 64, 980.998,  1 ! 116Te
     &, 54, 62, 967.333,  1 ! 116Xe
     &, 47, 70, 989.852,  2/! 117Ag
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 801, 810)/
     &  48, 69, 993.245,  2 ! 117Cd
     &, 55, 62, 968.202,  2 ! 117Cs
     &, 53, 64, 983.767,  2 ! 117I
     &, 49, 68, 994.991,  2 ! 117In
     &, 51, 66, 993.136,  2 ! 117Sb
     &, 50, 67, 995.663,  2 ! 117Sn
     &, 52, 65, 988.863,  2 ! 117Te
     &, 54, 63, 976.614,  2 ! 117Xe
     &, 48, 70,1001.608,  1 ! 118Cd
     &, 55, 63, 978.094,  1/! 118Cs
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 811, 820)/
     &  53, 65, 992.359,  1 ! 118I
     &, 49, 69,1001.568,  1 ! 118In
     &, 51, 67,1000.520,  1 ! 118Sb
     &, 50, 68,1004.990,  1 ! 118Sn
     &, 52, 66, 999.442,  1 ! 118Te
     &, 54, 64, 988.276,  1 ! 118Xe
     &, 56, 63, 981.243,  2 ! 119Ba
     &, 48, 71,1007.203,  2 ! 119Cd
     &, 55, 64, 990.025,  2 ! 119Cs
     &, 53, 66,1002.880,  2/! 119I
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 821, 830)/
     &  49, 70,1009.920,  2 ! 119In
     &, 51, 68,1010.108,  2 ! 119Sb
     &, 50, 69,1011.474,  2 ! 119Sn
     &, 52, 67,1007.032,  2 ! 119Te
     &, 54, 65, 997.108,  2 ! 119Xe
     &, 56, 64, 993.834,  1 ! 120Ba
     &, 48, 72,1015.025,  1 ! 120Cd
     &, 55, 65, 999.207,  1 ! 120Cs
     &, 53, 67,1011.132,  1 ! 120I
     &, 49, 71,1015.962,  1/! 120In
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 831, 840)/
     &  51, 69,1017.118,  1 ! 120Sb
     &, 50, 70,1020.581,  1 ! 120Sn
     &, 52, 68,1017.318,  1 ! 120Te
     &, 54, 66,1008.399,  1 ! 120Xe
     &, 56, 65,1003.426,  2 ! 121Ba
     &, 55, 66,1010.789,  2 ! 121Cs
     &, 53, 68,1021.344,  2 ! 121I
     &, 49, 72,1024.175,  2 ! 121In
     &, 51, 70,1026.357,  2 ! 121Sb
     &, 50, 71,1026.753,  2/! 121Sn
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 841, 850)/
     &  52, 69,1024.494,  2 ! 121Te
     &, 54, 67,1016.771,  2 ! 121Xe
     &, 56, 66,1015.188,  1 ! 122Ba
     &, 55, 67,1019.720,  1 ! 122Cs
     &, 53, 69,1029.435,  1 ! 122I
     &, 49, 73,1030.005,  1 ! 122In
     &, 51, 71,1033.164,  1 ! 122Sb
     &, 50, 72,1035.569,  1 ! 122Sn
     &, 52, 70,1034.362,  1 ! 122Te
     &, 54, 68,1027.653,  1/! 122Xe
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 851, 860)/
     &  56, 67,1024.390,  2 ! 123Ba
     &, 55, 68,1030.672,  2 ! 123Cs
     &, 53, 70,1039.317,  2 ! 123I
     &, 49, 74,1037.917,  2 ! 123In
     &, 51, 72,1042.129,  2 ! 123Sb
     &, 50, 73,1041.515,  2 ! 123Sn
     &, 52, 71,1041.295,  2 ! 123Te
     &, 54, 69,1035.854,  2 ! 123Xe
     &, 56, 68,1035.821,  1 ! 124Ba
     &, 55, 69,1039.384,  1/! 124Cs
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 861, 870)/
     &  53, 71,1046.780,  1 ! 124I
     &, 49, 75,1043.649,  1 ! 124In
     &, 51, 73,1048.597,  1 ! 124Sb
     &, 50, 74,1050.006,  1 ! 124Sn
     &, 52, 72,1050.719,  1 ! 124Te
     &, 54, 70,1046.086,  1 ! 124Xe
     &, 56, 69,1044.603,  2 ! 125Ba
     &, 55, 70,1049.965,  2 ! 125Cs
     &, 53, 72,1056.331,  2 ! 125I
     &, 49, 76,1051.120,  2/! 125In
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 871, 880)/
     &  51, 74,1057.307,  2 ! 125Sb
     &, 50, 75,1055.740,  2 ! 125Sn
     &, 52, 73,1057.292,  2 ! 125Te
     &, 54, 71,1053.818,  2 ! 125Xe
     &, 56, 70,1055.775,  1 ! 126Ba
     &, 55, 71,1058.327,  1 ! 126Cs
     &, 53, 73,1063.473,  1 ! 126I
     &, 49, 77,1056.592,  1 ! 126In
     &, 51, 75,1063.529,  1 ! 126Sb
     &, 50, 76,1063.933,  1/! 126Sn
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 881, 890)/
     &  52, 74,1066.411,  1 ! 126Te
     &, 54, 72,1063.942,  1 ! 126Xe
     &, 56, 71,1064.046,  2 ! 127Ba
     &, 55, 72,1068.275,  2 ! 127Cs
     &, 53, 74,1072.614,  2 ! 127I
     &, 49, 78,1063.934,  2 ! 127In
     &, 57, 70,1058.264,  2 ! 127La
     &, 51, 76,1071.903,  2 ! 127Sb
     &, 50, 77,1069.581,  2 ! 127Sn
     &, 52, 75,1072.701,  2/! 127Te
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 891, 900)/
     &  54, 73,1071.167,  2 ! 127Xe
     &, 56, 72,1074.840,  1 ! 128Ba
     &, 55, 73,1076.075,  1 ! 128Cs
     &, 53, 75,1079.439,  1 ! 128I
     &, 49, 79,1069.175,  1 ! 128In
     &, 57, 71,1067.256,  1 ! 128La
     &, 51, 77,1078.000,  1 ! 128Sb
     &, 50, 78,1077.493,  1 ! 128Sn
     &, 52, 76,1081.480,  1 ! 128Te
     &, 54, 74,1080.784,  1/! 128Xe
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 901, 910)/
     &  56, 73,1082.546,  2 ! 129Ba
     &, 55, 74,1085.775,  2 ! 129Cs
     &, 53, 76,1088.282,  2 ! 129I
     &, 49, 80,1076.027,  2 ! 129In
     &, 57, 72,1077.767,  2 ! 129La
     &, 51, 78,1085.972,  2 ! 129Sb
     &, 50, 79,1082.765,  2 ! 129Sn
     &, 52, 77,1087.567,  2 ! 129Te
     &, 54, 75,1087.692,  2 ! 129Xe
     &, 56, 74,1092.804,  1/! 130Ba
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 911, 920)/
     &  55, 75,1093.147,  1 ! 130Cs
     &, 53, 77,1094.746,  1 ! 130I
     &, 49, 81,1081.059,  1 ! 130In
     &, 57, 73,1086.319,  1 ! 130La
     &, 51, 79,1091.794,  1 ! 130Sb
     &, 50, 80,1090.576,  1 ! 130Sn
     &, 52, 78,1095.979,  1 ! 130Te
     &, 54, 76,1096.947,  1 ! 130Xe
     &, 56, 75,1100.299,  2 ! 131Ba
     &, 58, 73,1091.478,  2/! 131Ce
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 921, 930)/
     &  55, 76,1102.422,  2 ! 131Cs
     &, 53, 78,1103.371,  2 ! 131I
     &, 57, 74,1096.561,  2 ! 131La
     &, 51, 80,1099.585,  2 ! 131Sb
     &, 50, 81,1095.748,  2 ! 131Sn
     &, 52, 79,1101.904,  2 ! 131Te
     &, 54, 77,1103.559,  2 ! 131Xe
     &, 56, 76,1110.098,  1 ! 132Ba
     &, 58, 74,1102.420,  1 ! 132Ce
     &, 55, 77,1109.602,  1/! 132Cs
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 931, 940)/
     &  53, 79,1109.698,  1 ! 132I
     &, 57, 75,1104.602,  1 ! 132La
     &, 51, 81,1105.167,  1 ! 132Sb
     &, 50, 82,1102.730,  1 ! 132Sn
     &, 52, 80,1109.988,  1 ! 132Te
     &, 54, 78,1112.496,  1 ! 132Xe
     &, 56, 77,1117.285,  2 ! 133Ba
     &, 58, 75,1110.322,  2 ! 133Ce
     &, 55, 78,1118.588,  2 ! 133Cs
     &, 53, 80,1117.966,  2/! 133I
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 941, 950)/
     &  57, 76,1114.504,  2 ! 133La
     &, 59, 74,1105.339,  2 ! 133Pr
     &, 51, 82,1112.609,  2 ! 133Sb
     &, 52, 81,1115.776,  2 ! 133Te
     &, 54, 79,1118.943,  2 ! 133Xe
     &, 56, 78,1126.756,  1 ! 134Ba
     &, 58, 76,1120.993,  1 ! 134Ce
     &, 55, 79,1125.480,  1 ! 134Cs
     &, 53, 81,1124.106,  1 ! 134I
     &, 57, 77,1122.274,  1/! 134La
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 951, 960)/
     &  59, 75,1113.911,  1 ! 134Pr
     &, 51, 83,1115.570,  1 ! 134Sb
     &, 52, 82,1123.588,  1 ! 134Te
     &, 54, 80,1127.478,  1 ! 134Xe
     &, 56, 79,1133.730,  2 ! 135Ba
     &, 58, 77,1128.845,  2 ! 135Ce
     &, 55, 80,1134.307,  2 ! 135Cs
     &, 53, 82,1132.003,  2 ! 135I
     &, 57, 78,1131.747,  2 ! 135La
     &, 60, 75,1119.020,  2/! 135Nd
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 961, 970)/
     &  59, 76,1124.502,  2 ! 135Pr
     &, 52, 83,1126.590,  2 ! 135Te
     &, 54, 81,1133.931,  2 ! 135Xe
     &, 56, 80,1142.838,  1 ! 136Ba
     &, 58, 78,1138.867,  1 ! 136Ce
     &, 55, 81,1141.072,  1 ! 136Cs
     &, 53, 83,1135.709,  1 ! 136I
     &, 57, 79,1139.189,  1 ! 136La
     &, 60, 76,1129.992,  1 ! 136Nd
     &, 59, 77,1132.984,  1/! 136Pr
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 971, 980)/
     &  52, 84,1131.891,  1 ! 136Te
     &, 54, 82,1141.921,  1 ! 136Xe
     &, 56, 81,1149.736,  2 ! 137Ba
     &, 58, 79,1146.348,  2 ! 137Ce
     &, 55, 82,1149.346,  2 ! 137Cs
     &, 53, 84,1141.071,  2 ! 137I
     &, 57, 80,1148.351,  2 ! 137La
     &, 60, 77,1138.283,  2 ! 137Nd
     &, 61, 76,1132.301,  2 ! 137Pm
     &, 59, 78,1142.866,  2/! 137Pr
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 981, 990)/
     &  54, 83,1145.783,  2 ! 137Xe
     &, 56, 82,1158.348,  1 ! 138Ba
     &, 58, 80,1156.075,  1 ! 138Ce
     &, 55, 83,1153.627,  1 ! 138Cs
     &, 53, 85,1144.152,  1 ! 138I
     &, 57, 81,1155.816,  1 ! 138La
     &, 60, 78,1148.975,  1 ! 138Nd
     &, 61, 77,1141.193,  1 ! 138Pm
     &, 59, 79,1150.855,  1 ! 138Pr
     &, 54, 84,1151.670,  1/! 138Xe
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i= 991,1000)/
     &  56, 83,1163.072,  2 ! 139Ba
     &, 58, 81,1163.548,  2 ! 139Ce
     &, 55, 84,1159.559,  2 ! 139Cs
     &, 57, 82,1164.595,  2 ! 139La
     &, 60, 79,1157.067,  2 ! 139Nd
     &, 61, 78,1151.734,  2 ! 139Pm
     &, 59, 80,1160.653,  2 ! 139Pr
     &, 62, 77,1145.752,  2 ! 139Sm
     &, 54, 85,1155.462,  2 ! 139Xe
     &, 56, 84,1169.503,  1/! 140Ba
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1001,1010)/
     &  58, 82,1172.734,  1 ! 140Ce
     &, 55, 85,1164.241,  1 ! 140Cs
     &, 57, 83,1169.756,  1 ! 140La
     &, 60, 80,1167.308,  1 ! 140Nd
     &, 61, 79,1160.486,  1 ! 140Pm
     &, 59, 81,1168.564,  1 ! 140Pr
     &, 62, 78,1157.003,  1 ! 140Sm
     &, 54, 86,1160.963,  1 ! 140Xe
     &, 56, 85,1174.270,  2 ! 141Ba
     &, 58, 83,1178.163,  2/! 141Ce
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1011,1020)/
     &  55, 86,1170.072,  2 ! 141Cs
     &, 63, 78,1158.693,  2 ! 141Eu
     &, 57, 84,1176.515,  2 ! 141La
     &, 60, 81,1175.363,  2 ! 141Nd
     &, 61, 80,1170.848,  2 ! 141Pm
     &, 59, 82,1177.961,  2 ! 141Pr
     &, 62, 79,1165.505,  2 ! 141Sm
     &, 54, 87,1164.855,  2 ! 141Xe
     &, 56, 86,1180.182,  1 ! 142Ba
     &, 58, 84,1185.332,  1/! 142Ce
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1021,1030)/
     &  55, 87,1174.094,  1 ! 142Cs
     &, 63, 79,1168.364,  1 ! 142Eu
     &, 57, 85,1181.597,  1 ! 142La
     &, 60, 82,1185.181,  1 ! 142Nd
     &, 61, 81,1179.509,  1 ! 142Pm
     &, 59, 83,1183.804,  1 ! 142Pr
     &, 62, 80,1176.645,  1 ! 142Sm
     &, 54, 88,1169.977,  1 ! 142Xe
     &, 56, 87,1184.443,  2 ! 143Ba
     &, 58, 85,1190.478,  2/! 143Ce
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1031,1040)/
     &  55, 88,1179.576,  2 ! 143Cs
     &, 63, 80,1179.366,  2 ! 143Eu
     &, 64, 79,1172.684,  2 ! 143Gd
     &, 57, 86,1187.961,  2 ! 143La
     &, 60, 83,1191.303,  2 ! 143Nd
     &, 61, 82,1189.480,  2 ! 143Pm
     &, 59, 84,1191.151,  2 ! 143Pr
     &, 62, 81,1185.250,  2 ! 143Sm
     &, 56, 88,1190.535,  1 ! 144Ba
     &, 58, 86,1197.371,  1/! 144Ce
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1041,1050)/
     &  55, 89,1183.218,  1 ! 144Cs
     &, 63, 81,1188.664,  1 ! 144Eu
     &, 64, 80,1184.185,  1 ! 144Gd
     &, 57, 87,1192.653,  1 ! 144La
     &, 60, 84,1199.121,  1 ! 144Nd
     &, 61, 83,1196.009,  1 ! 144Pm
     &, 59, 85,1196.908,  1 ! 144Pr
     &, 62, 82,1195.774,  1 ! 144Sm
     &, 56, 89,1194.397,  2 ! 145Ba
     &, 58, 87,1202.132,  2/! 145Ce
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1051,1060)/
     &  55, 90,1189.079,  2 ! 145Cs
     &, 63, 82,1199.035,  2 ! 145Eu
     &, 64, 81,1193.257,  2 ! 145Gd
     &, 57, 88,1198.714,  2 ! 145La
     &, 60, 85,1204.877,  2 ! 145Nd
     &, 61, 84,1203.934,  2 ! 145Pm
     &, 59, 86,1203.854,  2 ! 145Pr
     &, 62, 83,1202.538,  2 ! 145Sm
     &, 65, 80,1185.775,  2 ! 145Tb
     &, 56, 90,1200.208,  1/! 146Ba
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1061,1070)/
     &  58, 88,1208.843,  1 ! 146Ce
     &, 63, 83,1206.282,  1 ! 146Eu
     &, 64, 82,1204.299,  1 ! 146Gd
     &, 57, 89,1203.326,  1 ! 146La
     &, 60, 86,1212.442,  1 ! 146Nd
     &, 61, 85,1210.178,  1 ! 146Pm
     &, 59, 87,1209.141,  1 ! 146Pr
     &, 62, 84,1210.938,  1 ! 146Sm
     &, 65, 81,1195.416,  1 ! 146Tb
     &, 58, 89,1213.395,  2/! 147Ce
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1071,1080)/
     &  66, 81,1199.105,  2 ! 147Dy
     &, 63, 84,1214.778,  2 ! 147Eu
     &, 64, 83,1211.667,  2 ! 147Gd
     &, 57, 90,1209.478,  2 ! 147La
     &, 60, 87,1217.734,  2 ! 147Nd
     &, 61, 86,1217.848,  2 ! 147Pm
     &, 59, 88,1215.813,  2 ! 147Pr
     &, 62, 85,1217.290,  2 ! 147Sm
     &, 65, 82,1206.188,  2 ! 147Tb
     &, 58, 90,1219.937,  1/! 148Ce
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1081,1090)/
     &  66, 82,1210.737,  1 ! 148Dy
     &, 63, 85,1221.550,  1 ! 148Eu
     &, 64, 84,1220.800,  1 ! 148Gd
     &, 60, 88,1225.069,  1 ! 148Nd
     &, 61, 87,1223.749,  1 ! 148Pm
     &, 59, 89,1220.954,  1 ! 148Pr
     &, 62, 86,1225.432,  1 ! 148Sm
     &, 65, 83,1214.390,  1 ! 148Tb
     &, 58, 91,1224.769,  2 ! 149Ce
     &, 66, 83,1218.569,  2/! 149Dy
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1091,1100)/
     &  63, 86,1229.825,  2 ! 149Eu
     &, 64, 85,1227.735,  2 ! 149Gd
     &, 67, 82,1211.786,  2 ! 149Ho
     &, 60, 89,1230.108,  2 ! 149Nd
     &, 61, 88,1231.014,  2 ! 149Pm
     &, 59, 90,1227.886,  2 ! 149Pr
     &, 62, 87,1231.304,  2 ! 149Sm
     &, 65, 84,1223.255,  2 ! 149Tb
     &, 66, 84,1228.250,  1 ! 150Dy
     &, 68, 82,1215.486,  1/! 150Er
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1101,1110)/
     &  63, 87,1236.214,  1 ! 150Eu
     &, 64, 86,1236.440,  1 ! 150Gd
     &, 67, 83,1220.368,  1 ! 150Ho
     &, 60, 90,1237.487,  1 ! 150Nd
     &, 61, 89,1236.573,  1 ! 150Pm
     &, 59, 91,1233.268,  1 ! 150Pr
     &, 62, 88,1239.289,  1 ! 150Sm
     &, 65, 85,1230.991,  1 ! 150Tb
     &, 66, 85,1235.783,  2 ! 151Dy
     &, 68, 83,1223.817,  2/! 151Er
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1111,1120)/
     &  63, 88,1244.180,  2 ! 151Eu
     &, 64, 87,1242.915,  2 ! 151Gd
     &, 67, 84,1229.840,  2 ! 151Ho
     &, 60, 91,1242.822,  2 ! 151Nd
     &, 61, 90,1244.481,  2 ! 151Pm
     &, 59, 92,1240.099,  2 ! 151Pr
     &, 62, 89,1244.886,  2 ! 151Sm
     &, 65, 86,1239.573,  2 ! 151Tb
     &, 69, 82,1215.635,  2 ! 151Tm
     &, 66, 86,1245.370,  1/! 152Dy
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1121,1130)/
     &  68, 84,1234.099,  1 ! 152Er
     &, 63, 89,1250.485,  1 ! 152Eu
     &, 64, 88,1251.522,  1 ! 152Gd
     &, 67, 85,1238.181,  1 ! 152Ho
     &, 60, 92,1250.095,  1 ! 152Nd
     &, 61, 91,1250.456,  1 ! 152Pm
     &, 62, 90,1253.145,  1 ! 152Sm
     &, 65, 87,1246.889,  1 ! 152Tb
     &, 69, 83,1224.716,  1 ! 152Tm
     &, 66, 87,1252.481,  2/! 153Dy
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1131,1140)/
     &  68, 85,1242.071,  2 ! 153Er
     &, 63, 90,1259.036,  2 ! 153Eu
     &, 64, 89,1258.009,  2 ! 153Gd
     &, 67, 86,1247.497,  2 ! 153Ho
     &, 60, 93,1255.380,  2 ! 153Nd
     &, 61, 92,1257.998,  2 ! 153Pm
     &, 62, 91,1259.012,  2 ! 153Sm
     &, 65, 88,1255.437,  2 ! 153Tb
     &, 69, 84,1234.848,  2 ! 153Tm
     &, 70, 83,1227.406,  2/! 153Yb
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1141,1150)/
     &  66, 88,1261.789,  1 ! 154Dy
     &, 68, 86,1252.272,  1 ! 154Er
     &, 63, 91,1265.471,  1 ! 154Eu
     &, 64, 90,1266.666,  1 ! 154Gd
     &, 67, 87,1255.250,  1 ! 154Ho
     &, 61, 93,1263.760,  1 ! 154Pm
     &, 62, 92,1266.981,  1 ! 154Sm
     &, 65, 89,1262.420,  1 ! 154Tb
     &, 69, 85,1243.580,  1 ! 154Tm
     &, 70, 84,1238.317,  1/! 154Yb
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1151,1160)/
     &  66, 89,1268.626,  2 ! 155Dy
     &, 68, 87,1259.961,  2 ! 155Er
     &, 63, 92,1273.641,  2 ! 155Eu
     &, 64, 91,1273.105,  2 ! 155Gd
     &, 67, 88,1264.741,  2 ! 155Ho
     &, 71, 84,1238.157,  2 ! 155Lu
     &, 61, 94,1270.481,  2 ! 155Pm
     &, 62, 93,1272.795,  2 ! 155Sm
     &, 65, 90,1271.507,  2 ! 155Tb
     &, 69, 86,1253.572,  2/! 155Tm
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1161,1170)/
     &  70, 85,1246.789,  2 ! 155Yb
     &, 66, 90,1278.068,  1 ! 156Dy
     &, 68, 88,1269.906,  1 ! 156Er
     &, 63, 93,1279.971,  1 ! 156Eu
     &, 64, 92,1281.642,  1 ! 156Gd
     &, 67, 89,1272.188,  1 ! 156Ho
     &, 71, 85,1247.438,  1 ! 156Lu
     &, 62, 94,1280.039,  1 ! 156Sm
     &, 65, 91,1278.421,  1 ! 156Tb
     &, 69, 87,1262.133,  1/! 156Tm
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1171,1180)/
     &  70, 86,1257.471,  1 ! 156Yb
     &, 66, 91,1285.037,  2 ! 157Dy
     &, 68, 89,1277.137,  2 ! 157Er
     &, 63, 94,1287.425,  2 ! 157Eu
     &, 64, 93,1288.002,  2 ! 157Gd
     &, 72, 85,1249.877,  2 ! 157Hf
     &, 67, 90,1281.720,  2 ! 157Ho
     &, 71, 86,1258.170,  2 ! 157Lu
     &, 62, 95,1285.602,  2 ! 157Sm
     &, 65, 92,1287.162,  2/! 157Tb
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1181,1190)/
     &  69, 88,1271.755,  2 ! 157Tm
     &, 70, 87,1265.752,  2 ! 157Yb
     &, 66, 92,1294.094,  1 ! 158Dy
     &, 68, 90,1287.149,  1 ! 158Er
     &, 63, 95,1293.271,  1 ! 158Eu
     &, 64, 94,1295.940,  1 ! 158Gd
     &, 72, 86,1261.209,  1 ! 158Hf
     &, 67, 91,1289.335,  1 ! 158Ho
     &, 71, 87,1267.002,  1 ! 158Lu
     &, 65, 93,1293.941,  1/! 158Tb
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1191,1200)/
     &  69, 89,1279.767,  1 ! 158Tm
     &, 70, 88,1276.084,  1 ! 158Yb
     &, 66, 93,1300.927,  2 ! 159Dy
     &, 68, 91,1294.581,  2 ! 159Er
     &, 63, 96,1300.033,  2 ! 159Eu
     &, 64, 95,1301.883,  2 ! 159Gd
     &, 72, 87,1269.861,  2 ! 159Hf
     &, 67, 92,1298.291,  2 ! 159Ho
     &, 71, 88,1277.333,  2 ! 159Lu
     &, 65, 94,1302.074,  2/! 159Tb
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1201,1210)/
     &  69, 90,1289.598,  2 ! 159Tm
     &, 70, 89,1283.916,  2 ! 159Yb
     &, 66, 94,1309.501,  1 ! 160Dy
     &, 68, 92,1304.314,  1 ! 160Er
     &, 63, 97,1305.715,  1 ! 160Eu
     &, 64, 96,1309.335,  1 ! 160Gd
     &, 72, 88,1280.883,  1 ! 160Hf
     &, 67, 93,1305.433,  1 ! 160Ho
     &, 71, 89,1285.845,  1 ! 160Lu
     &, 73, 87,1270.130,  1/! 160Ta
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1211,1220)/
     &  65, 95,1308.450,  1 ! 160Tb
     &, 69, 91,1297.610,  1 ! 160Tm
     &, 70, 90,1294.247,  1 ! 160Yb
     &, 66, 95,1315.955,  2 ! 161Dy
     &, 68, 93,1311.531,  2 ! 161Er
     &, 64, 97,1314.971,  2 ! 161Gd
     &, 72, 89,1289.334,  2 ! 161Hf
     &, 67, 94,1314.320,  2 ! 161Ho
     &, 71, 90,1296.067,  2 ! 161Lu
     &, 73, 88,1281.262,  2/! 161Ta
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1221,1230)/
     &  65, 96,1316.148,  2 ! 161Tb
     &, 69, 92,1307.232,  2 ! 161Tm
     &, 70, 91,1302.169,  2 ! 161Yb
     &, 66, 96,1324.152,  1 ! 162Dy
     &, 68, 94,1320.741,  1 ! 162Er
     &, 64, 98,1321.896,  1 ! 162Gd
     &, 72, 90,1300.036,  1 ! 162Hf
     &, 67, 95,1321.235,  1 ! 162Ho
     &, 71, 91,1304.398,  1 ! 162Lu
     &, 73, 89,1290.203,  1/! 162Ta
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1231,1240)/
     &  65, 97,1322.513,  1 ! 162Tb
     &, 69, 93,1315.163,  1 ! 162Tm
     &, 74, 88,1283.841,  1 ! 162W
     &, 70, 92,1312.181,  1 ! 162Yb
     &, 66, 97,1330.424,  2 ! 163Dy
     &, 68, 95,1327.646,  2 ! 163Er
     &, 72, 91,1308.118,  2 ! 163Hf
     &, 67, 96,1329.639,  2 ! 163Ho
     &, 71, 92,1314.500,  2 ! 163Lu
     &, 73, 90,1300.935,  2/! 163Ta
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1241,1250)/
     &  65, 98,1329.505,  2 ! 163Tb
     &, 69, 94,1324.465,  2 ! 163Tm
     &, 74, 89,1292.633,  2 ! 163W
     &, 70, 93,1320.083,  2 ! 163Yb
     &, 66, 98,1338.081,  1 ! 164Dy
     &, 68, 96,1336.489,  1 ! 164Er
     &, 72, 92,1318.699,  1 ! 164Hf
     &, 67, 97,1336.269,  1 ! 164Ho
     &, 71, 93,1322.782,  1 ! 164Lu
     &, 73, 91,1309.717,  1/! 164Ta
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1251,1260)/
     &  65, 99,1335.007,  1 ! 164Tb
     &, 69, 95,1331.745,  1 ! 164Tm
     &, 74, 90,1303.894,  1 ! 164W
     &, 70, 94,1329.864,  1 ! 164Yb
     &, 66, 99,1343.797,  2 ! 165Dy
     &, 68, 97,1343.139,  2 ! 165Er
     &, 72, 93,1326.751,  2 ! 165Hf
     &, 67, 98,1344.299,  2 ! 165Ho
     &, 71, 94,1332.433,  2 ! 165Lu
     &, 73, 92,1320.069,  2/! 165Ta
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1261,1270)/
     &  69, 96,1340.762,  2 ! 165Tm
     &, 74, 91,1312.596,  2 ! 165W
     &, 70, 95,1337.217,  2 ! 165Yb
     &, 66,100,1350.841,  1 ! 166Dy
     &, 68, 98,1351.614,  1 ! 166Er
     &, 72, 94,1337.043,  1 ! 166Hf
     &, 67, 99,1350.542,  1 ! 166Ho
     &, 71, 95,1340.445,  1 ! 166Lu
     &, 73, 93,1328.880,  1 ! 166Ta
     &, 69, 97,1347.784,  1/! 166Tm
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1271,1280)/
     &  74, 92,1323.478,  1 ! 166W
     &, 70, 96,1346.710,  1 ! 166Yb
     &, 68, 99,1358.050,  2 ! 167Er
     &, 72, 95,1344.784,  2 ! 167Hf
     &, 67,100,1357.863,  2 ! 167Ho
     &, 71, 96,1349.867,  2 ! 167Lu
     &, 75, 92,1323.937,  2 ! 167Re
     &, 73, 94,1338.802,  2 ! 167Ta
     &, 69, 98,1356.519,  2 ! 167Tm
     &, 74, 93,1332.019,  2/! 167W
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1281,1290)/
     &  70, 97,1353.782,  2 ! 167Yb
     &, 68,100,1365.821,  1 ! 168Er
     &, 72, 96,1354.806,  1 ! 168Hf
     &, 67,101,1363.888,  1 ! 168Ho
     &, 71, 97,1357.589,  1 ! 168Lu
     &, 75, 93,1333.059,  1 ! 168Re
     &, 73, 95,1347.324,  1 ! 168Ta
     &, 69, 99,1363.359,  1 ! 168Tm
     &, 74, 94,1342.641,  1 ! 168W
     &, 70, 98,1362.836,  1/! 168Yb
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1291,1300)/
     &  68,101,1371.825,  2 ! 169Er
     &, 72, 97,1362.308,  2 ! 169Hf
     &, 67,102,1370.483,  2 ! 169Ho
     &, 71, 98,1366.441,  2 ! 169Lu
     &, 76, 93,1335.198,  2 ! 169Os
     &, 75, 94,1343.560,  2 ! 169Re
     &, 73, 96,1357.025,  2 ! 169Ta
     &, 69,100,1371.394,  2 ! 169Tm
     &, 74, 95,1351.103,  2 ! 169W
     &, 70, 99,1369.704,  2/! 169Yb
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1301,1310)/
     &  68,102,1379.083,  1 ! 170Er
     &, 72, 98,1371.969,  1 ! 170Hf
     &, 67,103,1375.862,  1 ! 170Ho
     &, 71, 99,1373.951,  1 ! 170Lu
     &, 76, 94,1346.250,  1 ! 170Os
     &, 75, 95,1352.422,  1 ! 170Re
     &, 73, 97,1365.187,  1 ! 170Ta
     &, 69,101,1377.988,  1 ! 170Tm
     &, 74, 96,1361.205,  1 ! 170W
     &, 70,100,1378.173,  1/! 170Yb
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1311,1320)/
     &  68,103,1384.765,  2 ! 171Er
     &, 72, 99,1379.221,  2 ! 171Hf
     &, 77, 94,1346.189,  2 ! 171Ir
     &, 71,100,1382.525,  2 ! 171Lu
     &, 76, 95,1354.951,  2 ! 171Os
     &, 75, 96,1362.674,  2 ! 171Re
     &, 73, 98,1374.739,  2 ! 171Ta
     &, 69,102,1385.474,  2 ! 171Tm
     &, 74, 97,1369.256,  2 ! 171W
     &, 70,101,1384.788,  2/! 171Yb
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1321,1330)/
     &  68,104,1391.614,  1 ! 172Er
     &, 72,100,1388.323,  1 ! 172Hf
     &, 77, 95,1355.400,  1 ! 172Ir
     &, 71,101,1389.501,  1 ! 172Lu
     &, 76, 96,1365.703,  1 ! 172Os
     &, 75, 97,1371.155,  1 ! 172Re
     &, 73, 99,1382.620,  1 ! 172Ta
     &, 69,103,1391.720,  1 ! 172Tm
     &, 74, 98,1379.238,  1 ! 172W
     &, 70,102,1392.808,  1/! 172Yb
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1331,1340)/
     &  68,105,1396.924,  2 ! 173Er
     &, 72,101,1395.335,  2 ! 173Hf
     &, 77, 96,1366.062,  2 ! 173Ir
     &, 71,102,1397.718,  2 ! 173Lu
     &, 76, 97,1374.345,  2 ! 173Os
     &, 78, 95,1357.160,  2 ! 173Pt
     &, 75, 98,1381.087,  2 ! 173Re
     &, 73,100,1391.652,  2 ! 173Ta
     &, 69,104,1398.638,  2 ! 173Tm
     &, 74, 99,1386.970,  2/! 173W
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1341,1350)/
     &  70,103,1399.175,  2 ! 173Yb
     &, 72,102,1403.966,  1 ! 174Hf
     &, 77, 97,1375.114,  1 ! 174Ir
     &, 71,103,1404.481,  1 ! 174Lu
     &, 76, 98,1384.626,  1 ! 174Os
     &, 78, 96,1368.371,  1 ! 174Pt
     &, 75, 99,1389.369,  1 ! 174Re
     &, 73,101,1399.334,  1 ! 174Ta
     &, 69,105,1404.334,  1 ! 174Tm
     &, 74,100,1396.651,  1/! 174W
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1351,1360)/
     &  70,104,1406.641,  1 ! 174Yb
     &, 79, 96,1367.891,  2 ! 175Au
     &, 72,103,1410.756,  2 ! 175Hf
     &, 77, 98,1385.456,  2 ! 175Ir
     &, 71,104,1412.149,  2 ! 175Lu
     &, 76, 99,1392.788,  2 ! 175Os
     &, 78, 97,1377.153,  2 ! 175Pt
     &, 75,100,1399.011,  2 ! 175Re
     &, 73,102,1407.775,  2 ! 175Ta
     &, 69,106,1410.845,  2/! 175Tm
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1361,1370)/
     &  74,101,1404.093,  2 ! 175W
     &, 70,105,1412.464,  2 ! 175Yb
     &, 79, 97,1377.202,  1 ! 176Au
     &, 72,104,1418.847,  1 ! 176Hf
     &, 77, 99,1394.207,  1 ! 176Ir
     &, 71,105,1418.443,  1 ! 176Lu
     &, 76,100,1402.960,  1 ! 176Os
     &, 78, 98,1388.125,  1 ! 176Pt
     &, 75,101,1406.902,  1 ! 176Re
     &, 73,103,1414.967,  1/! 176Ta
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1371,1380)/
     &  69,107,1416.217,  1 ! 176Tm
     &, 74,102,1413.285,  1 ! 176W
     &, 70,106,1419.335,  1 ! 176Yb
     &, 79, 98,1388.064,  2 ! 177Au
     &, 72,105,1425.230,  2 ! 177Hf
     &, 80, 97,1378.742,  2 ! 177Hg
     &, 77,100,1404.259,  2 ! 177Ir
     &, 71,106,1425.516,  2 ! 177Lu
     &, 76,101,1410.841,  2 ! 177Os
     &, 78, 99,1397.006,  2/! 177Pt
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1381,1390)/
     &  75,102,1416.124,  2 ! 177Re
     &, 73,104,1423.290,  2 ! 177Ta
     &, 74,103,1420.506,  2 ! 177W
     &, 70,107,1424.902,  2 ! 177Yb
     &, 79, 99,1397.356,  1 ! 178Au
     &, 72,106,1432.857,  1 ! 178Hf
     &, 80, 98,1390.093,  1 ! 178Hg
     &, 77,101,1412.781,  1 ! 178Ir
     &, 71,107,1431.505,  1 ! 178Lu
     &, 76,102,1420.643,  1/! 178Os
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1391,1400)/
     &  78,100,1407.358,  1 ! 178Pt
     &, 75,103,1423.846,  1 ! 178Re
     &, 73,105,1430.161,  1 ! 178Ta
     &, 74,104,1429.288,  1 ! 178W
     &, 70,108,1431.648,  1 ! 178Yb
     &, 79,100,1407.767,  2 ! 179Au
     &, 72,107,1438.957,  2 ! 179Hf
     &, 80, 99,1399.035,  2 ! 179Hg
     &, 77,102,1422.472,  2 ! 179Ir
     &, 71,108,1438.387,  2/! 179Lu
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1401,1410)/
     &  76,103,1428.255,  2 ! 179Os
     &, 78,101,1415.810,  2 ! 179Pt
     &, 75,104,1432.737,  2 ! 179Re
     &, 73,106,1438.059,  2 ! 179Ta
     &, 74,105,1436.213,  2 ! 179W
     &, 79,101,1416.719,  1 ! 180Au
     &, 72,108,1446.345,  1 ! 180Hf
     &, 80,100,1410.167,  1 ! 180Hg
     &, 77,103,1430.584,  1 ! 180Ir
     &, 71,109,1444.029,  1/! 180Lu
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1411,1420)/
     &  76,104,1437.656,  1 ! 180Os
     &, 78,102,1425.992,  1 ! 180Pt
     &, 75,105,1440.048,  1 ! 180Re
     &, 73,107,1444.698,  1 ! 180Ta
     &, 74,106,1444.625,  1 ! 180W
     &, 79,102,1426.801,  2 ! 181Au
     &, 72,109,1452.041,  2 ! 181Hf
     &, 80,101,1419.168,  2 ! 181Hg
     &, 77,104,1440.066,  2 ! 181Ir
     &, 76,105,1444.918,  2/! 181Os
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1421,1430)/
     &  78,103,1434.003,  2 ! 181Pt
     &, 75,106,1448.731,  2 ! 181Re
     &, 73,108,1452.281,  2 ! 181Ta
     &, 74,107,1451.310,  2 ! 181W
     &, 79,103,1435.412,  1 ! 182Au
     &, 72,110,1458.700,  1 ! 182Hf
     &, 80,102,1429.660,  1 ! 182Hg
     &, 77,105,1447.777,  1 ! 182Ir
     &, 76,106,1454.160,  1 ! 182Os
     &, 78,104,1443.995,  1/! 182Pt
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1431,1440)/
     &  75,107,1455.792,  1 ! 182Re
     &, 73,109,1458.344,  1 ! 182Ta
     &, 74,108,1459.373,  1 ! 182W
     &, 79,104,1445.314,  2 ! 183Au
     &, 72,111,1464.050,  2 ! 183Hf
     &, 80,103,1438.212,  2 ! 183Hg
     &, 77,106,1456.959,  2 ! 183Ir
     &, 76,107,1461.142,  2 ! 183Os
     &, 78,105,1451.717,  2 ! 183Pt
     &, 75,108,1464.225,  2/! 183Re
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1441,1450)/
     &  73,110,1465.278,  2 ! 183Ta
     &, 81,102,1429.569,  2 ! 183Tl
     &, 74,109,1465.563,  2 ! 183W
     &, 79,105,1453.596,  1 ! 184Au
     &, 72,112,1470.333,  1 ! 184Hf
     &, 80,104,1448.633,  1 ! 184Hg
     &, 77,107,1464.451,  1 ! 184Ir
     &, 76,108,1469.956,  1 ! 184Os
     &, 78,106,1461.368,  1 ! 184Pt
     &, 75,109,1470.697,  1/! 184Re
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1451,1460)/
     &  73,111,1470.892,  1 ! 184Ta
     &, 81,103,1438.711,  1 ! 184Tl
     &, 74,110,1472.975,  1 ! 184W
     &, 79,106,1463.178,  2 ! 185Au
     &, 80,105,1456.805,  2 ! 185Hg
     &, 77,108,1473.302,  2 ! 185Ir
     &, 76,109,1476.582,  2 ! 185Os
     &, 82,103,1440.840,  2 ! 185Pb
     &, 78,107,1468.720,  2 ! 185Pt
     &, 75,110,1478.379,  2/! 185Re
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1461,1470)/
     &  73,112,1477.502,  2 ! 185Ta
     &, 81,104,1448.993,  2 ! 185Tl
     &, 74,111,1478.730,  2 ! 185W
     &, 79,107,1471.209,  1 ! 186Au
     &, 80,106,1467.087,  1 ! 186Hg
     &, 77,109,1480.240,  1 ! 186Ir
     &, 76,110,1484.854,  1 ! 186Os
     &, 82,104,1451.502,  1 ! 186Pb
     &, 78,108,1478.132,  1 ! 186Pt
     &, 75,111,1484.559,  1/! 186Re
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1471,1480)/
     &  73,113,1482.814,  1 ! 186Ta
     &, 81,105,1457.814,  1 ! 186Tl
     &, 74,112,1485.930,  1 ! 186W
     &, 79,108,1480.461,  2 ! 187Au
     &, 80,107,1474.868,  2 ! 187Hg
     &, 77,110,1488.866,  2 ! 187Ir
     &, 76,111,1491.146,  2 ! 187Os
     &, 82,105,1460.183,  2 ! 187Pb
     &, 78,109,1485.183,  2 ! 187Pt
     &, 75,112,1491.926,  2/! 187Re
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1481,1490)/
     &  81,106,1467.956,  2 ! 187Tl
     &, 74,113,1491.396,  2 ! 187W
     &, 79,109,1488.153,  1 ! 188Au
     &, 80,108,1484.760,  1 ! 188Hg
     &, 77,111,1495.551,  1 ! 188Ir
     &, 76,112,1499.135,  1 ! 188Os
     &, 82,106,1470.815,  1 ! 188Pb
     &, 78,110,1494.233,  1 ! 188Pt
     &, 75,113,1497.798,  1 ! 188Re
     &, 81,107,1476.388,  1/! 188Tl
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1491,1500)/
     &  74,114,1498.232,  1 ! 188W
     &, 79,110,1497.144,  2 ! 189Au
     &, 83,106,1470.474,  2 ! 189Bi
     &, 80,109,1492.162,  2 ! 189Hg
     &, 77,112,1503.779,  2 ! 189Ir
     &, 76,113,1505.060,  2 ! 189Os
     &, 82,107,1479.247,  2 ! 189Pb
     &, 78,111,1501.087,  2 ! 189Pt
     &, 75,114,1504.834,  2 ! 189Re
     &, 81,108,1486.189,  2/! 189Tl
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1501,1510)/
     &  74,115,1503.117,  2 ! 189W
     &, 79,111,1504.682,  1 ! 190Au
     &, 83,107,1479.526,  1 ! 190Bi
     &, 80,110,1501.984,  1 ! 190Hg
     &, 77,113,1510.071,  1 ! 190Ir
     &, 76,114,1512.852,  1 ! 190Os
     &, 82,108,1489.679,  1 ! 190Pb
     &, 78,112,1509.906,  1 ! 190Pt
     &, 75,115,1510.456,  1 ! 190Re
     &, 81,109,1494.401,  1/! 190Tl
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1511,1520)/
     &  74,116,1509.938,  1 ! 190W
     &, 79,112,1513.748,  2 ! 191Au
     &, 83,108,1489.798,  2 ! 191Bi
     &, 80,111,1509.575,  2 ! 191Hg
     &, 77,114,1518.141,  2 ! 191Ir
     &, 76,115,1518.613,  2 ! 191Os
     &, 82,109,1497.760,  2 ! 191Pb
     &, 78,113,1516.358,  2 ! 191Pt
     &, 75,116,1517.351,  2 ! 191Re
     &, 81,110,1503.983,  2/! 191Tl
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1521,1530)/
     &  79,113,1520.717,  1 ! 192Au
     &, 83,109,1498.489,  1 ! 192Bi
     &, 80,112,1519.137,  1 ! 192Hg
     &, 77,115,1524.340,  1 ! 192Ir
     &, 76,116,1526.172,  1 ! 192Os
     &, 82,110,1507.892,  1 ! 192Pb
     &, 78,114,1525.015,  1 ! 192Pt
     &, 81,111,1511.974,  1 ! 192Tl
     &, 79,114,1529.381,  2 ! 193Au
     &, 83,110,1508.451,  2/! 193Bi
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1531,1540)/
     &  80,113,1526.259,  2 ! 193Hg
     &, 77,116,1532.105,  2 ! 193Ir
     &, 76,117,1531.755,  2 ! 193Os
     &, 82,111,1515.744,  2 ! 193Pb
     &, 84,109,1500.419,  2 ! 193Po
     &, 78,115,1531.262,  2 ! 193Pt
     &, 81,112,1521.476,  2 ! 193Tl
     &, 79,115,1536.349,  1 ! 194Au
     &, 83,111,1516.943,  1 ! 194Bi
     &, 80,114,1535.516,  1/! 194Hg
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1541,1550)/
     &  77,117,1538.172,  1 ! 194Ir
     &, 76,118,1538.857,  1 ! 194Os
     &, 82,112,1525.555,  1 ! 194Pb
     &, 84,110,1510.990,  1 ! 194Po
     &, 78,116,1539.640,  1 ! 194Pt
     &, 81,113,1529.338,  1 ! 194Tl
     &, 79,116,1544.736,  2 ! 195Au
     &, 83,112,1526.715,  2 ! 195Bi
     &, 80,115,1542.432,  2 ! 195Hg
     &, 77,118,1545.421,  2/! 195Ir
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1551,1560)/
     &  76,119,1544.202,  2 ! 195Os
     &, 82,113,1533.367,  2 ! 195Pb
     &, 84,111,1519.312,  2 ! 195Po
     &, 78,117,1545.749,  2 ! 195Pt
     &, 81,114,1538.449,  2 ! 195Tl
     &, 85,111,1519.591,  1 ! 196At
     &, 79,117,1551.398,  1 ! 196Au
     &, 83,113,1534.866,  1 ! 196Bi
     &, 80,116,1551.300,  1 ! 196Hg
     &, 77,119,1551.241,  1/! 196Ir
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1561,1570)/
     &  82,114,1543.039,  1 ! 196Pb
     &, 84,112,1529.534,  1 ! 196Po
     &, 78,118,1553.671,  1 ! 196Pt
     &, 81,115,1546.021,  1 ! 196Tl
     &, 85,112,1529.643,  2 ! 197At
     &, 79,118,1559.458,  2 ! 197Au
     &, 83,114,1544.588,  2 ! 197Bi
     &, 80,117,1558.260,  2 ! 197Hg
     &, 77,120,1558.303,  2 ! 197Ir
     &, 82,115,1550.590,  2/! 197Pb
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1571,1580)/
     &  84,113,1537.625,  2 ! 197Po
     &, 78,119,1559.521,  2 ! 197Pt
     &, 81,116,1555.073,  2 ! 197Tl
     &, 85,113,1538.355,  1 ! 198At
     &, 79,119,1565.970,  1 ! 198Au
     &, 83,115,1552.550,  1 ! 198Bi
     &, 80,118,1566.561,  1 ! 198Hg
     &, 77,121,1563.464,  1 ! 198Ir
     &, 82,116,1559.932,  1 ! 198Pb
     &, 84,114,1547.537,  1/! 198Po
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1581,1590)/
     &  78,120,1567.083,  1 ! 198Pt
     &, 81,117,1562.315,  1 ! 198Tl
     &, 85,114,1548.226,  2 ! 199At
     &, 79,120,1573.555,  2 ! 199Au
     &, 83,116,1561.931,  2 ! 199Bi
     &, 80,119,1573.226,  2 ! 199Hg
     &, 82,117,1567.384,  2 ! 199Pb
     &, 84,115,1555.589,  2 ! 199Po
     &, 78,121,1572.654,  2 ! 199Pt
     &, 81,118,1570.966,  2/! 199Tl
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1591,1600)/
     &  85,115,1556.498,  1 ! 200At
     &, 79,121,1579.823,  1 ! 200Au
     &, 83,117,1569.853,  1 ! 200Bi
     &, 80,120,1581.254,  1 ! 200Hg
     &, 82,118,1576.335,  1 ! 200Pb
     &, 84,116,1565.351,  1 ! 200Po
     &, 78,122,1579.905,  1 ! 200Pt
     &, 86,114,1550.786,  1 ! 200Rn
     &, 81,119,1578.018,  1 ! 200Tl
     &, 85,116,1566.420,  2/! 201At
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1601,1610)/
     &  79,122,1586.995,  2 ! 201Au
     &, 83,118,1578.875,  2 ! 201Bi
     &, 80,121,1587.484,  2 ! 201Hg
     &, 82,119,1583.574,  2 ! 201Pb
     &, 84,117,1573.092,  2 ! 201Po
     &, 78,123,1585.117,  2 ! 201Pt
     &, 86,115,1559.067,  2 ! 201Rn
     &, 81,120,1586.215,  2 ! 201Tl
     &, 85,117,1574.491,  1 ! 202At
     &, 79,123,1592.526,  1/! 202Au
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1611,1620)/
     &  83,119,1586.576,  1 ! 202Bi
     &, 87,115,1559.246,  1 ! 202Fr
     &, 80,122,1595.240,  1 ! 202Hg
     &, 82,120,1592.261,  1 ! 202Pb
     &, 84,118,1582.534,  1 ! 202Po
     &, 86,116,1569.069,  1 ! 202Rn
     &, 81,121,1593.089,  1 ! 202Tl
     &, 85,118,1584.013,  2 ! 203At
     &, 79,124,1599.718,  2 ! 203Au
     &, 83,120,1595.208,  2/! 203Bi
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1621,1630)/
     &  87,116,1569.248,  2 ! 203Fr
     &, 80,123,1601.232,  2 ! 203Hg
     &, 82,121,1599.185,  2 ! 203Pb
     &, 84,119,1590.186,  2 ! 203Po
     &, 86,117,1577.261,  2 ! 203Rn
     &, 81,122,1600.942,  2 ! 203Tl
     &, 85,119,1592.085,  1 ! 204At
     &, 79,125,1605.010,  1 ! 204Au
     &, 83,121,1602.500,  1 ! 204Bi
     &, 87,117,1577.680,  1/! 204Fr
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1631,1640)/
     &  80,124,1608.730,  1 ! 204Hg
     &, 82,122,1607.579,  1 ! 204Pb
     &, 84,120,1599.147,  1 ! 204Po
     &, 88,116,1571.487,  1 ! 204Ra
     &, 86,118,1587.102,  1 ! 204Rn
     &, 81,123,1607.598,  1 ! 204Tl
     &, 85,120,1601.147,  2 ! 205At
     &, 83,122,1610.821,  2 ! 205Bi
     &, 87,118,1587.662,  2 ! 205Fr
     &, 80,125,1614.398,  2/! 205Hg
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1641,1650)/
     &  82,123,1614.311,  2 ! 205Pb
     &, 84,121,1606.545,  2 ! 205Po
     &, 88,117,1579.859,  2 ! 205Ra
     &, 86,119,1595.004,  2 ! 205Rn
     &, 81,124,1615.153,  2 ! 205Tl
     &, 85,121,1608.988,  1 ! 206At
     &, 83,123,1617.856,  1 ! 206Bi
     &, 87,119,1595.873,  1 ! 206Fr
     &, 80,126,1621.126,  1 ! 206Hg
     &, 82,124,1622.401,  1/! 206Pb
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1651,1660)/
     &  84,122,1615.231,  1 ! 206Po
     &, 88,118,1589.951,  1 ! 206Ra
     &, 86,120,1604.446,  1 ! 206Rn
     &, 81,125,1621.657,  1 ! 206Tl
     &, 85,122,1617.640,  2 ! 207At
     &, 83,124,1625.953,  2 ! 207Bi
     &, 87,120,1605.415,  2 ! 207Fr
     &, 82,125,1629.140,  2 ! 207Pb
     &, 84,123,1622.262,  2 ! 207Po
     &, 88,119,1598.282,  2/! 207Ra
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1661,1670)/
     &  86,121,1612.237,  2 ! 207Rn
     &, 81,126,1628.501,  2 ! 207Tl
     &, 85,123,1625.042,  1 ! 208At
     &, 83,125,1632.846,  1 ! 208Bi
     &, 87,121,1613.607,  1 ! 208Fr
     &, 82,126,1636.508,  1 ! 208Pb
     &, 84,124,1630.659,  1 ! 208Po
     &, 88,120,1608.124,  1 ! 208Ra
     &, 86,122,1621.179,  1 ! 208Rn
     &, 81,127,1632.299,  1/! 208Tl
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1671,1680)/
     &  89,120,1608.223,  2 ! 209Ac
     &, 85,124,1633.361,  2 ! 209At
     &, 83,126,1640.306,  2 ! 209Bi
     &, 87,122,1622.668,  2 ! 209Fr
     &, 82,127,1640.445,  2 ! 209Pb
     &, 84,125,1637.629,  2 ! 209Po
     &, 88,121,1616.156,  2 ! 209Ra
     &, 86,123,1628.685,  2 ! 209Rn
     &, 81,128,1637.253,  2 ! 209Tl
     &, 89,121,1616.555,  1/! 210Ac
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1681,1690)/
     &  85,125,1640.521,  1 ! 210At
     &, 83,127,1644.911,  1 ! 210Bi
     &, 87,123,1630.620,  1 ! 210Fr
     &, 82,128,1645.630,  1 ! 210Pb
     &, 84,126,1645.290,  1 ! 210Po
     &, 88,122,1625.588,  1 ! 210Ra
     &, 86,124,1637.370,  1 ! 210Rn
     &, 81,129,1640.926,  1 ! 210Tl
     &, 89,122,1626.087,  2 ! 211Ac
     &, 85,126,1648.270,  2/! 211At
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1691,1700)/
     &  83,128,1650.047,  2 ! 211Bi
     &, 87,124,1639.272,  2 ! 211Fr
     &, 82,129,1649.456,  2 ! 211Pb
     &, 84,127,1649.843,  2 ! 211Po
     &, 88,123,1633.489,  2 ! 211Ra
     &, 86,125,1644.595,  2 ! 211Rn
     &, 89,123,1634.378,  1 ! 212Ac
     &, 85,127,1653.313,  1 ! 212At
     &, 83,129,1654.388,  1 ! 212Bi
     &, 87,125,1646.813,  1/! 212Fr
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1701,1710)/
     &  82,130,1654.598,  1 ! 212Pb
     &, 84,128,1655.852,  1 ! 212Po
     &, 88,124,1642.451,  1 ! 212Ra
     &, 86,126,1652.572,  1 ! 212Rn
     &, 89,124,1643.460,  2 ! 213Ac
     &, 85,128,1659.349,  2 ! 213At
     &, 83,130,1659.568,  2 ! 213Bi
     &, 87,126,1654.751,  2 ! 213Fr
     &, 82,131,1658.247,  2 ! 213Pb
     &, 84,129,1660.205,  2/! 213Po
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1711,1720)/
     &  88,125,1650.123,  2 ! 213Ra
     &, 86,127,1657.684,  2 ! 213Rn
     &, 90,123,1636.608,  2 ! 213Th
     &, 89,125,1651.562,  1 ! 214Ac
     &, 85,129,1664.221,  1 ! 214At
     &, 83,131,1663.606,  1 ! 214Bi
     &, 87,127,1660.232,  1 ! 214Fr
     &, 82,132,1663.364,  1 ! 214Pb
     &, 84,130,1666.093,  1 ! 214Po
     &, 88,126,1658.394,  1/! 214Ra
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1721,1730)/
     &  86,128,1664.377,  1 ! 214Rn
     &, 90,124,1646.049,  1 ! 214Th
     &, 89,126,1659.824,  2 ! 215Ac
     &, 85,130,1670.165,  2 ! 215At
     &, 83,132,1668.758,  2 ! 215Bi
     &, 87,128,1667.029,  2 ! 215Fr
     &, 84,131,1670.226,  2 ! 215Po
     &, 88,127,1664.025,  2 ! 215Ra
     &, 86,129,1669.300,  2 ! 215Rn
     &, 90,125,1654.121,  2/! 215Th
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1731,1740)/
     &  89,127,1665.865,  1 ! 216Ac
     &, 85,131,1674.738,  1 ! 216At
     &, 83,133,1672.570,  1 ! 216Bi
     &, 87,129,1672.435,  1 ! 216Fr
     &, 84,132,1675.989,  1 ! 216Po
     &, 88,128,1671.343,  1 ! 216Ra
     &, 86,130,1675.948,  1 ! 216Rn
     &, 90,126,1662.673,  1 ! 216Th
     &, 89,128,1673.216,  2 ! 217Ac
     &, 85,132,1680.665,  2/! 217At
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1741,1750)/
     &  87,130,1679.175,  2 ! 217Fr
     &, 84,133,1679.869,  2 ! 217Po
     &, 88,129,1676.818,  2 ! 217Ra
     &, 86,131,1680.615,  2 ! 217Rn
     &, 90,127,1668.993,  2 ! 217Th
     &, 89,129,1679.152,  1 ! 218Ac
     &, 85,133,1685.019,  1 ! 218At
     &, 87,131,1684.504,  1 ! 218Fr
     &, 84,134,1685.546,  1 ! 218Po
     &, 88,130,1684.127,  1/! 218Ra
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1751,1760)/
     &  86,132,1687.124,  1 ! 218Rn
     &, 90,128,1676.844,  1 ! 218Th
     &, 89,130,1686.500,  2 ! 219Ac
     &, 85,134,1690.660,  2 ! 219At
     &, 87,132,1691.008,  2 ! 219Fr
     &, 88,131,1689.466,  2 ! 219Ra
     &, 86,133,1691.577,  2 ! 219Rn
     &, 90,129,1682.808,  2 ! 219Th
     &, 89,131,1692.385,  1 ! 220Ac
     &, 85,135,1695.062,  1/! 220At
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1761,1770)/
     &  87,133,1696.227,  1 ! 220Fr
     &, 88,132,1696.651,  1 ! 220Ra
     &, 86,134,1697.880,  1 ! 220Rn
     &, 90,130,1690.686,  1 ! 220Th
     &, 89,132,1699.686,  2 ! 221Ac
     &, 87,134,1702.504,  2 ! 221Fr
     &, 88,133,1702.029,  2 ! 221Ra
     &, 86,135,1702.171,  2 ! 221Rn
     &, 90,131,1696.487,  2 ! 221Th
     &, 89,133,1705.658,  1/! 222Ac
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1771,1780)/
     &  87,135,1707.502,  1 ! 222Fr
     &, 91,131,1698.751,  1 ! 222Pa
     &, 88,134,1708.746,  1 ! 222Ra
     &, 86,136,1708.253,  1 ! 222Rn
     &, 90,132,1704.296,  1 ! 222Th
     &, 89,134,1712.522,  2 ! 223Ac
     &, 87,136,1713.530,  2 ! 223Fr
     &, 91,132,1706.452,  2 ! 223Pa
     &, 88,135,1713.895,  2 ! 223Ra
     &, 90,133,1710.309,  2/! 223Th
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1781,1790)/
     &  89,135,1718.200,  1 ! 224Ac
     &, 87,137,1718.274,  1 ! 224Fr
     &, 91,133,1713.056,  1 ! 224Pa
     &, 88,136,1720.388,  1 ! 224Ra
     &, 90,134,1717.643,  1 ! 224Th
     &, 89,136,1724.864,  2 ! 225Ac
     &, 87,138,1724.265,  2 ! 225Fr
     &, 91,134,1720.605,  2 ! 225Pa
     &, 88,137,1725.286,  2 ! 225Ra
     &, 90,135,1723.405,  2/! 225Th
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1791,1800)/
     &  89,137,1730.261,  1 ! 226Ac
     &, 87,139,1728.667,  1 ! 226Fr
     &, 91,135,1726.968,  1 ! 226Pa
     &, 88,138,1731.679,  1 ! 226Ra
     &, 90,136,1730.591,  1 ! 226Th
     &, 92,134,1725.029,  1 ! 226U
     &, 89,138,1736.784,  2 ! 227Ac
     &, 87,140,1734.619,  2 ! 227Fr
     &, 91,136,1734.237,  2 ! 227Pa
     &, 88,139,1736.231,  2/! 227Ra
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1801,1810)/
     &  90,137,1736.045,  2 ! 227Th
     &, 92,135,1731.406,  2 ! 227U
     &, 89,139,1741.810,  1 ! 228Ac
     &, 91,137,1740.271,  1 ! 228Pa
     &, 88,140,1742.547,  1 ! 228Ra
     &, 90,138,1743.165,  1 ! 228Th
     &, 92,136,1739.137,  1 ! 228U
     &, 89,140,1748.057,  2 ! 229Ac
     &, 93,136,1741.889,  2 ! 229Np
     &, 91,138,1747.325,  2/! 229Pa
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1811,1820)/
     &  88,141,1746.840,  2 ! 229Ra
     &, 90,139,1748.414,  2 ! 229Th
     &, 92,137,1745.229,  2 ! 229U
     &, 89,141,1753.089,  1 ! 230Ac
     &, 93,137,1748.487,  1 ! 230Np
     &, 91,139,1753.118,  1 ! 230Pa
     &, 88,142,1753.071,  1 ! 230Ra
     &, 90,140,1755.205,  1 ! 230Th
     &, 92,138,1752.894,  1 ! 230U
     &, 89,142,1759.011,  2/! 231Ac
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1821,1830)/
     &  93,138,1756.165,  2 ! 231Np
     &, 91,140,1759.933,  2 ! 231Pa
     &, 90,141,1760.326,  2 ! 231Th
     &, 92,139,1758.793,  2 ! 231U
     &, 89,143,1763.842,  1 ! 232Ac
     &, 93,139,1762.572,  1 ! 232Np
     &, 91,141,1765.493,  1 ! 232Pa
     &, 94,138,1760.718,  1 ! 232Pu
     &, 90,142,1766.763,  1 ! 232Th
     &, 92,140,1766.048,  1/! 232U
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1831,1840)/
     &  93,140,1769.924,  2 ! 233Np
     &, 91,142,1772.012,  2 ! 233Pa
     &, 94,139,1767.110,  2 ! 233Pu
     &, 90,143,1771.549,  2 ! 233Th
     &, 92,141,1771.802,  2 ! 233U
     &, 95,139,1769.981,  1 ! 234Am
     &, 93,141,1776.055,  1 ! 234Np
     &, 91,143,1777.222,  1 ! 234Pa
     &, 94,140,1774.881,  1 ! 234Pu
     &, 90,144,1777.741,  1/! 234Th
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1841,1850)/
     &  92,142,1778.646,  1 ! 234U
     &, 95,140,1777.862,  2 ! 235Am
     &, 93,142,1783.038,  2 ! 235Np
     &, 91,144,1783.322,  2 ! 235Pa
     &, 94,141,1781.135,  2 ! 235Pu
     &, 90,145,1782.275,  2 ! 235Th
     &, 92,143,1783.944,  2 ! 235U
     &, 95,141,1784.564,  1 ! 236Am
     &, 96,140,1781.912,  1 ! 236Cm
     &, 93,143,1788.723,  1/! 236Np
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1851,1860)/
     &  91,145,1788.174,  1 ! 236Pa
     &, 94,142,1788.478,  1 ! 236Pu
     &, 92,144,1790.490,  1 ! 236U
     &, 95,142,1792.016,  2 ! 237Am
     &, 96,141,1788.703,  2 ! 237Cm
     &, 93,144,1795.352,  2 ! 237Np
     &, 91,146,1794.146,  2 ! 237Pa
     &, 94,143,1794.351,  2 ! 237Pu
     &, 92,145,1795.615,  2 ! 237U
     &, 95,143,1798.311,  1/! 238Am
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1861,1870)/
     &  97,141,1790.883,  1 ! 238Bk
     &, 96,142,1796.547,  1 ! 238Cm
     &, 93,145,1800.840,  1 ! 238Np
     &, 91,147,1798.587,  1 ! 238Pa
     &, 94,144,1801.349,  1 ! 238Pu
     &, 92,146,1801.768,  1 ! 238U
     &, 95,144,1805.410,  2 ! 239Am
     &, 97,142,1798.954,  2 ! 239Bk
     &, 96,143,1802.927,  2 ! 239Cm
     &, 93,146,1807.058,  2/! 239Np
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1871,1880)/
     &  94,145,1806.997,  2 ! 239Pu
     &, 92,147,1806.574,  2 ! 239U
     &, 95,145,1811.428,  1 ! 240Am
     &, 97,143,1805.596,  1 ! 240Bk
     &, 98,142,1802.494,  1 ! 240Cf
     &, 96,144,1810.376,  1 ! 240Cm
     &, 93,147,1812.226,  1 ! 240Np
     &, 94,146,1813.531,  1 ! 240Pu
     &, 92,148,1812.506,  1 ! 240U
     &, 95,146,1818.010,  2/! 241Am
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1881,1890)/
     &  97,144,1813.278,  2 ! 241Bk
     &, 98,143,1809.405,  2 ! 241Cf
     &, 96,145,1816.464,  2 ! 241Cm
     &, 93,148,1818.198,  2 ! 241Np
     &, 94,147,1818.772,  2 ! 241Pu
     &, 95,147,1823.552,  1 ! 242Am
     &, 97,145,1819.649,  1 ! 242Bk
     &, 98,144,1817.335,  1 ! 242Cf
     &, 96,146,1823.430,  1 ! 242Cm
     &, 93,149,1823.059,  1/! 242Np
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1891,1900)/
     &  94,148,1825.082,  1 ! 242Pu
     &, 95,148,1829.916,  2 ! 243Am
     &, 97,146,1826.836,  2 ! 243Bk
     &, 98,145,1823.829,  2 ! 243Cf
     &, 96,147,1829.126,  2 ! 243Cm
     &, 99,144,1819.156,  2 ! 243Es
     &, 94,149,1830.116,  2 ! 243Pu
     &, 95,149,1835.279,  1 ! 244Am
     &, 97,147,1832.947,  1 ! 244Bk
     &, 98,146,1831.345,  1/! 244Cf
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1901,1910)/
     &  96,148,1835.926,  1 ! 244Cm
     &, 99,145,1826.058,  1 ! 244Es
     &, 94,150,1836.137,  1 ! 244Pu
     &, 95,150,1841.332,  2 ! 245Am
     &, 97,148,1839.853,  2 ! 245Bk
     &, 98,147,1837.505,  2 ! 245Cf
     &, 96,149,1841.446,  2 ! 245Cm
     &, 99,146,1833.720,  2 ! 245Es
     &,100,145,1829.297,  2 ! 245Fm
     &, 94,151,1840.855,  2/! 245Pu
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1911,1920)/
     &  95,151,1846.381,  1 ! 246Am
     &, 97,149,1845.716,  1 ! 246Bk
     &, 98,148,1844.857,  1 ! 246Cf
     &, 96,150,1847.903,  1 ! 246Cm
     &, 99,147,1840.241,  1 ! 246Es
     &,100,146,1837.258,  1 ! 246Fm
     &, 94,152,1846.794,  1 ! 246Pu
     &, 95,152,1852.243,  2 ! 247Am
     &, 97,150,1852.324,  2 ! 247Bk
     &, 98,149,1850.875,  2/! 247Cf
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1921,1930)/
     &  96,151,1853.060,  2 ! 247Cm
     &, 99,148,1847.693,  2 ! 247Es
     &,100,147,1843.920,  2 ! 247Fm
     &, 95,153,1856.954,  1 ! 248Am
     &, 97,151,1857.890,  1 ! 248Bk
     &, 98,150,1857.854,  1 ! 248Cf
     &, 96,152,1859.273,  1 ! 248Cm
     &, 99,149,1854.095,  1 ! 248Es
     &,100,148,1851.641,  1 ! 248Fm
     &,101,147,1845.750,  1/! 248Md
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1931,1940)/
     &  97,152,1864.103,  2 ! 249Bk
     &, 98,151,1863.447,  2 ! 249Cf
     &, 96,153,1863.986,  2 ! 249Cm
     &, 99,150,1861.270,  2 ! 249Es
     &,100,149,1858.104,  2 ! 249Fm
     &,101,148,1853.561,  2 ! 249Md
     &, 97,153,1869.073,  1 ! 250Bk
     &, 98,152,1870.071,  1 ! 250Cf
     &, 96,154,1869.819,  1 ! 250Cm
     &, 99,151,1867.288,  1/! 250Es
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1941,1950)/
     & 100,150,1865.606,  1 ! 250Fm
     &,101,149,1860.293,  1 ! 250Md
     &, 97,154,1874.845,  2 ! 251Bk
     &, 98,153,1875.182,  2 ! 251Cf
     &, 99,152,1874.027,  2 ! 251Es
     &,100,151,1871.747,  2 ! 251Fm
     &,101,150,1867.935,  2 ! 251Md
     &, 97,155,1879.636,  1 ! 252Bk
     &, 98,154,1881.353,  1 ! 252Cf
     &, 99,153,1879.451,  1/! 252Es
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1951,1960)/
     & 100,152,1878.997,  1 ! 252Fm
     &,101,151,1874.536,  1 ! 252Md
     &,102,150,1871.392,  1 ! 252No
     &, 98,155,1886.157,  2 ! 253Cf
     &, 99,154,1885.661,  2 ! 253Es
     &,100,153,1884.545,  2 ! 253Fm
     &,101,152,1881.868,  2 ! 253Md
     &,102,151,1877.996,  2 ! 253No
     &, 98,156,1892.185,  1 ! 254Cf
     &, 99,155,1890.753,  1/! 254Es
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1961,1970)/
     & 100,154,1891.063,  1 ! 254Fm
     &,101,153,1887.790,  1 ! 254Md
     &,102,152,1885.668,  1 ! 254No
     &, 99,156,1896.736,  2 ! 255Es
     &,100,155,1896.238,  2 ! 255Fm
     &,103,152,1887.437,  2 ! 255Lr
     &,101,154,1894.371,  2 ! 255Md
     &,102,153,1891.599,  2 ! 255No
     &, 99,157,1901.628,  1 ! 256Es
     &,100,156,1902.625,  1/! 256Fm
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1971,1980)/
     & 103,153,1893.938,  1 ! 256Lr
     &,101,155,1899.903,  1 ! 256Md
     &,102,154,1898.740,  1 ! 256No
     &,100,157,1907.589,  2 ! 257Fm
     &,103,154,1900.860,  2 ! 257Lr
     &,101,156,1906.355,  2 ! 257Md
     &,102,155,1904.389,  2 ! 257No
     &,104,153,1897.097,  2 ! 257Rf
     &,103,155,1907.082,  1 ! 258Lr
     &,101,157,1911.647,  1/! 258Md
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1981,1990)/
     & 102,156,1911.164,  1 ! 258No
     &,104,154,1904.569,  1 ! 258Rf
     &,103,156,1913.973,  2 ! 259Lr
     &,102,157,1916.730,  2 ! 259No
     &,104,155,1910.691,  2 ! 259Rf
     &,105,155,1912.830,  1 ! 260Ha
     &,103,157,1919.905,  1 ! 260Lr
     &,104,156,1918.033,  1 ! 260Rf
     &,105,156,1920.092,  2 ! 261Ha
     &,104,157,1924.084,  2/! 261Rf
      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1991,1992)/
     & 105,157,1926.583,  1 ! 262Ha
     &,106,157,1929.603,  2/! 263Sg


*-----------------------------------------------------------------------
*     non exist nucleus
*-----------------------------------------------------------------------


      data (nzz0(i),nnn0(i),be0(i),jgs0(i),i=1993,2000)/
     &   1,  5,   0.000,  1 ! 6H
     &,  1,  6,   0.000,  1 ! 7H
     &,  1,  7,   0.000,  1 ! 8H
     &,  4,  1,   0.000,  1 ! 5Be
     &,  4,  2,   0.000,  1 ! 6Be
     &,  5,  1,   0.000,  1 ! 6B
     &,  6,  2,   0.000,  1 ! 7C
     &,  2,  7,   0.000,  1/! 9He


*-----------------------------------------------------------------------


      end




