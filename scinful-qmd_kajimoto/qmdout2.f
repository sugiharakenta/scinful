************************************************************************
*                                                                      *
*        PART 8: Display Utilities 1                                   *
*                                                                      *
*                                                                      *
*   List of subprograms in rough order of relevance with main purpose  *
*      ( s = subroutine, f = function, b = block data, e = entry )     *
*                                                                      *
*                                                                      *
*  s  disp01    to write all particle position and cluster contour     *
*               for ground state properties.                           *
*  s  disp02    to write all particle position and cluster contour     *
*               for reactions                                          *
*  s  disp03    to write particle position or momentum on file         *
*  s  disp04    to write color cluster plot of density for reactions   *
*  s  disp05    to write 2-d matrix on the file                        *
*  s  disp06    to display particle position on terminal               *
*  s  disp07    to write all particle momentum position and            *
*               Fermi sphere for ground state properties               *
*  s  disp08    to write all particle momentum position and            *
*               Fermi sphere for reactions                             *
*  s  densrr    to calculate spacial density                           *
*  s  cntdrw    to write contour line                                  *
*  s  cntsub    to make contour data and write them on file            *
*  s  cntout    to write contour data on file                          *
*                                                                      *
************************************************************************
*                                                                      *
************************************************************************
*                                                                      *
      subroutine disp01(mdim,scalo,scalp,isef,
     &                  time,itim,itdg,iefw)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to write all particle position and cluster contour      *
*              for ground state properties.                            *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              mdim        : size of display, -mdim to +mdim (fm)      *
*              scalo       : global scale                              *
*              scalp       : scale of the previous graph               *
*              isef        : =0; only position, =1; with momentum      *
*              time        : time (fm/c)                               *
*              itim        : =0; normal, =1; ending                    *
*              itdg        : =0; time is integer, =1; real             *
*              iefw        : =0; write data in one file                *
*                            =1; in each different file                *
*                                                                      *
*              xfx         : fixed x-position in a unit of x-axis      *
*              yfx         : fixed y-position in a unit of y-axis      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


      common /rannum/ iseed, iseed0, iseed1


      common /vriab0/ massal, massba, nmeson
      common /vriab1/ b, llnow, ntnow


      common /clusti/ itc(0:nnn,0:nnn)
      common /clustf/ nclst, iclust(nnn)


      common /summ00/ jdsp(50), idsp(50), fdsp(50), ifds(50)
      character    fdsp*80


      common /fnamid/ fname(mxfnm), idf(mxfnm), ifnl(mxfnm)
      character       fname*80


*-----------------------------------------------------------------------


      dimension       rho(mnrho,mnrho)
      dimension       dax(mnrho), day(mnrho)
      dimension       cuts(50), cols(50)


      dimension       it(0:nnn)


      dimension       icdp(10)


      character       ofnm*80


      data idp01 / 0 /


      save idp01


*-----------------------------------------------------------------------
*        multi events
*-----------------------------------------------------------------------


         if( time .lt. 0.0001 ) then


            ievnt = 1


         else


            ievnt = 0


         end if


*-----------------------------------------------------------------------
*        Out put file name and unit
*-----------------------------------------------------------------------


            io = idsp(1)


         if( iefw .eq. 1 ) then


            iq = idsp(1) + 1


            ofnm = fname(4)(1:ifnl(4))//'-gr'
            ifnm = ifnl(4) + 9


            ofnm(ifnm-3:ifnm) = '.ang'


         else


            iq = io


         end if


*-----------------------------------------------------------------------
*           Some constants
*-----------------------------------------------------------------------


               rmsh = 0.5
               imsh = 0


               afa0 = 1.0


*-----------------------------------------------------------------------


               maxx = nint( float(mdim) / rmsh )
               maxy = nint( float(mdim) / rmsh )


               mdimx = 2 * maxx
               mdimy = 2 * maxy


               xdd = rmsh
               ydd = rmsh


               ixnm = mdimx + 1
               iynm = mdimy + 1


               ipd2 = 0


               ncut = 2


               cuts(1) = rho0 * 1.0
               cuts(2) = rho0 * 0.1


               cols(1) = 3.0
               cols(2) = 2.0


*-----------------------------------------------------------------------


               do i = -maxx, maxx


                  ix = 1 + i + maxx


                  dax(ix) = float(i) * rmsh


               end do


               do i = -maxy, maxy


                  iy = 1 + i + maxy


                  day(iy) = float(i) * rmsh


               end do


*-----------------------------------------------------------------------
*        scal, form, and offset for ANGEL parameters
*-----------------------------------------------------------------------


               scal = 30.0 / float(mdimx) / rmsh
               form = 1.0


               sca0 = 0.45 * scalo


            if( isef .eq. 0 ) then


               xoff = -0.11
               yoff = -1.36


               xfx  =  0.5
               yfx  =  2.3


            else if( isef .eq. 1 ) then


               xoff = -0.66
               yoff = -1.36


               xfx  =  1.0
               yfx  =  2.3


            end if


               xorg =   4.0 / 14.0 * ( 1.0 / sca0 - scal )
     &                + xfx * ( 1.0 / sca0 - 1.0 )
     &                + 1.0 / sca0 * xoff


               yorg = ( 3.5 / 14.0 * ( 1.0 / sca0 - scal )
     &                + yfx * ( 1.0 / sca0 - form )
     &                + 1.0 / sca0 * yoff ) / form


               txps = 0.0
               typs = float(mdim)
     &              + float(mdim) * 0.2


*-----------------------------------------------------------------------
*        default values
*-----------------------------------------------------------------------


               sca00 =  1.00 * scalo


               xoff0 =  0.00
               yoff0 =  0.20


               xfx0  =  0.5
               yfx0  =  0.0


               xorg0 =   4.0 / 14.0 * ( 1.0 / sca00 - scal )
     &                 + xfx0 * ( 1.0 / sca00 - 1.0 )
     &                 + 1.0 / sca00 * xoff0


               yorg0 = ( 3.5 / 14.0 * ( 1.0 / sca00 - scal )
     &                + yfx0 * ( 1.0 / sca00 - form )
     &                + 1.0 / sca00 * yoff0 ) / form


*-----------------------------------------------------------------------
*     Write data on the file for ANGEL
*-----------------------------------------------------------------------


         if( ievnt .eq. 1 .and. llnow .gt. 1 ) then


            write(io,'( /''newpage:'')')


         end if


         if( isef .eq. 0 ) then


            write(io,'(/
     &          ''*'',71(''-'') /
     &          ''*     Snapshot of R-Space'',
     &               '' in Three Directions at T ='',f7.2,'' (fm/c)''/
     &          ''*     '',i2,'' th Event with iseed = '',i12/
     &          ''*'',71(''-''))') time, llnow, iseed1


            write(io,'(/
     &          ''msur:{\rm\it Snapshot of R-Space'',
     &          '' in Three Directions of \ ''/ i2,
     &          '' th Event'',i12,''}'')')
     &          llnow, iseed1


         else if( isef .eq. 1 ) then


            write(io,'(/
     &          ''*'',71(''-'') /
     &          ''*     Snapshot of R and P-Space'',
     &               '' in Three Directions at T ='',f7.2,'' (fm/c)''/
     &          ''*     '',i2,'' th Event with iseed = '',i12/
     &          ''*'',71(''-''))') time, llnow, iseed1


            write(io,'(/
     &          ''msur:{\rm\it Snapshot of R and P-Space'',
     &          '' in Three Directions of \ ''/ i2,
     &          '' th Event'',i12)')
     &          llnow, iseed1


            scalp = scal


         end if


            write(io,'(
     &          ''msul:{\rm\it JQMD  plotted by  \ANGEL}''/
     &          ''msdc:{\rm\Large \- \page  \-}'')')


*-----------------------------------------------------------------------
*     Three figures of z-x, z-y, x-y
*-----------------------------------------------------------------------


      do kk = 1, 3


         if( kk .ge. 2 ) then


            write(io,'(/''z:'')')


         end if


         if( iefw .eq. 1 ) then


            idp01 = idp01 + 1


            write(ofnm(ifnm-5:ifnm-4),'(i2.2)') idp01


            open( iq, file = ofnm, status='unknown')


            write(io,'( /''infl: {'',80a1)')
     &            ( ofnm(k:k), k = 1, ifnm ), '}'


         end if


*-----------------------------------------------------------------------


         if( kk .eq. 1 ) then


               kkx = 3
               kky = 1
               kkz = 2


            write(iq,'(/
     &          ''*'',71(''-'') /
     &          ''*     Snapshot of R-Space in z-x plane at T ='',f7.2,
     &          '' (fm/c)''/
     &          ''*     '',i2,'' th Event with iseed = '',i12/
     &          ''*'',71(''-'') )') time, llnow, iseed1


            write(iq,'(/''x: z (fm)'')')
            write(iq,'( ''y: x (fm)'')')


*-----------------------------------------------------------------------


         else if( kk .eq. 2 ) then


               kkx = 3
               kky = 2
               kkz = 1


            write(iq,'(/
     &          ''*'',71(''-'') /
     &          ''*     Snapshot of R-Space in z-y plane at T = '',f7.2,
     &          '' (fm/c)''/
     &          ''*     '',i2,'' th Event with iseed = '',i12/
     &          ''*'',71(''-'') )') time, llnow, iseed1


            write(iq,'(/''x: z (fm)'')')
            write(iq,'( ''y: y (fm)'')')


*-----------------------------------------------------------------------


         else if( kk .eq. 3 ) then


               kkx = 1
               kky = 2
               kkz = 3


            write(iq,'(/
     &          ''*'',71(''-'') /
     &          ''*     Snapshot of R-Space in x-y plane at T = '',f7.2,
     &          '' (fm/c)''/
     &          ''*     '',i2,'' th Event with iseed = '',i12/
     &          ''*'',71(''-'') )') time, llnow, iseed1


            write(iq,'(/''x: x (fm)'')')
            write(iq,'( ''y: y (fm)'')')


         end if


*-----------------------------------------------------------------------


            if( itdg .eq. 0 ) then


               write(iq,'( /''w: T = '',i3,''/ x('',f7.2,
     &                      '') y('',f7.2,'') ix(2) iy(1) s('',f7.4,
     &                      '')'')')
     &                   nint(time), txps, typs, afa0/scal*1.4


            else


               write(iq,'( /''w: T = '',f6.1,''/ x('',f7.2,
     &                      '') y('',f7.2,'') ix(2) iy(1) s('',f7.4,
     &                      '')'')')
     &                   time, txps, typs, afa0/scal*1.4


            end if


            write(iq,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       xorg0, yorg0


            write(iq,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                   '') xfac('',f7.4,'')'')')
     &                       scal * sca00, afa0 / scal, 1.0 / scal


            write(iq,'(  ''p: nofr noms nosp port form(1.0)'')')


            write(iq,
     &               '(  ''p: xmin('',i3,'') xmax('',i3,'') ymin('',i3,
     &                   '') ymax('',i3,'')'')')
     &                   -mdim, mdim, -mdim, mdim


*-----------------------------------------------------------------------
*        write particle data on file
*-----------------------------------------------------------------------


               do i = 1, 4


                  icdp(i) = 1


               end do


               call disp03(iq,kkx,kky,kkz,icdp,0)




*-----------------------------------------------------------------------
*        Calculate Density on grit points for clusters
*-----------------------------------------------------------------------


         do jj = 1, nclst


*-----------------------------------------------------------------------


               mclst = itc(jj,0)


*-----------------------------------------------------------------------


            if( mclst .gt. 1 ) then


*-----------------------------------------------------------------------
*           Calculate density on grit points
*-----------------------------------------------------------------------


                  it(0) = mclst


               do i = 1, mclst


                  it(i) = itc(jj,i)


               end do


                  call densrr(it,rho,kkx,kky,kkz,
     &                        -maxx,maxx,-maxy,maxy,imsh)


*-----------------------------------------------------------------------
*           write contour line data on file
*-----------------------------------------------------------------------


                  write(iq,'(/''c: Plot of density contour''/)')


               do i = 1, ncut


                  call cntdrw(iq,ixnm,iynm,xdd,ydd,dax,day,rho,cuts(i),
     &                        ipd2,cols(i))


               end do


*-----------------------------------------------------------------------


            end if


*-----------------------------------------------------------------------


         end do


*-----------------------------------------------------------------------


         if( iefw .eq. 1 ) then


            close(iq)


         end if


         if( io .eq. iq ) then


            write(io,'( /''*'',71(''-'')/
     &                   ''      Parameters for Multipage of Angel''/
     &                   ''*'',71(''-''))')


         end if


*-----------------------------------------------------------------------


         if( kk .eq. 1 ) then


            write(io,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       xorg, yorg


            write(io,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                   '') xfac('',f7.4,'')'')')
     &                       scal * sca0, afa0 / scal, 1.0 / scal


         else


            write(io,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       0.0, 1.3


            write(io,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                   '') xfac('',f7.4,'')'')')
     &                       1.0, afa0 / scal, 1.0 / scal


         end if


            if( kk .eq. 3 ) then


               write(io,'(  ''p: mssg'')')


            else


               write(io,'(  ''p: mssg nocm'')')


            end if


*-----------------------------------------------------------------------


      end do


*-----------------------------------------------------------------------


         if( isef .eq. 0 .and. itim .eq. 0 ) then


            write(io,'(/''newpage:'')')


         end if


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine disp02(nminx,nmaxx,nminy,nmaxy,
     &                  scalo,scalp,formp,ised,isef,mulg,
     &                  time,itim,itdg,iefw)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to write all particle position and cluster contour      *
*              for reactions                                           *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              nminx,nmaxx : size of display, -nmin to +nmax (fm)      *
*              nminx,nmaxx                                             *
*              scalo       : global scale                              *
*              scalp       : scale of the previous graph               *
*              formp       : form of the previous graph                *
*              ised        : =1; z-x, =2; z-y, =3; x-y                 *
*              isef        : =0; only this, =1; with momentum,         *
*                            =2; multi-graph in one page               *
*                            =3; multi-graph with momentum in one page *
*                            =4; with color plot                       *
*                            =5; multi-graph with color plot           *
*              mulg        : number of graphs in a page                *
*                            =< 0; default                             *
*              time        : time (fm/c)                               *
*              itim        : =0; normal, =1; ending                    *
*              itdg        : =0; time is integer, =1; real             *
*              iefw        : =0; write data in one file                *
*                            =1; in each different file                *
*                                                                      *
*              xfx         : fixed x-position in a unit of x-axis      *
*              yfx         : fixed y-position in a unit of y-axis      *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


      common /rannum/ iseed, iseed0, iseed1


      common /vriab0/ massal, massba, nmeson
      common /vriab1/ b, llnow, ntnow


      common /clusti/ itc(0:nnn,0:nnn)
      common /clustf/ nclst, iclust(nnn)


      common /summ00/ jdsp(50), idsp(50), fdsp(50), ifds(50)
      character    fdsp*80


      common /fnamid/ fname(mxfnm), idf(mxfnm), ifnl(mxfnm)
      character       fname*80


*-----------------------------------------------------------------------


      dimension       rho(mnrho,mnrho)
      dimension       dax(mnrho), day(mnrho)
      dimension       cuts(50), cols(50)


      dimension       it(0:nnn)


      dimension       icdp(10)


*-----------------------------------------------------------------------


      character       ofnm*80


      data idp04 / 0 /


      save idp04


*-----------------------------------------------------------------------


      save imult


*-----------------------------------------------------------------------
*        multi events
*-----------------------------------------------------------------------


         if( time .lt. 0.0001 ) then


            ievnt = 1
            imult = 0


         else


            ievnt = 0


         end if


*-----------------------------------------------------------------------
*        Out put file name and unit
*-----------------------------------------------------------------------


            io = idsp(4)


         if( iefw .eq. 1 ) then


            iq = idsp(4) + 1


            ofnm = fname(4)(1:ifnl(4))//'-r0'
            ifnm = ifnl(4) + 9


            ofnm(ifnm-3:ifnm) = '.ang'


         else


            iq = io


         end if


*-----------------------------------------------------------------------
*           Some constants
*-----------------------------------------------------------------------


               rmsh = 0.5
               imsh = 0


               mul  = mulg
               afa0 = 0.8


*-----------------------------------------------------------------------
*        size of display; form and scal
*-----------------------------------------------------------------------


               mminx = nint( float(nminx) / rmsh )
               mmaxx = nint( float(nmaxx) / rmsh )
               mminy = nint( float(nminy) / rmsh )
               mmaxy = nint( float(nmaxy) / rmsh )


               mdimx = mmaxx - mminx
               mdimy = mmaxy - mminy


               xdd = rmsh
               ydd = rmsh


               ixnm = mdimx + 1
               iynm = mdimy + 1


               ipd2 = 0


               ncut = 1


               cuts(1) = rho0 * 1.0
               cuts(2) = rho0 * 0.1


               cols(1) = 3.0
               cols(2) = 2.0




               form = float(mdimy) / float(mdimx)
               scal =  30.0 / float(mdimx) / rmsh


*-----------------------------------------------------------------------
*        offset and comment position
*-----------------------------------------------------------------------


            if( mul .le. 0 ) then


               mul = max( 2, int( 3.0 / form + 0.01 ) * 2 )


            end if


            if( isef .eq. 0 ) then


               sca0 =  1.0 * scalo


               xoff =  0.0
               yoff =  0.3


               xfx  =  0.5
               yfx  =  0.0


            else if( isef .eq. 1 .or. isef .eq. 4 ) then


               sca0 =  0.7 * scalo


               xoff = -0.1
               yoff = -0.37


               xfx  =  0.5
               yfx  =  1.0


            else if( isef .eq. 2 ) then


               sca0 =  0.52 * scalo


               xoff = -0.64
               yoff =  0.57


               xfx  =  1.0
               yfx  =  1.0


               xds0 =  0.0
               xds1 =  1.3
               yds1 = -1.1


               yofg =  1.1 * ( mul / 2 - 1 )


               imult = imult + 1


            else if( isef .eq. 3 .or. isef .eq. 5 ) then


               sca0 =  0.52 * scalo


               if( scalp .gt. 0.0 ) then


                  sca0 = 1.0 / scalp


               end if


               xoff = -0.64
               yoff =  0.57


               xfx  =  1.0
               yfx  =  1.0


               xds0 =  0.0
               xds1 = -1.3
               yds1 = -1.1


               yofg = 1.1 * ( mul / 2 - 1 )


               imult = imult + 1


            end if


               xorg =   4.0 / 14.0 * ( 1.0 / sca0 - scal )
     &                + xfx * ( 1.0 / sca0 - 1.0 )
     &                + 1.0 / sca0 * xoff


               yorg = ( 3.5 / 14.0 * ( 1.0 / sca0 - scal )
     &                + yfx * ( 1.0 / sca0 - form )
     &                + 1.0 / sca0 * yoff ) / form


               txps = float(nminx)
     &              + float( nmaxx - nminx ) * 0.05
               typs = float(nmaxy)
     &              - float( nmaxy - nminy ) * 0.05 / form


*-----------------------------------------------------------------------
*        default values
*-----------------------------------------------------------------------


               sca00 =  1.0 * scalo


               xoff0 =  0.0
               yoff0 =  0.3


               xfx0  =  0.5
               yfx0  =  0.0


               xorg0 =   4.0 / 14.0 * ( 1.0 / sca00 - scal )
     &                 + xfx0 * ( 1.0 / sca00 - 1.0 )
     &                 + 1.0 / sca00 * xoff0


               yorg0 = ( 3.5 / 14.0 * ( 1.0 / sca00 - scal )
     &                 + yfx0 * ( 1.0 / sca00 - form )
     &                 + 1.0 / sca00 * yoff0 ) / form


*-----------------------------------------------------------------------


               do i = mminx, mmaxx


                  ix = 1 + i - mminx


                  dax(ix) = float(i) * rmsh


               end do


               do i = mminy, mmaxy


                  iy = 1 + i - mminy


                  day(iy) = float(i) * rmsh


               end do


*-----------------------------------------------------------------------
*        write ANGEL parameters
*-----------------------------------------------------------------------


      if( ievnt .eq. 1 .and. llnow .gt. 1 ) then


            write(io,'( /''newpage:'')')


      end if


      if( isef .eq. 0 ) then


            write(io,'(/
     &          ''*'',71(''-'') /
     &          ''*     Snapshot of R-Space'',
     &               '' at T ='',f7.2,'' (fm/c)''/
     &          ''*     '',i2,'' th Event with iseed = '',i12/
     &          ''*'',71(''-''))') time, llnow, iseed1


            write(io,'(/
     &          ''msur:{\rm\it Snapshot of R-Space of \ ''/ i2,
     &          '' th Event'',i12,''}''/
     &          ''msul:{\rm\it JQMD  plotted by  \ANGEL}''/
     &          ''msdc:{\rm\Large \- \page  \-}'')')
     &          llnow, iseed1


      else if( isef .eq. 1 ) then


            write(io,'(/
     &          ''*'',71(''-'') /
     &          ''*     Snapshot of P and R-Space'',
     &               '' at T ='',f7.2,'' (fm/c)''/
     &          ''*     '',i2,'' th Event with iseed = '',i12/
     &          ''*'',71(''-''))') time, llnow, iseed1


            write(io,'(/
     &          ''msur:{\rm\it Snapshot of P and R-Space of \ ''/ i2,
     &          '' th Event'',i12,''}''/
     &          ''msul:{\rm\it JQMD  plotted by  \ANGEL}''/
     &          ''msdc:{\rm\Large \- \page  \-}'')')
     &          llnow, iseed1


      else if( isef .eq. 2 ) then


         if( ( imult - 1 ) / mul * mul .eq. imult - 1 ) then


            write(io,'(/
     &          ''*'',71(''-'') /
     &          ''*     Time Evolution of R-Space''/
     &          ''*     '',i2,'' th Event with iseed = '',i12/
     &          ''*'',71(''-''))') llnow, iseed1


            write(io,'(/
     &          ''msur:{\rm\it Time Evolution of R-Space of \ ''/ i2,
     &          '' th Event'',i12,''}''/
     &          ''msul:{\rm\it JQMD  plotted by  \ANGEL}''/
     &          ''msdc:{\rm\Large \- \page  \-}'')')
     &          llnow, iseed1


         else


            write(io,'( /''z:'')')


         end if


      else if( isef .eq. 3 ) then


         if( ( imult - 1 ) / ( mul / 2 ) * ( mul / 2 ) .eq.
     &         imult - 1 ) then


            write(io,'(/
     &          ''*'',71(''-'') /
     &          ''*     Time Evolution of R and P-Space''/
     &          ''*     '',i2,'' th Event with iseed = '',i12/
     &          ''*'',71(''-''))') llnow, iseed1


            write(io,'(/
     &          ''msur:{\rm\it Time Evolution of'',
     &          '' R and P-Space of \ ''/ i2,
     &          '' th Event'',i12,''}''/
     &          ''msul:{\rm\it JQMD  plotted by  \ANGEL}''/
     &          ''msdc:{\rm\Large \- \page  \-}'')')
     &          llnow, iseed1


         else


            write(io,'( /''z:'')')


         end if


      else if( isef .eq. 4 ) then


            write(io,'(/
     &          ''*'',71(''-'') /
     &          ''*     Snapshot of R-Space with'',
     &               '' Contour and Color'',
     &               '' at T ='',f7.2,'' (fm/c)''/
     &          ''*     '',i2,'' th Evnet with iseed = '',i12/
     &          ''*'',71(''-''))') time, llnow, iseed1


            write(io,'(/
     &          ''msur:{\rm\it Snapshot of R-Space'',
     &          '' with Contour and Color Plot of \ ''/ i2,
     &          '' th Event'',i12,''}''/
     &          ''msul:{\rm\it JQMD  plotted by  \ANGEL}''/
     &          ''msdc:{\rm\Large \- \page  \-}'')')
     &          llnow, iseed1


      else if( isef .eq. 5 ) then


         if( ( imult - 1 ) / ( mul / 2 ) * ( mul / 2 ) .eq.
     &         imult - 1 ) then


            write(io,'(/
     &          ''*'',71(''-'') /
     &          ''*     Time Evolution of R-Space'',
     &               '' with Contour and Color Plot''/
     &          ''*     '',i2,'' th Event with iseed = '',i12/
     &          ''*'',71(''-''))') llnow, iseed1


            write(io,'(/
     &          ''msur:{\rm\it Time Evolution of'',
     &          '' R-Space with Contour and Color Plot of \ ''/ i2,
     &          '' th Event'',i12,''}''/
     &          ''msul:{\rm\it JQMD  plotted by  \ANGEL}''/
     &          ''msdc:{\rm\Large \- \page  \-}'')')
     &          llnow, iseed1


         else


            write(io,'( /''z:'')')


         end if


      end if


*-----------------------------------------------------------------------


         if( iefw .eq. 1 ) then


            idp04 = idp04 + 1


            write(ofnm(ifnm-5:ifnm-4),'(i2.2)') idp04


            open( iq, file = ofnm, status='unknown')


            write(io,'( /''infl: {'',80a1)')
     &            ( ofnm(k:k), k = 1, ifnm ), '}'


         end if


*-----------------------------------------------------------------------


            write(iq,'(/
     &          ''*'',71(''-'') /
     &          ''*     Snapshot of R-Space'',
     &               '' at T ='',f7.2,'' (fm/c)''/
     &          ''*     '',i2,'' th Event with iseed = '',i12/
     &          ''*'',71(''-''))') time, llnow, iseed1


*-----------------------------------------------------------------------
*        X and Y axis
*-----------------------------------------------------------------------


            if( ised .eq. 1 ) then


               kkx = 3
               kky = 1
               kkz = 2


               write(iq,'(/''x: z (fm)'')')
               write(iq,'( ''y: x (fm)'')')


            else if( ised .eq. 2 ) then


               kkx = 3
               kky = 2
               kkz = 1


               write(iq,'(/''x: z (fm)'')')
               write(iq,'( ''y: y (fm)'')')


            else if( ised .eq. 3 ) then


               kkx = 1
               kky = 2
               kkz = 3


               write(iq,'(/''x: x (fm)'')')
               write(iq,'( ''y: y (fm)'')')


            end if


*-----------------------------------------------------------------------
*        Time and Frame
*-----------------------------------------------------------------------


            if( itdg .eq. 0 ) then


               write(iq,'( /''w: T = '',i3,''/ x('',f7.2,
     &                      '') y('',f7.2,'') ix(1) iy(3) s('',f7.4,
     &                      '')'')')
     &                   nint(time), txps, typs, afa0/scal*1.4


            else


               write(iq,'( /''w: T = '',f6.1,''/ x('',f7.2,
     &                      '') y('',f7.2,'') ix(1) iy(3) s('',f7.4,
     &                      '')'')')
     &                   time, txps, typs, afa0/scal*1.4


            end if


               write(iq,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       xorg0, yorg0


               write(iq,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                      '') xfac('',f7.4,'')'')')
     &                       scal * sca00, afa0 / scal, 1.0 / scal


               write(iq,'(  ''p: nofr noms nosp port form('',
     &                    f7.4,'')'')') form


               write(iq,
     &             '(''p: xmin('',i3,'') xmax('',i3,'') ymin('',i3,
     &             '') ymax('',i3,'')'')') nminx, nmaxx, nminy, nmaxy




*-----------------------------------------------------------------------
*        write particle data on file
*-----------------------------------------------------------------------


               do i = 1, 4


                  icdp(i) = 1


               end do


               call disp03(iq,kkx,kky,kkz,icdp,0)




*-----------------------------------------------------------------------
*        Calculate Density on grit points for clusters
*-----------------------------------------------------------------------


         do jj = 1, nclst


*-----------------------------------------------------------------------


               mclst = itc(jj,0)


*-----------------------------------------------------------------------


            if( mclst .gt. 1 ) then


*-----------------------------------------------------------------------
*           Calculate density on grit points
*-----------------------------------------------------------------------


                  it(0) = mclst


               do i = 1, mclst


                  it(i) = itc(jj,i)


               end do


                  call densrr(it,rho,kkx,kky,kkz,
     &                        mminx,mmaxx,mminy,mmaxy,imsh)


*-----------------------------------------------------------------------
*           write contour line data on file
*-----------------------------------------------------------------------


               do i = 1, ncut


                  call cntdrw(iq,ixnm,iynm,xdd,ydd,dax,day,rho,cuts(i),
     &                        ipd2,cols(i))


               end do


*-----------------------------------------------------------------------


            end if


*-----------------------------------------------------------------------


         end do


*-----------------------------------------------------------------------


         if( iefw .eq. 1 ) then


            close(iq)


         end if


         if( io .eq. iq ) then


            write(io,'( /''*'',71(''-'')/
     &                   ''      Parameters for Multipage of Angel''/
     &                   ''*'',71(''-''))')


         end if


*-----------------------------------------------------------------------


      if( isef .lt. 2 .or. isef .eq. 4 ) then


            write(io,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       xorg, yorg


            write(io,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                   '') xfac('',f7.4,'')'')')
     &                    scal * sca0, afa0 / scal, 1.0 / scal


            write(io,'(  ''p: mssg'')')


*-----------------------------------------------------------------------


      else if( isef .eq. 2 ) then


         if( ( imult - 1 ) / mul * mul .eq. imult - 1 ) then


            write(io,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       xorg, yorg


            write(io,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                   '') xfac('',f7.4,'')'')')
     &                    scal * sca0, afa0 / scal, 1.0 / scal


            write(io,'(  ''p: mssg'')')


         else if( ( imult - mul / 2 - 1 ) / mul * mul .eq.
     &              imult - mul / 2 - 1 ) then


            write(io,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       xds1, yofg


            write(io,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                   '') xfac('',f7.4,'')'')')
     &                    1.0, afa0 / scal, 1.0 / scal


         else


            write(io,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       xds0, yds1


            write(io,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                   '') xfac('',f7.4,'')'')')
     &                    1.0, afa0 / scal, 1.0 / scal


         end if


         if( imult / ( mul / 2 ) * ( mul / 2 ) .ne. imult .and.
     &       itim .eq. 0 ) then


            write(io,'(  ''p: noxn noxt'')')


         end if


*-----------------------------------------------------------------------


      else if( isef .eq. 3 .or. isef .eq. 5 ) then


         if( ( imult - 1 ) / ( mul / 2 ) * ( mul / 2 ) .eq.
     &         imult - 1 ) then


            write(io,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       xorg, yorg


            write(io,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                   '') xfac('',f7.4,'')'')')
     &                    scal * sca0, afa0 / scal, 1.0 / scal


            write(io,'(  ''p: mssg'')')


         else


            write(io,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       xds1, yds1


            write(io,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                   '') xfac('',f7.4,'')'')')
     &                    scal * sca0, afa0 / scal, 1.0 / scal


         end if


         if( imult / ( mul / 2 ) * ( mul / 2 ) .ne. imult .and.
     &       itim .eq. 0 ) then


            write(io,'(  ''p: noxn noxt'')')


         end if


         if( isef .eq. 5 ) then


            write(io,'(  ''p: clin(w)'')')


         end if


      end if


*-----------------------------------------------------------------------


            if( isef .eq. 1 .or. isef .eq. 3 .or.
     &          isef .eq. 4 .or. isef .eq. 5 ) then


               formp = form
               scalp = scal


            end if


*-----------------------------------------------------------------------


         if( isef .eq. 0 .and.
     &       itim .eq. 0 ) then


            write(io,'(/''newpage:'')')


         else if( isef .eq. 2 .and.
     &            imult / mul * mul .eq. imult .and.
     &            itim .eq. 0 ) then


            write(io,'(/''newpage:'')')


         end if


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine disp03(io,nx,ny,nz,icdp,irp)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to write particle position or momentum on file          *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              io          : io of the out put file                    *
*              nx, ny, nz  : x, y axis and integrated axis             *
*              icdp        : flag of display particles                 *
*              irp         : flag of coordinate or momentum            *
*                                    (0)           (1)                 *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


      common /vriab0/ massal, massba, nmeson


      common /coodrp/ r(3,nnn),  p(5,nnn)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)


      dimension       iren(nnn), irem(nnn)


      dimension       icdp(10)


*-----------------------------------------------------------------------
*        renumbering of the particle id according to the y-position
*-----------------------------------------------------------------------


            do i = 1, massal


               irem(i) = -1


            end do




         do j = 1, massal


                     posm = 1.0e30


            do i = 1, massal


               if( irem(i) .eq. -1 ) then


                  if( irp .eq. 0 ) then


                     rpv = r(nz,i)


                  else


                     rpv = p(nz,i)


                  end if


                  if( rpv .lt. posm ) then


                     irenp = i
                     posm  = rpv


                  end if


               end if


            end do


                     iren(j)     = irenp
                     irem(irenp) = j


         end do


*-----------------------------------------------------------------------
*        write particle position on file(io)
*-----------------------------------------------------------------------


         if( irp .eq. 0 ) then


            write(io,'(/''c: Plot of particle position''/)')


         else


            write(io,'(/''c: Plot of particle momentum''/)')


         end if


                  indp = -1
                  ichp = 10


         do 100 j = 1, massal


                  i = iren(j)


            if( icdp(inds(i)) .eq. 0 ) goto 100


            if( ( inds(i) .ne. 1 .and. inds(i) .ne. indp ) .or.
     &          ( inds(i) .eq. 1 .and. ichp .ne. ichg(i) ) ) then


               if( inds(i) .eq. 1 .and. ichg(i) .eq. 1 ) then


                  write(io,'(''h:  x     y,n3[o]aa'')')


                  ichp = ichg(i)


               else if( inds(i) .eq. 1 .and. ichg(i) .eq. 0 ) then


                  write(io,'(''h:  x     y,n3aa'')')


                  ichp = ichg(i)


               else if( inds(i) .eq. 2 ) then


                  write(io,'(''h:  x     y,n7[r]aa'')')


                  ichp = 10


               else if( inds(i) .eq. 3 ) then


                  write(io,'(''h:  x     y,n9[r]aa'')')


                  ichp = 10


               else if( inds(i) .eq. 4 ) then


                  write(io,'(''h:  x     y,n13[rrr]a'')')


                  ichp = 10


               end if


                  indp = inds(i)


            end if


                  if( irp .eq. 0 ) then


                     rpvx = r(nx,i)
                     rpvy = r(ny,i)


                  else


                     rpvx = p(nx,i) / hbc
                     rpvy = p(ny,i) / hbc


                  end if


                  write(io,'(2e12.4)') rpvx, rpvy


  100    continue


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine disp04(nminx,nmaxx,nminy,nmaxy,
     &                  scalo,scalp,formp,ised,isef,mulg,
     &                  time,itim,itdg,icdp,iefw)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to write color cluster plot of density  for reactions   *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              io          : io of the out put file                    *
*              nminx,nmaxx : size of display, -nmin to +nmax (fm)      *
*              nminx,nmaxx                                             *
*              scalo       : global scale                              *
*              scalp       : scale of the previous graph               *
*              formp       : form of the previous graph                *
*              ised        : =1; z-x, =2; z-y, =3; x-y                 *
*              isef        : =0; only this,                            *
*                            =2; multi-graph in one page               *
*                            =4; with contour plot                     *
*                            =5; multi-graph with contour plot         *
*              mulg        : number of graphs in a page                *
*                            =< 0; default                             *
*              time        : time (fm/c)                               *
*              itim        : =0; normal, =1; ending                    *
*              itdg        : =0; time is integer, =1; real             *
*                                                                      *
*              icdp        : display particle on color plot            *
*              iefw        : =0; write data in one file                *
*                            =1; in each different file                *
*                                                                      *
*              xfx         : fixed x-position in a unit of x-axis      *
*              yfx         : fixed y-position in a unit of y-axis      *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


      common /rannum/ iseed, iseed0, iseed1


      common /vriab0/ massal, massba, nmeson
      common /vriab1/ b, llnow, ntnow


      common /clusti/ itc(0:nnn,0:nnn)
      common /clustf/ nclst, iclust(nnn)


      common /summ00/ jdsp(50), idsp(50), fdsp(50), ifds(50)
      character    fdsp*80


      common /fnamid/ fname(mxfnm), idf(mxfnm), ifnl(mxfnm)
      character       fname*80


*-----------------------------------------------------------------------


      dimension       rho(mnrho,mnrho)


      dimension       it(0:nnn)


      dimension       icdp(10)


      data ihigh /0/
      save highi


      save imult


      character       ofnm*80


      data idp10 / 0 /


      save idp10


*-----------------------------------------------------------------------
*        multi events
*-----------------------------------------------------------------------


         if( time .lt. 0.0001 ) then


            ievnt = 1
            imult = 0


         else


            ievnt = 0


         end if


*-----------------------------------------------------------------------
*        Out put file name and unit
*-----------------------------------------------------------------------


            io = idsp(10)


         if( iefw .eq. 1 ) then


            iq = idsp(10) + 1


            ofnm = fname(4)(1:ifnl(4))//'-rc'
            ifnm = ifnl(4) + 9


            ofnm(ifnm-3:ifnm) = '.ang'


         else


            iq = io


         end if


*-----------------------------------------------------------------------
*           Some constants
*-----------------------------------------------------------------------


               rmsh = 0.5
               imsh = 0


               mul  = mulg
               afa0 = 0.8
               cma0 = 50.0


*-----------------------------------------------------------------------
*        size of display
*-----------------------------------------------------------------------


               mminx = nint( float(nminx) / rmsh )
               mmaxx = nint( float(nmaxx) / rmsh )
               mminy = nint( float(nminy) / rmsh )
               mmaxy = nint( float(nmaxy) / rmsh )


               mdimx = mmaxx - mminx
               mdimy = mmaxy - mminy


               xdd = rmsh
               ydd = rmsh


               ixnm = mdimx + 1
               iynm = mdimy + 1


               form = float(mdimy) / float(mdimx)
               scal =  30.0 / float(mdimx) / rmsh


*-----------------------------------------------------------------------
*        offset and comment position
*-----------------------------------------------------------------------


            if( mul .le. 0 ) then


               mul = max( 2, int( 3.0 / form + 0.01 ) * 2 )


            end if


            if( isef .eq. 0 ) then


               sca0 = 1.0 * scalo


               xoff = 0.0
               yoff = 0.3


               xfx  =  0.5
               yfx  =  0.0


            else if( isef .eq. 2 ) then


               sca0 =  0.52 * scalo


               xoff = -0.64
               yoff =  0.57


               xfx  =  1.0
               yfx  =  1.0


               xds0 =  0.0
               xds1 =  1.3
               yds1 = -1.1


               yofg = 1.1 * ( mul / 2 - 1 )


               imult = imult + 1


            else if( isef .eq. 4 ) then


               sca0 = 1.0 / scalp


               xoff = 0.0
               yoff = 1.0 + 0.3 / formp


            else if( isef .eq. 5 ) then


               sca0 = 1.0 / scalp


               xds1 =  1.3
               yds1 =  0.0


               imult = imult + 1


            end if


               xorg =   4.0 / 14.0 * ( 1.0 / sca0 - scal )
     &                + xfx * ( 1.0 / sca0 - 1.0 )
     &                + 1.0 / sca0 * xoff


               yorg = ( 3.5 / 14.0 * ( 1.0 / sca0 - scal )
     &                + yfx * ( 1.0 / sca0 - form )
     &                + 1.0 / sca0 * yoff ) / form


               txps = float(nminx)
     &              + float( nmaxx - nminx ) * 0.05
               typs = float(nmaxy)
     &              - float( nmaxy - nminy ) * 0.05 / form


*-----------------------------------------------------------------------
*        default values
*-----------------------------------------------------------------------


               sca00 =  1.0 * scalo


               xoff0 =  0.0
               yoff0 =  0.3


               xfx0  =  0.5
               yfx0  =  0.0


               xorg0 =   4.0 / 14.0 * ( 1.0 / sca00 - scal )
     &                 + xfx0 * ( 1.0 / sca00 - 1.0 )
     &                 + 1.0 / sca00 * xoff0


               yorg0 = ( 3.5 / 14.0 * ( 1.0 / sca00 - scal )
     &                 + yfx0 * ( 1.0 / sca00 - form )
     &                 + 1.0 / sca00 * yoff0 ) / form


*-----------------------------------------------------------------------
*        write ANGEL parameters
*-----------------------------------------------------------------------


      if( ievnt .eq. 1 .and. llnow .gt. 1 .and.
     &  ( isef .eq. 0 .or. isef .eq. 2 ) ) then


            write(io,'( /''newpage:'')')


      end if


      if( isef .eq. 0 ) then


            write(io,'(/
     &          ''*'',71(''-'') /
     &          ''*     Snapshot of R-Space with Color Plot'',
     &               '' at T ='',f7.2,'' (fm/c)''/
     &          ''*     '',i2,'' th Event with iseed = '',i12/
     &          ''*'',71(''-''))') time, llnow, iseed1


            write(io,'(/
     &          ''msur:{\rm\it Snapshot of R-Space of \ ''/ i2,
     &          '' th Event'',i12,''}''/
     &          ''msul:{\rm\it JQMD  plotted by  \ANGEL}''/
     &          ''msdc:{\rm\Large \- \page  \-}'')') llnow, iseed1


      else if( isef .eq. 2 ) then


         if( ( imult - 1 ) / mul * mul .eq. imult - 1 ) then


            write(io,'(/
     &          ''*'',71(''-'') /
     &          ''*     Time Evolution of R-Space with Color Plot''/
     &          ''*     '',i2,'' th Event with iseed = '',i12/
     &          ''*'',71(''-''))') llnow, iseed1


            write(io,'(/
     &          ''msur:{\rm\it Time Evolution of R-Space of \ ''/ i2,
     &          '' th Event'',i12,''}''/
     &          ''msul:{\rm\it JQMD  plotted by  \ANGEL}''/
     &          ''msdc:{\rm\Large \- \page  \-}'')') llnow, iseed1


         else


            write(io,'( /''z:'')')


         end if


      else if( isef .eq. 4 ) then


            write(io,'( /''z:'')')


      else if( isef .eq. 5 ) then


            write(io,'( /''z:'')')


      end if


*-----------------------------------------------------------------------


         if( iefw .eq. 1 ) then


            idp10 = idp10 + 1


            write(ofnm(ifnm-5:ifnm-4),'(i2.2)') idp10


            open( iq, file = ofnm, status='unknown')


            write(io,'( /''infl: {'',80a1)')
     &            ( ofnm(k:k), k = 1, ifnm ), '}'


         end if


*-----------------------------------------------------------------------


            write(iq,'(/
     &          ''*'',71(''-'') /
     &          ''*     Snapshot of R-Space with Color Plot'',
     &               '' at T ='',f7.2,'' (fm/c)''/
     &          ''*     '',i2,'' th Event with iseed = '',i12/
     &          ''*'',71(''-''))') time, llnow, iseed1


*-----------------------------------------------------------------------
*        X and Y axis
*-----------------------------------------------------------------------


            if( ised .eq. 1 ) then


               kkx = 3
               kky = 1
               kkz = 2


               write(iq,'(/''x: z (fm)'')')
               write(iq,'( ''y: x (fm)'')')


            else if( ised .eq. 2 ) then


               kkx = 3
               kky = 2
               kkz = 1


               write(iq,'(/''x: z (fm)'')')
               write(iq,'( ''y: y (fm)'')')


            else if( ised .eq. 3 ) then


               kkx = 1
               kky = 2
               kkz = 3


               write(iq,'(/''x: x (fm)'')')
               write(iq,'( ''y: y (fm)'')')


            end if


*-----------------------------------------------------------------------
*        Time and Frame
*-----------------------------------------------------------------------


            if( itdg .eq. 0 ) then


               write(iq,'( /''w: T = '',i3,''/ x('',f7.2,
     &                      '') y('',f7.2,'') ix(1) iy(3) s('',f7.4,
     &                      '') c(w)'')')
     &                   nint(time), txps, typs, afa0 / scal * 1.4


            else


               write(iq,'( /''w: T = '',f6.1,''/ x('',f7.2,
     &                      '') y('',f7.2,'') ix(1) iy(3) s('',f7.4,
     &                      '') c(w)'')')
     &                   time, txps, typs, afa0 / scal * 1.4


            end if


               write(iq,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       xorg0, yorg0


               write(iq,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                      '') xfac('',f7.4,'')'')')
     &                       scal * sca00, afa0 / scal, 1.0 / scal


               write(iq,'(  ''p: nofr noms nosp port form('',
     &                   f7.4,'')'')') form


               write(iq,
     &             '(''p: xmin('',i3,'') xmax('',i3,'') ymin('',i3,
     &             '') ymax('',i3,'')'')') nminx, nmaxx, nminy, nmaxy


*-----------------------------------------------------------------------
*           Calculate density on grit points for all barions
*-----------------------------------------------------------------------


                  it(0) = massba


               do i = 1, massba


                  it(i) = i


               end do


                  call densrr(it,rho,kkx,kky,kkz,
     &                        mminx,mmaxx,mminy,mmaxy,imsh)


*-----------------------------------------------------------------------
*           Maximum valu of the density rho(ix,iy)
*-----------------------------------------------------------------------


                  high = 0.0


               do ix = 1, ixnm
               do iy = 1, iynm


                  if( rho(ix,iy) .gt. high ) high = rho(ix,iy)


               end do
               end do


               if( ihigh .eq. 0 ) then


                  highi = high
                  ihigh = 1


                  cmax  = cma0


               else


                  cmax  = cma0 * highi / high


               end if


*-----------------------------------------------------------------------
*        write color culuster data on file
*-----------------------------------------------------------------------


               write(iq,'(/''p: ipdc clin(b) dmin(0.01)'')')
               write(iq,'( ''p: cmax('',f7.1,'')'')') cmax


               write(iq,
     &             '(/''hc: y = '',i3,'' to '',i3,'' by '',f7.4,
     &                '' ; x = '',i3,'' to '',i3,'' by '',f7.4,
     &              '' ;'')') nmaxy, nminy, rmsh, nminx, nmaxx, rmsh




               call disp05(iq,1,ixnm,1,iynm,rho,high)




*-----------------------------------------------------------------------
*        write pion position on file
*-----------------------------------------------------------------------




               call disp03(iq,kkx,kky,kkz,icdp,0)




*-----------------------------------------------------------------------


         if( iefw .eq. 1 ) then


            close(iq)


         end if


         if( io .eq. iq ) then


            write(io,'( /''*'',71(''-'')/
     &                   ''      Parameters for Multipage of Angel''/
     &                   ''*'',71(''-''))')


         end if


*-----------------------------------------------------------------------
*        write ANGEL parameters
*-----------------------------------------------------------------------


      if( isef .eq. 0 ) then


            write(io,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       xorg, yorg


            write(io,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                   '') xfac('',f7.4,'')'')')
     &                    scal * sca0, afa0 / scal, 1.0 / scal


            write(io,'(  ''p: mssg'')')


*-----------------------------------------------------------------------


      else if( isef .eq. 2 ) then


         if( ( imult - 1 ) / mul * mul .eq. imult - 1 ) then


            write(io,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       xorg, yorg


            write(io,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                   '') xfac('',f7.4,'')'')')
     &                    scal * sca0, afa0 / scal, 1.0 / scal


            write(io,'(  ''p: mssg'')')


         else if( ( imult - mul / 2 - 1 ) / mul * mul .eq.
     &              imult - mul / 2 - 1 ) then


            write(io,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       xds1, yofg


            write(io,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                   '') xfac('',f7.4,'')'')')
     &                    1.0, afa0 / scal, 1.0 / scal


         else


            write(io,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       xds0, yds1


            write(io,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                   '') xfac('',f7.4,'')'')')
     &                    1.0, afa0 / scal, 1.0 / scal


         end if


         if( imult / ( mul / 2 ) * ( mul / 2 ) .ne. imult .and.
     &       itim .eq. 0 ) then


            write(io,'(  ''p: noxn noxt'')')


         end if


*-----------------------------------------------------------------------


      else if( isef .eq. 4 ) then


            write(io,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       xoff, yoff


            write(io,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                   '') xfac('',f7.4,'')'')')
     &                       scal * sca0, afa0 / scal, 1.0 / scal


*-----------------------------------------------------------------------


      else if( isef .eq. 5 ) then


            write(io,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       xds1, yds1


            write(io,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                   '') xfac('',f7.4,'')'')')
     &                       scal * sca0, afa0 / scal, 1.0 / scal


         if( imult / ( mul / 2 ) * ( mul / 2 ) .ne. imult .and.
     &       itim .eq. 0 ) then


            write(io,'(  ''p: noxn noxt'')')


            scalp = scal


         else


            scalp = -1.0


         end if




      end if


*-----------------------------------------------------------------------


         if( ( isef .eq. 0 .or. isef .eq. 4 ) .and.
     &         itim .eq. 0 ) then


            write(io,'(/''newpage:'')')


         else if( isef .eq. 2 .and.
     &            imult / mul * mul .eq. imult .and.
     &            itim .eq. 0 ) then


            write(io,'(/''newpage:''/)')


         else if( isef .eq. 5 .and.
     &            imult / ( mul / 2 ) * ( mul / 2 ) .eq. imult .and.
     &            itim .eq. 0 ) then


            write(io,'(/''newpage:''/)')


         end if


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine disp05(io,minx,maxx,miny,maxy,rho,high)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to write 2-d matrix on the file                         *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              io          : io of the out put file                    *
*              minx, maxx                                              *
*              miny, maxy  : min. and max of the field                 *
*              rho         : 2-d data                                  *
*              high        : maximum value of the data                 *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


      character       zvalu(mnrho)*3
      dimension       rho(mnrho,mnrho)


*-----------------------------------------------------------------------


         do 210 j = maxy, miny, -1


            do 110 i = minx, maxx


                  help = rho(i,j) * 99.9999 / high


               if( help .ge. 1.0 ) then


                  write( zvalu(i), '(i3)') int( help )


               else if( help .ge. 0.9 ) then


                  zvalu(i) = ' .9'


               else if( help .ge. 0.8 ) then


                  zvalu(i) = ' .8'


               else if( help .ge. 0.7 ) then


                  zvalu(i) = ' .7'


               else if( help .ge. 0.6 ) then


                  zvalu(i) = ' .6'


               else if( help .ge. 0.5 ) then


                  zvalu(i) = ' .5'


               else if( help .ge. 0.4 ) then


                  zvalu(i) = ' .4'


               else if( help .ge. 0.3 ) then


                  zvalu(i) = ' .3'


               else if( help .ge. 0.2 ) then


                  zvalu(i) = ' .2'


               else if( help .ge. 0.1 ) then


                  zvalu(i) = ' .1'


               else


                  zvalu(i) = '  0'


               end if


  110       continue


                  write(io,'(60a3)') ( zvalu(i), i = minx, maxx )


  210    continue


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine disp06(ptime)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to display particle position on terminal                *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              ptime                                                   *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


      parameter     ( isx = 18 , isy = 10 )


*-----------------------------------------------------------------------


      common /vriab0/ massal, massba, nmeson


      common /vriab1/ b, llnow, ntnow


      common /rannum/ iseed, iseed0, iseed1


      common /coodrp/ r(3,nnn),  p(5,nnn)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)


      common /const0/ idnta, idnpr, massta, masspr, mstapr, msprpr


      common /const5/ pzpr, pxpr, rzpr, rxpr
      common /const6/ pzta, pxta, rzta, rxta
      common /const7/ betpr, gampr, prmas, radpr
      common /const8/ betta, gamta, tamas, radta


      common /clusti/ itc(0:nnn,0:nnn)
      common /clustf/ nclst, iclust(nnn)
      common /clustg/ jclust(0:7,nnn),  qclust(0:7,nnn)
      common /framtr/ betafr(0:2), gammfr(0:2)


      common /summ00/ jdsp(50), idsp(50), fdsp(50), ifds(50)
      character    fdsp*80


      common /summ01/ iefw, ised, mulg, icdp(4), scald, isef


*-----------------------------------------------------------------------


      dimension ncount(-isx:isx,-isy:isy)


*-----------------------------------------------------------------------


      save ifirst
      save etot0


*-----------------------------------------------------------------------
*           Some constants
*-----------------------------------------------------------------------


                  forms = float( isy ) / float( isx )


                  pmag = 2.5
                  pfm  = ( 3.0 / 2.0 * pi**2 * rho0 )**(1./3.)


                  pminx0 = pzta / hbc - pmag * pfm * sqrt(gamta)


               if( masspr .gt. 1 ) then


                  pmaxx0 = pzpr / hbc + pmag * pfm * sqrt(gampr)


               else


                  pmaxx0 = pzpr / hbc + pmag * pfm


               end if


                  pminx0 = real( nint( pminx0 ) ) * scald
                  pmaxx0 = real( nint( pmaxx0 ) ) * scald


               if( ised .eq. 3 ) then


                  pminx = - ( pmaxx0 - pminx0 ) / 2.0
                  pmaxx =   ( pmaxx0 - pminx0 ) / 2.0
                  pminy = - ( pmaxx0 - pminx0 ) / 2.0
                  pmaxy =   ( pmaxx0 - pminx0 ) / 2.0


               else


                  pminx =   pminx0
                  pmaxx =   pmaxx0
                  pminy = - forms * ( pmaxx0 - pminx0 ) / 2.0
                  pmaxy =   forms * ( pmaxx0 - pminx0 ) / 2.0


               end if




*-----------------------------------------------------------------------
*     write Reaction on 6 and idf(4)
*-----------------------------------------------------------------------


      do 5000 jj = 1, 2


*-----------------------------------------------------------------------
*     File unit of output
*-----------------------------------------------------------------------


         if( jj .eq. 1 ) then


            if( jdsp(25) .eq. 0 ) goto 5000


               io = 6


               idsrp = jdsp(25)


         else


            if( jdsp(26) .eq. 0 ) goto 5000


               io = idsp(20)


               idsrp = jdsp(26)


         end if


*-----------------------------------------------------------------------


         if( ptime .le. 0.0001 ) then


            write(io,'(//'' *'',69(''-''),''*'')')


            write(io,'(6x,''Event ='',i3,
     &               ''th   b = '',f7.3,'' (fm),   '',
     &               ''iseed = '',i12)')
     &                 llnow, b, iseed1


            write(io,'('' *'',69(''-''),''*'')')


            ifirst = 0


         end if


*-----------------------------------------------------------------------
*     Display
*-----------------------------------------------------------------------


         do i = -isx, isx
         do j = -isy, isy


            ncount(i,j)=0


         end do
         end do


            if( ised .eq. 1 ) then


               ix = 3
               iy = 1


            else if( ised .eq. 2 ) then


               ix = 3
               iy = 2


            else if( ised .eq. 3 ) then


               ix = 1
               iy = 2


            end if


            m1 = 1
            m2 = massal




         do 100 ii = m1, m2




            if( idsrp .eq. 1 ) then


              i = nint( r(ix,ii) * 0.72 / scald )
              j = nint( r(iy,ii) * 0.72 / scald )


            else


              i = nint( ( p(ix,ii) / hbc - pminx ) / ( pmaxx - pminx )
     &          * 2.0 * float(isx) - float(isx) )


              j = nint( ( p(iy,ii) / hbc - pminy ) / ( pmaxy - pminy )
     &          * 2.0 * float(isy) - float(isy) )


            end if




            if( abs(i) .le. isx .and.
     &          abs(j) .le. isy )    ncount(i,j) = ncount(i,j) + 1




 100     continue


            write(io,'()')
            write(io,'('' '',37i2.0)')
     &        ( ( ncount(i,j), i = -isx, isx ), j = isy, -isy, -1 )
            write(io,'()')




*-----------------------------------------------------------------------
*        Total energy and particle number
*-----------------------------------------------------------------------


                  nncl = 0
                  nres = 0
                  nmes = 0


                  call epotall(epot)


                  ekin = 0.0


            do i = 1, massal


                  ekin = ekin + p(4,i)


                  indd = inds(i)


               if( indd .eq. 1 ) then


                  nncl = nncl + 1


               else if( indd .eq. 2 .or. indd .eq. 3 ) then


                  nres = nres + 1


               else if( indd .eq. 4 ) then


                  nmes = nmes + 1


               end if


            end do


                  etot = ( epot + ekin ) / massba


*-----------------------------------------------------------------------
*     write information
*-----------------------------------------------------------------------


            if( ifirst .eq. 0 ) then


               etot0 = etot
               ifirst = 1


            end if


               ene0 = ( etot - etot0 ) * 1000.0


            write(io,'(/'' ** T ='',f7.2,'' (fm/c), Edif ='',f7.4,
     &               '' (MeV/A), [ N, R, P ] ('',i3,'','',i3,
     &               '',''i3,'')''/)')
     &                 ptime, ene0, nncl, nres, nmes




*-----------------------------------------------------------------------
*     summary of clusters
*-----------------------------------------------------------------------


                  nprt = 0
                  nnut = 0
                  nres = 0
                  ncls = 0


         do 1000 i = 1, nclst


                  mclst = itc(i,0)


                  nprot = jclust(1,i)
                  nnuet = jclust(2,i)
                  nreso = jclust(3,i) + jclust(3,i)
                  excit = qclust(6,i) / float(mclst)


            if( iclust(i) .eq. 0  ) then


               write(io,'(''  '',I3,'' ['',i3,'' ]''
     &                    '' ('',i3,'','',i3,'','',i3,'' ) {'',
     &                     f8.3, '' }'')')
     &               i, mclst, nprot, nnuet, nreso, excit


                  ncls = ncls + 1


            else


                  nprt = nprt + nprot
                  nnut = nnut + nnuet
                  nres = nres + nreso


            end if


 1000    continue


                  nfre = nprt + nnut + nres + nmes


               write(io,'(''  '',i3,
     &                 '' [  1 ] ('',i3,'','',i3,
     &                 '','',i3,'','',i3,'' ) Free Particles'')')
     &                 nfre, nprt, nnut, nres, nmes




*-----------------------------------------------------------------------


 5000 continue


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine disp07(pdim,scalo,scalp,isef,
     &                  time,itim,itdg,iefw)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to write all particle momentum position and             *
*              Fermi sphere for ground state properties                *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              pdim        : size of display, -pdim to +pdim (1/fm)    *
*              scalo       : global scale                              *
*              scalp       : scale of the previous graph               *
*              isef        : =0; only momentum, =1; with position      *
*              time        : time (fm/c)                               *
*              itim        : =0; normal, =1; ending                    *
*              itdg        : =0; time is integer, =1; real             *
*              iefw        : =0; write data in one file                *
*                            =1; in each different file                *
*                                                                      *
*              xfx         : fixed x-position in a unit of x-axis      *
*              yfx         : fixed y-position in a unit of y-axis      *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


      common /rannum/ iseed, iseed0, iseed1


      common /vriab0/ massal, massba, nmeson
      common /vriab1/ b, llnow, ntnow


*-----------------------------------------------------------------------


      common /summ00/ jdsp(50), idsp(50), fdsp(50), ifds(50)
      character    fdsp*80


      common /fnamid/ fname(mxfnm), idf(mxfnm), ifnl(mxfnm)
      character       fname*80


*-----------------------------------------------------------------------


      dimension       icdp(10)


      character       ofnm*80


      data idp02 / 0 /


      save idp02


*-----------------------------------------------------------------------
*        multi events
*-----------------------------------------------------------------------


         if( time .lt. 0.0001 ) then


            ievnt = 1


         else


            ievnt = 0


         end if


*-----------------------------------------------------------------------
*        Out put file name and unit
*-----------------------------------------------------------------------


            io = idsp(2)


         if( iefw .eq. 1 ) then


            iq = idsp(2) + 1


            ofnm = fname(4)(1:ifnl(4))//'-gp'
            ifnm = ifnl(4) + 9


            ofnm(ifnm-3:ifnm) = '.ang'


         else


            iq = io


         end if


*-----------------------------------------------------------------------
*        some constants
*-----------------------------------------------------------------------


               afa0 = 1.0


               sca0 = 0.45 * scalo


               scal = 3.0 / pdim
               form = 1.0


               txps = 0.0
               typs = pdim + pdim * 0.2


*-----------------------------------------------------------------------
*        ANGEL parameters
*-----------------------------------------------------------------------


         if( isef .eq. 0 ) then


               xoff = -0.11
               yoff = -1.36


               xfx  =  0.5
               yfx  =  2.3


               xorg =   4.0 / 14.0 * ( 1.0 / sca0 - scal )
     &                + xfx * ( 1.0 / sca0 - 1.0 )
     &                + 1.0 / sca0 * xoff


               yorg = ( 3.5 / 14.0 * ( 1.0 / sca0 - scal )
     &                + yfx * ( 1.0 / sca0 - form )
     &                + 1.0 / sca0 * yoff ) / form


         else if( isef .eq. 1 ) then


               sca0 = 1.0 / scalp


               xorg =  1.4
               yorg = -2.6


         end if


*-----------------------------------------------------------------------
*        default values
*-----------------------------------------------------------------------


               sca00 =  1.00 * scalo


               xoff0 =  0.00
               yoff0 =  0.20


               xfx0  =  0.5
               yfx0  =  0.0


               xorg0 =   4.0 / 14.0 * ( 1.0 / sca00 - scal )
     &                 + xfx0 * ( 1.0 / sca00 - 1.0 )
     &                 + 1.0 / sca00 * xoff0


               yorg0 = ( 3.5 / 14.0 * ( 1.0 / sca00 - scal )
     &                + yfx0 * ( 1.0 / sca00 - form )
     &                + 1.0 / sca00 * yoff0 ) / form


*-----------------------------------------------------------------------
*     Write data on the file for ANGEL
*-----------------------------------------------------------------------


      if( isef .eq. 0 ) then


         if( ievnt .eq. 1 .and. llnow .gt. 1 ) then


            write(io,'( /''newpage:'')')


         end if


            write(io,'(/
     &          ''*'',71(''-'') /
     &          ''*     Snapshot of P-Space'',
     &               '' in Three Directions at T ='',f7.2,'' (fm/c)''/
     &          ''*     '',i2,'' th Event with iseed = '',i12/
     &          ''*'',71(''-''))') time, llnow, iseed1


            write(io,'(/
     &          ''msur:{\rm\it Snapshot of P-Space'',
     &          '' in Three Directions of \ ''/ i2,
     &          '' th Event'',i12,''}'')') llnow, iseed1


            write(io,'(
     &          ''msul:{\rm\it JQMD  plotted by  \ANGEL}''/
     &          ''msdc:{\rm\Large \- \page  \-}'')')


      end if


*-----------------------------------------------------------------------
*     Three figures of z-x, z-y, x-y
*-----------------------------------------------------------------------


      do kk = 1, 3


         if( kk .gt. 1 .or. isef .ne. 0 ) then


            write(io,'(/''z:'')')


         end if


         if( iefw .eq. 1 ) then


            idp02 = idp02 + 1


            write(ofnm(ifnm-5:ifnm-4),'(i2.2)') idp02


            open(iq, file=ofnm, status='unknown')


            write(io,'( /''infl: {'',80a1)')
     &            ( ofnm(k:k), k = 1, ifnm ), '}'


         end if


*-----------------------------------------------------------------------


         if( kk .eq. 1 ) then


               kkx = 3
               kky = 1
               kkz = 2


            write(iq,'( /''*'',71(''-'') /
     &                   ''*     Snapshot of P-Space'',
     &                   '' in p_z-p_x plane at T ='',f7.2,
     &                   '' (fm/c)''/
     &                   ''*     '',i2,'' th Event with iseed = '',i12/
     &                   ''*'',71(''-'') )') time, llnow, iseed1


            write(iq,'(/''x: p_z (1/fm)'')')
            write(iq,'( ''y: p_x (1/fm)'')')


*-----------------------------------------------------------------------


         else if( kk .eq. 2 ) then


               kkx = 3
               kky = 2
               kkz = 1


            write(iq,'( /''*'',71(''-'') /
     &                   ''*     Snapshot of P-Space'',
     &                   '' in p_z-p_y plane at T = '',f7.2,
     &                   '' (fm/c)''/
     &                   ''*     '',i2,'' th Event with iseed = '',i12/
     &                   ''*'',71(''-'') )') time, llnow, iseed1


            write(iq,'(/''x: p_z (1/fm)'')')
            write(iq,'( ''y: p_y (1/fm)'')')


*-----------------------------------------------------------------------


         else if( kk .eq. 3 ) then


               kkx = 1
               kky = 2
               kkz = 3


            write(iq,'( /''*'',71(''-'') /
     &                   ''*     Snapshot of R-Space'',
     &                   '' in p_x-p_y plane at T = '',f7.2,
     &                   '' (fm/c)''/
     &                   ''*     '',i2,'' th Event with iseed1'',i12/
     &                   ''*'',71(''-'') )') time, llnow, iseed1


            write(iq,'(/''x: p_x (1/fm)'')')
            write(iq,'( ''y: p_y (1/fm)'')')


         end if


*-----------------------------------------------------------------------


            if( itdg .eq. 0 ) then


               write(iq,'( /''w: T = '',i3,''/ x('',f7.2,
     &                      '') y('',f7.2,'') ix(2) iy(1) s('',f7.4,
     &                      '')'')')
     &                   nint(time), txps, typs, afa0/scal*1.4


            else


               write(iq,'( /''w: T = '',f6.1,''/ x('',f7.2,
     &                      '') y('',f7.2,'') ix(2) iy(1) s('',f7.4,
     &                      '')'')')
     &                   time, txps, typs, afa0/scal*1.4


            end if


            write(iq,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       xorg0, yorg0


            write(iq,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                   '') xfac('',f7.4,'')'')')
     &                       scal * sca00, afa0 / scal, 1.0 / scal


            write(iq,'(  ''p: nofr noms nosp port form(1.0)'')')


            write(iq,
     &               '(''p: xmin('',f6.1,'') xmax('',f6.1,'') ymin('',
     &                 f6.1,'') ymax('',f6.1,'')'')')
     &                 -pdim, pdim, -pdim, pdim


            pfm = ( 3.0 / 2.0 * pi**2 * rho0 )**(1./3.)


            write(iq,'(/''c: Plot of Fermi Sphere''/)')
            write(iq,'(''set: c1['',f7.4,'']'')') pfm
            write(iq,'(''p: nocn'')')
            write(iq,
     &          '(''h: v=[0,2*pi,40] '',
     &            ''x=[c1*cos(v)] y=[c1*sin(v)],l0b'')')


*-----------------------------------------------------------------------
*        write particle data on file
*-----------------------------------------------------------------------


                  icdp(1) = 1


               do i = 2, 4


                  icdp(i) = 0


               end do


               call disp03(iq,kkx,kky,kkz,icdp,1)


*-----------------------------------------------------------------------


         if( iefw .eq. 1 ) then


            close(iq)


         end if


         if( io .eq. iq ) then


            write(io,'( /''*'',71(''-'')/
     &                   ''      Parameters for Multipage of Angel''/
     &                   ''*'',71(''-''))')


         end if


*-----------------------------------------------------------------------


         if( kk .eq. 1 ) then


            write(io,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                        xorg, yorg


            write(io,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                   '') xfac('',f7.4,'')'')')
     &                       scal * sca0, afa0 / scal, 1.0 / scal


         else


            write(io,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       0.0, 1.3


            write(io,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                   '') xfac('',f7.4,'')'')')
     &                       1.0, afa0 / scal, 1.0 / scal


         end if


            if( kk .eq. 3 ) then


               write(io,'(  ''p: mssg'')')


            else


               write(io,'(  ''p: mssg nocm'')')


            end if


*-----------------------------------------------------------------------


      end do


*-----------------------------------------------------------------------


         if( itim .eq. 0 ) then


            write(io,'(/''newpage:'')')


         end if


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine disp08(qminx,qmaxx,qminy,qmaxy,
     &                  scalo,scalp,formp,ised,isef,mulg,
     &                  time,itim,itdg,iefw)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to write all particle momentum position and             *
*              Fermi sphere for reactions                              *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              qminx,qmaxx : size of display (1/fm)                    *
*              qminy,qmaxy   if qmin .ge. 100.0,                       *
*                            size is default values                    *
*              scalo       : global scale                              *
*              scalp       : scale of the previous graph               *
*              formp       : form of the previous graph                *
*              ised        : =1; z-x, =2; z-y, =3; x-y                 *
*              isef        : =0; only momentum, =1; with position      *
*                            =2; multi-graph in one page               *
*                            =3; multi-graph with position in one page *
*              mulg        : number of graphs in a page                *
*                            =< 0; default                             *
*              time        : time (fm/c)                               *
*              itim        : =0; normal, =1; ending                    *
*              itdg        : =0; time is integer, =1; real             *
*              iefw        : =0; write data in one file                *
*                            =1; in each different file                *
*                                                                      *
*              xfx         : fixed x-position in a unit of x-axis      *
*              yfx         : fixed y-position in a unit of y-axis      *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


      common /rannum/ iseed, iseed0, iseed1


      common /vriab0/ massal, massba, nmeson
      common /vriab1/ b, llnow, ntnow


      common /const0/ idnta, idnpr, massta, masspr, mstapr, msprpr


      common /const5/ pzpr, pxpr, rzpr, rxpr
      common /const6/ pzta, pxta, rzta, rxta
      common /const7/ betpr, gampr, prmas, radpr
      common /const8/ betta, gamta, tamas, radta


      common /summ00/ jdsp(50), idsp(50), fdsp(50), ifds(50)
      character    fdsp*80


      common /fnamid/ fname(mxfnm), idf(mxfnm), ifnl(mxfnm)
      character       fname*80


*-----------------------------------------------------------------------


      dimension       icdp(10)


*-----------------------------------------------------------------------


      character       ofnm*80


      data idp06 / 0 /


      save idp06


      save imult


*-----------------------------------------------------------------------
*        multi events
*-----------------------------------------------------------------------


         if( time .lt. 0.0001 ) then


            ievnt = 1
            imult = 0


         else


            ievnt = 0


         end if


*-----------------------------------------------------------------------
*        Out put file name and unit
*-----------------------------------------------------------------------


            io = idsp(6)


         if( iefw .eq. 1 ) then


            iq = idsp(6) + 1


            ofnm = fname(4)(1:ifnl(4))//'-p0'
            ifnm = ifnl(4) + 9


            ofnm(ifnm-3:ifnm) = '.ang'


         else


            iq = io


         end if


*-----------------------------------------------------------------------
*           Some constants
*-----------------------------------------------------------------------


               pmag = 2.5
               pfm  = ( 3.0 / 2.0 * pi**2 * rho0 )**(1./3.)


               mul  = mulg
               afa0 = 0.8


*-----------------------------------------------------------------------
*        size of dispay; form and scal
*-----------------------------------------------------------------------


            if( qminx .ge. 100.0 .or.
     &          qmaxx .ge. 100.0 .or.
     &          qminy .ge. 100.0 .or.
     &          qmaxy .ge. 100.0 ) then


                  pminx0 = pzta / hbc - pmag * pfm * sqrt(gamta)


               if( masspr .gt. 1 ) then


                  pmaxx0 = pzpr / hbc + pmag * pfm * sqrt(gampr)


               else


                  pmaxx0 = pzpr / hbc + pmag * pfm


               end if


                  pminx0 = real( nint( pminx0 ) )
                  pmaxx0 = real( nint( pmaxx0 ) )


            end if


            if( qminx .ge. 100.0 ) then


               if( ised .eq. 3 ) then


                  pminx = - ( pmaxx0 - pminx0 ) / 2.0


               else


                  pminx = pminx0


               end if


            else


                  pminx = qminx


            end if


            if( qmaxx .ge. 100.0 ) then


               if( ised .eq. 3 ) then


                  pmaxx = ( pmaxx0 - pminx0 ) / 2.0


               else


                  pmaxx = pmaxx0


               end if


            else


                  pmaxx = qmaxx


            end if


            if( qminy .ge. 100.0 ) then


               if( ised .eq. 3 ) then


                  pminy = - ( pmaxx0 - pminx0 ) / 2.0


               else


                  pminy = - 0.75 * ( pmaxx0 - pminx0 ) / 2.0


               end if


            else


                  pminy = qminy


            end if


            if( qmaxy .ge. 100.0 ) then


               if( ised .eq. 3 ) then


                  pmaxy = ( pmaxx0 - pminx0 ) / 2.0


               else


                  pmaxy = 0.75 * ( pmaxx0 - pminx0 ) / 2.0


               end if


            else


                  pmaxy = qmaxy


            end if


            if( isef .eq. 3 ) then


                  pminy = - formp * ( pmaxx - pminx ) / 2.0
                  pmaxy =   formp * ( pmaxx - pminx ) / 2.0


            end if




         if( ised .eq. 1 ) then


                  pctax = pzta / hbc
                  pctay = pxta / hbc
                  pcprx = pzpr / hbc
                  pcpry = pxpr / hbc


                  pftax = pfm * gamta
                  pftay = pfm
                  pfprx = pfm * gampr
                  pfpry = pfm


         else if( ised .eq. 2 ) then


                  pctax = pzta / hbc
                  pctay = 0.0
                  pcprx = pzpr / hbc
                  pcpry = 0.0


                  pftax = pfm * gamta
                  pftay = pfm
                  pfprx = pfm * gampr
                  pfpry = pfm


         else if( ised .eq. 3 ) then


                  pctax = 0.0
                  pctay = pxta / hbc
                  pcprx = 0.0
                  pcpry = pxpr / hbc


                  pftax = pfm
                  pftay = pfm
                  pfprx = pfm
                  pfpry = pfm


         end if


                  pdimx = pmaxx - pminx
                  pdimy = pmaxy - pminy




                  form = pdimy / pdimx
                  scal = 6.0 / pdimx


*-----------------------------------------------------------------------
*        offset and comment position
*-----------------------------------------------------------------------


            if( mul .le. 0 ) then


               mul = max( 2, int( 3.0 / form + 0.01 ) * 2 )


            end if


            if( isef .eq. 0 ) then


               sca0 = 1.0 * scalo


               xoff = 0.0
               yoff = 0.3


               xfx  = 0.5
               yfx  = 0.0


            else if( isef .eq. 1 ) then


               sca0 = 1.0 / scalp


               xoff = 0.0
               yoff = 1.0 + 0.3 / formp


            else if( isef .eq. 2 ) then


               sca0 = 0.52 * scalo


               xoff = -0.64
               yoff =  0.57


               xfx  =  1.0
               yfx  =  1.0


               xds0 =  0.0
               xds1 =  1.3
               yds1 = -1.1


               yofg =  1.1 * ( mul / 2 - 1 )


               imult = imult + 1


            else if( isef .eq. 3 ) then


               sca0 = 1.0 / scalp


               xds1 =  1.3
               yds1 =  0.0


               imult = imult + 1


            end if


               xorg =   4.0 / 14.0 * ( 1.0 / sca0 - scal )
     &                + xfx * ( 1.0 / sca0 - 1.0 )
     &                + 1.0 / sca0 * xoff


               yorg = ( 3.5 / 14.0 * ( 1.0 / sca0 - scal )
     &                + yfx * ( 1.0 / sca0 - form )
     &                + 1.0 / sca0 * yoff ) / form


               txps = pminx + ( pmaxx - pminx ) * 0.05
               typs = pmaxy - ( pmaxy - pminy ) * 0.05 / form


*-----------------------------------------------------------------------
*        default values
*-----------------------------------------------------------------------


               sca00 = 1.0 * scalo


               xoff0 = 0.0
               yoff0 = 0.3


               xfx0  = 0.5
               yfx0  = 0.0


               xorg0 =   4.0 / 14.0 * ( 1.0 / sca00 - scal )
     &                 + xfx0 * ( 1.0 / sca00 - 1.0 )
     &                 + 1.0 / sca00 * xoff0


               yorg0 = ( 3.5 / 14.0 * ( 1.0 / sca00 - scal )
     &                 + yfx0 * ( 1.0 / sca00 - form )
     &                 + 1.0 / sca00 * yoff0 ) / form


*-----------------------------------------------------------------------
*        write ANGEL parameters
*-----------------------------------------------------------------------


      if( ievnt .eq. 1 .and. llnow .gt. 1 .and.
     &  ( isef .eq. 0 .or. isef .eq. 2 ) ) then


            write(io,'( /''newpage:'')')


      end if


      if( isef .eq. 0 ) then


            write(io,'(/
     &          ''*'',71(''-'') /
     &          ''*     Snapshot of P-Space'',
     &               '' at T ='',f7.2,'' (fm/c)''/
     &          ''*     '',i2,'' th Event with iseed = '',i12/
     &          ''*'',71(''-''))') time, llnow, iseed1


            write(io,'(/
     &          ''msur:{\rm\it Snapshot of P-Space of \ ''/ i2,
     &          '' th Event'',i12,''}''/
     &          ''msul:{\rm\it JQMD  plotted by  \ANGEL}''/
     &          ''msdc:{\rm\Large \- \page  \-}'')')
     &          llnow, iseed1


      else if( isef .eq. 1 ) then


            write(io,'( /''z:'')')


      else if( isef .eq. 2 ) then


         if( ( imult - 1 ) / mul * mul .eq. imult - 1 ) then


            write(io,'(/
     &          ''*'',71(''-'') /
     &          ''*     Time Evolution of P-Space''/
     &          ''*     '',i2,'' th Event with iseed = '',i12/
     &          ''*'',71(''-''))') llnow, iseed1


            write(io,'(/
     &          ''msur:{\rm\it Time Evolution of P-Space of \ ''/ i2,
     &          '' th Event'',i12,''}''/
     &          ''msul:{\rm\it JQMD  plotted by  \ANGEL}''/
     &          ''msdc:{\rm\Large \- \page  \-}'')')
     &          llnow, iseed1


         else


            write(io,'( /''z:'')')


         end if


      else if( isef .eq. 3 ) then


            write(io,'( /''z:'')')


      end if


*-----------------------------------------------------------------------


         if( iefw .eq. 1 ) then


            idp06 = idp06 + 1


            write(ofnm(ifnm-5:ifnm-4),'(i2.2)') idp06


            open( iq, file = ofnm, status='unknown')


            write(io,'( /''infl: {'',80a1)')
     &            ( ofnm(k:k), k = 1, ifnm ), '}'


         end if


*-----------------------------------------------------------------------


            write(iq,'(/
     &          ''*'',71(''-'') /
     &          ''*     Snapshot of P-Space'',
     &               '' at T ='',f7.2,'' (fm/c)''/
     &          ''*     '',i2,'' th Event with iseed = '',i12/
     &          ''*'',71(''-''))') time, llnow, iseed1


*-----------------------------------------------------------------------
*        X and Y axis
*-----------------------------------------------------------------------


            if( ised .eq. 1 ) then


               kkx = 3
               kky = 1
               kkz = 2


               write(iq,'(/''x: p_z (1/fm)'')')
               write(iq,'( ''y: p_x (1/fm)'')')


            else if( ised .eq. 2 ) then


               kkx = 3
               kky = 2
               kkz = 1


               write(iq,'(/''x: p_z (1/fm)'')')
               write(iq,'( ''y: p_y (1/fm)'')')


            else if( ised .eq. 3 ) then


               kkx = 1
               kky = 2
               kkz = 3


               write(iq,'(/''x: p_x (1/fm)'')')
               write(iq,'( ''y: p_y (1/fm)'')')


            end if


*-----------------------------------------------------------------------
*        Time and Frame
*-----------------------------------------------------------------------


            if( itdg .eq. 0 ) then


               write(iq,'( /''w: T = '',i3,''/ x('',f7.2,
     &                      '') y('',f7.2,'') ix(1) iy(3) s('',f7.4,
     &                      '')'')')
     &                   nint(time), txps, typs, afa0 / scal * 1.4


            else


               write(iq,'( /''w: T = '',f6.1,''/ x('',f7.2,
     &                      '') y('',f7.2,'') ix(1) iy(3) s('',f7.4,
     &                      '')'')')
     &                   time, txps, typs, afa0 / scal * 1.4


            end if




               write(iq,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       xorg0, yorg0


               write(iq,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                      '') xfac('',f7.4,'')'')')
     &                       scal * sca00, afa0 / scal, 1.0 / scal


               write(iq,'(  ''p: nofr noms nosp port form('',
     &                   f7.4,'')'')') form


               write(iq,
     &             '(''p: xmin('',f6.1,'') xmax('',f6.1,'') ymin('',
     &               f6.1,'') ymax('',f6.1,'')'')')
     &               pminx, pmaxx, pminy, pmaxy


*-----------------------------------------------------------------------
*        Plot of Fermi Spheres
*-----------------------------------------------------------------------


               write(iq,'(/''c: Plot of Fermi Spheres'')')
               write(iq,'(''p: nocn'')')


               write(iq,'(''set: c1['',f7.4,''] c3['',f7.4,'']'')')
     &               pftax, pftay
               write(iq,'(''set: c2['',f7.4,''] c4['',f7.4,'']'')')
     &               pctax, pctay
               write(iq,
     &             '(''h: v=[0,2*pi,40] '',
     &               ''x=[c1*cos(v)+c2] y=[c3*sin(v)+c4],l0b'')')


            if( masspr .gt. 1 ) then


               write(iq,'(''set: c5['',f7.4,''] c7['',f7.4,'']'')')
     &               pfprx, pfpry
               write(iq,'(''set: c6['',f7.4,''] c8['',f7.4,'']'')')
     &               pcprx, pcpry
               write(iq,
     &             '(''h: v=[0,2*pi,40] '',
     &               ''x=[c5*cos(v)+c6] y=[c7*sin(v)+c8],l0b'')')


            end if


*-----------------------------------------------------------------------
*        write particle data on file
*-----------------------------------------------------------------------


               do i = 1, 4


                  icdp(i) = 1


               end do


               call disp03(iq,kkx,kky,kkz,icdp,1)


*-----------------------------------------------------------------------


         if( iefw .eq. 1 ) then


            close(iq)


         end if


         if( io .eq. iq ) then


            write(io,'( /''*'',71(''-'')/
     &                   ''      Parameters for Multipage of Angel''/
     &                   ''*'',71(''-''))')


         end if


*-----------------------------------------------------------------------
*        write ANGEL parameters
*-----------------------------------------------------------------------


      if( isef .eq. 0 ) then


            write(io,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       xorg, yorg


            write(io,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                   '') xfac('',f7.4,'')'')')
     &                       scal * sca0, afa0 / scal, 1.0 / scal


            write(io,'(  ''p: mssg'')')


*-----------------------------------------------------------------------


      else if( isef .eq. 1 ) then


            write(io,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       xoff, yoff


            write(io,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                   '') xfac('',f7.4,'')'')')
     &                       scal * sca0, afa0 / scal, 1.0 / scal


*-----------------------------------------------------------------------


      else if( isef .eq. 2 ) then


         if( ( imult - 1 ) / mul * mul .eq. imult - 1 ) then


            write(io,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       xorg, yorg


            write(io,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                   '') xfac('',f7.4,'')'')')
     &                    scal * sca0, afa0 / scal, 1.0 / scal


            write(io,'(  ''p: mssg'')')


         else if( ( imult - mul / 2 - 1 ) / mul * mul .eq.
     &              imult - mul / 2 - 1 ) then


            write(io,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       xds1, yofg


            write(io,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                   '') xfac('',f7.4,'')'')')
     &                    1.0, afa0 / scal, 1.0 / scal


         else


            write(io,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       xds0, yds1


            write(io,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                   '') xfac('',f7.4,'')'')')
     &                    1.0, afa0 / scal, 1.0 / scal


         end if


         if( imult / ( mul / 2 ) * ( mul / 2 ) .ne. imult .and.
     &       itim .eq. 0 ) then


            write(io,'(  ''p: noxn noxt'')')


         end if


*-----------------------------------------------------------------------


      else if( isef .eq. 3 ) then


            write(io,'( /''p: xorg('',f7.4,'') yorg('',f7.4,'')'')')
     &                       xds1, yds1


            write(io,'(  ''p: scal('',f7.4,'') afac('',f7.4,
     &                   '') xfac('',f7.4,'')'')')
     &                       scal * sca0, afa0 / scal, 1.0 / scal


         if( imult / ( mul / 2 ) * ( mul / 2 ) .ne. imult .and.
     &       itim .eq. 0 ) then


            write(io,'(  ''p: noxn noxt'')')


            scalp = scal


         else


            scalp = -1.0


         end if


      end if


*-----------------------------------------------------------------------


         if( isef .le. 1 .and.
     &       itim .eq. 0 ) then


            write(io,'(/''newpage:'')')


         else if( isef .eq. 2 .and.
     &            imult / mul * mul .eq. imult .and.
     &            itim .eq. 0 ) then


            write(io,'(/''newpage:''/)')


         else if( isef .eq. 3 .and.
     &            imult / ( mul / 2 ) * ( mul / 2 ) .eq. imult .and.
     &            itim .eq. 0 ) then


            write(io,'(/''newpage:''/)')


         end if


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine densrr(it,rho,nx,ny,nz,
     &                  mminx,mmaxx,mminy,mmaxy,imesh)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate spacial density                            *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              it          : id for the particles                      *
*              rho         : output density                            *
*              nx, ny, nz  : x, y axis and integrated axis             *
*              mminx,mmaxx : min and max of x, y axis                  *
*              mminy,mmaxy                                             *
*              imesh       : =1; rmsh = 1.0, =0; rmsh = 0.5            *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


      parameter ( dcut = 5.0 )


*-----------------------------------------------------------------------


      common /vriab0/ massal, massba, nmeson
      common /poten1/ gamm, c0, c3, cs, cl, wl


      common /coodrp/ r(3,nnn),  p(5,nnn)


      dimension       rho(mnrho,mnrho)
      dimension       it(0:nnn)


      dimension       cm(24000)


      save            cm
      save            ip, iq, rmsh


      data ideni /0/


*-----------------------------------------------------------------------
*        Initialization of the weight
*-----------------------------------------------------------------------


         if( ideni .eq. 0 ) then


               ideni = 1


            if( imesh .eq. 1 ) then


               rmsh = 1.0
               rnrm = 1.0


               ip = 2
               iq = 2


            else if( imesh .eq. 0 ) then


               rmsh = 0.5
               rnrm = 8.0


               ip = 1
               iq = 4


            end if


               rd = rmsh / float( 2 * ip + 1 )


               do ix = -ip, ip
               do iy = -ip, ip
               do iz = -ip, ip


                     ic = 1 + ( iz + ip )
     &                      + ( iy + ip ) * ( 2 * ip + 1 )
     &                      + ( ix + ip ) * ( 2 * ip + 1 )**2


                     xi = float(ix) * rd
                     yi = float(iy) * rd
                     zi = float(iz) * rd


                     ib  = 0
                     sek = 0.0


                  do jx = -iq, iq
                  do jy = -iq, iq
                  do jz = -iq, iq


                     ib = ib + 1


                     id = ic + ( ib - 1 ) * ( 2 * ip + 1 )**3


                     xj = float(jx) * rmsh
                     yj = float(jy) * rmsh
                     zj = float(jz) * rmsh


                     rsqr = ( xi - xj )**2
     &                    + ( yi - yj )**2
     &                    + ( zi - zj )**2


                     if( rsqr .gt. dcut ) then


                        cm(id) = 0.0


                     else


                        cm(id) = exp( -rsqr / 2.0 / wl )


                     end if


                     sek = sek + cm(id)


                  end do
                  end do
                  end do


                     ib = 0


                  do jx = -iq, iq
                  do jy = -iq, iq
                  do jz = -iq, iq


                     ib = ib + 1


                     id = ic + ( ib - 1 ) * ( 2 * ip + 1 )**3


                     cm(id) = cm(id) / sek


                  end do
                  end do
                  end do


               end do
               end do
               end do


         end if


*-----------------------------------------------------------------------
*        Zero set of the density
*-----------------------------------------------------------------------


               do ix = 1, mmaxx - mminx + 1
               do iy = 1, mmaxy - mminy + 1


                  rho(ix,iy) = 0.0


               end do
               end do


*-----------------------------------------------------------------------
*        Sum up the density
*-----------------------------------------------------------------------


         do ii = 1, it(0)


               i = it(ii)


               ix = nint( r(nx,i) / rmsh )
               iy = nint( r(ny,i) / rmsh )
               iz = nint( r(nz,i) / rmsh )


               kx = nint( float( 2 * ip + 1 )
     &                  * ( r(nx,i) / rmsh - float(ix) )
     &                  + float( ix * ( ip + 1 ) ) )
     &            - ix * ( ip + 1 )


               ky = nint( float( 2 * ip + 1 )
     &                  * ( r(ny,i) / rmsh - float(iy) )
     &                  + float( iy * ( ip + 1 ) ) )
     &            - iy * ( ip + 1 )


               kz = nint( float( 2 * ip + 1 )
     &                  * ( r(nz,i) / rmsh - float(iz) )
     &                  + float( iz * ( ip + 1 ) ) )
     &            - iz * ( ip + 1 )


               ic = 1 + ( kz + ip )
     &                + ( ky + ip ) * ( 2 * ip + 1 )
     &                + ( kx + ip ) * ( 2 * ip + 1 )**2


            do jx = max(mminx,ix-iq), min(mmaxx,ix+iq)
            do jy = max(mminy,iy-iq), min(mmaxy,iy+iq)


               kx = 1 + jx - mminx
               ky = 1 + jy - mminy


               lx = jx - ix
               ly = jy - iy


            do jz = max(mminx,iz-iq), min(mmaxx,iz+iq)


               lz = jz - iz


               ib = 1 + ( lz + iq )
     &                + ( ly + iq ) * ( 2 * iq + 1 )
     &                + ( lx + iq ) * ( 2 * iq + 1 )**2


               id = ic + ( ib - 1 ) * ( 2 * ip + 1 )**3


               rho(kx,ky) = rho(kx,ky)
     &                       + cm(id) * p(5,i) / rmass * rnrm


            end do
            end do
            end do


         end do


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine cntdrw(jol,nx,ny,dx,dy,dax,day,p,z,
     &                  ipd2,col)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to write contour line                                   *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              jol         : io of the output file                     *
*              nx,ny       : dimension of matrix p( , )                *
*              dx,dy       : mesh width                                *
*              dax,day     : grit point coordinate                     *
*              p( , )      : data on grid                              *
*              z           : potential                                 *
*                                                                      *
*              ipd2        : 1-> spline 0-> no                         *
*              col         : color of line                             *
*                            0-> block, 3-> blue                       *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


      parameter     ( istx = 0, isty = 1, ists = 2 )


*-----------------------------------------------------------------------


      dimension     p(mnrho,mnrho), dax(mnrho), day(mnrho)


*-----------------------------------------------------------------------


      logical trace(mnrho,mnrho)
      logical ods


*-----------------------------------------------------------------------


         ods(i,j,ii,jj) = ( p(i,j)   .lt. z .and.
     &                      p(ii,jj) .gt. z ) .or.
     &                    ( p(i,j)   .gt. z .and.
     &                      p(ii,jj) .lt. z )


*-----------------------------------------------------------------------


            do i = 1, nx
            do j = 1, ny


               trace(i,j) = .false.


               if( p(i,j) .eq. z )
     &             p(i,j) = p(i,j) + abs( p(i,j) * 1.e-7 )


            end do
            end do


*-----------------------------------------------------------------------


           do 120 ii = 1, nx - 1


                  i = ii
                  j = 1


               if( trace(i,j) ) goto 120


               if( ods(i,j,i+1,j) ) then


                  istate = istx


               else


                  goto 120


               end if


                  call cntsub(p,nx,ny,dx,dy,dax,day,z,istate,30,i,j,
     &                        trace,jol,ipd2,col)


  120       continue


*-----------------------------------------------------------------------


            do 130 jj = 1, ny - 1


                  i = 1
                  j = jj


               if( trace(i,j) ) goto 130


               if( ods(i,j,i,j+1) ) then


                  istate = isty


               else


                  goto 130


               end if


                  call cntsub(p,nx,ny,dx,dy,dax,day,z,istate,30,i,j,
     &                        trace,jol,ipd2,col)


  130       continue


*-----------------------------------------------------------------------


            do 140 ii = 1, nx - 1


                  i = ii
                  j = ny


               if( trace(i,j) ) goto 140


               if( ods(i,j,i+1,j) ) then


                  istate = istx


               else


                  goto 140


               end if


                  call cntsub(p,nx,ny,dx,dy,dax,day,z,istate,50,i,j,
     &                        trace,jol,ipd2,col)


  140       continue


*-----------------------------------------------------------------------


            do 150 jj = 1, ny - 1


                  i = nx
                  j = jj


               if( trace(i,j) ) goto 150


               if( ods(i,j,i,j+1) ) then


                  istate = isty


               else


                  goto 150


               end if


                  call cntsub(p,nx,ny,dx,dy,dax,day,z,istate,50,i,j,
     &                        trace,jol,ipd2,col)


  150       continue


*-----------------------------------------------------------------------


            do 160 ii = 1, nx - 1
            do 160 jj = 1, ny - 1


                  i = ii
                  j = jj


               if( trace(i,j) ) goto 160


               if( ods(i,j,i+1,j) ) then


                  istate = istx


               else if( ods(i,j,i,j+1) ) then


                  istate = isty


               else


                  goto 160


               end if


                  call cntsub(p,nx,ny,dx,dy,dax,day,z,istate,30,i,j,
     &                        trace,jol,ipd2,col)


  160       continue


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine cntsub(p,nx,ny,dx,dy,dax,day,z,istate,ientry,i,j,
     &                  trace,jol,ipd2,col)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to make contour data and write them on file             *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              p( , )      : data on grid                              *
*              nx,ny       : dimension of matrix p( , )                *
*              dx,dy       : mesh width                                *
*              dax,day     : grit point coordinate                     *
*              z           : potential                                 *
*              istate      : = istx, isty, ists                        *
*              ientry      : 30 (boarder) 50 (normal)                  *
*                                                                      *
*              jol         : output file                               *
*              ipd2        : 1-> spline 0-> no                         *
*              col         : color of line                             *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


      parameter     ( istx = 0, isty = 1, ists = 2 )


*-----------------------------------------------------------------------


      dimension       p(mnrho,mnrho), dax(mnrho), day(mnrho)


*-----------------------------------------------------------------------


      logical         trace(mnrho,mnrho)
      logical         ods


*-----------------------------------------------------------------------


         ods(i,j,ii,jj) = ( p(i,j)   .lt. z .and.
     &                      p(ii,jj) .gt. z ) .or.
     &                    ( p(i,j)   .gt. z .and.
     &                      p(ii,jj) .lt. z )


*-----------------------------------------------------------------------


         divx(i,j) = dax(i)
     &             + ( z - p(i,j) ) / ( p(i+1,j) - p(i,j) ) * dx
         divy(i,j) = day(j)
     &             + ( z - p(i,j) ) / ( p(i,j+1) - p(i,j) ) * dy


*-----------------------------------------------------------------------


            ic = 0


            if( ientry .eq. 50 ) goto 50


*-----------------------------------------------------------------------


   30 continue


            if( istate .eq. istx ) then


               call cntout(ic,divx(i,j),day(j),jol,ipd2,col)


               ic = ic + 1


            end if


            if( istate .eq. isty ) then


               call cntout(ic,dax(i),divy(i,j),jol,ipd2,col)


               ic = ic + 1


            end if


            if( trace(i,j) ) goto 40


               trace(i,j) = .true.




            if( i .eq. nx .or. j .eq. ny ) goto 40




            if( istate .eq. istx ) then


               if( ods(i+1,j,i,j+1) ) then


                  istate = ists


               else


                  istate = isty


               end if


            else if( istate .eq. isty ) then


               if( ods(i+1,j,i,j+1) ) then


                  istate = ists


               else


                  istate = istx


               end if


            else


               if( ods(i,j,i+1,j) ) then


                  istate = istx


               else


                  istate = isty


               end if


            end if


*-----------------------------------------------------------------------


   50 continue




            if( istate .eq. istx ) then


               call cntout(ic,divx(i,j),day(j),jol,ipd2,col)


               ic = ic + 1


            else if( istate .eq. isty ) then


               call cntout(ic,dax(i),divy(i,j),jol,ipd2,col)


               ic = ic + 1


            end if


*-----------------------------------------------------------------------


            if( i .eq. 1 .and. istate .eq. isty ) goto 40
            if( j .eq. 1 .and. istate .eq. istx ) goto 40


*-----------------------------------------------------------------------


            if( istate .eq. istx ) then


               if( ods(i,j,i+1,j-1) ) then


                  istate = ists


                  j = j - 1


               else


                  istate = isty


                  i = i + 1
                  j = j - 1


               end if


            else if( istate .eq. isty ) then


               if( ods(i,j,i-1,j+1) ) then


                  istate = ists


                  i = i - 1


               else


                  istate = istx


                  i = i - 1
                  j = j + 1


               end if


            else


               if( ods(i+1,j,i+1,j+1) ) then


                  istate = isty


                  i = i + 1


               else


                  istate = istx


                  j = j + 1


               end if


            end if


                  goto 30


*-----------------------------------------------------------------------


   40 continue


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine cntout(ic,x,y,jol,ipd2,rcol)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to write contour data on file                           *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              ic          :  0  -> start of the line                  *
*                            >0  -> sequential lines                   *
*                            -1  -> end of line                        *
*              x, y        : coordinate of the points                  *
*              jol         : output file                               *
*              ipd2        : 1-> spline 0-> no                         *
*              rcol        : color of line                             *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      if( ic .eq. 0 ) then


         if( ipd2 .eq. 0 ) then


            if( nint(rcol*2) .eq. 2 ) then
               write(jol,'(''h: x y,l0r'')')
            else if( nint(rcol*2) .eq. 3 ) then
               write(jol,'(''h: x y,l0y'')')
            else if( nint(rcol*2) .eq. 4 ) then
               write(jol,'(''h: x y,l0g'')')
            else if( nint(rcol*2) .eq. 5 ) then
               write(jol,'(''h: x y,l0c'')')
            else if( nint(rcol*2) .eq. 6 ) then
               write(jol,'(''h: x y,l0b'')')
            else
               write(jol,'(''h: x y,l0'')')
            end if


         else if( ipd2 .ne. 0 ) then


            if( nint(rcol*2) .eq. 2 ) then
               write(jol,'(''h: x y,l0sr'')')
            else if( nint(rcol*2) .eq. 3 ) then
               write(jol,'(''h: x y,l0sy'')')
            else if( nint(rcol*2) .eq. 4 ) then
               write(jol,'(''h: x y,l0sg'')')
            else if( nint(rcol*2) .eq. 5 ) then
               write(jol,'(''h: x y,l0sc'')')
            else if( nint(rcol*2) .eq. 6 ) then
               write(jol,'(''h: x y,l0sb'')')
            else
               write(jol,'(''h: x y,l0s'')')
            end if


         end if


      end if


*-----------------------------------------------------------------------


               write(jol,'(2e12.4)') x, y


*-----------------------------------------------------------------------


      return
      end




