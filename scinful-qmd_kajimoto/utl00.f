************************************************************************
*                                                                      *
*        PART 10: Utilities                                            *
*                                                                      *
*                                                                      *
*   List of subprograms in rough order of relevance with main purpose  *
*      ( s = subroutine, f = function, b = block data, e = entry )     *
*                                                                      *
*                                                                      *
*  s  jbook1    to initialize one-dimentional histogram                *
*  e  jfill1    to fill one dimensional histogram                      *
*  e  jscale1   to scal one dimentional histgram data                  *
*  e  jprint1   to print one dimensional histogram                     *
*  e  jftot     to give total particle production cross section        *
*  e  jfddx     to give the double differential cross section          *
*               for one angel                                          *
*  s  ranint    to initialize the random number                        *
*  s  jqmdver   to store JQMD version and last reviced date            *
*  s  datetime  to detect date and time                                *
*  s  cputime   to detect cpu time                                     *
*  s  howmany   to control multi-event runs                            *
*  f  pcmsr     to determine the cm momentum from energy and mass      *
*  s  trfram    to determine energy and momentum by Lorentz transform  *
*  f  erf       to calculate the eror function                         *
*                                                                      *
************************************************************************
*                                                                      *
************************************************************************
*                                                                      *
      subroutine jbook1(id,title,
     &                  tfac,ilog,inum,ifac,nx,xmin,xmax,nw,wi,wt)
*                                                                      *
*                                                                      *
*        Last Revised:     1999 03 15                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to initialize one-dimentional histogram                 *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              id       : histogram id =< 100                          *
*              title    : title of the histgram ( character )          *
*                                                                      *
*              tfac     : total factor                                 *
*                                                                      *
*              ilog     : 0=> x-linear bin, 1=> x-log bin              *
*                                                                      *
*              inum     : =0 ; without event number                    *
*                         =1 ; with N event number                     *
*                         =2 ; with 1/sqrt(N)                          *
*                         =3 ; with Y/sqrt(N)                          *
*                                                                      *
*              ifac     : =0 ; without implicit scaling factor         *
*                         =1 ; with scaling factor 1/dx 1/dw           *
*                                                                      *
*              nx       : number of x-bins                             *
*              xmin     : minimum x-value                              *
*              xmax     : maximum x-value                              *
*                                                                      *
*              nw       : number of windows =< 10, or 20               *
*              wi       : edge values of window, wi(20)                *
*              wt       : edge values of window, wt(20)                *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      parameter ( idmax = 90000 )


*-----------------------------------------------------------------------


      common /startf/ iday0,imon0,iyer0,ihor0,imin0,isec0


*-----------------------------------------------------------------------


      dimension vmat(idmax)


      dimension wi(20), wt(20)


      dimension ind(100), ilg(100), inn(100), ifc(100)
      dimension inx(100), inw(100), ixp(100)
      dimension bmi(100), bma(100), bin(100)
      dimension wmi(100,20), wma(100,20)
      dimension gfc(100)


      character title*(*)
      character vtitle*80
      character stitle(100)*80
      dimension ititle(100)


      character br(20)*1


      dimension pxsdd(40)


*-----------------------------------------------------------------------


      data vmat /idmax*0.0/
      save vmat


      data ind / 100*0 /
      save ind, ilg, inn, ifc
      save inw, inx, ixp
      save bmi, bma, bin
      save wmi, wma


      data gfc / 100*1.0 /
      save gfc


      data indx / 0 /
      save indx


      data ixps / 0 /
      save ixps


      save stitle, ititle


      data br /20*' '/


*-----------------------------------------------------------------------
*        check
*-----------------------------------------------------------------------


            if( id .le. 0 .or. id .gt. 100 ) then


               write(*,*) ' **** Error at jbook1: id is out of range,',
     &                    ' 0 < id < 100'
               stop


            end if


            if( ilog .lt. 0 .or. ilog .gt. 1 ) then


               ilog = 0


            end if


            if( inum .lt. 0 .or. inum .gt. 3 ) then


               inum = 0


            end if


            if( ifac .lt. 0 .or. ifac .gt. 1 ) then


               ifac = 0


            end if


            if( xmin .gt. xmax ) then


               write(*,*) ' **** Error at jbook1: xmin and xmax ',
     &                    ' are wrong order'
               stop


            end if


            if( ilog .eq. 1 .and. xmin .le. 0.0 ) then


               write(*,*) ' **** Error at jbook1: x-range is wrong,',
     &                    ' as x-bin is log'
               stop


            end if


            if( nw .lt. 0 ) nw = 0


            if( inum .eq. 0 .and. nw .gt. 20 ) then


               write(*,*) ' **** Error at jbook1: nw should be less ',
     &                    ' than 20, when inum = 0'
               stop


            end if


            if( inum .ne. 0 .and. nw .gt. 10 ) then


               write(*,*) ' **** Error at jbook1: nw should be less ',
     &                    ' than 10, when inum > 0'
               stop


            end if


            if( ind(id) .ne. 0 ) then


               write(*,*) ' **** Error at jbook1: double booking'
               stop


            end if


*-----------------------------------------------------------------------
*     initialization of booking
*-----------------------------------------------------------------------
*        store configuration
*-----------------------------------------------------------------------


               indx = indx + 1


               ind(id)   = indx


               ilg(indx) = ilog
               inn(indx) = inum
               ifc(indx) = ifac
               inx(indx) = nx
               inw(indx) = nw


               ixp(indx) = ixps


               gfc(indx) = gfc(indx) * tfac


*-----------------------------------------------------------------------
*        store bin range
*-----------------------------------------------------------------------


               if( ilog .eq. 0 ) then


                     bmi(indx) = xmin
                     bma(indx) = xmax


               else


                     bmi(indx) = log(xmin)
                     bma(indx) = log(xmax)


               end if


                     bin(indx) = ( bma(indx) - bmi(indx) ) / nx


            do i = 1, nx + 1


               if( ilog .eq. 0 ) then


                     vmat(ixps+i) = xmin + bin(indx) * ( i - 1 )


               else


                     vmat(ixps+i) = xmin * exp( bin(indx) * ( i - 1 ) )


               end if


            end do


*-----------------------------------------------------------------------
*        store window range
*-----------------------------------------------------------------------


            if( nw .gt. 0 ) then


               do i = 1, nw


                  if( wi(i) .lt. wt(i) ) then


                     wmi(indx,i) = wi(i)
                     wma(indx,i) = wt(i)


                  else if( wi(i) .gt. wt(i) ) then


                     wmi(indx,i) = wt(i)
                     wma(indx,i) = wi(i)


                  else if( wi(i) .eq. wt(i) ) then


                     write(*,*) ' **** Error at jbook1:',
     &                          ' window size is zero'
                     stop


                  end if


               end do


            end if


*-----------------------------------------------------------------------
*        store column position
*-----------------------------------------------------------------------


               mx = nx + 1 + 4


            if( inum .eq. 0 ) then


               if( nw .eq. 0 ) then


                     ixps = ixps + 2 * mx


               else


                     ixps = ixps + ( nw + 1 ) * mx


               end if


            else


               if( nw .eq. 0 ) then


                     ixps = ixps + 3 * mx


               else


                     ixps = ixps + ( 2 * nw + 1 ) * mx


               end if


            end if




            if( ixps .gt. idmax) then


               write(*,*) ' **** Error at jbook1: over booking.',
     &                    ' Please increase idmax'
               stop


            end if


*-----------------------------------------------------------------------
*        store title of histgram
*-----------------------------------------------------------------------


         vtitle = title//' '


            do i = 80, 1, -1


               if( vtitle(i:i) .ne. ' ' ) goto 100


            end do


               i = 0


  100       continue


               ititle(indx) = i


            do i = 1, ititle(indx)


               stitle(indx)(i:i) = vtitle(i:i)


            end do


            do i = 1, 80


               vtitle(i:i) = ' '


            end do


*-----------------------------------------------------------------------


      return


*-----------------------------------------------------------------------




************************************************************************
*                                                                      *
      entry jfill1(id,x,v,wn,wf)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to fill one dimensional histogram                       *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              id       : histogram id =< 100                          *
*              x        : x-value                                      *
*              v        : window-value                                 *
*              wn       : number weight                                *
*              wf       : factor                                       *
*                                                                      *
*                         total weight = wn * wf                       *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------
*        check
*-----------------------------------------------------------------------


            ic = ind(id)


            if( ic .eq. 0 ) then


               write(*,*) ' **** Error at jfill1: unbooked id'
               stop


            end if


            if( ilg(ic) .eq. 1 .and. x .le. 0.0 ) return


*-----------------------------------------------------------------------
*        fill matrix
*-----------------------------------------------------------------------


                        ix = ixp(ic)
                        mx = inx(ic) + 1
                        lx = mx + 4




                     if( ilg(ic) .eq. 1 ) then


                        xx = log( x )


                     else


                        xx = x


                     end if


*-----------------------------------------------------------------------


                        vmat(ix+mx+3) = vmat(ix+mx+3) + wn * wf
                        vmat(ix+mx+4) = vmat(ix+mx+4) + wn


                     if( xx .lt. bmi(ic) ) then


                        vmat(ix+mx+1) = vmat(ix+mx+1) + wn * wf


                     else if( xx .ge. bma(ic) ) then


                        vmat(ix+mx+2) = vmat(ix+mx+2) + wn * wf


                     end if


*-----------------------------------------------------------------------


            if( inw(ic) .eq. 0 ) then


                  if( xx .lt. bmi(ic) ) then


                        vmat(ix+lx+mx+1) = vmat(ix+lx+mx+1) + wn


                  else if( xx .ge. bma(ic) ) then


                        vmat(ix+lx+mx+2) = vmat(ix+lx+mx+2) + wn


                  else


                        j = int( ( xx - bmi(ic) ) / bin(ic) ) + 1


                     if( ifc(ic) .eq. 1 ) then


                        fac = 1.0 / ( vmat(ix+j+1) - vmat(ix+j) )


                     else


                        fac = 1.0


                     end if


                        vmat(ix+lx+j) = vmat(ix+lx+j) + wn * wf * fac


                     if( inn(ic) .gt. 0 ) then


                        vmat(ix+2*lx+j) = vmat(ix+2*lx+j) + wn


                     end if


                  end if


*-----------------------------------------------------------------------


            else if( inw(ic) .gt. 0 ) then




               do ii = 1, inw(ic)
               if( v .ge. wmi(ic,ii) .and. v .lt. wma(ic,ii) ) then


                     if( inn(ic) .eq. 0 ) then


                        i = ii


                     else


                        i = 2 * ii - 1


                     end if


                     if( ifc(ic) .eq. 1 ) then


                        fac = 1.0 / ( wma(ic,ii) - wmi(ic,ii) )


                     else


                        fac = 1.0


                     end if


                        vmat(ix+i*lx+mx+3) = vmat(ix+i*lx+mx+3)
     &                                     + wn * wf * fac
                        vmat(ix+i*lx+mx+4) = vmat(ix+i*lx+mx+4)
     &                                     + wn


                  if( xx .lt. bmi(ic) ) then


                        vmat(ix+i*lx+mx+1) = vmat(ix+i*lx+mx+1)
     &                                     + wn * wf * fac


                     if( inn(ic) .gt. 0 ) then


                        vmat(ix+(i+1)*lx+mx+1) = vmat(ix+(i+1)*lx+mx+1)
     &                                         + wn


                     end if


                  else if( xx .ge. bma(ic) ) then


                        vmat(ix+i*lx+mx+2) = vmat(ix+i*lx+mx+2)
     &                                     + wn * wf * fac


                     if( inn(ic) .gt. 0 ) then


                        vmat(ix+(i+1)*lx+mx+2) = vmat(ix+(i+1)*lx+mx+2)
     &                                         + wn


                     end if


                  else


                        j = int( ( xx - bmi(ic) ) / bin(ic) ) + 1


                     if( ifc(ic) .eq. 1 ) then


                        fac = fac / ( vmat(ix+j+1) - vmat(ix+j) )


                     end if


                        vmat(ix+i*lx+j) = vmat(ix+i*lx+j)
     &                                  + wn * wf * fac


                     if( inn(ic) .gt. 0 ) then


                        vmat(ix+(i+1)*lx+j) = vmat(ix+(i+1)*lx+j)
     &                                      + wn


                     end if


                  end if


               end if
               end do


            end if




*-----------------------------------------------------------------------


      return




************************************************************************
*                                                                      *
      entry jscale1(id,scal)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to scal one dimentional histgram data                   *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              id       : histogram id =< 100                          *
*              scal     : global scaling factor                        *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------
*        check
*-----------------------------------------------------------------------


            ic = ind(id)


            if( ic .eq. 0 ) then


               write(*,*) ' **** Error at jscale1: unbooked id'
               return


            end if


*-----------------------------------------------------------------------


            gfc(ic) = gfc(ic) * scal


*-----------------------------------------------------------------------


      return




************************************************************************
*                                                                      *
      entry jprint1(id,io)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to print one dimensional histogram                      *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              id       : histogram id =< 100                          *
*              io       : output unit                                  *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------
*        check
*-----------------------------------------------------------------------


            ic = ind(id)


            if( ic .eq. 0 ) then


               write(*,*) ' **** Error at jprint1: unbooked id'
               return


            end if


*-----------------------------------------------------------------------
*        some constants
*-----------------------------------------------------------------------


                  gfac = gfc(ic)


                  ix = ixp(ic)
                  mx = inx(ic) + 1
                  lx = mx + 4


*-----------------------------------------------------------------------
*        write title
*-----------------------------------------------------------------------


               write(io,'(''#'')')


               write(io,'(''#'',5x,''JBOOK: '',80a1)')
     &              ( stitle(ic)(i:i), i = 1, ititle(ic) )


               write(io,'(''#'')')


               write(io,'(''#'',5x,''Histgram ID :'',i3)') id


               write(io,'(''#'',5x,''No. of entries :'',1pg15.7)')
     &               vmat(ix+mx+4)


               write(io,'(''#'',5x,''Date: '',
     &               i4.4,''-'',i2.2,''-'',i2.2,''   '',
     &               i2.2,'':'',i2.2,'':'',i2.2)')
     &               iyer0,imon0,iday0,ihor0,imin0,isec0


               write(io,'(''#'')')


*-----------------------------------------------------------------------
*        write header
*-----------------------------------------------------------------------


            if( inw(ic) .eq. 0 ) then


               if( inn(ic) .eq. 0 ) then


                     write(io,'(''#'',5x,''x  '',
     &                     9x,''y01'')')


               else


                     write(io,'(''#'',5x,''x  '',
     &                     9x,''y01'',9x,''n01'')')


               end if


            else if( inw(ic) .gt. 0 ) then


               if( inn(ic) .eq. 0 ) then


                     write(io,'(''#'',5x,''x  '',
     &                     20(a1,8x,''y'',i2.2))')
     &                     ( br(j), j, j=1,inw(ic) )


               else


                     write(io,'(''#'',5x,''x  '',
     &                     20(a1,8x,''y''i2.2,9x,''n'',i2.2))')
     &                     ( br(j), j, j, j=1,inw(ic) )


               end if


            end if


               write(io,'(''#'')')


*-----------------------------------------------------------------------
*        print histgrams
*-----------------------------------------------------------------------


         do k = 1, mx


            if( inw(ic) .eq. 0 ) then


               if( inn(ic) .eq. 0 ) then


                     write(io,'(1x,1p2g12.4)')
     &               vmat(ix+k), vmat(ix+lx+k) * gfac


               else


                  if( vmat(ix+2*lx+k) .gt. 0.0 ) then


                     if( inn(ic) .eq. 2 ) then


                        vmat(ix+2*lx+k) = 1.0 / sqrt( vmat(ix+2*lx+k) )


                     else if( inn(ic) .eq. 3 ) then


                        vmat(ix+2*lx+k) = vmat(ix+lx+k) * gfac
     &                                  / sqrt( vmat(ix+2*lx+k) )


                     end if


                  end if


                     write(io,'(1x,1p3g12.4)')
     &               vmat(ix+k), vmat(ix+lx+k) * gfac, vmat(ix+2*lx+k)


               end if


            else if( inw(ic) .gt. 0 ) then


               if( inn(ic) .eq. 0 ) then


                     write(io,'(1x,1p21g12.4)')
     &               vmat(ix+k), ( vmat(ix+j*lx+k) * gfac, j=1,inw(ic) )


               else


                  if( inn(ic) .eq. 2 .or. inn(ic) .eq. 3 ) then


                     do l = 2, inw(ic)*2, 2


                        if( vmat(ix+l*lx+k) .gt. 0.0 ) then


                           if( inn(ic) .eq. 2 ) then


                              vmat(ix+l*lx+k) = 1.0
     &                                        / sqrt( vmat(ix+l*lx+k) )


                           else if( inn(ic) .eq. 3 ) then


                              vmat(ix+l*lx+k) = vmat(ix+(l-1)*lx+k)
     &                                        * gfac
     &                                        / sqrt( vmat(ix+l*lx+k) )


                           end if


                        end if


                     end do


                  end if


                     write(io,'(1x,1p21g12.4)')
     &               vmat(ix+k), ( vmat(ix+(2*j-1)*lx+k) * gfac,
     &                             vmat(ix+2*j*lx+k), j=1,inw(ic) )


               end if


            end if


         end do


*-----------------------------------------------------------------------
*        write summary
*-----------------------------------------------------------------------


                     write(io,'(''#'')')


            if( inw(ic) .eq. 0 ) then


                     write(io,'(''#'',4x,''total  ='',1p1g12.4,
     &                                        ''# ='',1p1g15.7)')
     &                     vmat(ix+mx+3) * gfac, vmat(ix+mx+4)


                     write(io,'(''#'',4x,''under  ='',1p1g12.4,
     &                                        ''# ='',1p1g15.7)')
     &                     vmat(ix+mx+1) * gfac, vmat(ix+lx+mx+1)


                     write(io,'(''#'',4x,''over   ='',1p1g12.4,
     &                                        ''# ='',1p1g15.7)')
     &                     vmat(ix+mx+2) * gfac, vmat(ix+lx+mx+2)


            else if( inw(ic) .gt. 0 ) then


                     write(io,'(''#'',4x,''total  ='',1p1g12.4,
     &                                        ''# ='',1p1g15.7)')
     &                     vmat(ix+mx+3) * gfac, vmat(ix+mx+4)


                     write(io,'(''#'',4x,''under  ='',1p1g12.4)')
     &                     vmat(ix+mx+1) * gfac


                     write(io,'(''#'',4x,''over   ='',1p1g12.4)')
     &                     vmat(ix+mx+2) * gfac


                     write(io,'(''#'')')


               if( inn(ic) .eq. 0 ) then


                     write(io,'(''#'',4x,''win:'',
     &                     20(a1,8x,''y'',i2.2))')
     &                     ( br(j), j, j=1,inw(ic) )


                     write(io,'(''#'')')


                     write(io,'(''#'',4x,''total  ='',1p21g12.4)')
     &                     ( vmat(ix+j*lx+mx+3) * gfac, j=1,inw(ic) )


                     write(io,'(''#'',4x,''    #  ='',1p21g12.4)')
     &                     ( vmat(ix+j*lx+mx+4), j=1,inw(ic) )


                     write(io,'(''#'',4x,''under  ='',1p21g12.4)')
     &                     ( vmat(ix+j*lx+mx+1) * gfac, j=1,inw(ic) )


                     write(io,'(''#'',4x,''over   ='',1p21g12.4)')
     &                     ( vmat(ix+j*lx+mx+2) * gfac, j=1,inw(ic) )


               else


                     write(io,'(''#'',4x,''win:'',
     &                     20(a1,8x,''y''i2.2,9x,''#'',i2.2))')
     &                     ( br(j), j, j, j=1,inw(ic) )


                     write(io,'(''#'')')


                     write(io,'(''#'',4x,''total  ='',1p21g12.4)')
     &                     ( vmat(ix+(2*j-1)*lx+mx+3) * gfac,
     &                       vmat(ix+(2*j-1)*lx+mx+4), j=1,inw(ic) )


                     write(io,'(''#'',4x,''under  ='',1p21g12.4)')
     &                     ( vmat(ix+(2*j-1)*lx+mx+1) * gfac,
     &                       vmat(ix+(2*j  )*lx+mx+1), j=1,inw(ic) )


                     write(io,'(''#'',4x,''over   ='',1p21g12.4)')
     &                     ( vmat(ix+(2*j-1)*lx+mx+2) * gfac,
     &                       vmat(ix+(2*j  )*lx+mx+2), j=1,inw(ic) )


               end if


            end if


               write(io,'(''#'')')




*-----------------------------------------------------------------------


      return




************************************************************************
*                                                                      *
      entry jftot(id,pxst,pxsu,pxso,pfac)
*                                                                      *
*                                                                      *
*        Last Revised:     1999 03 15                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to give total particle production cross section         *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              id       : histogram id =< 100                          *
*              pxst     : total cross section within the range         *
*              pxsu     : total cross section under  the range         *
*              pxso     : total cross section over   the range         *
*              pfac     : normalization factor                         *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------
*        check
*-----------------------------------------------------------------------


            ic = ind(id)


            if( ic .eq. 0 ) then


               write(*,*) ' **** Error at jftot: unbooked id'
               return


            end if


*-----------------------------------------------------------------------
*        some constants
*-----------------------------------------------------------------------


                  gfac = gfc(ic) * pfac


                  ix = ixp(ic)
                  mx = inx(ic) + 1
                  lx = mx + 4


*-----------------------------------------------------------------------
*        write summary
*-----------------------------------------------------------------------


                  pxst = vmat(ix+mx+3) * gfac
                  pxsu = vmat(ix+mx+1) * gfac
                  pxso = vmat(ix+mx+2) * gfac


*-----------------------------------------------------------------------


      return




************************************************************************
*                                                                      *
      entry jfddx(id,jd,pxsdd)
*                                                                      *
*                                                                      *
*        Last Revised:     1999 03 16                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to give the double differential cross section           *
*              for one angel                                           *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              id       : histogram id =< 100                          *
*              jd       : angle id                                     *
*              pxsdd    : cross section                                *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------
*        check
*-----------------------------------------------------------------------


            ic = ind(id)


         if( ic .ne. 0 ) then


*-----------------------------------------------------------------------
*        some constants
*-----------------------------------------------------------------------


                  gfac = gfc(ic)


                  ix = ixp(ic)
                  mx = inx(ic) + 1
                  lx = mx + 4


*-----------------------------------------------------------------------
*        give ddx
*-----------------------------------------------------------------------


               j = jd


               pxsdd(1)   = vmat(ix+j*lx+mx+1) * gfac
     &                    / vmat(ix+1)


            do k = 1, mx - 1


               pxsdd(k+1) = vmat(ix+j*lx+k) * gfac


            end do


               pxsdd(mx+1)  = 0.0


*-----------------------------------------------------------------------


         else


            do k = 1, 40


               pxsdd(k) = 0.0


            end do


         end if


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine ranint
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to initialize the random number                         *
*                                                                      *
*                                                                      *
************************************************************************


      common /rannum/ iseed, iseed0, iseed1
      common /swich2/ icfg, imany, icpus, idatm


*-----------------------------------------------------------------------


            idatm0 = idatm
            idatm  = 1


               call datetime(iyer,imon,iday,ihor,imin,isec)


            idatm  = idatm0




         if( iseed .eq. 0 ) then


            iseed = isec*10000 + ihor*100 + imin


            if(iseed / 2 * 2 .eq. iseed ) iseed = iseed + 1


         end if


            iseed0 = iseed


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine jqmdver( Version, Last Revised )
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to store JQMD version and last reviced date             *
*                                                                      *
*                                                                      *
************************************************************************


      common /verqmd/ versn, lastr, iyeav, imonv, idayv


*-----------------------------------------------------------------------


               versn = Version
               lastr = Last Revised


               iyeav = lastr / 10000
               imonv = ( lastr - iyeav * 10000 ) / 100
               idayv = lastr - iyeav * 10000 - imonv * 100


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine datetime(iyer,imon,iday,ihor,imin,isec)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to detect date and time                                 *
*                                                                      *
*                                                                      *
************************************************************************


      common /swich2/ icfg, imany, icpus, idatm


*-----------------------------------------------------------------------


            if( idatm .ne. 1 ) return


*-----------------------------------------------------------------------


               call date_a_time(iyer,imon,iday,ihor,imin,isec)


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine cputime(ic)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to detect cpu time                                      *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      parameter ( icpum = 30 )


*-----------------------------------------------------------------------


      common /cputim/ stime(30), cputm(30)


      common /swich2/ icfg, imany, icpus, idatm


*-----------------------------------------------------------------------


      dimension estime(icpum), ecputm(icpum)
      dimension ietm(icpum)


      data  ietm  / icpum * 0 /
      data  stime / icpum * 0.0 /
      data  cputm / icpum * 0.0 /


      data  estime / icpum * 0.0 /
      data  ecputm / icpum * 0.0 /


      save ietm
      save estime, ecputm


*-----------------------------------------------------------------------


         if( icpus .ne. 1 ) return


         if( ic .lt. 1 .or. ic .gt. icpum ) then


            write(*,*) ' **** Error at cputime(ic) : ic is too large'
            write(*,*) ' ==========================================='
            stop


         end if


*-----------------------------------------------------------------------


         if( ietm(ic) .eq. 0 ) then


            estime(ic) = sect_a(0.0)
            ecputm(ic) = cput_a(0.0)


            ietm (ic) = 1


         else


            stime(ic) = stime(ic) + sect_a(estime(ic))
            cputm(ic) = cputm(ic) + cput_a(ecputm(ic))


            ietm (ic) = 0


         end if


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine howmany(icd)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to control multi-event runs                             *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              icd         : =0; initialization                        *
*                            =1; each event control                    *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


      common /const2/ dt, ntmax, iprun, iprun0


      common /vriab1/ b, llnow, ntnow


      common /swich2/ icfg, imany, icpus, idatm


      common /rannum/ iseed, iseed0, iseed1


      common /fnamid/ fname(mxfnm), idf(mxfnm), ifnl(mxfnm)
      character       fname*80


*-----------------------------------------------------------------------


         if( imany .eq. 0 ) return


*-----------------------------------------------------------------------
*        Open files for HOWMANY
*-----------------------------------------------------------------------


         if( icd .eq. 0 ) then


               open(idf(3),file=fname(3),status='unknown')


               write(idf(3),'('' 1 <-0: stop'')')
               write(idf(3),'('' event = '',i6,'' / '',i6)')
     &                                      1, iprun0
               write(idf(3),'('' iseed = '',i12)') iseed0


               close(idf(3))


*-----------------------------------------------------------------------
*        Multi Run Control
*-----------------------------------------------------------------------


         else if( icd .eq. 1 ) then


               open (idf(3),file=fname(3),status='old')
               read (idf(3),*) icont
               close(idf(3))


            if( icont .eq. 0 ) then


               open (idf(3),file=fname(3),status='OLD')
               write(idf(3),'('' 2 : next event is'')')
               write(idf(3),'('' event = '',i6,'' / '',i6)')
     &                                      llnow, iprun0
               write(idf(3),'('' iseed = '',i12)') iseed1
               close(idf(3))


               iprun = llnow - 1


               call finsumry


               write(6,'(/'' **** Run is Terminated by HOWMANY ****''/
     &                    ''      at Event = '',i6/)') iprun


               stop


            else if( icont .eq. 1 ) then


               open (idf(3),file=fname(3),status='OLD')
               write(idf(3),'('' 1 <-0: stop'')')
               write(idf(3),'('' event = '',i6,'' / '',i6)')
     &                                      llnow, iprun0
               write(idf(3),'('' iseed = '',i12)') iseed1
               close(idf(3))


            end if


         end if


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      function pcmsr(a,b,c)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to determine the cm momentum from energy and mass       *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              a           : sqrt of s                                 *
*              b           : mass b                                    *
*              c           : mass c                                    *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      if( a .ge. b + c ) then


         pcmsr = sqrt( ( a**2 - ( b + c )**2 )
     &               * ( a**2 - ( b - c )**2 ) )
     &               / ( 2.0 * a )


      else


         pcmsr = 0.0


      end if


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine trfram(px,py,pz,et,rm,ifrm)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to determine energy and momentum by Lorentz transform   *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              px, py, pz     : momentum                               *
*              et             : energy                                 *
*              rm             : rest mass                              *
*              ifrm           : 0-> to lab                             *
*                               1-> to cm                              *
*                               2-> to n-n                             *
*                                                                      *
*                                                                      *
************************************************************************


      common /framtr/ betafr(0:2), gammfr(0:2)


*-----------------------------------------------------------------------


         if( ifrm .lt. 0 .or. ifrm .gt. 2 ) then


            write(*,*) ' **** Error at [trfram], unrecognized frame '
            stop


         end if


*-----------------------------------------------------------------------


            gamm = gammfr(ifrm)
            beta = betafr(ifrm)


*-----------------------------------------------------------------------


            pz = pz * gamm - beta * gamm * et
            et = sqrt( px**2 + py**2 + pz**2 + rm**2 )


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      function erf(x)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate the eror function                          *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              x           : x value for error function                *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      parameter(  dmax = 5.0 )


*-----------------------------------------------------------------------


      real*8      a(6)


*-----------------------------------------------------------------------


      data a / 0.0705230784, 0.0422820123,
     &         0.0092705272, 0.0001520143,
     &         0.0002765672, 0.0000430638 /


*-----------------------------------------------------------------------


               dnm = 1.0
               xx  = 1.0


            do i = 1, 6


               xx  = xx * x
               dnm = dnm + xx * a(i)


            end do


            if( dnm .gt. dmax ) then


               erf = 1.0


            else


               erf = 1.0 - 1.0 / dnm**16


            end if


*-----------------------------------------------------------------------


      return
      end




