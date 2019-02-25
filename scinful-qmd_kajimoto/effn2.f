      subroutine effn
c       modified by d.satoh ('99.11.27)
c       last modification on '04.08.05 by d.satoh @ jaeri
c
      real bins(1000)
      common /effnaka/  bins, sboxt, maxl
c
      do 100 i=1,1000
            bins(i) = exp( alog(1000.) /1000. *i ) - 1.
100   continue
c
cd    bins(1)=0.0
cd    do 100 i=1,999
cd          bins(i+1)=bins(i)+ 5.0e-2
cd  100 continue
cd      bins(1000) = 1000
c
      return
      end


cccccccccccccccccccccccccccccccc
      subroutine efft
cccccccccccccccccccccccccccccccc
c_modified by d.satoh ('99.11.27)
c_modified at 2001.07.19
c ( add the routine to output the response function !!)
c_last modify at 2003.06.16 -> c_jaeri
c ( modified the format of the output file for response )
c
      integer sboxt(999),total
      real effic(5),count(5), bias2(5)
      common /effnaka/  bins(1000), sboxt, maxl
      common /INIT/ Nhist,Nhits
      common /MAXWEL/ Esourc,E1,T
      common /RES/ iswt, bias(5), icont, ictrl, iesc, ilight


      data count / 0.0, 0.0, 0.0, 0.0, 0.0/


      do 5000 ijk=1,5
      bias2(ijk)=0.0
      effic(ijk)=0.0
      count(ijk)=0.0
5000  continue


      total=0


      do 200 i=1,5
          bias2(i) = bias(i) / 1.25
      do 100 ii=1,maxl
        total=total+sboxt(ii)
          binw=bins(ii+1)-bins(ii)


        if(bias2(i).ge.bins(ii+1)) goto 100
        if(bias2(i).ge.bins(ii).and.bias2(i).lt.bins(ii+1)) then
            count(i)=count(i)+float(sboxt(ii))*(bins(ii+1)-bias2(i))/binw
        end if
        if(bias2(i).lt.bins(ii)) then
            count(i)=count(i)+float(sboxt(ii))
        end if


  100 continue


      effic(i)=count(i)/float(nhits)  ! Final result !


  200 continue


      write(*,20) (effic(i),i=1,5)


      if( icont .eq. 1 ) then
            write(31,20) (effic(i),i=1,5)
   20       format(1h ,'Efficiency ==> ',5(2x,f8.6))
            write(31,*) 'c ------------------------------------------------------'
      else
            write(31,21) Esourc, (effic(i),i=1,5)
   21       format(f7.4,5(3x,f8.6))
      end if


c+++++ response function +++++
c_2003/06/16@jaeri


      if( iswt .eq. 2 ) then


          write(33,1000) esourc
 1000       format(/,'++++++++++++++++++++++++++++++++++++++++++++++++',/
     +      26h Incident neutron energy =,f7.2,5h  MeV)
          write(33,1001) nhits, nhist
 1001       format(8h Nhits =,i10,/8h Nhist =,i10)
          write(33,*) '++++++++++++++++++++++++++++++++++++++++++++++++'
          write(33,*) 'LO_low(MeVee), LO_up(MeVee), Response(/MeVee/n),
     +Error(/MeVee/n)'


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
c -------------------------------------------------------------------------------
c added by kajimoto 10/08/06
              write(34,41) esourc,bins(kk)*1.25,bins(KK+1)*1.25,resp,err
   41         format(1pe11.3,1pe11.3,1x,1pe11.3,1x,1pe11.3,1x,1pe11.3)
          end do


      end if
c+++++ response function +++++


      return
      end
