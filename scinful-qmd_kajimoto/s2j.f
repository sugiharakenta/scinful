********************************************************************
* Subroutine jqmd                                                  *
*                                                                  *
*     August/2000                                                  *
********************************************************************
      subroutine jqmd (ip,itype)
*-----------------------------------------------------------------------


      include 'param01.inc'
      include 'param02.inc'


*-----------------------------------------------------------------------


      common /randm/ irx
      common /ncll/ ncll
      common /incidnt/ enp,u,v,w,x,y,z
      common /flag/ iflag


*-----------------------------------------------------------------------


      common /swich2/ icfg, imany, icpus, idatm
      common /input1/ mstq1(mxpa1), parq1(mxpa1)
      common /inecho/ icnum, indata(200), indlng(200)
      common /fnamid/ fname(mxfnm), idf(mxfnm), ifnl(mxfnm)


*-----------------------------------------------------------------------


      character       fname*80
      character       indata*80


*-----------------------------------------------------------------------
*     icfg = 4;  no input file and
*                take input data from main by variables
*-----------------------------------------------------------------------


               icfg = 4


*-----------------------------------------------------------------------


      enp = enp / 1000.0     ! (GeV)
      parq1(1) = enp         ! incident energy in lab [ GeV per nucleon ]
      mstq1(10) = irx        ! iseed : random seed.


*_Input data for JQMD----------------------------------------------------


*______________Identification of projectile!!___________________________
*-----------------------------------------------------------------------
*        neutron( ip = 1 ) and proton( ip = 2 )
*-----------------------------------------------------------------------


            if( ip .eq. 1 ) then


                    mstq1(1) = 0
                    mstq1(2) = 1
                    mstq1(3) = 0
                    goto 666


            end if


            if( ip .eq. 2 ) then


                    mstq1(1) = 0
                    mstq1(2) = 1
                    mstq1(3) = 1
                    goto 666


            end if


*-----------------------------------------------------------------------
*        pi^+( ip = 20 )  and  pi^-( ip = 21 )
*-----------------------------------------------------------------------


            if( ip .eq. 20 ) then


                    mstq1(1) = 211
                    mstq1(2) = 1
                    mstq1(3) = 1
                    goto 666


            end if


            if( ip .eq. 21 ) then


                    mstq1(1) = -211
                    mstq1(2) = 1
                    mstq1(3) = -1
                    goto 666


            end if


      write(*,*) 'Error at s_jqmd * can not identify that projectile!!'
            return


*-----------------------------------------------------------------------
*______________Identification of target!!_______________________________
  666 continue
            if( itype .eq. 1 ) then


                     mstq1(4) = 0    ! H target
                     mstq1(5) = 1
                     mstq1(6) = 1
                     mstq1(17) = 0   ! ielst : choice of elastic reaction
                     parq1(4) = 2.80 ! bmax  : maximum impact parameter (fm).


            else if( itype .eq. 2 ) then


                     mstq1(4) = 0    ! 12C target
                     mstq1(5) = 12
                     mstq1(6) = 6
                     mstq1(17) = 0   ! ielst : choice of elastic reaction
                     parq1(4) = 4.61 ! bmax  : maximum impact parameter (fm).


            else


                     write(*,*) 'Error at s_jqmd * wrong itype!!'
                     return


            end if
*-----------------------------------------------------------------------


      mstq1(7) = 1      ! iprun : total number of siumulation run
      mstq1(8) = 100    ! ntmax : total number of time step.
      mstq1(9) = 0      ! insys : (1=cm, 2=nn, 0=lab)
      parq1(3) = 0.0    ! bmin  : minimum impact parameter (fm).
      parq1(5) = 1.0    ! dt    : time step size(fm/c)


      mstq1(11) = 1     ! choice of impact parameter bin
      mstq1(12) = 20    ! number impact parameter bin
      mstq1(13) = 2     ! order of RKG
      mstq1(14) = 0     ! imany : multi-run control
      mstq1(15) = 0     ! output for QMD results
      mstq1(16) = 0     ! input of QMD results


      fname(4)  = 'jqmd'      ! Header of the Output File Name for ANGEL
      fname(10) = 'jqmd.dat'  ! output file for QMD results
      fname(11) = 'jqmd.dat'  ! input file of QMD results


      mstq1(120) = 1    ! issdm  : SDM included or not
      mstq1(121) = 0    ! iswids : choice of decay width
      mstq1(124) = 0    ! isgrnd : ground state decay


      mstq1(150) = 0    ! jdsp(20)   : Display Summary on file
      mstq1(151) = 0    ! jdsp(21)   : Display Header and Summary(1)
      mstq1(152) = 0    ! jdsp(22)   : Display Collision History 01
      mstq1(153) = 0    ! jdsp(23)   : Display Ground State Properties
      mstq1(154) = 0    ! jdsp(1-13) : Choice of out put graph
      mstq1(155) = 0    ! jdsp(24)   : Display JQMD Logo and CPU Time on Screen
      mstq1(156) = 0    ! jdsp(25)   : Display Reaction on Screen
      mstq1(157) = 0    ! jdsp(26)   : Display Reaction on file(4)
      mstq1(158) = 0    ! jdsp(27)   : Number of Events for Display or File
      mstq1(159) = 0    ! jdsp(30)   : Display Mass Distribution of QMD
      mstq1(160) = 0    ! jdsp(31)   : Display Mass Distribution of SDM
      mstq1(161) = 0    ! jdsp(32)   : Display Mass Distribution of Fission


      mstq1(169) = 0    ! icpus    : detect cpu time
      mstq1(170) = 0    ! idatm    : detect date and time
      mstq1(171) = 5    ! nfreq    : time interval for output : 0 => no
      mstq1(172) = 10   ! nfrec    : time interval for cluster output: 0 => no
      mstq1(173) = 1    ! nfred    : time interval for QMDDISP: 0 => no
      mstq1(174) = 0    ! jdsp(33) : QMDDISP
      mstq1(175) = 1    ! call user anal
*-----------------------------------------------------------------------


      ncll = ncll + 1


*-----------------------------------------------------------------------
*_satoh
*        if the reaction type is nucleon-nucleon,
*        avoid ichannel = 0 @ s_crosw.
*-----------------------------------------------------------------------
               iflag = 0
      if( ip .eq. 1 .or. ip .eq.2 ) then
          if( itype .eq. 1 ) iflag = 1   ! nucleon-nucleon
      end if


*-----------------------------------------------------------------------
*        call jqmd00
*-----------------------------------------------------------------------


               call jqmd00


*-----------------------------------------------------------------------


      end




************************************************************************
*                                                                      *
      subroutine anal_int
*                                                                      *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              user subroutine for analysis                            *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


      include 'param02.inc'


*-----------------------------------------------------------------------
      common /collis/ Nelm, Echrg(6)


      common /summ04/ ianal
      common /const1/ elab, rdist, bmin, bmax, ibch, ibin
      common /const2/ dt, ntmax, iprun, iprun0
      common /vriab3/ qmdfac, sdmfac
      common /swich3/ ielst, jelst, kelst


*-----------------------------------------------------------------------


      dimension wi(10), wt(10)


*-----------------------------------------------------------------------


      if( ianal .eq. 0 ) return


*-----------------------------------------------------------------------


      return




************************************************************************
*                                                                      *
      entry anal_qmd(ik,jj,iz,in,id,is,ic,iq,im,
     &               bi,px,py,pz,et,rm,ex)
*                                                                      *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              user entry for analysis of QMD                          *
*                                                                      *
*                                                                      *
************************************************************************


      if( ianal .eq. 0 ) return


*-----------------------------------------------------------------------
*        inelastic frag and detection frame
*-----------------------------------------------------------------------


            if( kelst .eq. 0 ) return


c              call trfram(px,py,pz,et,rm,0)


*-----------------------------------------------------------------------
*        nucleus( ik = 0 )
*-----------------------------------------------------------------------


            if( ik .eq. 0 ) then




            end if


*-----------------------------------------------------------------------
*        proton( ik = 1 ) and neutron( ik = 2 )
*-----------------------------------------------------------------------


            if( ik .eq. 1 ) then


            end if


            if( ik .eq. 2 ) then


            end if


*-----------------------------------------------------------------------
*        pions( ik = 5 )
*-----------------------------------------------------------------------


            if( ik .eq. 5 ) then


               if( ic .eq. 1 ) then


               else if( ic .eq. 0 ) then


               else if( ic .eq. -1 ) then


               end if


            end if




*-----------------------------------------------------------------------


      return




************************************************************************
*                                                                      *
      entry anal_sdm(ik,jj,iz,in,id,is,ic,iq,im,
     &               bi,px,py,pz,et,rm,ex)
*                                                                      *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              user entry for analysis of SDM                          *
*                                                                      *
*                                                                      *
************************************************************************
*-----------------------------------------------------------------------
      data c / 2.997925e+10 /
*-----------------------------------------------------------------------


      if( ianal .eq. 0 ) return


*-----------------------------------------------------------------------
*        inelastic frag and detection frame
*-----------------------------------------------------------------------


            if( kelst .eq. 0 ) return


            call trfram(px,py,pz,et,rm,0)


            ip=0
*-----------------------------------------------------------------------
*        nucleus( ik = 0 )
*-----------------------------------------------------------------------


            if( ik .eq. 0 ) then


                if( iz .eq. 1 ) then


                    if( in .eq. 1 ) then ! deuteron


                        ip = 10
                      Nelm = -3


                    else if( in .eq. 2 ) then ! triton


                        ip = 13
                      Nelm = -4


                    end if


                else if( iz .eq. 2 ) then


                    if( in .eq. 1 ) then ! 3He


                        ip = 19
                      Nelm = -5


                    else if( in .eq. 2 ) then ! alpha


                        ip = 3
                      Nelm = -2


                    end if


                 end if


            end if


*-----------------------------------------------------------------------
*        proton( ik = 1 ) and neutron( ik = 2 )
*-----------------------------------------------------------------------


            if( ik .eq. 1 ) then


                    ip = 2
                  Nelm = -6


            end if


            if( ik .eq. 2 ) then


                    ip = 1
                  Nelm = -10


            end if




*-----------------------------------------------------------------------
*        pions( ik = 5 )
*-----------------------------------------------------------------------


            if( ik .eq. 5 ) then


               if( ic .eq. 1 ) then


                       ip = 20
                     Nelm = -1


               else if( ic .eq. 0 ) then


               else if( ic .eq. -1 ) then


                        ip = 21
                      Nelm = -1


               end if


            end if


*-----------------------------------------------------------------------


      ekin = ( et - rm ) * 1000.0     ! (MeV)
      ett  = et * 1000.0


      betax = px/ett
      betay = py/ett
      betaz = pz/ett
      beta  = sqrt(betax**2+betay**2+betaz**2)


c_check
              if (beta .ge. 1.0 ) then
                  write(*,*) 'Error at s_jqmd : beta > 1 !!'
                  return
              end if


      vx = betax * c     ! (cm/s)
      vy = betay * c
      vz = betaz * c


      call transp(ip,ic,ekin,vx,vy,vz)


*-----------------------------------------------------------------------


      return




************************************************************************
*                                                                      *
      entry anal_fin
*                                                                      *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              user entry for analysis of weight                       *
*                                                                      *
*                                                                      *
************************************************************************


      if( ianal .eq. 0 ) return


*-----------------------------------------------------------------------




*-----------------------------------------------------------------------


      return
      end




