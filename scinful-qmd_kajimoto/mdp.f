************************************************************************
*                                                                      *
*        PART 11: Machine dependent routines                           *
*                                                                      *
*                          by d.satoh
*                                                                      *
*   List of subprograms in rough order of relevance with main purpose  *
*      ( s = subroutine, f = function, b = block data, e = entry )     *
*                                                                      *
*                                                                      *
*  f  rn()           to get random number                              *
*  s  date_a_time    to get current date                               *
*  f  sect_a         to get elapse time                                *
*  f  cput_a         to get cpu time                                   *
*                                                                      *
************************************************************************
*                                                                      *
************************************************************************
*                                                                      *
      function rn()
*                                                                      *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to get random number                                    *
*              use function : ran(iseed), iseed in common              *
*                                                                      *
*                                                                      *
************************************************************************


      common /rannum/ iseed, iseed0, iseed1


*-----------------------------------------------------------------------


            rn = ran(iseed)


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      subroutine date_a_time(iyer,imon,iday,ihor,imin,isec)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to get current date                                     *
*              use subroutine date_and_time                            *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              iyer     : 4-digit year                                 *
*              imon     : month                                        *
*              iday     : day                                          *
*              ihor     : hour                                         *
*              imin     : minite                                       *
*              isec     : secont                                       *
*                                                                      *
*                                                                      *
************************************************************************


      dimension ivalue(8)
      character*12 dum(3)


*-----------------------------------------------------------------------


            call date_and_time(dum(1),dum(2),dum(3),ivalue)


            iyer = ivalue(1)
            imon = ivalue(2)
            iday = ivalue(3)
            ihor = ivalue(5)
            imin = ivalue(6)
            isec = ivalue(7)


*-----------------------------------------------------------------------


      return
      end




************************************************************************
*                                                                      *
      function sect_a(stime)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to get elapse time                                      *
*              use function : seconds(stime)                           *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              stime    : starting time                                *
*              sect_a   : elapse time from starting time               *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------


            sect_a = secnds(stime)


*-----------------------------------------------------------------------


      return
      end


************************************************************************
*                                                                      *
      function cput_a(cputm)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to get cpu time                                         *
*              use function : dtime(tarray)                            *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              cputm    : starting time                                *
*              cput_a   : cpu time from starting time                  *
*                                                                      *
*                                                                      *
************************************************************************


      dimension tarray(2)


      data cpuall / 0.0 /


      save cpuall


*-----------------------------------------------------------------------


            cpuall = cpuall + dtime(tarray)


            cput_a = cpuall - cputm


*-----------------------------------------------------------------------


      return
      end


