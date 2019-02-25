****************************************************
*  SUBROUTINE TRANSP                               *
*  Purpose is to compute particle transport        *
*  for multiple-collision.                         *
*  (Considering energy loss of charged particle.)  *
*                                                  *
*             made by  d.satoh    may/2000         *
*             revised  at  August/2000 (c++)       *
****************************************************
c---------------------------------------------------
      subroutine transp(id,ic,en,vx,vy,vz)
c---------------------------------------------------
      common /incidnt/Enn,u,v,w,x,y,z
      common /clst/no
      common /vector/xpn,ypn,zpn,xn,yn,zn
      common /naid/ r,h,rcollim,rz
      common /cutoff/ ecutoff
      common /cascade/ ids(20),xs(20),ys(20),zs(20),ks(20),
     +           us(20),vs(20),ws(20),ens(20)
      common /collis/ nelm,echrg(6)
c
c      common /proton/ noprot, range(95), proten(95)
      COMMON/PROTON/NOPROT,RNGEN(144),PRAN(144),DRAN(144),ARAN(144),
     a              HRAN(144),TRAN(144)
      common /RES/ iswt, bias(5), icont, ictrl, iesc, ilight
c---------------------------------------------------
      if( id .eq. 0 ) return
c---------------------------------------------------
      en1 = 0.
      en2 = 0.
c---------------------------------------------------


      call dircos(vx,vy,vz,xpn,ypn,zpn)


c-direction cosine of incident particle
      zx=u
      zy=v
      zz=w


c-direction cosine of secondary particle
      call transvec(zx,zy,zz)
      uc=xn
      vc=yn
      wc=zn


      path=plngth(uc,vc,wc,x,y,z,r,h)


c-------------------------------------------------
      if( en .gt. 150 ) then


c_____________Above 150 MeV______________________


c------------------------------------------------
c     Considering only n(id=1),p(2),pi^+(20),pi^-(21)
c     as a projectile!!
c------------------------------------------------


      if( id .eq. 1 ) then
c3          totx=totalx3(id,en)   ! neutron
          totx=totalx2(en)   ! neutron


      else if( id .eq. 2 ) then
c3          totx=totalx3(id,en)   ! proton
          totx=totalx2(en)


      else if( id .eq. 20 ) then
          totx=totalx3(id,en)   ! pi^+


      else if( id .eq. 21 ) then
          totx=totalx3(id,en)   ! pi^-


      else
          goto 1919
      end if


      k=ibox2(path,totx,dist)


          if( k .lt. 1 ) then
             if( id .eq. 1 ) return   ! neutron
             goto 1919                ! other
          else


c_New coordinates
             x2=x+uc*dist
             y2=y+vc*dist
             z2=z+wc*dist
             if (z2.lt.0.0 .or. z2.gt.h) goto 1919


c_Estimate the energy-deposition along the path
             isw=2  ! energy-loss
             tls=0.


             call dedxl (isw,id,ic,en,dist,tls)


               en1 = en
             en2 = en1 - tls


             echrg(1)= en1
             echrg(2)= en2


             if(nelm .ne. -10) call bankr2


c_Infomation array as a projectile for JQMD
             no=no+1
             ids(no)=id
             xs(no)=x2
             ys(no)=y2
             zs(no)=z2
             ks(no)=k
             us(no)=uc
             vs(no)=vc
             ws(no)=wc
             ens(no)=en


          end if


      else


c_____________Below 150 MeV______________________
cs121201
c         if( en .le. 0.0 ) goto 1105
cs121201


         if( id .eq. 1 ) then


c------------------------------------------------
c           "neutron"
c------------------------------------------------
            totx=totalx(en)
            k=ibox(path,totx,dist)


            if( k .lt. 1 ) return
            if( en .lt. ecutoff ) return
c_New coordinates
            x2=x+uc*dist
            y2=y+vc*dist
            z2=z+wc*dist
            if( z2 .lt. 0.0 .or. z2 .gt. h) return


            call scin (k,en,uc,vc,wc,x2,y2,z2)


         else


c------------------------------------------------
c           "not neutron"
c------------------------------------------------
            goto 1919


         end if




      end if


      return


 1919 continue
c-------------------------------------------------
c    Caluculate the energy deposition,
c    and accumlate the light output!!
c-------------------------------------------------
      rg = 0.
      nterp = 4


      if( en .gt. 1.0 ) then
        if(nelm.eq.-6) then
          rg = EXTERP(RNGEN, PRAN, en, noprot, nterp) ! proton
ckajimoto            rg = EXTERP(proten, range, en, noprot, nterp)
        else if(nelm.eq.-3) then
          rg = EXTERP(RNGEN, DRAN, en, noprot, nterp) ! deuteron
        else if(nelm.eq.-4) then
          rg = EXTERP(RNGEN, TRAN, en, noprot, nterp) ! triton
        else if(nelm.eq.-5) then
          rg = EXTERP(RNGEN, HRAN, en, noprot, nterp) ! 3He
        else if(nelm.eq.-2) then
          rg = EXTERP(RNGEN, ARAN, en, noprot, nterp) ! alpha
        else
            isw=1  ! range
            call dedxl (isw,id,ic,en,rg,tls)
        end if
      end if

      rangeq = rg - path

      echrg(1) = en
      eprotf = 0.


      if( id .eq. 2 .and. en .le. 1e-1 ) goto 666     !Thd. proton ---> 0.1 MeV
      if( id .ne. 2 .and. en .le. 1e+1 ) goto 666     !Thd. other  ---> 10. MeV


      if( rangeq .gt. 0.1 ) then    ! escaping

c            if( iesc .eq. 2 ) goto 1105   ! elimination of escaping events.
c            if( id .eq. 2 .and. en .ge. 1e+1 ) then
        if(nelm.eq.-6) then
          reng = EXTERP(PRAN, RNGEN, rangeq, noprot, nterp) ! proton
          eprotf = en - reng
        else if(nelm.eq.-3) then
          reng = EXTERP(DRAN, RNGEN, rangeq, noprot, nterp) ! deuteron
          eprotf = en - reng
        else if(nelm.eq.-4) then
          reng = EXTERP(TRAN, RNGEN, rangeq, noprot, nterp) ! triton
          eprotf = en - reng
        else if(nelm.eq.-5) then
          reng = EXTERP(HRAN, RNGEN, rangeq, noprot, nterp) ! 3He
          eprotf = en - reng
        else if(nelm.eq.-2) then
          reng = EXTERP(ARAN, RNGEN, rangeq, noprot, nterp) ! alpha
          eprotf = en - reng
        else
          isw = 2     ! energy-loss
          tls = 0.
          call dedxl (isw, id, ic, en, path, tls)
          eprotf = en - tls
        end if

      end if


  666 continue
      echrg(2) = eprotf
      call bankr2


 1105 continue


      return
      end
