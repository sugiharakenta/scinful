**********************************************
*   SUBROUTINE DEDXL                         *
*   Purpose is to calculate raange or        *
*   energy-loss by Bethe formula.            *
*                                            *
*           August/2000                      *
**********************************************
c---------------------------------------------
      subroutine dedxl(isw,idn,ic,en,dist,tls)
c---------------------------------------------
       common /MASSES/ Emass(21)
       common /SCINT/ Scntyp
c---------------------------------------------
       real*8 ta,tz,fion,rat,dens,den,pm,dt,tta,
     + uu,range,pz
       integer*4 ntgz
       parameter( ntgz = 2 )
c---------------------------------------------
       dimension ta(2),tz(2),fion(2),
     + rat(2),den(2)
c---------------------------------------------
c      Configulation
c---------------------------------------------
       data( ta(i), i=1,ntgz ) /12.,1./
       data( tz(i), i=1,ntgz ) /6.,1./
       data( fion(i), i=1,ntgz ) /73.80,20.40/
       if( Scntyp .eq. 110 ) then
         data( rat(i), i=1,ntgz ) /1.,1.104/
         dens = 1.0320d0
       else
         data( rat(i), i=1,ntgz ) /1.,1.213/
         dens = 0.8740d0
       end if
c---------------------------------------------
       dt = 0.
       ma = Emass(idn)      ! mass of projectile (MeV)


       pm = dble(ma)
       ic = abs(ic)
       pz = dble(ic)        ! charge of projectile
       if( isw .eq. 1 ) then
       dt = 1.0d-1          ! thickness of each layers (cm)
       else if( isw .eq. 2 ) then
       dt = 1.0d-3          ! thickness of each layers (cm)
       end if
c---------------------------------------------
       dens = dens * 1.0d3  ! density (mg/cm3)
       tta = 0.


       do i =1, ntgz
          tta = tta + ta(i) * rat(i)
       end do


       do i =1, ntgz
          den(i) = dens * ta(i) * rat(i) / tta
       end do
c---------------------------------------------
       uu = dble (en)


       call rangecal( isw,ntgz,uu,pm,pz,den,dt,dist,range,
     + ta,tz,fion)


       if( isw.eq.1 ) then
           dist = sngl(range)


       else if( isw.eq.2 ) then
           tls = en - sngl(uu)


       end if


      if( tls .lt. 0.0 ) then
            write(*,*) 'Error @ s_dedxl: tls =',tls
      end if


       return
       end
*=========================================================*
c---------------------------------------
      subroutine rangecal( isw,ntgz,uu,pm,pz,den,dt,dist,range,
     +                     ta,tz,fion )
c---------------------------------------
      real*8 tdedx, beta2, a, dedx
      real*8 ta,tz,pz,fion,rat,dens,den,pm,dt,tta,
     +uu,range
      integer*4 ix,ntgz
c---------------------------------------
       dimension ta(2),tz(2),fion(2),
     + rat(2),den(2)
c---------------------------------------


      dd = dble(dist)


      do ix = 1,500000000
      tdedx = 0.


      do i = 1,ntgz
      beta2= ( uu**2.+2.0d0*pm*uu )/(uu+pm)**2.  ! square of beta


      a= log( (2.d0*511006.d0*beta2)/(fion(i)*(1.d0-beta2)) )-beta2


      dedx= 0.3070614557d-3 * pz**2 *  tz(i)*a/(beta2*ta(i))*den(i)


      tdedx= tdedx + dedx
      end do


      if( isw .eq. 1 ) then


          uu = uu - tdedx*dt
          if( uu .le. 0.2d-3 .or. a .le. 0.d0 ) then
          range = ix * dt


          go to 100
          end if


      else if( isw .eq. 2 ) then


          uu = uu - tdedx*dt
          range = ix * dt
          if( range .ge. dd-1.d-3 .and. range .le. dd+1.d-3 )goto 100
          if( uu .le. 0.2d-3 .or. a .le. 0.d0 ) then
          uu=0.
          goto 100
          end if


      else
          write(*,*) 'error at s_rangecal !!'
      end if


      end do


  100 return
      end
