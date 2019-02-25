      Function akalbach(Ein,Eout,Ac,Zc,Za,Na,Zb,Nb)
c
c
c
      Common /RANDM/  Irx
c
      real Na, Nb, akalbach, xxcos(180)
      real xcos(180), cos1(180), cos2(180), sigma(180), sigint(180)
c
      pi=3.14159265358979
c
      Aa=Za+Na
      Ab=Zb+Nb
c
      do i=1,180 ! initialized
        theta=i-0.5
        theta1=theta-0.5
        theta2=theta+0.5
        theta=theta*pi/180.0
        theta1=theta1*pi/180.0
        theta2=theta2*pi/180.0
        cos1(i)=cos(theta1)
        cos2(i)=cos(theta2)
        xcos(i)=cos(theta)
        xxcos(i)=cos(theta)/2.0+0.5
        sigma(i)=0.0
        sigint(i)=0.0
      enddo
c
      bind=BindEng(Ac,Zc,Zb,Nb) ! outgoing particle
      esys=Ein*(Ac-Aa)/Ac+BindEng(Ac,Zc,Za,Na) ! incident particle
c
      er=1.0
      if (Aa .gt. 3.0) er=2.0
c       er = Et1/130.
      esysr=esys/er
      e1=esys
      if (esysr .gt. 190.0) then
        e1=130.0*er
      else if (esysr .gt. 80.0) then
        x=(esysr-130.0)/12.0
        x=1.0+exp(x)
        e1=esys/x
        x=(136.0-esysr)/12.0
        x=1.0+ exp(x)
        e1=e1+130.0*er/x
      end if
      if (Aa .gt. 3.0) goto 20
      e3=esys
      if (esysr .gt. 51.0) then
        e3=35.0*er
      else if (esysr .gt. 21.0) then
        x=(esysr-35.0)/3.2
        x=1.0+exp(x)
        e3=esys/x
        x=(36.6-esysr)/3.2
        x=1.0+exp(x)
        e3=e3+35.0*er/x
      end if
c
      xmb=1.0
      if(Ab.eq.4)xmb=2.
      if(Zb.eq.0)xmb=0.5
   20 continue
c
c      if(fmsd(ne).le.0.)fmsd(ne)=1.0
      fmsd=1.0 ! kajimoto
      epscm=Eout+bind
      y=epscm*e1/(esys*130.0)
      a=5.2*y+4.0*y*y*y
      if (Aa .gt. 3) goto 22
      y=epscm*e3/(esys*35.0)
      a=a+1.9*y*y*y*y*xmb
   22 xnorm=a/(12.5664*sinh(a))
      slope=a
c
      do 24 i=1,90 ! forward angle (0 ~ 90 deg.)
        arg=a*xcos(i)
        sig=fmsd*sinh(arg)+cosh(arg)
        sigma(i)=sig*xnorm
        if(i .eq. 1) then
          sigint(i)=sigma(i)*2.0*pi*(cos1(i)-cos2(i))
        else
          sigint(i)=sigint(i-1)+sigma(i)*2.0*pi*(cos1(i)-cos2(i))
        endif
   24 continue
c
      do 34 i=91,180 ! backward angle (91 ~ 180 deg.)
        arg=slope*xcos(i)
        sig=fmsd*sinh(arg)+cosh(arg)
        sigma(i)=sig*xnorm
        sigint(i)=sigint(i-1)+sigma(i)*2.0*pi*(cos1(i)-cos2(i))
   34 continue
c
      samp=sigint(180)*Ran(Irx)
      akalbach=EXTERP(sigint,xxcos,samp,180,1)
      akalbach=(akalbach-0.5)*2.0
      if(akalbach .gt. 1.0) akalbach=1.0
      if(akalbach .lt. -1.0) akalbach=-1.0
c      if(Zb .eq. 2.0 .and. Nb .eq. 1.0) write(*,*) Ein, akalbach
      return
      end
c
c======================================================================
c
      function BindEng(acom,azcom,jpout,jnout)
c
c     calculate binding energies from semi-empirical masses
c     with pairing (and shell) effects not included
c
      real xs(5,5),jpout,jnout
c
      jjpout = int(jpout)
      jjnout = int(jnout)
c
      do i=1,5
        do j=1,5
          xs(i,j)=0.0
        enddo
      enddo
      xs(2,2)=2.22
      xs(2,3)=8.48
      xs(3,2)=7.72
      xs(3,3)=28.30
c
      ancom=acom-azcom
      xout=jpout+jnout
      xpout=jpout
      ares=acom-xout
      azres=azcom-xpout
      anres=ares-azres
      acom3=acom**0.33333
      ares3=ares**0.33333
      bin=15.68*xout
      bin=bin-18.56*(acom3*acom3-ares3*ares3)
      x=(ancom-azcom)*(ancom-azcom)/acom
      y=(anres-azres)*(anres-azres)/ares
      bin=bin-28.07*(x-y)
      bin=bin+33.228*(x/acom3-y/ares3)
      bin=bin-0.717*(azcom*azcom/acom3-azres*azres/ares3)
      bin=bin+1.211*(azcom*azcom/acom-azres*azres/ares)
      BindEng=bin-xs(jjpout+1,jjnout+1)
      return
      end
c
c======================================================================
c
