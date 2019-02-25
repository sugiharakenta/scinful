****************************************************
*  SUBROUTINE SICN                                 *
*  Poupose is to compute neutron transportation    *
*  and nuclear reactions below 150 MeV by SCINFUL. *
*                                                  *
*            by d.satoh   may/2000                 *
****************************************************
      subroutine scin (itype,enp,up,vp,wp,xp,yp,zp)
c--------------------------------------------------
      common /naid/ rdet,ht,rcollim,rz
      common /cutoff/ ecutoff
c++      common /incidnt/ enp,up,vp,wp,xp,yp,zp
      common /neutrn/ Eneut,u,v,w,x,y,z
      common /neutr2/ Eneut2,u2,v2,w2,x2,y2,z2
c++      common /particle/ ip,itype
      common /ncll/ Ncll
c008      common /ncolls/ Ndie
      Common /NCOLLS/ Nclsns(21),Ndie,Nfstcl,Ntcoll
c--------------------------------------------------
c_To  common /neutrn/
      Eneut=enp
      u=up
      v=vp
      w=wp
      x=xp
      y=yp
      z=zp
c--------------------------------------------------
cs121201
      if (Eneut .LE. 0.0) then
          write(*,*) '*** Error in Function SCIN; En = ',Eneut
          stop
      end if
c--------------------------------------------------
c_check
c      write(17,*) 'Check subroutine scin!!'
c      write(17,*) 'ncoll,itype',ncoll,itype
c      write(17,*) 'ip,itype',ip,itype
c      write(17,*) 'Eneut,u,v,w,x,y,z',Eneut,u,v,w,x,y,z
c--------------------------------------------------
      k=itype
   10 Ncll=Ncll+1
      goto (20,30,40,50,60,70,80,110,90,100) k
   20 call hydrog
      goto 200
   30 call pscat
      goto 200
   40 call inelas
      goto 200
   50 call nalpha
      goto 250
   60 call nn3alf
      goto 200
c--------------------------------------------------
c_Next four routines may or may not have
c a secondary output neutron.
   70 call n3he
          if(Eneut) 250,250,200
   80 call npx
          if(Eneut) 250,250,200
   90 call nd
          if(Eneut) 250,250,200
  100 call nt
          if(Eneut) 250,250,200
  110 call n2n
  200 if(Eneut.lt.ecutoff) goto 240
cs  200 write(17,*)'Eneut,ecutoff',Eneut,ecutoff
cs      if(Eneut.lt.ecutoff) goto 240
c--------------------------------------------------
c_Determine if the scattered neutron interacts again
c in the detector.
      path=plngth(u,v,w,x,y,z,rdet,ht)
      totx=totalx(Eneut)
  210 Itype=ibox(path,totx,dist)
c               write(17,*)'path,totx,dist,Itype',path,totx,dist,Itype
c               write(17,*)'pre/Eneut2=',Eneut2
c
c_If Itype = 0, the neutron has escaped without further collision.
      if(Itype.le.0) goto 250
c
c_Check to see if there is a Neutron No. 2 (from N2N or NPN routines)
      if(Eneut2.le.0.0) goto 220
c
c_If there is one, no sense in choosing the N2N routine to go
c through on this pass, so choose Itype again. (MODIFIED 8/88)
      if(k.eq.8) goto 210
c--------------------------------------------------
c_O.K! Neutron interacted again. Find new interaction spot!!
  220 x=x+u*dist
      y=y+v*dist
      z=z+w*dist
c      write(17,*)'(x,y,z),k',x,y,z,Itype
      if(z.lt.0. .or. z.gt.ht) goto 250


c_Last may be an unnecessary check on single-precision arithmetic.
      k=Itype
      goto 10
c--------------------------------------------------
c_That's for the case where E(neutron) < Ecutoff.
  240 Ndie=Ndie+1
c--------------------------------------------------
cs  250 write(17,*)'Eneut2=',Eneut2
cs      if(Eneut2.le.0.0) goto 280
  250 if(Eneut2.le.0.0) goto 280
c--------------------------------------------------
c_Someting in the 2nd neutron arry
c      write(17,*)'aft/Eneut2=',Eneut2
      Eneut=Eneut2
      u=u2
      v=v2
      w=w2
      x=x2
      y=y2
      z=z2
      Eneut2=0.
      u2=0.
      v2=0.
      w2=0.
      x2=0.
      y2=0.
      z2=0.
      goto 200
c--------------------------------------------------


cs  280 write(17,*)'O.K! subroutine scin finished!!'
  280 return
cs      return
      end


