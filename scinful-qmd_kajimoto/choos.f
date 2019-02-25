***************************************
*  SUBROUTINE                         *
*     choosa                          *
*     choosl                          *
*     choosm                          *
*     choosn                          *
*     choosp                          *
*     p3alph                          *
*                                     *
***************************************
c=============================================
c   This is file CHOOSA.FOR
c
c    Purpose is to choose an alpha energy from a simple formalism
c       modestly representing a distribution of the energy of one


c       of the alphas (the first in our case) from the 12-C --> 3 alphas


c       3-body breakup reaction.


c
      Function CHOOSA(Emax)
c
c         The formalism used is:   Phi = E**2 * (Emax - E)**2


c
      Common /RANDM/ Irx
      Dimension E(101), S(101)
      Data N/100/, M/101/, Nterp/1/, Eold/0.0/
c
      Em=Emax
      IF (Em .EQ. Eold)  goto 8
      Eold=Em
      D=0.01*Em
      Halfd=0.5*D
      E(1)=0.0
      S(1)=0.0
      Do 5 I=1,N
      Ek=E(I)+Halfd
      E(I+1)=Ek+Halfd
      Sqrtphi=(Em-Ek)*Ek
    5 S(I+1)=S(I) + Sqrtphi*Sqrtphi
c
      Sx=S(M)
c
    8 R=Sx*Ran(Irx)
      Chosen=EXTERP(S,E,R,M,Nterp)
      CHOOSA=Chosen
      Return
      END
c=============================================
c   This is file CHOOSL.FOR
c
c       Purpose to choose an azimuthal angle, or actually the cosine of the


c         chosen angle, from a distribution described by Legendre polynomial


c         expansion.  The technique is from an old "theorem" in which one


c         computes the maximum value the function can have within the


c         independent variable's domain of interest (in this case from -1


c         to +1) and then compares that to a value of the function for a


c         randomly chosen value of the independent variable.


c
      Function CHOOSL(Fi,N)
c
      Common /RANDM/ Irx
      Dimension P(10),Fi(1)
c       There are N (where N .LE. 10, so be careful) coefficients Fi at


c         entry.


c
C       Next portion of programming is from O5S code.


      A = 1.
      DO 19 J = 1, N
   19 A = A + FLOAT(J+J+1) * ABS (FI(J))
c               In this case -A- is .GE. the maximum value that can be


c               obtained from this polynomial expansion, and not just


c               equal to the maximum value.


   20 P(1) = 2. * RAN(IRX) - 1.
      P(2) = ( 3. * P(1)**2 - 1. ) / 2.
      NM = N - 1
      FL = 1.
      DO 21 L = 2, NM
      FL = FL + 1.
      P(L+1) = ( P(1) * ( FL + FL + 1. ) * P(L) - FL * P(L-1) ) /
     1  ( FL + 1. )
   21 CONTINUE
      X = 1.
      FL = 0.
      DO 22 L = 1, N
      FL = FL + 1.
      X = X + ( FL + FL +1. ) * FI(L) * P(L)
   22 CONTINUE
      R=RAN(Irx)
      IF ( X .LT. A * R )  GOTO 20
      CHOOSL = P(1)
      Return
      END
c=============================================
C   This is file CHOOSM.FOR
c
c     Purpose is to choose a neutron energy by random number generation from
c       a Maxwellian distribution of energies given by:


c
c               Phi(En) = SQRT(En) * EXP(-En/T),


c
c       where T is a "temperature".


c
      Function CHOOSM(Enmin,Enmax,Temp)
c
      Common /RANDM/ Irx
      Dimension En(401),Sn(401)
      Data Eold/0.0/, Elow/0.0/, Told/0.0/
c
      N=100
      Elo=Enmin
      T=Temp
      E=Enmax
      IF (E .EQ. Eold  .AND.  T .EQ. Told) goto 8
c       Set up new En, Sn arrays only for new value of Enmax and/or Temp


      Told=T
      Eold=E
      Edel=E
      IF (Edel .LE. 10.0) goto 2
      IF (Edel .GT. 20.0*T) Edel=20.0*T
      N=10*IFIX(Edel)
      IF (N .GT. 400) N=400
    2 D=Edel/FLOAT(N)
      En(1)=0.0
      Sn(1)=0.0
      Halfd=0.5*D
      Do 5 I=1,N
      Ek=En(I)+Halfd
      En(I+1)=Ek+Halfd
      Phi=SQRT(Ek)*EXP(-Ek/T)
    5 Sn(I+1)=Sn(I) + Phi
      M=N+1
      Snx=Sn(M)
c
    8 IF (Elo .EQ. Elow) goto 10
      Elow=Elo
      Sny=EXTERP(En,Sn,Elo,M,1)
   10 R=Sny + (Snx-Sny)*RAN(Irx)
      Chosen=EXTERP(Sn,En,R,M,1)
   15 CHOOSM=Chosen
      Return
      End
c=============================================
C   This is file CHOOSN.FOR
c
c     Purpose is to choose a neutron energy by random number generation from
c       a distribution of energies given by:


c
c               Phi(En) = EXP(-En/T) * En**1.5,


c
c       where T is a "temperature".  This functional dependence approximates


c        the relative "continuum" (of) 12-C excitations by neutrons having


c        20 MeV < En < 70 MeV as computed by the Hauser-Feshbach compound-


c        nucleus plus pre-compound formalism in the code TNG.


c
      Function CHOOSN(Enmin,Enmax,Temp)
c
      Common /RANDM/ Irx
      Dimension En(401),Sn(401)
      Data Eold/0.0/, Elow/0.0/, Told/0.0/
c
      N=100
      Elo=Enmin
      T=Temp
      E=Enmax
      IF (E .EQ. Eold  .AND.  T .EQ. Told) goto 8
c       Set up new En, Sn arrays only for new value of Enmax and/or Temp


      Told=T
      Eold=E
      Edel=E
      IF (Edel .LE. 10.0) goto 2
      IF (Edel .GT. 20.0*T) Edel=20.0*T
      N=10*IFIX(Edel)
      IF (N .GT. 400) N=400
    2 D=Edel/FLOAT(N)
      En(1)=0.0
      Sn(1)=0.0
      Halfd=0.5*D
      Do 5 I=1,N
      Ek=En(I)+Halfd
      En(I+1)=Ek+Halfd
      Phi=Ek**1.5 * EXP(-Ek/T)
    5 Sn(I+1)=Sn(I) + Phi
      M=N+1
      Snx=Sn(M)
c
    8 IF (Elo .EQ. Elow) goto 10
      Elow=Elo
      Sny=EXTERP(En,Sn,Elo,M,1)
   10 R=Sny + (Snx-Sny)*RAN(Irx)
      Chosen=EXTERP(Sn,En,R,M,1)
   15 CHOOSN=Chosen
      Return
      End
c=============================================
C   This is file CHOOSP.FOR
c
c    Purpose is to choose a proton energy from a simple formalism for the
c       (n,p) and (n,np) reactions.


c    The formulation used is:
c
c       Phi(Eproton) = Eprot*EXP(-Eprot/Temp)*(1.0-EXP(-K*Eprot/B))


c
c               where Temp = "Nuclear Temperature", B=barrier penetration


c               factor , and K=an empirical parameter selected to give the


c               best fit of (1.0-EXP(-K*Eprot/B)) to the barrier penetration


c               probability of 1.5 MeV.  The value for B was deduced to be


c               2.9 MeV for the O5S program, and it appears that the value


c               of K was taken to be 1.5 MeV.  That's why the "B" in the


c               Data statement below was set = (2.9/1.5).  The given


c               functional dependence reproduces adequately the shape of the


c               n + 12-C --> p + 12-B proton continuum (in 12-B) excitation


c               as computed by TNG.


c
c
      Function CHOOSP(Epmax,Temp)
c
      Common /RANDM/ Irx
      Dimension Epr(401),Sp(401)
      Data B/1.93333/, Eold/0.0/, Told/0.0/, Nterp/1/
c
      N=100
      E=Epmax
      IF (E .EQ. Eold .AND. Temp .EQ. Told) goto 10
c               Set up Epr and Sp arrays for new energy Epmax or Temp


      Told= Temp
      Eold=E
      IF (E .LE. 10.0) goto 2
      N=10*IFIX(E)
      IF (N .GT. 400) N=400
    2 D=E/FLOAT(N)
      Epr(1)=0.0
      Sp(1)=0.0
      Halfd=0.5*D
c
      Do 5 J=1,N
      Eprj=Epr(J)+Halfd
      Epr(J+1)=Eprj+Halfd
      Phi=Eprj*EXP(-Eprj/Temp) * (1.0 - EXP(-Eprj/B))
      Sp(J+1)=Sp(J)+Phi
    5 Continue
      M=N+1
      Spx=Sp(M)
c
   10 Rspx=Spx*RAN(Irx)
      Chosen=EXTERP(Sp,Epr,Rspx,M,Nterp)
   15 CHOOSP=Chosen
      Return
      END
c=============================================
C   This is file P3ALPH.FOR
c
c       Purpose to return array PBAR, the summed probabilities for the 14 (in


c       effect) sub-reactions leading to the (n,n' 3alpha) break-up reaction.


c       See Subroutine NN3ALF for the specific branching designations.


c
      SUBROUTINE P3ALPH(E)
c
      Common /P3ALF/ Pbar(13)
c
      Dimension Eofp(26),P1(26),P2(26),P3(26),P4(26),P5(26),P6(26),
     X P7(26),P8(26),P9(26),P10(26),P11(26),P12(26),P13(26)
c
      Data Eofp/9.7, 10.2,  10.5,  11.0,  11.85, 12.0,  12.5, 13.0,
     V 13.2,  13.7,  14.3,  15.3,  16.2,  17.5,  18.5,  19.6,  20.,
     W 22.,   23.2,  25.,   30.,   35.,   40.,   50.,   60.,   70./
      Data P1/1.0,   0.86,  0.672, 0.244, 0.105, 0.095, 0.074, 0.059,
     V 0.054, 0.045, 0.04,  0.032, 0.029, 0.031, 0.31,  0.034, 0.034,
     W 0.022, 0.017, 0.02,  0.026, 0.023, 0.016, 0.01,  0.0,   0.0/
      Data P2/1.0,   0.86,  0.672, 0.478, 0.476, 0.473, 0.493, 0.469,
     V 0.455, 0.43,  0.413, 0.369, 0.337, 0.322, 0.271, 0.238, 0.222,
     W 0.174, 0.180, 0.181, 0.183, 0.201, 0.196, 0.195, 0.195, 0.195/
      Data P3/1.0,   1.0,   0.943, 0.822, 0.748, 0.738, 0.729, 0.666,
     V 0.636, 0.577, 0.531, 0.449, 0.397, 0.369, 0.309, 0.272, 0.253,
     W 0.196, 0.196, 0.194, 0.191, 0.206, 0.196, 0.195, 0.195, 0.195/
      Data P4/1.0,   1.0,   1.0,   1.0,   1.0,   0.99,  0.955, 0.863,
     V 0.819, 0.73,  0.659, 0.536, 0.463, 0.422, 0.351, 0.309, 0.288,
     W 0.220, 0.214, 0.209, 0.199, 0.21,  0.196, 0.195, 0.195, 0.195/
      Data P5/1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   0.991, 0.921,
     V 0.88,  0.79,  0.729, 0.601, 0.52,  0.472, 0.391, 0.346, 0.322,
     W 0.242, 0.229, 0.218, 0.203, 0.212, 0.196, 0.195, 0.195, 0.195/
      Data P6/1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     W 0.994, 0.875, 0.938, 0.802, 0.696, 0.601, 0.508, 0.447, 0.416,
     W 0.305, 0.278, 0.254, 0.222, 0.222, 0.201, 0.195, 0.195, 0.195/
      Data P7/1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     W 1.0,   1.0,   0.997, 0.907, 0.827, 0.766, 0.662, 0.598, 0.563,
     W 0.407, 0.349, 0.301, 0.276, 0.256, 0.217, 0.197, 0.195, 0.195/
      Data P8/1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     W 1.0,   1.0,   1.0,   1.0,   0.986, 0.96,  0.833, 0.755, 0.709,
     W 0.505, 0.425, 0.358, 0.301, 0.269, 0.224, 0.197, 0.195, 0.195/
      Data P9/1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     W 1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   0.94,  0.792, 0.744,
     W 0.535, 0.452, 0.384, 0.317, 0.278, 0.229, 0.197, 0.195, 0.195/
      Data P10/1.0,  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     W 1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   0.985, 0.93,
     W 0.653, 0.533, 0.441, 0.373, 0.307, 0.247, 0.203, 0.195, 0.195/
      Data P11/1.0,  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     W 1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   0.985, 0.977,
     W 0.92,  0.894, 0.86,  0.81,  0.764, 0.745, 0.767, 0.765, 0.75/
      Data P12/1.0,  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     W 1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     X 1.0,   0.98,  0.94,  0.845, 0.780, 0.753, 0.767, 0.765, 0.75/
      Data P13/1.0,  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     W 1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     X 1.0,   1.0,   1.0,   0.96,  0.865, 0.809, 0.79,  0.765, 0.75/
c
      Data N/26/, Nterp/1/
c
      En=E
      Pbar(1) = EXTERP(Eofp, P1,En,N,Nterp)
      Pbar(2) = EXTERP(Eofp, P2,En,N,Nterp)
      Pbar(3) = EXTERP(Eofp, P3,En,N,Nterp)
      Pbar(4) = EXTERP(Eofp, P4,En,N,Nterp)
      Pbar(5) = EXTERP(Eofp, P5,En,N,Nterp)
      Pbar(6) = EXTERP(Eofp, P6,En,N,Nterp)
      Pbar(7) = EXTERP(Eofp, P7,En,N,Nterp)
      Pbar(8) = EXTERP(Eofp, P8,En,N,Nterp)
      Pbar(9) = EXTERP(Eofp, P9,En,N,Nterp)
      Pbar(10)= EXTERP(Eofp,P10,En,N,Nterp)
      Pbar(11)= EXTERP(Eofp,P11,En,N,Nterp)
      Pbar(12)= EXTERP(Eofp,P12,En,N,Nterp)
      Pbar(13)= EXTERP(Eofp,P13,En,N,Nterp)
c
      Return
      END
c=============================================
