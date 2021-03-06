c    This is file ETOL.FOR
c
c    Purpose is equivalence: light-unit value <---> energy (of charged
c    particle - proton, carbon ion, or alpha particle.)
c
      FUNCTION ETOL(N,A)
cccc modified by d.satoh (2006.07.07) cccc
c
c       N=1   Get proton energy (MeV) for A = proton light-unit value;
c       N=3   Get carbon energy (MeV) for A = carbon light-unit value;
c       N=5   Get alpha energy (MeV) for A = alpha light-unit value;
c       N=2   Get proton light-unit value for A = proton energy (MeV);
c       N=4   Get carbon light-unit value for A = carbon energy (MeV);
c       N=6   Get alpha light-unit value for A = alpha energy (MeV).
c       N=7   Get deuteron energy (MeV) for A = deuteron light unit value;
c       N=8   Get deuteron light-unit value for A = deuteron energy (MeV);
c
c       Tabulated data entered 6/87 -- estimated from measurements of
c       Bechetti, et al, Nucl. Instrum Methods 138 (1976) 93.
cccc
      Dimension X(127),Y(127)
      Common /LTABLE/ Itab,Ene(127),Hydl(127),Carl(127),Alpl(127),
     x Dlight(127)
      common /RES/ iswt, bias(5), icont, ictrl, iesc, ilight
c -------------------------------------------------------------------
      goto (987,999,987,987,987,999,987,999), N
c  999 if( ilight .eq. 1 .or. A .le. 1.e+1) goto 987  ! Database in SCINFUL
  999 if( ilight .eq. 1 ) goto 987  ! Database
c
c modified by d.satoh on June 2008
c
      En = A
      if( N .eq. 2 .and. ilight .eq. 3 .and. En .gt. 26.0 ) then
         a1 = 0.810                ! Satoh's parametrization for proton
         a2 = 2.430
         a3 = 0.290
      else if( N .eq. 2 .and. ilight .eq. 2 .and. En .gt. 10.0 ) then
         a1 = 0.810                ! Nakao's parametrization for proton
         a2 = 2.800
         a3 = 0.200
      else if( N .eq. 8 .and. ilight .eq. 3 .and. En .gt. 10.0 ) then
         a1 = 0.740                ! Satoh's parametrization for deuteron
         a2 = 3.450
         a3 = 0.200
      else if( N .eq. 8 .and. ilight .eq. 2 .and. En .gt. 1.0 ) then
         a1 = 0.750                ! Nakao's parametrization for deuteron
         a2 = 4.500
         a3 = 0.160
      else if( N .eq. 6 .and. ilight .eq. 3 .and. En .gt. 5.0 ) then
         a1 = 0.510                ! Satoh's parametrization for alpha
         a2 = 6.420
         a3 = 0.080
      else
              goto 987
      end if
c
      ETOL = a1*En - a2*(1.0-exp(-1*a3*En))     ! light output function (MeVee)
      ETOL = ETOL / 1.25      ! Na-unit
      return
c -------------------------------------------------------------------
  987 continue
c
      DATA ITAB/127/
      DATA (HYDL(I),I=1,127)/0.0,2.150E-4,1.030E-3,1.150E-3,1.270E-3,
     11.400E-3,1.460E-3,1.520E-3,1.580E-3,1.640E-3,1.760E-3,1.880E-3,
     22.000E-3,2.120E-3,2.240E-3,2.360E-3,2.480E-3,2.600E-3,2.720E-3,
     32.850E-3,3.090E-3,3.340E-3,3.570E-3,3.820E-3,4.070E-3,4.460E-3,
     44.840E-3,5.220E-3,5.620E-3,6.150E-3,6.710E-3,7.400E-3,8.120E-3,
     58.860E-3,9.620E-3,1.040E-2,1.122E-2,1.207E-2,1.292E-2,1.377E-2,
     61.465E-2,1.645E-2,1.838E-2,2.035E-2,2.245E-2,2.460E-2,2.680E-2,
     72.900E-2,3.130E-2,3.380E-2,3.650E-2,4.220E-2,4.830E-2,5.450E-2,
     86.090E-2,6.780E-2,7.860E-2,9.100E-2,1.040E-1,1.175E-1,1.367E-1,
     91.562E-1,1.825E-1,2.095E-1,2.385E-1,2.690E-1,2.995E-1,3.320E-1,
     A3.660E-1,4.000E-1,4.360E-1,4.725E-1,5.480E-1,6.250E-1,7.030E-1,
     B7.830E-1,8.660E-1,9.520E-1,1.042E+0,1.135E+0,1.230E+0,1.327E+0,
     C1.521E+0,1.718E+0,1.915E+0,2.112E+0,2.310E+0,2.630E+0,2.950E+0,
     D3.280E+0,3.620E+0,4.080E+0,4.550E+0,5.140E+0,5.750E+0,6.360E+0,
     E6.970E+0,7.580E+0,8.200E+0,8.830E+0,9.480E+0,1.014E+1,1.080E+1,
     F 12.15,   13.5,    14.9,    16.3,    17.7,    19.1,    20.5,
     G 21.9,    23.35,   24.8,    27.7,    30.75,   33.8,    36.9,
     H 40.0,    44.65,   49.4,    54.3,    59.2,    65.8,    72.4,
     S         397.03,          802.98,           965.36/
c       Light values for E(proton) between 0.1 and 40 MeV taken from the
c          O5S tabulation.  Same energy region for light values for
c          E(Carbon) also from O5S tabulation.  Alpha light, however,
c          was adjusted following comparisons with measured responses
c          during program development.
      DATA (CARL(I),I=1,127)/.0,6.0E-5,.000228,.00025,.000272,.000294,
     1.000305,.000316,.000327,.000338,.000360,.000381,.000402,.000422,
     2.000443,.000462,.000482,.000500,.000519,.000537,.000571,.000607,
     3.000640,.000673,.000707,.000756,.000806,.000856,.000905,.000972,
     4.001038,.001115,.001192,.001270,.001347,.001424,.001501,.001573,
     5.001645,.001717,.001788,.001932,.002076,.002219,.002363,.002506,
     6.002650,.002793,.002926,.003058,.003191,.003433,.003676,.003919,
     7.004140,.004361,.004692,.005023,.005354,.005686,.006127,.006569,
     8.007093,.007604,.008128,.008639,.009163,.009660,.010157,.010667,
     9.011150,.011647,.012641,.013634,.014628,.015622,.016615,.017664,
     A.018713,.019762,.020810,.021859,.023957,.026054,.028152,.030250,
     B.032347,.035549,.038750,.041952,.045154,.049570,.053986,.059616,
     C.065412,.071346,.077694,.084594,.091632,.098808,.105984,.113712,
     D.121440,.136896,.153456,.170016,.187680,.206448,.225216,.246192,
     E.268272,.290352,.312432,.357696,.406272,.460368,.516672,.575184,
     F.670128,.775008,.888720,1.00464,1.18680,1.36896,
     S        1.36869,1.36869,1.36896/
      DATA (ALPL(I),I=1,127)/0.0,6.500E-5,2.850E-4,3.210E-4,3.550E-4,
     13.880E-4,4.040E-4,4.200E-4,4.360E-4,4.520E-4,4.840E-4,5.160E-4,
     25.480E-4,5.800E-4,6.100E-4,6.400E-4,6.700E-4,7.000E-4,7.300E-4,
     37.580E-4,8.140E-4,8.700E-4,9.260E-4,9.820E-4,1.040E-3,1.130E-3,
     41.220E-3,1.310E-3,1.400E-3,1.520E-3,1.640E-3,1.790E-3,1.940E-3,
     52.090E-3,2.240E-3,2.400E-3,2.560E-3,2.720E-3,2.880E-3,3.040E-3,
     63.200E-3,3.520E-3,3.860E-3,4.200E-3,4.540E-3,4.900E-3,5.270E-3,
     75.640E-3,6.010E-3,6.380E-3,6.750E-3,7.510E-3,8.300E-3,9.120E-3,
     89.960E-3,1.080E-2,1.210E-2,1.350E-2,1.500E-2,1.656E-2,1.876E-2,
     92.100E-2,2.400E-2,2.700E-2,3.020E-2,3.430E-2,3.880E-2,4.360E-2,
     A4.860E-2,5.390E-2,5.940E-2,6.510E-2,7.730E-2,9.040E-2,1.044E-1,
     B1.193E-1,1.350E-1,1.520E-1,1.693E-1,1.876E-1,2.070E-1,2.270E-1,
     C2.690E-1,3.150E-1,3.640E-1,4.155E-1,4.705E-1,5.585E-1,6.530E-1,
     D7.545E-1,8.620E-1,1.016E+0,1.180E+0,1.410E+0,1.660E+0,1.910E+0,
     E2.200E+0,2.460E+0,2.770E+0,3.070E+0,3.390E+0,3.710E+0,4.000E+0,
     F4.700E+0,5.400E+0,6.140E+0,6.820E+0,7.600E+0,8.300E+0,9.100E+0,
     G1.000E+1,1.080E+1,1.160E+1,1.300E+1,1.470E+1,1.610E+1,1.790E+1,
     H 19.5,    21.8,    24.2,    26.5,    28.8,    31.9,    34.4,
     S        185.51,  374.02,  449.42/
cs      Data (Dlight(i),i=1,124)/0.000E+0,5.585E-5,6.144E-4,7.386E-4,8.627E-4,
c_4v
      Data Dlight/0.000E+0,5.585E-5,6.144E-4,7.386E-4,8.627E-4,
     &9.868E-4,1.049E-3,1.111E-3,1.173E-3,1.235E-3,1.359E-3,1.483E-3,
     &1.608E-3,1.732E-3,1.856E-3,1.980E-3,2.104E-3,2.228E-3,2.352E-3,
     &2.476E-3,2.725E-3,2.973E-3,3.221E-3,3.470E-3,3.718E-3,4.090E-3,
     &4.463E-3,4.835E-3,5.207E-3,5.704E-3,6.200E-3,6.821E-3,7.442E-3,
     &8.063E-3,8.684E-3,9.307E-3,9.934E-3,1.057E-2,1.120E-2,1.185E-2,
     &1.250E-2,1.384E-2,1.523E-2,1.665E-2,1.812E-2,1.964E-2,2.121E-2,
     &2.283E-2,2.451E-2,2.624E-2,2.802E-2,3.175E-2,3.569E-2,3.986E-2,
     &4.425E-2,4.885E-2,5.616E-2,6.395E-2,7.220E-2,8.090E-2,9.318E-2,
     &1.062E-1,1.235E-1,1.418E-1,1.611E-1,1.814E-1,2.024E-1,2.243E-1,
     &2.470E-1,2.704E-1,2.944E-1,3.192E-1,3.705E-1,4.243E-1,4.802E-1,
     &5.383E-1,5.984E-1,6.605E-1,7.245E-1,7.903E-1,8.577E-1,9.268E-1,
     &1.069E+0,1.218E+0,1.371E+0,1.530E+0,1.693E+0,1.945E+0,2.206E+0,
     &2.475E+0,2.751E+0,3.130E+0,3.520E+0,4.021E+0,4.536E+0,5.063E+0,
     &5.602E+0,6.151E+0,6.710E+0,7.278E+0,7.853E+0,8.437E+0,9.027E+0,
     &1.023E+1,1.145E+1,1.270E+1,1.396E+1,1.524E+1,1.654E+1,1.785E+1,
     &1.918E+1,2.051E+1,2.186E+1,2.459E+1,2.735E+1,3.015E+1,3.297E+1,
     &3.582E+1,4.014E+1,4.451E+1,4.892E+1,5.337E+1,5.935E+1,6.539E+1,
     S 355.28, 718.10, 863.22/
c
       DATA ENE /  0.000,  0.001,  0.010,  0.012,  0.014,
     1  0.016,  0.017,  0.018,  0.019,  0.020,  0.022,  0.024,  0.026,
     2  0.028,  0.030,  0.032,  0.034,  0.036,  0.038,  0.040,  0.044,
     3  0.048,  0.052,  0.056,  0.060,  0.066,  0.072,  0.078,  0.084,
     4  0.092,  0.100,  0.110,  0.120,  0.130,  0.140,  0.150,  0.160,
     5  0.170,  0.180,  0.190,  0.200,  0.220,  0.240,  0.260,  0.280,
     6  0.300,  0.320,  0.340,  0.360,  0.380,  0.400,  0.440,  0.480,
     7  0.520,  0.560,  0.600,  0.660,  0.720,  0.780,  0.840,  0.920,
     8  1.000,  1.100,  1.200,  1.300,  1.400,  1.500,  1.600,  1.700,
     9  1.800,  1.900,  2.000,  2.200,  2.400,  2.600,  2.800,  3.000,
     &  3.200,  3.400,  3.600,  3.800,  4.000,  4.400,  4.800,  5.200,
     &  5.600,  6.000,  6.600,  7.200,  7.800,  8.400,  9.200, 10.000,
     & 11.000, 12.000, 13.000, 14.000, 15.000, 16.000, 17.000, 18.000,
     & 19.000, 20.000, 22.000, 24.000, 26.000, 28.000, 30.000, 32.000,
     & 34.000, 36.000, 38.000, 40.000, 44.000, 48.000, 52.000, 56.000,
     & 60.000, 66.000, 72.000, 78.000, 84.000, 92.000,100.000,
     S 500.00,       1000.000,       1200.000/
c
      Aa=A
      M=N
      L=Itab
      IF (M .GT. 0) goto (1,9,1,9,1,9,1,9,25),M
   25 Etol = 0.0
      write(*,40)M
   40 Format(/'  *** Error Function ETOL, Input N ='I12//)
      Write (21,40) M
      Return
c
    1 Do 2 I=1,L
    2 Y(I)=Ene(I)
      IF (M - 3) 3,5,7
    3 Do 4 I=1,L
    4 X(I)=Hydl(I)
      Goto 19
    5 Do 6 I=1,L
    6 X(I)=Carl(I)
      Goto 19
    7 Do 8 I=1,L
      IF (M.EQ.5) X(I)=Alpl(I)
      IF (M.EQ.7) X(I)=Dlight(I)
    8 Continue
      Goto 19
    9 Do 10 I=1,L
   10 X(I)=Ene(I)
      IF (M - 4) 11,13,15
   11 Do 12 I=1,L
   12 Y(I)=Hydl(I)
      Goto 19
   13 Do 14 I=1,L
   14 Y(I)=Carl(I)
      Goto 19
   15 IF (M.EQ.8) goto 17
      Do 16 I=1,L
   16 Y(I)=Alpl(I)
      Goto 19
   17 Do 18 I=1,L
   18 Y(I)=Dlight(I)
c
   19 If (Aa .LT. X(3)) goto 20
      Etol=EXTERP(X,Y,Aa,L,4)
      Return
   20 Etol=EXTERP(X,Y,Aa,L,1)
      Return
      END
