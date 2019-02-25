C   This is file SIGCELAS.FOR
c     Purpose is to return a cross section (in barns) for elastic scattering
c      of Carbon by an incident neutron of energy En (in MeV).  These data are
c     taken from the ENDF/B-V Evaluation for En up to about 16 MeV.
c     Other data for En > 16 MeV were obtained by subtracting total
c       non-elastic cross sections from total cross sections.
c
      FUNCTION SIGCELAS(En)
      Dimension E(541),Sig(541)
      Data NN/541/
      Data (E(I),I=1,173)/
     W  2.530E-08, 1.0E-07, 1.0E-05, 0.001, 0.005, 0.01,
     X  0.015,  0.02,   0.04,   0.05,   0.075,  0.10,   0.125,  0.225,
     X  0.3250, 0.4250, 0.5250, 0.6250, 0.7250, 0.8500, 1.0000, 1.1000,
     Y  1.2000, 1.3000, 1.4000, 1.5000, 1.6000, 1.7000, 1.8000, 1.9000,
     Z  2.0000, 2.0250, 2.0500, 2.0520, 2.0540, 2.0560, 2.0580, 2.0600,
     X  2.0610, 2.0620, 2.0630, 2.0640, 2.0650, 2.0660, 2.0670, 2.0680,
     X  2.0690, 2.0700, 2.0710, 2.0720, 2.0730, 2.0740, 2.0750, 2.0760,
     X  2.0770, 2.0780, 2.0790, 2.0800, 2.0810, 2.0820, 2.0830, 2.0840,
     X  2.0850, 2.0860, 2.0870, 2.0880, 2.0890, 2.0900, 2.0910, 2.0920,
     X  2.0930, 2.0940, 2.0950, 2.0960, 2.0980, 2.1000, 2.1020, 2.1040,
     A  2.1060, 2.1080, 2.1100, 2.1140, 2.1180, 2.1200, 2.1400, 2.1600,
     B  2.1800, 2.2000, 2.2400, 2.2800, 2.3200, 2.36,   2.4000, 2.4400,
     C  2.4800, 2.5200, 2.5600, 2.6000, 2.6400, 2.6600, 2.6800, 2.7000,
     X  2.7200, 2.7400, 2.7600, 2.7800, 2.7900, 2.8000, 2.8040, 2.8060,
     P  2.8080, 2.8100, 2.8110, 2.8120, 2.8130, 2.8140, 2.8150, 2.8160,
     Q  2.8170, 2.8180, 2.8190, 2.8200, 2.8220, 2.8240, 2.8260, 2.8300,
     R  2.835,  2.8400, 2.8500, 2.8600, 2.8650, 2.8700, 2.8750, 2.8800,
     S  2.8850, 2.8900, 2.8950, 2.9000, 2.9050, 2.9100, 2.9150, 2.9200,
     T  2.9280, 2.9320, 2.9360, 2.9400, 2.9440, 2.9480, 2.9520, 2.9560,
     U  2.9600, 2.9640, 2.9680, 2.9720, 2.9760, 2.9800, 2.9840, 2.9880,
     V  2.9920, 2.9960, 3.0000, 3.0100, 3.0200, 3.0300, 3.0400, 3.0600,
     W  3.0800, 3.1000, 3.1200, 3.14,   3.1600, 3.1800, 3.20/
      Data (E(I), I=174,377)/
     X  3.2200, 3.2600, 3.3000, 3.3200, 3.3400, 3.36,   3.3800, 3.4000,
     X  3.4200, 3.4400, 3.4600, 3.4800, 3.5000, 3.52,   3.5400, 3.56,
     Y  3.5800, 3.6000, 3.6400, 3.66,   3.6800, 3.7000, 3.7400, 3.7800,
     Z  3.8000, 3.8400, 3.9200, 3.9400, 3.9800, 4.0000, 4.02,
     a  4.0400, 4.0600, 4.0800, 4.1000, 4.1100,
     b  4.1200, 4.1300, 4.1400, 4.1500, 4.1600, 4.1700, 4.1800, 4.1900,
     c  4.2000, 4.2200, 4.2300, 4.2400, 4.2500, 4.2600, 4.2700, 4.2800,
     d  4.2900, 4.3000, 4.3100, 4.3300, 4.3500, 4.3700, 4.3900, 4.4000,
     e  4.4200, 4.4400, 4.4600, 4.4800, 4.5000, 4.5400, 4.5600, 4.6000,
     f  4.6400, 4.6800, 4.7000, 4.7400, 4.8000, 4.8200, 4.8500, 4.9000,
     K  4.9194, 4.9300, 4.9345, 4.9355, 4.9367, 4.9385, 4.9395, 4.9432,
     L  4.9481, 4.9546, 4.9584, 4.9698, 4.9800, 5.0000, 5.0300, 5.1000,
     M  5.1200, 5.1800, 5.2000, 5.2300, 5.2800, 5.3000, 5.3300, 5.3350,
     N  5.3400, 5.3600, 5.3620, 5.3710, 5.3780, 5.3900, 5.4000, 5.4100,
     O  5.4200, 5.4400, 5.4600, 5.5000, 5.5500, 5.5530, 5.6000, 5.6500,
     P  5.7,    5.8000, 5.9000, 6.0000, 6.0500, 6.1250, 6.1600, 6.1800,
     Q  6.2000, 6.2100, 6.2200, 6.2300, 6.2400, 6.2500, 6.2850, 6.2950,
     R  6.3030, 6.3100, 6.3200, 6.3300, 6.3400, 6.3500, 6.3600, 6.3700,
     S  6.3900, 6.4000, 6.4100, 6.4200, 6.4300, 6.4400, 6.4500, 6.4700,
     T  6.4900, 6.5100, 6.5400, 6.5600, 6.5700, 6.5800, 6.5900, 6.6000,
     U  6.6200, 6.6400, 6.6575, 6.6650, 6.6800, 6.7000, 6.7500, 6.8100,
     V  6.9200, 7.0000, 7.1000, 7.1400, 7.1800, 7.2000, 7.2200, 7.2250,
     W  7.2500, 7.2700, 7.3400, 7.3500, 7.3700, 7.4000, 7.4200, 7.4700,
     X  7.5417, 7.5937, 7.6200, 7.6500, 7.6674, 7.6800, 7.7000, 7.7250,
     Y  7.7450, 7.7500, 7.7700, 7.7887, 7.8100, 7.8191, 7.8600, 7.8884,
     Z  7.8971, 7.9300, 8.0000, 8.0440, 8.0800, 8.1000, 8.1050, 8.120/
      Data (E(I), I=378,521)/
     A  8.1380, 8.1660, 8.2000, 8.2100, 8.2400, 8.2800, 8.2960, 8.3200,
     B  8.3300, 8.4000, 8.4260, 8.4500, 8.5000, 8.5200, 8.6000, 8.6641,
     C  8.7000, 8.7681, 8.8000, 8.8500, 8.9198, 8.9400, 8.9800, 9.0000,
     D  9.005,  9.0200, 9.0300, 9.0450, 9.0800, 9.1490, 9.1625, 9.1800,
     E  9.2189, 9.2500, 9.2535, 9.3000, 9.3600, 9.4500, 9.5000, 9.5220,
     F  9.5600, 9.5900, 9.63,   9.6400, 9.6800, 9.6920, 9.7000, 9.7259,
     G  9.7400, 9.7500, 9.8000, 9.8299, 9.9000, 9.9209, 10.000, 10.050,
     H  10.170, 10.250, 10.300, 10.372, 10.400, 10.500, 10.550, 10.620,
     I  10.690, 10.830, 10.940, 11.000, 11.053, 11.100, 11.170, 11.250,
     J  11.400, 11.500, 11.700, 11.800, 11.900, 11.917, 12.000, 12.050,
     K  12.100, 12.224, 12.250, 12.300, 12.400, 12.500, 12.599, 12.70,
     L  12.990, 13.000, 13.100, 13.12,  13.250, 13.300, 13.540, 13.587,
     M  13.700, 13.822, 13.830, 13.965, 14.000, 14.182, 14.250, 14.419,
     N  14.500, 14.566, 14.694, 14.750, 14.767, 14.812, 14.837, 14.863,
     O  14.888, 14.927, 14.962, 15.000, 15.045, 15.093, 15.250, 15.477,
     P  15.731, 15.970, 16.000, 16.068, 16.256, 16.440, 16.695, 16.820,
     Q  16.974, 17.138, 17.300, 17.467, 17.687, 17.900, 18.087, 18.273,
     R  18.632, 18.833, 19.030, 19.185, 19.346, 19.511, 19.660, 20.00/
      Data (Sig(I),I=1,173)/4.7392,4.7392,4.7391,4.7346,4.7161, 4.6991,
     X  4.6821, 4.6653, 4.5989, 4.5662, 4.4862, 4.4084, 4.3326, 4.0491,
     A  3.7937, 3.5626, 3.3527, 3.1615, 2.9868, 2.7888, 2.5774, 2.4503,
     B  2.3331, 2.2249, 2.1250, 2.0328, 1.9479, 1.8698, 1.7981, 1.7321,
     C  1.6704, 1.6591, 1.6849, 1.6951, 1.7089, 1.7278, 1.7542, 1.7918,
     X  1.8165, 1.8467, 1.8836, 1.9294, 1.9868, 2.0595, 2.1525, 2.2736,
     X  2.4328, 2.6449, 2.9297, 3.3122, 3.8176, 4.4511, 5.1565, 5.7664,
     X  6.0435, 5.8719, 5.3728, 4.7654, 4.1989, 3.7282, 3.3559, 3.0658,
     X  2.8399, 2.6626, 2.5218, 2.4085, 2.3166, 2.2407, 2.1776, 2.1245,
     X  2.0794, 2.0407, 2.0073, 1.9783, 1.9303, 1.8926, 1.8624, 1.8376,
     A  1.8170, 1.7996, 1.7848, 1.7609, 1.7424, 1.7347, 1.6863, 1.6611,
     B  1.6445, 1.6320, 1.6136, 1.6004, 1.5911, 1.5853, 1.5829, 1.5843,
     C  1.5899, 1.6008, 1.6165, 1.6390, 1.6709, 1.6905, 1.7132, 1.7396,
     X  1.7702, 1.8059, 1.8478, 1.8979, 1.9270, 1.9607, 1.9786, 1.9903,
     P  2.0070, 2.0389, 2.0718, 2.1400, 2.3280, 3.1605, 5.0709, 2.7762,
     Q  2.2895, 2.1589, 2.1090, 2.0863, 2.0697, 2.0673, 2.0702, 2.0825,
     R  2.1029, 2.1265, 2.1795, 2.2425, 2.2774, 2.3149, 2.3555, 2.3994,
     S  2.4470, 2.4987, 2.5548, 2.6157, 2.6818, 2.7532, 2.8296, 2.9103,
     T  3.0432, 3.1067, 3.1633, 3.2070, 3.2297, 3.2212, 3.1707, 3.0686,
     U  2.9106, 2.7013, 2.4557, 2.1960, 1.9458, 1.7235, 1.5395, 1.3962,
     V  1.2908, 1.2176, 1.1704, 1.1295, 1.1503, 1.1976, 1.2547, 1.3725,
     W  1.4825, 1.5826, 1.6744, 1.76,   1.8404, 1.9171, 1.9904/
      Data (Sig(I),I=174,377)/
     X  2.0609, 2.1936, 2.3132, 2.3678, 2.4183, 2.4644, 2.5058, 2.5420,
     X  2.5729, 2.5984, 2.6182, 2.6324, 2.6410, 2.6442, 2.6422, 2.6353,
     Y  2.6238, 2.6081, 2.5654, 2.5393, 2.5105, 2.4795, 2.4120, 2.3394,
     Z  2.3019, 2.2257, 2.0736, 2.0367, 1.9657, 1.9319, 1.9004,
     a  1.8711, 1.8456, 1.8255, 1.8126, 1.8098,
     b  1.8103, 1.8145, 1.8233, 1.8372, 1.8572, 1.8839, 1.9175, 1.9579,
     c  2.0043, 2.1066, 2.1559, 2.1988, 2.2319, 2.2528, 2.2609, 2.2568,
     d  2.2422, 2.2195, 2.1909, 2.1240, 2.0534, 1.9854, 1.9226, 1.8933,
     e  1.8390, 1.7899, 1.7454, 1.7049, 1.6677, 1.6015, 1.5717, 1.5175,
     f  1.4688, 1.4245, 1.4036, 1.3639, 1.3078, 1.2877, 1.2600, 1.2230,
     g  1.2387, 1.3107, 1.3415, 1.7193, 1.9671, 1.9407, 1.7372, 1.4172,
     L  1.2776, 1.2301, 1.2072, 1.1818, 1.1720, 1.1583, 1.1448, 1.1201,
     M  1.1125, 1.0746, 1.0580, 1.0350,0.99397,0.97397,0.96297,0.97447,
     N  1.0360, 1.5038, 1.5516, 1.7101, 1.5508, 1.3542, 1.2424, 1.1506,
     O  1.1088, 1.0483, 1.0158,0.99097,0.98096,0.98012, 0.971, 0.94796,
     P 0.92596,0.90196,0.89546,0.88695,0.88595,0.89595,0.91222,0.94508,
     Q  1.0129, 1.0430, 1.0830, 1.1529, 1.2430, 1.4230, 2.1414, 2.2325,
     R  2.1251, 1.8340, 1.5599, 1.3796, 1.2490, 1.1593, 1.0766, 1.0023,
     S 0.91364, 0.8635, 0.8280, 0.8120,0.82164, 0.8280, 0.8390, 0.8535,
     T 0.83905, 0.7746,0.63544, 0.5551, 0.5216,0.50227,0.49794, 0.5236,
     U 0.60744,0.68294,0.71994,0.71294,0.67155, 0.6448, 0.6387, 0.6270,
     V  0.5963, 0.5696,0.56836, 0.5822, 0.6180,0.64442, 0.6829, 0.6913,
     W  0.7734, 0.8646, 1.3667, 1.4166, 1.4483, 1.4353, 1.4247, 1.4085,
     X  1.4056, 1.3920, 1.3775, 1.3966, 1.4427, 1.4692, 1.5413, 1.7184,
     Y  1.9003, 1.9170, 1.7906, 1.6925, 1.5909, 1.5687, 1.4917, 1.4642,
     Z  1.4517, 1.4123, 1.3597, 1.2955, 1.27765,1.2459, 1.2464, 1.223/
      Data (Sig(I),I=378,541)/
     A 1.16900,1.09276,0.99394,0.95854,0.88064,0.82764,0.81564,0.79411,
     B 0.78703,0.77581,0.7821, 0.77684,0.77778,0.7679, 0.75438,0.74355,
     C 0.7327, 0.72074,0.70652,0.68353,0.65027,0.63938,0.63623,0.64866,
     D 0.64878,0.68113,0.68937,0.67623,0.66755,0.70213,0.71266,0.72048,
     E 0.69418,0.69261,0.69549,0.75252,0.7096, 0.67788,0.6738, 0.66661,
     F 0.6542, 0.66894,0.63361,0.63353,0.68819,0.69409,0.69591,0.66616,
     G 0.64781,0.64182,0.64082,0.63852,0.62361,0.61916,0.62124,0.61361,
     H 0.56131,0.55811,0.57448,0.61246,0.62573,0.67312,0.68875,0.75263,
     I 0.78372,0.82552,0.88923,0.85707,0.88929,0.91773,0.85653,0.86087,
     J 0.893,  0.89375,0.85125,0.86114,0.90718,0.91882,0.9785, 1.0188,
     K 1.01411,0.91976,0.89953,0.87388,0.8526, 0.84632,0.85985,0.87126,
     L 0.91038,0.9069, 0.873,  0.86616,0.8512, 0.8448, 0.86513,0.84551,
     M 0.79844,0.80251,0.80274,0.80207,0.79686,0.78104,0.78273,0.78597,
     N 0.7946, 0.80265,0.83275,0.8548, 0.85801,0.88437,0.91573,0.92163,
     O 0.90352,0.87374,0.88093,0.89,   0.89651,0.90612,0.9129, 0.934,
     P 0.94242,0.9473, 0.9483, 0.9514, 0.9448, 0.93634,0.92169,0.9023,
     Q 0.8771, 0.8562, 0.8483, 0.8527, 0.8711, 0.8894, 0.9005, 0.9073,
     R 0.9144, 0.921,  0.9339, 0.952,  0.9823, 1.0076, 1.0103, 1.0055,
c     (Next 3 rows for E(neutron) > 20 MeV, as given in next data
c       statement.)
     S 0.9475, 0.930,  0.9205, 0.913,  0.885,  0.91,   0.938,  0.964,
     T 0.884,  0.830,  0.780,  0.730,  0.65,   0.605,  0.545,  0.5,
     U 0.46,   0.389,  0.325,  0.214/
      Data (E(I), I=522,541)/
     S 20.8,   22.0,   24.0,   26.0,   28.6,   29.0,   29.25,  29.59,
     T 30.0,   35.0,   40.0,   45.0,   50.0,   55.0,   60.0,   65.0,
     U 70.0,   80.0,   90.0,  110.0/
c
      Enn=En
      Sigma=EXTERP(E,Sig,Enn,NN,1)
      IF (Enn .LE. E(10)) Sigma=EXTERP(E,Sig,Enn,NN,4)
      Sigcelas=Sigma
      Return
      END
c
C   This is file SIGCINEL.FOR
c
c     Purpose is to return a cross section (in barns) for inelastic
c       scattering to the 4.4-MeV level in 12-C by an incident neutron
c       of energy En (in MeV).  Data for En < 15 MeV from ENDF/B-V
c     evaluation.
c
c     For E(neutron) = 20.8, 22.0, 24.0 and 26.0 MeV see Phys. Med.
c      Biol. 29 (1984) 643 for measurements; for En > 26 MeV see same
c      reference for nuclear model predictions.
c
      FUNCTION SIGCINEL(En)
      Dimension E(122),Sig(122)
      Data NN/122/, Nterp/1/
      Data E/ 4.812, 4.850, 4.900, 4.920, 4.930, 4.940, 4.9500, 4.9800,
     X  5.0000, 5.0300, 5.1000, 5.1200, 5.1500, 5.1800, 5.2000, 5.2300,
     X  5.2800, 5.3600, 5.3700, 5.3800, 5.4300, 5.5000, 5.5500, 5.6000,
     X  5.6500, 5.9000, 6.0500, 6.2000, 6.2500, 6.3200, 6.3400, 6.3500,
     X  6.3600, 6.3900, 6.4100, 6.4300, 6.4500, 6.5400, 6.5600, 6.6200,
     X  6.6400, 6.6700, 6.7500, 6.8100, 6.9200, 7.1400, 7.1800, 7.2200,
     X  7.2500, 7.3600, 7.4200, 7.4700, 7.5417, 7.5937, 7.6674, 7.7887,
     X  7.8191, 7.8884, 7.9361, 8.0000, 8.0140, 8.0440, 8.1000, 8.1380,
     X  8.1660, 8.2000, 8.2400, 8.3200, 8.4260, 8.5000, 8.7500, 8.8330,
     Y  9.0000, 9.0450, 9.1490, 9.2500, 9.5000, 9.6920, 9.7500, 10.000,
     Z  10.250, 10.500, 10.690, 10.750, 10.830, 11.000, 11.250, 11.500,
     A  11.750, 11.909, 12.000, 12.224, 12.599, 13.000, 13.250, 13.500,
     B  13.748, 14.000, 14.500, 14.750, 14.807, 14.863, 14.909, 14.954,
     C  16.443, 18.6,   19.5,   20.000, 20.8,   22.0,   24.0,   26.0,
     D  30.0,   35.0,   40.0,   45.0,   50.0,   60.0,   65.0,   70.000,
     E  80.0,   110.0/
      Data Sig/ 0.0, 0.008, 0.022, 0.028, 0.032, 0.035, 0.037, 0.047,
     X 0.048,  0.047,  0.038,  0.036,  0.040,  0.045,  0.052,  0.066,
     X 0.092,  0.148,  0.150,  0.149,  0.133,  0.124,  0.124,  0.124,
     X 0.1370, 0.1970, 0.2360, 0.2520, 0.2770, 0.340,  0.3510, 0.3490,
     X 0.340,  0.2880, 0.2650, 0.2550, 0.2520, 0.2720, 0.2640, 0.200,
     X 0.1870, 0.1750, 0.1620, 0.1560, 0.1560, 0.1670, 0.1750, 0.1920,
     X 0.2120, 0.310,  0.3510, 0.3480, 0.32156,0.30943,0.31376,0.3857,
     X 0.38917,0.35883,0.3597, 0.3770, 0.400,  0.450,  0.490,  0.460,
     X 0.430,  0.400,  0.390,  0.3450, 0.2660, 0.2450, 0.260,  0.2650,
     A 0.280,  0.310,  0.3140, 0.290,  0.280,  0.3220, 0.350,  0.3250,
     B 0.3150, 0.310,  0.330,  0.3450, 0.3650, 0.360,  0.3350, 0.270,
     C 0.260,  0.250,  0.2390, 0.2270, 0.2140, 0.210,  0.2050, 0.2030,
     D 0.1950 ,0.190,  0.170,  0.1650, 0.190,  0.2150, 0.1983, 0.1817,
     E 0.1630, 0.126,  0.1192, 0.1305, 0.102,  0.093,  0.079,  0.07,
     F 0.0577, 0.0475, 0.0397, 0.0331, 0.0275, 0.0195, 0.0166, 0.015,
     G 0.0126, 0.01/
c
      Enn=En
      Sigma=EXTERP(E,Sig,Enn,NN,Nterp)
      Sigcinel=Sigma
      Return
      END
c
C    This is file SIGCN2N.FOR
c
c      Purpose is to return a cross-section value (in barns) for the reaction
c       n + 12-C --> n + n + 11-C  for an incident neutron having energy
c     En (MeV).  Data for En to 34 MeV from Zeit. fur Physik A301 (1981)
c     353. See also Phys. Rev. 73 (1947) 262; ibid 265.
c
c     4/87. Add in "cross sections" to account for added capability of
c      the N2N routine to compute some (n,2np) reactions.  Threshold for
c      the 12-C(n,2np)10-B reaction is about 34.5 MeV.
c
      FUNCTION SIGCN2N(En)
      Dimension E(21),Sig(21)
      Data NN/21/, Nterp/1/
      Data E/  20.3,   22.0,  22.8,  23.9,  25.0,  26.0,  26.7,  28.0,
     U  30.0,  32.0,   36.3,  37.5,  40.0,  42.0,  45.0,  50.0,  56.0,
     W  60.0,  75.0,   90.0, 110.0/
      Data Sig/ 0.0,  0.002, 0.003,0.0063,0.0114,0.0139,0.0163,0.0207,
     W0.0235,0.0256,  0.028, 0.029,0.0313,0.0323, 0.033,0.0333,0.0325,
     X0.0316,0.0283,  0.022, 0.0145/
c
      Enn=En
      Sigma=EXTERP(E,Sig,Enn,NN,Nterp)
      Sigcn2n=Sigma
      Return
      END
c
C   This is file SIGCN3HE.FOR
c
c     Purpose to return a cross section (in barns) for  n + 12-C -->
c        3-He + 10-Be  reactions.  Data for En = 39.7 and 60.7 MeV
c        extracted from 3-He spectra published in Phys. Rev. C28, 521
c        (1983).  Data at En = 90 MeV for all 3-He reactions yield a
c        cross section of about 6 mb.
c
      FUNCTION SIGCN3HE(En)
      Dimension E(22), Sig(22)
      Data E/22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 33.5, 39.7, 42.5, 45.0,
     x 47.0, 49.0, 51.0, 53.0,  55.0,  57.5,  60.7, 65.0, 70.0, 76.0,
     y 90.0, 110.0/
c     Next array -- cross sections in mb, not barns.
      Data Sig/0.0, 0.5, 1.2,  2.4,  3.3,  4.6,  5.8,  10.9, 12.7, 14.0,
     x 14.8, 15.5, 15.9, 16.15, 16.25, 16.25, 16.0, 15.3, 14.0, 11.7,
     y 8.0, 6.5/
ckaji     y 6.0, 2.5/
      Data N/22/,Nterp/1/
c
      Enn=En
      S=EXTERP(E,Sig, Enn, N,Nterp)
      Sigcn3He=0.001*S
      Return
      END
c
C   This is file SIGCNAL.FOR
c
c     Purpose is to return a cross-section value (in barns) for the reaction
c      n + 12-C --> alpha + 9-Be (ground state)  by an incident neutron
c      of energy En (MeV).
c
      FUNCTION SIGCNAL(En)
      Dimension E(45),Sig(45)
      Data NN/ 45/, Nterp/1/
      Data E/ 6.186, 6.340, 7.180, 7.280, 7.340, 7.400, 7.5287, 7.6977,
     X  7.8971, 8.0054, 8.1008, 8.2958, 8.4995, 8.6641, 8.7681, 8.9198,
     X  9.1625, 9.2189, 9.2535, 9.3099, 9.7259, 9.8212, 9.9209, 10.372,
     X  11.004, 11.4,   12.5,   13.5,   14.5,   15.0,   17.0,   19.0,
     Y  22.0,   23.0,   24.0,   27.0,   30.0,   32.5,   40.0,   50.,
     Z  55.3,   62.7,   72.8,   90.0,  110.0/
      Data Sig/ 0.0, 1.0E-05, 0.001,0.006,0.011, 0.018,0.04391,0.09476,
     X  0.1589, 0.171,  0.1589, 0.0965, 0.05952,0.05952,0.07338,0.135,
     X  0.248,  0.293,  0.299,  0.286,  0.186,  0.178,  0.189,  0.1306,
     X  0.08263,0.0785, 0.074,  0.0716, 0.066,  0.062,  0.0392, 0.026,
     Y  0.0175, 0.0155, 0.0142, 0.0112, 0.0085, 0.0060, 0.0130, 0.0130,
     Z  0.0130, 0.0090, 0.0080, 0.0010, 0.00007/
ckajiZ  0.0029, 0.0027, 0.0025, 0.0010, 0.00007/
c     Y  0.0175, 0.0155, 0.0142, 0.0112, 0.0095, 0.0084, 0.00605,0.004,
c     Z  0.0034, 0.0029, 0.00215,0.0019, 0.0017/
c
      Enn=En
      IF (Enn .LT. E(3)) goto 2
      Sigma = EXTERP(E,Sig,Enn,NN,Nterp)
      Sigcnal=Sigma
      Return
    2 Sigma = EXTERP(E,Sig,Enn,NN,2)
      Sigcnal=Sigma
      Return
      END
c
C   This is file SIGCND.FOR
c
c      Purpose is to return a cross-section value (in barns) for the reaction
c     n + 12-C --> d + 11-B for an incident neutron energy En (in MeV).
c     This reaction can occur for neutrons having energies too small for the
c     n + 12-C --> n + p + 11-B reaction; the ENDF/B-V evaluation gives
c     cross sections up to 20 MeV based on indirect evidence (it appears).
c
c      Also tried to fit energy distributions at 27.4, 39.7, and 60.7 MeV
c      of Subramanian et al. Phys. Rev. C28, 521 (1983).  The fits are
c      quite reasonable.
c
c     (2006.01.26 by d.satoh)
c     Cross sections between 110 and 150 MeV were added based on the calculation
c     results of TALYS code. (n, d) + (n, nd)
c
      FUNCTION SIGCND(En)
      Dimension E(25),Sig(25)
      Data NN/25/
      Data E/15.25, 15.48, 15.97, 16.44, 16.97, 17.9,  18.46, 18.75,
     P       19.0,  19.51, 20.0,  21.0,  23.0,  26.5,  33.5,
c     Q       42.7,  53.0,  61.0,  70.0,  90.0, 110.0/
     Q       42.7,  53.0,  61.0,  70.0,  90.0,
     +       110.0, 120.0, 130.0, 140.0, 150.0/
      Data Sig/0.0, 0.002, 0.02,  0.03,  0.04,  0.06,   0.07, 0.071,
     P       0.068, 0.06,  0.053, 0.0505,0.049, 0.0485, 0.047,
c     Q       0.0449,0.042, 0.0395,0.0356,0.025, 0.016/
ckajiQ       0.0449,0.042, 0.0395,0.0356,0.025,
     Q       0.043,0.040, 0.0365,0.0350,0.030,
c_true     +       0.013994, 0.011622, 0.010256, 0.0086828, 0.0075769/
     +       0.015501, 0.0085761, 0.0075167, 0.0062766, 0.0054025/
ckaji     +       0.010501, 0.0085761, 0.0075167, 0.0062766, 0.0054025/
c
      Enn=En
      Sigma=EXTERP(E,Sig,Enn,NN,1)
      Sigcnd=Sigma
      Return
      End
c
c
C   This is file SIGCNN3A.FOR
c     Purpose is to return a cross section (in barns) for the reaction
c      n + 12-C --> n + 3 alphas for incident neutron of energy
c      En (in MeV).  For En < 13 MeV the sum of partial (n,n') cross
c      sections are being used.
c     It appears that values for 11 to 35 MeV gotton from Nuclear
c     Physics A394 (1983) 87 are too large, especially for En > 17 MeV.
c
      FUNCTION SIGCNN3A(En)
      Dimension E(31),Sig(31)
      Data NN/31/, Nterp/1/
      Data E/  8.4,  8.7,   9.0,   9.5,   10.0,  10.5,  11.0,  12.0,
     W 13.0,  13.5,  14.0,  14.5,  15.0,  16.0,  17.0,  18.0,  20.0,
     Y 23.0,  25.0,  27.0,  29.0,  31.0,  35.0,  40.0,  45.0,  50.0,
     Z 55.0,  60.0,  70.0,  90.0, 110.0/
      Data Sig/0.0,  0.006, 0.013, 0.035, 0.05,  0.061, 0.09,  0.15,
     W 0.205, 0.235, 0.265, 0.292, 0.306, 0.314, 0.308, 0.292, 0.258,
     X 0.21,  0.179, 0.152, 0.134, 0.117, 0.075, 0.030, 0.022, 0.018,
     Z 0.017, 0.017, 0.016, 0.007, 0.005/
c     X 0.21,  0.179, 0.152, 0.134, 0.117, 0.094, 0.073, 0.06,  0.048,
c     Z 0.04,  0.034, 0.0267,0.019, 0.012/
c
      Enn=En
      IF (Enn .LT. E(4)) goto 2
      Sigma=EXTERP(E,Sig,Enn,NN,Nterp)
      Sigcnn3a=Sigma
      Return
    2 Sigcnn3a=EXTERP(E,Sig,Enn,NN,2)
      Return
      END
c
C   This is file SIGCNP.FOR
c
c      Purpose is to return a cross-section value (in barns) for the reaction
c       n + 12-C --> p + 12-B (the particle stable levels only) for incident
c     neutron of energy En (MeV).
c
c     Cross sections between 15 and 20 MeV are essentially ENDF/B-V values
c      which, in turn, were taken from Rimmer and Fisher, Nuclear Physics
c      A108 (1968) 567.  Values at 27.4, 39.7, and 60.7 MeV adjusted
c      to satisfy spectral data of Subramanian et al. Physical Review C28,
c      521 (1983).  These adjustments took into consideration protons
c      generated in the other subroutines.
c
      FUNCTION SIGCNP(En)
      Dimension E(29),Sig(29)
      Data NN/29/
      Data E/13.665, 14.0,   14.5,   15.0, 15.477, 15.966, 16.443,
     Y  16.974,17.467,17.901,18.46,  19.034, 19.522, 20.0, 21.0,
     Z  22.0,  23.0,  25.0,  27.0,   28.4,   30.0,  34.8,  41.8,
     A  45.5,  55.3,  62.7,  72.8,   95.6,  110.0/
c     A  45.0,  55.0,  74.0, 110.0/
      Data Sig/0.0, 0.00014, 0.0004, 0.001,  0.004, 0.008, 0.011,
     Y  0.013, 0.016, 0.019, 0.019,  0.018,  0.015, 0.013, 0.0088,
     Z  0.0088,0.01,  0.0121,0.0133, 0.0135, 0.0133,0.012, 0.0108,
     A  0.0085,0.0075,0.0055,0.0030, 0.0020, 0.0015/
c     A  0.0098,0.0085,0.0068,0.005/
c
      Enn=En
      Sigma=EXTERP(E,Sig,Enn,NN,3)
      Sigcnp=Sigma
      Return
      END
c
C    This is file SIGCNPN.FOR
c
c      Purpose is to return a cross-section value (in barns) for the generic
c      reaction n + 12-C --> p + n + 11-B for incident neutrons having
c      energy En (in MeV).  It includes reactions in the program
c      following breakup of the 11-B when energetically available.
c      It also includes other reactions following breakup of 12-B when
c      energetically available, e.g.  n + 12-C --> p + p + 11-Be.
c     Datum at 90 MeV derived from Kellogg, Phys. Rev. 90, 224 (1953).
c
c      modified by daiki satoh. (2005.Aug)
c      The values above 80 MeV are assumed constant. -->  restore the data table (2005.Nov.)
c
      FUNCTION SIGCNPN(En)
      Dimension E(28),Sig(28)
      Data NN/28/, Nterp/1/
cc      Dimension E(24),Sig(24)
cc      Data NN/24/, Nterp/1/
c      Data E/ 17.35, 18.0,  18.7,  19.3,  20.0,  21.0,  22.0, 23.0,
c     W 24.5,  26.0,  28.0,  30.0,  33.0,  36.0,  40.0,  42.0, 45.0,
cc     D 47.0,  50.0,  53.0,  58.0,  65.0,  70.0,  80.0/
c     X 47.0,  50.0,  53.0,  58.0,  65.0,  70.0,  80.0,  85.0, 90.0,
c     Y 100.0, 110.0/
c      Data Sig/ 0.0, 0.005, 0.009, 0.012, 0.016, 0.02,  0.0235,0.0275,
c     W 0.0355, 0.04, 0.0465,0.05,  0.0545,0.0575,0.06,  0.0615,0.0625,
cc     D 0.0627,0.0625,0.062, 0.062, 0.0635,0.065, 0.074/
c     X 0.0627,0.0625,0.062, 0.062, 0.0635,0.065, 0.074, 0.083, 0.095,
c     Y 0.117, 0.141/
      Data E/ 17.35, 18.0,  18.7,  19.3,  20.0,  21.0,  22.0, 23.0,
     W 24.5,  26.0,  28.0,  30.0,  35.5,  38.5,  40.0,  41.8, 45.5,
     X 47.0,  50.0,  55.3,  62.7,  72.8,  80.0,  85.0, 90.0, 100.0,
     Y 110.0, 140.0/
      Data Sig/ 0.0, 0.005, 0.009, 0.012, 0.016, 0.02,  0.0235,0.0275,
     W 0.0355, 0.04, 0.0465,0.05,  0.0645,0.0710,0.0705,0.0690,0.0730,
ckajiW 0.0355, 0.04, 0.0465,0.05,  0.0545,0.0720,0.0705,0.0710,0.0705,
     X 0.0745,0.0760,0.0760,0.080, 0.090, 0.084, 0.087, 0.086, 0.102,
     Y 0.113, 0.125/
ck   X 0.0745,0.0770,0.082, 0.0855,0.110, 0.098, 0.083, 0.090,
ck   Y 0.098, 0.106/
ckajiX 0.0627,0.0625,0.062, 0.062, 0.0635,0.065, 0.074, 0.083, 0.095,
ckajiY 0.117, 0.141/
c
      Enn=En
      Sigma=EXTERP(E,Sig,Enn,NN,Nterp)
      Sigcnpn=Sigma
      Return
      END
c
C   This is file SIGCNT.FOR
c
      FUNCTION SIGCNT(En)
      Dimension  E(12),Sig(12)
c     Purpose to return cross section for (N,T) reaction.  Data values
c      for En = 27.4, 39.7, and 60.7 MeV deduced from graphs shown for
c      double differential cross sections in paper of Subramanian, et al,
c      Physical Review C28, 521 (1983).  Datum at 90 MeV from Kellogg's
c      paper, Phys. Rev. 90, 224 (l953) which includes all reactions
c      resulting in a triton.  And there are quite a few of them.
c
      Data E/ 21.5,22.5,25.0,27.4,30.0, 39.7, 45., 50., 54., 60.7, 90.0,
     k  110.0/
ckaji Data Sig/0.0, 1.0, 3.8, 7.1,10.5, 23.3, 28.5,30.7,31.3,30.5, 22.0,
      Data Sig/0.0, 1.0, 3.8, 7.1,10.5, 15.3, 16.0,18.8,19.0,21.3, 22.0,
     k  13.0/
c     Sig values are in millibarns, not barns as in the other functions.
      Data Ndata/12/, Nterp/1/
c
      Enn=En
      Sigma=EXTERP(E, Sig, Enn, Ndata, Nterp)
      Sigcnt=0.001*Sigma
      Return
      END
c
C    This is file SIGHYD.FOR
c      Purpose is to obtain cross-section value (in barns) for n+p scattering
c      by neutron of energy En in MeV.
c_4v
      FUNCTION SIGHYD(En)
      Dimension Ep(41),Sp(41),Ep2(47),Sp2(47)
      Data Ep/1.0E-10,1.5E-10,2.0E-10,2.5E-10,3.E-10,4.E-10,5.E-10,
     A  6.E-10,8.E-10,1.0E-9,1.5E-9,2.0E-9,2.5E-9,3.E-9,4.E-9,5.E-9,
     B  6.E-9,8.E-9,1.0E-8,1.5E-8,2.0E-8,2.5E-8,3.E-8,4.E-8,5.E-8,
     C  6.E-8,8.E-8,1.0E-7,1.5E-7,2.0E-7,2.5E-7,3.E-7,4.E-7,5.E-7,
     D  6.E-7,8.E-7,1.0E-6,1.5E-6,2.0E-6,2.5E-6,2.8E-6/
      Data Sp/  67.48,  50.58,  41.63,  35.67, 32.60, 26.04, 22.57,
     A  20.09, 16.82, 14.74, 11.67, 10.09, 8.897, 8.106, 7.018,6.326,
     B  5.883,5.242, 4.838, 4.208, 3.813, 3.547, 3.37, 3.124, 2.977,
     C  2.839,2.592, 2.608, 2.476, 2.388, 2.340, 2.301, 2.232, 2.193,
     D  2.168,2.134, 2.049, 2.046, 2.041, 2.037, 2.034/
cs
C modified by meigo (Add a datum at 200MeV by satoh)
      Data Ep2/2.000E+1,2.100E+1,2.200E+1,2.300E+1,2.40E+1,
     &2.500E+1,2.600E+1,2.700E+1,2.800E+1,2.90E+1,3.000E+1,3.200E+1,
     &3.400E+1,3.600E+1,3.800E+1,4.000E+1,4.20E+1,4.400E+1,4.600E+1,
     &4.800E+1,5.000E+1,5.200E+1,5.400E+1,5.60E+1,5.800E+1,6.000E+1,
     &6.200E+1,6.400E+1,6.600E+1,6.800E+1,7.00E+1,7.200E+1,7.400E+1,
     &7.600E+1,7.800E+1,8.000E+1,8.200E+1,8.40E+1,8.600E+1,8.800E+1,
     &9.000E+1,9.200E+1,9.400E+1,9.600E+1,9.80E+1,1.000E+2,2.000E+2/
      Data Sp2/4.801E-1,4.548E-1,4.331E-1,4.145E-1,3.97E-01,
     &3.804E-1,3.648E-1,3.501E-1,3.362E-1,3.232E-1,3.109E-1,2.884E-1,
     &2.683E-1,2.503E-1,2.341E-1,2.195E-1,2.064E-1,1.948E-1,1.843E-1,
     &1.750E-1,1.666E-1,1.591E-1,1.522E-1,1.458E-1,1.398E-1,1.342E-1,
     &1.289E-1,1.240E-1,1.195E-1,1.152E-1,1.112E-1,1.075E-1,1.040E-1,
     &1.008E-1,9.782E-2,9.502E-2,9.241E-2,8.996E-2,8.766E-2,8.549E-2,
     &8.343E-2,8.147E-2,7.959E-2,7.779E-2,7.607E-2,7.442E-2,4.205E-2/
      Data NN/41/,NN2/47/
      Data Pi/3.14159265/
c
      E=En
      S=20.34
      IF (En - Ep(NN)) 1,3,2
    1 S=10.*EXTERP(Ep,Sp,En,Nn,4)
      Goto 3
    2 IF (En-Ep2(1)) 10,20,20
   10 A=(-1.86+0.09415*E+0.000136*E*E)**2
      B=(0.4223+0.13*E)**2
      C=1.206*E
      S=3.*Pi/(C+A) + Pi/(C+B)
      go to 3
   20 S=EXTERP(Ep2,Sp2,En,Nn2,4)
    3 Sighyd=S
      Return
      End
c=============================================================
c   This is file SIGHYD2.FOR
c
c      Purpose is to return a cross-section value (in barn) for n+p scattering
c       by neutron of energy En in MeV.


c       Cross sections between 20 and 1000 MeV were taken from S.Chiba et al.,
c        J.Nucl.Sci.Technol.33(1996)654.


c       Cross sections above 1000MeV from Monning/Schopper ,LANDOLT-BORNSTEIN.
c        Vol.7.
c
c modified by satoh  '99.11.11


      FUNCTION SIGHYD2(En)
      Dimension E(25),Sig(25)
      Data NN/25/
      Data E/20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0,
     +  100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0,
     +  900.0, 1000.0, 1047.0, 1135.0, 1225.0, 1315.0, 1400.0,
     +  1453.0, 1600.0/
c
      Data Sig/0.4801, 0.3109, 0.2195, 0.1666, 0.1342, 0.1112,
     +  0.09502, 0.08343, 0.07442, 0.04205, 0.03491, 0.03269,
     +  0.03299, 0.03467, 0.03709, 0.03871, 0.03752, 0.03768,
     +  0.03830, 0.03890, 0.03950, 0.03990, 0.04240, 0.04040,
     +  0.03830/
c
      Enn=En
      Sigma=EXTERP(E,Sig,Enn,NN,1)
      Sighyd2=Sigma
      return
      end
c -------------------------------------------------------------------------
c     This is file SIGTOT.for


c        Purpose is to return a total cross section (in barn) of
c       carbon.
c        These value are calculated by Parlstein.




c       Cross sections above 1000MeV from Monning/Schopper ,LANDOLT-BORNSTEIN.
c        Vol.7.
c
c  modified by satoh  '99.11.11
cs+  modified by satoh  '99.12.03  (changing cross sections to Monning/Schopper
c                                     above 100Mev   :ordered by KENJI)


      function sigctot(en)
      dimension e(38),sig(38)


      data E/78.,82.0,84.0,88.0,92.0,96.0,98.0,100.0,
     +  101., 106., 110., 111., 117., 120., 121., 126.,
     +  129., 139., 140., 141., 151., 156., 169.,
     +  180., 190., 216., 220., 270., 280., 351., 380.,
     +  410., 500., 590., 630., 765., 1400., 2975./
cs+     +  120.0,129.0,150.0,200.0,250.0,300.0,350.0,400.0,450.0,
c     +  500.0,550.0,600.0,650.0,700.0,750.0,800.0,
c     +  850.0,900.0,950.0,1000.0, 1400.0, 2975.0/
      data sig/0.6333, 0.60250,0.58830,0.56180,0.53790,
     +  0.51620,0.50610,0.49660,0.472,0.456,0.439,0.435,0.408,
     +  0.403,0.392,0.390,0.375,0.355,0.349,0.346,0.341,0.330,
     +  0.323,0.311,0.291,0.296,0.297,0.288,0.279,0.2853,
     +  0.286,0.297,0.306,0.319,0.338,0.3421,0.378,0.369/
cs+     +  0.51620,0.50610,0.49660,0.403,0.375,0.36060,0.31870,0.30720,
c     +  0.30640,0.30930,0.31330,0.31720,0.32070,0.32370,
c     +  0.32620,0.32820,0.32980,0.33110,0.33210,0.33290,
c     +  0.33350,0.33400,0.33440, 0.3780, 0.3690/
c
      nn=38
      enn=en
      sigma=EXTERP(E,sig,enn,nn,1)
      sigctot=sigma
      return
      end
c -------------------------------------------------------------------------
************************************************************************
*                                                                      *
      subroutine sigrc(incp,emev,ia,iz,sigt,sigr,sigs)
*                                                                      *
*        calculates total, nonelastic and elastic cross-sections       *
*                                                                      *
*     input:                                                           *
*        incp   : =1, proton, =2, neutron                              *
*        emev   : incident nucleon energy (MeV)                        *
*        ia     : mass number                                          *
*        iz     : charge number                                        *
*                                                                      *
*     output:                                                          *
*       sigt    : total cross-section (b)                              *
*       sigr    : nonelastic cross-section (b)                         *
*       sigs    : sigt-sigr=elastic scattering cross-section (b)       *
*                                                                      *
************************************************************************


      implicit real*8 (a-h,o-z)
      dimension       par(4),c(4,3)


      data (c(1,k),k=1,3)/1.353,0.4993,1.076/
      data (c(2,k),k=1,3)/0.6653,0.9900,1.040/
      data (c(3,k),k=1,3)/2.456,2.060,0.9314/
      data (c(4,k),k=1,3)/0.1167,1.623,0.9736/


*-----------------------------------------------------------------------


         rz = 0.14


         a = float(ia)
         z = float(iz)


         a13 = a**0.3333333
         a23 = a**0.6666667
         alg = log(a)


         ecst = 82.0
         acst = 238.0


*-----------------------------------------------------------------------


         c1 =  0.0825
         c2 = -0.0057
         c3 =  0.14
         c4 = -0.2
         c5 =  2.72
         c6 =  1.62
         c7 = -5.3


*-----------------------------------------------------------------------


         ap = 1.0
         pi = 3.1415927


*-----------------------------------------------------------------------
*     modified by Niita
*-----------------------------------------------------------------------


            cmev = 0.0575 * a + 12.31


            facp = min( 1.0d0, 0.684 + 1.327e-3 * a )
            facc = 1.0 - facp


*-----------------------------------------------------------------------


            g1 = 1.0 - 0.62 * exp( -cmev / 200. )
     &                      * sin( 10.9 / cmev**0.28 )
            g1 = g1 * 1.1
            g2 = 1.0 + 0.016 * sin( 5.3 - 2.63 * log(a) )


            sigpa = 0.045 * a**0.7 * g1 * g2


*-----------------------------------------------------------------------


            eroot = sqrt(cmev)
            rad   = rz * a**0.3333
            wvl   = 0.1 * 1.22 * ( a + ap ) / a * sqrt(14.1) / eroot


            sigcr = pi * ( rad + wvl )**2


*-----------------------------------------------------------------------


            sigpc = facp * sigcr + facc * sigpa


            ftpa = sigpc / sigpa
            ftcr = sigpc / sigcr


            enpa = ftpa * 1.1 - 1.0


*-----------------------------------------------------------------------


         do k=1,4


            par(k) = log(c(k,1))+log(c(k,2))*alg+log(c(k,3))*alg**2
            par(k) = exp(par(k))


         end do


            epk   = par(3) * a13
            fexp1 = par(2) * log(epk/cmev)


            sigtp  = sigpc * ( 1.0 + par(4) )
     &             + par(1) * a13 * exp( -fexp1**2 )


            esub = ecst * a13 / acst**0.3333
            epk2 = epk - esub


         if( epk2 .gt. 0.0 ) then


            fexp2 = par(2) * log( epk2 / cmev )


            sigtp = sigtp + par(1) * a13 * exp( -fexp2**2 )


         end if


*-----------------------------------------------------------------------


            if( eroot .lt. 3.0 ) as = 1.0 - c1 * ( 3.0 - eroot )**2
            if( eroot .ge. 3.0 ) as = 1.0


            p = c2 * eroot + c3
            q = c4 * eroot + c5


            if( eroot .ge. 4.0 ) r = 1.2
            if( eroot .lt. 4.0 ) r = c6 * eroot + c7


            sigtc = 2.0 * sigpc * ( as - p * cos( q * a**0.3333 - r ) )


*-----------------------------------------------------------------------


            facp = min( 1.0d0, 0.578 + 1.77e-3 * a )
            facc = 1.0 - facp


            sigtpc = facp * sigtc + facc * sigtp


            fttp = sigtpc / sigtp
            fttc = sigtpc / sigtc


*-----------------------------------------------------------------------
*     calculation of sigr as
*     J.Letaw etal., Astrophys. J. Supp. series 51,271-276(1983).
*     sigrb = asymptotic high energy reaction cross-section
*-----------------------------------------------------------------------


         sigrb = 0.045 * a**0.7


*-----------------------------------------------------------------------
*     energy dependent factor
*-----------------------------------------------------------------------


         f1 = 1.0 - 0.62 * exp( -emev / 200. )
     &                   * sin( 10.9 / emev**0.28 )


*-----------------------------------------------------------------------
*     Pearlstein added low energy enhancement factor of 10%
*-----------------------------------------------------------------------


         f1 = f1 * ( 1.0 + enpa
     &      * exp( -min( 50.d0, ( emev - cmev ) / 10. )) )




c        f1 = f1 * ( 1.0 + 0.1 * exp( -( emev - 20. ) / 10. ) )


*-----------------------------------------------------------------------
*     mass dependent factor
*-----------------------------------------------------------------------


         f2 = 1.0 + 0.016 * sin( 5.3 - 2.63 * log(a) )


*-----------------------------------------------------------------------


         sigr = sigrb * f1 * f2


*-----------------------------------------------------------------------
*     sigt according to fit by S. Pearlstein, June 86.
*-----------------------------------------------------------------------


      if( emev .ge. cmev ) then


         do 210 k=1,4


            par(k) = log(c(k,1))+log(c(k,2))*alg+log(c(k,3))*alg**2
            par(k) = exp(par(k))


  210    continue


            epk   = par(3) * a13
            fexp1 = par(2) * log(epk/emev)


            sigt  = sigr * ( 1.0 + par(4) )
     &            + par(1) * a13 * exp( -fexp1**2 )


            ecst = 82.0
            acst = 238.
            esub = ecst * a13 / acst**0.3333
            epk2 = epk - esub


         if( epk2 .gt. 0.0 ) then


            fexp2 = par(2) * log( epk2 / emev )


            sigt = sigt + par(1) * a13 * exp( -fexp2**2 )


         end if


            factp = 1.0 + ( fttp - 1.0 )
     &            * exp( - min( 50.d0,( emev - cmev ) / 10.0 ))


            sigt = sigt * factp


*-----------------------------------------------------------------------
*     calculates xsects according to Ramsauer effect
*     ref: Angeli and Csikai, NP A170,577-583(1971)
*     sigt/signe=as-p*cos(q*a**0.3333-r)
*-----------------------------------------------------------------------


      else


            eroot = sqrt(emev)


            if( eroot .lt. 3.0 ) as = 1.0 - c1 * ( 3.0 - eroot )**2
            if( eroot .ge. 3.0 ) as = 1.0


            p = c2 * eroot + c3
            q = c4 * eroot + c5


            if( eroot .ge. 4.0 ) r = 1.2
            if( eroot .lt. 4.0 ) r = c6 * eroot + c7


            rad = rz * a**0.3333
            wvl = 0.1 * 1.22 * ( a + ap ) / a * sqrt(14.1) / eroot


         facr = 1.0 + ( ftcr - 1.0 )
     &        * exp( - min( 50.d0,( cmev - emev ) / 10.0 ))


            sigr = pi * ( rad + wvl )**2 * facr
            sigt = 2.0 * sigr * ( as - p * cos( q * a**0.3333 - r ) )


            factc = 1.0 + ( fttc - 1.0 )
     &            * exp( - min( 50.d0, ( cmev - emev ) / 10.0 ))


            sigt = sigt * factc


      end if


*-----------------------------------------------------------------------
*     coulomb factor
*     ec, coulomb barrier half-height, similar to s. pearlstein,
*     j. nuc. energy 23,87(1975) using optical model inverse p cs's.
*-----------------------------------------------------------------------


      if( incp .eq. 1 .and. emev .lt. 200.0 ) then


            ec = 1.44 * z / ( 8.2 + 0.68 * a13 )


            w1    = 3.816 + 0.1974 * z
            fexc1 = exp( max( -50.d0, ( ec - emev ) / w1 ))
            qfac1 = 1.0 / ( 1.0 + fexc1 )


            w2    = 0.07246 * z + 6.058
            if( z .lt. 10.0 ) w2 = 12.0
            qfac2 =  1.0 - exp( - min( 50.d0, ( emev / w2 )**2 ))


            w3    = 2.0
            fexc3 = exp( max( -50.d0, ( ec - emev ) / w3 ))
            qfac3 = 1.0 / ( 1.0 + fexc3 )


            fcoul = qfac1 * qfac2 * qfac3


            sigt = sigt * fcoul
            sigr = sigr * fcoul


      end if


*-----------------------------------------------------------------------


            sigs = sigt - sigr


*-----------------------------------------------------------------------


      return
      end
