#!/usr/bin/env python3

from math import *

emin = 1  # MeV
emax = 1000  # MeV
nmax = 30

# log axis
# dltloge = (log10(emax) - log10(emin))/(nmax - 1)
dltloge = (log10(emax) - log10(emin))/(nmax)

energy0 = []
energy1 = []

f0 = open('bin.inp', 'w')
for n in range(nmax):
    energy0.append(10 ** (log10(emin) + dltloge*n))
    f0.write('{0}\n'.format(energy0[-1]))
    print('{0}'.format(energy0[-1]))
f0.close()

f1 = open('energy.inp', 'w')
for i in range(nmax-1):
    energy1.append(float((energy0[i+1]+energy0[i])/2))
    f1.write('{0}\n'.format(energy1[-1]))
f1.close()

'''
# liner axis
delte = float((emax-emin)/(nmax-1))
print(delte)

energy0 = []
energy1 = []

f0 = open('bin.inp', 'w')
for n in range(nmax):
    energy0.append(emin + n*delte)
    f0.write('{0}\n' .format(energy0[-1]))
    print('{0}' .format(energy0[-1]))
f0.close()

f1 = open('energy.inp', 'w')
for i in range(nmax-1):
    energy1.append(float((energy0[i+1]+energy0[i])/2))
    f1.write('{0}\n'.format(energy1[-1]))
f1.close()
'''
