#!/usr/bin/python

import math
import numpy as np

import boot
import fit

T = 64 # temporal extend of the lattice

################################################################################
# reading bootstrapped data
C2 = np.load('bootdata/C2_p0.npy')
C4 = np.load('bootdata/C4_p0.npy')

################################################################################
################ computing the fits ############################################

# array for timeslices
tt = []
for t in range(0, T):
	tt.append(float(t))

# fitting the correlation function directly ####################################
print '\n**********************************************************'
print ' Fitting the mass directly from the correlation function.'
print '**********************************************************'
fitfunc1 = lambda p,x: p[0]*np.cosh((T/2.0-t)*p[1])

for lo in range(12, 13):
  up = C2.shape[1] # fit includes last time slice
  t = np.asarray(tt[lo:up])
  print "\nResult from bootstrapsample fit in interval (2pt):", lo, up
  res = fit.fitting(fitfunc1, t, C2[:,lo:up], [10., 0.19])


# fitting the four-point function ##############################################
print '\n*******************************************'
print ' Fitting the four-point function directly.'
print '*******************************************'
fitfunc1 = lambda p,x: p[0]*np.cosh((T/2.0-t)*p[1]) + p[2]

for lo in range(12, 13):
  up = C4.shape[1]
  t = np.asarray(tt[lo:up])
  print "\nResult from bootstrapsample fit in interval (2pt):", lo, up
  fit.fitting(fitfunc1, t, C4[:,lo:up], [10., 0.19, 10.])

















