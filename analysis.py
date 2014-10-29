#!/usr/bin/python

import math
import numpy as np

import matplotlib
matplotlib.use('QT4Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.mlab as mlab

import boot
import fit
import IOcontraction
import plot

T = 64 # temporal extend of the lattice

################################################################################
# reading bootstrapped data
C2 = np.load('bootdata/C2_p0.npy')
C4 = np.load('bootdata/C4_p0.npy')
# computing mean values and errors
C2_mean, C2_error = boot.mean_error_print(C2)
C4_mean, C4_error = boot.mean_error_print(C4)

################################################################################
################ computing the fits ############################################

# array for timeslices
tt = []
for t in range(0, T):
	tt.append(float(t))
tt = np.asarray(tt)

# fitting the correlation function directly ####################################
#print '\n**********************************************************'
#print ' Fitting the mass directly from the correlation function.'
#print '**********************************************************'
plot_path = './plots/test.pdf'
IOcontraction.ensure_dir(plot_path)
pdfplot = PdfPages(plot_path)

fitfunc1 = lambda p,x: p[0]*np.cosh((T/2.0-x)*p[1])
for lo in range(12, 13):
  up = C2.shape[1] # fit includes last time slice
  t = np.asarray(tt[lo:up])
  print "\nResult from bootstrapsample fit in interval (2pt):", lo, up
  res = fit.fitting(fitfunc1, t, C2[:,lo:up], [10., 0.19])
  res_mean = np.mean(res, axis=0)
  label = ['t', 'C2(t)', 'C2', '%.3f*cosh((%d-t)%.3f)' \
           % (res_mean[0], T/2, res_mean[1])]
  plot.corr_fct_with_fit(tt, C2_mean, C2_error, fitfunc1, res_mean, \
                         [5., T/2+1], label, pdfplot, 1)

# fitting the four-point function ##############################################
print '\n*******************************************'
print ' Fitting the four-point function directly.'
print '*******************************************'

fitfunc1 = lambda p,x: p[0]*np.cosh((T/2.0-x)*p[1]) + p[2]
for lo in range(12, 13):
  up = C4.shape[1]
  t = np.asarray(tt[lo:up])
  print "\nResult from bootstrapsample fit in interval (2pt):", lo, up
  res = fit.fitting(fitfunc1, t, C4[:,lo:up], [10., 0.19, 10.])
  res_mean = np.mean(res, axis=0)
  label = ['t', 'C4(t)', 'C4', '%.3f*cosh((%d-t)%.3f) + %.3f' \
           % (res_mean[0], T/2, res_mean[1], res_mean[2])]
  plot.corr_fct_with_fit(tt, C4_mean, C4_error, fitfunc1, res_mean, \
                         [5., T/2+1], label, pdfplot, 1)
pdfplot.close()
















