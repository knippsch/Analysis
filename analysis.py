#!/usr/bin/python

import math
import numpy as np
import scipy

import matplotlib
matplotlib.use('QT4Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.mlab as mlab

import boot
import fit
import IOcontraction
import plot
import zeta

T = 48 # temporal extend of the lattice

################################################################################
################ computing the fits ############################################

# array for timeslices
tt = []
for t in range(0, T):
	tt.append(float(t))
tt = np.asarray(tt)

# fitting the correlation function directly ####################################
print '\n**********************************************************'
print ' Fitting the mass directly from the correlation function.'
print '**********************************************************'
plot_path = './plots/C2.pdf'
IOcontraction.ensure_dir(plot_path)
pdfplot = PdfPages(plot_path)

for p in range(0,1):
  filename = 'bootdata/C2_p%d.npy' % p
  print 'reading file:', filename
  C2 = np.load(filename)
  C2_mean, C2_error = boot.mean_error_print(C2)
  m, m_mean, m_error = boot.compute_mass(C2) 
  print 'The correlation function:\n'
  i = 0
  for a, b in zip(C2_mean, C2_error):
    print i, a, '+/-', b
    i += 1

  fitfunc1 = lambda p,t,tt: p[0]*(np.exp(-t*p[1]) + np.exp(-(T-t)*p[1]))
  fitfunc2 = lambda p,t: p[0]
  counter = 0
  for lo in range(13, C2.shape[1]):
    for up in range(C2.shape[1], 13, -1):
      if (up-lo-2) < 2: # only fits with dof>1 are ok
        continue
      t = np.asarray(tt[lo:up])
      print "\nResult from bootstrapsample fit in interval (2pt):", lo, up-1
      start_parm = [C2_mean[lo], m_mean[lo]]
      res_C2, chi2, pvalue = fit.fitting(fitfunc1, t, C2[:,lo:up], start_parm)
      res_C2_mean = np.mean(res_C2, axis=0)
      label = ['t', 'm_eff(t)', 'data', \
               'fit, tmin = %d, tmax = %d' % (lo, up-1)]
      plot.corr_fct_with_fit(tt, m_mean, m_error, fitfunc2, [res_C2_mean[1]], \
                             [7., m.shape[1]], label, pdfplot, 0)
      # collecting fitresults
      if counter is 0:
        fitresult = res_C2[:,1]
        fitdata = np.array([chi2, pvalue, lo, up-1])
        counter += 1
      else:
        fitresult = np.vstack((fitresult, res_C2[:,1]))
        fitdata = np.vstack((fitdata, np.array([chi2, pvalue, lo, up-1])))
        counter += 1
  # save fitted masses on disk
  path = 'bootdata/C2_massfit_p%d' % p
  IOcontraction.ensure_dir(path)
  np.save(path, fitresult)
  # save fitted masses on disk
  path = 'bootdata/C2_massfit_params_p%d' % p
  IOcontraction.ensure_dir(path)
  np.save(path, fitdata)

pdfplot.close()
################################################################################
################################################################################
################################################################################
#print '\n**********************************************'
#print ' Fitting the ratio of 4pt- and 2pt- functions.'
#print '**********************************************'
#plot_path = './plots/C4_ratio.pdf'
#IOcontraction.ensure_dir(plot_path)
#pdfplot = PdfPages(plot_path)
#
#tt1 = []
#for t in range(0, T-1):
#	tt1.append(float(t) + 0.5)
#
#def fitfunc(p, t, tt):
#  r = np.empty(len(t), dtype=float)
#  for i in range(0, len(t)):
#    if i < tt:
#      r[i] = p[0]*(np.cosh(p[1]*(t[i]-T*0.5)) + np.sinh(p[1]*(t[i]-T*0.5)) * \
#                   np.cosh(2*p[2]*(t[i]-T*0.5)) / \
#                   np.sinh(2*p[2]*(t[i]-T*0.5)))
#    else:
#      r[i] = p[3]*(np.exp(-t[i]*p[2]) + np.exp(-(T-t[i])*p[2]))
#  return r
#
#for p in range(0,1):
#  filename = 'bootdata/C4_p%d.npy' % p
#  C4 = np.load(filename)
#  C4, C4_mean, C4_error = boot.compute_ratio(C4, C2)
#  i = 0
#  for a, b in zip(C4_mean, C4_error):
#    print i, a, '+/-', b
#    i += 1
#
#  fitfunc1 = lambda p,t: p[0]*(np.cosh(p[1]*(t-T*0.5)) +\
#                               np.sinh(p[1]*(t-T*0.5)) * \
#                               np.cosh(2*p[2]*(t-T*0.5)) / \
#                               np.sinh(2*p[2]*(t-T*0.5)))
#  counter = 0
#  for loC2 in range(13, 15):
#    for upC2 in range(C2.shape[1], C2.shape[1]-1, -1):
#      if (upC2-loC2-2) < 2: # only fits with dof>1 are ok
#        continue
#      B = C2[:,loC2:upC2]
#      tC2 = np.asarray(tt[loC2:upC2])
#
#      for lo in range(13, 15):
#        for up in range(C4.shape[1], C4.shape[1]-1, -1):
#          if (up-lo-2) < 2: # only fits with dof>1 are ok
#            continue
#          t = np.asarray(tt1[lo:up])
#          A = C4[:,lo:up]
#          data = np.hstack((A, B))
#          ttt = np.hstack((t, tC2))
#          print "\nResult from bootstrapsample fit in interval (2pt):", lo, up-1
#          start_parm = [1., 0.006, m_mean[loC2], C4_mean[loC2]]
#          # FIT
#          res_C4, chi2, pvalue = fit.fitting(fitfunc, ttt, data, \
#                                             start_parm, correlated = 1, \
#                                             tcut = len(t))
#          res_C4_mean = np.mean(res_C4, axis=0)
#          label = ['t', 'ratio(t)', 'data', \
#                   'fit, \nC4:tmin = %d, tmax = %d \nC2: tmin = %d, tmax = %d' \
#                   % (lo, up-1, loC2, upC2-1)]
#          # PLOT
#          plot.corr_fct_with_fit(tt1, C4_mean, C4_error, fitfunc1, \
#                                 res_C4_mean[0:3], \
#                                 [7., C4.shape[1]], label, pdfplot, 0)
#          # collecting fitresults
#          if counter is 0:
#            fitresult = res_C4[:,1]
#            fitdata = np.array([chi2, pvalue, lo, up-1, loC2, upC2-1])
#            counter += 1
#          else:
#            fitresult = np.vstack((fitresult, res_C4[:,1]))
#            fitdata = np.vstack((fitdata, np.array([chi2, pvalue, lo, up-1, \
#                                 loC2, upC2-1])))
#            counter += 1
#      # save fitted masses on disk
#      path = 'bootdata/C4_ratio_massfit_p%d' % p
#      IOcontraction.ensure_dir(path)
#      np.save(path, fitresult)
#      # save fitted masses on disk
#      path = 'bootdata/C4_ratio_massfit_params_p%d' % p
#      IOcontraction.ensure_dir(path)
#      np.save(path, fitdata)
#
#pdfplot.close()
#
##################################################################################
##################################################################################
### fitting the four-point function ##############################################
#print '\n*******************************************'
#print ' Fitting the four-point function directly.'
#print '*******************************************'
#plot_path = './plots/C4_direct.pdf'
#IOcontraction.ensure_dir(plot_path)
#pdfplot = PdfPages(plot_path)
#
#for p in range(0,1):
#  filename = 'bootdata/C4_p%d.npy' % p
#  C4 = np.load(filename)
#  C4_mean, C4_error = boot.mean_error_print(C4)
#  m4, m4_mean, m4_error = boot.compute_mass(C4) 
#  i = 0
#  for a, b in zip(C4_mean, C4_error):
#    print i, a, '+/-', b
#    i += 1
#
#  fitfunc1 = lambda p,t,tt: p[0]*(np.exp(-t*p[1]) + np.exp(-(T-t)*p[1])) + p[2]
#  fitfunc2 = lambda p,t: p[0]*(np.exp(-t*p[1]) + np.exp(-(T-t)*p[1])) + p[2]
#  counter = 0
#  for lo in range(13, 15):
#    for up in range(C4.shape[1], C4.shape[1]-2, -1):
#      if (up-lo-3) < 2: # only fits with dof>1 are ok
#        continue
#      t = np.asarray(tt[lo:up])
#      print "\nResult from bootstrapsample fit in interval (2pt):", lo, up-1
#      start_parm = [C4_mean[lo]**2, m4_mean[lo], C4_mean[C4.shape[1]-1]]
#      res_C4, chi2, pvalue = fit.fitting(fitfunc1, t, C4[:,lo:up], \
#                                         start_parm, correlated = 1)
#      res_C4_mean = np.mean(res_C4, axis=0)
#      label = ['t', 'C4(t)', 'data', \
#               'fit, tmin = %d, tmax = %d' % (lo, up-1)]
#      plot.corr_fct_with_fit(tt, C4_mean, C4_error, fitfunc2, res_C4_mean, \
#                             [7., C4.shape[1]], label, pdfplot, 1)
#      # collecting fitresults
#      if counter is 0:
#        fitresult = res_C4[:,1]
#        fitdata = np.array([chi2, pvalue, lo, up-1])
#        counter += 1
#      else:
#        fitresult = np.vstack((fitresult, res_C4[:,1]))
#        fitdata = np.vstack((fitdata, np.array([chi2, pvalue, lo, up-1])))
#        counter += 1
#  # save fitted masses on disk
#  path = 'bootdata/C4_direkt_massfit_p%d' % p
#  IOcontraction.ensure_dir(path)
#  np.save(path, fitresult)
#  # save fitted masses on disk
#  path = 'bootdata/C4_direkt_massfit_params_p%d' % p
#  IOcontraction.ensure_dir(path)
#  np.save(path, fitdata)
#
#pdfplot.close()
#################################################################################
#################################################################################
#################################################################################
#print '\n***********************************************'
#print ' Fitting the derivative of four-point function.'
#print '***********************************************'
#plot_path = './plots/C4_derivative.pdf'
#IOcontraction.ensure_dir(plot_path)
#pdfplot = PdfPages(plot_path)
#
#for p in range(0,1):
#  filename = 'bootdata/C4_p%d.npy' % p
#  C4 = np.load(filename)
#  C4, C4_mean, C4_error = boot.compute_derivative(C4)
#  m4, m4_mean, m4_error = boot.compute_mass(C4) 
#  i = 0
#  for a, b in zip(C4_mean, C4_error):
#    print i, a, '+/-', b
#    i += 1
#
#  fitfunc1 = lambda p,t,tt: -p[0]*np.exp(-p[1]*(t+T+1.)) * \
#                         (-1.+np.exp(p[1])) * \
#                         (np.exp(p[1]+2.*p[1]*t)-np.exp(p[1]*T))
#  fitfunc2 = lambda p,t: p[0]
#  counter = 0
#  for lo in range(13, 15):
#    for up in range(C4.shape[1], C4.shape[1]-2, -1):
#      if (up-lo-2) < 2: # only fits with dof>1 are ok
#        continue
#      t = np.asarray(tt[lo:up])
#      print "\nResult from bootstrapsample fit in interval (2pt):", lo, up-1
#      start_parm = [C4_mean[lo], m4_mean[lo]]
#      res_C4, chi2, pvalue = fit.fitting(fitfunc1, t, C4[:,lo:up], \
#                                         start_parm, correlated = 1)
#      res_C4_mean = np.mean(res_C4, axis=0)
#      label = ['t', 'm_C4(t)', 'data', \
#               'fit, tmin = %d, tmax = %d' % (lo, up-1)]
#      plot.corr_fct_with_fit(tt, m4_mean, m4_error, fitfunc2, [res_C4_mean[1]], \
#                             [7., m4.shape[1]], label, pdfplot, 0)
#      # collecting fitresults
#      if counter is 0:
#        fitresult = res_C4[:,1]
#        fitdata = np.array([chi2, pvalue, lo, up-1])
#        counter += 1
#      else:
#        fitresult = np.vstack((fitresult, res_C4[:,1]))
#        fitdata = np.vstack((fitdata, np.array([chi2, pvalue, lo, up-1])))
#        counter += 1
#  # save fitted masses on disk
#  path = 'bootdata/C4_derivative_massfit_p%d' % p
#  IOcontraction.ensure_dir(path)
#  np.save(path, fitresult)
#  # save fitted masses on disk
#  path = 'bootdata/C4_derivative_massfit_params_p%d' % p
#  IOcontraction.ensure_dir(path)
#  np.save(path, fitdata)
#
#pdfplot.close()





