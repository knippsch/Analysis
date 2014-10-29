#!/usr/bin/python

from scipy.optimize import curve_fit
from scipy.optimize import leastsq
import numpy as np

###############################################################################
# Fit routine
def fitting(fitfunc, X, Y, start_parm,  correlated = 1, writescreen = 1):
  errfunc = lambda p, x, y, error: np.dot(error, (y-fitfunc(p,x)).T)
  # compute inverse, cholesky decomposed covariance matrix
  if correlated == 0:
    if writescreen:
      print '\tPerforming an uncorrelated fit!'
    cov = np.diag(np.diagonal(np.cov(Y.T)))
  else:
    if writescreen:
      print '\tPerforming a correlated fit!'
    cov = np.cov(Y.T)
  cov = np.linalg.inv(np.linalg.cholesky(cov))
  # degrees of freedom
  dof = float(Y.shape[1]-len(start_parm)) 
  # The FIT to the boostrap samples
  res = np.zeros((Y.shape[0], len(start_parm)))
  chisquare = np.zeros(Y.shape[0])
  for b in range(0, Y.shape[0]):
    p,cov1,infodict,mesg,ier = leastsq(errfunc, start_parm, \
                               args=(X, Y[b,:], cov), full_output=1)
    chisquare[b] = float(sum(infodict['fvec']**2.)/dof)
    res[b] = np.array(p)
  res_mean, res_std = np.mean(res, axis=0), np.std(res, axis=0)
  # The fit to the mean value
  y = np.mean(Y, axis=0)
  p,cov1,infodict,mesg,ier = leastsq(errfunc, start_parm, \
                             args=(X, y, cov), full_output=1)
  # writing results to screen
  if writescreen:
    print '\tDegrees of freedom:', dof
    print '\n\tFit results from bootstrap fit:'
    for rm, rs in zip(res_mean, res_std):
      print '\t%.6e' % rm, '+/-', '%.6e' % rs
    print '\tChi^2/dof: %.6e' % np.median(chisquare), '+/- %.6e' % np.std(chisquare)
    print '\n\tFit parameter from mean value fit:'
    for pp in p:
      print '\t%.6e' % pp
    print '\tChi^2/dof: %.6e' % float(sum(infodict['fvec']**2.)/dof)

  return res
