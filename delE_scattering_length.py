#!/usr/bin/python

import math
import numpy as np
import scipy.optimize
import scipy.stats
from itertools import chain

import matplotlib
matplotlib.use('QT4Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.mlab as mlab

import IOcontraction

# some global variables
################################################################################
p = 0 # momenta to analyse
L = 24 # lattice extend
pvalue_cut = 0.01 # minimal p-value
min_fitrange = 4 # minimal fit range

# normalisation of a Matrix
################################################################################
def normalise(COV):
  norm = np.empty([COV.shape[0], COV.shape[1]], dtype=float)
  for i in range(0, COV.shape[0]):
    for j in range(0, COV.shape[1]):
      norm[i, j] = COV[i, j]/(np.sqrt(COV[i,i]*COV[j,j]))
  return norm

# choose correct p-values and fit ranges
# TODO: it is a little bit hacked to ensure the correct dimensions of M_cut
################################################################################
def choose_pvalue_cut(M, params, pcut, min_fitrange):
  counter = 0
  for i in range(0, params.shape[0]):
    if (params[i,1] >= pcut) and (params[i,1] <= 1.-pcut):
      if params[i,3] - params[i,2] > min_fitrange:
        if counter is 0:
          M_cut = M[i,:]
          params_out = params[i]
          counter += 1
        else:
          M_cut = np.vstack((M_cut, M[i,:]))
          params_out = np.vstack((params_out, params[i]))
  if (counter > 0):
    if (M_cut.shape[0] != M.shape[1]):
      return M_cut, params_out
    else:
      return np.vstack((M_cut, np.empty([0, M.shape[1]], dtype=float))), \
             params_out
  else:
    return np.empty([0, M.shape[1]], dtype=float), []

# computation of the energy shift
################################################################################
def compute_energy_shift(C2, C2_params, C4, C4_params):
  count = 0
  if (C4.shape[0] != 0) or (C2.shape[0] != 0):
    for i in range(0, C4.shape[0]):
      for j in range(0, C2.shape[0]):
        V = np.add(C4[i,:], np.multiply(-2., C2[j,:]) )
        C2_h = C2[j,:]
        weight = 1./(abs(C2_params[j,1]-0.5) + abs(C4_params[i,1]-0.5))
        if count == 0:
          delE = V
          C2_out = C2_h
          weight_out = [weight]
          count = 1
        else:
          delE = np.vstack((delE, V))
          C2_out = np.vstack((C2_out, C2_h))
          weight_out.append(weight)
    return delE, weight_out, C2_out
  else:
    return np.empty([0, C2.shape[1]], dtype=float), [], \
           np.empty([0, C2.shape[1]], dtype=float)

# Luescher's function to third order - rewritten with less divisions
################################################################################
def Luscher_func_3rd(a, delE, L, m):
  return -delE + a * (-80.1129*a*a + 35.6535*a*L - 12.5664*L*L)/(L**5 * m)

# compute weights
################################################################################
def compute_weight(params):
  if len(params) != 0:
    weights = []
    for i in range(0, params.shape[0]):
      weights.append(1./(2.*abs(params[i,1]-0.5)))
    return weights
  else:
    return []

# compute the weighted quantile
################################################################################
def weighted_quantile(data, weights, quantile):
  ind_sorted = np.argsort(data)
  sorted_data = data[ind_sorted]
  sorted_weights = weights[ind_sorted]
  # Compute the auxiliary arrays
  Sn = np.cumsum(sorted_weights)
  # TODO: Check that the weights do not sum zero
  Pn = (Sn-0.5*sorted_weights)/np.sum(sorted_weights)
  # Get the value of the weighted median
  interpolated_quant = np.interp(quantile, Pn, sorted_data)
  for i in range(0, len(sorted_data)):
    a = sorted_data[i] - interpolated_quant
    if a > 0:
      b = abs(sorted_data[i-1] - interpolated_quant)
      if b > a:
        quant = sorted_data[i]
      else:
        quant = sorted_data[i-1]
      break
  return interpolated_quant, quant


################################################################################
################################################################################
################################################################################
# (2pt) reading data extracted from direct fit
filename = 'bootdata/C2_massfit_p%d.npy' % p
C2_read = np.load(filename)
filename = 'bootdata/C2_massfit_params_p%d.npy' % p
C2_params = np.load(filename)
# (4pt) reading data extracted from direct fit
filename = 'bootdata/C4_direkt_massfit_p%d.npy' % p
C4_direct_read = np.load(filename)
filename = 'bootdata/C4_direkt_massfit_params_p%d.npy' % p
C4_direct_params = np.load(filename)
# (4pt) reading data extracted from derivative fit
filename = 'bootdata/C4_derivative_massfit_p%d.npy' % p
C4_derv_read = np.load(filename)
filename = 'bootdata/C4_derivative_massfit_params_p%d.npy' % p
C4_derv_params = np.load(filename)
# (4pt) reading data extracted from direkt fit
filename = 'bootdata/C4_ratio_massfit_p%d.npy' % p
C4_ratio_read = np.load(filename)
filename = 'bootdata/C4_ratio_massfit_params_p%d.npy' % p
C4_ratio_params = np.load(filename)
# (4pt) reading data extracted from direkt fit
filename = 'bootdata/C2_ratio_massfit_p%d.npy' % p
C2_ratio_read = np.load(filename)
filename = 'bootdata/C2_ratio_massfit_params_p%d.npy' % p
C2_ratio_params = np.load(filename)

# choosing combinations with wanted p-values and fit ranges
################################################################################
C2_mass, C2_mass_params = \
         choose_pvalue_cut(C2_read, C2_params, pvalue_cut, min_fitrange)
C4_direct_mass, C4_direct_mass_params = \
         choose_pvalue_cut(C4_direct_read, C4_direct_params, pvalue_cut, \
                                                             min_fitrange)
C4_derv_mass, C4_derv_mass_params = \
         choose_pvalue_cut(C4_derv_read, C4_derv_params, pvalue_cut, \
                                                             min_fitrange)
C4_ratio_mass, C4_ratio_mass_params = \
         choose_pvalue_cut(C4_ratio_read, C4_ratio_params, pvalue_cut, \
                                                           min_fitrange)
C2_ratio_mass, C2_ratio_mass_params = \
         choose_pvalue_cut(C2_ratio_read, C2_ratio_params, pvalue_cut, \
                                                           min_fitrange)

# compute delta E and its weights
################################################################################
delE_direct, delE_direct_weight, C2_mass1 = \
    compute_energy_shift(C2_mass, C2_mass_params, \
                         C4_direct_mass, C4_direct_mass_params)
delE_derv, delE_derv_weight, C2_mass2 = \
    compute_energy_shift(C2_mass, C2_mass_params, \
                         C4_derv_mass, C4_derv_mass_params)
delE_ratio, delE_ratio_weight = \
    C4_ratio_mass, compute_weight(C4_ratio_mass_params)
delE = np.vstack((delE_direct, delE_derv, delE_ratio))
delE_weight = np.asarray(list(chain(
                delE_direct_weight, delE_derv_weight, delE_ratio_weight)), \
                                                                 dtype=float)
C2_mass_delE = np.vstack((C2_mass1, C2_mass2, C2_ratio_mass))

# computation of a_pipi
################################################################################
#a_pipi = np.empty([delE.shape[0], delE.shape[1]], dtype=float)
#for i in range(0, delE.shape[0]):
#  for j in range(0, delE.shape[1]):
#    start = -delE[i, j]*C2_mass_delE[i,j]*L*L*L/(4.*math.pi)
#    sol = scipy.optimize.root(Luscher_func_3rd, start, \
#                              args=(delE[i, j], L, C2_mass_delE[i,j]))
#    a_pipi[i, j] = sol.x
path='./bootdata/a_pipi'
#IOcontraction.ensure_dir(path)
#np.save(path, a_pipi)
a_pipi = np.load(path+'.npy')

# computation of quantiles and so on
################################################################################
a_pipi_median_inter, a_pipi_median = weighted_quantile(a_pipi[:,0], \
                                                       delE_weight, 0.5)
a_pipi_16quant_inter, a_pipi_16quant = weighted_quantile(a_pipi[:,0], \
                                                       delE_weight, 0.16)
a_pipi_84quant_inter, a_pipi_84quant = weighted_quantile(a_pipi[:,0], \
                                                       delE_weight, 0.84)

index = np.where(a_pipi == a_pipi_median)
a_pipi_std = np.std(a_pipi[index[0],:])
a_pipi_sys_lo = a_pipi_median - a_pipi_16quant_inter
a_pipi_sys_hi = a_pipi_84quant_inter - a_pipi_median
print 'median +/- stat + sys - sys'
print a_pipi_median, a_pipi_std, a_pipi_sys_hi, a_pipi_sys_lo
                                 
# creating a histogram with all data
################################################################################
hist, bins = np.histogram(a_pipi[:,0], 15, weights=delE_weight, density=True)
print hist
print bins
width = 0.7 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.xlabel('x')
plt.ylabel('distribution of a_pipi')
plt.grid(True)
x = np.linspace(center[0], center[-1], 1000)
plt.plot(x, scipy.stats.norm.pdf(x, loc=a_pipi_median, scale=a_pipi_std),\
         'r-', lw=3, alpha=1, label='median + std. error')
plt.plot(x, scipy.stats.norm.pdf(x, loc=a_pipi_median, \
         scale=0.5*(a_pipi_sys_lo+a_pipi_sys_hi)),\
         'g-', lw=3, alpha=1, label='median +  sys. error')
plt.legend()
plt.bar(center, hist, align='center', width=width, alpha=0.7)
plt.show()

# creating a histogramms for the three methods
################################################################################


