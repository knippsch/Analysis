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


def normalise(COV):
  norm = np.empty([COV.shape[0], COV.shape[1]], dtype=float)
  for i in range(0, COV.shape[0]):
    for j in range(0, COV.shape[1]):
      norm[i, j] = COV[i, j]/(np.sqrt(COV[i,i]*COV[j,j]))
  return norm

def choose_pvalue_cut(M, params, cut, min_fitrange):
  counter = 0
  for i in range(0, params.shape[0]):
    if (params[i,1] >= cut) and (params[i,1] <= 1.-cut):
      if params[i,3] - params[i,2] > min_fitrange:
        if counter is 0:
          M_cut = M[i,:]
          counter += 1
        else:
          M_cut = np.vstack((M_cut, M[i,:]))
  if (counter > 0):
    if (M_cut.shape[0] != M.shape[1]):
      return M_cut
    else:
      return np.vstack((M_cut, np.empty([0, M.shape[1]], dtype=float)))
  else:
    return np.empty([0, M.shape[1]], dtype=float)

def compute_energy_shift(C2, C4):
  count = 0
  if (C4.shape[0] is not 0) or (C2.shape[0] is not 0):
    for i in range(0, C4.shape[0]):
      for j in range(0, C2.shape[0]):
        V = np.add(C4[i,:], np.multiply(-2., C2[j,:]) )
        if count is 0:
          delE = V
          count = 1
        else:
          delE = np.vstack((delE, V))
    return delE
  else:
    return np.empty([0, C2.shape[1]], dtype=float)

p = 0

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

# choosing combinations with wanted p-values and fit ranges
pvalue_cut = 0.1
min_fitrange = 8
C2_mass = choose_pvalue_cut(C2_read, C2_params, pvalue_cut, min_fitrange)
C4_direct_mass = choose_pvalue_cut(C4_direct_read, C4_direct_params,\
                                   pvalue_cut, min_fitrange)
C4_derv_mass = choose_pvalue_cut(C4_derv_read, C4_derv_params, pvalue_cut,\
                                 min_fitrange)
C4_ratio_mass = choose_pvalue_cut(C4_ratio_read, C4_ratio_params, pvalue_cut,\
                                  min_fitrange)

# compute delta E
if(C4_ratio_mass.shape[0] is 0):
  delE = np.empty([0, C2_mass.shape[1]], dtype=float)
else:
  delE = C4_ratio_mass
delE = C4_ratio_mass
delE = np.vstack((delE, compute_energy_shift(C2_mass, C4_derv_mass)))
delE = np.vstack((delE, compute_energy_shift(C2_mass, C4_direct_mass)))

def Luscher_func_3rd(x, delE, L, m):
  return - delE - 4.*math.pi*x/(L*L*L*m) * \
         (1. - 2.837297*x/L + 6.375183*x*x/(L*L))
         

a_pipi = np.empty([delE.shape[0], delE.shape[1]], dtype=float)
L=24
THEMASS = C2_mass[0,0]
for i in range(0, delE.shape[0]):
  for j in range(0, delE.shape[1]):
    start = -delE[i, j]*THEMASS*L*L*L/(4.*math.pi)
    sol = scipy.optimize.root(Luscher_func_3rd, start, \
                              args=(delE[i, j], L, THEMASS))
    a_pipi[i, j] = sol.x*THEMASS
#    print a_pipi[i, j]

delE_mean = np.mean(delE, axis=1)
delE_err = np.std(delE, axis=1)
delE_mean_mean = np.mean(delE_mean)
delE_mean_err = np.std(delE_mean)
delE_err_mean = np.mean(delE_err)


print delE_mean_mean, delE_err_mean, delE_mean_err
#for a, b in zip(delE_mean, delE_err):
#  print a, ' +/- ', b

COV = np.cov(delE)
# devision by the number of elements comes from the mean value
delE_sys_err = np.sqrt(np.sum(COV))/len(delE)

print delE_sys_err

hist, bins = np.histogram(delE_mean, 15, density=True)
x = np.linspace(0.98*bins[1], 1.02*bins[-1], 1000)
width = 0.7 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.xlabel('x')
plt.ylabel('distribution of delE')
plt.grid(True)
plt.plot(x, scipy.stats.norm.pdf(x, loc=delE_mean_mean, scale=delE_err_mean),\
         'r-', lw=3, alpha=1, label='std. error of mean')
plt.plot(x, scipy.stats.norm.pdf(x, loc=delE_mean_mean, scale=delE_mean_err),\
         'c-', lw=3, alpha=1, label='mean std. error')
plt.plot(x, scipy.stats.norm.pdf(x, loc=delE_mean_mean, scale=delE_sys_err),\
         'y-', lw=3, alpha=1, label='sys. error from fit method')
plt.legend()
plt.bar(center, hist, align='center', width=width, alpha=0.7)
plt.show()

#print C4_ratio_mass.shape, C4_ratio_mass.shape
#
#COV = np.cov(delE)
#print COV
#
#COV_norm = normalise(COV)
#
#print COV_norm
#
#fig, ax = plt.subplots()
#heatmap = ax.pcolor(COV_norm, cmap=plt.cm.Blues)
#
## put the major ticks at the middle of each cell
##ax.set_xticks(np.arange(COV.shape[0])+0.5, minor=False)
##ax.set_yticks(np.arange(COV.shape[1])+0.5, minor=False)
#
## want a more natural, table-like display
#ax.invert_yaxis()
#ax.xaxis.tick_top()
#
##ax.set_xticklabels(row_labels, minor=False)
##ax.set_yticklabels(column_labels, minor=False)
#plt.show()


#
#L=32.
#m = res_C2[:,1]
#path = './bootdata/mass'
#IOcontraction.ensure_dir(path)
#np.save(path, m)
#
#E = res_C4[:,1]
#path = './bootdata/E'
#IOcontraction.ensure_dir(path)
#np.save(path, E)
#
#q2 = (E*E*0.25-m*m)*(L/(2.*math.pi))**2.
#path = './bootdata/q2'
#IOcontraction.ensure_dir(path)
#np.save(path, q2)
#
#Z = []
#for q, mm in zip(q2, m):
#  Z.append(zeta.Z(q))
#Z = np.asarray(Z)
#path = './bootdata/Z'
#IOcontraction.ensure_dir(path)
#np.save(path, Z)


