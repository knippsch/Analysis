#!/usr/bin/python

import os, errno, math, random, struct
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
import numpy as np

import IO
import boot
import fit

bootstrapsize = 1000
T = 64 # temporal extend of the lattice
start_cfg = 404
delta_cfg = 4
end_cnfg = 1200
nb_cfg = (end_cnfg-start_cfg)/delta_cfg + 1


################################################################################
# everthing else
name = '../4ptFunction_A40.32/data/C2_pi+-_conf'
corr = IO.extract_corr_fct(name, start_cfg, delta_cfg, nb_cfg, T, 0)
IO.write_corr_fct('./raw_data/pion_mom_0.dat', 'real', corr)
pi = []
for c in corr:
  pi.append(c.real)
print '\ncomputing the effective mass of the pion:\n'
mass_boot = boot.sym_and_boot(pi, T, nb_cfg, bootstrapsize)
boot.return_mean_corr(mass_boot)
massC2, massC2_mean, massC2_var = boot.compute_mass(mass_boot)
print '\n'

#name = './data/C4_1_conf'
#corr1 = extract_corr_fct(name, start_cfg, delta_cfg, nb_cfg, T, 0)
#write_corr_fct('./raw_data/C4_1_mom_0.dat', 'real', corr1)
#C4_1 = []
#for c in corr1:
#  C4_1.append(c.real)
#
#name = './data/C4_2_conf'
#corr2 = extract_corr_fct(name, start_cfg, delta_cfg, nb_cfg, T, 0)
#write_corr_fct('./raw_data/C4_2_mom_0.dat', 'real', corr2)
#C4_2 = []
#for c in corr2:
#  C4_2.append(c.real)
#
#name = './data/C4_3_conf'
#corr3 = extract_corr_fct(name, start_cfg, delta_cfg, nb_cfg, T, 0)
#write_corr_fct('./raw_data/C4_3_mom_0.dat', 'real', corr3)
#C4_3 = []
#for c in corr3:
#  C4_3.append(c.real)
#
#C4 = []
#for a, b, c in zip(C4_1, C4_2, C4_3):
#  C4.append(a+b-2.0*c)
#boot = sym_and_boot(C4)
#C4_mean, C4_var = return_mean_corr(boot)
#ratio, val, sigma = compute_derivative(boot)
#massC4, massC4_mean, massC4_var = compute_mass(ratio)


################################################################################
################ computing the fits ############################################

# array for timeslices
tt = []
for t in range(0, T):
	tt.append(float(t))

# fitting the mass #############################################################
print '\n************************************'
print ' Fitting the mass from the plateau.'
print '************************************'
fitfunc1 = lambda p,x: p[0]

for lo in range(9, 15):
  up = massC2.shape[1]
  t = np.asarray(tt[lo:up])
  print "\nResult from bootstrapsample fit in interval (2pt):", lo, up
  fit.fitting(fitfunc1, t, massC2[:,lo:up], [0.19])

# fitting the correlation function directly ####################################
print '\n**********************************************************'
print ' Fitting the mass directly from the correlation function.'
print '**********************************************************'
fitfunc1 = lambda p,x: p[0]*np.cosh((T/2.0-t)*p[1])

for lo in range(9, 15):
  up = mass_boot.shape[1]
  t = np.asarray(tt[lo:up])
  print "\nResult from bootstrapsample fit in interval (2pt):", lo, up
  res = fit.fitting(fitfunc1, t, mass_boot[:,lo:up], [10., 0.19])


















