#!/usr/bin/python

import os, errno, math, random, struct
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
import numpy as np

import IOcontraction
import boot
import fit

bootstrapsize = 2500
T = 48 # temporal extend of the lattice
start_cfg = 501
delta_cfg = 8
end_cnfg = 2997
nb_cfg = (end_cnfg-start_cfg)/delta_cfg + 1

################################################################################
# reading two-pint correlation function and bootstrapping it
name = '/data/LapHs/I_2_Scattering_Analysis/A100/data/C2_pi+-_conf'
corr = IOcontraction.extract_corr_fct(name, start_cfg, delta_cfg, nb_cfg, T, 0)
pi = []
for c in corr:
  pi.append(c.real)
C2 = boot.sym_and_boot(pi, T, nb_cfg, bootstrapsize, path='bootdata/C2_p0')

################################################################################
# reading four-pint correlation function and bootstrapping it
name = '/data/LapHs/I_2_Scattering_Analysis/A100/data/C4_1_conf'
corr1 = IOcontraction.extract_corr_fct(name, start_cfg, delta_cfg, nb_cfg, T, 0)
C4_1 = []
for c in corr1:
  C4_1.append(c.real)

name = '/data/LapHs/I_2_Scattering_Analysis/A100/data/C4_2_conf'
corr2 = IOcontraction.extract_corr_fct(name, start_cfg, delta_cfg, nb_cfg, T, 0)
C4_2 = []
for c in corr2:
  C4_2.append(c.real)

name = '/data/LapHs/I_2_Scattering_Analysis/A100/data/C4_3_conf'
corr3 = IOcontraction.extract_corr_fct(name, start_cfg, delta_cfg, nb_cfg, T, 0)
C4_3 = []
for c in corr3:
  C4_3.append(c.real)

pi = []
for a, b, c in zip(C4_1, C4_2, C4_3):
  pi.append(a+b-2.0*c)
C4 = boot.sym_and_boot(pi, T, nb_cfg, bootstrapsize, path='bootdata/C4_p0')




















