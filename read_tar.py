#!/usr/bin/python

#import os, math, tarfile
#import numpy as np

################################################################################

import boot, IOcontraction

T = 48
cnfg_min = 714 
cnfg_max = 2330
cnfg_del = 4

################################################################################
# reading 2pt functions
#print '************************'
#print ' reading 2pt functions:'
#print '************************'
#
#for p in range(0,5):
#  print 'reading data for momentum:', p
#  filename = './pi_corr_p%d.conf' % p
#  corr, nb_cnfg = IOcontraction.extract_corr_tar('../data_Liuming/pi.tar', \
#                                   filename, cnfg_min, cnfg_max, cnfg_del, T, 1)
#  pi = []
#  for c in corr:
#    pi.append(c.real)
#  filename = 'bootdata/C2_p%d' % p
#  C2 = boot.sym_and_boot(pi, T, nb_cnfg, 1000, path=filename)
#  C2_mean, C2_error = boot.mean_error_print(C2, 1)
#  print '\n'
#
#
################################################################################
# reading 4pt functions
print '************************'
print ' reading 4pt functions:'
print '************************'

for p1 in range(0,5):
#  for p2 in range(0,5):
  print 'reading data for momentum:', p1
  # diagram 1
  filename = './pipi_pipi_A1_corr1_p%d%d.conf' % (p1, p1)
  corr1, nb_cnfg = IOcontraction.extract_corr_tar('../data_Liuming/pipi.tar', \
                                   filename, cnfg_min, cnfg_max, cnfg_del, T, 1)
  C4_1 = []
  for c in corr1:
    C4_1.append(c.real)
  # diagram 2
  filename = './pipi_pipi_A1_corr2_p%d%d.conf' % (p1, p1)
  corr2, nb_cnfg = IOcontraction.extract_corr_tar('../data_Liuming/pipi.tar', \
                                   filename, cnfg_min, cnfg_max, cnfg_del, T, 1)
  C4_2 = []
  for c in corr2:
    C4_2.append(c.real)
  # diagram 3
  filename = './pipi_pipi_A1_corr3_p%d%d.conf' % (p1, p1)
  corr3, nb_cnfg = IOcontraction.extract_corr_tar('../data_Liuming/pipi.tar', \
                                   filename, cnfg_min, cnfg_max, cnfg_del, T, 1)
  C4_3 = []
  for c in corr3:
    C4_3.append(c.real)

  # writing C4 to disk
  pi = []
  for a, b, c in zip(C4_1, C4_2, C4_3):
    pi.append(a+b-2.0*c)
  filename = 'bootdata/C4_p%d%d' % (p1, p1)
  C4 = boot.sym_and_boot(pi, T, nb_cnfg, 1000, path=filename)
  C4_mean, C4_error = boot.mean_error_print(C4, 1)
  print '\n'








