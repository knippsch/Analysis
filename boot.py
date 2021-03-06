#!/usr/bin/python

import math, random
import numpy as np

# bootstrapping
def bootstrap(X, bootstrapsize):
  np.random.seed(1227)
  boot = np.empty(bootstrapsize, dtype=float)
  for i in range(0, bootstrapsize):
    rnd = np.random.random_integers(0, high=len(X)-1, size=len(X))
    boot_dummy = 0.0 
    for j in range(0, len(X)):
      boot_dummy += X[rnd[j]]  
    boot[i] = boot_dummy/len(X)
  return boot

# symmetrising and bootstrapping
def sym_and_boot(X, T, nb_cfg, bootstrapsize = 1000):
  boot = bootstrap(X[0:nb_cfg], bootstrapsize)
  for t in range(1, T/2):
    data = []
    for a, b in zip(X[t*nb_cfg:(t+1)*nb_cfg], X[(T-t)*nb_cfg:(T-t+1)*nb_cfg]):
      data.append( (a+b)/2.0)
    boot = np.c_[boot, bootstrap(data, bootstrapsize)]
  boot = np.c_[boot, bootstrap(X[(T/2)*nb_cfg:(T/2+1)*nb_cfg], bootstrapsize)]
  return boot

# ratio computation
def compute_ratio(C4, C2):
  print '\ncompute ratio:\n--------------\n' 
  ratio,  sigma, val = [], [], []
  for t in range(1, T/2-1):
    for b in range(0, bootstrapsize):
      a = (C4[t*bootstrapsize + b] - C4[(t+1)*bootstrapsize + b]) / \
          ((C2[t*bootstrapsize + b])**2 - (C2[(t+1)*bootstrapsize + b])**2)
      ratio.append(a)
    print t, mean(ratio[(t-1)*bootstrapsize:(t)*bootstrapsize]), \
                  std_error(ratio[(t-1)*bootstrapsize:(t)*bootstrapsize])
    sigma.append(std_error(ratio[(t-1)*bootstrapsize:(t)*bootstrapsize]))
    val.append(mean(ratio[(t-1)*bootstrapsize:(t)*bootstrapsize]))
  return ratio, sigma, val

# derivative
def compute_derivative(boot):
  print '\ncompute derivative:\n-------------------\n' 
  derv = np.empty([boot.shape[0], boot.shape[1]-1], dtype=float)
  # computing the derivative
  for b in range(0, boot.shape[0]):
    row = boot[b,:]
    for t in range(0, len(row)-1):
      derv[b, t] = row[t+1] - row[t]
  mean, err = mean_error_print(derv)
  return derv, mean, err

# computes the mean and the error, and writes both out
def mean_error_print(boot):
  mean = np.mean(boot, axis=0)
  err  = np.std(boot, axis=0)
  for t, m, e in zip(range(0, len(mean)), mean, err):
    print t, m, e
  return mean, err

# compute the mean correlator with error
def return_mean_corr(boot):
  print'\nmean correlator:\n----------------\n' 
  mean, err = mean_error_print(boot)
  return mean, err

# mass computation
def compute_mass(boot):
  print '\ncompute mass:\n-------------\n' 
  # creating mass array from boot array
  mass = np.empty([boot.shape[0], boot.shape[1]-2], dtype=float)
  # computing the mass via formula
  for b in range(0, boot.shape[0]):
    row = boot[b,:]
    for t in range(1, len(row)-1):
      mass[b, t-1] = (row[t-1] + row[t+1])/(2.0*row[t])
  mass = np.arccosh(mass)
  mean, err = mean_error_print(mass)
  return mass, mean, err


