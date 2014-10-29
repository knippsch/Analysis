#!/usr/bin/python

import os, errno, struct

################################################################################
# extract a correlation function and resorts it
def extract_corr_fct(name='', start_cfg=0, delta_cfg=0, nb_cfg=0, T=0, gamma=5,\
                     verbose = 0): 
  re = []
  im = []
  for x in range(start_cfg, start_cfg+delta_cfg*nb_cfg, delta_cfg):
    filename = name + "%04d" % x + '.dat'
    f = open(filename, "rb") # Open a file
    if verbose:
      print "reading from file: ", f.name
    f.seek(2*8*T*gamma)
    for t in range(0, T):
      re.insert(t, struct.unpack('d', f.read(8))) # returns a tuple -> convers.
      im.insert(t, struct.unpack('d', f.read(8))) # returns a tuple -> convers.
    f.close(); # close the file  
  # conversion of the tuple to list and reorganise
  corr = [complex(0.0, 0.0)]*nb_cfg*T
  t = 0
  for x, y in zip(re, im):
    corr[(t%T)*nb_cfg + t/T] = complex(x[0], y[0])
    t += 1
  return corr

################################################################################
# checks if the directory where the file will be written does exist
def ensure_dir(f):
  d = os.path.dirname(f)
  if not os.path.exists(d):
    os.makedirs(d)

################################################################################
# writing the correlation function in new order
def write_corr_fct(name, re_im, corr, verbose = 1):
  ensure_dir(name)
  f = open(name, "wb") # Open a file
  if verbose:
    print "writing to file: ", f.name
  if re_im == 'real':
    for x in corr:
      f.write(struct.pack('d', x.real))
  elif re_im == 'imag':
    for x in corr:
      f.write(struct.pack('d', x.imag))
  else:
    print 'wring re_im -> must be real or imag'
  f.close()

