#!/usr/bin/python

import os, math, random, struct
import numpy as np

import matplotlib
matplotlib.use('QT4Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.mlab as mlab

def corr_fct_with_fit(X, Y, dY, fitfunc, args, plotrange, label, pdfplot, logscale = 0):
  # plotting the data
  l = plotrange[0]
  u = plotrange[1]
  p1 = plt.errorbar(X[l:u], Y[l:u], dY[l:u], fmt='x' + 'b', label = label[3])
  # plotting the fit function
  x1 = np.linspace(l, u, 1000)
  y1 = []
  for i in x1:
    y1.append(fitfunc(args,i))
  y1 = np.asarray(y1)
  p2, = plt.plot(x1, y1, 'r', label = label[2])
  # adjusting the plot style
  plt.grid(True)
  plt.xlabel(label[0])
  plt.ylabel(label[1])
  plt.legend([p1, p2], [label[2], label[3]])
  if logscale:
    plt.yscale('log')
  # save pdf
  pdfplot.savefig()
  plt.clf()
