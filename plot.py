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
  l = int(plotrange[0])
  u = int(plotrange[1])
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


# this can be used to plot the chisquare distribution of the fits
#  x = np.linspace(scipy.stats.chi2.ppf(1e-6, dof), scipy.stats.chi2.ppf(1.-1e-6, dof), 1000)
#  hist, bins = np.histogram(chisquare, 50, density=True)
#  width = 0.7 * (bins[1] - bins[0])
#  center = (bins[:-1] + bins[1:]) / 2
#  plt.xlabel('x')
#  plt.ylabel('chi^2(x)')
#  plt.grid(True)
#  plt.plot(x, scipy.stats.chi2.pdf(x, dof), 'r-', lw=2, alpha=1, label='chi2 pdf')
#  plt.bar(center, hist, align='center', width=width)
#  plt.show()

