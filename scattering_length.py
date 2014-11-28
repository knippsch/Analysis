
      # choosing only fits with certain p values
      if (pvalue > pvalue_minmax) or (pvalue < 1.-pvalue_minmax):
        if counter is 0:
          fitresult = res_C2[:,1]
          fitdata = np.array([chi2, pvalue, lo, up-1])
          counter += 1
        else:
          fitresult = np.vstack((fitresult, res_C2[:,1]))
          fitdata = np.vstack((fitdata, np.array([chi2, pvalue, lo, up-1])))
          counter += 1

  print 'number of accepted fits: ', counter
  # the histogram of the mass fits
  meanvalues = np.mean(fitresult, axis=1)
  hist, bins = np.histogram(meanvalues, bins=counter/5+1, density=False)
  width = 0.7 * (bins[1] - bins[0])
  center = (bins[:-1] + bins[1:]) / 2
  plt.xlabel('m_eff')
  plt.ylabel('Probability')
  plt.grid(True)
  plt.bar(center, hist, align='center', width=width)
  pdfplot.savefig()
  plt.clf()
