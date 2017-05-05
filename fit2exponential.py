import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
from scipy.odr import *

import numbers

import warnings
warnings.filterwarnings(action="ignore", module="matplotlib", message="^tight_layout")

def fit2exponential (G, t, p, mu, plot = False):
  # find first zero in G and cut
  makeFinite = np.where(np.real(G) < 1e-100)# & np.abs(np.imag(np.log(G))) > 1e-2)
  if makeFinite[0].size:
    end = makeFinite[0][0]
  else:
    end = G.size - 1

  # removed zeros which will be troublesome when using log later
  t = t[:end]
  G = np.real(G[:end])

  upperBound = end - 1


  Z = 1
  E = 0
  lowerBound = 0
  fittingData = True
  while fittingData:

    # Define a function (quadratic in our case) to fit the data with.
    def exp_fit(params, t):
      Z, E = params
      return Z*np.exp(-(E - mu)*t)

    # Create a model for fitting.
    exp_model = Model(exp_fit)

    # Create a RealData object using our initiated data from above.
    data = RealData(t[lowerBound:], G[lowerBound:])

    # Set up ODR with the model and data.
    odr = ODR(data, exp_model, beta0=[Z, E])

    # Run the regression.
    out = odr.run()

    Z = out.beta[0]
    E = out.beta[1]

    g_fit = Z*np.exp(-(E - mu)*t)

    # find intersection
    signs = np.sign(G[lowerBound:] - g_fit[lowerBound:])
    firstSign = signs[0]

    for i, s in enumerate(signs):
      if s != firstSign:
        if i > 1:
          lowerBound += i;
          break
        else:
          fittingData = False
          break

      if i == signs.size:
        fittingData = False

  if plot:

    f, axarr = plt.subplots(3, sharex=True, figsize=(20, 10))

    axarr[0].set_title(r'$p = {0}, \, \mu = {1}, \, E_0 = {2}, \, Z_0 = {3}$'.format(p, mu, E, Z))

    axarr[0].plot(t, G)
    axarr[0].plot(t[lowerBound], G[lowerBound], 'o', color='green')
    axarr[0].plot(t, Z*np.exp(-(E - mu)*t), '--', color='red')
    axarr[0].plot(t[upperBound], G[upperBound], 'o', color='green')
    axarr[0].set_xlabel(r'$\tau$')
    axarr[0].set_ylabel(r'$G$')

    axarr[1].semilogy(t, G)
    axarr[1].semilogy(t[lowerBound], G[lowerBound], 'o', color='green')
    axarr[1].semilogy(t, Z*np.exp(-(E - mu)*t), '--', color='red')
    axarr[1].semilogy(t[upperBound], G[upperBound], 'o', color='green')
    axarr[1].set_xlabel(r'$\tau$')
    axarr[1].set_ylabel(r'$\log G$')

    axarr[2].plot(t[lowerBound:], G[lowerBound:] - g_fit[lowerBound:])
    axarr[2].plot(t[lowerBound:], np.zeros(t.size - lowerBound), color='red')
    axarr[2].set_xlabel(r'$\tau$')
    axarr[2].set_ylabel(r'$G - G_{\mathrm{fit}}$')

    plt.tight_layout()
    plt.savefig('plots/foo.pdf')
    plt.show()

  return E, Z
