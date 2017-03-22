import glob
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
from scipy.odr import *

import numbers

import warnings
warnings.filterwarnings(action="ignore", module="scipy", message="^internal gelsd")
warnings.filterwarnings(action="ignore", module="matplotlib", message="^tight_layout")

fileNames = ['a=1 N=30B t=80/*', 'a=1 N=30B t=50/*', 'a=1 N=50B/*']
# fileNames = ['a=1 N=50B/*']
for fileName in fileNames:
  paths = glob.glob('./data/' + fileName + '.txt', recursive=True)

  filters = [
    # {
    #   'key': 'a',
    #   'value': 1
    # }
    # ,{
    #   'key': 'p',
    #   'value': {
    #     'min': 0,
    #     'max': 2
    #   }
    # }
  ]

  Ps = []
  Es = []
  Zs = []
  As = []
  Es_err = []
  Zs_err = []


  for j, path in enumerate(paths):
    with open(path) as fp:

      T = [];
      G = [];

      counter = 0
      passed = True

      for i, line in enumerate(fp):

        # fetch data from first line
        if line.find('--------') == 0:
          data = line[9: -10].split(' ')
          parameters = {}
          for d in data:
            separated = d.split('=', 1)
            if (len(separated) == 1):
              parameters[key] += ' ' + separated[0]
            else:
              key = separated[0]
              parameters[key] = separated[-1]

          for key in parameters:
            if key != 'date':
              parameters[key] = float(parameters[key])

          # filter out
          for filter in filters:
            if isinstance(filter['value'], numbers.Number) and parameters[filter['key']] != filter['value']:
              passed = False
              break
            if isinstance(filter['value'], dict) and not (filter['value']['min'] <= parameters[filter['key']] <= filter['value']['max']):
              passed = False
              break
          if not passed:
            break

          counter = 1
          continue

        line = line.split('\n')[0]
        line = line.split(' ')

        if counter == 1:
          T = [float(part) for part in line]
          counter += 1
          continue

        if counter == 2:
          G = [float(part) for part in line]
          counter = 0

      if passed:



        # find first zero in G and cut
        try:
          makeFinite = G.index(0)
        except ValueError:
          makeFinite = len(G)


        T = np.array(T)
        G = np.array(G)


        t = np.array(T[:makeFinite])
        g = np.array(G[:makeFinite])


        logG = sig.savgol_filter(np.log(g), 51, 3)
        diffLogG = sig.savgol_filter(np.diff(logG), 51, 3)/(2*T[0])

        # lower bound
        lowerBound = 0
        interval = 200
        for i, val in enumerate(diffLogG):
          found = True
          if i + interval <= diffLogG.size - 1:
            for j in range(i + 1, i + 1 + interval):
              if diffLogG[i] < diffLogG[j]:
                found = False

          if found:
            lowerBound = i
            break

        
        upperBound = makeFinite - lowerBound - 250






        # Define a function (quadratic in our case) to fit the data with.
        def exp_fit(params, t):
          Z, E = params
          return Z*np.exp(-(E - parameters['mu'])*t)

        # Create a model for fitting.
        exp_model = Model(exp_fit)

        # Create a RealData object using our initiated data from above.
        data = RealData(T[lowerBound:lowerBound+upperBound], G[lowerBound:lowerBound+upperBound])

        # Set up ODR with the model and data.
        odr = ODR(data, exp_model, beta0=[1, 0])

        # Run the regression.
        out = odr.run()

        Z = out.beta[0]
        E = out.beta[1]

        Z_err = out.sd_beta[0]
        E_err = out.sd_beta[1]


        Ps.append(parameters['p'])
        Es.append(E)
        Zs.append(Z)
        As.append(parameters['a'])
        Es_err.append(E_err)
        Zs_err.append(Z_err)

        print('-------')
        print('a=' + str(parameters['a']))
        print('p=' + str(parameters['p']))
        print('E=' + str(E))
        print('Z=' + str(Z))


        if False:
          f, axarr = plt.subplots(3, sharex=True, figsize=(20, 10))

          axarr[0].set_title(
            r'$p = {0}, \, \alpha = {1}, \, \mu = {2}, \, N = {3} \, E = {4}, \, Z = {5}$'.format(
              parameters['p'],
              parameters['a'],
              parameters['mu'],
              parameters['N'],
              E,
              Z
            )
          )

          axarr[0].plot(T, G)
          axarr[0].plot(T[lowerBound], G[lowerBound], 'o', color='green')
          axarr[0].plot(T, Z*np.exp(-(E - parameters['mu'])*T), '--', color='#ff00ee')
          axarr[0].plot(T[upperBound + lowerBound], G[upperBound + lowerBound], 'o', color='green')


          axarr[1].plot(t, np.log(g))
          axarr[1].plot(t, logG, '--', color='orange')
          axarr[1].plot(T[lowerBound], logG[lowerBound], 'o', color='green')
          axarr[1].plot(T[upperBound + lowerBound], logG[upperBound + lowerBound], 'o', color='green')
          axarr[1].plot(T, np.log(Z) - (E - parameters['mu'])*T, '--', color='#ff00ee')


          axarr[2].plot(t[0: -1], np.diff(logG, n=1)/(2*T[0]))
          axarr[2].plot(t[0: -1], diffLogG, 'orange')
          axarr[2].plot(T[lowerBound], diffLogG[lowerBound], 'o', color='green')
          axarr[2].plot(T[upperBound + lowerBound], diffLogG[upperBound + lowerBound], 'o', color='green')
          axarr[2].plot(T, - (E - parameters['mu']) + T*0, '--', color='#ff00ee')

          plt.tight_layout()
          plt.show()

        # plt.plot(parameters['a'], E0, '*')


  # plt.plot(Ps, Es)

  plt.plot(Ps, Es, '-')
  plt.errorbar(Ps, Es, yerr=Es_err, linestyle='None', marker='x')


x = np.array([0, 20.657, 32.947, 46.544, 65.981, 92.46, 104, 131, 147, 161, 168, 174, 180, 186, 189, 192])/104
y = np.array([182, 180, 176, 168.9, 154, 125.8, 111, 72.3, 48.7, 30.3, 21.31, 13.23, 7.334, 3.234, 2.109, 1.5])/(-181)

plt.plot(x, y, ':')


#   plt.plot(Ps, Zs)

# print(Es_err)

plt.errorbar(Ps, Es, yerr=Es_err, linestyle='None', marker='x')

plt.savefig('plots/foo.pdf')
plt.show()

