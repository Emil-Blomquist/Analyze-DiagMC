import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
from scipy.odr import *

import numbers

import warnings
warnings.filterwarnings(action="ignore", module="matplotlib", message="^tight_layout")

# folder = 'p=0 a=0->6/'
folder = 'a=1 p=0->2.5/a=1 N=200B t=80 p=0->1.8'

# used to calculate mean and error bars
Eall = {'p': {}, 'a': {}}
Zall = {'p': {}, 'a': {}}

for root, dirs, files in os.walk(os.getcwd() + '/data/' + folder):
  # filter files to only accept data files (.txt)
  files = [fi for fi in files if fi.endswith(".txt")]

  filters = [
    # {
    #   'key': 'p',
    #   'value': 0
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


  for j, path in enumerate(files):
    with open(root + '/' + path) as fp:

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

        upperBound = makeFinite - 1

        Z = 1
        E = 0
        lowerBound = 0
        fittingData = True
        while fittingData:

          # Define a function (quadratic in our case) to fit the data with.
          def exp_fit(params, t):
            Z, E = params
            return Z*np.exp(-(E - parameters['mu'])*t)

          # Create a model for fitting.
          exp_model = Model(exp_fit)

          # Create a RealData object using our initiated data from above.
          data = RealData(t[lowerBound:], g[lowerBound:])

          # Set up ODR with the model and data.
          odr = ODR(data, exp_model, beta0=[Z, E])

          # Run the regression.
          out = odr.run()

          Z = out.beta[0]
          E = out.beta[1]

          g_fit = Z*np.exp(-(E - parameters['mu'])*t)

          # find intersection
          signs = np.sign(g[lowerBound:] - g_fit[lowerBound:])
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

        Z_err = out.sd_beta[0]
        E_err = out.sd_beta[1]

        Ps.append(parameters['p'])
        Es.append(E)
        Zs.append(Z)
        As.append(parameters['a'])
        Es_err.append(E_err)
        Zs_err.append(Z_err)

        # to calculate mean values and error bars
        if parameters['p'] in Eall['p'].keys():
          Eall['p'][parameters['p']] = np.append(Eall['p'][parameters['p']], E)
          Zall['p'][parameters['p']] = np.append(Zall['p'][parameters['p']], Z)
        else:
          Eall['p'][parameters['p']] = np.array([E])
          Zall['p'][parameters['p']] = np.array([Z])

        if parameters['a'] in Eall['a'].keys():
          Eall['a'][parameters['a']] = np.append(Eall['a'][parameters['a']], E)
          Zall['a'][parameters['a']] = np.append(Zall['a'][parameters['a']], Z)
        else:
          Eall['a'][parameters['a']] = np.array([E])
          Zall['a'][parameters['a']] = np.array([Z])

        print('-------')
        print('a=' + str(parameters['a']))
        print('p=' + str(parameters['p']))
        print('E=' + str(E))
        print('Z=' + str(Z))


        if False:
          f, axarr = plt.subplots(3, sharex=True, figsize=(20, 10))

          axarr[0].set_title(
            r'$p = {0}, \, \alpha = {1}, \, \mu = {2}, \, N = {3}, \, E_0 = {4}, \, Z_0 = {5}$'.format(
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
          axarr[0].plot(T, Z*np.exp(-(E - parameters['mu'])*T), '--', color='red')
          axarr[0].plot(T[upperBound], G[upperBound], 'o', color='green')
          axarr[0].set_xlabel(r'$\tau$')
          axarr[0].set_ylabel(r'$G$')

          axarr[1].semilogy(T, G)
          axarr[1].semilogy(T[lowerBound], G[lowerBound], 'o', color='green')
          axarr[1].semilogy(T, Z*np.exp(-(E - parameters['mu'])*T), '--', color='red')
          axarr[1].semilogy(T[upperBound], G[upperBound], 'o', color='green')
          axarr[1].set_xlabel(r'$\tau$')
          axarr[1].set_ylabel(r'$\log G$')

          axarr[2].plot(t[lowerBound:], g[lowerBound:] - g_fit[lowerBound:])
          axarr[2].plot(t[lowerBound:], np.zeros(t.size - lowerBound), color='red')
          axarr[2].set_xlabel(r'$\tau$')
          axarr[2].set_ylabel(r'$G - G_{\mathrm{fit}}$')

          plt.tight_layout()
          plt.savefig('plots/foo.pdf')
          plt.show()

  plt.errorbar(Ps, Es, yerr=Es_err, marker='x', label='My result', ls="none")
  # plt.errorbar(As, Es, yerr=Es_err, marker='x', label='My result')
  # plt.plot(Ps, Zs, '.')
  # plt.plot(As, Es, '.')
  # plt.plot(As, Zs, '.')


p = np.array([0, 20.657, 32.947, 46.544, 65.981, 92.46, 104, 131, 147, 161, 168, 174, 180, 186, 189, 192])/104
Ep = np.array([182.7, 180, 176, 168.9, 154, 125.8, 111, 72.3, 48.7, 30.3, 21.31, 13.23, 7.334, 3.234, 2.109, 1.5])/(-181)

a = np.array([0, 2.932, 8.945, 14.915, 29.804, 59.613, 89.362, 119.225, 134.287, 148.998, 178.837])/178.837*6
Ea = np.array([0, 2.697, 8.158, 13.656, 27.292, 55.575, 85.552, 116.807, 132.381, 149.587, 185.606])/188.337*(-7)


# plt.plot(p, Ep, '--', color='black', lw=2, label='Your result')
# plt.plot(a, Ea, ':', marker='.')

meff = 1/(1 - 1/6)
p = np.linspace(0.001, 1.4, 1000)
E0second = -1 - 1.26/100 + p**2/(2*meff)
E0first = 0.5*p**2 - 2**0.5/p * np.arcsin(p*(2**(-0.5)))
# plt.plot(p, E0second, ':', color='magenta', label=r'$E_0^{(2)}(p)$')
# plt.plot(p, E0first, ':', color='brown', label=r'$E_0^{(1)}(p)$')


# plt.legend(loc=4)
# plt.xlim(1.70, 1.86)
# plt.ylim(-0.07, 0)
# plt.title(r'$\alpha = {0}$'.format(parameters['a']))
# plt.xlabel(r'$p$')
# plt.ylabel(r'$E_0(p)$')

# calculate and plot mean and with bars

A = []
P = []
meanEofA = []
stdEofA = []
meanZofA = []
stdZofA = []
meanEofP = []
stdEofP = []
meanZofP = []
stdZofP = []


for a in sorted(Eall['a']):
  A += [a]
  meanEofA += [np.mean(Eall['a'][a])]
  stdEofA += [np.std(Eall['a'][a])]
  meanZofA += [np.mean(Zall['a'][a])]
  stdZofA += [np.std(Zall['a'][a])]


for p in sorted(Eall['p']):
  P += [p]
  meanEofP += [np.mean(Eall['p'][p])]
  stdEofP += [np.std(Eall['p'][p])]
  meanZofP += [np.mean(Zall['p'][p])]
  stdZofP += [np.std(Zall['p'][p])]


# plt.errorbar(A, meanEofA, yerr=stdEofA, marker='.', linestyle=':', color='magenta')
# plt.errorbar(A, meanZofA, yerr=stdZofA, marker='.', linestyle=':', color='magenta')
plt.errorbar(P, meanEofP, yerr=stdEofP, marker='.', linestyle=':', color='magenta')
# plt.errorbar(P, meanZofP, yerr=stdZofP, marker='.', linestyle=':', color='magenta')


print(P)
print(meanEofP)


plt.savefig('plots/foo.pdf')
plt.show()