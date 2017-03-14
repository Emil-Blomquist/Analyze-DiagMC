import glob
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig

import numbers

import warnings
warnings.filterwarnings(action="ignore", module="scipy", message="^internal gelsd")
warnings.filterwarnings(action="ignore", module="matplotlib", message="^tight_layout")

fileNames = ['a=1 mu=-0.02/*']
# fileNames = ['./p=0/*']
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
        diff2LogG = sig.savgol_filter(np.diff(diffLogG), 51, 3)/(2*T[0])

        lowerBound = np.argmax(diff2LogG < diff2LogG[0]*0.0000001)
        lowerBound = 250*3


        data = np.diff(diffLogG[lowerBound:])/(2*T[0])
        mean = np.mean(data)
        deviation = np.abs(data - mean)
        deviation = deviation/np.sum(deviation)
        deviation = deviation - deviation[0]
        maxmax = 0.01
        # maxmax = np.mean(deviation)
        upperBound = np.argmax(deviation > maxmax) - 1 #makeFinite - lowerBound - 3
        if upperBound == -1:
          upperBound = makeFinite - lowerBound - 3

        M1 = T[lowerBound:lowerBound+upperBound].reshape((T[lowerBound:lowerBound+upperBound].size, 1))
        M2 = np.ones((T[lowerBound:lowerBound+upperBound].size, 1))
        M = np.concatenate((M1, M2), axis=1)
        Y = np.log(G[lowerBound:lowerBound+upperBound]).reshape((G[lowerBound:lowerBound+upperBound].size, 1))

        k, m = tuple(np.linalg.inv(M.T.dot(M)).dot(M.T.dot(Y)))


        E0 = parameters['mu'] - k[0]
        Z0 = np.exp(m[0])

        Ps.append(parameters['p'])
        Es.append(E0)
        Zs.append(Z0)
        As.append(parameters['a'])

        print('-------')
        print('a=' + str(parameters['a']))
        print('p=' + str(parameters['p']))
        print('E0=' + str(E0))
        print('Z0=' + str(Z0))


        if False:
          f, axarr = plt.subplots(5, sharex=True, figsize=(20, 10))

          axarr[0].set_title(
            r'$p = {0}, \, \alpha = {1}, \, \mu = {2}, \, N = {3} \, E = {4}, \, Z = {5}$'.format(
              parameters['p'],
              parameters['a'],
              parameters['mu'],
              parameters['N'],
              E0,
              Z0
            )
          )

          axarr[0].plot(T, G)
          axarr[0].plot(T[lowerBound], G[lowerBound], 'o', color='green')
          axarr[0].plot(T, Z0*np.exp(-(E0 - parameters['mu'])*T), '--', color='#ff00ee')
          axarr[0].plot(T[upperBound + lowerBound], G[upperBound + lowerBound], 'o', color='green')


          axarr[1].plot(t, np.log(g))
          axarr[1].plot(t, logG, '--', color='orange')
          axarr[1].plot(T[lowerBound], logG[lowerBound], 'o', color='green')
          axarr[1].plot(T[upperBound + lowerBound], logG[upperBound + lowerBound], 'o', color='green')
          axarr[1].plot(T, m + k*T, '--', color='#ff00ee')


          axarr[2].plot(t[0: -1], np.diff(logG, n=1)/(2*T[0]))
          axarr[2].plot(t[0: -1], diffLogG, 'orange')
          axarr[2].plot(T[lowerBound], diffLogG[lowerBound], 'o', color='green')
          axarr[2].plot(T[upperBound + lowerBound], diffLogG[upperBound + lowerBound], 'o', color='green')
          axarr[2].plot(T, k + T*0, '--', color='#ff00ee')


          axarr[3].plot(t[0: -2], np.diff(diffLogG)/(2*T[0]))
          axarr[3].plot(t[0: -2], diff2LogG, 'orange')
          axarr[3].plot(T[lowerBound], diff2LogG[lowerBound], 'o', color='green')
          axarr[3].plot(T[upperBound + lowerBound], diff2LogG[upperBound + lowerBound], 'o', color='green')


          axarr[4].plot(t[lowerBound: -2], deviation)
          axarr[4].plot([0, t[-1]], [maxmax] * 2, '--', color='red')
          axarr[4].plot([0, t[-1]], [0] * 2, '--', color='red')
          axarr[4].plot(T[lowerBound], deviation[0], 'o', color='green')
          axarr[4].plot(T[upperBound + lowerBound], deviation[upperBound], 'o', color='green')

          plt.tight_layout()
          plt.show()

        # plt.plot(parameters['a'], E0, '*')


  # plt.plot(Ps, Es)

x = np.array([0, 20.657, 32.947, 46.544, 65.981, 92.46, 104, 131, 147, 161, 168, 174, 180, 186, 189, 192])/104
y = np.array([182, 180, 176, 168.9, 154, 125.8, 111, 72.3, 48.7, 30.3, 21.31, 13.23, 7.334, 3.234, 2.109, 1.5])/(-181)

plt.plot(x, y)
plt.plot(x, np.array(y) - 0.05)

print(x)
print(np.array(y) - 0.05)

# plt.show()


#   plt.plot(Ps, Zs)
plt.plot(Ps, Es, '-o')
plt.plot(Ps, Es, '-o')

plt.show()

