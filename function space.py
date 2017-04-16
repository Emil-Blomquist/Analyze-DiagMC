# from __future__ import division
import numpy as np
import matplotlib.pyplot as plt




L = 1.0
N = 100
dx = L/N



x0 = 0.3


index = int(x0/dx)
x0_discrete = (index + 0.5)*dx



x = np.linspace(0.5*dx, L - 0.5*dx, N)
# k1 = np.pi*1/L
# k2 = np.pi*10/L
# prod = np.sin(k1*x)*np.sin(k2*x)
# print(np.sum(prod))
# plt.plot(x, prod)



S = x*0


# plot()
for n in range(1, N + 1):

  kn = np.pi*n/L
  cn = np.sin(kn*x0_discrete)*1

  plt.plot(x, cn*np.sin(kn*x))

  S += cn*np.sin(kn*x)

plt.plot(x, S)


plt.show()

