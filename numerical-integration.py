import numpy as np
import matplotlib.pyplot as plt


N = 100000
tmax = 10
dt = tmax/N
t = np.linspace(0.5*dt, tmax - 0.5*dt, N)


def integrand (t, e, w):
  return np.exp(-(1j*w + e)*t)*t**-0.5



e = 2
w = 100


for w in np.linspace(0, 10, 100):

  I_approx = dt*np.sum(integrand(t, e, w))
  # I_exact = np.pi**0.5*(e + 1j*w)**-0.5

  # plt.plot(w, np.real(I_exact), '.r')
  plt.plot(w, I_approx, '.b')
  # plt.plot(w, np.imag(I_exact), '.r')
  # plt.plot(w, np.imag(I_approx), '.b')



  # print(I_approx, I_exact)

plt.show()


