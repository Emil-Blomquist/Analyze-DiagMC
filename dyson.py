import numpy as np
from numpy.fft import fft, ifft, fftshift, fftfreq
import matplotlib.pyplot as plt


def Dyson (SE_t, t, a, p, mu):
  EofP = 0.5*p**2

  # make time longer to increase resolution which is needed since we have an extremely sharply peaked distribution in w-space
  dN = 10000
  dt = t[1] - t[0]
  for i in range(0, dN):
    t = np.append(t, t[-1] + dt)
    SE_t = np.append(SE_t, 0)

  N = t.size
  L = t[-1] - t[0]
  dt = t[1] - t[0]

  # bare propagator
  G0_t = np.exp(-(EofP - mu)*t)

  w = fftshift(fftfreq(t.shape[-1], d=dt)) * 2*np.pi

  param = 1 # need to be larger than one!
  SE_t_sing = a*np.exp(-param*t)*(np.pi*t)**-0.5
  SE_w_sing = a*(param + 1j*w)**-0.5

  SE_t_reg = SE_t - SE_t_sing

  # w space
  G0_w = 1/(1j*w + (EofP - mu))
  SE_w_reg = fftshift(fft(SE_t_reg)) * L/N

  SE_w = SE_w_reg + SE_w_sing

  # self energy difference (to avoid numerical problems)
  dG_w = 1/(1/G0_w - SE_w) - G0_w

  # transform back
  dG_t = ifft(fftshift(dG_w)) * N/L

  # go back to original length of vectors
  G0 = G0_t[0: G0_t.size - dN]
  dG = dG_t[0: dG_t.size - dN]

  return np.array(G0 + dG), np.array(t[0: t.size - dN])
