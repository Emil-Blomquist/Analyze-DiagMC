import numpy as np
from numpy.fft import fft, ifft, fftshift, fftfreq
import matplotlib.pyplot as plt

def heaviside(x):
  return 0.5 * (np.sign(x) + 1)

N = 2000

t = np.linspace(-50, 50, N)
w = fftshift(fftfreq(t.shape[-1], d=t[1] - t[0])) * 2*np.pi

a = 2

# bare propagator
f_t = np.exp(-0.5*t**2*a**-2)
f_w_transformed = fftshift(fft(fftshift(f_t))) * (t[-1] - t[1])/N
f_w = a*(2*np.pi)**0.5*np.exp(-0.5*(a*w)**2)
f_t_transformed = fftshift(ifft(fftshift(f_w))) * N/(t[-1] - t[1])



plt.plot(t, f_t)
plt.plot(t, np.real(f_t_transformed))

# plt.plot(w, np.real(f_w))
# plt.plot(w, np.imag(f_w))
# plt.plot(w, np.real(f_w_transformed))
# plt.plot(w, np.imag(f_w_transformed))




plt.show()