from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import eval_sh_legendre



a = 1
omega = 1
mu = -0.217
p = 1.548
EofP = 0.5*p**2

def e_r (n, x, l):
  y = x*l

  if n == 0:
    v = 1
  elif n == 1:
    v = (-1 + 2*y)
  elif n == 2:
    v = (1 - 4*y +  2*y**2)
  elif n == 3:
    v = (-3 + 18*y - 18*y**2 + 4*y**3)/3
  elif n == 4:
    v = (3 - 24*y + 36*y**2 - 16*y**3 + 2*y**4)/3
  elif n == 5:
    v = (-15 + 150*y - 300*y**2 + 200*y**3 - 50*y**4 + 4*y**5)/15
  elif n == 6:
    v = (45 - 540*y + 1350*y**2 - 1200*y**3 + 450*y**4 - 72*y**5 + 4*y**6)/45
  elif n == 7:
    v = (-315 + 4410*y - 13230*y**2 + 14700*y**3 - 7350*y**4 + 1764*y**5 - 196*y**6 + 8*y**7)/315
  elif n == 8:
    v = (315 - 5040*y + 17640*y**2 - 23520*y**3 + 14700*y**4 - 4704*y**5 + 784*y**6 - 64*y**7 + 2*y**8)/315
  elif n == 9:
    v = (-2835 + 51030*y - 204120*y**2 + 317520*y**3 - 238140*y**4 + 95256*y**5 - 21168*y**6 + 2592*y**7 - 162*y**8 + 4*y**9)/2835
  elif n == 10:
    v = (14175 - 283500*y + 1275750*y**2 - 2268000*y**3 + 1984500*y**4 - 952560*y**5 + 264600*y**6 - 43200*y**7 + 4050*y**8 - 200*y**9 + 4*y**10)/14175
  elif n == 11:
    v = (-155925 + 3430350*y - 17151750*y**2 + 34303500*y**3 - 34303500*y**4 + 19209960*y**5 - 6403320*y**6 + 1306800*y**7 - 163350*y**8 + 12100*y**9 - 484*y**10 + 8*y**11)/155925
  elif n == 12:
    v = (467775 - 11226600*y + 61746300*y**2 - 137214000*y**3 + 154365750*y**4 - 98794080*y**5 + 38419920*y**6 - 9408960*y**7 + 1470150*y**8 - 145200*y**9 + 8712*y**10 - 288*y**11 + 4*y**12)/467775

  return v * (2*l)**0.5*np.exp(-y)

def e (x, n, l):
  return e_r(x, n, l)

def Cs_analytical (l, L, a):
  # l: is the decay constant of the basis
  # L: is the the decay constant in the singularity
  return np.array([
    2,
    -2*L * (l + L)**-1,
    (l**2 + 2*L**2) * (l + L)**-2,
    (-3*l**2*L - 2*L**3) * (l + L)**-3,
    (3*l**4 + 24*l**2*L**2 + 8*L**4)/4 * (l + L)**-4,
    -L*(15*l**4 + 40*l**2*L**2 + 8*L**4)/4 * (l + L)**-5,
    (5*l**6 + 90*l**4*L**2 + 120*l**2*L**4 + 16*L**6)/8 * (l + L)**-6,
    -L*(35*l**6 + 210*l**4*L**2 + 168*l**2*L**4 + 16*L**6)/8 * (l + L)**-7,
    (35*l**8 + 1120*l**6*L**2 + 3360*l**4*L**4 + 1792*l**2*L**6 + 128*L**8)/64 * (l + L)**-8,
    -L*(315*l**8 + 3360*l**6*L**2 + 6048*l**4*L**4 + 2304*l**2*L**6 + 128*L**8)/64 * (l + L)**-9,
    (63*l**10 + 3150*l**8*L**2 + 16800*l**6*L**4 + 20160*l**4*L**6 + 5760*l**2*L**8 + 256*L**10)/128 * (l + L)**-10,
    -L*(693*l**10 + 11550*l**8*L**2 + 36960*l**6*L**4 + 31680*l**4*L**6 + 7040*l**2*L**8 + 256*L**10)/128 * (l + L)**-11,
    (231*l**12 + 16632*l**10*L**2 + 138600*l**8*L**4 + 295680*l**6*L**6 + 190080*l**4*L**8 + 33792*l**2*L**10 + 1024*L**12)/512 * (l + L)**-12
  ]) * (0.5*np.pi*l/(l + L))**0.5 * a*np.pi**-0.5

Cs = np.zeros(13)


if False:
  # -------- p=0.0000000 a=1.0000000 mu=-1.0560000 t=10.0000000 N=1000000 date=2017-04-20 14:39:19 --------
  # n = 1
  p = 0
  mu = -1.056
  a = 1
  w = 1
  SE_t = a*np.exp(-(w - mu)*t)*(np.pi*t)**-0.5 - a*np.exp(-t)*(np.pi*t)**-0.5
  L = w - mu

if True:
  # -------- p=0.5000000 a=1.0000000 mu=-1.0560000 t=10.0000000 N=1000000 date=2017-04-21 11:00:26 --------
  # n = 1
  p = 0.5
  mu = -1.056
  a = 1
  w = 1
  l = w - mu

  t = np.array([0.0100000, 0.0300000, 0.0500000, 0.0700000, 0.0900000, 0.1100000, 0.1300000, 0.1500000, 0.1700000, 0.1900000, 0.2100000, 0.2300000, 0.2500000, 0.2700000, 0.2900000, 0.3100000, 0.3300000, 0.3500000, 0.3700000, 0.3900000, 0.4100000, 0.4300000, 0.4500000, 0.4700000, 0.4900000, 0.5100000, 0.5300000, 0.5500000, 0.5700000, 0.5900000, 0.6100000, 0.6300000, 0.6500000, 0.6700000, 0.6900000, 0.7100000, 0.7300000, 0.7500000, 0.7700000, 0.7900000, 0.8100000, 0.8300000, 0.8500000, 0.8700000, 0.8900000, 0.9100000, 0.9300000, 0.9500000, 0.9700000, 0.9900000, 1.0100000, 1.0300000, 1.0500000, 1.0700000, 1.0900000, 1.1100000, 1.1300000, 1.1500000, 1.1700000, 1.1900000, 1.2100000, 1.2300000, 1.2500000, 1.2700000, 1.2900000, 1.3100000, 1.3300000, 1.3500000, 1.3700000, 1.3900000, 1.4100000, 1.4300000, 1.4500000, 1.4700000, 1.4900000, 1.5100000, 1.5300000, 1.5500000, 1.5700000, 1.5900000, 1.6100000, 1.6300000, 1.6500000, 1.6700000, 1.6900000, 1.7100000, 1.7300000, 1.7500000, 1.7700000, 1.7900000, 1.8100000, 1.8300000, 1.8500000, 1.8700000, 1.8900000, 1.9100000, 1.9300000, 1.9500000, 1.9700000, 1.9900000, 2.0100000, 2.0300000, 2.0500000, 2.0700000, 2.0900000, 2.1100000, 2.1300000, 2.1500000, 2.1700000, 2.1900000, 2.2100000, 2.2300000, 2.2500000, 2.2700000, 2.2900000, 2.3100000, 2.3300000, 2.3500000, 2.3700000, 2.3900000, 2.4100000, 2.4300000, 2.4500000, 2.4700000, 2.4900000, 2.5100000, 2.5300000, 2.5500000, 2.5700000, 2.5900000, 2.6100000, 2.6300000, 2.6500000, 2.6700000, 2.6900000, 2.7100000, 2.7300000, 2.7500000, 2.7700000, 2.7900000, 2.8100000, 2.8300000, 2.8500000, 2.8700000, 2.8900000, 2.9100000, 2.9300000, 2.9500000, 2.9700000, 2.9900000, 3.0100000, 3.0300000, 3.0500000, 3.0700000, 3.0900000, 3.1100000, 3.1300000, 3.1500000, 3.1700000, 3.1900000, 3.2100000, 3.2300000, 3.2500000, 3.2700000, 3.2900000, 3.3100000, 3.3300000, 3.3500000, 3.3700000, 3.3900000, 3.4100000, 3.4300000, 3.4500000, 3.4700000, 3.4900000, 3.5100000, 3.5300000, 3.5500000, 3.5700000, 3.5900000, 3.6100000, 3.6300000, 3.6500000, 3.6700000, 3.6900000, 3.7100000, 3.7300000, 3.7500000, 3.7700000, 3.7900000, 3.8100000, 3.8300000, 3.8500000, 3.8700000, 3.8900000, 3.9100000, 3.9300000, 3.9500000, 3.9700000, 3.9900000, 4.0100000, 4.0300000, 4.0500000, 4.0700000, 4.0900000, 4.1100000, 4.1300000, 4.1500000, 4.1700000, 4.1900000, 4.2100000, 4.2300000, 4.2500000, 4.2700000, 4.2900000, 4.3100000, 4.3300000, 4.3500000, 4.3700000, 4.3900000, 4.4100000, 4.4300000, 4.4500000, 4.4700000, 4.4900000, 4.5100000, 4.5300000, 4.5500000, 4.5700000, 4.5900000, 4.6100000, 4.6300000, 4.6500000, 4.6700000, 4.6900000, 4.7100000, 4.7300000, 4.7500000, 4.7700000, 4.7900000, 4.8100000, 4.8300000, 4.8500000, 4.8700000, 4.8900000, 4.9100000, 4.9300000, 4.9500000, 4.9700000, 4.9900000, 5.0100000, 5.0300000, 5.0500000, 5.0700000, 5.0900000, 5.1100000, 5.1300000, 5.1500000, 5.1700000, 5.1900000, 5.2100000, 5.2300000, 5.2500000, 5.2700000, 5.2900000, 5.3100000, 5.3300000, 5.3500000, 5.3700000, 5.3900000, 5.4100000, 5.4300000, 5.4500000, 5.4700000, 5.4900000, 5.5100000, 5.5300000, 5.5500000, 5.5700000, 5.5900000, 5.6100000, 5.6300000, 5.6500000, 5.6700000, 5.6900000, 5.7100000, 5.7300000, 5.7500000, 5.7700000, 5.7900000, 5.8100000, 5.8300000, 5.8500000, 5.8700000, 5.8900000, 5.9100000, 5.9300000, 5.9500000, 5.9700000, 5.9900000, 6.0100000, 6.0300000, 6.0500000, 6.0700000, 6.0900000, 6.1100000, 6.1300000, 6.1500000, 6.1700000, 6.1900000, 6.2100000, 6.2300000, 6.2500000, 6.2700000, 6.2900000, 6.3100000, 6.3300000, 6.3500000, 6.3700000, 6.3900000, 6.4100000, 6.4300000, 6.4500000, 6.4700000, 6.4900000, 6.5100000, 6.5300000, 6.5500000, 6.5700000, 6.5900000, 6.6100000, 6.6300000, 6.6500000, 6.6700000, 6.6900000, 6.7100000, 6.7300000, 6.7500000, 6.7700000, 6.7900000, 6.8100000, 6.8300000, 6.8500000, 6.8700000, 6.8900000, 6.9100000, 6.9300000, 6.9500000, 6.9700000, 6.9900000, 7.0100000, 7.0300000, 7.0500000, 7.0700000, 7.0900000, 7.1100000, 7.1300000, 7.1500000, 7.1700000, 7.1900000, 7.2100000, 7.2300000, 7.2500000, 7.2700000, 7.2900000, 7.3100000, 7.3300000, 7.3500000, 7.3700000, 7.3900000, 7.4100000, 7.4300000, 7.4500000, 7.4700000, 7.4900000, 7.5100000, 7.5300000, 7.5500000, 7.5700000, 7.5900000, 7.6100000, 7.6300000, 7.6500000, 7.6700000, 7.6900000, 7.7100000, 7.7300000, 7.7500000, 7.7700000, 7.7900000, 7.8100000, 7.8300000, 7.8500000, 7.8700000, 7.8900000, 7.9100000, 7.9300000, 7.9500000, 7.9700000, 7.9900000, 8.0100000, 8.0300000, 8.0500000, 8.0700000, 8.0900000, 8.1100000, 8.1300000, 8.1500000, 8.1700000, 8.1900000, 8.2100000, 8.2300000, 8.2500000, 8.2700000, 8.2900000, 8.3100000, 8.3300000, 8.3500000, 8.3700000, 8.3900000, 8.4100000, 8.4300000, 8.4500000, 8.4700000, 8.4900000, 8.5100000, 8.5300000, 8.5500000, 8.5700000, 8.5900000, 8.6100000, 8.6300000, 8.6500000, 8.6700000, 8.6900000, 8.7100000, 8.7300000, 8.7500000, 8.7700000, 8.7900000, 8.8100000, 8.8300000, 8.8500000, 8.8700000, 8.8900000, 8.9100000, 8.9300000, 8.9500000, 8.9700000, 8.9900000, 9.0100000, 9.0300000, 9.0500000, 9.0700000, 9.0900000, 9.1100000, 9.1300000, 9.1500000, 9.1700000, 9.1900000, 9.2100000, 9.2300000, 9.2500000, 9.2700000, 9.2900000, 9.3100000, 9.3300000, 9.3500000, 9.3700000, 9.3900000, 9.4100000, 9.4300000, 9.4500000, 9.4700000, 9.4900000, 9.5100000, 9.5300000, 9.5500000, 9.5700000, 9.5900000, 9.6100000, 9.6300000, 9.6500000, 9.6700000, 9.6900000, 9.7100000, 9.7300000, 9.7500000, 9.7700000, 9.7900000, 9.8100000, 9.8300000, 9.8500000, 9.8700000, 9.8900000, 9.9100000, 9.9300000, 9.9500000, 9.9700000, 9.9900000])
  # SE_t = np.array([5.5202651, 3.0568192, 2.2659792, 1.8360280, 1.5509789, 1.3435837, 1.1851011, 1.0574761, 0.9512693, 0.8619317, 0.7853095, 0.7192153, 0.6608820, 0.6096077, 0.5634182, 0.5221652, 0.4846911, 0.4509921, 0.4204513, 0.3926222, 0.3666329, 0.3430010, 0.3212923, 0.3014545, 0.2826236, 0.2656413, 0.2493999, 0.2345191, 0.2208866, 0.2079489, 0.1958998, 0.1845900, 0.1741967, 0.1644579, 0.1554295, 0.1467549, 0.1386055, 0.1310163, 0.1237243, 0.1171964, 0.1107594, 0.1049353, 0.0993675, 0.0940865, 0.0892058, 0.0844640, 0.0799929, 0.0758424, 0.0719718, 0.0682141, 0.0647631, 0.0614221, 0.0582970, 0.0553217, 0.0525641, 0.0499019, 0.0473980, 0.0449691, 0.0426860, 0.0406201, 0.0385691, 0.0366567, 0.0348163, 0.0331066, 0.0314866, 0.0299059, 0.0284567, 0.0270738, 0.0257566, 0.0245303, 0.0233003, 0.0221791, 0.0210949, 0.0200807, 0.0191034, 0.0181889, 0.0173020, 0.0164860, 0.0156728, 0.0149274, 0.0142313, 0.0135314, 0.0129006, 0.0122865, 0.0116923, 0.0111403, 0.0106214, 0.0101173, 0.0096327, 0.0091767, 0.0087439, 0.0083393, 0.0079384, 0.0075693, 0.0072185, 0.0068802, 0.0065519, 0.0062469, 0.0059587, 0.0056785, 0.0054101, 0.0051598, 0.0049220, 0.0046895, 0.0044770, 0.0042626, 0.0040700, 0.0038816, 0.0037041, 0.0035354, 0.0033704, 0.0032154, 0.0030653, 0.0029241, 0.0027884, 0.0026600, 0.0025371, 0.0024241, 0.0023105, 0.0022030, 0.0021042, 0.0020026, 0.0019149, 0.0018270, 0.0017453, 0.0016665, 0.0015892, 0.0015163, 0.0014470, 0.0013834, 0.0013195, 0.0012601, 0.0012015, 0.0011476, 0.0010946, 0.0010458, 0.0009985, 0.0009530, 0.0009106, 0.0008683, 0.0008292, 0.0007919, 0.0007562, 0.0007222, 0.0006891, 0.0006590, 0.0006286, 0.0005998, 0.0005728, 0.0005470, 0.0005229, 0.0004990, 0.0004770, 0.0004554, 0.0004347, 0.0004156, 0.0003970, 0.0003791, 0.0003619, 0.0003458, 0.0003304, 0.0003156, 0.0003013, 0.0002875, 0.0002750, 0.0002628, 0.0002510, 0.0002394, 0.0002291, 0.0002187, 0.0002092, 0.0001997, 0.0001910, 0.0001824, 0.0001743, 0.0001665, 0.0001593, 0.0001523, 0.0001452, 0.0001388, 0.0001326, 0.0001268, 0.0001211, 0.0001157, 0.0001106, 0.0001057, 0.0001010, 0.0000965, 0.0000922, 0.0000883, 0.0000842, 0.0000805, 0.0000769, 0.0000735, 0.0000703, 0.0000672, 0.0000643, 0.0000613, 0.0000587, 0.0000561, 0.0000536, 0.0000513, 0.0000489, 0.0000469, 0.0000447, 0.0000427, 0.0000409, 0.0000391, 0.0000373, 0.0000357, 0.0000342, 0.0000326, 0.0000312, 0.0000299, 0.0000285, 0.0000273, 0.0000261, 0.0000249, 0.0000238, 0.0000228, 0.0000218, 0.0000208, 0.0000199, 0.0000190, 0.0000182, 0.0000174, 0.0000166, 0.0000159, 0.0000152, 0.0000145, 0.0000139, 0.0000133, 0.0000127, 0.0000122, 0.0000116, 0.0000111, 0.0000106, 0.0000102, 0.0000097, 0.0000093, 0.0000089, 0.0000085, 0.0000081, 0.0000078, 0.0000074, 0.0000071, 0.0000068, 0.0000065, 0.0000062, 0.0000059, 0.0000057, 0.0000054, 0.0000052, 0.0000050, 0.0000047, 0.0000045, 0.0000044, 0.0000042, 0.0000040, 0.0000038, 0.0000036, 0.0000035, 0.0000033, 0.0000032, 0.0000030, 0.0000029, 0.0000028, 0.0000027, 0.0000025, 0.0000024, 0.0000023, 0.0000022, 0.0000021, 0.0000020, 0.0000020, 0.0000019, 0.0000018, 0.0000017, 0.0000016, 0.0000016, 0.0000015, 0.0000014, 0.0000014, 0.0000013, 0.0000013, 0.0000012, 0.0000011, 0.0000011, 0.0000011, 0.0000010, 0.0000010, 0.0000009, 0.0000009, 0.0000008, 0.0000008, 0.0000008, 0.0000007, 0.0000007, 0.0000007, 0.0000006, 0.0000006, 0.0000006, 0.0000006, 0.0000005, 0.0000005, 0.0000005, 0.0000005, 0.0000005, 0.0000004, 0.0000004, 0.0000004, 0.0000004, 0.0000004, 0.0000003, 0.0000003, 0.0000003, 0.0000003, 0.0000003, 0.0000003, 0.0000003, 0.0000003, 0.0000002, 0.0000002, 0.0000002, 0.0000002, 0.0000002, 0.0000002, 0.0000002, 0.0000002, 0.0000002, 0.0000002, 0.0000002, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000])
  # SE_t -= a*np.exp(-l*t)*(np.pi*t)**-0.5

L = l
SE_t = a*np.exp(-L*t)*(np.pi*t)**-0.5

Cs_a = Cs_analytical(l, L, a)

# calculate weights
for i, val in enumerate(SE_t):
  for n in range(0, Cs.size):
    Cs[n] += val*e(n, t[i], l)

# approximation of histogram
S = t*0
for n, c in enumerate(Cs):
  S += c*e(n, t, l)

# analytical projection
S_a = t*0
for n, c in enumerate(Cs_a):
  S_a += c*e(n, t, l)

# just to match the amplitude
scaleFactor = np.mean(SE_t[:100]/S[:100])

###
### using exponential distribution
###

Cs_dist = np.zeros(Cs.size)
S_dist = t*0
for i in range(0, 10000):
  r = np.random.rand()
  T = - np.log(1 - r)/L
  for n in range(0, Cs_dist.size):
    Cs_dist[n] += e(n, T, l)

for n, c in enumerate(Cs_dist):
  S_dist += c*e(n, t, l)


print(Cs_a)
print(Cs*scaleFactor)

plt.plot(t, S*scaleFactor)
plt.plot()
plt.plot(t, S_a, color='magenta')
plt.plot(t, SE_t, ':r')
plt.plot(t, SE_t - S*scaleFactor, '--')
# plt.xlim(-0.1/l, 10/l)
plt.show()