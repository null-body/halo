import math
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from sympy import polylog

def sigma_r(x, c, k):
    return - k * x * (x + 1) **2/(x **2 + c **2) * (1/2 * (c **2/x + 1/(1 + x) **2 + c **2/(1 + x) **2 + 2/(1 + x) + (6 * c **2)/(1 + x) + c **2 * np.log(x) - c **2 * np.log(1 + x) - (c **2 * np.log(1 + x))/x **2 + (4 * c **2 * np.log(1 + x))/x + (2 * np.log(1 + x))/(1 + x) + (2 * c **2 * np.log(1 + x))/(1 + x) - np.log(1 + x) **2 - 3 * c **2 * np.log(1 + x) **2 - 2 * (1 + 3 * c **2) * polylog(2, -x)) - 1/2 * (c + 1/(1 + c) **2 + c **2/(1 + c) **2 + 2/(1 + c) + (6 * c **2)/(1 + c) + c **2 * np.log(c) - np.log(1 + c) + 4 * c * np.log(1 + c) - c **2 * np.log(1 + c) + (2 * np.log(1 + c))/(1 + c) + (2 * c **2 * np.log(1 + c))/(1 + c) - np.log(1 + c) **2 - 3 * c **2 * np.log(1 + c) **2 - 2 * (1 + 3 * c **2) * polylog(2, -c)))

def sigma_0(x, c, k):
    return - 1/2 * k * x * (x + 1) **2 * ((1/x + 1/(1 + x) **2 + 6/(1 + x) + np.log(x/c) - np.log(1 + x) - np.log(1 + x)/x **2 + (4 * np.log(1 + x))/x + (2 * np.log(1 + x))/(1 + x) - 3 * np.log(1 + x) **2 - 6 * polylog(2, -x)) - (1/c + 1/(1 + c) **2 + 6/(1 + c) - np.log(1 + c) - np.log(1 + c)/c **2 + (4 * np.log(1 + c))/c + (2 * np.log(1 + c))/(1 + c) - 3 * np.log(1 + c) **2 - 6 * polylog(2, -c)))

con_ = 5
spt_ = 200
gap_ = con_/spt_
kap_ = 1

x0, y0 = [], []
y_iso, y_p = [], []

for i in range(spt_):
    x0.append(0.0001 + gap_ * i)
    y0.append(float(sigma_r(0.0001 + gap_ * i, con_, kap_)))
    y_iso.append(float(sigma_0(0.0001 + gap_ * i, con_, kap_)))
    y_p.append(abs(1 - float(sigma_0(0.0001 + gap_ * i, con_, kap_)) / float(sigma_r(0.0001 + gap_ * i, con_, kap_))))

avg = sum(y_p)/len(y_p)
    
plt.figure(figsize = (16, 6))

rc = {"font.family" : "serif", "mathtext.fontset" : "stix"}
plt.rcParams.update(rc)
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]

plt.subplot(121)
plt.plot(x0, y0, "k-", linewidth = 1, label = '$Anisotropic$')
plt.plot(x0, y_iso, "k--", linewidth = 1, label = '$Isotropic$')

plt.title('$Radial$ $velocity$ $dispersion$ $at$ $c=$'+str(con_), fontsize = 18)
plt.xlabel('$r / r_c$', fontsize = 16)
plt.ylabel('$\sigma^2_r/\kappa$', fontsize = 16)
#plt.yscale('log')

plt.legend(loc = 'best')

plt.subplot(122)
plt.plot(x0, y_p, "k-", linewidth = 1)

plt.title('$\sigma^2_r - \sigma^2_0$', fontsize = 18)
plt.xlabel('$r / r_c$', fontsize = 16)
plt.ylabel('$\sigma^2/\kappa$', fontsize = 16)
plt.axhline(y = avg, color = 'k', linestyle = '-.', linewidth = 1, label = '{0:.3g}'.format(avg))

plt.legend(loc = 'best')

plt.show()
