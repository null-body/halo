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
spt_ = 150
gap_ = con_/spt_
kap_ = 1
x0, y0 = [], []

for i in range(spt_):
    x0.append(0.0001 + gap_ * i)
    y0.append(float(sigma_r(0.0001 + gap_ * i, con_, kap_)))

y, y_power = [], []

for i in range(len(x0)):
    y.append(7.5 * x0[i] **(-1.875))
    y_power.append(1/x0[i]/(x0[i] + 1) **2 * y0[i] **(-3/2))

plt.figure(figsize = (8, 6))

rc = {"font.family" : "serif", "mathtext.fontset" : "stix"}
plt.rcParams.update(rc)
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]

plt.plot(x0, y_power, "k-", linewidth = 1)
plt.plot(x0, y, "k--", linewidth = 1, label = "$7.5 \; x^{-1.875}$")

plt.title('$\\rho / \sigma^3_r$ at $c=$'+str(con_), fontsize = 20)
plt.xlabel('$r / r_c$', fontsize = 16)
plt.ylabel('$\\rho \sigma^{-3}_r / \kappa$', fontsize = 16)
plt.yscale('log')
plt.xscale('log')

plt.legend(loc = 'best')

plt.show()
