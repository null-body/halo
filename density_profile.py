import math
import matplotlib.pyplot as plt
import numpy as np
from pylab import *

rc('axes', linewidth = 1.4)
mpl.rcParams['xtick.major.width'] = 1.4
mpl.rcParams['xtick.minor.width'] = 1.4
mpl.rcParams['ytick.major.width'] = 1.4
mpl.rcParams['ytick.minor.width'] = 1.4
mpl.rcParams['xtick.labelsize'] = 15
mpl.rcParams['ytick.labelsize'] = 15

ein = []
m99 = []
nfw = []
rr = np.logspace(-3, 3, 1024)

for i in range(len(rr)):
    nfw.append(3000/rr[i]/(1 + rr[i]) **2)
    m99.append(1500/rr[i] **(1.5)/(1 + rr[i] **(1.5)))
    ein.append(750 * np.exp((-2/0.17) * (rr[i] **(0.17) - 1)))

fig = plt.figure(figsize = (9, 7.5))

plt.rcParams.update({"font.family" : "serif", "mathtext.fontset" : "stix"})
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]

plt.loglog(rr, nfw, 'k-', linewidth = 1, label = 'NFW')
plt.loglog(rr, m99, 'k:', linewidth = 1.5, label = 'M99')
plt.loglog(rr, ein, 'k--', linewidth = 1, label = 'Einasto')

plt.ylim(1.e-4, 1.e8)
plt.xlim(1.e-3, 1.e3)
plt.xlabel('$r$ / $r_c$', fontsize = 22.5)
plt.ylabel('$\\rho$ / $\\rho_b$', fontsize = 22.5)
plt.legend(loc = 'best', prop={'size': 15})
plt.show()
