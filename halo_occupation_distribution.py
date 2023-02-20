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

N_cen = []
N_sat = []
N_all = []
M = np.logspace(11, 15, 1024)

M_min = 10 **11.73
M_0 = 10 **12.09
M_prm = 10 **12.87
sig = 0.32
alp = 0.96

for i in range(len(M)):
    N_cen.append(1/2 * (1 + math.erf((np.log10(M[i]) - np.log10(M_min))/sig)))
    if M[i] > M_0:
        N_sat.append(((M[i] - M_0)/M_prm) **alp)
    else:
        N_sat.append(0)
    N_all.append(N_cen[i] + N_sat[i])

fig = plt.figure(figsize = (12, 9))

plt.rcParams.update({"font.family" : "serif", "mathtext.fontset" : "stix"})
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]

plt.loglog(M, N_cen, 'k--', linewidth = 1, label = 'Central')
plt.loglog(M, N_sat, 'k:', linewidth = 1.5, label = 'Satellite')
plt.loglog(M, N_all, 'k-', linewidth = 1, label = 'Combined')

plt.ylim(np.sqrt(10) * 1.e-2, np.sqrt(10) * 1.e1)
plt.xlim(1.e11, 1.e15)
plt.ylabel('$\langle N \\rangle$', fontsize = 25)
plt.xlabel('$M$ $[\\rm M_\odot]$', fontsize = 25)
plt.legend(loc = 'best', prop = {'size': 15})

plt.show()
