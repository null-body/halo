import math
import matplotlib.pyplot as plt
import numpy as np
from pylab import *

rc('axes', linewidth = 1.4)
mpl.rcParams['xtick.major.width'] = 1.4
mpl.rcParams['xtick.minor.width'] = 1.4
mpl.rcParams['ytick.major.width'] = 1.4
mpl.rcParams['ytick.minor.width'] = 1.4

rm = np.zeros(999)
rmult = np.zeros(999)
ff = np.zeros(999)
fm = np.zeros(999)
N = np.zeros(999)

om_m0 = 0.3
rho0 = 2.78e11 * om_m0

f = open('file_location')

fig = plt.figure(figsize = (9, 7.5))

plt.rcParams.update({"font.family" : "serif", "mathtext.fontset" : "stix"})
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]

for iz in range(6):
    f.readline()
    for i in range(999):
        line = f.readline()
        columns = line.split()
        rm[i] = float(columns[1])
        rmult[i] = float(columns[3])
    plt.loglog(rm, rmult, 'k-', linewidth = 1)
f.close()

plt.xscale('log')

plt.xlim((1.e9, 1.e16))
plt.ylim((3.e-4, 3.e-1))
plt.title('$z$ $=$ $0, 1, 2, 3, 4, 5$', fontsize = 22.5)
plt.xlabel('$M$ $(h^{-1} \, M_\odot)$', fontsize = 22.5)
plt.ylabel('$M^2$ $n(M)$ / $\\rho_0$ ', fontsize = 22.5)

plt.show()
