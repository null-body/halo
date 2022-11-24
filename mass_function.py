import matplotlib.pyplot as plt
import numpy as np
import math
from pylab import *
from matplotlib import cm
from sympy import polylog
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

rc('axes', linewidth=1.4)
mpl.rcParams['xtick.major.width'] = 1.4
mpl.rcParams['xtick.minor.width'] = 1.4
mpl.rcParams['ytick.major.width'] = 1.4
mpl.rcParams['ytick.minor.width'] = 1.4

rm = np.zeros(999)
rmult = np.zeros(999)

f = open('C:\\Users\\lexgu\\Desktop\\MPhys\\Data\\output.txt')

plt.figure(figsize = (9, 6))

rc = {"font.family" : "serif", "mathtext.fontset" : "stix"}
plt.rcParams.update(rc)
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]

for iz in range(11):
    f.readline()
    for i in range(999):
        line = f.readline()
        columns = line.split()
        rm[i] = float(columns[1])
        rmult[i] = float(columns[3])
    plt.loglog(rm, rmult, 'r-', linewidth = 2.5)
f.close()

plt.xlim((10 **9, 10 **16))
plt.ylim((10 **(-5), 0.1))

plt.xlabel('$M$  [$M_\odot$]')
plt.ylabel('$M^2f(M)/\\rho_0$ ')
plt.title('$z=0,1,\dots,10$', y = 1.03)

plt.show()
