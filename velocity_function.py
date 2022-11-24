import math
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from scipy.optimize import curve_fit
from sympy import polylog

def sigma_r(x, c, k):
    return - k * x * (x + 1) **2/(x **2 + c **2) * (1/2 * (c **2/x + 1/(1 + x) **2 + c **2/(1 + x) **2 + 2/(1 + x) + (6 * c **2)/(1 + x) + c **2 * np.log(x) - c **2 * np.log(1 + x) - (c **2 * np.log(1 + x))/x **2 + (4 * c **2 * np.log(1 + x))/x + (2 * np.log(1 + x))/(1 + x) + (2 * c **2 * np.log(1 + x))/(1 + x) - np.log(1 + x) **2 - 3 * c **2 * np.log(1 + x) **2 - 2 * (1 + 3 * c **2) * polylog(2, -x)) - 1/2 * (c + 1/(1 + c) **2 + c **2/(1 + c) **2 + 2/(1 + c) + (6 * c **2)/(1 + c) + c **2 * np.log(c) - np.log(1 + c) + 4 * c * np.log(1 + c) - c **2 * np.log(1 + c) + (2 * np.log(1 + c))/(1 + c) + (2 * c **2 * np.log(1 + c))/(1 + c) - np.log(1 + c) **2 - 3 * c **2 * np.log(1 + c) **2 - 2 * (1 + 3 * c **2) * polylog(2, -c)))

def sigma_0(x, c, k):
    return - 1/2 * k * x * (x + 1) **2 * ((1/x + 1/(1 + x) **2 + 6/(1 + x) + np.log(x/c) - np.log(1 + x) - np.log(1 + x)/x **2 + (4 * np.log(1 + x))/x + (2 * np.log(1 + x))/(1 + x) - 3 * np.log(1 + x) **2 - 6 * polylog(2, -x)) - (1/c + 1/(1 + c) **2 + 6/(1 + c) - np.log(1 + c) - np.log(1 + c)/c **2 + (4 * np.log(1 + c))/c + (2 * np.log(1 + c))/(1 + c) - 3 * np.log(1 + c) **2 - 6 * polylog(2, -c)))

def norm(x, mea, amp, std):
    return amp * np.exp(-(x - mea) **2 / (2 * std **2))

def norm_v(x, std):
    return 1 / np.sqrt(2 * np.pi * std) * np.exp(-1/2 * (x) **2 / std)

con_ = 5
kap_ = 1

x, x_v = np.arange(-con_, con_, con_/50), np.arange(-2, 2, 0.01)
x_len, x_v_len = len(x), len(x_v)
y, v = [[0 for i in range(x_len)] for j in range(9)], [[0 for i in range(x_v_len)] for j in range(9)]
#y_iso, v_iso = [[0 for i in range(x_len)] for j in range(9)], [[0 for i in range(x_v_len)] for j in range(9)]
rho = []

for j in range(9):
    r_p = 0.1 * (j + 1)
    for i in range(x_len):
        if (x[i] **2 + r_p **2 * con_ **2 <= con_ **2):
            rho.append(1/(np.sqrt(x[i] **2 + r_p **2 * con_ **2) * (np.sqrt(x[i] **2 + r_p **2 * con_ **2) + 1) **2))
            y[j][i] = float(((x[i] **2 + con_ **2)/(x[i] **2 + r_p **2 * con_ **2 + con_ **2) * sigma_r(np.sqrt(x[i] **2 + r_p **2 * con_ **2), con_, kap_)))
#            y_iso[j][i] = float(sigma_0(np.sqrt(x[i]**2 + r_p**2 * con_**2), con_, kap_))
        else:
            rho.append(0)

y[5][90], y[7][80] = 0, 0

for k in range(9):
    rho_sum = float(sum(rho[k * x_len:(k + 1) * x_len]))
    for i in range(x_v_len):
        for j in range(x_len):
            if y[k][j] != 0:
                v[k][i] = v[k][i] + (rho[j + k * x_len] * norm_v(x_v[i], y[k][j])) / rho_sum
#                v_iso[k][i] = v_iso[k][i] + (rho[j + k * x_len] * norm_v(x_v[i], y_iso[k][j])) / rho_sum

plt.figure(figsize = (16, 6))

plt.rcParams.update({"font.family" : "serif", "mathtext.fontset" : "stix"})
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]

plt.subplot(121)

for k in range(9):
    plt.plot(x_v, v[k], "r-", linewidth = 1.5)
    popt, popv = curve_fit(norm, x_v, v[k])
    y_gau = norm(x_v, popt[0], popt[1], popt[2])
    plt.plot(x_v, y_gau, 'g-', linewidth = 1)

#    plt.plot(x_v, v_iso[k], "r--", linewidth = 1.5)
#    popt, popv = curve_fit(norm, x_v, v_iso[k])
#    y_gau = norm(x_v, popt[0], popt[1], popt[2])
#    plt.plot(x_v, y_gau, 'g--', linewidth = 1)

plt.title('$R_p = 0.1,...,0.9 \; R_v \; at \; c=$'+str(con_), fontsize = 20)
plt.xlabel('$v$', fontsize = 16)
plt.ylabel('$f(v)$', fontsize = 16)
plt.xticks(fontsize = 10)
plt.yticks(fontsize = 10)

plt.xlim(-1, 1)
plt.ylim(0, 6)

plt.subplot(122)

for k in range(9):
    plt.plot(x_v, v[k], "r-", linewidth = 1.5)
    popt, popv = curve_fit(norm, x_v, v[k])
    y_gau = norm(x_v, popt[0], popt[1], popt[2])
    plt.plot(x_v, y_gau, 'g-', linewidth = 1)

#    plt.plot(x_v, v_iso[k], "r--", linewidth = 1.5)
#    popt, popv = curve_fit(norm, x_v, v_iso[k])
#    y_gau = norm(x_v, popt[0], popt[1], popt[2])
#    plt.plot(x_v, y_gau, 'g--', linewidth = 1)

plt.title('$R_p = 0.1,...,0.9 \; R_v \; at \; c=$'+str(con_), fontsize = 20)
plt.xlabel('$v$', fontsize = 16)
plt.ylabel('$f(v)$', fontsize = 16)
plt.xticks(fontsize = 10)
plt.yticks(fontsize = 10)

plt.xlim(-1, 1)
plt.ylim(1.e-5, 1.e+1)
plt.yscale('log')

plt.show()
