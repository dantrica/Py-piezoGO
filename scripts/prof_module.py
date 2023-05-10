# -*- coding: utf-8 -*-
"""
Created on Sat Apr 10 15:21:32 2021

@author: Silvi_Sami
"""

import numpy as np
import matplotlib.pyplot as plt 
from lmfit import minimize, Model, Parameters
from scipy.optimize import basinhopping 
from scipy.interpolate import interp1d
 
data = np.loadtxt('G:/Mi unidad/Colab_DanielTriana/data/20200622/Z_p59_1.txt', skiprows=1)
f = np.array(data[:, 1])
Z = np.array(data[:, 2]) - 1j*np.array(data[:, 3])

# Zi = interp1d(f, Z)
# f = np.linspace(0.1, 1e6, int(1e6))
# Z = Zi(f)
# # plt.plot(Zi.real, -Zi.imag)
# plt.plot(Z.real, -Z.imag)

# s = I == min(I)
# ss = R <= R[s]
# f = f[ss]
# R = R[ss]
# I = I[ss]

rho_0 = 14000*(np.pi*(15e-3)**2)/20e-3
M_f_1 = 1e-3
tau_1 = 1e-1
c_1 = 0.1
M_f_2 = 1e-3
tau_2 = 1e-3 
c_2 = 0.5
M_f_3 = 1e-3
tau_3 = 1e-5 
c_3 = 0.9

# def gemtip_sh(f, rho_0, M_f_1, tau_1, c_1, M_f_2, tau_2, c_2, M_f_3, tau_3, c_3):
#     w = 2*np.pi*np.array(f)
#     p = np.array([M_f_1, tau_1, c_1, M_f_2, tau_2, c_2, M_f_3, tau_3, c_3])
#     params = p.reshape((3, 3))
#     sum_elements = 0
#     for k in params:
#         M_f_l = k[0]; tau_l = k[1]; c_l = k[2]
#         sum_elements += M_f_l * (1 - 1 / ( 1 + (1j * w * tau_l) ** c_l ))
#     sigma_e = (1 + sum_elements)/rho_0
#     Z = 20e-3/sigma_e/(np.pi*(15e-3)**2)
#     return Z

def gemtip_sh(f, rho_0, M_f_1, tau_1, c_1):
    w = 2*np.pi*np.array(f)
    sum_elements = M_f_1 * (1 - 1 / ( 1 + (1j * w * tau_1) ** c_1 ))
    sigma_e = (1 + sum_elements)/rho_0
    Z = 20e-3/sigma_e/(np.pi*(15e-3)**2)
    return Z

gmodel = Model(gemtip_sh)
params =Parameters()
params.add(gmodel.param_names[0], value=rho_0, vary=False)
params.add(gmodel.param_names[1], value=M_f_1, vary=False)
params.add(gmodel.param_names[2], value=tau_1, vary=True)
params.add(gmodel.param_names[3], value=c_1, vary=False)
# params.add(gmodel.param_names[4], value=M_f_2, vary=True)
# params.add(gmodel.param_names[5], value=tau_2, vary=True)
# params.add(gmodel.param_names[6], value=c_2, vary=True)
# params.add(gmodel.param_names[7], value=M_f_3, vary=True)
# params.add(gmodel.param_names[8], value=tau_3, vary=True)
# params.add(gmodel.param_names[9], value=c_3, vary=True)

result = gmodel.fit(data=Z, f=f, nan_policy='omit')


# plt.figure(dpi=120)
# plt.plot(Z.real, -Z.imag)