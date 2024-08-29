#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 11:14:05 2023

@author: md1621
"""

"##############################################################################"
"######################## Importing important packages ########################"
"##############################################################################"

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.cm as cm
from scipy.optimize import minimize

np.random.seed(1998)


"##############################################################################"
"###################### Model-Based Design of Experiments #####################"
"##############################################################################"

def kinetic_model(x, init):
    CA,CB,CC,CD,CE = init
    k1, k2, k3 = np.array([1.56885688, 3.70554929, 1.95498632])
    dAdt = - k1*CA
    dBdt = k1*CA + k2*CD + k3*CE
    dCdt = k3*CE
    dDdt = k1*CA - k2*CD
    dEdt = k2*CD - k3*CE
    return dAdt,dBdt,dCdt,dDdt,dEdt


def kinetic_model_2(x, init):
    CA,CB,CC,CD,CE,CF = init
    k1, k2, k3, k4 = np.array([1.514, 5.259, 9.352, 2.359])
    dAdt = - k1*CA
    dBdt = k1*CA + k2*CD + k4*CF
    dCdt = k4*CF
    dDdt = k1*CA - k2*CD
    dEdt = k2*CD - k3*CE
    dFdt = k3*CE - k4*CF
    return dAdt,dBdt,dCdt,dDdt,dEdt,dFdt


def MBDoE(ic, time):
    t = [0, np.max(time)]
    t_eval = list(time)
    ic_1 = np.append(ic, [0] * 2)
    ic_2 = np.append(ic, [0] * 3)
    solution_1 = solve_ivp(kinetic_model, t, ic_1, t_eval = t_eval, method = "RK45")
    solution_2 = solve_ivp(kinetic_model_2, t, ic_2, t_eval = t_eval, method = "RK45")
    diff_1 = np.array([solution_1.y.T[:, 0], solution_1.y.T[:, 1], solution_1.y.T[:, 2]]).T
    diff_2 = np.array([solution_2.y.T[:, 0], solution_2.y.T[:, 1], solution_2.y.T[:, 2]]).T
    difference = -np.sum((diff_1 - diff_2)**2)
    return difference


def Opt_Rout(multistart, number_parameters, lower_bound, upper_bound, to_opt, \
    time):
    localsol = np.empty([multistart, number_parameters])
    localval = np.empty([multistart, 1])
    boundss = tuple([(lower_bound[i], upper_bound[i]) for i in range(len(lower_bound))])
    
    for i in range(multistart):
        x0 = np.random.uniform(lower_bound, upper_bound, size = number_parameters)
        res = minimize(to_opt, x0, args = (time), \
                        method = 'L-BFGS-B', bounds = boundss)
        localsol[i] = res.x
        localval[i] = res.fun

    minindex = np.argmin(localval)
    opt_val = localval[minindex]
    opt_param = localsol[minindex]
    
    return opt_val, opt_param


"##############################################################################"
"########################## MBDoE on Competing Models #########################"
"##############################################################################"

multistart = 10
number_parameters = 3
lower_bound = np.array([4, 0, 0])
upper_bound = np.array([6, 2, 1])
to_opt = MBDoE
timesteps = 30
time = np.linspace(0, 2, timesteps)


# THESE MODELS AND THE REAL ONE ARE INDISTINGUISHABLE WITHIN BOUNDS
a, b = Opt_Rout(multistart, number_parameters, lower_bound, upper_bound, to_opt, \
    time)

print('Optimal difference: ', -1 * a)
print('Optimal experiment: ', b)


"##############################################################################"
"########################### Plot MBDoE Experiment ############################"
"##############################################################################"

color_1 = cm.plasma(np.linspace(0, 1, 3))
marker = ['o', 'o', 'o', 'o', 'o', 'o']
species = ['A', 'B', 'C', 'D', 'E', 'F']


# Plotting the in-silico data for visualisation
fig, ax = plt.subplots()
# ax.set_title("Experiment " + str(i + 1))
ax.set_ylabel("Concentrations $(M)$", fontsize = 18)
ax.set_xlabel("Time $(h)$", fontsize = 18)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.tick_params(axis = 'both', which = 'major', labelsize = 18)

timesteps = 30
time = np.linspace(0, 2, timesteps)
t = [0, np.max(time)]
t_eval = list(time)


ic = b
ic_1 = np.append(ic, [0] * 2)
ic_2 = np.append(ic, [0] * 3)
solution_1 = solve_ivp(kinetic_model, t, ic_1, t_eval = t_eval, method = "RK45")
solution_2 = solve_ivp(kinetic_model_2, t, ic_2, t_eval = t_eval, method = "RK45")


for j in np.array([0, 1, 2]):
    y_1 = solution_1.y[j]
    y_2 = solution_2.y[j]
    ax.plot(time, y_1, marker[j], linestyle = '--', markersize = 4, label = species[j], color = color_1[j])
    ax.plot(time, y_2, marker[j], linestyle = '-', markersize = 4, label = species[j], color = color_1[j])

ax.grid(alpha = 0.5)
ax.legend(loc='upper right', fontsize = 15)


plt.show()