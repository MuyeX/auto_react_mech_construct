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
    CA,CB,CD,CE,CC = init
    k1, k2, k3 = np.array([1.56944133, 4.24332355, 1.91130801])
    dAdt = - k1*CA
    dBdt = k1*CA + k2*CD + k3*CE
    dDdt = k1*CA - k2*CD
    dEdt = k2*CD - k3*CE
    dCdt = k3*CE
    return dAdt,dBdt,dDdt,dEdt,dCdt

def kinetic_model_2(x, init):
    CA,CB,CD,CE,CF,CC = init
    k1, k2, k3, k4 = np.array([1.514, 5.259, 2.359, 9.352])
    dAdt = - k1*CA
    dBdt = k1*CA + k2*CD + k4*CF
    dDdt = k1*CA - k2*CD
    dEdt = k2*CD - k3*CE
    dFdt = k3*CE - k4*CF
    dCdt = k4*CF
    return dAdt,dBdt,dDdt,dEdt,dFdt,dCdt

def MBDoE(ic, time):
    t = [0, np.max(time)]
    t_eval = list(time)
    ic_1 = np.insert(ic, 2, [0] * 2)
    ic_2 = np.insert(ic, 2, [0] * 3)
    solution_1 = solve_ivp(kinetic_model, t, ic_1, t_eval = t_eval, method = "RK45")
    solution_2 = solve_ivp(kinetic_model_2, t, ic_2, t_eval = t_eval, method = "RK45")
    diff_1 = np.array([solution_1.y.T[:, 0], solution_1.y.T[:, 1], solution_1.y.T[:, -1]]).T
    diff_2 = np.array([solution_2.y.T[:, 0], solution_2.y.T[:, 1], solution_2.y.T[:, -1]]).T
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

multistart = 100
number_parameters = 3
lower_bound = np.array([1, 0, 0])
upper_bound = np.array([10, 10, 10])
to_opt = MBDoE
timesteps = 30
time = np.linspace(0, 2, timesteps)


# THESE MODELS AND THE REAL ONE ARE INDISTINGUISHABLE WITHIN BOUNDS
a, b = Opt_Rout(multistart, number_parameters, lower_bound, upper_bound, to_opt, \
    time)

print('Optimal experiment: ', b)