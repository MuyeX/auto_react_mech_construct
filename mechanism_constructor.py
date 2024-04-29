#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 09:37:56 2023

@author: md1621
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from ODE_generator import make_system
import matplotlib.cm as cm
import pandas as pd
import os 

# Function that reads files from a directory and returns a dictionary
def read_files(directory):
    files = os.listdir(directory)
    data = {}

    for file in files:
        if not file.startswith('.'):
            data[file[:-4]] = pd.read_csv(os.path.join(directory, file), 
                                          header = 0, index_col = None,
                                          encoding = 'latin1').values
    return data

in_silico_data = read_files("exp_data")

r1 = 'A -> B + I'
reactions = [r1]
mechanism = make_system(reactions)
# The function executed below is called kinetic_model
exec(mechanism)

multistart = 10
number_parameters = 1
lower_bound = 0.0001
upper_bound = 10

# abc_obj = abc(sse, [(lower_bound, upper_bound) for i in range(number_parameters)])
# abc_obj.fit() 

# solution = abc_obj.get_solution()
# print('Initial guess = ', solution)

solution = np.array([0.1])

opt_val, opt_param = Opt_Rout(multistart, number_parameters, solution, lower_bound, \
    upper_bound, sse, kinetic_model)

print('MSE = ', opt_val)
print('Optimal parameters = ', opt_param)

# timesteps = 30
# time = np.linspace(0, 2, timesteps)
# t = [0, np.max(time)]
# t_eval = list(time)
# rate_constants = np.array([1.514])
# ic = np.array([4, 0, 0])
# solution = solve_ivp(kinetic_model, t, ic, t_eval = t_eval, method = "RK45", 
#                      args = rate_constants)
# print(solution.y.T)



