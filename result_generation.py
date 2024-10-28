#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 08:54:25 2024

@author: md1621
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from ODE_generator import make_system
import matplotlib.cm as cm
import pandas as pd
from parameter_estimation import sse, callback, timeout_handler, Opt_Rout, evaluate
from read_data import read_files, reverse_dict
np.random.seed(1998)


# name_file = "exp_data_fruc_HMF"
# name_file = "exp_data_hypoth"
name_file = "exp_data_aldol_condensation"
place_holder = read_files(name_file)
in_silico_data = reverse_dict(place_holder)


solution_1 = np.array( [[-1, -1,  1,  1]]) 
solution_2 = np.array( [[-1,  0,  0,  0,  1],[ 0, -1,  1,  1, -1]])
solution_3 = np.array( [[-1,  0,  0,  0,  1,  0],[ 0, -1,  0,  0, -1,  1],[ 0,  0,  1,  1,  0, -1]])
solution_4 = np.array( [[-1,  0,  0,  0,  0,  0,  1],[ 0, -1,  0,  0,  1,  1, -1],[ 0,  0,  0,  1, -1,  0,  0],[ 0,  0,  1,  0,  0, -1,  0]])
model_pred, opt_param, nll, aic = evaluate(solution_4)

num_observable_species = 4

color_1 = ['salmon', 'royalblue', 'darkviolet', 'limegreen']
marker = ['o' for i in range(num_observable_species)]

num_exp = 5
timesteps = 30
time = np.linspace(0, 10, timesteps)
species = ['A', 'B', 'C', 'D', 'E']

    
# Plotting the in-silico data for visualisation
for i in range(num_exp):
    fig, ax = plt.subplots()
    # ax.set_title("Experiment " + str(i + 1))
    ax.set_ylabel("Concentrations $(M)$", fontsize = 18)
    ax.set_xlabel("Time $(h)$", fontsize = 18)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 18)
    
    for j in range(num_observable_species):
        y = model_pred["exp_" + str(i + 1)][:,j]
        ax.plot(time, y, label = species[j], color = color_1[j])
        yy = in_silico_data["exp_" + str(i + 1)][:,j]
        ax.plot(time, yy, marker[j], markersize = 4, color = color_1[j])

    ax.grid(alpha = 0.5)
    ax.legend(loc='upper right', fontsize = 15)
    
    if i == 1:
        file_path = 'Aldol_Results_Iteration_4.png'
        plt.savefig(file_path, dpi = 600, bbox_inches = "tight")

plt.show()

print(aic)
print(opt_param)