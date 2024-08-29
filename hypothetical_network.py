#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 11:39:51 2024

@author: md1621
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from ODE_generator import make_system
import matplotlib.cm as cm
import pandas as pd
np.random.seed(1998)

r1 = 'A + A -> B'
r2 = 'A -> D'
r3 = 'D -> E'
r4 = 'A + E -> C'
reactions = [r1, r2, r3, r4]
mechanism = make_system(reactions)
print(mechanism)

def kinetic_model(x, init, k1, k2, k3, k4):
    CA,CB,CC,CD,CE = init
    dAdt = - k1*CA*CA - k2*CA - k4*CA*CE
    dBdt = k1*CA*CA
    dCdt = k4*CA*CE
    dDdt = k2*CA - k3*CD
    dEdt = k3*CD - k4*CA*CE
    return dAdt,dBdt,dCdt,dDdt,dEdt

# Plotting the data given
num_observable_species = 3
species = ['A', 'B', 'C', 'D', 'E']

initial_conditions = {
    "ic_1": np.array([10, 0, 2, 0, 0]),
    "ic_2": np.array([10, 2, 0, 0, 0]),
    "ic_3": np.array([10, 2, 2, 0, 0]),
    "ic_4": np.array([5 , 0, 0, 0, 0]),
    "ic_5": np.array([10, 0, 0, 0, 0])
    }

rate_constants = np.array([0.1, 0.2, 0.13, 0.25])
    
num_exp = len(initial_conditions)
num_species = len(species)

timesteps = 30
time = np.linspace(0, 10, timesteps)
t = [0, np.max(time)]
t_eval = list(time)
STD = 0.15
noise = [np.random.normal(0, STD, size = (num_species, timesteps)) for i in range(num_exp)]
in_silico_data = {}
no_noise_data = {}
obs_data = {}

for i in range(num_exp):
    ic = initial_conditions["ic_" + str(i + 1)]
    solution = solve_ivp(kinetic_model, t, ic, t_eval = t_eval, method = "RK45", 
                          args = rate_constants)
    in_silico_data["exp_" + str(i + 1)] = np.clip(solution.y + noise[i], 0, 1e99)
    no_noise_data["exp_" + str(i + 1)] = solution.y
    obs_data["exp_" + str(i + 1)] = in_silico_data["exp_" + str(i + 1)][[j for j in range(num_observable_species)]]


def dict_to_csv(input_dict, filename):
    # For CSV, each key will be saved as a separate file, hence adding the key to the filename.
    for key, value in input_dict.items():
        if isinstance(value, np.ndarray):
            value = value.T
            value = pd.DataFrame(value)
        
        csv_filename = filename + key + '.csv'
        value.to_csv(csv_filename, index=False)


dict_to_csv(obs_data, 'exp_data_hypoth/')

color_1 = cm.plasma(np.linspace(0, 1, num_observable_species))
marker = ['o' for i in range(num_observable_species)]

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
        y = in_silico_data["exp_" + str(i + 1)][j]
        ax.plot(time, y, marker[j], markersize = 4, label = species[j], color = color_1[j])
        # y = no_noise_data["exp_" + str(i + 1)][j]
        # ax.plot(time, y, label = species[j], color = color_1[j])
    
    ax.grid(alpha = 0.5)
    ax.legend(loc='upper right', fontsize = 15)
    
    file_path = 'Hypoth_Experiment_' + str(i + 1) +'.png'
    plt.savefig(file_path, dpi = 600, bbox_inches = "tight")

plt.show()