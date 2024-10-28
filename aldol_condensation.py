#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 14:54:22 2024

@author: md1621
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from ODE_generator import make_system
import matplotlib.cm as cm
import pandas as pd
np.random.seed(1998)

r1 = 'A -> E'
r2 = 'B + E -> F'
r3 = 'F -> C + D'
reactions = [r1, r2, r3]
mechanism = make_system(reactions)
print(mechanism)

def kinetic_model(x, init, k1, k2, k3):
    CA,CB,CC,CD,CE,CF = init
    dAdt = - k1*CA
    dBdt = - k2*CE*CB
    dCdt = k3*CF
    dDdt = k3*CF
    dEdt = k1*CA - k2*CE*CB
    dFdt = k2*CE*CB - k3*CF
    return dAdt,dBdt,dCdt,dDdt,dEdt,dFdt

# Plotting the data given
num_observable_species = 4
species = ['A', 'B', 'C', 'D', 'E', 'F']

initial_conditions = {
    "ic_1": np.array([5 , 10, 0, 0, 0, 0]),
    "ic_2": np.array([5 , 5 , 2, 0, 0, 0]),
    "ic_3": np.array([5 , 10, 0, 2, 0, 0]),
    "ic_4": np.array([10, 10, 0, 2, 0, 0]),
    "ic_5": np.array([10, 10, 2, 2, 0, 0])
    }

rate_constants = np.array([0.7593, 0.2928, 0.6813])
    
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
    obs_data["exp_" + str(i + 1)] = in_silico_data["exp_" + str(i + 1)][[i for i in range(num_observable_species)]]


def dict_to_csv(input_dict, filename):
    # For CSV, each key will be saved as a separate file, hence adding the key to the filename.
    for key, value in input_dict.items():
        if isinstance(value, np.ndarray):
            value = value.T
            value = pd.DataFrame(value)
        
        csv_filename = filename + key + '.csv'
        value.to_csv(csv_filename, index=False)


dict_to_csv(obs_data, 'exp_data_aldol_condensation/')

color_1 = ['salmon', 'royalblue', 'darkviolet', 'limegreen']
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
        # yy = no_noise_data["exp_" + str(i + 1)][j]
        # ax.plot(time, yy, color = color_1[j])
    
    ax.grid(alpha = 0.5)
    ax.legend(loc=(0.78,0.18), fontsize = 15)
    
    file_path = 'Aldol_Condensation_Experiment_' + str(i + 1) +'.png'
    plt.savefig(file_path, dpi = 600, bbox_inches = "tight")

plt.show()


# from metrics import NLL_mechanism, Akaike

# exp_data = np.vstack(list(in_silico_data.values()))
# model_pred = np.vstack(list(no_noise_data.values()))

# a = NLL_mechanism(exp_data, model_pred)
# b = Akaike(a, np.array([1,2,3]))
# print(b)
