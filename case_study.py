#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 13:22:34 2023

@author: md1621
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from mechanism_generator import make_system
import matplotlib.cm as cm

# r1 = 'A -> B + C'
# r2 = 'C -> D'
# r3 = 'D -> E'
# r4 = 'E -> B + F'
# r5 = 'F -> G'
# r6 = 'G -> B + H'
# r7 = 'H -> I'
# reactions = [r1, r2, r3, r4, r5, r6, r7]
# mechanism = make_system(reactions)
# print(mechanism)

def kinetic_model(x, init, k1, k2, k3, k4, k5, k6, k7):
    CA,CB,CC,CD,CE,CF,CG,CH,CI = init
    dAdt = - k1*CA
    dBdt = k1*CA + k4*CE + k6*CG
    dCdt = k1*CA - k2*CC
    dDdt = k2*CC - k3*CD
    dEdt = k3*CD - k4*CE
    dFdt = k4*CE - k5*CF
    dGdt = k5*CF - k6*CG
    dHdt = k6*CG - k7*CH
    dIdt = k7*CH
    return dAdt,dBdt,dCdt,dDdt,dEdt,dFdt,dGdt,dHdt,dIdt

# Plotting the data given
species = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
initial_conditions = {
    "ic_1": np.array([4, 0, 0, 0, 0, 0, 0, 0, 0])
    }

rate_constants = np.array([1.514, 8.259, 8.359, 9.352, 7.001, 7.621 , 6.493])
    
num_exp = len(initial_conditions)
num_species = len(species)

timesteps = 30
time = np.linspace(0, 2, timesteps)
t = [0, np.max(time)]
t_eval = list(time)
STD = 0.0
noise = [np.random.normal(0, STD, size = (num_species, timesteps)) for i in range(num_exp)]
in_silico_data = {}
no_noise_data = {}

for i in range(num_exp):
    ic = initial_conditions["ic_" + str(i + 1)]
    solution = solve_ivp(kinetic_model, t, ic, t_eval = t_eval, method = "RK45", 
                          args = rate_constants)
    in_silico_data["exp_" + str(i + 1)] = np.clip(solution.y + noise[i], 0, 1e99)
    no_noise_data["exp_" + str(i + 1)] = solution.y

color_1 = cm.plasma(np.linspace(0, 1, 9))
marker = ['o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o']

# Plotting the in-silico data for visualisation
for i in range(num_exp):
    fig, ax = plt.subplots()
    # ax.set_title("Experiment " + str(i + 1))
    ax.set_ylabel("Concentrations $(M)$", fontsize = 18)
    ax.set_xlabel("Time $(h)$", fontsize = 18)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 18)
    
    for j in np.array([0, 1, -1]):
        y = in_silico_data["exp_" + str(i + 1)][j]
        ax.plot(time, y, marker[j], markersize = 4, label = species[j], color = color_1[j])
    
    ax.grid(alpha = 0.5)
    ax.legend(loc='upper right', fontsize = 15)
    
    file_path = 'Experiment_' + str(i + 1) +'.png'
    plt.savefig(file_path, dpi = 600, bbox_inches = "tight")

plt.show()