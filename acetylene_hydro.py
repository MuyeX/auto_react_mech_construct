#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 14:27:34 2024

@author: md1621
"""

import numpy as np
import matplotlib.pyplot as plt
from read_data import read_files, reverse_dict
np.random.seed(1998)

num_observable_species = 4
num_exp = 9
timesteps = 4
time = np.linspace(0, 1.5, timesteps)
species = ['A', 'B', 'C', 'D']

name_file = "exp_data_acetylene_hydro"
place_holder = read_files(name_file)
in_silico_data = reverse_dict(place_holder)


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
        y = in_silico_data["exp_" + str(i + 1)][:,j]
        ax.plot(time, y, marker[j], markersize = 4, label = species[j], color = color_1[j])
        # y = no_noise_data["exp_" + str(i + 1)][j]
        # ax.plot(time, y, label = species[j], color = color_1[j])
    
    ax.grid(alpha = 0.5)
    ax.legend(loc='upper right', fontsize = 15)
    
    file_path = 'Acetylene_Experiment_' + str(i + 1) +'.png'
    plt.savefig(file_path, dpi = 600, bbox_inches = "tight")

plt.show()