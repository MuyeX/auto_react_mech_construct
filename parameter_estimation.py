#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 13:18:02 2023

@author: md1621
"""

"##############################################################################"
"######################## Importing important packages ########################"
"##############################################################################"

import numpy as np
from scipy.integrate import solve_ivp
import pandas as pd
from scipy.optimize import minimize
from beecolpy import abc
import os
from ODE_generator import make_system
import matplotlib.pyplot as plt
import matplotlib.cm as cm
np.random.seed(1998)


"##############################################################################"
"############################ Optimise Rate Model #############################"
"##############################################################################"

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

# This takes the first column from each entry of the dictionary and puts it into another dictionary
initial_conditions = {}
for key, value in in_silico_data.items():
    aa = "ic_" + key[-1]
    initial_conditions[aa] = value[0]

num_exp = len(initial_conditions)
timesteps = 30
time = np.linspace(0, 2, timesteps)
t = [0, np.max(time)]
t_eval = list(time)


def sse(params, num_species):
    num_exp = len(initial_conditions)
    total = np.zeros((num_exp, 1))

    for i in range(num_exp):
        aa = initial_conditions["ic_" + str(i+1)]
        ic = np.zeros(num_species)
        ic[0], ic[1], ic[-1] = aa[0], aa[1], aa[2]
        observations = in_silico_data["exp_" + str(i + 1)]
        solution = solve_ivp(kinetic_model, t, ic, t_eval = t_eval, method = "RK45",
                             args = params)
        model_response = np.array([solution.y.T[:, 0], solution.y.T[:, 1], solution.y.T[:, -1]]).T

        SSE = (observations - model_response)**2
        total[i] = np.sum(SSE)

    return np.sum(total)


def callback(xk):
    # Print out the current solution
    print(f"Current solution: {xk}")


def Opt_Rout(multistart, number_parameters, x0, lower_bound, upper_bound, to_opt, num_species):
    localsol = np.empty([multistart, number_parameters])
    localval = np.empty([multistart, 1])
    boundss = tuple([(lower_bound, upper_bound) for i in range(number_parameters)])
    
    for i in range(multistart):
        res = minimize(to_opt, x0, method = 'L-BFGS-B', \
                       bounds = boundss, callback = callback, args = num_species)
        localsol[i] = res.x
        localval[i] = res.fun

    minindex = np.argmin(localval)
    opt_val = localval[minindex]
    opt_param = localsol[minindex]
    
    return opt_val, opt_param

# r1 = 'A -> B + C'
# r2 = 'C -> D'
# r3 = 'D -> E'
# r4 = 'E -> B + F'
# r5 = 'F -> G'
# r6 = 'G -> B + H'
# r7 = 'H -> I'

r1 = 'A -> B + C'
r2 = 'C -> B + D'
r3 = 'D -> B + E'
r4 = 'E -> I'
reactions = [r1, r2, r3, r4]
mechanism = make_system(reactions)
# The function executed below is called kinetic_model
exec(mechanism)

# Create an empty set to store the unique letters
unique_letters = set()

# Loop through each reaction and add the unique letters to the set
for reaction in reactions:
    for letter in reaction:
        if letter.isalpha():
            unique_letters.add(letter)

num_species = len(unique_letters)
multistart = 10
number_parameters = len(reactions)
lower_bound = 0.0001
upper_bound = 10

# abc_obj = abc(sse, [(lower_bound, upper_bound) for i in range(number_parameters)])
# abc_obj.fit() 

# solution = abc_obj.get_solution()
# print('Initial guess = ', solution)

solution = np.array([0.1 for i in range(number_parameters)])

opt_val, opt_param = Opt_Rout(multistart, number_parameters, solution, lower_bound, \
    upper_bound, sse, num_species)

print('MSE = ', opt_val)
print('Optimal parameters = ', opt_param)


def NLL_mechanism(exp_data, model_pred):
    number_species = np.shape(exp_data)[1]
    number_datapoints = np.shape(exp_data)[0]
    output = np.zeros(number_species)
    mse = np.zeros(number_species)
    variance = np.zeros(number_species)
    
    for i in range(number_species):
        a = ((exp_data[:, i] - model_pred[:, i])**2)
        mse[i] = np.sum(a)
        variance[i] = mse[i] / (number_datapoints)
    
    for i in range(number_species):
        likelihood = ((exp_data[:, i] - model_pred[:, i])**2) / (2 * (variance[i])) \
            - np.log(1 / (np.sqrt(2 * np.pi * (variance[i]))))
        output[i] = np.sum(likelihood)
    
    return np.sum(output)

def Akaike(nll, params):
    num_param = len(params)
    aic = 2 * nll + 2 * num_param
    
    return aic

model_predictions = {}
for i in range(num_exp):
    aa = initial_conditions["ic_" + str(i+1)]
    ic = np.zeros(num_species)
    ic[0], ic[1], ic[-1] = aa[0], aa[1], aa[2]
    solution = solve_ivp(kinetic_model, t, ic, t_eval = t_eval, method = "RK45",
                         args = opt_param)
    model_predictions["exp_" + str(i + 1)] = np.array([solution.y.T[:, 0], 
                                                       solution.y.T[:, 1], 
                                                       solution.y.T[:, -1]]).T

exp_data = np.vstack(list(in_silico_data.values()))
model_pred = np.vstack(list(model_predictions.values()))

nll = NLL_mechanism(exp_data, model_pred)
AIC = Akaike(nll, opt_param)
print('NLL value: ', nll)
print('AIC value: ', AIC)


species = list(sorted(unique_letters))
color_1 = cm.plasma(np.linspace(0, 1, 3))
marker = ['o' for i in range(3)]

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
        y = in_silico_data["exp_" + str(i + 1)][:, j]
        yy = model_predictions["exp_" + str(i + 1)][:, j]
        ax.plot(time, y, marker[j], markersize = 4, label = species[j], color = color_1[j])
        ax.plot(time, yy, '-', color = color_1[j])
    
    ax.grid(alpha = 0.5)
    ax.legend(loc='upper right', fontsize = 15)
    
plt.show()