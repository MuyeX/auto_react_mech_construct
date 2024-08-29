#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 13:18:02 2023

@author: md1621
"""
import numpy as np
np.random.seed(1998)
from scipy.integrate import solve_ivp
import pandas as pd
from scipy.optimize import minimize
import os
from ODE_generator import make_system
from metrics import NLL_mechanism, Akaike
from matrix_to_reaction_string import format_matrix
from read_data import *
import signal
from read_data import read_files, reverse_dict


"##############################################################################"
"############################ Optimise Rate Model #############################"
"##############################################################################"

name_file = "exp_data_fruc_HMF"
# name_file = "exp_data_hypoth"
# name_file = "exp_data_aldol_condensation"

place_holder = read_files(name_file)

in_silico_data = reverse_dict(place_holder)

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


def sse(kinetic_model, params, num_species):
    num_observable_species = 3
    num_exp = len(initial_conditions)
    total = np.zeros((num_exp, 1))

    for i in range(num_exp):
        aa = initial_conditions["ic_" + str(i+1)]
        ic = np.zeros(num_species)
        
        for j in range(num_observable_species):
            ic[j] = aa[j]
        
        sorted_data = {ii: in_silico_data[ii] for ii in sorted(in_silico_data.keys(), key=lambda x: int(x.split('_')[1]))}
        observations = sorted_data["exp_" + str(i + 1)]
        solution = solve_ivp(kinetic_model, t, ic, t_eval = t_eval, method = "RK45",
                             args = params)
        model_response = np.array([solution.y.T[:, k] for k in range(num_observable_species)]).T

        SSE = (observations - model_response)**2
        total[i] = np.sum(SSE)

    return np.sum(total)


def callback(xk):
    # Print out the current solution
    pass
    # print(f"Current solution: {xk}")


# def Opt_Rout(multistart, number_parameters, x0, lower_bound, upper_bound, to_opt, num_species):
#     localsol = np.empty([multistart, number_parameters])
#     localval = np.empty([multistart, 1])
#     boundss = tuple([(lower_bound, upper_bound) for i in range(number_parameters)])
    
#     for i in range(multistart):
#         res = minimize(to_opt, x0, method='L-BFGS-B', bounds=boundss, callback=callback, args=(num_species,))
#         localsol[i] = res.x
#         localval[i] = res.fun

#     minindex = np.argmin(localval)
#     opt_val = localval[minindex]
#     opt_param = localsol[minindex]
    
#     return opt_val, opt_param

def timeout_handler(signum, frame):
    raise TimeoutError

def Opt_Rout(multistart, number_parameters, x0, lower_bound, upper_bound, to_opt, num_species):
    localsol = np.empty([multistart, number_parameters])
    localval = np.empty([multistart, 1])
    boundss = tuple([(lower_bound, upper_bound) for i in range(number_parameters)])
    
    signal.signal(signal.SIGALRM, timeout_handler)

    try:
        for i in range(multistart):
            signal.alarm(1)  # Set the timeout to x seconds
            res = minimize(to_opt, x0, method='L-BFGS-B', bounds=boundss, args=(num_species,))
            localsol[i] = res.x
            localval[i] = res.fun
            signal.alarm(0)  # Disable the alarm if optimization completes in time
    
    except TimeoutError:
        # print("Optimization took longer than 1 seconds.")
        opt_param = np.full(number_parameters, upper_bound)
        opt_val = 1e99
        return opt_val, opt_param

    minindex = np.argmin(localval)
    opt_val = localval[minindex]
    opt_param = localsol[minindex]
    
    return opt_val, opt_param


def evaluate(reaction_matrix):
    
    num_observable_species = 3
    reactions = format_matrix(reaction_matrix)
    mechanism = make_system(reactions)
    # The function executed below is called kinetic_model
    exec(mechanism, globals())

    number_parameters, num_species = np.shape(reaction_matrix)
    multistart = 2
    lower_bound = 0.0001
    upper_bound = 10

    solution = np.array([np.random.uniform(lower_bound, upper_bound) for i in range(number_parameters)])
    # solution = np.array([0.1 for i in range(number_parameters)])
    
    # aaaa = np.array([[-1, 1, 0, 0, 0, 1],
    #   [0, 1, 0, 0 , 1 , -1],
    #   [0, 0, 0, 1 , -1, 0 ],
    #   [0, 1, 1, -1, 0 , 0 ]])
    
    # if np.array_equal(aaaa, reaction_matrix):
    #     lower_bound = 0.0001
    #     upper_bound = 10
    #     solution = np.array([1.514, 5.259, 9.352, 2.359])
    #     print('SUIIIIIIIIIIIIIII')

    opt_val, opt_param = Opt_Rout(multistart, number_parameters, solution, lower_bound, 
        upper_bound, lambda params, ns: sse(kinetic_model, params, ns), num_species)

    # print('MSE = ', opt_val)
    # print('Optimal parameters = ', opt_param)

    model_predictions = {}
    for i in range(num_exp):
        aa = initial_conditions["ic_" + str(i+1)]
        ic = np.zeros(num_species)
        
        for j in range(num_observable_species):
            ic[j] = aa[j]
                
        # ic[0], ic[1], ic[2], ic[3] = aa[0], aa[1], aa[2], aa[3]
        
        try:
            # Set the signal handler and a 1-second alarm
            signal.signal(signal.SIGALRM, timeout_handler)
            signal.alarm(1)  # Set the timeout to x second
            
            solution = solve_ivp(kinetic_model, t, ic, t_eval=t_eval, method="RK45", args=opt_param)

            # Disable the alarm after successful completion
            signal.alarm(0)
            
            model_predictions["exp_" + str(i + 1)] = np.array([solution.y.T[:, k] for k in range(num_observable_species)]).T
            
        except TimeoutError:
            # print(f"Integration for experiment {i + 1} took longer than 1 second.")
            model_predictions["exp_" + str(i + 1)] = np.full((len(t_eval), num_species), 1e99)
    
    sorted_data = {ii: in_silico_data[ii] for ii in sorted(in_silico_data.keys(), key=lambda x: int(x.split('_')[1]))}
    exp_data = np.vstack(list(sorted_data.values()))
    model_pred = np.vstack(list(model_predictions.values()))

    nll = NLL_mechanism(exp_data, model_pred)
    AIC = Akaike(nll, opt_param)
    # print('NLL value: ', nll)
    # print('AIC value: ', AIC)
    
    return model_predictions, opt_param, nll, AIC
