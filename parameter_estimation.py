#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 13:40:46 2024

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
from functools import partial


"##############################################################################"
"############################ Optimise Rate Model #############################"
"##############################################################################"

# name_file = "exp_data_fruc_HMF"
# name_file = "exp_data_hypoth"
# name_file = "exp_data_aldol_condensation"
name_file = "exp_data_fruc_HMF2"

place_holder = read_files(name_file)

in_silico_data = reverse_dict(place_holder)

# This takes the first column from each entry of the dictionary and puts it into another dictionary
# initial_conditions = {}
# for key, value in in_silico_data.items():
#     aa = "ic_" + key[-1]
#     initial_conditions[aa] = value[0]
    
if name_file == "exp_data_fruc_HMF":
    initial_conditions = {
        "ic_1": np.array([4, 0, 0, 0, 0, 0]),
        "ic_2": np.array([6, 2, 1, 0, 0, 0]),
        "ic_3": np.array([4, 2, 0, 0, 0, 0]),
        "ic_4": np.array([6, 0, 0, 0, 0, 0]),
        "ic_5": np.array([6, 2, 0, 0, 0, 0])
        }

if name_file == "exp_data_aldol_condensation":
    initial_conditions = {
        "ic_1": np.array([5 , 10, 0, 0, 0, 0]),
        "ic_2": np.array([5 , 5 , 2, 0, 0, 0]),
        "ic_3": np.array([5 , 10, 0, 2, 0, 0]),
        "ic_4": np.array([10, 10, 0, 2, 0, 0]),
        "ic_5": np.array([10, 10, 2, 2, 0, 0])
        }
    
if name_file == "exp_data_hypoth":
    initial_conditions = {
        "ic_1": np.array([10, 0, 2, 0, 0]),
        "ic_2": np.array([10, 2, 0, 0, 0]),
        "ic_3": np.array([10, 2, 2, 0, 0]),
        "ic_4": np.array([5 , 0, 0, 0, 0]),
        "ic_5": np.array([10, 0, 0, 0, 0])
        }
    
if name_file == "exp_data_fruc_HMF2":
    initial_conditions = {
        "ic_1": np.array([4, 0, 0]),
        "ic_2": np.array([6, 2, 1]),
        "ic_3": np.array([4, 2, 0]),
        "ic_4": np.array([6, 0, 0]),
        "ic_5": np.array([6, 2, 0])
        }


num_exp = len(initial_conditions)
timesteps = 30
time = np.linspace(0, 90, timesteps)
t = [0, np.max(time)]
t_eval = list(time)


def adjust_ic_length(ic, num_species):
    """
    Adjusts the length of the initial condition array to match the number of species.
    If too long, truncate from the end. If too short, append zeros at the end.

    Parameters:
    - ic (array-like): Initial condition array.
    - num_species (int): Desired length of the array.

    Returns:
    - adjusted_ic (np.array): Adjusted initial condition array.
    """
    ic = np.array(ic)
    if len(ic) > num_species:
        return ic[:num_species]
    elif len(ic) < num_species:
        return np.pad(ic, (0, num_species - len(ic)), 'constant')
    return ic

def sse(kinetic_model, params, num_species):
    """
    Calculates the sum of squared errors (SSE) for fitting an ODE system to experimental data.

    Parameters:
    - kinetic_model (function): The ODE system to be fitted.
    - params (list or np.array): Initial guess for the kinetic parameters.
    - num_species (int): Number of species in the ODE model.

    Returns:
    - sse (float): The sum of squared errors for the given model and parameters.
    """

    def simulate_experiment(ic, params):
        """Simulates the ODE system for given initial conditions and parameters."""
        solution = solve_ivp(
            lambda t, y: kinetic_model(t, y, *params),
            [time[0], time[-1]],
            ic,
            t_eval=time,
            method="RK45"
        )
        return solution.y

    # Initialize SSE
    sse = 0.0
    num_observable_species = 3

    # Iterate over all experiments
    for i, (exp, ic) in enumerate(initial_conditions.items()):
        # Adjust the initial conditions to match the number of species
        adjusted_ic = adjust_ic_length(ic, num_species)

        # Simulate the ODE system
        simulated_data = simulate_experiment(adjusted_ic, params)

        # Only consider the first `num_species` (observable species)
        simulated_observable = simulated_data[:num_observable_species, :]

        # Get observed data for the current experiment
        observed_data = in_silico_data[f"exp_{i+1}"]
        
        # Transpose observed data if necessary to match shapes
        if observed_data.shape != simulated_observable.shape:
            observed_data = observed_data.T

        # Compute squared errors
        squared_errors = (simulated_observable - observed_data) ** 2

        # Accumulate SSE
        sse += np.sum(squared_errors)

    return sse


class TimeoutException(Exception):
    """Custom exception for handling timeouts."""
    pass

def timeout_handler(signum, frame):
    """Signal handler for the timeout."""
    raise TimeoutException()


def Opt_Rout(multistart, number_parameters, x0, lower_bound, upper_bound, to_opt, num_species):
    """
    Estimates the parameters of a given kinetic_model using optimization.

    Parameters:
    - multistart (int): Number of optimization runs with different initial guesses.
    - number_parameters (int): Number of tunable parameters in the kinetic model.
    - x0 (array-like): Initial guess for the parameters.
    - lower_bound (array-like): Lower bounds for the parameters.
    - upper_bound (array-like): Upper bounds for the parameters.
    - to_opt (function): Objective function to minimize (e.g., SSE).
    - kinetic_model (function): The ODE system to be fitted.
    - num_species (int): Number of species in the kinetic model.

    Returns:
    - best_params (np.array): Best-fit parameters from optimization.
    - best_sse (float): SSE corresponding to the best-fit parameters.
    """
    best_sse = np.inf
    best_params = None

    bounds = [(lower_bound, upper_bound) for _ in range(number_parameters)]

    # Partially apply to_opt with fixed kinetic_model and num_species arguments
    def partial_to_opt(params):
        return to_opt(kinetic_model, params, num_species)

    for _ in range(multistart):
        # Generate random initial guess within bounds
        random_x0 = np.random.uniform(lower_bound, upper_bound, number_parameters)

        try:
            # Set up the timeout
            signal.signal(signal.SIGALRM, timeout_handler)
            signal.alarm(10)  # Set the timeout to 1 second

            # Run optimization
            result = minimize(
                partial_to_opt,  # Pass the custom wrapper function
                random_x0,
                method='L-BFGS-B',
                bounds=bounds
            )

            # Cancel the alarm after successful optimization
            signal.alarm(0)

            # Update best parameters if current run is better
            if result.fun < best_sse:
                best_sse = result.fun
                best_params = result.x

        except TimeoutException:
            # Optimization timed out
            print("Optimization timed out. Skipping this run...")
            best_params = random_x0  # Random parameters
            best_sse = np.inf  # Infinite SSE

    return best_params, best_sse



def evaluate(reaction_matrix):
    
    num_observable_species = 3
    reactions = format_matrix(reaction_matrix)
    mechanism = make_system(reactions)
    # The function executed below is called kinetic_model
    exec(mechanism, globals())

    number_parameters, num_species = np.shape(reaction_matrix)
    multistart = 2
    lower_bound = 0.0001
    upper_bound = 3

    # solution = np.array([np.random.uniform(lower_bound, upper_bound) for i in range(number_parameters)])
    solution = np.array([1 for i in range(number_parameters)])
    
    # aaaa = np.array([[-1,  0,  0,  1,  0,  0],[ 0,  1,  0,  -1,  1,  0],[ 0,  1,  0,  0,  -1,  1],[ 0,  1,  1, 0,  0,  -1]])
    
    # if np.array_equal(aaaa, reaction_matrix):
    #     lower_bound = 1
    #     upper_bound = 10
    #     solution = np.array([0.10159522, 0.35146577, 1.71775928, 0.04800828])
    #     print('SUIIIIIIIIIIIIIII')

    # opt_val, opt_param = Opt_Rout(multistart, number_parameters, solution, lower_bound, 
    #     upper_bound, lambda params, ns: sse(kinetic_model, params, ns), num_species)
    
    opt_param, opt_val = Opt_Rout(multistart, number_parameters, solution, lower_bound, 
        upper_bound, sse, num_species)


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
