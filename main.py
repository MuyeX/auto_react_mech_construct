#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 11:39:13 2024

@author: md1621
"""

import numpy as np
np.random.seed(1998)
import time
import multiprocessing
from multiprocessing import Pool, Manager
from parallel_backtracking import make_matrix, find_empty, sum_pos_neg_excluding_nine, cumulative_sum, is_valid, parallel_solve, solve, generate_initial_matrices
from matrix_to_reaction_string import format_matrix
from ODE_generator import make_system
from parameter_estimation import sse, callback, timeout_handler, Opt_Rout, evaluate
import logging

# Configure logging to write to a file in append mode
# name_file = "output_fruc_to_hmf.log"
# name_file = "output_hypoth.log"
name_file = "testing.log"
# name_file = "output_aldol_condensation.log"
logging.basicConfig(filename = name_file, level = logging.INFO, \
                    format = '%(message)s', filemode = 'a')


def bob_the_mechanism_builder(elementary_reactions, number_species, stoichiometry, intermediate, product, reactant, time_budget):
    
    iteration_counter = 0

    # Initialize flag variables
    last_AIC = 1e99
    min_AIC_value = 1e99
    min_AIC_solution = ['Nothing']
    opt_solution = {}

    while min_AIC_value <= last_AIC:
        # Keep track of AIC values and last optimal solution for breaking the loop 
        last_AIC = min_AIC_value
        last_mech = min_AIC_solution
        
        if iteration_counter > 0:
            elementary_reactions += 1
            number_species += 1
            stoichiometry.append(0)
            model_pred, opt_param, nll, aic = evaluate(min_AIC_solution)
            opt_solution["model_predictions"] = model_pred
            opt_solution["opt_param"] = opt_param
            opt_solution["reaction_chain"] = format_matrix(min_AIC_solution)
            opt_solution["reaction_matrix"] = min_AIC_solution
            opt_solution["nll"] = nll
            opt_solution["AIC"] = aic

        matrix = make_matrix(elementary_reactions, number_species)
    
        start = time.time()
        find = find_empty(matrix)
        row, col = find
    
        # tasks = [(matrix.copy(), stoichiometry, intermediate, product, reactant, time_budget, start, row, col, i) for i in range(-2, 3)]

        tasks = []
        initial_matrices = generate_initial_matrices(elementary_reactions, number_species)

        for matrix in initial_matrices:
            find = find_empty(matrix)
            if find:
                row, col = find
                num_tasks = range(-2, 3)
                for i in num_tasks:
                    tasks.append((matrix.copy(), stoichiometry, intermediate, product, reactant, time_budget, start, row, col, i))

    
        with Pool(processes=multiprocessing.cpu_count()) as pool:
            results = pool.starmap(parallel_solve, tasks)
    
        solutions = []
        count = 0
        for result in results:
            _solutions, _count = result
            solutions.extend(_solutions)
            count += _count
    
        all_AIC = []
                    
        for i, solution in enumerate(solutions):
            model_pred, opt_param, nll, aic = evaluate(solution)
            
            # Store the AIC value, solution, and position in a dictionary
            all_AIC.append({'aic': aic, 'solution': solution, 'position': i})
            
            print(i + 1, '/', len(solutions))
            print(solution)
            print(aic)
            print('\n')
            
        if len(solutions) == 0:
            all_AIC.append({'aic': 1e99, 'solution': 1e99, 'position': 1e99})
            
        # Find the dictionary with the smallest AIC value
        min_aic_value = min(x['aic'] for x in all_AIC)

        # Collect all entries with the minimum AIC value
        min_AIC_entries = [entry for entry in all_AIC if entry['aic'] == min_aic_value]

        # Choose one entry from the minimum AIC entries to update the loop variables
        min_AIC_entry = min_AIC_entries[0]
        min_AIC_value = min_AIC_entry['aic']
        min_AIC_position = min_AIC_entry['position']
        min_AIC_solution = min_AIC_entry['solution']
        
        # min_AIC_entry = min(all_AIC, key=lambda x: x['aic'])
        # min_AIC_value = min_AIC_entry['aic']
        # min_AIC_position = min_AIC_entry['position']
        # min_AIC_solution = min_AIC_entry['solution']
        
        iteration_counter += 1 
        print('ITERATION NUMBER:', iteration_counter)
        print('Current best mechanism: \n', min_AIC_solution)
        print('Current AIC value:', min_AIC_value)
        print('Previous best mechanism: \n', last_mech)
        print('Previous AIC value:', last_AIC)
        
        logging.info(f'ITERATION NUMBER: {iteration_counter}')
        logging.info(f'Current best mechanism: \n {min_AIC_solution}')
        logging.info(f'Current AIC value: {min_AIC_value}')
        logging.info(f'Previous best mechanism: \n {last_mech}')
        logging.info(f'Previous AIC value: {last_AIC}')
        
        if len(min_AIC_entries) > 1:
            print('All possible solutions:', min_AIC_entries)
        
    return opt_solution



if __name__ == '__main__':
    # Example usage:
    elementary_reactions = 2
    number_species = 3
    stoichiometry = [-4, 1, 1]
    intermediate = 3
    product = 1
    reactant = 0
    time_budget = 60 * 60 * 2
    found_mechanism = bob_the_mechanism_builder(elementary_reactions, \
                                                number_species, stoichiometry, \
                                                intermediate, product, reactant, \
                                                time_budget)
        
    print("\n", "#"*80, "\n")
    print("Solution found!")
    print(f"Optimal reaction matrix: \n {found_mechanism['reaction_matrix']}")
    print(f"Optimal reaction chain: {found_mechanism['reaction_chain']}")
    print(f"Optimal reaction parameters: {found_mechanism['opt_param']}")
    print(f"Optimal NLL: {found_mechanism['nll']}")
    print(f"Optimal AIC: {found_mechanism['AIC']}")    

            
        
        
        
        
        
        
        
    

