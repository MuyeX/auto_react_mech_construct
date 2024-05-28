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
from parallel_backtracking import make_matrix, find_empty, sum_pos_neg_excluding_nine, cumulative_sum, is_valid, parallel_solve, solve
from matrix_to_reaction_string import format_matrix
from ODE_generator import make_system
from parameter_estimation import sse, callback, Opt_Rout, evaluate


def bob_the_mechanism_builder(elementary_reactions, number_species, stoichiometry, intermediate, product, reactant, time_budget):
    
    iteration_counter = 0

    # Initialize flag variables
    last_AIC = 1e99
    min_AIC_value = 1e99
    opt_solution = {}

    while min_AIC_value <= last_AIC:
        # Keep track of AIC values and last optimal solution for breaking the loop 
        last_AIC = min_AIC_value
        
        if iteration_counter > 0:
            elementary_reactions += 1
            number_species += 1
            stoichiometry.append(0)
            model_pred, nll, aic = evaluate(min_AIC_solution)
            opt_solution["model_predictions"] = model_pred
            opt_solution["reaction_chain"] = format_matrix(min_AIC_solution)
            opt_solution["nll"] = nll
            opt_solution["AIC"] = aic

        matrix = make_matrix(elementary_reactions, number_species)
    
        start = time.time()
        find = find_empty(matrix)
        row, col = find
    
        tasks = [(matrix.copy(), stoichiometry, intermediate, product, reactant, time_budget, start, row, col, i) for i in range(-2, 3)]
    
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
            #TODO: if no solutions are found, output the best previous solution
            model_pred, nll, aic = evaluate(solution)
            
            # Store the AIC value, solution, and position in a dictionary
            all_AIC.append({'aic': aic, 'solution': solution, 'position': i})
            
            print(aic)
            print('----------------------------')
            
        # Find the dictionary with the smallest AIC value
        min_AIC_entry = min(all_AIC, key=lambda x: x['aic'])
        min_AIC_value = min_AIC_entry['aic']
        min_AIC_position = min_AIC_entry['position']
        min_AIC_solution = min_AIC_entry['solution']
        
        iteration_counter += 1 
        print('ITERATION NUMBER:', iteration_counter)
    
    # Print important information of the chosen solution
    print("#"*50, "\n")
    print("Solution found!")
    print(f"Optimal reaction chain: {opt_solution['reaction_chain']}")
    print(f"Optimal NLL: {opt_solution['nll']}")
    print(f"Optimal AIC: {opt_solution['AIC']}")



if __name__ == '__main__':
    # Example usage:
    elementary_reactions = 2
    number_species = 4
    stoichiometry = [-1, 2, 1, 0]
    intermediate = 3
    product = 1
    reactant = 0
    time_budget = 1200
    bob_the_mechanism_builder(elementary_reactions, number_species, stoichiometry, intermediate, product, reactant, time_budget)
            
        
        
        
        
        
        
        
    

