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
from parameter_estimation import adjust_ic_length, sse, TimeoutException, timeout_handler, Opt_Rout, evaluate
import logging

# Configure logging to write to a file in append mode
# name_file = "output_fruc_to_hmf.log"
# name_file = "output_hypoth.log"
# name_file = "output_aldol_condensation.log"
# name_file = "output_fruc_to_hmf2.log"
name_file = "aldol_test.log"
logging.basicConfig(filename = name_file, level = logging.INFO, \
                    format = '%(message)s', filemode = 'a')


def evaluate_solution_parallel(solution):
    sol, idx, len_sol = solution
    model_pred, opt_param, nll, aic = evaluate(sol)
    print(idx + 1, '/', len_sol)
    print(sol)
    print(aic)
    print('\n')
    return {'aic': aic, 'solution': sol}


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

        # solutions, count = solve_dfs(matrix.copy(), stoichiometry, intermediate, product, reactant, time_budget)

        for i in range(len(solutions)):
            solutions[i] = (solutions[i], i, len(solutions))

        all_AIC = []
        # print(solutions)
        # Parallelize the evaluation of solutions
        with Pool(processes=multiprocessing.cpu_count()) as pool:
            all_AIC = pool.map(evaluate_solution_parallel, solutions)

        for i, result in enumerate(all_AIC):
            print(i + 1, '/', len(solutions))
            # print(result['solution'])
            print(result['aic'])
            print('\n')
            
        if len(solutions) == 0:
            all_AIC.append({'aic': 1e99, 'solution': 1e99})

        for i in range(len(all_AIC)):
            all_AIC[i]['position'] = i
            
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
        # print('Current best mechanism: \n', min_AIC_solution)
        print('Current AIC value:', min_AIC_value)
        print('Previous best mechanism: \n', last_mech)
        print('Previous AIC value:', last_AIC)

        logging.info('\n' + f'ITERATION NUMBER: {iteration_counter}')
        logging.info("#" * 80 + "\n")
        logging.info(f'Current best mechanism: \n {min_AIC_solution}')
        logging.info(f'Current AIC value: {min_AIC_value}')
        logging.info(f'Previous best mechanism: \n {last_mech}')
        logging.info(f'Previous AIC value: {last_AIC}')
        logging.info( "#" * 80 + "\n")
        
        if len(min_AIC_entries) > 1:
            print('All possible solutions:', min_AIC_entries)
        
    return opt_solution



if __name__ == '__main__':
    # Example usage:
    elementary_reactions = 1
    number_species = 4
    stoichiometry = [-1, -1, 1, 1]
    intermediate = 4
    product = 2
    reactant = 0
    time_budget = 10
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

    logging.info("\n" + "#"*80 + "\n")
    logging.info("Solution found!")
    logging.info(f"Optimal reaction matrix: \n {found_mechanism['reaction_matrix']}")
    logging.info(f"Optimal reaction chain: {found_mechanism['reaction_chain']}")
    logging.info(f"Optimal reaction parameters: {found_mechanism['opt_param']}")
    logging.info(f"Optimal NLL: {found_mechanism['nll']}")
    logging.info(f"Optimal AIC: {found_mechanism['AIC']}")
    logging.info("\n" + "#"*80 + "\n")

        
