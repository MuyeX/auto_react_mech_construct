#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 17:33:44 2024

@author: md1621
"""

import numpy as np
import time
import multiprocessing
from multiprocessing import Pool, Manager

def make_matrix(rows, cols):
    starting_matrix = [[9 for _ in range(cols)] for _ in range(rows)]
    return np.array(starting_matrix)

def find_empty(matrix):
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if matrix[i][j] == 9:
                return (i, j)  # row, col
    return None

def sum_pos_neg_excluding_nine(arr):
    sum_negative = sum(value for value in arr if value < 0)
    sum_positive = sum(value for value in arr if value > 0 and value != 9)
    return sum_negative, sum_positive

def cumulative_sum(arr):
    cum_sum = []
    cum_sum.append(arr[0])
    total = 0

    # Iterate over the original array, summing as we go
    for num in arr:
        total += num  
        cum_sum.append(total) 

    return np.array(cum_sum[1:])

def is_valid(matrix, stoichiometry, intermediate, product, reactant):
    if not len(matrix[0]) == len(stoichiometry):
        raise AssertionError('Stoichiometry and matrix generated do not match.')

    for i in range(len(matrix)):
        sum_negative, sum_positive = sum_pos_neg_excluding_nine(matrix[i])
        if sum_negative < -2 or sum_positive > 2:
            return False
    
    if not any(9 in row for row in matrix):
        for i in range(len(matrix)):
            sum_negative, sum_positive = sum_pos_neg_excluding_nine(matrix[i])
            if sum_negative == 0 or sum_positive == 0:
                return False
    
    if not any(9 in row for row in matrix):
        matrix_stoich = np.sum(matrix, axis=0)
        if np.all(matrix_stoich == stoichiometry) == False:
            return False
        
    for i in range(intermediate, np.shape(matrix)[1]):
        array = matrix[:, i]
        cumulative = cumulative_sum(array)
        if np.all(cumulative >= 0) == False or np.all(cumulative == 0) == True:
            return False
        
    for i in range(product, intermediate):
        array = matrix[:, i]
        cumulative = cumulative_sum(array)
        if np.all(cumulative >= 0) == False:
            return False

    for i in range(reactant, product):
        array = matrix[:, i]
        cumulative = cumulative_sum(array)
        if np.all(cumulative >= stoichiometry[i]) == False:
            return False

    return True

def parallel_solve(matrix, stoichiometry, intermediate, product, reactant, time_budget, start, row, col, i):
    matrix[row][col] = i
    solutions = []
    count = [0]
    _solutions, _count = solve(matrix, stoichiometry, intermediate, product, reactant, time_budget, solutions, count, start)
    return _solutions, _count

def solve(matrix, stoichiometry, intermediate, product, reactant, time_budget, solutions=None, count=[0], start=None):
    if start is None:
        start = time.time()
    
    if solutions is None:
        solutions = []

    find = find_empty(matrix)
    
    if not find:
        if is_valid(matrix, stoichiometry, intermediate, product, reactant):
            solutions.append(np.copy(matrix))
        return solutions, count[0]

    row, col = find
    
    for i in range(-2, 3):  # From -2 to 2 inclusive.
        matrix[row][col] = i
        end = time.time()
        
        if end - start > time_budget:
            print('Time budget reached!')
            return solutions, count[0]
        
        if is_valid(matrix, stoichiometry, intermediate, product, reactant):
            count[0] += 1  # Increment the counter when is_valid is called.
            _solutions, _ = solve(matrix, stoichiometry, intermediate, product, reactant, time_budget, solutions, count, start)
        
        matrix[row][col] = 9  # Backtrack.
    
    return solutions, count[0]

if __name__ == '__main__':
    # Example usage:
    hi = time.time()
    elementary_reactions = 4
    number_species = 6
    matrix = make_matrix(elementary_reactions, number_species)
    stoichiometry = [-1, 2, 1, 0, 0, 0]
    intermediate = 3
    product = 1
    reactant = 0
    time_budget = 60

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


    for solution in solutions:
        print(solution)
        print('----------------------------')

    print('Number of solutions found: ', len(solutions))
    print('Total number of possible matrices: ', 5**(elementary_reactions * number_species))
    print('Number of matrices checked: ', count)
    print('Percentage of space checked: ', (count * 100) / 5**(elementary_reactions * number_species), '%')
    bye = time.time()
    print('Time spent: ', bye - hi, 's')
        
