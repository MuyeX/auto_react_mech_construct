#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 15:52:32 2024

@author: md1621
"""

import numpy as np
import time


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
    
    if not any(9 in row for row in matrix) == True:
        for i in range(len(matrix)):
            sum_negative, sum_positive = sum_pos_neg_excluding_nine(matrix[i])
            if sum_negative == 0 or sum_positive == 0:
                return False
    
    if not any(9 in row for row in matrix) == True:
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


def solve(matrix, stoichiometry, intermediate, product, reactant, time_budget, solutions=None, count=[0]):
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
            _solutions, _ = solve(matrix, stoichiometry, intermediate, product, reactant, time_budget, solutions, count)
        
        matrix[row][col] = 9  # Backtrack.
    
    return solutions, count[0]


# Example usage:
matrix = make_matrix(3, 5)
stoichiometry = [-1, 2, 1, 0, 0]
intermediate = 3
product = 1
reactant = 0
time_budget = 300
solutions, count = solve(matrix, stoichiometry, intermediate, product, reactant, time_budget)
print('Number of solutions found:', len(solutions))
print('Number of matrices checked:', count)

for solution in solutions:
    print(solution)
    print('----------------------------')