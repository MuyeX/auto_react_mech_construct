#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 16:36:35 2024

@author: md1621
"""

import itertools
import numpy as np
import time

start = time.time()

def generate_matrices_iteratively(rows, cols):
    # Define the range of values each entry can take
    value_range = range(-2, 3)
    
    # Calculate the total number of elements in the matrix
    total_elements = rows * cols
    
    # Use itertools.product to generate all possible combinations of values
    # for the matrix entries
    all_combinations = itertools.product(value_range, repeat=total_elements)
    
    # Yield each combination as a matrix one by one
    for combination in all_combinations:
        yield np.array(combination).reshape(rows, cols)

# Example usage
rows, cols = 2, 4  # For demonstration purposes; large values will result in huge output
matrix_generator = generate_matrices_iteratively(rows, cols)


def find_good_matrices(rows, cols, matrix_generator, stoichiometry, position_intermediates):
    
    total_elements = rows * cols
    number_matrices = 5**total_elements
    good_matrices = []
    
    for i in range(number_matrices):
        
        current_matrix = next(matrix_generator)
        
        if not np.shape(current_matrix[0]) == np.shape(stoichiometry):
            raise AssertionError('Stoichiometry and matrix generated do not match.')
        
        matrix_stoich = np.sum(current_matrix, axis = 0)
        negative_sums = np.sum(current_matrix * (current_matrix < 0), axis=1)
        positive_sums = np.sum(current_matrix * (current_matrix > 0), axis=1)
        intermediates = current_matrix[:, position_intermediates:]
        
        if np.all(matrix_stoich == stoichiometry) == True:
            #the mechanism needs to add up to the total stoichiometry
            if np.all(negative_sums >= -2) == True:
                #each elementary step cannot have more than 2 molecules reacting
                if np.all(negative_sums <= -1) == True:
                    #each elementary step must have at least 1 molecule reacting
                    if np.all(positive_sums <= 2) == True:
                        #each elementary step cannot generate more than 2 molecules
                        if np.all(positive_sums >= 1) == True:
                            #each elementary step must generate at least 1 molecule
                            inter_row, inter_col = np.shape(intermediates)
                            sum_inter = np.zeros((inter_row, inter_col))
                            sum_inter[0] = intermediates[0]
                            
                            for i in range(inter_row - 1):
                                sum_inter[i + 1] = sum_inter[i] + intermediates[i + 1]
                            
                            if np.all(sum_inter >= 0) == True:
                                #before reacting an intermediate, the intermediate must be generated
                                if not np.all(sum_inter == 0) == True:
                                    #before reacting an intermediate, the intermediate must be generated
                                    good_matrices.append(current_matrix)
    return good_matrices

stoichiometry = np.array([-1, 2, 1, 0])
position_intermediates = 3
tests = find_good_matrices(rows, cols, matrix_generator, stoichiometry, position_intermediates)

end = time.time()
# print('Time for execution:', end - start, 's')

for test in tests:
    print(test, '\n')
        