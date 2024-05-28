#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 10:38:38 2024

@author: md1621
"""

def format_matrix(matrix):
    def index_to_letter(index):
        # Convert numerical index to alphabetical character (0 is A, 1 is B, etc.)
        return chr(index + 65)

    formatted_rows = []

    for row in matrix:
        # Collect letters with negative and positive values
        left_side = []
        right_side = []

        for idx, value in enumerate(row):
            if value < 0:
                # Repeat the letter for the absolute value of the entry
                left_side.extend([index_to_letter(idx)] * abs(value))
            elif value > 0:
                # Repeat the letter for the value of the entry
                right_side.extend([index_to_letter(idx)] * value)

        # Join the collected letters with '+', and form the final string for this row
        left_part = ' + '.join(left_side)
        right_part = ' + '.join(right_side)
        
        # Combine both parts with '->' in between
        formatted_row = f"{left_part} -> {right_part}" if left_side else f" -> {right_part}"
        formatted_rows.append(formatted_row)

    return formatted_rows

# # Example usage:
# matrix = [[-1, 0, 1, 1], [0, 2, 0, -1]]
# formatted = format_matrix(matrix)
# print(formatted)
