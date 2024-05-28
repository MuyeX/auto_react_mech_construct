#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 16:49:45 2023

@author: md1621
"""

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, Parameter, report_fit
from sklearn.metrics import r2_score
import random
import progressbar
import regex

def make_system(reactions):
    unique = []
    for reaction in reactions:
        terms = regex.findall("(?i)[a-z]+\d*", reaction)
        unique.append(terms)
    unique = [item for sublist in unique for item in sublist]
    unique = sorted(set(unique), key=unique.index)
    unique = sorted(unique)  # Sort the species alphabetically
    
    firstline = ""
    lastline = ""
    for item in unique:
        firstline += f'C{item},'
        lastline += f'd{item}dt,'
    firstline = firstline[:-1]
    firstline += " = init"
    lastline = lastline[:-1]
    lastline = f'    return {lastline}'
    
    before_equals = ""
    after_equals = ""
    error_sol = ""
    for i in range(len(reactions)):
        before_equals += f"k{i+1},"
        after_equals += f"rate_const[{i}],"
        error_sol += f" k{i+1},"
    before_equals = before_equals[:-1]
    after_equals = after_equals[:-1]
    error_sol = error_sol[:-1]
    second_line = before_equals + " = " + after_equals
    
    differentials = []
    for item in unique:
        item = f'd{item}dt'
        differentials.append(f'{item} =')
    
    diff_dict = {}
    i = 1
    for differential in differentials:
        diff_dict[differential] = i
        i += 1
    
    reactdict = {}
    i = 1
    arrows = {}
    for reaction in reactions:
        matches = regex.finditer("\S+", reaction)
        sublist = []
        for match in matches:
            if match.group() == "->":
                arrows[i] = match.span()[0]
            if match.group() != "+":
                if match.group() != "->":
                    sublist.append([match.group(), match.span()])
        reactdict[f"r{i}"] = sublist
        i += 1
    
    reaction_num = 1    
    for reaction in reactdict.items():
        reaction_string_before = f"k{reaction_num}"
        arrow_digit = arrows[reaction_num]
        to_add_before = {}
        to_add_after = {}
        for pair in reaction[1]:
            for component in diff_dict.items():
                pair_to_compare = pair[0]
                component_to_compare = component[0][1:-4]
                digit = 1
                if pair_to_compare[0].isdigit():
                    digits_start = 0
                    digits_end = 1
                    i = 1
                    while pair_to_compare[i].isdigit():
                        digits_end += 1
                        i += 1
                    digit = pair_to_compare[digits_start:digits_end]
                    pair_to_compare = pair_to_compare[digits_end:]
                if pair_to_compare == component_to_compare:
                    if pair[1][0] < arrow_digit:
                        to_add_before[component[1]] = int(digit)
                        reaction_string_before = reaction_string_before + f"*C{pair_to_compare}"
                    else:
                        to_add_after[component[1]] = int(digit)
        reaction_num += 1
        for addition in to_add_before.items():
            differentials[addition[0]-1] += " "
            differentials[addition[0]-1] += "- "
            if addition[1] != 1:
                differentials[addition[0]-1] += f'{str(addition[1])}*'
            differentials[addition[0]-1] += f"{reaction_string_before}"
        for addition in to_add_after.items():
            differentials[addition[0]-1] += " "
            if differentials[addition[0]-1][-2] != '=':
                differentials[addition[0]-1] += '+ '
            if addition[1] != 1:
                differentials[addition[0]-1] += f'{str(addition[1])}*'
            differentials[addition[0]-1] += f"{reaction_string_before}"
    
    function_string = "def kinetic_model(x, init," + str(error_sol) + "):\n"
    function_string += f"    {firstline}\n"
    for differential in differentials:
        function_string += f"    {differential}\n"
    function_string += lastline
    
    return function_string


# def make_system(reactions):
#     unique = []
#     for reaction in reactions:
#         terms = regex.findall("(?i)[a-z]+\d*", reaction)
#         unique.append(terms)
#     unique = [item for sublist in unique for item in sublist]
#     unique = sorted(set(unique), key=unique.index)
#     unique = list(unique)
#     firstline = ""
#     lastline = ""
#     for item in unique:
#         firstline += f'C{item},'
#         lastline += f'd{item}dt,'
#     firstline = firstline[:-1]
#     firstline += " = init"
#     lastline = lastline [:-1]
#     lastline = f'    return {lastline}'
    
#     before_equals = ""
#     after_equals = ""
#     error_sol = ""
#     for i in range(len(reactions)):
#         before_equals += f"k{i+1},"
#         # after_equals += f"rate_const['k{i+1}'],"
#         after_equals += f"rate_const[{i}],"
#         error_sol += f" k{i+1},"
#     before_equals = before_equals[:-1]
#     after_equals = after_equals[:-1]
#     error_sol = error_sol[:-1]
#     second_line = before_equals + " = " + after_equals
    
#     differentials = []
#     for item in unique:
#         item = f'd{item}dt'
#         differentials.append(f'{item} =')
    
#     diff_dict = {}
#     i = 1
#     for differential in differentials:
#         diff_dict[differential] = i
#         i+=1
    
#     reactdict = {}
#     i =1
#     arrows = {}
#     for reaction in reactions:
#         matches = regex.finditer("\S+", reaction)
#         sublist = []
#         for match in matches:
#             if match.group() == "->":
#                 arrows[i] = match.span()[0]
#             if match.group() != "+":
#                 if match.group() != "->":
#                     sublist.append([match.group(),match.span()])
#         reactdict[f"r{i}"] = sublist
#         i+=1
    
#     reaction_num = 1    
#     for reaction in reactdict.items():
#         reaction_string_before = f"k{reaction_num}"
#         arrow_digit = arrows[reaction_num]
#         to_add_before = {}
#         to_add_after = {}
#         for pair in reaction[1]:
#             for component in diff_dict.items():
#                 pair_to_compare = pair[0]
#                 component_to_compare = component[0][1:-4]
#                 digit = 1
#                 if pair_to_compare[0].isdigit():
#                     digits_start = 0
#                     digits_end = 1
#                     i=1
#                     while pair_to_compare[i].isdigit():
#                         digits_end +=1
#                         i+=1
#                     digit = pair_to_compare[digits_start:digits_end]
#                     pair_to_compare = pair_to_compare[digits_end:]
#                 if pair_to_compare == component_to_compare:
#                     if pair[1][0] < arrow_digit:
#                         to_add_before[component[1]] = int(digit)
#                         reaction_string_before= reaction_string_before + f"*C{pair_to_compare}"
#                     else:
#                         to_add_after[component[1]] = int(digit)
#         reaction_num += 1
#         for addition in to_add_before.items():
#             differentials[addition[0]-1] += " "
#             differentials[addition[0]-1]+= "- "
#             if addition[1] != 1:
#                 differentials[addition[0]-1] += f'{str(addition[1])}*'
#             differentials[addition[0]-1] += f"{reaction_string_before}"
#         for addition in to_add_after.items():
#             differentials[addition[0]-1] += " "
#             if differentials[addition[0]-1][-2] != '=':
#                 differentials[addition[0]-1] += '+ '
#             if addition[1] != 1:
#                 differentials[addition[0]-1] += f'{str(addition[1])}*'
#             differentials[addition[0]-1] += f"{reaction_string_before}"
    
#     function_string = "def kinetic_model(x, init," + str(error_sol) + "):\n"
#     function_string += f"    {firstline}\n"
#     # function_string += f"    {second_line}\n"
#     for differential in differentials:
#         function_string += f"    {differential}\n"
#     function_string += lastline
    
#     return function_string 



