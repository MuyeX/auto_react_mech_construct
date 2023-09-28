#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 09:37:56 2023

@author: md1621
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from ODE_generator import make_system
import matplotlib.cm as cm
import pandas as pd

df = pd.read_excel("kinetic_data.xlsx")
kinetic_data = df.to_numpy()

r1 = 'A -> B + I'
reactions = [r1]
mechanism = make_system(reactions)
exec(mechanism)

timesteps = 30
time = np.linspace(0, 2, timesteps)
t = [0, np.max(time)]
t_eval = list(time)
rate_constants = np.array([1.514])
ic = np.array([4, 0, 0])
solution = solve_ivp(kinetic_model, t, ic, t_eval = t_eval, method = "RK45", 
                     args = rate_constants)
print(solution.y.T)






