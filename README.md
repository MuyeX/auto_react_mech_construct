<<<<<<< HEAD
# Simplest Mechanism Builder Algorithm (SiMBA)

Welcome to **SiMBA (Simplest Mechanism Builder Algorithm)**, an automated microkinetic model discovery tool. SiMBA is designed to facilitate the generation of robust, accurate, and computationally efficient microkinetic models from kinetic data, bridging the gap between theoretical exploration and practical applicability in chemical engineering.

## Overview

**SiMBA** is an open-source Python tool that helps automate the discovery of microkinetic models for chemical reactions. Manual construction of these models is often time-consuming and error-prone. SiMBA addresses this challenge by automating the generation, refinement, and evaluation of reaction mechanisms.

The key phases of the algorithm are:
1. **Mechanism Generation**: Uses a parallelized backtracking algorithm to generate feasible reaction pathways, ensuring that only physically sensible mechanisms are proposed.
2. **Mechanism Translation**: Converts generated mechanisms into systems of ordinary differential equations (ODEs), making them executable models.
3. **Parameter Estimation**: Estimates kinetic parameters by minimizing the error between model predictions and experimental data using optimization techniques.
4. **Model Comparison**: Compares generated models using the Akaike Information Criterion (AIC) to balance model accuracy and complexity.

## Features
- Automated mechanism generation based on available kinetic data.
- Parallelized backtracking to efficiently explore possible reaction pathways.
- Model translation from reaction matrices to ODEs for simulation.
- Iterative optimization and model comparison for mechanism refinement.
- Suitable for complex chemical systems such as aldol condensations and biomass dehydration reactions.

## Installation

To install SiMBA, clone the repository from GitHub:

```sh
$ git clone https://github.com/OptiMaL-PSE-Lab/auto_react_mech_construct
$ cd auto_react_mech_construct

=======
Simplest Mechanism Builder Algorithm (SiMBA): An Automated Microkinetic Model Discovery Tool
>>>>>>> 469968f54c02f02f926ab9d0be0ff260c81b0eee
