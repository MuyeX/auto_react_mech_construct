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
```

Ensure that you have the required Python libraries installed. You can install them via the provided requirements.txt file:

```sh
$ pip install -r requirements.txt
```

## Usage

The kinetic data available should be loaded into the directory. This should be added as a single folder, where each file is a different experiment (each experiment should be named as "exp_{number of experiment}.csv").

To record the output from the algorithm, define an output file (e.g., name_file = "output_hypoth.log") in main.py.

In parameter_estimation.py define the name of the folder within the directory that has the kinetic data (e.g., name_file = "exp_data_hypoth").

In parameter_estimation.py ensure that num_observable_species is defined according to the system (line 75 and 153).

At the bottom of main.py, define elementary_reactions (i.e., the number of elementary steps of the smallest possible mechanism given the stoichiometry of the reaction), number_species (i.e.,the number of species present in the smallest possible mechanism given the stocihiometry of the reaction), stoichiometry (i.e., the stoichiometric coefficients of the main reactants and products), intermediate (i.e., this is the position within the vector from where the intermediates will start showing; this number is just the same as the number_species), product (i.e., this is the position within the vector from where the products will start showing), reactant (i.e., this is the position within the vector from where the reactants will start showing; this is always 0), and finally time_budget (i.e., the number of seconds SiMBA will explore a specific layer of complexity).

Once all this is done, run the main.py file.

## Example Case Studies

SiMBA has been tested on several case studies, including:

- Hypothetical Reaction: A system involving five species and four elementary steps, demonstrating the algorithm's versatility in handling unknown intermediates, as well as first and second order elementary steps.
- Aldol Condensation Reaction: The condensation between benzaldehyde and acetophenone, showcasing SiMBA's capability in practical contexts.
- Dehydration of Fructose: Conversion to 5-hydroxymethylfurfural, illustrating SiMBA's application in complex chemical systems where only half of the species involved in the mechanism are observable.

## Limitations and Future Work

While SiMBA significantly enhances the exploration of microkinetic models, it does not inherently chemically identify intermediates, requiring expert input for complex systems. Future improvements will focus on incorporating chemical identification and uncertainty quantification to make the models even more robust and autonomous.

## Acknowledgments

This work was supported by the Engineering and Physical Sciences Research Council (EPSRC) under grant EP/S023232/1.
