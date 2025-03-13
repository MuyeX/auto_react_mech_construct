# read a json file with the following content:
# elementary_reactions = 1
# number_species = 4
# stoichiometry = [-1, -1, 1, 1]
# intermediate = 4
# product = 2
# reactant = 0
# time_budget = 10
# name_file = "aldol_test.log"
# num_observable_species = 4
# name_file = "exp_data_aldol_condensation"
# initial_conditions = {
#     "ic_1": np.array([4, 0, 0, 0, 0, 0]),
#     "ic_2": np.array([6, 2, 1, 0, 0, 0]),
#     "ic_3": np.array([4, 2, 0, 0, 0, 0]),
#     "ic_4": np.array([6, 0, 0, 0, 0, 0]),
#     "ic_5": np.array([6, 2, 0, 0, 0, 0])
# }


import json
import numpy as np


def load_config_file(file_name: str):
    """
    Load a JSON file with the configuration parameters for the model.
    Args:
        file_name: file name of the JSON file

    Returns:

    """

    # Load the JSON file
    with open(file_name, 'r') as file:
        config_data = json.load(file)

    # Extract parameters
    elementary_reactions = config_data['elementary_reactions']
    number_species = config_data['number_species']
    stoichiometry = config_data['stoichiometry']
    intermediate = config_data['intermediate']
    product = config_data['product']
    reactant = config_data['reactant']
    time_budget = config_data['time_budget']
    log_file = config_data['log_file']
    input_dir = config_data['input_dir']
    num_observable_species = config_data['num_observable_species']
    initial_conditions= {key: np.array(value) for key, value in config_data['initial_conditions'].items()}
    config_data['initial_conditions'] = initial_conditions

    # Print the loaded parameters to verify
    print("Elementary Reactions:", elementary_reactions)
    print("Number of Species:", number_species)
    print("Stoichiometry:", stoichiometry)
    print("Intermediate:", intermediate)
    print("Product:", product)
    print("Reactant:", reactant)
    print("Time Budget:", time_budget)
    print("Name File:", log_file)
    print("Input Directory:", input_dir)
    print("Number of Observable Species:", num_observable_species)
    print("Initial Conditions:", initial_conditions)

    return config_data