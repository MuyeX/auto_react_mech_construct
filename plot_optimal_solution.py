import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
np.random.seed(1998)
from parameter_estimation import read_files

def plot(unique_letters, model_predictions):

    in_silico_data = read_files("auto_react_mech_construct/exp_data")

    # This takes the first column from each entry of the dictionary and puts it into another dictionary
    initial_conditions = {}
    for key, value in in_silico_data.items():
        aa = "ic_" + key[-1]
        initial_conditions[aa] = value[0]

    num_exp = len(initial_conditions)
    timesteps = 30
    time = np.linspace(0, 2, timesteps)
    t = [0, np.max(time)]
    t_eval = list(time)


    species = list(sorted(unique_letters))
    color_1 = cm.plasma(np.linspace(0, 1, 3))
    marker = ['o' for i in range(3)]

    # Plotting the in-silico data for visualisation
    for i in range(num_exp):
        fig, ax = plt.subplots()
        # ax.set_title("Experiment " + str(i + 1))
        ax.set_ylabel("Concentrations $(M)$", fontsize = 18)
        ax.set_xlabel("Time $(h)$", fontsize = 18)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.tick_params(axis = 'both', which = 'major', labelsize = 18)
        
        
        for j in np.array([0, 1, -1]):
            y = in_silico_data["exp_" + str(i + 1)][:, j]
            yy = model_predictions["exp_" + str(i + 1)][:, j]
            ax.plot(time, y, marker[j], markersize = 4, label = species[j], color = color_1[j])
            ax.plot(time, yy, '-', color = color_1[j])
        
        ax.grid(alpha = 0.5)
        ax.legend(loc='upper right', fontsize = 15)
        
    plt.show()