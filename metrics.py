import numpy as np
np.random.seed(1998)

def NLL_mechanism(exp_data, model_pred):
    number_species = np.shape(exp_data)[1]
    number_datapoints = np.shape(exp_data)[0]
    output = np.zeros(number_species)
    mse = np.zeros(number_species)
    variance = np.zeros(number_species)
    
    for i in range(number_species):
        a = ((exp_data[:, i] - model_pred[:, i])**2)
        mse[i] = np.sum(a)
        variance[i] = mse[i] / (number_datapoints)
    
    for i in range(number_species):
        likelihood = ((exp_data[:, i] - model_pred[:, i])**2) / (2 * (variance[i])) \
            - np.log(1 / (np.sqrt(2 * np.pi * (variance[i]))))
        output[i] = np.sum(likelihood)
    
    return np.sum(output)

def Akaike(nll, params):
    num_param = len(params)
    aic = 2 * nll + 2 * num_param
    
    return aic