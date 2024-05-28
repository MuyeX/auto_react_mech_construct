#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 12:29:09 2024

@author: md1621
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import mean_squared_error

# Generate synthetic data
np.random.seed(42)
X = np.linspace(0, 10, 50)
y_true = 2 * X**3 - 5 * X**2 + X + 3
y = y_true + np.random.normal(scale=75, size=X.shape)

# Function to fit and plot polynomial regression models
def plot_polynomial_fit(X, y, degrees):
    plt.figure(figsize=(32, 16))

    for i, degree in enumerate(degrees, start=1):
        plt.subplot(1, 2, i)

        # Transform input data for polynomial regression
        poly = PolynomialFeatures(degree)
        X_poly = poly.fit_transform(X.reshape(-1, 1))

        # Fit the polynomial regression model
        model = LinearRegression()
        model.fit(X_poly, y)

        # Predict using the model
        y_pred = model.predict(X_poly)

        # Calculate mean squared error
        mse = mean_squared_error(y, y_pred)
        aic = mse + 2 * (degree + 1)

        # Plot the data and the model
        plt.scatter(X, y, color='blue', label='Data', s=300, alpha=0.7)
        plt.plot(X, y_true, color='green', linestyle='--', label='True Function', linewidth=8)
        plt.plot(X, y_pred, color='red', label=f'AIC: {aic:.2f}', linewidth=8)
        plt.xlabel('X', fontsize=52)
        plt.ylabel('y', fontsize=52)
        plt.legend(fontsize=50)
        plt.title(f'Polynomial Regression Degree {degree}', fontsize=56)
        plt.grid(True, linestyle='-', alpha=0.6)
        
        # Adjust tick parameters
        plt.tick_params(axis='both', which='major', labelsize=50)  # Set major tick label size
        
        # Remove the top and right spines
        ax = plt.gca()  # Get current axis
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig('model_comparison_poster', dpi = 600, bbox_inches = "tight")
    plt.show()

# Define degrees of polynomials to fit
degrees = [2, 3]

# Plot polynomial fits
plot_polynomial_fit(X, y, degrees)


# # Generate data for model complexity
# model_complexity = np.linspace(1, 10, 100)

# # Generate arbitrary training and test error data
# # Training error asymptotically approaches zero
# training_error = 1 / (model_complexity * 0.5)

# # Test error forms a convex curve (U-shape) with a minimum
# test_error = (model_complexity - 5)**2 / 10 + 0.6
# AIC_value = (model_complexity - 5)**2 / 10 + 0.8

# # Plotting the training and test error
# plt.figure(figsize=(32, 16))
# plt.plot(model_complexity, training_error, label='Training Error', color='blue', linewidth=8)
# plt.plot(model_complexity, test_error, label='Test Error', color='red', linestyle='--', linewidth=8)
# plt.plot(model_complexity, AIC_value, label='AIC Value', color='green', linestyle='-.', linewidth=8)

# # Adding titles and labels
# plt.title('Training Error vs Test Error', fontsize=56)
# plt.xlabel('Model Complexity', fontsize=52)
# plt.ylabel('Error', fontsize=52)
# plt.legend(fontsize=50)
# plt.grid(True, linestyle='-', alpha=0.6)

# # Adjust tick parameters
# plt.tick_params(axis='both', which='both', length=0, labelsize=0)

# # Remove the top and right spines
# ax = plt.gca()  # Get current axis
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)

# # Save the plot as SVG
# plt.savefig('training_test_error', dpi = 600, bbox_inches = "tight")

# # Show the plot
# plt.show()
