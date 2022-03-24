#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# System integration
import os
import secrets
import sys
import abc

# Standard python libraries
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ML libraries
import tensorflow as tf
from tensorflow.keras import Model, Input, regularizers, layers, models
from tensorflow.keras.models import Sequential, load_model


# Impact conditions
scenario = {'Diameter [m]':             [70],
            'Density [kg/m^3]':         [3500],
            'Strength [Pa]':            [1e5],
            'Alpha coefficient [-]':    [0.2],
            'Velocity [m/s]':           [1e4],
            'Incidence angle [degree]': [45],
            'Azimuth [degrees]':        [180],
            'LumEff [-]':               [3e-3],
            'Ablation [kg/J]':          [1e-9]}

# Level of damage considered
damage = 'BlastRad1'

# Load the scaling parameters
inputs = ['Diameter', 'Density', 'Strength', 'Velocity', 'Angle', 'Azimuth', 'Alpha', 'LumEff', 'Ablation']
x_scalings = pd.read_csv('./models/' + damage + '/Scaling_parameters.csv')[inputs].to_numpy()
y_scalings = pd.read_csv('./models/' + damage + '/Scaling_parameters.csv')[damage].to_numpy()

# Load the machine learning models already trained
class_model = models.load_model('./models/' + damage + '/Classification_model')
regr_model = models.load_model('./models/' + damage + '/Regression_model')

# Rescale the input data and give it the appropriate format
scenario = pd.DataFrame.from_dict(scenario).to_numpy()
scenario = np.divide(np.subtract(scenario, x_scalings[0:1,:]), x_scalings[1:2,:])

print(" The normalized inputs are: ", scenario)


# # FORWARD PASS

# # Predict the probability of hazard
# probability_of_threat = class_model.predict(scenario)[0,0]

# if probability_of_threat > 0.5:
#     # Predict the size of the affected region
#     dummy = np.zeros((1,1))
#     radius_normalized = regr_model.predict([scenario, dummy, dummy])
#     print("Radius unscaled: ", radius_normalized)
#     radius_of_damage = np.multiply(radius_normalized[0], y_scalings[1:2]) + y_scalings[0:1]
#     print('The area of damage has a radius estimated to', round(radius_of_damage[0,0],1), 'm')

# else:
#     # Print info message
#     print('The entry conditions specified are not expected to produce this level of damage')
    

# BACKWARD PASS

# Create a new Tensorflow model with the same weights
x=Input(shape=(9,), name='input')
y=layers.Dense(64, activation='relu', name='dense1')(x)
y=layers.Dense(128, activation='relu', name='dense2')(y)
y=layers.Dense(256, activation='relu', name='dense3')(y)
y=layers.Dense(1, activation='linear', name='output')(y)
newModel=Model(inputs=x,outputs=y)
newModel.compile(loss='mse',optimizer='adam')

# Prepare the weights for the new model
layer1 = []; layer2 = []; layer3 = []; layer4 = []
layer1.append(regr_model.layers[1].weights[0].numpy())
layer1.append(regr_model.layers[1].bias.numpy().reshape((64,)))
layer2.append(regr_model.layers[3].weights[0].numpy())
layer2.append(regr_model.layers[3].bias.numpy().reshape((128,)))
layer3.append(regr_model.layers[5].weights[0].numpy())
layer3.append(regr_model.layers[5].bias.numpy().reshape((256,)))
layer4.append(regr_model.layers[9].weights[0].numpy())
layer4.append(regr_model.layers[9].bias.numpy().reshape((1,)))

newModel.layers[1].set_weights(layer1)
newModel.layers[2].set_weights(layer2)
newModel.layers[3].set_weights(layer3)
newModel.layers[4].set_weights(layer4)

# Backward pass
x_input = tf.Variable(scenario, dtype='float32')
with tf.GradientTape() as tape:
    preds = newModel(x_input)
grads = tape.gradient(preds, x_input)


dummy = np.zeros((1,1))
radius_normalized = regr_model.predict([scenario, dummy, dummy])
print("Radius unscaled: ", radius_normalized)
radius_of_damage = np.multiply(radius_normalized[0], y_scalings[1:2]) + y_scalings[0:1]

print("The predictions of the previous model are: ", np.multiply(regr_model.predict([scenario, np.zeros((1,1)), np.zeros((1,1))])[0], y_scalings[1:2]) + y_scalings[0:1])

print("The predictions of the new model are: ", np.multiply(newModel.predict(scenario)[0], y_scalings[1:2]) + y_scalings[0:1])
print("The derivatives of the new model are: ", grads)

