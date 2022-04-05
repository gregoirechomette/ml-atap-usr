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
x_scalings = pd.read_csv('../models/' + damage + '/Scaling_parameters.csv')[inputs].to_numpy()
y_scalings = pd.read_csv('../models/' + damage + '/Scaling_parameters.csv')[damage].to_numpy()

# Load the machine learning models already trained
class_model = models.load_model('../models/' + damage + '/Classification_model')
regr_model = models.load_model('../models/' + damage + '/Regression_model')

# Rescale the input data and give it the appropriate format
scenario = pd.DataFrame.from_dict(scenario).to_numpy()
scenario = np.divide(np.subtract(scenario, x_scalings[0:1,:]), x_scalings[1:2,:])

print(" The normalized inputs are: ", scenario)


""" ================== FORWARD PASS =================="""

# With loaded TF model
dummy = np.zeros((1,1))
radius_normalized = regr_model.predict([scenario, dummy, dummy])
radius_of_damage = np.multiply(radius_normalized[0], y_scalings[1:2]) + y_scalings[0:1]
print('The radius according to loaded TF model is: ', max(round(radius_of_damage[0,0],0),0), 'm')

# With recreated TF model
x=Input(shape=(9,), name='input')
y=layers.Dense(64, activation='relu', name='dense1')(x)
y=layers.Dense(128, activation='relu', name='dense2')(y)
y=layers.Dense(256, activation='relu', name='dense3')(y)
y=layers.Dense(1, activation='linear', name='output')(y)
newModel=Model(inputs=x,outputs=y)
newModel.compile(loss='mse',optimizer='adam')

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

new_damage_radius = np.multiply(newModel.predict(scenario)[0], y_scalings[1:2]) + y_scalings[0:1]
print('The radius according to recreated TF model is: ', max(round(new_damage_radius[0],0),0), 'm')

# With manual method
def preActiv(x,w,b):
    return np.matmul(w.T,x) + b

def activ(z):
    return np.maximum(z,0)

def dAdZ(z):
    dadz = np.copy(z)
    for i in range(z.shape[0]):
        if z[i,0] < 0:
            dadz[i,0] = 0
        else:
            dadz[i,0] = 1
            
    dadz = dadz.reshape((dadz.shape[0]))
    return np.diag(dadz)

w1 = regr_model.layers[1].weights[0].numpy()
w2 = regr_model.layers[3].weights[0].numpy()
w3 = regr_model.layers[5].weights[0].numpy()
w4 = regr_model.layers[9].weights[0].numpy()

b1 = regr_model.layers[1].bias.numpy().reshape((64,1))
b2 = regr_model.layers[3].bias.numpy().reshape((128,1))
b3 = regr_model.layers[5].bias.numpy().reshape((256,1))
b4 = regr_model.layers[9].bias.numpy().reshape((1,1))

# First layer
z1 = preActiv(scenario.T, w1, b1)
a1 = activ(z1)

# Second layer
z2 = preActiv(a1, w2, b2)
a2 = activ(z2)

# Third layer
z3 = preActiv(a2, w3, b3)
a3 = activ(z3)

# Fourth layer
z4 = preActiv(a3, w4, b4)

manual_damage_radius = np.multiply(z4[0], y_scalings[1:2]) + y_scalings[0:1]
print("The radius according to the manual method is: ", max(round(manual_damage_radius[0], 0),0), 'm')



""" ================== BACKWARD PASS =================="""

# With recreated TF model
x_input = tf.Variable(scenario, dtype='float32')
with tf.GradientTape() as tape:
    preds = newModel(x_input)
grads = tape.gradient(preds, x_input)
print('The gradients with the TF method are: ', grads.numpy().reshape((9,1)))

# With manual method
der = np.matmul(dAdZ(z3),w4)
der = np.matmul(w3,der)
der = np.matmul(dAdZ(z2), der)
der = np.matmul(w2, der)
der = np.matmul(dAdZ(z1), der)
der = np.matmul(w1, der)
print("The gradients with the manual method are:", der)
