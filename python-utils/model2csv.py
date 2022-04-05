#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# System integration
import os
import sys
import abc

# Standard python libraries
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ML libraries
import tensorflow as tf
    
# Load the model
layer_units = list([64,128,256])
model = tf.keras.models.load_model('../models/BlastRad1/Regression_model')
print(model.summary())

# Function to write the weights to csv
def writeWeights(f, w, b):
    for i in range(w.shape[0]):
        for j in range(w.shape[1]):
            f.write(str(w[i,j]) + ',')
    f.write('\n')
    for i in range(b.shape[0]):
        f.write(str(b[i,0]) + ',')
    f.write('\n')

# Retrieve the weights
w1 = model.layers[1].weights[0].numpy()
w2 = model.layers[3].weights[0].numpy()
w3 = model.layers[5].weights[0].numpy()
w4 = model.layers[9].weights[0].numpy()
 
b1 = model.layers[1].bias.numpy().reshape((layer_units[0],1))
b2 = model.layers[3].bias.numpy().reshape((layer_units[1],1))
b3 = model.layers[5].bias.numpy().reshape((layer_units[2],1))
b4 = model.layers[9].bias.numpy().reshape((1,1))

# Write the radii in a .csv file
fileName = 'weights-2.csv'
os.system('rm ' + fileName)
f = open(fileName, 'a+')
writeWeights(f, w1, b1)
writeWeights(f, w2, b2)
writeWeights(f, w3, b3)
writeWeights(f, w4, b4)
f.close()
