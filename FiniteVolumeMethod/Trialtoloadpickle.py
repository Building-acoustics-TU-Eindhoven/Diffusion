# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 09:46:05 2024

@author: 20225533
"""
#import numpy as np
#import os
import pickle

def load(filename):
    with open(filename, 'rb') as f:
        variables = pickle.load(f)
        globals().update(variables)

# Example usage
load('results.pkl')

#np.load(os.path.join('variables.pkl'))