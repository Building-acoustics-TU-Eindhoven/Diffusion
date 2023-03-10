# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 15:25:01 2023

@author: 20225533
"""

import numpy as np

def clarity(time_div, decay, fspatial):
    """
    Clarity `C_i` determined from a Schoeder decay.
    :param time: Time in miliseconds (e.g.: 50, 80).
    :param decay: name of the decay calculated with the diffusion equation.
    :param f_spatial: Frequency sample for the calculation of the x range in the graph.
    """
    
    c = 10*np.log10(np.sum(decay[:time_div]) / np.sum(decay[time_div:]))
    print(c)
    return c

#dt = 0.0001 #distance between grid points on the time discretization [s]
#fspatial = 1/dt #frequency spatial resolution (sampling period)
#bands = np.array([63,125,250,500,1000,2000,4000]) #array of frequency bands
#decay = np.load('w_rec.npy')  

#c50 = clarity(50, decay, fspatial)
