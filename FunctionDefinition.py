# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 15:25:01 2023

@author: 20225533
"""

import numpy as np

def definition(time, decay, fspatial):
    """
    Definition `D_i` determined from a Schoeder decay.
    :param time: Time in miliseconds (e.g.: 50, 80).
    :param decay: name of the decay calculated with the diffusion equation.
    :param f_spatial: Frequency sample for the calculation of the x range in the graph.
    """

    t = int((time / 1000.0) * fspatial + 1)
    d = 100*((np.sum(decay[:t]) / np.sum(decay[:])))
    return d

