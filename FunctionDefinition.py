# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 15:25:01 2023

@author: 20225533
"""
import math
import numpy as np

def definition(t60, V, Eq_A, S, c0, dist):
    """
    Definition determined from a SPL decay using Barron's formula
    :param t60: Reverberation time in s
    :param V: volume of the room
    :param Eq_A: Equivalent absoption area of the room
    :param S: total surface area of the rom
    :param c0: sound speed
    :param dist: distance between source and receiver
    """

    d50 = 100*(((t60/13.8*V)*(math.exp(-(Eq_A/S))-math.exp(-((0.69/t60)+(Eq_A/S)))) + (1/(4*math.pi*c0*dist**2)))/
                      (((t60/13.8*V)*(math.exp(-(Eq_A/S)))) + (1/(4*math.pi*c0*dist**2))))
    return d50

