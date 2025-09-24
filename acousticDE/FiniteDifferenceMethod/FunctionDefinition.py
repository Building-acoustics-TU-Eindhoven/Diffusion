# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 15:25:01 2023

@author: Ilaria Fichera
"""
import math

def definition(t60, V, Eq_A, S, c0, dist):
    """
    Definition determined from a SPL decay using Barron's revised formula [Vorlander 2008]
    :param t60: Reverberation time in s
    :param V: volume of the room
    :param Eq_A: Equivalent absoption area of the room
    :param S: total surface area of the rom
    :param c0: sound speed
    :param dist: distance between source and receiver
    """

    d50 = 100*(((t60/13.8*V)*(math.exp(-(Eq_A/S))-math.exp(-((0.69/t60)+(Eq_A/S)))) + (1/(4*math.pi*c0*dist**2)))/
                      (((t60/13.8*V)*(math.exp(-(Eq_A/S)))) + (1/(4*math.pi*c0*dist**2))))
    
#    d50 = 100*(((t60/13.8*V)*(1-math.exp(-(276/t60)))+ (1/(4*math.pi*c0*dist**2)))/
#               ((t60/13.8*V) + (1/(4*math.pi*c0*dist**2))))
    
    return d50


