# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 15:25:01 2023

@author: Ilaria Fichera
"""
import math

def definition(t60, V, Eq_A, S, c0, dist):
    """
    Calculation of definition D50 from a SPL decay using Barron's revised formula [Vorlander 2008]

    Parameters
    ----------
        t60 : float
            Reverberation time T_{60}
        V : float
            Volume of the room
        Eq_A : array of floats
            Equivalent absoption area of the room
        S : float
            Total surface area of the rom
        c0 : int 
            Speed of sound
        dist: float
            Distance between source and receiver

    Returns
    -------
        d50 : float
            Definition D_{50}
    """  
    d50 = 100*(((t60/13.8*V)*(math.exp(-(Eq_A/S))-math.exp(-((0.69/t60)+(Eq_A/S)))) + (1/(4*math.pi*c0*dist**2)))/
                      (((t60/13.8*V)*(math.exp(-(Eq_A/S)))) + (1/(4*math.pi*c0*dist**2))))
    
#    d50 = 100*(((t60/13.8*V)*(1-math.exp(-(276/t60)))+ (1/(4*math.pi*c0*dist**2)))/
#               ((t60/13.8*V) + (1/(4*math.pi*c0*dist**2))))
    
    return d50


