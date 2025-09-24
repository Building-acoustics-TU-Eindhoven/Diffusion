# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 15:25:01 2023

@author: Ilaria Fichera
"""
import math
import numpy as np

def clarity(t60, V, Eq_A, S, c0, dist):
    """
    Calculation of clarity C80 from a SPL decay using Barron's revised formula [Vorlander 2008]

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
        c80 : float
            Clarity C_{80}
    """  
    
    c80 = 10*np.log10(((t60/13.8*V)*(math.exp(-(Eq_A/S))-math.exp(-((1.104/t60)+(Eq_A/S)))) + (1/(4*math.pi*c0*dist**2)))/
                      ((t60/13.8*V)*math.exp(-((1.104/t60)+(Eq_A/S))))) 
    
#    c80 = 10*np.log10(((math.exp((13.8*t60)/V))*(1+((13.8*V)/(4*math.pi*c0*(dist**2)*t60))))-1) #equation for A<<S or alpha<<1
    
    return c80

#dt = 0.0001 #distance between grid points on the time discretization [s]
#fspatial = 1/dt #frequency spatial resolution (sampling period)
#bands = np.array([63,125,250,500,1000,2000,4000]) #array of frequency bands
#decay = np.load('w_rec.npy')  

#c50 = clarity(50, decay, fspatial)

