# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 15:25:01 2023

@author: Ilaria Fichera
"""

def centretime(t60, Eq_A, S):
    """
    Calculation of center time Ts from a SPL decay using Barron's revised formula [Vorlander 2008]

    Parameters
    ----------
        t60 : float
            Reverberation time T_{60}
        Eq_A : array of floats
            Equivalent absoption area of the room
        S : float
            Total surface area of the rom

    Returns
    -------
        c80 : float
            Clarity C_{80}
    """  

    ts = 1000*((t60/13.8) * ((Eq_A/S)+1)) #with direct and reverberant field
    
    #ts = 1000*((t60/13.8)) #with direct field only
    return ts

