# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 15:25:01 2023

@author: Ilaria Fichera
"""

def centretime(t60, Eq_A, S):
    """
    Clarity determined from a SPL decay using Barron's revised formula [Vorlander 2008]
    :param t60: Reverberation time in s
    :param Eq_A: Equivalent absoption area of the room
    :param S: total surface area of the rom
    """
#    ts = 1000*((t60/13.8) * ((Eq_A/S)+1)) #with direct and reverberant field
    
    ts = 1000*((t60/13.8)) #with direct field only
    return ts

