# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 15:25:01 2023

@author: 20225533
"""
import math
import numpy as np

def centretime(t60, Eq_A, S):
    """
    Clarity determined from a SPL decay using Barron's formula
    :param t60: Reverberation time in s
    :param Eq_A: Equivalent absoption area of the room
    :param S: total surface area of the rom
    """
    ts = (t60/13.8) * ((Eq_A/S)+1)
    return ts
