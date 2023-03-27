# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 08:56:15 2023

@author: 20225533
"""
import matplotlib.pyplot as plt #import matplotlib as mpl
import numpy as np
from scipy import stats

def t60_decay(t, w_rec, fsample, rt='t30'):
    """
    Reverberation time from a SPL decay - Schroeder integration (Schroeder 1965)
    :param t: time axis
    :param decay: name of the SPL decay calculated with the diffusion equation.
    :param fsample: Frequency sample for the calculation of the x range in the graph.
    :param rt: Reverberation time estimator - 't30', 't20', 't10', 'edt'
    :returns: Reverberation time T_{60}
    """        
    if rt == 't30':
        init = -5.0
        end = -35.0
        factor = 2.0
    elif rt == 't20':
        init = -5.0
        end = -25.0
        factor = 3.0
    elif rt == 't10':
        init = -5.0
        end = -15.0
        factor = 6.0
    elif rt == 'edt':
        init = 0.0
        end = -10.0
        factor = 6.0

    #t60 = np.zeros(bands.size) #this will need to be included if there is a frequency dependent t60 to calculate.
    
    sch = np.cumsum(w_rec[::-1]**2)[::-1]
    sch_db = 10.0 * np.log10(sch / np.max(sch))
    
    #Linear regression
    zerodecay = np.where(sch_db == 0)[0][0] #index at which the normalised spl is zero
    rt_decay = sch_db[zerodecay:] #decay of spl removing all the part before the zerodecay
    array = np.abs(rt_decay - init) #array of the decay modified by the starting point of the rt calculation
    array_min_idx = np.abs(rt_decay - init).argmin() #index at which there is the minimum value for the start of the linear decay
    array_min = rt_decay[array_min_idx] #minimum value (should correspond more or less with the "init" value)

    array_max_idx = np.abs(rt_decay - end).argmin() #index at which there is the minimum value for the end of the linear decay
    array_max = rt_decay[array_max_idx] #minimum value (should correspond more or less with the "end" value)
    
    init_sample = np.where(rt_decay == array_min)[0][0] #initial value of sch_db in x axis; it is the index of the position of sch_init when it is equal to decay; twice the [0] because before it give the first element of a tuple and then the first element of the list
    end_sample = np.where(rt_decay == array_max)[0][0] #end value of sch_db in x axis; it is the index of the position of sch_init when it is equal to decay; twice the [0] because before it give the first element of a tuple and then the first element of the list
    x = np.arange(init_sample, end_sample + 1) / fsample #spaced x value divided the sample frequency
    y = rt_decay[init_sample:end_sample + 1] #y values correspondent to the x values of sch_dB
    slope, intercept = stats.linregress(x, y)[0:2] #calculates the slope and intercept of a linear least-squares regression from the two sets of measurements given
    res = stats.linregress(x, y) #linear regression
    #plt.plot(t[zerodecay:],rt_decay)
    #plt.plot(x,res.intercept + res.slope*x, 'r', label='fitted line') 
    

    #Reverberation time (T30, T20, T10 or EDT)
    db_regress_init = (init - intercept) / slope 
    db_regress_end = (end - intercept) / slope
    t60 = factor * (db_regress_end - db_regress_init)
    return t60