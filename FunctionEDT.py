# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 08:56:15 2023

@author: 20225533
"""
import matplotlib.pyplot as plt #import matplotlib as mpl
import numpy as np
from scipy import stats

def edt_decay(t, sch_db, idx_w_rec):
    """
    Reverberation time from a Schroeder decay (Schroeder 1965)
    :param t: time axis
    :param sch_db: name of the Schroeder decay calculated with the diffusion equation.
    :param idx_w_rec: Index at which the t array is equal to the sourceon_time - calculation of EDT from when the source stops.
    :returns: Early Decay Time edt
    """    
    init = 0.0 #because I want the edt, I need to start at 0
    end = -10.0 #because I want the edt, I need to finish at -10
    factor = 6.0 #factor of 2 since I need the edt
    
    #Linear regression
    idxL1 = np.where(sch_db <= init)[0][0] #index at which the rtdecay is equal to -5
    idxL2 = np.where(sch_db <= end)[0][0] #index at which the rtdecay is equal to -35
       
    timeL1 = t[idxL1] #index at which the time vector is equal to the idxL1
    timeL2 = t[idxL2] #index at which the time vector is equal to the idxL2
    
    # Classical Approach (T30 = 2*(t[-35dB]-t[-5dB])
    RTCalc = factor*(timeL2 - timeL1)
    
    # Linregress approach
    slope,intercept = stats.linregress(t[idxL1:idxL2],sch_db[idxL1:idxL2])[0:2] #calculating the slope and the interception of the line connecting the two points
    db_regress_init = (init - intercept) / slope #dB initial
    db_regress_end = (end - intercept) / slope #dB End
    edtI = factor * (db_regress_end - db_regress_init) #t60 according to linregress approach
    
    # Poly-based Approach y = Ax + B
    CoefAlpha = np.polyfit(t[idxL1:idxL2], sch_db[idxL1:idxL2], 1) ##calculating the slope and the interception of the line connecting the two points
    edt = (-60/CoefAlpha[0]) #t60 according to polyfit approach
        
    y_axis = (slope*t[idx_w_rec:] + intercept) + slope

    plt.figure()
    plt.plot(t[idx_w_rec:],sch_db, color ='b', linewidth = 1.8)
    plt.plot(t[idx_w_rec:],y_axis,color='r',linewidth=2)
    plt.plot(t[idx_w_rec:][idxL1],np.real(sch_db[idxL1]),'o',linewidth=2)
    plt.plot(t[idx_w_rec:][idxL2],np.real(sch_db[idxL2]),'o',linewidth=2)
    plt.axvline(x=edt,ymin=-100,ymax=0,linestyle='--',linewidth=2)
    plt.ylabel('Normalized Magnitude (dB)')
    plt.xlabel('Time (s)')
    plt.legend(['EDC','Line Fitting','Upper Point','Lower Point','Estimate EDT'])
    plt.title('EDT = ' + str(round(edt,2)) + ' s.')
    plt.grid(True)
    plt.ylim([-100,0])
    plt.show()
    
    return edt