# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 08:56:15 2023

@author: 20225533
"""
#import matplotlib.pyplot as plt #import matplotlib as mpl
import numpy as np
from scipy import stats

def t60_decay(t, sch_db, idx_w_rec, rt='t30'):
    """
    Reverberation time from a Schroeder decay (Schroeder 1965)
    :param t: time axis
    :param sch_db: name of the Schroeder decay calculated with the diffusion equation.
    :param idx_w_rec: Index at which the t array is equal to the sourceon_time - calculation of RT from when the source stops.
    :returns: Reverberation time T_{60}
    """    
    if rt == 't30':
        init = -5.0 #because I want the T30, I need to start at -5
        end = -35.0 #because I want the T30, I need to finish at -35
        factor = 2.0 #factor of 2 since I need the T30
    elif rt == 't20':
        init = -5.0 
        end = -25.0 
        factor = 3.0 
    elif rt == 't60':
        init = -5.0 
        end = -65.0 
        factor = 1.0 
    elif rt == 't10':
        init = -5.0
        end = -15.0
        factor = 6.0
    elif rt == 'edt':
        init = 0.0
        end = -10.0
        factor = 6.0
    
    #rt_decay = sch_db #decay to be used to calculate the RT
    
    #Linear regression
    idxL1 = np.argmin(np.abs(sch_db - init)) #np.where(sch_db <= init)[0][0] #index at which the rtdecay is equal to -5
    idxL2 = np.argmin(np.abs(sch_db - end)) #np.where(sch_db <= end)[0][0] #index at which the rtdecay is equal to -35
       
    timeL1 = t[idxL1] #index at which the time vector is equal to the idxL1
    timeL2 = t[idxL2] #index at which the time vector is equal to the idxL2
    
    # Classical Approach (T30 = 2*(t[-35dB]-t[-5dB])
    RTCalc = factor*(timeL2 - timeL1)
    
    # Linregress approach
    slope,intercept = stats.linregress(t[idxL1:idxL2],sch_db[idxL1:idxL2])[0:2] #calculating the slope and the interception of the line connecting the two points
    db_regress_init = (init - intercept) / slope #dB initial
    db_regress_end = (end - intercept) / slope #dB End
    t60I = factor * (db_regress_end - db_regress_init) #t60 according to linregress approach
    
    # Poly-based Approach y = Ax + B
    CoefAlpha = np.polyfit(t[idxL1:idxL2], sch_db[idxL1:idxL2], 1) ##calculating the slope and the interception of the line connecting the two points
    t60 = (-60/CoefAlpha[0]) #t60 according to polyfit approach
    
    y_axis = (slope*t[idx_w_rec:] + intercept) + slope

    # plt.figure(3)
    # plt.plot(t[idx_w_rec:],sch_db, color ='b', linewidth = 1.8)
    # plt.plot(t[idx_w_rec:],y_axis,color='r',linewidth=2)
    # plt.plot(t[idx_w_rec:][idxL1],np.real(sch_db[idxL1]),'o',linewidth=2)
    # plt.plot(t[idx_w_rec:][idxL2],np.real(sch_db[idxL2]),'o',linewidth=2)
    # plt.axvline(x=t60,ymin=-100,ymax=0,linestyle='--',linewidth=2)
    # plt.ylabel('Normalized Magnitude (dB)')
    # plt.xlabel('Time (s)')
    # plt.legend(['EDC','Line Fitting','Upper Point','Lower Point','Estimate T_{30}'])
    # plt.title("Figure 3: 'T_{30} = ' + str(round(t60,2)) + ' s.'")
    # plt.grid(True)
    # plt.ylim([-100,0])
    # #plt.show()
    
    return t60
