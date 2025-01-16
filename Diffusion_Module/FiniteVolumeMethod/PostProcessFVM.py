# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 11:18:22 2024

@author: 20225533
"""
#%%
###############################################################################
#IMPORT RESULTS FROM FVM.py
###############################################################################
#import numpy as np
#import os
import pickle

def load(filename):
    with open(filename, 'rb') as f:
        variables = pickle.load(f)
        globals().update(variables)

# Example usage
load('results.pkl')
#np.load(os.path.join('variables.pkl'))

#%%
###############################################################################
#POST-PROCESS AND CREATION OF IMPULSE RESPONSE FROM EDC
###############################################################################

#To add 

#%%
###############################################################################
#POST-PROCESS AND FIGURES
###############################################################################
import matplotlib.pyplot as plt
import numpy as np

if tcalc == "decay":
    t_off = t_off - t_off[0] #removing the t_off[0] to make the vector start from zero.
    
    #Decay of SPL in the recording_time
    plt.figure(figsize=(12, 8))
    plt.title('SPL over time at the receiver per octave band')
    for fi in range(nBands):
        plt.subplot(nBands, 1, fi+1)
        plt.plot(t_off, spl_r_off_band[fi], label=f'{center_freq[fi]} Hz')
        plt.xlabel('Time [s]')
        plt.ylabel('SPL [dB]')
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.axhline(0, color='black', linewidth=0.5)
        plt.legend(loc='best')
        plt.show()
    
    #Decay of SPL in the recording_time normalised to maximum 0dB
    plt.figure(figsize=(12, 8))
    plt.title('Normalised SPL over time at the receiver per octave band')
    for fi in range(nBands):
        plt.subplot(nBands, 1, fi+1)
        plt.plot(t, spl_r_norm_band[fi], label=f'{center_freq[fi]} Hz')
        plt.xlabel('Time [s]')
        plt.ylabel('SPL [dB]')
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.axhline(0, color='black', linewidth=0.5)
        plt.legend(loc='best')
        plt.show()
    
    #Schroeder decay
    plt.figure(figsize=(12, 8))
    plt.title('Schroeder decay (Energy Decay Curve) at the receiver per octave band')
    for fi in range(nBands):
        plt.subplot(nBands, 1, fi+1)
        plt.plot(t_off, sch_db_band[fi], label=f'{center_freq[fi]} Hz')
        plt.xlabel('Time [s]')
        plt.ylabel('Energy decay [dB]')
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.axhline(0, color='black', linewidth=0.5)
        plt.legend(loc='best')
        plt.show()
    
    #T30 at the receiver over frequency
    plt.plot(center_freq,t30_band)
    plt.title("T30 over frequency at the receiver")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("T30 [s]")
    #plt.xlim()
    plt.ylim(0, max(t30_band)+1)
    plt.xticks(center_freq)
    
    #edt at the receiver over frequency
    plt.plot(center_freq,edt_band)
    plt.title("Early Decay Time over frequency at the receiver")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("EDT [s]")
    #plt.xlim()
    plt.ylim(0, max(edt_band)+1)
    plt.xticks(center_freq)
    
    #c80 at the receiver over frequency
    plt.plot(center_freq,c80_band)
    plt.title("Clarity over frequency at the receiver")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("C80 [s]")
    #plt.xlim()
    plt.ylim(0, max(c80_band)+1)
    plt.xticks(center_freq)
    
    
    #for iBand in range(nBands):     
        # #Energy density at the receiver over time
        # plt.figure(3)
        # plt.plot(t,w_rec_band[iBand])
        # plt.title("Energy density over time at the receiver")
        # plt.xlabel("t [s]")
        # plt.ylabel("Energy density [kg m^-1 s^-2]")
        # plt.xlim()
        # plt.ylim()
        # plt.xticks(np.arange(0, recording_time +0.1, 0.1))
    
if tcalc == "stationarysource":
    for iBand in range(nBands):

        #Decay of SPL in the recording_time at the receiver
        plt.figure(1)
        plt.plot(t,spl_r_band[iBand]) #plot sound pressure level with Pref = (2e-5)**5
        plt.title("SPL over time at the receiver")
        plt.xlabel("t [s]")
        plt.ylabel("SPL [dB]")
        plt.xlim()
        plt.ylim()
        plt.xticks(np.arange(0, recording_time +0.1, 0.5))
        #plt.yticks(np.arange(0, 120, 20))
    
        #Decay of SPL in the recording_time normalised to maximum 0dB
        plt.figure(2)
        plt.title("Normalised SPL over time at the receiver")
        plt.plot(t,spl_r_norm_band[iBand])
        plt.xlabel("t [s]")
        plt.ylabel("SPL [dB]")
        plt.xlim()
        plt.ylim()
        plt.xticks(np.arange(0, recording_time +0.1, 0.1))
        plt.yticks(np.arange(0, -60, -10))
    
        #Energy density over time at the receiver
        plt.figure(3)
        plt.title("Energy density over time at the receiver")
        plt.plot(t,w_rec_band[iBand])
        plt.ylabel('$\mathrm{Energy \ Density \ [kg/ms^2]}$')
        plt.xlabel("t [s]")
        
        #Sound pressure level stationary over the space y.
        plt.figure(4)
        t_dim = len(t)
        last_time_index = t_dim-1
        plt.title("SPL over the y axis")
        plt.plot(y_axis,spl_stat_y_band[iBand])
        plt.xticks(np.arange(0, 20, 5))
        plt.yticks(np.arange(75, 105, 5))
        plt.ylabel('$\mathrm{Sound \ Pressure\ Level \ [dB]}$')
        plt.xlabel('$\mathrm{Distance \ along \ y \ axis \ [m]}$')
        
        #Sound pressure level stationary over the space x.
        plt.figure(5)
        t_dim = len(t)
        last_time_index = t_dim-1
        plt.title("SPL over the x axis")
        plt.plot(x_axis,spl_stat_x_band[iBand])
        #plt.xticks(np.arange(0, 35, 5))
        plt.yticks(np.arange(90, 97, 1))
        plt.ylabel('$\mathrm{Sound \ Pressure \ Level \ [dB]}$')
        plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$')
        
        #Energy density at t=recording_time over the space x.
        plt.figure(6)
        plt.title("Energy density over the x axis at t=recording_time")
        plt.plot(x_axis,w_rec_x_band[iBand])
        plt.ylabel('$\mathrm{Energy \ Density \ [kg/ms^2]}$')
        plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$')