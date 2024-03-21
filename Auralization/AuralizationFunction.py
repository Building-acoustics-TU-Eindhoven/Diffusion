# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 14:33:58 2024

@author: 20225533
"""

import numpy as np
#import matplotlib
import matplotlib.pyplot as plt
from FunctionRT import *
import sounddevice as sd
import soundfile as sf
import scipy
from scipy.io import wavfile

#Import anechoic signal
filename = 'C:/Users/20225533/Diffusion/Auralization/Frequency(english).wav' #name of the anechoic signal file

# Extract data and sampling rate from file
data_signal, fs = sf.read(filename) #this returns "data_signal", which is the 
#audiodata (one_dimentional array) of the anechoic signal. It returns also the
#"fs" sample frequency of the signal
 
#sd.play(data, fs) #this line allows to listen to the anechoic signal as it is.
#status = sd.wait()  #Wait until file is done playing

#Import the energy decay curve
edc = np.load('C:/Users/20225533/Diffusion/Auralization/w_rec_off.npy') #energy decay curve taken from the results of the diffusion equation model
t = np.load('C:/Users/20225533/Diffusion/Auralization/t.npy') #time steps array
edc_deriv = np.load('C:/Users/20225533/Diffusion/Auralization/w_rec_off_deriv.npy') #energy decay curve differentiated (or also impulse response of the room)
t_off = np.load('C:/Users/20225533/Diffusion/Auralization/t_off.npy') #decay time of the energy decay curve
t_off = t_off - t_off[0] #removing the t_off[0] to make the vector start from zero.

sch_db = 10.0 * np.log10(edc / max(edc)) #level of the array: schroeder decay
idx_w_rec = 17000
t60 = t60_decay(t, sch_db, idx_w_rec) #called function for calculation of t60 [s] before auralization


edc_deriv = np.load('C:/Users/20225533/Diffusion/Auralization/w_rec_off_deriv.npy') #energy decay curve differentiated (or also envelope of the impulse response of the room)
t_off = np.load('C:/Users/20225533/Diffusion/Auralization/t_off.npy') #decay time of the energy decay curve
t_off = t_off - t_off[0] #removing the t_off[0] to make the vector start from zero.

#Random noise creation
random_array = np.random.rand(1, len(edc_deriv))*2 - 1 #random noise vector with numbers between -1 and 1
random_array = sum(random_array) #this line of code is used for passing from a row vector to a column vector

#From the envelope of the impulse response, we need to get the impulse response
square_root = np.sqrt(edc_deriv)

#Multiplication with random noise
imp_rand = square_root*random_array 

#Create a file wav for impulse response
scipy.io.wavfile.write("imp_resp.wav", fs, imp_rand)

#Play the impulse response
sd.play(imp_rand, fs)

#Convolution of the impulse_rand with the anechoic signal
st = np.arange(0,(len(data_signal))/fs,1/fs) #Time vector of the speech signal
ht = np.arange(0,(len(imp_rand))/fs,1/fs)  #Time vector of the room impulse response

#Create impulse response
sh_conv = np.convolve(imp_rand,data_signal) #convolution of the impulse response with the anechoic signal
sh_conv = sh_conv/max(abs(sh_conv)) #normalized to the maximum value of the convolved signal

t_conv = np.arange(0,(len(sh_conv))/fs,1/fs) #Time vector of the convolved signal

plt.plot(st,data_signal) #plot the anechoic signal
 
plt.plot(ht,imp_rand) #plot the impulse response

plt.plot(t_conv,sh_conv) #plot the convolved signal

#Play the convolved signal
#sd.play(sh_conv, fs)

###############################################################################
################################### CHECKS ####################################
#Backward integration of the imp_rand to see if the RT matches the one calculated from the diffusion equation model

env_imp_resp = imp_rand**2 #envelope of the impulse response created by squaring the impulse response

# Schroeder integration
sch = np.cumsum(env_imp_resp[::-1])[::-1] #backward integration
sch_db = 10.0 * np.log10(sch / np.max(sch))

#Reverberation time calculation
init = -5.0 #because I want the T30, I need to start at -5
end = -35.0 #because I want the T30, I need to finish at -35
factor = 2.0 #factor of 2 since I need the T30

#Linear regression
idxL1 = np.where(sch_db <= init)[0][0] #index at which the rtdecay is equal to -5
idxL2 = np.where(sch_db <= end)[0][0] #index at which the rtdecay is equal to -35
   
timeL1 = t[idxL1] #index at which the time vector is equal to the idxL1
timeL2 = t[idxL2] #index at which the time vector is equal to the idxL2

# Classical Approach (T30 = 2*(t[-35dB]-t[-5dB])
RTCalc = factor*(timeL2 - timeL1)

x = t[idxL1:idxL2]
y = sch_db[idxL1:idxL2]

# Linregress approach
from scipy import stats
slope,intercept = stats.linregress(t[idxL1:idxL2],sch_db[idxL1:idxL2])[0:2] #calculating the slope and the interception of the line connecting the two points
db_regress_init = (init - intercept) / slope #dB initial
db_regress_end = (end - intercept) / slope #dB End
t60_after = factor * (db_regress_end - db_regress_init) #t60 according to linregress approach

# Poly-based Approach y = Ax + B
CoefAlpha = np.polyfit(t[idxL1:idxL2], sch_db[idxL1:idxL2], 1) ##calculating the slope and the interception of the line connecting the two points
t60_after = (-60/CoefAlpha[0]) #t60 according to polyfit approach






