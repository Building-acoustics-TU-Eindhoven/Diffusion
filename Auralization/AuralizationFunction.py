# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 14:33:58 2024

@author: 20225533
"""

import numpy as np
#import matplotlib
import matplotlib.pyplot as plt

#import soundfile as sf
#data, samplerate = sf.read('Frequency (english).wav')
#import scipy.io.wavfile as wav

filename = 'C:/Users/20225533/Diffusion/Auralization/Frequency (english).wav' #name of the anechoic signal file
#fs, signal = wav.read(filename)
#signal = signal / 32767 # 2**15 - 1

import sounddevice as sd
import soundfile as sf
# Extract data and sampling rate from file
#with open('Frequency (english).wav', 'rb') as f:
#    data, samplerate = sf.read(f)

data_signal, fs = sf.read(filename) #this returns "data_signal", which is the 
#audiodata (one_dimentional array) of the anechoic signal. It returns also the
#"fs" sample frequency of the signal
 
#sd.play(data, fs) #this line allows to listen to the anechoic signal as it is.
#status = sd.wait()  #Wait until file is done playing

#Import the energy decay curve
edc = np.load('C:/Users/20225533/Diffusion/Auralization/w_rec_off.npy') #energy decay curve taken from the results of the diffusion equation model
edc_deriv = np.load('C:/Users/20225533/Diffusion/Auralization/w_rec_off_deriv.npy') #energy decay curve differentiated (or also impulse response of the room)
t_off = np.load('C:/Users/20225533/Diffusion/Auralization/t_off.npy') #decay time of the energy decay curve
t_off = t_off - t_off[0] #removing the t_off[0] to make the vector start from zero.

###############################################################################
################# METHOD 1 (DO NOT CONSIDER)
###############################################################################
#Frissen et al. applies energy decays to a normally distributed random number sequence.
# Is this correct??? Do I need to multiply the energy decay or the derivation of the energy decay?
#I need the pressure to multiply with a random numbers because the data from the anechoic file is pressure
 
#Random noise creation
# random_array = np.random.rand(1, len(edc_deriv))*2 - 1 #random noise vector with numbers between -1 and 1
# random_array = sum(random_array) #this line of code is used for passing from a row vector to a column vector

#Create impulse response
# imp_rand = edc_deriv*random_array

#Play the impulse response
#sd.play(imp_rand, fs)

# #ip = sum(imp_rand)
# ip = imp_rand/max(abs(imp_rand))

# plt.plot(t_off,ip)

# #Signal out
# from scipy.signal import lfilter
# audioOut = lfilter(ip, 1, data) ###Check what this thing does???? 
# audioOut = audioOut/max(audioOut)
#sd.play(audioOut, fs)

#TO DO
#check fft(ip) and fft(data) and multiply it then ifft(multiplication)
#make sure to zero pad both signals to the len(ip)+len(data)-1, do it before fft.

###############################################################################
################# METHOD 2
###############################################################################

#Random noise creation
# random_array = np.random.rand(1, len(edc_deriv))*2 - 1 #random noise vector with numbers between -1 and 1
# random_array = sum(random_array) #this line of code is used for passing from a row vector to a column vector

#Multiplication with random noise
# imp_rand = edc_deriv*random_array 

# #Play the impulse response
# #sd.play(imp_rand, fs)

# #Convolution of the impulse_rand with the anechoic signal
# st = np.arange(0,(len(data_signal))/fs,1/fs) #Time vector of the speech signal
# ht = np.arange(0,(len(imp_rand))/fs,1/fs)  #Time vector of the room impulse response

# #Create impulse response
# sh_conv = np.convolve(imp_rand,data_signal) #convolution of the impulse response with the anechoic signal
# sh_conv = sh_conv/max(abs(sh_conv)) #normalized to the maximum value of the convolved signal

# t_conv = np.arange(0,(len(sh_conv))/fs,1/fs) #Time vector of the convolved signal

# plt.plot(st,data_signal) #plot the anechoic signal
 
# plt.plot(ht,imp_rand) #plot the impulse response

# plt.plot(t_conv,sh_conv) #plot the convolved signal

#Play the convolved signal
#sd.play(sh_conv, fs)


###############################################################################
################# METHOD 3 (DO NOT CONSIDER)
###############################################################################
#Random noise creation
# random_array = np.random.rand(1, len(edc_deriv))*2 - 1 #random noise vector with numbers between -1 and 1
# random_array = sum(random_array) #this line of code is used for passing from a row vector to a column vector

#Multiplication with random noise
# imp_rand = edc_deriv*random_array 

#Luizard2019: ifft of the decaying noise #Wouter is checking it!
# imp_resp = np.fft.ifft(imp_rand)

# plt.plot(t_off, imp_resp.real, "b-", t_off, imp_resp.imag, "r--")# label='real')


###############################################################################
################# METHOD 4
###############################################################################
#Use a Gaussian instead of a random noise

#Gaussian distribution
mean = 0
sigmax_gau = 1 #??
gaussian_noise = (np.random.normal(mean, sigmax_gau, len(edc_deriv)))*2 - 1

#Multiplication with random noise
imp_rand = edc_deriv*gaussian_noise 

#Play the impulse response
#sd.play(imp_rand, fs)

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

# ###############################################################################
# ###############################################################################
# ######### CHECKS ####################################
#Backward integration of the imp_rand to see if the RT matches the one calculated from the diffusion equation model
R=np.zeros(len(t_off))

ir_squared = abs(imp_rand**2)
Ro=sum(ir_squared)
for q in range(1,len(t_off)):
    R[q]=sum(ir_squared[q:len(t_off)])

sch_db = 10.0 * np.log10(R / max(R))

sch_db = np.delete(sch_db,0)
sch_db = np.delete(sch_db,-1)

#Reverberation time calculation
init = -5.0 #because I want the T30, I need to start at -5
end = -35.0 #because I want the T30, I need to finish at -35
factor = 2.0 #factor of 2 since I need the T30

#rt_decay = sch_db #decay to be used to calculate the RT

#Linear regression
idxL1 = np.where(sch_db <= init)[0][0] #index at which the rtdecay is equal to -5
idxL2 = np.where(sch_db <= end)[0][0] #index at which the rtdecay is equal to -35
   
timeL1 = t_off[idxL1] #index at which the time vector is equal to the idxL1
timeL2 = t_off[idxL2] #index at which the time vector is equal to the idxL2

# Classical Approach (T30 = 2*(t[-35dB]-t[-5dB])
RTCalc = factor*(timeL2 - timeL1)

# Linregress approach
from scipy import stats
slope,intercept = stats.linregress(t_off[idxL1:idxL2],sch_db[idxL1:idxL2])[0:2] #calculating the slope and the interception of the line connecting the two points
db_regress_init = (init - intercept) / slope #dB initial
db_regress_end = (end - intercept) / slope #dB End
t60I = factor * (db_regress_end - db_regress_init) #t60 according to linregress approach

# Poly-based Approach y = Ax + B
CoefAlpha = np.polyfit(t_off[idxL1:idxL2], sch_db[idxL1:idxL2], 1) ##calculating the slope and the interception of the line connecting the two points
t60 = (-60/CoefAlpha[0]) #t60 according to polyfit approach

y_axis = (slope*t_off + intercept) + slope






