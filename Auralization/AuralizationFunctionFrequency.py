# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 14:33:58 2024

@author: 20225533
"""

#%%
###############################################################################
#IMPORT PACKAGES
###############################################################################

import numpy as np
#import matplotlib
import matplotlib.pyplot as plt
from FunctionRT import *
import sounddevice as sd
import soundfile as sf
import scipy
from scipy.io import wavfile
from scipy import signal

#%%
###############################################################################
#IMPORT ANECHOIC SIGNAL
###############################################################################
#Import anechoic signal
filename = 'C:/Users/20225533/Diffusion/Auralization/Frequency(english).wav' #name of the anechoic signal file

# Extract data and sampling rate from file
data_signal, fs = sf.read(filename) #this returns "data_signal", which is the 
#audiodata (one_dimentional array) of the anechoic signal. It returns also the
#"fs" sample frequency of the signal
 
#sd.play(data, fs) #this line allows to listen to the anechoic signal as it is.
#status = sd.wait()  #Wait until file is done playing

#%%
###############################################################################
#IMPORT ENERGY DECAY CURVES
###############################################################################
#Import the energy decay curve
edc_band = np.load('C:/Users/20225533/Diffusion/Auralization/w_rec_band.npy')
#edc = np.load('C:/Users/20225533/Diffusion/Auralization/w_rec_off.npy') #energy decay curve taken from the results of the diffusion equation model
t = np.load('C:/Users/20225533/Diffusion/Auralization/t.npy') #time steps array

edc_deriv_band = np.load('C:/Users/20225533/Diffusion/Auralization/w_rec_band_deriv.npy') #energy decay curve differentiated (or also impulse response of the room)
t_off = np.load('C:/Users/20225533/Diffusion/Auralization/t_off.npy') #decay time of the energy decay curve
t_off = t_off - t_off[0] #removing the t_off[0] to make the vector start from zero.

#%%
###############################################################################
#SQUARE-ROOT of ENVELOPE
###############################################################################
#From the envelope of the impulse response, we need to get the impulse response
square_root = np.sqrt(edc_deriv_band) #this gives the impulse response

#%%
###############################################################################
#CREATION OF RANDOM NOISE
###############################################################################
#Random noise creation
random_array = np.random.rand(1, edc_deriv_band.shape[1])*2 - 1 #random noise vector with numbers between -1 and 1
random_array = sum(random_array) #this line of code is used for passing from a row vector to a column vector


#%%
###############################################################################
#CREATION OF FILTER
###############################################################################

#Import the frequency from the main FVM calculation
center_freq = np.load('C:/Users/20225533/Diffusion/Auralization/center_freq.npy')
center_freq = center_freq.astype(np.int32)
nBands = len(center_freq) #np.load('C:/Users/20225533/Diffusion/Auralization/nBands.npy')

Nyquist_freq = int(fs/2) 
Nb = 5 #number of biquad sections of the desired system

h_all=np.zeros(shape=(nBands,Nyquist_freq)) #matrix to store the band-pass filter responses for each band

#Create octave band filters
h_low=np.zeros(Nyquist_freq)
h_low[0] = 1
h_low_all=np.zeros(shape=(nBands,Nyquist_freq))
for fi in reversed(range(nBands)):
    wn_low = np.power(2,1./2)*center_freq[fi]/Nyquist_freq #normalized critical frequency for the low-pass filter
    b_low, a_low  = signal.butter(Nb,wn_low,'low') #coefficients of the Butterworth low-pass filter
    h_low = signal.lfilter(b_low,a_low,h_low) #applies the filter to the impulse response -> creates the filter impulse response
    h_low_all[fi,:]=h_low

h_high=np.zeros(Nyquist_freq)
h_high[0] = 1
h_high_all=np.zeros(shape=(nBands,Nyquist_freq))

for fi in range(nBands):
    wn_high =np.power(2,-1./2)*center_freq[fi]/Nyquist_freq #normalized critical frequency for the high-pass filter
    b_high, a_high = signal.butter(Nb, wn_high,'high') #coefficients of the Butterworth high-pass filter
    h_high = signal.lfilter(b_high,a_high,h_high)
    h_high_all[fi,:]=h_high
    
    h=np.convolve(h_high_all[fi,:],h_low_all[fi,:]) #band-pass filter is the result of convolving the corresponding low-pass and high-pass filters. It allows frequencies within a specific range (octave band) to pass while attenuating frequencies outside this range
    h_all[fi,:]=h[0:Nyquist_freq]


# Plot the frequency responses
freqs = np.linspace(0, Nyquist_freq, Nyquist_freq)
plt.figure(figsize=(14, 10))

for i in range(nBands):
    plt.subplot(nBands, 1, i + 1)
    plt.plot(freqs, 20 * np.log10(np.abs(h_all[i, :])))
    plt.title(f'Band-pass Filter for Center Frequency {center_freq[i]} Hz')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Magnitude (dB)')
    plt.grid(True)

plt.tight_layout()
plt.show()

# Frequency response
w, h_low = signal.freqz(b_low, a_low, worN=fs)
w, h_high = signal.freqz(b_high, a_high, worN=fs)

# Plotting
plt.figure(figsize=(14, 6))

# Low-pass filter response
plt.subplot(2, 1, 1)
plt.plot(w / np.pi * Nyquist_freq, 20 * np.log10(np.abs(h_low)), 'b')
plt.title('Low-Pass Filter Frequency Response')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Magnitude (dB)')
plt.grid()

# High-pass filter response
plt.subplot(2, 1, 2)
plt.plot(w / np.pi * Nyquist_freq, 20 * np.log10(np.abs(h_high)), 'r')
plt.title('High-Pass Filter Frequency Response')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Magnitude (dB)')
plt.grid()

plt.tight_layout()
plt.show()

#%%
###############################################################################
#CREATION OF FILTERED RANDOM NOISE
###############################################################################

# Determine the length of the convolved result between filter and random noise
conv_length = len(np.convolve(h_all[0, :], random_array))

#Convolution of random noise with filter


filt_noise_band = np.empty((nBands, conv_length), dtype=float)
for fi in range(nBands):
    filt_noise=np.convolve(h_all[fi,:],random_array)
    filt_noise_band[fi, :] = filt_noise



#%%
###############################################################################
#MULTIPLICATION OF SQUARE-ROOT with FILTERED RANDOM NOISE
###############################################################################
#Multiplication of SQUARE-ROOT of envelope with random noise only (UNFILTERED)
imp_unfilt_band = []
for fi in range(nBands):
    imp_unfilt = square_root[fi,:] * random_array
    imp_unfilt_band.append(imp_unfilt)


#Padding the square-root to the same length as the filtered random noise
pad_length = filt_noise_band.shape[1]-edc_deriv_band.shape[1]
square_root_padded = np.pad(square_root, ((0,0),(0,pad_length)) ,mode='constant' )

#Multiplication of SQUARE-ROOT of envelope with filtered random noise (FILTERED)
imp_filt_band = []
for fi in range(nBands):
    imp_filt = square_root_padded[fi,:]*filt_noise_band[fi,:]
    #imp_filt=np.convolve(h_all[fi,:],square_root[fi,:]*random_array)
    imp_filt_band.append(imp_filt)


#Sum of the bands
imp_tot = [sum(imp_filt_band[i][j] for i in range(len(imp_filt_band))) for j in range(len(imp_filt_band[0]))]
imp_tot = np.array(imp_tot, dtype=float)


#Create a file wav for impulse response
scipy.io.wavfile.write("imp_resp.wav", fs, imp_tot)

#Play the impulse response
sd.play(imp_tot, fs)

#%%
###############################################################################
#CONVOLUTION FOR AURALIZATION
###############################################################################

#Convolution of the impulse_rand with the anechoic signal
st = np.arange(0,(len(data_signal))/fs,1/fs) #Time vector of the speech signal
ht = np.arange(0,(len(imp_tot))/fs,1/fs)  #Time vector of the room impulse response

#Create impulse response
sh_conv = np.convolve(imp_tot,data_signal) #convolution of the impulse response with the anechoic signal
sh_conv = sh_conv/max(abs(sh_conv)) #normalized to the maximum value of the convolved signal

t_conv = np.arange(0,(len(sh_conv))/fs,1/fs) #Time vector of the convolved signal

plt.plot(st,data_signal) #plot the anechoic signal
 
plt.plot(ht,imp_tot) #plot the impulse response

plt.plot(t_conv,sh_conv) #plot the convolved signal

#Play the convolved signal
sd.play(sh_conv, fs)


#%%
###############################################################################
#CHECKS
###############################################################################
#Expected energy ratio
n = 1 #for octave band
dt = 1/fs #???? is this correct?
Ci_exp = []
for fi in range(len(center_freq)):
    Ci = center_freq[fi]*((2**(1/n)-1)/(2**(1/2*n)))*2*dt
    Ci_exp.append(Ci)
    

#Predicted energy ratio
n = 1 #for octave band
dt = 1/fs #???? is this correct?
Ci_pred = []
for fi in range(len(center_freq)):
    Ci = sum(abs(imp_filt_band[fi])**2*dt)/(sum(abs(imp_unfilt_band[fi])**2*dt))
    Ci_pred.append(Ci)    