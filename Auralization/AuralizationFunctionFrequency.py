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
dt_sim = np.load('C:/Users/20225533/Diffusion/Auralization/dt.npy')
original_fs = 1/dt_sim
edc_band = np.load('C:/Users/20225533/Diffusion/Auralization/w_rec_off_band.npy')
#edc_band = signal.resample(edc_band, int(len(edc_band) * fs / original_fs))

edc_deriv_band = np.load('C:/Users/20225533/Diffusion/Auralization/w_rec_off_deriv_band.npy') #energy decay curve differentiated (or also impulse response of the room)
#edc_deriv_band = signal.resample(edc_deriv_band, int(edc_deriv_band.shape[1] * fs / original_fs)))
num_samples = int(edc_deriv_band.shape[1] * fs / original_fs)
edc_deriv_band_resampled = np.zeros((edc_deriv_band.shape[0], num_samples))

for i in range(edc_deriv_band.shape[0]):
    edc_deriv_band_resampled[i, :] = signal.resample(edc_deriv_band[i, :], num_samples)

t_off = np.load('C:/Users/20225533/Diffusion/Auralization/t_off.npy') #decay time of the energy decay curve
t_off = t_off - t_off[0] #removing the t_off[0] to make the vector start from zero.
t_off = np.linspace(0, len(t_off) / fs, len(t_off))

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
#noise = np.random.random(edc_deriv_band.shape[1]) #This is random from 0 to 1 but NOT uniform distribution
#noise1= np.random.normal(0,1,edc_deriv_band.shape[1]) #This is how Gerd does it -> This is not from -1 to 1 but from -4 to 4
#random = np.random.rand(1, edc_deriv_band.shape[1]) #random noise vector with unifrom distribution and with numbers between 0 and 1
#random = sum(random) #this line of code is used for passing from a row vector to a column vector

#FIRST ATTEMPT of noise creation
# noise = np.random.rand(1, edc_deriv_band.shape[1])*2 - 1 #random noise vector with unifrom distribution and with numbers between -1 and 1
# noise = sum(random_noise_array) #this line of code is used for passing from a row vector to a column vector
# mean_value = np.mean(noise)
# difference_squared = (noise - mean_value)**2
# variance = np.mean(difference_squared) 
#In this way the variance is 0.33 and not 1

#SECOND ATTEMPT of noise creation: after speaking with Wouter; the noise needs to be normal distribution and the variance needs to be one
# noise = np.random.random(edc_deriv_band.shape[1]) #This is random from 0 to 1 but NOT uniform distribution
# noise = noise*2 -1
# mean_value = np.mean(noise)
# difference_squared = (noise - mean_value)**2
# variance = np.mean(difference_squared)
#In this way the variance is 0.33 and not 1

#THIRD ATTEMPT of noise creation: after speaking with Wouter; the noise needs to be normalized and the variance needs to be one
noise = np.random.normal(0,1,edc_deriv_band.shape[1]) #This is how Gerd does it -> This is not from -1 to 1 but from -4 to 4
mean_value = np.mean(noise)
difference_squared = (noise - mean_value)**2
variance = np.mean(difference_squared)
#Variance here is 0.997 so this could be good -> the variance is 1 but the noise is not between -1 and 1

#FOURTH ATTEMPT of noise creation: after speaking with Wouter; the noise needs to be normalized and the variance needs to be one
# noise = np.random.rand(1, edc_deriv_band.shape[1]) #random noise vector with unifrom distribution and with numbers between 0 and 1
# noise = sum(noise) #this line of code is used for passing from a row vector to a column vector
# mean_value = np.mean(noise)
# difference_squared = (noise - mean_value)**2
# variance = np.mean(difference_squared)
#In this way the variance is 0.08 and not 1

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
    plt.xlabel('Time (s)')
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
conv_length = len(np.convolve(h_all[0, :], noise))

#Convolution of random noise with filter
filt_noise_band = np.empty((nBands, conv_length), dtype=float)
for fi in range(nBands):
    filt_noise=np.convolve(h_all[fi,:],noise)
    filt_noise_band[fi, :] = filt_noise

#%%
###############################################################################
#MULTIPLICATION OF SQUARE-ROOT with FILTERED RANDOM NOISE
###############################################################################
#Multiplication of SQUARE-ROOT of envelope with random noise only (UNFILTERED)
imp_unfilt_band = []
for fi in range(nBands):
    imp_unfilt = square_root[fi,:] * noise
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

#%%
###############################################################################
#ALL FREQUENCY IMPULSE RESPONSE
###############################################################################
#Sum of the bands
imp_tot = [sum(imp_filt_band[i][j] for i in range(len(imp_filt_band))) for j in range(len(imp_filt_band[0]))]
imp_tot = np.array(imp_tot, dtype=float)

#Create a file wav for impulse response
scipy.io.wavfile.write("imp_resp.wav", fs, imp_tot)

#Play the impulse response
#sd.play(imp_tot, fs)

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
#sd.play(sh_conv, fs)

#Create a file wav for impulse response
scipy.io.wavfile.write("auralization.wav", fs, sh_conv)

#%%
###############################################################################
#FROM FLOATING POINT FORMAT TO standard integer format such as 16-bit
###############################################################################
# Normalize the floating-point data to the range of int16
sh_conv_normalized = np.int16(sh_conv / np.max(np.abs(sh_conv)) * 32767) #32767 scales the normalized data to the range of 16-bit integers (-32768 to 32767).

# Write the normalized data to a WAV file
scipy.io.wavfile.write("auralization.wav", fs, sh_conv_normalized)

#%%
###############################################################################
#CHECKS
###############################################################################
#Expected energy ratio
n = 1 #for octave band
dt = 1/fs 
Ci_exp = []
for fi in range(len(center_freq)):
    Ci = center_freq[fi]*((2**(1/n)-1)/(2**(1/2*n)))*2*dt
    Ci_exp.append(Ci)
    

#Predicted energy ratio
n = 1 #for octave band
dt = 1/fs
Ci_pred = []
for fi in range(len(center_freq)):
    Ci = sum(abs(imp_filt_band[fi])**2*dt)/(sum(abs(imp_unfilt_band[fi])**2*dt))
    Ci_pred.append(Ci)    