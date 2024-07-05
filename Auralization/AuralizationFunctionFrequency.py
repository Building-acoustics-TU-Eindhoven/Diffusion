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
from scipy import fft
from scipy.signal import butter, sosfilt, sosfreqz

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
#IMPORT DATA FROM DIFFUSION EQUATION CODE: 
#ENERGY DECAY CURVES/PRESSURE CURVES, dt, frequencies etc...
###############################################################################
#Import the energy decay curve
dt_sim = np.load('C:/Users/20225533/Diffusion/Auralization/dt.npy')
original_fs = 1/dt_sim
t_off = np.load('C:/Users/20225533/Diffusion/Auralization/t_off.npy') #decay time of the energy decay curve
t_off = t_off - t_off[0] #removing the t_off[0] to make the vector start from zero.
#edc_band = np.load('C:/Users/20225533/Diffusion/Auralization/w_rec_off_band.npy') #energy decay in terms of energy density
#edc_band = signal.resample(edc_band, int(len(edc_band) * fs / original_fs))

######
#ENERGY
######
edc_deriv_band = np.load('C:/Users/20225533/Diffusion/Auralization/w_rec_off_deriv_band.npy') #energy decay curve in terms of energy density differentiated
#num_samples = int(edc_deriv_band.shape[1] * fs / original_fs)
#edc_deriv_band_resampled = np.zeros((edc_deriv_band.shape[0], num_samples))

#for i in range(edc_deriv_band.shape[0]):
#    edc_deriv_band_resampled[i, :] = signal.resample_poly(edc_deriv_band[i, :], up=int(fs), down=int(original_fs))

#t_off = np.load('C:/Users/20225533/Diffusion/Auralization/t_off.npy') #decay time of the energy decay curve
#plt.plot(t_off,edc_deriv_band[0])

#t_off = t_off - t_off[0] #removing the t_off[0] to make the vector start from zero.
#t_off = np.linspace(0, len(t_off), num_samples)
#plt.plot(t_off,edc_deriv_band_resampled[4])


######
#PRESSURE
######
p_rec_off_deriv_band = np.load('C:/Users/20225533/Diffusion/Auralization/p_rec_off_deriv_band.npy') #energy decay curve in terms of pressure differentiated
p_rec_off_deriv_band_n4 = p_rec_off_deriv_band[4]

#Import the frequency from the main FVM calculation
center_freq = np.load('C:/Users/20225533/Diffusion/Auralization/center_freq.npy')
center_freq = center_freq.astype(np.int32)
nBands = len(center_freq) #np.load('C:/Users/20225533/Diffusion/Auralization/nBands.npy')

#%%
###############################################################################
#RESAMPLING PRESSURE ENVELOPE
###############################################################################
num_samples = int(p_rec_off_deriv_band.shape[1] * fs / original_fs)
p_rec_off_deriv_band_resampled = np.zeros((p_rec_off_deriv_band.shape[0], num_samples))

for i in range(p_rec_off_deriv_band.shape[0]):

    p_rec_off_deriv_band_resampled[i, :] = signal.resample_poly(p_rec_off_deriv_band[i, :], up=int(fs), down=int(original_fs))
    

#Clip negative values to zero
p_rec_off_deriv_band_resampled = np.clip(p_rec_off_deriv_band_resampled, a_min=0, a_max=None)

t_off_resampled = np.linspace(0, t_off[-1], num_samples)

plt.plot(t_off,p_rec_off_deriv_band[4])
plt.plot(t_off_resampled,p_rec_off_deriv_band_resampled[4])

#%%
###############################################################################
#SQUARE-ROOT of ENVELOPE
###############################################################################
#From the envelope of the impulse response, we need to get the impulse response
square_root = np.sqrt(p_rec_off_deriv_band_resampled) #this gives the impulse response at each frequency
plt.plot(t_off_resampled,square_root[4])

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
# noise = np.random.rand(1, p_rec_off_deriv_band_resampled.shape[1])*2 - 1 #random noise vector with unifrom distribution and with numbers between -1 and 1
# noise = sum(noise) #this line of code is used for passing from a row vector to a column vector
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
# noise = np.random.normal(0,1,p_rec_off_deriv_band_resampled.shape[1]) #This is how Gerd does it -> This is not from -1 to 1 but from -4 to 4
# mean_value = np.mean(noise)
# difference_squared = (noise - mean_value)**2
# variance = np.mean(difference_squared)
#Variance here is 0.997 so this could be good -> the variance is 1 but the noise is not between -1 and 1

#FOURTH ATTEMPT of noise creation
noise = np.random.rand(1, p_rec_off_deriv_band_resampled.shape[1])*2*(np.sqrt(3)) - (np.sqrt(3)) #random noise vector with unifrom distribution and with numbers between -1 and 1
noise = sum(noise) #this line of code is used for passing from a row vector to a column vector
mean_value = np.mean(noise)
difference_squared = (noise - mean_value)**2
variance = np.mean(difference_squared) 
#In this way the variance is 0.33 and not 1

#plt.plot(t_off_resampled,noise)

#FIFTH ATTEMPT of noise creation: after speaking with Wouter; the noise needs to be normalized and the variance needs to be one
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

######
#BUTTER FILTER
######
Nyquist_freq = int(fs/2) 
filter_order = 8 #number of biquad sections of the desired system
nth_octave = 1  # e.g., 3 for third-octave

# Create filter
filter_tot = []
for fc in center_freq:
    # Calculate low and high cutoff frequencies for each band
    lowcut = fc / (2 ** (1 / (2 * nth_octave)))
    highcut = fc * (2 ** (1 / (2 * nth_octave)))
    
    # Normalize the cutoff frequencies by the Nyquist frequency
    low = lowcut / Nyquist_freq
    high = highcut / Nyquist_freq
    
    # Design Butterworth bandpass filter
    butter_band = butter(filter_order, [low, high], btype='band', output='sos') # butter_band contains the second-order sections representation of the Butterworth filter
    
    # Append filter coefficients to filter tot
    filter_tot.append(butter_band)

# Plot frequency responses
plt.figure(figsize=(12, 8))
for band in range(len(filter_tot)):
    print(band)
    # Compute the frequency response of each filter
    w, h = sosfreqz(filter_tot[band], worN=2000, fs=fs)
    # w is the array of frequencies at which the response is computed
    # h is the frequency response of the filter
    
    # Plot the magnitude response in decibels
    plt.semilogx(w, 20 * np.log10(abs(h)), label=f'{center_freq[band]}')

plt.title('Frequency response of octave band filters')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Gain [dB]')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.axhline(0, color='black', linewidth=0.5)
plt.legend(loc='best')
plt.show()

######
#PREVIOUS FILTER
######
# Nyquist_freq = int(fs/2) 
# Nb = 7 #number of biquad sections of the desired system

# h_all=np.zeros(shape=(nBands,Nyquist_freq)) #matrix to store the band-pass filter responses for each band

# #Create octave band filters
# h_low=np.zeros(Nyquist_freq)
# h_low[0] = 1
# h_low_all=np.zeros(shape=(nBands,Nyquist_freq))
# for fi in reversed(range(nBands)):
#     wn_low = np.power(2,1./2)*center_freq[fi]/Nyquist_freq #normalized critical frequency for the low-pass filter
#     b_low, a_low  = signal.butter(Nb,wn_low,'low') #coefficients of the Butterworth low-pass filter
#     h_low = signal.lfilter(b_low,a_low,h_low) #applies the filter to the impulse response -> creates the filter impulse response
#     h_low_all[fi,:]=h_low

# h_high=np.zeros(Nyquist_freq)
# h_high[0] = 1
# h_high_all=np.zeros(shape=(nBands,Nyquist_freq))

# for fi in range(nBands):
#     wn_high =np.power(2,-1./2)*center_freq[fi]/Nyquist_freq #normalized critical frequency for the high-pass filter
#     b_high, a_high = signal.butter(Nb, wn_high,'high') #coefficients of the Butterworth high-pass filter
#     h_high = signal.lfilter(b_high,a_high,h_high) #-> creates the filter impulse response
#     h_high_all[fi,:]=h_high
    
#     h=np.convolve(h_high_all[fi,:],h_low_all[fi,:])[0:Nyquist_freq] #band-pass filter is the result of convolving the corresponding low-pass and high-pass filters. It allows frequencies within a specific range (octave band) to pass while attenuating frequencies outside this range
#     h_all[fi,:]=h
    
#     # Frequency response
#     H = np.fft.fft(h)
#     H = H[:Nyquist_freq]
    
#     # Plot the filter response for this band
#     plt.figure(figsize=(10, 6))
#     freqs = np.linspace(0, Nyquist_freq, len(H))
#     plt.plot(freqs, 20 * np.log10(abs(H)))
#     plt.title(f'Band-pass Filter {fi + 1} ({center_freq[fi]} Hz)')
#     plt.xlabel('Frequency (Hz)')
#     plt.ylabel('Magnitude (dB)')
#     plt.grid()
#     plt.show()

# # Plot all filter responses combined
# plt.figure(figsize=(10, 6))
# for fi in range(nBands):
#     H = np.fft.fft(h_all[fi, :])
#     H = H[:Nyquist_freq]
#     freqs = np.linspace(0, Nyquist_freq, len(H))
#     plt.plot(freqs, 20 * np.log10(abs(H)), label=f'{center_freq[fi]} Hz')
# plt.title('All Band-pass Filters Combined')
# plt.xlabel('Frequency (Hz)')
# plt.ylabel('Magnitude (dB)')
# plt.legend()
# plt.grid()
# plt.show()

#%%
###############################################################################
#CREATION OF FILTERED RANDOM NOISE
###############################################################################

#TIME DOMAIN OF THE FILTERED RANDOM NOISE: for each band the sosfilt creates a time domain convolution of the noise with the filter
filt_noise_band = [sosfilt(band, noise) for band in filter_tot] #this is in the time domain because the sosfilt gives the time domain

for fi in range(nBands):
    plt.plot(t_off_resampled, 20 * np.log10(abs(filt_noise_band[fi])), label=f'{center_freq[fi]} Hz')

#FREQUENCY DOMAIN OF THE FILTERED RANDOM NOISE: Frequency response of the filtered random noise
filt_noise_band_freq = np.fft.fft(filt_noise_band)

# fv is the freqeuncy vector for the x axis
nSamples = len(noise)
fv = np.arange(nSamples) * (fs/nSamples) #This can also be written with linspace as #np.linspace((0, (nSamples-1)))*fs/nSamples;

#Plot the frequency response
plt.figure(figsize=(12, 8))
plt.title('Frequency response of octave band filtered random noise')
for fi in range(nBands):
    plt.subplot(nBands, 1, fi+1)
    plt.semilogx(fv, 20 * np.log10(abs(filt_noise_band_freq[fi, :])), label=f'{center_freq[fi]} Hz')

    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Gain [dB]')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    #plt.axhline(0, color='black', linewidth=0.1)
    plt.legend(loc='best')
    plt.show()

#Make the filt_noise_band list into an array
filt_noise_band = np.array(filt_noise_band)

# # Determine the length of the convolved result between filter and random noise
# conv_length = len(np.convolve(h_all[0, :], noise))

# #Convolution of random noise with filter
# filt_noise_band = np.empty((nBands, conv_length), dtype=float)
# for fi in range(nBands):
#     filt_noise=np.convolve(h_all[fi,:],noise)
#     filt_noise_band[fi, :] = filt_noise

#%%
###############################################################################
#MULTIPLICATION OF SQUARE-ROOT with FILTERED RANDOM NOISE
###############################################################################
#Multiplication of SQUARE-ROOT of envelope with random noise only (UNFILTERED)
imp_unfilt_band = []
for fi in range(nBands):
    imp_unfilt = square_root[fi,:] * noise
    imp_unfilt_band.append(imp_unfilt)

plt.figure(figsize=(12, 8))
plt.title('Time domain of UNfiltered impulse response per frequency band')
for fi in range(nBands):
    plt.subplot(nBands, 1, fi+1)
    plt.plot(t_off_resampled, imp_unfilt_band[fi], label=f'{center_freq[fi]} Hz')
    
    plt.xlabel('Time [s]')
    plt.ylabel('Magnitude [dB]')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    #plt.axhline(0, color='black', linewidth=0.1)
    plt.legend(loc='best')
    plt.show()

#Padding the square-root to the same length as the filtered random noise
pad_length = filt_noise_band.shape[1]-p_rec_off_deriv_band_resampled.shape[1]
square_root_padded = np.pad(square_root, ((0,0),(0,pad_length)) ,mode='constant' )
t_off_padded = np.pad(t_off_resampled, ((0,pad_length)) ,mode='constant' )

#Multiplication of SQUARE-ROOT of envelope with filtered random noise (FILTERED)
imp_filt_band = []
plt.figure(figsize=(12, 8))
plt.title('Time domain of filtered impulse response per frequency band')
for fi in range(nBands):
    plt.subplot(nBands, 1, fi+1)
    imp_filt = square_root_padded[fi,:]*filt_noise_band[fi,:]
    imp_filt_band.append(imp_filt)
    
    plt.plot(t_off_padded, imp_filt_band[fi], label=f'{center_freq[fi]} Hz')
    
    plt.xlabel('Time [s]')
    plt.ylabel('Magnitude [dB]')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    #plt.axhline(0, color='black', linewidth=0.1)
    plt.legend(loc='best')
    plt.show()

#%%
###############################################################################
#ALL FREQUENCY IMPULSE RESPONSE
###############################################################################
#Sum of the bands in the time domain
imp_tot = [sum(imp_filt_band[i][j] for i in range(len(imp_filt_band))) for j in range(len(imp_filt_band[0]))]
imp_tot = np.array(imp_tot, dtype=float)

# Plot impulse response in the time domain
plt.figure(figsize=(12, 8))
plt.plot(t_off_padded, imp_tot)
plt.xlabel('Time [s]')
plt.ylabel('Magnitude [dB]')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
#plt.axhline(0, color='black', linewidth=0.1)
#plt.legend(loc='best')
plt.show()


#Frequency domain of the total impulse response
freq_spectrum = 20*np.log10(abs(np.fft.fft(imp_tot)))

# Plot impulse response in the frequency domain
plt.figure(figsize=(12, 8))
plt.semilogx(fv, freq_spectrum)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Magnitude [dB]')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
#plt.axhline(0, color='black', linewidth=0.1)
#plt.legend(loc='best')
plt.show()

# #%%
# ###############################################################################
# #ADDING DIRECT SOUND & SHIFT EVERYTHING TO THE DIRECT SOUND ARRIVAL TIME STEP
# ###############################################################################
# W = 0.01 #the power is in Watts
# dist_sr = 1.5
# rho = 1.21
# c0 = 343

# #POSSIBLE SOLUTION FOR HAVING A DIFFERENT POWER PER BANDWIDTH??????? i am not sure
# #???
# # # Calculate bandwidths for each octave band
# # bandwidths = center_freq * (2**(1/2) - 2**(-1/2))

# # # Normalize power per band by bandwidth
# # power_per_band = W * (bandwidths / np.sum(bandwidths))

# # # Calculate direct component for each band
# # press_dir_sound = np.sqrt(power_per_band / (4 * np.pi * dist_sr**2) * rho * c0)
# #???


# press_dir_sound = np.sqrt((W/(4*np.pi*dist_sr**2))*rho*c0)

# time_dir_sound = dist_sr/c0
# time_dir_sound_step = int(time_dir_sound/dt_sim)
# time_dir_sound_step_resampled = int(time_dir_sound_step * fs/ original_fs)

# #Shift everythingto the direct sound arrival
# zero_array = np.zeros([time_dir_sound_step_resampled]) #array of zeros in front of the impulse response
# imp_tot = np.hstack([zero_array,imp_tot])

# #Add the direct component to the impulse response at the right time step
# imp_tot[time_dir_sound_step_resampled] += press_dir_sound

# pad_length_direct = imp_tot.shape[0]-t_off_padded.shape[0]
# t_off_padded_direct = np.pad(t_off_padded, ((0,pad_length_direct)) ,mode='constant' )

# plt.plot(t_off_padded_direct,imp_tot)



#%%
###############################################################################
#ADDING DIRECT SOUND AND FILTER IT & SHIFT EVERYTHING TO THE DIRECT SOUND ARRIVAL TIME STEP
###############################################################################
# W = 0.01 #the power is in Watts
# dist_sr = 1.5
# rho = 1.21
# c0 = 343

# #POSSIBLE SOLUTION FOR HAVING A DIFFERENT POWER PER BANDWIDTH??????? i am not sure
# #???
# # # Calculate bandwidths for each octave band
# # bandwidths = center_freq * (2**(1/2) - 2**(-1/2))

# # # Normalize power per band by bandwidth
# # power_per_band = W * (bandwidths / np.sum(bandwidths))

# # # Calculate direct component for each band
# # press_dir_sound = np.sqrt(power_per_band / (4 * np.pi * dist_sr**2) * rho * c0)
# #???


# press_dir_sound = np.sqrt((W/(4*np.pi*dist_sr**2))*rho*c0)

# #Multiplication of press_dir_sound with filtered random noise (FILTERED)
# press_dir_sound_band = []
# for fi in range(nBands):
#     dir_filt = press_dir_sound*h_all[fi,:]
#     press_dir_sound_band.append(dir_filt)


# time_dir_sound = dist_sr/c0
# time_dir_sound_step = int(time_dir_sound/dt_sim)
# time_dir_sound_step_resampled = int(time_dir_sound_step * fs/ original_fs)

# #Shift everything to the direct sound arrival
# zero_array = np.zeros([time_dir_sound_step_resampled]) #array of zeros in front of the impulse response
# imp_tot = np.hstack([zero_array,imp_tot])

# #Add the direct component to the impulse response at the right time step
# imp_tot[time_dir_sound_step_resampled] += press_dir_sound

# pad_length_direct = imp_tot.shape[0]-t_off_padded.shape[0]
# t_off_padded_direct = np.pad(t_off_padded, ((0,pad_length_direct)) ,mode='constant' )

# plt.plot(t_off_padded_direct,imp_tot)


#%%
###############################################################################
#ALL FREQUENCY IMPULSE RESPONSE
###############################################################################
#Sum of the bands


# imp_filt_band_direct = [imp_filt_band[i]+press_dir_sound_band[i] for i in range(len(imp_filt_band))]


# imp_tot = [sum(imp_filt_band_direct[i][j] for i in range(len(imp_filt_band_direct))) for j in range(len(imp_filt_band_direct[0]))]
# imp_tot = np.array(imp_tot, dtype=float)

# plt.plot(t_off_padded,imp_tot)


# #Frequency spectrum
# freq_spectrum = 20*np.log10(abs(np.fft.fft(imp_tot)))









#%%
###############################################################################
#CREATE IMPULSE RESPONSE FILE
###############################################################################
#Create a file wav for impulse response
imp_resp_norm = np.int16(imp_tot / np.max(np.abs(imp_tot)) * 32767) 

scipy.io.wavfile.write("imp_resp.wav", fs, imp_resp_norm)

#Play the impulse response
#sd.play(imp_tot, fs)


#%%
###############################################################################
#TRIAL IMPULSE RESPONSE ONLY FOR ONE BAND
###############################################################################
#Create a file wav for impulse response

#125Hz
imp_resp_band0 = np.int16(imp_filt_band[0] / np.max(np.abs(imp_filt_band[0])) * 32767) 

scipy.io.wavfile.write("imp_resp_band0.wav", fs, imp_resp_band0)

#250Hz
imp_resp_band1 = np.int16(imp_filt_band[1] / np.max(np.abs(imp_filt_band[1])) * 32767) 

scipy.io.wavfile.write("imp_resp_band1.wav", fs, imp_resp_band1)

#500Hz
imp_resp_band2 = np.int16(imp_filt_band[2] / np.max(np.abs(imp_filt_band[2])) * 32767) 

scipy.io.wavfile.write("imp_resp_band2.wav", fs, imp_resp_band2)

#1000Hz
imp_resp_band3 = np.int16(imp_filt_band[3] / np.max(np.abs(imp_filt_band[3])) * 32767) 

scipy.io.wavfile.write("imp_resp_band3.wav", fs, imp_resp_band3)

#2000Hz
imp_resp_band4 = np.int16(imp_filt_band[4] / np.max(np.abs(imp_filt_band[4])) * 32767) 

scipy.io.wavfile.write("imp_resp_band4.wav", fs, imp_resp_band4)

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

#Create a file wav for auralization
#scipy.io.wavfile.write("auralization.wav", fs, sh_conv)

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