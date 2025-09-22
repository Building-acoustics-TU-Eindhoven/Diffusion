# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 14:33:58 2024

@author: Ilaria Fichera
"""
#%%
###############################################################################
#IMPORT PACKAGES
###############################################################################

import numpy as np
#import matplotlib.pyplot as plt
# import sounddevice as sd
import soundfile as sf
import scipy
from scipy import signal
from scipy.signal import butter, sosfilt, sosfreqz
import pickle
import os


def load(filename):
    """
    Loading variables

    Parameters
    ----------
        filename : str
            Name of the file to load the variables from
    Returns
    -------
        variables : dict
            Variables from the FVM results
    """
    import pickle
    with open(filename, 'rb') as f:
        variables = pickle.load(f, encoding="latin1")
        globals().update(variables)
    
    return variables


def extract_anechoic(anechoic_signal):
    """
    Extracting information form the anechoic signal

    Parameters
    ----------
        anechoic_signal : wav
            wav file of the anechoic sound
    Returns
    -------
        data_signal : array of floats 
            Signal data of the anechoic file
        fs : float
            Sample frequency of the anechoic signal
    """
    # Extract data and sampling rate from file
    data_signal, fs = sf.read(anechoic_signal) #this returns "data_signal", which is the 
    #audiodata (one_dimentional array) of the anechoic signal. It returns also the
    #"fs" sample frequency of the signal

    return data_signal, fs

def resample_pressure(p_rec_off_deriv_band, fs, fs_sim, t_off_sim):
    """
    Resampling the pressure over time to the dimension of the anechoic signal

    Parameters
    ----------
        p_rec_off_deriv_band : list of arrays
            Pressure over time after the source is switched off per each frequency
        fs : float
            Sample frequency of the anechoic signal
        fs_sim : float
            Sample frequency of the FVM simulation
        t_off_sim : array of float
            Time array after the source is switched off
            
    Returns
    -------
        p_rec_off_deriv_band_resampled : list of arrays
            Pressure over time after the source is switched off per each frequency resampled
        t_off_resampled : array of float
            Time array after the source is switched off resampled
    """
    num_samples = int(p_rec_off_deriv_band.shape[1] * fs / fs_sim)
    p_rec_off_deriv_band_resampled = np.zeros((p_rec_off_deriv_band.shape[0], num_samples))
    
    for i in range(p_rec_off_deriv_band.shape[0]):
    
        p_rec_off_deriv_band_resampled[i, :] = signal.resample_poly(p_rec_off_deriv_band[i, :], up=int(fs), down=int(fs_sim))
        
    #Clip negative values to zero
    p_rec_off_deriv_band_resampled = np.clip(p_rec_off_deriv_band_resampled, a_min=0, a_max=None)
    
    t_off_resampled = np.linspace(0, t_off_sim[-1], num_samples)
    
    return p_rec_off_deriv_band_resampled, t_off_resampled


def square_root_fun(p_rec_off_deriv_band_resampled):
    """
    Resampling the pressure over time to the dimension of the anechoic signal

    Parameters
    ----------
        p_rec_off_deriv_band_resampled : list of arrays
            Pressure over time after the source is switched off per each frequency resampled
            
    Returns
    -------
        square_root : list of arrays
            Square rooted pressure over time after the source is switched off per each frequency resampled
    """
    #From the envelope of the impulse response, we need to get the impulse response
    square_root = np.sqrt(p_rec_off_deriv_band_resampled) #this gives the impulse response at each frequency
    return square_root


def create_noise(p_rec_off_deriv_band_resampled):
    """
    Creation of random noise

    Parameters
    ----------
        p_rec_off_deriv_band_resampled : list of arrays
            Pressure over time after the source is switched off per each frequency resampled
            
    Returns
    -------
        noise : array of float
            Random noise vector with unifrom distribution over time
    """
    noise = np.random.rand(1, p_rec_off_deriv_band_resampled.shape[1])*2*(np.sqrt(3)) - (np.sqrt(3)) #random noise vector with unifrom distribution and with numbers between -1 and 1
    noise = sum(noise) #this line of code is used for passing from a row vector to a column vector
    mean_value = np.mean(noise)
    difference_squared = (noise - mean_value)**2
    variance = np.mean(difference_squared) 
    return noise


def butter_filter(fs, center_freq):
    """
    Creation of Butterworth bandpass filter

    Parameters
    ----------
        fs : float
            Sample frequency of the anechoic signal
        center_freq : array of int
            Array of all the frequencies to calculate.
            
    Returns
    -------
        filter_tot : list of arrays
            Filter coefficient per each frequency
    """
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
    return filter_tot


def filtered_noise(noise, filter_tot, fs):
    """
    Calculation of filtered noise

    Parameters
    ----------
        noise : array of float
            Random noise vector with unifrom distribution over time
        filter_tot : list of arrays
            Filter coefficient per each frequency
        fs : float
            Sample frequency of the anechoic signal
        
    Returns
    -------
        filt_noise_band : array of floats
            Results of the time domain convolution of the noise with the filter per each frequency
        filt_noise_band_freq : array of complex
            Frequency response of the filtered random noise
        fv : array of float
            Frequency vector
    """
    
    #TIME DOMAIN OF THE FILTERED RANDOM NOISE: for each band the sosfilt creates a time domain convolution of the noise with the filter
    filt_noise_band = [sosfilt(band, noise) for band in filter_tot] #this is in the time domain because the sosfilt gives the time domain     
        
    #FREQUENCY DOMAIN OF THE FILTERED RANDOM NOISE: Frequency response of the filtered random noise
    filt_noise_band_freq = np.fft.fft(filt_noise_band)
    
    # fv is the freqeuncy vector for the x axis
    nSamples = len(noise)
    fv = np.arange(nSamples) * (fs/nSamples) #This can also be written with linspace as #np.linspace((0, (nSamples-1)))*fs/nSamples;
       
    #Make the filt_noise_band list into an array
    filt_noise_band = np.array(filt_noise_band)
    return filt_noise_band, filt_noise_band_freq, fv


def unfiltered_envelope(nBands, square_root, noise):
    """
    Calculation of unfiltered envelope

    Parameters
    ----------
        nBands : int
            Number of frequency bands.
        square_root : list of arrays
            Square rooted pressure over time after the source is switched off per each frequency resampled
        noise : array of float
            Random noise vector with unifrom distribution over time
        
    Returns
    -------
        imp_unfilt_band : list of arrays
            Multiplication of the envelope with random noise (and not filter)
            
    """
    #Multiplication of SQUARE-ROOT of envelope with random noise only (UNFILTERED)
    imp_unfilt_band = []
    for fi in range(nBands):
        imp_unfilt = square_root[fi,:] * noise
        imp_unfilt_band.append(imp_unfilt)
    return imp_unfilt_band


def filtered_envelope(filt_noise_band,p_rec_off_deriv_band_resampled, square_root, t_off_resampled, nBands):
    """
    Calculation of filtered envelope

    Parameters
    ----------
        filt_noise_band : array of floats
            Results of the time domain convolution of the noise with the filter per each frequency
        p_rec_off_deriv_band_resampled : list of arrays
            Pressure over time after the source is switched off per each frequency resampled
        square_root : list of arrays
            Square rooted pressure over time after the source is switched off per each frequency resampled
        t_off_resampled : array of float
            Time array after the source is switched off resampled
        nBands : int
            Number of frequency bands.
        
    Returns
    -------
        imp_filt_band : list of arrays
            Multiplication of the envelope with filtered random noise
        t_off_padded : 
            Time array after the source is switched off padded to the legnth of the filt_noise_band   
    """
    #Padding the square-root to the same length as the filtered random noise
    pad_length = filt_noise_band.shape[1]-p_rec_off_deriv_band_resampled.shape[1]
    square_root_padded = np.pad(square_root, ((0,0),(0,pad_length)) ,mode='constant' )
    t_off_padded = np.pad(t_off_resampled, ((0,pad_length)) ,mode='constant' )
    
    #Multiplication of SQUARE-ROOT of envelope with filtered random noise (FILTERED)
    imp_filt_band = []
    for fi in range(nBands):
        imp_filt = square_root_padded[fi,:]*filt_noise_band[fi,:]
        imp_filt_band.append(imp_filt)
    return imp_filt_band, t_off_padded


def generate_imp_resp(imp_filt_band):
    """
    Generation of impulse response

    Parameters
    ----------
        imp_filt_band : list of arrays
            Multiplication of the envelope with filtered random noise
        
    Returns
    -------
        imp_tot : array of floats
            Total impulse response in the time domain
        freq_spectrum : array of floats
            Total impulse response in the frequency domain
        imp_resp_norm : array of int
            Total impulse response in the time domain normalised to its maximum value
    """
    #Sum of the bands in the time domain
    imp_tot = [sum(imp_filt_band[i][j] for i in range(len(imp_filt_band))) for j in range(len(imp_filt_band[0]))]
    imp_tot = np.array(imp_tot, dtype=float)
    
    #Frequency domain of the total impulse response
    freq_spectrum = 20*np.log10(abs(np.fft.fft(imp_tot)))
    
    imp_resp_norm = np.int16(imp_tot / np.max(np.abs(imp_tot)) * 32767) 
    
    return imp_tot, freq_spectrum, imp_resp_norm


def convolution_fun(data_signal,fs,imp_tot):
    """
    Convolution

    Parameters
    ----------
        data_signal : array of floats 
            Signal data of the anechoic file
        fs : float
            Sample frequency of the anechoic signal
        imp_tot : array of floats
            Total impulse response in the time domain
        
    Returns
    -------
        st : array of floats
            Time vector of the speech signal
        ht : array of floats
            Time vector of the room impulse response
        sh_conv : array of floats
            Convolution of the impulse response with the anechoic signal
        t_conv : array of floats
            Time vector of the convolved signal
        sh_conv_normalized : array of floats
            Convolution of the impulse response with the anechoic signal normalised to its maximum values
        
    """
    #Convolution of the impulse_rand with the anechoic signal
    st = np.arange(0,(len(data_signal))/fs,1/fs) #Time vector of the speech signal
    ht = np.arange(0,(len(imp_tot))/fs,1/fs)  #Time vector of the room impulse response
    
    #Create impulse response
    sh_conv = np.convolve(imp_tot,data_signal) #convolution of the impulse response with the anechoic signal
    sh_conv = sh_conv/max(abs(sh_conv)) #normalized to the maximum value of the convolved signal
    
    t_conv = np.arange(0,(len(sh_conv))/fs,1/fs) #Time vector of the convolved signal
    
    # Normalize the floating-point data to the range of int16
    sh_conv_normalized = np.int16(sh_conv / np.max(np.abs(sh_conv)) * 32767) #32767 scales the normalized data to the range of 16-bit integers (-32768 to 32767).

    return st, ht, sh_conv, t_conv, sh_conv_normalized
