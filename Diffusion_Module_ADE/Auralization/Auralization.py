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
import scipy
import pickle
import os
from Diffusion_Module_ADE.Auralization.Auralizationfunctions import *

def run_auralization(anechoic_signal_path,resultsFVM_path):
    """
    Function for running the full auralization calculation. It will use all the functions defined in Auralizationfunctions.py file

    Parameters
    ----------
        anechoic_signal_path : str
            String of the wav file of the anechoic sound
        resultsFVM_path : str
            String of the results from the FVM simulation
        
    Returns
    -------
        results : dict
            Dictionary of all the variable calculated including auralization and impulse response wav files.
    """
  
    with open(resultsFVM_path, "rb") as f:
        resultsFVM = pickle.load(f)
   
    #Import data needed from the resultsFVM pickle file
    dt_sim = resultsFVM["dt"]   #Import delta t (time step) from the simulation calc
    t_off_sim = resultsFVM["t_off"]   #Import the time t array since the source has been switched off from the simulation calc
    p_rec_off_deriv_band = np.array(resultsFVM["p_rec_off_deriv_band"]) #import pressure curve from the simulation calc
    center_freq = np.array(resultsFVM["center_freq"]) #import frequency bands from the simulation calc
    
    #%%
    ###############################################################################
    #GETTING FIXED DATA
    ###############################################################################
    fs_sim = 1/dt_sim #sampling frequency
    t_off_sim = t_off_sim - t_off_sim[0] #removing the t_off[0] to make the vector start from zero.
    center_freq = center_freq.astype(np.int32) #centre frequency considered as integers
    nBands = len(center_freq) #number of frequency bands
    
    data_signal, fs = extract_anechoic(anechoic_signal_path)
    #sd.play(data, fs) #this line allows to listen to the anechoic signal as it is.
    #status = sd.wait()  #Wait until file is done playing
    
    p_rec_off_deriv_band_resampled, t_off_resampled = resample_pressure(p_rec_off_deriv_band, fs, fs_sim, t_off_sim)
    # #FIGURE 1
    # plt.figure(figsize=(12, 8))
    # plt.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)  #This line disables ticks on the main figure
    # plt.title('Time domain envelope of squared impulse response per octave band filters')
    # for fi in range(nBands):
    #     plt.subplot(nBands, 1, fi+1)
    #     plt.plot(t_off, p_rec_off_deriv_band[fi], label=f'{center_freq[fi]} Hz ORIGINAL')
    #     plt.plot(t_off_resampled,p_rec_off_deriv_band_resampled[fi], label=f'{center_freq[fi]} Hz RESAMPLED')
    
    #     plt.xlabel('Time [s]')
    #     plt.ylabel('Magnitude [dB]')
    #     plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    #     plt.axhline(0, color='black', linewidth=0.5)
    #     plt.legend(loc='best')
    # plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust layout to fit suptitle
    # plt.show()

    square_root = square_root_fun(p_rec_off_deriv_band_resampled)
    
    noise = create_noise(p_rec_off_deriv_band_resampled)
    # #FIGURE 2
    # plt.figure(figsize=(12, 8))
    # plt.plot(t_off_resampled,noise)
    

    filter_tot = butter_filter(fs, center_freq)
    # #FIGURE 3
    # # Plot frequency responses
    # plt.figure(figsize=(12, 8))
    # for band in range(len(filter_tot)):
    #     #print(band)
    #     # Compute the frequency response of each filter
    #     w, h = sosfreqz(filter_tot[band], worN=2000, fs=fs)
    #     # w is the array of frequencies at which the response is computed
    #     # h is the frequency response of the filter

    #     # Plot the magnitude response in decibels
    #     plt.semilogx(w, 20 * np.log10(abs(h)), label=f'{center_freq[band]}')
    # plt.title('Frequency response of octave band filters')
    # plt.xlabel('Frequency [Hz]')
    # plt.ylabel('Gain [dB]')
    # plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    # plt.axhline(0, color='black', linewidth=0.5)
    # plt.legend(loc='best')
    # plt.show()
    
    
    
    filt_noise_band, filt_noise_band_freq, fv = filtered_noise(noise, filter_tot, fs)
    # #FIGURE 4
    # #Plot the time domain of the filtered noise
    # plt.figure(figsize=(12, 8))
    # plt.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)  #This line disables ticks on the main figure
    # plt.title('Time domain response of octave band filtered random noise')
    # for fi in range(nBands):
    #     plt.subplot(nBands, 1, fi+1)
    #     plt.plot(t_off_resampled, 20 * np.log10(abs(filt_noise_band[fi])), label=f'{center_freq[fi]} Hz')
    
    #     plt.xlabel('Time [s]')
    #     plt.ylabel('Magnitude [dB]')
    #     plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    #     plt.axhline(0, color='black', linewidth=0.5)
    #     plt.legend(loc='best')
    # plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust layout to fit suptitle
    # plt.show()
        
    # #FIGURE 5
    # #Plot the frequency response
    # plt.figure(figsize=(12, 8))
    # plt.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)  #This line disables ticks on the main figure
    # plt.title('Frequency response of octave band filtered random noise')
    # for fi in range(nBands):
    #     plt.subplot(nBands, 1, fi+1)
    #     plt.semilogx(fv, 20 * np.log10(abs(filt_noise_band_freq[fi, :])), label=f'{center_freq[fi]} Hz')
    
    #     plt.xlabel('Frequency [Hz]')
    #     plt.ylabel('Gain [dB]')
    #     plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    #     #plt.axhline(0, color='black', linewidth=0.1)
    #     plt.legend(loc='best')
    # plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust layout to fit suptitle
    # plt.show()
    

    imp_unfilt_band = unfiltered_envelope(nBands, square_root, noise)
    # #FIGURE 6
    # plt.figure(figsize=(12, 8))
    # plt.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)  #This line disables ticks on the main figure
    # plt.title('Time domain of UNfiltered impulse response per frequency band')
    # for fi in range(nBands):
    #     plt.subplot(nBands, 1, fi+1)
    #     plt.plot(t_off_resampled, imp_unfilt_band[fi], label=f'{center_freq[fi]} Hz')
        
    #     plt.xlabel('Time [s]')
    #     plt.ylabel('Magnitude [dB]')
    #     plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    #     #plt.axhline(0, color='black', linewidth=0.1)
    #     plt.legend(loc='best')
    # plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust layout to fit suptitle
    # plt.show()
    
    imp_filt_band, t_off_padded = filtered_envelope(filt_noise_band,p_rec_off_deriv_band_resampled, square_root, t_off_resampled, nBands) 
    # #FIGURE 7
    # plt.figure(figsize=(12, 8))
    # plt.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)  #This line disables ticks on the main figure
    # plt.title('Time domain of filtered impulse response per frequency band')
    # for fi in range(nBands):
    #     plt.subplot(nBands, 1, fi+1)    
    #     plt.plot(t_off_padded, imp_filt_band[fi], label=f'{center_freq[fi]} Hz')
        
    #     plt.xlabel('Time [s]')
    #     plt.ylabel('Magnitude [dB]')
    #     plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    #     #plt.axhline(0, color='black', linewidth=0.1)
    #     plt.legend(loc='best')
    # plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust layout to fit suptitle
    # plt.show()
    
    imp_tot, freq_spectrum, imp_resp_norm = generate_imp_resp(imp_filt_band)  
    # #FIGURE 8
    # # Plot impulse response in the time domain
    # plt.figure(figsize=(12, 8))
    # plt.plot(t_off_padded, imp_tot)
    # plt.xlabel('Time [s]')
    # plt.ylabel('Magnitude [dB]')
    # plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    # #plt.axhline(0, color='black', linewidth=0.1)
    # #plt.legend(loc='best')
    # plt.show()
       
    # #FIGURE 9
    # # Plot impulse response in the frequency domain
    # plt.figure(figsize=(12, 8))
    # plt.semilogx(fv, freq_spectrum)
    # plt.xlabel('Frequency [Hz]')
    # plt.ylabel('Magnitude [dB]')
    # plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    # #plt.axhline(0, color='black', linewidth=0.1)
    # #plt.legend(loc='best')
    # plt.show()
    
    scipy.io.wavfile.write("imp_resp.wav", fs, imp_resp_norm)
    #Play the impulse response
    #sd.play(imp_tot, fs)
    
    print("Starting convolution...")
    st, ht, sh_conv, t_conv, sh_conv_normalized = convolution_fun(data_signal,fs,imp_tot)
    #Play the convolved signal
    #sd.play(sh_conv, fs)
    
    # #FIGURE 10
    # plt.figure(figsize=(12, 8))
    # plt.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)  #This line disables ticks on the main figure
    # plt.title('Convolution')
    # signals = [
    #     (st, data_signal, 'anechoic_signal'),
    #     (ht, imp_tot, 'impulse_response'),
    #     (t_conv, sh_conv, 'convolved_signal')
    # ]
    # for fi, (x, y, label) in enumerate(signals):
    #     plt.subplot(3, 1, fi + 1)
    #     plt.plot(x, y, label=label)
    #     plt.xlabel('Time [s]')
    #     plt.ylabel('Magnitude [dB]')
    #     plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    #     plt.legend(loc='best')
    # plt.tight_layout(rect=[0, 0, 1, 0.96])
    # plt.show()

    # Write the normalized data to a WAV file
    scipy.io.wavfile.write("auralization.wav", fs, sh_conv_normalized)
    
    results = locals()
    
    # remove any un-pickleable objects
    for k, v in list(results.items()):
        if hasattr(v, "read"):   # catches file-like objects
            print(f"Removing unpickleable object: {k}")
            del results[k]
    
    with open('resultsAuralization.pkl', 'wb') as f:
        pickle.dump(results, f)
        
    print("Simulation finished successfully! Results in resultsAuralization.pkl file")
    
    return results

#%%
###############################################################################
#INPUT VARIABLES
###############################################################################
script_dir = os.path.dirname(os.path.abspath(__file__))

# Load input data
if __name__ == "__main__":
    
    #Calling function %run_auralization_sim%
    results = run_auralization('Frequency(english).wav','resultsFVM.pkl')
