# Auralization Documentation

An auralization, is an audible file that is created by the convolution of an anaechoic audio signal with a simulated room impulse response. It shows how a specific sound (the sound produced by the anechoic file) would sound in the room specifically simulated (room impulse response). The diffusion eqaution method is the method used to produce the energy densities inside the room, which needs to be converted into a time-domain impulse reponse usable for auralization.

In order to create an auralization, the first step is to get the results from the diffusion equation method. 
In particular, the diffusion equation simulator provides the energy density at the receiver position. To make it audible, acoustic pressure is needed which is related to the energy density by the formula $p_{rec}^2 = w_{rec} \rho c_0^2$, where w_{rec} is the energy density at the receiver, $\rho$ is the air density in $\frac{kg}{m^3}$ and $c_0$ is the speed of sound.   

By using the _FVM.py_ script, the acoustic pressure is already converted into an Envelope of the square impulse response at the receiver position for each frequency band, called *p_\{rec\}_\{off\}_\{deriv\}_\{band\}*. This is calculated by doing the opposite of a backward integration of the pressure (basically differenziation).

After getting the file signal depending on the characteristics of the room from the _FVM.py_ script, the envelope of the square impulse response at the receiver position for each frequency band is resampled to the sampling frequency of the anechoic signal chosen. This is to have consistent temporal resolutions between the two signals.

This consists only of square rooting the resampled envelope of the squared impulse response, since energy is proportional to pressure squared. From the squared envelope of the impulse response, the envelope of the impulse response is calculated. This step converts energy information into amplitude information.

This step include the creation of random noise, to represent the random phase fluctuations of the diffuse sound field. In particular a sequence of uniformly distributed random numbers with a mean of zero and a variance of one is introduced, ensuring that the noise is 'neutral' before being shaped by filters.

Since each octave band is modeled separately in the diffusion equation model, a Butterwoth bandpass filter ensures that the random noise occupies the correct frequency range. This step include the creation of the filter for the random noise. A order 8 for the bandpass filter ensures that the transition band is shapr enough to isolate the frequency regions while keeping stability.

The filtered random noise is created by convolving the time domain of the Butterworth bandpass filter and the random noise for each frequency band. This gives the noise in each frequency the correct spectral charateristic. 

The square-root of the envelope of the impulse response (calculated at step 2) is to be padded to get the same length of the filtered random noise (calculated at step 5). Padding ensures temporal alignment without distorting the envelope and allows for element-wise multiplication.

This step involves the multiplication between the filtered random noise and the square root of the resampled envolope of the squared impulse response. It shapes the noise in time according to how energy decays in the room.

The time vectors for each frequency band are summed element-wise, creating the total broadband impulse response without direct sound.

The final step is the convolution of the total broadband simulated impulse response to a dry recording (anechoic file). This allows the listener to enter the simulated room acoustics environment and hear how the soundwould propagate in that specific room.
