# Auralization Use

## Requirements
1. Download and install Anaconda or download and install any Python software IDE you prefer;
2. Clone/Fork this repository to a folder of your preference;
3. Run the "FVM.py" file through the preferred IDE and get the results;
4. Open the "Auralization.py" file through the preferred IDE and run it.

## Libraries
To properly run the software, the following libraries are needed:
- Python version 3.10.9 or above

Libraries for python:
- math
- matplotlib
- soundfile
- scipy
- pickle
- os

## Python running files 
The main files for users to run the software is _Auralization.py_ to run the auralization. The _Auralization.py_ can be run **only if** the results from the _FVM.py_ have been obtained. 

## Algorithm
The software is organised in multiple sections:
- Input variables:
    The inputs are the final results of the _FVM.py_ simulation and the anechoic signal.
- Creation of impulse response:
    From the energy density gotten from the _FVM.py_ calculation, a full impulse response is generated. This impulse response include only the reverberant part of the decay.
- Auralization:
    From the impulse response, a convolution with the anechoic signal is made to generate the .wav file.

## Inputs

### Anechoic signal
The anechoic signal is a wav file needed for the convolution with the impulse response of the room in question. An anechoic file is a sound/speech recorder inside a fully absorbing chamber and therefore without any reflections. In the repository only one anechoic file is present; the user could decide to use their own. This anechoic file has been provided by ODEON software.

### Energy density/Pressure
The energy density and pressure curves are obtained simulating the room in question by the _FVM.py_ script. it include the pressure curve over time at each frequency band per receiver. 

## Outputs

### Impulse response file
An impulse response wav file representing the reverberant fingerprint of the room is the outoput of the calculation.

### Auralization file
The calculation include also the creating of a auralization wav file.
