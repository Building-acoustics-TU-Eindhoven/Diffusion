# Auralization Use

## Requirements
1. Set up acousticDE following the instructions in the installation section. 
3. Run the "FVM.py" file and get the results;
4. Use these results for the auralization.

## Usage & files
To use the software, the following files are to be used:
- _Auralization.py_: it contains the main function run_auralization_sim to run the full simulation and generate an auralization file.

The main software works with the following associated functions:
+ _Auralizationfunctions.py_ include all the main functions that are used in the full simulation;


## Inputs

### Anechoic signal path
The anechoic signal is a wav file needed for the convolution with the impulse response of the room in question. An anechoic file is a sound/speech recorder inside a fully absorbing chamber and therefore without any reflections. In the repository only one anechoic file is present; the user could decide to use their own. This anechoic file has been provided by ODEON software.

### resultsFVM.pkl path
The result file from the FVM calculation is needed to run the script. 

## Outputs

### Impulse response file
An impulse response wav file representing the reverberant fingerprint of the room is the outoput of the calculation.

### Auralization file
The calculation include also the creating of a auralization wav file.
