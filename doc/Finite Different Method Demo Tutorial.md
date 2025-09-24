# Finite Different Method Demo Tutorial

Test the software with the following inputs (Navarro et al., 2012):

```
input_data = {
    "room_dim": [8.0, 8.0, 8.0], #dimension of the room x,y,z
    "coord_source": [4.0, 4.0, 4.0], #source coordinates x,y,z
    "coord_rec": [2.0, 2.0, 2.0], #rec coordinates x,y,z
    "alpha_1": [1/6, 1/6, 1/6, 1/6, 1/6, 1/6], #Absorption coefficient for Surface1 - Floor
    "alpha_2": [1/6, 1/6, 1/6, 1/6, 1/6, 1/6], #Absorption coefficient for Surface2 - Ceiling
    "alpha_3": [1/6, 1/6, 1/6, 1/6, 1/6, 1/6], #Absorption coefficient for Surface3 - Wall Front
    "alpha_4": [1/6, 1/6, 1/6, 1/6, 1/6, 1/6], #Absorption coefficient for Surface4 - Wall Back
    "alpha_5": [1/6, 1/6, 1/6, 1/6, 1/6, 1/6], #Absorption coefficient for Surface5 - Wall Left
    "alpha_6": [1/6, 1/6, 1/6, 1/6, 1/6, 1/6], #Absorption coefficient for Surface6 - Wall Right
    "fc_low": 125, #lowest frequency
    "fc_high": 4000, #highest frequency
    "num_octave": 1, # 1 or 3 depending on how many octave you want
    "dx": 0.5,
    "dt": 1/8000, #time discretization
    "m_atm": 0, #air absorption coefficient [1/m]
    "th": 3, #int(input("Enter type Absortion conditions (option 1,2,3):")) # options Sabine (th=1), Eyring (th=2) and modified by Xiang (th=3)
    "tcalc": "decay" #Choose "decay" if the objective is to calculate the energy decay of the room with all its energetic parameters; Choose "stationarysource" if the aim is to understand the behaviour of a room subject to a stationary source
}
```

Test if the software provides the following results (Navarro et al., 2012):

- Reverberation time (RT): 1.18 s; 
- Early Decay Time (EDT): 1.18 s;
- Clarity ($C_{80}$): 2.08 dB;
- Definition ($D_{50}$): 45.63 %
- Centre Time ($T_{s}$): 85.33 ms

The result file is a pickle file called _resultsFDM.pkl_. All the results are included in this file.

## References
- J. M. Navarro, J. Escolano, J. J. Lopez, Implementation and evaluation of a diffusion equation model based on finite difference schemes for sound field prediction in rooms, Applied Acoustics 73 (6-7) (2012) 659–665.
