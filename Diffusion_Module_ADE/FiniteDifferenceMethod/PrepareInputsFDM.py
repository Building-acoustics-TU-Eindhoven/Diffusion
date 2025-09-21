# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 11:52:04 2025

@author: 20225533
"""

import os
import json
#%%
###############################################################################
#GENERAL INPUT VARIABLES
###############################################################################
input_data = {
    "room_dim": [8.0, 8.0, 8.0],
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

script_dir = os.path.dirname(os.path.abspath(__file__)) #script directory

# Save to JSON
with open(os.path.join(script_dir,"simulation_fdm_inputs.json"), "w") as f:
    json.dump(input_data, f, indent=4)

print("Input file successfully created: simulation_inputs.json")
