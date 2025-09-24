"""
Simulate the energy decay in a cuboid room using the Finite Volume diffusion
equation.
"""
# %%
import os
import pandas as pd
import json
import numpy as np
import matplotlib.pyplot as plt
from Diffusion_Module_ADE.FiniteDifferenceMethod.FDM import run_fdm_sim

script_dir = os.path.dirname(os.path.abspath(__file__))

#%%
###############################################################################
#GENERAL INPUT VARIABLES
###############################################################################
input_data = {
    "room_dim": [3.0, 3.0, 3.0],
    "coord_source": [1.5, 1.5, 1.5], #source coordinates x,y,z
    "coord_rec": [2.0, 1.5, 1.5], #rec coordinates x,y,z
    "alpha_1": [0.10, 0.15, 0.20, 0.25, 0.25, 0.30], #Absorption coefficient for Surface1 - Floor
    "alpha_2": [0.07, 0.10, 0.13, 0.15, 0.15, 0.16], #Absorption coefficient for Surface2 - Ceiling
    "alpha_3": [0.08, 0.09, 0.11, 0.15, 0.14, 0.14], #Absorption coefficient for Surface3 - Wall Front
    "alpha_4": [0.08, 0.09, 0.11, 0.15, 0.14, 0.14], #Absorption coefficient for Surface4 - Wall Back
    "alpha_5": [0.08, 0.09, 0.11, 0.15, 0.14, 0.14], #Absorption coefficient for Surface5 - Wall Left
    "alpha_6": [0.08, 0.09, 0.11, 0.15, 0.14, 0.14], #Absorption coefficient for Surface6 - Wall Right
    "fc_low": 125, #lowest frequency
    "fc_high": 4000, #highest frequency
    "num_octave": 1, # 1 or 3 depending on how many octave you want
    "dx": 0.5,
    "dt": 1/8000, #time discretization
    "m_atm": 0, #air absorption coefficient [1/m]
    "th": 3, #int(input("Enter type Absortion conditions (option 1,2,3):")) # options Sabine (th=1), Eyring (th=2) and modified by Xiang (th=3)
    "tcalc": "decay" #Choose "decay" if the objective is to calculate the energy decay of the room with all its energetic parameters; Choose "stationarysource" if the aim is to understand the behaviour of a room subject to a stationary source
}

#%%
###############################################################################
#SAVE TO JSON
###############################################################################
fname_input_configuration = "cube_input_fdm.json"
with open(os.path.join(script_dir, fname_input_configuration), "w") as f:
    json.dump(input_data, f, indent=4)

print("Input file successfully created: cube_input_fdm.json")

# %%
result = run_fdm_sim(os.path.join(script_dir, fname_input_configuration))
# %%
times = result['t'][:len(result['t'])//2]
energy_decay_curve = np.array(result['w_rec_off_band'])
# %%
