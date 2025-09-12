# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 14:37:26 2025

@author: 20225533
"""

import gmsh
import os
import pandas as pd
import sys
import json
from Diffusion_Module.FiniteVolumeMethod.FVMfunctions import *

#%%
###############################################################################
#GENERAL INPUT VARIABLES
###############################################################################
input_data = {
    "file_name": "3x3x3.msh", #name of the mesh file
    "coord_source": [0.5, 0.5, 1.0], #source coordinates x,y,z
    "coord_rec": [2.0, 0.5, 1.0], #rec coordinates x,y,z
    "fc_low": 125, #lowest frequency
    "fc_high": 2000, #highest frequency
    "num_octave": 1, # 1 or 3 depending on how many octave you want
    "dt": 1/20000, #time discretization
    "m_atm": 0, #air absorption coefficient [1/m]
    "th": 3, #int(input("Enter type Absortion conditions (option 1,2,3):")) # options Sabine (th=1), Eyring (th=2) and modified by Xiang (th=3)
    "tcalc": "decay" #Choose "decay" if the objective is to calculate the energy decay of the room with all its energetic parameters; Choose "stationarysource" if the aim is to understand the behaviour of a room subject to a stationary source
}

script_dir = os.path.dirname(os.path.abspath(__file__))
file_path = os.path.join(script_dir, input_data["file_name"]) # Full path to the file

#Calling function %create_vgroups_names%
vGroupsNames = create_vgroups_names(file_path)
input_data["vGroupsNames"] = vGroupsNames

#Calling function %number_freq%
x_frequencies, nBands, center_freq = number_freq(input_data["num_octave"], input_data["fc_high"], input_data["fc_low"])
input_data["center_freq"] = center_freq.tolist()
input_data["nBands"] = nBands
input_data["x_frequencies"] = x_frequencies

# Save to JSON
with open(os.path.join(script_dir,"simulation_inputs.json"), "w") as f:
    json.dump(input_data, f, indent=4)

print("Input file successfully created: simulation_inputs.json")
    
gmsh.finalize()

#%%
###############################################################################
#CREATION OF CSV FOR ABSORPTION
###############################################################################

#Creation of csv
column_names = ["Material"] + [f"{int(fc)}Hz" for fc in center_freq]
csv_path = os.path.join(script_dir, "absorption_coefficients.csv")

# Check if CSV exists, otherwise create template
if not os.path.exists(csv_path):
    surface_names = [group[2] for group in vGroupsNames if group[0] == 2]
    df_template = pd.DataFrame(columns=column_names)
    df_template["Material"] = surface_names
    df_template.to_csv(csv_path, index=False)
    print(f"Template created: {csv_path}. Please fill in the absorption coefficients.")
    sys.exit()