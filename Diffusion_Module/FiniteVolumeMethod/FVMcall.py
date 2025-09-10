# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 16:12:40 2023

@author: 20225533
"""

#%%
###############################################################################
#IMPORT LIBRARIES
###############################################################################
#Code developed by Ilaria Fichera for the analysis of the FVM method adapted solving the 3D diffusion equation with one intermittent omnidirectional sound source
#Import modules

import time
import os
from math import log
import gmsh
import numpy as np
import json
import pandas as pd
from Diffusion_Module.FiniteVolumeMethod.FunctionClarity import clarity
from Diffusion_Module.FiniteVolumeMethod.FunctionDefinition import definition
from Diffusion_Module.FiniteVolumeMethod.FunctionCentreTime import centretime
from Diffusion_Module.FiniteVolumeMethod.FunctionRT import t60_decay
from Diffusion_Module.FiniteVolumeMethod.FVMfunctions import *


#%%
###############################################################################
#INPUT VARIABLES
###############################################################################

# # Source position
# coord_source = [0.5,0.5,1.0] #coordinates of the receiver position in an list: x , y and z direction in meters [m]

# # Receiver position
# coord_rec = [2.0,0.5,1.0] #coordinates of the receiver position in an list: x , y and z direction in meters [m]

# # Type of Calculation
# #Choose "decay" if the objective is to calculate the energy decay of the room with all its energetic parameters; 
# #Choose "stationarysource" if the aim is to understand the behaviour of a room subject to a stationary source
# tcalc = "decay"

# # Frequency range
# fc_low = 125
# fc_high = 2000
# num_octave = 1    

# # Time discretization
# dt = 1/20000 #time discretization

# # Air absorption coefficient
# m_atm = 0 #air absorption coefficient [1/m]

# # Absorption term and Absorption coefficients
# th = 3 #int(input("Enter type Absortion conditions (option 1,2,3):")) 
# # options Sabine (th=1), Eyring (th=2) and modified by Xiang (th=3)

# file_name = "3x3x3.msh" #Insert file name, msh file created from sketchUp and then gmsh

# #%%
# ###############################################################################
# #FIXED INPUTS AND INITIALISE GMSH
# ###############################################################################
# #Folder where to run
# script_dir = os.path.dirname(os.path.abspath(__file__))
# file_path = os.path.join(script_dir, file_name) # Full path to the file

# # ####?????? Ilaria: here it initialise the file. This is also important because the gmsh.exe file needs to be in the same folder as the code!!!
# # gmsh.initialize() #Initialize msh file
# # mesh = gmsh.open(file_path) #open the file

# # dim = -1 #dimensions of the entities, 0 for points, 1 for curves/edge/lines, 2 for surfaces, 3 for volumes, -1 for all the entities 
# # tag = -1 #all the nodes of the room

# #Calling function %create_vgroups_names% (This needs to stay outside the run_sim function because then it is used for the input of the absorption coefficient)
# vGroupsNames = create_vgroups_names(file_path)

# #Calling function %number_freq%
# x_frequencies, nBands, center_freq = number_freq(num_octave,fc_high,fc_low)


# #%%
# ###############################################################################
# #INPUTS ABSORPTION COEFFICIENT VIA CSV
# ###############################################################################

# import pandas as pd
# import sys

# # Define frequency bands
# freq_bands = center_freq  # e.g., [125, 250, 500, 1000, 2000]
# column_names = ["Surface"] + [f"{int(fc)}Hz" for fc in freq_bands]
# csv_path = os.path.join(script_dir, "absorption_coefficients.csv")

# # Check if CSV exists, otherwise create template
# if not os.path.exists(csv_path):
#     surface_names = [group[2] for group in vGroupsNames if group[0] == 2]
#     df_template = pd.DataFrame(columns=column_names)
#     df_template["Material"] = surface_names
#     df_template.to_csv(csv_path, index=False)
#     print(f"Template created: {csv_path}. Please fill in the absorption coefficients.")
#     sys.exit()


# #General settings
# c0= 343 #adiabatic speed of sound [m.s^-1]

# #Set initial condition - Source Info (interrupted method)
# Ws = 0.01 #Source point power [Watts] interrupted after "sourceon_time" seconds; 10^-2 W => correspondent to 100dB

# # Reference pressure in Pa
# pRef = 2 * (10**-5) #Reference pressure in Pa

# # Air density
# rho = 1.21 #air density [kg.m^-3] at 20°C

# #%%
# ###############################################################################
# #ABSORPTION COEFFICIENT INPUT
# ###############################################################################

# # Load filled CSV
# df_abs = pd.read_csv(csv_path)


# # Initialize a list to store surface tags and their absorption coefficients
# surface_absorption = [] #initialization absorption term (alpha*surfaceofwall) for each wall of the room
# triangle_face_absorption = [] #initialization absorption term for each triangle face at the boundary and per each wall
# absorption_coefficient_dict = {}

# for group in vGroupsNames:
#     if group[0] != 2:
#         continue

#     name_group = group[2]
#     row = df_abs[df_abs["Material"] == name_group]
#     if row.empty:
#         raise ValueError(f"No absorption data found for surface: {name_group}")

#     abscoeff = row.iloc[0, 1:].values.astype(float).tolist() #input(f"Enter absorption coefficient for frequency {fc_low} to {fc_high} for {name_abs_coeff}:") 
#     surface_materials (group, abscoeff, surface_absorption, absorption_coefficient_dict,nBands,th,c0)

# for entity, Abs_term in surface_absorption:
#     triangle_faces, _ = gmsh.model.mesh.getElementsByType(2, entity) #Get all the triangle faces for the current surface
#     triangle_face_absorption.extend([Abs_term] * len(triangle_faces)) #Append the Abs_term value for each triangle face

# print("Correctly inputted surface materials. Starting initial geometry calculations...")


# Load input data
with open("simulation_inputs.json", "r") as f:
    inputs = json.load(f)


#%%
###############################################################################
#RUN SIMULATION FUNCTION
###############################################################################
#def run_sim(coord_source, coord_rec,fc_low,fc_high,num_octave, dt,m_atm , c0, Ws, th, pRef, rho, file_name,  dim, tag, center_freq, tcalc = "decay"):
def run_sim(inputs, abs_coeff_path):
    
    st = time.time() #start time of calculation
         
    # Access input variables
    coord_source = inputs["coord_source"]
    coord_rec = inputs["coord_rec"]
    fc_low = inputs["fc_low"]
    fc_high = inputs["fc_high"]
    num_octave = inputs["num_octave"]
    dt = inputs["dt"]
    m_atm = inputs["m_atm"]
    th = inputs["th"]
    file_name = inputs["file_name"]
    center_freq = inputs["center_freq"]
    nBands = inputs["nBands"]
    x_frequencies = inputs["x_frequencies"]
    vGroupsNames = inputs["vGroupsNames"]
    tcalc = inputs.get("tcalc", "decay")  # default fallback
    
    df_abs = pd.read_csv(abs_coeff_path)
    # validate and convert to dictionary
    if df_abs.isnull().values.any():
        raise ValueError("Absorption coefficient file contains missing values.")
        
    c0= 343 #adiabatic speed of sound [m.s^-1]
    Ws = 0.01 #Source point power [Watts] interrupted after "sourceon_time" seconds; 10^-2 W => correspondent to 100dB
    pRef = 2 * (10**-5) #Reference pressure in Pa
    rho = 1.21 #air density [kg.m^-3] at 20°C
    
    gmsh.initialize() #Initialize msh file
    mesh = gmsh.open(file_name) #open the file
    
    dim = -1
    tag = -1
    
    #Calling function %get_nodes_elem%
    nodecoords, node_indices, bounEl, bounNode, voluEl, voluNode, belemNodes, velemNodes, boundaryEl_dict, volumeEl_dict = get_nodes_elem(dim,tag)
    
    #Calling function %velem_volume_centre%
    cell_center, cell_volume = velem_volume_centre(volumeEl_dict, nodecoords, node_indices)

    #Calling function %belem_area_centre%
    barea_dict, centre_area = belem_area_centre(boundaryEl_dict, nodecoords, node_indices)
    
    #Calling function %get_neighbour_faces%
    fxt, txt, neighbourVolume = get_neighbour_faces(voluEl)
    print("Completed initial geometry calculation. Starting internal tetrahedrons calculations...")
    
    #Calling function %interior_tetra%
    interior_tet, interior_tet_sum = interior_tetra(voluEl, cell_center, velemNodes, nodecoords, node_indices)
    print("Completed internal tetrahedrons calculation. Starting boundary tetrahedrons calculations...")
    
    #Calling function %surface_absorption%
    surface_absorption, triangle_face_absorption, absorption_coefficient_dict = surface_absorption_fun(vGroupsNames,df_abs,nBands,th,c0)
    
    #Calling function %surface_area%
    surface_areas = surface_area(surface_absorption, nodecoords, node_indices)
    
    #Calling function %boundary_triang%
    boundary_areas, total_boundArea = boundary_triang(velemNodes, nBands, bounNode, nodecoords, node_indices, triangle_face_absorption)
    print("Completed boundary tetrahedrons calculation. Starting main diffusion equation calculations over time and frequency...")
    
    gmsh.finalize()
    
    #Calling function %equiv_absorp%
    V,S,Eq_A = equiv_absorp(cell_volume, total_boundArea, surface_areas, absorption_coefficient_dict)
    
    #Calling function %calculate_sourceon_time%
    sourceon_time = calculate_sourceon_time(nBands, V, Eq_A)
    
    #Calling function %rec_time%
    recording_time, t, recording_steps = rec_time(sourceon_time, dt, edt=-1, ir_length=-1)
    
    #Calling function %diff_coeff%
    Dx, Dy, Dz = diff_coeff(V, S, c0)
    
    #Calling function %dist_source_receiver%
    dist_sr = dist_source_receiver(coord_rec, coord_source)
    
    #Calling function %source_interp%
    cl_tet_s_keys, total_weights_s = source_interp(cell_center, coord_source)
    
    #Calling function %source_volume%
    Vs = source_volume(velemNodes, nodecoords, coord_source, cell_volume)
    
    #Calling function %initial_cond%
    source1, sourceon_steps = initial_cond(Ws, Vs, sourceon_time, dt, recording_steps)
    
    #Calling function %source_matrix%
    s = source_matrix(voluEl,cl_tet_s_keys, source1, total_weights_s)
    
    #Calling function %receiver_interp%
    cl_tet_r_keys, total_weights_r = receiver_interp(cell_center, coord_rec)
    
    #Calling function %room_dimensions%
    room_length, room_width, room_height = room_dimensions(nodecoords)
    
    #Calling function %line_receivers%
    x_axis, y_axis, line_rec_x_idx_list, dist_x, line_rec_y_idx_list, dist_y = line_receivers(room_length, room_width, coord_rec, coord_source, cell_center)
    
    #Calling function %beta_zero%
    beta_zero_freq = beta_zero(boundary_areas, dt, Dx, interior_tet_sum, cell_volume)
    
    #Calling function %computing_energy_density%
    w_new_band, w_rec_band, w_rec_off_band, w_rec_off_deriv_band, p_rec_off_deriv_band, idx_w_rec, t_off = computing_energy_density(nBands, voluEl, recording_steps, beta_zero_freq, dt, c0, m_atm, Dx, interior_tet, cell_volume, s, cl_tet_r_keys, total_weights_r, tcalc, cl_tet_s_keys, source1, total_weights_s, t, sourceon_time, rho)
    print("100% of main calculation completed")
    
    #Calling function %freq_parameters%
    w_rec_x_band, w_rec_y_band, spl_stat_x_band, spl_stat_y_band, spl_r_band, spl_r_off_band, spl_r_norm_band, sch_db_band, t30_band, edt_band, c80_band, d50_band, ts_band = freq_parameters(nBands, line_rec_x_idx_list, w_new_band, line_rec_y_idx_list, rho, c0, Ws, dist_x, dist_y, pRef, w_rec_band, w_rec_off_band, tcalc, t, idx_w_rec, V, Eq_A, S, dist_sr)
    
    et = time.time() #end time
    elapsed_time = et - st
    
    results = locals()
    
    #Calling function %save%
    save('resultsFVM.pkl', results)
    
    return results

#Calling function %run_sim%
#results = run_sim(coord_source, coord_rec,fc_low,fc_high,num_octave, dt,m_atm , c0, Ws, th, pRef, rho, file_name, dim, tag, center_freq,tcalc = "decay")
results = run_sim(inputs,'absorption_coefficients.csv')