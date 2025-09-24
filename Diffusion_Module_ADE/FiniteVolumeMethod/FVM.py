# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 16:12:40 2023

@author: Ilaria Fichera
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
from Diffusion_Module_ADE.FiniteVolumeMethod.FunctionClarity import clarity
from Diffusion_Module_ADE.FiniteVolumeMethod.FunctionDefinition import definition
from Diffusion_Module_ADE.FiniteVolumeMethod.FunctionCentreTime import centretime
from Diffusion_Module_ADE.FiniteVolumeMethod.FunctionRT import t60_decay
from Diffusion_Module_ADE.FiniteVolumeMethod.FVMfunctions import *

#%%
###############################################################################
#RUN SIMULATION FUNCTION
###############################################################################
#def run_sim(coord_source, coord_rec,fc_low,fc_high,num_octave, dt,m_atm , c0, Ws, th, pRef, rho, file_name,  dim, tag, center_freq, tcalc = "decay"):
def run_fvm_sim(mesh_file_path, inputs_path, abs_coeff_path):
    """Function for running the full calculation. It will use all the functions defined in FVMfunctions.py file

    Args:
        mesh_file_path : str
            String of the name of the mesh file
        inputs_path : str
            String of the json file with all the inputs for the simulation.
        abs_coeff_path : str
            String of the csv file of the absorption coefficient inputs

    Raises:
        ValueError: If the values in the absorption coefficient csv file are missing, it will flag an error.

    Returns:
        results : dict
            Dictionary of all the variable calculated during the simulation
    """
    
    st = time.time() #start time of calculation
    
    #file_path = os.path.join(script_dir, mesh_file) # Full path to the file

    #Calling function %create_vgroups_names%
    vGroupsNames = create_vgroups_names(mesh_file_path)
    
    with open(os.path.join(script_dir,inputs_path), "r") as f:
        inputs = json.load(f)
         
    # Access input variables
    coord_source = inputs["coord_source"]
    coord_rec = inputs["coord_rec"]
    fc_low = inputs["fc_low"]
    fc_high = inputs["fc_high"]
    num_octave = inputs["num_octave"]
    dt = inputs["dt"]
    m_atm = inputs["m_atm"]
    th = inputs["th"]
    #file_name = inputs["file_name"]
    #center_freq = inputs["center_freq"]
    #nBands = inputs["nBands"]
    #x_frequencies = inputs["x_frequencies"]
    #vGroupsNames = inputs["vGroupsNames"]
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
    mesh = gmsh.open(mesh_file_path) #open the file
    
    dim = -1
    tag = -1
    
    #Calling function %create_vgroups_names%
    vGroupsNames = create_vgroups_names(mesh_file_path)
    
    #Calling function %number_freq%
    nBands, center_freq = number_freq(num_octave, fc_high, fc_low)
    
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
    V,S,Eq_A = equiv_absorp_area(cell_volume, total_boundArea, surface_areas, absorption_coefficient_dict)
    
    #Calling function %calculate_sourceon_time%
    sourceon_time = calculation_sourceon_time(nBands, V, Eq_A)
    
    #Calling function %rec_time%
    recording_time, t, recording_steps = calculation_rec_time(sourceon_time, dt, edt=-1, ir_length=-1)
    
    #Calling function %diff_coeff%
    Dx, Dy, Dz = diffusion_coeff(V, S, c0)
    
    #Calling function %dist_source_receiver%
    dist_sr = distance_source_receiver(coord_rec, coord_source)
    
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
    beta_zero_freq = beta_zero_freq_fun(boundary_areas, dt, Dx, interior_tet_sum, cell_volume)
    
    #Calling function %computing_energy_density%
    w_new_band, w_rec_band, w_rec_off_band, w_rec_off_deriv_band, p_rec_off_deriv_band, idx_w_rec, t_off = computing_energy_density(nBands, voluEl, recording_steps, beta_zero_freq, dt, c0, m_atm, Dx, interior_tet, cell_volume, s, cl_tet_r_keys, total_weights_r, tcalc, cl_tet_s_keys, source1, total_weights_s, t, sourceon_time, rho)
    print("100% of main calculation completed")
    
    #Calling function %freq_parameters%
    w_rec_x_band, w_rec_y_band, spl_stat_x_band, spl_stat_y_band, spl_r_band, spl_r_off_band, spl_r_norm_band, sch_db_band, t30_band, edt_band, c80_band, d50_band, ts_band = freq_parameters(nBands, line_rec_x_idx_list, w_new_band, line_rec_y_idx_list, rho, c0, Ws, dist_x, dist_y, pRef, w_rec_band, w_rec_off_band, tcalc, t, idx_w_rec, V, Eq_A, S, dist_sr)
    
    et = time.time() #end time
    elapsed_time = et - st
    
    #results = locals()
    
    results = {"coord_source" : coord_source,
               "coord_rec" : coord_rec,
               "fc_low" : fc_low,
               "fc_high" : fc_high,
               "num_octave" : num_octave,
               "dt" : dt,
               "m_atm" : m_atm,
               "th" : th,
               "tcalc" : tcalc,
               "df_abs" : df_abs,
               "center_freq" : center_freq,
               "t": t,
               "w_new_band" : w_new_band, 
               "w_rec_band" : w_rec_band, 
               "w_rec_off_band" : w_rec_off_band, 
               "w_rec_off_deriv_band" : w_rec_off_deriv_band, 
               "p_rec_off_deriv_band" : p_rec_off_deriv_band, 
               "idx_w_rec" : idx_w_rec, 
               "t_off" : t_off,
               "w_rec_x_band" : w_rec_x_band, 
               "w_rec_y_band" : w_rec_y_band, 
               "spl_stat_x_band" : spl_stat_x_band, 
               "spl_stat_y_band" : spl_stat_y_band, 
               "spl_r_band" : spl_r_band, 
               "spl_r_off_band" : spl_r_off_band, 
               "spl_r_norm_band" : spl_r_norm_band, 
               "sch_db_band" : sch_db_band, 
               "t30_band" : t30_band, 
               "edt_band" : edt_band, 
               "c80_band" : c80_band, 
               "d50_band" : d50_band, 
               "ts_band" : ts_band,
               }
    
    save_fvm('resultsFVM.pkl',results)
        
    print("Simulation finished successfully! Results in resultsFVM.pkl file")
    
    return results