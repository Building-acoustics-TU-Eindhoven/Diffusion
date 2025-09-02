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
from Diffusion_Module.FiniteVolumeMethod.FunctionClarity import clarity
from Diffusion_Module.FiniteVolumeMethod.FunctionDefinition import definition
from Diffusion_Module.FiniteVolumeMethod.FunctionCentreTime import centretime
from Diffusion_Module.FiniteVolumeMethod.FunctionRT import t60_decay
from Diffusion_Module.FiniteVolumeMethod.FVMfunctions import *
import logging


####?????? Ilaria: Do we need this logger??
# Create logger for this module
logger = logging.getLogger(__name__)
####??????

#%%
###############################################################################
#FVM
###############################################################################

####?????? Ilaria: Do we need this main etc..??
if __name__ == '__main__':
    st = time.time() #start time of calculation
####??????

    #%%
    ####?????? Ilaria: Now, these for me are the things that needs to stay outside the run_sim function and these are the inputs! 
    
    ###############################################################################
    #INPUT VARIABLES
    ###############################################################################

    # Source position
    coord_source = [0.5,0.5,1.0] #coordinates of the receiver position in an list: x , y and z direction in meters [m]

    # Receiver position
    coord_rec = [2.0,0.5,1.0] #coordinates of the receiver position in an list: x , y and z direction in meters [m]

    # Type of Calculation
    #Choose "decay" if the objective is to calculate the energy decay of the room with all its energetic parameters; 
    #Choose "stationarysource" if the aim is to understand the behaviour of a room subject to a stationary source
    tcalc = "decay"

    # Frequency range
    fc_low = 125
    fc_high = 2000
    num_octave = 1    

    # Time discretization
    dt = 1/20000 #time discretization

    # Air absorption coefficient
    m_atm = 0 #air absorption coefficient [1/m]
    
    # Absorption term and Absorption coefficients
    th = 3 #int(input("Enter type Absortion conditions (option 1,2,3):")) 
    # options Sabine (th=1), Eyring (th=2) and modified by Xiang (th=3)

    ###############################################################################
    #FIXED INPUTS
    ###############################################################################
    #General settings
    c0= 343 #adiabatic speed of sound [m.s^-1]

    #Set initial condition - Source Info (interrupted method)
    Ws = 0.01 #Source point power [Watts] interrupted after "sourceon_time" seconds; 10^-2 W => correspondent to 100dB

    # Reference pressure in Pa
    pRef = 2 * (10**-5) #Reference pressure in Pa

    # Air density
    rho = 1.21 #air density [kg.m^-3] at 20°C

    #%%
    ###############################################################################
    #INITIALISE GMSH
    ###############################################################################
    
    
    ####?????? Ilaria: Now, this is a bit more complicated. The file name is an input for sure but then gmsh needs to open it and initialise it with the lines
    
    file_name = "3x3x3.msh" #Insert file name, msh file created from sketchUp and then gmsh

    ####?????? Ilaria: What is this?
    script_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(script_dir, file_name) # Full path to the file
    ####??????

    ####?????? Ilaria: here it initialise the file. This is also important because the gmsh.exe file needs to be in the same folder as the code!!!
    gmsh.initialize() #Initialize msh file
    mesh = gmsh.open(file_path) #open the file
    #gmsh.fltk.run() #run the file to see it in gmsh
    dim = -1 #dimensions of the entities, 0 for points, 1 for curves/edge/lines, 2 for surfaces, 3 for volumes, -1 for all the entities 
    tag = -1 #all the nodes of the room
    ####??????
    
    ####?????? Ilaria: Call of the first function. This needs to stay outside the run_sim function because then it is used for the input of the absorption coefficient
    #Calling function %create_vgroups_names%
    vGroupsNames = create_vgroups_names()

    # Initialize a list to store surface tags and their absorption coefficients
    surface_absorption = [] #initialization absorption term (alpha*surfaceofwall) for each wall of the room
    triangle_face_absorption = [] #initialization absorption term for each triangle face at the boundary and per each wall
    absorption_coefficient_dict = {}

    for group in vGroupsNames:
        if group[0] != 2:
            continue

        name_group = group[2]
        name_split = name_group.split("$")
        name_abs_coeff = name_split[0]

        abscoeff = input(f"Enter absorption coefficient for frequency {fc_low} to {fc_high} for {name_abs_coeff}:") 
        surface_materials (group, abscoeff, surface_absorption, absorption_coefficient_dict)

    for entity, Abs_term in surface_absorption:
        triangle_faces, _ = gmsh.model.mesh.getElementsByType(2, entity) #Get all the triangle faces for the current surface
        triangle_face_absorption.extend([Abs_term] * len(triangle_faces)) #Append the Abs_term value for each triangle face

    print("Correctly inputted surface materials. Starting initial geometry calculations...")


    def run_sim(coord_source, coord_rec,fc_low,fc_high,num_octave, dt,m_atm , c0, Ws, th, pRef, rho, file_name, tcalc = "decay"):
        
        #Calling function %number_freq%
        x_frequencies, nBands,center_freq = number_freq()
        
        #Calling function %get_nodes_elem%
        nodecoords, node_indices, bounEl, bounNode, voluEl, voluNode, belemNodes, velemNodes, boundaryEl_dict, volumeEl_dict = get_nodes_elem()
        
        #Calling function %velem_volume_centre%
        cell_center, cell_volume = velem_volume_centre()
    
        #Calling function %belem_area_centre%
        barea_dict, centre_area = belem_area_centre()
        
        #Calling function %get_neighbour_faces%
        fxt, txt, neighbourVolume = get_neighbour_faces()
        print("Completed initial geometry calculation. Starting internal tetrahedrons calculations...")
        
        #Calling function %interior_tetra%
        interior_tet, interior_tet_sum = interior_tetra()
        print("Completed internal tetrahedrons calculation. Starting boundary tetrahedrons calculations...")
        
        #Calling function %surface_absorption%
        surface_areas = surface_area(surface_absorption)
        
        #Calling function %boundary_triang%
        boundary_areas, total_boundArea = boundary_triang(triangle_face_absorption)
        print("Completed boundary tetrahedrons calculation. Starting main diffusion equation calculations over time and frequency...")
        
        ####?????? Ilaria: Here the gmsh needs to be finialised. The operations with the mesh have finished and from now we pass to the main calcualtions of the diffusion equation method.
        gmsh.finalize()
        ####??????
        
        #Calling function %equiv_absorp%
        V,S,Eq_A = equiv_absorp(surface_areas)
        
        #Calling function %calculate_sourceon_time%
        sourceon_time = calculate_sourceon_time()
        
        #Calling function %rec_time%
        recording_time, t, recording_steps = rec_time(sourceon_time, dt)
        
        #Calling function %diff_coeff%
        Dx, Dy, Dz = diff_coeff()
        
        #Calling function %dist_source_receiver%
        dist_sr = dist_source_receiver()
        
        #Calling function %source_interp%
        cl_tet_s_keys, total_weights_s = source_interp()
        
        #Calling function %source_volume%
        Vs = source_volume()
        
        #Calling function %initial_cond%
        source1, sourceon_steps = initial_cond()
        
        #Calling function %source_matrix%
        s = source_matrix()
        
        #Calling function %receiver_interp%
        cl_tet_r_keys, total_weights_r = receiver_interp()
        
        #Calling function %room_dimensions%
        room_length, room_width, room_height = room_dimensions()
        
        #Calling function %line_receivers%
        x_axis, y_axis, line_rec_x_idx_list, dist_x, line_rec_y_idx_list, dist_y = line_receivers()
        
        #Calling function %beta_zero%
        beta_zero_freq = beta_zero()
        
        #Calling function %computing_energy_density%
        w_new_band, w_rec_band, w_rec_off_band, w_rec_off_deriv_band, p_rec_off_deriv_band, idx_w_rec, t_off = computing_energy_density()
        print("100% of main calculation completed")
        
        #Calling function %freq_parameters%
        w_rec_x_band, w_rec_y_band, spl_stat_x_band, spl_stat_y_band, spl_r_band, spl_r_off_band, spl_r_norm_band, sch_db_band, t30_band, edt_band, c80_band, d50_band, ts_band = freq_parameters()
        
        ####?????? Ilaria: What do we do with this?
        et = time.time() #end time
        elapsed_time = et - st
        ####??????
    
        #Calling function %save%
        save('resultsFVM.pkl')
