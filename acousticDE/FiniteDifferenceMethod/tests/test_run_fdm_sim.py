# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 16:12:40 2023

@author: Ilaria Fichera
"""

import time

from acousticDE.FiniteDifferenceMethod.FunctionClarity import clarity
from acousticDE.FiniteDifferenceMethod.FunctionDefinition import definition
from acousticDE.FiniteDifferenceMethod.FunctionCentreTime import centretime
from acousticDE.FiniteDifferenceMethod.FunctionRT import t60_decay
from acousticDE.FiniteDifferenceMethod.FDMfunctions import *

#%%
###############################################################################
#RUN SIMULATION FUNCTION
###############################################################################

def test_run_fdm_sim():
    """
    Function for running the full calculation. It will use all the functions defined in FDMfunctions.py file

    Parameters
    ----------
    inputs : dict
        Dictionary of all the inputs for the simulation (taken from the json file).
    
    Returns
    -------
    results : dict
        Dictionary of all the variable calculated during the simulation
    """
    
    inputs = {
        "room_dim": [3.0, 3.0, 3.0],
        "coord_source": [1.0, 1.0, 1.0], #source coordinates x,y,z
        "coord_rec": [2.0, 2.0, 2.0], #rec coordinates x,y,z
        "alpha_1": [0.80, 0.82, 0.84, 0.86, 0.88, 0.90], #Absorption coefficient for Surface1 - Floor
        "alpha_2": [0.80, 0.82, 0.84, 0.86, 0.88, 0.90], #Absorption coefficient for Surface2 - Ceiling
        "alpha_3": [0.80, 0.82, 0.84, 0.86, 0.88, 0.90], #Absorption coefficient for Surface3 - Wall Front
        "alpha_4": [0.80, 0.82, 0.84, 0.86, 0.88, 0.90], #Absorption coefficient for Surface4 - Wall Back
        "alpha_5": [0.80, 0.82, 0.84, 0.86, 0.88, 0.90], #Absorption coefficient for Surface5 - Wall Left
        "alpha_6": [0.80, 0.82, 0.84, 0.86, 0.88, 0.90], #Absorption coefficient for Surface6 - Wall Right
        "fc_low": 125, #lowest frequency
        "fc_high": 4000, #highest frequency
        "num_octave": 1, # 1 or 3 depending on how many octave you want
        "dx": 0.5,
        "dt": 1/8000, #time discretization
        "m_atm": 0, #air absorption coefficient [1/m]
        "th": 3, #int(input("Enter type Absortion conditions (option 1,2,3):")) # options Sabine (th=1), Eyring (th=2) and modified by Xiang (th=3)
        "tcalc": "decay" #Choose "decay" if the objective is to calculate the energy decay of the room with all its energetic parameters; Choose "stationarysource" if the aim is to understand the behaviour of a room subject to a stationary source
    }
    
    st = time.time() #start time of calculation
    
    # Access input variables
    length = inputs["room_dim"][0]
    width = inputs["room_dim"][1]
    height = inputs["room_dim"][2]
    coord_source = inputs["coord_source"]
    coord_rec = inputs["coord_rec"]
    alpha_1 = inputs["alpha_1"]
    alpha_2 = inputs["alpha_2"]
    alpha_3 = inputs["alpha_3"]
    alpha_4 = inputs["alpha_4"]
    alpha_5 = inputs["alpha_5"]
    alpha_6 = inputs["alpha_6"]
    fc_low = inputs["fc_low"]
    fc_high = inputs["fc_high"]
    num_octave = inputs["num_octave"]
    dx = inputs["dx"]
    dy = inputs["dx"]
    dz = inputs["dx"]
    dt = inputs["dt"]
    m_atm = inputs["m_atm"]
    th = inputs["th"]
    tcalc = inputs.get("tcalc", "decay")  # default fallback
    
    #Fixed inputs
    #Set initial condition - Source Info (interrupted method)
    c0= 343 #adiabatic speed of sound [m.s^-1]
    Ws = 0.01 #Source point power [Watts] interrupted after "sourceon_time" seconds; 10^-2 W => correspondent to 100dB
    pRef = 2 * (10**-5) #Reference pressure in Pa
    rho = 1.21 #air density [kg.m^-3] at 20Â°C

    nBands, center_freq = number_freq(num_octave,fc_high,fc_low)
    
    #Calling function %room_charact%
    S1, S2, S3, S4, S5, S6, S, V = room_charact(length,width,height)

    #Calling function %create_mesh%
    x,y,z,Nx,Ny,Nz,xx,yy,zz = create_mesh(length,width,height, dx, dy, dz)

    #Calling function %abs_term%
    Abs_1 = abs_term(th,alpha_1,c0) #absorption term for S1
    Abs_2 = abs_term(th,alpha_2,c0) #absorption term for S2
    Abs_3 = abs_term(th,alpha_3,c0) #absorption term for S3
    Abs_4 = abs_term(th,alpha_4,c0) #absorption term for S4
    Abs_5 = abs_term(th,alpha_5,c0) #absorption term for S5
    Abs_6 = abs_term(th,alpha_6,c0) #absorption term for S6
    
    #Calling function %equiv_absorp%
    alpha_average, Eq_A = equiv_absorp(alpha_1, alpha_2, alpha_3, alpha_4, alpha_5, alpha_6, S1, S2, S3, S4, S5, S6, S)

    #Calling function %rec_sourceon_time%
    RT_Sabine_band, sourceon_time, recording_time, t, recording_steps, sourceon_steps = rec_sourceon_time(nBands, V, Eq_A, dt)

    #Calling function %diff_coeff%
    D_th, Dx, Dy, Dz = diff_coeff(V, S, c0)
    
    #Calling function %beta_zero%
    beta_zero_x, beta_zero_y, beta_zero_z, beta_zero, beta_zero_condition = beta_zero_fun(Dx, Dy, Dz, dx, dy, dz, dt)
        
    #Calling function %initial_cond%
    Vs, w1, s1, source1 = initial_condition(dx, dy, dz, Ws, sourceon_steps, recording_steps)
    
    #Calling function %source_interp%
    s, row_lr_s, row_up_s, col_lr_s, col_up_s, dep_lr_s, dep_up_s, weight_row_lr_s, weight_row_up_s, weight_col_lr_s, weight_col_up_s, weight_dep_lr_s, weight_dep_up_s = source_interpolation(coord_source, dx, dy, dz, source1, Nx,Ny,Nz)
    
    #Calling function %rec_interp%
    row_lr_r, row_up_r, col_lr_r, col_up_r, dep_lr_r, dep_up_r, weight_row_lr_r, weight_row_up_r, weight_col_lr_r, weight_col_up_r, weight_dep_lr_r, weight_dep_up_r = receiver_interpolation(coord_rec, dx, dy, dz)
    
    #Calling function %dist_source_receiver%
    dist_sr = dist_source_receiver(coord_rec, coord_source)
    
    #Calling function %dist_source_x_y%
    dist_x, dist_y = dist_source_x_y(xx, yy, zz, row_lr_r, row_up_r, col_lr_r, col_up_r, dep_lr_r, dep_up_r, weight_row_lr_r, weight_row_up_r, weight_col_lr_r, weight_col_up_r, weight_dep_lr_r, weight_dep_up_r, coord_source)
    
    #Calling function %computing_energy_density%
    w_new_band, w_t0_band, w_rec_band, w_rec_off_band, w_rec_x_t0_band, idx_w_rec = comput_energy_density(nBands, c0, m_atm, Nx, Ny, Nz, recording_steps, x, y, z, recording_time, dt, beta_zero, 
                                 beta_zero_x , beta_zero_y, beta_zero_z, dx, dy, dz, Dx, Dy, Dz, t, sourceon_time, sourceon_steps,
                                 Abs_1, Abs_2, Abs_3, Abs_4, Abs_5, Abs_6, s, source1, 
                                 row_lr_s, row_up_s, col_lr_s, col_up_s, dep_lr_s, dep_up_s, 
                                 weight_row_lr_s, weight_row_up_s, weight_col_lr_s, weight_col_up_s, weight_dep_lr_s, weight_dep_up_s, 
                                 row_lr_r, row_up_r, col_lr_r, col_up_r, dep_lr_r, dep_up_r, 
                                 weight_row_lr_r, weight_row_up_r, weight_col_lr_r, weight_col_up_r, weight_dep_lr_r, weight_dep_up_r, tcalc)
    
    print("Post-processing calculations...")
    #Calling function %freq_parameters%
    spl_r_band, spl_r_off_band, spl_r_norm_band, sch_db_band, spl_t0_band, spl_new_band, spl_rec_x_t0_band, t30_band, edt_band, c80_band, d50_band, ts_band = freq_param(nBands, 
                                                                                                                                                                              rho, pRef, c0, Ws, 
                                                                                                                                                                              w_new_band, w_t0_band, w_rec_band, 
                                                                                                                                                                              w_rec_off_band, w_rec_x_t0_band, idx_w_rec, 
                                                                                                                                                                              t, dist_sr, m_atm, V, S, Eq_A, tcalc)
        
    et = time.time() #end time
    elapsed_time = et - st
    
    results = locals()
    
    return results

results = test_run_fdm_sim()  

