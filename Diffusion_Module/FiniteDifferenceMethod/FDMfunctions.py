# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 10:39:42 2023

@author: 20225533
"""
#Code developed by Ilaria Fichera for the analysis of the FDM method Du Fort & Frankel solving the 3D diffusion equation with one intermittent omnidirectional sound source
import math
import time
from math import ceil
from math import log

import numpy as np
#uncomment this if you need drawnow
#from drawnow import drawnow

#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import pickle
import types

from Diffusion_Module.FiniteDifferentMethod.FunctionClarity import clarity
from Diffusion_Module.FiniteDifferentMethod.FunctionDefinition import definition
from Diffusion_Module.FiniteDifferentMethod.FunctionCentreTime import centretime
from Diffusion_Module.FiniteDifferentMethod.FunctionRT import t60_decay

#%%
###############################################################################
#CALCULATION SECTION
###############################################################################

def number_freq(num_octave,fc_high,fc_low):
    """
    Calculate the frequency array and the number of frequency bands.

    Parameters
    ----------
    num_octave : int
        Number of octaves to calculate (1 or 3 octaves).
    fc_high : int
        The highest frequency in the calculation.
    fc_low : int
        The lowest frequency in the calculation.

    Returns
    -------
    nBands : int
        Number of frequency bands.
    center_freq : list of float
        Array of all the frequencies to calculate.
    """
    x_frequencies  = num_octave * log(fc_high/fc_low) / log(2)
    nBands = int(num_octave * log(fc_high/fc_low) / log(2) + 1)
    center_freq = fc_low * np.power(2,((np.arange(0,x_frequencies+1) / num_octave)))
    return nBands, center_freq



#Room characteristics
def room_charact(length,width,height):
    S1,S2 = length*width, length*width #xy planes
    S3,S4 = length*height, length*height #xz planes
    S5,S6 = width*height, width*height #yz planes
    S = length*width*2 + length*height*2 + width*height*2 #Total Surface Area[m2]
    V = length*width*height #Volume of the room [m^3]
    return S1, S2, S3, S4, S5, S6, S, V


#Creating of meshgrid
def create_mesh(length,width,height, dx, dy, dz):
    x = np.arange(0, length+dx, dx) #mesh points in space x direction
    y = np.arange(0, width+dy, dy) #mesh points in space y direction
    z = np.arange(0, height+dz, dz) #mesh points in space z direction
    Nx = len(x) #number of point in the x direction
    Ny = len(y) #number of point in the y direction
    Nz = len(z) #number of point in the z direction
    yy, xx , zz = np.meshgrid(y,x,z) #Return coordinate matrices from coordinate vectors; create the 3D grid
    return x,y,z,Nx,Ny,Nz,xx,yy,zz

#Absorption term for boundary conditions 
def abs_term(th,alpha,c0):
    Absx_array = np.array([])
    for abs_coeff in alpha:
        #print(abs_coeff)
        if th == 1:
            Absx = (c0*abs_coeff)/4 #Sabine
        elif th == 2:
            Absx = (c0*(-log(1-abs_coeff)))/4 #Eyring
        elif th == 3:
            Absx = (c0*abs_coeff)/(2*(2-abs_coeff)) #Modified by Xiang
        Absx_array = np.append(Absx_array, Absx)
    return Absx_array


#%%
###############################################################################
#CALCULATION OF EQUIVALENT ABSORPTION AREA
###############################################################################
#Absorption parameters for room
def equiv_absorp(alpha_1, alpha_2, alpha_3, alpha_4, alpha_5, alpha_6, S1, S2, S3, S4, S5, S6, S):
    alpha_average = (np.multiply(alpha_1,S1) + np.multiply(alpha_2,S2) + np.multiply(alpha_3,S3) + np.multiply(alpha_4,S4) + np.multiply(alpha_5,S5) + np.multiply(alpha_6,S6))/S #average absorption
    Eq_A = np.multiply(alpha_1,S1) + np.multiply(alpha_2,S2) + np.multiply(alpha_3,S3) + np.multiply(alpha_4,S4) + np.multiply(alpha_5,S5) + np.multiply(alpha_6,S6) #equivalent absorption area of the room
    return alpha_average, Eq_A


def rec_sourceon_time(nBands, center_freq, V, Eq_A, dt):
    RT_Sabine_band = []
    for iBand in range(nBands):
        #freq = center_freq[iBand]
        RT_Sabine = 0.16*V/Eq_A[iBand]
        RT_Sabine_band.append(RT_Sabine)

    sourceon_time = round(max(RT_Sabine_band),1)#time that the source is ON before interrupting [s]
    recording_time = 2*sourceon_time #total time recorded for the calculation [s]
    #Time resolution
    t = np.arange(0, recording_time, dt) #mesh point in time
    recording_steps = ceil(recording_time/dt) #number of time steps to consider in the calculation
    sourceon_steps = ceil(sourceon_time/dt) #time steps at which the source is calculated/considered in the calculation
    
    return RT_Sabine_band, sourceon_time, recording_time, t, recording_steps, sourceon_steps


def diff_coeff(V, S, c0):
    """
    Calculation of the theoretical diffusion coefficient

    Parameters
    ----------
        V : float
            Volume of the room
        S : float
            Total surface are of the room
        c0 : int
            Speed of sound 

    Returns
    -------
        Dx : float 
            Diffusion coefficient in the x direction (equal to the theoretical diffusion coefficient)
        Dy : float 
            Diffusion coefficient in the y direction (equal to the theoretical diffusion coefficient)
        Dz : float
            Diffusion coefficient in the z direction (equal to the theoretical diffusion coefficient)
    """
    mean_free_path = (4*V)/S #mean free path for 3D
    #mean_free_time = mean_free_path/c0 #mean free time for 3D
    #mean_free_time_step = int(mean_free_time/dt)
    D_th = (mean_free_path*c0)/3
    Dx = (mean_free_path*c0)/3 #diffusion coefficient for proportionate rooms x direction
    Dy = (mean_free_path*c0)/3 #diffusion coefficient for proportionate rooms y direction
    Dz = (mean_free_path*c0)/3 #diffusion coefficient for proportionate rooms z direction
    return D_th, Dx, Dy, Dz


#Mesh numbers
def beta_zero(Dx, Dy, Dz, dx, dy, dz, dt):
    beta_zero_x = (2*Dx*dt)/(dx**2) #mesh number in x direction
    beta_zero_y = (2*Dy*dt)/(dy**2) #mesh number in x direction
    beta_zero_z = (2*Dz*dt)/(dz**2) #mesh number in x direction
    beta_zero = beta_zero_x + beta_zero_y + beta_zero_z #beta_zero is the condition for all the directions deltax, deltay and deltaz.
    #Condition for the model to be unconditionally stable
    beta_zero_condition = ((beta_zero**2) - 1)/(1+(beta_zero**2)+(2*beta_zero)) #from Navarro 2012 paper
    if beta_zero_condition >1:
        print("aa! error! Check beta condition")
    return beta_zero_x, beta_zero_y, beta_zero_z, beta_zero, beta_zero_condition
 
    
#Initial condition - Source Info (interrupted method)
def initial_cond(dx, dy, dz, Ws, sourceon_steps, recording_steps):
    Vs=dx*dy*dz  #Volume of the source
    w1=Ws/Vs #w1 = round(Ws/Vs,4) #power density of the source [Watts/(m^3))]
    s1 = np.multiply(w1,np.ones(sourceon_steps)) #energy density of source number 1 at each time step position
    source1 = np.append(s1, np.zeros(recording_steps-sourceon_steps)) #This would be equal to s1 if and only if recoding_steps = sourceon_steps
    return Vs, w1, s1, source1


###############################################################################
#SOURCE INTERPOLATION
###############################################################################
def source_interp(coord_source, dx, dy, dz, source1, Nx,Ny,Nz):
    # Calculate the fractional indices
    row_lr_s = int(np.floor(coord_source[0] / dx))
    row_up_s = row_lr_s + 1
    col_lr_s = int(np.floor(coord_source[1] / dy))
    col_up_s = col_lr_s + 1
    dep_lr_s = int(np.floor(coord_source[2] / dz))
    dep_up_s = dep_lr_s + 1
    
    # Calculate the interpolation weights
    weight_row_up_s = (coord_source[0] / dx) - row_lr_s
    weight_row_lr_s = 1 - weight_row_up_s
    weight_col_up_s = (coord_source[1] / dy) - col_lr_s
    weight_col_lr_s = 1 - weight_col_up_s
    weight_dep_up_s = (coord_source[2] / dz) - dep_lr_s
    weight_dep_lr_s = 1 - weight_dep_up_s
    
    s = np.zeros((Nx,Ny,Nz)) #matrix of zeros for source
    
    # Perform linear interpolation
    s[row_lr_s, col_lr_s, dep_lr_s] += source1[0] * weight_row_lr_s * weight_col_lr_s * weight_dep_lr_s
    s[row_lr_s, col_lr_s, dep_up_s] += source1[0] * weight_row_lr_s * weight_col_lr_s * weight_dep_up_s
    s[row_lr_s, col_up_s, dep_lr_s] += source1[0] * weight_row_lr_s * weight_col_up_s * weight_dep_lr_s
    s[row_lr_s, col_up_s, dep_up_s] += source1[0] * weight_row_lr_s * weight_col_up_s * weight_dep_up_s
    s[row_up_s, col_lr_s, dep_lr_s] += source1[0] * weight_row_up_s * weight_col_lr_s * weight_dep_lr_s
    s[row_up_s, col_lr_s, dep_up_s] += source1[0] * weight_row_up_s * weight_col_lr_s * weight_dep_up_s
    s[row_up_s, col_up_s, dep_lr_s] += source1[0] * weight_row_up_s * weight_col_up_s * weight_dep_lr_s
    s[row_up_s, col_up_s, dep_up_s] += source1[0] * weight_row_up_s * weight_col_up_s * weight_dep_up_s

    return s, row_lr_s, row_up_s, col_lr_s, col_up_s, dep_lr_s, dep_up_s, weight_row_lr_s, weight_row_up_s, weight_col_lr_s, weight_col_up_s, weight_dep_lr_s, weight_dep_up_s

s, row_lr_s, row_up_s, col_lr_s, col_up_s, dep_lr_s, dep_up_s, weight_row_lr_s, weight_row_up_s, weight_col_lr_s, weight_col_up_s, weight_dep_lr_s, weight_dep_up_s = source_interp(coord_source, dx, dy, dz, source1)


###############################################################################
#RECEIVER INTERPOLATION
###############################################################################

def rec_interp(coord_rec, dx, dy, dz):
    #Calculate the fractional indices for receiver
    row_lr_r = int(np.floor(coord_rec[0] / dx))
    row_up_r = row_lr_r + 1
    col_lr_r = int(np.floor(coord_rec[1] / dy))
    col_up_r = col_lr_r + 1
    dep_lr_r = int(np.floor(coord_rec[2] / dz))
    dep_up_r = dep_lr_r + 1
       
    #Calculate the interpolation weights for receiver
    weight_row_up_r = (coord_rec[0] / dx) - row_lr_r #weight x upper
    weight_row_lr_r = 1 - weight_row_up_r #weight x lower
    weight_col_up_r = (coord_rec[1] / dy) - col_lr_r #weight y upper
    weight_col_lr_r = 1 - weight_col_up_r #weight y lower
    weight_dep_up_r = (coord_rec[2] / dz) - dep_lr_r #weight z upper
    weight_dep_lr_r = 1 - weight_dep_up_r #weight z lower
    
    return row_lr_r, row_up_r, col_lr_r, col_up_r, dep_lr_r, dep_up_r, weight_row_lr_r, weight_row_up_r, weight_col_lr_r, weight_col_up_r, weight_dep_lr_r, weight_dep_up_r

row_lr_r, row_up_r, col_lr_r, col_up_r, dep_lr_r, dep_up_r, weight_row_lr_r, weight_row_up_r, weight_col_lr_r, weight_col_up_r, weight_dep_lr_r, weight_dep_up_r = rec_interp(coord_rec, dx, dy, dz)


def dist_source_receiver(coord_rec, coord_source):
    """
    Calculation of the distance between source and receiver positions

    Parameters
    ----------
        coord_rec : list 
            Coordinates of the receiver position
        coord_source : list
            Coordinates of the source position

    Returns
    -------
        dist_sr : float
            Distance between source and receiver position
    """
    #distance between source and receiver
    dist_sr = math.sqrt((abs(coord_rec[0] - coord_source[0]))**2 + (abs(coord_rec[1] - coord_source[1]))**2 + (abs(coord_rec[2] - coord_source[2]))**2) #distance between source and receiver
    return dist_sr

dist_sr = dist_source_receiver(coord_rec, coord_source)



def dist_source_x_y(xx, yy, zz, row_lr_r, row_up_r, col_lr_r, col_up_r, dep_lr_r, dep_up_r, weight_row_lr_r, weight_row_up_r, weight_col_lr_r, weight_col_up_r, weight_dep_lr_r, weight_dep_up_r, coord_source):
    #distance between source and each mesh point in the x direction 
    dist_x = np.sqrt((((xx[:, col_lr_r, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_lr_r))+\
        (xx[:, col_lr_r, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
            (xx[:, col_up_r, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
                (xx[:, col_up_r, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                    (xx[:, col_lr_r, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                        (xx[:, col_lr_r, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                            (xx[:, col_up_r, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                                (xx[:, col_up_r, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r))) - coord_source[0])**2 +\
                     (((yy[row_lr_r, col_lr_r, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r* weight_dep_lr_r))+\
                         (yy[row_lr_r, col_lr_r, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
                             (yy[row_lr_r, col_up_r, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
                                 (yy[row_lr_r, col_up_r, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                                     (yy[row_up_r, col_lr_r, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                                         (yy[row_up_r, col_lr_r, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                                             (yy[row_up_r, col_up_r, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                                                 (yy[row_up_r, col_up_r, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r))) - coord_source[1])**2 +\
                         (((zz[row_lr_r, col_lr_r, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_lr_r))+\
                             (zz[row_lr_r, col_lr_r, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
                                 (zz[row_lr_r, col_up_r, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
                                     (zz[row_lr_r, col_up_r, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                                         (zz[row_up_r, col_lr_r, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                                             (zz[row_up_r, col_lr_r, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                                                 (zz[row_up_r, col_up_r, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                                                     (zz[row_up_r, col_up_r, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r))) - coord_source[2])**2)
    
    
    #distance between source and each mesh point in the y direction  
    dist_y = np.sqrt((((xx[row_lr_r, col_lr_r, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_lr_r))+\
        (xx[row_lr_r, col_lr_r, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
            (xx[row_lr_r, col_up_r, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
                (xx[row_lr_r, col_up_r, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                    (xx[row_up_r, col_lr_r, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                        (xx[row_up_r, col_lr_r, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                            (xx[row_up_r, col_up_r, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                                (xx[row_up_r, col_up_r, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r))) - coord_source[0])**2 +\
                     (((yy[row_lr_r, :, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r* weight_dep_lr_r))+\
                         (yy[row_lr_r, :, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
                             (yy[row_lr_r, :, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
                                 (yy[row_lr_r, :, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                                     (yy[row_up_r, :, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                                         (yy[row_up_r, :, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                                             (yy[row_up_r, :, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                                                 (yy[row_up_r, :, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r))) - coord_source[1])**2 +\
                         (((zz[row_lr_r, col_lr_r, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_lr_r))+\
                             (zz[row_lr_r, col_lr_r, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
                                 (zz[row_lr_r, col_up_r, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
                                     (zz[row_lr_r, col_up_r, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                                         (zz[row_up_r, col_lr_r, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                                             (zz[row_up_r, col_lr_r, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                                                 (zz[row_up_r, col_up_r, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                                                     (zz[row_up_r, col_up_r, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r))) - coord_source[2])**2)

    return dist_x, dist_y

dist_x, dist_y = dist_source_x_y(xx, yy, zz, row_lr_r, row_up_r, col_lr_r, col_up_r, dep_lr_r, dep_up_r, weight_row_lr_r, weight_row_up_r, weight_col_lr_r, weight_col_up_r, weight_dep_lr_r, weight_dep_up_r, coord_source)

#%%
###############################################################################
#MAIN CALCULATION - COMPUTING ENERGY DENSITY
############################################################################### 

def computing_energy_density(nBands, Nx, Ny, Nz, recording_steps, x, y, z, center_freq, recording_time, dt, beta_zero, 
                             Abs_1, Abs_2, Abs_3, Abs_4, Abs_5, Abs_6, s, source1, 
                             row_lr_s, row_up_s, col_lr_s, col_up_s, dep_lr_s, dep_up_s, 
                             weight_row_lr_s, weight_row_up_s, weight_col_lr_s, weight_col_up_s, weight_dep_lr_s, weight_dep_up_s, 
                             row_lr_r, row_up_r, col_lr_r, col_up_r, dep_lr_r, dep_up_r, 
                             weight_row_lr_r, weight_row_up_r, weight_col_lr_r, weight_col_up_r, weight_dep_lr_r, weight_dep_up_r, 
                             tcalc = "decay"):
    w_new_band = np.zeros((nBands, Nx, Ny, Nz))
    w_t0_band = np.zeros((nBands, Nx, Ny, Nz))
    w_rec_band = np.zeros((nBands, recording_steps))
    w_rec_off_band = []
    w_rec_x_t0_band = np.zeros((nBands, len(x)))
    
    curPercent = 0
    
    #Computing w;
    for iBand in range(nBands):
        #freq = center_freq[iBand]  
        
        w_new = np.zeros((Nx,Ny,Nz)) #unknown w at new time level (n+1)
        w = w_new.copy() #w at n level
        w_old = w.copy() #w_old at n-1 level
        # w = w_new #w at n level
        # w_old = w #w_old at n-1 level
    
        w_rec = np.arange(0,recording_time,dt) #energy density at the receiver
        w_rec_x_t0 = np.zeros((1,len(x))) 
        
        for steps in range(0, recording_steps):
            #Compute w at inner mesh points
            time_steps = steps*dt #total time for the calculation
            
            #In the x direction
            w_iminus1 = w[0:Nx-1, 0:Ny, 0:Nz]
            w_m_i = (w_iminus1[0,:,:])
            w_m_i = np.expand_dims(w_m_i, axis=0) #Expand the dimensions of w_m_i to match the shape of w_iminus1
            w_iminus1 = np.concatenate((w_m_i,w_iminus1),axis = 0)
        
            w_iplus1 = w[1:Nx, 0:Ny, 0:Nz]
            w_p_i = (w_iplus1[-1,:,:])
            w_p_i = np.expand_dims(w_p_i, axis=0) #Expand the dimensions of w_p_i to match the shape of w_iplus1
            w_iplus1 = np.concatenate((w_iplus1,w_p_i), axis=0)
            
            #In the y direction
            w_jminus1 = w[0:Nx, 0:Ny-1, 0:Nz]
            w_m_j = (w_jminus1[:,0,:])
            w_m_j = np.expand_dims(w_m_j, axis=1) #Expand the dimensions of w_m_j to match the shape of w_jminus1
            w_jminus1 = np.concatenate((w_m_j, w_jminus1), axis=1)
            
            w_jplus1 = w[0:Nx, 1:Ny, 0:Nz]
            w_p_j = (w_jplus1[:,-1,:])
            w_p_j = np.expand_dims(w_p_j, axis=1) #Expand the dimensions of w_p_j to match the shape of w_jplus1
            w_jplus1 = np.concatenate((w_jplus1,w_p_j), axis=1)
            
            #In the z direction
            w_kminus1 = w[0:Nx, 0:Ny, 0:Nz-1]
            w_m_k = (w_kminus1[:,:,0])
            w_m_k = np.expand_dims(w_m_k, axis=2) # Expand the dimensions of w_m_k to match the shape of w_kminus1
            w_kminus1 = np.concatenate((w_m_k, w_kminus1), axis=2)
            
            w_kplus1 = w[0:Nx, 0:Ny, 1:Nz]
            w_p_k = (w_kplus1[:,:,-1])
            w_p_k = np.expand_dims(w_p_k, axis=2) #Expand the dimensions of w_p_k to match the shape of w_kplus1
            w_kplus1 = np.concatenate((w_kplus1,w_p_k), axis=2)
            
            #Computing w_new (w at n+1 time step)
            w_new = np.divide((np.multiply(w_old,(1-beta_zero))),(1+beta_zero)) - \
                np.divide((2*dt*c0*m_atm*w),(1+beta_zero)) + \
                    np.divide((2*dt*s),(1+beta_zero)) + \
                        np.divide((np.multiply(beta_zero_x,(w_iplus1+w_iminus1))),(1+beta_zero)) + \
                            np.divide((np.multiply(beta_zero_y,(w_jplus1+w_jminus1))),(1+beta_zero)) + \
                                np.divide((np.multiply(beta_zero_z,(w_kplus1+w_kminus1))),(1+beta_zero))
            
            #Insert boundary conditions  
            w_new[0,:,:] = np.divide((4*w_new[1,:,:] - w_new[2,:,:]),(3+((2*Abs_5[iBand]*dx)/Dx))) #boundary condition at x=0, any y, any z
            w_new[-1,:,:] = np.divide((4*w_new[-2,:,:] - w_new[-3,:,:]),(3+((2*Abs_6[iBand]*dx)/Dx))) #boundary condition at lx=length, any y, any z
        
        
            w_new[:,0,:] = np.divide((4*w_new[:,1,:] - w_new[:,2,:]),(3+((2*Abs_3[iBand]*dx)/Dy))) #boundary condition at y=0, any x, any z
            w_new[:,-1,:] = np.divide((4*w_new[:,-2,:] - w_new[:,-3,:]),(3+((2*Abs_4[iBand]*dx)/Dy))) #boundary condition at at ly=width, any x, any z
         
            
            w_new[:,:,0] = np.divide((4*w_new[:,:,1] - w_new[:,:,2]),(3+((2*Abs_1[iBand]*dx)/Dz))) #boundary condition at z=0, any x, any y
            w_new[:,:,-1] = np.divide((4*w_new[:,:,-2] - w_new[:,:,-3]),(3+((2*Abs_2[iBand]*dx)/Dz))) #boundary condition at at lz=height, any x, any y
            
            #sdl = 10*np.log10(abs(w_new),where=abs(w_new)>0) #sound density level
            #spl = 10*np.log10(((abs(w_new))*rho*(c0**2))/(pRef**2)) #,where=press_r>0, sound pressure level in the 3D space
                  
            #Update w before next step
            w_old = w.copy() #The w at n step becomes the w at n-1 step
            w = w_new.copy() #The w at n+1 step becomes the w at n step
    
            # w_old = w #The w at n step becomes the w at n-1 step
            # w = w_new #The w at n+1 step becomes the w at n step
            
            #w_rec is the energy density at the specific receiver
            w_rec[steps] = ((w_new[row_lr_r, col_lr_r, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_lr_r))+\
                (w_new[row_lr_r, col_lr_r, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
                    (w_new[row_lr_r, col_up_r, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
                        (w_new[row_lr_r, col_up_r, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                            (w_new[row_up_r, col_lr_r, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                                (w_new[row_up_r, col_lr_r, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                                    (w_new[row_up_r, col_up_r, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                                        (w_new[row_up_r, col_up_r, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r)))
            
            idx_w_rec = np.argmin(np.abs(t - sourceon_time))
            w_rec_off = w_rec[idx_w_rec:]
            
            if steps == sourceon_steps:
                #print("Steps for source:",steps)
                w_t0 = w_new 
            
            #Updating the source term
            if tcalc == "decay":
                s[row_lr_s, col_lr_s, dep_lr_s] = source1[steps] * (weight_row_lr_s * weight_col_lr_s * weight_dep_lr_s)
                s[row_lr_s, col_lr_s, dep_up_s] = source1[steps] * (weight_row_lr_s * weight_col_lr_s * weight_dep_up_s)
                s[row_lr_s, col_up_s, dep_lr_s] = source1[steps] * (weight_row_lr_s * weight_col_up_s * weight_dep_lr_s)
                s[row_lr_s, col_up_s, dep_up_s] = source1[steps] * (weight_row_lr_s * weight_col_up_s * weight_dep_up_s)
                s[row_up_s, col_lr_s, dep_lr_s] = source1[steps] * (weight_row_up_s * weight_col_lr_s * weight_dep_lr_s)
                s[row_up_s, col_lr_s, dep_up_s] = source1[steps] * (weight_row_up_s * weight_col_lr_s * weight_dep_up_s)
                s[row_up_s, col_up_s, dep_lr_s] = source1[steps] * (weight_row_up_s * weight_col_up_s * weight_dep_lr_s)
                s[row_up_s, col_up_s, dep_up_s] = source1[steps] * (weight_row_up_s * weight_col_up_s * weight_dep_up_s)
            if tcalc == "stationarysource":
                s[row_lr_s, col_lr_s, dep_lr_s] = source1[0] * (weight_row_lr_s * weight_col_lr_s * weight_dep_lr_s)
                s[row_lr_s, col_lr_s, dep_up_s] = source1[0] * (weight_row_lr_s * weight_col_lr_s * weight_dep_up_s)
                s[row_lr_s, col_up_s, dep_lr_s] = source1[0] * (weight_row_lr_s * weight_col_up_s * weight_dep_lr_s)
                s[row_lr_s, col_up_s, dep_up_s] = source1[0] * (weight_row_lr_s * weight_col_up_s * weight_dep_up_s)
                s[row_up_s, col_lr_s, dep_lr_s] = source1[0] * (weight_row_up_s * weight_col_lr_s * weight_dep_lr_s)
                s[row_up_s, col_lr_s, dep_up_s] = source1[0] * (weight_row_up_s * weight_col_lr_s * weight_dep_up_s)
                s[row_up_s, col_up_s, dep_lr_s] = source1[0] * (weight_row_up_s * weight_col_up_s * weight_dep_lr_s)
                s[row_up_s, col_up_s, dep_up_s] = source1[0] * (weight_row_up_s * weight_col_up_s * weight_dep_up_s)
            
            
            #print(time_steps)
            percentDone = round(100*time_steps/recording_time);
            if (percentDone > curPercent):
                print(str(curPercent + 1) + "% done")
                curPercent += 1;
        
                    
        w_rec_x_t0 = ((w_t0[:, col_lr_r, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_lr_r))+\
            (w_t0[:, col_lr_r, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
                (w_t0[:, col_up_r, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
                    (w_t0[:, col_up_r, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                        (w_t0[:, col_lr_r, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                            (w_t0[:, col_lr_r, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                                (w_t0[:, col_up_r, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                                    (w_t0[:, col_up_r, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r)))         
        
        w_new_band[iBand] = w_new
        w_t0_band[iBand] = w_t0
        w_rec_band[iBand] = w_rec
        w_rec_off_band.append(w_rec_off) 
        w_rec_x_t0_band[iBand] = w_rec_x_t0
    
    plt.show()
    
    w_rec_x_end = ((w_new[:, col_lr_r, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_lr_r))+\
        (w_new[:, col_lr_r, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
            (w_new[:, col_up_r, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
                (w_new[:, col_up_r, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                    (w_new[:, col_lr_r, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                        (w_new[:, col_lr_r, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                            (w_new[:, col_up_r, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                                (w_new[:, col_up_r, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r)))
        
    w_rec_y_end = ((w_new[row_lr_r, :, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_lr_r))+\
        (w_new[row_lr_r, :, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
            (w_new[row_lr_r, :, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
                (w_new[row_lr_r, :, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                    (w_new[row_up_r, :, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                        (w_new[row_up_r, :, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                            (w_new[row_up_r, :, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                                (w_new[row_up_r, :, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r)))

    return w_new_band, w_t0_band, w_rec_band, w_rec_off_band, w_rec_x_t0_band, idx_w_rec

w_new_band, w_t0_band, w_rec_band, w_rec_off_band, w_rec_x_t0_band, idx_w_rec = computing_energy_density(nBands, Nx, Ny, Nz, recording_steps, x, y, z, center_freq, recording_time, dt, beta_zero, 
                             Abs_1, Abs_2, Abs_3, Abs_4, Abs_5, Abs_6, s, source1, 
                             row_lr_s, row_up_s, col_lr_s, col_up_s, dep_lr_s, dep_up_s, 
                             weight_row_lr_s, weight_row_up_s, weight_col_lr_s, weight_col_up_s, weight_dep_lr_s, weight_dep_up_s, 
                             row_lr_r, row_up_r, col_lr_r, col_up_r, dep_lr_r, dep_up_r, 
                             weight_row_lr_r, weight_row_up_r, weight_col_lr_r, weight_col_up_r, weight_dep_lr_r, weight_dep_up_r, 
                             tcalc = "decay")
#%%

###############################################################################
#RESULTS
###############################################################################
def freq_parameters(nBands, rho, c0, w_new_band, w_t0_band, w_rec_band, w_rec_off_band, w_rec_x_t0_band, idx_w_rec, t, dist_sr, m_atm, V, S, Eq_A):

    spl_r_band = []
    spl_r_off_band = []
    spl_r_norm_band = []
    t30_band = []
    edt_band = []
    c80_band = []
    d50_band = []
    ts_band = []
    sch_db_band = []
    spl_new_band = []
    spl_t0_band = []
    spl_rec_x_t0_band = []
    
    
    for iBand in range(nBands):
     
        spl_t0 = 10*np.log10(((abs(w_t0_band[iBand]))*rho*(c0**2))/(pRef**2))
        spl_new = 10*np.log10(((abs(w_new_band[iBand]))*rho*(c0**2))/(pRef**2))
        
        press_r = ((abs(w_rec_band[iBand]))*rho*(c0**2)) #pressure at the receiver
        spl_r = 10*np.log10(((abs(w_rec_band[iBand]))*rho*(c0**2))/(pRef**2)) #,where=press_r>0, sound pressure level at the receiver
        spl_r_off = 10*np.log10(((abs(w_rec_off_band[iBand]))*rho*(c0**2))/(pRef**2))
        
        spl_r_norm = 10*np.log10((((abs(w_rec_band[iBand]))*rho*(c0**2))/(pRef**2)) / np.max(((abs(w_rec_band[iBand]))*rho*(c0**2))/(pRef**2))) #normalised to maximum to 0dB
        spl_r_tot = 10*np.log10(rho*c0*((Ws/(4*math.pi*dist_sr**2))*np.exp(-m_atm*dist_sr) + ((abs(w_rec_band[iBand]))*c0))/(pRef**2)) #spl total (including direct field) at the receiver position????? but it will need to be calculated for a stationary source 100dB
        
        spl_rec_x_t0 = 10*np.log10(rho*c0**2*abs(w_rec_x_t0_band[iBand])/pRef**2)
        
        #Schroeder integration
        schroeder = w_rec_off_band[iBand] #energy_r_rev_cum[::-1] #reverting the array again -> creating the schroder decay
        sch_db = 10.0 * np.log10(schroeder / max(schroeder)) #level of the array: schroeder decay
        
        if tcalc == "decay":
            t30 = t60_decay(t, sch_db, idx_w_rec, rt='t30') #called function for calculation of t60 [s]
            edt = t60_decay(t, sch_db, idx_w_rec, rt='edt') #called function for calculation of edt [s]
            c80 = clarity(t30, V, Eq_A[iBand], S, c0, dist_sr) #called function for calculation of c80 [dB]
            d50 = definition(t30, V, Eq_A[iBand], S, c0, dist_sr) #called function for calculation of d50 [%]
            ts = centretime(t30, Eq_A[iBand], S) #called function for calculation of ts [ms]
        
            t30_band.append(t30)
            edt_band.append(edt)
            c80_band.append(c80)
            d50_band.append(d50)
            ts_band.append(ts)
            
        spl_r_band.append(spl_r)
        spl_r_off_band.append(spl_r_off)
        spl_r_norm_band.append(spl_r_norm)
        sch_db_band.append(sch_db)
        spl_t0_band.append(spl_t0)
        spl_new_band.append(spl_new)
        spl_rec_x_t0_band.append(spl_rec_x_t0)
        
    spl_r_off_band = np.array(spl_r_off_band)
    t30_band = np.array(t30_band)
    edt_band = np.array(edt_band)
    c80_band = np.array(c80_band)
    d50_band = np.array(d50_band)
    ts_band = np.array(ts_band)
    
    return spl_r_band, spl_r_off_band, spl_r_norm_band, sch_db_band, spl_t0_band, spl_new_band, spl_rec_x_t0_band, t30_band, edt_band, c80_band, d50_band, ts_band

spl_r_band, spl_r_off_band, spl_r_norm_band, sch_db_band, spl_t0_band, spl_new_band, spl_rec_x_t0_band, t30_band, edt_band, c80_band, d50_band, ts_band = freq_parameters(nBands, rho, c0, w_new_band, w_t0_band, w_rec_band, w_rec_off_band, w_rec_x_t0_band, idx_w_rec, t, dist_sr, m_atm, V, S, Eq_A)
    
et = time.time() #end time
elapsed_time = et - st

#%%
###############################################################################
#FIGURES & POST-PROCESSING
###############################################################################

# if tcalc == "decay":
#     #Figure 5: Decay of SPL in the recording_time
#     plt.figure(5)
#     plt.plot(t, spl_r)  # plot sound pressure level with Pref = (2e-5)**5
#     plt.title("Figure 5 :SPL over time at the receiver")
#     plt.xlabel("t [s]")
#     plt.ylabel("SPL [dB]")
#     plt.xlim()
#     plt.ylim()
#     plt.xticks(np.arange(0, recording_time + 0.1, 0.5))
#     plt.yticks(np.arange(0, 120, 20))

#     #Figure 6: Decay of SPL in the recording_time normalised to maximum 0dB
#     plt.figure(6)
#     plt.plot(t,spl_r_norm)
#     plt.title("Figure 6: Normalised SPL over time at the receiver")
#     plt.xlabel("t [s]")
#     plt.ylabel("SPL [dB]")
#     plt.xlim()
#     plt.ylim()
#     plt.xticks(np.arange(0, recording_time +0.1, 0.1))
#     plt.yticks(np.arange(0, -60, -10))
    
#     #Figure 7: Energy density at the receiver over time
#     plt.figure(7)
#     plt.plot(t,w_rec)
#     plt.title("Figure 7: Energy density over time at the receiver")
#     plt.xlabel("t [s]")
#     plt.ylabel("Energy density [kg m^-1 s^-2]")
#     plt.xlim()
#     plt.ylim()
#     plt.xticks(np.arange(0, recording_time +0.1, 0.1))
    
#     #Figure 8: Schroeder decay
#     plt.figure(8)
#     plt.plot(t[idx_w_rec:],sch_db)
#     plt.title("Figure 8: Schroeder decay (Energy Decay Curve)")
#     plt.xlabel("t [s]")
#     plt.ylabel("Energy decay [dB]")
#     plt.xlim()
#     plt.ylim()
#     #plt.xticks(np.arange(t, recording_time +0.1, 0.1))
    
#     #Figure 9: 2D image of the energy density in the room
#     w_new_2d = w_new[:,:,dep_up_r] #The 3D w_new array is slised at the the desired z level
#     plt.figure(9)
#     plt.imshow(w_new_2d, origin='lower', extent=[x[0], x[-1], y[0], y[-1]], aspect='equal') #plot with the extent being the room dimension x and y 
#     plt.colorbar(label='Energy Density [kg/ms^2]')
#     plt.xlabel('X [m]')
#     plt.ylabel('Y [m]')
#     plt.title('Figure 9: Energy Density at Z = z_rec and t = recording_time')
#     plt.show()
    
#     #Figure 10: 2D image of the SDL in the room
#     sdl_2d = sdl[:,:,dep_up_r] #The 3D w_new array is slised at the the desired z level
#     plt.figure(10)
#     plt.imshow(sdl_2d, origin='lower', extent=[x[0], x[-1], y[0], y[-1]], aspect='equal') #plot with the extent being the room dimension x and y 
#     plt.colorbar(label='Sound Density Level [dB]')
#     plt.xlabel('X [m]')
#     plt.ylabel('Y [m]')
#     plt.title('Figure 10: Sound Density level at Z = z_rec and t = recording_time')
#     plt.show()
    
#     #Figure 11: 2D image of the SPL in the room
#     spl_2d = spl[:,:,dep_up_r] #The 3D w_new array is slised at the the desired z level
#     plt.figure(11)
#     plt.imshow(spl_2d, origin='lower', extent=[x[0], x[-1], y[0], y[-1]], aspect='equal') #plot with the extent being the room dimension x and y 
#     plt.colorbar(label='Sound Pressure Level [dB]')
#     plt.xlabel('X [m]')
#     plt.ylabel('Y [m]')
#     plt.title('Figure 11: Sound Pressure level at Z = z_rec and t = recording_time')
#     plt.show()
    
#     #Figure 12: Energy density at t=recording_time over the space x.
#     plt.figure(12)
#     plt.title("Figure 12: Energy density over the x axis at t=recording_time")
#     plt.plot(x,w_rec_x_end)
#     plt.ylabel('$\mathrm{Energy \ Density \ [kg/ms^2]}$')
#     plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$')
    
#     #Figure 13: Energy density at t=sourceoff_step over the space x.
#     plt.figure(13)
#     plt.title("Figure 13: Energy density over the x axis at t=0")
#     plt.plot(x,w_rec_x_t0)
#     plt.ylabel('$\mathrm{Energy \ Density \ [kg/ms^2]}$')
#     plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$') 
    
#     #Figure 14: SPL at t=sourceoff_step over the space x. reverb sound only
#     plt.figure(14)
#     plt.title("Figure 14: SPL REVERB over the x axis at t=0")
#     plt.plot(x,spl_rec_x_t0)
#     plt.ylabel('$\mathrm{SPL \ [dB]}$')
#     plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$') 
       
#     #Figure 15: Energy density at t=1*mean_free over the space x.
#     plt.figure(15)
#     plt.title("Figure 15: Energy density over the x axis at t=1*mean_free_time")
#     plt.plot(x,w_rec_x_1l)
#     plt.ylabel('$\mathrm{Energy \ Density \ [kg/ms^2]}$')
#     plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$') 
    
#     #Figure 16: SPL at  t=1*mean_free over the space x. reverb sound only
#     plt.figure(16)
#     plt.title("Figure 16: SPL REVERB over the x axis at t=1*mean_free")
#     plt.plot(x,spl_rec_x_1l)
#     plt.ylabel('$\mathrm{SPL \ [dB]}$')
#     plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$')    
    
#     #Figure 17: Energy density at t=2*mean_free over the space x.
#     plt.figure(17)
#     plt.title("Figure 17: Energy density over the x axis at t=2*mean_free_time")
#     plt.plot(x,w_rec_x_2l)
#     plt.ylabel('$\mathrm{Energy \ Density \ [kg/ms^2]}$')
#     plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$') 
    
#     #Figure 18: SPL at  t=2*mean_free over the space x. reverb sound only
#     plt.figure(18)
#     plt.title("Figure 18: SPL REVERB over the x axis at t=2*mean_free")
#     plt.plot(x,spl_rec_x_2l)
#     plt.ylabel('$\mathrm{SPL \ [dB]}$')
#     plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$')      
    
#     #Figure 19: Energy density at t=3*mean_free over the space x.
#     plt.figure(19)
#     plt.title("Figure 19: Energy density over the x axis at t=3*mean_free_time")
#     plt.plot(x,w_rec_x_3l)
#     plt.ylabel('$\mathrm{Energy \ Density \ [kg/ms^2]}$')
#     plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$') 
    
#     #Figure 20: SPL at  t=3*mean_free over the space x. reverb sound only
#     plt.figure(20)
#     plt.title("Figure 20: SPL REVERB over the x axis at t=3*mean_free")
#     plt.plot(x,spl_rec_x_3l)
#     plt.ylabel('$\mathrm{SPL \ [dB]}$')
#     plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$') 
    
#     #Figure 21: Energy density at t=5*mean_free over the space x.
#     plt.figure(21)
#     plt.title("Figure 21: Energy density over the x axis at t=5*mean_free_time")
#     plt.plot(x,w_rec_x_5l)
#     plt.ylabel('$\mathrm{Energy \ Density \ [kg/ms^2]}$')
#     plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$') 
    
#     #Figure 22: SPL at  t=5*mean_free over the space x. reverb sound only
#     plt.figure(22)
#     plt.title("Figure 22: SPL REVERB over the x axis at t=5*mean_free")
#     plt.plot(x,spl_rec_x_5l)
#     plt.ylabel('$\mathrm{SPL \ [dB]}$')
#     plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$') 
    
# if tcalc == "stationarysource":

#     #Figure 3: Decay of SPL in the recording_time at the receiver
#     plt.figure(3)
#     plt.plot(t,spl_r) #plot sound pressure level with Pref = (2e-5)**5
#     plt.title("Figure 3: SPL over time at the receiver")
#     plt.xlabel("t [s]")
#     plt.ylabel("SPL [dB]")
#     plt.xlim()
#     plt.ylim()
#     plt.xticks(np.arange(0, recording_time +0.1, 0.5))
#     #plt.yticks(np.arange(0, 120, 20))

#     #Figure 4: Decay of SPL in the recording_time normalised to maximum 0dB
#     plt.figure(4)
#     plt.title("Figure 4: Normalised SPL over time at the receiver")
#     plt.plot(t,spl_r_norm)
#     plt.xlabel("t [s]")
#     plt.ylabel("SPL [dB]")
#     plt.xlim()
#     plt.ylim()
#     plt.xticks(np.arange(0, recording_time +0.1, 0.1))
#     plt.yticks(np.arange(0, -60, -10))

#     #Figure 5: Energy density over time at the receiver
#     plt.figure(5)
#     plt.title("Figure 5: Energy density over time at the receiver")
#     plt.plot(t,w_rec)
#     plt.ylabel('$\mathrm{Energy \ Density \ [kg/ms^2]}$')
#     plt.xlabel("t [s]")

    
#     #Figure 8: Energy density at t=recording_time over the space x.
#     plt.figure(8)
#     plt.title("Figure 8: Energy density over the x axis at t=recording_time")
#     plt.plot(x,w_rec_x_end)
#     plt.ylabel('$\mathrm{Energy \ Density \ [kg/ms^2]}$')
#     plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$')
    
#%%
###############################################################################
#SAVING
###############################################################################
# Save all variables to a file
def save(filename):
    """
    Saving of variables
    
    Parameters
    ----------
        filename : str
            Name of the file to save the results
        variables : dict
            Compilation of all the variables of the overal simulation
    """
    with open(filename, 'wb') as f:
        # Filter out modules, functions, and other unsupported types
        filtered_variables = {}
        for k, v in globals().items():
            try:
                # Check if the object can be pickled
                pickle.dumps(v)
                # Exclude some types explicitly known to cause issues
                if not k.startswith('__') and not isinstance(v, (types.ModuleType, types.FunctionType, types.BuiltinFunctionType, types.LambdaType, types.MethodType, types.MappingProxyType)):
                    filtered_variables[k] = v
            except Exception as e:
                print(f"Could not pickle {k}: {str(e)}")

        pickle.dump(filtered_variables, f)

# To save all current variables
save('resultsFDM.pkl')
