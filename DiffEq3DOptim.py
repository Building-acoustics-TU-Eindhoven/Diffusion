# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 10:39:42 2023

@author: 20225533
"""
#Code developed by Ilaria Fichera for the analysis of the FDM method Du Fort & Frankel solving the 3D diffusion equation with one intermittent omnidirectional sound source
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simps
from scipy import linalg
import sys
#uncomment this if you need drawnow
from drawnow import drawnow
from math import ceil
from math import log
from FunctionRT import *
from FunctionEDT import *
from FunctionClarity import *
from FunctionDefinition import *
from FunctionCentreTime import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import time as time
from scipy import stats
from scipy.interpolate import griddata
from matplotlib.animation import FuncAnimation

st = time.time() #start time of calculation


def calculate_energy_density_diffusion(D):

    #%%
    ###############################################################################
    #INPUT VARIABLES
    ###############################################################################
    
    #General settings
    c0= 343 #adiabatic speed of sound [m.s^-1]
    m_atm = 0 #air absorption coefficient [1/m] from Billon 2008 paper and Navarro paper 2012
    
    #Room dimensions
    length = 3.0 #point x finish at the length of the room in the x direction [m] %Length
    width = 3.0 #point y finish at the length of the room in the y direction [m] %Width
    height = 3.0 #point z finish at the length of the room in the x direction [m] %Height
    
    # Source position
    x_source = 0.5  #position of the source in the x direction [m]
    y_source = 0.5  #position of the source in the y direction [m]
    z_source = 1.0  #position of the source in the z direction [m]
    
    # Receiver position
    x_rec = 2.0 #position of the receiver in the x direction [m]
    y_rec = 0.5 #position of the receiver in the y direction [m]
    z_rec = 1.0 #position of the receiver in the z direction [m]
    
    #Spatial discretization
    dx = 0.5 #distance between grid points x direction [m] #See Documentation for more insight about dt and dx
    dy = dx #distance between grid points y direction [m]
    dz = dx #distance between grid points z direction [m]
    
    #Time discretization
    dt = 1/8000 #distance between grid points on the time discretization [s] #See Documentation for more insight about dt and dx
    
    #Absorption term and Absorption coefficients
    th = 3 #int(input("Enter type Absortion conditions (option 1,2,3):")) 
    # options Sabine (th=1), Eyring (th=2) and modified by Xiang (th=3)
    alpha_1 = 1/6 #Absorption coefficient for Surface1 - Floor
    alpha_2 = 1/6 #Absorption coefficient for Surface2 - Ceiling
    alpha_3 = 1/6 #Absorption coefficient for Surface3 - Wall Front
    alpha_4 = 1/6 #Absorption coefficient for Surface4 - Wall Back
    alpha_5 = 1/6 #Absorption coefficient for Surface5 - Wall Left
    alpha_6 = 1/6 #Absorption coefficient for Surface6 - Wall Right
    
    #Type of Calculation
    #Choose "decay" if the objective is to calculate the energy decay of the room with all its energetic parameters; 
    #Choose "stationarysource" if the aim is to understand the behaviour of a room subject to a stationary source
    tcalc = "decay"
    
    #Set initial condition - Source Info (interrupted method)
    Ws = 0.01 #Source point power [Watts] interrupted after "sourceon_time" seconds; 10^-2 W => correspondent to 100dB
    sourceon_time =  0.50 #time that the source is ON before interrupting [s]
    recording_time = 2.00 #total time recorded for the calculation [s]
    
    #%%
    ###############################################################################
    #CALCULATION SECTION
    ###############################################################################
    
    #Fixed inputs
    pRef = 2 * (10**-5) #Reference pressure in Pa
    rho = 1.21 #air density [kg.m^-3] at 20°C
    
    #Frequency resolution & spatial parameters
    fsample = 1/dt #frequency spatial resolution (sampling period)
    
    #Time resolution
    t = np.arange(0, recording_time, dt) #mesh point in time
    recording_steps = ceil(recording_time/dt) #number of time steps to consider in the calculation
    
    #Room characteristics
    S1,S2 = length*width, length*width #xy planes
    S3,S4 = length*height, length*height #xz planes
    S5,S6 = width*height, width*height #yz planes
    S = length*width*2 + length*height*2 + width*height*2 #Total Surface Area[m2]
    V = length*width*height #Volume of the room [m^3]
    
    #Creating of meshgrid
    x = np.arange(0, length+dx, dx) #mesh points in space x direction
    y = np.arange(0, width+dy, dy) #mesh points in space y direction
    z = np.arange(0, height+dz, dz) #mesh points in space z direction
    Nx = len(x) #number of point in the x direction
    Ny = len(y) #number of point in the y direction
    Nz = len(z) #number of point in the z direction
    
    yy, xx , zz = np.meshgrid(y,x,z) #Return coordinate matrices from coordinate vectors; create the 3D grid
    
    #uncoment this when using drawnow
    #Figure 1: Visualization of 3D meshgrid
    fig = plt.figure(1)
    ax = plt.axes(projection ="3d")
    ax.scatter(xx, yy, zz, c=zz, cmap='Greens')
    plt.title("Figure 1: Visualization of 3D meshgrid")
    
    #Absorption term for boundary conditions 
    def abs_term(th,alpha):
        if th == 1:
            Absx = (c0*alpha)/4 #Sabine
        elif th == 2:
            Absx = (c0*(-log(1-alpha)))/4 #Eyring
        elif th == 3:
            Absx = (c0*alpha)/(2*(2-alpha)) #Modified by Xiang
        return Absx
    
    Abs_1 = abs_term(th,alpha_1) #absorption term for S1
    Abs_2 = abs_term(th,alpha_2) #absorption term for S2
    Abs_3 = abs_term(th,alpha_3) #absorption term for S3
    Abs_4 = abs_term(th,alpha_4) #absorption term for S4
    Abs_5 = abs_term(th,alpha_5) #absorption term for S5
    Abs_6 = abs_term(th,alpha_6) #absorption term for S6
    
    #Absorption parameters for room
    alpha_average = (alpha_1*S1 + alpha_2*S2 + alpha_3*S3 + alpha_4*S4 + alpha_5*S5 + alpha_6*S6)/S #average absorption
    Eq_A = alpha_1*S1 + alpha_2*S2 + alpha_3*S3 + alpha_4*S4 + alpha_5*S5 + alpha_6*S6 #equivalent absorption area of the room
    
    #Diffusion parameters
    lambda_path = (4*V)/S #mean free path for 3D
    lambda_time= lambda_path/c0 #mean free time for 3D
    lambda_time_step = int(lambda_time/dt)
    Dx = D
    Dy = D
    Dz = D
    
    #Mesh numbers
    beta_zero_x = (2*Dx*dt)/(dx**2) #mesh number in x direction
    beta_zero_y = (2*Dy*dt)/(dy**2) #mesh number in x direction
    beta_zero_z = (2*Dz*dt)/(dz**2) #mesh number in x direction
    beta_zero = beta_zero_x + beta_zero_y + beta_zero_z #beta_zero is the condition for all the directions deltax, deltay and deltaz.
     
    #Condition for the model to be unconditionally stable
    beta_zero_condition = ((beta_zero**2) - 1)/(1+(beta_zero**2)+(2*beta_zero)) #from Navarro 2012 paper
    if beta_zero_condition >1:
        print("aa! error! Check beta condition")
    
    #Initial condition - Source Info (interrupted method)
    Vs=dx*dy*dz  #Volume of the source
    w1=Ws/Vs #w1 = round(Ws/Vs,4) #power density of the source [Watts/(m^3))]
    sourceon_steps = ceil(sourceon_time/dt) #time steps at which the source is calculated/considered in the calculation
    s1 = np.multiply(w1,np.ones(sourceon_steps)) #energy density of source number 1 at each time step position
    source1 = np.append(s1, np.zeros(recording_steps-sourceon_steps)) #This would be equal to s1 if and only if recoding_steps = sourceon_steps
    
    ###############################################################################
    #SOURCE INTERPOLATION
    ###############################################################################
    #Finding index in meshgrid of the source position
    coord_source = [x_source,y_source,z_source] #coordinates of the receiver position in an list
    
    # Calculate the fractional indices
    row_lower = int(np.floor(x_source / dx))
    row_upper = row_lower + 1
    col_lower = int(np.floor(y_source / dy))
    col_upper = col_lower + 1
    depth_lower = int(np.floor(z_source / dz))
    depth_upper = depth_lower + 1
    
    # Calculate the interpolation weights
    weight_row_upper = (x_source / dx) - row_lower
    weight_row_lower = 1 - weight_row_upper
    weight_col_upper = (y_source / dy) - col_lower
    weight_col_lower = 1 - weight_col_upper
    weight_depth_upper = (z_source / dz) - depth_lower
    weight_depth_lower = 1 - weight_depth_upper
    
    s = np.zeros((Nx,Ny,Nz)) #matrix of zeros for source
    
    # Perform linear interpolation
    s[row_lower, col_lower, depth_lower] += source1[1] * weight_row_lower * weight_col_lower * weight_depth_lower
    s[row_lower, col_lower, depth_upper] += source1[1] * weight_row_lower * weight_col_lower * weight_depth_upper
    s[row_lower, col_upper, depth_lower] += source1[1] * weight_row_lower * weight_col_upper * weight_depth_lower
    s[row_lower, col_upper, depth_upper] += source1[1] * weight_row_lower * weight_col_upper * weight_depth_upper
    s[row_upper, col_lower, depth_lower] += source1[1] * weight_row_upper * weight_col_lower * weight_depth_lower
    s[row_upper, col_lower, depth_upper] += source1[1] * weight_row_upper * weight_col_lower * weight_depth_upper
    s[row_upper, col_upper, depth_lower] += source1[1] * weight_row_upper * weight_col_upper * weight_depth_lower
    s[row_upper, col_upper, depth_upper] += source1[1] * weight_row_upper * weight_col_upper * weight_depth_upper
    
    ###############################################################################
    #RECEIVER INTERPOLATION
    ###############################################################################
    #Finding index in meshgrid of the receiver position
    coord_receiver = [x_rec,y_rec,z_rec] #coordinates of the receiver position in an list
    
    #Calculate the fractional indices for receiver
    row_lr = int(np.floor(x_rec / dx))
    row_ur = row_lr + 1
    col_lr = int(np.floor(y_rec / dy))
    col_ur = col_lr + 1
    depth_lr = int(np.floor(z_rec / dz))
    depth_ur = depth_lr + 1
       
    #Calculate the interpolation weights for receiver
    weight_row_ur = (x_rec / dx) - row_lr #weight x upper
    weight_row_lr = 1 - weight_row_ur #weight x lower
    weight_col_ur = (y_rec / dy) - col_lr #weight y upper
    weight_col_lr = 1 - weight_col_ur #weight y lower
    weight_depth_ur = (z_rec / dz) - depth_lr #weight z upper
    weight_depth_lr = 1 - weight_depth_ur #weight z lower
    
    
    #distance between source and receiver
    dist_sr = math.sqrt((abs(x_rec - x_source))**2 + (abs(y_rec - y_source))**2 + (abs(z_rec - z_source))**2) #distance between source and receiver
    
    #distance between source and each mesh point in the x direction 
    dist_x = np.sqrt((((xx[:, col_lr, depth_lr]*(weight_row_lr * weight_col_lr * weight_depth_lr))+\
        (xx[:, col_lr, depth_ur]*(weight_row_lr * weight_col_lr * weight_depth_ur))+\
            (xx[:, col_ur, depth_lr]*(weight_row_lr * weight_col_ur * weight_depth_lr))+\
                (xx[:, col_ur, depth_ur]*(weight_row_lr * weight_col_ur * weight_depth_ur))+\
                    (xx[:, col_lr, depth_lr]*(weight_row_ur * weight_col_lr * weight_depth_lr))+\
                        (xx[:, col_lr, depth_ur]*(weight_row_ur * weight_col_lr * weight_depth_ur))+\
                            (xx[:, col_ur, depth_lr]*(weight_row_ur * weight_col_ur * weight_depth_lr))+\
                                (xx[:, col_ur, depth_ur]*(weight_row_ur * weight_col_ur * weight_depth_ur))) - x_source)**2 +\
                     (((yy[row_lr, col_lr, depth_lr]*(weight_row_lr * weight_col_lr * weight_depth_lr))+\
                         (yy[row_lr, col_lr, depth_ur]*(weight_row_lr * weight_col_lr * weight_depth_ur))+\
                             (yy[row_lr, col_ur, depth_lr]*(weight_row_lr * weight_col_ur * weight_depth_lr))+\
                                 (yy[row_lr, col_ur, depth_ur]*(weight_row_lr * weight_col_ur * weight_depth_ur))+\
                                     (yy[row_ur, col_lr, depth_lr]*(weight_row_ur * weight_col_lr * weight_depth_lr))+\
                                         (yy[row_ur, col_lr, depth_ur]*(weight_row_ur * weight_col_lr * weight_depth_ur))+\
                                             (yy[row_ur, col_ur, depth_lr]*(weight_row_ur * weight_col_ur * weight_depth_lr))+\
                                                 (yy[row_ur, col_ur, depth_ur]*(weight_row_ur * weight_col_ur * weight_depth_ur))) - y_source)**2 +\
                         (((zz[row_lr, col_lr, depth_lr]*(weight_row_lr * weight_col_lr * weight_depth_lr))+\
                             (zz[row_lr, col_lr, depth_ur]*(weight_row_lr * weight_col_lr * weight_depth_ur))+\
                                 (zz[row_lr, col_ur, depth_lr]*(weight_row_lr * weight_col_ur * weight_depth_lr))+\
                                     (zz[row_lr, col_ur, depth_ur]*(weight_row_lr * weight_col_ur * weight_depth_ur))+\
                                         (zz[row_ur, col_lr, depth_lr]*(weight_row_ur * weight_col_lr * weight_depth_lr))+\
                                             (zz[row_ur, col_lr, depth_ur]*(weight_row_ur * weight_col_lr * weight_depth_ur))+\
                                                 (zz[row_ur, col_ur, depth_lr]*(weight_row_ur * weight_col_ur * weight_depth_lr))+\
                                                     (zz[row_ur, col_ur, depth_ur]*(weight_row_ur * weight_col_ur * weight_depth_ur))) - z_source)**2)
    
    
    #distance between source and each mesh point in the y direction 
    dist_y = np.sqrt((((xx[row_lr, col_lr, depth_lr]*(weight_row_lr * weight_col_lr * weight_depth_lr))+\
        (xx[row_lr, col_lr, depth_ur]*(weight_row_lr * weight_col_lr * weight_depth_ur))+\
            (xx[row_lr, col_ur, depth_lr]*(weight_row_lr * weight_col_ur * weight_depth_lr))+\
                (xx[row_lr, col_ur, depth_ur]*(weight_row_lr * weight_col_ur * weight_depth_ur))+\
                    (xx[row_ur, col_lr, depth_lr]*(weight_row_ur * weight_col_lr * weight_depth_lr))+\
                        (xx[row_ur, col_lr, depth_ur]*(weight_row_ur * weight_col_lr * weight_depth_ur))+\
                            (xx[row_ur, col_ur, depth_lr]*(weight_row_ur * weight_col_ur * weight_depth_lr))+\
                                (xx[row_ur, col_ur, depth_ur]*(weight_row_ur * weight_col_ur * weight_depth_ur))) - x_source)**2 +\
                     (((yy[row_lr, :, depth_lr]*(weight_row_lr * weight_col_lr * weight_depth_lr))+\
                         (yy[row_lr, :, depth_ur]*(weight_row_lr * weight_col_lr * weight_depth_ur))+\
                             (yy[row_lr, :, depth_lr]*(weight_row_lr * weight_col_ur * weight_depth_lr))+\
                                 (yy[row_lr, :, depth_ur]*(weight_row_lr * weight_col_ur * weight_depth_ur))+\
                                     (yy[row_ur, :, depth_lr]*(weight_row_ur * weight_col_lr * weight_depth_lr))+\
                                         (yy[row_ur, :, depth_ur]*(weight_row_ur * weight_col_lr * weight_depth_ur))+\
                                             (yy[row_ur, :, depth_lr]*(weight_row_ur * weight_col_ur * weight_depth_lr))+\
                                                 (yy[row_ur, :, depth_ur]*(weight_row_ur * weight_col_ur * weight_depth_ur))) - y_source)**2 +\
                         (((zz[row_lr, col_lr, depth_lr]*(weight_row_lr * weight_col_lr * weight_depth_lr))+\
                             (zz[row_lr, col_lr, depth_ur]*(weight_row_lr * weight_col_lr * weight_depth_ur))+\
                                 (zz[row_lr, col_ur, depth_lr]*(weight_row_lr * weight_col_ur * weight_depth_lr))+\
                                     (zz[row_lr, col_ur, depth_ur]*(weight_row_lr * weight_col_ur * weight_depth_ur))+\
                                         (zz[row_ur, col_lr, depth_lr]*(weight_row_ur * weight_col_lr * weight_depth_lr))+\
                                             (zz[row_ur, col_lr, depth_ur]*(weight_row_ur * weight_col_lr * weight_depth_ur))+\
                                                 (zz[row_ur, col_ur, depth_lr]*(weight_row_ur * weight_col_ur * weight_depth_lr))+\
                                                     (zz[row_ur, col_ur, depth_ur]*(weight_row_ur * weight_col_ur * weight_depth_ur))) - z_source)**2)

    #%%
    ###############################################################################
    #MAIN CALCULATION - COMPUTING ENERGY DENSITY
    ############################################################################### 

    w_new = np.zeros((Nx,Ny,Nz)) #unknown w at new time level (n+1)
    w = w_new #w at n level
    w_old = w #w_old at n-1 level
    
    w_rec = np.arange(0,recording_time,dt) #energy density at the receiver
    
    #Computing w;
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
        w_new[0,:,:] = np.divide((4*w_new[1,:,:] - w_new[2,:,:]),(3+((2*Abs_5*dx)/Dx))) #boundary condition at x=0, any y, any z
        w_new[-1,:,:] = np.divide((4*w_new[-2,:,:] - w_new[-3,:,:]),(3+((2*Abs_6*dx)/Dx))) #boundary condition at lx=length, any y, any z
    
    
        w_new[:,0,:] = np.divide((4*w_new[:,1,:] - w_new[:,2,:]),(3+((2*Abs_3*dx)/Dy))) #boundary condition at y=0, any x, any z
        w_new[:,-1,:] = np.divide((4*w_new[:,-2,:] - w_new[:,-3,:]),(3+((2*Abs_4*dx)/Dy))) #boundary condition at at ly=width, any x, any z
     
        
        w_new[:,:,0] = np.divide((4*w_new[:,:,1] - w_new[:,:,2]),(3+((2*Abs_1*dx)/Dz))) #boundary condition at z=0, any x, any y
        w_new[:,:,-1] = np.divide((4*w_new[:,:,-2] - w_new[:,:,-3]),(3+((2*Abs_2*dx)/Dz))) #boundary condition at at lz=height, any x, any y
        
        sdl = 10*np.log10(abs(w_new),where=abs(w_new)>0) #sound density level
        spl = 10*np.log10(((abs(w_new))*rho*(c0**2))/(pRef**2)) #,where=press_r>0, sound pressure level in the 3D space
        
        #Update w before next step
        w_old = w #The w at n step becomes the w at n-1 step
        w = w_new #The w at n+1 step becomes the w at n step
        
        #w_rec is the energy density at the specific receiver
        w_rec[steps] = ((w_new[row_lr, col_lr, depth_lr]*(weight_row_lr * weight_col_lr * weight_depth_lr))+\
            (w_new[row_lr, col_lr, depth_ur]*(weight_row_lr * weight_col_lr * weight_depth_ur))+\
                (w_new[row_lr, col_ur, depth_lr]*(weight_row_lr * weight_col_ur * weight_depth_lr))+\
                    (w_new[row_lr, col_ur, depth_ur]*(weight_row_lr * weight_col_ur * weight_depth_ur))+\
                        (w_new[row_ur, col_lr, depth_lr]*(weight_row_ur * weight_col_lr * weight_depth_lr))+\
                            (w_new[row_ur, col_lr, depth_ur]*(weight_row_ur * weight_col_lr * weight_depth_ur))+\
                                (w_new[row_ur, col_ur, depth_lr]*(weight_row_ur * weight_col_ur * weight_depth_lr))+\
                                    (w_new[row_ur, col_ur, depth_ur]*(weight_row_ur * weight_col_ur * weight_depth_ur)))
    
        if steps == sourceon_steps:
            #print("Steps for source:",steps)
            w_t0 = w_new
    
        if steps == 5*lambda_time_step + sourceon_steps:
            w_5l = w_new
        
        #Updating the source term
        if tcalc == "decay":
            s[row_lower, col_lower, depth_lower] = source1[steps] * weight_row_lower * weight_col_lower * weight_depth_lower
            s[row_lower, col_lower, depth_upper] = source1[steps] * weight_row_lower * weight_col_lower * weight_depth_upper
            s[row_lower, col_upper, depth_lower] = source1[steps] * weight_row_lower * weight_col_upper * weight_depth_lower
            s[row_lower, col_upper, depth_upper] = source1[steps] * weight_row_lower * weight_col_upper * weight_depth_upper
            s[row_upper, col_lower, depth_lower] = source1[steps] * weight_row_upper * weight_col_lower * weight_depth_lower
            s[row_upper, col_lower, depth_upper] = source1[steps] * weight_row_upper * weight_col_lower * weight_depth_upper
            s[row_upper, col_upper, depth_lower] = source1[steps] * weight_row_upper * weight_col_upper * weight_depth_lower
            s[row_upper, col_upper, depth_upper] = source1[steps] * weight_row_upper * weight_col_upper * weight_depth_upper
        if tcalc == "stationarysource":
            s[row_lower, col_lower, depth_lower] = source1[0] * weight_row_lower * weight_col_lower * weight_depth_lower
            s[row_lower, col_lower, depth_upper] = source1[0] * weight_row_lower * weight_col_lower * weight_depth_upper
            s[row_lower, col_upper, depth_lower] = source1[0] * weight_row_lower * weight_col_upper * weight_depth_lower
            s[row_lower, col_upper, depth_upper] = source1[0] * weight_row_lower * weight_col_upper * weight_depth_upper
            s[row_upper, col_lower, depth_lower] = source1[0] * weight_row_upper * weight_col_lower * weight_depth_lower
            s[row_upper, col_lower, depth_upper] = source1[0] * weight_row_upper * weight_col_lower * weight_depth_upper
            s[row_upper, col_upper, depth_lower] = source1[0] * weight_row_upper * weight_col_upper * weight_depth_lower
            s[row_upper, col_upper, depth_upper] = source1[0] * weight_row_upper * weight_col_upper * weight_depth_upper
            
        
        #print(time_steps)
    
    plt.show()
    print(D)
    
    w_rec_x = ((w_new[:, col_lr, depth_lr]*(weight_row_lr * weight_col_lr * weight_depth_lr))+\
            (w_new[:, col_lr, depth_ur]*(weight_row_lr * weight_col_lr * weight_depth_ur))+\
                (w_new[:, col_ur, depth_lr]*(weight_row_lr * weight_col_ur * weight_depth_lr))+\
                    (w_new[:, col_ur, depth_ur]*(weight_row_lr * weight_col_ur * weight_depth_ur))+\
                        (w_new[:, col_lr, depth_lr]*(weight_row_ur * weight_col_lr * weight_depth_lr))+\
                            (w_new[:, col_lr, depth_ur]*(weight_row_ur * weight_col_lr * weight_depth_ur))+\
                                (w_new[:, col_ur, depth_lr]*(weight_row_ur * weight_col_ur * weight_depth_lr))+\
                                    (w_new[:, col_ur, depth_ur]*(weight_row_ur * weight_col_ur * weight_depth_ur)))
            
    w_rec_x_5l = ((w_5l[:, col_lr, depth_lr]*(weight_row_lr * weight_col_lr * weight_depth_lr))+\
                (w_5l[:, col_lr, depth_ur]*(weight_row_lr * weight_col_lr * weight_depth_ur))+\
                    (w_5l[:, col_ur, depth_lr]*(weight_row_lr * weight_col_ur * weight_depth_lr))+\
                        (w_5l[:, col_ur, depth_ur]*(weight_row_lr * weight_col_ur * weight_depth_ur))+\
                            (w_5l[:, col_lr, depth_lr]*(weight_row_ur * weight_col_lr * weight_depth_lr))+\
                                (w_5l[:, col_lr, depth_ur]*(weight_row_ur * weight_col_lr * weight_depth_ur))+\
                                    (w_5l[:, col_ur, depth_lr]*(weight_row_ur * weight_col_ur * weight_depth_lr))+\
                                        (w_5l[:, col_ur, depth_ur]*(weight_row_ur * weight_col_ur * weight_depth_ur)))
        
    w_rec_y = ((w_new[row_lr, :, depth_lr]*(weight_row_lr * weight_col_lr * weight_depth_lr))+\
            (w_new[row_lr, :, depth_ur]*(weight_row_lr * weight_col_lr * weight_depth_ur))+\
                (w_new[row_lr, :, depth_lr]*(weight_row_lr * weight_col_ur * weight_depth_lr))+\
                    (w_new[row_lr, :, depth_ur]*(weight_row_lr * weight_col_ur * weight_depth_ur))+\
                        (w_new[row_ur, :, depth_lr]*(weight_row_ur * weight_col_lr * weight_depth_lr))+\
                            (w_new[row_ur, :, depth_ur]*(weight_row_ur * weight_col_lr * weight_depth_ur))+\
                                (w_new[row_ur, :, depth_lr]*(weight_row_ur * weight_col_ur * weight_depth_lr))+\
                                    (w_new[row_ur, :, depth_ur]*(weight_row_ur * weight_col_ur * weight_depth_ur)))
    return w_rec_x

#Dx = 228

#w_rec_x = calculate_energy_density_diffusion(Dx)