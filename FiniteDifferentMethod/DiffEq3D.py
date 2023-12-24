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
#from drawnow import drawnow
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

#%%
###############################################################################
#INPUT VARIABLES
###############################################################################

#General settings
c0= 343 #adiabatic speed of sound [m.s^-1]
m_atm = 0 #air absorption coefficient [1/m] from Billon 2008 paper and Navarro paper 2012

#Room dimensions
length = 21.0 #point x finish at the length of the room in the x direction [m] %Length
width = 3.0 #point y finish at the length of the room in the y direction [m] %Width
height = 3.0 #point z finish at the length of the room in the x direction [m] %Height

# Source position
x_source = 0.5  #position of the source in the x direction [m]
y_source = 0.5  #position of the source in the y direction [m]
z_source = 1.0  #position of the source in the z direction [m]

# Receiver position
x_rec = 20.0 #position of the receiver in the x direction [m]
y_rec = 0.5 #position of the receiver in the y direction [m]
z_rec = 1.0 #position of the receiver in the z direction [m]

#Spatial discretization
dx = 0.50 #distance between grid points x direction [m] #See Documentation for more insight about dt and dx
dy = dx #distance between grid points y direction [m]
dz = dx #distance between grid points z direction [m]

#Time discretization
dt = 1/2000 #distance between grid points on the time discretization [s] #See Documentation for more insight about dt and dx

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
sourceon_time =  0.5 #time that the source is ON before interrupting [s]
recording_time = 1.0 #total time recorded for the calculation [s]

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
mean_free_path = (4*V)/S #mean free path for 3D
mean_free_time= mean_free_path/c0 #mean free time for 3D
mean_free_time_step = int(mean_free_time/dt)
Dx = (mean_free_path*c0)/3 #diffusion coefficient for proportionate rooms x direction
Dy = (mean_free_path*c0)/3 #diffusion coefficient for proportionate rooms y direction
Dz = (mean_free_path*c0)/3 #diffusion coefficient for proportionate rooms z direction


#Longest dimension in the room
longest_dimension = math.sqrt(length**2+width**2)
longest_dimension_time = longest_dimension/c0
longest_dimension_step = int(longest_dimension_time/dt)

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
row_lr_s = int(np.floor(x_source / dx))
row_up_s = row_lr_s + 1
col_lr_s = int(np.floor(y_source / dy))
col_up_s = col_lr_s + 1
dep_lr_s = int(np.floor(z_source / dz))
dep_up_s = dep_lr_s + 1

# Calculate the interpolation weights
weight_row_up_s = (x_source / dx) - row_lr_s
weight_row_lr_s = 1 - weight_row_up_s
weight_col_up_s = (y_source / dy) - col_lr_s
weight_col_lr_s = 1 - weight_col_up_s
weight_dep_up_s = (z_source / dz) - dep_lr_s
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

###############################################################################
#RECEIVER INTERPOLATION
###############################################################################
#Finding index in meshgrid of the receiver position
coord_receiver = [x_rec,y_rec,z_rec] #coordinates of the receiver position in an list

#Calculate the fractional indices for receiver
row_lr_r = int(np.floor(x_rec / dx))
row_up_r = row_lr_r + 1
col_lr_r = int(np.floor(y_rec / dy))
col_up_r = col_lr_r + 1
dep_lr_r = int(np.floor(z_rec / dz))
dep_up_r = dep_lr_r + 1
   
#Calculate the interpolation weights for receiver
weight_row_up_r = (x_rec / dx) - row_lr_r #weight x upper
weight_row_lr_r = 1 - weight_row_up_r #weight x lower
weight_col_up_r = (y_rec / dy) - col_lr_r #weight y upper
weight_col_lr_r = 1 - weight_col_up_r #weight y lower
weight_dep_up_r = (z_rec / dz) - dep_lr_r #weight z upper
weight_dep_lr_r = 1 - weight_dep_up_r #weight z lower


#distance between source and receiver
dist_sr = math.sqrt((abs(x_rec - x_source))**2 + (abs(y_rec - y_source))**2 + (abs(z_rec - z_source))**2) #distance between source and receiver

#distance between source and each mesh point in the x direction 
dist_x = np.sqrt((((xx[:, col_lr_r, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_lr_r))+\
    (xx[:, col_lr_r, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
        (xx[:, col_up_r, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
            (xx[:, col_up_r, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                (xx[:, col_lr_r, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                    (xx[:, col_lr_r, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                        (xx[:, col_up_r, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                            (xx[:, col_up_r, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r))) - x_source)**2 +\
                 (((yy[row_lr_r, col_lr_r, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r* weight_dep_lr_r))+\
                     (yy[row_lr_r, col_lr_r, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
                         (yy[row_lr_r, col_up_r, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
                             (yy[row_lr_r, col_up_r, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                                 (yy[row_up_r, col_lr_r, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                                     (yy[row_up_r, col_lr_r, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                                         (yy[row_up_r, col_up_r, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                                             (yy[row_up_r, col_up_r, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r))) - y_source)**2 +\
                     (((zz[row_lr_r, col_lr_r, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_lr_r))+\
                         (zz[row_lr_r, col_lr_r, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
                             (zz[row_lr_r, col_up_r, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
                                 (zz[row_lr_r, col_up_r, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                                     (zz[row_up_r, col_lr_r, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                                         (zz[row_up_r, col_lr_r, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                                             (zz[row_up_r, col_up_r, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                                                 (zz[row_up_r, col_up_r, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r))) - z_source)**2)


#distance between source and each mesh point in the y direction  
dist_y = np.sqrt((((xx[row_lr_r, col_lr_r, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_lr_r))+\
    (xx[row_lr_r, col_lr_r, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
        (xx[row_lr_r, col_up_r, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
            (xx[row_lr_r, col_up_r, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                (xx[row_up_r, col_lr_r, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                    (xx[row_up_r, col_lr_r, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                        (xx[row_up_r, col_up_r, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                            (xx[row_up_r, col_up_r, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r))) - x_source)**2 +\
                 (((yy[row_lr_r, :, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r* weight_dep_lr_r))+\
                     (yy[row_lr_r, :, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
                         (yy[row_lr_r, :, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
                             (yy[row_lr_r, :, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                                 (yy[row_up_r, :, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                                     (yy[row_up_r, :, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                                         (yy[row_up_r, :, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                                             (yy[row_up_r, :, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r))) - y_source)**2 +\
                     (((zz[row_lr_r, col_lr_r, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_lr_r))+\
                         (zz[row_lr_r, col_lr_r, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
                             (zz[row_lr_r, col_up_r, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
                                 (zz[row_lr_r, col_up_r, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                                     (zz[row_up_r, col_lr_r, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                                         (zz[row_up_r, col_lr_r, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                                             (zz[row_up_r, col_up_r, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                                                 (zz[row_up_r, col_up_r, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r))) - z_source)**2)


#Function to draw figure
##def draw_fig():
##    for i in range(w_new.shape[2]):
##        print(i)
##        plt.figure(2)
##        plt.imshow(w_new[:, :, i], cmap='hot', vmin=w_new.min(), vmax=w_new.max())

############################################################################
#TRIAL 4D PLOT
############################################################################
#Initialize the figure and axes
#if tcalc == "decay":
#    plt.figure()
#    ax = fig.add_subplot(111, projection='3d')
#    plt.title("4D plot")
#    ax.view_init(elev=30, azim=45) #Set the initial view angle

#%%
###############################################################################
#MAIN CALCULATION - COMPUTING ENERGY DENSITY
############################################################################### 

w_new = np.zeros((Nx,Ny,Nz)) #unknown w at new time level (n+1)
w = w_new #w at n level
w_old = w #w_old at n-1 level

w_rec = np.arange(0,recording_time,dt) #energy density at the receiver
w_rec_all = np.zeros((1,len(x))) 

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
    
   #uncomment when activating the drawnow library   
   #Visualization of the energy density changes while the calculation is progressing
    ##if (steps % 100 == 0): #draw only on certain steps and not all the steps
    ##    print("A")
    ##    drawnow(draw_fig)
    
    #Update w before next step
    w_old = w #The w at n step becomes the w at n-1 step
    w = w_new #The w at n+1 step becomes the w at n step
    
    #w_rec is the energy density at the specific receiver
    w_rec[steps] = ((w_new[row_lr_r, col_lr_r, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_lr_r))+\
        (w_new[row_lr_r, col_lr_r, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
            (w_new[row_lr_r, col_up_r, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
                (w_new[row_lr_r, col_up_r, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                    (w_new[row_up_r, col_lr_r, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                        (w_new[row_up_r, col_lr_r, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                            (w_new[row_up_r, col_up_r, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                                (w_new[row_up_r, col_up_r, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r)))
            
    if steps == sourceon_steps:
        #print("Steps for source:",steps)
        w_t0 = w_new

    if steps == round(1*mean_free_time_step + sourceon_steps + (dist_sr/c0)):
        index_1l = steps
        w_1l = w_new
        
    if steps == round(2*mean_free_time_step + sourceon_steps + (dist_sr/c0)):
        index_2l = steps
        w_2l = w_new
        
    if steps == round(3*mean_free_time_step + sourceon_steps + (dist_sr/c0)):
        index_3l = steps
        w_3l = w_new

    if steps == round(5*mean_free_time_step + sourceon_steps + (dist_sr/c0)):
        index_5l = steps
        w_5l = w_new
    
    if steps == round(2*longest_dimension_step): #+ sourceon_steps + (dist_sr/c0)):
        index_2ld = steps
        w_2ld = w_new
        
    if steps == round(4*longest_dimension_step): #+ sourceon_steps + (dist_sr/c0)):
        index_4ld = steps
        w_4ld = w_new
    #4D Visualization????
    #Flatten the coordinates and w_new values for scatter plot
    ##coords = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])
    ##w_new_flat = w_new.ravel()
    
    #Normalize the w_new values to [0, 1] range
    ##norm = plt.Normalize(vmin=np.min(w_new_flat), vmax=np.max(w_new_flat))
    ##colors = cm.jet(norm(w_new_flat))  # Use 'jet' colormap for a range of red to yellow colors
    
    #Update the scatter plot with the new w_new values
    ##if steps == 0:
    ##    scatter = ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=w_new_flat, cmap='hot')
    ##else:
    ##    scatter.set_array(w_new_flat)
        
    #Set a suitable title and labels
    ##ax.set_title("Time Step: 'time_steps'")
    ##ax.set_xlabel('X')
    ##ax.set_ylabel('Y')
    ##ax.set_zlabel('Z')

    #Adjust the plot limits
    ##ax.set_xlim(0, length)
    ##ax.set_ylim(0, width)
    ##ax.set_zlim(0, height)

    #Add a colorbar and legend
    ##cbar = fig.colorbar(scatter, ax=ax,fraction=0.04, pad=0.1)
    ##cbar.set_label('Energy Density')
    ##cbar.ax.yaxis.set_ticks_position('left')

    #Pause to create an animated effect
    ##plt.pause(0.01)  # Adjust the pause duration as needed
    
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
    
    
    print(time_steps)

plt.show()

w_rec_x_end = ((w_new[:, col_lr_r, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_lr_r))+\
    (w_new[:, col_lr_r, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
        (w_new[:, col_up_r, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
            (w_new[:, col_up_r, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                (w_new[:, col_lr_r, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                    (w_new[:, col_lr_r, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                        (w_new[:, col_up_r, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                            (w_new[:, col_up_r, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r)))

spl_rec_x_end = 10*np.log10(rho*c0**2*w_rec_x_end/pRef**2)
    
w_rec_x_1l = ((w_1l[:, col_lr_r, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_lr_r))+\
     (w_1l[:, col_lr_r, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
         (w_1l[:, col_up_r, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
             (w_1l[:, col_up_r, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                 (w_1l[:, col_lr_r, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                     (w_1l[:, col_lr_r, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                         (w_1l[:, col_up_r, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                             (w_1l[:, col_up_r, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r)))

spl_rec_x_1l = 10*np.log10(rho*c0**2*w_rec_x_1l/pRef**2)

w_rec_x_2l = ((w_2l[:, col_lr_r, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_lr_r))+\
       (w_2l[:, col_lr_r, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
           (w_2l[:, col_up_r, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
               (w_2l[:, col_up_r, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                   (w_2l[:, col_lr_r, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                       (w_2l[:, col_lr_r, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                           (w_2l[:, col_up_r, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                               (w_2l[:, col_up_r, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r)))

spl_rec_x_2l = 10*np.log10(rho*c0**2*w_rec_x_2l/pRef**2)
    
w_rec_x_3l = ((w_3l[:, col_lr_r, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_lr_r))+\
       (w_3l[:, col_lr_r, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
           (w_3l[:, col_up_r, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
               (w_3l[:, col_up_r, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                   (w_3l[:, col_lr_r, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                       (w_3l[:, col_lr_r, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                           (w_3l[:, col_up_r, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                               (w_3l[:, col_up_r, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r)))

spl_rec_x_3l = 10*np.log10(rho*c0**2*w_rec_x_3l/pRef**2)

w_rec_x_5l = ((w_5l[:, col_lr_r, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_lr_r))+\
    (w_5l[:, col_lr_r, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
        (w_5l[:, col_up_r, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
            (w_5l[:, col_up_r, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                (w_5l[:, col_lr_r, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                    (w_5l[:, col_lr_r, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                        (w_5l[:, col_up_r, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                            (w_5l[:, col_up_r, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r)))    

spl_rec_x_5l = 10*np.log10(rho*c0**2*w_rec_x_5l/pRef**2)

w_rec_x_2ld = ((w_2ld[:, col_lr_r, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_lr_r))+\
    (w_2ld[:, col_lr_r, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
        (w_2ld[:, col_up_r, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
            (w_2ld[:, col_up_r, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                (w_2ld[:, col_lr_r, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                    (w_2ld[:, col_lr_r, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                        (w_2ld[:, col_up_r, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                            (w_2ld[:, col_up_r, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r)))    

spl_rec_x_2ld = 10*np.log10(rho*c0**2*w_rec_x_2ld/pRef**2)

w_rec_x_4ld = ((w_4ld[:, col_lr_r, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_lr_r))+\
    (w_4ld[:, col_lr_r, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
        (w_4ld[:, col_up_r, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
            (w_4ld[:, col_up_r, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                (w_4ld[:, col_lr_r, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                    (w_4ld[:, col_lr_r, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                        (w_4ld[:, col_up_r, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                            (w_4ld[:, col_up_r, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r)))    

spl_rec_x_4ld = 10*np.log10(rho*c0**2*w_rec_x_4ld/pRef**2)

    
w_rec_x_t0 = ((w_t0[:, col_lr_r, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_lr_r))+\
    (w_t0[:, col_lr_r, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
        (w_t0[:, col_up_r, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
            (w_t0[:, col_up_r, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                (w_t0[:, col_lr_r, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                    (w_t0[:, col_lr_r, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                        (w_t0[:, col_up_r, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                            (w_t0[:, col_up_r, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r)))        

spl_rec_x_t0 = 10*np.log10(rho*c0**2*w_rec_x_t0/pRef**2)
    
w_rec_y_end = ((w_new[row_lr_r, :, dep_lr_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_lr_r))+\
    (w_new[row_lr_r, :, dep_up_r]*(weight_row_lr_r * weight_col_lr_r * weight_dep_up_r))+\
        (w_new[row_lr_r, :, dep_lr_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_lr_r))+\
            (w_new[row_lr_r, :, dep_up_r]*(weight_row_lr_r * weight_col_up_r * weight_dep_up_r))+\
                (w_new[row_up_r, :, dep_lr_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_lr_r))+\
                    (w_new[row_up_r, :, dep_up_r]*(weight_row_up_r * weight_col_lr_r * weight_dep_up_r))+\
                        (w_new[row_up_r, :, dep_lr_r]*(weight_row_up_r * weight_col_up_r * weight_dep_lr_r))+\
                            (w_new[row_up_r, :, dep_up_r]*(weight_row_up_r * weight_col_up_r * weight_dep_up_r)))

#%%

###############################################################################
#RESULTS
###############################################################################


spl_stat_x_t0 = 10*np.log10(rho*c0*(((Ws)/(4*math.pi*(dist_x**2))) + ((abs(w_rec_x_t0)*c0)))/(pRef**2)) #with direct sound
spl_stat_x_5l = 10*np.log10(rho*c0**2*w_rec_x_5l/pRef**2)

spl_stat_x_rev = 10*np.log10(rho*c0*(((abs(w_rec_x_end)*c0)))/(pRef**2))

spl_stat_x = 10*np.log10(rho*c0*(((Ws)/(4*math.pi*(dist_x**2))) + ((abs(w_rec_x_end)*c0)))/(pRef**2))
spl_stat_y = 10*np.log10(rho*c0*(((Ws)/(4*math.pi*(dist_y**2))) + ((abs(w_rec_y_end)*c0)))/(pRef**2)) #It should be the spl stationary

press_r = ((abs(w_rec))*rho*(c0**2)) #pressure at the receiver
spl_r = 10*np.log10(((abs(w_rec))*rho*(c0**2))/(pRef**2)) #,where=press_r>0, sound pressure level at the receiver
spl_r_norm = 10*np.log10((((abs(w_rec))*rho*(c0**2))/(pRef**2)) / np.max(((abs(w_rec))*rho*(c0**2))/(pRef**2))) #normalised to maximum to 0dB
spl_r_tot = 10*np.log10(rho*c0*((Ws/(4*math.pi*dist_sr**2))*np.exp(-m_atm*dist_sr) + ((abs(w_rec))*c0)/(pRef**2))) #spl total (including direct field) at the receiver position????? but it will need to be calculated for a stationary source 100dB

#Find the energy decay part of the overal calculation
idx_w_rec = np.where(t == sourceon_time)[0][0] #index at which the t array is equal to the sourceon_time; I want the RT to calculate from when the source stops.
w_rec_off = w_rec[idx_w_rec:] #cutting the energy density array at the receiver from the idx_w_rec to the end

#Impulse response from the energy density
w_rec_off_deriv = w_rec_off #initialising an array of derivative equal to the w_rec_off -> this will be the impulse response after modifying it
w_rec_off_deriv = np.delete(w_rec_off_deriv, 0) #delete the first element of the array -> this means shifting the array one step before and therefore do a derivation
w_rec_off_deriv = np.append(w_rec_off_deriv,0) #add a zero in the last element of the array -> for derivation and to have the same length as previously
impulse = ((w_rec_off - w_rec_off_deriv))/dt#/(rho*c0**2) #This is the difference between the the energy density and the impulse response 

#Schroeder integration
#energy_r_rev = (w_rec_off)[::-1] #reverting the array
#The energy density is related to the pressure with the following relation: w = p^2
#energy_r_rev_cum = np.cumsum(energy_r_rev) #cumulative summation of all the item in the array
schroeder = w_rec_off #energy_r_rev_cum[::-1] #reverting the array again -> creating the schroder decay
sch_db = 10.0 * np.log10(schroeder / max(schroeder)) #level of the array: schroeder decay

if tcalc == "decay":
    t60 = t60_decay(t, sch_db, idx_w_rec) #called function for calculation of t60 [s]
    edt = edt_decay(t, sch_db, idx_w_rec) #called function for calculation of edt [s]
    c80 = clarity(t60, V, Eq_A, S, c0, dist_sr) #called function for calculation of c80 [dB]
    d50 = definition(t60, V, Eq_A, S, c0, dist_sr) #called function for calculation of d50 [%]
    ts = centretime(t60, Eq_A, S) #called function for calculation of ts [ms]

et = time.time() #end time
elapsed_time = et - st

#%%
###############################################################################
#FIGURES & POST-PROCESSING
###############################################################################

if tcalc == "decay":
    #Figure 5: Decay of SPL in the recording_time
    plt.figure(5)
    plt.plot(t, spl_r)  # plot sound pressure level with Pref = (2e-5)**5
    plt.title("Figure 5 :SPL over time at the receiver")
    plt.xlabel("t [s]")
    plt.ylabel("SPL [dB]")
    plt.xlim()
    plt.ylim()
    plt.xticks(np.arange(0, recording_time + 0.1, 0.5))
    plt.yticks(np.arange(0, 120, 20))

    #Figure 6: Decay of SPL in the recording_time normalised to maximum 0dB
    plt.figure(6)
    plt.plot(t,spl_r_norm)
    plt.title("Figure 6: Normalised SPL over time at the receiver")
    plt.xlabel("t [s]")
    plt.ylabel("SPL [dB]")
    plt.xlim()
    plt.ylim()
    plt.xticks(np.arange(0, recording_time +0.1, 0.1))
    plt.yticks(np.arange(0, -60, -10))
    
    #Figure 7: Energy density at the receiver over time
    plt.figure(7)
    plt.plot(t,w_rec)
    plt.title("Figure 7: Energy density over time at the receiver")
    plt.xlabel("t [s]")
    plt.ylabel("Energy density [kg m^-1 s^-2]")
    plt.xlim()
    plt.ylim()
    plt.xticks(np.arange(0, recording_time +0.1, 0.1))
    
    #Figure 8: Schroeder decay
    plt.figure(8)
    plt.plot(t[idx_w_rec:],sch_db)
    plt.title("Figure 8: Schroeder decay (Energy Decay Curve)")
    plt.xlabel("t [s]")
    plt.ylabel("Energy decay [dB]")
    plt.xlim()
    plt.ylim()
    #plt.xticks(np.arange(t, recording_time +0.1, 0.1))
    
    #Figure 9: 2D image of the energy density in the room
    w_new_2d = w_new[:,:,dep_up_r] #The 3D w_new array is slised at the the desired z level
    plt.figure(9)
    plt.imshow(w_new_2d, origin='lower', extent=[x[0], x[-1], y[0], y[-1]], aspect='equal') #plot with the extent being the room dimension x and y 
    plt.colorbar(label='Energy Density [kg/ms^2]')
    plt.xlabel('X [m]')
    plt.ylabel('Y [m]')
    plt.title('Figure 9: Energy Density at Z = z_rec and t = recording_time')
    plt.show()
    
    #Figure 10: 2D image of the SDL in the room
    sdl_2d = sdl[:,:,dep_up_r] #The 3D w_new array is slised at the the desired z level
    plt.figure(10)
    plt.imshow(sdl_2d, origin='lower', extent=[x[0], x[-1], y[0], y[-1]], aspect='equal') #plot with the extent being the room dimension x and y 
    plt.colorbar(label='Sound Density Level [dB]')
    plt.xlabel('X [m]')
    plt.ylabel('Y [m]')
    plt.title('Figure 10: Sound Density level at Z = z_rec and t = recording_time')
    plt.show()
    
    #Figure 11: 2D image of the SPL in the room
    spl_2d = spl[:,:,dep_up_r] #The 3D w_new array is slised at the the desired z level
    plt.figure(11)
    plt.imshow(spl_2d, origin='lower', extent=[x[0], x[-1], y[0], y[-1]], aspect='equal') #plot with the extent being the room dimension x and y 
    plt.colorbar(label='Sound Pressure Level [dB]')
    plt.xlabel('X [m]')
    plt.ylabel('Y [m]')
    plt.title('Figure 11: Sound Pressure level at Z = z_rec and t = recording_time')
    plt.show()
    
    #Figure 12: Energy density at t=recording_time over the space x.
    plt.figure(12)
    plt.title("Figure 12: Energy density over the x axis at t=recording_time")
    plt.plot(x,w_rec_x_end)
    plt.ylabel('$\mathrm{Energy \ Density \ [kg/ms^2]}$')
    plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$')
    
    #Figure 13: Energy density at t=sourceoff_step over the space x.
    plt.figure(13)
    plt.title("Figure 13: Energy density over the x axis at t=0")
    plt.plot(x,w_rec_x_t0)
    plt.ylabel('$\mathrm{Energy \ Density \ [kg/ms^2]}$')
    plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$') 
    
    #Figure 14: SPL at t=sourceoff_step over the space x. reverb sound only
    plt.figure(14)
    plt.title("Figure 14: SPL REVERB over the x axis at t=0")
    plt.plot(x,spl_rec_x_t0)
    plt.ylabel('$\mathrm{SPL \ [dB]}$')
    plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$') 
       
    #Figure 15: Energy density at t=1*mean_free over the space x.
    plt.figure(15)
    plt.title("Figure 15: Energy density over the x axis at t=1*mean_free_time")
    plt.plot(x,w_rec_x_1l)
    plt.ylabel('$\mathrm{Energy \ Density \ [kg/ms^2]}$')
    plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$') 
    
    #Figure 16: SPL at  t=1*mean_free over the space x. reverb sound only
    plt.figure(16)
    plt.title("Figure 16: SPL REVERB over the x axis at t=1*mean_free")
    plt.plot(x,spl_rec_x_1l)
    plt.ylabel('$\mathrm{SPL \ [dB]}$')
    plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$')    
    
    #Figure 17: Energy density at t=2*mean_free over the space x.
    plt.figure(17)
    plt.title("Figure 17: Energy density over the x axis at t=2*mean_free_time")
    plt.plot(x,w_rec_x_2l)
    plt.ylabel('$\mathrm{Energy \ Density \ [kg/ms^2]}$')
    plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$') 
    
    #Figure 18: SPL at  t=2*mean_free over the space x. reverb sound only
    plt.figure(18)
    plt.title("Figure 18: SPL REVERB over the x axis at t=2*mean_free")
    plt.plot(x,spl_rec_x_2l)
    plt.ylabel('$\mathrm{SPL \ [dB]}$')
    plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$')      
    
    #Figure 19: Energy density at t=3*mean_free over the space x.
    plt.figure(19)
    plt.title("Figure 19: Energy density over the x axis at t=3*mean_free_time")
    plt.plot(x,w_rec_x_3l)
    plt.ylabel('$\mathrm{Energy \ Density \ [kg/ms^2]}$')
    plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$') 
    
    #Figure 20: SPL at  t=3*mean_free over the space x. reverb sound only
    plt.figure(20)
    plt.title("Figure 20: SPL REVERB over the x axis at t=3*mean_free")
    plt.plot(x,spl_rec_x_3l)
    plt.ylabel('$\mathrm{SPL \ [dB]}$')
    plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$') 
    
    #Figure 21: Energy density at t=5*mean_free over the space x.
    plt.figure(21)
    plt.title("Figure 21: Energy density over the x axis at t=5*mean_free_time")
    plt.plot(x,w_rec_x_5l)
    plt.ylabel('$\mathrm{Energy \ Density \ [kg/ms^2]}$')
    plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$') 
    
    #Figure 22: SPL at  t=5*mean_free over the space x. reverb sound only
    plt.figure(22)
    plt.title("Figure 22: SPL REVERB over the x axis at t=5*mean_free")
    plt.plot(x,spl_rec_x_5l)
    plt.ylabel('$\mathrm{SPL \ [dB]}$')
    plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$') 
    
    #Figure 23: Energy density at t=2*ld over the space x.
    plt.figure(23)
    plt.title("Figure 23: Energy density over the x axis at t=2*longest_dimension_time")
    plt.plot(x,w_rec_x_2ld)
    plt.ylabel('$\mathrm{Energy \ Density \ [kg/ms^2]}$')
    plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$') 
    
    #Figure 24: SPL at  t=2*ld over the space x. reverb sound only
    plt.figure(24)
    plt.title("Figure 24: SPL REVERB over the x axis at t=2*longest_dimension_time")
    plt.plot(x,spl_rec_x_2ld)
    plt.ylabel('$\mathrm{SPL \ [dB]}$')
    plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$') 
    
    #Figure 25: Energy density at t=4*ld over the space x.
    plt.figure(25)
    plt.title("Figure 25: Energy density over the x axis at t=4*longest_dimension_time")
    plt.plot(x,w_rec_x_4ld)
    plt.ylabel('$\mathrm{Energy \ Density \ [kg/ms^2]}$')
    plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$') 
    
    #Figure 26: SPL at  t=4*ld over the space x. reverb sound only
    plt.figure(26)
    plt.title("Figure 26: SPL REVERB over the x axis at t=4*longest_dimension_time")
    plt.plot(x,spl_rec_x_4ld)
    plt.ylabel('$\mathrm{SPL \ [dB]}$')
    plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$') 
    
if tcalc == "stationarysource":

    #Figure 3: Decay of SPL in the recording_time at the receiver
    plt.figure(3)
    plt.plot(t,spl_r) #plot sound pressure level with Pref = (2e-5)**5
    plt.title("Figure 3: SPL over time at the receiver")
    plt.xlabel("t [s]")
    plt.ylabel("SPL [dB]")
    plt.xlim()
    plt.ylim()
    plt.xticks(np.arange(0, recording_time +0.1, 0.5))
    #plt.yticks(np.arange(0, 120, 20))

    #Figure 4: Decay of SPL in the recording_time normalised to maximum 0dB
    plt.figure(4)
    plt.title("Figure 4: Normalised SPL over time at the receiver")
    plt.plot(t,spl_r_norm)
    plt.xlabel("t [s]")
    plt.ylabel("SPL [dB]")
    plt.xlim()
    plt.ylim()
    plt.xticks(np.arange(0, recording_time +0.1, 0.1))
    plt.yticks(np.arange(0, -60, -10))

    #Figure 5: Energy density over time at the receiver
    plt.figure(5)
    plt.title("Figure 5: Energy density over time at the receiver")
    plt.plot(t,w_rec)
    plt.ylabel('$\mathrm{Energy \ Density \ [kg/ms^2]}$')
    plt.xlabel("t [s]")

    #Figure 6: Sound pressure level stationary over the space y.
    plt.figure(6)
    t_dim = len(t)
    last_time_index = t_dim-1
    #spl_y = spl_stat[rows_r,:,dept_r]
    spl_y = spl_stat_y
    data_y = spl_y
    plt.title("Figure 6: SPL over the y axis")
    plt.plot(y,data_y)
    #plt.xticks(np.arange(0, 20, 5))
    #plt.yticks(np.arange(75, 105, 5))
    plt.ylabel('$\mathrm{Sound \ Pressure\ Level \ [dB]}$')
    plt.xlabel('$\mathrm{Distance \ along \ y \ axis \ [m]}$')

    #Figure 7: Sound pressure level stationary over the space x.
    plt.figure(7)
    t_dim = len(t)
    last_time_index = t_dim-1
    spl_x = spl_stat_x
    data_x = spl_x
    plt.title("Figure 7: SPL over the x axis")
    plt.plot(x,data_x)
    #plt.xticks(np.arange(0, 35, 5))
    #plt.yticks(np.arange(90, 97, 1))
    plt.ylabel('$\mathrm{Sound \ Pressure \ Level \ [dB]}$')
    plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$')
    
    #Figure 8: Energy density at t=recording_time over the space x.
    plt.figure(8)
    plt.title("Figure 8: Energy density over the x axis at t=recording_time")
    plt.plot(x,w_rec_x_end)
    plt.ylabel('$\mathrm{Energy \ Density \ [kg/ms^2]}$')
    plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$')
    
#%%
###############################################################################
#SAVING
###############################################################################
np.save('w_rec',w_rec)
np.save('w_rec_off',w_rec_off)
np.save('w_rec_x_t0',w_rec_x_t0)
np.save('w_rec_x_1l',w_rec_x_1l)
np.save('w_rec_x_2l',w_rec_x_2l)
np.save('w_rec_x_3l',w_rec_x_3l)
np.save('w_rec_x_5l',w_rec_x_5l)

np.save('spl_rec_x_t0',spl_rec_x_t0)
np.save('spl_rec_x_1l',spl_rec_x_1l)
np.save('spl_rec_x_2l',spl_rec_x_2l)
np.save('spl_rec_x_3l',spl_rec_x_3l)
np.save('spl_rec_x_5l',spl_rec_x_5l)

np.save('w_rec_x',w_rec_x_end)
np.save('w_rec_y',w_rec_y_end)
np.save('D0',Dx)

# Collect all variables in the workspace using locals()
FDM_8x8x8_S444_R222_AM0_alpha016_dx05 = {name: value for name, value in locals().items()}

import pickle
import dill
# Save the dictionary to a file using pickle
with open('FDM_8x8x8_S444_R222_AM0_alpha016_dx05.pkl', 'wb') as f:
    dill.dump(FDM_8x8x8_S444_R222_AM0_alpha016_dx05, f)