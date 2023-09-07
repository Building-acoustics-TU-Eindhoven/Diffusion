# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 10:39:42 2023

@author: 20225533
"""
#Code developed by Ilaria Fichera for the analysis of the FDM method Du Fort & Frankel on the 3D diffusion equation with one impulse source
import math
import matplotlib.pyplot as plt #import matplotlib as mpl
import numpy as np
from scipy.integrate import simps
from scipy import linalg
import sys
from drawnow import drawnow
from math import ceil
from math import log
#from FunctionRT import *
#from FunctionRT1 import *
from FunctionClarity import *
from FunctionDefinition import *
from FunctionCentreTime import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
import time as time
from scipy import stats
#from numpy import inf

st = time.time() #start time

#General settings
c0= 343 #sound particle velocity [m.s^-1]
rho = 1.21 #air density [Kg.m^-3] at 20°C
m_atm = 0 #air absorption coefficient [1/m] from Billon 2008 paper and Navarro paper 2012
pRef = 2 * (10**-5) #Reference pressure

#Spatial discretization
dx = 0.5 #distance between grid points x direction [m]
dy = dx #distance between grid points y direction [m]
dz = dx #distance between grid points z direction [m]

#Room dimensions
lxmin = 0 #point x starts at zero [m]
lxmax = 5.0 #point x finish at the length of the room in the x direction [m] %Length
lymin = 0 #point y starts at zero [m]
lymax = 5.0 #point y finish at the length of the room in the y direction [m] %Width
lzmin = 0 #point z starts at zero [m]
lzmax = 5.0 #point z finish at the length of the room in the x direction [m] %Height

S1,S2 = lxmax*lymax, lxmax*lymax #xy planes
S3,S4 = lxmax*lzmax, lxmax*lzmax #xz planes
S5,S6 = lymax*lzmax, lymax*lzmax #yz planes

S = lxmax*lymax*2 + lxmax*lzmax*2 + lymax*lzmax*2 #Surface [m2]
V = lxmax*lymax*lzmax #volume of the room [m^3]

x = np.arange(lxmin, lxmax+dx, dx) #mesh points in space x direction
y = np.arange(lymin, lymax+dy, dy) #mesh points in space y direction
z = np.arange(lzmin, lzmax+dz, dz) #mesh points in space z direction

Nx = len(x) #number of point in the x direction
Ny = len(y) #number of point in the y direction
Nz = len(z) #number of point in the z direction

yy, xx , zz = np.meshgrid(y,x,z) #Return coordinate matrices from coordinate vectors; create the 3D grid
#fig = plt.figure()
#ax = plt.axes(projection='3d')
#ax = ax.scatter3D(xx, yy, zz, c=zz, cmap='Greens')

#Absorption term for boundary conditions - options Sabine, Eyring and modified by Xiang
def abs_term(th,alpha):
    if th == 1:
        Absx = (c0*alpha)/4 #Sabine
    elif th == 2:
        Absx = (c0*(-log(1-alpha)))/4 #Eyring
    elif th == 3:
        Absx = (c0*alpha)/(2*(2-alpha)) #Modified by Xiang
    return Absx

th = 3 #int(input("Enter type Asbortion conditions (option 1,2,3):")) #input 1,2,3 just to understand the type of boundary chosen
alpha_1 = 0.5 #Absorption coefficient for Surface1 - Floor
alpha_2 = 0.5 #Absorption coefficient for Surface2 - Ceiling
alpha_3 = 0.5 #Absorption coefficient for Surface3 - Wall Front
alpha_4 = 0.5 #Absorption coefficient for Surface4 - Wall Back
alpha_5 = 1.0 #Absorption coefficient for Surface5 - Wall Left
alpha_6 = 0.5 #Absorption coefficient for Surface6 - Wall Right

Abs_1 = abs_term(th,alpha_1) #absorption term for S1
Abs_2 = abs_term(th,alpha_2) #absorption term for S2
Abs_3 = abs_term(th,alpha_3) #absorption term for S3
Abs_4 = abs_term(th,alpha_4) #absorption term for S4
Abs_5 = abs_term(th,alpha_5) #absorption term for S5
Abs_6 = abs_term(th,alpha_6) #absorption term for S6

alpha_average = (alpha_1*S1 + alpha_2*S2 + alpha_3*S3 + alpha_4*S4 + alpha_5*S5 + alpha_6*S6)/S #average absorption
Eq_A = alpha_1*S1 + alpha_2*S2 + alpha_3*S3 + alpha_4*S4 + alpha_5*S5 + alpha_6*S6 #equivalent absorption area of the room

#Diffusion parameters
lambda_path = (4*V)/S #mean free path for 3D
Dx = (lambda_path*c0)/3 #diffusion coefficient for proportionate rooms x direction
Dy = (lambda_path*c0)/3 #diffusion coefficient for proportionate rooms y direction
Dz = (lambda_path*c0)/3 #diffusion coefficient for proportionate rooms z direction

#Time discretization
dt = 1/32000 #distance between grid points on the time discretization [s]
recording_time = 4 #time recorded for the source [s]
recording_steps = ceil(recording_time/dt) #number of time steps to consider in the calculation
t = np.arange(0, recording_time, dt) #mesh point in time

beta_zero_x = (2*Dx*dt)/(dx**2) #mesh number in x direction
beta_zero_y = (2*Dy*dt)/(dy**2) #mesh number in x direction
beta_zero_z = (2*Dz*dt)/(dz**2) #mesh number in x direction
beta_zero = beta_zero_x + beta_zero_y + beta_zero_z #beta_zero is the condition for all the directions deltax, deltay and deltaz.
 
#Condition for the model to be unconditionally stable
beta_zero_condition = ((beta_zero**2) - 1)/(1+(beta_zero**2)+(2*beta_zero)) #from Navarro 2012 paper
if beta_zero_condition >1:
    print("aa! errors! Check beta condition")

#Frequency resolution & spatial parameters
fsample = 1/dt #frequency spatial resolution (sampling period)

#Finding index in meshgrid of the source position
x_source = 2.5 #int(ceil(Nx/2))#4 #position of the source in the x direction [m]
y_source = 2.5 #int(ceil(Ny/2))#4 #position of the source in the y direction [m]
z_source = 2.5 #int(ceil(Nz/2))#4 #position of the source in the z direction [m]
coord_source = [x_source , y_source, z_source] #coordinates of the source position in an list
rows_s = np.argmin(abs(xx[:,0,0] - coord_source[0])) #Find index of grid point with minimum distance from source along x direction
cols_s = np.argmin(abs(yy[0,:,0] - coord_source[1])) #Find index of grid point with minimum distance from source along y direction
dept_s = np.argmin(abs(zz[0,0,:] - coord_source[2])) #Find index of grid point with minimum distance from source along z direction

#Finding index in meshgrid of the receiver position
x_rec = 2.5#int(ceil(Nx/4)) #position of the receiver in the x direction [m]
y_rec = 1.0#int(ceil(Nx/4)) #position of the receiver in the y direction [m]
z_rec = 1.0#int(ceil(Nx/4)) #position of the receiver in the z direction [m]
coord_receiver = [x_rec,y_rec,z_rec] #coordinates of the receiver position in an list
rows_r = np.argmin(abs(xx[:,0,0] - coord_receiver[0])) #Find index of grid point with minimum distance from receiver along x direction
cols_r = np.argmin(abs(yy[0,:,0] - coord_receiver[1])) #Find index of grid point with minimum distance from receiver along y direction
dept_r = np.argmin(abs(zz[0,0,:] - coord_receiver[2])) #Find index of grid point with minimum distance from receiver along z direction

#Gaussian distribution
#n_gau = 1 #??
#sigmax_gau = n_gau*dx #??
#Ax_gau = 5e-3 #??
#ax_gau = 1/(2*(sigmax_gau**2)) #??
#power_x = (np.subtract(xx,x[x_source-1]))**2 #?? 
#power_y = (np.subtract(yy,y[y_source-1]))**2 #??
#power_z = (np.subtract(zz,z[z_source-1]))**2 #??
#P = np.multiply(Ax_gau, np.exp(np.multiply(-ax_gau,(power_x+power_y+power_z)))) #??
P = np.zeros((Nx,Ny,Nz)) #matrix of zeros for source
radius_s = 0.62 #[m]
Vs = dx*dy*dz #4/3*math.pi*(radius_s**3) #????????????????????????????
Ws = 0.01 # power of the source in [W]
P[rows_s, cols_s, dept_s] = Ws/Vs

dist = math.sqrt((abs(x_rec - x_source))**2 + (abs(y_rec - y_source))**2 + (abs(z_rec - z_source))**2) #distance between source and receiver

dist_x = np.sqrt((xx[:,cols_r,dept_r] - x_source)**2 + (yy[rows_r,cols_r,dept_r] - y_source)**2 + (zz[rows_r,cols_r,dept_r] - z_source)**2)
dist_y = np.sqrt((xx[rows_r,cols_r,dept_r] - x_source)**2 + (yy[rows_r,:,dept_r] - y_source)**2 + (zz[rows_r,cols_r,dept_r] - z_source)**2)

#dist_y = math.sqrt(np.add((abs(x_rec - x_source))**2, np.power((abs(y[:] - y_source)),2) , (abs(z_rec - z_source))**2)) 

w_new = np.zeros((Nx,Ny,Nz)) #unknown w at new time level (n+1)
w = w_new #w at n level
w_old = w #w_old at n-1 level

w_rec = np.arange(0,recording_time,dt) #energy density at the receiver; mesh for the graph later on

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
            np.divide((2*dt*P),(1+beta_zero)) + \
                np.divide((np.multiply(beta_zero_x,(w_iplus1+w_iminus1))),(1+beta_zero)) + \
                    np.divide((np.multiply(beta_zero_y,(w_jplus1+w_jminus1))),(1+beta_zero)) + \
                        np.divide((np.multiply(beta_zero_z,(w_kplus1+w_kminus1))),(1+beta_zero))
    
     
    #Insert boundary conditions  
    w_new[0,:,:] = np.divide((4*w_new[1,:,:] - w_new[2,:,:]),(3+((2*Abs_1*dx)/Dx))) #boundary condition at x=0, any y, any z
    w_new[-1,:,:] = np.divide((4*w_new[-2,:,:] - w_new[-3,:,:]),(3+((2*Abs_2*dx)/Dx))) #boundary condition at lx=lxmax, any y, any z


    w_new[:,0,:] = np.divide((4*w_new[:,1,:] - w_new[:,2,:]),(3+((2*Abs_3*dx)/Dy))) #boundary condition at y=0, any x, any z
    w_new[:,-1,:] = np.divide((4*w_new[:,-2,:] - w_new[:,-3,:]),(3+((2*Abs_4*dx)/Dy))) #boundary condition at at ly=lymax, any x, any z
 
    
    w_new[:,:,0] = np.divide((4*w_new[:,:,1] - w_new[:,:,2]),(3+((2*Abs_5*dx)/Dz))) #boundary condition at z=0, any x, any y
    w_new[:,:,-1] = np.divide((4*w_new[:,:,-2] - w_new[:,:,-3]),(3+((2*Abs_6*dx)/Dz))) #boundary condition at at lz=lzmax, any x, any y
    
    sdl = 10*np.log10(abs(w_new),where=abs(w_new)>0) #sound density level

    #Update w before next step
    w_old = w #The w at n step becomes the w at n-1 step
    w = w_new #The w at n+1 step becomes the w at n step

    #w_rec is the energy density at the receiver specifically
    w_rec[steps] = w_new[rows_r, cols_r, dept_r] #energy density at the receiver is equal to the energy density new calcuated in time
    
    print(time_steps)

w_rec_x = w_new[:, cols_r, dept_r]  
w_rec_y = w_new[rows_r, :, dept_r]

spl_stat_x = 10*np.log10(rho*c0*(((Ws)/(4*math.pi*(dist_x**2))) + ((abs(w_rec_x)*c0)))/(pRef**2))
spl_stat_y = 10*np.log10(rho*c0*(((Ws)/(4*math.pi*(dist_y**2))) + ((abs(w_rec_y)*c0)))/(pRef**2)) #It should be the spl stationary
#spl_stat_norm = 10*np.log10(rho*c0*((P/(4*math.pi*dist**2)) + ((abs(w_new))*c0)/(pRef**2))/ np.max(rho*c0*((P/(4*math.pi*dist**2)) + ((abs(w_new))*c0)/(pRef**2)))) #It should be the spl stationary

spl_norm = 10*np.log10((((abs(w_rec))*rho*(c0**2))/(pRef**2)) / np.max(((abs(w_rec))*rho*(c0**2))/(pRef**2))) #normalised to maximum to 0dB

#Figure 1: Decay of SPL in the recording_time at the receiver
plt.figure(1) 
press_r = ((abs(w_rec))*rho*(c0**2))
spl = 10*np.log10(((abs(w_rec))*rho*(c0**2))/(pRef**2)) #,where=press_r>0

plt.plot(t,spl) #plot sound pressure level with Pref = (2e-5)**5
plt.title("SPL over time at the receiver")
plt.xlabel("t")
plt.ylabel("SPL")
plt.xlim()
plt.ylim()
plt.xticks(np.arange(0, recording_time +0.1, 0.5))
#plt.yticks(np.arange(0, 120, 20))

#Figure 2: Decay of SPL in the recording_time normalised to maximum 0dB
plt.figure(2)
plt.title("Normalised SPL over time at the receiver")
plt.plot(t,spl_norm)
plt.xlabel("t")
plt.ylabel("SPL")
plt.xlim()
plt.ylim()
plt.xticks(np.arange(0, recording_time +0.1, 0.1))
plt.yticks(np.arange(0, -60, -10))

#Figure 3: Energy density over time at the receiver
plt.figure(3)
plt.title("Energy density over time at the receiver")
plt.plot(t,w_rec)

#Figure 4: Sound pressure level stationary over the space y.
plt.figure(4)
t_dim = len(t)
last_time_index = t_dim-1
#spl_y = spl_stat[rows_r,:,dept_r]
spl_y = spl_stat_y
data_y = spl_y
plt.title("SPL over the y axis")
plt.plot(y,data_y)
plt.xticks(np.arange(0, 20, 5))
plt.yticks(np.arange(75, 105, 5))
plt.ylabel('$\mathrm{Sound \ Pressure\ Level \ [dB]}$')
plt.xlabel('$\mathrm{Distance \ along \ y \ axis \ [m]}$')

#Figure 5: Sound pressure level stationary over the space x.
plt.figure(5)
t_dim = len(t)
last_time_index = t_dim-1
spl_x = spl_stat_x
data_x = spl_x
plt.title("SPL over the x axis")
plt.plot(x,data_x)
plt.xticks(np.arange(0, 35, 5))
plt.yticks(np.arange(65, 105, 5))
plt.ylabel('$\mathrm{Sound \ Pressure \ Level \ [dB]}$')
plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$')

et = time.time() #end time
elapsed_time = et - st