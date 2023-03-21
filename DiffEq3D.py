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
from FunctionRT import *
from FunctionClarity import *
from FunctionDefinition import *
from FunctionCentreTime import *

#General settings
c0= 343 #sound particle velocity [m.s^-1]
rho = 1.21 #air density [Kg.m^-3] at 20°C
m_atm = 0 #air absorption coefficient [1/m] from Billon 2008 paper and Navarro paper 2012
pRef = 2 * (10**-5) #Reference pressure

#Spatial discretization
dx = 0.5 #distance between grid points x direction [m]
dy = 0.5 #distance between grid points y direction [m]
dz = 0.5 #distance between grid points z direction [m]

#Time discretization
dt = 0.0001 #distance between grid points on the time discretization [s]
recording_time = 4 #time recorded for the source [s]
recording_steps = ceil(recording_time/dt) #number of time steps to consider in the calculation
t = np.arange(0, recording_time, dt) #mesh point in time

#Frequency resolution & spatial parameters
fspatial = 1/dt #frequency spatial resolution (sampling period)

#Room dimensions
lxmin = 0 #point x starts at zero [m]
lxmax = 10.0 #point x finish at the length of the room in the x direction [m] %Length
lymin = 0 #point y starts at zero [m]
lymax = 10.0 #point y finish at the length of the room in the y direction [m] %Width
lzmin = 0 #point z starts at zero [m]
lzmax = 10.0 #point z finish at the length of the room in the x direction [m] %Height

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

th = 2 #int(input("Enter type Asbortion conditions (option 1,2,3):")) #input 1,2,3 just to understand the type of boundary chosen
alpha_1 = 0.1#1/6 #Absorption coefficient for Surface1
alpha_2 = 0.1 #Absorption coefficient for Surface2
alpha_3 = 0.1 #Absorption coefficient for Surface3
alpha_4 = 0.1 #Absorption coefficient for Surface4
alpha_5 = 0.1 #Absorption coefficient for Surface5
alpha_6 = 0.1 #Absorption coefficient for Surface6

Abs_1 = round(abs_term(th,alpha_1),4) #absorption term for S1
Abs_2 = round(abs_term(th,alpha_2),4) #absorption term for S2
Abs_3 = round(abs_term(th,alpha_3),4) #absorption term for S3
Abs_4 = round(abs_term(th,alpha_4),4) #absorption term for S4
Abs_5 = round(abs_term(th,alpha_5),4) #absorption term for S5
Abs_6 = round(abs_term(th,alpha_6),4) #absorption term for S6

alpha_average = (alpha_1*S1 + alpha_2*S2 + alpha_3*S3 + alpha_4*S4 + alpha_5*S5 + alpha_6*S6)/S #average absorption
Eq_A = alpha_1*S1 + alpha_2*S2 + alpha_3*S3 + alpha_4*S4 + alpha_5*S5 + alpha_6*S6 #equivalent absorption area of the room

#Diffusion parameters
lambda_path = round(4*V/S,4) #mean free path for 3D
Dx = round((lambda_path*c0)/3,4) #diffusion coefficient for proportionate rooms x direction
Dy = round((lambda_path*c0)/3,4) #diffusion coefficient for proportionate rooms y direction
Dz = round((lambda_path*c0)/3,4) #diffusion coefficient for proportionate rooms z direction

beta_zero_x = round((2*Dx*dt)/(dx**2),4) #mesh number in x direction
beta_zero_y = round((2*Dy*dt)/(dy**2),4) #mesh number in x direction
beta_zero_z = round((2*Dz*dt)/(dz**2),4) #mesh number in x direction
beta_zero = beta_zero_x + beta_zero_y + beta_zero_z #beta_zero is the condition for all the directions deltax, deltay and deltaz.
 
#Condition for the model to be unconditionally stable
beta_zero_condition = round(((beta_zero**2) - 1)/(1+(beta_zero**2)+(2*beta_zero)),4) #from Navarro 2012 paper
if beta_zero_condition >1:
    print("aa! errors! Check beta condition")

#Set initial condition - Source Info (excitation with Gaussian) 
Ws=10**-2#0.005 #Source point power [Watts] interrupted after 2seconds; 10^-2 value taken from Jing 2007; correspondent to a SWL of 100dB
#Vs=0.0001
#Vs=round(4/3*round(pi,4)*(dx**3),10) #Source volume
#w1 = round(Ws/Vs,4) #power density of the source [Watts/(m^3))]

sourceon_time =  1 #time that the source is on before interrupting [s]
sourceon_steps = ceil(sourceon_time/dt) #time steps at which the source is calculated/considered in the calculation
s1 = np.multiply(Ws,np.ones(sourceon_steps)) #energy density of source number 1 at each time step position #does the source not need to be only at the time 0 to 2seconds and after that there should not be any source term? Yes
source1 = np.append(s1, np.zeros(recording_steps-sourceon_steps)) #This would be equal to s1 if and only if recoding_steps = sourceon_steps

#np.around(s1, 4, s1) #evenly round to the given number of decimals

x_source = 9#ceil(Nx/2)#4.0 #position of the source in the x direction [m]
y_source = 9#ceil(Nx/2) #position of the source in the y direction [m]
z_source = 9#ceil(Nx/2) #position of the source in the z direction [m]

coord_source = [x_source , y_source, z_source] #coordinates of the source position in an list

for i in range(0,Nx):
    np.around(xx[i], 2, xx[i]) #evenly round to the given number of decimals
    np.around(yy[i], 2, yy[i]) #evenly round to the given number of decimals
    np.around(zz[i], 2, zz[i]) #evenly round to the given number of decimals
coord_sourceRound0 = round(coord_source[0], 2)
coord_sourceRound1 = round(coord_source[1], 2)
coord_sourceRound2 = round(coord_source[2], 2)
index_source = (np.argwhere((xx == coord_sourceRound0) & (yy == coord_sourceRound1) & (zz == coord_sourceRound2)))[0] #finding the index of the source in the meshgrid
rows_s, cols_s, dept_s = index_source[0], index_source[1], index_source[2] #the row index is the first item in the list; the col index is the second item in the list, the dept is the third item in the list

#Set initial condition - Receiver Info
x_rec = ceil(Nx/4) #2.0 #position of the receiver in the x direction [m]
y_rec = ceil(Nx/4) #position of the receiver in the y direction [m]
z_rec = ceil(Nx/4) #position of the receiver in the z direction [m]

dist = math.sqrt((abs(x_rec - x_source))**2 + (abs(y_rec - y_source))**2 + (abs(z_rec - z_source))**2) #distance between source and receiver

coord_receiver = [x_rec,y_rec,z_rec] #coordinates of the receiver position in an list
index_receiver = (np.argwhere((xx==coord_receiver[0]) & (yy==coord_receiver[1]) & (zz==coord_receiver[2])))[0] #finding the index of the receiver in the meshgrid
rows_r, cols_r, dept_r = index_receiver[0], index_receiver[1], index_receiver[2] #the row index is the first item in the list; the col index is the second item in the list,the dept is the third item in the list

s = np.zeros((Nx,Ny,Nz)) #matrix of zeros for source 
s[rows_s, cols_s, dept_s] = source1[1] #at the index where the different between the source and x is zero, the source value is the energy density of the source, for all the other values it is zero.

w_new = np.zeros((Nx,Ny,Nz)) #unknown w at new time level (n+1)
w = w_new #w at n level
w_old = w #w_old at n-1 level

w_rec = np.arange(0,recording_time,dt) #energy density at the receiver; mesh for the graph later on

#Function to draw figure1
def draw_fig1():
    plt.plot(x, w_new) #plot the figure energy density in space x
    plt.grid()
    plt.xlabel("x")
    plt.ylabel("w")
    plt.xlim(0, 1)
    plt.ylim(1e-10, 1e-1)
    plt.xticks(np.arange(0, 1+0.1, 0.1))
    plt.yticks([1e-10, 1e-8, 1e-6, 1e-4, 1e-2])
  
#Function to draw figure2
def draw_fig2():
    plt.plot(x, sdl) #plot the figure energy density in decibel in space x
    plt.grid()
    plt.xlabel("x")
    plt.ylabel("w")
    plt.xlim(0, 1)
    plt.ylim(1e-10, 1e-1)
    plt.xticks(np.arange(0, 1+0.1, 0.1))
    plt.yticks([1e-10, 1e-8, 1e-6, 1e-4, 1e-2])

def draw_fig():
    plt.imshow(w_new.transpose(), vmin=w_new.min(), vmax=w_new.max());

#Computing w;
for steps in range(0, recording_steps):
    #Compute w at inner mesh points
    time = steps*dt #total time for the calculation
    s[rows_s, cols_s, dept_s] = source1[steps] #array of zero of the source apart from the index_dist_source = energy density of the source at each step position
    #w_trans = np.transpose(w) #transpose of w could be the trans???
    
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
    w_new[0,:,:] = np.divide((4*w_new[1,:,:] - w_new[2,:,:]),(3+((2*Abs_1*dx)/Dx))) #boundary condition at x=0, any y, any z
    w_new[-1,:,:] = np.divide((4*w_new[-2,:,:] - w_new[-3,:,:]),(3+((2*Abs_2*dx)/Dx))) #boundary condition at lx=lxmax, any y, any z
    
    w_new[:,0,:] = np.divide((4*w_new[:,1,:] - w_new[:,2,:]),(3+((2*Abs_3*dx)/Dy))) #boundary condition at y=0, any x, any z
    w_new[:,-1,:] = np.divide((4*w_new[:,-2,:] - w_new[:,-3,:]),(3+((2*Abs_4*dx)/Dy))) #boundary condition at at ly=lymax, any x, any z
 
    w_new[:,:,0] = np.divide((4*w_new[:,:,1] - w_new[:,:,2]),(3+((2*Abs_5*dx)/Dz))) #boundary condition at z=0, any x, any y
    w_new[:,:,-1] = np.divide((4*w_new[:,:,-2] - w_new[:,:,-3]),(3+((2*Abs_6*dx)/Dz))) #boundary condition at at lz=lzmax, any x, any y
 
    sdl = 10*np.log10(abs(w_new),where=abs(w_new)>0) #sound density level
    #if (steps % 100 == 0): #draw only on certain steps and not all the steps
    #    drawnow(draw_fig)
    #Update w before next step
    w_old = w #The w at n step becomes the w at n-1 step
    w = w_new #The w at n+1 step becomes the w at n step

    #w_rec is the energy density at the receiver specifically
    w_rec[steps] = w_new[rows_r, cols_r,dept_r] #energy density at the receiver is equal to the energy density new calcuated in time
    #drawnow(draw_fig1)
    #drawnow(draw_fig2)

#Figure 3: Decay of SPL in the recording_time
plt.figure(3) 
press_r = ((abs(w_rec))*rho*(c0**2))
#max_press_r = np.max(((abs(w_rec))*rho*(c0**2))/(pRef**2))
spl = 10*np.log10(((abs(w_rec))*rho*(c0**2))/(pRef**2)) #,where=press_r>0
plt.plot(t,spl) #plot sound pressure level with Pref = (2e-5)**5
plt.xlabel("t")
plt.ylabel("SPL")
plt.xlim()
plt.ylim()
plt.xticks(np.arange(0, recording_time +0.1, 0.5))
#plt.yticks(np.arange(0, 120, 20))

sch_db = 10*np.log10((((abs(w_rec))*rho*(c0**2))/(pRef**2)) / np.max(((abs(w_rec))*rho*(c0**2))/(pRef**2))) #normalised to maximum to 0dB

#Figure 4: Decay of SPL in the recording_time normalised to maximum 0dB
plt.figure(4)
plt.plot(t,sch_db)
plt.xlabel("t")
plt.ylabel("SPL")
plt.xlim()
plt.ylim()
plt.xticks(np.arange(0, recording_time +0.1, 0.5))
#plt.yticks(np.arange(0, -120, -20))


error = (dt**2) * (dx**(-2)) * (10**(-8))

t60 = t60_decay(t, sch_db, fspatial, rt='t30') #called function for calculation of t60 [s]
edt = t60_decay(t, sch_db, fspatial, rt='edt') #called function for calculation of edt [s]
c80 = clarity(t60, V, Eq_A, S, c0, dist) #called function for calculation of c80 [dB]
d50 = definition(t60, V, Eq_A, S, c0, dist) #called function for calculation of d50 [%]
ts = centretime(t60, Eq_A, S) #called function for calculation of ts [ms]
