# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 10:39:42 2023

@author: 20225533
"""
#Code developed by Ilaria Fichera for the analysis of the FDM method Du Fort & Frankel on the 1D diffusion equation with one impulse source
import math
import matplotlib.pyplot as plt #import matplotlib as mpl
import numpy as np
from scipy.integrate import simps
from scipy import linalg
import sys
from drawnow import drawnow
from math import ceil
from math import log
from math import pi
from scipy.io import wavfile
from scipy import stats
from acoustics.utils import _is_1d
from acoustics.signal import bandpass
from acoustics.bands import (_check_band_type, octave_low, octave_high, third_low, third_high)
from FunctionRT import *
from scipy.spatial import KDTree

#General settings
c0= 343 #sound particle velocity [m.s^-1]
rho = 1.21 #air density [Kg.m^-3] at 20°C
m_atm = 0 #air absorption coefficient [1/m] from Billon 2008 paper
pRef = 2 * (10**-5) #Reference pressure

#Spatial discretization
dx = 0.1 #distance between grid points x direction [m]
dy = 0.1 #distance between grid points y direction [m]

#Time discretization
dt = 0.0001 #distance between grid points on the time discretization [s]
recording_time = 3 #time recorded for the source [s]
recording_steps = ceil(recording_time/dt) #number of time steps to consider in the calculation
t = np.arange(0, recording_time, dt) #mesh point in time

#Frequency resolution & spatial parameters
fspatial = 1/dt #frequency spatial resolution (sampling period)

#Room dimensions
lxmin = 0 #point x starts at zero [m]
lxmax = 39.7 #point x finish at the length of the room in the x direction [m] %Length
lymin = 0 #point y starts at zero [m]
lymax = 2.8 #point y finish at the length of the room in the y direction [m] %Width

S_2D = lxmax*lymax #Surface [m2]
P_2D = 2*(lxmax+lymax) #Perimeter [m]

x = np.arange(lxmin, lxmax+dx, dx) #mesh points in space x direction
y = np.arange(lymin, lymax+dy, dy) #mesh points in space y direction

Nx = len(x) #number of point in the x direction
Ny = len(y) #number of point in the y direction

yy , xx = np.meshgrid(y,x) #Return coordinate matrices from coordinate vectors; create the 2D grid
plt.plot(xx, yy, marker='.', color='k', linestyle='none') #plots the 2D grid

#Absorption term for boundary conditions - options Sabine, Eyring and modified by Xiang
def abs_term(th,alpha):
    if th == 1:
        Absx = (c0*alpha)/4 #Sabine
    elif th == 2:
        Absx = (c0*(-log(1-alpha)))/4 #Eyring
    elif th == 3:
        Absx = (c0*alpha)/(2*(2-alpha)) #Modified by Xiang
    return Absx

th = 1#int(input("Enter type Asbortion conditions (option 1,2,3):")) #input 1,2,3 just to understand the type of boundary chosen
alpha_x = 0.99#float(input("Enter absorption coefficient:")) #enter absorption coefficient uniform for all the surfaces in the x direction and for one frequency only
alpha_y = 0.1#float(input("Enter absorption coefficient:")) #enter absorption coefficient uniform for all the surfaces in the y direction and for one frequency only
Abs_x = round(abs_term(th,alpha_x),4) #absorption term for x
Abs_y = round(abs_term(th,alpha_y),4) #absorption term for y

#Diffusion parameters
lambda_path_2D = round(pi*S_2D/P_2D,4) #mean free path for 2D
Dx = round((lambda_path_2D*c0)/3,4) #diffusion coefficient for proportionate rooms x direction
Dy = round((lambda_path_2D*c0)/3,4) #diffusion coefficient for proportionate rooms y direction

beta_zero_x = round((2*Dx*dt)/(dx**2),4) #mesh number in x direction
beta_zero_y = round((2*Dy*dt)/(dy**2),4) #mesh number in x direction
beta_zero = beta_zero_x + beta_zero_y #beta_zero is the condition for all the directions deltax, deltay and deltaz.
 
#Condition for the model to be unconditionally stable
beta_zero_condition = round(((beta_zero**2) - 1)/(1+(beta_zero**2)+(2*beta_zero)),4) #from Navarro 2012 paper
if beta_zero_condition >1:
    print("aa! errors! Check beta condition")

#Set initial condition - Source Info (excitation with Gaussian) 
Ws=10**-2 #Source point power [Watts] interrupted after 2seconds; value taken from Jing 2007; correspondent to a SWL of 100dB
Vs=0.0001
#Vs=round(4/3*round(pi,4)*(dx**3),10) #Source volume
w1 = round(Ws/Vs,4) #power density of the source [Watts/(m^3))]

sourceon_time =  0.5 #time that the source is on before interrupting [s]
sourceon_steps = ceil(sourceon_time/dt) #time steps at which the source is calculated/considered in the calculation
s1 = np.multiply(w1,np.ones(sourceon_steps)) #energy density of source number 1 at each time step position
#####does the source not need to be only at the time 0 to 2seconds and after that there should not be any source term? Yes
source1 = np.append(s1, np.zeros(recording_steps-sourceon_steps)) #This would be equal to s1 if and only if recoding_steps = sourceon_steps

np.around(s1, 4, s1) #evenly round to the given number of decimals

x_source = 10.2 #position of the source in the x direction [m]
y_source = 1.2 #position of the source in the x direction [m]

coord_source = [x_source , y_source] #coordinates of the source position in an list

for i in range(0,Nx):
    np.around(xx[i], 2, xx[i]) #evenly round to the given number of decimals
    np.around(yy[i], 2, yy[i]) #evenly round to the given number of decimals
coord_sourceRound0 = round(coord_source[0], 2)
coord_sourceRound1 = round(coord_source[1], 2)
index_source = (np.argwhere((xx == coord_sourceRound0) & (yy == coord_sourceRound1)))[0] #finding the index of the source in the meshgrid
rows_s, cols_s = index_source[0], index_source[1] #the row index is the first item in the list; the col index is the second item in the list

#Set initial condition - Receiver Info
x_rec = 20.2
y_rec = 2.2 #lymax/2
#index_receiver = kdt.query([[x_rec,y_rec]])[1] #index of the receiver position
coord_receiver = [x_rec,y_rec] #coordinates of the receiver position in an list
index_receiver = (np.argwhere((xx==coord_receiver[0]) & (yy==coord_receiver[1])))[0] #finding the index of the receiver in the meshgrid
rows_r, cols_r = index_receiver[0], index_receiver[1] #the row index is the first item in the list; the col index is the second item in the list

s = np.zeros((Nx,Ny)) #matrix of zeros for source 
s[rows_s, cols_s] = source1[1] #at the index where the different between the source and x is zero, the source value is the energy density of the source, for all the other values it is zero.

w_new = np.zeros((Nx,Ny)) #unknown w at new time level (n+1)
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
    s[rows_s, cols_s] = source1[steps] #array of zero of the source apart from the index_dist_source = energy density of the source at each step position
    #w_trans = np.transpose(w) #transpose of w could be the trans???
    
    #In the x direction
    w_iminus1 = w[0:Nx-1, 0:Ny] 
    w_iminus1 = np.vstack([np.ones(Ny),w_iminus1])
    w_iplus1 = w[1:Nx, 0:Ny]
    w_iplus1 = np.vstack([w_iplus1,np.ones(Ny)])
    
    #In the y direction
    w_jminus1 = w[0:Nx, 0:Ny-1] 
    w_jminus1 = np.hstack([np.ones((Nx,1)),w_jminus1])
    w_jplus1 = w[0:Nx, 1:Ny]
    w_jplus1 = np.hstack([w_jplus1,np.ones((Nx,1))])
    
    #Compututing w_new (w at n+1 time step)
    w_new = np.divide((np.multiply(w_old,(1-beta_zero))),(1+beta_zero)) - np.divide((2*dt*c0*m_atm*w),(1+beta_zero)) + np.divide((2*dt*s),(1+beta_zero)) + np.divide((np.multiply(beta_zero_x,(w_iplus1+w_iminus1))),(1+beta_zero)) + np.divide((np.multiply(beta_zero_y,(w_jplus1+w_jminus1))),(1+beta_zero))
    
    #Insert boundary conditions
    w_new[0,0:Ny] = np.divide((4*w_new[1,0:Ny] - w_new[2,0:Ny]),(3+((2*Abs_x*dx)/Dx))) #boundary condition at x=0, any y
    w_new[Nx-1,0:Ny] = np.divide((4*w_new[Nx-2,0:Ny] - w_new[Nx-3,0:Ny]),(3+((2*Abs_x*dx)/Dx))) #boundary condition at lx=lxmax, any y
    
    w_new[0:Nx,0] = np.divide((4*w_new[0:Nx,1] - w_new[0:Nx,2]),(3+((2*Abs_y*dx)/Dy))) #boundary condition at y=0, any x
    w_new[0:Nx, Ny-1] = np.divide((4*w_new[0:Nx,Ny-2] - w_new[0:Nx,Ny-3]),(3+((2*Abs_y*dx)/Dy))) #boundary condition at at ly=lymax, any x
 
    sdl = 10*np.log10(abs(w_new),where=abs(w_new)>0) #sound density level
    #if (steps % 100 == 0): #draw only on certain steps and not all the steps
    #    drawnow(draw_fig)
    #Update w before next step
    w_old = w #The w at n step becomes the w at n-1 step
    w = w_new #The w at n+1 step becomes the w at n step

    #w_rec is the energy density at the receiver specifically
    w_rec[steps] = w_new[rows_r, cols_r] #energy density at the receiver is equal to the energy density new calcuated in time
    #drawnow(draw_fig1)
    #drawnow(draw_fig2)

#Figure 3: Decay of SPL in the recording_time
plt.figure(3) 
press_r = ((abs(w_rec))*rho*(c0**2))/(pRef**2)
#max_press_r = np.max(((abs(w_rec))*rho*(c0**2))/(pRef**2))
spl = 10*np.log10(press_r) #,where=press_r>0
plt.plot(t,spl) #plot sound pressure level with Pref = (2e-5)**5
plt.xlabel("t")
plt.ylabel("SPL")
plt.xlim()
plt.ylim()
plt.xticks(np.arange(0, recording_time +0.1, 0.5))
plt.yticks(np.arange(0, 120, 20))

sch_db = 10*np.log10((((abs(w_rec))*rho*(c0**2))/(pRef**2)) / np.max(((abs(w_rec))*rho*(c0**2))/(pRef**2))) #normalised to maximum to 0dB

#Figure 4: Decay of SPL in the recording_time normalised to maximum 0dB
plt.figure(4)
plt.plot(t,sch_db)
plt.xlabel("t")
plt.ylabel("SPL")
plt.xlim()
plt.ylim()
plt.xticks(np.arange(0, recording_time +0.1, 0.5))
plt.yticks(np.arange(0, -120, -20))

sch_db[sch_db == -np.inf] = 0 #replace the -inf with zero in the decay

bands = np.array([63,125,250,500,1000,2000,4000]) #array of frequency bands 
t60 = t60_decay(t, sch_db, fspatial, rt='t30') #called function for calculation of t60 









