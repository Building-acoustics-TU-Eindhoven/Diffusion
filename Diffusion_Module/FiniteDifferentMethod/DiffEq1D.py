# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 09:26:51 2023

@author: 20225533
"""
#Code developed by Ilaria Fichera for the analysis of the FDM method Du Fort & Frankel on the 1D diffusion equation with one impulse source
from math import ceil
from math import log
from math import pi

import numpy as np
# from drawnow import drawnow

from scipy import stats
from scipy import linalg
from scipy.io import wavfile
from scipy.integrate import simps

import matplotlib.pyplot as plt  # import matplotlib as mpl

from acoustics.room import t60_impulse
# from acoustics.signal import bandpass
# from acoustics.utils import _is_1d
# from acoustics.bands import (_check_band_type, octave_low, octave_high, third_low, third_high)

from FunctionRT import t60_decay

#General settings
c0= 343 #sound particle velocity [m.s^-1]
rho = 1.21 #air density [Kg.m^-3] at 20°C
m_atm = 0 #air absorption coefficient [1/m] from Billon 2008 paper
pRef = 2 * (10**-5) #Reference pressure

#Spatial discretization
dx = 0.1 #distance between grid points [m]

#Time discretization
dt = 0.0001 #distance between grid points on the time discretization [s]
recording_time = 3 #time recorded for the source [s]
recording_steps = ceil(recording_time/dt) #number of time steps to consider in the calculation
t = np.arange(0, recording_time, dt) #mesh point in time

#Frequency resolution & spatial parameters
fspatial = 1/dt #frequency spatial resolution (sampling period)

#Room dimensions
lxmin = 0 #point x starts at zero [m]
lxmax = 39.7 #point x finish at the length of the room in the x direction [m]
x = np.arange(lxmin, lxmax+dx, dx) #mesh points in space
N = len(x) #number of point in the x direction

#Absorption term for boundary conditions - options Sabine, Eyring and modified by Xiang
def abs_term(th,alpha):
    if th == 1:
        Absx = (c0*alpha)/4 #Sabine
    elif th == 2:
        Absx = (c0*(-log(1-alpha)))/4 #Eyring
    elif th == 3:
        Absx = (c0*alpha)/(2*(2-alpha)) #Modified by Xiang
    return Absx

th = int(input("Enter type Asbortion conditions (option 1,2,3):")) #input 1,2,3 just to understand the type of boundary chosen
alpha = float(input("Enter absorption coefficient:")) #enter absorption coefficient uniform for all the surfaces and for one frequency only
Abs = round(abs_term(th,alpha),4) #absorption term

#Diffusion parameters
lambda_path_1D = lxmax #mean free path for 1D
D = round((lambda_path_1D*c0)/3,4) #diffusion coefficient for proportionate rooms

beta_zero_x = round((2*D*dt)/(dx**2),4) #mesh number
beta_zero = beta_zero_x #beta_zero is the condition for all the directions deltax, deltay and deltaz.
 
#Condition for the model to be unconditionally stable
beta_zero_condition = round(((beta_zero**2) - 1)/(1+(beta_zero**2)+(2*beta_zero)),4) #from Navarro 2012 paper
if beta_zero_condition >1:
    print("aa! errors! Check beta condition")

#Set initial condition - Source Info (excitation with Gaussian) 
Ws=10**-2   #Source point power [Watts] interrupted after 2seconds
Vs=round(4/3*round(pi,4)*(dx**3),4) #Source volume
w1 = round(Ws/Vs,4) #energy density of the source [Watts/(m∙(s^3))]

sourceon_time =  1 #time that the source is on before interrupting [s]
sourceon_steps = ceil(sourceon_time/dt) #time steps at which the source is calculated/considered in the calculation
s1 = np.multiply(w1,np.ones(sourceon_steps)) #energy density of source number 1 at each time step position
#####does the source not need to be only at the time 0 to 1seconds and after that there should not be any source term? Yes
source1 = np.append(s1, np.zeros(recording_steps-sourceon_steps)) #This would be equal to s1 if and only if recoding_steps = sourceon_steps

x_source = 0.5 #position of the source in the x direction [m]
source_distance = x-x_source #distance between the x value (where the calculation is happening) and the source position
index_dist_source = [i for i in range(source_distance.size) if source_distance[i] == 0] #shorten of the for loop just below
#for i in range(source_distance.size):
#    if source_distance[i] == 0:
#        index_dist_source = i #get the index at which the source is on

#Set initial condition - Receiver Info
x_rec = 10
rec_distance = x-x_rec #distance between the x value (where the calculation is happening) and the receiver position
index_dist_rec = [j for j in range(rec_distance.size) if rec_distance[j] == 0] #shorten of the for loop just below
#for j in range(rec_distance.size):
#    if rec_distance[j] == 0:
#        index_dist_rec = j

s = np.zeros(N) #array of zeros for source 
s[index_dist_source] = source1[index_dist_source] #at the index where the different between the source and x is zero, the source value is the energy density of the source, for all the other values it is zero.

w_new = np.zeros(N) #unknown w at new time level (n+1)
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
    plt.plot(x, SDL) #plot the figure energy density in decibel in space x
    plt.grid()
    plt.xlabel("x")
    plt.ylabel("w")
    plt.xlim(0, 1)
    plt.ylim(1e-10, 1e-1)
    plt.xticks(np.arange(0, 1+0.1, 0.1))
    plt.yticks([1e-10, 1e-8, 1e-6, 1e-4, 1e-2])

#Computing w;
for steps in range(0, recording_steps):
    #Compute w at inner mesh points
    time = steps*dt #total time for the calculation
    s[index_dist_source] = source1[steps] #array of zero of the source apart from the index_dist_source = energy density of the source at each step position
    #w_trans = np.transpose(w) #transpose of w could be the trans???
    
    w_iminus1 = w[0:N-1] 
    w_iminus1 = np.insert(w_iminus1,0 ,0)
    w_iplus1 = w[1:N]
    w_iplus1 = np.append(w_iplus1,0)
    #w_iminus1 = w_trans[0:N-1] #transpose of w could be the trans???
    #w_iminus1 = np.insert(w_iminus1,0 ,1)
    #w_iplus1 = w_trans[1:N] #transpose of w could be the trans???
    #w_iplus1 = np.append(w_iplus1,0)
    
    w_new = np.divide((np.multiply(w_old,(1-beta_zero))),(1+beta_zero)) - np.divide((2*dt*c0*m_atm*w),(1+beta_zero)) + np.divide((2*dt*s),(1+beta_zero)) + np.divide((np.multiply(beta_zero_x,(w_iplus1+w_iminus1))),(1+beta_zero))
    
    #Insert boundary conditions
    w_new[0] = np.divide((4*w_new[1] - w_new[2]),(3+((2*Abs*dx)/D))) #boundary condition at x=0
    w_new[N-1] = np.divide((4*w_new[N-2] - w_new[N-3]),(3+((2*Abs*dx)/D))) #boundary condition at x=xmax
 
    SDL = 20*np.log10(abs(w_new)) #sound density level
    
    #Update w before next step
    w_old = w #The w at n step becomes the w at n-1 step
    w = w_new #The w at n+1 step becomes the w at n step

    #w_rec is the energy density at the receiver specifically
    w_rec[steps] = w_new[index_dist_rec] #energy density at the receiver is equal to the energy density new calcuated in time
#    drawnow(draw_fig1)
#    drawnow(draw_fig2)

#Figure 3: Decay of SPL in the recording_time
plt.figure(3) 
SPL = 10*np.log10(((abs(w_rec))*rho*(c0**2))/(pRef**2))
plt.plot(t,SPL) #plot sound pressure level with Pref = (2e-5)**5
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
t60 = t60_impulse(t, sch_db, fspatial, rt='t30') #called function for calculation of t60 
