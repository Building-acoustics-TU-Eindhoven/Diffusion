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
#from drawnow import drawnow
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

#Time discretization
dt = 1/32000 #distance between grid points on the time discretization [s]
recording_time = 4 #time recorded for the source [s]
recording_steps = ceil(recording_time/dt) #number of time steps to consider in the calculation
t = np.arange(0, recording_time, dt) #mesh point in time

#Frequency resolution & spatial parameters
fsample = 1/dt #frequency spatial resolution (sampling period)

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

th = 2 #int(input("Enter type Asbortion conditions (option 1,2,3):")) #input 1,2,3 just to understand the type of boundary chosen
alpha_1 = 0.2 #Absorption coefficient for Surface1
alpha_2 = 0.2 #Absorption coefficient for Surface2
alpha_3 = 0.2 #Absorption coefficient for Surface3
alpha_4 = 0.5 #Absorption coefficient for Surface4
alpha_5 = 0.2 #Absorption coefficient for Surface5
alpha_6 = 0.2 #Absorption coefficient for Surface6

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

beta_zero_x = (2*Dx*dt)/(dx**2) #mesh number in x direction
beta_zero_y = (2*Dy*dt)/(dy**2) #mesh number in x direction
beta_zero_z = (2*Dz*dt)/(dz**2) #mesh number in x direction
beta_zero = beta_zero_x + beta_zero_y + beta_zero_z #beta_zero is the condition for all the directions deltax, deltay and deltaz.
 
#Condition for the model to be unconditionally stable
beta_zero_condition = ((beta_zero**2) - 1)/(1+(beta_zero**2)+(2*beta_zero)) #from Navarro 2012 paper
if beta_zero_condition >1:
    print("aa! errors! Check beta condition")

#Set initial condition - Source Info (interrupted method)
Ws=0.01 #Source point power [Watts] interrupted after 2seconds; 10^-2 value taken from Jing 2007; correspondent to a SWL of 100dB
Vs=0.2
w1=Ws
#w1 = round(Ws/Vs,4) #power density of the source [Watts/(m^3))]
sourceon_time =  2 #time that the source is on before interrupting [s]
sourceon_steps = ceil(sourceon_time/dt) #time steps at which the source is calculated/considered in the calculation
s1 = np.multiply(w1,np.ones(sourceon_steps)) #energy density of source number 1 at each time step position #does the source not need to be only at the time 0 to 2seconds and after that there should not be any source term? Yes
source1 = np.append(s1, np.zeros(recording_steps-sourceon_steps)) #This would be equal to s1 if and only if recoding_steps = sourceon_steps

#Finding index in meshgrid of the source position
x_source = 2.5 #int(ceil(Nx/2))#4 #position of the source in the x direction [m]
y_source = 2.5 #int(ceil(Ny/2))#4 #position of the source in the y direction [m]
z_source = 2.5 #int(ceil(Nz/2))#4 #position of the source in the z direction [m]
coord_source = [x_source , y_source, z_source] #coordinates of the source position in an list
rows_s = np.argmin(abs(xx[:,0,0] - coord_source[0])) #Find index of grid point with minimum distance from source along x direction
cols_s = np.argmin(abs(yy[0,:,0] - coord_source[1])) #Find index of grid point with minimum distance from source along y direction
dept_s = np.argmin(abs(zz[0,0,:] - coord_source[2])) #Find index of grid point with minimum distance from source along z direction

#Finding index in meshgrid of the receiver position
x_rec = 1.0 #int(ceil(Nx/4)) #position of the receiver in the x direction [m]
y_rec = 1.0 #int(ceil(Nx/4)) #position of the receiver in the y direction [m]
z_rec = 1.0 #int(ceil(Nx/4)) #position of the receiver in the z direction [m]
coord_receiver = [x_rec,y_rec,z_rec] #coordinates of the receiver position in an list
rows_r = np.argmin(abs(xx[:,0,0] - coord_receiver[0])) #Find index of grid point with minimum distance from receiver along x direction
cols_r = np.argmin(abs(yy[0,:,0] - coord_receiver[1])) #Find index of grid point with minimum distance from receiver along y direction
dept_r = np.argmin(abs(zz[0,0,:] - coord_receiver[2])) #Find index of grid point with minimum distance from receiver along z direction

dist_sr = math.sqrt((abs(x_rec - x_source))**2 + (abs(y_rec - y_source))**2 + (abs(z_rec - z_source))**2) #distance between source and receiver

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

#Computing w;
for steps in range(0, recording_steps):
    #Compute w at inner mesh points
    time_steps = steps*dt #total time for the calculation
    #s[rows_s, cols_s, dept_s] = source1[steps] #array of zero of the source apart from the index_dist_source = energy density of the source at each step position
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
    w_new[0,:,:] = np.divide((4*w_new[1,:,:] - w_new[2,:,:]),(3+((2*Abs_5*dx)/Dx))) #boundary condition at x=0, any y, any z
    w_new[-1,:,:] = np.divide((4*w_new[-2,:,:] - w_new[-3,:,:]),(3+((2*Abs_6*dx)/Dx))) #boundary condition at lx=lxmax, any y, any z


    w_new[:,0,:] = np.divide((4*w_new[:,1,:] - w_new[:,2,:]),(3+((2*Abs_3*dx)/Dy))) #boundary condition at y=0, any x, any z
    w_new[:,-1,:] = np.divide((4*w_new[:,-2,:] - w_new[:,-3,:]),(3+((2*Abs_4*dx)/Dy))) #boundary condition at at ly=lymax, any x, any z
 
    
    w_new[:,:,0] = np.divide((4*w_new[:,:,1] - w_new[:,:,2]),(3+((2*Abs_1*dx)/Dz))) #boundary condition at z=0, any x, any y
    w_new[:,:,-1] = np.divide((4*w_new[:,:,-2] - w_new[:,:,-3]),(3+((2*Abs_2*dx)/Dz))) #boundary condition at at lz=lzmax, any x, any y
    
    sdl = 10*np.log10(abs(w_new),where=abs(w_new)>0) #sound density level
    #if (steps % 100 == 0): #draw only on certain steps and not all the steps
    #    drawnow(draw_fig)
    #Update w before next step
    w_old = w #The w at n step becomes the w at n-1 step
    w = w_new #The w at n+1 step becomes the w at n step

    #w_rec is the energy density at the receiver specifically
    w_rec[steps] = w_new[rows_r, cols_r, dept_r] #energy density at the receiver is equal to the energy density new calcuated in time
    
    s[rows_s, cols_s, dept_s] = source1[steps] #array of zero of the source apart from the index_dist_source = energy density of the source at each step position
    
    print(time_steps)
    #drawnow(draw_fig1)
    #drawnow(draw_fig2)

spl_tot = 10*np.log10(rho*c0*((Ws/(4*math.pi*dist_sr**2))*np.exp(-m_atm*dist_sr) + ((abs(w_rec))*c0)/(pRef**2))) #It should be the spl total (including direct field) at the receiver position????? but it will need to be calculated for astationary source 100dB

#Figure 3: Decay of SPL in the recording_time
plt.figure(3) 
press_r = ((abs(w_rec))*rho*(c0**2))
#max_press_r = np.max(((abs(w_rec))*rho*(c0**2))/(pRef**2))
spl = 10*np.log10(((abs(w_rec))*rho*(c0**2))/(pRef**2)) #,where=press_r>0
plt.plot(t,spl) #plot sound pressure level with Pref = (2e-5)**5
plt.title("SPL over time at the receiver")
plt.xlabel("t")
plt.ylabel("SPL")
plt.xlim()
plt.ylim()
plt.xticks(np.arange(0, recording_time +0.1, 0.5))
#plt.yticks(np.arange(0, 120, 20))

spl_norm = 10*np.log10((((abs(w_rec))*rho*(c0**2))/(pRef**2)) / np.max(((abs(w_rec))*rho*(c0**2))/(pRef**2))) #normalised to maximum to 0dB

#Figure 4: Decay of SPL in the recording_time normalised to maximum 0dB
plt.figure(4)
plt.plot(t,spl_norm)
plt.title("Normalised SPL over time at the receiver")
plt.xlabel("t")
plt.ylabel("SPL")
plt.xlim()
plt.ylim()
plt.xticks(np.arange(0, recording_time +0.1, 0.1))
plt.yticks(np.arange(0, -60, -10))

plt.figure(5)
plt.title("Energy density over time at the receiver")
plt.plot(t,w_rec)

#Failed trial for graph of sound pressure level stationary over the space.
#plt.figure(6)
#t_dim = len(t)
#last_time_index = t_dim-1
#spl_y = spl_tot[index_receiver[0],:,index_receiver[2]]
#data_y = spl_y
#plt.plot(y,data_y)
#plt.title("SPL over the y axis")

#plt.figure(7)
#plt.plot(t,spl_tot) #plot sound pressure level with Pref = (2e-5)**5
#plt.title("SPL_tot over time at the receiver")
#plt.xlabel("t")
#plt.ylabel("SPL")
#plt.xlim()
#plt.ylim()
#plt.xticks(np.arange(0, recording_time +0.1, 0.5))

#3D image of the energy density in the room
#fig = plt.figure(1)
#ax = fig.add_subplot(111, projection='3d')
#ax.set_box_aspect([0.5, 0.5, 0.5])  # Set the aspect ratio of the plot
#X, Y= np.meshgrid(x, y)
#ax.contourf(X, Y, w_new[:, :, 4], cmap='hot', alpha=0.8)
#ax.view_init(azim=-120, elev=30)  # Set the viewing angle
#plt.show()

init = -5.0 #because I want the T30, I need to start at -5
end = -35.0 #because I want the T30, I need to finish at -35
factor = 2.0 #factor of 2 since I need the T30

#Schroeder integration
idx_w_rec = np.where(t == sourceon_time)[0][0] #index at which the t array is equal to the sourceon_time; I want the RT to calculate from when the source stops.
w_rec = w_rec[idx_w_rec:] #cutting the energy density array at the receiver from the idx_w_rec to the end   
press_r_rev = (w_rec)[::-1] #reverting the array
#The energy density is related to the pressure with the following relation: w = p^2
press_r_rev_cum = np.cumsum(press_r_rev) #sumulative summation of all the item in the array
schroeder = press_r_rev_cum[::-1] #reverting the array again -> creating the schroder decay
sch_db = 10.0 * np.log10(schroeder / np.max(schroeder)) #level of the array: schroeder decay
#plt.plot(t,sch_db)
    
#Linear regression
rt_decay = sch_db #decay to be used to calcuate the RT
 
timeVector = np.arange(0, recording_time, dt) #equal to t
idxL1 = np.where(rt_decay <= init)[0][0] #index at which the rtdecay is equal to -5
idxL2 = np.where(rt_decay <= end)[0][0] #index at which the rtdecay is equal to -35
   
timeL1 = timeVector[idxL1] #index at which the time vector is equal to the idxL1
timeL2 = timeVector[idxL2] #index at which the time vector is equal to the idxL2
EDCTimeVec = np.arange(0, (len(rt_decay)-1)/fsample+1/fsample, 1/fsample) #equal to t and timeVector

# Classical Approach (T30 = 2*(t[-35dB]-t[-5dB])
RTCalc = factor*(timeL2 - timeL1)

# Linregress approach
slope,intercept = stats.linregress(EDCTimeVec[idxL1:idxL2],rt_decay[idxL1:idxL2])[0:2] #calculating the slope and the interception of the line connecting the two points
db_regress_init = (init - intercept) / slope #dB initial
db_regress_end = (end - intercept) / slope #dB End
t60I = factor * (db_regress_end - db_regress_init) #t60 according to linregress approach

# Poly-based Approach y = Ax + B
CoefAlpha = np.polyfit(EDCTimeVec[idxL1:idxL2], rt_decay[idxL1:idxL2], 1) ##calculating the slope and the interception of the line connecting the two points
t60 = (-60/CoefAlpha[0]) #t60 according to polyfit approach
   
EDCL1 = rt_decay[idxL1]
EDCL2 = rt_decay[idxL2]
    
y_axis = (slope*timeVector[idx_w_rec:] + intercept) + slope

plt.figure()
plt.plot(timeVector[idx_w_rec:],rt_decay, color ='b', linewidth = 1.8)
plt.plot(timeVector[idx_w_rec:],y_axis,color='r',linewidth=2)
plt.plot(timeVector[idx_w_rec:][idxL1],np.real(rt_decay[idxL1]),'o',linewidth=2)
plt.plot(timeVector[idx_w_rec:][idxL2],np.real(rt_decay[idxL2]),'o',linewidth=2)
plt.axvline(x=t60I,ymin=-100,ymax=0,linestyle='--',linewidth=2)
plt.ylabel('Normalized Magnitude (dB)')
plt.xlabel('Time (s)')
plt.legend(['EDC','Line Fitting','Upper Point','Lower Point','Estimate T_{30}'])
plt.title('T_{30} = ' + str(round(t60I,2)) + ' s.')
plt.grid(True)
plt.ylim([-100,0])
plt.show()

#t60 = t60_decay(t, press_r, fsample, rt='t30') #called function for calculation of t60 [s]
#t60M = t60_decayM(t, spl_norm, fsample, rt='t30') #called function for calculation of t60 [s]
#edt = t60_decay(t, spl_norm, fsample, rt='edt') #called function for calculation of edt [s]
c80 = clarity(t60, V, Eq_A, S, c0, dist_sr) #called function for calculation of c80 [dB]
d50 = definition(t60, V, Eq_A, S, c0, dist_sr) #called function for calculation of d50 [%]
ts = centretime(t60, Eq_A, S) #called function for calculation of ts [ms]

#EDT Calc
initEDT = 0.0 #because I want the T30, I need to start at -5
endEDT = -10.0 #because I want the T30, I need to finish at -35
factorEDT = 6.0 #factor of 2 since I need the T30

idxL1EDT = np.where(rt_decay <= initEDT)[0][0] #index at which the rtdecay is equal to -5
idxL2EDT = np.where(rt_decay <= endEDT)[0][0] #index at which the rtdecay is equal to -35
   
timeL1EDT = timeVector[idxL1EDT] #index at which the time vector is equal to the idxL1
timeL2EDT = timeVector[idxL2EDT] #index at which the time vector is equal to the idxL2
EDCTimeVec = np.arange(0, (len(rt_decay)-1)/fsample+1/fsample, 1/fsample) #equal to t and timeVector

# Classical Approach (T30 = 2*(t[-35dB]-t[-5dB])
EDTCalc = factorEDT*(timeL2EDT - timeL1EDT)

# Linregress approach
slopeEDT,interceptEDT = stats.linregress(EDCTimeVec[idxL1EDT:idxL2EDT],rt_decay[idxL1EDT:idxL2EDT])[0:2] #calculating the slope and the interception of the line connecting the two points
db_regress_initEDT = (initEDT - interceptEDT) / slopeEDT #dB initial
db_regress_endEDT = (endEDT - interceptEDT) / slopeEDT #dB End
EDTI = factorEDT * (db_regress_endEDT - db_regress_initEDT) #t60 according to linregress approach

# Poly-based Approach y = Ax + B
CoefAlpha = np.polyfit(EDCTimeVec[idxL1EDT:idxL2EDT], rt_decay[idxL1EDT:idxL2EDT], 1) ##calculating the slope and the interception of the line connecting the two points
EDT = (-60/CoefAlpha[0]) #t60 according to polyfit approach




et = time.time() #end time
elapsed_time = et - st