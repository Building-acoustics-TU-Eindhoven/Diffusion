# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 10:39:42 2023

@author: 20225533
"""
#Code developed by Ilaria Fichera for the analysis of the FDM method Du Fort & Frankel solving the 3D diffusion equation with one intermittent omnidirectional sound source
import math
import matplotlib
import matplotlib.pyplot as plt #import matplotlib as mpl
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
import numpy as np
import time as time
from scipy import stats
from scipy.interpolate import griddata
from matplotlib.animation import FuncAnimation

st = time.time() #start time

#%%
#Input variables

# I advice to make a section of input variables, the things that the user can change.
# You may then decide to make that a separate .py file, like I did in PSTDbox.
# All operations that are working wit the input, or are fixed values (like the reference pressure value), I would keep in this file.

#General settings
# MH: you could also ask users to set a temperature, and then calculate the speed of sound based on in
c0= 343 #adiabatic speed of sound [m.s^-1]
m_atm = 0.1 #air absorption coefficient [1/m] from Billon 2008 paper and Navarro paper 2012


#Room dimensions
# MH: Your coordinate system starts in the corner, so basically users only need to enter the width, length and height of the room
lxmin = 0 #point x starts at zero [m]
lxmax = 30.0 #point x finish at the length of the room in the x direction [m] %Length
lymin = 0 #point y starts at zero [m]
lymax = 40.0 #point y finish at the length of the room in the y direction [m] %Width
lzmin = 0 #point z starts at zero [m]
lzmax = 4.0 #point z finish at the length of the room in the x direction [m] %Height

# Source position
x_source = 20.0  #position of the source in the x direction [m]
y_source = 1.0  #position of the source in the y direction [m]
z_source = 3.0  #position of the source in the z direction [m]

# Receier position
x_rec = 2.0 #position of the receiver in the x direction [m]
y_rec = 2.0 #position of the receiver in the y direction [m]
z_rec = 2.0 #position of the receiver in the z direction [m]

#Spatial discretization
dx = 0.5 #distance between grid points x direction [m]
dy = dx #distance between grid points y direction [m]
dz = dx #distance between grid points z direction [m]

#Time discretization
# MH: Do you have directions for the user on how to chose this dt value?
dt = 1/8000 #distance between grid points on the time discretization [s]
recording_time = 1.00 #time recorded for the source [s]

th = 3 #int(input("Enter type Absortion conditions (option 1,2,3):")) 
# options Sabine (th=1), Eyring (th=2) and modified by Xiang (th=3)
alpha_1 = 0.17 #Absorption coefficient for Surface1 - Floor
alpha_2 = 0.17 #Absorption coefficient for Surface2 - Ceiling
alpha_3 = 0.17 #Absorption coefficient for Surface3 - Wall Front
alpha_4 = 0.17 #Absorption coefficient for Surface4 - Wall Back
alpha_5 = 0.17 #Absorption coefficient for Surface5 - Wall Left
alpha_6 = 0.17 #Absorption coefficient for Surface6 - Wall Right

#Set initial condition - Source Info (interrupted method)
Ws=0.005 #Source point power [Watts] interrupted after "sourceon_time" seconds; 10^-2 W => correspondent to 100dB
Vs=0.2  #MH WHAT IS Vs?
sourceon_time =  0.1 #time that the source is on before interrupting [s]

## MH: here the fixed input section starts
#Frequency resolution & spatial parameters
fsample = 1/dt #frequency spatial resolution (sampling period)

pRef = 2 * (10**-5) #Reference pressure in Pa
rho = 1.21 #air density [kg.m^-3] at 20°C
t = np.arange(0, recording_time, dt) #mesh point in time
recording_steps = ceil(recording_time/dt) #number of time steps to consider in the calculation

S1,S2 = lxmax*lymax, lxmax*lymax #xy planes
S3,S4 = lxmax*lzmax, lxmax*lzmax #xz planes
S5,S6 = lymax*lzmax, lymax*lzmax #yz planes

S = lxmax*lymax*2 + lxmax*lzmax*2 + lymax*lzmax*2 #Total Surface Area[m2]
V = lxmax*lymax*lzmax #Volume of the room [m^3]

x = np.arange(lxmin, lxmax+dx, dx) #mesh points in space x direction
y = np.arange(lymin, lymax+dy, dy) #mesh points in space y direction
z = np.arange(lzmin, lzmax+dz, dz) #mesh points in space z direction

Nx = len(x) #number of point in the x direction
Ny = len(y) #number of point in the y direction
Nz = len(z) #number of point in the z direction

yy, xx , zz = np.meshgrid(y,x,z) #Return coordinate matrices from coordinate vectors; create the 3D grid

# MH: uncoment this when using drawnow
##Visualization of 3D meshgrid
##fig = plt.figure()
##ax = plt.axes(projection='3d')
##ax = ax.scatter3D(xx, yy, zz, c=zz, cmap='Greens')
##plt.title("Visualization of 3D meshgrid")

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
    print("aa! error! Check beta condition")

#Set initial condition - Source Info (interrupted method)
w1=Ws #w1 = round(Ws/Vs,4) #power density of the source [Watts/(m^3))]
sourceon_steps = ceil(sourceon_time/dt) #time steps at which the source is calculated/considered in the calculation
s1 = np.multiply(w1,np.ones(sourceon_steps)) #energy density of source number 1 at each time step position
source1 = np.append(s1, np.zeros(recording_steps-sourceon_steps)) #This would be equal to s1 if and only if recoding_steps = sourceon_steps

#############################################################################
#SOURCE INTERPOLATION
#############################################################################
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

#############################################################################
#RECEIVER INTERPOLATION
#############################################################################
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

dist_sr = math.sqrt((abs(x_rec - x_source))**2 + (abs(y_rec - y_source))**2 + (abs(z_rec - z_source))**2) #distance between source and receiver

w_new = np.zeros((Nx,Ny,Nz)) #unknown w at new time level (n+1)
w = w_new #w at n level
w_old = w #w_old at n-1 level

w_rec = np.arange(0,recording_time,dt) #energy density at the receiver

#Function to draw figure
def draw_fig():
    for i in range(w_new.shape[2]):
        print(i)
        plt.imshow(w_new[:, :, i], cmap='hot', vmin=w_new.min(), vmax=w_new.max())
        plt.show()  # Display each slice separately

############################################################################
#TRIAL 4D PLOT - INITIALIZATION (main part within the for loop)
############################################################################
#Initialize the figure and axes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Set the initial view angle
ax.view_init(elev=30, azim=45)

############################################################################
#MAIN CALCULATION - COMPUTING ENERGY DENSITY
############################################################################ 
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
    w_new[-1,:,:] = np.divide((4*w_new[-2,:,:] - w_new[-3,:,:]),(3+((2*Abs_6*dx)/Dx))) #boundary condition at lx=lxmax, any y, any z


    w_new[:,0,:] = np.divide((4*w_new[:,1,:] - w_new[:,2,:]),(3+((2*Abs_3*dx)/Dy))) #boundary condition at y=0, any x, any z
    w_new[:,-1,:] = np.divide((4*w_new[:,-2,:] - w_new[:,-3,:]),(3+((2*Abs_4*dx)/Dy))) #boundary condition at at ly=lymax, any x, any z
 
    
    w_new[:,:,0] = np.divide((4*w_new[:,:,1] - w_new[:,:,2]),(3+((2*Abs_1*dx)/Dz))) #boundary condition at z=0, any x, any y
    w_new[:,:,-1] = np.divide((4*w_new[:,:,-2] - w_new[:,:,-3]),(3+((2*Abs_2*dx)/Dz))) #boundary condition at at lz=lzmax, any x, any y
    
    sdl = 10*np.log10(abs(w_new),where=abs(w_new)>0) #sound density level
    
 # MH uncomment when activating the drawnow library   
    ##Visualization of the energy density changes while the calculation is progressing
   # if (steps % 100 == 0): #draw only on certain steps and not all the steps
   #     print("A")
   #     drawnow(draw_fig)
    
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
        print("Steps for source:",steps)
        w_t0 = w_new
    
    # MH if the commented secionts are not needed, remove them from this script. 
    
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
    ##ax.set_title('Time Step: ', time_steps)
    ##ax.set_xlabel('X')
    ##ax.set_ylabel('Y')
    ##ax.set_zlabel('Z')

    #Adjust the plot limits
    ##ax.set_xlim(lxmin, lxmax)
    ##ax.set_ylim(lymin, lymax)
    ##ax.set_zlim(lzmin, lzmax)

    # Add a colorbar and legend
    #cbar = fig.colorbar(scatter, ax=ax,fraction=0.04, pad=0.1)
    #cbar.set_label('Energy Density')
    #cbar.ax.yaxis.set_ticks_position('left')

    # Pause to create an animated effect
    ##plt.pause(0.01)  # Adjust the pause duration as needed
    
    #Updating the source term
    s[row_lower, col_lower, depth_lower] = source1[steps] * weight_row_lower * weight_col_lower * weight_depth_lower
    s[row_lower, col_lower, depth_upper] = source1[steps] * weight_row_lower * weight_col_lower * weight_depth_upper
    s[row_lower, col_upper, depth_lower] = source1[steps] * weight_row_lower * weight_col_upper * weight_depth_lower
    s[row_lower, col_upper, depth_upper] = source1[steps] * weight_row_lower * weight_col_upper * weight_depth_upper
    s[row_upper, col_lower, depth_lower] = source1[steps] * weight_row_upper * weight_col_lower * weight_depth_lower
    s[row_upper, col_lower, depth_upper] = source1[steps] * weight_row_upper * weight_col_lower * weight_depth_upper
    s[row_upper, col_upper, depth_lower] = source1[steps] * weight_row_upper * weight_col_upper * weight_depth_lower
    s[row_upper, col_upper, depth_upper] = source1[steps] * weight_row_upper * weight_col_upper * weight_depth_upper
     
    print(time_steps)

plt.show()

press_r = ((abs(w_rec))*rho*(c0**2)) #pressure
spl = 10*np.log10(((abs(w_rec))*rho*(c0**2))/(pRef**2)) #,where=press_r>0
spl_tot = 10*np.log10(rho*c0*((Ws/(4*math.pi*dist_sr**2))*np.exp(-m_atm*dist_sr) + ((abs(w_rec))*c0)/(pRef**2))) #spl total (including direct field) at the receiver position????? but it will need to be calculated for astationary source 100dB
spl_norm = 10*np.log10((((abs(w_rec))*rho*(c0**2))/(pRef**2)) / np.max(((abs(w_rec))*rho*(c0**2))/(pRef**2))) #normalised to maximum to 0dB

#Schroeder integration
idx_w_rec = np.where(t == sourceon_time)[0][0] #index at which the t array is equal to the sourceon_time; I want the RT to calculate from when the source stops.
w_rec_off = w_rec[idx_w_rec:] #cutting the energy density array at the receiver from the idx_w_rec to the end   
energy_r_rev = (w_rec_off)[::-1] #reverting the array

#The energy density is related to the pressure with the following relation: w = p^2
energy_r_rev_cum = np.cumsum(energy_r_rev) #cumulative summation of all the item in the array
schroeder = energy_r_rev_cum[::-1] #reverting the array again -> creating the schroder decay
sch_db = 10.0 * np.log10(schroeder / np.max(schroeder)) #level of the array: schroeder decay

t60 = t60_decay(t, sch_db, idx_w_rec) #called function for calculation of t60 [s]
edt = edt_decay(t, sch_db, idx_w_rec) #called function for calculation of edt [s]
c80 = clarity(t60, V, Eq_A, S, c0, dist_sr) #called function for calculation of c80 [dB]
d50 = definition(t60, V, Eq_A, S, c0, dist_sr) #called function for calculation of d50 [%]
ts = centretime(t60, Eq_A, S) #called function for calculation of ts [ms]

et = time.time() #end time
elapsed_time = et - st

#%%
###Figures###

#Figure 3: Decay of SPL in the recording_time
plt.figure(3) 
plt.plot(t,spl) #plot sound pressure level with Pref = (2e-5)**5
plt.title("SPL over time at the receiver")
plt.xlabel("t [s]")
plt.ylabel("SPL [dB]")
plt.xlim()
plt.ylim()
plt.xticks(np.arange(0, recording_time +0.1, 0.5))
plt.yticks(np.arange(0, 120, 20))

#Figure 4: Decay of SPL in the recording_time normalised to maximum 0dB
plt.figure(4)
plt.plot(t,spl_norm)
plt.title("Normalised SPL over time at the receiver")
plt.xlabel("t [s]")
plt.ylabel("SPL [dB]")
plt.xlim()
plt.ylim()
plt.xticks(np.arange(0, recording_time +0.1, 0.1))
plt.yticks(np.arange(0, -60, -10))

#Figure 5: Energy density at the receiver over time
plt.figure(5)
plt.plot(t,w_rec)
plt.title("Energy density over time at the receiver")
plt.xlabel("t [s]")
plt.ylabel("Energy density [kg m^-1 s^-2]")
plt.xlim()
plt.ylim()
plt.xticks(np.arange(0, recording_time +0.1, 0.1))

#Figure 6: Schroeder decay
plt.figure(6)
plt.plot(t[idx_w_rec:],sch_db)
plt.title("Schroeder decay")
plt.xlabel("t [s]")
plt.ylabel("Energy decay [dB]")
plt.xlim()
plt.ylim()
plt.xticks(np.arange(t[idx_w_rec], recording_time +0.1, 0.1))

# MH please write a caption of this figure to help the reader what they see. I suggest to make a 2D colorpolot here, 
# I think that it gives a better picture on the distribution of the energy over space. Also, at what time is w shown?
#Figure 7: 3D image of the energy density in the room
plt.figure(7)
fig = plt.figure(7)
ax = fig.add_subplot(111, projection='3d')
ax.set_box_aspect([0.5, 0.5, 0.5])  # Set the aspect ratio of the plot
X, Y= np.meshgrid(x, y)
X = X.T
Y = Y.T
ax.contourf(X, Y, w_new[:, :, 4], cmap='hot', alpha=0.8)
ax.view_init(azim=-120, elev=30)  # Set the viewing angle
plt.show()

##############################################################################
#TRIAL CALCULATION PICAUT 1997
##############################################################################

sigma = c0*alpha_1/lambda_path
mean_w_t0 = np.mean(w_t0)
summation=np.sum(w_t0)
#w_inf = np.sum(w_t0)*(dx*dy*dz)/V
w_inf = np.sum(w_t0)/((Nx-1)*(Ny-1)*(Nz-1))
Q = abs((w_rec/w_inf)*np.exp(sigma*t) - 1)
phi = np.log10(Q)
A_x = np.log10(abs((2*np.cos(np.pi*x/lxmax)))) 

#phi_norm = phi / np.max(phi)
two_lambda_time = 2*lambda_path/c0

difference = np.abs(t-(sourceon_time+two_lambda_time))
index_diff = difference.argmin()

idx_phi = np.where(t == t[index_diff])[0][0] 
phi_off = phi[idx_phi:]  

# Linregress approach
slope,intercept = stats.linregress(t[idx_phi:],phi_off)[0:2] 

Diff = ((slope*lxmax**2)/(np.pi**2))

# Poly-based Approach y = Ax + B
CoefAlpha = np.polyfit(t[idx_phi:], phi_off, 1) ##calculating the slope and the interception of the line connecting the two points
slope2 = CoefAlpha[0]

#Diffusion = (np.subtract(A_x[70], phi)*lxmax**2)/(np.pi**2*t)
#Aminusphi =np.subtract(A_x[70], phi)

D_3D = ((6*math.log(10)/t60) * ((4*lxmax*lzmax)/(np.pi**2*(lxmax+4*lzmax))))*lxmax #not applicable for rooms; only for streets
D_2D = ((6*np.log10(10)/t60) * ((lxmax)/(np.pi**2))) #not applicable for rooms; only for streets

plt.figure(9)
plt.plot(x/lxmax,A_x)
plt.title("Representation of A{x}")
plt.xlabel("x/L")
plt.ylabel("A{x}")
plt.ylim([-4, 1])
plt.yticks(np.arange(-4, 1, 1))

#plt.figure(10)
#plt.plot(t,Diffusion)
#plt.title("Diffusion coefficient over time")
#plt.xlabel("t [s]")
#plt.ylabel("Diffusion coefficient [m^2/s]")
#plt.ylim([0,20000])
#plt.yticks(np.arange(-4, 1, 1))
