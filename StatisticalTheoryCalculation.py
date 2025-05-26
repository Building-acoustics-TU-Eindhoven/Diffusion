# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 09:38:17 2023

@author: 20225533
"""

import math
import matplotlib.pyplot as plt #import matplotlib as mpl
import numpy as np

from FunctionClarity import clarity
from FunctionDefinition import definition
from FunctionCentreTime import centretime

#Values calculated with the Barron revised theory
c0= 343 #sound particle velocity [m.s^-1]

#Room dimensions
lxmin = 0 #point x starts at zero [m]
lxmax = 8.0 #point x finish at the length of the room in the x direction [m] %Length
lymin = 0 #point y starts at zero [m]
lymax = 8.0 #point y finish at the length of the room in the y direction [m] %Width
lzmin = 0 #point z starts at zero [m]
lzmax = 8.0 #point z finish at the length of the room in the x direction [m] %Height

S1,S2 = lxmax*lymax, lxmax*lymax #xy planes
S3,S4 = lxmax*lzmax, lxmax*lzmax #xz planes
S5,S6 = lymax*lzmax, lymax*lzmax #yz planes

S = lxmax*lymax*2 + lxmax*lzmax*2 + lymax*lzmax*2 #Surface [m2]
V = lxmax*lymax*lzmax #volume of the room [m^3]

#Absorption terms
alpha_1 = 1/6 #Absorption coefficient for Surface1
alpha_2 = 1/6 #Absorption coefficient for Surface2
alpha_3 = 1/6 #Absorption coefficient for Surface3
alpha_4 = 1/6 #Absorption coefficient for Surface4
alpha_5 = 1/6 #Absorption coefficient for Surface5
alpha_6 = 1/6 #Absorption coefficient for Surface6

alpha_average = (alpha_1*S1 + alpha_2*S2 + alpha_3*S3 + alpha_4*S4 + alpha_5*S5 + alpha_6*S6)/S #average absorption
Eq_A = alpha_1*S1 + alpha_2*S2 + alpha_3*S3 + alpha_4*S4 + alpha_5*S5 + alpha_6*S6 #equivalent absorption area of the room

#Source Position
x_source = 4.0 #position of the source in the x direction [m]
y_source = 4.0 #position of the source in the y direction [m]
z_source = 4.0 #position of the source in the z direction [m]

#Receiver Position
x_rec = 6.0 #position of the receiver in the x direction [m]
y_rec = 6.0 #position of the receiver in the y direction [m]
z_rec = 6.0 #position of the receiver in the z direction [m]

dist = math.sqrt((abs(x_rec - x_source))**2 + (abs(y_rec - y_source))**2 + (abs(z_rec - z_source))**2) #distance between source and receiver

#Calculated t60 based on the Eyring formula
t60E = round((0.16*V)/(-S*np.log(1-alpha_average)),1)

#Calculated t60 based on the Sabine formula
t60S = (0.16*V)/Eq_A

d = 100/(dist**2) #direct field
er = ((31200*t60E)/V)*(math.exp(-(0.04*dist)/t60E))*(1-(math.exp(-(1.11/t60E)))) #early reflected
l = ((31200*t60E)/V)*(math.exp(-(0.04*dist)/t60E))*(math.exp(-(1.11/t60E))) #late reflected

c80dr = 10*math.log10((d+er)/l) #with direct and reverberant field
c80l = 10*math.log10((er)/l) #with direct field only

c80 = clarity(t60E, V, Eq_A, S, c0, dist) #called function for calculation of c80 [dB]
d50 = definition(t60E, V, Eq_A, S, c0, dist) #called function for calculation of d50 [%]
ts = centretime(t60E, Eq_A, S) #called function for calculation of ts [ms]

