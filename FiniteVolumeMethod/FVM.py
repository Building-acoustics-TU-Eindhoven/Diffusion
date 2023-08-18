# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 16:12:40 2023

@author: 20225533
"""

#Code developed by Ilaria Fichera for the analysis of the FV method with Du Fort & Frankel method adapted solving the 3D diffusion equation with one intermittent omnidirectional sound source
#Import modules
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

import gmsh
import sys
 
st = time.time() #start time of calculation

#%%
###############################################################################
#INPUT VARIABLES
###############################################################################

#General settings
c0= 343 #adiabatic speed of sound [m.s^-1]
m_atm = 0 #air absorption coefficient [1/m] from Billon 2008 paper and Navarro paper 2012

length_mesh = 2 ###???We do not need it!
dt = np.sqrt(length_mesh*(10**-8)/2) #time discretizatione

# Source position
x_source = 0.5  #position of the source in the x direction [m]
y_source = 0.5  #position of the source in the y direction [m]
z_source = 1.0  #position of the source in the z direction [m]

# Receiver position
x_rec = 2.0 #position of the receiver in the x direction [m]
y_rec = 0.5 #position of the receiver in the y direction [m]
z_rec = 1.0 #position of the receiver in the z direction [m]

#Absorption term and Absorption coefficients
th = 3 #int(input("Enter type Absortion conditions (option 1,2,3):")) 
# options Sabine (th=1), Eyring (th=2) and modified by Xiang (th=3)
alpha_1 = 1/6 #Absorption coefficient for Surface1 - Floor
alpha_2 = 1/6 #Absorption coefficient for Surface2 - Ceiling
alpha_3 = 1/6 #Absorption coefficient for Surface3 - Wall Front
alpha_4 = 1/6 #Absorption coefficient for Surface4 - Wall Back
alpha_5 = 1/6 #Absorption coefficient for Surface5 - Wall Left
alpha_6 = 1/6 #Absorption coefficient for Surface6 - Wall Right

#%%
###############################################################################
#INITIALISE GMSH
###############################################################################
    
file_name = "rect.msh" #Insert file name, msh file created from sketchUp and then gmsh
gmsh.initialize() #Initialize msh file
mesh = gmsh.open(file_name) #open the file

#gmsh.fltk.run() #run the file to see it in gmsh

dim = -1 #dimensions of the entities, 0 for points, 1 for curves/edge/lines, 2 for surfaces, 3 for volumes, -1 for all the entities 
tag = -1 #all the nodes of the room

#Nodes
nodeTags, coords, parametricCoord = gmsh.model.mesh.getNodes(dim,tag) #gets the tags for each node and the coordinates of each node
nodecoords = coords.reshape((-1,3)) #coordinates reshaped in a matrix 3xnumber of nodes


#Element Types
elemTypes,elemTags,elemNodeTags = gmsh.model.mesh.getElements(dim,tag)
#elemTypes = 1 for lines, 2 for surfaces, 4 for tetrahedron
#elemTags =  list of list of lines, boundary elements (surfaces) and volume elements (tetrahedron)
#elemNodeTags = did not understand this yet, probably a tag gives to each line, surface and tetrahedron
for e_type in elemTypes:
    if e_type == 1: #if the e_type = 1, then get all the elements (lines) with that e_type
        tag = -1 
        edgeEl, edgeNodeTagstype = gmsh.model.mesh.getElementsByType(e_type,tag)
        #edgeEl = numbered edge elements (lines)
    elif e_type == 2: #if the e_type = 2, then get all the elements (surfaces) with that e_type
        tag = -1 
        bounEl, bounNodeTagstype = gmsh.model.mesh.getElementsByType(e_type,tag)
        bounNode = bounNodeTagstype.reshape((-1,3))
        #boundEl = numbered boundary elements (surfaces)
        #bounNode = nodes of the surfaces
    elif e_type == 4:
        tag = -1#if the e_type = 4, then get all the elements (tetrahedron) with that e_type
        voluEl, voluNodeTagstype = gmsh.model.mesh.getElementsByType(e_type,tag)
        voluNode = voluNodeTagstype.reshape((-1,4))
        #voluEl = numbered volume elements (tetrahedron)
        #voluNode = nodes of the tetrahedon
    else:
        print("The procedure is not possible")

        
velemNodes = elemNodeTags[2].reshape((-1,4)) #nodes per each tetrahedron 
belemNodes = elemNodeTags[1].reshape((-1,3)) #nodes per each surface boundary 

eelement = 0
for elem in range(len(edgeEl)):
    eelement = eelement + 1 #scalar number of the edge elements

belement = 0
for elem in range(len(bounEl)):
    belement = belement + 1 #scalar number of the boundary elements

velement = 0
for elem in range(len(voluEl)):
    velement = velement + 1 #scalar number of the volume elements


#Volume Element dictionary + nodes of each volume elements (4 nodes per element)
volumeEl_dict = {} #Dictionary of volumelements + its nodes
for i in range(len(voluEl)):
    #print(i)
    volumeEl_dict[voluEl[i]] = velemNodes[i]
  
#Each element has 4 nodes, each node has three coordinates; here we store the coordinates per each node    
vnodeCoord_dict = {}
for i in range(len(elemNodeTags[2])):
    if not elemNodeTags[2][i] in vnodeCoord_dict.keys(): 
        vnodeCoord_dict[elemNodeTags[2][i]] = nodecoords[int(elemNodeTags[2][i])-1]
    #print(nodes)
    
#Boundary Element dictionary + node per each surface elements (3 nodes per element)
boundaryEl_dict = {} #Dictionary of boundary elements + its nodes
for i in range(len(bounEl)):
    boundaryEl_dict[bounEl[i]] = belemNodes[i]    

bnodeCoord_dict = {}
for i in range(len(elemNodeTags[1])):
    if not elemNodeTags[1][i] in bnodeCoord_dict.keys(): 
        bnodeCoord_dict[elemNodeTags[1][i]] = nodecoords[int(elemNodeTags[1][i])-1]
    #print(nodes)

#volumeEl_dict = {} #Dictionary of volumelements + its nodes
#for i in range(velement):
#    volumeEl_dict[i] = velemNodes[i]
    
#boundaryEl_dict = {} #Dictionary of boundary elements + its nodes
#for i in range(belement):
#    boundaryEl_dict[i] = belemNodes[i]


#Calculation of volume cells and centre of volume    
vcell_dict = {} #volume of each element tetrahedron initialization
centre_cell = {} #centre of the element tetrahedron initialization
for i in volumeEl_dict.keys():
    coord_centre_cell = np.zeros(3)
    centre_cell[i] = []
    #print(i)
    vc0 = vnodeCoord_dict[volumeEl_dict[i][0]] #Coordinates of the node number zero of the volume element i
    #print(nc0)
    vc1 = vnodeCoord_dict[volumeEl_dict[i][1]] #Coordinates of the node number one of the volume element i
    #print(nc1)
    vc2 = vnodeCoord_dict[volumeEl_dict[i][2]] #Coordinates of the node number two of the volume element i
    #print(nc2)
    vc3 = vnodeCoord_dict[volumeEl_dict[i][3]] #Coordinates of the node number three of the volume element i
    #print(nc3)
    for j in range(3):
        coord_centre_cell[j] = (vc0[j]+vc1[j]+vc2[j]+vc3[j])/4
        centre_cell[i].append(coord_centre_cell[j])
        #centre_cell[i,j] = (vc0[j]+vc1[j]+vc2[j]+vc3[j])/4; #coordinates of the centre of each volume element
    vcell_dict[i] = abs(np.dot(np.cross(vc1-vc3,vc2-vc3),vc0-vc3))/6 #volume of each volume element


#Calculation of boundary elements area and centre   
barea_dict = {} #surface of each element boundary initialization
centre_area = {} #centre of the element tetrahedron initialization
for i in boundaryEl_dict.keys():
    coord_centre_area = np.zeros(3)
    centre_area[i] = []
    #print(i)
    bc0 = bnodeCoord_dict[boundaryEl_dict[i][0]] #Coordinates of the node number zero of the volume element i
    #print(nc0)
    bc1 = bnodeCoord_dict[boundaryEl_dict[i][1]] #Coordinates of the node number one of the volume element i
    #print(nc1)
    bc2 = bnodeCoord_dict[boundaryEl_dict[i][2]] #Coordinates of the node number two of the volume element i
    #print(nc2)
    for j in range(3):
        coord_centre_area[j] = (bc0[j]+bc1[j]+bc2[j])/3 #coordinates of the centre of each volume element
        centre_area[i].append(coord_centre_area[j])
    barea_dict[i] = abs(sum(np.cross(bc2-bc1,bc1-bc0)))/2 #volume of each volume element


#distance between volume elements
#for j in range(velement):
#    for k in range(velement):
#        if j != k:
#            dist_jk = math.sqrt((abs(centre_cell[i][0] - centre_cell[j][0]))**2 + (abs(centre_cell[i][1] - centre_cell[j][1]))**2 + (abs(centre_cell[i][2] - centre_cell[j][2]))**2) #distance between volume elements

##################
#This is for each volume element with its own numbering
##################
#Neighbours calculation; What are the neighbours faces of each volume? 3 per each minimum?
facenodes = gmsh.model.mesh.getElementFaceNodes(4, 3) #4 is the element type (tetrahedron) and three are the nodes per each face

#Computing face x tetrahedon incidence
faces = []
fxt = {} #dictionary with keys as the nodes of each face and values the volume elements of which this face is neighbour
for i in range(0, len(facenodes), 3): # per ecah element basically, goes trhough the nodes of each face 3by3
    #print(i)
    f = tuple(sorted(facenodes[i:i + 3])) #nodes of each face put in a tuple from node i to node i plus 3
    faces.append(f)
    t = voluEl[i // 12] #volume element number at which the faces are associated?
    if not f in fxt: #if the face f (with its node) is alrady in the dictionary, just append the volume element neighbour to
        fxt[f] = [t]
    else:
        fxt[f].append(t)

#Computing neighbors by face
txt = {}
for i in range(0, len(faces)):
    #print(i)
    f = faces[i]
    t = voluEl[i // 4]
    if not t in txt:
        txt[t] = []
    for tt in fxt[f]:
        if tt != t:
            txt[t].append(tt) #volumes neighbours to each volume

#print("--- done: neighbors by face =", txt)

dist_dict = {}
#distance between neighbour volume elements
for j in txt.keys():
    #print(j)
    centre_cell_j = centre_cell[j]
    distances = []
    for k in txt[j]:
        #print(k)
        centre_cell_k = centre_cell[k]
        if j != k:
            dist_jk = math.sqrt((abs(centre_cell[j][0] - centre_cell[k][0]))**2 + (abs(centre_cell[j][1] - centre_cell[k][1]))**2 + (abs(centre_cell[j][2] - centre_cell[k][2]))**2) #distance between volume elements
        distances.append(dist_jk)
    dist_dict[j] = distances

#    for k in range(velement):
#        if j != k:
#            dist_jk = math.sqrt((abs(centre_cell[i][0] - centre_cell[j][0]))**2 + (abs(centre_cell[i][1] - centre_cell[j][1]))**2 + (abs(centre_cell[i][2] - centre_cell[j][2]))**2) #distance between volume elements


#namGroup = gmsh.model.getPhysicalName(dimGroup, tagGroup)


#Get Physical tags, groups to put boundary conditions

vGroups = gmsh.model.getPhysicalGroups()
vGroupsNames = []
for iGroup in vGroups:
    dimGroup = iGroup[0]  #entity tag: 1 lines, 2 surfaces, 3 volumes (1D, 2D or 3D)
    tagGroup = iGroup[1]  #physical tag group 
    namGroup = gmsh.model.getPhysicalName(dimGroup, tagGroup)
    alist = [dimGroup,tagGroup,namGroup]
    vGroupsNames.append(alist)

surfEl = {}
for el in vGroupsNames:
    if el[0] == 2:
        surfaces, bounNodeTagstype = gmsh.model.mesh.getElementsByType(el[0],el[1])
        abscoeff = input(f"Enter absorption coefficient input for {el[2]}:")
        surfaces = np.append(surfaces,abscoeff)
        surfEl[el[2]] = surfaces

        #vGroupsNames[el].append(abscoeff)





farea_dict = {}

for tetrahedron, neighbors in txt.items():
    tetrahedron_faces = fxt.keys()
    fareas = []

    for neighbor in neighbors:
        shared_faces = [face for face, tetrahedron_list in fxt.items() if tetrahedron in tetrahedron_list and neighbor in tetrahedron_list]

        for shared_face in shared_faces:
            fc0 = vnodeCoord_dict[shared_face[0]]
            fc1 = vnodeCoord_dict[shared_face[1]]
            fc2 = vnodeCoord_dict[shared_face[2]]
            farea = abs(np.dot(np.cross(fc2 - fc0, fc1 - fc0), fc2 - fc0)) / 2.0
            fareas.append(farea)

    farea_dict[tetrahedron] = fareas


















#farea_dict = {}
#for i in txt.keys():
#    print(i)
#    fareas = []
#    for value in txt[i]:
#        print(value)
#        thelist = [i,value]
#        for j in fxt.values():
#            if j == thelist:
#                face = #this needs to be the key correspondent to the value j
#                fc0 = vnodeCoord_dict[face[0]] #Coordinates of the node number zero of the volume element i
#                #print(nc0)
#                fc1 = vnodeCoord_dict[face[1]] #Coordinates of the node number one of the volume element i
#                #print(nc1)
#                fc2 = vnodeCoord_dict[face[2]] #Coordinates of the node number two of the volume element i
#                #print(nc2)
#                farea = abs(sum(np.cross(fc2-fc1,fc1-fc0)))/2 #volume of each volume element
#        fareas.append(farea)
#    farea_dict[i] = fareas        
            
            
            
            
#        for key in fxt.keys():
#            print(key)
#            if fxt[key] == thelist:
#                face = key
#                fc0 = vnodeCoord_dict[face[0]] #Coordinates of the node number zero of the volume element i
#                #print(nc0)
#                fc1 = vnodeCoord_dict[face[1]] #Coordinates of the node number one of the volume element i
#                #print(nc1)
#                fc2 = vnodeCoord_dict[face[2]] #Coordinates of the node number two of the volume element i
#                #print(nc2)
#                farea = abs(sum(np.cross(fc2-fc1,fc1-fc0)))/2 #volume of each volume element
#            flist.append(farea)
#    farea_dict[i] = flist

       # face = {i for i in fxt if fxt[i]==thelist}
        
        

gmsh.finalize()


#%%

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

Abs_terms = 10

#Set initial condition - Source Info (interrupted method)
Ws = 0.01 #Source point power [Watts] interrupted after "sourceon_time" seconds; 10^-2 W => correspondent to 100dB
Vs = 4 #Volume of the cell #volume of the source = to volume of cells ###?????

sourceon_time =  0.50 #time that the source is ON before interrupting [s]
recording_time = 2.00 #total time recorded for the calculation [s]

V = sum(vcell_dict.values())
S = 100 # surface area of the room ###???

#%%
###############################################################################
#CALCULATION SECTION
###############################################################################

#Fixed inputs
pRef = 2 * (10**-5) #Reference pressure in Pa
rho = 1.21 #air density [kg.m^-3] at 20°C

#Time resolution
t = np.arange(0, recording_time, dt) #mesh point in time
recording_steps = ceil(recording_time/dt) #number of time steps to consider in the calculation

#Initial condition - Source Info (interrupted method)
w1=Ws/Vs #w1 = round(Ws/Vs,4) #power density of the source [Watts/(m^3))]
sourceon_steps = ceil(sourceon_time/dt) #time steps at which the source is calculated/considered in the calculation
s1 = np.multiply(w1,np.ones(sourceon_steps)) #energy density of source number 1 at each time step position
source1 = np.append(s1, np.zeros(recording_steps-sourceon_steps)) #This would be equal to s1 if and only if recoding_steps = sourceon_steps

#Diffusion parameters
lambda_path = (4*V)/S #mean free path for 3D
Dx = (lambda_path*c0)/3 #diffusion coefficient for proportionate rooms x direction
Dy = (lambda_path*c0)/3 #diffusion coefficient for proportionate rooms y direction
Dz = (lambda_path*c0)/3 #diffusion coefficient for proportionate rooms z direction

Sjk =  10#area face between two adjacent tetrahedron shapes ###?????
hjk = 10 #distance between the central nodes of two adjacent tetrahedron shapes ###?????

beta_zero = np.divide((dt*(Dx*Sjk/hjk - hjk*Abs_terms)),vcell_dict.values()) ###????? my interpretation of the beta_zero

#Finding index in meshgrid of the source position
coord_source = [x_source,y_source,z_source] #coordinates of the receiver position in an list

s = np.zeros((velement,velement)) #matrix of zeros for source

#Finding index in meshgrid of the receiver position
coord_receiver = [x_rec,y_rec,z_rec] #coordinates of the receiver position in an list


#%%
###############################################################################
#MAIN CALCULATION - COMPUTING ENERGY DENSITY
############################################################################### 

w_new = np.zeros((velement,velement)) #unknown w at new time level (n+1)
w = w_new #w at n level
w_old = w #w_old at n-1 level

w_rec = np.arange(0,recording_time,dt) #energy density at the receiver

#Computing w;
for steps in range(0, recording_steps):
    #Compute w at inner mesh points
    time_steps = steps*dt #total time for the calculation
    
    #Computing w_new (w at n+1 time step)
    w_new = np.divide((np.multiply(w_old,(1-beta_zero))),(1+beta_zero)) - \
        np.divide((2*dt*c0*m_atm*w),(1+beta_zero)) + \
            np.divide((np.divide((2*dt*Dx*(Sjk/hjk)*w),vcell_dict.values())),(1+beta_zero)) + \
                np.divide((2*dt*s),(1+beta_zero)) #The absorption term is part of beta_zero
                 ###?????
                 
    #Update w before next step
    w_old = w #The w at n step becomes the w at n-1 step
    w = w_new #The w at n+1 step becomes the w at n step
    print(time_steps)
    #w_rec is the energy density at the specific receiver
#    w_rec[steps] = w_new[row_lr, col_lr, depth_lr]












