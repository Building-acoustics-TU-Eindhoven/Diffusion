# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 16:12:40 2023

@author: 20225533
"""

#%%
###############################################################################
#IMPORT LIBRARIES
###############################################################################
#Code developed by Ilaria Fichera for the analysis of the FVM method adapted solving the 3D diffusion equation with one intermittent omnidirectional sound source
#Import modules
import itertools
import math
import pickle
import os
import time
import types

from math import ceil
from math import log

import gmsh
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from mpl_toolkits.mplot3d import Axes3D
# from matplotlib.animation import FuncAnimation

# from scipy import linalg
# from scipy.integrate import simps
# from scipy import stats
# from scipy.interpolate import griddata
# from scipy.sparse import lil_matrix
from scipy.spatial import cKDTree

from FunctionRT import t60_decay
from FunctionClarity import clarity
from FunctionDefinition import definition
from FunctionCentreTime import centretime

st = time.time() #start time of calculation

#%%
###############################################################################
#INPUT VARIABLES
###############################################################################

# Source position
x_source = 0.5  #position of the source in the x direction [m]
y_source = 0.5  #position of the source in the y direction [m]
z_source = 1.0  #position of the source in the z direction [m]

# Receiver position
x_rec = 2.0 #position of the receiver in the x direction [m]
y_rec = 0.5 #position of the receiver in the y direction [m]
z_rec = 1.0 #position of the receiver in the z direction [m]

# Type of Calculation
#Choose "decay" if the objective is to calculate the energy decay of the room with all its energetic parameters; 
#Choose "stationarysource" if the aim is to understand the behaviour of a room subject to a stationary source
tcalc = "decay"

# Frequency range
fc_low = 125
fc_high = 2000
num_octave = 1

x_frequencies  = num_octave * log(fc_high/fc_low) / log(2)
nBands = int(num_octave * log(fc_high/fc_low) / log(2) + 1)
center_freq = fc_low * np.power(2,((np.arange(0,x_frequencies+1) / num_octave)))

# Time discretization
dt = 1/20000 #time discretization

# Air absorption coefficient
m_atm = 0 #air absorption coefficient [1/m]

#%%
###############################################################################
#FIXED INPUTS
###############################################################################
#General settings
c0= 343 #adiabatic speed of sound [m.s^-1]

#Set initial condition - Source Info (interrupted method)
Ws = 0.01 #Source point power [Watts] interrupted after "sourceon_time" seconds; 10^-2 W => correspondent to 100dB

# Absorption term and Absorption coefficients
th = 3 #int(input("Enter type Absortion conditions (option 1,2,3):")) 
# options Sabine (th=1), Eyring (th=2) and modified by Xiang (th=3)

# Reference pressure in Pa
pRef = 2 * (10**-5) #Reference pressure in Pa

# Air density
rho = 1.21 #air density [kg.m^-3] at 20°C

#%%
###############################################################################
#INITIALISE GMSH
###############################################################################
file_name = "3x3x3.msh" #Insert file name, msh file created from sketchUp and then gmsh
gmsh.initialize() #Initialize msh file
mesh = gmsh.open(file_name) #open the file

#gmsh.fltk.run() #run the file to see it in gmsh
dim = -1 #dimensions of the entities, 0 for points, 1 for curves/edge/lines, 2 for surfaces, 3 for volumes, -1 for all the entities 
tag = -1 #all the nodes of the room

#%%
###############################################################################
#SURFACE MATERIALS
###############################################################################
vGroups = gmsh.model.getPhysicalGroups(-1) #these are the entity tag and physical groups in the msh file. 
vGroupsNames = [] #these are the entity tag and physical groups in the msh file + their names
for iGroup in vGroups:
    dimGroup = iGroup[0]  #entity tag: 1 lines, 2 surfaces, 3 volumes (1D, 2D or 3D)
    tagGroup = iGroup[1]  #physical tag group (depending on material properties defined in SketchUp)
    namGroup = gmsh.model.getPhysicalName(dimGroup, tagGroup) #names of the physical groups defined in SketchUp   
    alist = [dimGroup,tagGroup,namGroup] #creates a list of the entity tag, physical tag group and name
    #print(alist)
    vGroupsNames.append(alist)
    
# Initialize a list to store surface tags and their absorption coefficients
surface_absorption = [] #initialization absorption term (alpha*surfaceofwall) for each wall of the room
triangle_face_absorption = [] #initialization absorption term for each triangle face at the boundary and per each wall
absorption_coefficient = {}

# Absorption term for boundary conditions 
def abs_term(th,abscoeff_list):
    Absx_array = np.array([])
    for abs_coeff in abscoeff_list:
        #print(abs_coeff)
        if th == 1:
            Absx = (c0*abs_coeff)/4 #Sabine
        elif th == 2:
            Absx = (c0*(-log(1-abs_coeff)))/4 #Eyring
        elif th == 3:
            Absx = (c0*abs_coeff)/(2*(2-abs_coeff)) #Modified by Xiang
        Absx_array = np.append(Absx_array, Absx)
    return Absx_array

for group in vGroupsNames:
    if group[0] != 2:
        continue
    name_group = group[2]
    name_split = name_group.split("$")
    name_abs_coeff = name_split[0]
    abscoeff = input(f"Enter absorption coefficient for frequency {fc_low} to {fc_high} for {name_abs_coeff}:") 
    abscoeff = abscoeff.split(",")
    #abscoeff = [float(i) for i in abscoeff][-1] #for one frequency
    abscoeff_list = [float(i) for i in abscoeff] #for multiple frequencies
    
    physical_tag = group[1] #Get the physical group tag
    entities = gmsh.model.getEntitiesForPhysicalGroup(2, physical_tag) #Retrieve all the entities in this physical group (the entities are the number of walls in the physical group)

    Abs_term = abs_term(th, abscoeff_list) #calculates the absorption term based on the type of boundary condition th
    for entity in entities:
        absorption_coefficient[entity] = abscoeff_list
        surface_absorption.append((entity, Abs_term)) #absorption term (alpha*surfaceofwall) for each wall of the room
        surface_absorption = sorted(surface_absorption, key=lambda x: x[0])

for entity, Abs_term in surface_absorption:
    triangle_faces, _ = gmsh.model.mesh.getElementsByType(2, entity) #Get all the triangle faces for the current surface
    triangle_face_absorption.extend([Abs_term] * len(triangle_faces)) #Append the Abs_term value for each triangle face

print("Correctly inputted surface materials. Starting initial geometry calculations...")
#%%
###############################################################################
#GMSH GEOMETRY CALCULATION
###############################################################################
#Nodes
nodeTags, coords, parametricCoord = gmsh.model.mesh.getNodes(dim,tag) #gets the tags for each node and the coordinates of each node
nodecoords = coords.reshape((-1,3)) #coordinates reshaped in a matrix 3xnumber of nodes

node_indices = {tag: index for (index, tag) in enumerate(nodeTags)}

#Element Types
elemTypes,elemTags,elemNodeTags = gmsh.model.mesh.getElements(dim,tag)
#elemTypes = 1 for lines, 2 for surfaces, 4 for tetrahedron
#elemTags =  list of list of lines, boundary elements (surfaces) and volume elements (tetrahedron)
#elemNodeTags = did not understand this yet, probably a tag gives to each line, surface and tetrahedron
for e_type in elemTypes:
    if e_type == 1: #if the e_type = 1, then get all the elements (lines) with that e_type
        tag = -1 
        edgeEl, edgeNodeTagstype = gmsh.model.mesh.getElementsByType(e_type,tag) #get the edge element tags and the nodes; edgeEl are lines
        #edgeEl = numbered edge elements (lines)
    elif e_type == 2: #if the e_type = 2, then get all the elements (surfaces) with that e_type
        tag = -1 
        bounEl, bounNodeTagstype = gmsh.model.mesh.getElementsByType(e_type,tag) #get the boundary/surface element tags and the nodes; bounEl are surfaces
        bounNode = bounNodeTagstype.reshape((-1,3))
        #boundEl = numbered boundary elements (surfaces)
        #bounNode = nodes of the surfaces
    elif e_type == 4:
        tag = -1#if the e_type = 4, then get all the elements (tetrahedron) with that e_type
        voluEl, voluNodeTagstype = gmsh.model.mesh.getElementsByType(e_type,tag) #get the volume element tags and the nodes; voluEl are tetrahedrons
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
volumeEl_dict = {} #Initialization Dictionary of volumelements + its nodes
for i in range(len(voluEl)):
    #print(i)
    volumeEl_dict[voluEl[i]] = velemNodes[i] #Dictionary of volumelements + its nodes
    
#Boundary Element dictionary + node per each surface elements (3 nodes per element)
boundaryEl_dict = {} #Initialization Dictionary of boundary elements + its nodes
for i in range(len(bounEl)):
    boundaryEl_dict[bounEl[i]] = belemNodes[i] #Dictionary of boundary elements + its nodes

#Calculation of volume cells and centre of volume    
vcell_dict = {} #volume of each element tetrahedron initialization
centre_cell = {} #centre of the element tetrahedron initialization
for i in volumeEl_dict.keys():
    coord_centre_cell = np.zeros(3)
    centre_cell[i] = []
    #print(i)
    vc0 = nodecoords[node_indices[volumeEl_dict[i][0]],:] #coordinates of node 0
    #vc0 = gmsh.model.mesh.getNode(volumeEl_dict[i][0])[0] #Coordinates of the node number zero of the volume element i
    #print(nc0)
    vc1 = nodecoords[node_indices[volumeEl_dict[i][1]],:]
    vc2 = nodecoords[node_indices[volumeEl_dict[i][2]],:]
    vc3 = nodecoords[node_indices[volumeEl_dict[i][3]],:]
    #vc1 = gmsh.model.mesh.getNode(volumeEl_dict[i][1])[0] #Coordinates of the node number one of the volume element i
    #print(nc1)
    #vc2 = gmsh.model.mesh.getNode(volumeEl_dict[i][2])[0] #Coordinates of the node number two of the volume element i
    #print(nc2)
    #vc3 = gmsh.model.mesh.getNode(volumeEl_dict[i][3])[0] #Coordinates of the node number three of the volume element i
    #print(nc3)
    for j in range(3): #three coordinates per each node
        coord_centre_cell[j] = (vc0[j]+vc1[j]+vc2[j]+vc3[j])/4 #coordinates of the centre of each volume element
        centre_cell[i].append(coord_centre_cell[j])
    vcell_dict[i] = abs(np.dot(np.cross(vc1-vc3,vc2-vc3),vc0-vc3))/6 #volume of each volume element

#The dictionary centre cell is modified into an array of floats (list)
cell_center = np.array(()) #initialization of an array of centre cell coordinates from the centre_cell dictionary
for key in centre_cell:
    cell_center = np.append(cell_center,centre_cell[key])
cell_center = cell_center.reshape((-1,3))

cell_volume = np.array(()) #initialization of an array of cell volumes from the vcell_dict dictionary
for key in vcell_dict:
    cell_volume = np.append(cell_volume,vcell_dict[key])

#Calculation of boundary elements area and centre   
barea_dict = {} #surface of each element boundary initialization
centre_area = {} #centre of the element tetrahedron initialization
for i in boundaryEl_dict.keys():
    coord_centre_area = np.zeros(3)
    centre_area[i] = []
    #print(i)
    bc0 = nodecoords[node_indices[boundaryEl_dict[i][0]],:]
    
    #bc0 = gmsh.model.mesh.getNode(boundaryEl_dict[i][0])[0]
    #bnodeCoord_dict[boundaryEl_dict[i][0]] #Coordinates of the node number zero of the volume element i
    #print(nc0)
    bc1 = nodecoords[node_indices[boundaryEl_dict[i][1]],:]
    #gmsh.model.mesh.getNode(boundaryEl_dict[i][1])[0] #Coordinates of the node number one of the volume element i
    #print(nc1)
    bc2 = nodecoords[node_indices[boundaryEl_dict[i][2]],:]
    #gmsh.model.mesh.getNode(boundaryEl_dict[i][2])[0] #Coordinates of the node number two of the volume element i
    #print(nc2)
    for j in range(3):
        coord_centre_area[j] = (bc0[j]+bc1[j]+bc2[j])/3 #coordinates of the centre of each volume element
        centre_area[i].append(coord_centre_area[j])
    barea_dict[i] = abs(sum(np.cross(bc2-bc1,bc1-bc0)))/2 #volume of each volume element


#Neighbours calculation; What are the neighbours faces of each volume? 3 per each minimum?
facenodes = gmsh.model.mesh.getElementFaceNodes(4, 3) #4 is the element type (tetrahedron) and three are the nodes per each face #get all the face tags of all the faces of the tetrahedrons

#Computing face x tetrahedon incidence
faces = [] #initialization of list of tuples of the faces nodes
fxt = {} #dictionary with keys as the nodes of each face and values the volume elements of which this face is neighbour
for i in range(0, len(facenodes), 3): # per each element basically, goes trhough the nodes of each face 3by3
    #print(i)
    f = tuple(sorted(facenodes[i:i + 3])) #nodes of each face put in a tuple from node i to node i plus 3
    faces.append(f)
    tet = voluEl[i // 12] #volume element number at which the faces are associated?
    if not f in fxt: #if the face f (with its node) is already in the dictionary, just append the volume element neighbour to
        fxt[f] = [tet]
    else:
        fxt[f].append(tet)

#Computing neighbors by face
txt = {} #dictionary with keys as the tetrahedron tag and values as the tet that are neighbours (the tet at the boundary are not counted)
for i in range(0, len(faces)):
    #print(i)
    f = faces[i] #f is a tuple of the nodes of the face into consideration
    tet = voluEl[i // 4] #tetrahedron at which the face is neighbour
    if not tet in txt: #if the tet is not in the dictionary, add it, otherwise append the new tt
        txt[tet] = []
    for tt in fxt[f]:
        if tt != tet:
            txt[tet].append(int(tt-(voluEl[0]-1))) #volumes neighbours to each volume
for values in txt.values(): 
    if len(values)==2: #if there are only two tetrahedrons neighbours it means that the other two are boundary tetrahedrons, therefore add a zero for each tetrahedron missing
        values.append(0)
        values.append(0)
    if len(values) == 3: #if there are only three tetrahedrons neighbours it means that the other one is a boundary tetrahedron, therefore add a zero
        values.append(0)

neighbourVolume = np.array(()) #initialization of an array with the neighbours tetrahedron per each tetrahedron in order from 0 to the number of tetrahedrons
for key in txt:
    neighbourVolume = np.append(neighbourVolume,txt[key])
for item in neighbourVolume:
        item = int(item)   
    
neighbourVolume = neighbourVolume.reshape((-1,4)) #reshape the array so that we have the 4 tetrahedrons neighbours of the tetrahedron in consideration

print("Completed initial geometry calculation. Starting internal tetrahedrons calculations...")

# Interior Tetrahedrons calculations
interior_tet = np.zeros((velement, velement)) #initialization matrix of tetrahedron per tetrahedron

for i in range(velement): #for each tetrahedron, take its centre
    #print(i)
    cell_center_i = cell_center[i]
    for j in range(velement): #for each tetrahedron, take its centre
        #cell_center_j = cell_center[j]
        #print(j)
        if i != j: #if the tetrahedrons are not the same one, then check if there are shared nodes in between the two tetrahedron i and j
            shared_nodes = np.intersect1d(velemNodes[i], velemNodes[j])
            #shared_nodes = []
            #count = 0
            #for node in velemNodes[i]: #for each node in tetrahedron i
            #    print(node)
            #    if node in velemNodes[j]: #if each node of the tetrahedron i is in nodelist of tetrahedron j
            #        count += 1
            #        shared_nodes.append(node) #append the node that it is in common
            if len(shared_nodes) == 3: #after have done this for all the nodes, if the cound is 3 then calculate the shared area between the tetrahedrons
                sc0 = nodecoords[node_indices[shared_nodes[0]],:]
                #sc0 = gmsh.model.mesh.getNode(shared_nodes[0])[0] #coordinates of node 0
                sc1 = nodecoords[node_indices[shared_nodes[1]],:]
                #sc1 = gmsh.model.mesh.getNode(shared_nodes[1])[0] #coordinates of node 1
                sc2 = nodecoords[node_indices[shared_nodes[2]],:]
                #sc2 = gmsh.model.mesh.getNode(shared_nodes[2])[0] #coordinates of node 2
                shared_area = np.linalg.norm(np.cross(sc2-sc0,sc1-sc0))/2 #compute shared area
                shared_distance = np.linalg.norm(cell_center_i - cell_center[j])
                    #sqrt((abs(cell_center_i[0] - cell_center_j[0]))**2 + (abs(cell_center_i[1] - cell_center_j[1]))**2 + (abs(cell_center_i[2] - cell_center_j[2]))**2) #distance between volume elements
                interior_tet[i, j] = shared_area/shared_distance #division between shared area and shared distance
            else:
                shared_area = 0
                interior_tet[i, j] = shared_area

interior_tet_sum = np.sum(interior_tet, axis=1) #sum of interior_tet per columns (so per i element)
print("Completed internal tetrahedrons calculation. Starting boundary tetrahedrons calculations...")

# Calculation surface areas
surface_areas = {}   
for entity, Abs_term in surface_absorption:   
    face_nodes_per_entity= gmsh.model.mesh.getElementFaceNodes(2, 3, tag=entity)
    surf_area_tot = 0
    for i in range(0, len(face_nodes_per_entity), 3): # per each element basically, goes trhough the nodes of each face 3by3
        #print(i)
        f = tuple(sorted(face_nodes_per_entity[i:i + 3])) 
        fc0 = nodecoords[node_indices[f[0]],:]
        #fc0 = gmsh.model.mesh.getNode(f[0])[0] #coordinates of vertix 0
        fc1 = nodecoords[node_indices[f[1]],:]
        fc2 = nodecoords[node_indices[f[2]],:]
        #fc1 = gmsh.model.mesh.getNode(f[1])[0] #coordinates of vertix 1
        #fc2 = gmsh.model.mesh.getNode(f[2])[0] #coordinates of vertix 2
        face_area = 0.5 * np.linalg.norm(np.cross(fc1 - fc0, fc2 - fc0)) #Compute the area using half of the cross product's magnitude
        surf_area_tot += face_area
        surface_areas[entity] = surf_area_tot


#FACE AREA & boundary_areas
total_boundArea = 0 #initialization of total surface area of the room
boundary_areas = []  #Initialize a list to store boundary_areas values for each tetrahedron

face_areas = np.zeros(len(velemNodes)) #Per each tetrahedron, if there is a face that is on the boundary, include the area, otehrwise zero
for idx, element in enumerate(velemNodes): #for index and element in the number of tetrahedrons
    #if idx == 491:
        tetrahedron_boundary_areas = 0 #initialization tetrahedron face on boundary*its absorption term
        total_tetrahedron_boundary_areas = np.zeros(nBands) #initialization total tetrahedron face on boundary*its absorption term if there are more than one face in the tetrahedron that is on the boundary
        #print(idx)
        node_combinations = [list(nodes) for nodes in itertools.combinations(element, 3)] #all possible combinations of the nodes of the tetrahedrons (it checks also for the order of the nodes in the same combination)
        # Check if the nodes are in any order in bounNode
        is_boundary = False #variable to say that at the beginning the face in not on a boundary
        for nodes in node_combinations: #for each node in each combination
            for surface_idx, surface in enumerate(bounNode): #for index and surface in the number of nodes
                surface_set = sorted(set(surface)) #creates a set of the surface nodes
                surface_set_idx = surface_idx
                nodes_set = sorted(set(nodes)) #create a set of the node combination of the tetrahedron into consideration
                surface_list = list(surface)
                if nodes_set == surface_set: #if these are equal, it means that the tetrahedron into consideration has a surface in the boundary and therefore is_boundary gets the value of True.
                    #print(surface_set)
                    #print(surface_list)
                    is_boundary = True
                    if is_boundary: #if the surface is at the boundary, then take the coordinates of each vertix
                        #Convert the vertices to NumPy arrays for vector operations
                        bc0 = nodecoords[node_indices[nodes[0]],:]
                        bc1 = nodecoords[node_indices[nodes[1]],:]
                        bc2 = nodecoords[node_indices[nodes[2]],:]
                        
                        
                        #bc0 = gmsh.model.mesh.getNode(nodes[0])[0] #coordinates of vertix 0
                        #bc1 = gmsh.model.mesh.getNode(nodes[1])[0] #coordinates of vertix 1
                        #bc2 = gmsh.model.mesh.getNode(nodes[2])[0] #coordinates of vertix 2
                        
                        face_area = 0.5 * np.linalg.norm(np.cross(bc1 - bc0, bc2 - bc0)) #Compute the area using half of the cross product's magnitude
                        #print(face_area)
                        
                        face_areas[idx] = face_area #area of the surface that is on boundary per each tetrahedron
                        total_boundArea += face_area #add to the total boundary area
                        
                        if face_area > 0:
                            # Use the index to access the corresponding absorption area
                            face_absorption_product = face_area * triangle_face_absorption[surface_set_idx] #calculate the product between the area*the correspondent absorption term
                            #print(face_absorption_product)
                            
                            tetrahedron_boundary_areas += face_absorption_product #add the calculation to the tetrahedron correspondent
                            
                            total_tetrahedron_boundary_areas = tetrahedron_boundary_areas #if there are multiple surfaces on the boundary per each tetrahedron, then add also the second and the third one
                            
        boundary_areas.append(np.array(total_tetrahedron_boundary_areas)) #Append the total boundary_areas for the tetrahedron to the list
        #print(total_tetrahedron_boundary_areas)

print("Completed boundary tetrahedrons calculation. Starting main diffusion equation calculations over time and frequency...")

gmsh.finalize()

#%%
###############################################################################
#CALCULATION SECTION
###############################################################################

V = sum(cell_volume)
S = total_boundArea #surface area of the room

sum_alpha_average = 0
Eq_A = 0
#Absorption parameters for room
for entity in surface_areas:
    #print(entity)
    sum_alpha_average += np.multiply(absorption_coefficient[entity],surface_areas[entity])
    Eq_A += np.multiply(absorption_coefficient[entity],surface_areas[entity])
alpha_average = sum_alpha_average/S #average absorption

#%%
###############################################################################
#TOTAL RECORDING TIME AND SOURCE ON TIME
###############################################################################

RT_Sabine_band = []
for iBand in range(nBands):
    freq = center_freq[iBand]
    RT_Sabine = 0.16*V/Eq_A[iBand]
    RT_Sabine_band.append(RT_Sabine)

sourceon_time = round(max(RT_Sabine_band),1)#time that the source is ON before interrupting [s]
recording_time = 2*sourceon_time #total time recorded for the calculation [s]

#%%
###############################################################################
#TIME RESOLUTION
###############################################################################
#Time resolution
t = np.arange(0, recording_time, dt) #mesh point in time
recording_steps = ceil(recording_time/dt) #number of time steps to consider in the calculation

#%%
###############################################################################
#DIFFUSION PARAMETERS
###############################################################################
#Diffusion parameters
mean_free_path = (4*V)/S #mean free path for 3D
mean_free_time= mean_free_path/c0 #mean free time for 3D
mean_free_time_step = int(mean_free_time/dt)
Dx = (mean_free_path*c0)/3 #diffusion coefficient for proportionate rooms x direction
Dy = (mean_free_path*c0)/3 #diffusion coefficient for proportionate rooms y direction
Dz = (mean_free_path*c0)/3 #diffusion coefficient for proportionate rooms z direction

#%%
###############################################################################
#SOURCE & RECEIVER
###############################################################################

#distance between source and receiver
dist_sr = math.sqrt((abs(x_rec - x_source))**2 + (abs(y_rec - y_source))**2 + (abs(z_rec - z_source))**2) #distance between source and receiver

coord_source = [x_source,y_source,z_source] #coordinates of the receiver position in an list
coord_rec = [x_rec,y_rec,z_rec] #coordinates of the receiver position in an list

#%%
###############################################################################
#SOURCE INTERPOLATION
###############################################################################

#ORIGINAL
#Position of source is the centre of a cell so the minimum distance with the centre of a cell has been calculated to understand which cell is the closest
# dist_source_cc_list = []
# for i in range(len(cell_center)):
#     dist_source_cc = math.sqrt(np.sum((cell_center[i] - coord_source)**2))
#     dist_source_cc_list.append(dist_source_cc)
# source_idx = np.argmin(dist_source_cc_list)

#SOURCE INTERPOLATION CALCULATED WITHIN N CENTRE CELL SELECTED WITHIN THE LENGTH OF MESH
#Position of source is the centre of a cell so the minimum distance with the centre of a cell has been calculated to understand which cell is the closest
# dist_source_cc_list = []
# for i in range(len(cell_center)):
#     dist_source_cc = math.sqrt(np.sum((cell_center[i] - coord_source)**2))
#     dist_source_cc_list.append(dist_source_cc)
# source_idx = np.argmin(dist_source_cc_list)
# import math
# #cl_tet_s stands for cl=closest, tet=tetrahedron, s=to the source
# cl_tet_s = {} #initialise dictionary fro closest tetrahedrons to the source
# for i in range(len(dist_source_cc_list)):
#     if dist_source_cc_list[i] < length_of_mesh: #the number of closest tetrahedron to the source depends on the mesh and I would say that 0.5 should be equal to the length of mesh
#         cl_tet_s[i] = dist_source_cc_list[i]
# total_weights_s = {} #initialise weights for each tetrahedron around the actual source position
# sum_weights_s = 0
# Vs = 0
# for i, dist in cl_tet_s.items(): #for each key and value in the dictionary (so for each closest tetrahedron to the source)
#     weights = np.divide(1.0 , dist)  #calculate the inverse distance weights, so closer to the point means higher weight
#     #print(weights)
#     sum_weights_s += weights
#     #weights /= np.sum(weights)  # Normalize weights to sum to 1
#     total_weights_s[i] = weights #put the wweigths (values) to the correspondent closest tetrahedron (keys)
# #    Vs += cell_volume[i] #volume of the source calculated summing the volumes of all the tetrahedrons involved
# #total_weights_s_values = total_weights_s.values()
# for i,weight in total_weights_s.items():
#     total_weights_s[i] = weight/sum_weights_s if sum_weights_s != 0 else 0
# cl_tet_s_keys = cl_tet_s.keys() #take only the keys of the cl_tet_s dictionary (so basically the indexes of the tetrahedrons)

#SOURCE INTERPOLATION CALCULATED WITHIN 4 CENTRE CELL SELECTED (TETRAHEDRON)
#Position of source is the centre of a cell so the minimum distance with the centre of a cell has been calculated to understand which cell is the closest
dist_source_cc_list = [] #initialise the list for all the distances between each cell centre and the source
for i in range(len(cell_center)): #for each tetra
    dist_source_cc = math.sqrt(np.sum((cell_center[i] - coord_source)**2)) #calculate the distance between its centre cell and the source coordinate
    dist_source_cc_list.append(dist_source_cc) #append the distance in a list
source_idx = np.argmin(dist_source_cc_list) #take the minimum distance index; this is where the source will be positioned

dist_source_cc_list_sorted = sorted(dist_source_cc_list) #sorted from the minimum to the maximum distance
selected_source_cc_list = dist_source_cc_list_sorted[:4] #take only the first four element of the sorted list (take the first 4 cell centres closest to the source)

dist_source_cc_list_sorted_indices = np.argsort(dist_source_cc_list)[:4] #takes the indeces of the minimum distances
selected_source_cc_list_indices = dist_source_cc_list_sorted_indices[:4] #does exactly the same as the previous line

#cl_tet_s stands for cl=closest, tet=tetrahedron, s=to the source
cl_tet_s = {} #initialise dictionary for closest tetrahedrons to the source
for i in range(len(dist_source_cc_list_sorted_indices)): #for each of the 4 tetra closest to the source
    idx = selected_source_cc_list[i] #take its index
    cl_tet_s[dist_source_cc_list_sorted_indices[i]] = idx #put it in a dictionary

total_weights_s = {} #initialise weights for each tetrahedron around the actual source position
sum_weights_s = 0
#Vs = 0
for i, dist in cl_tet_s.items(): #for each key and value in the dictionary (so for each closest tetrahedron to the source)
    weights = np.divide(1.0 , dist)  #calculate the inverse distance weights, so closer to the point means higher weight
    #print(weights)
    sum_weights_s += weights
    #weights /= np.sum(weights)  # Normalize weights to sum to 1
    total_weights_s[i] = weights #put the wweigths (values) to the correspondent closest tetrahedron (keys)
    #Vs += cell_volume[i] #volume of the source calculated summing the volumes of all the tetrahedrons involved

#total_weights_s_values = total_weights_s.values()
for i,weight in total_weights_s.items():
    total_weights_s[i] = weight/sum_weights_s if sum_weights_s != 0 else 0

cl_tet_s_keys = cl_tet_s.keys() #take only the keys of the cl_tet_s dictionary (so basically the indexes of the tetrahedrons)

#%%
###############################################################################
#SOURCE VOLUME DEFINITION
###############################################################################
#To make sure that the source is in the correct tetrahedron position
node_ids = velemNodes.T
ori=nodecoords[node_ids[0,:]-1,:]
#ori = nodecoords[node_indices[node_ids[:,0],:]]
v_tet_s1=nodecoords[node_ids[1,:]-1,:]-ori
v_tet_s2=nodecoords[node_ids[2,:]-1,:]-ori
v_tet_s3=nodecoords[node_ids[3,:]-1,:]-ori
n_tet=len(node_ids.T)
v1s = v_tet_s1.T.reshape((3,1,n_tet))
v2s = v_tet_s2.T.reshape((3,1,n_tet))
v3s = v_tet_s3.T.reshape((3,1,n_tet))
mat = np.concatenate((v1s,v2s,v3s), axis=1)
inv_mat = np.linalg.inv(mat.T).T
#if coord_source.size==3:  # to make rec has a dimension of (N_rec,3)
#    rec=rec.reshape((1,3))
coord_source_array = np.array(coord_source)
if coord_source_array.size==3:  # to make rec has a dimension of (N_rec,3)
    coord_source_array=coord_source_array.reshape((1,3))
N_sou=coord_source_array.shape[0]
orir=np.repeat(ori[:,:,np.newaxis], N_sou, axis=2)
newp=np.einsum('imk,kmj->kij',inv_mat,coord_source_array.T-orir)
val=np.all(newp>=0, axis=1) & np.all(newp <=1, axis=1) & (np.sum(newp, axis=1)<=1)
id_tet, id_p = np.nonzero(val)
res = -np.ones(N_sou, dtype=id_tet.dtype) # Sentinel value
res[id_p]=id_tet

###############################################################################

# #VOLUME CALCULATED WITHIN 4 CENTRE CELL SELECTED (TETRAHEDRON)
#Calculate volume of the source with the 4 cell centres as the vertices of the tetrahedron
# vertices_source = np.array([]).reshape(0, 3)  # Initialize as an empty 2D array with 3 columns
# for tet in cl_tet_s_keys:
#     vertices = cell_center[tet]
#     vertices_source = np.vstack((vertices_source, vertices))

# # Calculate edge vectors
# AB = vertices_source[0] - vertices_source[3]
# CD = vertices_source[1] - vertices_source[3]
# EF = vertices_source[2] - vertices_source[3]

# Calculate volume using scalar triple product
# Vs = np.abs(np.dot(AB, np.cross(CD, EF))) / 6.0


#VOLUME CALCULATED WITHIN N CENTRE CELL SELECTED
# vertices_source = np.array([]).reshape(0, 3)  # Initialize as an empty 2D array with 3 columns
# for tet in cl_tet_s_keys:
#     vertices = cell_center[tet]
#     vertices_source = np.vstack((vertices_source, vertices)) 
# from scipy.spatial import ConvexHull

# def calculate_volume(vertices_source):
#     # Create a ConvexHull object
#     hull = ConvexHull(vertices_source)

#     # Calculate the volume of the convex hull
#     volume = hull.volume

#     return volume
# Vs = calculate_volume(vertices_source)

#VOLUME ORIGINAL
Vs = cell_volume[res[0]] #volume of the source = to volume of cell where the source is 
# Vs = 1

################SOURCE INTERPOLATION WITH VERTICES OF TETRAHEDRON IN WHICH SOURCE IS IN
#Position of source is the centre of a cell so the minimum distance with the centre of a cell has been calculated to understand which cell is the closest
# dist_source_cc_list = []
# for i in range(len(cell_center)):
#     dist_source_cc = math.sqrt(np.sum((cell_center[i] - coord_source)**2))
#     dist_source_cc_list.append(dist_source_cc)
# source_idx = np.argmin(dist_source_cc_list)

# node_source_in = voluNode[source_idx]
# vertices_source_in = np.array([]).reshape(0, 3) 
# for vx in node_source_in:
#     vertices = nodecoords[vx]
#     vertices_source_in = np.vstack((vertices_source_in, vertices))
    
# dist_source_vx_list = []
# for i in range(len(vertices_source_in)):
#     dist_source_vx = math.sqrt(np.sum((vertices_source_in[i] - coord_source)**2))
#     dist_source_vx_list.append(dist_source_vx)


# total_weights_s = {} #initialise weights for the tetrahedron around the actual source position
# sum_weights_s = 0
# Vs = 0
# for i in range(len(dist_source_vx_list)): #for each key and value in the dictionary (so for each closest tetrahedron to the source)
#     weights = np.divide(1.0 , dist_source_vx_list[i])  #calculate the inverse distance weights, so closer to the point means higher weight
#     #print(weights)
#     sum_weights_s += weights
#     #weights /= np.sum(weights)  # Normalize weights to sum to 1
#     total_weights_s[i] = weights #put the wweigths (values) to the correspondent closest tetrahedron (keys)
#     #Vs += cell_volume[i] #volume of the source calculated summing the volumes of all the tetrahedrons involved

# #total_weights_s_values = total_weights_s.values()
# for i,weight in total_weights_s.items():
#     total_weights_s[i] = weight/sum_weights_s if sum_weights_s != 0 else 0

# Vs = cell_volume[source_idx] #volume of the source = to volume of cells where the volume is 

#Initial condition - Source Info (interrupted method)
w1=Ws/Vs #w1 = round(Ws/Vs,4) #power density of the source [Watts/(m^3))]
sourceon_steps = ceil(sourceon_time/dt) #time steps at which the source is calculated/considered in the calculation
s1 = np.multiply(w1,np.ones(sourceon_steps)) #energy density of source number 1 at each time step position
source1 = np.append(s1, np.zeros(recording_steps-sourceon_steps)) #This would be equal to s1 if and only if recoding_steps = sourceon_steps

#ORIGINAL SOURCE MATRIX
# s = np.zeros((velement)) #matrix of zeros for source
# s[source_idx] = source1[0]

#INTERPOLATION WITH CELL CENTRES - SOURCE MATRIX
s = np.zeros((velement)) #matrix of zeros for source
for tet_s in cl_tet_s_keys:
    s[tet_s] = source1[0] *total_weights_s[tet_s]

#INTERPOLATION WITH VERTICES - SOURCE MATRIX
#s = np.zeros((velement)) #matrix of zeros for source
#for i in total_weights_s:
#    s[source_idx] += source1[0] *total_weights_s[i]

#%%
###############################################################################
#RECEIVER INTERPOLATION
###############################################################################
#ORIGINAL
# dist_rec_cc_list = []
# for i in range(len(cell_center)):
#     dist_rec_cc = math.sqrt(np.sum((cell_center[i] - coord_rec)**2))
#     dist_rec_cc_list.append(dist_rec_cc)
# rec_idx = np.argmin(dist_rec_cc_list)

#INTERPOLATION WITH CLOSEST CENTRES WITHIN THE LENGTH OF MESH
#Position of receiver is the centre of a cell so the minimum distance with the centre of a cell has been calculated to understand which cell is the closest
# dist_rec_cc_list = []
# for i in range(len(cell_center)):
#     dist_rec_cc = math.sqrt(np.sum((cell_center[i] - coord_rec)**2))
#     dist_rec_cc_list.append(dist_rec_cc)
# rec_idx = np.argmin(dist_rec_cc_list)


# #cl_tet_r stands for cl=closest, tet=tetrahedron, r=to the receiver
# cl_tet_r = {} #initialise dictionary fro closest tetrahedrons to the source
# for i in range(len(dist_rec_cc_list)):
#     if dist_rec_cc_list[i] < length_of_mesh: #the number of closest tetrahedron to the source depends on the mesh and I would say that 0.5 should be equal to the length of mesh
#         cl_tet_r[i] = dist_rec_cc_list[i]

# total_weights_r = {} #initialise weights for each tetrahedron around the actual source position
# sum_weights_r = 0
# for i, dist in cl_tet_r.items(): #for each key and value in the dictionary (so for each closest tetrahedron to the source)
#     weights = np.divide(1.0 , dist)  #calculate the inverse distance weights, so closer to the point means higher weight
#     #print(weights)
#     sum_weights_r += weights
#     #weights /= np.sum(weights)  # Normalize weights to sum to 1
#     total_weights_r[i] = weights #put the wweigths (values) to the correspondent closest tetrahedron (keys)

# #total_weights_s_values = total_weights_s.values()
# for i,weight in total_weights_r.items():
#     total_weights_r[i] = weight/sum_weights_r if sum_weights_r != 0 else 0

# cl_tet_r_keys = cl_tet_r.keys() #take only the keys of the cl_tet_s dictionary (so basically the indexes of the tetrahedrons)

#RECEIVER INTERPOLATION CALCULATED WITHIN 4 CENTRE CELL SELECTED (TETRAHEDRON)
#Position of receiver is the centre of a cell so the minimum distance with the centre of a cell has been calculated to understand which cell is the closest
dist_rec_cc_list = [] #initialise the list for all the distances between each cell centre and the source
for i in range(len(cell_center)): #for each tetra
    dist_rec_cc = math.sqrt(np.sum((cell_center[i] - coord_rec)**2)) #calculate the distance between its centre cell and the source coordinate
    dist_rec_cc_list.append(dist_rec_cc) #append the distance in a list
rec_idx = np.argmin(dist_rec_cc_list) #take the minimum distance index; this is where the source will be positioned

dist_rec_cc_list_sorted = sorted(dist_rec_cc_list) #sorted from the minimum to the maximum distance
selected_rec_cc_list = dist_rec_cc_list_sorted[:4] #take only the first four element of the sorted list (take the first 4 cell centres closest to the source)

dist_rec_cc_list_sorted_indices = np.argsort(dist_rec_cc_list)[:4] #takes the indeces of the minimum distances
selected_rec_cc_list_indices = dist_rec_cc_list_sorted_indices[:4] #does exactly the same as the previous line

#cl_tet_s stands for cl=closest, tet=tetrahedron, s=to the source
cl_tet_r = {} #initialise dictionary for closest tetrahedrons to the source
for i in range(len(dist_rec_cc_list_sorted_indices)): #for each of the 4 tetra closest to the source
    idx = selected_rec_cc_list[i] #take its index
    cl_tet_r[dist_rec_cc_list_sorted_indices[i]] = idx #put it in a dictionary

total_weights_r = {} #initialise weights for each tetrahedron around the actual source position
sum_weights_r = 0
#Vs = 0
for i, dist in cl_tet_r.items(): #for each key and value in the dictionary (so for each closest tetrahedron to the source)
    weights = np.divide(1.0 , dist)  #calculate the inverse distance weights, so closer to the point means higher weight
    #print(weights)
    sum_weights_r += weights
    #weights /= np.sum(weights)  # Normalize weights to sum to 1
    total_weights_r[i] = weights #put the wweigths (values) to the correspondent closest tetrahedron (keys)
    #Vs += cell_volume[i] #volume of the source calculated summing the volumes of all the tetrahedrons involved

#total_weights_s_values = total_weights_s.values()
for i,weight in total_weights_r.items():
    total_weights_r[i] = weight/sum_weights_r if sum_weights_r != 0 else 0

cl_tet_r_keys = cl_tet_r.keys() #take only the keys of the cl_tet_s dictionary (so basically the indexes of the tetrahedrons)

#%%
###############################################################################
#CALCULATION OF LENGTH OF ROOM
###############################################################################
#Extract x-coordinates of all nodes
x_coordinates = nodecoords[:, 0]

#Find the minimum and maximum x-coordinates to determine the length of the room
min_x = np.min(x_coordinates)
max_x = np.max(x_coordinates)

#Calculate the length of the room
room_length = max_x - min_x

#%%
###############################################################################
#CALCULATION OF WIDTH OF ROOM
###############################################################################
#Extract y-coordinates of all nodes
y_coordinates = nodecoords[:, 1]

#Find the minimum and maximum x-coordinates to determine the width of the room
min_y = np.min(y_coordinates)
max_y = np.max(y_coordinates)

#Calculate the width of the room
room_width = max_y - min_y

#%%
###############################################################################
#CALCULATION OF HEIGHT OF ROOM
###############################################################################
#Extract y-coordinates of all nodes
z_coordinates = nodecoords[:, 2]

#Find the minimum and maximum x-coordinates to determine the width of the room
min_z = np.min(z_coordinates)
max_z = np.max(z_coordinates)

#Calculate the height of the room
room_height = max_z - min_z

#%%
###############################################################################
#LINESPACE LINES OVER THE X AND Y SPACE AND Z AXIS
###############################################################################
#Arange linespace lines
dx = 0.5
x_axis = np.arange(0,room_length+dx,dx) #lispace on x_axis with distance dx
y_axis = np.arange(0,room_width+dx,dx)
z_axis = np.arange(0,room_height+dx,dx)

#%%
###############################################################################
#CALCULATION OF RECEIVERS IN A X LINE
###############################################################################
line_rec_x_idx_list = []
dist_x = np.array([])
for x_chang in x_axis:
    line_rec = [x_chang, y_rec, z_rec]
    #Position of line_receiver is the centre of a cell
    dist_line_rec_x =  math.sqrt((abs(line_rec[0] - x_source))**2 + (abs(line_rec[1] - y_source))**2 + (abs(line_rec[2] - z_source))**2) #distance between source and line_receiver
    dist_x = np.append(dist_x, dist_line_rec_x)  # Append to the NumPy array
    dist_line_rec_x_cc_list = []
    for i in range(len(cell_center)):
        dist_line_rec_x_cc = math.sqrt(np.sum((cell_center[i] - line_rec)**2))
        dist_line_rec_x_cc_list.append(dist_line_rec_x_cc)
    line_rec_x_idx = np.argmin(dist_line_rec_x_cc_list)
    line_rec_x_idx_list.append(line_rec_x_idx)   

#%%
###############################################################################
#CALCULATION OF RECEIVERS IN A Y LINE
###############################################################################
line_rec_y_idx_list = []
dist_y = np.array([])
for y_chang in y_axis:
    line_rec = [x_rec, y_chang, z_rec]
    #Position of line_receiver is the centre of a cell
    dist_line_rec_y =  math.sqrt((abs(line_rec[0] - x_source))**2 + (abs(line_rec[1] - y_source))**2 + (abs(line_rec[2] - z_source))**2) #distance between source and line_receiver
    dist_y = np.append(dist_y, dist_line_rec_y)  # Append to the NumPy array
    dist_line_rec_y_cc_list = []
    for i in range(len(cell_center)):
        dist_line_rec_y_cc = math.sqrt(np.sum((cell_center[i] - line_rec)**2))
        dist_line_rec_y_cc_list.append(dist_line_rec_y_cc)
    line_rec_y_idx = np.argmin(dist_line_rec_y_cc_list)
    line_rec_y_idx_list.append(line_rec_y_idx)  


#def interpolate_receiver_position(interior_tet_tet, cell_centers, receiver_position):
#    distances = np.sqrt(np.sum((cell_centers - receiver_position)**2, axis=1))
#    weights = interior_tet_tet[rec_idx, :] / distances  # Use the row corresponding to the receiver index
#    interpolated_position = np.dot(weights, cell_centers) / np.sum(weights)
#    return interpolated_position

#interpolated_receiver_position = interpolate_receiver_position(interior_tet_tet, cell_center, coord_rec)

#dist_sr_interpolated = np.linalg.norm(interpolated_receiver_position - coord_source)

#%%
###############################################################################
#CALCULATION OF RECEIVERS IN A GRID LINE (HEATMAP)
###############################################################################
target_height = 1.5  # height in meters where the heatmap is desired
grid_x, grid_y = np.meshgrid(x_axis,y_axis) #Definition of the grid at the target height

#Find the closest tetrahedron center point for each grid point
tetra_points_2d = cell_center[:, :2]  # Only x, y coordinates
tree = cKDTree(tetra_points_2d)
grid_points_2d = np.c_[grid_x.ravel(), grid_y.ravel()]
distances, indices = tree.query(grid_points_2d)

#%%
###############################################################################
#CALCULATION OF BETA_ZERO
###############################################################################

boundary_areas = np.array(boundary_areas)
boundary_areas = boundary_areas.T
beta_zero_freq = []
for iBand in range(len(boundary_areas)):
    #print(iBand)
    #freq = center_freq[iBand]
    #print(boundary_areas[iBand])
    beta_zero_element = np.divide(dt*((Dx *interior_tet_sum) + boundary_areas[iBand]),cell_volume) #my interpretation of the beta_zero
    beta_zero_freq.append(beta_zero_element)

#%%
###############################################################################
#MAIN CALCULATION - COMPUTING ENERGY DENSITY
############################################################################### 

w_new_band = []
w_rec_band = []
w_rec_off_band = []
w_rec_off_deriv_band = []
p_rec_off_deriv_band = []
w_t0_band = []
    
curPercent = 0

for iBand in range(nBands):
    freq = center_freq[iBand]    
    
    w_new = np.zeros(velement) #unknown w at new time level (n+1)
    #w_old = np.zeros(velement) 
    w = w_new #w at n level
    w_old = w #w_old at n-1 level
    #w[source_idx] = w1 #w (m time step) at source position -> impulse source
    
    w_rec = np.zeros(recording_steps) #energy density at the receiver
    
    #Computing w;
    for steps in range(0, recording_steps):
        #Compute w at inner mesh points
        time_steps = steps*dt #total time for the calculation
        
        #Computing w_new (w at n+1 time step)
                    
        w_new = np.divide((np.multiply(w_old,(1-beta_zero_freq[iBand]))),(1+beta_zero_freq[iBand])) - \
            np.divide((2*dt*c0*m_atm*w),(1+beta_zero_freq[iBand])) + \
                np.divide(np.divide((2*dt*Dx*(interior_tet@w)),cell_volume),(1+beta_zero_freq[iBand])) + \
                    np.divide((2*dt*s),(1+beta_zero_freq[iBand])) #The absorption term is part of beta_zero
                     
        #Update w before next step
        w_old = w #The w at n step becomes the w at n-1 step
        w = w_new #The w at n+1 step becomes the w at n step
    
        #w_new_total = 0
        #w_rec is the energy density at the specific receiver
        #ORIGINAL
        #w_rec[steps] = w_new[rec_idx]
        
        #INTERPOLATION WITH N CELL CENTRES OR 4 CELL CENTRES
        for tet_r in cl_tet_r_keys:
            w_rec[steps] += w_new[tet_r] *total_weights_r[tet_r]   
        
        if steps == sourceon_steps:
            #print("Steps for source:",steps)
            w_t0 = w_new
    
        if steps == round(1*mean_free_time_step + sourceon_steps + (dist_sr/c0)):
            w_1l = w_new
            
        if steps == round(2*mean_free_time_step + sourceon_steps + (dist_sr/c0)):
            w_2l = w_new
            
        if steps == round(3*mean_free_time_step + sourceon_steps + (dist_sr/c0)):
            w_3l = w_new
    
        if steps == round(5*mean_free_time_step + sourceon_steps + (dist_sr/c0)):
            w_5l = w_new
        
        if tcalc == "decay":
            #INTERPOLATION WITH N CELL CENTRES OR 4 CELL CENTRES
            for tet_s in cl_tet_s_keys:
                 s[tet_s] = source1[steps] *total_weights_s[tet_s]
            #ORIGINAL
            #s[source_idx] = source1[steps]
            #INTERPOLATION WITH VERTICES
            #for i in total_weights_s:
            #    s[source_idx] = source1[steps] *total_weights_s[i]
            
        
        if tcalc == "stationarysource":
            #ORIGINAL
            #s[source_idx] = source1[0]
            
            #INTERPOLATION SOURCE
            for tet_s in cl_tet_s_keys:
                 s[tet_s] = source1[0] *total_weights_s[tet_s]
        
        idx_w_rec = np.argmin(np.abs(t - sourceon_time)) #index at which the t array is equal to the sourceon_time; I want the RT to calculate from when the source stops.
        w_rec_off = w_rec[idx_w_rec:]
        p_rec_off = (w_rec[idx_w_rec:])*rho*c0**2
        t_off = t[idx_w_rec:]
        
        #Envelope of Impulse response from the energy density
        w_rec_off_deriv = w_rec_off #initialising an array of derivative equal to the w_rec_off -> this will be the impulse response after modifying it
        w_rec_off_deriv = np.delete(w_rec_off_deriv, 0) #delete the first element of the array -> this means shifting the array one step before and therefore do a derivation
        w_rec_off_deriv = np.append(w_rec_off_deriv,0) #add a zero in the last element of the array -> for derivation and to have the same length as previously
        impulse = ((w_rec_off - w_rec_off_deriv))/dt#/(rho*c0**2) #This is the difference between the the energy density and the impulse response 
        
        #Envelope of Impulse response from the pressure
        p_rec_off_deriv = p_rec_off #initialising an array of derivative equal to the w_rec_off -> this will be the impulse response after modifying it
        p_rec_off_deriv = np.delete(p_rec_off_deriv, 0) #delete the first element of the array -> this means shifting the array one step before and therefore do a derivation
        p_rec_off_deriv = np.append(p_rec_off_deriv,0) #add a zero in the last element of the array -> for derivation and to have the same length as previously
        
        #print(time_steps)
        
    percentDone = round(100*iBand/nBands);
    #if (percentDone > curPercent):
    print(str(percentDone) + "% of main calculation completed")
        #curPercent += 1;

    w_t0_band.append(w_t0)
    w_new_band.append(w_new)
    w_rec_band.append(w_rec)
    w_rec_off_band.append(w_rec_off)
    w_rec_off_deriv_band.append(w_rec_off_deriv)
    p_rec_off_deriv_band.append(p_rec_off_deriv)

plt.show()

print("100% of main calculation completed")

#%%
###############################################################################
#RESULTS
###############################################################################

w_rec_x_band = []
w_rec_y_band = []
spl_stat_x_band = []
spl_stat_y_band = []
spl_r_band = []
spl_r_off_band = []
spl_r_norm_band = []
t30_band = []
edt_band = []
sch_db_band = []
spl_new_band = []
spl_t0_band = []
w_rec_x_t0_band = []
w_rec_y_t0_band = []
spl_rec_x_t0_band = []
spl_rec_y_t0_band = []

for iBand in range(nBands):
    
    w_rec_x_end = np.array([])
    w_rec_x_t0 = np.array([])
    for xr in line_rec_x_idx_list:
        w_rec_x_new = w_new_band[iBand][xr]
        w_rec_x_end = np.append(w_rec_x_end, w_rec_x_new)
        w_rec_x0 = w_t0_band[iBand][xr]
        w_rec_x_t0 = np.append(w_rec_x_t0, w_rec_x0)
        
    w_rec_y_end = np.array([])
    w_rec_y_t0 = np.array([])
    for yr in line_rec_y_idx_list:
        w_rec_y = w_new_band[iBand][yr]
        w_rec_y_end = np.append(w_rec_y_end, w_rec_y)
        w_rec_y0 = w_t0_band[iBand][yr]
        w_rec_y_t0 = np.append(w_rec_y_t0, w_rec_y0)
        
    spl_t0 = 10*np.log10(((abs(w_t0_band[iBand]))*rho*(c0**2))/(pRef**2))

    spl_new = 10*np.log10(((abs(w_new_band[iBand]))*rho*(c0**2))/(pRef**2))

    spl_stat_x = 10*np.log10(rho*c0*(((Ws)/(4*math.pi*(dist_x**2))) + ((abs(w_rec_x_end)*c0)))/(pRef**2))
    spl_stat_y = 10*np.log10(rho*c0*(((Ws)/(4*math.pi*(dist_y**2))) + ((abs(w_rec_y_end)*c0)))/(pRef**2)) #It should be the spl stationary 

    press_r = ((abs(w_rec_band[iBand]))*rho*(c0**2)) #pressure at the receiver
    spl_r = 10*np.log10(((abs(w_rec_band[iBand]))*rho*(c0**2))/(pRef**2)) #,where=press_r>0, sound pressure level at the receiver
    spl_r_off = 10*np.log10(((abs(w_rec_off_band[iBand]))*rho*(c0**2))/(pRef**2))
    
    spl_r_norm = 10*np.log10((((abs(w_rec_band[iBand]))*rho*(c0**2))/(pRef**2)) / np.max(((abs(w_rec_band[iBand]))*rho*(c0**2))/(pRef**2))) #normalised to maximum to 0dB
    spl_r_tot = 10*np.log10(rho*c0*((Ws/(4*math.pi*dist_sr**2))*np.exp(-m_atm*dist_sr) + ((abs(w_rec_band[iBand]))*c0))/(pRef**2)) #spl total (including direct field) at the receiver position????? but it will need to be calculated for a stationary source 100dB
    
    spl_rec_x_t0 = 10*np.log10(((abs(w_rec_x_t0))*rho*(c0**2))/(pRef**2))
    spl_rec_y_t0 = 10*np.log10(((abs(w_rec_y_t0))*rho*(c0**2))/(pRef**2))
    #Find the energy decay part of the overal calculation
    #idx_w_rec = np.argmin(np.abs(t - sourceon_time)) #index at which the t array is equal to the sourceon_time; I want the RT to calculate from when the source stops.
    #idx_w_rec = np.where(t == sourceon_time)[0][0] #index at which the t array is equal to the sourceon_time; I want the RT to calculate from when the source stops.
    #w_rec_off = w_rec_band[iBand][idx_w_rec:] #cutting the energy density array at the receiver from the idx_w_rec to the end
    
    #Schroeder integration
    #energy_r_rev = (w_rec_off)[::-1] #reverting the array
    #The energy density is related to the pressure with the following relation: w = p^2
    #energy_r_rev_cum = np.cumsum(energy_r_rev) #cumulative summation of all the item in the array
    schroeder = w_rec_off_band[iBand] #energy_r_rev_cum[::-1] #reverting the array again -> creating the schroder decay
    sch_db = 10.0 * np.log10(schroeder / max(schroeder)) #level of the array: schroeder decay
    
    if tcalc == "decay":
        t30 = t60_decay(t, sch_db, idx_w_rec, rt='t30') #called function for calculation of t60 [s]
        edt = t60_decay(t, sch_db, idx_w_rec, rt='edt') #called function for calculation of edt [s]
        #Eq_A = 0.16*V/t60 #equivalent absorption area defined from the RT 
        c80 = clarity(t30, V, Eq_A[iBand], S, c0, dist_sr) #called function for calculation of c80 [dB]
        d50 = definition(t30, V, Eq_A[iBand], S, c0, dist_sr) #called function for calculation of d50 [%]
        ts = centretime(t30, Eq_A[iBand], S) #called function for calculation of ts [ms]
        
        t30_band.append(t30)
        edt_band.append(edt)
        
    w_rec_x_band.append(w_rec_x_end)
    w_rec_y_band.append(w_rec_y_end)
    spl_stat_x_band.append(spl_stat_x)
    spl_stat_y_band.append(spl_stat_y)
    spl_r_band.append(spl_r)
    spl_r_off_band.append(spl_r_off)
    spl_r_norm_band.append(spl_r_norm)
    sch_db_band.append(sch_db)
    spl_t0_band.append(spl_t0)
    spl_new_band.append(spl_new)
    w_rec_x_t0_band.append(w_rec_x_t0)
    w_rec_y_t0_band.append(w_rec_y_t0)
    spl_rec_x_t0_band.append(spl_rec_x_t0)
    spl_rec_y_t0_band.append(spl_rec_y_t0)
    
    
spl_r_off_band = np.array(spl_r_off_band)
t30_band = np.array(t30_band)

et = time.time() #end time
elapsed_time = et - st

#%%
###############################################################################
#FIGURES & POST-PROCESSING
###############################################################################

if tcalc == "decay":
    for iBand in range(nBands):
        #Figure 5: Decay of SPL in the recording_time
        plt.figure(5)
        plt.plot(t, spl_r_band[iBand])  # plot sound pressure level with Pref = (2e-5)**5
        plt.title("Figure 5 :SPL over time at the receiver")
        plt.xlabel("t [s]")
        plt.ylabel("SPL [dB]")
        plt.xlim()
        plt.ylim()
        plt.xticks(np.arange(0, recording_time + 0.1, 0.5))
        plt.yticks(np.arange(0, 120, 20))
    
        #Figure 6: Decay of SPL in the recording_time normalised to maximum 0dB
        plt.figure(6)
        plt.plot(t,spl_r_norm_band[iBand])
        plt.title("Figure 6: Normalised SPL over time at the receiver")
        plt.xlabel("t [s]")
        plt.ylabel("SPL [dB]")
        plt.xlim()
        plt.ylim()
        plt.xticks(np.arange(0, recording_time +0.1, 0.1))
        plt.yticks(np.arange(0, -60, -10))
        
        #Figure 7: Energy density at the receiver over time
        plt.figure(7)
        plt.plot(t,w_rec_band[iBand])
        plt.title("Figure 7: Energy density over time at the receiver")
        plt.xlabel("t [s]")
        plt.ylabel("Energy density [kg m^-1 s^-2]")
        plt.xlim()
        plt.ylim()
        plt.xticks(np.arange(0, recording_time +0.1, 0.1))
        
        #Figure 8: Schroeder decay
        plt.figure(8)
        plt.plot(t[idx_w_rec:],sch_db_band[iBand])
        plt.title("Figure 8: Schroeder decay (Energy Decay Curve)")
        plt.xlabel("t [s]")
        plt.ylabel("Energy decay [dB]")
        plt.xlim()
        plt.ylim()
        plt.xticks(np.arange(t[idx_w_rec], recording_time +0.1, 0.1))
        
        #Figure 9: SPL at t=sourceoff_step over the space x. reverb sound only
        plt.figure(9)
        plt.title("Figure 9: SPL REVERB over the x axis at t=0")
        plt.plot(x_axis,spl_rec_x_t0_band[iBand])
        plt.ylabel('$\mathrm{SPL \ [dB]}$')
        plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$') 
        
        #Figure 10: SPL at t=sourceoff_step over the space y. reverb sound only
        plt.figure(10)
        plt.title("Figure 10: SPL REVERB over the y axis at t=0")
        plt.plot(y_axis,spl_rec_y_t0_band[iBand])
        plt.ylabel('$\mathrm{SPL \ [dB]}$')
        plt.xlabel('$\mathrm{Distance \ along \ y \ axis \ [m]}$')
        
        #Plot the heatmap
        grid_values = spl_t0_band[iBand][indices].reshape(grid_x.shape) #map the spl to the grid
        plt.figure()
        plt.contourf(grid_x, grid_y, grid_values, levels=50, cmap='viridis')
        plt.colorbar(label="Sound Pressure Level (SPL)")
        plt.title(f"Sound Pressure Level Heatmap at {target_height} m")
        plt.xlabel("X (m)")
        plt.ylabel("Y (m)")
        #plt.show()
    
if tcalc == "stationarysource":
    for iBand in range(nBands):

        #Figure 3: Decay of SPL in the recording_time at the receiver
        plt.figure(3)
        plt.plot(t,spl_r_band[iBand]) #plot sound pressure level with Pref = (2e-5)**5
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
        plt.plot(t,spl_r_norm_band[iBand])
        plt.xlabel("t [s]")
        plt.ylabel("SPL [dB]")
        plt.xlim()
        plt.ylim()
        plt.xticks(np.arange(0, recording_time +0.1, 0.1))
        plt.yticks(np.arange(0, -60, -10))
    
        #Figure 5: Energy density over time at the receiver
        plt.figure(5)
        plt.title("Figure 5: Energy density over time at the receiver")
        plt.plot(t,w_rec_band[iBand])
        plt.ylabel('$\mathrm{Energy \ Density \ [kg/ms^2]}$')
        plt.xlabel("t [s]")
        
        #Figure 6: Sound pressure level stationary over the space y.
        plt.figure(6)
        t_dim = len(t)
        last_time_index = t_dim-1
        #spl_y = spl_stat_y
        #data_y = spl_y
        plt.title("Figure 6: SPL over the y axis")
        plt.plot(y_axis,spl_stat_y_band[iBand])
        plt.xticks(np.arange(0, 20, 5))
        plt.yticks(np.arange(75, 105, 5))
        plt.ylabel('$\mathrm{Sound \ Pressure\ Level \ [dB]}$')
        plt.xlabel('$\mathrm{Distance \ along \ y \ axis \ [m]}$')
        
        #Figure 7: Sound pressure level stationary over the space x.
        plt.figure(7)
        t_dim = len(t)
        last_time_index = t_dim-1
        #spl_x = spl_stat_x
        #data_x = spl_x
        plt.title("Figure 7: SPL over the x axis")
        plt.plot(x_axis,spl_stat_x_band[iBand])
        #plt.xticks(np.arange(0, 35, 5))
        plt.yticks(np.arange(90, 97, 1))
        plt.ylabel('$\mathrm{Sound \ Pressure \ Level \ [dB]}$')
        plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$')
        
        #Figure 8: Energy density at t=recording_time over the space x.
        plt.figure(8)
        plt.title("Figure 8: Energy density over the x axis at t=recording_time")
        plt.plot(x_axis,w_rec_x_band[iBand])
        plt.ylabel('$\mathrm{Energy \ Density \ [kg/ms^2]}$')
        plt.xlabel('$\mathrm{Distance \ along \ x \ axis \ [m]}$')

#%%
###############################################################################
#SAVING
###############################################################################
np.save(os.path.join('dt'),dt)
np.save(os.path.join('w_rec_off_band'),w_rec_off_band)
np.save(os.path.join('w_rec_off_deriv_band'),w_rec_off_deriv_band)
np.save(os.path.join('p_rec_off_deriv_band'),p_rec_off_deriv_band)
np.save(os.path.join('t_off'),t_off)


# Save all variables to a file
def save(filename):
    with open(filename, 'wb') as f:
        # Filter out modules, functions, and other unsupported types
        filtered_variables = {}
        for k, v in globals().items():
            try:
                # Check if the object can be pickled
                pickle.dumps(v)
                # Exclude some types explicitly known to cause issues
                if not k.startswith('__') and not isinstance(v, (types.ModuleType, types.FunctionType, types.BuiltinFunctionType, types.LambdaType, types.MethodType, types.MappingProxyType)):
                    filtered_variables[k] = v
            except Exception as e:
                print(f"Could not pickle {k}: {str(e)}")

        pickle.dump(filtered_variables, f)

# To save all current variables
save('results.pkl')
