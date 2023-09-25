# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 16:26:36 2023

@author: 20225533
"""
import gmsh
import sys
from math import ceil
from math import log
from math import sqrt
import numpy as np

file_name = "rect.msh" #Insert file name, msh file created from sketchUp and then gmsh
gmsh.initialize() #Initialize msh file
mesh = gmsh.open(file_name) #open the file

c0= 343 #adiabatic speed of sound [m.s^-1]
th = 3

dim = -1 #dimensions of the entities, 0 for points, 1 for curves/edge/lines, 2 for surfaces, 3 for volumes, -1 for all the entities 
tag = -1 #all the nodes of the room

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

total_boundArea = 0
nodeTags, coords, parametricCoord = gmsh.model.mesh.getNodes(dim,tag) #gets the tags for each node and the coordinates of each node
nodecoords = coords.reshape((-1,3)) #coordinates reshaped in a matrix 3xnumber of nodes
velemNodes = elemNodeTags[2].reshape((-1,4)) #nodes per each tetrahedron 
belemNodes = elemNodeTags[1].reshape((-1,3)) #nodes per each surface boundary      


belement = 0
for elem in range(len(bounEl)):
    belement = belement + 1 #scalar number of the boundary elements

velement = 0
for elem in range(len(voluEl)):
    velement = velement + 1 #scalar number of the volume elements
 
######################################TXT###################################
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
            txt[t].append(int(tt-(voluEl[0]-1))) #volumes neighbours to each volume
for values in txt.values():
    if len(values)==2:
        values.append(0)
        values.append(0)
    if len(values) == 3:
        values.append(0)



########################neighbourVolume#######################################


neighbourVolume = np.array(())
for key in txt:
    neighbourVolume = np.append(neighbourVolume,txt[key])
for item in neighbourVolume:
        item = int(item)   
    
neighbourVolume = neighbourVolume.reshape((-1,4))

###############################################################################
#Absorption term
###############################################################################
#Absorption term for boundary conditions 
def abs_term(th,alpha):
    if th == 1:
        Absx = (c0*alpha)/4 #Sabine
    elif th == 2:
        Absx = (c0*(-log(1-alpha)))/4 #Eyring
    elif th == 3:
        Absx = (c0*alpha)/(2*(2-alpha)) #Modified by Xiang
    return Absx

Abs_terms = []
vGroups = gmsh.model.getPhysicalGroups(-1)
vGroupsNames = []
for iGroup in vGroups:
    dimGroup = iGroup[0]  #entity tag: 1 lines, 2 surfaces, 3 volumes (1D, 2D or 3D)
    tagGroup = iGroup[1]  #physical tag group 
    namGroup = gmsh.model.getPhysicalName(dimGroup, tagGroup)         
    alist = [dimGroup,tagGroup,namGroup]
    print(alist)
    vGroupsNames.append(alist)
    vEntities = gmsh.model.getEntitiesForPhysicalGroup(dimGroup, tagGroup)
    print(vEntities)

# Initialize a list to store surface tags and their absorption coefficients
surface_absorption = []
triangle_face_absorption = []

for group in vGroupsNames:
    if group[0] != 2:
        continue

    # Get the physical group tag
    physical_tag = group[1]

    # Retrieve all the entities in this physical group
    entities = gmsh.model.getEntitiesForPhysicalGroup(2, physical_tag)

    abscoeff = float(input(f"Enter absorption coefficient input for {group[2]}:"))
    Abs_term = abs_term(th, abscoeff)

    for entity in entities:
        surface_absorption.append((entity, Abs_term))
        # Get all the triangle faces for the current surface
        triangle_faces, _ = gmsh.model.mesh.getElementsByType(2, entity)
        
        # Append the Abs_term value for each triangle face
        triangle_face_absorption.extend([Abs_term] * len(triangle_faces))

###############################################################################
#Face areas
###############################################################################
       
total_boundArea = 0
boundaryV = []  # Initialize a list to store boundaryV values for each tetrahedron
import itertools
face_areas = np.zeros(len(velemNodes))
for idx, element in enumerate(velemNodes):
    tetrahedron_boundaryV = 0
    total_tetrahedron_boundaryV = 0
    print(idx)
    node_combinations = [list(nodes) for nodes in itertools.combinations(element, 3)]
    #nodes = element[:3]  # Take the first 3 nodes for a face
        # Check if the nodes are in any order in bounNode
    is_boundary = False
    for nodes in node_combinations:
        for surface in bounNode:
            surface_set = set(surface)
            nodes_set = set(nodes)
            if nodes_set.issubset(surface_set):
                is_boundary = True
                if is_boundary:
                    # Convert the vertices to NumPy arrays for vector operations
                    bc0 = gmsh.model.mesh.getNode(nodes[0])[0]
                    bc1 = gmsh.model.mesh.getNode(nodes[1])[0]
                    bc2 = gmsh.model.mesh.getNode(nodes[2])[0]
                    
                    # Calculate the area using half of the cross product's magnitude
                    face_area = 0.5 * np.linalg.norm(np.cross(bc1 - bc0, bc2 - bc0))
                    
                    face_areas[idx] = face_area
                    total_boundArea += face_area
                    
                    if face_area > 0:
                        #Calculate the index of the current face within the tetrahedron
                        face_index = np.where((bounNode == nodes))[0][0]
                        
                        # Use the index to access the corresponding absorption area
                        face_absorption_product = face_area * triangle_face_absorption[face_index]
                        
                        tetrahedron_boundaryV += face_absorption_product
                        
                    total_tetrahedron_boundaryV += tetrahedron_boundaryV
                        
    # Append the total boundaryV for the tetrahedron to the list
    boundaryV.append(total_tetrahedron_boundaryV)


######PREVIOUS########################################################
total_boundArea = 0

import itertools
face_areas = np.zeros(len(velemNodes))
for idx, element in enumerate(velemNodes):
    print(idx)
    node_combinations = [list(nodes) for nodes in itertools.combinations(element, 3)]
    #nodes = element[:3]  # Take the first 3 nodes for a face
        # Check if the nodes are in any order in bounNode
    is_boundary = False
    for nodes in node_combinations:
        for surface in bounNode:
            surface_set = set(surface)
            nodes_set = set(nodes)
            if nodes_set.issubset(surface_set):
                is_boundary = True
                if is_boundary:
                    # Convert the vertices to NumPy arrays for vector operations
                    bc0 = gmsh.model.mesh.getNode(nodes[0])[0]
                    bc1 = gmsh.model.mesh.getNode(nodes[1])[0]
                    bc2 = gmsh.model.mesh.getNode(nodes[2])[0]
                    
                    # Calculate the area using half of the cross product's magnitude
                    face_area = 0.5 * np.linalg.norm(np.cross(bc1 - bc0, bc2 - bc0))

                    face_areas[idx] = face_area
                    total_boundArea += face_area




###############################################################################
###############################################################################

faceNodes = gmsh.model.mesh.getElementFaceNodes(2, 3)
#faceTags, faceOrientations = gmsh.model.mesh.getFaces(3, faceNodes)
#elementTags, elementNodeTags = gmsh.model.mesh.getElementsByType(2)
#faces2Elements = {}
#for i in range(len(faceTags)): # 4 faces per tetrahedron
#    if not faceTags[i] in faces2Elements:
#        faces2Elements[faceTags[i]] = [elementTags[i // 4]]
#    else:
#        faces2Elements[faceTags[i]].append(elementTags[i // 4])
        
        
#for i in range(len(neighbourVolume)):
#    facebyvolume = gmsh.model.mesh.getFaces(3,neighbourVolume[i])
    
    
#boundaryAreaTotal = 0
#boundaryV = []
#Computing list of faces per tetrahedron
#faces = []
#facenodes2 = facenodes.reshape((-1,3)) #gmsh.model.mesh.getElementFaceNodes(4, 3) #4 is the element type (tetrahedron) and three are the nodes per each face

#faces_per_tet = [] #dictionary with keys as the nodes of each face and values the volume elements of which this face is neighbour

#for i in range(0, len(facenodes2)): # per ecah element basically, goes trhough the nodes of each face 3by3
#    print(i)
#    f = list(sorted(facenodes2[i])) #nodes of each face put in a tuple from node i to node i plus 3
#    #faces.append(f)
#    #t = voluEl[i // 12] #volume element number at which the faces are associated?
#    faces_per_tet.append(f)
#    #if not f in faces_per_tet: #if the face f (with its node) is alrady in the dictionary, just append the volume element neighbour to
    #    faces_per_tet[i] = [f]
    #else:
    #    faces_per_tet[i].append(f)

#faces_per_tet = np.reshape(faces_per_tet, (len(voluEl),4))

#total_faces = []
#list_facesxt = []
#for i in range(0,len(faces_per_tet),4):
#    list_facesxt.append(list(faces_per_tet[i:i+4]))

    
#Computing boundary area*abs_term
#boundaryAreaTotal = 0
#boundaryV = np.zeros(velement)
            
# Create a set of faces that are on the boundary for faster lookup

#for i in range(len(list_facesxt)):
    #print(i)
    #if i == 159:
    #    break
    #    is_boundary = False  # Flag to check if any face of this tetrahedron is on the boundary
#    for j in range(len(belemNodes)):
#        if (list_facesxt[i][0:4] == belemNodes[j]).any():
#            print(j)
#            print(belemNodes[j])
#            print(list_facesxt[i])
            #break
                #j = int(j)
                #face_in_consideration = gmsh.model.mesh.getFaces(3, j)
#            a_surf = face_areas[i] * triangle_face_absorption[j]
#    boundaryV[i] = a_surf




#face_per_tet = {}
#for i in range(len(neighbourVolume)):
#    print(i)
#new = []
#for key,value in fxt.items():
#flipped = {}
 
#for key, value in fxt.items():
#    if value[0] not in flipped:
#        flipped[value[0]] = [list(key)]
#    else:
#        flipped[value[0]].append(list(key))

    


#face_per_tet = {}
#for i in range(len(neighbourVolume)):
#    print(i)
#for key,value in fxt.items():
#    print(key)
#    a_list_key = list(key)
#    if len(fxt[key]) <2:
#        tet = fxt[key][0]
#        if tet not in face_per_tet:
#            face_per_tet[tet] = a_list_key
#        else:
#            face_per_tet[tet].append(a_list_key)
    
"""

boundaryAreaTotal = 0
boundaryV = []
            
# Create a set of faces that are on the boundary for faster lookup
boundary_faces = set(bounEl)

for i in range(len(neighbourVolume)):
    print(i)
    is_boundary = False  # Flag to check if any face of this tetrahedron is on the boundary
    for j in bounEl:
        print(j)
        if j in range(len(boundary_faces)):
            j = int(j)
            a_surf = face_areas[i] * triangle_face_absorption[j]
            boundaryV.append(a_surf)
            is_boundary = True
            break  # No need to check other faces if one is on the boundary
    if not is_boundary:
        boundaryV.append(0)
        
"""