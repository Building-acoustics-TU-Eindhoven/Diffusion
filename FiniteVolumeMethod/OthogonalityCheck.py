# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 13:00:03 2024

@author: 20225533
"""

import numpy as np
import gmsh

file_name = "scenario2.msh" #Insert file name, msh file created from sketchUp and then gmsh
gmsh.initialize() #Initialize msh file
mesh = gmsh.open(file_name) #open the file

dim = -1 #dimensions of the entities, 0 for points, 1 for curves/edge/lines, 2 for surfaces, 3 for volumes, -1 for all the entities 
tag = -1 #all the nodes of the room

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
    
    
interior_tet = np.zeros((velement, velement)) #initialization matrix of tetrahedron per tetrahedron
angles_ortho = np.zeros((velement, velement))
shared_area_vector = np.zeros((3))
angles = np.zeros((velement))

for i in range(velement): #for each tetrahedron, take its centre
    print(i)
    cell_center_i = cell_center[i]
    angle_non_ortho = []
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
                face_centre = (sc0+sc1+sc2)/3
                face_to_cell = face_centre - cell_center_i
                shared_area_vector[0] = (sc2[1]-sc0[1])*(sc1[2]-sc0[2]) - (sc1[1]-sc0[1])*(sc2[2]-sc0[2])
                shared_area_vector[1] = (sc2[2]-sc0[2])*(sc1[0]-sc0[0]) - (sc1[2]-sc0[2])*(sc2[0]-sc0[0])
                shared_area_vector[2] = (sc2[0]-sc0[0])*(sc1[1]-sc0[1]) - (sc1[0]-sc0[0])*(sc2[1]-sc0[1])
                shared_area = np.linalg.norm(np.cross(sc2-sc0,sc1-sc0))/2 #compute shared area
                Area_vector = (shared_area*shared_area_vector)/np.linalg.norm(shared_area_vector) 
                #The "area vector" here is a mathematical construct used to 
                #represent the signed area of a triangle in three-dimensional space. 
                #It's not a traditional vector in the sense of a direction and 
                #magnitude, but rather a mathematical abstraction that encodes 
                #information about the orientation and magnitude of the area.
                shared_distance_vector = cell_center[j] - cell_center_i #vector distance between the two cell centre of the tetrahedrons
                shared_distance_magnitude = np.linalg.norm(cell_center_i - cell_center[j]) #magnitude of distance between the two cell centre of the tetrahedrons
                fv1_n = sc1 - sc0 #vector1 for normal of the shared face
                fv2_n = sc2 - sc0 #vector2 for normal of the shared face
                f_normal = np.cross(fv1_n, fv2_n) #normal vector of shared face
                
                #Max of the arcocosine #????
                #ang = max((np.degrees(np.arccos(np.dot(face_to_cell, f_normal) / (np.linalg.norm(face_to_cell) * np.linalg.norm(f_normal)))),
                #    np.degrees(np.arccos(np.dot(shared_distance_vector, f_normal) / (np.linalg.norm(shared_distance_vector) * np.linalg.norm(f_normal))))))
                
                #ang = min((np.dot(face_to_cell, f_norm) / (np.linalg.norm(face_to_cell) * np.linalg.norm(f_norm))),
                #    (np.dot(shared_distance_vector, f_norm) / (np.linalg.norm(shared_distance_vector) * np.linalg.norm(f_norm))))
                
                #Based on ANSYN
                # #This is with the area of the face/vector face to cell and/or area of the face/vector distance between two cells
                # This is an orthogonality inex for all the elements
                ang = min(abs(np.dot(Area_vector, face_to_cell) / (np.linalg.norm(Area_vector) * np.linalg.norm(face_to_cell))),
                   abs(np.dot(Area_vector, shared_distance_vector) / (np.linalg.norm(Area_vector) * np.linalg.norm(shared_distance_vector))))
                
                interior_tet[i, j] = shared_area/shared_distance_magnitude #division between shared area and shared distance
                angle_non_ortho.append(ang)
                if shared_area == 0:
                    angles_ortho[i,j] = 0 
                else:
                    angles_ortho[i,j] = ang
            else:
                shared_area = 0
                interior_tet[i, j] = shared_area
    angles[i] = min(angle_non_ortho)
    
#Normalise if it is an angle
#normalised_data = (angles - min(angles))/(max(angles) - min(angles)) #????



























# def coms_line_dist(volume_coords, face_nodes):
#     b_c_1 = np.mean([volume_coords[node - 1] for node in face_nodes], axis=0)
#     b_c_2 = np.mean([volume_coords[node - 1] for node in face_nodes], axis=0)
#     coms_line = b_c_2 - b_c_1
#     center_distance = np.linalg.norm(coms_line)
#     return coms_line, center_distance

# def face_normal_dir(face_coords):
#     v1 = face_coords[1] - face_coords[0]
#     v2 = face_coords[2] - face_coords[1]
#     f_norm = np.cross(v1, v2)
#     return f_norm

# def check_non_orth(coms_line, f_norm, avNonOrth, non_orthogonal_threshold):
#     ang = np.arccos(np.dot(coms_line, f_norm) / (np.linalg.norm(coms_line) * np.linalg.norm(f_norm)))
#     avNonOrth += ang
#     maxNonOrth = 0
#     if ang > non_orthogonal_threshold:
#         maxNonOrth = max(maxNonOrth, ang)
#     return avNonOrth, maxNonOrth

# def com_of_face(face_coords):
#     return np.mean(face_coords, axis=0)

# def check_skewness(p_f_com, coms_line, center_distance, avSkew, skewness_threshold):
#     skewness = np.linalg.norm(np.cross(p_f_com - coms_line[0], p_f_com - coms_line[1])) / center_distance
#     avSkew += skewness
#     maxSkew = 0
#     if skewness > skewness_threshold:
#         maxSkew = max(maxSkew, skewness)
#     return avSkew, maxSkew

# def calc_averages(avNonOrth, avSkew, total_faces):
#     avNonOrth /= total_faces if total_faces > 0 else 1
#     avSkew /= total_faces if total_faces > 0 else 1
#     return avNonOrth, avSkew

# def print_stats(avNonOrth, maxNonOrth, avSkew, maxSkew):
#     print("MeshQualityCheck:")
#     print("Average non-orthogonality:", avNonOrth)
#     print("Maximum non-orthogonality:", maxNonOrth)
#     print("Average skewness:", avSkew)
#     print("Maximum skewness:", maxSkew)


# # Example usage
# file_name = "3x3x3.msh" #Insert file name, msh file created from sketchUp and then gmsh

# non_orthogonal_threshold=65
# skewness_threshold=0.5

# gmsh.initialize()
# gmsh.open(file_name)

# volume_tags, volume_coords, _ = gmsh.model.mesh.getNodes(dim=3)

# avNonOrth = 0
# maxNonOrth = 0
# avSkew = 0
# maxSkew = 0
# total_faces = 0

# for vol_tag in volume_tags:
#     volume_faces = gmsh.model.getBoundary([(3, vol_tag)])[0]
#     for face_tag in volume_faces:
#         total_faces += 1
#         face_nodes = gmsh.model.mesh.getElements(dim=2, tag=face_tag)[0][1]
#         face_coords = np.array([volume_coords[node - 1] for node in face_nodes])
#         coms_line, center_distance = coms_line_dist(volume_coords, face_nodes)

#         f_norm = face_normal_dir(face_coords)
#         avNonOrth, maxNonOrth = check_non_orth(coms_line, f_norm, avNonOrth, non_orthogonal_threshold)

#         p_f_com = com_of_face(face_coords)
#         avSkew, maxSkew = check_skewness(p_f_com, coms_line, center_distance, avSkew, skewness_threshold)

# avNonOrth, avSkew = calc_averages(avNonOrth, avSkew, total_faces)
# print_stats(avNonOrth, maxNonOrth, avSkew, maxSkew)

# gmsh.finalize()