# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 08:56:01 2023

@author: 20225533
"""

#Code developed by Ilaria Fichera for the analysis of the FVM method adapted solving the 3D diffusion equation with one intermittent omnidirectional sound source
#Import modules

#import py2gmsh as p2g
#from py2gmsh import (Mesh, Entity, Field)

# Initialize py2gmsh
#geom = p2g.Geometry()

#Load the Geo file
#geom.load_file("8x8x8.geo")

#Set mesh sizes
#lc = 0.1  # Set your desired mesh size here
#geom.set_element_size(lc)

# Generate the mesh
#geom.generate_mesh()

# Save the mesh to an MSH file
#geom.write_mesh("output.msh")

import gmsh
import pygmsh as pg
import py2gmsh as p2g
from py2gmsh import (Mesh, Entity, Field)

import subprocess

#Set your desired mesh size
#lc = 10

#Define the path to your Gmsh executable
#gmsh_path = "C:\\Users\\20225533\\Diffusion\\FiniteVolumeMethod\\gmsh.exe"

#Define the path to your Geo file
#geo_file = "C:\\Users\\20225533\\Diffusion\\FiniteVolumeMethod\\8x8x8.geo"

#gmsh.model.mesh.generate()
#8x8x8.geo -3 -clmax lc -o my_mesh.msh

#mesh = geo_file.generate_mesh(dim=2,algorithm = 6)

#Define the path for the output mesh file (MSH format)
#output_mesh_file = "C:\\Users\\20225533\\Diffusion\\FiniteVolumeMethod\\output.msh"

#Construct the Gmsh command to generate the mesh
#command = [gmsh_path, geo_file, "-3", "-clmax", str(lc), "-o", output_mesh_file] #path, file, -3 = 3D mesh, clmax=caracteristic length maximum

#Run Gmsh to generate the mesh
#subprocess.run(command)



#%%
###############################################################################
#INPUT VARIABLES
###############################################################################
# Specify the path to your original Geo file
geo_file_path = '5x5x5.geo'
max_mesh_size = 1
name_gmsh_file = '5x5x5.msh'
length_of_mesh = 0.5

#%%
#import os
#import subprocess

# Assuming lengthOfMesh is a variable in Python
# and filename is the name of the file without extension
#if lengthOfMesh == 'H':
#    strimp = f'gmsh -3 -format msh4 {os.getcwd()}/"8x8x8".geo'
#    subprocess.run(strimp, shell=True)
#    lengthOfMesh = 2
#else:

#strimp = f'gmsh -3 -format msh4 -clscale {lengthOfMesh} "{os.getcwd()}/"8x8x8".geo"'
#subprocess.run(strimp, shell=True)
#lengthOfMesh = float(lengthOfMesh)

# Assuming msh_file is the path to the generated mesh file
#msh_file = f'{os.getcwd()}/GeoModels/"8x8x8".msh'


#%%
###############################################################################
#INITIALISE GMSH
###############################################################################
# Read the content of the Geo file
with open(geo_file_path, 'r') as file:
    geo_content = file.readlines()

# Specify the line to be removed
line_to_remove = 'Mesh.RemeshAlgorithm = 1; // automatic\n'

# Remove the specified line from the content
if line_to_remove in geo_content:
    geo_content.remove(line_to_remove)

# Write the modified content back to the Geo file
with open(geo_file_path, 'w') as file:
    file.writelines(geo_content)

gmsh.initialize()
gmsh.open(geo_file_path)

#gmsh.option.setNumber('Mesh.MeshSizeMin', 1)
#gmsh.option.setNumber('Mesh.MeshSizeMax', max_mesh_size)

# Set the mesh size factor to control the characteristic length
gmsh.option.setNumber('Mesh.MeshSizeFactor', length_of_mesh) 

gmsh.model.mesh.generate(3)
print(gmsh.logger.get())
gmsh.write(name_gmsh_file)
gmsh.finalize()

gmsh.initialize() #Initialize msh file
mesh = gmsh.open(name_gmsh_file) #open the file

#gmsh.fltk.run() #run the file to see it in gmsh
