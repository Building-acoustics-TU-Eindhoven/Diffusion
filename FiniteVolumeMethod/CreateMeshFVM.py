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

gmsh.initialize()
gmsh.open('8x8x8.geo')
gmsh.option.setNumber('Mesh.MeshSizeMin', 1)
gmsh.option.setNumber('Mesh.MeshSizeMax', 1)
gmsh.model.mesh.generate(3)
gmsh.write('mesh.msh')
gmsh.finalize()
