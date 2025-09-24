# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 08:56:01 2023

@author: Ilaria Fichera
"""

# Code developed by Ilaria Fichera for the analysis of the FVM method adapted solving the 3D diffusion equation with one intermittent omnidirectional sound source
# Import modules

import gmsh
import argparse
import os

# %%
###############################################################################
# INITIALISE AND GENERATE GMSH
###############################################################################

def generate_mesh(geo_file_path, name_gmsh_file, characteristic_length):
    
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

    # If an lc is given in the geo file, we want to compensate for this
    lc_value = 1 # set to 1 by default
    for line in geo_content:
        if "lc =" in line:
            lc_value = float(line.split('=')[1].strip().strip(';'))
            print("Extracted value:", lc_value)
            break
    gmsh.initialize()
    gmsh.open(geo_file_path)

    # Divide by the lc value given in the .geo file to get a correct characteristic length 
    gmsh.option.setNumber('Mesh.MeshSizeFactor', characteristic_length / lc_value)

    gmsh.model.mesh.generate(3)
    print(gmsh.logger.get())
    gmsh.write(name_gmsh_file)
    
    gmsh.finalize()
    # mesh = gmsh.open(name_gmsh_file)  # open the file

    # gmsh.fltk.run()  # run the file to see it in gmsh


# %%
###############################################################################
# INPUT VARIABLES
###############################################################################
# The variables to be assigned by the user directly
# Only name of the file
name_of_geo_file = '3x3x3.geo'
name_of_gmsh_file = "3x3x3.msh"
length_of_mesh = 1


###############################################################################

base_dir = os.path.dirname(os.path.abspath(__file__))  # Directory of the current script

def main():
    #gmsh.initialize()
    parser = argparse.ArgumentParser(description='process the input path, taking dynamically.')

    parser.add_argument(
        'geo_file_path',
        type=str,
        help='Path to the original Geo file',
        nargs='?',  # This makes the argument optional
        default=os.path.join(
            base_dir, name_of_geo_file
        ),
    )
    parser.add_argument(
        'name_gmsh_file',
        type=str,
        nargs='?',  # This makes the argument optional
        help='Name for the output Gmsh file',
        default=os.path.join(
            base_dir, name_of_gmsh_file
        )
    )
    parser.add_argument(
        'length_of_mesh',
        type=int,
        nargs='?',  # This makes the argument optional
        help='Length of the mesh',
        default=length_of_mesh
    )

    args = parser.parse_args()
    generate_mesh(args.geo_file_path, args.name_gmsh_file, args.length_of_mesh)
    #gmsh.finalize()

if __name__ == '__main__':
    main()

