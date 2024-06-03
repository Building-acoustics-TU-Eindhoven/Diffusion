# Finite Volume Method Tutorial

After installing all the softwares and library and reading the documentation in [Finite Volume Method Use Documentation](https://building-acoustics-tu-eindhoven.github.io/Diffusion/Finite%20Volume%20Method%20Use.html), the software can be tested with the following files:
- 3x3x3.skp
- 3x3x3.geo
- 3x3x3.msh

The file "3x3x3.skp" is a SketchUp file with a room of volume 3x3x3 m^3. The SketchUp file has been already created following the steps in the [Finite Volume Method Use Documentation, Paragraph 'Geometry'](https://building-acoustics-tu-eindhoven.github.io/Diffusion/Finite%20Volume%20Method%20Use.html#geometry).

The file "3x3x3.geo" is the .geo file created and saved according to the instructions in [Finite Volume Method Use Documentation, Paragraph 'Geometry'](https://building-acoustics-tu-eindhoven.github.io/Diffusion/Finite%20Volume%20Method%20Use.html#geometry).

The file "3x3x3.msh" is the mesh file created by running the python script "CreateMeshFVM.py" following the instructions in [Finite Volume Method Use Documentation, Paragraph 'Mesh creation'](https://building-acoustics-tu-eindhoven.github.io/Diffusion/Finite%20Volume%20Method%20Use.html#mesh-creation).

Currently in the repository, the following files "3x3x3.skp", "3x3x3.geo", "3x3x3.msh" are already created.

Therefore, the only thing to do is to test the software by opening the script "FVM.py", inputting the mesh file in the "file_name" variable and press "Run". 

During the calculation (quite soon after you pressed "Run"), the python script will prompt you with a question regarding the absorption coefficients of the surface of the room. It is important to include the absorption coefficient per each frequency (maximum 5 frequencies) separated by commas as per this example:
"`Enter absorption coefficient for frequency {fc_low} to {fc_high} for Wall1:`0.15,0.36,0.64,0.83,0.88".

Once you have done that for all the different surface material types, the software should continue running without any prompt. The simulation time depends on the room dimensions but mostly on the mesh dimension.

The software should provide results of Sound Pressure Levels at the receiver position, Reverberation time, Clarity and other energetic parameters.
