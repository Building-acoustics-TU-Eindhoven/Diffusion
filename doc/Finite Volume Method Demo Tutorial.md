# Finite Volume Method Demo Tutorial

After installing all the softwares/libraries and reading the documentation in [Finite Volume Method Use Documentation](https://building-acoustics-tu-eindhoven.github.io/Diffusion/Finite%20Volume%20Method%20Use.html), the software can be tested with the following files:
- _3x3x3.skp_
- _3x3x3.geo_
- _3x3x3.msh_

The file _3x3x3.skp_ is a SketchUp file with a room of volume 3x3x3 $m^3$. The SketchUp file has been already created following the steps in the [Finite Volume Method Use Documentation, Paragraph 'Geometry'](https://building-acoustics-tu-eindhoven.github.io/Diffusion/Finite%20Volume%20Method%20Use.html#geometry).

The file _3x3x3.geo_ is the .geo file created and saved according to the instructions in [Finite Volume Method Use Documentation, Paragraph 'Geometry'](https://building-acoustics-tu-eindhoven.github.io/Diffusion/Finite%20Volume%20Method%20Use.html#geometry).

The file _3x3x3.msh_ is the mesh file created by running the python script _CreateMeshFVM.py_ following the instructions in [Finite Volume Method Use Documentation, Paragraph 'Mesh creation'](https://building-acoustics-tu-eindhoven.github.io/Diffusion/Finite%20Volume%20Method%20Use.html#mesh-creation).

Currently in the repository, the following files _3x3x3.skp_, _3x3x3.geo_, _3x3x3.msh_ are already created.

Therefore, to test the software:
1. Open the script _FVM.py_;
2. Input the mesh file in the "file_name" variable;
3. Input all the other parameters in the input section based on [Finite Volume Method Use Documentation, Paragraph 'Inputs'](https://building-acoustics-tu-eindhoven.github.io/Diffusion/Finite%20Volume%20Method%20Use.html#inputs)
4. Press "Run";
5. During the calculation (quite soon after you pressed "Run"), the python script will prompt you with a question regarding the absorption coefficients of the surface of the room. It is important to include the absorption coefficient per each frequency separated by commas as per this example:
"`Enter absorption coefficient for frequency {fc_low} to {fc_high} for Wall1:`0.15,0.36,0.64,0.83,0.88".
6. Once you have done that for all the different surface material types, the software should continue running without any prompt. The simulation time depends on the room dimensions but mostly on the mesh dimension.

The software should provide results of Sound Pressure Levels at the receiver position, Reverberation time, Clarity and other energetic parameters in a pickle file called resultsFVM.pkl.