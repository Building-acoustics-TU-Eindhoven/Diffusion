DIFFUSION EQUATION CODE FOR ROOM ACOUSTICS MODELING
Code is developed by Ilaria Fichera with Eindhoven University of Technology.

The diffusion equation, generally used in heat conduction, can be used also for acoustics purposes.
The acoustic modelling method of the diffusion equation allows to study the acoustics properties of a room and to obtain spatial distribution of acoustic energies over time in a specific room. 
The code is a currently a PROTOTYPE. The code is written in Python language.


HOW TO RUN THE CODE
The code is based on the finite different method Du Fort&Frankel.
It can be used for cube/rectangular parallelepiped shapes.
The code is developed for a 1D, 2D and 3D situation, respectively DiffEq1D, DiffEq2D, DiffEq3D.
THe "DiffEq3D" is the main code and it works with the following associated functions: FunctionRT, FunctionDefinition, FunctionClarity, FunctionCentreTime.
The Function RT calculate the reverberation time of the room in question
The Function Clarity calculate the clarity80ms of the room in question based on Barron's revised theory formulas.
The Function Definition calculate the definition50ms of the room in question based on Barron's revised theory formulas.
The Function Centre Time calculate the centre time of the room in question based on Barron's revised theory formulas.


REFERENCES
 - 
 -
 -


LICENCE
Diffusion Code is licenced under MIT licence. See LICENCE.md for more details.