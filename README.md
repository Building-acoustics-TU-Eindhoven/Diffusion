## Diffusion Equation Code for Room Acoustics Modeling
Code is developed by Ilaria Fichera at Eindhoven University of Technology.

The diffusion equation, generally used in heat conduction, can be used also for acoustics purposes.
The acoustic modelling method of the diffusion equation allows to study the acoustics properties of a room and to obtain spatial distribution of acoustic energies over time in a specific room.
The code is a currently a **PROTOTYPE**. The code is written in Python language.

## Release version

## Documentation

### How to run the code
The code is based on the finite different method Du Fort&Frankel.
It can be used for cube/rectangular parallelepiped shapes.
The code is developed for a 1D, 2D and 3D situation, respectively DiffEq1D, DiffEq2D, DiffEq3D.

The "DiffEq3D" is the main code and it works with the following associated functions:
+ The Function RT calculate the reverberation time of the room in question
+ The Function Clarity calculate the clarity80ms of the room in question based on Barron's revised theory formulas.
+ The Function Definition calculate the definition50ms of the room in question based on Barron's revised theory formulas.
+ The Function Centre Time calculate the centre time of the room in question based on Barron's revised theory formulas.

## Authors

## Publications

### References
+ ....
+ ....
+ ....
+ ....

## Fundings

### Licence
Diffusion Code is licenced under MIT licence. See LICENCE.md for more details.
