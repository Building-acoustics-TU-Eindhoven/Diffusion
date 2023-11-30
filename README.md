## Diffusion Equation Software for Room Acoustics Modeling
The software implementation is part of an ongoing research in the Building Acoustics Group at the Built Environment Department of Eindhoven University of Technology.
The software is currectly **UNDER DEVELOPMENT** and it is being implemented by Ilaria Fichera in Python language.

The goal of the project is to specifically implement, develop and validate the diffusion equation modelling method.
Solving the diffusion equation allows to study the acoustics properties of a room and to obtain the distribution of acoustic energy over space and time in a specific room.

## Release version
Version 1.0

## Setup and Usage instructions
1. Download and install Anaconda or download and install any Python software IDE you prefer
2. Clone/Fork this repository to a folder of your preference
3. Open the Main files through the preferred IDE and test the software

### Repository structure
The solution of the Diffusion Equation is currently investigated with two different numerical methods: the Finite Different Method (FDM) by Du Fort&Frankel (Navarro et al., 2012) and the Finite VOlume Method (FVM) (Munoz, 2019). 
The FDM can be used for cuboid shapes while the FVM can be used for any shape/geometry.
The FDM script is developed for a 1D, 2D and 3D spaces, respectively reflected in the following Python files DiffEq1D.py, DiffEq2D.py, DiffEq3D.py.
The FVM script is developed for 3D spaces in the Python file FVM.py. 

For both the numerical methods, the main software works with the following associated functions:
+ FunctionRT.py calculates the reverberation time of the room in question
+ FunctionEDT.py calculates the early decay time of the room in question
+ FunctionClarity.py calculates the clarity $C_{80}$ of the room in question based on Barron's revised theory formula.
+ FunctionDefinition.py calculates the definition $D_{50}$ of the room in question based on Barron's revised theory formula.
+ FunctionCentreTime.py calculates the centre time $T_s$ of the room in question based on Barron's revised theory formula.

## Documentation
Documentation documents can be found below:
- [Documentation 1D](https://ilariafichera.github.io/Diffusion/Documentation1D.html)
- [Documentation 2D](https://ilariafichera.github.io/Diffusion/Documentation2D.html)
- [Documentation 3D](https://ilariafichera.github.io/Diffusion/Documentation3D.html)
- [Documentation FVM](https://ilariafichera.github.io/Diffusion/DocumentationFVM.html)

## Sample Calculation FDM

Test the software with the following inputs (Navarro et al., 2012):
- Spatial discretization:             dx = 0.5 m
- Time discretization:                dt = 1/8000 s
- Recording total time:               T = 2.0 s
- Room dimension:                     8.0 m x 8.0 m x 8.0 m 
- Absorption Coefficient surface i:   1/6 ≈ 0.17
- Air absorption:                     0.0 1/m
- Absorption term:                    $A_{M}$
- Source On time:                     S = 1.0 s 
- Source position coordinates:        4.0, 4.0, 4.0 
- Receiver position coordinates:      2.0, 2.0, 2.0

Test if the software provides the following results (Navarro et al., 2012):
- Reverberation time (RT):            1.18 s
- Early Decay Time (RDT):             1.18 s
- Clarity ($C_{80}$):                 2.08 dB
- Definition ($D_{50}$):              45.63 %
- Centre Time ($T_s$):                85.33 ms

## Authors
Software is being developed by Ilaria Fichera at Eindhoven University of Technology.

## References
+ J. M. Navarro, J. Escolano and J. J. Lopez, Implementation and evaluation of a diffusion equation model based on finite difference schemes for sound field prediction in rooms, Appl. Acoust.73 (2012) 659–665.
+ R. P. Muñoz, "Numerical modeling for urban sound propagation: developments in wave-based and energy based methods," PhD Thesis, Technische Universiteit Eindhoven, 2019.

## Funding
The project is funded by <u>[NWO](https://www.nwo.nl/projecten/19430), in the Netherlands.

## Licence
Diffusion Code is licenced under GNU General Public License v2.0. See LICENCE.md for more details.
