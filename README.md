## Diffusion Equation Code for Room Acoustics Modeling
Code has been developed by Ilaria Fichera at Eindhoven University of Technology.

The diffusion equation, generally used for heat conduction problems, can also be used for acoustics purposes.
Solving the diffusion equation allows to study the acoustics properties of a room and to obtain the distribution of acoustic energy over space and time in a specific room.
The code is a currently a **PROTOTYPE**. The code is written in Python language.

## Release version
Version 1.0

## Documentation
Documentation documents can be found below:
- [Documentation 1D](https://ilariafichera.github.io/Diffusion/Documentation1D.html)
- [Documentation 2D](https://ilariafichera.github.io/Diffusion/Documentation2D.html)
- [Documentation 3D](https://ilariafichera.github.io/Diffusion/Documentation3D.html)

### How to run the code
The numerical solver of the Diffusion Equation is based on the finite differences method by Du Fort&Frankel.
It can be used for cuboid shapes.
The code is developed for a 1D, 2D and 3D spaces, respectively reflected in the following Python codes DiffEq1D.py, DiffEq2D.py, DiffEq3D.py

DiffEq3D.py is the main code and it works with the following associated functions:
+ FunctionRT.py calculates the reverberation time of the room in question
+ FunctionClarity.py calculates the clarity $C_{80}$ of the room in question based on Barron's revised theory formula.
+ FunctionDefinition.py calculates the definition $D_{50}$ of the room in question based on Barron's revised theory formula.
+ FunctionCentreTime.py calculates the centre time $T_s$ of the room in question based on Barron's revised theory formula.

## Authors
Code has been developed by Ilaria Fichera at Eindhoven University of Technology.

## References
+ J. M. Navarro, J. Escolano and J. J. Lopez, Implementation and evaluation of a diffusion equa-tion model based on finite difference schemes for sound field prediction in rooms,Appl. Acoust.73(2012) 659â€“665.
+ Billon A, Picaut J, Foy C, Valeau V, Sakout A. Introducing atmospheric attenuation within a diffusion model for room-acoustic predictions. J Acoust Soc Am. 2008 Jun;123(6):4040-3. doi: 10.1121/1.2903872. PMID: 18537354.
+ Vorländer M. Auralization: fundamentals of acoustics, modelling, simulation,algorithms and acoustic virtual reality. Springer 2008

## Fundings
The project is funded by <u>[NWO](https://www.nwo.nl/projecten/19430), in the Netherlands.

## Licence
Diffusion Code is licenced under MIT licence. See LICENCE.md for more details.
