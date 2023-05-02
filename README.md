## Diffusion Equation Code for Room Acoustics Modeling
Code is developed by Ilaria Fichera at Eindhoven University of Technology.

The diffusion equation, generally used in heat conduction, can be used also for acoustics purposes.
The acoustic modelling method of the diffusion equation allows to study the acoustics properties of a room and to obtain spatial distribution of acoustic energies over time in a specific room.
The code is a currently a **PROTOTYPE**. The code is written in Python language.

## Release version
Version 1.0

## Documentation
Documentation documents can be found below:
- [Documentation 1D](https://ilariafichera.github.io/Diffusion/Documentation1D.html)
- [Documentation 2D](https://ilariafichera.github.io/Diffusion/Documentation2D.html)
- [Documentation 3D](https://ilariafichera.github.io/Diffusion/Documentation3D.html)

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
Code is developed by Ilaria Fichera at Eindhoven University of Technology.

## References
+ J. M. Navarro, J. Escolano and J. J. Lopez, Implementation and evaluation of a diffusion equa-tion model based on finite difference schemes for sound field prediction in rooms,Appl. Acoust.73(2012) 659–665.
+ Billon A, Picaut J, Foy C, Valeau V, Sakout A. Introducing atmospheric attenuation within a diffusion model for room-acoustic predictions. J Acoust Soc Am. 2008 Jun;123(6):4040-3. doi: 10.1121/1.2903872. PMID: 18537354.
+ Vorländer M. Auralization: fundamentals of acoustics, modelling, simulation,algorithms and acoustic virtual reality. Springer 2008

## Fundings
The project is funded by NWO, in the Netherlands.

## Licence
Diffusion Code is licenced under MIT licence. See LICENCE.md for more details.
