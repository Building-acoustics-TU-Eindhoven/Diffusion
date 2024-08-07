# Finite Different Method Use

## Requirements
1. Download and install Anaconda or download and install any Python software IDE you prefer;
2. Clone/Fork this repository to a folder of your preference;
3. Open the "DiffEq3D.py" file through the preferred IDE;
4. After reading the "Libraries" and "Input" sections below, test the software with your preferred inputs.

## Libraries
To properly run the software, the following libraries are needed:
- Python version 3.10.9 or above

Libraries for python:
- math
- matplotlib
- numpy
- scipy
- sys
- drawnow
- time

## Python running files
Currently, in the repository there are a lot of python files. This is because the software is still in development. 
Currently, the updated files for users to run the software are:
- DiffEq1D.py: to simulate a line geometry;
- DiffEq2D.py: to simulate a surface geometry;
- DiffEq3D.py: to simulate a 3D geometry.

## Algorithm
The software is organised in three sections:
- Input variables:
    The inputs regarding the room dimensions, source and receiver positions along with other are to be inserted for the specific room in question.
- Calculation loop:
    The "for loop" would loop over the time to calculate the energy density at each position in the mesh grid and at each time step.
- Results and Post-processing
    Results of SPL and other are included in this section together with graphs for the analysis.

## Inputs

### Geometry
The geometry for this method is defined within the inputs variable section of the main python file "DiffEq3D.py".
The method is suitable for parallelepiped and cuboid spaces and therefore, the only inputs are the length, width and height of the room.

### Surface materials
The surface materials are defined within the inputs variable section of the main python file "DiffEq3D.py".
Currently, the surface material would require only one value of absorption coefficient per surface.
It is not possible to include doors, windows or other items within a surface. 

### Sound source
The model allows for the insertion of only one source position per calculation. 
The sound source is defined as an omnidirectional source. The users can input the sound power of the source *W_{s}* in Watts and its position in the room in the following variables *x_\{source\}*,*y_\{source\}*,*z_\{source\}* in m in the x,y,z directions.

The source can be defined as an interrupted noise source or an impulse source. 
- If an interrupted noise source is chosen, then the time within which the source stays on before getting switch off need to be defined. The variable is *sourceon_\{time\}*. The time is defined in seconds and it will need to be long enough for the room to be filled with sound before switching it off. The source on time needs to be inputted as the 2/3 of the theoretical calculated Sabine reverberation time of the room.
In this case, there is the option to to toggle the sound source on and off or keep it stationary. This can be done by changing the variable "tcalc" from "decay" to "stationarysource".
- If an impulse source is chosen, the source will be automatically defined.

### Receiver
The model allows for the insertion of only one acoustics receiver position per calculation. These are defined as point omnidirectional receivers. The users will input the position of the receiver in the room in the following vaariables *x_\{rec\}*,*y_\{rec\}*,*z_\{rec\}* in m in the x,y,z directions.

### Discretization variables
#### Spatial discretization dv
The Finite Different method works with a spatial discretization. The space is defined by a mesh grid of points at a distance dv between each other. The distance dv is equal for each dimension x,y,z, therefore dv = dx = dy = dz, and it is defined in meters. 
It is important to choose an appropriate dv for the precision of the calculation. 
A dv of 0.5 m normally would suffice for a correct calculation, although the choice is contingent upon user preferences for details in parameters values and room dimensions. It is strongly advised that dv is a factor of the dimensions of the room (e.g. for a room of 3x3x3, dv could be 0.5 m or 0.2 m or 0.1 m but it should not be 0.4 m). 

#### Time discretization dt
According the Navarro 2012, to get good converged results, the time discretization dt will need to be defined depending on the dv chosen. 
To make sure that the predictions converge to a fixed value with a very low error, the following empirical cretirion will need to apply.
```{math}
10^{-8} = dt^2 dv^{-2}
```
The time discretization is defined in seconds. 

### Total Recording time
The total recording time is the amount of time for the calculation to run. It is the sum of the source on time and the time for the decay. 
The decay time needs to be inputted as the 2/3 of the theoretical calculated Sabine reverberation time of the room.
The total recording time is defined in seconds.

### Air absorption
The absorption of the air will need to be inputted. The air absorption is defined as *m_\{atm\}* and in 1/meters. 

## Fixed inputs
Within the fixed inputs, there are:
- Adiabatic speed of sound defined as 343 m/s;
- Absorption conditions term (Option 1 Sabine, Option 2 Eyring, Option 3 Modified). These are absorption factors for the boundary conditions. Currently, the most accurate absorption factor is set as the Option 3 Modified as it has bee demostrated that this is accurate for low and high absorption coefficients;
- Reference pressure defined as $'2*10^{-5}'$ Pascal;
- Air density at 20°C defined as 1.21 [$'kg/m^{-3}'$].

## Acoustics Calculation
The acoustic calculation is based on the Du Fort and Frankel method (explicit unconditionally stable Finite difference method) solving the diffusion equation (Navarro et al., 2012). More information regarding the Finite Different Method in the paragraph below.

## Acoustics parameters

The diffusion equation method predicts the time-dependent propagation of the sound energy density w(r, t) in the evaluated frequency band. 

#### Sound Density Level
The sound density level can be expressed as function of sound energy density w(r, t) as:
```{math}
SDL = 10 log_{10}⁡(w(r,t))
```

#### Sound Pressure Level
After predicting the time-dependent sound energy density in the room, the sound pressure level decay curve can be expressed as function of sound energy density w(r, t) as:
```{math}
SPL = 20 log_{10}⁡((w(r,t)ρc^2)/p_{ref}^2) 
```
where *p_{ref}* is $'2*10^{-5}'$ Pa and ρ is the air density.

#### Reverberation time and Early Decay Time (EDT)
From the sound pressure level decay curve, the Reverberation time can be estimated. The RT is defined by the time that it takes for the sound pressure level to decay of 60 dB. Depending on the room geometry, occasionally it is difficult to evaluate 60 dB of decay and therefore, the T30 is evaluated. This is obtained from the slope between -5 and -35 dB of the maximum starting level.  

The Early Decay time is defined by the time that it takes for the sound pressure level to decay of 10 dB and it is calculated in a similar way. This is obtained from the slope between 0 and -10 dB of the maximum starting level. 

#### Clarity, Definition and Centre Time

The Clarity (C80) parameter is the early to late arriving sound energy ratio. Clarity refers to how clear the sound quality is and it is calculated from the impulse response with the following relation:

```{math}
C_{80}=10 log⁡(\int_0^{80ms} p^2(t) \, dt/\int_{80ms}^\infty p^2 (t) \, dt)   \,  [dB]
```

The Definition (D50) parameter is the ratio of the early received sound energy (0-50ms after direct sound arrival) to the total received energy. It referres only to the speech and it is defined as: 

```{math}
D_{50}=10 log⁡(\int_0^{50ms} p^2(t) \, dt/\int_0^\infty p^2(t) \, dt)  \,  [%]
```

The Centr Time (Ts) parameter is the center of gravity of the squared impulse response. Centre Time avoids the discrete division of the impulse response into early and late periods. 

```{math}
T_{s}=10 log⁡(\int_0^\infty tp^2(t) \, dt/\int_0^\infty p^2(t) \, dt)   \,  [s]
```

A low value indicate that most of the energy arrives early, a high value reveals that there is much reverberance.

The values for all these parameters are calculated from the Barron’s revisited theory formulas (Vorlander, 2008) with the influence of the direct field neglected.

## References
- J. M. Navarro, J. Escolano and J. J. Lopez, Implementation and evaluation of a diffusion equa-tion model based on finite difference schemes for sound field prediction in rooms,Appl. Acoust.73(2012) 659–665.
- Vorländer M. Auralization: fundamentals of acoustics, modelling, simulation,algorithms and acoustic virtual reality. Springer 2008