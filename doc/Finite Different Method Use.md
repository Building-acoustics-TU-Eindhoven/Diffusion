# Finite Different Method Use

## Requirements
Set up acousticDE following the instructions in the installation section. 

## Usage & files
To use the software, the following files are to be used:
- _PrepareInputsFDM.py_: to create the json file with the inputs of the room;
- _FDM.py_: it contains the main function run_fdm_sim to run the full simulation and calculate the acoustics parameters in the room.

The main software works with the following associated functions:
+ _FDMfunctions.py_ include all the main functions that are used in the full simulation;
+ _FunctionRT.py_ calculates the reverberation time of the room in question;
+ _FunctionClarity.py_ calculates the clarity $C_{80}$ of the room in question based on Barron's revised theory formula;
+ _FunctionDefinition.py_ calculates the definition $D_{50}$ of the room in question based on Barron's revised theory formula;
+ _FunctionCentreTime.py_ calculates the centre time $T_s$ of the room in question based on Barron's revised theory formula.

## Inputs

The general inputs needs to be set by using the script <a href="https://raw.githubusercontent.com/Building-acoustics-TU-Eindhoven/Diffusion/refs/heads/master/Diffusion_Module_ADE/FiniteDifferenceMethod/PrepareInputsFDM.py" download>⬇ Download PrepareInputsFDM.py</a>.

Please download the script and define the following inputs:
```
input_data = {
    "room_dim": [8.0, 8.0, 8.0], #dimension of the room x,y,z
    "coord_source": [4.0, 4.0, 4.0], #source coordinates x,y,z
    "coord_rec": [2.0, 2.0, 2.0], #rec coordinates x,y,z
    "alpha_1": [1/6, 1/6, 1/6, 1/6, 1/6, 1/6], #Absorption coefficient for Surface1 - Floor
    "alpha_2": [1/6, 1/6, 1/6, 1/6, 1/6, 1/6], #Absorption coefficient for Surface2 - Ceiling
    "alpha_3": [1/6, 1/6, 1/6, 1/6, 1/6, 1/6], #Absorption coefficient for Surface3 - Wall Front
    "alpha_4": [1/6, 1/6, 1/6, 1/6, 1/6, 1/6], #Absorption coefficient for Surface4 - Wall Back
    "alpha_5": [1/6, 1/6, 1/6, 1/6, 1/6, 1/6], #Absorption coefficient for Surface5 - Wall Left
    "alpha_6": [1/6, 1/6, 1/6, 1/6, 1/6, 1/6], #Absorption coefficient for Surface6 - Wall Right
    "fc_low": 125, #lowest frequency
    "fc_high": 4000, #highest frequency
    "num_octave": 1, # 1 or 3 depending on how many octave you want
    "dx": 0.5,
    "dt": 1/8000, #time discretization
    "m_atm": 0, #air absorption coefficient [1/m]
    "th": 3, #int(input("Enter type Absortion conditions (option 1,2,3):")) # options Sabine (th=1), Eyring (th=2) and modified by Xiang (th=3)
    "tcalc": "decay" #Choose "decay" if the objective is to calculate the energy decay of the room with all its energetic parameters; Choose "stationarysource" if the aim is to understand the behaviour of a room subject to a stationary source
}
```

### Room dimensions
The geometry for this method is defined within the inputs variable section of the python file _PrepareInputsFDM.py_.
The method is suitable for parallelepiped and cuboid spaces and therefore, the only inputs are the length, width and height of the room.

### Sound source
The model allows for the insertion of only one source position per calculation. 
The sound source is defined as an omnidirectional source. The users can input the sound power of the source *W_{s}* in Watts and its position in the room in the following variables *x_{source\}*,*y_\{source\}*,*z_\{source\}* in meters in the x,y,z directions.

### Receiver
The model allows for the insertion of only one acoustics receiver position per calculation. These are defined as point omnidirectional receivers. The users will input the position of the receiver in the room in the following variables *x_\{rec\}*,*y_\{rec\}*,*z_\{rec\}* in meters in the x,y,z directions. The receiver position cannot be equal to the source position, otherwise it will give an error.

### Surface material
The surface materials are defined within the inputs variable section of the python file _PrepareInputsFDM.py. The surface material requires a value of absorption coefficient per frequency per surface. It is not possible to include doors, windows or other items within a surface.

### Frequency range
The frequency range for this method is defined within the _PrepareInputsFVM.py_ python script. 
The frequency resolution should be included as inputs variables *fc_\{low\}* and *fc_\{high\}*; these should be the middle frequency of a band. The maximum number of frequencies is set in octave bands and can be choosen by the user. Normally, *fc_\{low\}* is set to 125 Hz and *fc_\{high\}* is set to 2000 Hz. The number of octave need to be also set: 1 for one octave and 3 for third octave bands.

### Spatial discretization $\Delta x$
The Finite Different method works with a spatial discretization. The space is defined by a mesh grid of points at a distance $\Delta v$ between each other. The distance $\Delta v$ is equal for each dimension $x,y,z$, therefore $\Delta v = \Delta x = \Delta y = \Delta z$, and it is defined in meters. 
It is important to choose an appropriate $\Delta v$ for the precision of the calculation. 
A $\Delta v$ of 0.5 meters normally would suffice for a correct calculation, although the choice is contingent upon user preferences for details in parameters values and room dimensions. It is strongly advised that $\Delta v$ is a factor of the dimensions of the room (e.g. for a room of 3 x 3 x 3 $m^3$, $\Delta v$ could be 0.5 meters or 0.2 meters or 0.1 meters but it should not be 0.4 meters). 

### Time discretization $\Delta t$
According the Navarro 2012, to get good converged results, the time discretization $\Delta t$ will need to be defined depending on the $\Delta v$ chosen. 
To make sure that the predictions converge to a fixed value with a very low error, the following empirical cretirion will need to apply.
```{math}
10^{-8} = \Delta t^2 \Delta v^{-2}
```
The time discretization is defined in seconds. 

### Air absorption
The absorption of the air will need to be inputted. The air absorption is defined as *m_\{atm\}* and in 1/meters. 

### Absorption Term
In the diffusion equation model, three different absorption terms conditions can be used. More info about these terms are shown in the software theory documentation of the FVM. The user can choose the preferred one. The options are Option 1 for the Sabine absorption term, Option 2 for Eyring absorption term and Option 3 for the Modified absorption term. These are absorption factors for the boundary conditions. Currently, the most accurate absorption factor is set as the Option 3 Modified as it has been demostrated that this is accurate for low and high absorption coefficients.

### Type of calculation
The type of calculation for this method is defined within the _PrepareInputsFDM.py_ python script. 
The source is defined as an interrupted noise source. The time within which the source stays on before getting switch off is predefined depending on the theoretical calculated Sabine reverberation time of the room in the variable *sourceon_\{time\}*. There is the option to toggle the sound source on and off or keep it stationary. This can be done by changing the variable "tcalc" from "decay" to "stationarysource".


The file _PrepareInputsFDM.py_ will create a .json file with the general inputs. This will be used in the main running function. 


## Acoustics Calculation
Once all the inputs have been set, the main calculation can be run using the function run_fvm_sim as described in the API references. To run, this function, a .json file path is needed. The acoustic calculation is based on the Du Fort and Frankel method (explicit unconditionally stable Finite difference method) solving the diffusion equation (Navarro et al., 2012). More information regarding the Finite Different Method in the paragraph below.

### Acoustics parameters

The diffusion equation method predicts the time-dependent propagation of the sound energy density $w(\mathbf{r}, t)$ in the evaluated frequency band. 

#### Sound Density Level
The sound density level can be expressed as function of sound energy density $w(\mathbf{r}, t)$ as:
```{math}
SDL = 10 log_{10}⁡(w(\mathbf{r}, t))
```

#### Sound Pressure Level
After predicting the time-dependent sound energy density in the room, the sound pressure level decay curve can be expressed as function of sound energy density w(r, t) as:
```{math}
SPL = 20 log_{10} \left( \frac{w(\mathbf{r}, t) \rho c^2}{p_{ref}^2} \right )
```
where $p_{ref}$ is $2 \cdot 10^{-5}$ Pa and $\rho $ is the air density.

#### Reverberation time (RT) and Early Decay Time (EDT)
From the sound pressure level decay curve, the Reverberation time can be estimated. The RT is defined by the time that it takes for the sound pressure level to decay of 60 dB. Depending on the room geometry, occasionally it is difficult to evaluate 60 dB of decay and therefore, the $T_{30}$ is evaluated. This is obtained from the slope between -5 and -35 dB of the maximum starting level.  

The Early Decay time is defined by the time that it takes for the sound pressure level to decay of 10 dB and it is calculated in a similar way. This is obtained from the slope between 0 and -10 dB of the maximum starting level. 

#### Clarity, Definition and Centre Time

The Clarity ($C_{80}$) parameter is the early to late arriving sound energy ratio. Clarity refers to how clear the sound quality is and it is calculated from the impulse response with the following relation:

```{math}
C_{80} = 10 log⁡ \left( \frac{\int_0^{80ms} p^2(t) \, dt}{\int_{80ms}^\infty p^2 (t) \, dt} \right ) \,  [dB]
```

The Definition ($D_{50}$) parameter is the ratio of the early received sound energy (0-50ms after direct sound arrival) to the total received energy. It referres only to the speech and it is defined as: 

```{math}
D_{50} = 10 log \left( \frac{\int_0^{50ms} p^2(t) \, dt}{\int_0^\infty p^2 (t) \, dt} \right ) \,  [\%]
```

The Centre Time (T_s) parameter is the center of gravity of the squared impulse response. Centre Time avoids the discrete division of the impulse response into early and late periods. 

```{math}
T_{s}=10 log⁡ \left( \frac{\int_0^\infty tp^2(t) \, dt}{\int_0^\infty p^2(t) \, dt} \right)   \,  [s]
```

A low value indicate that most of the energy arrives early, a high value reveals that there is much reverberance.

The values for all these parameters are calculated from the Barron’s revisited theory formulas (Vorlander, 2008) with the influence of the direct field neglected.

## References
- J. M. Navarro, J. Escolano, J. J. Lopez, Implementation and evaluation of a diffusion equation model based on finite difference schemes for sound field prediction in rooms, Applied Acoustics 73 (6-7) (2012) 659–665.

- M. Vorländer, Auralization: fundamentals of acoustics, modelling, simulation, algorithms and acoustic virtual reality,  Springer 2008.
