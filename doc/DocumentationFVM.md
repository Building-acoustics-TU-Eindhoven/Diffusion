# Diffusion Equation Code Finite Volume Method Documentation

The software is based on the Finite Volume method (FVM) solving the diffusion equation (Munoz, 2019).
It is suitable for cuboid spaces.

## Objectives
The software provides the energy density and sound pressure level over 3D space and time and the energy parameters such as Reverberation time (RT), Early Decay Time (EDT), Clarity, Definition and Centre time.

## Requirements
1. Download and install Anaconda or download and install any Python software IDE you prefer
2. Clone/Fork this repository to a folder of your preference
3. Open the Main files through the preferred IDE and test the software

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

## Inputs
The inputs are:
- SketchUp model
- G mesh creation
- Air absorption $m_\{atm\}$ in 1/m;
- Position of the source $x_\{source\},y_\{source\},z_\{source\}$ in m in the x,y,z directions;
- Position of the receiver $x_\{rec\},y_\{rec\},z_\{rec\}$ in m in the x,y,z directions;
- Time discretization dt in s;
- Recording time of the calculation in s;
- Absorption conditions term (option 1 Sabine, 2 Eyring, 3 Modified);
- Absorption coefficient $\\alpha_i$ of each surface $$ for one frequency only;
- Source point power $W_s$ in Watts;
- Volume of the source $V_s$ in $m^3$;
- Source time in s (time of the source being on before interrupting).

## Theory of Diffusion Equation model

The model for the sound energy density w(r,t) at position r and at time t on a domain V is based on the following partial differential equation:
```{math}
∂w(r,t)/∂t- D((∂^2 w(r,t))/(∂x^2 )+(∂^2 w(r,t))/(∂y^2 )+(∂^2 w(r,t))/(∂z^2 ))+ cmw(r,t)=P(t)δ(r-r_s ) in V
```
where ∂^2/dv is the Laplace operator and D = (λc)/3 is the so-called diffusion coefficient with c being the speed of sound. The diffusion coefficient is a constant value that takes into account the room geometry and volume trough the mean free path defined for proportionate rooms as 4*V/S with volume V of the room and S the total surface area. The term P(t) indicates a sound source term at position r_s. The term cmw(r, t) accounts for the atmospheric attenuation within the room, where m is the absorption coefficient of air (Billon et al., 2008).

The main partial differential equation is associated with mixed boundary conditions on a domain ∂V, as follows:
```{math}
-D*∂w(r,t)/∂n = A_{x}(r,α)*cw(r,t) on ∂V
```
This models the effect of the absorption at the surfaces in the sound field.
The term n indicates the vector normal to the surface and the term A_{x} is dependent on the absorption coefficient α. Depending on the scope and on the absorption coefficient to use, the absorption term A_{x} is defined differently. There are three different absorption factors: the Sabine (Picaut at al., 1999; Valeau at al., 2006), the Eyring (Jing et al., 2007; Billon et al., 2008) and the modified by Xiang (Jing et al., 2008).

## Finite Volume Scheme

### Diffusion Equation


![Grid 1D](images/Surfaces.png)

The full discretised partial differential equation is:


### Boundary Conditions

Additional equations are needed to describe the sound field at the boundaries. 


## Results

The diffusion equation method predicts the time-dependent propagation of the sound energy density w(r, t) in the evaluated frequency band. 

#### Sound Density Level
The sound density level can be expressed as:
```{math}
SDL = 10 log_10⁡(w(r,t))
```
which discretised becomes:
```{math}
SDL = 10 log_10⁡(w_{i,j,k}^{n+1})
```

#### Sound Pressure Level
After predicting the time-dependent sound energy density in the room, the sound pressure level decay curve can be expressed as:
```{math}
SPL = 20 log_10⁡((w(r,t)*ρ*c^2)/p_{ref}^2) 
```
which discretised becomes:
```{math}
SPL = 20 log_10⁡(((w_{i_r,j_r,k_r}^{n+1}ρc_{0}^2 ))/p_{ref}^2) 
```
where p_{ref} is 2 × 10−5 Pa and ρ is the air density.

#### Reverberation time and EDT
From the sound pressure level decay curve, the Reverberation time can be estimated using the Schroeder method to the decay curve(Schroeder, 1965). The method consists on performing a backward integration of the slope of the decay. The RT is defined by the time that it takes for the sound pressure level to decay of 60 dB. Depending on the room geometry, occasionally it is difficult to evaluate 60 dB of decay and therefore, the T30 is evaluated. This is obtained from the slope between -5 and -35 dB of the maximum starting level.  

The Early Decay time is defined by the time that it takes for the sound pressure level to decay of 10 dB and it is calculated in a similar way. This is obtained from the slope between 0 and -10 dB of the maximum starting level. 

#### Clarity, Definition and Centre Time

The Clarity (C80) parameter is the early to late arriving sound energy ratio. Clarity refers to how clear the sound quality is and it is calculated from the impulse response with the following relation:
C_80=10 log⁡〖(∫_0^80ms▒〖p^2 (t)〗 dt)/(∫_80ms^∞▒〖p^2 (t)〗 dt)〗    [dB]

The Definition (D50) parameter is the ratio of the early received sound energy (0-50ms after direct sound arrival) to the total received energy. It referres only to the speech and it is defined as: 

D_50=10 log⁡〖(∫_0^50ms▒〖p^2 (t)  dt〗)/(∫_50ms^∞▒〖p^2 (t) 〗 dt)〗    [%]

The Centr Time (Ts) parameter is the center of gravity of the squared impulse response. Centre Time avoids the discrete division of the impulse response into early and late periods. 

T_s=10 log⁡〖(∫_0^∞▒〖τ∙p^2 (t)  dt〗)/(∫_0^∞▒〖p^2 (t) 〗  dt)〗    [s]
A low value indicate that most of the energy arrives early, a high value reveals that there is much reverberance.

The values for all these parameters are calculated from the Barron’s revisited theory formulas (Vorlander, 2008) with the influence of the direct field neglected.

## Algorithm
The software is organised in three sections:
- Input variables:
    The inputs regarding the room dimensions, source and receiver positions along with other are to be inserted for the specific room in question.
- Calculation loop:
    The for loop would loop over the time to calculate the energy density at each position in the mesh grid and at each time step.
- Results and Post-processing
    Results of SPL and other are included in this section together with graphs for the analysis.

## References
- Vorländer M. Auralization: fundamentals of acoustics, modelling, simulation,algorithms and acoustic virtual reality. Springer 2008