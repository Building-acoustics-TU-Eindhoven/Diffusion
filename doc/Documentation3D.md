# Diffusion Equation Code 3D documentation

The code is based on the Du Fort and Frankel method (explicit unconditionally stable Finite difference method) solving the diffusion equation.

## Objectives/what it can do
This code can calculate the acoustic energy density over space and time caused by an omnidirectional sound source located at an arbitrary point in space. It will return the acoustic energy distribution over a 3D room.

## Requirements
Install Spyder 5.4.2

## Libraries
To properly run the code, the following libraries are needed:
- Python version 3.10.9 installed through Anaconda
- Libraries for python:
- import math
- import matplotlib.pyplot as plt #import matplotlib as mpl
- import numpy as np
- from math import ceil
- from math import log
- from math import pi
- from FunctionRT import *

## Inputs
The inputs are:
- Length, width, height of the room (lxmax, lumax, lzmax in m);
- Air absorption $m_{atm}$ in 1/m;
- Distance between grid points dx, dy and dz in m;
- Time discretization dt in s;
- Recording time of the calculation in s;
- Source time in s (time of the source being on before interrupting);
- Absorption conditions term (option 1 Sabine, 2 Eyring, 3 Modified); 
- Absorption coefficient $\alpha_i$ of each surface $$ for one frequency only;
- Source point power $W_s$ in Watts;
- Volume of the source $V_s$ in $m^3$;
- Position of the source $x_{source},y_{source},z_{source}$ in m in the x,y,z directions;
- Position of the receiver $x_{rec},y_{rec},z_{rec}$ in m in the x,y,z directions.

## Implementation
1. Geometry model - Set up the dimension of the room 
2. Set up the spatial distribution and temporal distribution (dx, dy, dz and dt) 
3. Calculate Absorption term based on one absorption coefficient and the type of absorption condition (Sabine, Eyring or Modified) \
4. Define Diffusion coefficient Dx and Dy as constants:
```{math}
Dx = Dy = Dz = (?c)/3
```
Where:
? = mean free path

5. Define stability criterion (Eq. 18 of Navarro 2012) 
```{math}
(?_0^{-1}/(1+?_0^2+2?_0)? 1
```
Where $\beta_{0}$ comes from the discretization of the diffusion equation.

6. Define source:
The source is positioned at x,y,z (to define) in the room and it changes over time depending on the recording time and the source time. The source is defined with intermittent noise level.

7. Define receiver: 
The receiver position will need to be defined.

8. Define the energy density temporal discretization:
- $w_{new}$ = unknown energy density at new time level (n+1)
- $w$ = energy density at n level
- $w_{old}$ = energy density at n-1 level

9. According to the Dufort and Frankel method, the discretaization is based on the following image – based on Navarro 2012

![Grid 1D](images/Surfaces.png)

- i is the spatial element along x direction
- j is the spatial element along y direction
- k is the spatial element along z direction
- n is the temporal element

The partial differential equation is:
```{math}
?w/?t- D(?^2 w/?x^2+?^2 w/?y^2+?^2 w/?z^2)+ cmw=P?(r-r_s ) in V
```
Each term of the equation is discretized as follows:
```{math}
?w/?t=(w_{i,j,k}^{n+1}- w_{i,j,k}^{n-1})/(2?t)
```
```{math}
?^2 w/?x^2=(w_{i+1,j,k}^n - 2((w_{i,j,k}^{n+1}+w_{i,j,k}^{n-1})/2)+w_{i-1,j,k}^n)/(?x)^2
```
```{math}
?^2 w/?y^2=(w_{i,j+1,k}^n - 2((w_{i,j,k}^{n+1}+w_{i,j,k}^{n-1})/2)+w_{i,j-1,k}^n)/(?y)^2
```
```{math}
?^2 w/?z^2=(w_{i,j,k+1}^n - 2((w_{i,j,k}^{n+1}+w_{i,j,k}^{n-1})/2)+w_{i,j,k-1}^n)/(?z)^2
```
```{math}
cmw=cmw_{i,j,k}^n
```
```{math}
P?(r-r_s )=P_{i,j,k}^n
```
The full discretised equation is:
```{math}
w_{i,j,k}^{n+1}=   (w_{i,j,k}^{n-1}(1-?_{0} )-2?tc_0mw_{i,j,k}^n - 2?tP_{i_s,j_s,k_s}^n + ?_{0{x}}(w_{i+1,j,k}^n+ w_{i-1,j,k}^n )+ ?_{0{y}}(w_{i,j+1,k}^n+ w_{i,j-1,k}^n )+ ?_{0{z}}(w_{i,j,k+1}^n+ w_{i,j,k-1}^n ))/(1+ ?_{0}),
```
where:
```{math}
?_{0{x}}=?_{0{y}}=?_{0{z}}=(2D?t)/(?x)^2, 
```
and
```{math}
?_{0}=?_{0{x}}+?_{0{y}}+?_{0{z}}, 
```
and
```{math}
c_0=343 m/s.
```
$P_(i_s,j_s,k_s)^n$=source term (soft source) at position $i_{s},j_{s},k_{s}$ and at the time step of n; this comes from the source term function (Navarro 2012) 

{math}`w_{i_s,j_s,k_s}^{n+1}= w_{i_s,j_s,k_s}^{n+1}+2∆tP_{i_s,j_s,k_s}^n`

m = air absorption coefficient = 0 from Billon paper 2008

10. Define boundary conditions:
:::{table} Boundary conditions
:widths: auto
:align: center

| Boundaries         | Formula                  |
| ------------------ | ------------------------ |
| Boundary at x = 0  | {math}`D ?w/?x-cA_x w=0` |
| Boundary at y = 0  | {math}`D ?w/?y-cA_y w=0` |
| Boundary at x = Lx | {math}`D ?w/?x+cA_x w=0` |
| Boundary at y = Lx | {math}`D ?w/?y+cA_y w=0` |
:::

The discretized equations are based on the forward and backward three points formula (Necati book). Only the one directional formulas have been included as examples.

Forward Difference Approximation (first derivative - three points formula) for x=0
```{math}
?w/?x=(-3 w_{0}^{n+1}+4w_{1}^{n+1}- w_{2}^{n+1})/(2?x)
```
Backward Difference Approximation (first derivative - three points formula) for x=Lx
```{math}
?w/?x=(-3 w_{0}^{n+1}+4w_{1}^{n+1}- w_{2}^{n+1})/(2?x)
```
And the discretized boundary reshaped are:
Boundary at x = 0
```{math}
w_{0}^{n+1}=   (4w_{1}^{n+1}-2w_{2}^{n+1})/(3+(2 A_{x_{0}}?x)/D_{x})= \mbox{boundary at} x=0
```
Boundary at x = Lx
```{math}
w_{L_{x}}^{n+1}=   (4w_{L_{x-1}}^{n+1}-2w_{L_{x-2}}^{n+1})/(3+(2 A_{x_{L_x}}?x)/D_{x})=\mbox{boundary at} x=L_x
```
11. Calculate Sound Density Level and Sound Pressure Level:
```{math}
SDL = 10 log_{10}?(w_{i,j,k}^{n+1})
```
```{math}
SPL = 20 log_{10?}(\frac{w_{i_r,j_r,k_r}^{n+1}?c_{0}^2}{p_{ref}^2}) 
```
	
12. Calculate Reverberation time:
The Reverberation time is calculated using backward integration and in between -5 and -35 dB difference. 

13. Calculate Clarity, Definition and Centre Time
The values are calculate from the Barron’s revisited theory formulas (Vorlander 2008) with the influence of the direct field neglected.

## References
- J. M. Navarro, J. Escolano and J. J. Lopez, Implementation and evaluation of a diffusion equa-tion model based on finite difference schemes for sound field prediction in rooms,Appl. Acoust.73(2012) 659–665.
- Billon A, Picaut J, Foy C, Valeau V, Sakout A. Introducing atmospheric attenuation within a diffusion model for room-acoustic predictions. J Acoust Soc Am. 2008 Jun;123(6):4040-3. doi: 10.1121/1.2903872. PMID: 18537354.
- Vorländer M. Auralization: fundamentals of acoustics, modelling, simulation,algorithms and acoustic virtual reality. Springer 2008