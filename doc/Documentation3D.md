# Diffusion Equation Code 3D documentation

The software is based on the Du Fort and Frankel method (explicit unconditionally stable Finite difference method) solving the diffusion equation (Navarro et al., 2012).
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
- Air absorption $m_\{atm\}$ in 1/m;
- Length, width, height of the room (lxmax, lumax, lzmax in m);
- Position of the source $x_\{source\},y_\{source\},z_\{source\}$ in m in the x,y,z directions;
- Position of the receiver $x_\{rec\},y_\{rec\},z_\{rec\}$ in m in the x,y,z directions;
- Distance between grid points dx, dy and dz in m;
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

## Finite Different Scheme

### Diffusion Equation
The software is implemented using the explicit unconditionally stable numerical Finite Different Method called Du Fort and Frankel to solve the diffusion equation (Navarro et al., 2012).
The Dufort and Frankel method is an expansion of the FTCS method with a truncation order of O[(∆t)^2] + O[(∆v)^2] + O[(∆t)^2 x (∆v)^2]. It is an explicit and unconditionally stable method with ∆t as time step discretization and ∆v=∆x=∆y=∆z as spatial discretization. The ∆v and ∆t are to define depending on the room in order for the results to converge to the exact one (Navarro et al. 2012).

According to the Dufort and Frankel method, the discretization is based on the following:

- i is the spatial element along x direction\
- j is the spatial element along y direction\
- k is the spatial element along z direction\
- n is the temporal element\

with the following surface definition:

![Grid 1D](images/Surfaces.png)

The full discretised partial differential equation is:
```{math}
w_{i,j,k}^{n+1}=  (w_{i,j,k}^{n-1}(1-β_{0} )-2∆tcmw_{i,j,k}^n - 2∆tP_{i_s,j_s,k_s}^n + β_{0{x}}(w_{i+1,j,k}^n+ w_{i-1,j,k}^n )+
+ β_{0{y}}(w_{i,j+1,k}^n+ w_{i,j-1,k}^n )+ β_{0{z}}(w_{i,j,k+1}^n+ w_{i,j,k-1}^n ))/(1+ β_{0})
```

Where:
```{math}
β_{0{x}}=β_{0{y}}=β_{0{z}}=(2D∆t)/(∆x)^2 
```
And
```{math}
β_{0}=β_{0{x}}+β_{0{y}}+β_{0{z}} 
```

Each term of the equation is discretized as follows:
```{math}
∂w(r,t)/∂t => w_{i,j,k}^{n} = (w_{i,j,k}^{n+1} - w_{i,j,k}^{n-1})/(2∆t)
```
```{math}
∂^2 w(r,t)/∂x^2 => w_{i,j,k}^n = (w_{i+1,j,k}^n - 2((w_{i,j,k}^{n+1}+w_{i,j,k}^{n-1})/2)+w_{i-1,j,k}^n)/(∆x)^2
```
```{math}
∂^2 w(r,t)/∂y^2 => w_{i,j,k}^n = (w_{i,j+1,k}^n - 2((w_{i,j,k}^{n+1}+w_{i,j,k}^{n-1})/2)+w_{i,j-1,k}^n)/(∆y)^2
```
```{math}
∂^2 w(r,t)/∂z^2 => w_{i,j,k}^n = (w_{i,j,k+1}^n - 2((w_{i,j,k}^{n+1}+w_{i,j,k}^{n-1})/2)+w_{i,j,k-1}^n)/(∆z)^2
```
```{math}
cmw(r,t) => cmw_{i,j,k}^n
```
```{math}
P(t)δ(r-r_s ) => P_{i_{s},j_{s},k_{s}}^n
```
This is a soft source term at position $i_\{s\},j_\{s\},k_\{s\}$ and at the time step of n.
According to Navarro et al. 2012, the source term function is: \
{math}`w_{i_s,j_s,k_s}^{n+1}= w_{i_s,j_s,k_s}^{n+1}+2∆tP_{i_s,j_s,k_s}^n`

The stability criterion is the following (Navarro et al., 2012):
```{math}
(β_0^{-1}/(1+β_0^2+2β_0)≤ 1
```
Where $\\beta_\{0\}$ is a term of constants in the discretization of the diffusion equation defined below.

### Boundary Conditions

Additional equations are needed to describe the sound field at the boundaries. The discretized equations are based on the forward and backward three points formula (Necati et al., 2017). 
Only the one directional formulas have been included as examples.

1. Forward Difference Approximation (first derivative - three points formula) for x=0, any j and any k
```{math}
∂w(r,t)/∂n => w_{0,j,k}^{n+1} = (-3 w_{0}^{n+1}+4w_{1}^{n+1}- w_{2}^{n+1})/(2∆x)
```
2. Backward Difference Approximation (first derivative - three points formula) for x=L_{x}, any j and any k
```{math}
∂w(r,t)/∂n => w_{L_{x}}^{n+1} = (3 w_{L_{x}}^{n+1}-4w_{L_{x-1}}^{n+1}+ w_{L_{x-2}}^{n+1})/(2∆x)
```
And the discretized boundary reshaped are:

1.1. Forward Difference Approximation for x=0, any j and any k
```{math}
w_{0}^{n+1}=   (4w_{1}^{n+1}-2w_{2}^{n+1})/(3+(2 A_{x_{0}}∆x)/D_{x})=boundary at x=0
```
2.1. Backward Difference Approximation for x=L_{x}, any j and any k
```{math}
w_{L_{x}}^{n+1}=   (4w_{L_{x-1}}^{n+1}-2w_{L_{x-2}}^{n+1})/(3+(2 A_{x_{L_x}}∆x)/D_{x})=boundary at x=Lx
```

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
- J. M. Navarro, J. Escolano and J. J. Lopez, Implementation and evaluation of a diffusion equa-tion model based on finite difference schemes for sound field prediction in rooms,Appl. Acoust.73(2012) 659–665.
- Billon A, Picaut J, Foy C, Valeau V, Sakout A. Introducing atmospheric attenuation within a diffusion model for room-acoustic predictions. J Acoust Soc Am. 2008 Jun;123(6):4040-3. doi: 10.1121/1.2903872. PMID: 18537354.
- Picaut, J., L. Simon, and J. D. Polack. 1997. “A Mathematical Model of Diffuse Sound Field Based on a Diffusion Equation.” Acta Acustica United with Acustica 83 (4): 614–621.
- Valeau, V., J. Picaut, and M. Hodgson. 2006. “On the Use of a Diffusion Equation for Room-Acoustic Prediction.” Journal of the Acoustical Society of America 119 (3): 1504–1513.
- Jing, Y., and N. Xiang. 2007. “A Modified Diffusion Equation for Room-Acoustic Predication (L).” Journal of the Acoustical Society of America 121 (6): 3284–3287.
- Billon, A., J. Picaut, and A. Sakout. 2008. “Prediction of the Reverberation Time in High Absorbent Room Using a Modified-Diffusion Model.” Applied Acoustics 69 (1): 68–74.
- Jing, Y., and N. Xiang. 2008. “On Boundary Conditions for the Diffusion Equation in Room Acoustic Predictions: Theory, Simulations, and Experiments.” Journal of the Acoustical Society of America 123 (1): 145–153.
- Ozisik, M. Necati.  (1977).  Basic heat transfer.  New York :  McGraw-Hill
- Vorländer M. Auralization: fundamentals of acoustics, modelling, simulation,algorithms and acoustic virtual reality. Springer 2008