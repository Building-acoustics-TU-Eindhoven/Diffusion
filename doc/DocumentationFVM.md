# Diffusion Equation Finite Volume Method (FVM) Documentation

The software is based on the Finite Volume method (FVM) solving the diffusion equation (Munoz, 2019).

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
The software is implemented using the numerical Finite Volume Method to solve the partial differential diffusion equation in the form of algebraic equations (Munoz, 2019).
The Finite Volume method is based on the discretization of the room in discrete nodes surrounded by finite volumes within the volume domain.

As for the finite element method, the finite volume method consists on the creation of a mesh of the room. The mesh elements are called control volumes. The solution is obtained from the fluxes of energy from one volume element to the other, at each time step, until the recording time. The results from the diffusion equation and boundary condition are calculated in the centre of each control volume.

For the discrtization, the following are important:
- j control volume
- k control volume adjacent to j
- wj energy density computed at the centre of control volume j
- wk energy density computed at the centre of control volume k
- Sjk common area between control volume j and control volume k
- djk distance between centre of control volumej and centre of control volume k

The diffusion equation need to be integrated over one element j via the following equation:

```{math}
$$ ∂/dt \int_{\Omega_j}  w(r,t)dr - D \int_{\∂Omega} \nabla  w(r,t)ndr = \int_{\Omega_j} P(t)δ(r-r_s )
```

### Boundary Conditions
Additional equations are needed to describe the sound field at the boundaries. 
The boundary condition above can be discretised for every face of a control volume being part of the boundary of the domain. 
```{math}
- D \nabla  w_j * n = - D ∂/∂n  w_j(r,t) = - h_{(b)j,k}  w_j
```

### Discretization
The full discretised partial differential diffusion equation is:
```{math}
V_j * (w_{j}^{n+1} - w_{j}^{n-1})/2*dt - D \sum_{k=1}^{Nf} S_i,k/hjk ( w_{k}^{n} - (w_{j}^{n+1} - w_{j}^{n-1})/2) - \sum_{k=1}^{Nbf} S_i,k*  h_{(b)j,k} (w_{j}^{n+1} - w_{j}^{n-1})/2) = V_j P_{j}^{n}
```

## References
- R. P. Muñoz, "Numerical modeling for urban sound propagation: developments in wave-based and energy based methods," PhD Thesis, Technische Universiteit Eindhoven, 2019.
- Vorländer M. Auralization: fundamentals of acoustics, modelling, simulation,algorithms and acoustic virtual reality. Springer 2008