# Diffusion Equation Finite Volume Method (FVM) Documentation

The software is based on the Finite Volume method (FVM) solving the diffusion equation (Munoz, 2019).

## Theory of Diffusion Equation model

The model for the sound energy density $w(\mathbf{r}, t)$ at position $\mathbf{r} = [x,y,z]$ and at time $t$ on a domain $V$ is based on the following partial differential equation:

```{math}
\frac{\partial w(\mathbf{r}, t)}{\partial t}
- D \left( \frac{\partial^2 w(\mathbf{r}, t)}{\partial x^2}
+ \frac{\partial^2 w(\mathbf{r}, t)}{\partial y^2}
+ \frac{\partial^2 w(\mathbf{r}, t)}{\partial z^2} \right)
+ c m w(\mathbf{r}, t) = P(t)\delta(\mathbf{r} - \mathbf{r}_s)
```

where $\frac{\partial^2}{\partial x^2}$, $\frac{\partial^2}{\partial y^2}$, $\frac{\partial^2}{\partial z^2}$ are the Laplace operators and $D = \frac{\lambda c}{3}$ is the so-called theoretical diffusion coefficient with $c$ being the speed of sound. The diffusion coefficient is a constant value that takes into account the room geometry and volume trough the mean free path defined for proportionate rooms as $\lambda = \frac{4 V}{S}$ with volume $V$ of the room and $S$ the total surface area. The term $P(t)$ indicates a sound source term at position $r_s$. The term $c m w(\mathbf{r}, t)$ accounts for the atmospheric attenuation within the room, where $m$ is the absorption coefficient of air (Billon et al., 2008).

The main partial differential equation is associated with mixed boundary conditions on a domain $\partial V$, as follows:


```{math}
-D \frac{\partial w(\mathbf{r}, t)} {\partial n} = A_{r}(\mathbf{r},\alpha ) c w(\mathbf{r}, t) \text{ on } \partial V
```
This models the effect of the absorption at the surfaces in the sound field.
The term $n$ indicates the vector normal to the surface and the term $A_{r}$ is dependent on the absorption coefficient $\alpha$. Depending on the scope and on the absorption coefficient to use, the absorption term $A_{r}$ is defined differently. There are three different absorption factors $A_{r}$: the Sabine (Picaut at al., 1999; Valeau at al., 2006), the Eyring (Jing et al., 2007; Billon et al., 2008) and the modified by Xiang (Jing et al., 2008).

## Finite Volume Scheme

### Diffusion Equation
The software is implemented using the numerical Finite Volume Method to solve the partial differential diffusion equation in the form of algebraic equations (Munoz, 2019).
The Finite Volume method is based on the discretization of the room in discrete nodes surrounded by finite volumes within the volume domain.

As for the finite element method, the finite volume method consists on the creation of a mesh of the room. The mesh elements are called control volumes. The solution is obtained from the fluxes of energy from one volume element to the other, at each time step, until the total recording time. The results from the diffusion equation and boundary condition are calculated in the centre of each control volume.

For the discretization, the following are important:
- $j$ control volume;
- $k \text{control volume adjacent to } j$;
- $w_j \text{energy density computed at the centre of control volume } j$;
- $w_k$ energy density computed at the centre of control volume $k$;
- $S_{j,k} \text{common area between control volume} $j$ \text{and control volume} k$;
- $d_{j,k} \text{distance between centre of control volume} $j$ \text{and centre of control volume} k$.

The diffusion equation need to be integrated over one element $j$ via the following equation:

```{math}
\frac{\partial}{\partial t} \int_{\Omega_j}  w(\mathbf{r}, t) dr - D \int_{\partial Omega} \nabla  w(\mathbf{r}, t) n dr = \int_{\Omega_j} P(t) delta(r-r_s )
```

### Boundary Conditions
Additional equations are needed to describe the sound field at the boundaries. 
The boundary condition above can be discretised for every face of a control volume being part of the boundary of the domain. 
```{math}
- D \frac{\partial}{\partial n}  w_j(\mathbf{r}, t) = - h_{(b)j,k}  w_j
```
The term $n$ indicates the vector normal to the surface and the term $h_{(b)j,k}$ is dependent on the absorption coefficient $\alpha$ of the surface are of the $j$ control volume.

### Discretization
The full discretised partial differential diffusion equation is:
```{math}
\frac{V_j (w_{j}^{n+1} - w_{j}^{n-1})}{2 \Delta t} - D \sum_{k=1}^{N_f} \frac{S_{i,k}}{h_{j,k}} \frac{w_{k}^{n} - (w_{j}^{n+1} - w_{j}^{n-1}}{2}) - \sum_{k=1}^{N_{bf}} S_{i,k}  h_{(b)j,k} \frac{(w_{j}^{n+1} - w_{j}^{n-1}}{2} = V_j P_{j}^{n}
```

## References
- R. P. Muñoz, "Numerical modeling for urban sound propagation: developments in wave-based and energy based methods," PhD Thesis, Technische Universiteit Eindhoven, 2019.
- Vorländer M. Auralization: fundamentals of acoustics, modelling, simulation,algorithms and acoustic virtual reality. Springer 2008