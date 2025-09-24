---
title: '`acousticDE`: A diffusion equation model package for room acoustics simulations'
tags:
  - python
  - diffusion equation
  - room acoustics
  - simulations
  - diffusion process
  - building physics
authors:
  - name: Ilaria Fichera
    orcid: 0000-0002-0097-1486
    equal-contrib: true
    corresponding: true
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: CÃ©dric Van hoorickx
    orcid: 0000-0002-9671-5558
    equal-contrib: false # (This is how you can denote equal contributions between multiple authors)
    affiliation: "1" 
  - name: Maarten Hornikx
    orcid: 0000-0002-8343-6613
    equal-contrib: false # (This is how you can denote equal contributions between multiple authors)
    affiliation: "1" 
affiliations:
 - name: Department of Built Environment, Eindhoven University of Technology
   index: 1
   ror: 02c2kyt77
date: 24 September 2025
bibliography: paper.bib

---

# Summary
<!-- A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience. -->

The acoustic properties of enclosed spaces are important for improving the quality of environments and reduce noise induced health effects. With the rise of numerical room acoustics modelling [@Kuttruff2019RoomAcoustics], tools are needed that are both computationally efficient and versatile in capturing sound behaviour in enclosed spaces. `acousticDE` is an open-source computer program using the diffusion equation approach, a modelling method that has gained considerable interest because it strikes an optimal balance between efficiency and flexibility. To simulate the sound field in a room, `acousticDE` requires a digital 3D model of the geometry, as well as material properties describing the acoustic boundaries. It can incorporate a sound sources, compute the spatial and temporal distribution of sound energy, and output key room acoustics parameters such as reverberation time. The diffusion equation method, originally introduced by Ollendorff [Ollendorff1969AbsorptionWaves] and later refined [Picaut1997AEquation], overcomes some of the limitations of statistical room acoustics, which is grounded in the diffuse field assumption and the classical Sabine reverberation formula [Sabine1922CollectedAcoustics]. Compared to wave-based approaches [Hamilton2017FDTDTime], [Wang2019RoomMethod], ray- and beam-tracing algorithms [@Kuttruff2019RoomAcoustics], or the image source method [@Kuttruff2019RoomAcoustics], the DE method provides a computationally lightweight yet physically informed framework for room acoustics simulations, research and practice [Valeau2006OnPrediction].

# Statement of need
<!-- A Statement of need section that clearly illustrates the research purpose of the software and places it in the context of related work. -->

Room acoustics simulation is essential for research and for consultancies to accurately predict and design the acoustics inside a room. However, existing methods often face a trade off between accuracy, transparency in the calculation method and codes and computational efficiency. Traditional statistical methods, such as Sabine and Eyring calculation, provides a fast global estimate of the acoustic parameters, however these are based on unreal assumptions (diffuse field environment) [@Kuttruff2019RoomAcoustics]. Some more accurate approaches are not commercially available and therefore it is very difficult to see how the calculation is actually done insider the software; this reduces transparency in the method of calculation. Additionally, full wave-based approaches offer higher accuracy at the cost of significant computational resources, making them impractical for larger spaces. 

The diffusion equation model has grown to be an effective compromise, providing physically meaningful spatial and temporal distributions of acoustic energy while remaining computationally efficient in the high-frequency range [Foy2016IncludingApproach], [Mou2023AnSpaces], [SuGul2020ComparativeStations]. However, despite its demonstrated usefulness in the literature, accessible and user-friendly software implementations of diffusion equation model have been limited.

`acousticDE` is an open source software package designed for the simulation of room acoustics using the diffusion equation model (DE). 
The diffusion equation software is investigated with two different numerical methods: the Finite Different Method (FDM) by Du Fort&Frankel (Navarro et al., 2012) and the Finite Volume Method (FVM) (Munoz, 2019). The FDM can be used for parallelepiped shapes and the room is discretised into a 3D grid of points where the energy density is calculated. The FVM can be used for any shape/geometry, since the discretization is volumetric and tetrahedrical. Beyond numerical simulation, `acousticDE` enables auralization, allowing users to generate perceptual audio demonstrations of acoustic conditions. The Diffusion Equation method software excels for its easeness of use and understanding and for its computational speed, allowing researchers and engineers to simulate quickly a wide range of room acoustics problems with accuracy and computational efficiency. The software is written in Python language. In its simpler settings, the `acousticDE` package is available as a back-end acoustic simulator of the Community Hub for Open-source Room Acoustics Software [CHORAS](https://github.com/choras-org/CHORAS) [Willemsen2025].

The software has already been applied in a number of scientific publications [Fichera2024DeterminationApproach], [Fichera2025] and graduate-level teaching at the Eindhoven University of Technology (TU/e), where it has proven valuable for understanding the impact of room geometry and material absorption on acoustic behavior. By combining ease of use, computational efficiency, and research-grade accuracy, `acousticDE` fills a critical need for researchers, consultants, and educators seeking a practical and transparent tool for room acoustics simulation.

# Method

The diffusion equation method aims to find the temporal and spatial distribution of acoustic energy within a specific room. The modelling method is based on solving the partial differential Diffusion Equation, which is applicable in the high-frequency range and it assumes only diffusely reflecting boundaries following the  Lambertâ€™s law [Valeau2006OnPrediction]. The diffusion equation together with its boundary condition is as follow:

\begin{equation}
    \begin{cases}
      \frac{\partial{w(\mathbf{r}, t)}} {\partial{t}} = - D \mathbf{\nabla}^2 w(\mathbf{r}, t) - m c w(\mathbf{r}, t) + q(\mathbf{r}, t) \\
      -D \dfrac{\partial w}{\partial n} = h w \\
    \end{cases}  
    \label{eq:diffusionequationwithRobin}
\end{equation}

where $w(\mathbf{r}, t)$ is the energy density at each position $\mathbf{r}$ and at time $t$, $m$ is the atmospheric attenuation coefficient, $c$ is the speed of sound, $q(\mathbf{r}, t)$ is the source term and $h$ is the absorption term. Several methodologies propose expressions for $h$ based on different assumptions on how to treat the suface absorption coefficients of the room [Picaut2002NumericalProcess], [Jing2007APredication], [Jing2008OnExperiments]. The diffusion coefficient $D$ is a term in $m^2/s$ that scales the energy density and indicates how quickly the sound diffuses around the room. This depends on the mean free path of the room.  

This method allows for an efficient way to calculate acoustics parameters in the room (e.g. sound pressure level and reverberation time) and for more accurate solutions of these parameters compared to the statistical methods of Sabine and Eyring, with the advantage of spatial and temporal distributions of the energies in the rooms while keeping the computational cost reasonable. 

# Acknowledgements and Fundings

Ilaria Fichera acknowledges contributions from Silvin Willemsen, Marco Berzborn, Hassan Teymoori, Felipe Raymann, and Radovan Bast for their useful tips. 

This research is founded by the Dutch Research Council (NWO),Applied and Engineering Sciences (AES) under grant agreement No. 19430, with project title â€A new era of room acoustics simulation software: from academic advances to a sustainable open-source project and communityâ€.

# References
<!-- A list of key references, including to other software addressing related needs. Note that the references should include full names of venues, e.g., journals and conferences, not abbreviations only understood in the context of a specific discipline. -->
