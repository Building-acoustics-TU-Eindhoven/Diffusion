## Diffusion Equation Software for Room Acoustics Modeling
acousticDE is an open-source software package developed for the simulation of room acoustics based on the diffusion equation method. The package provides two complementary numerical implementation of this method: the Finite Volume Method (FVM) and the Finite Difference Method (FDM). In addition, it includes a tool for auralization enabling users to listen to acoustic renderings generated from the results of the Finite Volume Method. The main strenght of acousticDE is in its ability to predict the distribution of acoustic energy over space and time within a given room. With only a minimal set of input parameters, the software delivers accurate and computationally efficient estimates of both the energy propagation and room acoustics properties. This makes it a valuable tool for researchers and practitioners interested in architectural and room acoustics. The software is being developed as part of an ongoing research within the Building Acoustics Group at the Department of Built Environment, Eindhoven University of Technology. It is currently under active development and implemented in Python by Ilaria Fichera. 

## Release version
Version 0.1.0

## Repository structure
The package acousticDE is formed by three subfolders: 
+ _Auralization_ generates the auralization wav file for the room in question;
+ _FiniteDifferenceMethod_ computes the room acoustics based on the diffusion equation using the finite difference method (Navarro et al., 2012);
+ _FiniteVolumeMethod_ computes the room acoustics based on the diffusion equation using the finite volume method (Munoz, 2019);

## Installation
Use pip to install acousticDE

```bash 
pip install acousticDE
```

To run the codes/functions, check the [documentation](https://building-acoustics-tu-eindhoven.github.io/Diffusion/) below depending on the method you want to use and check the Tutorial sections of the documentation.

## Usage & Documentation
The [documentation](https://building-acoustics-tu-eindhoven.github.io/Diffusion/) is created to help to use and develop acousticDE effectively. To use acousticDE, please refer to the Tutorial section of the [documentation](https://building-acoustics-tu-eindhoven.github.io/Diffusion/). In addition, the documentation gives an introduction of the package for both FDM and FVM.

## Authors
Software is being developed by Ilaria Fichera at Eindhoven University of Technology (TU/e). 

## Funding
This research is founded by the Dutch Research Council (<u>[NWO](https://www.nwo.nl/projecten/19430)), Applied and Engineering Sciences (AES) under grant agreement No. 19430, with project title â€A new era of room acoustics simulation software: from academic advances to a sustainable open-source project and communityâ€.

## License
Diffusion is under copyright of Building Acoustics Group at the Eindhoven University of Technology and is licensed under GNU General Public License v2.0. See LICENSE.md for more details.

## References
J. M. Navarro, J. Escolano and J. J. Lopez, Implementation and evaluation of a diffusion equation model based on finite difference schemes for sound field prediction in rooms, Applied Acoustics 73 (2012) 659â€“665.

R. P. MuÃ±oz, Numerical modeling for urban sound propagation: developments in wave-based and energy based methods, PhD Thesis, Technische Universiteit Eindhoven, 2019.

