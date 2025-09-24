# Finite Volume Method Demo Tutorial

After installing the package and reading the documentation in [Finite Volume Method Use Documentation](https://building-acoustics-tu-eindhoven.github.io/Diffusion/Finite%20Volume%20Method%20Use.html), the software can be tested with the following files:
- _cube.skp_
- _cube.geo_
- _cube.msh_

The following files can be found in the example folder.

The file _cube.skp_ is a SketchUp file with a room of volume 3x3x3 $m^3$. The files _cube.geo_ and _cube.msh_ are both important files created and saved according to the instructions in [Finite Volume Method Use Documentation, Paragraph 'Geometry & Mesh'](https://building-acoustics-tu-eindhoven.github.io/Diffusion/Finite%20Volume%20Method%20Use.html#geometry-mesh). 

The inputs needs to be prepared. For that, follow the instruction in [Finite Volume Method Use Documentation, Paragraph 'General Inputs'](https://building-acoustics-tu-eindhoven.github.io/Diffusion/Finite%20Volume%20Method%20Use.html#general-inputs). This will create a csv file that you need to fill with the absorption coefficients.

Once all the files are created (msh file, json file and csv file), the main acoustics simulation can be run as below:

```
results = run_fvm_sim('C:\....\cube.msh', 'C:\....\cube_input_fvm.json', 'C:\....\absorption_coefficients.csv')
```

The software should provide results of Sound Pressure Levels at the receiver position, Reverberation time, Clarity and other energetic parameters in a pickle file called _resultsFVM.pkl_.
