# Finite Different Method Demo Tutorial

Test the software with the following inputs (Navarro et al., 2012):

- Spatial discretization: $\Delta x = 0.5 m$;
- Time discretization: $\Delta t = 1/8000 s$;
- Recording total time: $T = 2.0 s$;
- Room dimension: 8.0 m x 8.0 m x 8.0 m;
- Absorption Coefficient surface $i$: $\alpha = 1/6 \approx 0.17$;
- Air absorption: $m = 0.0$ 1/m;
- Absorption term: Modified Absorption Term $A_{M}$;
- Source position coordinates: 4.0 m, 4.0 m, 4.0 m;
- Receiver position coordinates: 2.0 m, 2.0 m, 2.0 m;

Test if the software provides the following results (Navarro et al., 2012):

- Reverberation time (RT): 1.18 s; 
- Early Decay Time (EDT): 1.18 s;
- Clarity ($C_{80}$): 2.08 dB;
- Definition ($D_{50}$): 45.63 %
- Centre Time ($T_{s}$): 85.33 ms

The result file is a pickle file called _resultsFDM.pkl_. All the results are included in this file.

## References
- J. M. Navarro, J. Escolano, J. J. Lopez, Implementation and evaluation of a diffusion equation model based on finite difference schemes for sound field prediction in rooms, Applied Acoustics 73 (6-7) (2012) 659–665.
