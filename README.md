# cerebral_fluid_dynamics
Developed code for the analysis and presentation of neural artery blood flow. 

This code uses two toolkit packages to aid in the fluid dynamics computations. I have included the entirely of their toolkits in this upload for simplicity of future users if they would like to do further analysis.

The input pressure data in the base neuralarteries script relies on computationally produced values from the Anatomically Detailed Arterial Network by Blanco et al. (2015). A limitation of this source is that the pressure values are only 1-Dimensional. 

The included function package applying fourier computations is also 1-Dimensional. The included vorticity package calculating vorticity is 2-Dimensional and required further work with 2-Dimensional data to finish connection.
