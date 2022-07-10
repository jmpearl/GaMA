# GaMA
The Gravitational and Mesh Adaptation library is a set of optimized matlab classes used to model the gravitational fields of asteroids and comets. Implemented gravity models include: the spherical harmonic model [[1]](https://ntrs.nasa.gov/citations/19940025085), analytic polyhedral model [[2]](https://doi.org/10.1007/BF00053511), mascon model, approximate polyhedral model [[3]](https://doi.org/10.1093/mnras/stz3461), and a curvilinear surface model.  It contains a number of supporting classes and functions required to initialize the gravitational models, integrate trajectories, and post-process results. The gravitational models are easy to initialize and contain methods to calculate the potential, acceleration, laplacian, and gravitational gradient given a position. It also gives the user the ability to create composite models -- superimposing as many base models as desired. The current version is serial. 

### Getting Started
To access the classes and functions of GaMA, add the directory tree to the run path. This can be done by running the following command:

```
addpath(genpath('/path/to/GaMA')) 
```
A number of example scripts highlight different functionality within the code and can be found in the Examples subdirectory. Included is a [script](scripts/Examples/GravityModels/TrajectoryIntegration.m) that integrates a trajectory around Asteroid Eros and is a good starting point. Another [script](scripts/Examples/GravityModels/ErosSurfaceSlopes.m) shows how to calculate the gravity slopes on the surface of Eros and export the results to a vtk file for viewing with Paraview or Visit.

![Eros Slopes](ErosSlopes.png)


### Version Requirements
GaMA can be run with Matlab or GNU Octave* and does not require pay-for toolboxes. The current version of GaMA has been tested on Matlab versions R2018b-R2021b and GNU Octave 5.2.0. Portions of the code make use of [implicit expansion](https://blogs.mathworks.com/loren/2016/10/24/matlab-arithmetic-expands-in-r2016b/) and
this will cause errors when used with Matlab R2016b or earlier. The vecnorm function is also used which was not introduced until R2017b. Note with GNU Octave performance is notably worse and
```
warning('off','Octave:data-file-in-path')
```
can be used to suppress the fopen warning associated with searching the load path for file names.

### Contributing
Contribution are welcome and those interested in contributing should see the [CONTRIBUTING](CONTRIBUTING)  file for details.

### Release and License
Copyright (c) 2022, Lawrence Livermore National Security, LLC. 

This work was produced under the auspices of the U.S. Department of Energy by Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344, with additional funded by the Vermont Space Grant Consortium under NASA Cooperative Agreement NNX15AP86H.

LLNL-CODE-830936

Written by Jason M. Pearl [pearl3@llnl.gov](mailto:pearl3@llnl.gov)


see the [LISCENSE](LICENSE) and [NOTICE](NOTICE) files for details.
