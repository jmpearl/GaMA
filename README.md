# GaMA
The Gravitational and Mesh Adaptation library is a set of optimized matlab classes used to model the gravitational fields of asteroids and comets. It contains a number of supporting classes and functions required to initialize gravitational models, integrate trajectories, and post-process results. The gravitational models are easy to initialize and contain methods to calculate the potential, acceleration, laplacian, and gravitational gradient given a position. It also gives the user the ability to create composite models -- superimposing as many base models as desired. The current version is serial. 

###### Getting Started
To get access to the classes and functions that make up GaMA add the directory tree to the run path. This can be done by running the following command:

```
addpath(genpath('/path/to/GaMA')) 
```

A number of example scripts highlight different functionality within the code and can be found in the Examples subdirectory. Included is a [script]<Examples/GravityModels/IntegrateTrajectory.m> that integrates a trajectory around Asteroid Eros and is a good starting point. 

GaMA can be run with Matlab or GNU Octave and does not require pay-for toolboxes. The current version of GaMA has been tested on Matlab versions R2018b-R2021b and GNU Octave 5.2.0. Portions of the code make use of [implicit expansion](https://blogs.mathworks.com/loren/2016/10/24/matlab-arithmetic-expands-in-r2016b/) and
this will cause errors when used with Matlab R2016b or earlier. The vecnorm function is also used which was not introduced until R2017b. 

###### Contributing
Contribution are welcome and those interested in contributing should see the [CONTRIBUTING]<CONTRIBUTING>  file for details.

###### Release and License
Copyright (c) 2022, Lawrence Livermore National Security, LLC. 
Produced at Lawrence Livermore National Laboratory. 
LLNL-CODE-XXXXXXXX

Written by Jason M. Pearl [pearl3@llnl.gov](mailto:pearl3@llnl.gov)


see the [LISCENSE]<LICENSE>  and [NOTICE]<NOTICE> files for details.
