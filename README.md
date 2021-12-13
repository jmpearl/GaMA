# GaMA
The Gravitational and Mesh Adapatation library is a set of optimized matlab classes used to model the gravitational fields of asteroids and comets. It contains a number of supporting classes and functions required to initialize gravitational models, integrate trajectories, post process results, and make relevant plots.  

###### Gravitational Model
- analytic polyhedral model [1]
- spherical harmonic model [2]
- mascon model 
- approximate polyhedral model [3]
- curvilinear surface model

###### Meshing
- surface mesh class with:
  - array-base half-edge data structure
  - mesh repair
  - smoothing
  - coarsening and refining
  - feature-based refinement
  - quadrature rule generation
  - plotting utilities
  -
-volume mesh class with:
- - graded mesh generation
- - quadrature rule generation

- [1] R.A. Werner and D.J. Scheeres. Exterior gravitation of a polyhedral derived and compared with harmonic and mascon gravitation representations of asteroid 4769 castalia. Celestial Mechanics and Dynamical Astronomy, 65:313,344, 1997.
- [2] Gottlieb, R.G., "Fast Gravity, Gravity Partials, Normalized Gravity, Gravity Gradient Torque and Magnetic Field: Derivation, Code and Data," NASA CR-188243, 1993.
- [3] Pearl, J.M., Hitt, D.L., "A fast quadrature-based gravity model for the homogeneous polyhedron," MNRAS, Vol. 492, pp. 420-430, 2020.
