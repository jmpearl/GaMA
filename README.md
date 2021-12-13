# GaMA
The Gravitational and Mesh Adaptation library is a set of optimized matlab classes used to model the gravitational fields of asteroids and comets. It contains a number of supporting classes and functions required to initialize gravitational models, integrate trajectories, and post-process results. The gravitational models are easy to initialize and contain methods to calculate the potential, acceleration, laplacian, and gravitational gradient given a position. The current version is serial. 

###### Gravitational Model
- analytic polyhedral model [1]
- spherical harmonic model [2]
- mascon model 
- approximate polyhedral model [3]
- curvilinear surface model
- third body model

###### Meshing
- surface mesh class with:
  - an efficient array-base half-edge data structure, curvilinear faces up to degree 4, mesh repair, smoothing, coarsening, refinement, projection, quadrature rule generation, and plotting utilities.
- volume mesh class with:
  - array based face-vertex data structure, curvilinear cells up to degree 2, graded mesh generation, mesh smoothing, quadrature rule generation


###### References
- [1] R.A. Werner and D.J. Scheeres. Exterior gravitation of a polyhedral derived and compared with harmonic and mascon gravitation representations of asteroid 4769 castalia. Celestial Mechanics and Dynamical Astronomy, 65:313,344, 1997.
- [2] Gottlieb, R.G., "Fast gravity, gravity partials, normalized gravity, gravity gradient torque and magnetic field: derivation, code and data," NASA CR-188243, 1993.
- [3] Pearl, J.M., Hitt, D.L., "A fast quadrature-based gravity model for the homogeneous polyhedron," MNRAS, Vol. 492, pp. 420-430, 2020.
