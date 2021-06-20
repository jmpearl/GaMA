clear all; clc; close all;

% body parameters
G = 6.67*10^-11;                       % gravitational constant (N m2/kg2)
load('Eros.mat');                      % load stored Eros properties                                  
Mu = bodyProperties.mass*G;            % set stand grav parameter
Omega = bodyProperties.rotationRate;   % 

% some model parameters
meshfile = 'Eros_46906.obj';               % mesh file to load 
mesh = SurfaceMesh(meshfile);              % create the surface mesh object
meshCoarse = mesh.coarsen(8000);           % create coarsened version
meshCoarse = meshCoarse.smooth(5);         % smooth things out

meshCoarseP2 = meshCoarse.coarsen(2000);
meshCoarseP2 = meshCoarseP2.smooth(5);
meshCoarseP2 = meshCoarseP2.setDegree(2);   % make a P2 mesh
meshCoarseP2 = meshCoarseP2.curve(mesh);    % curve mesh through projection

% initialize our gravity models on the coarsened mesh
quadratureModelP2 = ApproximatePolyhedralModel(meshCoarseP2,Mu);


analyticModelTruth = AnalyticPolyhedralModel(mesh,Mu);

% calculate grav field and error on surface
pts = mesh.coordinates; 
truthAcceleration = analyticModelTruth.acceleration(pts);
approximateAccelerationP2 = quadratureModelP2.acceleration(pts);
totApproxErrorPotentialP2 = 100*vecnorm(approximateAccelerationP2-truthAcceleration,2,2)./vecnorm(truthAcceleration,2,2);


meshOffset = mesh.offsetSurfaceMesh(meshCoarse.resolution,mesh.numVertices); 
pts=meshOffset.coordinates;
truthAcceleration = analyticModelTruth.acceleration(pts);
approximateAccelerationP2 = quadratureModelP2.acceleration(pts);
totApproxErrorPotentialP2Offset = 100*vecnorm(approximateAccelerationP2-truthAcceleration,2,2)./vecnorm(truthAcceleration,2,2);


% load the fields in and export for viewing w/ paraview
mesh = mesh.clearFields();
mesh = mesh.addNodeField(totApproxErrorPotentialP2,'TotalApproxErrorAccelerationP2');

mesh.writeVTK('Eros_SurfaceGravFields.vtk');

meshOffset = meshOffset.clearFields();
meshOffset = meshOffset.addNodeField(totApproxErrorPotentialP2Offset,'TotalApproxErrorAccelerationP2Offset');
meshOffset.writeVTK('Eros_OffsetGravMesh.vtk');
