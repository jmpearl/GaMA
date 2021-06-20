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

meshCoarseP3 = meshCoarse.coarsen(888);
meshCoarseP3 = meshCoarseP3.smooth(5);
meshCoarseP3 = meshCoarseP3.setDegree(3);   % make a P2 mesh
meshCoarseP3 = meshCoarseP3.curve(mesh);    % curve mesh through projection

meshCoarseP4 = meshCoarse.coarsen(500);
meshCoarseP4 = meshCoarseP4.smooth(5);
meshCoarseP4 = meshCoarseP4.setDegree(4);   % make a P2 mesh
meshCoarseP4 = meshCoarseP4.curve(mesh);    % curve mesh through projection

% initialize our gravity models on the coarsened mesh
quadratureModel   = ApproximatePolyhedralModel(meshCoarse,Mu,'L1');
quadratureModelP2 = ApproximatePolyhedralModel(meshCoarseP2,Mu);
quadratureModelP3 = ApproximatePolyhedralModel(meshCoarseP3,Mu);
quadratureModelP4 = ApproximatePolyhedralModel(meshCoarseP4,Mu);

analyticModelTruth = AnalyticPolyhedralModel(mesh,Mu);

% calculate grav field and error on surface
pts = mesh.coordinates; 

truthAcceleration = analyticModelTruth.acceleration(pts);
truthPotential = analyticModelTruth.potential(pts);

approximatePotential = quadratureModel.potential(pts);
approximatePotentialP2 = quadratureModelP2.potential(pts);
approximatePotentialP3 = quadratureModelP3.potential(pts);
approximatePotentialP4 = quadratureModelP4.potential(pts);

                    
totApproxErrorPotential = 100*(approximatePotential-truthPotential)./truthPotential;         
totApproxErrorPotentialP2 = 100*(approximatePotentialP2-truthPotential)./truthPotential;
totApproxErrorPotentialP3 = 100*(approximatePotentialP3-truthPotential)./truthPotential;
totApproxErrorPotentialP4 = 100*(approximatePotentialP4-truthPotential)./truthPotential;

% load the fields in and export for viewing w/ paraview
mesh = mesh.clearFields();
mesh = mesh.addNodeField(totApproxErrorPotential,'TotalApproxErrorPotential');
mesh = mesh.addNodeField(totApproxErrorPotentialP2,'TotalApproxErrorPotentialP2');
mesh = mesh.addNodeField(totApproxErrorPotentialP3,'TotalApproxErrorPotentialP3');
mesh = mesh.addNodeField(totApproxErrorPotentialP4,'TotalApproxErrorPotentialP4');

mesh.writeVTK('Eros_GravFields.vtk');
meshCoarse.writeVTK('Eros_CoarseMesh.vtk');

