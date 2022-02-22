%==========================================================================
% Script used to create Figure 5
% Pearl J.M., Hitt, D.L., "Cutting Corners: Curvilinear Surface-Based
% Gravity Models for Asteroids and Comets"
% 2022 (submitted).
%
% ** Note GaMA is under development small variations from published 
% ** results may occur.
%
% writes vtk file with percent error of the gravitational acceleration
% field with different models on the surface of eros.
%==========================================================================

clear all; clc; close all;

% body parameters
G = 6.67*10^-11;                       % gravitational constant (N m2/kg2)
load('Eros.mat');                      % load stored Eros properties                                  
Mu = bodyProperties.mass*G;            % set stand grav parameter
Omega = bodyProperties.rotationRate;   % 

% Surface Meshes
%--------------------------------------------------------------------------
meshfile = 'Eros_46906.obj';               % mesh file to load 
mesh = SurfaceMesh(meshfile);              % create the surface mesh object

meshCoarseP1 = SurfaceMesh(mesh);          % copy construct
meshCoarseP1.coarsen(8000);                % create coarsened version
meshCoarseP1.smooth(5);                    % smooth things out
meshCoarseP1.projectOnTo(mesh);

meshCoarseP2 = SurfaceMesh(meshCoarseP1);
meshCoarseP2.coarsen(2000);
meshCoarseP2.smooth(5);
meshCoarseP2.setDegree(2);                 % make a P2 mesh
meshCoarseP2.curve(mesh);                  % curve mesh through projection

% Gravity Models
%--------------------------------------------------------------------------
quadratureModelP1 = ApproximatePolyhedralModel(meshCoarseP2,Mu);
quadratureModelP2 = ApproximatePolyhedralModel(meshCoarseP2,Mu);

analyticModelTruth = AnalyticPolyhedralModel(mesh,Mu);

% Calculate Acceleration Error on the Surface
%--------------------------------------------------------------------------
pts = mesh.coordinates;

approximateAccelerationP1 = quadratureModelP1.acceleration(pts);
approximateAccelerationP2 = quadratureModelP2.acceleration(pts);
truthAcceleration = analyticModelTruth.acceleration(pts);

totErrorP1 = 100*vecnorm(approximateAccelerationP1-truthAcceleration,2,2)./vecnorm(truthAcceleration,2,2);
totErrorP2 = 100*vecnorm(approximateAccelerationP2-truthAcceleration,2,2)./vecnorm(truthAcceleration,2,2);


% load the fields in and export for viewing w/ paraview
%--------------------------------------------------------------------------
mesh.clearFields();
mesh.addNodeField(totErrorP1,'TotalApproxErrorAccelerationP2');
mesh.addNodeField(totErrorP2,'TotalApproxErrorAccelerationP2');
mesh.writeVTK('ErosSurfaceAccelerationError.vtk');


