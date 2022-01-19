%==========================================================================
% Exports gravity slopes on the surface of eros in degrees for viewing
% with paraview or visit.
%==========================================================================
clear all; close all; clc;


G = 6.67e-11;

% load body data
%--------------------------------------------------------------------------
load('Eros.mat');                                 
Mu = bodyProperties.mass*G;
Omega = bodyProperties.rotationRate;

% set up the surface mesh
%--------------------------------------------------------------------------
meshFile = 'Eros_46906.obj';    % define surface mesh file
mesh = SurfaceMesh(meshFile);   % create tri surface mesh
meshCoarse = SurfaceMesh(mesh); % copy construct 
meshCoarse.coarsen(2000);       % coursen to 2000 faces

% set up the gravity models
%--------------------------------------------------------------------------
gravityModel = AnalyticPolyhedralModel(mesh,Mu);     % Werner 1994

% calc at fine mesh vertices
%--------------------------------------------------------------------------
acceleration = gravityModel.acceleration(mesh.coordinates);
accelerationDirection = acceleration./vecnorm(acceleration,2,2);
normals = mesh.nodeNormals();

slope = acosd(dot(accelerationDirection,normals,2));

% add to the mesh and export
%--------------------------------------------------------------------------
mesh.addNodeField(slope,'Slopes');
mesh.writeVTK('ErosSlopesExample.vtk');