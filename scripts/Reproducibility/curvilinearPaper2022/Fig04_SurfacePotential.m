%==========================================================================
% Script used to create Figure 4
% Pearl J.M., Hitt, D.L., "Cutting Corners: Curvilinear Surface-Based
% Gravity Models for Asteroids and Comets"
% 2022 (submitted).
%
% ** Note GaMA is under development small variations from published 
% ** results may occur.
%
% Writes vtk file with percent error of the gravitational acceleration
% field with different models on the surface of eros. 
%==========================================================================

clear all; clc; close all;

% body parameters
%--------------------------------------------------------------------------
G = 6.67*10^-11;                       % gravitational constant (N m2/kg2)
load('Eros.mat');                      % load stored Eros properties                                  
Mu = bodyProperties.mass*G;            % set stand grav parameter
Omega = bodyProperties.rotationRate;   % 

% Create Our Meshes
%--------------------------------------------------------------------------
meshfile = 'Eros_46906.obj';               % mesh file to load 
mesh = SurfaceMesh(meshfile);              % create the surface mesh object

meshCoarseP1 = SurfaceMesh(mesh);            % copy constructr
meshCoarseP1.coarsen(8000);                  % coarsen to 8000 faces
meshCoarseP1.smooth(5);                      % smooth things out
meshCoarseP1.projectOnTo(mesh);              % project back onto original mesh

meshCoarseP2 = SurfaceMesh(meshCoarseP1);    % copy construct
meshCoarseP2.coarsen(2000);                % coarsen to 2000 faces
meshCoarseP2.smooth(5);                    % smooth 5 interations
meshCoarseP2.setDegree(2);                 % make a P2 mesh
meshCoarseP2.curve(mesh);                  % curve mesh through projection

meshCoarseP3 = SurfaceMesh(meshCoarseP1);
meshCoarseP3.coarsen(888);
meshCoarseP3.smooth(5);
meshCoarseP3.setDegree(3);   % make a P3 mesh
meshCoarseP3.curve(mesh);    % curve mesh through projection

meshCoarseP4 = SurfaceMesh(meshCoarseP1);
meshCoarseP4.coarsen(500);
meshCoarseP4.smooth(5);
meshCoarseP4.setDegree(4);   % make a P4 mesh
meshCoarseP4.curve(mesh);    % curve mesh through projection

% gravity models
%--------------------------------------------------------------------------
quadratureModelP1 = ApproximatePolyhedralModel(meshCoarseP1,Mu);
quadratureModelP2 = ApproximatePolyhedralModel(meshCoarseP2,Mu);
quadratureModelP3 = ApproximatePolyhedralModel(meshCoarseP3,Mu);
quadratureModelP4 = ApproximatePolyhedralModel(meshCoarseP4,Mu);

analyticModelTruth = AnalyticPolyhedralModel(mesh,Mu);


% calculate grav field and error on surface
%--------------------------------------------------------------------------
pts = mesh.coordinates; 

truthAcceleration = analyticModelTruth.acceleration(pts);
truthPotential = analyticModelTruth.potential(pts);

approximatePotentialP1 = quadratureModelP1.potential(pts);
approximatePotentialP2 = quadratureModelP2.potential(pts);
approximatePotentialP3 = quadratureModelP3.potential(pts);
approximatePotentialP4 = quadratureModelP4.potential(pts);
           
totApproxErrorPotentialP1 = 100*(approximatePotentialP1-truthPotential)./truthPotential;         
totApproxErrorPotentialP2 = 100*(approximatePotentialP2-truthPotential)./truthPotential;
totApproxErrorPotentialP3 = 100*(approximatePotentialP3-truthPotential)./truthPotential;
totApproxErrorPotentialP4 = 100*(approximatePotentialP4-truthPotential)./truthPotential;


% load the fields in and export for viewing w/ paraview
%--------------------------------------------------------------------------
mesh.clearFields();
mesh.addNodeField(totApproxErrorPotentialP1,'TotalApproxErrorPotentialP1');
mesh.addNodeField(totApproxErrorPotentialP2,'TotalApproxErrorPotentialP2');
mesh.addNodeField(totApproxErrorPotentialP3,'TotalApproxErrorPotentialP3');
mesh.addNodeField(totApproxErrorPotentialP4,'TotalApproxErrorPotentialP4');

mesh.writeVTK('ErosSurfacePotentialError.vtk');
meshCoarseP1.writeVTK('ErosCoarseMesh.vtk');

