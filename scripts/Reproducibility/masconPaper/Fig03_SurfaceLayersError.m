%==========================================================================
% Script used to create Figure 9 and 10
% Pearl J.M., Hitt, D.L., "Cutting Corners: A Quadrature-based Gravity
% Model for Comets and Asteroids using Curvilinear Surface Definitions"
% MNRAS, 2022 (submitted).
%
% ** Note GaMA is under development small variations from published 
% ** results may occur.
%
% plots acceleration error of quadrature and analytic models as a function
% of normalized altitude
%==========================================================================

clear all; clc; close all;

% Constants
%--------------------------------------------------------------------------
% constant altitude surface parameters
numPointsCAS = 500;       % number of points on CAS
numPointsBaseCAS = 4000;  % number of points on origin mesh
numAltitudes = 5;
altitudes = logspace(-3,4,numAltitudes);

% meshes
numFacesCoarse    =  5000;  % coarsest test mesh
numMascons        =  5000;  % number of mascons
numFacesTruthMesh = 50000;  % reference mesh

% names for our fields we'll dump to VTK
fieldnames = {'AnalyticPolyhedron',...
              'MasconDegree0Lattice',...
              'MasconDegree1LatticeVertex',...
              'MasconDegree2LatticeVertex',...
              'MasconDegree1LatticeCell',...
              'MasconDegree1LatticeNode',...
              'MasconDegree1OctreeVertex'};

% body parameters
%--------------------------------------------------------------------------
G = 6.67*10^-11;                       % gravitational constant (N m2/kg2)
load('Eros.mat');                      % load stored Eros properties                                  
Mu = bodyProperties.mass*G;            % set stand grav parameter
Omega = bodyProperties.rotationRate;   % 

% meshes 
%--------------------------------------------------------------------------
meshfile = 'Eros_46906.obj';           % mesh file to load 
mesh = SurfaceMesh(meshfile);          % create the surface mesh object
mesh.setNumFaces(numFacesTruthMesh);

% initialize our course mesh we'll manipulate
meshCoarse = SurfaceMesh(mesh);           % copy  construct
meshCoarse.setNumFaces(numFacesCoarse);   % coarsen to
sm = SurfaceMesh(meshCoarse);
volMesh = VolumeMesh(meshCoarse);
volMesh.initializeFromSimpleLattice(numMascons);
volMesh.smooth(5)


% Gravity Models
%-------------------------------------------------------------------------
% fine analytic polyhedron
truthGravityModel = AnalyticPolyhedralModel(mesh,Mu);

% coarse analytic polyhedron
polyhedralModel = AnalyticPolyhedralModel(meshCoarse,Mu);

% mascon degree 2 lattice
masconModel{1} = MasconModel(meshCoarse,Mu,numMascons);

% mascon degree 1 lattice
masconModel{2} = MasconModel();
masconModel{2}.initializeFromVolumeMesh(volMesh,Mu,'vertex');

% mascon degree 2 lattice
volMesh.setDegree(2);
volMesh.curve(mesh);

masconModel{3} = MasconModel();
masconModel{3}.initializeFromVolumeMesh(volMesh,Mu,'vertex');

masconModel{4} = MasconModel();
masconModel{4}.initializeFromVolumeMesh(volMesh,Mu,'cell');

masconModel{5} = MasconModel();
masconModel{5}.initializeFromVolumeMesh(volMesh,Mu,'node');

% mascon degree 2 octree

volMesh2 = VolumeMesh(sm);
volMesh2.initializeFromOctree(1000,2,2);
volMesh2.smooth(1);
volMesh2.setDegree(2);
volMesh2.curve(mesh);

masconModel{6} = MasconModel();
masconModel{6}.initializeFromVolumeMesh(volMesh2,Mu,'vertex');


% Acceleration and Error
%--------------------------------------------------------------------------
pts = mesh.coordinates; 

% "true" potential/acceleration
potOriginal = truthGravityModel.potential(pts);
accOriginal = truthGravityModel.acceleration(pts);
accOgMag = vecnorm(accOriginal,2,2);

% analytic coarse polyhedron potential/acceleration 
potTempPoly = polyhedralModel.potential(pts);
accTempPoly = polyhedralModel.acceleration(pts);
potError = 100*abs((potTempPoly - potOriginal)./potOriginal);
accError = 100*vecnorm((accTempPoly - accOriginal),2,2)./accOgMag;

% mascon models 
for j=1:length(masconModel)
    potTemp = masconModel{j}.potential(pts);
    accTemp = masconModel{j}.acceleration(pts);
    potError = [potError,100*abs((potTemp - potOriginal)./potOriginal)];
    accError = [accError,100*vecnorm((accTemp - accOriginal),2,2)./accOgMag];
end

% add them to our surface mesh 
for j = 1:size(potError,2)

    mesh.addNodeField(potError(:,j),['Potential',fieldnames{j}])
    mesh.addNodeField(accError(:,j),['Acceleration',fieldnames{j}])

end

% export our vtk file to visualize with paraview
mesh.writeVTK('masconSurfaceError.vtk')

