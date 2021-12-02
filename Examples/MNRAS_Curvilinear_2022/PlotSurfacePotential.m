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
%analyticModel = AnalyticPolyhedralModel(meshCoarse,Mu);
analyticModelTruth = AnalyticPolyhedralModel(mesh,Mu);

% calculate grav field and error on surface
pts = mesh.coordinates; 

truthAcceleration = analyticModelTruth.acceleration(pts);
truthPotential = analyticModelTruth.potential(pts);

%analyticAcceleration = analyticModel.acceleration(pts);
%analyticPotential = analyticModel.potential(pts);

%approximateAcceleration = quadratureModel.acceleration(pts);
approximatePotential = quadratureModel.potential(pts);

%approximateAccelerationP2 = quadratureModelP2.acceleration(pts);
approximatePotentialP2 = quadratureModelP2.potential(pts);

%approximateAccelerationP3 = quadratureModelP3.acceleration(pts);
approximatePotentialP3 = quadratureModelP3.potential(pts);

%approximateAccelerationP4 = quadratureModelP4.acceleration(pts);
approximatePotentialP4 = quadratureModelP4.potential(pts);

%totAnalyticErrorAcceleration = 100*vecnorm(analyticAcceleration-truthAcceleration,2,2)./...
%                        vecnorm(truthAcceleration,2,2);
                    
%totAnalyticErrorPotential = 100*(analyticPotential-truthPotential)./truthPotential;


%totApproxErrorAcceleration = 100*vecnorm(approximateAcceleration-truthAcceleration,2,2)./...
%                        vecnorm(truthAcceleration,2,2);
                    
totApproxErrorPotential = 100*(approximatePotential-truthPotential)./truthPotential;

%totApproxErrorAccelerationP2 = 100*vecnorm(approximateAccelerationP2-truthAcceleration,2,2)./...
%                        vecnorm(truthAcceleration,2,2);
                    
totApproxErrorPotentialP2 = 100*(approximatePotentialP2-truthPotential)./truthPotential;

totApproxErrorPotentialP3 = 100*(approximatePotentialP3-truthPotential)./truthPotential;
totApproxErrorPotentialP4 = 100*(approximatePotentialP4-truthPotential)./truthPotential;
%numErrorAcceleration = 100*vecnorm(approximateAcceleration-analyticAcceleration,2,2)./...
%                        vecnorm(analyticAcceleration,2,2);
                    
%numErrorPotential = 100*(approximatePotential-analyticPotential)./analyticPotential;

% load the fields in and export for viewing w/ paraview
mesh = mesh.clearFields();
%mesh = mesh.addNodeField(truthAcceleration,'AnalyticAcceleration46906');
%mesh = mesh.addNodeField(truthPotential,'AnalyticPotential46906');
%mesh = mesh.addNodeField(analyticAcceleration,'AnalyticAcceleration1996');
%mesh = mesh.addNodeField(analyticPotential,'AnalyticPotential1996');
%mesh = mesh.addNodeField(approximateAcceleration,'ApproximateAcceleration');
%mesh = mesh.addNodeField(approximatePotential,'ApproximatePotential');
%mesh = mesh.addNodeField(numErrorAcceleration,'NumericalErrorAcceleration');
%mesh = mesh.addNodeField(numErrorPotential,'NumericalErrorPotential');
%mesh = mesh.addNodeField(totApproxErrorAcceleration,'TotalApproxErrorAcceleration');
mesh = mesh.addNodeField(totApproxErrorPotential,'TotalApproxErrorPotential');
%mesh = mesh.addNodeField(totAnalyticErrorAcceleration,'TotalAnalyticErrorAcceleration');
%mesh = mesh.addNodeField(totAnalyticErrorPotential,'TotalAnalyticErrorPotential');
%mesh = mesh.addNodeField(totApproxErrorAccelerationP2,'TotalApproxErrorAccelerationP2');
mesh = mesh.addNodeField(totApproxErrorPotentialP2,'TotalApproxErrorPotentialP2');
mesh = mesh.addNodeField(totApproxErrorPotentialP3,'TotalApproxErrorPotentialP3');
mesh = mesh.addNodeField(totApproxErrorPotentialP4,'TotalApproxErrorPotentialP4');

mesh.writeVTK('Eros_GravFields.vtk');
meshCoarse.writeVTK('Eros_CoarseMesh.vtk');

