%--------------------------------------------------------------------------
% initializes the 4 gravity models, calculates field values and prints 
% out the results for the first field point
%--------------------------------------------------------------------------

clear all; close all; clc


% body parameters
G = 6.67*10^-11;                       % gravitational constant (N m2/kg2)
load('Eros.mat');                      % load stored Eros properties                                  
Mu = bodyProperties.mass*G;            % set stand grav parameter
Omega = bodyProperties.rotationRate;   % 

% some model parameters
meshfile = 'Eros_7624.obj';    % mesh file to load 
mesh = SurfaceMesh(meshfile);  % create the surface mesh object
mesh.coarsen(3000);            % coarsen it to 3k faces
N = 6;                         % Order of SH expansion
numMascons = 4000;             % number of mascons 

% location of calculation point
Nsample = 1;
P = [rand(Nsample,1),rand(Nsample,1),rand(Nsample,1)];
P = P./vecnorm(P,2,2)*20000;

% initialize models 
harmonicModel = SphericalHarmonicModel(mesh,Mu,N);
quadratureModel = ApproximatePolyhedralModel(mesh,Mu,'B2');
polyhedralModel = AnalyticPolyhedralModel(mesh,Mu);
masconModel = MasconModel(mesh, Mu, numMascons);

% potential
potH = harmonicModel.potential(P);
potQ = quadratureModel.potential(P);
potP = polyhedralModel.potential(P);
potM = masconModel.potential(P);

% acceleration 
accH = harmonicModel.acceleration(P);
accQ = quadratureModel.acceleration(P);
accP = polyhedralModel.acceleration(P);
accM = masconModel.acceleration(P);

% Laplacian
%LapH = harmonicModel.Laplacian(P); % not implemented yet
LapQ = quadratureModel.Laplacian(P);
LapP = polyhedralModel.Laplacian(P);
LapM = masconModel.Laplacian(P);

% Laplacian @ COM
%Lap0H = harmonicModel.Laplacian([0,0,0]); % not implemented yet
Lap0Q = quadratureModel.Laplacian([0,0,0]);
Lap0P = polyhedralModel.Laplacian([0,0,0]);
Lap0M = masconModel.Laplacian([0,0,0]);

% Grav Gradient
%Lap0H = harmonicModel.Laplacian([0,0,0]);  % not implemented yet
gravGradQ = quadratureModel.gravityGradient(P);
gravGradP = polyhedralModel.gravityGradient(P);
gravGradM = masconModel.gravityGradient(P);

disp('================================================================')
disp([' Testing GravityModels -- Geometry: ',meshfile])
disp('----------------------------------------------------------------')
disp('  SphericalHarmonicModel')
disp('  AnalyticPolyhedralModel')
disp('  ApproximatePolyhedralModel')
disp('  MasconModel')
disp('================================================================')
disp(' ')
disp('Approximate Polyhedral Model')
disp(' ')
disp(['Potential    =  ',num2str(potQ(1,1),'%.8e')])
disp(' ')
disp(['Acceleration =  ',num2str(accQ(1,1:3),'%.8e')])
disp(' ')
disp(['GravGradient = |',num2str(gravGradQ(1,[1,2,3]),'%.8e'),'|'])
disp(['               |',num2str(gravGradQ(1,[2,4,5]),'%.8e'),'|'])
disp(['               |',num2str(gravGradQ(1,[3,5,6]),'%.8e'),'|'])
disp(' ')
disp(['Laplacian(P) =  ',num2str(LapQ(1,1),'%.8e')])
disp(['Laplacian(0) =  ',num2str(Lap0Q(1,1),'%.8e')])
disp(' ')
disp('--------------------------------------------------------------------')
disp(' ')
disp('Analytic Polyhedral Model')
disp(' ')
disp(['Potential    =  ',num2str(potP(1,1),'%.8e')])
disp(' ')
disp(['Acceleration =  ',num2str(accP(1,1:3),'%.8e')])
disp(' ')
disp(['GravGradient = |',num2str(gravGradP(1,[1,2,3]),'%.8e'),'|'])
disp(['               |',num2str(gravGradP(1,[2,4,5]),'%.8e'),'|'])
disp(['               |',num2str(gravGradP(1,[3,5,6]),'%.8e'),'|'])
disp(' ')
disp(['Laplacian(P) =  ',num2str(LapP(1,1),'%.8e')])
disp(['Laplacian(0) =  ',num2str(Lap0P(1,1),'%.8e')])
disp(' ')
disp('--------------------------------------------------------------------')
disp(' ')
disp('Mascon Model')
disp(' ')
disp(['Potential    =  ',num2str(potM(1,1),'%.8e')])
disp(' ')
disp(['Acceleration =  ',num2str(accM(1,1:3),'%.8e')])
disp(' ')
disp(['GravGradient = |',num2str(gravGradM(1,[1,2,3]),'%.8e'),'|'])
disp(['               |',num2str(gravGradM(1,[2,4,5]),'%.8e'),'|'])
disp(['               |',num2str(gravGradM(1,[3,5,6]),'%.8e'),'|'])
disp(' ')
disp(['Laplacian(P) =  ',num2str(LapM(1,1),'%.8e')])
disp(['Laplacian(0) =  ',num2str(Lap0M(1,1),'%.8e')])
disp(' ')
disp('--------------------------------------------------------------------')
disp(' ')
disp('Spherical Harmonic Model')
disp(' ')
disp(['Potential    =  ',num2str(potH(1,1),'%.8e')])
disp(' ')
disp(['Acceleration =  ',num2str(accH(1,1:3),'%.8e')])
disp(' ')
% disp(['GravGradient = |',num2str(gravGradM(1,[1,2,3]),'%.8e'),'|'])
% disp(['               |',num2str(gravGradM(1,[2,4,5]),'%.8e'),'|'])
% disp(['               |',num2str(gravGradM(1,[3,5,6]),'%.8e'),'|'])
% disp(' ')
%disp(['Laplacian(P) =  ',num2str(LapH(1,1),'%.8e')])
%disp(['Laplacian(0) =  ',num2str(Lap0H(1,1),'%.8e')])
disp(' ')
disp('--------------------------------------------------------------------')

