clear all; close all; clc;

fileName = 'Didymos_6m.obj';%
%fileName = 'Didymos_25cm_to_4m_GradedMesh.obj';%

% constants (units kg m s)
G = 6.67e-11;
M = 4.8e9;
T = 11.9 * 3600;

omega = -2*pi/T;
Mu = G*M;

intergratorTolerance = 1e-7;

didyM = 540e9;
didyR = 1.2e3; 
didyA = G*didyM/didyR^2;

Rcom = didyR*didyM / (didyM + M);

omega = -sqrt(G*(M+didyM)/Rcom)/Rcom;
% create our surface mesh
surf = SurfaceMesh(fileName);
%surf.setNumFaces(1000)
surf.coordinates = surf.coordinates * 1000;
surf.coordinates=surf.coordinates -surf.centroid + [Rcom,0.0,0.0];
surf.resetBulkProperties();

% create the gravity model
%--------------------------------------------------------------------------

dimorphosGrav = ApproximatePolyhedralModel(surf,Mu);

didymosGrav = MasconModel(surf, G*didyM, 1);
didymosGrav.coordinates = [-didyR+Rcom,0,0];

grav = CompositeModel(dimorphosGrav,...
                      didymosGrav);


x = surf.coordinates;

rotTermx = omega^2.*x(:,1);
rotTermy = omega^2.*x(:,2);
rotTermPot = 1/2*omega^2*vecnorm(x,2,2).^2;

RotAcc = [rotTermx,rotTermy,rotTermx*0];

potential = grav.potential(x)+rotTermPot;
acceleration = grav.acceleration(x);

TotAcc = acceleration + RotAcc;
accMag = vecnorm(TotAcc,2,2);
accelerationDirection = TotAcc./accMag;
normals = surf.nodeNormals();

slope = acosd(dot(-accelerationDirection,normals,2));


accDidy = didymosGrav.acceleration(x);
accDimo = dimorphosGrav.acceleration(x);

altitude = vecnorm(surf.coordinates-surf.centroid,2,2);
surf.addNodeField(altitude,'altitude');
surf.addNodeField(accMag,'accMag');
surf.addNodeField(potential,'pot');
surf.addNodeField(slope,'slope');
surf.addNodeField(TotAcc,'acceleration');
surf.addNodeField(accDidy,'acceleration_didymos');
surf.addNodeField(accDidy+RotAcc,'acceleration_NOTdimorphos');
surf.addNodeField(accDimo,'acceleration_dimorphos');
surf.addNodeField(RotAcc,'acceleration_RotPsuedo');
surf.writeVTK('dimporphos.vtk');