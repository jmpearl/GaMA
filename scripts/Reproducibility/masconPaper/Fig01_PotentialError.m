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
numAltitudes = 20;
altitudes = logspace(-1,2,numAltitudes);

% meshes
numFacesCoarse    =  5000;  % coarsest test mesh
numMascons        =  2500;  % number of mascons
numFacesTruthMesh = 50000;  % reference mesh

% body parameters
%--------------------------------------------------------------------------
G = 6.67*10^-11;                       % gravitational constant (N m2/kg2)
load('Eros.mat');                      % load stored Eros properties                                  
Mu = bodyProperties.mass*G;            % set stand grav parameter

% meshes 
%--------------------------------------------------------------------------
meshfile = 'Eros_46906.obj';           % mesh file to load 
mesh = SurfaceMesh(meshfile);          % create the surface mesh object
mesh.setNumFaces(numFacesTruthMesh);

% initialize our course mesh we'll manipulate
meshCoarse = SurfaceMesh(mesh);           % copy  construct
meshCoarse.setNumFaces(numFacesCoarse);   % coarsen to

insetSurfaceMesh = meshCoarse.offsetSurfaceMesh(-meshCoarse.resolution/2, ...
                                                 meshCoarse.numVertices);

meshCAS = SurfaceMesh(mesh);
meshCAS.coarsen(numPointsBaseCAS);

vmCoarse = VolumeMesh(meshCoarse);
vmCoarse.initializeFromSimpleLattice(numMascons);
vmCoarse.smooth(5)

vmFine = VolumeMesh(mesh);
vmFine.initializeFromSimpleLattice(numFacesCoarse);
vmFine.smooth(5);

vmOctreeDegree2 = VolumeMesh(meshCoarse);
vmOctreeDegree2.initializeFromOctree(200,2,2);
vmOctreeDegree2.smooth(5);
vmOctreeDegree2.setDegree(2);
vmOctreeDegree2.curve(mesh);

% Gravity Models
%-------------------------------------------------------------------------

% polyhedral truth
truthGravityModel = AnalyticPolyhedralModel(mesh,Mu);

% polyhedral test mesh
polyhedralModel{1} = AnalyticPolyhedralModel(meshCoarse,Mu);

% mascon - packing
masconModel{1} = MasconModel(insetSurfaceMesh,Mu,numFacesCoarse);

% mascon - volume mesh rectilinear d1 quad
masconModel{2} = MasconModel();
masconModel{2}.initializeFromVolumeMesh(vmCoarse,Mu,'vertex');

% mascon - volume mesh P2 mesh d1, d1, d2, quad
vmCoarse.setDegree(2);
vmCoarse.curve(mesh);
masconModel{3} = MasconModel();
masconModel{3}.initializeFromVolumeMesh(vmCoarse,Mu,'vertex');
masconModel{4} = MasconModel();
masconModel{4}.initializeFromVolumeMesh(vmCoarse,Mu,'cell');
masconModel{5} = MasconModel();
masconModel{5}.initializeFromVolumeMesh(vmCoarse,Mu,'node');

% mascon - degree 2 octree
masconModel{6} = MasconModel();
masconModel{6}.initializeFromVolumeMesh(vmOctreeDegree2,Mu,'vertex');

% based on true mesh excluding surface points
masconModel{7} = MasconModel();
masconModel{7}.initializeFromVolumeMesh(vmFine,Mu,'excludesurface');

% Acceleration and Error
%--------------------------------------------------------------------------
iter = 1;
for i = 1:numAltitudes

    % get pts on constant alt surface
    constAltSurface = meshCAS.offsetSurfaceMesh(altitudes(i)*mesh.resolution, numPointsCAS);
    pts = constAltSurface.coordinates;

    % "true" acceleration
    accOriginal = truthGravityModel.acceleration(pts);
    accOgMag = vecnorm(accOriginal,2,2);
    
    % loop through each mesh resolution calculating the error
    for k = 1:length(polyhedralModel)
        accTempPoly = polyhedralModel{k}.acceleration(pts);
        error(iter,1) = 100*mean(vecnorm((accTempPoly - accOriginal),2,2)./accOgMag);
        maxerror(iter,1) = 100*max(vecnorm((accTempPoly - accOriginal),2,2)./accOgMag);

        for j=1:length(masconModel)
            accTemp = masconModel{j}.acceleration(pts);
            error(iter,j+1) = 100*mean(vecnorm((accTemp - accOriginal),2,2)./accOgMag);
            maxerror(iter,j+1) = 100*max(vecnorm((accTemp - accOriginal),2,2)./accOgMag);

        end
        h(iter) = altitudes(i);
        iter = iter+1
    end
end

stop

FS = 12;
LW = 1;
MS = 8;
figure(1)
hold on
plot(h,error(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
plot(h,error(:,2),'ko-','MarkerFaceColor','r','LineWidth',LW,'MarkerSize',MS)
plot(h,error(:,3),'ks-','MarkerFaceColor','c','LineWidth',LW,'MarkerSize',MS)
plot(h,error(:,4),'ks-','MarkerFaceColor','b','LineWidth',LW,'MarkerSize',MS)
plot(h,error(:,5),'ks-','MarkerFaceColor','g','LineWidth',LW,'MarkerSize',MS)
plot(h,error(:,6),'ks-','MarkerFaceColor','m','LineWidth',LW,'MarkerSize',MS)
plot(h,error(:,8),'ks-','MarkerFaceColor','y','LineWidth',LW,'MarkerSize',MS)
legend({['Analytic Polyhedron $N_f=',num2str(meshCoarse.numFaces),'$'],...
        ['Simple Packing $N_i=',num2str(masconModel{1}.numElements),'$']...
        ['FVM P1Q1 $N_i=',num2str(masconModel{2}.numElements),'$']...
        ['FVM P2Q1 $N_i=',num2str(masconModel{3}.numElements),'$']...
        ['FVM P2Q1 $N_i=',num2str(masconModel{4}.numElements),'$']...
        ['FVM P2Q2 $N_i=',num2str(masconModel{5}.numElements),'$']...
        ['FVM P2Q1 $N_i=',num2str(masconModel{7}.numElements),'$']...
        },'interpreter','latex','FontSize',FS)
ylabel('error $\%$','interpreter','latex','FontSize',FS)
xlabel('altitude','interpreter','latex','FontSize',FS)
set(gcf,'Color',[1,1,1]);
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')
set(gca,'xscale','log'); set(gca,'yscale','log');
%axis([0.1,10000,0.1,1000])
legend boxoff
box on

