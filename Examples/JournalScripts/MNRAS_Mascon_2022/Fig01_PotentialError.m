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
altitudes = logspace(-1,4,numAltitudes);

% meshes
numFacesCoarse    =  2100;  % coarsest test mesh
numMascons        = 10000;  % number of mascons
numFacesTruthMesh = 50000;  % reference mesh

% degrees of quad rules
quadRules = [1,2,3,4];   % mesh degrees we'll test

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

% Gravity Models
%-------------------------------------------------------------------------
truthGravityModel = AnalyticPolyhedralModel(mesh,Mu);

polyhedralModel{1} = AnalyticPolyhedralModel(meshCoarse,Mu);
masconModel{1} = MasconModel(meshCoarse,Mu,numMascons);

% Acceleration and Error
%--------------------------------------------------------------------------
% create dummy mesh for our constant altitude surfaces
meshCoarse = SurfaceMesh(mesh);
meshCoarse.coarsen(numPointsBaseCAS);

iter = 1;
%progressMeter = progressbar('loop over altitudes :');
for i = 1:numAltitudes
    %progressMeter.update(100*i/numAltitudes)

    % get pts on constant alt surface
    constAltSurface = meshCoarse.offsetSurfaceMesh(altitudes(i)*mesh.resolution, numPointsCAS);
    pts = constAltSurface.coordinates;

    % "true" acceleration
    accOriginal = truthGravityModel.acceleration(pts);
    accOgMag = vecnorm(accOriginal,2,2);
    
    % loop through each mesh resolution calculating the error
    for k = 1:length(polyhedralModel)
        accTempPoly = polyhedralModel{k}.acceleration(pts);
        error(iter,1) = 100*mean(vecnorm((accTempPoly - accOriginal),2,2)./accOgMag);
        stderror(iter,1) = std(vecnorm((accTempPoly - accOriginal),2,2)./accOgMag);
        
        for j=1:length(masconModel)
            accTemp = masconModel{j}.acceleration(pts);
            stderror(iter,j+1) = 100*std(vecnorm((accTemp - accOriginal),2,2)./accOgMag);
            stdnumError(iter,j+1) = 100*std(vecnorm((accTemp - accTempPoly),2,2)./vecnorm(accTempPoly,2,2));
            error(iter,j+1) = 100*mean(vecnorm((accTemp - accOriginal),2,2)./accOgMag);
            numError(iter,j+1) = 100*mean(vecnorm((accTemp - accTempPoly),2,2)./vecnorm(accTempPoly,2,2));
            numNodes(iter,j+1) = masconModel{j}.numElements;
        end
        h(iter) = altitudes(i);
        %delta(iter) = res(k);
        %Nf(iter) = numFacesFinest(k);
        iter = iter+1
    end
end

stop
%progressMeter.close();


% Process Data
%--------------------------------------------------------------------------
Nres=length(res);

errorSquareAnalytic = [];
errorSquareP1 = [];
errorSquareP2 = [];
errorSquareP3 = [];
errorSquareP4 = [];

for i=1:length(altitudes)
    i1 = Nres*(i-1)+1;
    i2 = Nres*i;

    errorSquareAnalytic = [errorSquareAnalytic,error(i1:i2,1).*Nf(i1:i2)'];
    errorSquareP1 = [errorSquareP1,error(i1:i2,2).*Nf(i1:i2)'];
    errorSquareP2 = [errorSquareP2,error(i1:i2,3).*Nf(i1:i2)'];
    errorSquareP3 = [errorSquareP3,error(i1:i2,4).*Nf(i1:i2)'];
    errorSquareP4 = [errorSquareP4,error(i1:i2,5).*Nf(i1:i2)'];
end


% Plot
%--------------------------------------------------------------------------

% Fig 9
%-------
FS=16; LW=1; MS=8
figure(1)
hold on
plot(altitudes'/mesh.resolution,mean(errorSquareP1,1),'ko','MarkerFaceColor','m','lineWidth',LW,"markerSize",MS)
plot(altitudes'/mesh.resolution,mean(errorSquareP2,1),'k^','MarkerFaceColor','b','lineWidth',LW,"markerSize",MS)
plot(altitudes'/mesh.resolution,mean(errorSquareP3,1),'ks','MarkerFaceColor','r','lineWidth',LW,"markerSize",MS)
plot(altitudes'/mesh.resolution,mean(errorSquareP4,1),'kp','MarkerFaceColor','y','lineWidth',LW,"markerSize",MS)
plot(altitudes'/mesh.resolution,mean(errorSquareAnalytic,1),'kd','MarkerFaceColor','k','lineWidth',LW,"markerSize",MS)
legend({'P1','P2','P3','P4','Analytic'},'interpreter','latex','FontSize',FS)
xlabel('$h/\delta_0$','interpreter','latex','FontSize',FS)
ylabel('$\epsilon\cdot N_f$','interpreter','latex','FontSize',FS)
set(gcf,'Color',[1,1,1]);
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')
set(gca,'xscale','log'); set(gca,'yscale','log');
xticks([0.1,1,10,100,1000,10000])
yticks(10.^[-2,0,2,4,5])
%axis([0.1 10000,0.99e-2,1e5])
legend boxoff
box on

% Fig 10 (entry for one body)
%-----------------------------
FS=16; LW=1; MS=8;
figure(4)
hold on
i1 = i;
Nf(i1)
plot(altitudes'/mesh.resolution,mean(errorSquareAnalytic./errorSquareP2,1),'ko','MarkerFaceColor','r','lineWidth',LW,"markerSize",MS)
legend({'Eros','67P','Phobos','Bennu'},'interpreter','latex','FontSize',FS)
xlabel('$h/\delta_0$','interpreter','latex','FontSize',FS)
ylabel('$\mathrm{error\ ratio\ }\frac{\mathrm{analytic}}{\mathrm{quadrature}}$','interpreter','latex','FontSize',FS)
set(gcf,'Color',[1,1,1]);
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')
set(gca,'xscale','log'); set(gca,'yscale','log');
%axis([0.1,10000,0.1,1000])
legend boxoff
box on


