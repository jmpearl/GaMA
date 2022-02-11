%==========================================================================
% Script used to create Figure 9 and 10
% Pearl J.M., Hitt, D.L.
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
altitudes = logspace(-1,2,numAltitudes);
insetDistances = linspace(0,1,20);

% meshes
numFacesCoarse    =  5000;  % coarsest test mesh
numMascons        =  5000;  % number of mascons
numFacesTruthMesh = 50000;  % reference mesh

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

ds = (mesh.volume/numMascons)^(1/3)
for i=1:length(insetDistances)
    i
    insetSurfaceMesh = meshCoarse.offsetSurfaceMesh(-insetDistances(i)*ds, ...
                                                    meshCoarse.numVertices);
    masconModel{i} = MasconModel(insetSurfaceMesh,Mu,numMascons);

end

% Acceleration and Error
%--------------------------------------------------------------------------
% create dummy mesh for our constant altitude surfaces
meshCoarse = SurfaceMesh(mesh);
meshCoarse.coarsen(numPointsBaseCAS);

iter = 1;
%progressMeter = progressbar('loop over altitudes :');
for i = 1:numAltitudes


    % get pts on constant alt surface
    constAltSurface = meshCoarse.offsetSurfaceMesh(altitudes(i)*mesh.resolution, numPointsCAS);
    pts = constAltSurface.coordinates;

    % "true" acceleration
    accOriginal = truthGravityModel.acceleration(pts);
    accOgMag = vecnorm(accOriginal,2,2);
    
    % loop through each mesh resolution calculating the error
    for k = 1:length(polyhedralModel)
        accTempPoly = polyhedralModel{k}.acceleration(pts);
        maxerror(iter,1) = 100*max(vecnorm((accTempPoly - accOriginal),2,2)./accOgMag);
        error(iter,1) = 100*mean(vecnorm((accTempPoly - accOriginal),2,2)./accOgMag);
        stderror(iter,1) = std(vecnorm((accTempPoly - accOriginal),2,2)./accOgMag);
        
        for j=1:length(masconModel)
            accTemp = masconModel{j}.acceleration(pts);
            stderror(iter,j+1) = 100*std(vecnorm((accTemp - accOriginal),2,2)./accOgMag);
            stdnumError(iter,j+1) = 100*std(vecnorm((accTemp - accTempPoly),2,2)./vecnorm(accTempPoly,2,2));
            error(iter,j+1) = 100*mean(vecnorm((accTemp - accOriginal),2,2)./accOgMag);
            maxerror(iter,j+1) = 100*max(vecnorm((accTemp - accOriginal),2,2)./accOgMag);
            numError(iter,j+1) = 100*mean(vecnorm((accTemp - accTempPoly),2,2)./vecnorm(accTempPoly,2,2));
            numNodes(iter,j+1) = masconModel{j}.numElements;
        end
        h(iter) = altitudes(i);

        iter = iter+1
    end
end



figure(1)
LW = 2;
FS = 12;
hold on
plot(insetDistances,error(5,2:end),'k.-','LineWidth',LW)
xlabel('inset depth / lattice spacing','interpreter','latex','FontSize',FS)
ylabel('acceleration error $\%$','interpreter','latex','FontSize',FS)
set(gcf,'Color',[1,1,1]);
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')
box on
