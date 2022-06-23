%==========================================================================
% Pearl J.M., Hitt, D.L.
%
% ** Note GaMA is under development small variations from published 
% ** results may occur if using the latest version
%==========================================================================

clear all; clc; close all;

% Constants
%--------------------------------------------------------------------------
% meshes
numFacesSampleMesh = 20000;
numRefinements = 15;
numFacesCoarse = round(logspace(2,4.3,numRefinements));  % coarsest test mesh

numAltitudes = 10;
altitudes = logspace(1.2460,4.2460,numAltitudes);

% body parameters
%--------------------------------------------------------------------------
G = 6.67*10^-11;                       % gravitational constant (N m2/kg2)
load('Eros.mat');                      % load stored Eros properties
Mu = bodyProperties.mass*G;            % set stand grav parameter

% meshes
%--------------------------------------------------------------------------
meshfile = 'Eros_46906.obj';           % mesh file to load
mesh = SurfaceMesh(meshfile);          % create the surface mesh object

meshSample = SurfaceMesh(mesh);
meshSample.setNumFaces(numFacesSampleMesh)

truthGravityModel = AnalyticPolyhedralModel(mesh,Mu);

% data sets and truth accelerations
%--------------------------------------------------------------------------
for j = 1:numAltitudes
    if j == 1
        pts{j} = meshSample.coordinates;
    else
        meshSampleNew = meshSample.offsetSurfaceMesh(altitudes(j),numFacesSampleMesh);
        pts{j} = meshSampleNew.coordinates;
    end
    tic;
    accOriginal{j} = truthGravityModel.acceleration(pts{j}); toc
    accOgMag{j} = vecnorm(accOriginal{j},2,2);
    numSamplesActual(j) = size(pts{j},1);
end


% do the dead
%--------------------------------------------------------------------------
for i = 1:numRefinements

    % set up our models
    %----------------------------------------------------------------------
    Ni = numFacesCoarse(i);            % number of mascons/faces
    meshCoarse = SurfaceMesh(mesh);    % copy  construct
    meshCoarse.setNumFaces(Ni);        % coarsen to

    polyhedralModel = AnalyticPolyhedralModel(meshCoarse,Mu);

    % storage function for use across scripts (outputs struct)
    masconModel = generateMasconModels(mesh,meshCoarse,Mu,Ni);

    % get our data points
    %----------------------------------------------------------------------
    for k = 1:numAltitudes

        % Acceleration coarse polyhedral model
        %------------------------------------------------------------------
        tic
        accTempPoly = polyhedralModel.acceleration(pts{k});
        time{k}(i,1)=toc
        numElements{k}(i,1)=meshCoarse.numFaces;

        % polyhedral model error norms
        L1{k}(i,1) = 100*mean(vecnorm((accTempPoly - accOriginal{k}),2,2)./accOgMag{k});
        L2{k}(i,1) = 100/numSamplesActual(k)*sqrt(sum((vecnorm((accTempPoly - accOriginal{k}),2,2)./accOgMag{k}).^2));
        Linf{k}(i,1) = 100*max(vecnorm((accTempPoly - accOriginal{k}),2,2)./accOgMag{k});

        % mascon models
        %------------------------------------------------------------------
        for j=1:length(masconModel)
            [i,k,j]

            % mascon model acceleration
            tic
            accTemp = masconModel{j}.acceleration(pts{k});
            time{k}(i,j+1)=toc

            numElements{k}(i,j+1)=masconModel{j}.numElements;

            % mascon model error norms
            L1{k}(i,j+1) = 100*mean(vecnorm((accTemp - accOriginal{k}),2,2)./accOgMag{k});
            L2{k}(i,j+1) = 100/numSamplesActual(k)*sqrt(sum((vecnorm((accTemp - accOriginal{k}),2,2)./accOgMag{k}).^2));
            Linf{k}(i,j+1) = 100*max(vecnorm((accTemp - accOriginal{k}),2,2)./accOgMag{k});

        end
    end

end

if strcmp(body,'Eros')
    save('Eros_46906_Mascon.mat');
elseif strcmp(body,'Bennu')
    save('Bennu_200000_Mascon.mat');
else
    error('invalid body')
end

