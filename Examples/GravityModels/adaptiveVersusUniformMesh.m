%==========================================================================
% Gravitation and Mesh Adaptation Library
% J.M. Pearl
% Dec 4 2021
%==========================================================================
% compares the accuracy of the polyhedral gravity model using feature 
% adapted and unform meshes. The feature-based is actaully
% worse.
%--------------------------------------------------------------------------
clear all; close all; clc;
addpath(genpath('/home/jmpearl/GravityModels/CLEO'))

% standard gravitational parameter
Mu=1;

% options
isCentered=true;

% plot options
markers = ["ko","ks","k^","k>","kv","k<","kp","kh"];
mfc = ['m','r','y','g','c','b','k'];

% base mesh file
surfaceMeshName = 'Eros_46906.obj';
%surfaceMeshName = 'Itokawa_50000.obj';

% Set our meshes
smOriginal = SurfaceMesh(surfaceMeshName);
smFeature  = SurfaceMesh(smOriginal);
smUniform  = SurfaceMesh(smOriginal);
smUniform.coarsenOptions.method = 'uniform';

% point sets to assess error 
ptsSurface = smOriginal.coordinates(randi(smOriginal.numVertices,100,1));
offsetDistance = max(vecnorm(ptsSurface,2,2));
offsetMesh = smOriginal.offsetSurfaceMesh(offsetDistance, 200);
ptsBrillouinSphere = offsetMesh.coordinates;

% loop through and deterine accuracy
i=1;
numFaces = round(smOriginal.numFaces*1/2); % start at 1/2 the original res
while numFaces(end) > 10000
    
    % coarsen
    numFaces(i) = numFaces(end) - 2000;
    smFeature.setNumFaces(numFaces(i));
    smUniform.setNumFaces(numFaces(i));

    % create our curved meshes
    smFeatureP2 = SurfaceMesh(smFeature);
    smFeatureP2.setDegree(2)
    smFeatureP2.curve(smOriginal);
    smUniformP2 = SurfaceMesh(smUniform);
    smUniformP2.setDegree(2)
    smUniformP2.curve(smOriginal);

    % center the meshes if requested
    if isCentered
        smOriginal.center();
        smFeature.center();
        smUniform.center();
        smFeatureP2.center();
        smUniformP2.center();
    end

    % store our refs to the surface meshes
    sm = [smFeature,smUniform,smFeatureP2,smUniformP2];

    % volume & center of mass error
    for j = 1:length(sm)
        volumeError(i,j) = 100*(sm(j).volume-smOriginal.volume)/smOriginal.volume;
        comError(i,j) = 100*norm(sm(j).centroid-smOriginal.centroid);
    end

    % gravity models 
    gm{1} = AnalyticPolyhedralModel(smOriginal,Mu);
    gm{2} = AnalyticPolyhedralModel(smFeature,Mu);
    gm{3} = AnalyticPolyhedralModel(smUniform,Mu);
    gm{4} = ApproximatePolyhedralModel(smFeature,Mu);
    gm{5} = ApproximatePolyhedralModel(smUniform,Mu);
    gm{6} = ApproximatePolyhedralModel(smFeatureP2,Mu);
    gm{7} = ApproximatePolyhedralModel(smUniformP2,Mu);

    % acceleration from diff grav models
    for j = 1:length(gm)
        accBrillouin{j} = gm{j}.acceleration(ptsBrillouinSphere);
        accSurface{j} = gm{j}.acceleration(ptsSurface);
    end

    % error rel to original mesh analytic polyhedral model
    for j = 1:length(gm)
        accSurfaceError(i,j) = 100*mean(vecnorm(accSurface{j}-accSurface{1},2,2)./...
                                        vecnorm(accSurface{1},2,2));
        accBrillouinError(i,j) = 100*mean(vecnorm(accBrillouin{j}-accBrillouin{1},2,2)./...
                                          vecnorm(accBrillouin{1},2,2));
    end

    % store our mesh resolution
    resolution(i) = smUniform.resolution;
    i=i+1;

end


figure(1)
hold on
for i = 1:length(sm)
    plot(numFaces,volumeError(:,i),markers(i),'MarkerFaceColor',mfc(i))
end
set(gcf,'Color',[1,1,1])
xlabel('number of faces','Interpreter','latex')
ylabel('volume error','Interpreter','latex')
legend('feature-p1','uniform-p1','feature-p2','uniform-p2')
figure(4)
hold on
for i = 1:length(sm)
    plot(numFaces,comError(:,i),markers(i),'MarkerFaceColor',mfc(i))
end
set(gcf,'Color',[1,1,1])
xlabel('number of faces','Interpreter','latex')
ylabel('com error(m)','Interpreter','latex')
legend('feature-p1','uniform-p1','feature-p2','uniform-p2')
figure(2)
hold on
for i = 1:size(accSurfaceError,2)
    plot(numFaces,accSurfaceError(:,i),markers(i),'MarkerFaceColor',mfc(i))
end
set(gcf,'Color',[1,1,1])
xlabel('number of faces','Interpreter','latex')
ylabel('acceleration error %','Interpreter','latex')
legend('Original','Analytic-feature','Analytic-uniform','ApproxP1-feature','ApproxP1-uniform','ApproxP2-feature','ApproxP2-uniform')

figure(3)
hold on
for i = 1:size(accSurfaceError,2)
    plot(numFaces,accBrillouinError(:,i),markers(i),'MarkerFaceColor',mfc(i))
end
set(gcf,'Color',[1,1,1])
title('at altitude equal to')
xlabel('number of faces','Interpreter','latex')
ylabel('acceleration error %','Interpreter','latex')
legend('Original','Analytic-feature','Analytic-uniform','ApproxP1-feature','ApproxP1-uniform','ApproxP2-feature','ApproxP2-uniform')