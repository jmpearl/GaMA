%==========================================================================
% Pearl J.M., Hitt, D.L.
%
% ** Note GaMA is under development small variations from published 
% ** results may occur.
%
% plots our distributions of mascons. Uncomment the distribution you want
% to get it to plot. For the extended tetrahedral distribution, ln. 24 
% setNumFaces should be 500 instead of 2000.
%==========================================================================

clear all; clc; close all;

% Constants
%--------------------------------------------------------------------------
Mu = 1;                     % dummy value for standard grav parameter
numMascons =  5000;         % number of mascons
numFacesTruthMesh = 50000;  % reference mesh
                            

% meshes 
%--------------------------------------------------------------------------
meshfile = 'Bennu_2692.obj';           % mesh file to load 
mesh = SurfaceMesh(meshfile);          % create the surface mesh object
mesh.setNumFaces(500);

meshInset = mesh.offsetSurfaceMesh(-0.02,2000);


% packing simple lattice
%---------------------------------------------------------------------------
 markerScaleFactor = 16;
 fileName = 'distSimple.png';
masconModel = MasconModel();
masconModel.initializeSimplePacking(meshInset,0.04,Mu)


% Plot
%--------------------------------------------------------------------------

% Fig 10 (entry for one body)
%-----------------------------
FS=16; LW=1; MS=8;
thresh = 0.05;
figure(4)
hold on

% plot faces in certain range x
c = mesh.faceCentroids();
for i = 1:mesh.numFaces
    if (c(i,1) < thresh && c(i,1) > -thresh)
        pts = [mesh.coordinates(mesh.faces(i,:),:);...
               mesh.coordinates(mesh.faces(i,1),:)];
        plot3(pts(:,1),pts(:,2),pts(:,3),'k-')
    end
end

% plot insetfaces in certain range x
c = meshInset.faceCentroids();
for i = 1:meshInset.numFaces
    if (c(i,1) < thresh && c(i,1) > -thresh)
        pts = [meshInset.coordinates(mesh.faces(i,:),:);...
               meshInset.coordinates(mesh.faces(i,1),:)];
        plot3(pts(:,1),pts(:,2),pts(:,3),'r-')
    end
end

c = masconModel.coordinates;
mu = masconModel.mu;

resolution = mu.^(1/3);
MS = resolution./max(resolution)*markerScaleFactor;
color = 1.0-(MS-min(MS))/max(max(MS)-min(MS),mean(MS)/10000);
colorVecMax = [1,1,0];
colorVecMin = [0,0,1];
for i = 1:masconModel.numElements
    if (c(i,1) < thresh && c(i,1) > -thresh)
        colorVeci = (1-color(i))*colorVecMax + color(i)*colorVecMin;
        MSi = MS(i);
        plot3(masconModel.coordinates(i,1),...
              masconModel.coordinates(i,2),...
              masconModel.coordinates(i,3),'ko','MarkerFaceColor',colorVeci,'MarkerSize',MSi)
    end
end
set(gcf,'Color',[1,1,1]);
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')
daspect([1,1,1])
view([1,0,0])
box off
axis off

[f1, cdata] = myaa([4, 2]); imwrite(cdata, 'insetMesh.png', 'png');
