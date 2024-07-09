clear all; close all; clc;
addpath(genpath('/home/jmpearl/GravityModels/CLEO'))

Mu = 1.0;

meshfile = 'Eros_46906.obj'; % mesh file to load
sm = SurfaceMesh(meshfile);
sm.coarsen(15000);

vm = VolumeMesh(sm);
vm.initializeFromSurfaceIteration(1/3)

masconModel = MasconModel();
masconModel.initializeFromVolumeMesh(vm,Mu,"excludesurface")

vm.addCellField(vm.cellVolumes(),'cellVolumes')
vm.writeVTK("Fig4_volumeMesh.vtk")


% Fig 10 (entry for one body)
%-----------------------------
FS=16; LW=1; MS=8;
thresh = 2000;
markerScaleFactor = 24;
figure(4)
hold on

% plot faces in certain range x
c = sm.faceCentroids();
for i = 1:sm.numFaces
    if (c(i,3) < thresh && c(i,3) > -thresh)
        pts = [sm.coordinates(sm.faces(i,:),:);...
               sm.coordinates(sm.faces(i,1),:)];
        plot3(pts(:,1),pts(:,2),pts(:,3),'k-')
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
    if (c(i,3) < thresh && c(i,3) > -thresh)
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
view([0,0,1])
box off
axis off

%[f1, cdata] = myaa([4, 2]); imwrite(cdata, [fileName], 'png');