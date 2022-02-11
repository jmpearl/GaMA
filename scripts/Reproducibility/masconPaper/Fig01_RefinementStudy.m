% ==========================================================================
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
% meshes
numFacesSampleMesh = 20000;
numFacesCoarse = round(logspace(2,4.3,15));  % coarsest test mesh
numRefinements = length(numFacesCoarse);

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

pts = meshSample.coordinates;
tic
accOriginal = truthGravityModel.acceleration(pts);
toc
accOgMag = vecnorm(accOriginal,2,2);


for i = 1:numRefinements

    Ni = numFacesCoarse(i);  % number of mascons


    % initialize our course mesh we'll manipulate
    meshCoarse = SurfaceMesh(mesh);           % copy  construct
    meshCoarse.setNumFaces(Ni);   % coarsen to

    insetSurfaceMesh = meshCoarse.offsetSurfaceMesh(-meshCoarse.resolution/2, ...
                                                     meshCoarse.numVertices);

    vmCoarse = VolumeMesh(meshCoarse);
    vmCoarse.initializeFromSimpleLattice(Ni-meshCoarse.numVertices);
    vmCoarse.smooth(5)

    vmFine = VolumeMesh(mesh);
    vmFine.initializeFromSimpleLattice(Ni);
    vmFine.smooth(5);

    vmOctreeDegree2 = VolumeMesh(mesh);
    vmOctreeDegree2.initializeFromOctree(round(Ni/4.0),2,2);
    vmOctreeDegree2.smooth(5);
    vmOctreeDegree2.setDegree(2);
    vmOctreeDegree2.curve(mesh);

    vmIter = VolumeMesh(meshCoarse);
    vmIter.initializeFromSurfaceIteration(1/2);
    vmIter.smooth(2)
    vmIter.setDegree(2);
    vmIter.curve(mesh);

    vmIter2 = VolumeMesh(meshCoarse);
    vmIter2.initializeFromSurfaceIteration(1/3);
    vmIter2.smooth(2)
    vmIter2.setDegree(2);
    vmIter2.curve(mesh);

    vmIter4 = VolumeMesh(mesh);
    vmIter4.initializeFromSurfaceIteration(1/2,meshCoarse.numVertices);
    vmIter4.smooth(2)

    % Gravity Models
    %-------------------------------------------------------------------------
    % polyhedral test mesh
    polyhedralModel = AnalyticPolyhedralModel(meshCoarse,Mu);

    % mascon - packing
    masconModel{1} = MasconModel(insetSurfaceMesh,Mu,Ni);

    % mascon - volume mesh P1 vertex quads
    masconModel{2} = MasconModel();
    masconModel{2}.initializeFromVolumeMesh(vmCoarse,Mu,'vertex');

    % mascon - volume mesh P2 mesh vertex quad
    vmCoarse.setDegree(2);
    vmCoarse.curve(mesh);
    masconModel{3} = MasconModel();
    masconModel{3}.initializeFromVolumeMesh(vmCoarse,Mu,'vertex');
    masconModel{4} = MasconModel();
    masconModel{4}.initializeFromVolumeMesh(vmCoarse,Mu,'cell');
    masconModel{5} = MasconModel();
    masconModel{5}.initializeFromVolumeMesh(vmCoarse,Mu,'node');

    % mascon - degree 2 octree vertex quad
    masconModel{6} = MasconModel();
    masconModel{6}.initializeFromVolumeMesh(vmOctreeDegree2,Mu,'excludesurface');

    % based on true mesh excluding surface points
    masconModel{7} = MasconModel();
    masconModel{7}.initializeFromVolumeMesh(vmFine,Mu,'excludesurface');

    % based on true mesh excluding surface points
    masconModel{8} = MasconModel();
    masconModel{8}.initializeFromVolumeMesh(vmIter,Mu,'excludesurface');

    masconModel{9} = MasconModel();
    masconModel{9}.initializeFromVolumeMesh(vmIter2,Mu,'excludesurface');

    masconModel{10} = MasconModel();
    masconModel{10}.initializeFromVolumeMesh(vmIter4,Mu,'excludesurface');

    % Acceleration and Error
    %--------------------------------------------------------------------------
    tic
    accTempPoly = polyhedralModel.acceleration(pts);
    time(i,1)=toc
    numElements(i,1)=meshCoarse.numFaces;
    % loop through each mesh resolution calculating the error
    error(i,1) = 100*mean(vecnorm((accTempPoly - accOriginal),2,2)./accOgMag);
    maxerror(i,1) = 100*max(vecnorm((accTempPoly - accOriginal),2,2)./accOgMag);

    for j=1:length(masconModel)
        [i,j]
        tic
        accTemp = masconModel{j}.acceleration(pts);
        time(i,j+1)=toc
        numElements(i,j+1)=masconModel{j}.numElements;
        error(i,j+1) = 100*mean(vecnorm((accTemp - accOriginal),2,2)./accOgMag);
        maxerror(i,j+1) = 100*max(vecnorm((accTemp - accOriginal),2,2)./accOgMag);

    end

end

stop

FS = 12;
LW = 1;
MS = 8;
figure(1)
hold on
plot(numElements(:,1),error(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,2),error(:,2),'ko-','MarkerFaceColor','r','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,3),error(:,3),'ks-','MarkerFaceColor','c','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,4),error(:,4),'ks-','MarkerFaceColor','b','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,5),error(:,5),'ks-','MarkerFaceColor','g','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,6),error(:,6),'ks-','MarkerFaceColor','m','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,7),error(:,7),'ks-','MarkerFaceColor',[1,0.5,0.0],'LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,8),error(:,8),'ks-','MarkerFaceColor',[1.0,0.5,1.0],'LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,9),error(:,9),'ks-','MarkerFaceColor',[0.5,0.0,1.0],'LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,10),error(:,10),'ks-','MarkerFaceColor',[1.0,1.0,1.0],'LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,11),error(:,11),'ks-','MarkerFaceColor',[0.9,0.9,0.9],'LineWidth',LW,'MarkerSize',MS)

legend({['Analytic Polyhedron $N_f=',num2str(meshCoarse.numFaces),'$'],...
    ['Simple Packing $N_i=',num2str(masconModel{1}.numElements),'$'],...
    ['FVM P1Q1 $N_i=',num2str(masconModel{2}.numElements),'$'],...
    ['FVM P2Q1 $N_i=',num2str(masconModel{3}.numElements),'$'],...
    ['FVM P2QC $N_i=',num2str(masconModel{4}.numElements),'$'],...
    ['FVM P2Q2 $N_i=',num2str(masconModel{5}.numElements),'$'],...
    ['FVM P2Q1Octree $N_i=',num2str(masconModel{6}.numElements),'$'],...
    ['FVM P2Q1NoSurf $N_i=',num2str(masconModel{7}.numElements),'$']...
    ['FVM P2Q1NoSurf layers $N_i=',num2str(masconModel{8}.numElements),'$']...
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



FS = 12;
LW = 1;
MS = 8;
figure(2)
hold on
plot(numElements(:,1),maxerror(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,2),maxerror(:,2),'ko-','MarkerFaceColor','r','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,3),maxerror(:,3),'ks-','MarkerFaceColor','c','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,4),maxerror(:,4),'ks-','MarkerFaceColor','b','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,5),maxerror(:,5),'ks-','MarkerFaceColor','g','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,6),maxerror(:,6),'ks-','MarkerFaceColor','m','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,7),maxerror(:,7),'ks-','MarkerFaceColor',[1,0.5,0.0],'LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,8),maxerror(:,8),'ks-','MarkerFaceColor',[1.0,0.5,1.0],'LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,9),maxerror(:,9),'ks-','MarkerFaceColor',[0.5,0.0,1.0],'LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,10),maxerror(:,10),'ks-','MarkerFaceColor',[1.0,1.0,1.0],'LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,11),maxerror(:,11),'ks-','MarkerFaceColor',[0.9,0.9,0.9],'LineWidth',LW,'MarkerSize',MS)

legend({['Analytic Polyhedron $N_f=',num2str(meshCoarse.numFaces),'$'],...
    ['Simple Packing $N_i=',num2str(masconModel{1}.numElements),'$'],...
    ['FVM P1Q1 $N_i=',num2str(masconModel{2}.numElements),'$'],...
    ['FVM P2Q1 $N_i=',num2str(masconModel{3}.numElements),'$'],...
    ['FVM P2QC $N_i=',num2str(masconModel{4}.numElements),'$'],...
    ['FVM P2Q2 $N_i=',num2str(masconModel{5}.numElements),'$'],...
    ['FVM P2Q1Octree $N_i=',num2str(masconModel{6}.numElements),'$'],...
    ['FVM P2Q1NoSurf $N_i=',num2str(masconModel{7}.numElements),'$']...
    ['FVM P2Q1NoSurf layers $N_i=',num2str(masconModel{8}.numElements),'$']...
    },'interpreter','latex','FontSize',FS)
ylabel('max error $\%$','interpreter','latex','FontSize',FS)
xlabel('altitude','interpreter','latex','FontSize',FS)
set(gcf,'Color',[1,1,1]);
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')
set(gca,'xscale','log'); set(gca,'yscale','log');
%axis([0.1,10000,0.1,1000])
legend boxoff
box on


% location of quadrature points
FS = 12;
LW = 1;
MS = 8;
figure(3)
hold on
plot(numElements(:,1),error(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,2),error(:,2),'ko-','MarkerFaceColor','r','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,4),error(:,4),'ks-','MarkerFaceColor','b','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,5),error(:,5),'ks-','MarkerFaceColor','g','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,6),error(:,6),'ks-','MarkerFaceColor','m','LineWidth',LW,'MarkerSize',MS)

legend({['Analytic Polyhedron'],...
    ['Simple Packing'],...
    ['P2 vertex'],...
    ['P2 centroid'],...
    ['P2 degree-2'],...
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



FS = 12;
LW = 1;
MS = 8;
figure(4)
hold on
plot(numElements(:,1),maxerror(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,2),maxerror(:,2),'ko-','MarkerFaceColor','r','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,4),maxerror(:,4),'ks-','MarkerFaceColor','b','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,5),maxerror(:,5),'ks-','MarkerFaceColor','g','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,6),maxerror(:,6),'ks-','MarkerFaceColor','m','LineWidth',LW,'MarkerSize',MS)

legend({['Analytic Polyhedron'],...
    ['Simple Packing'],...
    ['P2 vertex'],...
    ['P2 centroid'],...
    ['P2 degree-2'],...
    },'interpreter','latex','FontSize',FS)
ylabel('max error $\%$','interpreter','latex','FontSize',FS)
xlabel('altitude','interpreter','latex','FontSize',FS)
set(gcf,'Color',[1,1,1]);
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')
set(gca,'xscale','log'); set(gca,'yscale','log');
%axis([0.1,10000,0.1,1000])
legend boxoff
box on



% degree of mesh
FS = 12;
LW = 1;
MS = 8;
figure(5)
hold on
plot(numElements(:,1),error(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,2),error(:,2),'ko-','MarkerFaceColor','r','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,3),error(:,3),'ks-','MarkerFaceColor','c','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,4),error(:,4),'ks-','MarkerFaceColor','b','LineWidth',LW,'MarkerSize',MS)
legend({['Analytic Polyhedron'],...
    ['Simple Packing'],...
    ['P1 vertex'],...
    ['P2 vertex'],...
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



FS = 12;
LW = 1;
MS = 8;
figure(6)
hold on
plot(numElements(:,1),maxerror(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,2),maxerror(:,2),'ko-','MarkerFaceColor','r','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,3),maxerror(:,3),'ks-','MarkerFaceColor','c','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,4),maxerror(:,4),'ks-','MarkerFaceColor','b','LineWidth',LW,'MarkerSize',MS)
legend({['Analytic Polyhedron'],...
    ['Simple Packing'],...
    ['P1 vertex'],...
    ['P2 vertex'],...
    ['P2 degree-2'],...
    },'interpreter','latex','FontSize',FS)
ylabel('max error $\%$','interpreter','latex','FontSize',FS)
xlabel('altitude','interpreter','latex','FontSize',FS)
set(gcf,'Color',[1,1,1]);
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')
set(gca,'xscale','log'); set(gca,'yscale','log');
%axis([0.1,10000,0.1,1000])
legend boxoff
box on





% degree of mesh
FS = 12;
LW = 1;
MS = 8;
figure(7)
hold on
plot(numElements(:,1),error(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,2),error(:,2),'ko-','MarkerFaceColor','r','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,9),error(:,9),'ko-','MarkerFaceColor','m','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,10),error(:,10),'ks-','MarkerFaceColor','c','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,10),error(:,11),'ks-','MarkerFaceColor','b','LineWidth',LW,'MarkerSize',MS)
legend({['Analytic Polyhedron'],...
    ['Simple Packing'],...
    ['P2 1/2'],...
    ['P2 1/3'],...
    ['Piecewise 1/2'],...
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



FS = 12;
LW = 1;
MS = 8;
figure(8)
hold on
plot(numElements(:,1),maxerror(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,2),maxerror(:,2),'ko-','MarkerFaceColor','r','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,9),maxerror(:,9),'ko-','MarkerFaceColor','m','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,10),maxerror(:,10),'ks-','MarkerFaceColor','c','LineWidth',LW,'MarkerSize',MS)
plot(numElements(:,11),maxerror(:,11),'ks-','MarkerFaceColor','b','LineWidth',LW,'MarkerSize',MS)
legend({['Analytic Polyhedron'],...
    ['Simple Packing'],...
    ['P2 1/2'],...
    ['P2 1/3'],...
    ['Piecewise 1/2'],...
    },'interpreter','latex','FontSize',FS)
ylabel('max error $\%$','interpreter','latex','FontSize',FS)
xlabel('altitude','interpreter','latex','FontSize',FS)
set(gcf,'Color',[1,1,1]);
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')
set(gca,'xscale','log'); set(gca,'yscale','log');
%axis([0.1,10000,0.1,1000])
legend boxoff
box on




