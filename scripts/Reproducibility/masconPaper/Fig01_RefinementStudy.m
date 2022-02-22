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




