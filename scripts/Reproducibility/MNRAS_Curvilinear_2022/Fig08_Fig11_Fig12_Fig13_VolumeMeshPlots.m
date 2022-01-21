%==========================================================================
% Script used to create Figure 8, 11, 12, 13
% Pearl J.M., Hitt, D.L., "Cutting Corners: A Quadrature-based Gravity
% Model for Comets and Asteroids using Curvilinear Surface Definitions"
% MNRAS, 2022 (submitted).
% 
% ** Note GaMA is under development small variations from published 
% ** results may occur.
%
% writes a volume mesh vtk file with acceleration, laplacian and model
% error information for asteroid Eros.
%==========================================================================

clear all; close all; clc

G = 6.67e-11; 

% some model parameters
meshfile = 'Eros_46906.obj';           % mesh file to load 
load('Eros.mat');                      % load stored Eros properties                                  
Mu = bodyProperties.mass*G;            % set stand grav parameter
Omega = bodyProperties.rotationRate;   % 
numIters = 4;

% set up meshes
%--------------------------------------------------------------------------
mesh = SurfaceMesh(meshfile);              % create the surface mesh object
meshCoarse = SurfaceMesh(mesh);            % copy construct
meshCoarse.coarsen(8000);                  % create coarsened version
meshCoarse.smooth(1);                      % smooth things out

meshCoarseP2 = SurfaceMesh(mesh);
meshCoarseP2.coarsen(2000);
meshCoarseP2.smooth(1);
meshCoarseP2.setDegree(2);   % make a P2 mesh
meshCoarseP2.curve(mesh);    % curve mesh through projection

% set up gravity models 
quadratureModelFullRes = ApproximatePolyhedralModel(mesh,Mu);
quadratureModel     = ApproximatePolyhedralModel(meshCoarse,Mu);
quadratureModelP2   = ApproximatePolyhedralModel(meshCoarseP2,Mu);
analyticModelTruth  = AnalyticPolyhedralModel(mesh,Mu);
analyticModelCoarse = AnalyticPolyhedralModel(meshCoarse,Mu);

% create volume mesh around body. We're manually creating this for now so
% that we can track the altitude of each vertex
%--------------------------------------------------------------------------

% offseting outward
volumePoints = [mesh.coordinates];
sm2 = SurfaceMesh(mesh);
newVertexCount = floor(mesh.numVertices);
h = zeros(size(volumePoints,1),1);
hnew=0.0;
for i = 1:numIters
    sm2 = sm2.offsetSurfaceMesh(sm2.resolution/3.0,newVertexCount);
    hnew = hnew + sm2.resolution/3.0;
    h = [h;hnew*ones(newVertexCount,1)];
    volumePoints = [volumePoints;sm2.coordinates];
    newVertexCount = newVertexCount-1000;
end

% offsetting inward
sm3 = SurfaceMesh(mesh);
newVertexCount = floor(mesh.numVertices);
hnew=0.0;
for i = 1:numIters
    sm3 = sm3.offsetSurfaceMesh(-sm3.resolution/3.0,newVertexCount);
    hnew = hnew - sm3.resolution/3.0;
    h = [h;hnew*ones(newVertexCount,1)];
    volumePoints = [volumePoints;sm3.coordinates];
    newVertexCount = newVertexCount-1000;
end

volumeMesh = VolumeMesh();                          % construct w/out surface mesh
volumeMesh.initializeFromPointCloud(volumePoints);  % delaunay tet mesh from points

% calculation our field values at the nodes (this may take a while)
%--------------------------------------------------------------------------
pts = volumeMesh.coordinates;
truthAcceleration = analyticModelTruth.acceleration(pts);
polyhedralAcceleration    = analyticModelCoarse.acceleration(pts);
approximateAcceleration   = quadratureModel.acceleration(pts);
approximateAccelerationP2 = quadratureModelP2.acceleration(pts);

lap_fR  = quadratureModelFullRes.Laplacian(pts)/Mu*meshCoarse.volume;
lap_QM  = quadratureModel.Laplacian(pts)/Mu*meshCoarse.volume;
lap_QMP2 = quadratureModelP2.Laplacian(pts)/Mu*meshCoarseP2.volume;

epsTotP2 = 100*vecnorm(approximateAccelerationP2-truthAcceleration,2,2)./vecnorm(truthAcceleration,2,2);
epsNumP1 = 100*vecnorm(approximateAcceleration-polyhedralAcceleration,2,2)./vecnorm(polyhedralAcceleration,2,2);
epsPoly = 100*vecnorm(truthAcceleration-polyhedralAcceleration,2,2)./vecnorm(truthAcceleration,2,2);


% load the fields in and export for viewing w/ paraview
%--------------------------------------------------------------------------

% Fig 8, 11, 12 created in paraview
%----------------------------------
volumeMesh.clearFields();
volumeMesh.addNodeField(h,'altitude');
volumeMesh.addNodeField(lap_fR,'lap_fR2');
volumeMesh.addNodeField(lap_QM,'lap_QM');
volumeMesh.addNodeField(lap_QMP2,'lap_QMP2');
volumeMesh.addNodeField(epsTotP2,'epsP2');
volumeMesh.addNodeField(epsNumP1,'epsP1');
volumeMesh.addNodeField(epsPoly,'epsPoly');
volumeMesh.addNodeField(approximateAcceleration,'acc_P1');
volumeMesh.addNodeField(approximateAccelerationP2,'acc_P2');
volumeMesh.addNodeField(polyhedralAcceleration,'acc_poly');
volumeMesh.addNodeField(truthAcceleration,'acc');
volumeMesh.writeVTK('testExternalMeshOut.vtk')



% Fig 13
%--------
figure(1)
hold on 
plot(h(lap_fR<-6.28)/mesh.resolution,lap_fR(lap_fR<-6.28),'ko','MarkerFaceColor','c')
plot(h(lap_fR>-6.28)/mesh.resolution,lap_fR(lap_fR>-6.28),'ko','MarkerFaceColor','r')
box on
set(gca,'YTick',-8*pi:2*pi:4*pi)
set(gca,'YTickLabel',{'-8\pi','-6\pi','-4\pi','-2\pi','0','2\pi','4\pi'})
set(gca,'FontSize',11)
xlabel('$h/\delta$','interpreter','latex','FontSize',15)
ylabel('$\frac{\nabla^2 U}{G\rho} $','interpreter','latex','FontSize',20)
set(gcf,'Color',[1,1,1])