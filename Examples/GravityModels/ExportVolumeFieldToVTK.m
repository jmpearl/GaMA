clear all; close all; clc

G = 6.67e-11; 

% some model parameters
meshfile = 'Eros_46906.obj';           % mesh file to load 
load('Eros.mat');                      % load stored Eros properties                                  
Mu = bodyProperties.mass*G;            % set stand grav parameter
Omega = bodyProperties.rotationRate;   % 
numIters = 4;

% set up meshes
mesh = SurfaceMesh(meshfile);              % create the surface mesh object
meshCoarse = mesh.coarsen(8000);           % create coarsened version
meshCoarse = meshCoarse.smooth(5);         % smooth things out

meshCoarseP2 = mesh.coarsen(2000);
%meshCoarseP2 = meshCoarseP2.refineFeatures(mesh,2,9,2);
meshCoarseP2 = meshCoarseP2.smooth(5);
meshCoarseP2 = meshCoarseP2.setDegree(2);   % make a P2 mesh
meshCoarseP2 = meshCoarseP2.curve(mesh);    % curve mesh through projection

% set up gravity models 
quadratureModelFullRes = ApproximatePolyhedralModel(mesh,Mu);
quadratureModel     = ApproximatePolyhedralModel(meshCoarse,Mu);
quadratureModelP2   = ApproximatePolyhedralModel(meshCoarseP2,Mu);
analyticModelTruth  = AnalyticPolyhedralModel(mesh,Mu);
analyticModelCoarse = AnalyticPolyhedralModel(meshCoarse,Mu);

% create volume mesh around body
volumePoints = [mesh.coordinates];
sm2 = mesh;
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


sm3 = mesh;
newVertexCount = floor(mesh.numVertices);
hnew=0.0;
for i = 1:numIters
    sm3 = sm3.offsetSurfaceMesh(-sm3.resolution/3.0,newVertexCount);
    hnew = hnew - sm3.resolution/3.0;
    h = [h;hnew*ones(newVertexCount,1)];
    volumePoints = [volumePoints;sm3.coordinates];
    newVertexCount = newVertexCount-1000;
end

disp('here we are')
volumeMesh = VolumeMesh(mesh);
volumeMesh.coordinates = volumePoints;
volumeMesh.numVertices = size(volumePoints,1);
volumeMesh = volumeMesh.delaunayTriangulation();
volumeMesh.numVertices
volumeMesh.surfaceMesh = sm2;
disp('here we are')
volumeMesh = volumeMesh.clipExternalCells();


disp('here we are')
% calculation our field values at the nodes
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
volumeMesh = volumeMesh.clearFields();
volumeMesh = volumeMesh.addNodeField(h,'altitude');
volumeMesh = volumeMesh.addNodeField(lap_fR,'lap_fR2');
volumeMesh = volumeMesh.addNodeField(lap_QM,'lap_QM');
volumeMesh = volumeMesh.addNodeField(lap_QMP2,'lap_QMP2');
volumeMesh = volumeMesh.addNodeField(epsTotP2,'epsP2');
volumeMesh = volumeMesh.addNodeField(epsNumP1,'epsP1');
volumeMesh = volumeMesh.addNodeField(epsPoly,'epsPoly');
volumeMesh = volumeMesh.addNodeField(approximateAcceleration,'acc_P1');
volumeMesh = volumeMesh.addNodeField(approximateAccelerationP2,'acc_P2');
volumeMesh = volumeMesh.addNodeField(polyhedralAcceleration,'acc_poly');
volumeMesh = volumeMesh.addNodeField(truthAcceleration,'acc');
volumeMesh.writeVTK('testExternalMeshOut.vtk')


%sMeshCoarse = meshCoarse;
%submergedQuadratureModel = ApproximatePolyhedralModel(sMeshCoarse,Mu);
pts = mesh.coordinates;
surfAccTrue = analyticModelTruth.acceleration(pts);
surfAcc = analyticModelCoarse.acceleration(pts);
surfSQM = quadratureModelP2.acceleration(pts);
surfQM = quadratureModel.acceleration(pts);
epsApprox = 100*vecnorm(surfQM-surfAccTrue,2,2)./vecnorm(surfAccTrue,2,2);
epsApproxSub = 100*vecnorm(surfSQM-surfAccTrue,2,2)./vecnorm(surfAccTrue,2,2);
epsPoly = 100*vecnorm(surfAcc-surfAccTrue,2,2)./vecnorm(surfAccTrue,2,2);


stop
surfQM2 = quadratureModel.acceleration(pts);
epsApprox2 = 100*vecnorm(surfQM2-surfAccTrue,2,2)./vecnorm(surfAccTrue,2,2);
mesh = mesh.clearFields();
mesh = mesh.addNodeField(epsApprox,'epsapproxP1');
mesh = mesh.addNodeField(epsApprox2,'epsapproxNoMax');
mesh = mesh.addNodeField(epsApproxSub,'epsapproxP2');
mesh = mesh.addNodeField(epsPoly,'epspoly');

mesh.writeVTK('surfaceMesh.vtk')
meshCoarse.writeVTK('surfaceMeshCoarse.vtk')
% sMeshCoarse = meshCoarseP2;
% 
% Factor = 0.00;
% sMeshCoarse.coordinates = sMeshCoarse.coordinates - sMeshCoarse.nodeNormals()*Factor*sMeshCoarse.resolution;
% sMeshCoarse.volume = sMeshCoarse.calculateVolume();
% submergedQuadratureModel = ApproximatePolyhedralModel(sMeshCoarse,Mu);
% surfSQM = submergedQuadratureModel.acceleration(pts);
% epsApproxSub = 100*vecnorm(surfSQM-surfAccTrue,2,2)./vecnorm(surfAccTrue,2,2);
% 
% mesh = mesh.clearFields();
% mesh = mesh.addNodeField(epsApprox,'epsapprox');
% mesh = mesh.addNodeField(epsApproxSub,'epsapproxSUB2');
% mesh = mesh.addNodeField(epsPoly,'epspoly');
% mesh.writeVTK('surfaceMesh.vtk')
% meshCoarseP2.writeVTK('featureRefine.vtk')
% sMeshCoarse.writeVTK('subfeatureRefine.vtk')

figure(2)
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