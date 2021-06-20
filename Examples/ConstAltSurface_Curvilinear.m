clear all; clc; close all;

% body parameters
G = 6.67*10^-11;                       % gravitational constant (N m2/kg2)
load('Eros.mat');                      % load stored Eros properties                                  
Mu = bodyProperties.mass*G;            % set stand grav parameter
Omega = bodyProperties.rotationRate;   % 

% some model parameters
meshfile = 'Eros_46906.obj';    % mesh file to load 
mesh = SurfaceMesh(meshfile);   % create the surface mesh object

meshOriginal = mesh;
i=1;
quadRules = [1,2,3,4];
baseLine = AnalyticPolyhedralModel(mesh,Mu);

while mesh.numVertices>2000
    mesh=mesh.coarsen(mesh.numFaces-2000);
    mesh=mesh.smooth(5);
    mesh = mesh.center();
    res(i) = mesh.resolution;
    numFaces(i) = mesh.numFaces;
    polyhedralModel{i} = AnalyticPolyhedralModel(mesh,Mu);
    for j = 1:length(quadRules)
        meshCurve = mesh.setDegree(quadRules(j));   % make a P2 mesh
        meshCurve = meshCurve.curve(meshOriginal);
        meshCurve = meshCurve.center();
        quadratureModel{i}{j} = ApproximatePolyhedralModel(meshCurve,Mu);
    end
    i=i+1
end


altitudes = logspace(-1,4,20)*meshOriginal.resolution;
meshCoarse = meshOriginal.coarsen(4000);
iter = 1
for i = 1:length(altitudes)
    i
    constAltSurface = meshCoarse.offsetSurfaceMesh(altitudes(i), 300);
    pts = constAltSurface.coordinates;
    accOriginal = baseLine.acceleration(pts);
    accOgMag = vecnorm(accOriginal,2,2);
    
    for k = 1:length(polyhedralModel)
        accTempPoly = polyhedralModel{k}.acceleration(pts);
        error(iter,1) = 100*mean(vecnorm((accTempPoly - accOriginal),2,2)./accOgMag);
        stderror(iter,1) = std(vecnorm((accTempPoly - accOriginal),2,2)./accOgMag);
        for j=1:length(quadRules)
            accTemp = quadratureModel{k}{j}.acceleration(pts);
            stderror(iter,j+1) = 100*std(vecnorm((accTemp - accOriginal),2,2)./accOgMag);
            stdnumError(iter,j+1) = 100*std(vecnorm((accTemp - accTempPoly),2,2)./vecnorm(accTempPoly,2,2));
            error(iter,j+1) = 100*mean(vecnorm((accTemp - accOriginal),2,2)./accOgMag);
            numError(iter,j+1) = 100*mean(vecnorm((accTemp - accTempPoly),2,2)./vecnorm(accTempPoly,2,2));
        end
        h(iter) = altitudes(i);
        delta(iter) = res(k);
        Nf(iter) = numFaces(k);
        iter = iter+1
    end
end

stop

FS=16; LW=1
figure(3)
hold on
    %plot(altitudes(j)/res(i),error{1}(i,j)/res(i)*altitudes(j),'kd','MarkerFaceColor','k','lineWidth',LW)
plot(h(22:23:end)',100*error(22:23:end,1),'k+','MarkerFaceColor','m','lineWidth',LW)
plot(h(22:23:end)',error(22:23:end,2),'ko','MarkerFaceColor','m','lineWidth',LW)
plot(h(22:23:end)',error(22:23:end,3),'k^','MarkerFaceColor','b','lineWidth',LW)
plot(h(22:23:end)',error(22:23:end,4),'ks','MarkerFaceColor','r','lineWidth',LW)
plot(h(22:23:end)',error(22:23:end,5),'k^','MarkerFaceColor','y','lineWidth',LW)
legend({'P1','P2','P3','P4'},'interpreter','latex','FontSize',FS)
xlabel('$h/\delta$','interpreter','latex','FontSize',FS)
ylabel('$\epsilon\cdot N_f^d\cdot \frac{h}{R}$','interpreter','latex','FontSize',FS)
set(gcf,'Color',[1,1,1]);
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')
set(gca,'xscale','log'); set(gca,'yscale','log');
legend boxoff
box on


FS=16; LW=1
figure(3)
hold on
    %plot(altitudes(j)/res(i),error{1}(i,j)/res(i)*altitudes(j),'kd','MarkerFaceColor','k','lineWidth',LW)
plot(h'./delta',numError(:,2).*Nf','ko','MarkerFaceColor','m','lineWidth',LW)
plot(h'./delta',numError(:,3).*Nf','k^','MarkerFaceColor','b','lineWidth',LW)
plot(h'./delta',numError(:,4).*Nf','ks','MarkerFaceColor','r','lineWidth',LW)
plot(h'./delta',numError(:,5).*Nf','k^','MarkerFaceColor','y','lineWidth',LW)
plot(h'./delta',numError(:,6).*Nf','ks','MarkerFaceColor','g','lineWidth',LW)
plot(h'./delta',numError(:,7).*Nf','k^','MarkerFaceColor','c','lineWidth',LW)
plot(h'./delta',numError(:,8).*Nf','ks','MarkerFaceColor',[0.5,0,1],'lineWidth',LW)
plot(h'./delta',numError(:,9).*Nf','ks','MarkerFaceColor',[1.0,0.5,0],'lineWidth',LW)

legend({'L1','L2','L3','G1','G2','B2','B3','O4'},'interpreter','latex','FontSize',FS)
xlabel('$h/\delta$','interpreter','latex','FontSize',FS)
ylabel('$\epsilon\cdot N_f^d\cdot \frac{h}{R}$','interpreter','latex','FontSize',FS)
set(gcf,'Color',[1,1,1]);
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')
set(gca,'xscale','log'); set(gca,'yscale','log');
legend boxoff
box on




R = mean(vecnorm(mesh.coordinates,2,2));
FS=16; LW=1
figure(4)
hold on

 plot(h./R,error(:,1),'kd','MarkerFaceColor','k','lineWidth',LW)

legend({'Analytic'},'interpreter','latex','FontSize',FS)
xlabel('h/R','interpreter','latex','FontSize',FS)
ylabel('$\epsilon \cdot N_f$','interpreter','latex','FontSize',FS)
set(gcf,'Color',[1,1,1]);
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')
set(gca,'xscale','log'); set(gca,'yscale','log');
legend boxoff
box on

