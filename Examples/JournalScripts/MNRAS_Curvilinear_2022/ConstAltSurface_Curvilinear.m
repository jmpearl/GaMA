clear all; clc; close all;

% body parameters
G = 6.67*10^-11;                       % gravitational constant (N m2/kg2)
load('Bennu.mat');                      % load stored Eros properties                                  
Mu = bodyProperties.mass*G;            % set stand grav parameter
Omega = bodyProperties.rotationRate;   % 

% some model parameters
meshfile = 'Bennu_200000.obj';    % mesh file to load 
mesh = SurfaceMesh(meshfile);   % create the surface mesh object
mesh.setNumFaces(50000);
%mesh = mesh.center();

meshOriginal = mesh;

i=1;
quadRules = [1,2,3,4];
baseLine = AnalyticPolyhedralModel(mesh,Mu);
meshCoarse = SurfaceMesh(mesh);
meshCoarse.setNumFaces(23500);
while mesh.numVertices>1050
    meshCoarse.coarsen(mesh.numFaces-2000);
    meshCoarse.smooth(10,'uniform',false);
    meshCoarse.coordinates=meshOriginal.project(mesh.coordinates,mesh.vertexNormals());
    meshCoarse.coarsen(mesh.numFaces-1);
    %mesh = mesh.center();
    res(i) = mesh.resolution;
    numFaces(i) = mesh.numFaces;
    polyhedralModel{i} = AnalyticPolyhedralModel(mesh,Mu);
    for j = 1:length(quadRules)
        meshCoarse.setDegree(quadRules(j));   % make a P2 mesh
        meshCoarse.curve(meshOriginal);
        %meshCurve = meshCurve.center()
        quadratureModel{i}{j} = ApproximatePolyhedralModel(meshCoarse,Mu);
        numNodes2(i,j) = quadratureModel{i}{j}.numElements;
    end
    i=i+1;
end


altitudes = logspace(-1,4,20)*meshOriginal.resolution;
meshCoarse = meshOriginal.coarsen(4000);
iter = 1
for i = 1:length(altitudes)
    i
    constAltSurface = meshCoarse.offsetSurfaceMesh(altitudes(i), 500);
    %qm = ApproximatePolyhedralModel(meshCoarse,1);
    %weights = vecnorm(qm.vertexAreaVectors(),2,2);
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
            numNodes(iter,j+1) = quadratureModel{k}{j}.numElements;
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
Nres=43
i1 = i;
Nf(i1)
plot(h(i1:Nres:end)'/meshOriginal.resolution,error(i1:Nres:end,1),'kd','MarkerFaceColor','k','lineWidth',LW)
plot(h(i1:Nres:end)'/meshOriginal.resolution,error(i1:Nres:end,2),'ko','MarkerFaceColor','m','lineWidth',LW)
plot(h(i1:Nres:end)'/meshOriginal.resolution,error(i1:Nres:end,3),'k^','MarkerFaceColor','b','lineWidth',LW)
plot(h(i1:Nres:end)'/meshOriginal.resolution,error(i1:Nres:end,4),'ks','MarkerFaceColor','r','lineWidth',LW)
plot(h(i1:Nres:end)'/meshOriginal.resolution,error(i1:Nres:end,5),'k^','MarkerFaceColor','y','lineWidth',LW)
plot(h(i1:Nres:end)'/meshOriginal.resolution,error(i1:Nres:end,1),'kd','MarkerFaceColor','k','lineWidth',LW)
legend({'P1','P2','P3','P4','Analytic'},'interpreter','latex','FontSize',FS)
xlabel('$h/\delta$','interpreter','latex','FontSize',FS)
ylabel('$\epsilon\cdot N_f^d\cdot \frac{h}{R}$','interpreter','latex','FontSize',FS)
set(gcf,'Color',[1,1,1]);
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')
set(gca,'xscale','log'); set(gca,'yscale','log');
legend boxoff
box on

Nres=43
errorSquare1 = [];
errorSquare2 = [];
errorSquare3 = [];
errorSquare4 = [];
errorSquar5 = [];
for i=1:length(altitudes)
    i1 = Nres*(i-1)+1
    i2 = Nres*i
    errorSquare1 = [errorSquare1,error(i1:i2,1).*Nf(i1:i2)'];
    errorSquare2 = [errorSquare2,error(i1:i2,2).*Nf(i1:i2)'];
    errorSquare3 = [errorSquare3,error(i1:i2,3).*Nf(i1:i2)'];
    errorSquare4 = [errorSquare4,error(i1:i2,4).*Nf(i1:i2)'];
    errorSquare5 = [errorSquare5,error(i1:i2,5).*Nf(i1:i2)'];
end

FS=16; LW=1; MS=8
figure(3)
hold on
i1 = i;
factor67P = 1.00;
factorEros =   1.00;
factorPhobos = 1.00;
factorBennu=1.00
Nf(i1)
%plot(altitudes'/meshOriginal.resolution,mean(errorSquare2,1)*factorEros,'ko','MarkerFaceColor','m','lineWidth',LW,"markerSize",MS)
plot(altitudes'/meshOriginal.resolution,mean(errorSquare3,1)*factorEros,'ko','MarkerFaceColor','b','lineWidth',LW,"markerSize",MS)
%plot(altitudes'/meshOriginal.resolution,mean(errorSquare4,1)*factorEros,'ko','MarkerFaceColor','r','lineWidth',LW,"markerSize",MS)
%plot(altitudes'/meshOriginal.resolution,mean(errorSquare5,1)*factorEros,'ko','MarkerFaceColor','y','lineWidth',LW,"markerSize",MS)
plot(altitudes'/meshOriginal.resolution,mean(errorSquare1,1)*factorEros,'ko','MarkerFaceColor','k','lineWidth',LW,"markerSize",MS)
%plot(altitudes'/meshOriginal.resolution,mean(errorSquare67P2,1)*factor67P,'ks','MarkerFaceColor','m','lineWidth',LW,"markerSize",MS)
plot(altitudes'/meshOriginal.resolution,mean(errorSquare67P3,1)*factor67P,'ks','MarkerFaceColor','b','lineWidth',LW,"markerSize",MS)
%plot(altitudes'/meshOriginal.resolution,mean(errorSquare67P4,1)*factor67P,'ks','MarkerFaceColor','r','lineWidth',LW,"markerSize",MS)
%plot(altitudes'/meshOriginal.resolution,mean(errorSquare67P5,1)*factor67P,'ks','MarkerFaceColor','y','lineWidth',LW,"markerSize",MS)
plot(altitudes'/meshOriginal.resolution,mean(errorSquare67P1,1)*factor67P,'ks','MarkerFaceColor','k','lineWidth',LW,"markerSize",MS)
%plot(altitudes'/meshOriginal.resolution,mean(errorSquarePhobos2,1)*factorPhobos,'kd','MarkerFaceColor','m','lineWidth',LW,"markerSize",MS)
plot(altitudes'/meshOriginal.resolution,mean(errorSquarePhobos3,1)*factorPhobos,'kd','MarkerFaceColor','b','lineWidth',LW,"markerSize",MS)
%plot(altitudes'/meshOriginal.resolution,mean(errorSquarePhobos4,1)*factorPhobos,'kd','MarkerFaceColor','r','lineWidth',LW,"markerSize",MS)
%plot(altitudes'/meshOriginal.resolution,mean(errorSquarePhobos5,1)*factorPhobos,'kd','MarkerFaceColor','y','lineWidth',LW,"markerSize",MS)
plot(altitudes'/meshOriginal.resolution,mean(errorSquarePhobos1,1)*factorPhobos,'kd','MarkerFaceColor','k','lineWidth',LW,"markerSize",MS)
%plot(altitudes'/meshOriginal.resolution,mean(errorSquareBennu2,1)*factorBennu,'k^','MarkerFaceColor','m','lineWidth',LW,"markerSize",MS)
plot(altitudes'/meshOriginal.resolution,mean(errorSquareBennu3,1)*factorBennu,'k^','MarkerFaceColor','b','lineWidth',LW,"markerSize",MS)
%plot(altitudes'/meshOriginal.resolution,mean(errorSquareBennu4,1)*factorBennu,'k^','MarkerFaceColor','r','lineWidth',LW,"markerSize",MS)
%plot(altitudes'/meshOriginal.resolution,mean(errorSquareBennu5,1)*factorBennu,'k^','MarkerFaceColor','y','lineWidth',LW,"markerSize",MS)
plot(altitudes'/meshOriginal.resolution,mean(errorSquareBennu1,1)*factorBennu,'k^','MarkerFaceColor','k','lineWidth',LW,"markerSize",MS)

%plot(h'/meshOriginal.resolution,error(:,4).*Nf','rs','MarkerFaceColor','r','lineWidth',LW)
legend({'P1','P2','P3','P4','Analytic'},'interpreter','latex','FontSize',FS)
xlabel('$h/\delta_0$','interpreter','latex','FontSize',FS)
ylabel('$\epsilon\cdot N_f$','interpreter','latex','FontSize',FS)
set(gcf,'Color',[1,1,1]);
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')
set(gca,'xscale','log'); set(gca,'yscale','log');
xticks([0.1,1,10,100,1000,10000])
yticks(10.^[-2,0,2,4,5])
axis([0.1 10000,0.99e-2,1e5])
legend boxoff
box on

%1.2e4 10000, 4600, 3800, 3400
FS=16; LW=1; MS=8;
figure(4)
hold on
i1 = i;
Nf(i1)
%plot(altitudes'/meshOriginal.resolution,mean(errorSquare1./errorSquare2,1),'ko','MarkerFaceColor','m','lineWidth',LW,"markerSize",MS)
%plot(altitudes'/meshOriginal.resolution,mean(errorSquare1./errorSquare3,1),'ko','MarkerFaceColor','r','lineWidth',LW,"markerSize",MS)
plot(altitudes'/meshOriginal.resolution,mean(errorSquare1./errorSquare4,1),'ko','MarkerFaceColor','r','lineWidth',LW,"markerSize",MS)
%plot(altitudes'/meshOriginal.resolution,mean(errorSquare1./errorSquare5,1),'ko','MarkerFaceColor','y','lineWidth',LW,"markerSize",MS)
%plot(altitudes'/meshOriginal.resolution,mean(errorSquare67P1./errorSquare67P2,1),'ks','MarkerFaceColor','m','lineWidth',LW,"markerSize",MS)
%plot(altitudes'/meshOriginal.resolution,mean(errorSquare67P1./errorSquare67P3,1),'ks','MarkerFaceColor','c','lineWidth',LW,"markerSize",MS)
plot(altitudes'/meshOriginal.resolution,mean(errorSquare67P1./errorSquare67P4,1),'ks','MarkerFaceColor','c','lineWidth',LW,"markerSize",MS)
%plot(altitudes'/meshOriginal.resolution,mean(errorSquare67P1./errorSquare67P5,1),'ks','MarkerFaceColor','y','lineWidth',LW,"markerSize",MS)
%plot(altitudes'/meshOriginal.resolution,mean(errorSquarePhobos1./errorSquarePhobos2,1),'kd','MarkerFaceColor','m','lineWidth',LW,"markerSize",MS)
%plot(altitudes'/meshOriginal.resolution,mean(errorSquarePhobos1./errorSquarePhobos3,1),'kd','MarkerFaceColor','b','lineWidth',LW,"markerSize",MS)
plot(altitudes'/meshOriginal.resolution,mean(errorSquarePhobos1./errorSquarePhobos4,1),'kd','MarkerFaceColor','b','lineWidth',LW,"markerSize",MS)
%plot(altitudes'/meshOriginal.resolution,mean(errorSquarePhobos1./errorSquarePhobos5,1),'kd','MarkerFaceColor','y','lineWidth',LW,"markerSize",MS)
%plot(altitudes'/meshOriginal.resolution,mean(errorSquareBennu1./errorSquareBennu2,1),'k^','MarkerFaceColor','m','lineWidth',LW,"markerSize",MS)
%plot(altitudes'/meshOriginal.resolution,mean(errorSquareBennu1./errorSquareBennu3,1),'k^','MarkerFaceColor','m','lineWidth',LW,"markerSize",MS)
plot(altitudes'/meshOriginal.resolution,mean(errorSquareBennu1./errorSquareBennu4,1),'k^','MarkerFaceColor','m','lineWidth',LW,"markerSize",MS)
%plot(altitudes'/meshOriginal.resolution,mean(errorSquareBennu1./errorSquareBennu5,1),'k^','MarkerFaceColor','y','lineWidth',LW,"markerSize",MS)

%plot(h'/meshOriginal.resolution,error(:,4).*Nf','rs','MarkerFaceColor','r','lineWidth',LW)
legend({'Eros','67P','Phobos','Bennu'},'interpreter','latex','FontSize',FS)
xlabel('$h/\delta_0$','interpreter','latex','FontSize',FS)
ylabel('$\mathrm{error\ ratio\ }\frac{\mathrm{analytic}}{\mathrm{quadrature}}$','interpreter','latex','FontSize',FS)
set(gcf,'Color',[1,1,1]);
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')
set(gca,'xscale','log'); set(gca,'yscale','log');
%xticks([0.1,1,10,100,1000,10000])
%yticks(10.^[-4,-2,0,2,4,6])
axis([0.1,10000,0.1,1000])
legend boxoff
box on




FS=16; LW=1
figure(3)
hold on
i=10
altitudes(i)/meshOriginal.resolution
i1 = 22*(i-1)+1
i2 = 22*i
plot(Nf(i1:i2),error(i1:i2,1),'k+','MarkerFaceColor','m','lineWidth',LW)
plot(Nf(i1:i2),error(i1:i2,2),'ko','MarkerFaceColor','m','lineWidth',LW)
plot(Nf(i1:i2),error(i1:i2,3),'k^','MarkerFaceColor','b','lineWidth',LW)
plot(Nf(i1:i2),error(i1:i2,4),'ks','MarkerFaceColor','r','lineWidth',LW)
plot(Nf(i1:i2),error(i1:i2,5),'k^','MarkerFaceColor','y','lineWidth',LW)
legend({'Analytic','P1','P2','P3','P4'},'interpreter','latex','FontSize',FS)
xlabel('$N_f$','interpreter','latex','FontSize',FS)
ylabel('$\epsilon$','interpreter','latex','FontSize',FS)
set(gcf,'Color',[1,1,1]);
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')
set(gca,'xscale','log'); set(gca,'yscale','log');
legend boxoff
box on


FS=16; LW=1
figure(4)
hold on
plot(h(i:23:end)'/meshOriginal.resolution,error(i:23:end,1),'k+','MarkerFaceColor','m','lineWidth',LW)
plot(h(i:23:end)'/meshOriginal.resolution,error(i:23:end,2),'ko','MarkerFaceColor','m','lineWidth',LW)
plot(h(i:23:end)'/meshOriginal.resolution,error(i:23:end,3),'k^','MarkerFaceColor','b','lineWidth',LW)
plot(h(i:23:end)'/meshOriginal.resolution,error(i:23:end,4),'ks','MarkerFaceColor','r','lineWidth',LW)
plot(h(i:23:end)'/meshOriginal.resolution,error(i:23:end,5),'k^','MarkerFaceColor','y','lineWidth',LW)
legend({'Analytic','P1','P2','P3','P4'},'interpreter','latex','FontSize',FS)
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

