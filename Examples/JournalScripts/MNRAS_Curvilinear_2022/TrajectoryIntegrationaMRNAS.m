%==========================================================================
% Integrates trajectory in the Body-fixed frame and Inertial frame and then
% plots the result in both frames.
%==========================================================================
clear all; close all; clc;


G = 6.67e-11;

% load body data
%--------------------------------------------------------------------------
load('Eros.mat');                                 
Mu = bodyProperties.mass*G;
Omega = bodyProperties.rotationRate;
Omega=0.0;
% set up the surface meshes we'll use
%--------------------------------------------------------------------------
meshFile = 'Eros_46906.obj';           % define surface mesh file
mesh = SurfaceMesh(meshFile);          % create tri surface mesh
mesh = mesh.center();

meshCoarse = mesh.setNumFaces(8000);   % coursen to 8000 faces
meshCoarse = meshCoarse.smooth(5);
meshCoarse = meshCoarse.center();

meshP2 = mesh.setNumFaces(2000);       % coursen to 2000 faces
meshP2 = meshP2.smooth(5);             % smooth mesh
meshP2 = meshP2.setDegree(2);          % quadratic faces
meshP2 = meshP2.curve(mesh);           % curve through projection
meshP2 = meshP2.center();

meshP3 = mesh.setNumFaces(889);
meshP3 = meshP3.smooth(5);
meshP3 = meshP3.setDegree(3);
meshP3 = meshP3.curve(mesh);
meshP3 = meshP3.center();

meshP4 = mesh.setNumFaces(500);
meshP4 = meshP4.smooth(5);
meshP4 = meshP4.setDegree(4);
meshP4 = meshP4.curve(mesh);
meshP4 = meshP4.center();

% set up the gravity models
%--------------------------------------------------------------------------
gravityModel{1} = AnalyticPolyhedralModel(mesh,Mu);          % Werner 1994 - truth
gravityModel{2} = AnalyticPolyhedralModel(meshCoarse,Mu);    % Werner 1994 - 8000-faces
gravityModel{3} = ApproximatePolyhedralModel(meshCoarse,Mu); % P1 - 8000-faces

meshCoarse=meshCoarse.setDegree(2);
meshCoarse=meshCoarse.curve(mesh);
meshCoarse=meshCoarse.center();
gravityModel{4} = ApproximatePolyhedralModel(meshCoarse,Mu); % P2 - 8000-faces

meshCoarse=meshCoarse.flatten();
meshCoarse=meshCoarse.setDegree(3);
meshCoarse=meshCoarse.curve(mesh);
meshCoarse=meshCoarse.center();
gravityModel{5} = ApproximatePolyhedralModel(meshCoarse,Mu); % P3 - 8000-faces

meshCoarse=meshCoarse.flatten();
meshCoarse=meshCoarse.setDegree(4);
meshCoarse=meshCoarse.curve(mesh);
meshCoarse=meshCoarse.center();
gravityModel{6} = ApproximatePolyhedralModel(meshCoarse,Mu); % P4 - 8000-faces

gravityModel{7} = ApproximatePolyhedralModel(meshP2,Mu);     % P2 - 2000-faces
gravityModel{8} = ApproximatePolyhedralModel(meshP3,Mu);     % P3 -  890-faces
gravityModel{9} = ApproximatePolyhedralModel(meshP4,Mu);     % P4 -  500-faces

% Integrate in the Body Fixed Frame
%--------------------------------------------------------------------------
disp('Body-Fixed Frame Integration')


% initial conditions
Ro = 45000;                 % initial radius
Vcirc = sqrt(Mu/Ro);        % initially circular
Xo = [Ro,0,0,0,0,0];        % I.C.s [rx,vx,ry,vy,rz,vz]
tol = 1e-6;                 % integration tolerances

% time period
nPeriods = 1000;
nSamplesPeriod = 20;
Tperiod = 2*pi*sqrt((Ro)^3/Mu);
T = [0,nPeriods*Tperiod];   
t_interp = linspace(T(1),T(2),nPeriods*nSamplesPeriod);

% integrator
odeOptions = odeset('RelTol',tol,'AbsTol',[tol tol tol tol tol tol]);
integrator = Integrator(gravityModel{1},@ode45);
integrator = integrator.setOdeOptions(odeOptions);
integrator = integrator.setFrame('inertial');
integrator.omega = Omega;

for i =1:length(gravityModel)
    
    U = gravityModel{i}.potential([Xo(1),Xo(3),Xo(5)]);
    Xo(4) = -sqrt(-U);
    integrator.gravityModel = gravityModel{i};         % reset 2 next model
    
    tic
    [tout{i},xout{i}] = integrator.integrate(T,Xo);    % integrate
    toc

    % calc our output metrics
    p{i} = [xout{i}(:,1),xout{i}(:,3),xout{i}(:,5)];
    v{i} = [xout{i}(:,2),xout{i}(:,4),xout{i}(:,6)];

    xoutBFF = inertial2BFF(tout{i},p{i},Omega);     
    U = gravityModel{i}.potential(xoutBFF);
    energy{i} = 0.5*sum(v{i}.*v{i},2)+U;
    
    pInterp{i} = interp1(tout{i},p{i},t_interp,'spline');
    eInterp{i} = interp1(tout{i},energy{i},t_interp,'spline');
 
end

meshCoarse.coordinates=meshCoarse.coordinates/1000;
meshCoarse.resolution=meshCoarse.calculateResolution();

t_interp = t_interp/Tperiod;

figure(2)
FS = 12
MS=8
hold on
n=1000
plot3(p{1}(end,1)/1000,p{1}(end,2)/1000,p{1}(end,3)/1000,'kd','MarkerFaceColor','k','MarkerSize',MS)
plot3(p{2}(end,1)/1000,p{2}(end,2)/1000,p{2}(end,3)/1000,'k>','MarkerFaceColor','k','MarkerSize',MS)
plot3(p{3}(end,1)/1000,p{3}(end,2)/1000,p{3}(end,3)/1000,'ks','MarkerFaceColor','m','MarkerSize',MS)
plot3(p{4}(end,1)/1000,p{4}(end,2)/1000,p{4}(end,3)/1000,'ko','MarkerFaceColor','b','MarkerSize',MS)
plot3(p{5}(end,1)/1000,p{5}(end,2)/1000,p{5}(end,3)/1000,'ko','MarkerFaceColor','r','MarkerSize',MS)
plot3(p{6}(end,1)/1000,p{6}(end,2)/1000,p{6}(end,3)/1000,'ko','MarkerFaceColor','y','MarkerSize',MS)

nPeriods2Plot = 0.65;
nDataPoints = size(p{1},1);
start = floor(nDataPoints-nDataPoints/nPeriods*nPeriods2Plot);
for i = start:nDataPoints
    colorVec = (nDataPoints-i)/(nDataPoints-start)*[1,1,1];
    plot3(p{1}(i-1:i,1)/1000,p{1}(i-1:i,2)/1000,p{1}(i-1:i,3)/1000,'k-','Color',colorVec,'LineWidth',2)
end

plot3(p{1}(end,1)/1000,p{1}(end,2)/1000,p{1}(end,3)/1000,'kd','MarkerFaceColor','k','MarkerSize',MS)
plot3(p{2}(end,1)/1000,p{2}(end,2)/1000,p{2}(end,3)/1000,'k>','MarkerFaceColor','k','MarkerSize',MS)
plot3(p{3}(end,1)/1000,p{3}(end,2)/1000,p{3}(end,3)/1000,'ks','MarkerFaceColor','m','MarkerSize',MS)
plot3(p{4}(end,1)/1000,p{4}(end,2)/1000,p{4}(end,3)/1000,'ko','MarkerFaceColor','b','MarkerSize',MS)
plot3(p{5}(end,1)/1000,p{5}(end,2)/1000,p{5}(end,3)/1000,'ko','MarkerFaceColor','r','MarkerSize',MS)
plot3(p{6}(end,1)/1000,p{6}(end,2)/1000,p{6}(end,3)/1000,'ko','MarkerFaceColor','y','MarkerSize',MS)


meshCoarse = meshCoarse.flatten();
meshCoarse.plot('plotType','userDefined',... 
                'edgeColor','None',...             
                'faceColor','w',...         
                'faceAlpha',0.9,...
                'lighting',[0,1,-0.1],... % specify vec of incoming light ray
                'shadowing','on')      % ray tracing on (default)
%daspect([1,1,1])
xlabel('x (km)','interpreter','latex'); 
ylabel('y (km)','interpreter','latex'); 
zlabel('z (km)','interpreter','latex');
%xticks([-50,0,50])
%yticks([-50,0,50])
set(gcf,'Color',[1,1,1]);
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')
box on
legend({'Analytic-46906','Analytic-8000','P1-8000','P2-8000','P3-8000','P4-8000'},'interpreter','latex')
legend boxoff


% filter out our periodic oscilations
t_interplog = [ones(length(t_interp),1),log10(t_interp')];
for i = 2:length(eInterp)
    errorE{i-1} = 100*abs((eInterp{i}-eInterp{1})./eInterp{1})';
    errorP{i-1} = 100*vecnorm(pInterp{i}-pInterp{1},2,2)/Ro;
    
    
    for j = 2:length(eInterp{i})
        
        errorE{i-1}(j) = max(errorE{i-1}(j),errorE{i-1}(j-1));
        errorP{i-1}(j) = max(errorP{i-1}(j),errorP{i-1}(j-1));
        
    end
    
    paddedE = log10(errorE{i-1});
    paddedP = log10(errorP{i-1});
    %betaE{i-1} = t_interplog(1000:5000,:)\paddedE(1000:5000,:);
    %betaP{i-1} = t_interplog(1000:5000,:)\paddedP(1000:5000,:);
end
%figure(1)
%plot(t_interplog(10000:end),log10(errorE{i-1}(10000:end,:)));

figure(3)
LW=2
LWA=3
MS=8
YC1=0.9
YC2=0.85
hold on
plot(t_interp,errorE{1},'k-','LineWidth',LWA)
plot(t_interp,errorE{2},'m-','LineWidth',LW)
plot(t_interp,errorE{3},'b--','LineWidth',LW)
plot(t_interp,errorE{4},'r:','LineWidth',LW)
plot(t_interp,errorE{5},'y-.','Color',[YC1,YC2,0],'LineWidth',LW)
%plot(t_interp,100*abs((eInterp{7}-eInterp{1})./eInterp{1}),'b--','LineWidth',LW)
%plot(t_interp,100*abs((eInterp{8}-eInterp{1})./eInterp{1}),'r--','LineWidth',LW)
%plot(t_interp,100*abs((eInterp{9}-eInterp{1})./eInterp{1}),'r--','Color',[1,0.5,0],'LineWidth',LW)

box on

axis([1,1000,10^-4,1])
xlabel('time (periods)','interpreter','latex'); 
ylabel('energy error ($\%$)','interpreter','latex'); 
set(gcf,'Color',[1,1,1]);
set(gca,'yscale','log');set(gca,'xscale','log');
legend({'Analytic-8000','P1-8000','P2-8000','P3-8000','P4-8000'},'interpreter','latex','location','NorthWest')
legend boxoff 
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')



figure(6)
LW=2
FS=12
hold on
plot(t_interp,errorP{1},'k-','LineWidth',LWA)
plot(t_interp,errorP{2},'m-','LineWidth',LW)
plot(t_interp,errorP{3},'b--','LineWidth',LW)
plot(t_interp,errorP{4},'r:','LineWidth',LW)
plot(t_interp,errorP{5},'r-.','Color',[YC1,YC2,0],'LineWidth',LW)
%plot(t_interp,100*vecnorm(pInterp{7}-pInterp{1},2,2)/Ro,'b--','LineWidth',LW)
%plot(t_interp,100*vecnorm(pInterp{8}-pInterp{1},2,2)/Ro,'r--','LineWidth',LW)
%plot(t_interp,100*vecnorm(pInterp{9}-pInterp{1},2,2)/Ro,'r--','Color',[1,0.5,0],'LineWidth',LW)

xlabel('time (periods)','interpreter','latex'); 
ylabel('position error ($\%$ of orbital radius)','interpreter','latex'); 
set(gcf,'Color',[1,1,1]);
set(gca,'yscale','log');set(gca,'xscale','log');
axis([1,1000,10^-3,200])
legend({'Analytic-8000','P1-8000','P2-8000','P3-8000','P4-8000'},'interpreter','latex','location','NorthWest')
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')
box on
legend boxoff

figure(7)
LW=2
FS=12
hold on
plot(t_interp,errorP{1},'k-','LineWidth',LWA)
plot(t_interp,errorP{2},'m-','LineWidth',LW)%
%plot(t_interp,100*vecnorm(pInterp{4}-pInterp{1},2,2)/Ro,'b-','LineWidth',LW)
%plot(t_interp,100*vecnorm(pInterp{5}-pInterp{1},2,2)/Ro,'r-','LineWidth',LW)
%plot(t_interp,100*vecnorm(pInterp{6}-pInterp{1},2,2)/Ro,'r-','Color',[1,0.5,0],'LineWidth',LW)
plot(t_interp,errorP{6},'b--','LineWidth',LW)
plot(t_interp,errorP{7},'r:','LineWidth',LW)
plot(t_interp,errorP{8},'y-.','Color',[YC1,YC2,0],'LineWidth',LW)

xlabel('time (periods)','interpreter','latex'); 
ylabel('position error ($\%$ of orbital radius)','interpreter','latex'); 
set(gcf,'Color',[1,1,1]);
set(gca,'yscale','log');set(gca,'xscale','log');
axis([1,1000,10^-3,200])
legend({'Analytic-8000','P1-8000','P2-2000','P3-890','P4-500'},'interpreter','latex','location','NorthWest')
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')
box on
legend boxoff

figure(8)
LW=2
FS=12
hold on
plot(t_interp,errorE{1},'k-','LineWidth',LWA)
plot(t_interp,errorE{2},'m-','LineWidth',LW)%
plot(t_interp,errorE{6},'b--','LineWidth',LW)
plot(t_interp,errorE{7},'r:','LineWidth',LW)
plot(t_interp,errorE{8},'r-.','Color',[YC1,YC2,0],'LineWidth',LW)

xlabel('time (periods)','interpreter','latex'); 
ylabel('energy error ($\%$)','interpreter','latex'); 
set(gcf,'Color',[1,1,1]);
set(gca,'yscale','log');set(gca,'xscale','log');
axis([1,1000,10^-4,1])
legend({'Analytic-8000','P1-8000','P2-2000','P3-890','P4-500'},'interpreter','latex','location','NorthWest')
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')
box on
legend boxoff


