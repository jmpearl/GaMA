%==========================================================================
% Script used for Figure 10 (a) (b) (c) of
% Pearl J.M., Hitt, D.L., "Cutting Corners: Curvilinear Surface-Based
% Gravity Models for Asteroids and Comets"
% 2022 (submitted).
%
% ** Note GaMA is under development small variations from published
% ** results may occur.
%
% Integrates trajectories around Eros for 1000 periods using different
% gravity models. Position and energy error is assessed relative to the
% polyhedral model on a fine, 46906-facet, mesh.
%==========================================================================
clear all; close all; clc;


G = 6.67e-11;

minRo = 35000;       % min orbital radius
maxRo = 50000;       % max orbital radius
numTrajectories=25;  % number of trajectories we're integrating
nPeriods = 100;      % integration interval in periods (circular)
nSamplesPeriod = 20; % number of data points / period in interpolation

% load body data
%--------------------------------------------------------------------------
load('Eros.mat');                     % load in on Eros
Mu = bodyProperties.mass*G;           % standard gravitational parameter
Omega = bodyProperties.rotationRate;  % rotational velocity
Omega=0.0;                            % set rotation to zero for this test

% set up the surface meshes we'll use
%--------------------------------------------------------------------------
meshFile = 'Eros_46906.obj';      % define surface mesh file
mesh = SurfaceMesh(meshFile);     % create tri surface mesh
mesh.center();                    % move centroid to (0,0,0)

meshP1 = SurfaceMesh(mesh);       % copy construct
meshP1.setNumFaces(8000);         % coursen to 8000 faces
meshP1.smooth(1);                 % smooth 1 iteration
meshP1.projectOnTo(mesh);         % project back onto mesh

meshP2 = SurfaceMesh(meshP1);     % copy construct
meshP2.setDegree(2);              % 2nd order faces
meshP2.curve(mesh);               % curve through projection

meshP3 = SurfaceMesh(meshP1);
meshP3.setDegree(3);
meshP3.curve(mesh);

meshP4 = SurfaceMesh(meshP1);
meshP4.setDegree(4);
meshP4.curve(mesh);

% set up the gravity models
%--------------------------------------------------------------------------
gravityModel{1} = AnalyticPolyhedralModel(mesh,Mu);      % Werner 1994 - truth
gravityModel{2} = AnalyticPolyhedralModel(meshP1,Mu);    % Werner 1994 - 8000-faces
gravityModel{3} = ApproximatePolyhedralModel(meshP1,Mu); % P1 - 8000-faces
gravityModel{4} = ApproximatePolyhedralModel(meshP2,Mu); % P2 - 8000-faces
gravityModel{5} = ApproximatePolyhedralModel(meshP3,Mu); % P3 - 8000-faces
gravityModel{6} = ApproximatePolyhedralModel(meshP4,Mu); % P4 - 8000-faces


% Integrate in the Body Fixed Frame
%--------------------------------------------------------------------------

for kk=1:numTrajectories

    % initial conditions
    Ro = minRo + rand * (maxRo-minRo);           % initial radius
    Vcirc = sqrt(Mu/Ro);                         % for calc initial energy
    xo = [rand-0.5,rand-0.5,rand-0.5];           % iinitial position vector dir
    xo = xo/norm(xo)*Ro;                         % initial position vector
    x1 = [rand-0.5,rand-0.5,rand-0.5];           % dummy vector
    vo = cross(x1,xo);                           % perpendicular vector
    vo = vo/norm(vo)*Vcirc;                      % initial velocity
    Xo = [xo(1),vo(1),xo(2),vo(2),xo(3),vo(3)];  % initial state


    % time period
    Tperiod = 2*pi*sqrt((Ro)^3/Mu);
    T = [0,nPeriods*Tperiod];
    t_interp = linspace(T(1),T(2),nPeriods*nSamplesPeriod);

    % integrator
    tol = 1e-8;
    odeOptions = odeset('RelTol',tol,'AbsTol',[tol tol tol tol tol tol]);
    integrator = Integrator(gravityModel{1},@ode45);
    integrator.setOdeOptions(odeOptions);
    integrator.setFrame('inertial');
    integrator.omega = Omega;

    for i =1:length(gravityModel)
        
        integrator.gravityModel = gravityModel{i};         % reset 2 next model

        tic
        [tout{kk}{i},xout{kk}{i}] = integrator.integrate(T,Xo);    % integrate
        toc

        % calc our output metrics
        p{kk}{i} = [xout{kk}{i}(:,1),xout{kk}{i}(:,3),xout{kk}{i}(:,5)];
        v{kk}{i} = [xout{kk}{i}(:,2),xout{kk}{i}(:,4),xout{kk}{i}(:,6)];

        xoutBFF = inertial2BFF(tout{kk}{i},p{kk}{i},Omega);
        U = gravityModel{i}.potential(xoutBFF);
        energy{kk}{i} = 0.5*sum(v{kk}{i}.*v{kk}{i},2)+U;

        pInterp{kk}{i} = interp1(tout{kk}{i},p{kk}{i},t_interp,'spline');
        eInterp{kk}{i} = interp1(tout{kk}{i},energy{kk}{i},t_interp,'spline');

    end

end

% Process data, remove oscillations be taking integrated max
%--------------------------------------------------------------------------
t_interplog = [ones(length(t_interp),1),log10(t_interp')];
for kk=1:numTrajectories

    for i = 2:length(eInterp{kk})
        errorE{kk}{i-1} = 100*abs((eInterp{kk}{i}-eInterp{kk}{1})./eInterp{kk}{1})';
        errorP{kk}{i-1} = 100*vecnorm(pInterp{kk}{i}-pInterp{kk}{1},2,2)/Ro;


        for j = 2:length(eInterp{kk}{i})

            errorE{kk}{i-1}(j) = max(errorE{kk}{i-1}(j),errorE{kk}{i-1}(j-1));
            errorP{kk}{i-1}(j) = max(errorP{kk}{i-1}(j),errorP{kk}{i-1}(j-1));

        end

        paddedE = log10(errorE{kk}{i-1});
        paddedP = log10(errorP{kk}{i-1});

    end
end

stop
meshP1.coordinates = meshP1.coordinates/1000;
meshP1.resolution  = meshP1.calculateResolution();
t_interp = t_interp/Tperiod;


% Plot our results
%--------------------------------------------------------------------------
% Fig 10a
%----------
figure(1)
LW=1; LWA=3; MS=8; YC1=0.9; YC2=0.85; FS = 13; n=1000;
hold on
nPeriods2Plot = 0.65;
for kk=1:numTrajectories
    nDataPoints = size(p{kk}{1},1);
    start = floor(nDataPoints-nDataPoints/nPeriods*nPeriods2Plot);
for i = start:nDataPoints
    colorVec = (nDataPoints-i)/(nDataPoints-start)*[1,1,1];
    plot3(p{kk}{1}(i-1:i,1)/1000,p{kk}{1}(i-1:i,2)/1000,p{kk}{1}(i-1:i,3)/1000,'k-','Color',colorVec,'LineWidth',1)
end
end

meshP1.flatten();
meshP1.plot('plotType','userDefined',...
    'edgeColor','None',...
    'faceColor','w',...
    'faceAlpha',0.9,...
    'lighting',[0,1,-0.1],... % specify vec of incoming light ray
    'shadowing','on')         % ray tracing on (default)
daspect([1,1,1])
xlabel('x (km)','interpreter','latex','FontSize',FS);
ylabel('y (km)','interpreter','latex','FontSize',FS);
zlabel('z (km)','interpreter','latex','FontSize',FS);
xticks([-50,0,50])
yticks([-50,0,50])
set(gcf,'Color',[1,1,1]);
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')
box on


% Figure 10 b
%---------------
figure(2)
LW=2; FS=14;
hold on
for kk=1:numTrajectories
    plot(t_interp,errorE{kk}{1},'k-','LineWidth',LWA)
    plot(t_interp,errorE{kk}{2},'m-','LineWidth',LW)
    plot(t_interp,errorE{kk}{3},'b-','LineWidth',LW)
    plot(t_interp,errorE{kk}{4},'r-','LineWidth',LW)
    plot(t_interp,errorE{kk}{5},'y-','Color',[YC1,YC2,0],'LineWidth',LW)
end
box on
axis([1,100,0.8e-5,.1])
xlabel('time (periods)','interpreter','latex','FontSize',FS);
ylabel('energy error ($\%$)','interpreter','latex','FontSize',FS);
set(gcf,'Color',[1,1,1]);
set(gca,'yscale','log');set(gca,'xscale','log');
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')

% Figure 10 c
%----------------
stepSize=60
alpha = 0.6
figure(3)
hold on
for i =1:5
averageErrorP{i} = errorP{1}{i};
for kk = 2:numTrajectories
    averageErrorP{i} = [averageErrorP{i},errorP{kk}{i}];
end
meanErrorP{i} = mean(averageErrorP{i},2);

end
plot(t_interp(1:stepSize:end)',meanErrorP{1}(1:stepSize:end),'kd-','MarkerFaceColor','k','MarkerSize',MS)
plot(t_interp(1:stepSize:end)',meanErrorP{2}(1:stepSize:end),'ko-','MarkerFaceColor','m','MarkerSize',MS)
plot(t_interp(1:stepSize:end)',meanErrorP{3}(1:stepSize:end),'k^-','MarkerFaceColor','b','MarkerSize',MS)
plot(t_interp(1:stepSize:end)',meanErrorP{4}(1:stepSize:end),'ks-','MarkerFaceColor','r','MarkerSize',MS)
plot(t_interp(1:stepSize:end)',meanErrorP{5}(1:stepSize:end),'kp-','MarkerFaceColor',[1,1,0],'MarkerSize',MS)
for kk=1:numTrajectories
    plot(t_interp,errorP{kk}{1},'k-','LineWidth',LWA,'Color',alpha*[1,1,1]+(1-alpha)*[0,0,0])
    plot(t_interp,errorP{kk}{2},'m-','LineWidth',LW,'Color',alpha*[1,1,1]+(1-alpha)*[1,0,1])
    plot(t_interp,errorP{kk}{3},'b-','LineWidth',LW,'Color',alpha*[1,1,1]+(1-alpha)*[0,0,1])
    plot(t_interp,errorP{kk}{4},'r-','LineWidth',LW,'Color',alpha*[1,1,1]+(1-alpha)*[1,0,0])
    plot(t_interp,errorP{kk}{5},'r-','LineWidth',LW,'Color',alpha*[1,1,1]+(1-alpha)*[YC1,YC2,0])
end
plot(t_interp(1:stepSize:end)',meanErrorP{1}(1:stepSize:end),'kd-','MarkerFaceColor','k','MarkerSize',MS)
plot(t_interp(1:stepSize:end)',meanErrorP{2}(1:stepSize:end),'ko-','MarkerFaceColor','m','MarkerSize',MS)
plot(t_interp(1:stepSize:end)',meanErrorP{3}(1:stepSize:end),'k^-','MarkerFaceColor','b','MarkerSize',MS)
plot(t_interp(1:stepSize:end)',meanErrorP{4}(1:stepSize:end),'ks-','MarkerFaceColor','r','MarkerSize',MS)
plot(t_interp(1:stepSize:end)',meanErrorP{5}(1:stepSize:end),'kp-','MarkerFaceColor',[1,1,0],'MarkerSize',MS)
xlabel('time (periods)','interpreter','latex','FontSize',FS);
ylabel('position error ($\%$ of orbital radius)','interpreter','latex','FontSize',FS);
set(gcf,'Color',[1,1,1]);
set(gca,'yscale','log');set(gca,'xscale','log');
axis([1,100,0.7*10^-4,100])
legend({'Analytic-7996','P1-7996','P2-7996','P3-7996','P4-7996'},'interpreter','latex','location','NorthWest','FontSize',FS)
legend boxoff
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')
box on




