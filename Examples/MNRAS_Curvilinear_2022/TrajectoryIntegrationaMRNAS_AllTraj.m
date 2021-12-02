%==========================================================================
% Integrates trajectory in the Body-fixed frame and Inertial frame and then
% plots the result in both frames.
%==========================================================================
clear all; close all; clc;


G = 6.67e-11;

%meshFile = ['Eros_46906.mat','
%bodies = ['Eros.mat','Bennu.mat','Phobos.mat','67P.mat'];
aOverRmaxRange = [2.0,10.0];
eccentricityRange = [0.0,0.0];
numFacesRange = [500,12000];
numFacesTruth = 25000;
numOrbits = 10;
nPeriods = 100;

% load body data
%--------------------------------------------------------------------------
load('Eros.mat');                                 
Mu = bodyProperties.mass*G;
Omega = bodyProperties.rotationRate;
Omega=0.0;


for j = 1:numOrbits
% our random orbit variables
%--------------------------------------------------------------------------
numFacesi = floor((numFacesRange(2)-numFacesRange(1))*rand + numFacesRange(1));
ecci = (eccentricityRange(2)-eccentricityRange(1))*rand + eccentricityRange(1);
aRi = (aOverRmaxRange(2)-aOverRmaxRange(1))*rand + aOverRmaxRange(1);

% set up the surface meshes we'll use
%--------------------------------------------------------------------------
meshFile = 'Eros_46906.obj';            % define surface mesh file
mesh = SurfaceMesh(meshFile);           % create tri surface mesh
mesh = mesh.setNumFaces(numFacesTruth);
mesh = mesh.center();

meshCoarse = mesh.setNumFaces(numFacesi);   
meshCoarse = meshCoarse.smooth(5);
meshCoarse = meshCoarse.center();

meshP2 = meshCoarse.setDegree(2);        
meshP2 = meshP2.curve(mesh);           
meshP2 = meshP2.center();

meshP3 = meshCoarse.setDegree(3);
meshP3 = meshP3.curve(mesh);
meshP3 = meshP3.center();

meshP4 = meshCoarse.setDegree(4);
meshP4 = meshP4.curve(mesh);
meshP4 = meshP4.center();

% set up the gravity models
%--------------------------------------------------------------------------
gravityModel{1} = AnalyticPolyhedralModel(mesh,Mu);          % Werner 1994 - truth
gravityModel{2} = AnalyticPolyhedralModel(meshCoarse,Mu);    % Werner 1994 - test
gravityModel{3} = ApproximatePolyhedralModel(meshCoarse,Mu); % P1 - test
gravityModel{4} = ApproximatePolyhedralModel(meshP2,Mu);     % P2 - test
gravityModel{5} = ApproximatePolyhedralModel(meshP3,Mu);     % P3 - test
gravityModel{6} = ApproximatePolyhedralModel(meshP4,Mu);     % P4 - test

% Integrate in the Body Fixed Frame
%--------------------------------------------------------------------------
disp('Body-Fixed Frame Integration')

% initial conditions
Rmax = max(vecnorm(mesh.coordinates,2,2));
ai = aRi*Rmax;
xo = [rand,rand,rand]'-0.5;
xo = xo/norm(xo)*ai/(1+ecci);

Uo  = Mu*(2/norm(xo)-1/ai);

vo = [rand,rand,rand]'-0.5;
vo = cross(xo,vo);
vo = vo/norm(vo)*sqrt(Uo);

Xo = [xo(1);vo(1);xo(2);vo(2);xo(3);vo(3)];

% time period
nSamplesPeriod = 20;
Tperiod = 2*pi*sqrt((ai)^3/Mu);
T = [0,nPeriods*Tperiod];   
t_interp = linspace(T(1),T(2),nPeriods*nSamplesPeriod);

% integrator
tol = 1e-6;                 % integration tolerances
odeOptions = odeset('RelTol',tol,'AbsTol',[tol tol tol tol tol tol]);
integrator = Integrator(gravityModel{1},@ode45);
integrator = integrator.setOdeOptions(odeOptions);
integrator = integrator.setFrame('inertial');
integrator.omega = Omega;

for i =1:length(gravityModel)
    
    
    integrator.gravityModel = gravityModel{i};         
    
    tic
    [tout{i}{j},xout{i}{j}] = integrator.integrate(T,Xo);  
    toc

    % calc our output metrics
    p{i}{j} = [xout{i}{j}(:,1),xout{i}{j}(:,3),xout{i}{j}(:,5)];
    v{i}{j} = [xout{i}{j}(:,2),xout{i}{j}(:,4),xout{i}{j}(:,6)];

    xoutBFF = inertial2BFF(tout{i}{j},p{i}{j},Omega); 
    
    U = gravityModel{i}.potential(xoutBFF);
    energy{i}{j} = 0.5*sum(v{i}{j}.*v{i}{j},2)+U;
    
    pInterp{i}{j} = interp1(tout{i}{j},p{i}{j},t_interp,'spline');
    eInterp{i}{j} = interp1(tout{i}{j},energy{i}{j},t_interp,'spline');
    

    errorE{i}{j} = 100*abs((eInterp{i}{j}-eInterp{1}{j})./eInterp{1}{j})';
    errorP{i}{j} = 100*vecnorm(pInterp{i}{j}-pInterp{1}{j},2,2)/aRi;
        
        
    for k = 2:length(eInterp{i}{j})
            
        errorE{i}{j}(k) = max(errorE{i}{j}(k),errorE{i}{j}(k-1));
        errorP{i}{j}(k) = max(errorP{i}{j}(k),errorP{i}{j}(k-1));
            
    end

    numFacesStore(i,j) = numFacesi;
    aOverRStore(i,j) = aRi;
    eccStore(i,j) = ecci;
 
end

end

stop
meshCoarse.coordinates=meshCoarse.coordinates/1000;
meshCoarse.resolution=meshCoarse.calculateResolution();

figure(2)
FS = 12
MS=8
hold on
n=1000
plot3(p{1}{j}(end,1)/1000,p{1}{j}(end,2)/1000,p{1}{j}(end,3)/1000,'kd','MarkerFaceColor','k','MarkerSize',MS)
% plot3(p{2}(end,1)/1000,p{2}(end,2)/1000,p{2}(end,3)/1000,'k>','MarkerFaceColor','k','MarkerSize',MS)
% plot3(p{3}(end,1)/1000,p{3}(end,2)/1000,p{3}(end,3)/1000,'ks','MarkerFaceColor','m','MarkerSize',MS)
% plot3(p{4}(end,1)/1000,p{4}(end,2)/1000,p{4}(end,3)/1000,'ko','MarkerFaceColor','b','MarkerSize',MS)
% plot3(p{5}(end,1)/1000,p{5}(end,2)/1000,p{5}(end,3)/1000,'ko','MarkerFaceColor','r','MarkerSize',MS)
% plot3(p{6}(end,1)/1000,p{6}(end,2)/1000,p{6}(end,3)/1000,'ko','MarkerFaceColor','y','MarkerSize',MS)

nPeriods2Plot = 2.0;
nDataPoints = size(p{1}{j},1);
start = floor(nDataPoints-nDataPoints/nPeriods*nPeriods2Plot);
for i = start:nDataPoints
    colorVec = (nDataPoints-i)/(nDataPoints-start)*[1,1,1];
    plot3(p{1}{j}(i-1:i,1)/1000,p{1}{j}(i-1:i,2)/1000,p{1}{j}(i-1:i,3)/1000,'k-','Color',colorVec,'LineWidth',2)
end

plot3(p{1}{j}(end,1)/1000,p{1}{j}(end,2)/1000,p{1}{j}(end,3)/1000,'kd','MarkerFaceColor','k','MarkerSize',MS)
% plot3(p{2}(end,1)/1000,p{2}(end,2)/1000,p{2}(end,3)/1000,'k>','MarkerFaceColor','k','MarkerSize',MS)
% plot3(p{3}(end,1)/1000,p{3}(end,2)/1000,p{3}(end,3)/1000,'ks','MarkerFaceColor','m','MarkerSize',MS)
% plot3(p{4}(end,1)/1000,p{4}(end,2)/1000,p{4}(end,3)/1000,'ko','MarkerFaceColor','b','MarkerSize',MS)
% plot3(p{5}(end,1)/1000,p{5}(end,2)/1000,p{5}(end,3)/1000,'ko','MarkerFaceColor','r','MarkerSize',MS)
% plot3(p{6}(end,1)/1000,p{6}(end,2)/1000,p{6}(end,3)/1000,'ko','MarkerFaceColor','y','MarkerSize',MS)


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


% 
% % filter out our periodic oscilations
% t_interplog = [ones(length(t_interp),1),log10(t_interp')];
% for i = 2:length(eInterp)
%     errorE{i-1} = 100*abs((eInterp{i}-eInterp{1})./eInterp{1})';
%     errorP{i-1} = 100*vecnorm(pInterp{i}-pInterp{1},2,2)/Ro;
%     
%     
%     for j = 2:length(eInterp{i})
%         
%         errorE{i-1}(j) = max(errorE{i-1}(j),errorE{i-1}(j-1));
%         errorP{i-1}(j) = max(errorP{i-1}(j),errorP{i-1}(j-1));
%         
%     end
%     
%     %paddedE = log10(errorE{i-1});
%     %paddedP = log10(errorP{i-1});
%     %betaE{i-1} = t_interplog(1000:5000,:)\paddedE(1000:5000,:);
%     %betaP{i-1} = t_interplog(1000:5000,:)\paddedP(1000:5000,:);
% end
% %figure(1)
% %plot(t_interplog(10000:end),log10(errorE{i-1}(10000:end,:)));
t_interp = t_interp/Tperiod;

figure(3)
LW=2
LWA=3
MS=8
FS=12
YC1=0.9
YC2=0.85
hold on
for i = 1:length(errorE)
    for j = 1:length(errorE{i})
        avgErrorE{i}(:,j) = errorE{i}{j}
        avgErrorP{i}(:,j) = errorP{i}{j}
    end
end
    Color = numFacesStore(2,i)/max(numFacesStore(2,:));
    %Color = aOverRStore(2,i)/max(aOverRStore(2,:));
    %plot(t_interp,errorP{4}{i}./errorP{2}{i},'-','LineWidth',LWA,'Color',[1,Color,0])
    plot(t_interp,mean(avgErrorP{3}./avgErrorP{2},2),'m-','LineWidth',LW)
    plot(t_interp,mean(avgErrorP{4}./avgErrorP{2},2),'b--','LineWidth',LW)
    plot(t_interp,mean(avgErrorP{5}./avgErrorP{2},2),'r:','LineWidth',LW)
    plot(t_interp,mean(avgErrorP{6}./avgErrorP{2},2),'y-.','Color',[YC1,YC2,0],'LineWidth',LW)

box on

%axis([1,1000,1,1000])
xlabel('time (periods)','interpreter','latex'); 
ylabel('position error ratio','interpreter','latex'); 
set(gcf,'Color',[1,1,1]);
set(gca,'yscale','log');set(gca,'xscale','log');
legend({'Analytic-8000','P1-8000','P2-8000','P3-8000','P4-8000'},'interpreter','latex','location','NorthWest')
legend boxoff 
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')

figure(4)
LW=2
LWA=3
MS=8
FS=12
YC1=0.9
YC2=0.85
hold on
plot(t_interp,mean(avgErrorE{3}./avgErrorE{2},2),'m-','LineWidth',LW)
plot(t_interp,mean(avgErrorE{4}./avgErrorE{2},2),'b--','LineWidth',LW)
plot(t_interp,mean(avgErrorE{5}./avgErrorE{2},2),'r:','LineWidth',LW)
plot(t_interp,mean(avgErrorE{6}./avgErrorE{2},2),'y-.','Color',[YC1,YC2,0],'LineWidth',LW)

box on

%axis([1,1000,1,1000])
xlabel('time (periods)','interpreter','latex'); 
ylabel('energy error ratio','interpreter','latex'); 
set(gcf,'Color',[1,1,1]);
set(gca,'yscale','log');set(gca,'xscale','log');
legend({'Analytic-8000','P1-8000','P2-8000','P3-8000','P4-8000'},'interpreter','latex','location','NorthWest')
legend boxoff 
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')


