%==========================================================================
% Integrates trajectory in the Body-fixed frame and Inertial frame and then
% plots the result in both frames.
%==========================================================================
clear all; close all; clc;

% trajectory I.C.s
%--------------------------------------------------------------------------
G = 6.67e-11;              % gravitational constant
T = [0,50]*24*3600;        % ten day integration period
Ro = 35000;                % initial radius

% load body data
%--------------------------------------------------------------------------
load('Eros.mat');                                 
Mu = bodyProperties.mass*G;
Omega = bodyProperties.rotationRate;

% set up the surface mesh
%--------------------------------------------------------------------------
meshFile = 'Eros_7624.obj';    % define surface mesh file
mesh = SurfaceMesh(meshFile);  % create tri surface mesh
mesh.coarsen(2000);            % coursen to 2000 faces

% set up the gravity models
%--------------------------------------------------------------------------
gravityModel1 = ApproximatePolyhedralModel(mesh,Mu/2);
gravityModel2 = ApproximatePolyhedralModel(mesh,Mu/2);

gravityModel = CompositeModel(gravityModel1,...
                              gravityModel2);       % superpose models

a_nominal = gravityModel.acceleration([0,Ro,0]);
thirdBodyModel = ConstantAccelerationModel(0.1*a_nominal);

% Integrate in the Body Fixed Frame
%--------------------------------------------------------------------------
disp('Body-Fixed Frame Integration')

Vcirc = sqrt(Mu/Ro);       % initially circular
Xo = [Ro,0,0,0,0,Vcirc];   % I.C.s [rx,vx,ry,vy,rz,vz]
tol = 1e-8;                % integration tolerances

Xo = convertState2BFF(0,Xo,Omega);

odeOptions = odeset('RelTol',tol,'AbsTol',[tol tol tol tol tol tol]);

integrator = Integrator(gravityModel,@ode45);
integrator.thirdBodyModel = thirdBodyModel;
integrator.setOdeOptions(odeOptions);
integrator.setFrame('BFF'); % defaults to inertial
integrator.omega = Omega;

tic
[toutBFF,xoutBFF] = integrator.integrate(T,Xo);
toc

disp(' ')

% Integrate in the Inertial Frame
%--------------------------------------------------------------------------
disp('Inertial Frame Integration')

integrator.setFrame('inertial');

Xo = [Ro,0,0,0,0,Vcirc];   

tic
[toutInertial,xoutInertial] = integrator.integrate(T,Xo);
toc

% plot results in BFF frame
%--------------------------------------------------------------------------
posInertial = inertial2BFF(toutInertial,...
                           [xoutInertial(:,1),xoutInertial(:,3),xoutInertial(:,5)],...
                            Omega);
posBFF = [xoutBFF(:,1),xoutBFF(:,3),xoutBFF(:,5)];

figure(1)
hold on
% plot trajectories
plot3(posBFF(:,1),posBFF(:,2),posBFF(:,3),'k-')
plot3(posInertial(:,1),posInertial(:,2),posInertial(:,3),'r-')

% plot mesh smoothness metric
mesh.plot('alpha') 

% plot final positions
plot3(posBFF(end,1),posBFF(end,2),posBFF(end,3),'kx')
plot3(posInertial(end,1),posInertial(end,2),posInertial(end,3),'rs')

daspect([1,1,1])
title('Body-Fxied Frame')
xlabel('x'); 
ylabel('y'); 
zlabel('z');
set(gcf,'Color',[1,1,1]);
legend('BFF','Inertial')

% plot results in Inertial frame
%--------------------------------------------------------------------------
posBFF = BFF2Inertial(toutBFF,...
                      [xoutBFF(:,1),xoutBFF(:,3),xoutBFF(:,5)],...
                       Omega);
posInertial = [xoutInertial(:,1),xoutInertial(:,3),xoutInertial(:,5)];

figure(2)
hold on
% plot trajectories
plot3(posBFF(:,1),posBFF(:,2),posBFF(:,3),'k-')
plot3(posInertial(:,1),posInertial(:,2),posInertial(:,3),'r-')

% plot mesh smoothness metric
mesh.plot('alpha') 

% plot final positions
plot3(posBFF(end,1),posBFF(end,2),posBFF(end,3),'kx')
plot3(posInertial(end,1),posInertial(end,2),posInertial(end,3),'rs')

daspect([1,1,1])
title('Inertial')
xlabel('x'); 
ylabel('y'); 
zlabel('z');
set(gcf,'Color',[1,1,1]);
legend('BFF','Inertial')
