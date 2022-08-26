function failed = Test_Integrator(failed)
disp('====================================================================')
disp('Testing Integrator')
disp('====================================================================')
disp('Comparing trajectories integrated around Eros in the BFF and ')
disp('intertial frames. Approximate polyhedral, third body, and composite')
disp('gravity models are all implicitly tested.')
disp('--------------------------------------------------------------------')
disp(' ')

format long


% trajectory I.C.s & tolerances
%--------------------------------------------------------------------------
G = 6.67e-11;              % gravitational constant
T = [0,10]*24*3600;        % ten day integration period
Ro = 35000;                % initial radius

errorThreshold = 0.0025;
intergratorTolerance = 1e-8;

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

% set up the gravity models -- here I'm splitting the gravity model into
% two as an example and test of the composite model class.
%--------------------------------------------------------------------------
gravityModel1 = ApproximatePolyhedralModel(mesh,Mu/2);
gravityModel2 = ApproximatePolyhedralModel(mesh,Mu/2);

gravityModel = CompositeModel(gravityModel1,...
    gravityModel2);

a_nominal = gravityModel.acceleration([0,Ro,0]);
thirdBodyModel = ConstantAccelerationModel(0.1*a_nominal);

% Integrate in the Body Fixed Frame
%--------------------------------------------------------------------------
disp(' ')
disp('Body-Fixed Frame Integration')

Vcirc = sqrt(Mu/Ro);       % initially circular
Xo = [Ro,0,0,0,0,Vcirc];   % I.C.s [rx,vx,ry,vy,rz,vz]

Xo = convertState2BFF(0,Xo,Omega);

odeOptions = odeset('RelTol',intergratorTolerance, ...
                    'AbsTol',[intergratorTolerance, ...
                              intergratorTolerance, ...
                              intergratorTolerance, ...
                              intergratorTolerance, ...
                              intergratorTolerance, ...
                              intergratorTolerance]);

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

% convert to BFF frame
%--------------------------------------------------------------------------
posInertial = inertial2BFF(toutInertial,...
    [xoutInertial(:,1),xoutInertial(:,3),xoutInertial(:,5)],...
    Omega);
posBFF = [xoutBFF(:,1),xoutBFF(:,3),xoutBFF(:,5)];


posErrorBFF = 100*norm(posBFF(end,:)-posInertial(end,:))/norm(posBFF(end,:));
disp(' ')
disp(['position error BFF frame      : ',num2str(posErrorBFF),'%'])

% convert to BFF inertial frame and make sure we get similar error
%--------------------------------------------------------------------------
posInertial = BFF2Inertial(toutInertial,posInertial,Omega);
posBFF = BFF2Inertial(toutBFF,posBFF,Omega);

posErrorInertial = 100*norm(posBFF(end,:)-posInertial(end,:))/norm(posBFF(end,:));

disp(['position error intertial frame: ',num2str(posErrorBFF),'%'])


% check to see if we passed
%--------------------------------------------------------------------------
disp(' ')
if posErrorBFF < errorThreshold && posErrorInertial < errorThreshold
    disp('    PASSED: integrator test')
else
    disp(' ')
    disp('    FAILED: integrator test')
    disp(' ')
    failed=true;
end


end