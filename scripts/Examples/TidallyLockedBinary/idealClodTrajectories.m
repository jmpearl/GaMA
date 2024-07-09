clear all; close all; clc;

addpath(genpath('/home/jmpearl/GravityModels/CLEO'))

fileName = 'Didymos_6m.obj';%
%fileName = 'Didymos_25cm_to_4m_GradedMesh.obj';%

% constants (units kg m s)
G = 6.67e-11;
M = 4.8e9;
T = 11.9 * 3600;

omega = -2*pi/T;
Mu = G*M;

intergratorTolerance = 1e-6;

didyM = 540e9;
didyR = 1.2e3; 
didyA = G*didyM/didyR^2;

Rcom = didyR*didyM / (didyM + M);

% create our surface mesh
surf = SurfaceMesh(fileName);
surf.setNumFaces(1000)
surf.coordinates = surf.coordinates * 1000;
surf.coordinates=surf.coordinates -surf.centroid + [Rcom,0.0,0.0];
surf.resetBulkProperties();

% create the gravity model
%--------------------------------------------------------------------------

dimorphosGrav = ApproximatePolyhedralModel(surf,Mu);

%grav = ApproximatePolyhedralModel(surf,Mu);
%didymosGrav = ConstantAccelerationModel([-didyA,0.0,0.0]);
didymosGrav = MasconModel(surf, G*didyM, 1);
didymosGrav.coordinates = [-didyR+Rcom,0,0];

grav = CompositeModel(dimorphosGrav,...
                      didymosGrav);


% integration I.C. and time interval
reductionFactor = linspace(0.1,1.8,50);
alpha = 25;
phi = linspace(0,360,15);

Ro = 1.01*abs(max(surf.coordinates(:,1)-surf.centroid(:,1)));
Vcirc = sqrt(Mu/Ro);
To = [0,15] * 3600;

for j=1:length(phi)
for i=1:length(reductionFactor)

    Vo = Vcirc*reductionFactor(i);
    Rparticle = Rcom+Ro;
    Vparticle = sqrt((M+didyM)*G/Rparticle);
    
    Xo = [Rparticle,Vo*sind(alpha),0,-Vo*cosd(alpha)*cosd(phi(j))-Vparticle,0,Vo*cosd(alpha)*sind(phi(j))];
    %Xo = [Rparticle,Vo*sind(alpha),0,-Vo*cosd(alpha)*cosd(phi(j))+omega*Rparticle,0,Vo*cosd(alpha)*sind(phi(j))];
    %Xo = [Ro,0,0,0,0,0];
    disp([reductionFactor(i),phi(j),Vo])
    disp(Xo)
    % set up the integrator
    odeOptions = odeset('RelTol', intergratorTolerance,  ...
                        'AbsTol',[intergratorTolerance, ...
                                  intergratorTolerance, ...
                                  intergratorTolerance, ...
                                  intergratorTolerance, ...
                                  intergratorTolerance, ...
                                  intergratorTolerance]);
    inte = Integrator(grav);
    %inte.thirdBodyModel = didymosGrav;
    inte.inertialThirdBodyFrame=false;
    inte.odeOptions = odeOptions;
    inte.omega = omega;

    %tic
    [tout,xout] = inte.integrate(To,Xo);
    %toc

    %pos = [xout(:,1),xout(:,3),xout(:,5)];
    %omega2 = sqrt(didyM*G/Rparticle)/Rparticle
    pos = inertial2BFF(tout,[xout(:,1),xout(:,3),xout(:,5)],omega);
    storage{i}{j} = pos;
end
end


figure(1)
hold on

colors = ["m-","r-","y-","g-","c-","b-","k-"];
% plot trajectories
for i=1:length(storage)
    for j = 1:length(storage{i})
    %disp([i,j])
    pos = storage{i}{j};
    isInside = surf.isInside(pos);
    impactId = find(isInside, 1, 'first')-1;
%     plotColor = "r-";
%     if isempty(impactId)
%         disp(Vcirc*reductionFactor(j))
%         impactId = length(pos(:,1));
%         plotColor = "b-";
%     end
    disp(colors(mod(j,length(colors))+1))
    plot3(pos(1:impactId,1),pos(1:impactId,2),pos(1:impactId,3),colors(mod(j,length(colors))+1),"LineWidth",2)
    %plot3(pos(:,1),pos(:,2),pos(:,3),colors(mod(j,length(colors))+1),"LineWidth",2)
    end
end
surf.plot('plotType','userDefined',... 
                'edgeColor','k',...             
                'faceColor','w',...         
                'faceAlpha',0.9,...
                'lineStyle','-',...
                'lighting',[-1,-1,-1],... % specify vec of incoming light ray
                'shadowing','on')      % ray tracing on (default)
plot3(-didyR+Rcom,0,0,'ko')
daspect([1,1,1])
title('impact trajectory')
xlabel('x'); 
ylabel('y'); 
zlabel('z');
ext = 200
axis([-ext+1200,ext+1200,-ext,ext,-ext,ext])

set(gcf,'Color',[1,1,1]);

figure(2)
hold on
% plot trajectories
colors = ["m","r","y","g","c","b","k"];
for i=1:length(storage)
    for j = 1:length(storage{i})

        

    pos = storage{i}{j};
    isInside = surf.isInside(pos);
    impactId = find(isInside, 1, 'first')-1;
    %disp(impactId);
    plot3(pos(impactId,1),pos(impactId,2),pos(impactId,3),'ko','MarkerFaceColor',colors(mod(j,length(colors))+1),"LineWidth",2)
    end
end
surf.plot('plotType','userDefined',... 
                'edgeColor','k',...             
                'faceColor','w',...         
                'faceAlpha',0.9,...
                'lineStyle','-',...
                'lighting',[-1,0.35,0],... % specify vec of incoming light ray
                'shadowing','on')      % ray tracing on (default)
daspect([1,1,1])
title('impact locations')
xlabel('x'); 
ylabel('y'); 
zlabel('z');
set(gcf,'Color',[1,1,1]);
view([0,-1,0])
%[f1, cdata] = myaa([4, 2]); imwrite(cdata, '42_refined.png', 'png');
%legend(')
