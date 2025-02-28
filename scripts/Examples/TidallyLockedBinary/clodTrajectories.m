clear all; close all; clc;

fileName = 'Didymos_6m.obj';
runname = "run_20241102_3";
Tf = 15;

% figure set 1
%datFileName = 'BustedClods/1xVesc-CLOD-pos-vel-55.091767.dat';

% figure set 2
datFileName = 'BustedClods/2xVesc-CLOD-pos-vel-50.519090.dat';

% load in the clod and rotate to the anti-didymos octant
stride = 1;
unitConversion = 0.01;

headerLength = 2;
theta = 90; 
rotMat = [ cosd(theta), sind(theta),0;...
          -sind(theta), cosd(theta),0;...
           0,           0,          1];

[r,drdt] = readDat(datFileName,headerLength);

r = r(drdt(:,2)<0,:);
drdt = drdt(drdt(:,2)<0,:);

%r = r - center;
r = r * unitConversion;
drdt = drdt*unitConversion;

% rotate the clod
r = r*rotMat;
drdt = drdt*rotMat;
mean(drdt)
min(drdt)
max(drdt)

r = r(1:stride:end,:);
drdt = drdt(1:stride:end,:);

angle = asind(drdt(:,1)./vecnorm(drdt,2,2));

% read in positions (r) and velocities (v)
%[r,drdt] = readDat(datFileName,headerLength);

% rotate the clod
% r = r*rotMat;
% drdt = drdt*rotMat;


% constants (units kg m s)
G = 6.67e-11;           % grav constant
M = 4.8e9;              % mass dimorphos
T = 11.921473 * 3600;   % orbital period

omega = -2*pi/T;   % rotation rate (rad/sec)
Mu = G*M;          % gravitational parameter

intergratorTolerance = 1e-5;

didyM = 540e9;            % mass of didymos
didyR = 1.2e3;            % orbital radius of dimorphos
didyA = G*didyM/didyR^2;  % acceleration didymos on dimorphos

Rcom = didyR*didyM / (didyM + M);   % center of mass location
Ro = 86;
% create our surface mesh -- first we load in the obj, then coarsen to 1000
% faces, then scale km->m, then position it relative to the COM of the
% system, then we recalculate the bulk properties
surf = SurfaceMesh(fileName);
disp(surf.numFaces)
surf.setNumFaces(1000)
surf.coordinates = surf.coordinates * 1000;
surf.coordinates=surf.coordinates -surf.centroid + [Rcom,0.0,0.0]; 
surf.resetBulkProperties();  

% bump out the busted clod positions so they are outside the surface model
r = r + [Rcom+2.5,0.0,0.0];


% create the gravity model
%--------------------------------------------------------------------------
% dimorphos is a polyhedron and didymos a point mass. We stitch these
% together in a composite gravity model which is acceptable b/c things are
% tidally locked.

dimorphosGrav = ApproximatePolyhedralModel(surf,Mu);

didymosGrav = MasconModel(surf, G*didyM, 1);
didymosGrav.coordinates = [-didyR+Rcom,0,0];

grav = CompositeModel(dimorphosGrav,...
                      didymosGrav);


% integration I.C. and time interval
To = [0,Tf] * 3600;
odeOptions = odeset('RelTol', intergratorTolerance,  ...
                        'AbsTol',[intergratorTolerance, ...
                                  intergratorTolerance, ...
                                  intergratorTolerance, ...
                                  intergratorTolerance, ...
                                  intergratorTolerance, ...
                                  intergratorTolerance]);
inte = Integrator(grav);
inte.odeOptions = odeOptions;
inte.omega = omega;

Rparticle = Rcom+Ro;
Vparticle = abs(Rparticle*omega);%sqrt(abs(grav.potential([0,Rparticle,0])));

for i = 1:size(r,1)
    
    
    Xo = [r(i,1), drdt(i,1), r(i,2), drdt(i,2)-1.1*Vparticle, r(i,3), drdt(i,3)];

    tic
    [tout,xout] = inte.integrate(To,Xo);
    toc

    % convert to the body-fixed frame and store in a cell array
    pos = inertial2BFF(tout,[xout(:,1),xout(:,3),xout(:,5)],omega);
    storage{i} = pos;
end


figure(1)
hold on
colors = ["m-","r-","y-","g-","c-","b-","k-"];
% plot trajectories
for i=1:length(storage)

    pos = storage{i};
    isInside = surf.isInside(pos);
    impactId = find(isInside, 1, 'first')-1;
    disp(impactId)
    plotColor = "r-";
    if isempty(impactId)
        
        impactId = length(pos(:,1));
        plotColor = "b-";
    end
    disp(colors(mod(i,length(colors))+1))
    plot3(pos(1:impactId,1),pos(1:impactId,2),pos(1:impactId,3),colors(mod(i,length(colors))+1),"LineWidth",2)
    %plot3(pos(:,1),pos(:,2),pos(:,3),colors(mod(j,length(colors))+1),"LineWidth",2)

end
surf.plot('plotType','userDefined',... 
                'edgeColor','k',...             
                'faceColor','w',...         
                'faceAlpha',0.75,...
                'lineStyle','-',...
                'lighting',[-1,-1,-1],... % specify vec of incoming light ray
                'shadowing','on')      % ray tracing on (default)
plot3(-didyR+Rcom,0,0,'ko')
daspect([1,1,1])
title('impact trajectory')
xlabel('x'); 
ylabel('y'); 
zlabel('z');
ext = 130;
axis([-ext+1280,ext+1280,-ext,ext,-ext,ext])
set(gcf,'Color',[1,1,1]);
[f1, cdata] = myaa([4, 2]); imwrite(cdata, runname+'_impact_trajectories.png', 'png');
close all

figure(2)
hold on
% plot trajectories
colors = ["m","r","y","g","c","b","k"];
for i=1:length(storage)
    pos = storage{i};
    isInside = surf.isInside(pos);
    impactId = find(isInside, 1, 'first')-1;
    %disp(impactId);
    if impactId > 0
        plot3(pos(impactId,1),pos(impactId,2),pos(impactId,3),'ko','MarkerFaceColor','c',"LineWidth",2)
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
[f1, cdata] = myaa([4, 2]); imwrite(cdata, runname+'_impact_locations_view2.png', 'png');
close all

figure(3)
histogram(vecnorm(drdt,2,2))
axis([0,0.1,0,1000])
set(gcf,'Color',[1,1,1]);
xlabel("speed (m/s)")
ylabel("number of particles")
[f1, cdata] = myaa([4, 2]); imwrite(cdata, runname+'_velocity_histogram.png', 'png');
close all

figure(4)
histogram(angle)
%axis([0,60,0,])
set(gcf,'Color',[1,1,1]);
xlabel("ejection angle (degrees)")
ylabel("number of particles")
[f1, cdata] = myaa([4, 2]); imwrite(cdata, runname+'_angle_histogram.png', 'png');