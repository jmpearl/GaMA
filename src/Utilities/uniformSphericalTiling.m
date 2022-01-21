function [ coordinates, area] = uniformSphericalTiling(numSteps)
% Tiles a sphere with approx equal area elements
%--------------------------------------------------------------------------
% Inputs:
%   numSteps ----- number of latitude bands
%--------------------------------------------------------------------------
% Outputs:
%   coordinates - xyz coordinates
%   areas ------- surface area
%--------------------------------------------------------------------------

    % phi spacing on 1/4 great circle
    ds = pi/2/numSteps;
    phi = linspace(0,pi/2,numSteps+1);

    % preprocess more than is needed
    Aq = zeros(8*numSteps^2,1);
    phiq = zeros(8*numSteps^2,1);
    thetaq = zeros(8*numSteps^2,1);
    N0 = 1;
    N1 = 0;

    % loop through and create tiles of approximate 1-1 aspect ratio for each
    % phi band, spacing based on phi value at centroid.
    for i = 1:numSteps
        phi2 = phi(i+1);
        phi1 = phi(i);
        phic = (phi2*sin(phi2)-phi1*sin(phi1)+cos(phi2)-cos(phi1))/(sin(phi2)-sin(phi1));
        N_step = round(2*pi*cos(phic)/ds);
        N1 = N1+N_step;

        dtheta = 2*pi/N_step;
        A = dtheta*(sin(phi2)-sin(phi1));

        thetaq(N0:N1,1) = linspace(dtheta,2*pi,N_step)'-dtheta/2;
        phiq(N0:N1,1) = phic;
        Aq(N0:N1,1) = A;

        N0 = N1+1;
    end

    % flip for lower half of unit sphere
    thetaq(N0:2*N1,1) = thetaq(1:N1,1);
    phiq(N0:2*N1,1) = -phiq(1:N1,1);
    Aq(N0:2*N1,1) = Aq(1:N1,1);

    % delete unused memory slots
    thetaq(2*N1+1:8*numSteps^2,:)=[];
    phiq(2*N1+1:8*numSteps^2,:)=[];
    Aq(2*N1+1:8*numSteps^2,:)=[];

    coordinates = [cos(phiq).*cos(thetaq),cos(phiq).*sin(thetaq),sin(phiq)];
    area = Aq;
end

