function [x_s] = cartesian2spherical(x)
%==========================================================================
%   Jason M. Pearl
%   University of Vermont
%   College of Engineering and Mathematical Sciences
%
%   jmpearl@uvm.edu | jasonmpearl91@gmail.com
%==========================================================================
% list of vectors (Nx3) is converted into list of unit vectors
%--------------------------------------------------------------------------
% Inputs:
%   x --------- (Nx3)  position in cartesian coordinates
%--------------------------------------------------------------------------
% Outputs:
%	a --------- (Nx3) position in spherical coordinates (r,theta,phi)
%               theta - longitude, phi - latitude
%==========================================================================

r = sqrt(x(:,1).^2+x(:,2).^2+x(:,3).^2);
theta = atan2(x(:,2),x(:,1));
phi = atan2(x(:,3),sqrt(x(:,1).^2+x(:,2).^2));

x_s=[r,theta,phi];

end

