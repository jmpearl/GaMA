function [ phi,phi_u,phi_v ] = phi2_new( u,v )
% quadratic interpolation functions for the unit triangle 

phi   = zeros(length(u),6);
phi_u = zeros(length(u),6);
phi_v = zeros(length(u),6);

phi(:,1) = 2*v.^2-v;
phi(:,2) = 2*u.^2+2*v.^2+4*u.*v-3*u-3*v+1;
phi(:,3) = 2*u.^2-u;
phi(:,4) = 4*(-v.^2-u.*v+v);
phi(:,5) = 4*(-u.^2-u.*v+u);
phi(:,6) = 4*u.*v;

phi_v(:,1) = 4*v-1;
phi_v(:,2) = 4*v+4*u-3;
phi_v(:,3) = 0;
phi_v(:,4) = -8*v-4*u+4;
phi_v(:,5) = -4*u;
phi_v(:,6) = 4*u;

phi_u(:,1) = 0;
phi_u(:,2) = 4*u+4*v-3;
phi_u(:,3) = 4*u-1;
phi_u(:,4) = -4*v;
phi_u(:,5) = -8*u-4*v+4;
phi_u(:,6) = 4*v;

end

