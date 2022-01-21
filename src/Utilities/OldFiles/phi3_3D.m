function [phi,dphidu,dphidv,dphidw] = phi3_3D(u,v,w)



Q = [3,0,0,0;...
    0,3,0,0;...
    0,0,3,0;...
    0,0,0,3;...
    2,1,0,0;...
    0,2,1,0;...
    1,0,2,0;...
    2,0,0,1;...
    0,2,0,1;...
    0,0,2,1;...
    1,2,0,0;...
    0,1,2,0;...
    2,0,1,0;...
    1,0,0,2;...
    0,1,0,2;...
    0,0,1,2;...
    1,1,1,0;...
    1,1,0,1;...
    0,1,1,1;...
    1,0,1,1];

X  = Q(:,1)/3*[1,0,0]+Q(:,2)/3*[0,1,0]+Q(:,3)/3*[0,0,1]+Q(:,4)/3*[0,0,0];

u_i = X(:,1); % u coordinates
v_i = X(:,2); % v coordinates
w_i = X(:,3); % v coordinates

N_q = length(u_i); % number of interpolation points

A_i = [u_i.^3,v_i.^3,w_i.^3,u_i.^2.*v_i,v_i.^2.*u_i,...
    v_i.^2.*w_i,w_i.^2.*v_i,u_i.^2.*w_i,w_i.^2.*u_i,w_i.*v_i.*u_i,...
    u_i.^2,v_i.^2,w_i.^2,u_i.*v_i,w_i.*v_i,u_i.*w_i,...
    u_i,v_i,w_i,0*u_i+1];

A_s = [u.^3,v.^3,w.^3,u.^2.*v,v.^2.*u,...
    v.^2.*w,w.^2.*v,u.^2.*w,w.^2.*u,w.*v.*u,...
    u.^2,v.^2,w.^2,u.*v,w.*v,u.*w,...
    u,v,w,0*u+1];



for i = 1:N_q
    b = zeros(N_q,1);
    b(i) = 1;
    c(1:20,i) = A_i\b; % solve for polynomial coefficents
    
end

phi = A_s*c;

dAdu_s = [3*u.^2,0*u,0*u,2*u.*v,v.^2,...
    0*u,0*u,2*u.*w,w.^2,w.*v,...
    2*u,0*u,0*u,v,0*u,w,...
    0*u+1,0*u,0*u,0*u];

dAdv_s = [0*u,3*v.^2,0*u,u.^2,2*v.*u,...
    2*v.*w,w.^2,0*u,0*u,w.*u,...
    0*u,2*v,0*u,u,w,0*u,...
    0*u,0*u+1,0*u,0*u];

dAdw_s = [0*u,0*u,3*w.^2,0*u,0*u,...
    v.^2,2*w.*v,u.^2,2*w.*u,v.*u,...
    0*u,0*u,2*w,0*u,v,u,...
    0*u,0*u,0*u+1,0*u];

dphidu = dAdu_s*c;
dphidv = dAdv_s*c;
dphidw = dAdw_s*c;

end
