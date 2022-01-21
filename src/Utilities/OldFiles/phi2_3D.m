function [phi,phiu,phiv,phiw] = phi2_3D(u,v,w)

    phi = zeros(length(u),10);
    phiu = zeros(length(u),10);
    phiv = zeros(length(u),10);
    phiw = zeros(length(u),10);
    
    c1 = u;
    c2 = v;
    c3 = w;
    c4 = 1-u-v-w;
    
    dc1du = 1;
    %dc2du = 0;
    %dc3du = 0;
    dc4du = -1;
    
    %dc1dv = 0;
    dc2dv = 1;
    %dc3dv = 0;
    dc4dv = -1;
    
    %dc1dw = 0;
    %dc2dw = 0;
    dc3dw = 1;
    dc4dw = -1;
    
    phi(:,1) = 2*c1.^2-c1;
    phi(:,2) = 2*c2.^2-c2;
    phi(:,3) = 2*c3.^2-c3;
    phi(:,4) = 2*c4.^2-c4;
    phi(:,5) = 4*c1.*c2;
    phi(:,6) = 4*c2.*c3;
    phi(:,7) = 4*c3.*c1;
    phi(:,8) = 4*c1.*c4;
    phi(:,9) = 4*c2.*c4;
    phi(:,10) = 4*c3.*c4;
    
    phiu(:,1) = (4*c1-1)*dc1du;
    phiu(:,2) = 0;
    phiu(:,3) = 0;
    phiu(:,4) = (4*c4-1)*dc4du;
    phiu(:,5) = 4*c2*dc1du;
    phiu(:,6) = 0;
    phiu(:,7) = 4*c3*dc1du;
    phiu(:,8) = 4*c4*dc1du + 4*c1*dc4du;
    phiu(:,9) = 4*c2*dc4du;
    phiu(:,10) = 4*c3*dc4du;
    
    phiv(:,1) = 0;
    phiv(:,2) = (4*c2-1)*dc2dv;
    phiv(:,3) = 0;
    phiv(:,4) = (4*c4-1)*dc4dv;
    phiv(:,5) = 4*c1*dc2dv;
    phiv(:,6) = 4*c3*dc2dv;
    phiv(:,7) = 0;
    phiv(:,8) = 4*c1*dc4dv;
    phiv(:,9) = 4*c2*dc4dv+4*c4*dc2dv;
    phiv(:,10) = 4*c3*dc4dv;
    
    phiw(:,1) = 0;
    phiw(:,2) = 0;
    phiw(:,3) = (4*c3-1)*dc3dw;
    phiw(:,4) = (4*c4-1)*dc4dw;
    phiw(:,5) = 0;
    phiw(:,6) = 4*c2*dc3dw;
    phiw(:,7) = 4*c1*dc3dw;
    phiw(:,8) = 4*c1*dc4dw;
    phiw(:,9) = 4*c2*dc4dw;
    phiw(:,10) = 4*c3*dc4dw+4*c4*dc3dw;
    
end

