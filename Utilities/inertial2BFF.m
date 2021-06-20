function [xout] = inertial2BFF(t,x,Omega)
    
    c = cos(Omega*t); 
    s = sin(Omega*t);
    
    xout = zeros(size(x,1),3);
    xout(:,1) = c.*x(:,1)+s.*x(:,2);
    xout(:,2) = -s.*x(:,1)+c.*x(:,2);
    xout(:,3) = x(:,3);
    
end

