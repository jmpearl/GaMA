function [xout] = convertState2BFF(t,x,Omega)
    
    if size(x,2) == 6
        
        c = cos(Omega*t);
        s = sin(Omega*t);
        
        Rot = [c, s, 0;
              -s, c, 0;
               0, 0, 1]; 
           
        dRotdt = [-Omega*s, Omega*c, 0;
            -Omega*c, -Omega*s, 0;
            0,        0,       0];
        
        r = [x(:,1),x(:,3),x(:,5)];
        v = [x(:,2),x(:,4),x(:,6)];
        
        r = r*Rot';
        v = (v*Rot')+(r*dRotdt');
        
        xout=[r(:,1),v(:,1),r(:,2),v(:,2),r(:,3),v(:,3)];
    end
end

