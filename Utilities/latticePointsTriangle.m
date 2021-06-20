function [u,v] = latticePointsTriangle(d)
    
    numPts = (d+2)*(d+1)/2;
    
    u = zeros(numPts,1);
    v = zeros(numPts,1);
    
    i1 = 1;
    for i = 0:d
        i2 = i1+d-i;
        u(i1:i2) =  i/d*ones(1,d-i+1);
        v(i1:i2) = linspace(d-i,0,d-i+1)/d;
        i1 = i2+1;
    end
    
end

