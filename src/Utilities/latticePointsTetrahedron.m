function [u,v,w] = latticePointsTetrahedron(d)
 
% centroid model for d=0
if d == 0
    u = 0.25;
    v = 0.25;
    w = 0.25;
    
% all other degrees
else
    numPts = (d+3)*(d+2)*(d+1)/6;
    
    u = zeros(numPts,1);
    v = zeros(numPts,1);
    w = zeros(numPts,1);
    
    % initialize
    i = 0;
    j = 0;
    k = 0;
    l = d;
    iter = 1;
    
    while i <= d
        
        u(iter) = i/d;
        v(iter) = j/d;
        w(iter) = k/d;
        
        % iteration logic
        if l>0
            l=l-1;
            k=k+1;
        elseif i+j+k==d && i+j<d
            k=0;
            j=j+1;
            l=d-(i+j+k);
        elseif i+j==d && i<d
            i=i+1;
            j=0;
            l=d-(i+j+k);
        else
            break
        end
        iter = iter+1;
        
    end
    
end

