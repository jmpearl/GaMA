function [phi,dphidu,dphidv,dphidw] =  LagrangeInterpolantsTetrahedron( u, v, w, d)
% lagrange interp polynomials and derviatives for unit right triangle
%==========================================================================
% Inputs:
%   u ----- parametric coordinate 1
%   v ----- parametric coordinate 2
%   w ----- parametric coordinate 3
%   d ----- degree of interpolating polynomial functions  
%--------------------------------------------------------------------------
% Outputs:
%   phi ----- values of interpolants at specified coordinates (u,v,w)
%   dphidu -- partial derivatives with respect to u
%   dphidv -- partial derivatives with respect to v
%   dphidw -- partial derivatives with respect to w
%=========================================================================%
% I wouldn't push this too far past degree 6

tiny = 1e-6;

% check bounds
if any(u>1+tiny) || any(v>1+1e-6) || any(w>1+1e-6) || any((1-u-v-w)>1+tiny) ...
        || any((1-u-v-w)<-tiny) || any(u<-tiny) || any(v<-tiny) || any(w<-tiny)
    error('bounds error: atleast one input outside unit right tetrahedron')
end

numPts = (d+3)*(d+2)*(d+1)/6;

if d == 0
    
    phi = ones(length(u),1);
    dphidu = zeros(length(u),1);
    dphidv = zeros(length(u),1);
    dphidw = zeros(length(u),1);
    
else
     
    [ui,vi,wi] = latticePointsTetrahedron(d);
    
    % convert to pt index 
    ui = ui*d;
    vi = vi*d;
    wi = wi*d;
    
    % monomials and their derivatives
    As = (u.^(ui')) .*  (v.^(vi')) .*  (w.^(wi'));
    dAsdu = ( (v).^(vi') ) .*  (w.^(wi')) .* ((ui').*(( u).^max(ui'-1,0)));
    dAsdv = ( (u).^(ui') ) .*  (w.^(wi')) .* ((vi').*(( v).^max(vi'-1,0)));
    dAsdw = ( (v).^(vi') ) .*  (u.^(ui')) .* ((wi').*(( w).^max(wi'-1,0)));
      
    % monomial coeffs
    b = zeros(numPts,1);
    coeffMatrix = zeros(numPts,numPts);
    Ai = ((ui/d).^(ui')) .* ((vi/d).^(vi')) .*  ((wi/d).^(wi'));
    
    for i = 1:numPts
        
        if i > 1
            b(i-1)=0;
        end
        
        b(i) = 1;
        coeffMatrix(1:numPts,i) = Ai\b;
        
    end
    
    phi = As*coeffMatrix;
    dphidu = dAsdu*coeffMatrix;
    dphidv = dAsdv*coeffMatrix;
    dphidw = dAsdw*coeffMatrix;
   
end
 
% slight corrections
phi = phi ./ sum(phi,2);

weightdu = abs(dphidu)./sum(abs(dphidu),2);
sumdu = sum(dphidu,2);
dphidu = dphidu - sumdu.*weightdu;

weightdv = abs(dphidv)./sum(abs(dphidv),2);
sumdv = sum(dphidv,2);
dphidv = dphidv - sumdv.*weightdv;

weightdw = abs(dphidw)./sum(abs(dphidw),2);
sumdw = sum(dphidw,2);
dphidw = dphidw - sumdw.*weightdw;

end



