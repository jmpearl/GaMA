function [phi,dphidu,dphidv] =  LagrangeInterpolantsTriangle( u, v, d)
% lagrange interp polynomials and derviatives for unit right triangle
%==========================================================================
%
%     o  (0,1)              Provides the values of the lagrange 
%     |\                    interpolating polynomials and derivatives for
%     |  \                  a set of user specified coordinates on the 
%   v |    \                unit triangle and a user specified order
%     |      \              of the interpolants. Outputs are organized
%     o--------o            in columns e.g. each interpolant is given a 
%  (0,0)  u    (1,0)        column with its values at u0....un the rows
%
%--------------------------------------------------------------------------
%
%    1 o
%      | \ 4      Ordering of interpolation points
%    2 o   o
%      |   5 \ 6
%    3 o---o---o
%
%--------------------------------------------------------------------------
% Inputs:
%   u ----- parametric coordinate 1
%   v ----- parametric coordinate 2
%   d ----- degree of interpolating polynomial functions  
%--------------------------------------------------------------------------
% Outputs:
%   phi ----- values of interpolants at specified coordinates (u,v)
%   dphidu -- partial derivatives with respect to u
%   dphidv -- partial derivatives with respect to v
%=========================================================================%
% the matrix inversion craps out at degree 15 ish

tiny = 1e-6;

% check bounds
if any(u>1+tiny) || any(v>1+1e-6) || any((1-u-v)>1+tiny) ...
        || any((1-u-v)<-tiny) || any(u<-tiny) || any(v<-tiny)
    error('bounds error: atleast one input outside unit right triangle')
end

numPts = (d+2)*(d+1)/2;

if d == 0
    
    phi = ones(length(u),1);
    dphidu = zeros(length(u),1);
    dphidv = zeros(length(u),1);
    
elseif d == 1
    
    phi   = zeros(length(u),3);
    dphidu = zeros(length(u),3);
    dphidv = zeros(length(u),3);
    
    phi(:,1) = v;
    phi(:,2) = 1-u-v;
    phi(:,3) = u;
    
    dphidv(:,1) = 1;
    dphidv(:,2) = -1;
    dphidv(:,3) = 0;
    
    dphidu(:,1) = 0;
    dphidu(:,2) = -1;
    dphidu(:,3) = 1;
    
elseif d == 2
    
    phi   = zeros(length(u),6);
    dphidu = zeros(length(u),6);
    dphidv = zeros(length(u),6);
    
    phi(:,1) = 2*v.^2-v;
    phi(:,2) = 4*(-v.^2-u.*v+v);
    phi(:,3) = 2*u.^2+2*v.^2+4*u.*v-3*u-3*v+1;
    phi(:,4) = 4*u.*v;
    phi(:,5) = 4*(-u.^2-u.*v+u);
    phi(:,6) = 2*u.^2-u;
    
    dphidv(:,1) = 4*v-1;
    dphidv(:,2) = -8*v-4*u+4;
    dphidv(:,3) = 4*v+4*u-3;
    dphidv(:,4) = 4*u;
    dphidv(:,5) = -4*u;
    dphidv(:,6) = 0;
    
    dphidu(:,1) = 0;
    dphidu(:,2) = -4*v;
    dphidu(:,3) = 4*u+4*v-3;
    dphidu(:,4) = 4*v;
    dphidu(:,5) = -8*u-4*v+4;
    dphidu(:,6) = 4*u-1; 
    
else
    
    [ui,vi] = latticePointsTriangle(d);
    
    % convert to pt index 
    ui = ui*d;
    vi = vi*d;
    
    % monomials and their derivatives
    As = (u.^(ui')) .*  (v.^(vi'));
    dAsdv = ( (u).^(ui') ) .* ((vi').*(( v).^max(vi'-1,0)));
    dAsdu = ( (v).^(vi') ) .* ((ui').*(( u).^max(ui'-1,0)));
     
    % monomial coeffs
    b = zeros(numPts,1);
    coeffMatrix = zeros(numPts,numPts);
    Ai = ((ui/d).^(ui')) .* ((vi/d).^(vi'));
    
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

end

