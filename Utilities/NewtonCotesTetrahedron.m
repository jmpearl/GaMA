function [u,v,w,weights] =  NewtonCotesTetrahedron( Order )
% Newton Cotes quadrature rules for tetrahedra
%=========================================================================%
% Function that returns the computational nodes and integration weights for
% an degree-d langrange polynomial interpolation scheme.
%-------------------------------------------------------------------------%
% Inputs:
%   Order - Order of interpolating polynomial functions
%
% Outputs:
%   u ------ natural coordinate 1 of computational nodes
%   v ------ natural coordinate 2 of computational nodes
%   w ------ natural coordinate 3 of compuational nodes
%   weight - weights of quadrature points
%
%=========================================================================%


[u,v,w]=latticePointsTetrahedron(Order);

if Order == 0
    
    weights = 1;
    
elseif Order == 1
    
    weights = [1/4,1/4,1/4,1,4]';
    
elseif Order == 2
    w1 = -1/20;
    w2 = 1/5;
    
    weights = [w1,w2,w1,w2,w2,w1,w2,w2,w2,w1]';
    
elseif Order == 3
    w1 = 1/40;
    w2 = 0;
    w3 = 9/40;
    
    weights = [w1,w2,w2,w1,w2,w3,w2,w2,w2,w1,w2,w3,w2,w3,w3,w2,w2,w2,w2,w1]';
    
elseif Order == 4
    w1 = -1/84;
    w2 =  4/105;
    w3 = -1/35;
    w4 =  32/105;
    
    weights= [w1,w2,w3,w2,w1,w2,w2,w2,w2,w3,w2,w3,w2,w2,w1,w2,w2,...
              w2,w2,w2,w4,w2,w2,w2,w2,w3,w2,w3,w2,w2,w3,w2,w2,w2,w1]';
          
elseif Order == 5
    
    w1 =  11/1344;
    w2 = -5/576;
    w3 =  5/576;
    w4 =  275/4032;
    w5 = -25/1344;
    w6 = 125/1344;
    
    weights= [w1,w2,w3,w3,w2,w1,w2,w4,w5,w4,w2,w3,w5,w5,w3,w3,w4,w3,...
              w2,w2,w1,w2,w4,w5,w4,w2,w4,w6,w6,w4,w5,w6,w5,w4,w4,w2,...
              w3,w5,w5,w3,w5,w6,w5,w5,w5,w3,w3,w4,w3,w4,w4,w3,w2,w2,w2,w1]';

elseif Order == 6
    
    w1 = -1/200;
    w2 =  3/175;
    w3 = -3/140;
    w4 =  1/35;
    w5 =  0;
    w6 =  3/140;
    w7 = -9/280;
    w8 =  9/70;
    
    weights= [w1,w2,w3,w4,w3,w2,w1,w2,w5,w6,w6,w5,w2,w3,w6,w7,w6,w3,...
              w4,w6,w6,w4,w3,w5,w3,w2,w2,w1,w2,w5,w6,w6,w5,w2,w5,w8,...
              w5,w8,w5,w6,w5,w5,w6,w6,w8,w6,w5,w5,w2,w3,w6,w7,w6,w3,...
              w6,w5,w5,w6,w7,w5,w7,w6,w6,w3,w4,w6,w6,w4,w6,w8,w6,w6,...
              w6,w4,w3,w5,w3,w5,w5,w3,w2,w2,w2,w1]';
else
    disp('Polynomial Order not supported')
    disp('    Valid Options: 0 - 6')
    disp('    ps 0 is really a midpoint rule')
end

    
end

