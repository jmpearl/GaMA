function [coordinates] = createLattice(ds,minExtent,maxExtent,latticeType)
% general generator for lattice distributions
%--------------------------------------------------------------------------
% Inputs:
%   ds ---------- spacing of the lattice
%   minExtent --- corner of box with min value in all cartesian directions
%   maxExtent --- corner of box with max value in all cartesian directions
%   latticeType - type of distribution 
%                 pc - primitive cubic lattice
%                 fcc - face centered cubic lattice
%                 bcc - body center cubic lattice
%--------------------------------------------------------------------------
% Outputs:
%   coordinates - Nx3 x,y,z values for each point
%--------------------------------------------------------------------------

    % process inputs
    if nargin < 4 
        latticeType = "pc";
    end
    latticeType = lower(latticeType);
    assert(length(maxExtent)==3,"max extent must be a 1x3")
    assert(length(minExtent)==3,"min extent must be a 1x3")
    assert(all(maxExtent>minExtent),"min extent must be smaller than max")
    assert(any(strcmp(latticeType,["pc","fcc","bcc"])),"lattice type must be: simple, fcc, or bcc") %#ok<STCI> 
    
    % Number of elements per dimension
    numStepsx = round((maxExtent(1)-minExtent(1))/ds);
    numStepsy = round((maxExtent(2)-minExtent(2))/ds);
    numStepsz = round((maxExtent(3)-minExtent(3))/ds);

    % Correct so that ds is constant
    maxExtent = minExtent+[numStepsx+1,numStepsy+1,numStepsz+1]*ds;
    
    % 1D x-y-z coordinates
    x = linspace(minExtent(1),maxExtent(1),numStepsx+2);
    y = linspace(minExtent(2),maxExtent(2),numStepsy+2);
    z = linspace(minExtent(3),maxExtent(3),numStepsz+2);
    
    % get the pc grid
    [X,Y,Z] = meshgrid(x,y,z);
    
    % if we have a bcc or fcc lattice we need to copy and shift
    coordinates = [X(:),Y(:),Z(:)];
    if strcmp(latticeType,"bcc")
        coordinates = [coordinates; c
                       coordinates + [0.5,   0.5, 0.5]*ds];
    elseif strcmpi(latticeType,"fcc")
        coordinates = [coordinates;
                       coordinates + [0.5, 0.5,   0]*ds;
                       coordinates + [0,   0.5, 0.5]*ds;
                       coordinates + [0.5, 0,   0.5]*ds];
    end
end