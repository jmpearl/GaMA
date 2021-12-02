function [X,Y,Z] = createSimpleLattice(minCoordinates,maxCoordinates,numDivisions)
% creates a 3d cartesian grid
%--------------------------------------------------------------------------
% Inputs:
%   minCoordinates -- [minx,miny,minz]
%   maxCoordinates -- [maxx,maxy,maxz]
%   numDivisions ---- number of cells in x,y,z
%--------------------------------------------------------------------------

    if length(minCoordinates)~=3
        error('minCoordinates must be 1x3')
    end
    if length(maxCoordinates)~=3
        error('maxCoordinates must be 1x3')
    end
    if length(numDivisions)~=3
        error('numDivisions must be 1x3')
    end

    x = linspace(minCoordinates(1),maxCoordinates(1),numDivisions(1));
    y = linspace(minCoordinates(2),maxCoordinates(2),numDivisions(2));
    z = linspace(minCoordinates(3),maxCoordinates(3),numDivisions(3));

    [X,Y,Z] = meshgrid(x,y,z);

end