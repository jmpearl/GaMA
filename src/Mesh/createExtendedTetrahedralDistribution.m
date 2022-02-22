function [ centroids,volumes] = createExtendedTetrahedralDistribution( mesh, numLayers, method)
% generates an extended tet distribution from surface mesh of arb. degree
%=========================================================================%
% Mascon distribution method proposed by T.G.G Chanut, S. Aljbaae, and V.
% Carruba in "Mascon Gravitation Model Using A Shaped Polyhedral Source"
% MNRAS 450,3742-3749 (2015)
%
% triangular facets of a closed surface mesh are used to create tetrahedra
% that extend down to the center of the body. The tetrahedra can then be
% sliced to create layers of prisms surrounding a tetrahedral core.
%
% Their approach has been generalized for curvilinear meshes and allows
% the inner most layer to be sliced at a constant radius and lumped into
% a single central mascon (method='lumpcore').
%--------------------------------------------------------------------------
% Inputs:
%   mesh ------- Surface Mesh or .obj file
%   numLayers -- 1-pure tetrahedra
%                2-tetrahedra surrounded by prisms
%   method ----- "lumpcore" will combine all tets in inner most layer into
%                into a single large mascon
%--------------------------------------------------------------------------
% Outputs:
%   centroids  - centroids of tetrahedra (Nx3)
%   volumes  --- centroids of tetrahedra (Nx3)
%=========================================================================%

    % process our inputs
    assert(nargin==2 || nargin==3,'incorrect number of inputs')
    assert(isa(mesh,"SurfaceMesh"),"first entry must be SurfaceMesh");
    assert(numLayers > 0,"numLayers must be positive integer");
    if nargin == 2
        method = 'standard';
    end
    if strcmpi(method,"lumpcore") && numLayers < 2
        warning("increasing numLayers to 2 for lumpedcore method")
        numLayers=2;
    end

    % were going to need to modify the mesh so make a local copy
    smLocal = SurfaceMesh(mesh);
    numMasconsPerLayer = smLocal.numFaces;

    % get ext. tet. volume for each layer and subtract to get prism volumes
    for i = 1:numLayers
        [comi,voli] = smLocal.extendedTetrahedaCentroids();

        % subtract inner tet from outer tet to get the prism
        if i > 1
            centroids(layerIndices,1:3)=centroids(layerIndices,1:3).*volumes(layerIndices,1) - comi.*voli;
            volumes(layerIndices,1)=volumes(layerIndices,1)-voli;
            centroids(layerIndices,1:3)=centroids(layerIndices,1:3)./volumes(layerIndices,1);
        end

        % add in new layer
        if strcmpi(method,"lumpcore") && i == numLayers
            centroids(end+1,1:3) = voli'*comi/sum(voli);
            volumes(end+1,1) = sum(voli);
        else
            layerIndices = numMasconsPerLayer*(i-1)+1 : numMasconsPerLayer*i;
            centroids(layerIndices,1:3) = comi;
            volumes(layerIndices,1) = voli;
        end

        % surface of next layer
        if strcmpi(method,"lumpcore") && i == numLayers-1
            posMag = vecnorm(smLocal.coordinates,2,2);
            radius = min(posMag)-smLocal.resolution;
            if radius<0.0
                warning("radius < 0 skipping last lumped layer")
                break
            end
            smLocal.coordinates = smLocal.coordinates./posMag * radius;
        else
            smLocal.coordinates = smLocal.coordinates * (numLayers-i)/numLayers;
        end
        
    end

end

