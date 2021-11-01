classdef VolumeMesh < handle
%==========================================================================
% GaMA - Gravity and Mesh Adaptation Libray for Asteroids and Comets
% J.M.Pearl 
% Oct 2021
%--------------------------------------------------------------------------
% Simple array based tetrahedral mesh class. This as designed as a mascon
% distribution generator -- using numerical quadrature to get the locations
% and masses. 
% 
% it is assumed that the surface mesh vertices are the first set of nodes 
% followed by internal vertices and then any higher degree nodes.
%--------------------------------------------------------------------------
% Abreviations:
%       e ----- unique edge
%       f ----- face
%       p ----- index of a vertex
%       coord - vertex coordinates 
%       v ----- vector (typically defining an edge)
%       c ----- centroid (face or volume)
%       A ----- area vector
%       n ----- unit normal
%--------------------------------------------------------------------------
    properties (GetAccess=public)
        coordinates         % coordinates of nodes
        cells               % cells (node indices)
        
        surfaceMesh         % surfaceMesh definition
        isBoundaryNode      % bool array indicating if its a BC node
        
        cellFields          % cell array of stored fields
        nodeFields          % cell array of stored fields 
        
        hasSurfaceMesh      % bool to track if we have a surf def
        isCurved            % logical to track if mesh is curvilinear
        degree              % polynomial degree of mesh
        volume              % total volume
        surfaceArea         % total surface area
        centroid            % centroid of enclosed volume
        resolution          % mean resolution length scale of mesh
        
        numBoundaryVertices % number of vertices on surface
        numBoundaryNodes    % number of nodes on surface
        numNodes            % number of nodes for degree>1 meshes
        numVertices         % number of vertices
        numCells            % number of half edges
        numNodeFields       % number of node fields
        numCellFields       % number of face fields
        
    end
    methods
        function obj = VolumeMesh(mesh)
        % constructor to set things up from surface or vol mesh
        %------------------------------------------------------------------
            if nargin == 1
                if isa(mesh,'SurfaceMesh')
                    obj.surfaceMesh=mesh;
                    obj.hasSurfaceMesh=true;

                    obj.numNodeFields = 0;
                    obj.numCellFields = 0;
                    obj.numCells = 0;
                    obj.numVertices = 0;
                    obj.numNodes = 0;
                    obj.numBoundaryVertices = 0;
                    obj.numBoundaryNodes = 0;

                    obj.isCurved=false;
                    obj.degree = 1;

                elseif isa(mesh,'VolumeMesh')
                    obj.copy(mesh);
                end
            elseif nargin ~= 0
                error("incorrection number of inputs")
            end
        end
        function copy(obj,mesh)
        % deep copy from another volume mesh
        %------------------------------------------------------------------
            
            obj.coordinates = mesh.coordinates;
            obj.cells = mesh.cells;
            obj.surfaceMesh = mesh.surfaceMesh;
            obj.isBoundaryNode = mesh.isBoundaryNode;
            
            obj.cellFields = mesh.cellFields;
            obj.nodeFields = mesh.nodeFields;

            obj.hasSurfaceMesh = mesh.hasSurfaceMesh;
            obj.isCurved = mesh.isCurved;
            obj.degree = mesh.degree;
            obj.volume = mesh.volume;
            obj.surfaceArea = mesh.surfaceArea;
            obj.centroid = mesh.centroid;
            obj.resolution = mesh.resolution;

            obj.numBoundaryVertices = mesh.numBoundaryVertices;
            obj.numBoundaryNodes = mesh.numBoundaryNodes;
            obj.numNodes = mesh.numNodes;
            obj.numVertices = mesh.numVertices;
            obj.numCells = mesh.numCells;
            obj.numNodeFields = mesh.numNodeFields;
            obj.numCellFields = mesh.numCellFields;

        end
        function clear(obj)
        % resets the volume mesh
        %------------------------------------------------------------------
            clear obj.coordinates obj.cells obj.surfaceMesh
            clear obj.isBoundaryNode

            obj.clearFields()

            obj.hasSurfaceMesh=false;
            obj.isCurved=false;
            obj.degree=1
            obj.volume=0;
            obj.surfaceArea=0;
            obj.centroid=0;
            obj.resolution=0;

            obj.numBoundaryVertices=0;
            obj.numBoundaryNodes=0;
            obj.numNodes=0;
            obj.numVertices=0;
            obj.numCells=0;
            obj.numNodeFields=0;
            obj.numCellFields=0;

        end

        function initializeFromPointCloud(obj,internalCoordinates)
        % constructions mesh from set of internal vertices
        %------------------------------------------------------------------
        % Inputs:
        %   internalCoordinates -- internal mesh vertices
        %------------------------------------------------------------------

            obj.coordinates = [obj.surfaceMesh.coordinates;internalCoordinates];
            obj.numVertices = size(obj.coordinates,1);
            obj.numNodes = size(obj.coordinates,1);
            obj.delaunayTriangulation();
            obj.clipExternalCells();

            obj.isBoundaryNode = zeros(obj.numNodes,1);
            obj.isBoundaryNode(1:obj.surfaceMesh.numVertices) = 1;
            obj.numBoundaryVertices = sum(obj.isBoundaryNode);
            obj.numBoundaryNodes = obj.numBoundaryVertices;

            obj.initializeBulkProperties;

        end
        function initializeFromSimpleLattice(obj,numInternalVertices)
        % initializes a volume mesh w/ cubic lattice point distribution
        %------------------------------------------------------------------
        % Inputs:
        %   numInternalVertices -- approximate number of internal vertices
        %------------------------------------------------------------------

            % approximate lattice spacing
            ds = (obj.surfaceMesh.volume/numInternalVertices)^(1/3);

            maxExtent=max(obj.surfaceMesh.coordinates,[],1) + 0.5*ds;
            minExtent=min(obj.surfaceMesh.coordinates,[],1) - 0.5*ds; 

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

            % get the grid
            [X,Y,Z] = meshgrid(x,y,z);

            % create a buffer 
            insetSurfaceMesh = obj.surfaceMesh.offsetSurfaceMesh(-obj.surfaceMesh.resolution/2, ...
                                                                  obj.surfaceMesh.numVertices);
            % clip to internal
            candidates = [X(:),Y(:),Z(:)];
            internalNodes = insetSurfaceMesh.isInside(candidates)==1;
            candidates = candidates(internalNodes,:);

            obj.initializeFromPointCloud(candidates);

        end
        function initializeFromBCCLattice(obj,numInternalVertices)
        % initializes a volume mesh w/ BCC lattice point distribution
        %------------------------------------------------------------------
        % Inputs:
        %   numInternalVertices -- approximate number of internal vertices
        %------------------------------------------------------------------

            % approximate lattice spacing
            ds = (2*obj.surfaceMesh.volume/numInternalVertices)^(1/3);

            maxExtent=max(obj.surfaceMesh.coordinates,[],1) + 0.5*ds;
            minExtent=min(obj.surfaceMesh.coordinates,[],1) - 0.5*ds; 

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

            % get the grid
            [X,Y,Z] = meshgrid(x,y,z);

            % create a buffer 
            insetSurfaceMesh = obj.surfaceMesh.offsetSurfaceMesh(-obj.surfaceMesh.resolution/2, ...
                                                                  obj.surfaceMesh.numVertices);
            % clip to internal
            candidates = [X(:),Y(:),Z(:)];
            candidates = [candidates; candidates + [1,1,1]*0.5*ds];
            internalNodes = insetSurfaceMesh.isInside(candidates)==1;
            candidates = candidates(internalNodes,:);

            obj.initializeFromPointCloud(candidates);
        end
        function initializeFromFCCLattice(obj,numInternalVertices)
        % initializes a volume mesh w/ BCC lattice point distribution
        %------------------------------------------------------------------
        % Inputs:
        %   numInternalVertices -- approximate number of internal vertices
        %------------------------------------------------------------------

            % approximate lattice spacing
            ds = (4*obj.surfaceMesh.volume/numInternalVertices)^(1/3);

            maxExtent=max(obj.surfaceMesh.coordinates,[],1) + 0.5*ds;
            minExtent=min(obj.surfaceMesh.coordinates,[],1) - 0.5*ds; 

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

            % get the grid
            [X,Y,Z] = meshgrid(x,y,z);

            % create a buffer 
            insetSurfaceMesh = obj.surfaceMesh.offsetSurfaceMesh(-obj.surfaceMesh.resolution/2, ...
                                                                  obj.surfaceMesh.numVertices);
            % clip to internal
            candidates = [X(:),Y(:),Z(:)];
            candidates = [candidates; 
                          candidates + [0.5, 0,   0]  *ds;
                          candidates + [0,   0.5, 0]  *ds;
                          candidates + [0,   0,   0.5]*ds];

            internalNodes = insetSurfaceMesh.isInside(candidates)==1;
            candidates = candidates(internalNodes,:);

            obj.initializeFromPointCloud(candidates);
        end
        function initializeFromSurfaceIteration(obj,ratio)
        % initialize volume mesh by iteratively offseting the surface 
        %------------------------------------------------------------------
        % Inputs:
        %   ratio --- ratio numVertices 
        %------------------------------------------------------------------

            %assert (obj.hasSurfaceMesh==1) "no surface mesh loaded"

            % default ratio 
            if nargin == 1 || ratio > 0.95
                ratio = 1/2;
            end

            % copy our surface mesh
            tempMesh = SurfaceMesh(obj.surfaceMesh);
            tempMesh.coarsenOptions.reproject=false;
            tempMesh.coarsenOptions.method = 'uniform';
            tempMesh.smoothOptions.reproject=false;

            % add the origin
            newCoordinates = [0,0,0];

            
            % offset the surface mesh inward coarsening along the way
            while tempMesh.volume > 0.1*obj.surfaceMesh.volume

                % reduce vertex count 
                newVertexCount = round(tempMesh.numVertices*ratio);
                tempMesh.setNumVertices(newVertexCount);

                % offset mesh
                tempMesh = tempMesh.offsetSurfaceMesh(-tempMesh.resolution);

                % we don't want to reproject things
                tempMesh.coarsenOptions.reproject=false;
                tempMesh.coarsenOptions.method = 'uniform';
                tempMesh.smoothOptions.reproject=false;

                % smooth it out 
                tempMesh.smooth(10,'cotan',false);

                % prevent slivers meshes
                if tempMesh.volume< 0.05*obj.surfaceMesh.volume ||...
                   obj.surfaceMesh.numFaces<100
                    break
                end

                % add our points
                newCoordinates = [newCoordinates;tempMesh.coordinates];
            end
            obj.initializeFromPointCloud(newCoordinates);
        end
        function initializeFromOctree(obj,numInitialInternalVertices,...
                                          numLayers,...
                                          smoothingFactor)
        % octree distributed internal points
        %------------------------------------------------------------------
        % Inputs:
        %   numLayers ----------------- number of refinement steps
        %   numInitialIntenalVertices - initial grid size
        %   smoothingFactor ----------- in [1,2] controls slip proximity
        %------------------------------------------------------------------
            
            % defaults
            if nargin <= 2
                numLayers = 2;
            end
            if nargin <= 3
                smoothingFactor = 1;
            end

            % approximate lattice spacing
            ds = (obj.surfaceMesh.volume/numInitialInternalVertices)^(1/3);

            maxExtent=max(obj.surfaceMesh.coordinates,[],1) + 0.5*ds;
            minExtent=min(obj.surfaceMesh.coordinates,[],1) - 0.5*ds; 

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

            % get the grid
            [X,Y,Z] = meshgrid(x,y,z);

            % 1d-ify
            X=X(:);
            Y=Y(:);
            Z=Z(:);

            % Begin: Octree Refinement and Mascon Distribution Generation
            %-------------------------------------------------------------------------%

            %Initialize
            %Octree_Level =[];
            cm=[];

            for j = 1:numLayers
                % Only keep mascons located within the body
                % location of vertices relative to centroid dependent upon otree level
                Offset = 1/2^j*ds*smoothingFactor;

                if j<=numLayers

                    InOrOutVar = zeros(length(X),8);

                    % returns 0 for outside 1 for inside
                    [ InOrOutVar(:,1) ] = obj.surfaceMesh.isInside( [X+Offset,Y+Offset,Z+Offset] );
                    [ InOrOutVar(:,2) ] = obj.surfaceMesh.isInside( [X-Offset,Y+Offset,Z+Offset] );
                    [ InOrOutVar(:,3) ] = obj.surfaceMesh.isInside( [X+Offset,Y-Offset,Z+Offset] );
                    [ InOrOutVar(:,4) ] = obj.surfaceMesh.isInside( [X-Offset,Y-Offset,Z+Offset] );
                    [ InOrOutVar(:,5) ] = obj.surfaceMesh.isInside( [X+Offset,Y+Offset,Z-Offset] );
                    [ InOrOutVar(:,6) ] = obj.surfaceMesh.isInside( [X-Offset,Y+Offset,Z-Offset] );
                    [ InOrOutVar(:,7) ] = obj.surfaceMesh.isInside( [X+Offset,Y-Offset,Z-Offset] );
                    [ InOrOutVar(:,8) ] = obj.surfaceMesh.isInside( [X-Offset,Y-Offset,Z-Offset] );
                    
                    % 8 -- fully internal
                    % 0 -- fully external
                    % inbetween is a bc cell
                    IOVar = sum(InOrOutVar,2);
                    
                    % Add cubes to the octree distribution if they lie 
                    % completely inside the body of interest
                    cm = [cm;X(IOVar>7.5),Y(IOVar>7.5),Z(IOVar>7.5)];

                    if j < numLayers
                        %Cubes that the domain boundary passes through
                        X = X(IOVar>.5 & IOVar<7.5);
                        Y = Y(IOVar>.5 & IOVar<7.5);
                        Z = Z(IOVar>.5 & IOVar<7.5);
    
                        %Subdivide boundary cube into 8 smaller cubes
                        C=1/2^(j+1);
                        cm_next = [X+C*ds,Y+C*ds,Z+C*ds;...
                                   X-C*ds,Y+C*ds,Z+C*ds;...
                                   X+C*ds,Y-C*ds,Z+C*ds;...
                                   X-C*ds,Y-C*ds,Z+C*ds;...
                                   X+C*ds,Y+C*ds,Z-C*ds;...
                                   X-C*ds,Y+C*ds,Z-C*ds;...
                                   X+C*ds,Y-C*ds,Z-C*ds;...
                                   X-C*ds,Y-C*ds,Z-C*ds];
    
                        %Redefine x,y,z for next iteration of loop
                        X = cm_next(:,1);
                        Y = cm_next(:,2);
                        Z = cm_next(:,3);
                    end

                else
                    % For the last refinement step determine if centroid 
                    % of cube is internal or external
                    %InOrOutVar = obj.surfaceMesh.isInside( [X,Y,Z] );

                    %IOVar = sum(InOrOutVar,2);

                    %cm = [cm;X(IOVar>0.5),Y(IOVar>0.5),Z(IOVar>0.5)];
                end

            end
            obj.initializeFromPointCloud(cm);
            obj.initializeBulkProperties();
        end
        function initializeBulkProperties(obj)
        % recalculate stored bulk properties
        %------------------------------------------------------------------

            obj.volume = obj.calculateVolume;
            obj.centroid = obj.calculateCentroid;
            obj.resolution = obj.calculateResolution;

            if obj.hasSurfaceMesh
                obj.surfaceArea = obj.surfaceMesh.surfaceArea;
            end

        end

        function createGrid(obj,...
                                      minCoordinates,...
                                      maxCoordinates,...
                                      numDivisions)
        % creates a 3d cartesian grid 
        %------------------------------------------------------------------
        % Inputs:
        %   minCoordinates -- [minx,miny,minz]
        %   maxCoordinates -- [maxx,maxy,maxz]
        %   numDivisions ---- number of cells in x,y,z
        %------------------------------------------------------------------
        
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
            
            obj.coordinates = [obj.coordinates;X(:),Y(:),Z(:)];
            obj.numVertices = size(obj.coordinates,1);
                              
        end
        function delaunayTriangulation(obj)
        % creates DT from the stored vertices
        %------------------------------------------------------------------
            DT = delaunayTriangulation(obj.coordinates);
            obj.cells = DT.ConnectivityList;
            obj.numCells = size(obj.cells,1);
            obj.numVertices = size(obj.coordinates,1);
            obj.numNodes = obj.numVertices;
                              
        end
        function clipExternalCells(obj)
        % if cell centroid is external to surface def cull it
        %------------------------------------------------------------------
            if ~ obj.hasSurfaceMesh
                error("requires surface definition to clip")
            end
            
            centroids = obj.cellCentroids();
            isInside = obj.surfaceMesh.isInside(centroids);
            obj.cells = obj.cells(isInside>0.5,:);
            
            obj.numCells = size(obj.cells,1);
            obj.numVertices = size(obj.coordinates,1);
            obj.numNodes = obj.numVertices;

        end
        function clipInternalCells(obj)
        % if cell centroid is internal to surface def cull it
        %------------------------------------------------------------------
               
            if ~ obj.hasSurfaceMesh
                error("requires surface definition to clip")
            end
            
            centroids = obj.cellCentroids;
            isInside = obj.surfaceMesh.isInside(centroids);
            obj.cells = obj.cells(isInside<0.5,:);
            
            obj.numCells = size(obj.cells,1);
            obj.numVertices = size(obj.coordinates,1);
            obj.numNodes = obj.numVertices;
            
        end
        
        function setDegree(obj,degree)
        % set the degree of the mesh (1 or 2)
        %------------------------------------------------------------------
        % for degree 2 the cell ordering is 
        % [v1, e12, v2, e13, e23, v3, e14, e24, e34, v4]
        % w/ v being vertex and e edges
        %------------------------------------------------------------------
            if degree ~= obj.degree
                switch degree
                    case 1
                        obj.flatten();
                    case 2
                        [edges, ~, ~, cellEdges] = obj.edgesAndRelations();
                        size(obj.coordinates);

                        % add midpoints to our coords
                        obj.coordinates = [obj.coordinates;...
                            0.5*(obj.coordinates(edges(:,1),:) +...
                            obj.coordinates(edges(:,2),:))];
                        size(obj.coordinates);

                        % midpoints on boundary
                        isBoundaryEdge = 0.5*(obj.isBoundaryNode(edges(:,1),:) +...
                            obj.isBoundaryNode(edges(:,2),:));
                        isBoundaryEdge(isBoundaryEdge<0.75) = 0;
                        obj.isBoundaryNode=[obj.isBoundaryNode;...
                            isBoundaryEdge];

                        % increase index so cells point to right location
                        cellEdges = cellEdges + obj.numVertices;

                        % set up our 10 nodes in correct local orientation
                        obj.cells = [obj.cells(:,1),...
                            cellEdges(:,1),...
                            obj.cells(:,2),...
                            cellEdges(:,2),...
                            cellEdges(:,3),...
                            obj.cells(:,3),...
                            cellEdges(:,4),...
                            cellEdges(:,5),...
                            cellEdges(:,6),...
                            obj.cells(:,4)];

                        obj.numBoundaryNodes = sum(obj.isBoundaryNode);
                        obj.numNodes = size(obj.coordinates,1);
                        obj.isCurved=true;
                        obj.degree = 2;
                        obj.initializeBulkProperties();
                    otherwise
                        error('specified degree not supported')
                end
            end
        end
        function flatten(obj)
        % force degree to 1
        %------------------------------------------------------------------
            if obj.degree==2
                obj.coordinates = obj.coordinates(1:obj.numVertices,:);
                obj.isBoundaryNode = obj.isBoundaryNode(1:obj.numVertices,:);
                obj.cells = [obj.cells(:,1),...
                             obj.cells(:,3),...
                             obj.cells(:,6),...
                             obj.cells(:,10)];

                obj.isCurved=false;
                obj.degree = 1;
                obj.initializeBulkProperties();
                obj.numBoundaryNodes = sum(obj.isBoundaryNode);
                obj.numNodes = size(obj.coordinates,1);
            end
        end
        function smooth(obj)
        % simple smoothing with uniform weights
        %------------------------------------------------------------------

            delta = zeros(obj.numVertices,3);
            wsum = zeros(obj.numVertices,1);
            for i = 1:obj.numCells
                vertices = obj.cells(i,:);
                coords = obj.coordinates(vertices,:);
                v1 = coords(2,:)-coords(1,:);
                v2 = coords(3,:)-coords(1,:);
                v3 = coords(4,:)-coords(1,:);
                v4 = coords(3,:)-coords(2,:);
                v5 = coords(4,:)-coords(2,:);
                v6 = coords(4,:)-coords(3,:);

                delta(vertices(1),:) = delta(vertices(1),:) +v1+v2+v3;
                delta(vertices(2),:) = delta(vertices(2),:) -v1+v4+v5;
                delta(vertices(3),:) = delta(vertices(3),:) -v2-v4+v6;
                delta(vertices(4),:) = delta(vertices(4),:) -v3-v5-v6;
                
                wsum(vertices(1)) = wsum(vertices(1))+3;
                wsum(vertices(2)) = wsum(vertices(2))+3;
                wsum(vertices(3)) = wsum(vertices(3))+3;
                wsum(vertices(4)) = wsum(vertices(4))+3;


            end
            delta(obj.isBoundaryNode==1,:) = 0;
            obj.coordinates = obj.coordinates + delta./wsum;
        end
        function edges = edges(obj)
        % unique edges within mesh
        %------------------------------------------------------------------
        % we make a map using a cell array. the min vertex index in a
        % vertex pair is used as the cell array index and the max
        % vertex index stored if it is unique.
        %------------------------------------------------------------------
        % Outputs:
        %   edges ----- vertex ids bounding each unique edge
        %------------------------------------------------------------------
            
            % initialize our maps
            cellEdges{obj.numVertices} = [];
            edges = zeros(6*obj.numVertices,2);
            numUniqueEdges=0;
            
            % local vertex indices for 6 edges
            nodei = [1,2,3,4,4,1];
            nodej = [2,3,4,1,2,3];
            
            % for each cell loop the local edges
            for i = 1:obj.numCells
                
                for j = 1:6
                    
                    vertex1 = obj.cells(i,nodei(j));
                    vertex2 = obj.cells(i,nodej(j));
                    vmin = min(vertex1,vertex2);
                    vmax = max(vertex1,vertex2);
                    
                    % if we haven't seen this unique edge before add it
                    if ~any(cellEdges{vmin}==vmax)
                        
                        numUniqueEdges=numUniqueEdges+1;
                        
                        cellEdges{vmin}=[cellEdges{vmin},vmax];
                        edges(numUniqueEdges,1:2) = [vmin,vmax];
                    end
                end
            end
            
            % clip if we over alloted
            edges = edges(1:numUniqueEdges,1:2);
            
        end
        function [edges, neighborVertices, edgeCells, cellEdges] = edgesAndRelations(obj)
        % unique edges within mesh and cell associativity
        %------------------------------------------------------------------
        % we make a map using a cell array. the min vertex index in a
        % vertex pair is used as the cell array index and the max
        % vertex index stored if it is unique.
        %------------------------------------------------------------------
        % Outputs:
        %   edges ----- vertex ids bounding each unique edge
        %   cellEdges - edges assoc w/ each cell
        %   edgeCells - cells assoc w/ each ege
        %------------------------------------------------------------------
            
            % initialize our maps
            neighborVertices{obj.numVertices} = [];
            edgeCells{obj.numVertices} = [];
            cellEdges = zeros(obj.numCells,6);
            uniqueEdgeIndices{obj.numVertices} = [];
            edges = zeros(6*obj.numVertices,2);
            numUniqueEdges=0;
            
            % local vertex indices for 6 edges
            nodei = [1,2,3,4,4,1];
            nodej = [2,3,4,1,2,3];
            
            % map  local vertex (index1,index2) -> (local edge index)
            vertexEdgeMap = [0, 1, 2, 4;...
                             1, 0, 3, 5;...
                             2, 3, 0, 6;...
                             4, 5, 6, 0];

            % for each cell loop the local edges
            for i = 1:obj.numCells
                for j = 1:6
                    
                    % local indices of edge nodes
                    p1 = nodei(j);
                    p2 = nodej(j);

                    localEdgeIndex = vertexEdgeMap(p1,p2);

                    % order our vertices by their index
                    vertex1 = obj.cells(i,p1);
                    vertex2 = obj.cells(i,p2);
                    vmin = min(vertex1,vertex2);
                    vmax = max(vertex1,vertex2);
                    
                    % if we haven't seen this unique edge before add it
                    if ~any(neighborVertices{vmin}==vmax)
                        
                        numUniqueEdges=numUniqueEdges+1;
                        
                        % neighbor vertices to the min index vertex of edge
                        neighborVertices{vmin}=[neighborVertices{vmin},vmax];

                        % unique edges attached to vertex vmin
                        uniqueEdgeIndices{vmin} = [uniqueEdgeIndices{vmin},numUniqueEdges];

                        % nodes of edges
                        edges(numUniqueEdges,1:2) = [vmin,vmax];
                    end

                    % unique edge
                    uniqueEdgeij = uniqueEdgeIndices{vmin}(neighborVertices{vmin}==vmax);

                    % map cell - > edges w/ appropriate ordering
                    cellEdges(i,localEdgeIndex) = uniqueEdgeij;
                    
                    % map edges - > cells
                    if uniqueEdgeij > length(edgeCells)
                        edgeCells{uniqueEdgeij}=[i];
                    else
                        edgeCells{uniqueEdgeij} = [edgeCells{uniqueEdgeij},i];
                    end
                    
                    
                end
            end
            
            % clip if we over alloted
            edges = edges(1:numUniqueEdges,1:2);
            
        end
        
        function [vertices] = cellVertices(obj)
        % gets the vertices for the given cell
        %------------------------------------------------------------------
        % with degree > 1 cells, the vertices are interwoven with the other 
        % geometric nodes. This function with separate them out.
        %------------------------------------------------------------------

            switch obj.degree
                case 1
                    vertices = obj.cells;
                case 2
                    vertices = [obj.cells(:,1),...
                                obj.cells(:,3),...
                                obj.cells(:,6),...
                                obj.cells(:,10)];
            end

        end

        function [centroids,volumes] = cellCentroids(obj)
        % centroids of cells
        %------------------------------------------------------------------
        % Outputs:
        %   centroids - coordinates of centroids
        %------------------------------------------------------------------
            if ~obj.isCurved
                p1 = obj.coordinates(obj.cells(:,1),:);
                p2 = obj.coordinates(obj.cells(:,2),:); 
                p3 = obj.coordinates(obj.cells(:,3),:); 
                p4 = obj.coordinates(obj.cells(:,4),:); 
                centroids=1/4*(p1+p2+p3+p4); 
            else
                % high degree quadrature rule (Nix1)
                quadratureOrder = 6;
                [u,v,w,weights] =  NewtonCotesTetrahedron(quadratureOrder);

                % Ng basis functions assessed at the Ni quad points (NixNg)
                [phi,dphidu,dphidv,dphidw] =  LagrangeInterpolantsTetrahedron( u, v, w, obj.degree);

                % x,y,z coords of interpolation points for cells (Nc x Ng)
                numNodesPerCell = size(obj.cells,2);
                Xx = zeros(obj.numCells,numNodesPerCell);
                Xy = zeros(obj.numCells,numNodesPerCell);
                Xz = zeros(obj.numCells,numNodesPerCell);
                
                for i = 1:numNodesPerCell
                    Xx(:,i) = obj.coordinates(obj.cells(:,i),1);
                    Xy(:,i) = obj.coordinates(obj.cells(:,i),2);
                    Xz(:,i) = obj.coordinates(obj.cells(:,i),3);
                end

                % coordinates of all quadrature points ( NixNc )
                x = phi*Xx';
                y = phi*Xy';
                z = phi*Xz';

                % components of jacobian at all quad points (NixNc)
                dXdx11 = dphidu*Xx';
                dXdx12 = dphidv*Xx';
                dXdx13 = dphidw*Xx';
                dXdx21 = dphidu*Xy';
                dXdx22 = dphidv*Xy';
                dXdx23 = dphidw*Xy';
                dXdx31 = dphidu*Xz';
                dXdx32 = dphidv*Xz';
                dXdx33 = dphidw*Xz';

                % det of jacobian at all quad points ( NixNc )
                J = dXdx11.*(dXdx22.*dXdx33 - dXdx23.*dXdx32) - ...
                    dXdx12.*(dXdx21.*dXdx33 - dXdx23.*dXdx31) + ...
                    dXdx13.*(dXdx21.*dXdx32 - dXdx22.*dXdx31);

                J = abs(J);

                % position times det(J) (NixNc) -- for centroid calc
                xJ = x.*J;
                yJ = y.*J;
                zJ = z.*J;

                volumes =  J'*weights;
                centroids =  [xJ'*weights,yJ'*weights,zJ'*weights];

                centroids = centroids./volumes;
                volumes = volumes/6.0;
            end

            
        end
        function [centroids,volumes] = vertexCentroids(obj)
        % centroids & volumes of nodes
        %------------------------------------------------------------------
        % after vectorizing the general part got a little hard to read.
        % Dimensions of matrices are listed with the following definitions:
        % Ng -- number of local geometric quadrature points defining cell
        % Ni -- number of high res quadrature points used to integrate over
        %       the cell
        % Nc -- number of cells
        %------------------------------------------------------------------
        % Outputs:
        %   centroids - coordinates of node centroids
        %   volumes --- volumes of nodes
        %------------------------------------------------------------------
            
            centroids = zeros(obj.numVertices,3);
            volumes = zeros(obj.numVertices,1);

            if ~ obj.isCurved

                cellCentroids = obj.cellCentroids();
                cellVolumes = obj.cellVolumes();

                for i = 1:obj.numCells
                    ids = obj.cells(i,:);
                    ci = cellVolumes(i).*cellCentroids(i,:);

                    volumes(ids) = volumes(ids)+cellVolumes(i);
                    centroids(ids,:) = centroids(ids,:) + [ci;ci;ci;ci];
                end
                centroids = centroids./volumes;
                volumes = volumes/4;

            else
                
                % high degree quadrature rule (Nix1)
                quadratureOrder = 6;
                [u,v,w,weights] =  NewtonCotesTetrahedron(quadratureOrder);

                % Ng basis functions assessed at the Ni quad points (NixNg)
                [phi,dphidu,dphidv,dphidw] =  LagrangeInterpolantsTetrahedron( u, v, w, obj.degree);

                % vertex basis functions (NixNg)
                [phiq,~,~,~] =  LagrangeInterpolantsTetrahedron( u, v, w, 1);
                rectilinearCells = obj.cellVertices();

                % x,y,z coords of interpolation points for cells (Nc x Ng)
                numNodesPerCell = size(obj.cells,2);
                Xx = zeros(obj.numCells,numNodesPerCell);
                Xy = zeros(obj.numCells,numNodesPerCell);
                Xz = zeros(obj.numCells,numNodesPerCell);
                
                for i = 1:numNodesPerCell
                    Xx(:,i) = obj.coordinates(obj.cells(:,i),1);
                    Xy(:,i) = obj.coordinates(obj.cells(:,i),2);
                    Xz(:,i) = obj.coordinates(obj.cells(:,i),3);
                end

                % coordinates of all quadrature points ( NixNc )
                x = phi*Xx';
                y = phi*Xy';
                z = phi*Xz';

                % components of jacobian at all quad points (NixNc)
                dXdx11 = dphidu*Xx';
                dXdx12 = dphidv*Xx';
                dXdx13 = dphidw*Xx';
                dXdx21 = dphidu*Xy';
                dXdx22 = dphidv*Xy';
                dXdx23 = dphidw*Xy';
                dXdx31 = dphidu*Xz';
                dXdx32 = dphidv*Xz';
                dXdx33 = dphidw*Xz';

                % det of jacobian at all quad points ( NixNc )
                J = dXdx11.*(dXdx22.*dXdx33 - dXdx23.*dXdx32) - ...
                    dXdx12.*(dXdx21.*dXdx33 - dXdx23.*dXdx31) + ...
                    dXdx13.*(dXdx21.*dXdx32 - dXdx22.*dXdx31);

                J = abs(J);

                % position times det(J) (NixNc) -- for centroid calc
                xJ = x.*J;
                yJ = y.*J;
                zJ = z.*J;

                % (Nix1)' - ( (NixNg) .* Nix1)
                for i = 1:obj.numCells
                    vertexIndices = rectilinearCells(i,:); 
                    volumes(vertexIndices) = volumes(vertexIndices,1) + (phiq.*J(:,i))'*weights;
                    centroids(vertexIndices,1) = centroids(vertexIndices,1) + (phiq.*xJ(:,i))'*weights;
                    centroids(vertexIndices,2) = centroids(vertexIndices,2) + (phiq.*yJ(:,i))'*weights;
                    centroids(vertexIndices,3) = centroids(vertexIndices,3) + (phiq.*zJ(:,i))'*weights;

                end
                centroids = centroids./volumes;
                volumes = volumes/6.0;
            end
        end
        function [centroids,volumes] = nodeCentroids(obj)
        % centroids & volumes of nodes
        %------------------------------------------------------------------
        % after vectorizing the general part got a little hard to read.
        % Dimensions of matrices are listed with the following definitions:
        % Ng -- number of local geometric quadrature points defining cell
        % Ni -- number of high res quadrature points used to integrate over
        %       the cell
        % Nc -- number of cells
        %------------------------------------------------------------------
        % Outputs:
        %   centroids - coordinates of node centroids
        %   volumes --- volumes of nodes
        %------------------------------------------------------------------
            
            centroids = zeros(obj.numNodes,3);
            volumes = zeros(obj.numNodes,1);

            if ~ obj.isCurved

                cellCentroids = obj.cellCentroids();
                cellVolumes = obj.cellVolumes();

                for i = 1:obj.numCells
                    ids = obj.cells(i,:);
                    ci = cellVolumes(i).*cellCentroids(i,:);

                    volumes(ids) = volumes(ids)+cellVolumes(i);
                    centroids(ids,:) = centroids(ids,:) + [ci;ci;ci;ci];
                end
                centroids = centroids./volumes;
                volumes = volumes/4;

            else
                
                % high degree quadrature rule (Nix1)
                quadratureOrder = 6;
                [u,v,w,weights] =  NewtonCotesTetrahedron(quadratureOrder);

                % Ng basis functions assessed at the Ni quad points (NixNg)
                [phi,dphidu,dphidv,dphidw] =  LagrangeInterpolantsTetrahedron( u, v, w, obj.degree);

                % x,y,z coords of interpolation points for cells (Nc x Ng)
                numNodesPerCell = size(obj.cells,2);
                Xx = zeros(obj.numCells,numNodesPerCell);
                Xy = zeros(obj.numCells,numNodesPerCell);
                Xz = zeros(obj.numCells,numNodesPerCell);
                
                for i = 1:numNodesPerCell
                    Xx(:,i) = obj.coordinates(obj.cells(:,i),1);
                    Xy(:,i) = obj.coordinates(obj.cells(:,i),2);
                    Xz(:,i) = obj.coordinates(obj.cells(:,i),3);
                end

                % coordinates of all quadrature points ( NixNc )
                x = phi*Xx';
                y = phi*Xy';
                z = phi*Xz';

                % components of jacobian at all quad points (NixNc)
                dXdx11 = dphidu*Xx';
                dXdx12 = dphidv*Xx';
                dXdx13 = dphidw*Xx';
                dXdx21 = dphidu*Xy';
                dXdx22 = dphidv*Xy';
                dXdx23 = dphidw*Xy';
                dXdx31 = dphidu*Xz';
                dXdx32 = dphidv*Xz';
                dXdx33 = dphidw*Xz';

                % det of jacobian at all quad points ( NixNc )
                J = dXdx11.*(dXdx22.*dXdx33 - dXdx23.*dXdx32) - ...
                    dXdx12.*(dXdx21.*dXdx33 - dXdx23.*dXdx31) + ...
                    dXdx13.*(dXdx21.*dXdx32 - dXdx22.*dXdx31);

                J = abs(J);

                % position times det(J) (NixNc) -- for centroid calc
                xJ = x.*J;
                yJ = y.*J;
                zJ = z.*J;

                % (Nix1)' - ( (NixNg) .* Nix1)
                for i = 1:obj.numCells
                    nodeIndices = obj.cells(i,:); 
                    volumes(nodeIndices) = volumes(nodeIndices,1) + (phi.*J(:,i))'*weights;
                    centroids(nodeIndices,1) = centroids(nodeIndices,1) + (phi.*xJ(:,i))'*weights;
                    centroids(nodeIndices,2) = centroids(nodeIndices,2) + (phi.*yJ(:,i))'*weights;
                    centroids(nodeIndices,3) = centroids(nodeIndices,3) + (phi.*zJ(:,i))'*weights;

                end
                centroids = centroids./volumes;
                volumes = volumes/6.0;
            end
        end
        function volumes = cellVolumes(obj)
        % centroids of cells
        %------------------------------------------------------------------
        % Outputs:
        %   centroids - coordinates of centroids
        %------------------------------------------------------------------
            
            if ~obj.isCurved

                p1 = obj.coordinates(obj.cells(:,1),:);
                p2 = obj.coordinates(obj.cells(:,2),:)-p1; 
                p3 = obj.coordinates(obj.cells(:,3),:)-p1; 
                p4 = obj.coordinates(obj.cells(:,4),:)-p1; 
                volumes = 1/6*abs(dot(p4,cross(p2,p3,2),2));

            else
                    
                % high degree quadrature rule (Nix1)
                quadratureOrder = obj.degree;
                [u,v,w,weights] =  NewtonCotesTetrahedron(quadratureOrder);

                % Ng basis functions assessed at the Ni quad points (NixNg)
                [~,dphidu,dphidv,dphidw] =  LagrangeInterpolantsTetrahedron( u, v, w, obj.degree);

                % x,y,z coords of interpolation points for cells (Nc x Ng)
                numNodesPerCell = size(obj.cells,2);
                Xx = zeros(obj.numCells,numNodesPerCell);
                Xy = zeros(obj.numCells,numNodesPerCell);
                Xz = zeros(obj.numCells,numNodesPerCell);
                
                for i = 1:numNodesPerCell
                    Xx(:,i) = obj.coordinates(obj.cells(:,i),1);
                    Xy(:,i) = obj.coordinates(obj.cells(:,i),2);
                    Xz(:,i) = obj.coordinates(obj.cells(:,i),3);
                end

                % components of jacobian at all quad points (NixNc)
                dXdx11 = dphidu*Xx';
                dXdx12 = dphidv*Xx';
                dXdx13 = dphidw*Xx';
                dXdx21 = dphidu*Xy';
                dXdx22 = dphidv*Xy';
                dXdx23 = dphidw*Xy';
                dXdx31 = dphidu*Xz';
                dXdx32 = dphidv*Xz';
                dXdx33 = dphidw*Xz';

                % det of jacobian at all quad points ( NixNc )
                J = dXdx11.*(dXdx22.*dXdx33 - dXdx23.*dXdx32) - ...
                    dXdx12.*(dXdx21.*dXdx33 - dXdx23.*dXdx31) + ...
                    dXdx13.*(dXdx21.*dXdx32 - dXdx22.*dXdx31);

                J = abs(J);

                volumes = J'*weights/6.0;
            end
            
        end
        function volumes = vertexVolumes(obj)
        % calculate volume assoc. w/ each node
        %------------------------------------------------------------------
            
            volumes = zeros(obj.numVertices,1);
        
            if ~obj.isCurved
                cellVolumes = obj.cellVolumes();
                for i = 1:obj.numCells
                    ids = obj.cells(i,:);
                    volumes(ids) = volumes(ids)+0.25*cellVolumes(i);
                end
            else
                
                % high degree quadrature rule (Nix1)
                quadratureOrder = obj.degree+2;
                [u,v,w,weights] =  NewtonCotesTetrahedron(quadratureOrder);

                % Ng basis functions assessed at the Ni quad points (NixNg)
                [~,dphidu,dphidv,dphidw] =  LagrangeInterpolantsTetrahedron( u, v, w, obj.degree);

                % vertex basis functions
                [phiq,~,~,~] =  LagrangeInterpolantsTetrahedron( u, v, w, 1);
                rectilinearCells = obj.cellVertices();

                % x,y,z coords of interpolation points for cells (Nc x Ng)
                numNodesPerCell = size(obj.cells,2);
                Xx = zeros(obj.numCells,numNodesPerCell);
                Xy = zeros(obj.numCells,numNodesPerCell);
                Xz = zeros(obj.numCells,numNodesPerCell);
                
                for i = 1:numNodesPerCell
                    Xx(:,i) = obj.coordinates(obj.cells(:,i),1);
                    Xy(:,i) = obj.coordinates(obj.cells(:,i),2);
                    Xz(:,i) = obj.coordinates(obj.cells(:,i),3);
                end

                % components of jacobian at all quad points (NixNc)
                dXdx11 = dphidu*Xx';
                dXdx12 = dphidv*Xx';
                dXdx13 = dphidw*Xx';
                dXdx21 = dphidu*Xy';
                dXdx22 = dphidv*Xy';
                dXdx23 = dphidw*Xy';
                dXdx31 = dphidu*Xz';
                dXdx32 = dphidv*Xz';
                dXdx33 = dphidw*Xz';

                % det of jacobian at all quad points ( NixNc )
                J = dXdx11.*(dXdx22.*dXdx33 - dXdx23.*dXdx32) - ...
                    dXdx12.*(dXdx21.*dXdx33 - dXdx23.*dXdx31) + ...
                    dXdx13.*(dXdx21.*dXdx32 - dXdx22.*dXdx31);

                J = abs(J);

                % (Nix1)' - ( (NixNg) .* Nix1)
                for i = 1:obj.numCells
                    vertexIndices = rectilinearCells(i,:); 
                    volumes(vertexIndices) = volumes(vertexIndices,1) + (phiq.*J(:,i))'*weights;
                end
                volumes = volumes/6.0;
            end
        end
        function volumes = nodeVolumes(obj)
        % calculate volume assoc. w/ each node
        %------------------------------------------------------------------
            
            volumes = zeros(obj.numNodes,1);
        
            if ~obj.isCurved
                cellVolumes = obj.cellVolumes();
                for i = 1:obj.numCells
                    ids = obj.cells(i,:);
                    volumes(ids) = volumes(ids)+0.25*cellVolumes(i);
                end
            else
                
                % high degree quadrature rule (Nix1)
                quadratureOrder = 6;
                [u,v,w,weights] =  NewtonCotesTetrahedron(quadratureOrder);

                % Ng basis functions assessed at the Ni quad points (NixNg)
                [phi,dphidu,dphidv,dphidw] =  LagrangeInterpolantsTetrahedron( u, v, w, obj.degree);

                % x,y,z coords of interpolation points for cells (Nc x Ng)
                numNodesPerCell = size(obj.cells,2);
                Xx = zeros(obj.numCells,numNodesPerCell);
                Xy = zeros(obj.numCells,numNodesPerCell);
                Xz = zeros(obj.numCells,numNodesPerCell);
                
                for i = 1:numNodesPerCell
                    Xx(:,i) = obj.coordinates(obj.cells(:,i),1);
                    Xy(:,i) = obj.coordinates(obj.cells(:,i),2);
                    Xz(:,i) = obj.coordinates(obj.cells(:,i),3);
                end

                % components of jacobian at all quad points (NixNc)
                dXdx11 = dphidu*Xx';
                dXdx12 = dphidv*Xx';
                dXdx13 = dphidw*Xx';
                dXdx21 = dphidu*Xy';
                dXdx22 = dphidv*Xy';
                dXdx23 = dphidw*Xy';
                dXdx31 = dphidu*Xz';
                dXdx32 = dphidv*Xz';
                dXdx33 = dphidw*Xz';

                % det of jacobian at all quad points ( NixNc )
                J = dXdx11.*(dXdx22.*dXdx33 - dXdx23.*dXdx32) - ...
                    dXdx12.*(dXdx21.*dXdx33 - dXdx23.*dXdx31) + ...
                    dXdx13.*(dXdx21.*dXdx32 - dXdx22.*dXdx31);

                J = abs(J);

                % (Nix1)' - ( (NixNg) .* Nix1)
                for i = 1:obj.numCells
                    nodeIndices = obj.cells(i,:); 
                    volumes(nodeIndices) = volumes(nodeIndices,1) + (phi.*J(:,i))'*weights;
                end
                volumes = volumes/6.0;
            end
        end
        
        function volume = calculateVolume(obj)
        % volume contained in mesh
        %------------------------------------------------------------------
        % generalized volume calc for curvilinear meshes. 
        %------------------------------------------------------------------
        % Outputs:
        %   vol - volume
        %------------------------------------------------------------------
            volume = sum(obj.nodeVolumes());
        end
        function centroid = calculateCentroid(obj)
        % centroid of body
        %------------------------------------------------------------------
        % Outputs:
        %   c - 1x3 centroid
        %------------------------------------------------------------------
            [c,v] = obj.nodeCentroids;
            centroid = sum(v.*c,1)/sum(v); 
        end
        function resolution = calculateResolution(obj)
            V = obj.cellVolumes();
            resolution =(sum(V)/size(V,1))^(1/3);
        end
        
        function addCellField(obj,fieldData,fieldName)
        % adds a face field that can be export to vtk
        %------------------------------------------------------------------
        % Inputs:
        %   fieldi ---- field data or function handle to calc field data
        %   fieldname - string giving it a name
        %-----------------------------------------------------------------
            if isa(fieldData,'function_handle')
                fieldName = func2str(fieldData);
                fieldData = fieldData(obj.faceCentroids());
            end
            if size(fieldData,1)~=obj.numCells
                error('faceFieldData must have num of row equal to numFaces')
            end
            obj.numCellFields = obj.numCellFields + 1;
            if nargin == 2
                fieldName = ['field_',num2str(obj.numCellFields)];
            elseif nargin ~=3
                error('incorrect number of inputs')
            end
            obj.cellFields{obj.numCellFields}.data = fieldData;
            obj.cellFields{obj.numCellFields}.name = fieldName;
        end
        function addNodeField(obj,fieldData,fieldName)
        % adds a node field that can be export to vtk
        %------------------------------------------------------------------
        % Inputs:
        %   fieldi ---- field data or function handle to calc field data
        %   fieldname - string giving it a name
        %-----------------------------------------------------------------
            if isa(fieldData,'function_handle')
                fieldName = func2str(fieldData);
                fieldData = fieldData(obj.coordinates);
            end
            if size(fieldData,1)~=obj.numNodes
                error('nodeFieldData must have num of row equal to numNodes')
            end
            obj.numNodeFields = obj.numNodeFields + 1;
            if nargin == 2
                fieldName = ['field_',num2str(obj.numNodeFields)];
            elseif nargin ~=3
                error('incorrect number of inputs')
            end
            obj.nodeFields{obj.numNodeFields}.data = fieldData;
            obj.nodeFields{obj.numNodeFields}.name = fieldName;
        end
        function clearFields(obj)
        % clears out the field data 
        %------------------------------------------------------------------
        
            obj.numNodeFields = 0;
            obj.numCellFields = 0;
            
            obj.nodeFields = [];
            obj.cellFields = [];
        end
        
        function cellValues = cellValuesFromVertices(obj,vertexValues)
        % get cell value from vertices
        %------------------------------------------------------------------
            rectilinearCells = obj.cellVertices();

            p1 = rectilinearCells(:,1);
            p2 = rectilinearCells(:,2);
            p3 = rectilinearCells(:,3);
            p4 = rectilinearCells(:,4);
           
            cellValues = 1/4*(vertexValues(p1,:)+...
                              vertexValues(p2,:)+...
                              vertexValues(p3,:)+...
                              vertexValues(p4,:)); 

        end
        function nodeValues = nodeValuesFromCellValues(obj,cellValues)
        % gets node values from vol weight cell values
        %------------------------------------------------------------------
            nodeValues = zeros(obj.numNodes,size(cellValues,2));
            nodeVolumes = zeros(obj.numNodes,1);

            cellVolumes = obj.cellVolumes;
            for i = 1:obj.numCells
                ids = obj.cells(i,:);
                nodeVolumes(ids,:) = nodeVolumes(ids,:) + cellVolumes(i);
                nodeValues(ids,:) = nodeValues(ids,:) + cellVolumes(i) * cellValues(i,:);
            end
            nodeValues = nodeValues./nodeVolumes;
        end
        function vertexValues = vertexValuesFromCellValues(obj,cellValues)
        % gets vertex values from vol weight cell values
        %------------------------------------------------------------------
            vertexValues = zeros(obj.numNodes,size(cellValues,2));
            nodeVolumes = zeros(obj.numNodes,1);

            rectilinearCells = obj.cellVertices();
            cellVolumes = obj.cellVolumes;
            for i = 1:obj.numCells
                ids = rectilinearCells(i,:);
                nodeVolumes(ids,:) = nodeVolumes(ids,:) + cellVolumes(i);
                vertexValues(ids,:) = vertexValues(ids,:) + cellVolumes(i) * cellValues(i,:);
            end
            vertexValues = vertexValues./nodeVolumes;
        end
        
        function coords = boundaryCoordinates(obj)
        % get the unit normals of the boundary nodes
        %------------------------------------------------------------------
            coords = obj.coordinates(obj.isBoundaryNode==1,:);
        end
        function normals = boundaryNodeNormals(obj)
        % get the unit normals of the boundary nodes
        %------------------------------------------------------------------
            normals = zeros(obj.numVertices,3);
            normals(1:obj.surfaceMesh.numVertices,:) = obj.surfaceMesh.vertexNormals();

            % were going map to the cell and back to the nodes to get
            % normals for our other boundary nodes and smooth things out
            normals = obj.cellValuesFromVertices(normals);
            normals = obj.nodeValuesFromCellValues(normals);
            normals = normals(obj.isBoundaryNode==1,:);
            normals = normals./sqrt(normals(:,1).^2 + ...
                                    normals(:,2).^2 + ...
                                    normals(:,3).^2);
        end
        function curve(obj,surfaceMesh)
        % curve 
        %------------------------------------------------------------------
            if isa(surfaceMesh,"SurfaceMesh")
                boundaryNodeCoords = obj.boundaryCoordinates();
                boundaryNormals = obj.boundaryNodeNormals();
                newCoordinates = surfaceMesh.project(boundaryNodeCoords,boundaryNormals);
                obj.coordinates(obj.isBoundaryNode==1,:) = newCoordinates;
                
            else
                error('requires surface mesh input')
            end
        end
        
        function writeVTK(obj,fileName,precision)
        % writes mesh to vtk file for visualization w/ paraview etc
        %------------------------------------------------------------------
        % Inputs:
        %   fileName ----- name of obj file (including extension)
        %   percision ---- number of sigfigs for coordinate ouput
        %==================================================================
            
            % handle optional inputs
            if nargin<2
                fileName = strcat('UnnamedAsteroid_',num2str(obj.numFaces),'.vtk');
            end
            if nargin <3
                precision = 8;
            end
            
            %vtk flags for rectilinear and quadratic tetrahedra -- need to
            %figure out how to get multiple cell types
            if obj.degree == 1
                cellType = '10\n';
                numLocalNodes = 4;
            elseif obj.degree == 2 
                cellType = '24\n';
                numLocalNodes = 10;
            else 
                error('Degree greater than 2 not supported')
            end
            
            fileID = fopen(fileName,'w');
            
            % Header
            %--------------------------------------------------------------------------
            fprintf(fileID,'# vtk DataFile Version 3.0\n');
            fprintf(fileID,'vtk output\n');
            fprintf(fileID,'ASCII\n');
            fprintf(fileID,'DATASET UNSTRUCTURED_GRID\n');
            
            % Write Node Data
            %--------------------------------------------------------------
            fprintf(fileID,['POINTS ',num2str(obj.numNodes),' float\n']);
            for i = 1:obj.numNodes
                outputi = [num2str(obj.coordinates(i,:),precision),'\n'];
                fprintf(fileID, outputi );
            end
            
            % Write Cells --- needs fixing for tetra10
            %--------------------------------------------------------------
            fprintf(fileID,['CELLS ',num2str(obj.numCells),' ',num2str(obj.numCells*(numLocalNodes+1)),'\n']);
            for i = 1:obj.numCells
                if obj.degree == 1
                    outputi = [num2str([length(obj.cells(i,:)),obj.cells(i,:)-1]),'\n'];
                elseif obj.degree == 2
                    teti = obj.cells(i,[1,3,6,10,2,5,4,7,8,9])-1;
                    outputi = [num2str([length(obj.cells(i,:)),teti]),'\n'];
                end
                fprintf(fileID,outputi);
            end
            
            % Write CellTypes
            %--------------------------------------------------------------
            fprintf(fileID,['CELL_TYPES ',num2str(obj.numCells),'\n']);
            for i = 1:obj.numCells
                fprintf(fileID,cellType);
                
            end
            
            % Write Cell Data
            %--------------------------------------------------------------
            fprintf(fileID,['CELL_DATA ',num2str(obj.numCells),'\n']);
            for i = 1:size(obj.cellFields,2)
                if size(obj.cellFields{i}.data,2)==1
                    vtkDataHeader = ['SCALARS ',obj.cellFields{i}.name,' float 1\nLOOKUP_TABLE default\n'];
                elseif size(obj.cellFields{i}.data,2)==3
                    vtkDataHeader = ['VECTORS ',obj.cellFields{i}.name,' float\n'];
                else
                    warning(['skipping faceData{',num2str(i),'} not scalar or vector'])
                end
                fprintf(fileID,vtkDataHeader);
                if size(obj.cellFields{i}.data,1)~=obj.numCells
                    error(['faceData{',num2str(i),'} inconsistent w/ number of mesh faces'])
                end
                for j = 1:obj.numCells
                    outputi = [num2str(obj.cellFields{i}.data(j,:),precision),'\n'];
                    fprintf(fileID,outputi);
                end
                
            end
            
            %Write Nodal Data
            %--------------------------------------------------------------
            fprintf(fileID,['POINT_DATA ',num2str(obj.numNodes),'\n']);
            for i = 1:size(obj.nodeFields,2)
                if size(obj.nodeFields{i}.data,2)==1
                    vtkDataHeader = ['SCALARS ',obj.nodeFields{i}.name,' float 1\nLOOKUP_TABLE default\n'];
                elseif size(obj.nodeFields{i}.data,2)==3
                    vtkDataHeader = ['VECTORS ',obj.nodeFields{i}.name,' float\n'];
                else
                    warning(['skipping faceData{',num2str(i),'} not scalar or vector'])
                end
                fprintf(fileID,vtkDataHeader);
                if size(obj.nodeFields{i}.data,1)~=obj.numNodes
                    error(['faceData{',num2str(i),'} inconsistent w/ number of mesh faces'])
                end
                for j = 1:obj.numNodes
                    outputi = [num2str(obj.nodeFields{i}.data(j,:),precision),'\n'];
                    fprintf(fileID,outputi);
                end
                
            end
            fclose(fileID);
            
        end
        
    end
    
    methods(Access=private)
        
    end
    methods(Static,Access=public)
        
    end
end



