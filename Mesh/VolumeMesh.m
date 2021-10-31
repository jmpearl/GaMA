classdef VolumeMesh < handle
%==========================================================================
% CLEO - Gravity and Mesh Adaptation Libray for Asteroids and Comets
% J.M.Pearl 
% Oct 2021
%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Abreviations:
%   I attempted to stick to a common convention but of course that didn't
%   always work out. In general for variables: 
%
%       he ---- halfedge
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
        isBoundaryNode      % bool array indicating if its a BC face
        
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
        numEdges            % number of edges
        numFaces            % number of faces
        numCells            % number of half edges
        numNodeFields       % number of node fields
        numCellFields       % number of face fields
        
    end
    methods
        function obj = VolumeMesh(surfaceMesh)
            obj.hasSurfaceMesh=false;
            if nargin == 1
                if isa(surfaceMesh,'SurfaceMesh')
                    obj.surfaceMesh=surfaceMesh;
                    obj.hasSurfaceMesh=true;
                end
            end

            obj.numNodeFields = 0;
            obj.numCellFields = 0;
            obj.numCells = 0;
            obj.numVertices = 0;
            obj.numNodes = 0;
            obj.numBoundaryVertices = 0;
            obj.numBoundaryNodes = 0;

            obj.isCurved=false;
            obj.degree = 1;
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
        
        function centroids = cellCentroids(obj)
        % centroids of cells
        %------------------------------------------------------------------
        % Outputs:
        %   centroids - coordinates of centroids
        %------------------------------------------------------------------
            if ~ obj.isCurved
                p1 = obj.coordinates(obj.cells(:,1),:);
                p2 = obj.coordinates(obj.cells(:,2),:); 
                p3 = obj.coordinates(obj.cells(:,3),:); 
                p4 = obj.coordinates(obj.cells(:,4),:); 
                centroids=1/4*(p1+p2+p3+p4); 
            end
            
        end
        function volumes = cellVolumes(obj)
        % centroids of cells
        %------------------------------------------------------------------
        % Outputs:
        %   centroids - coordinates of centroids
        %------------------------------------------------------------------
            if ~ obj.isCurved
                p1 = obj.coordinates(obj.cells(:,1),:);
                p2 = obj.coordinates(obj.cells(:,2),:)-p1; 
                p3 = obj.coordinates(obj.cells(:,3),:)-p1; 
                p4 = obj.coordinates(obj.cells(:,4),:)-p1; 
                volumes = 1/6*dot(p4,cross(p2,p3,2),2);
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
            volume = sum(obj.cellVolumes);
        end
        function centroid = calculateCentroid(obj)
        % centroid of body
        %------------------------------------------------------------------
        % Outputs:
        %   c - 1x3 centroid
        %------------------------------------------------------------------
            vf = obj.cellVolumes;
            cf = obj.cellCentroids;
            centroid = cf/sum(vf); 
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
                    teti = obj.cells(i,[1,3,10,2,4,5,6,7,8,9])-1;
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



