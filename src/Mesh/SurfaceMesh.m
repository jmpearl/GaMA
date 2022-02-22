classdef SurfaceMesh < handle
%==========================================================================
% GaMA - Gravity and Mesh Adaptation 
% J.M.Pearl 
% Feb 2021
%--------------------------------------------------------------------------
% Array-based implementation of half-edge data structure for tri meshes
%--------------------------------------------------------------------------
% half-edges are ordered in sets of three where each set defines a face 
% making the assoc. face index implicit. 
%
% This class support curvilinear surface definitions up to degree 4 in a
% limited capacity. The curvilinear feature is primarily designed for the
% generation of efficient integration rules and would typically be the last
% step in any mesh manipulation process. For example, the coarsen, refine,
% and smooth methods do not support curvilinear surface definitions and
% will either return an error or flatten the mesh. Additionally, methods
% prefixed with face... vertex... or halfEdge... as well as the edge and
% isInsideRigorous methods all treat curvilinear meshes as if they are
% rectilinear -- i.e. they return results as if the faces were flat w/ the
% same vertices.
%--------------------------------------------------------------------------
% Quick note on nomenclature:
%   In this class, node is used to refer to all the unique computational
%   nodes. For example a T6 triangle node would be used to refer to the
%   edge midpoints and vertices.
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
        next              % index of next half-edge in face cycle
        pair              % index of pair half-edge
        ends              % index of vertex half-edge points to
        vertexHalfEdges   % index of half-edge vertex points to
        coordinates       % coordinates of vertices
        faces             % nodes of faces
        
        faceFields        % cell array of stored fields
        nodeFields        % cell array of stored fields 
        
        isCurved          % logical to track if mesh is curvilinear
        degree            % polynomial degree of mesh
        volume            % total volume
        surfaceArea       % total surface area
        centroid          % centroid of enclosed volume
        resolution        % mean resolution length scale of mesh
        
        numNodes          % number of nodes for degree>1 meshes
        numVertices       % number of vertices
        numFixedVertices  % number of vertices that can't be altered
        numEdges          % number of edges
        numHalfEdges      % number of half edges
        numFixedHalfEdges % number of half edges that can't be altered
        numFaces          % number of faces
        numNodeFields     % number of node fields
        numFaceFields     % number of face fields
        
        coarsenOptions    % options for edge collapse algorithm
        refineOptions     % options for refinement algorithm
        smoothOptions     % options for smoothing algorithm

        isFixedVertex     % vertices that won't be modified
        isFixedHalfEdge   % edges that won't be modified
    end

    methods
        function obj = SurfaceMesh(mesh)
        % constructor
        %------------------------------------------------------------------
        % Input:
        %   mesh --- can be either a obj file name or another mesh object
        %------------------------------------------------------------------   
            if nargin == 1
                if isa(mesh,'SurfaceMesh')
                    
                    obj.copy(mesh);
                    
                elseif isstring(mesh) || ischar(mesh)
                    
                    if contains(mesh,'obj')
                        [tempVertices,tempFaces] = obj.readOBJ(mesh);
                        obj.initializeFromFaceData(tempVertices,tempFaces);
                    else
                        error('currently only supports obj files')
                    end
                    
                end
            elseif nargin~=0
                error('incorrect number of inputs')
            end
            
            obj.setDefaultOptions()
            
        end
        function copy(obj,mesh)
        % implemention of deep copy constructor
        %------------------------------------------------------------------
        
            obj.next = mesh.next;
            obj.pair = mesh.pair;
            obj.ends = mesh.ends;
            obj.vertexHalfEdges = mesh.vertexHalfEdges;
            obj.coordinates = mesh.coordinates;
            obj.faces = mesh.faces;      
            
            obj.faceFields = mesh.faceFields;
            obj.nodeFields = mesh.nodeFields;
            
            obj.isCurved = mesh.isCurved;
            obj.degree = mesh.degree;
            
            obj.volume = mesh.volume;
            obj.surfaceArea = mesh.surfaceArea;
            obj.centroid = mesh.centroid;
            obj.resolution = mesh.resolution;
            
            obj.numNodes = mesh.numNodes;
            obj.numVertices = mesh.numVertices;
            obj.numEdges = mesh.numEdges;
            obj.numHalfEdges = mesh.numHalfEdges;
            obj.numFaces  = mesh.numFaces;
            obj.numNodeFields = mesh.numNodeFields;
            obj.numFaceFields = mesh.numFaceFields;

            obj.coarsenOptions = mesh.coarsenOptions;
            obj.refineOptions  = mesh.refineOptions;
            obj.smoothOptions  = mesh.smoothOptions;

            obj.isFixedVertex = mesh.isFixedVertex;
            obj.isFixedHalfEdge = mesh.isFixedHalfEdge;
            
        end
        function initializeFromFaceData(obj,tempVertices,tempFaces)
        % reads in vertex-face representation and initialized he mesh
        %------------------------------------------------------------------
        % Inputs:
        %   vertices -- vertices of mesh
        %   faces ----- connectivity of triangles [v1,v2,v3]
        %------------------------------------------------------------------
        
                obj.numVertices = size(tempVertices,1);
                obj.numFaces = size(tempFaces,1);
                obj.numEdges = 3/2*obj.numFaces;
                obj.numHalfEdges = 3*obj.numFaces;
                obj.numNodes = obj.numVertices;
                obj.numFaceFields = 0;
                obj.numNodeFields = 0;
                
                obj.coordinates = tempVertices;
                obj.vertexHalfEdges = zeros(obj.numVertices,1);
                obj.next = zeros(obj.numHalfEdges,1);
                obj.pair = zeros(obj.numHalfEdges,1);
                obj.ends = zeros(obj.numHalfEdges,1);
                
                for i = 1:obj.numFaces
                    obj.next(3*i-2:3*i) = [3*i-1;3*i;3*i-2];
                    obj.ends(3*i-2:3*i) = tempFaces(i,:)';
                    obj.vertexHalfEdges(tempFaces(i,:)) = [3*i-2,3*i-1,3*i];
                end
                obj.faces = [obj.ends(1:3:end),...
                             obj.ends(2:3:end),...
                             obj.ends(3:3:end)];
                         
                % The rest of this is going to find the pairs
                %--------------------------------------------
                % half edges pointing to each vertex
                allVertexHalfEdges{obj.numVertices}=[];
                for i = 1:obj.numHalfEdges
                    vertexIndex = obj.ends(i);
                    allVertexHalfEdges{vertexIndex} = [allVertexHalfEdges{vertexIndex};i];
                end
                
                % for each half edge we're going to find the previous half
                % edge in the face loop which points to the pair vertex.
                % Then we'll loop over the half edges that point to that
                % vertex to find the correct pair edge
                for i = 1:obj.numHalfEdges
                    
                    if obj.pair(i)==0
                        
                        % starting vertex for half edge i
                        nexti = obj.next(i);
                        previ = obj.next(nexti);
                        prevVertex = obj.ends(previ);
                        
                        % half edges pointing to start vertex 
                        candidates = allVertexHalfEdges{prevVertex};
                        for j = 1:length(candidates)
                            
                            candidatej=candidates(j);
                            nextj = obj.next(candidatej);
                            prevj = obj.next(nextj);
                            pairVertex = obj.ends(prevj);
                            
                            if pairVertex == obj.ends(i)
                                obj.pair(i) = candidatej;
                                obj.pair(candidatej) = i;
                                break
                            end
                            
                        end
                    end
                end
                
                obj.isFixedVertex   = zeros(obj.numVertices,1);
                obj.isFixedHalfEdge = zeros(obj.numHalfEdges,1);

                obj.degree = 1;
                obj.isCurved = false;
                
                obj.volume = obj.calculateVolume();
                obj.centroid = obj.calculateCentroid();
                obj.surfaceArea = obj.calculateSurfaceArea();
                obj.resolution = obj.calculateResolution(); 
                
                obj.cleanDeadVertices();
                obj.cleanInvertedFaces();

                obj.isValid()
                
        end
        function cleanDeadVertices(obj)
        % cleans out vertices that aren't included in connectivity
        %------------------------------------------------------------------
            
            if any(obj.vertexHalfEdges==0)
                disp("-----------------------------------------------------------------")
                warning('dead vertices detected ... attempting to cull and recover the mesh')
                disp("-----------------------------------------------------------------")
                vertexIndices = 1:obj.numVertices;
                oldVertexIndices = vertexIndices(obj.vertexHalfEdges~=0);
                newVertexIndices = 1:length(oldVertexIndices);
                
                disp(" Cleaning out dead vertices ")
                disp(' ')
                disp([' total vertices : ',num2str(obj.numVertices)])
                disp([' live vertices  : ',num2str(newVertexIndices(end))])
                disp(' ')
                
                % correct the vertex index half-edges point to
                obj.replaceHalfEdgeVertex(oldVertexIndices,...
                    newVertexIndices);
                
                % clip coordinates to living set
                obj.coordinates = obj.coordinates(oldVertexIndices,:);
                obj.vertexHalfEdges = obj.vertexHalfEdges(oldVertexIndices,:);
                
                % update vertex ids in our face data structure
                obj.faces = [obj.ends(1:3:end),...
                    obj.ends(2:3:end),...
                    obj.ends(3:3:end)];
                
                % correct our numVertices and num nodes
                obj.numVertices = length(obj.vertexHalfEdges);
                obj.numNodes = obj.numVertices;
            end
        end
        function cleanInvertedFaces(obj)
        % If there are any folds in the mesh this will smooth them out
        %------------------------------------------------------------------
        % we want to mostly leave the imported mesh alone so we find any 
        % locations where the mesh is potentially inverted and smooth in 
        % that region
        %------------------------------------------------------------------
            
            vertexMaxDihedrals = obj.vertexMaxDihedralAngles();
            if max(vertexMaxDihedrals)>160
                disp("-----------------------------------------------------------------")
                warning('detected inverted face(s) attempting to smooth and recover')
                disp("-----------------------------------------------------------------")
         
                % isolate our metric around inverted faces
                vertexMaxDihedrals(vertexMaxDihedrals<90)=0;
                vertexMaxDihedrals = vertexMaxDihedrals/max(vertexMaxDihedrals);

                % diffuse metric around those points
                for i = 1:10
                    faceMaxDihedrals = obj.faceValuesFromVertices(vertexMaxDihedrals);
                    vertexMaxDihedrals = obj.vertexValuesFromFaces(faceMaxDihedrals);
                    vertexMaxDihedrals = vertexMaxDihedrals/max(vertexMaxDihedrals);
                end
                
                % apply smoothing in these zones
                iter = 1;
                vertexMaxDihedralsTrigger = obj.vertexMaxDihedralAngles();
                while max(vertexMaxDihedralsTrigger)>90 && iter < 10
                    obj.smooth(1,'cotan',false,vertexMaxDihedrals);
                    vertexMaxDihedralsTrigger = obj.vertexMaxDihedralAngles();
                    disp(['smoothing iteration: ',num2str(iter),'  number of inverted faces: ',num2str(length(find(vertexMaxDihedralsTrigger>90)))]);
                    iter = iter +1;
                end
                
                % throw a warning if we dont reach our quality metric
                if iter ==10
                    warning('Max smoothing iters exceeded. Inverted faces may still exist')
                end
                
            end
        end
        function setDefaultOptions(obj)
        % sets our default options for coarsening, refining, smoothing
        %------------------------------------------------------------------
            
            obj.coarsenOptions.method = 'uniform';
            obj.coarsenOptions.averageCollapseVertexPosition = true;
            obj.coarsenOptions.reproject = true;
            obj.coarsenOptions.progressThreshold = 1000;
            
            obj.refineOptions.dihedralAngleThreshold = 9;
            obj.refineOptions.numIterations = 1;
            obj.refineOptions.maxLevel = 1;
            
            obj.smoothOptions.method = 'cotangent';
            obj.smoothOptions.reproject = false;
        end

        function setFixedVertices(obj,vertices)
        % set vertices that can't be changed
        %------------------------------------------------------------------

            obj.isFixedVertex(vertices) = 1;
            for i = 1:length(vertices)
                vertexi = vertices(i);
                fixedHalfEdgesi = obj.spokeHalfEdges(vertexi);
                obj.isFixedHalfEdge(fixedHalfEdgesi) = 1;
            end

            obj.numFixedVertices  = sum(obj.isFixedVertex);
            obj.numFixedHalfEdges = sum(obj.isFixedHalfEdge);
        end
        function clearFixedVertices(obj)
        % set vertices that can't be changed
        %------------------------------------------------------------------
            obj.isFixedVertex(:) = 0;
            obj.isFixedHalfEdge(:) = 0;

            obj.numFixedVertices  = 0;
            obj.numFixedHalfEdges = 0;
        end

        function addFaceField(obj,fieldData,fieldName)
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
            if size(fieldData,1)~=obj.numFaces
                error('faceFieldData must have num of row equal to numFaces')
            end
            obj.numFaceFields = obj.numFaceFields + 1;
            if nargin == 2
                fieldName = ['field_',num2str(obj.numFaceFields)];
            elseif nargin ~=3
                error('incorrect number of inputs')
            end
            obj.faceFields{obj.numFaceFields}.data = fieldData;
            obj.faceFields{obj.numFaceFields}.name = fieldName;
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
            obj.numFaceFields = 0;
            
            obj.nodeFields = [];
            obj.faceFields = [];
        end
        
        function [edges,edgeHalfEdges,halfEdgeEdges] = edges(obj)
        % unique edges within mesh
        %------------------------------------------------------------------
        % if the pair index is greater than the half edge index we take it
        % as a unique edge
        %------------------------------------------------------------------
        % Outputs:
        %   edges --------- vertex ids bounding each unique edge
        %   edgeHalfEdges - halfEdge ids assoc w/ each unique edge
        %   halfEdgeEdges - unique edge ids assoc w/ each half edge
        %------------------------------------------------------------------
        
            halfEdgeEdges = zeros(obj.numHalfEdges,1);
            
            
            indices = linspace(1,obj.numHalfEdges,obj.numHalfEdges)';
            rootEdges = obj.pair(obj.pair>indices);
            pairEdges = obj.pair(rootEdges);
            
            edges  = [obj.ends(rootEdges),...
                      obj.ends(pairEdges)];
             
            edgeHalfEdges = [rootEdges,pairEdges];  
            
            uniqueEdgeIndices = linspace(1,obj.numEdges,obj.numEdges)';
            halfEdgeEdges(rootEdges) = uniqueEdgeIndices;
            halfEdgeEdges(pairEdges) = uniqueEdgeIndices;
            
        end
        
        function faceValues = faceValuesFromVertices(obj,vertexValues)
        % convert vertex field to face field
        %------------------------------------------------------------------
        % Inputs:
        %   vertexValues -- field values at vertices
        %------------------------------------------------------------------
        % Outputs:
        %   faceValues --- field values at faces
        %------------------------------------------------------------------
            
            p1 = obj.ends(1:3:end);
            p2 = obj.ends(2:3:end); 
            p3 = obj.ends(3:3:end); 
           
            faceValues = 1/3*(vertexValues(p1,:)+...
                              vertexValues(p2,:)+...
                              vertexValues(p3,:)); 
        end
        function areaVectors = faceAreaVectors(obj)
        % outward facing area vectors of faces
        %------------------------------------------------------------------
        % Outputs:
        %   areaVectors - unit normals scaled by face area
        %------------------------------------------------------------------
            
            % points defining each triangle
            p1 = obj.coordinates(obj.ends(1:3:end),:);
            p2 = obj.coordinates(obj.ends(2:3:end),:); 
            p3 = obj.coordinates(obj.ends(3:3:end),:); 
           
            areaVectors = 1/2*cross(p2-p1,p3-p2,2); 
        end
        function normals = faceNormals(obj)
        % outward facing unit normals of faces
        %------------------------------------------------------------------
        % Outputs:
        %   normals - unit normal vectors for each face
        %------------------------------------------------------------------
            a = obj.faceAreaVectors();
            normals = a./sqrt(a(:,1).^2+a(:,2).^2+a(:,3).^2);
        end
        function area = faceAreas(obj)
        % areas of faces
        %------------------------------------------------------------------
        % Outputs:
        %   area - areas of faces
        %------------------------------------------------------------------
            
            a = obj.faceAreaVectors();
            area = sqrt(a(:,1).^2+a(:,2).^2+a(:,3).^2);
            
        end
        function centroids = faceCentroids(obj)
        % centroids of faces -- degree 1
        %------------------------------------------------------------------
        % Outputs:
        %   centroids - coordinates of centroids
        %------------------------------------------------------------------
            p1 = obj.coordinates(obj.ends(1:3:end),:);
            p2 = obj.coordinates(obj.ends(2:3:end),:); 
            p3 = obj.coordinates(obj.ends(3:3:end),:); 
           
            centroids=1/3*(p1+p2+p3); 
            
        end
        function neighbors = faceNeighbors(obj)
        % indices of 3 neighbor face for each face
        %------------------------------------------------------------------
        % Half edges are in order for each face e.g. half edges 1-3 belong
        % for face 1 and so on. Here we round the indices of the pair half
        % edges for each face to get the indices of the neighbor faces.
        %------------------------------------------------------------------
        % Outputs:
        %   neighbors - neighbor indices
        %------------------------------------------------------------------
            
            he1 = obj.pair(1:3:end);
            he2 = obj.pair(2:3:end);
            he3 = obj.pair(3:3:end);
            
            neighbors = ceil([he1/3,he2/3,he3/3]);
            
        end
        function alpha = faceAverageMeshAngles(obj)
        % calculates average angle between each face w/ neighbors (degrees)
        %------------------------------------------------------------------
        % Outputs:
        %   alpha - average angle in degrees between face normal and
        %           neighboring face normals (degrees)
        %------------------------------------------------------------------
        
            neighbors  = obj.faceNeighbors();
            nf = obj.faceNormals();
            
            n1 = nf(neighbors(:,1),:); 
            n2 = nf(neighbors(:,2),:); 
            n3 = nf(neighbors(:,3),:);
            
            alpha0 = sum(nf.*(n1+n2+n3),2)/3;
            alpha = acosd(min(alpha0,1)); 
            
        end
        function angle = faceMinIncludedAngles(obj)
        % calculates minimum included ange for all triangles
        %------------------------------------------------------------------
        % Ouputs:
        %   angle - angle in degrees
        %------------------------------------------------------------------
        
            p1 = obj.coordinates(obj.ends(1:3:end),:);
            p2 = obj.coordinates(obj.ends(2:3:end),:);
            p3 = obj.coordinates(obj.ends(3:3:end),:);
        
            v1 = p2-p1;
            v2 = p3-p2;
            v3 = p1-p3;
        
            angles = zeros(obj.numFaces,3);
            angles(:,1) = acosd(dot(v2,-v1,2)./...
                               (vecnorm(v2,2,2).*vecnorm(v1,2,2)) );
            angles(:,2) = acosd(dot(v2,-v3,2)./...
                               (vecnorm(v2,2,2).*vecnorm(v3,2,2)) );
        
            angles(:,3) = 180-(angles(:,1)+angles(:,2));
           
           
            angle = min(angles,[],2);              
       
        end
        function angle = faceMaxIncludedAngles(obj)
        % calculates maximum included ange for all triangles
        %------------------------------------------------------------------
        % Ouputs:
        %   angle - angle in degrees
        %------------------------------------------------------------------
        
            p1 = obj.coordinates(obj.ends(1:3:end),:);
            p2 = obj.coordinates(obj.ends(2:3:end),:);
            p3 = obj.coordinates(obj.ends(3:3:end),:);
        
            v1 = p2-p1;
            v2 = p3-p2;
            v3 = p1-p3;
        
            angles = zeros(obj.numFaces,3);
            
            angles(:,1) = acosd(dot(v2,-v1,2)./...
                               (vecnorm(v2,2,2).*vecnorm(v1,2,2)) );
            angles(:,2) = acosd(dot(v2,-v3,2)./...
                               (vecnorm(v2,2,2).*vecnorm(v3,2,2)) );
        
            angles(:,3) = 180-(angles(:,1)+angles(:,2));
           
            angle = max(angles,[],2);              
       
        end
        
        function angles = halfEdgeAngles(obj)
        % angle between adjacent faces for each half edge
        %------------------------------------------------------------------
        % Outputs:
        %    angles - in degrees
        %------------------------------------------------------------------
            nf = obj.faceNormals();
            faceA = ceil(1/3*linspace(1,obj.numHalfEdges,obj.numHalfEdges));
            faceB = ceil(obj.pair/3);
            
            angles = acosd(dot(nf(faceA,:),nf(faceB,:),2));
            
            
        end
        function lengths = halfEdgeLengths(obj,ids)
        % angle between adjacent faces for each half edge
        %------------------------------------------------------------------
        % Inputs:
        %   ids ----- optional one can specify a subset of half edges
        %------------------------------------------------------------------
        % Outputs:
        %    angles - in degrees
        %------------------------------------------------------------------
            if nargin==2
                p1 = obj.coordinates(obj.ends(ids),:);
                p2 = obj.coordinates(obj.ends(obj.pair(ids)),:);
                v = p2-p1;
                lengths = sqrt(v(:,1).^2+v(:,2).^2+v(:,3).^2);
            else
                p1 = obj.coordinates(obj.ends,:);
                p2 = obj.coordinates(obj.ends(obj.pair),:);
                v = p2-p1;
                lengths = sqrt(v(:,1).^2+v(:,2).^2+v(:,3).^2);
            end
            
        end
        
        function areaVectors = edgeAreaVectors(obj)
        % returns areaVectors for each edge, avg over adjacent faces
        %------------------------------------------------------------------
        % We take the indices of the half edges associated w/ each unique
        % edge convert that to the face indices. Left and right really have
        % no directional meaning here just faces on opposite side of edges
        %------------------------------------------------------------------
        % Outputs:
        %   AreaVectors - area vectors
        %------------------------------------------------------------------
            
            a = obj.faceAreaVectors();  
            [~,edgeHalfEdges,~] = obj.edges();
            
            left  = ceil(edgeHalfEdges(:,1)/3); 
            right = ceil(edgeHalfEdges(:,2)/3);
 
            areaVectors = (a(left,:)+a(right,:))/3;

        end
        function normals = edgeNormals(obj)
        % returns normals for each edge, avg over adjacent faces
        %------------------------------------------------------------------
        % Outputs:
        %   normals - unit normal area-weight averaged over adj faces
        %------------------------------------------------------------------
            a = obj.edgeAreaVectors();  % area of faces
            normals = a./sqrt(a(:,1).^2+a(:,2).^2+a(:,3).^2);
        end
        
        function vertexValues = vertexValuesFromFaces(obj,faceValues)
        % area weight vertex values from face centroid values
        %------------------------------------------------------------------
        % Inputs:
        %   faceValues -- field values at faces 
        %------------------------------------------------------------------
        % Outputs:
        %   AreaVectors - fieldValues at 
        %------------------------------------------------------------------
        
            areas = obj.faceAreas;
            
            startingEdge = obj.vertexHalfEdges;
            nexti = obj.next(startingEdge);
            pairi = obj.pair(nexti);
            
            areaSum = areas(ceil(startingEdge/3));
            vertexValues = faceValues(ceil(startingEdge/3),:).*areas(ceil(startingEdge/3));
            
            ids = linspace(1,obj.numVertices,obj.numVertices)';
   
            while  any(pairi ~= startingEdge)
                    keepLoops = pairi ~= startingEdge;
                    pairi = pairi(keepLoops); 
                    startingEdge = startingEdge(keepLoops); 
                    ids = ids(keepLoops);
                    
                    areaSum(ids) = areaSum(ids) +...
                                     areas(ceil(pairi/3));
                                 
                    vertexValues(ids,:) = vertexValues(ids,:) +...
                                          faceValues(ceil(pairi/3),:).*areas(ceil(pairi/3));
                                      
                    nexti = obj.next(pairi);
                    pairi = obj.pair(nexti);
            end   
            vertexValues = vertexValues ./ areaSum;  
        end
        function vertexAreaVectors = vertexAreaVectors(obj)
        % returns areaVectors for each vertex, avg over adjacent faces
        %------------------------------------------------------------------
        % Vectorization of vertex loops to calculate the area vector
        % associated with each vertex. Steps to next face through half edge
        % relations and culls the vertex loop if it reaches the starting
        % point
        %------------------------------------------------------------------
        % Outputs:
        %   AreaVectors - area vectors
        %------------------------------------------------------------------
        
            faceAreasVectors = obj.faceAreaVectors()/3; % 1/3 vertex weighting here
               
            startingEdge = obj.vertexHalfEdges;
            nexti = obj.next(startingEdge);
            pairi = obj.pair(nexti);
            
            vertexAreaVectors = faceAreasVectors(ceil(startingEdge/3),:);
            
            ids = linspace(1,obj.numVertices,obj.numVertices)';
            
            while  any(pairi ~= startingEdge)
                    keepLoops = pairi ~= startingEdge;
                    pairi = pairi(keepLoops); 
                    startingEdge = startingEdge(keepLoops); 
                    ids = ids(keepLoops);

                    vertexAreaVectors(ids,:) = vertexAreaVectors(ids,:) +...
                                               faceAreasVectors(ceil(pairi/3),:);
                    nexti = obj.next(pairi);
                    pairi = obj.pair(nexti);
            end
                

        end
        function normals = vertexNormals(obj)
        % returns normals for each vertex, avg over adjacent faces
        %------------------------------------------------------------------
        % Outputs:
        %   normals - unit normal area-weight averaged over adj faces
        %------------------------------------------------------------------
            faceNormals = obj.faceNormals(); 
               
            startingEdge = obj.vertexHalfEdges;
            nexti = obj.next(startingEdge);
            pairi = obj.pair(nexti);
            
            normals = faceNormals(ceil(startingEdge/3),:);
            numSum = ones(obj.numVertices,1);
            
            ids = linspace(1,obj.numVertices,obj.numVertices)';
            
            while  any(pairi ~= startingEdge)
                    keepLoops = pairi ~= startingEdge;
                    pairi = pairi(keepLoops); 
                    startingEdge = startingEdge(keepLoops); 
                    ids = ids(keepLoops);

                    normals(ids,:) = normals(ids,:) +...
                                               faceNormals(ceil(pairi/3),:);
                                           
                    numSum(ids,:) = numSum(ids,:) + 1;
                    
                    nexti = obj.next(pairi);
                    pairi = obj.pair(nexti);
            end
            
            normals = normals./vecnorm(normals,2,2);
        end
        function numSpokes = vertexSpokeCounts(obj,vertexSubset)
        % get the number of unique edge attached to each vertex
        %------------------------------------------------------------------
        % Inputs:
        %   activeVertices -- allows calculation on subset of vertices
        %------------------------------------------------------------------
        % Outputs:
        %   numSpokes (numVertices x 1) integer w/ number of spokes
        %------------------------------------------------------------------
        
            startingEdge = obj.vertexHalfEdges; % first edge in loop
            nexti = obj.next(startingEdge);     % next in face loop
            pairi = obj.pair(nexti);            % next in vertex loop
            
            % decide if we're calcing for all or subset
            %--------------------------------------------------------------
            if nargin>1
                numActiveVertices = length(vertexSubset);
                numSpokes = ones(numActiveVertices,1);
                ids = linspace(1,numActiveVertices,numActiveVertices)';
                pairi = pairi(vertexSubset);
                startingEdge = startingEdge(vertexSubset);
            else
                ids = linspace(1,obj.numVertices,obj.numVertices)';
                numSpokes = ones(obj.numVertices,1);
            end

            % iterate the vertex loops
            %--------------------------------------------------------------
            while  any(pairi ~= startingEdge)
                    keepLoops = pairi ~= startingEdge;
                    pairi = pairi(keepLoops); 
                    startingEdge = startingEdge(keepLoops); 
                    ids = ids(keepLoops);
                    numSpokes(ids,1) = numSpokes(ids,1) + 1;
                    nexti = obj.next(pairi);
                    pairi = obj.pair(nexti);
            end
            
        end 
        function angles = vertexAverageAngles(obj,vertexSubset)
        % average of average angles for adjacent faces
        %------------------------------------------------------------------
        % Uses a vectorized version of a standard half-edge mesh vertex
        % cycle updating the ids tracking vertex loops that haven't been
        % completed yet.
        %------------------------------------------------------------------
        %             o p3             | Look through and calcuate the 
        %              ^               | dihedral angle of each spoke edge
        %               \              | and each "rim" edge. Spoke edges
        %           n2   \ v32         | are weighted twice as much as rim
        %                 \            | edges.
        %           v24    \           | 
        %    p4 o---------->o p2       | p1 is the central vertex we are
        %      / ^                     | calculating the normal for.
        %     /   \                    |
        %    /     \ v41   n1          | 
        %   / v45   \                  |
        %  v         \                 |
        % o p5   n3   o p1             | 
        %------------------------------------------------------------------
        % Outputs:
        %   angles -- average of average mesh angle for surrounding faces
        %-----------------------------------------------------------------
                        
            %calculate subset if indices are input
            if nargin == 1
                startingEdge = obj.vertexHalfEdges();
                coord1 = obj.coordinates;
            elseif nargin == 2
                startingEdge = obj.vertexHalfEdges(vertexSubset);
                coord1 = obj.coordinates(vertexSubset,:);
            end
            numInitialVertices = length(startingEdge);
            
            % initialize i-1 cycle things
            pairi = startingEdge;
            
            p2 = obj.ends(obj.pair(startingEdge));
            p4 = obj.ends(obj.next(obj.pair(startingEdge)));
            coord2 = obj.coordinates(p2,:);
            coord4 = obj.coordinates(p4,:);
            v41 = coord4-coord1;
            v24 = coord2-coord4;
            n1 = obj.orthogonalNormal(v24,v41);
            
            % initialize next step
            coord4 = coord2; 
            n3 = n1;
            
            % initialize our outputs
            angles = ones(numInitialVertices,1);
            weight = ones(numInitialVertices,1);
            
            ids = linspace(1,numInitialVertices,numInitialVertices)';
            
            count = 0;
            while  any(pairi ~= startingEdge) || count == 0
                
                if count>0
                    % clip to incomplete cycles
                    keepLoops = pairi ~= startingEdge;
                    pairi = pairi(keepLoops);
                    startingEdge = startingEdge(keepLoops);
                    ids = ids(keepLoops);
                
                    % from prev cycle
                    coord1 = coord1(keepLoops,:);
                    coord4 = coord2(keepLoops,:);
                    n3 = n1(keepLoops,:);
                end
                
                % spokes
                nexti = obj.next(pairi);
                pairi = obj.pair(nexti);
                
                % relavant coordinates
                p2 = obj.ends(nexti);
                p3 = obj.ends(obj.next(obj.pair(obj.next(nexti))));
                coord2 = obj.coordinates(p2,:);
                coord3 = obj.coordinates(p3,:);

                % edge vectors
                v41 = coord4 - coord1;
                v24 = coord2 - coord4;
                v32 = coord3 - coord2;
                
                % face normals
                n1 = obj.orthogonalNormal(v24,v41);
                n2 = obj.orthogonalNormal(v24,v32);

                % 2 x spoke dihedral 1x rim dihedral angle
                angles(ids,:) = angles(ids,:) + sum(n1.*(2*n3+n2),2);
                weight(ids,:) = weight(ids,:) + 3;
                
                count = count + 1;
            end
            
            angles = acosd(min(angles./weight,1));
        end
        function angles = vertexMaxDihedralAngles(obj,vertexSubset)
        % max dihedral angle of spokes around vertices
        %------------------------------------------------------------------
        % Uses a vectorized version of a standard half-edge mesh vertex
        % cycle updating the ids tracking vertex loops that haven't been
        % completed yet.
        %------------------------------------------------------------------
        % Outputs:
        %   angles -- average of average mesh angle for surrounding faces
        %-----------------------------------------------------------------
                        
            %calculate subset if indices are input
            if nargin == 1
                startingEdge = obj.vertexHalfEdges();
                coord1 = obj.coordinates;
            elseif nargin == 2
                startingEdge = obj.vertexHalfEdges(vertexSubset);
                coord1 = obj.coordinates(vertexSubset,:);
            end
            numInitialVertices = length(startingEdge);
            
            % initialize i-1 cycle things
            pairi = startingEdge;
            
            p2 = obj.ends(obj.pair(startingEdge));
            p4 = obj.ends(obj.next(obj.pair(startingEdge)));
            coord2 = obj.coordinates(p2,:);
            coord4 = obj.coordinates(p4,:);
            v41 = coord4-coord1;
            v24 = coord2-coord4;
            n1 = obj.orthogonalNormal(v24,v41);
            
            % initialize next step
            coord4 = coord2; 
            n3 = n1;
            
            % initialize our outputs
            angles = ones(numInitialVertices,1);
            
            ids = linspace(1,numInitialVertices,numInitialVertices)';
            
            count = 0;
            while  any(pairi ~= startingEdge) || count == 0
                
                if count>0
                    % clip to incomplete cycles
                    keepLoops = pairi ~= startingEdge;
                    pairi = pairi(keepLoops);
                    startingEdge = startingEdge(keepLoops);
                    ids = ids(keepLoops);
                
                    % from prev cycle
                    coord1 = coord1(keepLoops,:);
                    coord4 = coord2(keepLoops,:);
                    n3 = n1(keepLoops,:);
                end
                
                % spokes
                nexti = obj.next(pairi);
                pairi = obj.pair(nexti);
                
                % relavant coordinates
                p2 = obj.ends(nexti);
                coord2 = obj.coordinates(p2,:);

                % edge vectors
                v41 = coord4 - coord1;
                v24 = coord2 - coord4;
                
                % face normals
                n1 = obj.orthogonalNormal(v24,v41);

                % 2 x spoke dihedral 1x rim dihedral angle
                angles(ids) = min(angles(ids),sum(n1.*n3,2));
                count = count + 1;
            end
            
            angles = acosd(min(angles,1));
        end
        
        function areaVectors = nodeAreaVectors(obj)
        % method to return of area vectors assoc with all mesh nodes 
        %------------------------------------------------------------------
        % Generalization to calculate the area vectors of nodes defining
        % curvilinear triangular mesh assumed nodes are equidistant with
        % configuration consistent w/ that of the NewtonCotes and
        % LangrangeInterpolants functions.
        %
        % for the curvilinear case the product of the basis function
        % assocated with each node and the differential surface vector must
        % be integrated over all faces the nodes is associated with. This
        % is accomplished by using a higher degree of exactness quadrature
        % rule that can integrate the product exactly.
        %
        % for the rectilinear case, the differential surface vector is
        % constant for each face and the integration can be performed with
        % an NC rule of equivalent degree to that of the mesh.
        %------------------------------------------------------------------
        % Outputs:
        %   areaVectors - area vector assoc w/ each node
        %------------------------------------------------------------------
        
            areaVectors = zeros(obj.numNodes,3);
            
            % Rectilinear Meshes
            if ~obj.isCurved
                Af = obj.faceAreaVectors();
                [~,~,w] = NewtonCotesTriangle(obj.degree);
            
                for i=1:obj.numFaces
                    areaVectors(obj.faces(i,:),:) = areaVectors(obj.faces(i,:),:) + w*Af(i,:);
                end
                
            % Curvilinear Meshes
            else
                
                % quadrature rule for unit triangle (Nix1)
                d = 3*obj.degree - 1;
                [u,v,w] = NewtonCotesTriangle(d);
                
                % Ng interpolants/derivatives at Ni quadrature points (NixNg)
                [ phi, dphidu, dphidv ] = LagrangeInterpolantsTriangle( u,v,obj.degree );
                
                %-- x,y,z coords of nodes separated in (Nf x Ng) format
                numNodesPerFace = size(obj.faces,2);
                Xx = zeros(obj.numFaces,numNodesPerFace);
                Xy = zeros(obj.numFaces,numNodesPerFace);
                Xz = zeros(obj.numFaces,numNodesPerFace);
                
                for i = 1:numNodesPerFace
                    Xx(:,i) = obj.coordinates(obj.faces(:,i),1);
                    Xy(:,i) = obj.coordinates(obj.faces(:,i),2);
                    Xz(:,i) = obj.coordinates(obj.faces(:,i),3);
                end
                
                %--  (Ni x Nf) = (NixNg) * (NgxNf)
                Ax = (dphidu*Xy').*(dphidv*Xz') - (dphidu*Xz').*(dphidv*Xy');
                Ay = (dphidu*Xz').*(dphidv*Xx') - (dphidu*Xx').*(dphidv*Xz');
                Az = (dphidu*Xx').*(dphidv*Xy') - (dphidu*Xy').*(dphidv*Xx');
                
                
                % nodes indices are not uniquely valued column wise i.e.
                % each node might be the vertex-1 for multiple faces so we
                % gotta loop over the faces.
                %-------------------------------------------
                % (Nix1)' + ( (NixNg) .* Nix1)
                for i = 1:obj.numFaces
                    nodeIndices = obj.faces(i,:); 
                    areaVectors(nodeIndices,1) = areaVectors(nodeIndices,1) + (phi.*Ax(:,i))'*w;
                    areaVectors(nodeIndices,2) = areaVectors(nodeIndices,2) + (phi.*Ay(:,i))'*w;
                    areaVectors(nodeIndices,3) = areaVectors(nodeIndices,3) + (phi.*Az(:,i))'*w;
                end
                
                areaVectors = areaVectors/2.0;
            end
            
        end
        function normals = nodeNormals(obj)
        % unit normals for the mesh nodes
        %------------------------------------------------------------------
        % the nodeAreaVector method does the heavy lifting, this just wraps
        % and normalizes the result.
        %------------------------------------------------------------------
        % Outputs:
        %   normals -- outward facing surface normals
        %------------------------------------------------------------------
            
            % Rectilinear
            if ~obj.isCurved
                nf = obj.faceNormals();
                normals = zeros(obj.numNodes,3);
                for i=1:obj.numFaces
                    normals(obj.faces(i,:),:) = normals(obj.faces(i,:),:) + nf(i,:);
                end
            
            % Curvilinear
            else
                normals = obj.nodeAreaVectors(); 
            end
            normals = normals./...
                 sqrt(normals(:,1).^2+normals(:,2).^2+normals(:,3).^2);
        end
        
        function resetBulkProperties(obj)
        % resets stored bulk properties
        %------------------------------------------------------------------
            obj.volume = obj.calculateVolume();
            obj.surfaceArea = obj.calculateSurfaceArea();
            obj.centroid = obj.calculateCentroid();
            obj.resolution = obj.calculateResolution();
        end
        function volume = calculateVolume(obj)
        % volume contained in mesh
        %------------------------------------------------------------------
        % generalized volume calc for curvilinear meshes. 
        %------------------------------------------------------------------
        % Outputs:
        %   vol - volume
        %------------------------------------------------------------------
            volume = 1/3*sum(dot(obj.nodeAreaVectors(),obj.coordinates,2));
        end
        function surfaceArea = calculateSurfaceArea(obj)
        % total surface area of mesh
        %------------------------------------------------------------------
        % Outputs:
        %   surfaceArea - surface area
        %------------------------------------------------------------------
            
            % Rectilinear Meshes
            if ~obj.isCurved
                 
                 a = obj.faceAreaVectors();
                 surfaceArea = sum(sqrt(a(:,1).^2 + a(:,2).^2 + a(:,3).^2));

            % Curvilinear Meshes
            else
                
                % quadrature rule for unit triangle (Nix1)
                d = 3*obj.degree - 1;
                [u,v,w] = NewtonCotesTriangle(d);
                
                % Ng interpolants/derivatives at Ni quadrature points (NixNg)
                [ ~, dphidu, dphidv ] = LagrangeInterpolantsTriangle( u,v,obj.degree );
                
                %-- x,y,z coords of nodes separated in (Nf x Ng) format
                numNodesPerFace = size(obj.faces,2);
                Xx = zeros(obj.numFaces,numNodesPerFace);
                Xy = zeros(obj.numFaces,numNodesPerFace);
                Xz = zeros(obj.numFaces,numNodesPerFace);
                
                for i = 1:numNodesPerFace
                    Xx(:,i) = obj.coordinates(obj.faces(:,i),1);
                    Xy(:,i) = obj.coordinates(obj.faces(:,i),2);
                    Xz(:,i) = obj.coordinates(obj.faces(:,i),3);
                end
                
                %--  (Ni x Nf) = (NixNg) * (NgxNf)
                Ax = (dphidu*Xy').*(dphidv*Xz') - (dphidu*Xz').*(dphidv*Xy');
                Ay = (dphidu*Xz').*(dphidv*Xx') - (dphidu*Xx').*(dphidv*Xz');
                Az = (dphidu*Xx').*(dphidv*Xy') - (dphidu*Xy').*(dphidv*Xx');
                
                % volume of each quadrature point (Ni x Nf)
                surfaceArea = 0.5*sum(w'*sqrt(Ax.^2 + Ay.^2 + Az.^2));
                
            end
        end
        function centroid = calculateCentroid(obj)
        % centroid of body
        %------------------------------------------------------------------
        % Outputs:
        %   c - 1x3 centroid
        %------------------------------------------------------------------
            
            % Rectilinear Meshes
            if ~obj.isCurved
                vf = dot(obj.nodeAreaVectors(),obj.coordinates,2);
                cf = 0.75*vf'*obj.coordinates;
                centroid = cf/sum(vf);
                
            % Curvilinear Meshes
            else
                
                % quadrature rule for unit triangle (Nix1)
                d = 3*obj.degree - 1;
                [u,v,w] = NewtonCotesTriangle(d);
                
                % Ng interpolants/derivatives at Ni quadrature points (NixNg)
                [ phi, dphidu, dphidv ] = LagrangeInterpolantsTriangle( u,v,obj.degree );
                
                %-- x,y,z coords of nodes separated in (Nf x Ng) format
                numNodesPerFace = size(obj.faces,2);
                Xx = zeros(obj.numFaces,numNodesPerFace);
                Xy = zeros(obj.numFaces,numNodesPerFace);
                Xz = zeros(obj.numFaces,numNodesPerFace);
                
                for i = 1:numNodesPerFace
                    Xx(:,i) = obj.coordinates(obj.faces(:,i),1);
                    Xy(:,i) = obj.coordinates(obj.faces(:,i),2);
                    Xz(:,i) = obj.coordinates(obj.faces(:,i),3);
                end
                
                % coordinates of quadrature points ( NixNf )
                xx = phi*Xx';
                xy = phi*Xy';
                xz = phi*Xz';
                
                %--  (Ni x Nf) = (NixNg) * (NgxNf)
                Ax = (dphidu*Xy').*(dphidv*Xz') - (dphidu*Xz').*(dphidv*Xy');
                Ay = (dphidu*Xz').*(dphidv*Xx') - (dphidu*Xx').*(dphidv*Xz');
                Az = (dphidu*Xx').*(dphidv*Xy') - (dphidu*Xy').*(dphidv*Xx');
                
                % volume of each quadrature point (Ni x Nf)
                voli = xx.*Ax + xy.*Ay + xz.*Az;
                
                % centroids of quadrature points
                cx = xx.*voli;
                cy = xy.*voli;
                cz = xz.*voli;
                
                %  face volume (Nix1)' x (NixNf)
                faceVolume = w'*voli;
                
                % body centroid 
                % 1xNf x (NfxNi) x (Nix1)  x 3 / (1x1) ---> (1x3)
                centroid = 3/4*sum([cx'*w,cy'*w,cz'*w],1)/sum(faceVolume);
                
            end
            
        end
        function resolution = calculateResolution(obj)
            A = obj.faceAreas();
            resolution = sqrt(sum(A)/size(A,1));
            
        end
        function center(obj)
        % resets coordinate origin to centroid of body
        %------------------------------------------------------------------
            obj.coordinates = obj.coordinates - obj.centroid;
            obj.centroid = [0,0,0];
        end
        
        function [centroids,volumes] = extendedTetrahedaCentroids(obj)
        % returns centroids and volumes of tets (tri->mesh centroid)
        %------------------------------------------------------------------
        % this is an extension of Chanut et. al. 2015 's mascon generation
        % approach generalized for curvilinear meshes
        %------------------------------------------------------------------
        % Outputs:
        %   centroids -- centroids of tetrahedra
        %   volumes ---- assoc. volumes
        %------------------------------------------------------------------

         % Rectilinear Meshes
            if ~obj.isCurved
                volumes = 1/3*dot(obj.faceAreaVectors(),obj.faceCentroids(),2);
                centroids = 3/4*obj.faceCentroids();
                
            % Curvilinear Meshes
            else
                
                % quadrature rule for unit triangle (Nix1)
                d = 3*obj.degree - 1;
                [u,v,w] = NewtonCotesTriangle(d);
                
                % Ng interpolants/derivatives at Ni quadrature points (NixNg)
                [ phi, dphidu, dphidv ] = LagrangeInterpolantsTriangle( u,v, obj.degree );
                
                %-- x,y,z coords of nodes separated in (Nf x Ng) format
                numNodesPerFace = size(obj.faces,2);
                Xx = zeros(obj.numFaces,numNodesPerFace);
                Xy = zeros(obj.numFaces,numNodesPerFace);
                Xz = zeros(obj.numFaces,numNodesPerFace);
                
                for i = 1:numNodesPerFace
                    Xx(:,i) = obj.coordinates(obj.faces(:,i),1);
                    Xy(:,i) = obj.coordinates(obj.faces(:,i),2);
                    Xz(:,i) = obj.coordinates(obj.faces(:,i),3);
                end
                
                % coordinates of quadrature points ( NixNf )
                xx = phi*Xx';
                xy = phi*Xy';
                xz = phi*Xz';
                
                %--  (Ni x Nf) = (NixNg) * (NgxNf)
                Ax = (dphidu*Xy').*(dphidv*Xz') - (dphidu*Xz').*(dphidv*Xy');
                Ay = (dphidu*Xz').*(dphidv*Xx') - (dphidu*Xx').*(dphidv*Xz');
                Az = (dphidu*Xx').*(dphidv*Xy') - (dphidu*Xy').*(dphidv*Xx');
                
                % volume of each quadrature point (Ni x Nf)
                voli = xx.*Ax + xy.*Ay + xz.*Az;
                
                % centroids of quadrature points (Ni x Nf)
                cx = xx.*voli;
                cy = xy.*voli;
                cz = xz.*voli;
                
                %  face volume (Nix1)' x (NixNf) ---> Nfx1
                volumes = (1/6*w'*voli)';
                
                % body centroid (NfxNi)*(Nix1) ./ (Nfx1)  ---> Nfx1
                centroids = 1/8*([cx'*w,cy'*w,cz'*w])./volumes;
                
            end
        end
        
        function [ptsq,Aq] = createQuadrature(obj,degreeOfExactness)
        % returns quadrature points and assoc area vectors for different
        % rules. displayQuadratureRules() for options.
        %------------------------------------------------------------------
        % Inputs:
        %   degreeOfExactness - degree of exactness of the quadrature rule
        %                       must be an integer greater than or equal to
        %                       the degree of the mesh. If unspecified then
        %                       the degree of the mesh will be used as the
        %                       degree of exactness or quadrature.
        %  
        %                       also support string inputs if specific
        %                       rules for rectilinear meshes from Pearl and
        %                       Hitt 2020 are desired.                    
        %------------------------------------------------------------------
        % Outputs:
        %   ptsq - quadrature point coordinates
        %   Aq --- area vectors associated w/ each quadrature point
        %------------------------------------------------------------------
            
            % if the mesh is curved or no degree of exactness is specificed
            % were going to base it of the degree of the mesh.
            if obj.isCurved || nargin==1
                
                ptsq = obj.coordinates;
                Aq = obj.nodeAreaVectors();
                
            else
                
                % defaults for specified degree of exactness
                switch degreeOfExactness
                    case 1 % degree 1 default
                        degreeOfExactness='L1';
                    case 2 % degree 2 default
                        degreeOfExactness='B2';
                    case 3 % degree 3 default
                        degreeOfExactness='B3';
                    case 4 % degree 4 default
                        degreeOfExactness='O3';
                end
                
                
                A = obj.faceAreaVectors();
            
                switch degreeOfExactness
                
                    case 'G1' % Guassian degree-1  [Stroud 1971]
                        Aq = A;
                        ptsq = obj.faceCentroids();
                        
                    case 'L1' % Newton-Cotes degree-1 [Lauffer 1955]
                        Aq = obj.vertexAreaVectors();
                        ptsq = obj.coordinates;
                        
                    case 'L2' % Newton-Cotes degree-2 [Stroud 1971]
                        [e,~,~] = obj.edges();
                        Aq = obj.edgeAreaVectors();
                        ptsq = (obj.coordinates(e(:,1),:)+obj.coordinates(e(:,2),:))/2;
                        
                    case 'B2' % Blended degree-2 [Lyness & Jespersen 1975]
                        Av = obj.vertexAreaVectors();
                        Aq = [3/4*A;1/4*Av];
                        ptsq = [obj.faceCentroids();obj.coordinates];
                        
                    case 'G2' % Guassian degree-2 [Stroud 1971]
                        Aq = 1/3*[A;A;A];
                        c = obj.faceCentroids();
                        ptsq = [c+1/2*(obj.coordinates(obj.faces(:,1),:)-c);...
                                c+1/2*(obj.coordinates(obj.faces(:,2),:)-c);...
                                c+1/2*(obj.coordinates(obj.faces(:,3),:)-c)];
                        
                    case 'L3' % Newton-Cotes degree-3 [Lauffer 1955]
                        c = obj.faceCentroids();
                        [e,~,~] = obj.edges();
                        Av = obj.vertexAreaVectors();
                        Ae = obj.edgeAreaVectors();
                        
                        Ae = [Ae;Ae];
                        ptse1 = (2*obj.coordinates(e(:,1),:)+obj.coordinates(e(:,2),:))/3;
                        ptse2 = (obj.coordinates(e(:,1),:)+2*obj.coordinates(e(:,2),:))/3;
                        
                        Aq = [1/10*Av;9/40*Ae;9/20*A];
                        ptsq = [obj.coordinates;ptse1;ptse2;c];
                        
                    case 'B3' % degree-3 [Stroud 1971]
                        c = obj.faceCentroids();
                        [e,~,~] = obj.edges();
                        Av = obj.vertexAreaVectors();
                        Ae = obj.edgeAreaVectors();
                        
                        ptse = (obj.coordinates(e(:,1),:)+obj.coordinates(e(:,2),:))/2;
                        Aq = [3/20*Av;2/5*Ae;9/20*A];
                        ptsq = [obj.coordinates;ptse;c];
                        
                    case 'O4' % degree-4 formula [Lyness & Jespersen 1975]
                        c = obj.faceCentroids();
                        [e,~,~] = obj.edges();
                        Av = obj.vertexAreaVectors();
                        Ae = obj.edgeAreaVectors();
                        
                        d=(3-sqrt(3))/6; % normalized distance vertex->edge pt
                        Ae = [Ae;Ae];
                        ptse1 = (d*obj.coordinates(e(:,1),:)+(1-d)*obj.coordinates(e(:,2),:));
                        ptse2 = ((1-d)*obj.coordinates(e(:,1),:)+d*obj.coordinates(e(:,2),:));
                        Aq = [-1/20*Av;3/10*Ae;9/20*A];
                        ptsq = [obj.coordinates;ptse1;ptse2;c];
                        
                    otherwise
                        disp(' ')
                        disp('Warning: missing or invalid quadrature rule selected.')
                        disp('         use displayQuadratureRules() method')
                        disp('         to see available options:')
                        disp(' ')
                end
        
            end
            
        end
        function inside = isInside(obj,coords)
        % returns true if points are inside mesh (analytic)
        %------------------------------------------------------------------
        % exact up to surface w/out singularities
        %------------------------------------------------------------------
        % Inputs:
        %   p ------ Mx3 array of M sample coordinates
        %------------------------------------------------------------------
        % Outputs:
        %   inside - Mx1 bool, true for internal points
        %------------------------------------------------------------------
                        
            inside = zeros(size(coords,1),1);
            
            coord1 = obj.coordinates(obj.ends(1:3:end),:);
            coord2 = obj.coordinates(obj.ends(2:3:end),:); 
            coord3 = obj.coordinates(obj.ends(3:3:end),:); 
           
            for i=1:size(coords,1)
                
                r1 = coord1-coords(i,:); 
                r2 = coord2-coords(i,:); 
                r3 = coord3-coords(i,:);
                
                r1 = r1./sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2);
                r2 = r2./sqrt(r2(:,1).^2+r2(:,2).^2+r2(:,3).^2);
                r3 = r3./sqrt(r3(:,1).^2+r3(:,2).^2+r3(:,3).^2);
                
                % explicitly writing box product faster Matlab R2018a
                omega = atan2(r1(:,1).*(r2(:,2).*r3(:,3)-r2(:,3).*r3(:,2))+...
                    r1(:,2).*(r2(:,3).*r3(:,1)-r2(:,1).*r3(:,3))+...
                    r1(:,3).*(r2(:,1).*r3(:,2)-r2(:,2).*r3(:,1)),...
                    1+sum(r1.*r2+r3.*(r2+r1),2));
                omega(isinf(omega)|isnan(omega)) = 0;
                
                inside(i) = (sum(omega)-pi)>0.0;
            end
        end
        
        function meshTemp = offsetSurfaceMesh(obj,altitude,numSamplePoints)
        % creates an offset surface mesh at a given altitude
        %------------------------------------------------------------------
        
            if nargin < 3
                numSamplePoints = obj.numVertices;
            end

            altitudei=0.0;
            numWorkingFaces = 2000;
            tol = obj.resolution/100;

            % make a copy
            meshTemp = SurfaceMesh(obj);
            meshTemp.flatten();
            
            % fix this to take care of the case when we want more elements
            if (altitude < obj.resolution*5)
                if (numSamplePoints<obj.numVertices)
                    meshTemp.coarsen(numSamplePoints*2);
                end
                meshTemp.coordinates = meshTemp.coordinates + altitude*meshTemp.vertexNormals();
            else  
                
                if meshTemp.numFaces > numWorkingFaces 
                    meshTemp.coarsen(numWorkingFaces);
                    meshTemp.smooth(5,'uniform',false);
                end
            
            
                % step the mesh out at the resolution/2 until altitude
                %--------------------------------------------------------------
                while altitudei<=altitude-tol
                    
                    stepSize = min(meshTemp.resolution/2,altitude-altitudei);
                    pointNormals = meshTemp.vertexNormals;
                    
                    % need to correct drift away from const alt surf
                    minDist = zeros(meshTemp.numVertices,1);
                    if altitudei > obj.resolution*5
                        for i = 1:meshTemp.numVertices
                            r = meshTemp.coordinates(i,:)-obj.coordinates;
                            rmag = vecnorm(r,2,2);
                            minDist(i,1) = (altitudei-min(rmag))/meshTemp.resolution;
                        end
                    end
                    minDist = max(min(minDist,1.0),-1.0);
                
                    % make sure things behave properly
                    meshTemp.coordinates = meshTemp.coordinates+pointNormals.*(stepSize.*(1+minDist));
                    meshTemp.coarsen(meshTemp.numFaces-24);
                    meshTemp.smooth(5,'uniform',false);
                
                    % keep numFaces within range
                    if meshTemp.numFaces<numWorkingFaces/4
                        meshTemp.refine();
                    end
                    altitudei= altitudei+stepSize;
                
                end
            
                % now were going to step things up to final number of pts
                %----------------------------------------------------------
                meshTemp.setNumVertices(numSamplePoints);
                meshTemp.smooth(10,'uniform',false);
                
                minDist = zeros(meshTemp.numVertices,1);
                for i = 1:meshTemp.numVertices
                    r = meshTemp.coordinates(i,:)-obj.coordinates;
                    rmag = vecnorm(r,2,2);
                    minDist(i,1) = (altitude-min(rmag));
                end
                
                % now make sure were right on the const surf
                %----------------------------------------------------------
                minDistTolerance = 0.01*meshTemp.resolution;
                iter = 0;
                while max(abs(minDist))>minDistTolerance && iter < 20
                    
                    pointNormals = meshTemp.vertexNormals;
                    pointNormals = meshTemp.faceValuesFromVertices(pointNormals);
                    pointNormals = meshTemp.vertexValuesFromFaces(pointNormals);
                    pointNormals = pointNormals./vecnorm(pointNormals,2,2);
                    
                    %stepSize = obj.resolution;
                    
                    minDist = zeros(meshTemp.numVertices,1);
                    
                    for i = 1:meshTemp.numVertices
                        r = meshTemp.coordinates(i,:)-obj.coordinates;
                        rmag = vecnorm(r,2,2);
                        minDist(i,1) = (altitude-min(rmag));
                    end
                    %max(abs(minDist))
                    %minDist = max(min(minDist,1.0),-1.0);
                    meshTemp.coordinates = meshTemp.coordinates+pointNormals.*(minDist);
                    
                    iter = iter+1;
                end
            end
            
        end
        
        function setResolution(obj,newResolution)
        % coarsen/refines mesh to a spec. resolution (length scale)
        %------------------------------------------------------------------
        % Inputs:
        %   newResolution -- length scale of mesh
        %------------------------------------------------------------------
            
           newNumFaces = obj.surfaceArea/newResolution^2;
            
            if newNuFaces < 20
                warning('very low face count might fail')
            end
            warni
            obj.setNumFaces(newNumFaces);

        end
        function setNumFaces(obj,numFaces)
        % wrapper to refine/coarsen to a desired number of faces
        %------------------------------------------------------------------
        % Inputs:
        %   numFaces ------ number of desired faces
        %   fixedVertices - vertices that won't be altered
        %------------------------------------------------------------------

            if round(numFaces/2.0) < obj.numFixedVertices
                error("too many fixed vertices to coarsen that much")
            end

            % if we request more faces refine
            while obj.numFaces<numFaces
                obj.refine();
            end

             % if projection is desired, store our original mesh
            if obj.coarsenOptions.reproject
                originalMesh = SurfaceMesh(obj);
            end

            % if we're are really dropping the node count use the faster 
            % uniform algorthim for the bulk of it.
            if strcmp(obj.coarsenOptions.method,'feature') &&...
                numFaces*10 < obj.numFaces
                obj.coarsen(numFaces*10,'uniform',false)
            end

            % finish things up
            obj.coarsen(numFaces);

            % if projection is desired, do the deed
            if obj.coarsenOptions.reproject
                obj.projectOnTo(originalMesh);
            end

        end
        function setNumVertices(obj,numVertices)
        % wrapper to refine/coarsen to a desired number of Vertices
        %------------------------------------------------------------------
        % Inputs:
        %   numVertices --- number of desired vertices
        %   fixedVertices - vertices that won't be altered
        %------------------------------------------------------------------

            while obj.numVertices < numVertices
                obj.refine();
            end
            obj.setNumFaces(2*numVertices);
        end
        function edgeFlipAll(obj)
        % checks all edges for flippage
        %------------------------------------------------------------------
        %             o p2             | To avoid mem allocations, edges
        %            / ^               | that don't need to be flipped
        %           /   \              | are left in the stack and we 
        %      he3 /     \ he2         | just work our was down the list.
        %         /       \            | When modified edges need to be
        %        v         \           | added to the stack they overwrite
        %    p3 o----he1--->o p1       | entries on the stack that are 
        %        \<--he4---^           | just above our "iterator". 
        %         \       /            |
        %      he5 \     / he6         |
        %           \   /              | 
        %            v /               |
        %             o p4             |
        %------------------------------------------------------------------
         
            flipStack = linspace(1,obj.numHalfEdges,obj.numHalfEdges)';
            flipStack = flipStack(flipStack < obj.pair);
            
            % things needed to track progress
            progressMeter = progressbar('Flipping non-Delaunay edges : ');
            stackLength = length(flipStack);
            stackThresh=length(flipStack);
            stackThreshStep=length(flipStack)/100;
            initialStackLength = stackLength;
            
            % work our way down flipping when we have to
            while stackLength>=1
                
                he1 = flipStack(stackLength);
                flip = obj.flipCriterion(he1);
                
                % track our progress
                if stackLength < stackThresh
                    progressMeter.update(100*(1-stackLength/initialStackLength));
                    stackThresh = stackThresh-stackThreshStep;
                end
                
                if flip
                    he2 = obj.next(he1);
                    he3 = obj.next(he2);
                    he4 = obj.pair(he1);
                    he5 = obj.next(he4);
                    he6 = obj.next(he5);

                    obj.flipEdge(he1);
                    
                    flipStack(stackLength+1:stackLength+3)=[he4,he3,he6];
                    stackLength=stackLength+3;
                else
                    stackLength=stackLength-1;
                end
            end
            
            % reset our connectivity
            obj.faces = [obj.ends(1:3:end),...
                         obj.ends(2:3:end),...
                         obj.ends(3:3:end)];
            
            % clear out invalidated face fields 
            obj.numFaceFields = 0;
            obj.faceFields = [];

            % recalc bulk volume, surface area, etc
            obj.resetBulkProperties();
            progressMeter.close();

        end  
        function curve(obj,mesh)
        % projects nodes of  a degree > 1 mesh onto a finer mesh
        %------------------------------------------------------------------
        % Curves mesh through projection. Requires that the degree of the
        % degree of the mesh be set higher than 1 using setDegree method.
        % mesh volume, centroid and surface area are updated to reflect the
        % change in topology.
        %------------------------------------------------------------------
        % Inputs:
        %   mesh -- SurfaceMesh object more defining more refined mesh 
        %------------------------------------------------------------------
            assert(isa(mesh,'SurfaceMesh'),'must input a SurfaceMesh to curve onto')
            assert(mesh.degree==1,'projection/curving is not set up for curvilinear meshes')
            assert(~mesh.isCurved,'projection/curving is not set up for curvilinear meshes')
            if mesh.numFaces<obj.numFaces
                warning('Input mesh for projection is coarser than original mesh')
            end
            
            normals = obj.nodeNormals();
            obj.coordinates = mesh.project(obj.coordinates, normals);
            obj.isCurved = true;
            
            obj.resetBulkProperties();
            obj.clearFields();
        end
        function setDegree(obj,degree)
        % set polynomial degree of mesh faces
        %------------------------------------------------------------------
        % creates the node and face set-up for curvilinear surface
        % definitions. Mesh will be flattened so any previous curvilinear
        % info will be lost. Resulting mesh will be rectilinear with 
        % degree-d faces.
        %
        % Nodes are in evenly space lattices organized into descending 
        % columns left to right on the unit triangle.
        %
        %------------------------------------------------------------------
        % Inputs:
        %   degree -- degree of mesh (1-4)
        %------------------------------------------------------------------
            if degree ~= obj.degree
                
                obj.flatten();
                obj.degree = degree;

                if degree == 2
                    
                    obj.numNodes = obj.numVertices + obj.numEdges;
                    [edges,~,halfEdgeEdges] = obj.edges();
                    
                    obj.faces = zeros(obj.numFaces,6);
                    
                    obj.faces(:,1) = obj.ends(1:3:end);  
                    obj.faces(:,3) = obj.ends(2:3:end);
                    obj.faces(:,6) = obj.ends(3:3:end);
                    obj.faces(:,2) = halfEdgeEdges(2:3:obj.numHalfEdges)+obj.numVertices;
                    obj.faces(:,4) = halfEdgeEdges(1:3:obj.numHalfEdges)+obj.numVertices;
                    obj.faces(:,5) = halfEdgeEdges(3:3:obj.numHalfEdges)+obj.numVertices;
  
                    obj.coordinates = [obj.coordinates;...
                                    (obj.coordinates(edges(:,1),:)+...
                                     obj.coordinates(edges(:,2),:))/2];
                             
                    
                elseif degree == 3
                    
                    obj.numNodes = obj.numVertices + 2*obj.numEdges + obj.numFaces;
                    
                    [edges,edgeHalfEdges,~] = obj.edges();
                    
                    halfEdgePoints = zeros(obj.numHalfEdges,2);
                    uniqueEdgeIndices = linspace(1,obj.numEdges,obj.numEdges)';
                    rootEdges = edgeHalfEdges(:,1);
                    pairEdges = edgeHalfEdges(:,2);
                    
                    halfEdgePoints(rootEdges,1) = uniqueEdgeIndices+obj.numEdges+obj.numVertices;
                    halfEdgePoints(rootEdges,2) = uniqueEdgeIndices+obj.numVertices;
                    halfEdgePoints(pairEdges,1) = uniqueEdgeIndices+obj.numVertices;
                    halfEdgePoints(pairEdges,2) = uniqueEdgeIndices+obj.numEdges+obj.numVertices;
            

                    obj.faces = zeros(obj.numFaces,10);
                    
                    obj.faces(:,1)  = obj.ends(1:3:end); 
                    obj.faces(:,4)  = obj.ends(2:3:end);      
                    obj.faces(:,10) = obj.ends(3:3:end);
                    
                    obj.faces(:,6) = ((obj.numNodes-obj.numFaces+1):obj.numNodes)';
                    
                    obj.faces(:,2) = halfEdgePoints(2:3:obj.numHalfEdges,1);
                    obj.faces(:,3) = halfEdgePoints(2:3:obj.numHalfEdges,2);
                    obj.faces(:,5) = halfEdgePoints(1:3:obj.numHalfEdges,2);
                    obj.faces(:,7) = halfEdgePoints(3:3:obj.numHalfEdges,1);
                    obj.faces(:,8) = halfEdgePoints(1:3:obj.numHalfEdges,1);
                    obj.faces(:,9) = halfEdgePoints(3:3:obj.numHalfEdges,2);
                    
                    obj.coordinates = [obj.coordinates;...
                                    (2*obj.coordinates(edges(:,1),:)+...
                                     obj.coordinates(edges(:,2),:))/3;...
                                    (obj.coordinates(edges(:,1),:)+...
                                     2*obj.coordinates(edges(:,2),:))/3;...
                                    (obj.coordinates(obj.ends(1:3:end),:)+...
                                     obj.coordinates(obj.ends(2:3:end),:)+...
                                     obj.coordinates(obj.ends(3:3:end),:))/3];
                                 
                elseif degree == 4 
                    
                    obj.numNodes = obj.numVertices + 3*obj.numEdges + 3*obj.numFaces;
                    
                    [edges,edgeHalfEdges,~] = obj.edges();
                    
                    halfEdgePoints = zeros(obj.numHalfEdges,3);
                    
                    uniqueEdgeIndices = linspace(1,obj.numEdges,obj.numEdges)';
                    rootEdges = edgeHalfEdges(:,1);
                    pairEdges = edgeHalfEdges(:,2);
                    
                    halfEdgePoints(rootEdges,1) = uniqueEdgeIndices+2*obj.numEdges+obj.numVertices;
                    halfEdgePoints(rootEdges,2) = uniqueEdgeIndices+obj.numEdges+obj.numVertices;
                    halfEdgePoints(rootEdges,3) = uniqueEdgeIndices+obj.numVertices;
                    halfEdgePoints(pairEdges,1) = uniqueEdgeIndices+obj.numVertices;
                    halfEdgePoints(pairEdges,2) = uniqueEdgeIndices+obj.numEdges+obj.numVertices;
                    halfEdgePoints(pairEdges,3) = uniqueEdgeIndices+2*obj.numEdges+obj.numVertices;

                    obj.faces = zeros(obj.numFaces,15);
                    
                    obj.faces(:,1)  = obj.ends(1:3:end); 
                    obj.faces(:,5)  = obj.ends(2:3:end);      
                    obj.faces(:,15) = obj.ends(3:3:end);
                    
                    startIndex = (obj.numVertices+3*obj.numEdges+1);
                    endIndex = (obj.numVertices+3*obj.numEdges+obj.numFaces);
                    obj.faces(:,7) =  (startIndex:endIndex)';
                    
                    startIndex = endIndex+1;
                    endIndex = endIndex+obj.numFaces;
                    obj.faces(:,8) =  (startIndex:endIndex)';
                    
                    startIndex = endIndex+1;
                    endIndex = endIndex+obj.numFaces;
                    obj.faces(:,11) = (startIndex:endIndex)';
                    
                    obj.faces(:,2) =  halfEdgePoints(2:3:obj.numHalfEdges,1);
                    obj.faces(:,3) =  halfEdgePoints(2:3:obj.numHalfEdges,2);
                    obj.faces(:,4) =  halfEdgePoints(2:3:obj.numHalfEdges,3);
                    obj.faces(:,9) =  halfEdgePoints(3:3:obj.numHalfEdges,1);
                    obj.faces(:,12) = halfEdgePoints(3:3:obj.numHalfEdges,2);
                    obj.faces(:,14) = halfEdgePoints(3:3:obj.numHalfEdges,3);
                    obj.faces(:,13) = halfEdgePoints(1:3:obj.numHalfEdges,1);
                    obj.faces(:,10) = halfEdgePoints(1:3:obj.numHalfEdges,2);
                    obj.faces(:,6) =  halfEdgePoints(1:3:obj.numHalfEdges,3);
                    
                    obj.coordinates = [obj.coordinates;...
                              (3*obj.coordinates(edges(:,1),:)+...
                                 obj.coordinates(edges(:,2),:))/4;...
                                 ...
                              (2*obj.coordinates(edges(:,1),:)+...
                               2*obj.coordinates(edges(:,2),:))/4;...
                               ...
                                (obj.coordinates(edges(:,1),:)+...
                               3*obj.coordinates(edges(:,2),:))/4;...
                                     ...
                              (2*obj.coordinates(obj.ends(1:3:end),:)+...
                                 obj.coordinates(obj.ends(2:3:end),:)+...
                                 obj.coordinates(obj.ends(3:3:end),:))/4;...
                                     ...
                                (obj.coordinates(obj.ends(1:3:end),:)+...
                               2*obj.coordinates(obj.ends(2:3:end),:)+...
                                 obj.coordinates(obj.ends(3:3:end),:))/4;...
                                     ...
                                (obj.coordinates(obj.ends(1:3:end),:)+...
                                 obj.coordinates(obj.ends(2:3:end),:)+...
                               2*obj.coordinates(obj.ends(3:3:end),:))/4;];
                elseif degree ~= 1
                    error('invalid mesh degree specificied, must be 1-4')
                end

            end
            
        end
        function flatten(obj)
        % curvilinear mesh --> rectilinear mesh
        %------------------------------------------------------------------
        % flattens a curvilinear mesh by removing the excess nodes and 
        % converting the face definition to be 3 pts.
        %------------------------------------------------------------------
            if obj.degree > 1
                obj.faces = [obj.ends(1:3:end),obj.ends(2:3:end),obj.ends(3:3:end)];
                obj.coordinates(obj.numVertices+1:obj.numNodes,:) = [];
                obj.numNodes = obj.numVertices;
                obj.isCurved=false;
                obj.degree = 1;
                
                obj.resetBulkProperties();
                obj.clearFields();
            end
            

                         
        end
        
        function coarsen(obj,numFacesCoarse,method,averageCollapse)
        % coarsens mesh attempting to maintain uniform mesh resolution
        %------------------------------------------------------------------
        % Edge collapse algorithm used to remove N vertices maintaining a
        % constant mesh resolution independent of features. Half-edge he1
        % is collapsed towards P3 removing half-edges he1-he6, vertex P1, 
        % and faces fa fb
        %------------------------------------------------------------------
        %            
        %                  o p2                  To make sure the collapse
        %                 / ^                   is valid we'll need to
        %                /   \                  check the four adjacent
        %           he3 /     \ he2             angles around  at p2 and p4
        %              /   fa  \                to ensure the collapse  
        %             v         \               doesn't invert any of the 
        %         p3 o----he1--->o p1           face adjacent to fa and fb.
        %             \<--he4---^               
        %              \       /                The number of spokes is 
        %           he5 \ fb  / he6             checked for p2 and p4 as 
        %                \   /                  well to prevent low spoke
        %                 v /                   count vertices that can
        %                  o p4                 cause trouble later.
        %      
        %------------------------------------------------------------------
        % Inputs:
        %   numFacesCoarse --- number of faces in our post-coarsened mesh
        %   method ----------- override our default coarsening method set
        %                      in coarsenOptions
        %   averageCollapse -- override our default to turn on/off
        %                      averaging of the collapse vertices
        %------------------------------------------------------------------

            % check validity of numFaces
            if obj.numFaces < numFacesCoarse
                error('requested a finer mesh from the coarsen method')
            elseif numFacesCoarse < 50
                warning('requested very coarse mesh... algorithm may fail')
            end
            
            % default inputs if not specified
            if nargin == 2
                method = obj.coarsenOptions.method;
                averageCollapse = obj.coarsenOptions.averageCollapseVertexPosition;
            elseif nargin <= 3
                averageCollapse = obj.coarsenOptions.averageCollapseVertexPosition;
            elseif nargin~=4
                error('Incorrect number of inputs. coarsen(numCoarseFaces,method,averageCollapseVertices)')
            end
            
            % make sure we got a valid string input for our method
            method = lower(method);
            isFeatureBased = strcmp(method,'feature');
            if ~(strcmp(method,'uniform') || strcmp(method,'feature'))
                error('coarsening method must be "uniform" or "feature"')
            end
            
            % set our vertex counts 
            numVerticesCoarse = round(numFacesCoarse/2);
            numDeletedVertices = obj.numVertices-numVerticesCoarse;
            if numVerticesCoarse < obj.numFixedVertices
                error("too many fixed vertices to coarsen that much")
            end
            
            % process our mesh 
            obj.flatten();      % only valid for rectilinear meshes
            obj.edgeFlipAll();  % make sure we're delaunay
            
            % prealloc our flagging arrays
            vertexFlags        = zeros(obj.numVertices,1);
            halfEdgeFlags      = zeros(obj.numHalfEdges,1);
            isModifiedHalfEdge = zeros(obj.numHalfEdges,1);
            
            numSpokes = obj.vertexSpokeCounts();
            lowValenceCandidates = obj.vertexHalfEdges(find(numSpokes==3));
                
            % metrics for edge collapse
            selectionMetric = obj.halfEdgeLengths();
            if isFeatureBased
                angleVertices = obj.vertexAverageAngles();
                edgeAngles = (angleVertices(obj.ends) +...
                              angleVertices(obj.ends(obj.pair)))/2;
                selectionMetric=selectionMetric.*edgeAngles;
            end
            selectionMetric(obj.isFixedHalfEdge==1) = nan;
            [~,sortIndex] = sort(selectionMetric);
            
            % initialize iterators
            iSortVector=1; % track how deep into sorted metric vec we are
            i=1;           % tracks number of collapses
            
            % track our progress
            progressMeter = progressbar('Collapsing Edges            : ');
            progressThreshold = obj.coarsenOptions.progressThreshold;
            lastProgressOutput = 1;

            thresholdResort = min(1000,numDeletedVertices/10);

            while i <= numDeletedVertices
                
                % print our progress 
                if mod(i,progressThreshold)==0 && lastProgressOutput<i
                    lastProgressOutput = i;
                    progressMeter.update(100*i/numDeletedVertices);
                end
                   
                % resort our metric every so often
                if mod(iSortVector,thresholdResort)==0
                    isModifiedHalfEdge(isModifiedHalfEdge==1) = 0;
                    [~,sortIndex] = sort(selectionMetric);
                    iSortVector=1;
                    thresholdResort = floor(max(min(2500,(numDeletedVertices-i)/10),100));
                end

                % if we have three spoke nodes eliminate them else
                % do the normal thing based on length/ mesh angle
                if isempty(lowValenceCandidates) 
                    he1 = sortIndex(iSortVector);
                else 
                    he1 = lowValenceCandidates(end);
                    lowValenceCandidates(end)=[];
                end
                he4=obj.pair(he1);

                % step our sort index
                iSortVector = iSortVector+1;

                % if things checkout do the collapse
                if (isModifiedHalfEdge(he1)==0) && obj.isValidCollapse(he1)

                    % select collapse direction
                    %------------------------------------------------------
                    p1 = obj.ends(he1);
                    p3 = obj.ends(he4);
                    
                    % flip collapse direction if better quality
                    p1Spokes =  obj.spokeHalfEdges(p1);
                    p3Spokes =  obj.spokeHalfEdges(p3);
                
                    p1MeanSpokeLength = mean(selectionMetric(p1Spokes));
                    p3MeanSpokeLength = mean(selectionMetric(p3Spokes));
                
                    if p3MeanSpokeLength < p1MeanSpokeLength && length(p1Spokes)>4
                        temp = he4;
                        he4 = he1;
                        he1 = temp;
                    
                        temp = p3;
                        p3 = p1;
                        p1 = temp;
                    end
                    
                    % now p1 and p3 are finalized we need to find the two
                    % nodes that will have reduced spokes counts after the
                    % collapse. We'll come back to these at the end.
                    reducedSpokeVertices = obj.ends(obj.next([he1;he4]));
                    
                    % flag half edges of collapsed faces for removal
                    %------------------------------------------------------
                    fa = ceil(he1/3);
                    fb = ceil(he4/3);
                
                    removedHalfEdges = [3*fa-2, 3*fa-1, 3*fa,...
                                        3*fb-2, 3*fb-1, 3*fb];
                
                    if any(halfEdgeFlags(removedHalfEdges)==1)
                        error('half edge already deleted')
                    end
                
                    % if any vert points to a flagged half edge fix it
                    %------------------------------------------------------
                    obj.correctVertexHalfEdgeRelations(removedHalfEdges);
                
                    % change vertices of relevant half edges
                    obj.replaceHalfEdgeVertex(p1, p3);
                    
                    % average coordinate
                    if averageCollapse
                        if isFeatureBased
                            % weight towards features
                            w = angleVertices([p1;p3]);
                            w1 = max(w(1),0.01);
                            w2 = max(w(2),0.01);
                            obj.coordinates(p3,:) = (w1*obj.coordinates(p1,:)+w2*obj.coordinates(p3,:))/max(w1+w2,0.01);
                        else
                            obj.coordinates(p3,:) = 0.5*(obj.coordinates(p1,:)+obj.coordinates(p3,:));
                        end
                    end
                    
                    % fix pair connectivity around collapse faces
                    obj.fixCollapsedFacePairs(he1);

                    % make sure things are delaunay
                    [modifiedEdges,modifiedVertices] = obj.edgeFlipModifiedRegion(p3);

                    % update our metric for collpase
                    modifiedEdges = unique([modifiedEdges;obj.spokeHalfEdges(p3)]);
                    modifiedVertices = unique([modifiedVertices;obj.collapseModifiedVertices(p3)]);

                    modifiedFixedEdges = obj.isFixedHalfEdge(modifiedEdges);
                    modifiedEdges = modifiedEdges(modifiedFixedEdges==0);
                    modifiedMetric = obj.halfEdgeLengths(modifiedEdges);

                    if isFeatureBased
                        angleVertices(modifiedVertices) = obj.vertexAverageAngles(modifiedVertices);
                        modifiedAngles = (angleVertices(obj.ends(modifiedEdges)) +...
                            angleVertices(obj.ends(obj.pair(modifiedEdges))))/2;
                        modifiedMetric=modifiedMetric.*modifiedAngles;
                    end

                    selectionMetric(modifiedEdges) = modifiedMetric;

                    % check if we create any questionable cycles
                    % if so, push them to the front of the line
                    newLowValenceVertices = obj.updateLowValenceVertices(reducedSpokeVertices);
                    lowValenceCandidates = [lowValenceCandidates;...
                                            newLowValenceVertices];

                    % track removed entities
                    vertexFlags(p1) = 1;
                    halfEdgeFlags(removedHalfEdges) = 1;
                    isModifiedHalfEdge(modifiedEdges) = 1;
                    isModifiedHalfEdge(removedHalfEdges) = 2;
                    selectionMetric(removedHalfEdges)=nan;
                    
                    i = i+1;
                end
                
            end
            
            % clean things up
            obj.cleanDeletedEntities(vertexFlags,halfEdgeFlags) 
            obj.resetBulkProperties();
            obj.clearFields();
            
            progressMeter.close();
        end
        function smooth(obj,N,method,projectPoints,additionalWeights)
        % mesh smoothing
        %------------------------------------------------------------------
        % Inputs:
        %   N ------------------ number of smoothing iterations
        %   method ------------- ('uniform','cotangent',or 'area')
        %   projectPoints ------ project back onto original mesh when done?
        %   additionalWeights -- specify additional weight for vertices
        %------------------------------------------------------------------
             
            % process our inputs
            if nargin == 2
                method = obj.smoothOptions.method;
            end
            if nargin >= 3
                method = lower(method);
                if strcmp(method,"cotan") || strcmp(method,"cot")
                    method = "cotangent";
                end
                if ~strcmp(method,"cotangent")  && ~strcmp(method,"uniform") && ~strcmp(method,"area")
                    error('incorrect method specified must be: "cotangent", "uniform" or "area"')
                end
            end
            if nargin == 1 || nargin > 5
                error('incorrect number of inputs')
            end
            if nargin < 4
                projectPoints=obj.smoothOptions.reproject;
            end
            if nargin == 5 
                if max(additionalWeights) > 1 || min(additionalWeights) < 0
                    error('additional weights must between 1 and 0')
                end
                if length(additionalWeights) == obj.numFaces
                    additionalWeights = obj.vertexValuesFromFaces(additionalWeights);
                elseif length(additionalWeights) ~= obj.numVertices
                    error('additional smoothing weights must be for vertices or faces')
                end
            else
                additionalWeights = ones(obj.numVertices,1);
            end

            % if its a fixed vertex don't let it move
            additionalWeights(obj.isFixedVertex==1) = 0;
            

            % smooth rectilinear
            if obj.degree == 1

                progressMeter = progressbar('Smoothing Mesh              : ');

                for j = 1:N
                    progressMeter.update(100*j/N);

                    delta = zeros(obj.numVertices,3);
                    wsum = zeros(obj.numVertices,1);
                    
                    % for area weighting
                    A1 = ones(obj.numFaces,1);
                    A2 = ones(obj.numFaces,1);
                    A3 = ones(obj.numFaces,1);

                    tmpCoords = obj.coordinates;
                    
                    % vertex indices defining faces
                    p1 = obj.ends(1:3:end);
                    p2 = obj.ends(2:3:end);
                    p3 = obj.ends(3:3:end);

                    if strcmp(method,'area')
                        pointAreas = obj.vertexAreaVectors();
                        pointAreas = sqrt(pointAreas(:,1).^2+...
                                          pointAreas(:,2).^2+...
                                          pointAreas(:,3).^2);
                        A1 = pointAreas(p1);
                        A2 = pointAreas(p2);
                        A3 = pointAreas(p3);
                    end
                    
                    % coordinates of vertices
                    coord1 = tmpCoords(p1,:);
                    coord2 = tmpCoords(p2,:);
                    coord3 = tmpCoords(p3,:);
                    
                    % edge vectors
                    v1 = coord2-coord1;
                    v2 = coord3-coord2;
                    v3 = coord1-coord3;
                    
                    if strcmp(method,'cotangent')
                        n1 = v1./sqrt(v1(:,1).^2+v1(:,2).^2+v1(:,3).^2);
                        n2 = v2./sqrt(v2(:,1).^2+v2(:,2).^2+v2(:,3).^2);
                        n3 = v3./sqrt(v3(:,1).^2+v3(:,2).^2+v3(:,3).^2);
                    
                        w1 = cot(cos(-dot(n3,n1,2)));
                        w2 = cot(cos(-dot(n2,n1,2)));
                        w3 = cot(cos(-dot(n3,n2,2)));

                    else % uniform weights
                        w1 = ones(obj.numFaces,1);
                        w2 = ones(obj.numFaces,1);
                        w3 = ones(obj.numFaces,1);
                    end

                    % get our offsets
                    d1 = -w2.*v3.*A3 + w3.*v1.*A2;
                    d2 =  w1.*v2.*A3 - w3.*v1.*A1;
                    d3 =  w2.*v3.*A1 - w1.*v2.*A2;
                    
                    % need to loop to deal with repeat entries
                    for i = 1:obj.numFaces
                        
                        delta(p1(i),:) = delta(p1(i),:) + d1(i,:);
                        delta(p2(i),:) = delta(p2(i),:) + d2(i,:);
                        delta(p3(i),:) = delta(p3(i),:) + d3(i,:);

                        wsum(p1(i)) = wsum(p1(i)) + w3(i)*A3(i)+w2(i)*A2(i);
                        wsum(p2(i)) = wsum(p2(i)) + w3(i)*A3(i)+w1(i)*A1(i);
                        wsum(p3(i)) = wsum(p3(i)) + w2(i)*A1(i)+w1(i)*A2(i);
                    end
                    
                    tmpCoords = tmpCoords + 1./wsum.*delta.*additionalWeights;
                    
                end

                progressMeter.close();

                % project back onto original mesh if need be
                if projectPoints
                    oldCoords = obj.coordinates;
                    obj.coordinates = tmpCoords;
                    pointNormals = obj.vertexNormals;
                    obj.coordinates = oldCoords;
                    tmpCoords = obj.project(tmpCoords,pointNormals);
                end
                
                obj.coordinates = tmpCoords;
                
                % recalc bulk stuff
                obj.resetBulkProperties()
                obj.clearFields();
            else
                error("Smoothing is not set up for degree>1 meshes, you'll need to flatten it")
            end

        end
        function refine(obj,refineFaces)
        % converts each tri to 4. Allows subset of faces to be refined.
        %------------------------------------------------------------------
        %                 coord2                | generates 4 triangles 
        %                o                      | with he ordering as 
        %  he2        / / ^ ^                   | indicated.
        %  (coarse)  / /   \ \                   
        %           / /5   4\ \
        %          / /       \ \
        %         / v    6    \ \ he1 (coarse)
        %        / o<---------o  \
        %       / / \   11   ^ ^  \
        %      / /  ^\12  10//  \  \           type-3 refinement
        %     / /7  9\\    //2  1\  \
        %    / /      \\  //      \  \
        %   v v   8    \v/v   3    \  \           
        %    o--------->o---------->o coord1
        %     --------------------->
        %              he3 (coarse)
        %------------------------------------------------------------------       
            
            % make sure we're rectilinear and faces connectivity is correct
            obj.flatten();
            
            if nargin==1
            
                % get unique edge defs and halfEdge --> uniqueEdge transform
                [edges,~,halfEdgeEdges] = obj.edges();
            
                % add edge midpoints to coordinates
                tempVertices = [obj.coordinates;...
                               (obj.coordinates(edges(:,1),:)+...
                                obj.coordinates(edges(:,2),:))/2];
            
                % split into 4 faces
                tempFaces = zeros(4*obj.numFaces,3);
                tempFaces(1:4:end,1) = obj.numVertices + halfEdgeEdges(1:3:end,1);
                tempFaces(1:4:end,2) = obj.numVertices + halfEdgeEdges(3:3:end,1);
                tempFaces(1:4:end,3) = obj.faces(:,3);
                tempFaces(2:4:end,1) = obj.faces(:,1);
                tempFaces(2:4:end,2) = obj.numVertices + halfEdgeEdges(2:3:end,1);
                tempFaces(2:4:end,3) = obj.numVertices + halfEdgeEdges(1:3:end,1);
                tempFaces(3:4:end,1) = obj.faces(:,2);
                tempFaces(3:4:end,2) = obj.numVertices + halfEdgeEdges(3:3:end,1);
                tempFaces(3:4:end,3) = obj.numVertices + halfEdgeEdges(2:3:end,1);
                tempFaces(4:4:end,1) = obj.numVertices + halfEdgeEdges(1:3:end,1);
                tempFaces(4:4:end,2) = obj.numVertices + halfEdgeEdges(2:3:end,1);
                tempFaces(4:4:end,3) = obj.numVertices + halfEdgeEdges(3:3:end,1);
                
            elseif nargin == 2
                
                if max(refineFaces) > obj.numFaces ||...
                   any(refineFaces <= 0)  ||...
                   min(size(refineFaces)) ~= 1
                    error('something went wrong with refineFaces input')
                end
                   
                % flag half edges of refine faces and their pairs
                halfEdgeIsSplit = zeros(obj.numHalfEdges,1);
                halfEdgeIsSplit(3*refineFaces-2) = 1;
                halfEdgeIsSplit(3*refineFaces-1) = 1;
                halfEdgeIsSplit(3*refineFaces)   = 1;
                halfEdgeIsSplit(obj.pair(3*refineFaces-2)) = 1;
                halfEdgeIsSplit(obj.pair(3*refineFaces-1)) = 1;
                halfEdgeIsSplit(obj.pair(3*refineFaces))   = 1;
                
                % determine split type
                splitType = halfEdgeIsSplit(1:3:end) + ...
                            halfEdgeIsSplit(2:3:end) + ...
                            halfEdgeIsSplit(3:3:end);
                        
                numNewFaces = sum(splitType)+obj.numFaces;
                
                % unique edge that split
                [edges,edgeHalfEdges,halfEdgeEdges] = obj.edges();
                edgeIsSplit = zeros(obj.numEdges,1);
                for i = 1:obj.numEdges
                    if halfEdgeIsSplit(edgeHalfEdges(i,1)) || halfEdgeIsSplit(edgeHalfEdges(i,2))
                        edgeIsSplit(i)=1;
                    end
                end
                splitEdges = 1:obj.numEdges;
                splitEdges = splitEdges(edgeIsSplit==1);
                edges2SplitEdges = zeros(obj.numEdges,1);
                edges2SplitEdges(splitEdges) = 1:length(splitEdges);

                % add all midpoints to tempVertices we'll cull later
                tempVertices = [obj.coordinates;...
                               (obj.coordinates(edges(splitEdges,1),:)+...
                                obj.coordinates(edges(splitEdges,2),:))/2];
                            
                tempFaces = zeros(numNewFaces,3);
                newFaceIndex = 1;
                for i = 1:obj.numFaces

                    if splitType(i) == 0
                        tempFaces(newFaceIndex,:)=obj.faces(i,:);
                        newFaceIndex = newFaceIndex+1;
                        
                    elseif splitType(i) == 1
                        
                        if halfEdgeIsSplit(3*i-2) == 1
                            he1 = 3*i-2;
                        elseif halfEdgeIsSplit(3*i-1) == 1
                            he1 = 3*i-1;
                        elseif halfEdgeIsSplit(3*i) == 1
                            he1 = 3*i;
                        else
                            error('something went wrong here')
                        end
                       
                        he2 = obj.next(he1);
                        he3 = obj.next(he2);
                        
                        p1 = obj.ends(he1);
                        p2 = obj.ends(he2);
                        p3 = obj.ends(he3);
                        
                        uniqueEdgeId = halfEdgeEdges(he1);
                        uniqueSplitEdgeId = edges2SplitEdges(uniqueEdgeId);
                        p4 = obj.numVertices+uniqueSplitEdgeId;
                        
                        tempFaces(newFaceIndex,:)   = [p1,p2,p4];
                        tempFaces(newFaceIndex+1,:) = [p4,p2,p3];
                        newFaceIndex = newFaceIndex+2;
                        
                    elseif splitType(i) == 2
                        
                        if halfEdgeIsSplit(3*i-2) == 0
                            he1 = 3*i-2;
                        elseif halfEdgeIsSplit(3*i-1) == 0
                            he1 = 3*i-1;
                        elseif halfEdgeIsSplit(3*i) == 0
                            he1 = 3*i;
                        else
                            error('something went wrong here')
                        end
                       
                        he2 = obj.next(he1);
                        he3 = obj.next(he2);
                        
                        p1 = obj.ends(he1);
                        p2 = obj.ends(he2);
                        p3 = obj.ends(he3);
                        
                        uniqueEdgeId = halfEdgeEdges(he2);
                        uniqueSplitEdgeId = edges2SplitEdges(uniqueEdgeId);
                        p4 = obj.numVertices+uniqueSplitEdgeId;
                        
                        uniqueEdgeId = halfEdgeEdges(he3);
                        uniqueSplitEdgeId = edges2SplitEdges(uniqueEdgeId);
                        p5 = obj.numVertices+uniqueSplitEdgeId;
                        
                        tempFaces(newFaceIndex,:)   = [p1,p4,p3];
                        tempFaces(newFaceIndex+1,:) = [p3,p4,p5];
                        tempFaces(newFaceIndex+2,:) = [p5,p4,p2];
                        
                        newFaceIndex = newFaceIndex+3;
                        
                    elseif splitType(i) == 3

                        he1 = 3*i-2;
                        he2 = obj.next(he1);
                        he3 = obj.next(he2);
                        
                        p1 = obj.ends(he1);
                        p2 = obj.ends(he2);
                        p3 = obj.ends(he3);
                        
                        uniqueEdgeId = halfEdgeEdges(he1);
                        uniqueSplitEdgeId = edges2SplitEdges(uniqueEdgeId);
                        p4 = obj.numVertices+uniqueSplitEdgeId;
                        
                        uniqueEdgeId = halfEdgeEdges(he2);
                        uniqueSplitEdgeId = edges2SplitEdges(uniqueEdgeId);
                        p5 = obj.numVertices+uniqueSplitEdgeId;
                        
                        uniqueEdgeId = halfEdgeEdges(he3);
                        uniqueSplitEdgeId = edges2SplitEdges(uniqueEdgeId);
                        p6 = obj.numVertices+uniqueSplitEdgeId;
                        
                        tempFaces(newFaceIndex,:)   = [p4,p5,p6];
                        tempFaces(newFaceIndex+1,:) = [p4,p1,p5];
                        tempFaces(newFaceIndex+2,:) = [p2,p6,p5];
                        tempFaces(newFaceIndex+3,:) = [p4,p6,p3];
                        
                        newFaceIndex = newFaceIndex+4;
                    end  
                end               
            end
            obj.initializeFromFaceData(tempVertices,tempFaces);
        end
        function refineFeatures(obj,fineMesh, numIters, alphaLimit, maxLevel)
        % wrapper method to allow for easier feature based-refinement
        %------------------------------------------------------------------
        % This attempts to refine the mesh in areas with large dihedral
        % angles. Because we're starting w/ a course mesh and refining,
        % it won't necessarily get everything. The refined mesh is then
        % projected onto another ideally more refined mesh.
        %------------------------------------------------------------------
        % Inputs:
        %   fineMesh --- mesh new pts are projected onto 
        %   numIters --- number of refinement iterations
        %   alphaLimit - criterion for refinement, smoothed avg angle 
        %                between neighboring faces (7-10 reccomended) 
        %   maxLevel --- max refinement level for a single face. 
        %------------------------------------------------------------------
            disp('--------------------------------------------------')
            disp(['Refining features, initial numFaces: ',num2str(obj.numFaces)])
            if nargin == 2
                disp('nargin 1')
                alphaLimit = 9;
                numIters = 1;
                maxLevel = 1;
            elseif nargin == 3
                disp('nargin 2')
                alphaLimit = 9;
                maxLevel = 4;
            elseif nargin == 4
                disp('nargin 3')
                maxLevel = 4;
            elseif nargin > 5
                error('incorrect number of inputs')
            end
            
            
            obj.flatten;
            
            refineLevels = zeros(obj.numFaces,1);
            for i =1:numIters
                alpha = obj.faceAverageMeshAngles();
                alphav = obj.vertexValuesFromFaces(alpha);
                alpha = obj.faceValuesFromVertices(alphav);
                alphav = obj.vertexValuesFromFaces(alpha);
                alpha = obj.faceValuesFromVertices(alphav);
                
                faceIndices = 1:obj.numFaces;
                candidates = faceIndices(alpha>alphaLimit);
          
                candidateLevels = refineLevels(candidates);
                levelLimitedCandidates = candidates(candidateLevels<maxLevel);
          
                if ~isempty(levelLimitedCandidates)
                    [obj,refineLevels] = obj.refineAndTrackLevel(levelLimitedCandidates,refineLevels);
                end
                obj.coordinates = fineMesh.project(obj.coordinates,obj.vertexNormals());

                disp(['   refinementIter: ',num2str(i),'   numFaces: ',num2str(obj.numFaces)])
            end
            disp('--------------------------------------------------')
            disp(' ')
            
            obj.clearFields();
            obj.resolution = obj.calculateResolution();
        end
        
        function projectOnTo(obj, mesh)
        % projects onto another mesh
        %------------------------------------------------------------------
        % wrapper for the project points method that makes UI easier
        %------------------------------------------------------------------
        % Inputs:
        %   mesh ---------- SurfaceMesh object, must be of similar topology
        %                   for the projection routine to work.
        %------------------------------------------------------------------
            
            obj.coordinates = mesh.project(obj.coordinates,obj.nodeNormals());
            
            obj.volume = obj.calculateVolume;
            obj.surfaceArea = obj.calculateSurfaceArea;
            obj.centroid = obj.calculateCentroid;
            obj.clearFields();
                
        end
        function projectedPoints = project(obj, points, pointNormals)
        % projects set of near-surface points onto the mesh.
        %------------------------------------------------------------------
        % for each point, begin by testing the nearest n-facets in the
        % local octant. If no appropriate facet is found the next nearest
        % n-facets is searched and this is repeated several times. If still
        % no facet for projection is found reset the search domain to all
        % facets (i.e. all 8 octants) since we're probably dealing with a
        % octant boundary case. 
        %------------------------------------------------------------------
        % Inputs:
        %   points -------- points to be projected onto the mesh
        %   pointNormals -- projection direction for each point 
        %------------------------------------------------------------------
        % Outputs:
        %   projectPoints - new post-projection coordinates
        %------------------------------------------------------------------
            threshold = 1e-4*obj.resolution;
            numProjPoints = size(points,1);
            numCandidates = 12;
            maxProjectionDistance = min(max(obj.coordinates,[],1) - ...
                                        min(obj.coordinates,[],1))/3;
            projectedPoints = zeros(numProjPoints,3);
            
            % if no projection direction specified use radial
            if nargin==2
                pointNormals = points...
                    ./sqrt(points(:,1).^2+points(:,2).^2+points(:,3).^2);
            end

            % We're going to project onto this mesh providing a theshold
            % for whether a line intersects a triangle
            faceCentroids = obj.faceCentroids();
            faceNormals = obj.faceNormals();
            faceCoordinates1 = obj.coordinates(obj.faces(:,1),:);
            faceCoordinates2 = obj.coordinates(obj.faces(:,2),:);
            faceCoordinates3 = obj.coordinates(obj.faces(:,3),:);
            
            % parse receiving faces into respective octants
            %--------------------------------------------------------------
            Nsteps = 2;
            octantFaceIndices{8} = [];
            octantFaceCentroids{8} = [];
            
            indices = 1:obj.numFaces;
            
            % create some overlap between octants so that we don't really
            % have to do special stuff for bc cases.
            offset = ones(8,3);
            offset(2:2:end,1) = -1;
            offset([3,4,7,8],2) = -1;
            offset(5:end,3) = -1;
            for i = 1:8
                nodeOctants = 1 +(faceCentroids(:,1) > (2*obj.resolution * offset(i,1)))...
                        + Nsteps*(faceCentroids(:,2) > (2*obj.resolution * offset(i,2)))...
                      + Nsteps^2*(faceCentroids(:,3) > (2*obj.resolution * offset(i,3)));
                  
                octantFaceIndices{i}=indices(nodeOctants==i);
                octantFaceCentroids{i}=faceCentroids(octantFaceIndices{i},:); 
            end
            progressMeter = progressbar('Projecting Points           : ');
            
            % Project
            %--------------------------------------------------------------
            for i=1:numProjPoints

                progressMeter.update(100*i/numProjPoints);

                ni = pointNormals(i,:);    
                coordi = points(i,:); 
                
                octIndex = 1 +(coordi(1,1)>0)+...
                            2*(coordi(1,2)>0)+...
                            4*(coordi(1,3)>0);
                        
                dist = octantFaceCentroids{octIndex}-coordi;             
                dist = dist(:,1).^2 + dist(:,2).^2 + dist(:,3).^2;
                
                allFaces=false;   
                projFaceFound=false;
                numCandidatesi = numCandidates;
                
                while ~projFaceFound
                    
                    % select the closest n faces
                    [~,candidates] = mink(dist,numCandidatesi); 
                    if ~allFaces
                        candidates = octantFaceIndices{octIndex}(candidates);
                    end
                    
                    j0 = numCandidatesi-numCandidates+1;
                    j1 = numCandidatesi;
                    
                    for j = j0:j1
                    
                        % things for face-j
                        nf = faceNormals(candidates(j),:);      
                        coord1 = faceCoordinates1(candidates(j),:); 
                        coord2 = faceCoordinates2(candidates(j),:);
                        coord3 = faceCoordinates3(candidates(j),:);

                        % intersection algo -
                        %--------------------------------------------------
                        % there was some brutal overhead using the func
                        
                        distance = - (nf(1)*(coordi(1)-coord1(1))+...
                                      nf(2)*(coordi(2)-coord1(2))+...
                                      nf(3)*(coordi(3)-coord1(3)))/...
                                  (nf(1)*ni(1)+nf(2)*ni(2)+nf(3)*ni(3));
                        
                        % point of intersection
                        intersectionCoord = coordi + distance*ni;
                        
                        v1 = coord2-coord1;            % edge vector 1
                        v2 = coord3-coord1;            % edge vector 2
                        v0 = intersectionCoord-coord1; % vector to intersection point
                        
                        % calc parametric coordinates of intersection point
                        dotv2v2 = v2(1)*v2(1) + v2(2)*v2(2) + v2(3)*v2(3);
                        dotv1v2 = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3);
                        dotv1v1 = v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3);
                        dotv0v2 = v0(1)*v2(1) + v0(2)*v2(2) + v0(3)*v2(3);
                        dotv0v1 = v0(1)*v1(1) + v0(2)*v1(2) + v0(3)*v1(3);
                        denom = dotv1v2 ^ 2 - dotv1v1 * dotv2v2;
                        
                        xj = (dotv1v2 * dotv0v2 - dotv2v2 * dotv0v1)/denom;
                        yj = (dotv1v2 * dotv0v1 - dotv1v1 * dotv0v2)/denom;
                        
                        if        xj > -threshold && ...
                                  yj > -threshold && ...
                           1-(xj+yj) > -threshold && ...
                           distance > - maxProjectionDistance
                       
                            projectedPoints(i,1:3) = intersectionCoord;
                            projFaceFound = true;
                            break
                        end
                        
                        % expand search if we fail
                        %--------------------------------------------------
                        if ~projFaceFound && j == j1
                            numCandidatesi = numCandidatesi+numCandidates;
                            
                            % expand search if we fail a lot
                            %----------------------------------------------
                            if (numCandidatesi > 4.5 * numCandidates) && allFaces==false
                                %disp('failed ... seaching all...')
                                allFaces=true;
                                numCandidatesi = numCandidates;
                                dist = faceCentroids-coordi;             
                                dist = dist(:,1).^2 + dist(:,2).^2 + dist(:,3).^2;
                            end
                            
                            % if we keep failing reduce tolerance, its
                            % probably a collocated node or something like
                            % that 
                            if (numCandidatesi > 4.5 * numCandidates) && ...
                                allFaces==true && ...
                                threshold < 0.11*obj.resolution
                                
                                threshold = min(threshold*10,0.11*obj.resolution);
                                numCandidatesi = numCandidates;
                                %dist = faceCentroids-coordi;             
                                %dist = dist(:,1).^2 + dist(:,2).^2 + dist(:,3).^2;
                            end
                        end
                        if allFaces == true && numCandidatesi > (obj.numFaces-numCandidates)
                            error('candidate face for projection could not be found')
                        end
                    end
                    
                end
            end
            progressMeter.close();    
        end
        function projectedPoints = projectRobust(obj, points, pointNormals)
        % projects set of near-surface points onto the mesh.
        %------------------------------------------------------------------
        % for each point, begin by testing the nearest n-facets in the
        % local octant. If no appropriate facet is found the next nearest
        % n-facets is searched and this is repeated several times. If still
        % no facet for projection is found reset the search domain to all
        % facets (i.e. all 8 octants) since we're probably dealing with a
        % octant boundary case. 
        %------------------------------------------------------------------
        % Inputs:
        %   points -------- points to be projected onto the mesh
        %   pointNormals -- projection direction for each point 
        %------------------------------------------------------------------
        % Outputs:
        %   projectPoints - new post-projection coordinates
        %------------------------------------------------------------------
            
            threshold = 1e-4;
            numProjPoints = size(points,1);
            projectedPoints = zeros(numProjPoints,3);
            
            % if no projection direction specified use radial
            if nargin==2
                pointNormals = points...
                    ./sqrt(points(:,1).^2+points(:,2).^2+points(:,3).^2);
            end

            if obj.degree>1
                warning('projection not set up for degree>1 ... flattening')
                obj.flatten();
            end
            
            % We're going to project onto this mesh providing a theshold
            % for whether a line intersects a triangle
            faceNormals = obj.faceNormals();
            faceCoordinates1 = obj.coordinates(obj.faces(:,1),:);
            faceCoordinates2 = obj.coordinates(obj.faces(:,2),:);
            faceCoordinates3 = obj.coordinates(obj.faces(:,3),:);

             progressMeter = progressbar('Projecting Points           : ');

            % Project
            %--------------------------------------------------------------
            for i=1:numProjPoints

                progressMeter.update(100*i/numProjPoints);

                ni = pointNormals(i,:);   
                coordi = points(i,:);
                        
                distancej = -obj.resolution*1e10;
                
                for j = 1:obj.numFaces
                    
                    % things for face-j
                    nf = faceNormals(j,:);
                    coord1 = faceCoordinates1(j,:);
                    coord2 = faceCoordinates2(j,:);
                    coord3 = faceCoordinates3(j,:);
                    
                    % intersection algo -
                    %--------------------------------------------------
                    % there was some brutal overhead using the func
                    
                    distance = - (nf(1)*(coordi(1)-coord1(1))+...
                                  nf(2)*(coordi(2)-coord1(2))+...
                                  nf(3)*(coordi(3)-coord1(3)))/...
                               (nf(1)*ni(1)+nf(2)*ni(2)+nf(3)*ni(3));
                    
                    % point of intersection
                    intersectionCoord = coordi + distance*ni;
                    
                    v1 = coord2-coord1;            % edge vector 1
                    v2 = coord3-coord1;            % edge vector 2
                    v0 = intersectionCoord-coord1; % vector to intersection point
                    
                    % calc parametric coordinates of intersection point
                    dotv2v2 = v2(1)*v2(1) + v2(2)*v2(2) + v2(3)*v2(3);
                    dotv1v2 = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3);
                    dotv1v1 = v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3);
                    dotv0v2 = v0(1)*v2(1) + v0(2)*v2(2) + v0(3)*v2(3);
                    dotv0v1 = v0(1)*v1(1) + v0(2)*v1(2) + v0(3)*v1(3);
                    denom = dotv1v2 ^ 2 - dotv1v1 * dotv2v2;
                    
                    xj = (dotv1v2 * dotv0v2 - dotv2v2 * dotv0v1)/denom;
                    yj = (dotv1v2 * dotv0v1 - dotv1v1 * dotv0v2)/denom;
                    
                    if      xj > -threshold && ...
                            yj > -threshold && ...
                     1-(xj+yj) > -threshold && ...
                      distance > distancej
                      if i ==1
                         disp([xj,yj,1-xj-yj]) 
                      end
                        projectedPoints(i,1:3) = intersectionCoord;
                        distancej = distance;
                    end
                    
                end
            end
            progressMeter.close();  
        end
        
        function  isValid(obj)
        % alias for check valid 
        %------------------------------------------------------------------
            obj.checkValid();
        end  
        function  checkValid(obj)
        % performs some simple tests to see if our mesh is invalid
        %------------------------------------------------------------------
        % If this passes is doesn't necessarily mean the mesh is valid, but
        % it will tell you if somethings went wrong
        %------------------------------------------------------------------
        
            % test pairs
            pairs = obj.pair;
            indices = linspace(1,obj.numHalfEdges,obj.numHalfEdges)';
            if ~all(obj.pair(pairs) == indices)
                error('invalid mesh topology: pairs of pairs invalid')
            end
            
            
            % test face loops
            nexti = indices(3:3:obj.numHalfEdges);
            for i = 1:3
                nexti = obj.next(nexti);
                if ~all(nexti == indices(i:3:obj.numHalfEdges))
                    error(['invalid mesh topology: next invalid ', i])
                end
            end
            
            
            % test vertex loops
            maxIters = 20;          
            startingEdge = obj.vertexHalfEdges;
            nexti = obj.next(startingEdge);
            pairi = obj.pair(nexti);
            iters = 0;
            while  any(pairi ~= startingEdge)
                iters = iters+1;
                keepLoops = pairi ~= startingEdge;
                pairi = pairi(keepLoops); 
                startingEdge = startingEdge(keepLoops); 

                nexti = obj.next(pairi);
                pairi = obj.pair(nexti);
                if iters==maxIters 
                    break
                end
                    
            end
            if iters == maxIters
                error('invalid mesh topology: infinite vertex loop')
            end
            
            % test vertices point to right half edges
            indices = linspace(1,obj.numVertices,obj.numVertices)';
            newIndices = obj.ends(obj.vertexHalfEdges);
            if ~all(newIndices == indices)
                error(['invalid mesh topology: vertex not pointing to correct half edge vertex numbers:',num2str(indices(newIndices ~= indices))])
            end
            
        end  
        
        function [] = writeOBJ(obj,fileName,precision)
        % writes mesh to obj file
        %------------------------------------------------------------------
        % Inputs:
        %   fileName ----- name of obj file (including extension)
        %   AsteroidName - name of asteroid
        %   DataSource --- where original shape model was retrieved
        %   precision ---- precision of the write (number of sig-figs)
        %==================================================================
            
            % handle optional inputs
            if nargin<2
                fileName = strcat('unnamedMesh_',num2str(obj.numFaces),'.obj');
            end
            if nargin <3
                precision = 8;
            end
            
            
            
            fileID = fopen(fileName,'w');
            
            % Header
            %--------------------------------------------------------------------------
            fprintf(fileID, '##########################################################\n');
            fprintf(fileID, '#\n');
            fprintf(fileID, ['# Vertices: ',num2str(obj.numVertices),'\n']);
            fprintf(fileID, ['# Faces: ',num2str(obj.numFaces),'\n']);
            fprintf(fileID, '#\n');
            fprintf(fileID, '##########################################################\n');
            
            
            % Write Data
            %--------------------------------------------------------------
            
            for i = 1:obj.numVertices
                fprintf(fileID, ['v ',num2str(obj.coordinates(i,1),precision)...
                    ,' ',num2str(obj.coordinates(i,2),precision)...
                    ,' ',num2str(obj.coordinates(i,3),precision),'\n'] );
            end
            
            
            for i = 1:obj.numFaces
                fprintf(fileID,'f');
                for j =1:length(obj.faces(i,:))
                    fprintf(fileID,[' ',num2str(obj.faces(i,j))]);
                end
                fprintf(fileID,'\n');
            end
            
            fclose(fileID);
        end
        function [] = writeVTK(obj,fileName,precision)
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
            precisionStr = num2str(precision);
            
            %vtk flags for rectilinear and quadratic triangles
            if obj.degree == 1
                cellType = 5*ones(obj.numFaces,1);
                numLocalNodes = 3;
            elseif obj.degree == 2 
                cellType = 22*ones(obj.numFaces,1);
                numLocalNodes = 6;
            else 
                error('Degree greater than 2 not supported')
            end
            
            fileID = fopen(fileName,'w');
            
            % Header
            %--------------------------------------------------------------
            fprintf(fileID,'# vtk DataFile Version 3.0\n');
            fprintf(fileID,'vtk output\n');
            fprintf(fileID,'ASCII\n');
            fprintf(fileID,'DATASET UNSTRUCTURED_GRID\n');
            
            % formation specifications
            %--------------------------------------------------------------
            formatSpecVectorFloat = ['%.',precisionStr,'g %.',precisionStr,'g %.',precisionStr,'g\n'];
            formatSpecScalarFloat = ['%.',precisionStr,'g\n'];
            formatSpecFacesLinear = '%.0f %.0f %.0f %.0f\n';
            formatSpecFacesQuadratic = '%.0f %.0f %.0f %.0f %.0f %.0f %.0f\n';
            
            % Write Node Coordinates
            %--------------------------------------------------------------
            fprintf(fileID,['POINTS ',num2str(obj.numNodes),' float\n']);
            fprintf(fileID,formatSpecVectorFloat,obj.coordinates');
            
            % Write Cells
            %--------------------------------------------------------------
            fprintf(fileID,['CELLS ',num2str(obj.numFaces),' ',num2str(obj.numFaces*(numLocalNodes+1)),'\n']);
            if obj.degree == 1
                facesOut = [3*ones(obj.numFaces,1),obj.faces-1]';
                fprintf(fileID,formatSpecFacesLinear,facesOut);
            elseif obj.degree == 2
                facesOut = [6*ones(obj.numFaces,1),obj.faces(:,[1,3,6,2,5,4])-1]';
                fprintf(fileID,formatSpecFacesQuadratic,facesOut);
                
            end
            
            % Write CellTypes
            %--------------------------------------------------------------
            fprintf(fileID,['CELL_TYPES ',num2str(obj.numFaces),'\n']);
            fprintf(fileID,'%.0f\n',cellType);
            
            % Write Face Normals
            %--------------------------------------------------------------
            fprintf(fileID,['CELL_DATA ',num2str(obj.numFaces),'\n']);
           
            vtkDataHeader = 'VECTORS faceNormals float\n';
            fprintf(fileID,vtkDataHeader);
            normals = obj.faceNormals();
            fprintf(fileID,formatSpecVectorFloat,normals');
            
            % Write Stored Face Fields
            %--------------------------------------------------------------
            for i = 1:size(obj.faceFields,2)
                if size(obj.faceFields{i}.data,2)==1
                    vtkDataHeader = ['SCALARS ',obj.faceFields{i}.name,' float 1\nLOOKUP_TABLE default\n'];
                    fprintf(fileID,vtkDataHeader);
                    fprintf(fileID,formatSpecScalarFloat,obj.faceFields{i}.data);
                elseif size(obj.faceFields{i}.data,2)==3
                    vtkDataHeader = ['VECTORS ',obj.faceFields{i}.name,' float\n'];
                    fprintf(fileID,vtkDataHeader);
                    fprintf(fileID,formatSpecVectorFloat,obj.faceFields{i}.data);
                else
                    warning(['skipping faceData{',num2str(i),'} not scalar or vector'])
                end
                
                if size(obj.faceFields{i}.data,1)~=obj.numFaces
                    error(['faceData{',num2str(i),'} inconsistent w/ number of mesh faces'])
                end
%                 for j = 1:obj.numFaces
%                     outputi = [num2str(obj.faceFields{i}.data(j,:),precision),'\n'];
%                     fprintf(fileID,outputi);
%                 end
                
            end
            
            %Write Stored Node Fields
            %--------------------------------------------------------------
            fprintf(fileID,['POINT_DATA ',num2str(obj.numNodes),'\n']);
            for i = 1:size(obj.nodeFields,2)
                if size(obj.nodeFields{i}.data,2)==1
                    vtkDataHeader = ['SCALARS ',obj.nodeFields{i}.name,' float 1\nLOOKUP_TABLE default\n'];
                    fprintf(fileID,vtkDataHeader);
                    fprintf(fileID,formatSpecScalarFloat,obj.nodeFields{i}.data);
                elseif size(obj.nodeFields{i}.data,2)==3
                    vtkDataHeader = ['VECTORS ',obj.nodeFields{i}.name,' float\n'];
                    fprintf(fileID,vtkDataHeader);
                    fprintf(fileID,formatSpecVectorFloat,obj.nodeFields{i}.data);
                else
                    warning(['skipping faceData{',num2str(i),'} not scalar or vector'])
                end
                if size(obj.nodeFields{i}.data,1)~=obj.numNodes
                    error(['faceData{',num2str(i),'} inconsistent w/ number of mesh faces'])
                end
%                 for j = 1:obj.numNodes
%                     outputi = [num2str(obj.nodeFields{i}.data(j,:),precision),'\n'];
%                     fprintf(fileID,outputi);
%                 end
                
            end
            fclose(fileID);
            
        end
        
        function faceShading = shadowedFaces(obj, lightDirection)
        % finds shadowed faces based on lighting direction
        %------------------------------------------------------------------
        % uses a 2d background mesh to search for possible shadowing. The 
        % search is limited to faces who's normals have a negative dot 
        % product with the light direction.
        %
        %------------------------------------------------------------------
        % Inputs:
        %   lightDirection - direction light comes in at
        %------------------------------------------------------------------
        % Outputs:
        %   faceShading ----- scale 0 - 1 giving face shading
        %------------------------------------------------------------------
         % reused x1 and x2 badddness
            if obj.degree>1
                error('shading only implemented for rectilinear grids')
            end
            if ~size(lightDirection,2)==3 || ~size(lightDirection,1)==1
                error('lightDirection must be 1x3 vector')
            end
            lightDirection = lightDirection/norm(lightDirection);
            
            threshold = 1e-4;
            
            % initialize faces on backside as shadowed
            faceShading = obj.faceNormals()*lightDirection';
            faceShading(faceShading<0) = 0;
            faceShading(faceShading>0) = 1;
            
            
            % some things we'll need and easy aliases
            faceVertices1 = obj.coordinates(obj.faces(:,1),:); 
            faceVertices2 = obj.coordinates(obj.faces(:,2),:); 
            faceVertices3 = obj.coordinates(obj.faces(:,3),:); 
            faceCentroids = obj.faceCentroids();
            faceNormals = obj.faceNormals();
      
            % construct orthnormal basis for lighting plane
            %--------------------------------------------------------------
            z0 = -lightDirection/norm(lightDirection);
            x0 = [-z0(2), z0(1), 0;...
                  -z0(3), 0,     z0(1);...
                   0,    -z0(3), z0(2)];
            [~,i] = max(dot(x0,x0,2));
            x0 = x0(i,:);
            x0 = x0/norm(x0);
            
            y0 = cross(z0,x0);
            y0 = y0/norm(y0);
            
            % coordinates in new basis
            x0 = faceCentroids*x0'; % parametric coord 1
            y0 = faceCentroids*y0'; % parametric coord 2
            z0 = faceCentroids*z0'; % parametric coord 3 (height off plane)

            % construct background mesh, group face-centroids by mesh-cell
            %--------------------------------------------------------------
            Nsteps = ceil(min(max(x0)-min(x0),max(y0)-min(y0))/obj.resolution/2.0);

            x1 = x0-min(x0)+obj.resolution;
            x1 = x1/(max(x1)+obj.resolution)*Nsteps;
            
            y1 = y0-min(y0)+obj.resolution;
            y1 = y1/(max(y1)+obj.resolution)*Nsteps;
            
            faceIndices = 1:obj.numFaces;
          
            meshIndex = Nsteps*floor(y1)+ceil(x1);
            
            cellNodes{Nsteps^2}=[];
            for i = 1:Nsteps^2
                cellNodes{i} = faceIndices(meshIndex==i);
            end
            tempCellNodes = cellNodes;
            
            for i = 1:Nsteps^2
                if i > Nsteps
                    if rem(i,Nsteps)~=1
                        cellNodes{i} = [cellNodes{i},tempCellNodes{i-Nsteps-1}];
                    end
                    cellNodes{i} = [cellNodes{i},tempCellNodes{i-Nsteps}];
                    if rem(i,Nsteps)~=0
                        cellNodes{i} = [cellNodes{i},tempCellNodes{i-Nsteps+1}];
                    end
                end
                if i <= Nsteps^2-Nsteps
                    if rem(i,Nsteps)~=1
                        cellNodes{i} = [cellNodes{i},tempCellNodes{i+Nsteps-1}];
                    end
                    cellNodes{i} = [cellNodes{i},tempCellNodes{i+Nsteps}];
                    if rem(i,Nsteps)~=0
                        cellNodes{i} = [cellNodes{i},tempCellNodes{i+Nsteps+1}];
                    end
                end
                if rem(i,Nsteps)~=1
                    cellNodes{i} = [cellNodes{i},tempCellNodes{i-1}];
                end
                if rem(i,Nsteps)~=0
                    cellNodes{i} = [cellNodes{i},tempCellNodes{i+1}];
                end
            end
            
         
            % check for shadowing
            %--------------------------------------------------------------
            for i=1:obj.numFaces
                
                if faceShading(i)==0
                    
                    coordi = faceCentroids(i,:);
                    ni  = lightDirection;
                    
                    candidates = cellNodes{meshIndex(i)};
                    numCandidates = length(candidates);
                    
                    candidateNormals   = faceNormals(candidates,:);        
                    candidateVertices1 = faceVertices1(candidates,:); 
                    candidateVertices2 = faceVertices2(candidates,:); 
                    candidateVertices3 = faceVertices3(candidates,:); 
                    
                    for j = 2:numCandidates
                        
                        nf     = candidateNormals(j,:);    
                        coord1 = candidateVertices1(j,:);  
                        coord2 = candidateVertices2(j,:); 
                        coord3 = candidateVertices3(j,:); 
                        
                        % intersection algo - "inlined"
                        %--------------------------------------------------
                        distance = - (nf(1)*(coordi(1)-coord1(1))+...
                            nf(2)*(coordi(2)-coord1(2))+...
                            nf(3)*(coordi(3)-coord1(3)))/...
                            (nf(1)*ni(1)+nf(2)*ni(2)+nf(3)*ni(3));
                        
                        % point of intersection
                        intersectionCoord = coordi + distance*ni;
                        
                        v1 = coord2-coord1;            % edge vector 1
                        v2 = coord3-coord1;            % edge vector 2
                        v0 = intersectionCoord-coord1; % vector to intersection point
                        
                        % calc parametric coordinates of intersection point
                        dotv2v2 = v2(1)*v2(1) + v2(2)*v2(2) + v2(3)*v2(3);
                        dotv1v2 = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3);
                        dotv1v1 = v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3);
                        dotv0v2  =  v0(1)*v2(1) + v0(2)*v2(2) +  v0(3)*v2(3);
                        dotv0v1  =  v0(1)*v1(1) + v0(2)*v1(2) +  v0(3)*v1(3);
                        denom = dotv1v2 ^ 2 - dotv1v1 * dotv2v2;
                        
                        xj = (dotv1v2 * dotv0v2 - dotv2v2 * dotv0v1)/denom;
                        yj = (dotv1v2 * dotv0v1 - dotv1v1 * dotv0v2)/denom;
                        
                        if       xj  > -threshold && ...
                                 yj  > -threshold && ...
                           1-(xj+yj) > -threshold && ...
                           z0(candidates(j))>z0(i)
                       
                            faceShading(i)=1;
                            break
                            
                        end
                        
                        
                    end
                end
            end 
        end
        function plot(obj,varargin)
        % plot utility for the surfaceMesh
        %------------------------------------------------------------------
        % Inputs:
        %   plotType -- how are we doing this? 
        %               'wireframe'
        %------------------------------------------------------------------
            hold on
        
            % defaults
            %--------------------------------------------------------------
            edgeColor   = 'k';
            plotType    = 'wireframe';
            faceColor   = [0.75,0.75,0.75];
            faceAlpha   = 1;
            isBounded   = 0;
            customColor = 0;
            lineStyle   = '-';
            color0 = [0,0,1];
            color1 = [1,1,0];
            lighting = 1;
            shadowingOn = 1;
            
            % parse keywords
            %--------------------------------------------------------------
            if length(varargin)==1
                plotType = varargin{1};
            elseif rem(length(varargin),2)~=0
                error('each arg requires a keyword beforehand')
            end
            
            for i = 1:length(varargin)/2
                keyword = lower(varargin{2*i-1});
                if strcmp(keyword,'plottype')
                    plotType = lower(varargin{2*i});
                elseif strcmp(keyword,'edgecolor')
                    edgeColor = varargin{2*i};
                elseif strcmp(keyword,'facecolor')
                    faceColor = varargin{2*i};
                elseif strcmp(keyword,'facealpha')
                    faceAlpha = varargin{2*i};
                elseif strcmp(keyword,'bounds')
                    bounds = varargin{2*i};
                    isBounded = 1;
                elseif strcmp(keyword,'maxcolor')
                    color1 = varargin{2*i};
                elseif strcmp(keyword,'mincolor')
                    color0 = varargin{2*i};
                elseif strcmp(keyword,'lighting')
                    lighting = varargin{2*i};
                elseif strcmp(keyword,'linestyle')
                    lineStyle = varargin{2*i};
                elseif strcmp(keyword,'shadow') || strcmp(keyword,'shadowing')
                    if strcmp(lower(varargin{2*i}),'off') || varargin{2*i}(1)==0
                        shadowingOn=0;
                    end
                end
            end
            
            % process entries
            %--------------------------------------------------------------
            
            % allow for transposes
            if size(faceColor,2) == obj.numFaces
                faceColor=faceColor';
            end
            
            if size(lighting,2) == obj.numFaces
                lighting=lighting';
            end
            
            % convert matlab color keys to rgb
            if isa(faceColor,'char')
                switch faceColor
                    case 'k'
                        faceColor = [0,0,0];
                    case 'w'
                        faceColor = [1,1,1];
                    case 'm'
                        faceColor = [1,0,1];
                    case 'r'
                        faceColor = [1,0,0];
                    case 'y'
                        faceColor = [1,1,0];
                    case 'g'
                        faceColor = [0,1,0];
                    case 'c'
                        faceColor = [0,1,1];
                    case 'b'
                        faceColor = [0,0,1];
                    otherwise
                        error('invalid faceColor char entry, must be matlab standard color char')
                end
            end
            
         
            % make sure were getting a valid faceColor
            if strcmp(plotType, 'userdefined') &&...
                 ~(size(faceColor,2) == 3 && size(faceColor,1) == 1) &&...
                 ~(size(faceColor,2) == 3 && size(faceColor,1) == obj.numFaces) &&...
                 ~(size(faceColor,2) == 1 && size(faceColor,1) == obj.numFaces)
                error('incorrect keyword input for faceColor. Valid options: char, 1x3 vector, Nfx1 vector, Nfx3 vector')
            end
            
            % if a direction is specified for illumination
            if size(lighting,1)==1 && size(lighting,2)==3
                lighting = lighting/norm(lighting);
                lightingVector = lighting;
                lighting = obj.faceNormals*lighting';
                
                
                if shadowingOn == 1  % raytrace
                    faceShading = obj.shadowedFaces(lightingVector);
                    
                    vertexShading=zeros(obj.numVertices,1);
                    vertexFrequency = zeros(obj.numVertices,1);
                    for i = 1:obj.numFaces
                        vertexShading(obj.faces(i,:)) = vertexShading(obj.faces(i,:))+faceShading(i); 
                        vertexFrequency(obj.faces(i,:)) = vertexFrequency(obj.faces(i,:))+1;
                    end
                    vertexShading = vertexShading./vertexFrequency;
                    faceShading = 1/3*(vertexShading(obj.faces(:,1)) + vertexShading(obj.faces(:,2)) + vertexShading(obj.faces(:,3)));
                    lighting=lighting.*(1-faceShading);
                    lighting(lighting>0) = 0;
                    lighting = lighting * -1;
                    
                else % simple normal orientation shading
                    lighting(lighting>0) = 0;
                    lighting = lighting * -1;
                end
              
            elseif (size(lighting,1) == 1 && size(lighting,2) == 1)
                lighting = lighting*ones(obj.numFaces,1);
            end
             
            % Now plot the stuff
            %--------------------------------------------------------------
            pts = obj.coordinates(obj.ends,:);
   
            % wireFrame
            if strcmp(plotType, 'wireframe')
                for i =1:obj.numHalfEdges
                    if obj.pair(i)>i
                        indices = [i,obj.pair(i)];
                        plot3(pts(indices,1),pts(indices,2),pts(indices,3),'k','color',edgeColor,'lineStyle',lineStyle)
                    end
                end
                
            % singleFaceColorPlots   
            elseif strcmp(plotType, 'userdefined') &&...
                    (size(faceColor,2) == 3 && size(faceColor,1) == 1)
                
                for i =1:obj.numFaces
                    colori = faceColor*lighting(i);
                    indices = [3*i-2,3*i-1,3*i,3*i-2];
                    fill3(pts(indices,1),pts(indices,2),pts(indices,3),...
                        colori,'EdgeColor',edgeColor,'faceAlpha',faceAlpha,'lineStyle',lineStyle)
                end
                
            % Variable faceColor plots   
            else
                % mesh quality plots
                if strcmp(plotType, 'maxincludedangle')
                    color = obj.faceMaxIncludedAngles();
                elseif strcmp(plotType, 'minincludedangle')
                    color = obj.faceMinIncludedAngles();
                elseif strcmp(plotType, 'alpha')
                    color = obj.faceAverageMeshAngles;
                    
                % user defined face color by scalar variable for each face
                elseif strcmp(plotType, 'userdefined') && ...
                       size(faceColor,1) == obj.numFaces &&...
                       size(faceColor,2) == 1
                    color = faceColor;
                    
                % user defined face color by rgb vector for each face
                elseif strcmp(plotType, 'userdefined') && ...
                       size(faceColor,1) == obj.numFaces &&...
                       size(faceColor,2) == 3
                    color = faceColor;
                    customColor = 1;
                else
                    error('incorrect plotType, available options: maxIncludedAngle, minIncludedAngle, alpha, userDefined, wireFrame')
                end
                
                color01 = (color-min(color))/(max(color)-min(color));
                for i =1:obj.numFaces
                    if customColor==0
                        colori = color01(i)*color0 + (1-color01(i))*color1;
                    else
                        colori=color(i,:);
                    end
                    colori=colori*lighting(i);
                    indices = [3*i-2,3*i-1,3*i,3*i-2];
                    fill3(pts(indices,1),pts(indices,2),pts(indices,3),...
                        colori,'EdgeColor',edgeColor,'faceAlpha',faceAlpha,'lineStyle',lineStyle)
                    
                end
            end
            
            if isBounded == 1
                xlim([bounds(1),bounds(2)])
                ylim([bounds(3),bounds(4)])
                zlim([bounds(5),bounds(6)])
            end
            daspect([1,1,1])
            set(gcf,'Color',[1,1,1])  
            xlabel('x')
            ylabel('y')
            zlabel('z')
        end
        function plotFace(obj,faceIndex)
        % plot a single face (for diagnostics)
        %------------------------------------------------------------------
        % Inputs:
        %   faceIndex -- id of face to plot
        %------------------------------------------------------------------
            hold on
            % Now plot the stuff
            %--------------------------------------------------------------
            pts = obj.coordinates(obj.ends,:);
            for i =faceIndex
                indices = [3*i-2,3*i-1,3*i,3*i-2];
                fill3(pts(indices,1),pts(indices,2),pts(indices,3),...
                    [1,1,1],'EdgeColor','k')
            end
            daspect([1,1,1])
            set(gcf,'Color',[1,1,1])
            xlabel('x')
            ylabel('y')
            zlabel('z')
        end
    end
    
    methods(Access=private)
        function validity = isValidCollapse(obj,he1)
        % add check on spoke count for vertices P2 and P4
        %------------------------------------------------------------------
        %             o p2             | The collapse could still mess
        %            / ^               | things up, this is just a quick
        %           /   \              | and easy check to help prevent
        %      he3 /     \ he2         | that from happening in most cases
        %         /       \            | 
        %        v         \           | 
        %    p3 o----he1--->o p1       | 
        %        \<--he4---^           | 
        %         \       /            |
        %      he5 \     / he6         | 
        %           \   /              |
        %            v /               |
        %             o p4             | 
        %------------------------------------------------------------------
        % Inputs:
        %   he1 ------ collapse half edge pointing to collapse vertex
        %------------------------------------------------------------------
        % Outputs:
        %   validity - True if you "can't see it from my house"
        %------------------------------------------------------------------
       
            validity = false;
            
            he2 = obj.next(he1);
            he4 = obj.pair(he1);
            he5 = obj.next(he4);
            
            %p1 = obj.ends(he1);
            p2 = obj.ends(he2);
            %p3 = obj.ends(he4);
            p4 = obj.ends(he5);
            
            p2NumSpokes = obj.numSpokes(p2);
            p4NumSpokes = obj.numSpokes(p4);
            
            minSpokeCount = min(p2NumSpokes,p4NumSpokes);
            
            if minSpokeCount > 3
                validity = true;
            end
        end
        function [modifiedEdges,modifiedVertices] = edgeFlipModifiedRegion(obj,iVertex)
        % edge flip algo for ring around collapsed vertex
        %------------------------------------------------------------------
        %             o             
        %              ^\            
        %               \\           
        %           nexti\\ pairi     
        %                 \\          
        %                  \v          
        %       o----hei--->o -----------> o   
        %         <--------^/  <----------        
        %                 //          
        %                //      
        %               //             
        %              /v             
        %             o               
        %------------------------------------------------------------------
        %             o p2             |             o p2             
        %            / ^               |            / ^               
        %           /   \              |           / ^ \             
        %      he3 /     \ he2         |      he3 /  || \ he4         
        %         /       \            |         /   ||  \            
        %        v         \                    v    ||   \           
        %    p3 o----he1--->o p1      ==>   p3 o  he2||he5 o p1       
        %        \<--he4---^                    \    ||   ^  
        %         \       /            |         \   ||  /            
        %      he5 \     / he6         |      he1 \  || / he6          
        %           \   /              |           \  v/               
        %            v /               |            v /                           
        %             o p4             |             o p4             
        %------------------------------------------------------------------

        
            modifiedEdges = [];
            modifiedVertices = [];
            flipStack = [];
            
            % add in all incoming spokes and outer ring
            hei = obj.vertexHalfEdges(iVertex);
            pairi = hei;
            count = 0;
            while pairi ~= hei || count == 0
                count = count + 1;
                nexti = obj.next(pairi);
                nextnexti = obj.next(nexti);
                flipStack = [flipStack,pairi,nextnexti];
                pairi = obj.pair(nexti);
            end
            
            % test them all for edge flip and track modified edges
            stackLength = length(flipStack);
            while stackLength>=1
                he1 = flipStack(stackLength);
                flip = obj.flipCriterion(he1);
                if flip
                    he2 = obj.next(he1);
                    he3 = obj.next(he2);
                    he4 = obj.pair(he1);
                    he5 = obj.next(he4);
                    he6 = obj.next(he5);
                    
                    modifiedEdges=[modifiedEdges;he1;he2;he4;he5];
                    
                    p1 = obj.ends(he1);
                    p2 = obj.ends(he2);
                    p3 = obj.ends(he4);
                    p4 = obj.ends(he5);
                    
                    modifiedVertices=[modifiedVertices;p1;p2;p3;p4];
                    
                    obj.flipEdge(he1);
                    
                    flipStack=[flipStack(1:stackLength),he4,he3,he6];
                    stackLength=stackLength+3;
                else
                    stackLength=stackLength-1;
                end
            end
            modifiedEdges=unique(modifiedEdges);
            modifiedVertices=unique(modifiedVertices);
        end  
        function replaceHalfEdgeVertex(obj, iVertexOld , iVertexNew)
        % replaces vertex ids in half edge data structure with a new one
        %------------------------------------------------------------------
        % walks the vertex loop for the oldIndex vertex replacing it with
        % the new vertex for halfEdges that pointed to oldIndex vertex.
        %------------------------------------------------------------------
        % Inputs:
        %   iVertexOld -- id of vertex being replaced
        %   iVertexNew -- id of vertex being sub-ed in
        %------------------------------------------------------------------
            
            % original code
            if length(iVertexOld)==1
                startingEdge = obj.vertexHalfEdges(iVertexOld);
                obj.ends(startingEdge) = iVertexNew;
            
                nexti = obj.next(startingEdge);
                pairi = obj.pair(nexti);
            
                while  pairi ~= startingEdge
                    obj.ends(pairi) = iVertexNew;
                    nexti = obj.next(pairi);
                    pairi = obj.pair(nexti);
                end
              
            % for big data we want to vectorize   
            else
                startingEdge = obj.vertexHalfEdges(iVertexOld);
                obj.ends(startingEdge) = iVertexNew;
                
                nexti = obj.next(startingEdge);    
                pairi = obj.pair(nexti);
   
                i=1;
                while  any(pairi ~= startingEdge)
                    keepLoops = pairi ~= startingEdge;  
                    %disp(['Spoke: ',num2str(i)])
                    %disp(['  numEdgeLoops: ',num2str(length(keepLoops))])
                    %max(iVertexNew)
                    pairi = pairi(keepLoops); 
                    startingEdge = startingEdge(keepLoops); 
                    iVertexNew = iVertexNew(keepLoops);
                    %ids = ids(keepLoops);
                    
                    obj.ends(pairi) = iVertexNew;
                                      
                    nexti = obj.next(pairi);
                    pairi = obj.pair(nexti);
                    i=i+1;
                end   
            
            end
                
            
        end
        function fixCollapsedFacePairs(obj, he1)
        % fixing the pairing of half edges associated w/ a collapsed face
        %------------------------------------------------------------------
        % The half edges of adjacent faces to the collapsed faces will now
        % be paired with each other to restore the appropriate connectivity
        %------------------------------------------------------------------
        %             o p2             | 
        %            / ^               | 
        %           /   \              |
        %      he3 /     \ he2         |   
        %         /       \            | 
        %        v         \           | 
        %    p3 o----he1--->o p1       | 
        %        \<--he4---^           |
        %         \       /            |
        %      he5 \     / he6         |
        %           \   /              | 
        %            v /               |
        %             o p4             |  
        %------------------------------------------------------------------
        % Inputs:
        %   i-- index of half edge being collapsed
        %------------------------------------------------------------------
            
            
            he2 = obj.next(he1);
            he3 = obj.next(he2);
            
            % half edges of adjacent faces
            he2pair = obj.pair(he2);
            he3pair = obj.pair(he3);
            
            % set the adjacent pairs 
            obj.pair(he2pair) = he3pair;
            obj.pair(he3pair) = he2pair;
            
            % repeat for pair
            he4 = obj.pair(he1);
            
            he5 = obj.next(he4);
            he6 = obj.next(he5);
            
            % half edges of adjacent faces
            he5pair = obj.pair(he5);
            he6pair = obj.pair(he6);
            
            % set the adjacent pairs 
            obj.pair(he5pair) = he6pair;
            obj.pair(he6pair) = he5pair;
                
        end
        function correctVertexHalfEdgeRelations(obj,removedHalfEdges) 
        % vertices pointing to specified he's switched to next he
        %------------------------------------------------------------------
        % loops through a list of half edge indices. For each specified
        % index tests if the half edge's vertex points to the specified
        % half edge. If so, the vertex is switched to point to another
        % one of its spoke half edges. This method is used in the edge
        % collapse algorithm to fix mesh connectivity after half edge
        % removal.
        %------------------------------------------------------------------
        % Inputs:
        %   removedHalfEdges -- list of he indices
        %------------------------------------------------------------------
            
            for j = 1:length(removedHalfEdges)
                iVertex = obj.ends(removedHalfEdges(j));

                % take one
                if obj.vertexHalfEdges(iVertex) == removedHalfEdges(j)
                    nexti = obj.next(removedHalfEdges(j));
                    pairi = obj.pair(nexti);

                    % take two
                    if any(removedHalfEdges == pairi)
                        nexti = obj.next(pairi);
                        pairi = obj.pair(nexti);
                        
                        % it still no good we have a flap face
                        if any(removedHalfEdges == pairi)
                            error('degenerate flap face: 2 faces w/ same vertices')
                        end
                    end

                    obj.vertexHalfEdges(iVertex) = pairi;
                end
            end
        end
        function selectionMetric = initializeSelectionMetric(obj,isFeatureBased)
        % get our selection metric
        %------------------------------------------------------------------
            modifiedEdges = 1:obj.numHalfEdges;
            modifiedVertices = 1:obj.numVertices;
            selectionMetric = obj.calculateSelectionMetric(modifiedEdges,modifiedVertices,isFeatureBased);
            
        end
        function modifiedMetric = calculateSelectionMetric(obj,modifiedEdges,modifiedVertices,isFeatureBased)
        % update selection metric post collapse
        %------------------------------------------------------------------

            % update lengths of modified edges from collapse
            %------------------------------------------------------
            modifiedFixedEdges = obj.isFixedHalfEdge(modifiedEdges);
            modifiedEdges = modifiedEdges(modifiedFixedEdges==0);
            modifiedMetric = obj.halfEdgeLengths(modifiedEdges);

            if isFeatureBased
                angleVertices(modifiedVertices) = obj.vertexAverageAngles(modifiedVertices);
                modifiedAngles = (angleVertices(obj.ends(modifiedEdges)) +...
                    angleVertices(obj.ends(obj.pair(modifiedEdges))))/2;
                modifiedMetric=modifiedMetric.*modifiedAngles;
            end

        end
        function newLowValenceVertices = updateLowValenceVertices(obj,vertexIndices)
        % adds vertices to lowValenceCandidates if 3 spokes or fewer
        %------------------------------------------------------------------
        % Inputs:
        %   vertexIndicies -- vertex indices to test for low spoke count
        %------------------------------------------------------------------
        % Outputs:
        %   newLowValencyVertices -- indices of vertices with few spokes
        %------------------------------------------------------------------
            newLowValenceVertices = [];
            for j=1:length(vertexIndices)
                pj = vertexIndices(j);
                numSpokes = obj.numSpokes(pj);
                if numSpokes == 3
                    hej = obj.spokeHalfEdges(pj);
                    newLowValenceVertices=[newLowValenceVertices;hej];
                elseif numSpokes == 2
                    disp(pj)
                    disp(obj.coordinates(pj,:))
                    error("collapse produced a degenerate flap face")
                end
            end
        end
        function cleanDeletedEntities(obj,vertexFlags,halfEdgeFlags)
        % cleans after edge collapse algorithm
        %------------------------------------------------------------------
        % Inputs:
        %   vertexFlags ---- vertices removed from mesh
        %   halfEdgeFlags -- half edges removed from mesh
        %------------------------------------------------------------------

            % map old indices to new 
            vertexNewIndices   = zeros(obj.numVertices,1);
            halfEdgeNewIndices = zeros(obj.numHalfEdges,1);

            iter = 0;
            for i = 1:obj.numHalfEdges
                if halfEdgeFlags(i) == 0
                    iter = iter+1;
                    halfEdgeNewIndices(i) = iter;
                end
            end
            
            iter = 0;
            numDeletedVertices=0;
            for i = 1:obj.numVertices
                if vertexFlags(i) == 0
                    iter = iter+1;
                    vertexNewIndices(i) = iter;
                else
                    numDeletedVertices=numDeletedVertices+1;
                end
            end
            
            % read our fixed flags
            obj.isFixedHalfEdge = obj.isFixedHalfEdge(halfEdgeFlags==0);
            obj.isFixedVertex = obj.isFixedVertex(vertexFlags==0);

            % set relations to new indices
            obj.ends = vertexNewIndices(obj.ends);
            obj.pair = halfEdgeNewIndices(obj.pair);
            obj.next = halfEdgeNewIndices(obj.next);
            obj.vertexHalfEdges = halfEdgeNewIndices(obj.vertexHalfEdges);
            
            % remove the inactive elements
            obj.ends = obj.ends(halfEdgeFlags==0);
            obj.next = obj.next(halfEdgeFlags==0);
            obj.pair = obj.pair(halfEdgeFlags==0);

            obj.coordinates = obj.coordinates(vertexFlags==0,:);
            obj.vertexHalfEdges = obj.vertexHalfEdges(vertexFlags==0,:);
            
            
            obj.faces = [obj.ends(1:3:end),...
                         obj.ends(2:3:end),...
                         obj.ends(3:3:end)];
            
            % correct our counter variables
            obj.numHalfEdges = obj.numHalfEdges - 6*numDeletedVertices;
            obj.numEdges = obj.numEdges - 3*numDeletedVertices;
            obj.numVertices = obj.numVertices - numDeletedVertices;
            obj.numFaces = obj.numFaces - 2*numDeletedVertices;
            obj.numNodes = obj.numVertices;
        end
    end

    methods(Access=public)

        function [halfEdges,halfEdgeAngles] = outerRingAngles(obj,iVertex)
        % iterate over the spoke and calc outter angle attached to each
        %------------------------------------------------------------------
        % Inputs:
        %   iVertex -------- index of central vertex
        %------------------------------------------------------------------
        % Outputs:
        %   halfEdges ------ half edges spokes for vertex i
        %   halfEdgeAngles - angle of outer ring at junction w/ half edge
        %------------------------------------------------------------------
            halfEdges = [];
            halfEdgeAngles = [];
            
    
            startingEdge = obj.vertexHalfEdges(iVertex);
            startingPair = obj.pair(startingEdge);
            startingNext = obj.next(startingPair);
            
            p0 = obj.coordinates(obj.ends(startingNext),:);
            p1 = obj.coordinates(obj.ends(startingPair),:);
            
            pairi = startingEdge;
            nexti = obj.next(startingEdge);
            
            i=1;
            
            while  pairi ~= startingEdge || i==1
                
                p2 = obj.coordinates(obj.ends(nexti),:);
                
                v1 = p0-p1;
                v2 = p2-p1;
                
                halfEdges(i) =  pairi;
                halfEdgeAngles(i) = acosd(v2*v1'/(sqrt(v1*v1')*sqrt(v2*v2')));
                i = i+1;
                
                p0 = p1;
                p1 = p2;
                
                pairi = obj.pair(nexti);
                nexti = obj.next(pairi);
                
            end
        end
        function count = numSpokes(obj,iVertex)
        % counts number of unique edges connected to vertex 
        %------------------------------------------------------------------
        
            count = 0;
            
            startingEdge = obj.vertexHalfEdges(iVertex);

            pairi = startingEdge;
            nexti = obj.next(startingEdge);
            
            while  pairi ~= startingEdge || count==0
                
                count = count+1;
                
                pairi = obj.pair(nexti);
                nexti = obj.next(pairi);
                if count > 20
                    warning(['infinite vertex loop: ',num2str(iVertex)])
                    count = 1;
                    break
                end
                
            end
            
        end
        function halfEdges = spokeHalfEdges(obj,iVertex)
            
            halfEdges = zeros(20,1);
            
            % initialze the iteration
            startingEdge = obj.vertexHalfEdges(iVertex);

            pairi = startingEdge;
            nexti = obj.next(startingEdge);
            
            i=0;
            
            while  pairi ~= startingEdge || i==0
                
                i = i+1;
                halfEdges(2*i) =  pairi;
                halfEdges(2*i-1) =  nexti;
                
                pairi = obj.pair(nexti);
                nexti = obj.next(pairi);
                
            end
            
            halfEdges = halfEdges(1:2*i);
            
        end
        function vertices = neighborVertices(obj,iVertex)
        % returns all vertices connected to iVertex
        %------------------------------------------------------------------
        % Inputs:
        %   iVertex --- vertex index
        %------------------------------------------------------------------
        % Outputs:
        %   vertices -- neighbor vertex indices
        %------------------------------------------------------------------
        
            vertices = zeros(15,1);
            
            % initialze the iteration
            startingEdge = obj.vertexHalfEdges(iVertex);

            pairi = startingEdge;
            nexti = obj.next(startingEdge);
            
            i=0;
            
            while  pairi ~= startingEdge || i==0
                
                i = i+1;
                vertices(i) =  obj.ends(nexti);
                
                pairi = obj.pair(nexti);
                nexti = obj.next(pairi);
                
            end
            
            vertices = vertices(1:i);
            
        end
        function vertices = collapseModifiedVertices(obj,iVertex)
        % return vertices in ring around iVertex including iVertex
        %------------------------------------------------------------------
        
            %preallocate
            vertices = zeros(20,1);
           
            
            % seed halfedges
            startingEdge = obj.vertexHalfEdges(iVertex);
            pairi = startingEdge;
            nexti = obj.next(startingEdge);
            
            % seed vertices
            vertices(1) = iVertex;
            vertices(2) = obj.ends(nexti);
            
            % start index
            i=2;
            
            while  pairi ~= startingEdge || i==2
                
                i = i+1;
                
                % add our ring indices
                pairi = obj.pair(nexti);
                nexti = obj.next(pairi);
                vertices(i) = obj.ends(nexti);
                
            end
            
            vertices=vertices(1:i);
            
        end
        function flip = flipCriterion(obj,he1)
        % returns true if min included angle is increased by edge flipping
        %------------------------------------------------------------------
        %             o p2             | Returns true if the min included
        %            / ^               | angle is increased by a
        %           /   \              | multiplicative threshold.
        %      he3 /     \ he2         |   
        %         /       \            | Max included angle of the  
        %        v         \           | quadralateral is used to ensure
        %    p3 o----he1--->o p1       | edge flipping won't result in
        %        \<--he4---^           | face inversion.
        %         \       /            |
        %      he5 \     / he6         | This gets called a lot so the
        %           \   /              | angle calculation are explicitly
        %            v /               | written out for speed
        %             o p4             |
        %------------------------------------------------------------------
        %   Original Configuration     |     Flipped Configuration
        %------------------------------------------------------------------          
        %             o p2             |              o p2 
        %            / ^               |             /|^               
        %           /   \              |            / | \          
        %      v2  /     \ v1          |       v2  /2b|2a\ v1         
        %         /       \            |          /   |   \        
        %        v 3a   1b \           |         v    |    \    
        %    p3 o<---v5-----o p1       |     p3 o    v6     o p1   
        %        \ 3b   1a ^           |         \    |    ^            
        %         \       /            |          \   |   /   
        %      v3  \     / v4          |       v3  \4a|4b/ v4       
        %           \   /              |            \ | /             
        %            v /               |             vv/             
        %             o p4             |              o p4 
        %------------------------------------------------------------------
        % Inputs:
        %   he1 -- collapsed half edge pointing towards collapsed nodes
        %------------------------------------------------------------------
        % Outputs:
        %   flip - boolean true if edge should be fliped
        %------------------------------------------------------------------
            
            flip = false;
            
            % relevant halfEdges
            he2 = obj.next(he1);
            he4 = obj.pair(he1);
            he5 = obj.next(he4);

            % points defining quadralateral 
            coord1 = obj.coordinates(obj.ends(he1),:);
            coord2 = obj.coordinates(obj.ends(he2),:);
            coord3 = obj.coordinates(obj.ends(he4),:);
            coord4 = obj.coordinates(obj.ends(he5),:);
            
            % vectors defining quadralateral
            v1 = coord2-coord1;
            v2 = coord3-coord2;
            v3 = coord4-coord3;
            v4 = coord1-coord4;
           
            % face normals (inlined critical code)
            n1 = [v1(2).*v2(3) - v1(3).*v2(2),...
                -v1(1).*v2(3) + v1(3).*v2(1),...
                v1(1).*v2(2) - v1(2).*v2(1)];
            n2 = [v3(2).*v4(3) - v3(3).*v4(2),...
                -v3(1).*v4(3) + v3(3).*v4(1),...
                v3(1).*v4(2) - v3(2).*v4(1)];
            n3 = [v4(2).*v1(3) - v4(3).*v1(2),...
                -v4(1).*v1(3) + v4(3).*v1(1),...
                v4(1).*v1(2) - v4(2).*v1(1)];
            n4 = [v2(2).*v3(3) - v2(3).*v3(2),...
                -v2(1).*v3(3) + v2(3).*v3(1),...
                v2(1).*v3(2) - v2(2).*v3(1)];
            n1Mag = sqrt(n1(1).*n1(1) + n1(2).*n1(2) + n1(3).*n1(3));
            n2Mag = sqrt(n2(1).*n2(1) + n2(2).*n2(2) + n2(3).*n2(3));
            n3Mag = sqrt(n3(1).*n3(1) + n3(2).*n3(2) + n3(3).*n3(3));
            n4Mag = sqrt(n4(1).*n4(1) + n4(2).*n4(2) + n4(3).*n4(3));
            n1 = n1./n1Mag;
            n2 = n2./n2Mag;
            n3 = n3./n3Mag;
            n4 = n4./n4Mag;
            
            % dihedral angles
            dihedral1 = sum(n1.*n2);
            dihedral2 = sum(n3.*n4);

            % split vectors
            v5 = coord3-coord1;
            v6 = coord4-coord2;

            % lengths
            d1 = sqrt(v1(1).*v1(1)+v1(2).*v1(2)+v1(3).*v1(3));
            d2 = sqrt(v2(1).*v2(1)+v2(2).*v2(2)+v2(3).*v2(3));
            d3 = sqrt(v3(1).*v3(1)+v3(2).*v3(2)+v3(3).*v3(3));
            d4 = sqrt(v4(1).*v4(1)+v4(2).*v4(2)+v4(3).*v4(3));
            d5 = sqrt(v5(1).*v5(1)+v5(2).*v5(2)+v5(3).*v5(3));
            d6 = sqrt(v6(1).*v6(1)+v6(2).*v6(2)+v6(3).*v6(3));

            % angles of triangles along quadralateral split
            angle1A = acosd(-(v5(1).*v4(1)+v5(2).*v4(2)+v5(3).*v4(3))/(d4*d5));
            angle1B = acosd((v5(1).*v1(1)+v5(2).*v1(2)+v5(3).*v1(3))/(d1*d5));
            angle2A = acosd(-(v6(1).*v1(1)+v6(2).*v1(2)+v6(3).*v1(3))/(d1*d6));
            angle2B = acosd((v6(1).*v2(1)+v6(2).*v2(2)+v6(3).*v2(3))/(d2*d6));
            angle3A = acosd((v2(1).*v5(1)+v2(2).*v5(2)+v2(3).*v5(3))/(d2*d5));
            angle3B = acosd(-(v3(1).*v5(1)+v3(2).*v5(2)+v3(3).*v5(3))/(d3*d5));
            angle4A = 180-angle2B-angle3A-angle3B;
            angle4B = 180-angle2A-angle1A-angle1B;

            currentMinAngle = min([angle1A,angle1B,angle3A,angle3B]);
            flippedMinAngle = min([angle2A,angle2B,angle4A,angle4B]);

            maxQuadAngles = max(angle1A+angle1B, angle3A+angle3B);
            
            % conditions for flip
            isFlippable1 = flippedMinAngle > 1.1*currentMinAngle && maxQuadAngles<175;
            isFlippable2 = flippedMinAngle > 10 && dihedral2 > dihedral1+0.2;
            if isFlippable2 || isFlippable1
                flip=true;
            end
           
        end
        function flipEdge(obj,he1)
        % flips connectivity of face pair
        %------------------------------------------------------------------
        % split half-edges (he1,he4) switched to point to p4 and p2 
        % respectively. This alters the pairing of he2,he5, he1 and he4 so 
        % that all has o be adjusted. If the vertex points to a he that
        % was altered the vertex's he needs to be switched to another that 
        % still points to it.
        %------------------------------------------------------------------
        %             o p2             |             o p2             
        %            / ^               |            / ^               
        %           /   \              |           / ^ \             
        %      he3 /     \ he2         |      he3 /  || \ he4         
        %         /       \            |         /   ||  \            
        %        v         \                    v    ||   \           
        %    p3 o----he1--->o p1      ==>   p3 o  he2||he5 o p1       
        %        \<--he4---^                    \    ||   ^  
        %         \       /            |         \   ||  /            
        %      he5 \     / he6         |      he1 \  || / he6          
        %           \   /              |           \  v/               
        %            v /               |            v /                           
        %             o p4             |             o p4             
        %------------------------------------------------------------------
        % Inputs:
        %   he1 -- index of bissecting half edge
        %------------------------------------------------------------------
        
            %face 1
            he2 = obj.next(he1);
            he3 = obj.next(he2);
            
            % face 2
            he4 = obj.pair(he1);
            he5 = obj.next(he4);
            he6 = obj.next(he5);

            % check if we're going to mess up vertex pairing
            if obj.vertexHalfEdges(obj.ends(he1))==he1
                obj.vertexHalfEdges(obj.ends(he1)) = he6;
            end
            if obj.vertexHalfEdges(obj.ends(he4))==he4
                obj.vertexHalfEdges(obj.ends(he4)) = he3;
            end
            
            % switch which vertices split edges point to
            obj.ends(he1) = obj.ends(he5);
            obj.ends(he4) = obj.ends(he2);
            
            % fix the pairing 
            he5pair = obj.pair(he5);
            he2pair = obj.pair(he2);
           
            obj.pair(he5) = he2;
            obj.pair(he2) = he5;
            
            obj.pair(he5pair) = he1;
            obj.pair(he2pair) = he4;
            
            obj.pair(he1) = he5pair;
            obj.pair(he4) = he2pair;
            
        end
        
        function [obj,newRefineLevel] = refineAndTrackLevel(obj,refineFaces,refineLevel)
        % converts each tri to 4. Allows subset of faces to be refined.
        %------------------------------------------------------------------
        %                 coord2                | generates 4 triangles 
        %                o                      | with he ordering as 
        %  he2        / / ^ ^                   | indicated. for select
        %  (coarse)  / /   \ \                  | faces w/ proper connect 
        %           / /5   4\ \                 | splits.
        %          / /       \ \
        %         / /    6    \ \ he1 (coarse)
        %        / o<---------o  \
        %       / / \   11   ^ ^  \
        %      / /  ^\12  10//  \  \           type-3 refinement
        %     / /7  9\\    //2  1\  \
        %    / /      \\  //      \  \
        %   v v   8    \v/v   3    \  \           
        %    o--------->o---------->o coord1
        %     --------------------->
        %              he3 (coarse)
        %------------------------------------------------------------------       
        % Inputs:
        %   refineFaces -- indices of faces to be refined
        %   refineLevel -- split level of incoming faces
        %------------------------------------------------------------------
        % Outputs:
        %   newRefineLevels -- new levels for refined mesh.
        %------------------------------------------------------------------
        
            if nargin ==2 
                refineLevel = zeros(obj.numFaces,1);
            end
            
            % make sure we're rectilinear and faces connectivity is correct
            obj.flatten();
            
            if nargin==1
            
                % get unique edge defs and halfEdge --> uniqueEdge transform
                [edges,~,halfEdgeEdges] = obj.edges();
            
                % add edge midpoints to coordinates
                tempVertices = [obj.coordinates;...
                               (obj.coordinates(edges(:,1),:)+...
                                obj.coordinates(edges(:,2),:))/2];
            
                % split into 4 faces
                tempFaces = zeros(4*obj.numFaces,3);
                tempFaces(1:4:end,1) = obj.numVertices + halfEdgeEdges(1:3:end,1);
                tempFaces(1:4:end,2) = obj.numVertices + halfEdgeEdges(3:3:end,1);
                tempFaces(1:4:end,3) = obj.faces(:,3);
                tempFaces(2:4:end,1) = obj.faces(:,1);
                tempFaces(2:4:end,2) = obj.numVertices + halfEdgeEdges(2:3:end,1);
                tempFaces(2:4:end,3) = obj.numVertices + halfEdgeEdges(1:3:end,1);
                tempFaces(3:4:end,1) = obj.faces(:,2);
                tempFaces(3:4:end,2) = obj.numVertices + halfEdgeEdges(3:3:end,1);
                tempFaces(3:4:end,3) = obj.numVertices + halfEdgeEdges(2:3:end,1);
                tempFaces(4:4:end,1) = obj.numVertices + halfEdgeEdges(1:3:end,1);
                tempFaces(4:4:end,2) = obj.numVertices + halfEdgeEdges(2:3:end,1);
                tempFaces(4:4:end,3) = obj.numVertices + halfEdgeEdges(3:3:end,1);
                
                newRefineLevel(1:4:end) = refineLevel+1;
                newRefineLevel(2:4:end) = refineLevel+1;
                newRefineLevel(3:4:end) = refineLevel+1;
                newRefineLevel(4:4:end) = refineLevel+1;
                
            elseif nargin == 2 || nargin == 3
                
                if max(refineFaces) > obj.numFaces ||...
                   any(refineFaces) <= 0 ||...
                   min(size(refineFaces)) ~= 1
                    error('something went wrong with refineFaces input')
                end
                   
                % flag half edges of refine faces and their pairs
                halfEdgeIsSplit = zeros(obj.numHalfEdges,1);
                halfEdgeIsSplit(3*refineFaces-2) = 1;
                halfEdgeIsSplit(3*refineFaces-1) = 1;
                halfEdgeIsSplit(3*refineFaces)   = 1;
                halfEdgeIsSplit(obj.pair(3*refineFaces-2)) = 1;
                halfEdgeIsSplit(obj.pair(3*refineFaces-1)) = 1;
                halfEdgeIsSplit(obj.pair(3*refineFaces))   = 1;
                
                % determine split type
                splitType = halfEdgeIsSplit(1:3:end) + ...
                            halfEdgeIsSplit(2:3:end) + ...
                            halfEdgeIsSplit(3:3:end);
                        
                numNewFaces = sum(splitType)+obj.numFaces;

                % unique edge that split
                [edges,edgeHalfEdges,halfEdgeEdges] = obj.edges();
                edgeIsSplit = zeros(obj.numEdges,1);
                for i = 1:obj.numEdges
                    if halfEdgeIsSplit(edgeHalfEdges(i,1)) || halfEdgeIsSplit(edgeHalfEdges(i,2))
                        edgeIsSplit(i)=1;
                    end
                end
                splitEdges = 1:obj.numEdges;
                splitEdges = splitEdges(edgeIsSplit==1);
                edges2SplitEdges = zeros(obj.numEdges,1);
                edges2SplitEdges(splitEdges) = 1:length(splitEdges);

                % add all midpoints to tempVertices we'll cull later
                tempVertices = [obj.coordinates;...
                               (obj.coordinates(edges(splitEdges,1),:)+...
                                obj.coordinates(edges(splitEdges,2),:))/2];
                            
                tempFaces = zeros(numNewFaces,3);
                newRefineLevel = zeros(numNewFaces,1);
                newFaceIndex = 1;
                
                % loop through old face and make new faces based on the
                % number of edges split.
                for i = 1:obj.numFaces
                    
                    % no edges split
                    if splitType(i) == 0
                        tempFaces(newFaceIndex,:)=obj.faces(i,:);
                        newFaceIndex = newFaceIndex+1;
                    
                    % one edge split
                    elseif splitType(i) == 1
                        
                        if halfEdgeIsSplit(3*i-2) == 1
                            he1 = 3*i-2;
                        elseif halfEdgeIsSplit(3*i-1) == 1
                            he1 = 3*i-1;
                        elseif halfEdgeIsSplit(3*i) == 1
                            he1 = 3*i;
                        else
                            error('something went wrong here')
                        end
                       
                        he2 = obj.next(he1);
                        he3 = obj.next(he2);
                        
                        p1 = obj.ends(he1);
                        p2 = obj.ends(he2);
                        p3 = obj.ends(he3);
                        
                        uniqueEdgeId = halfEdgeEdges(he1);
                        uniqueSplitEdgeId = edges2SplitEdges(uniqueEdgeId);
                        p4 = obj.numVertices+uniqueSplitEdgeId;
                        
                        tempFaces(newFaceIndex,:)   = [p1,p2,p4];
                        tempFaces(newFaceIndex+1,:) = [p4,p2,p3];
                        
                                                
                        newRefineLevel(newFaceIndex)   = refineLevel(i)+1;
                        newRefineLevel(newFaceIndex+1) = refineLevel(i)+1;
                        
                        newFaceIndex = newFaceIndex+2;
                    
                    % two edges split
                    elseif splitType(i) == 2
                        
                        if halfEdgeIsSplit(3*i-2) == 0
                            he1 = 3*i-2;
                        elseif halfEdgeIsSplit(3*i-1) == 0
                            he1 = 3*i-1;
                        elseif halfEdgeIsSplit(3*i) == 0
                            he1 = 3*i;
                        else
                            error('something went wrong here')
                        end
                       
                        he2 = obj.next(he1);
                        he3 = obj.next(he2);
                        
                        p1 = obj.ends(he1);
                        p2 = obj.ends(he2);
                        p3 = obj.ends(he3);
                        
                        uniqueEdgeId = halfEdgeEdges(he2);
                        uniqueSplitEdgeId = edges2SplitEdges(uniqueEdgeId);
                        p4 = obj.numVertices+uniqueSplitEdgeId;
                        
                        uniqueEdgeId = halfEdgeEdges(he3);
                        uniqueSplitEdgeId = edges2SplitEdges(uniqueEdgeId);
                        p5 = obj.numVertices+uniqueSplitEdgeId;
                        
                        tempFaces(newFaceIndex,:)   = [p1,p4,p3];
                        tempFaces(newFaceIndex+1,:) = [p3,p4,p5];
                        tempFaces(newFaceIndex+2,:) = [p5,p4,p2];
                        
                        newRefineLevel(newFaceIndex)   = refineLevel(i)+1;
                        newRefineLevel(newFaceIndex+1) = refineLevel(i)+1;
                        newRefineLevel(newFaceIndex+2) = refineLevel(i)+1;
                        
                        newFaceIndex = newFaceIndex+3;
                    
                    % three edges split
                    elseif splitType(i) == 3

                        he1 = 3*i-2;
                        he2 = obj.next(he1);
                        he3 = obj.next(he2);
                        
                        p1 = obj.ends(he1);
                        p2 = obj.ends(he2);
                        p3 = obj.ends(he3);
                        
                        uniqueEdgeId = halfEdgeEdges(he1);
                        uniqueSplitEdgeId = edges2SplitEdges(uniqueEdgeId);
                        p4 = obj.numVertices+uniqueSplitEdgeId;
                        
                        uniqueEdgeId = halfEdgeEdges(he2);
                        uniqueSplitEdgeId = edges2SplitEdges(uniqueEdgeId);
                        p5 = obj.numVertices+uniqueSplitEdgeId;
                        
                        uniqueEdgeId = halfEdgeEdges(he3);
                        uniqueSplitEdgeId = edges2SplitEdges(uniqueEdgeId);
                        p6 = obj.numVertices+uniqueSplitEdgeId;
                        
                        tempFaces(newFaceIndex,:)   = [p4,p5,p6];
                        tempFaces(newFaceIndex+1,:) = [p4,p1,p5];
                        tempFaces(newFaceIndex+2,:) = [p2,p6,p5];
                        tempFaces(newFaceIndex+3,:) = [p4,p6,p3];
                        
                        newRefineLevel(newFaceIndex)   = refineLevel(i)+1;
                        newRefineLevel(newFaceIndex+1) = refineLevel(i)+1;
                        newRefineLevel(newFaceIndex+2) = refineLevel(i)+1;
                        newRefineLevel(newFaceIndex+3) = refineLevel(i)+1;
                        
                        newFaceIndex = newFaceIndex+4;
                    end  
                end   
            else
                error('incorrect number of inputs')
            end
            
            obj.initializeFromFaceData(tempVertices,tempFaces);

        end
    end
    methods(Static,Access=public)
        function [isIntersecting,intersectionCoord]...
                                = isVectorIntersectingTriangle(...
                                     coord0,n0,...
                                  coord1,coord2,coord3,nf)
        % determines if vector intersects triangle in 3d
        %------------------------------------------------------------------
        %                                
        %                     o coord3  
        %                    / \                          (0,1)
        %              vec0 /   \                        o 
        %    -----o -------->o   \  triNormal            | \
        % coord0 ^        /       \----->                |   \
        %       /        /         \            ==>    t | o   \
        %      /        o ----------o coord1             |       \
        %     /         coord2                           |         \ (1,0)
        %    /                                     (0,0) o-----------o
        %                                                     s
        %------------------------------------------------------------------
        % 1) find intersection point of vector in face plane. 
        % 2) covert to parametric coords of triangle (s,t)
        % 3) test in/out
        %------------------------------------------------------------------
        % Inputs: 
        %   coord0 -- anchor for intersecting vector
        %   n0 ------ vector direction
        %   coord1 -- coordinates of vertex 1
        %   coord2 -- coordinates of vertex 2
        %   coord3 -- coordinates of vertex 3
        %   nf ------ outward faceing normal CCW
        %------------------------------------------------------------------
        % Outputs:
        %   isIntersecting ---- bool true if intersection point is 
        %                       contained within the triangle in question
        %   intersectionCoord - coord of project point in face plane
        %------------------------------------------------------------------    

            isIntersecting = false;
            
            % t - distance from coord0 to face plane along vec0
            distance = - (nf(1)*(coord0(1)-coord1(1))+...
                          nf(2)*(coord0(2)-coord1(2))+...
                          nf(3)*(coord0(3)-coord1(3)))/...
                (nf(1)*n0(1)+nf(2)*n0(2)+nf(3)*n0(3));
            
            % point of intersection
            intersectionCoord = coord0 + distance*n0;
            
            v1 = coord2-coord1;            % edge vector 1
            v2 = coord3-coord1;            % edge vector 2
            v0 = intersectionCoord-coord1; % vector to intersection point
            
            % calc parametric coordinates of intersection point
            dotv2v2 = v2(1)*v2(1) + v2(2)*v2(2) + v2(3)*v2(3);
            dotv1v2 = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3);
            dotv1v1 = v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3);
            dotv0v2  =  v0(1)*v2(1) + v0(2)*v2(2) +  v0(3)*v2(3);
            dotv0v1  =  v0(1)*v1(1) + v0(2)*v1(2) +  v0(3)*v1(3);
            denom = dotv1v2 ^ 2 - dotv1v1 * dotv2v2;
            
            s = (dotv1v2 * dotv0v2 - dotv2v2 * dotv0v1)/denom;
            t = (dotv1v2 * dotv0v1 - dotv1v1 * dotv0v2)/denom;
            
            threshold = 1e-6*sqrt(dotv2v2);
            
            if s > -threshold && t > -threshold && s+t<1+threshold
                isIntersecting = true;
            end
        end
        function n = orthogonalNormal(a,b)
        % returns the othogonal normal of Nx3 Nx3 arrays of vectors
        %-----------------------------------------------------------------
        % Inputs:
        %   a - Nx3 array containing coords of N vectors
        %   b - same 
        %------------------------------------------------------------------
        % Output:
        %   n - orthogonal unit vector
        %------------------------------------------------------------------
            n = [ a(:,2).*b(:,3) - a(:,3).*b(:,2),...
                 -a(:,1).*b(:,3) + a(:,3).*b(:,1),...
                  a(:,1).*b(:,2) - a(:,2).*b(:,1)];
            nMag = sqrt(n(:,1).*n(:,1) + n(:,2).*n(:,2) + n(:,3).*n(:,3));
            n = n./nMag;
        end
        function [pts,f] = readOBJ(fileName)
        % read in surface mesh from .obj file
        %------------------------------------------------------------------
        % Inputs:
        %   filename - string specifing file name
        %------------------------------------------------------------------
         
            assert(contains(fileName,'.obj'),'unsupported file type')
            
            vertexchar = [''];
            facechar = [''];
            
            fid = fopen(fileName,'r');
            tline = fgetl(fid);
            while ischar(tline)
                str = strtrim(tline);
                if isempty(tline)
                    %disp('empty')
                elseif str(1) == 'v'
                    vertexchar = [vertexchar,str(2:end)];
                elseif str(1) =='f'
                    facechar = [facechar,str(2:end)];
                end
                tline = fgetl(fid);
            end
            fclose(fid);
            pts = cell2mat(textscan(vertexchar,'%f%f%f'));
            f = round(cell2mat(textscan(facechar,'%u%u%u')));
        end
        function [] = displayQuadratureRules()
        % lists the available quadrature rules and reccomended settings
        %------------------------------------------------------------------
        % R. Lauffer. Interpolation mehrfacher intergrale. Archiv der 
        % Mathematik, 6:159-164, 1955.
        %
        % A.H. Stroud. Approximate Calculation of Multiple Integrals.
        % Prentice Hall, Inc, Englewood Cliffs, NJ, 1971.
        %
        % J.N. Lyness and D. Jespersen. "Moderate degree symmetric 
        % quadrature rules for the triangle," Journal of the Institute of 
        % Mathematics and its Applications, Vol.15, pp. 19-32, 1975.
        %
        % R. Cools and P. Rabinowitz. Monomial cubature rules since 
        % “stroud”: a compilation. Journal of Computational and Applied 
        % Mathematics, Vol. 48, pp. 309-326, 1993.
        %------------------------------------------------------------------
            disp(' ')
            disp('Available Quadrature Rules for degree-1 meshes:')
            disp(' ')
            disp('  G1 - degree-1, 1-pt Gaussian formula')
            disp('       one quadrature point at each facet centroid')
            disp('       total unique q-points = Nf')
            disp(' ')
            disp('  L1 - degree-1, 3-pt Newton-Cotes formula')
            disp('       total unique q-points = 0.5 * Nf')
            disp(' ')
            disp('  L2 - degree-2, 3-pt Newton-Cotes (Gaussian) formula')
            disp('       total unique q-points = 1.5 * Nf')
            disp(' ')
            disp('  B2 - degree-2 4-pt formula')
            disp('       total unique q-points = 1.5 * Nf')
            disp(' ')
            disp('  G2 - degree-2 3-pt Gaussian formula')
            disp('       total unique q-points = 3 * Nf')
            disp(' ')
            disp('  L3 - degree-3 10-pt Newton-Cotes formula')
            disp('       total unique q-points = 4.5 * Nf')
            disp(' ')
            disp('  B3 - degree-3 7-pt formula')
            disp('       total unique q-points = 3 * Nf')
            disp(' ')
            disp('  O4 - degree-4 10-pt formula')
            disp('       total unique q-points = 4.5 * Nf')
            disp(' ')
            disp('Reccomended Settings for ApproximatePolyhedralModel:')
            disp(' ')
            disp('  Potential field calculation........... L1')
            disp('  Determing In/Out w/ Laplacian......... L1')
            disp('  Acceleration field calculation........ L1')
            disp('  Trajectory Integration < 50 periods... L2')
            disp('  Trajectory Integration < 500 periods.. O4')
            disp(' ')
            disp('  Acceleration calculation becomes invalid for altitudes')
            disp('  at or below the surface mesh resolution.')
        end
    end
end



