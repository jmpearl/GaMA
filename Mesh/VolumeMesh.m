classdef VolumeMesh
%==========================================================================
% CLEO - Gravity and Mesh Adaptation Libray for Asteroids and Comets
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
        coordinates      % coordinates of vertices
        cells            % cells (indices of faces)
        
        faces2Cells      % cells associated w/ faces
        isBoundaryFace   % bool array indicating if its a BC face
        
        faceFields       % cell array of stored fields
        nodeFields       % cell array of stored fields 
        
        isCurved         % logical to track if mesh is curvilinear
        degree           % polynomial degree of mesh
        volume           % total volume
        surfaceArea      % total surface area
        centroid         % centroid of enclosed volume
        resolution       % mean resolution length scale of mesh
        
        numNodes         % number of nodes for degree>1 meshes
        numVertices      % number of vertices
        numEdges         % number of edges
        numFaces         % number of faces
        numCells         % number of half edges
        numNodeFields    % number of node fields
        numCellFields    % number of face fields
        
    end
    methods
        function obj = VolumeMesh(fileName)
            
            if nargin==1 && contains(fileName,'obj')
                [tempVertices,tempFaces] = obj.readOBJ(fileName);
                obj = obj.initializeFromFaceData(tempVertices,tempFaces);
            elseif nargin > 1
                error('incorrect number of inputs')
            elseif nargin == 1 
                error('unsupported mesh file format')
            end
            
        end
        function obj = initializeGrid(obj,...
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
        
            if length(minCoordinates)==3
                error('minCoordinates must be 1x3')
            end
            if length(maxCoordinates)==3
                error('maxCoordinates must be 1x3')
            end
            if length(numDivisions)==3
                error('numDivisions must be 1x3')
            end
            
            x = linspace(minCoordinates(1),maxCoordinates(1),numDivisions(1));
            y = linspace(minCoordinates(2),maxCoordinates(2),numDivisions(2));
            z = linspace(minCoordinates(3),maxCoordinates(3),numDivisions(3));
            
            [X,Y,Z] = meshgrid(x,y,z);
            
            obj.coordinates = [X(:),Y(:),Z(:)];
            obj.numVertices = size(obj.coordinates,1);
                              
        end
        function obj = delaunayTriangulation(obj)
        % creates DT from the stored vertices
        %------------------------------------------------------------------
            DT = delaunayTriangulation(obj.coordinates(:,1),...
                                       obj.coordinates(:,2),...
                                       obj.coordinates(:,3));
            obj.numCells = size(obj.cells,1);
                              
        end
        function obj = initializeFromFaceData(obj,tempVertices,tempFaces)
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
                        next = obj.next(i);
                        prev = obj.next(next);
                        prevVertex = obj.ends(prev);
                        
                        % half edges pointing to start vertex 
                        candidates = allVertexHalfEdges{prevVertex};
                        for j = 1:length(candidates)
                            
                            candidatej=candidates(j);
                            next = obj.next(candidatej);
                            prev = obj.next(next);
                            pairVertex = obj.ends(prev);
                            
                            if pairVertex == obj.ends(i)
                                obj.pair(i) = candidatej;
                                obj.pair(candidatej) = i;
                                break
                            end
                            
                        end
                    end
                end
                
                obj.degree = 1;
                obj.isCurved = false;
                
                obj.volume = obj.calculateVolume();
                obj.centroid = obj.calculateCentroid();
                obj.surfaceArea = obj.calculateSurfaceArea();
                obj.resolution = obj.calculateResolution(); 
                
                if any(obj.vertexHalfEdges==0)
                    warning('dead vertices detected ... attempting to cull and recover the mesh')
                    obj = obj.cleanDeadVertices();
                end
            
            
        end
        function obj = cleanDeadVertices(obj)
        % cleans out vertices not part of the mesh
        %------------------------------------------------------------------
            vertexIndices = 1:obj.numVertices;
            oldVertexIndices = vertexIndices(obj.vertexHalfEdges~=0);
            newVertexIndices = 1:length(oldVertexIndices);
            for i = 1:length(oldVertexIndices)
                obj = obj.replaceHalfEdgeVertex(oldVertexIndices(i),...
                                                newVertexIndices(i));
            end
            obj.coordinates = obj.coordinates(oldVertexIndices,:);
            obj.vertexHalfEdges = obj.vertexHalfEdges(oldVertexIndices,:);
            
            obj.faces = [obj.ends(1:3:end),...
                         obj.ends(2:3:end),...
                         obj.ends(3:3:end)];
                     
            obj.numVertices = length(obj.vertexHalfEdges);
            obj.numNodes = obj.numVertices;
        end
        
        function centroids = cellCentroids(obj)
        % centroids of cells
        %------------------------------------------------------------------
        % Outputs:
        %   centroids - coordinates of centroids
        %------------------------------------------------------------------
            p1 = obj.coordinates(obj.cells(:,1),:);
            p2 = obj.coordinates(obj.cells(:,2),:); 
            p3 = obj.coordinates(obj.cells(:,3),:); 
            p4 = obj.coordinates(obj.cells(:,4),:); 
            centroids=1/4*(p1+p2+p3+p4); 
            
        end
        
        function obj = addFaceField(obj,fieldData,fieldName)
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
        function obj = addNodeField(obj,fieldData,fieldName)
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
        function obj = clearFields(obj)
        % clears out the field data 
        %------------------------------------------------------------------
        
            obj.numNodeFields = 0;
            obj.numFaceFields = 0;
            
            obj.nodeFields = [];
            obj.faceFields = [];
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
            
            %vtk flags for rectilinear and quadratic triangles
            if obj.degree == 1
                cellType = '5\n';
                numLocalNodes = 3;
            elseif obj.degree == 2 
                cellType = '22\n';
                numLocalNodes = 6;
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
            
            % Write Cells
            %--------------------------------------------------------------
            fprintf(fileID,['CELLS ',num2str(obj.numFaces),' ',num2str(obj.numFaces*(numLocalNodes+1)),'\n']);
            for i = 1:obj.numFaces
                if obj.degree == 1
                    outputi = [num2str([length(obj.faces(i,:)),obj.faces(i,:)-1]),'\n'];
                elseif obj.degree == 2
                    facei = obj.faces(i,[1,3,6,2,5,4])-1;
                    outputi = [num2str([length(obj.faces(i,:)),facei]),'\n'];
                end
                fprintf(fileID,outputi);
            end
            
            % Write CellTypes
            %--------------------------------------------------------------
            fprintf(fileID,['CELL_TYPES ',num2str(obj.numFaces),'\n']);
            for i = 1:obj.numFaces
                fprintf(fileID,cellType);
                
            end
            
            % Write Cell Data
            %--------------------------------------------------------------
            fprintf(fileID,['CELL_DATA ',num2str(obj.numFaces),'\n']);
           
            vtkDataHeader = 'VECTORS faceNormals float\n';
            fprintf(fileID,vtkDataHeader);
            normals = obj.faceNormals();
            for j = 1:obj.numFaces
                outputi = [num2str(normals(j,:),precision),'\n'];
                fprintf(fileID,outputi);
            end
            
            for i = 1:size(obj.faceFields,2)
                if size(obj.faceFields{i}.data,2)==1
                    vtkDataHeader = ['SCALARS ',obj.faceFields{i}.name,' float 1\nLOOKUP_TABLE default\n'];
                elseif size(obj.faceFields{i}.data,2)==3
                    vtkDataHeader = ['VECTORS ',obj.faceFields{i}.name,' float\n'];
                else
                    warning(['skipping faceData{',num2str(i),'} not scalar or vector'])
                end
                fprintf(fileID,vtkDataHeader);
                if size(obj.faceFields{i}.data,1)~=obj.numFaces
                    error(['faceData{',num2str(i),'} inconsistent w/ number of mesh faces'])
                end
                for j = 1:obj.numFaces
                    outputi = [num2str(obj.faceFields{i}.data(j,:),precision),'\n'];
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



