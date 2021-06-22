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
        
        surfaceMesh      % surfaceMesh definition
        
        faces2Cells      % cells associated w/ faces
        isBoundaryFace   % bool array indicating if its a BC face
        
        faceFields       % cell array of stored fields
        nodeFields       % cell array of stored fields 
        
        hasSurfaceMesh   % bool to track if we have a surf def
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
        function obj = VolumeMesh(surfaceMesh)
            if nargin == 1
                if isa(surfaceMesh,SurfaceMesh)
                    obj.surfaceMesh=surfaceMesh;
                end
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
            
            obj.coordinates = [obj.coordinates;X(:),Y(:),Z(:)];
            obj.numVertices = size(obj.coordinates,1);
                              
        end
        function obj = delaunayTriangulation(obj)
        % creates DT from the stored vertices
        %------------------------------------------------------------------
            DT = delaunayTriangulation(obj.coordinates);
            obj.cells = DT.ConnectivityMap;
            obj.numCells = size(obj.cells,1);
                              
        end
        function obj = clipExternalCells(obj)
        % if cell centroid is external to surface def cull it
        %------------------------------------------------------------------
            if ~ obj.hasSurfaceMesh
                error("requires surface definition to clip")
            end
            
            centroids = obj.cellCentroids;
            isInside = obj.surfaceMesh.isInsideRigorous(centroids);
            obj.cells = obj.cells(isInside,:);

        end
        function obj = clipInternalCells(obj)
        % if cell centroid is internal to surface def cull it
        %------------------------------------------------------------------
               
            if ~ obj.hasSurfaceMesh
                error("requires surface definition to clip")
            end
            
            centroids = obj.cellCentroids;
            isInside = obj.surfaceMesh.isInsideRigorous(centroids);
            obj.cells = obj.cells(~isInside,:);
            
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



