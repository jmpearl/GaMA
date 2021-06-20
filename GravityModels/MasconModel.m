classdef MasconModel
% Mass-concentration model 
%--------------------------------------------------------------------------
%
    properties(GetAccess=public)
        mu;          % G x mass of mascon
        positions;   % coordinates of mascons
        
        numElements  % number of computational elements
    end    
    methods
        function [obj] = MasconModel(mesh, Mu, numMascons, distributionType)
        % Provides several methods to distribute mascon throughout interior
        %------------------------------------------------------------------
        % Inputs:
        %   mesh -------------- SurfaceMesh object or string name of obj
        %   distributionType -- type of distribution
        %   N ----------------- number of mascons (rough)
        %   Mu ---------------- standard gravitational parameter for object
        %------------------------------------------------------------------
           
            if nargin==3
                distributionType = 'simplePacking';
            elseif nargin==2
                distributionType = 'extendedTetrahedra';
                numMascons = 0;
            elseif nargin~=0 && nargin~=4
                error('incorrect number of arguments: (SurfaceMesh,Mu,numMascons,distributionType)')
            end
            
            if nargin~=0
                
                if isstring(mesh) || ischar(mesh)
                    mesh = SurfaceMesh(mesh); % we need this to be consistent
                elseif ~isa(mesh,'SurfaceMesh')
                    error('mesh input must be SurfaceMesh object or string specifing .obj mesh file')
                end
                
                if strcmp(distributionType,'simplePacking')
                    
                    ds = (mesh.volume/numMascons)^(1/3);
                    obj = obj.initializeSimplePacking( mesh, ds, Mu);
                    obj.mu = Mu*ones(size(obj.positions,1),1)/size(obj.positions,1);
                    
                elseif strcmp(distributionType,'extendedTetrahedra')
                    
                    numLayers = max(1, round(numMascons/mesh.numFaces));
                    obj = obj.initializeExtendedTetrahedra(mesh,numLayers,Mu);
                    
                else
                    
                    error('invalid distributionType - available options: simplePacking, extendedTetrahedra')
                    
                end
                
                obj.numElements = size(obj.mu,1);
            end
            
        end
        function obj = initializeSimplePacking(obj, mesh, spacing, Mu)
        % mascons packed in a primitive cubic lattice.
        %------------------------------------------------------------------
        % Inputs:
        %   mesh ----- SurfaceMesh object or string name of obj file
        %   spacing -- spacing between mascons
        %   Mu ------- standard gravitational parameter for object
        %------------------------------------------------------------------
            
            if isstring(mesh) || ischar(mesh)
                mesh = SurfaceMesh(mesh); % we need this to be consistent
            elseif ~isa(mesh,'SurfaceMesh')
                error('Mesh input must be SurfaceMesh object or string specifing .obj mesh file')
            end
            if nargin == 3
                Mu = mesh.volume;
            elseif nargin < 3 || nargin > 4
                error('incorrect number of inputs')
            end
            
            bodyMax = max(mesh.coordinates,[],1); % max in each dimension
            bodyMin = min(mesh.coordinates,[],1); % min in each dimension
            
            bodyRange=(bodyMax-bodyMin); % ranges in each direction
            
            %Extend slightly past ranges
            bodyMin = bodyMin - 0.05*bodyRange;
            bodyMax = bodyMax + 0.05*bodyRange;
            
            % Number of elements per dimension
            numSteps = floor((bodyMax-bodyMin)/spacing);
            
            % Correct so that ds is constant
            bodyMax = bodyMin + (numSteps+1)*spacing;
            
            % unique x-y-z coordinates
            x = linspace(bodyMin(1),bodyMax(1),numSteps(1)+2); %#ok<*PROPLC>
            y = linspace(bodyMin(2),bodyMax(2),numSteps(2)+2);
            z = linspace(bodyMin(3),bodyMax(3),numSteps(3)+2);
            
            % Grid up the domain
            [X,Y,Z] = meshgrid(x,y,z);
            candidates = [X(:),Y(:),Z(:)];
            
            % Only keep mascons located within the body
            isKeeper = mesh.isInside(candidates);
            
            % package as Nx3 matrix
            obj.positions = candidates(isKeeper>0.5,:);
            obj.mu = Mu/size(obj.positions,1);
            obj.numElements = size(obj.mu,1);
            
        end
        function obj = initializeExtendedTetrahedra(obj, mesh, numLayers, Mu)
        % Mascons distribution technique of Chanut et al. 2015
        %------------------------------------------------------------------
        % Mascon distribution method proposed by T.G.G Chanut, S. Aljbaae,
        % and V. Carruba in "Mascon Gravitation Model Using A Shaped
        % Polyhedral Source" Monthly Notices of the Royal Astronomical
        % Society 450,3742-3749 (2015)
        %
        % triangular facets of a closed surface mesh are used to create
        % tetrahedra that extend down to the center of the body. Each
        % tetrahedron is sliced and the prisms/tets are approximated as
        % point-masses
        %------------------------------------------------------------------
        % Inputs:
        %   mesh ------ SurfaceMesh object or string name of obj file
        %   numLayers - number of layers
        %   Mu -------- standard gravitational parameter
        %------------------------------------------------------------------
            
            if isstring(mesh) || ischar(mesh)
                mesh = SurfaceMesh(mesh); % we need this to be consistent
            elseif ~isa(mesh,'SurfaceMesh')
                error('Mesh input must be SurfaceMesh object or string specifing .obj mesh file')
            end
            nargin
            if nargin == 3
                Mu = mesh.volume;
            elseif nargin < 3 || nargin > 4
                error('incorrect number of inputs')
            end
            
            % set centroid as origin
            vertices = mesh.coordinates-mesh.centroid;
            
            % three vertices defining each exterior facet
            p1 = vertices(mesh.faces(:,1),:);
            p2 = vertices(mesh.faces(:,2),:);
            p3 = vertices(mesh.faces(:,3),:);
            
            
            obj.positions = zeros(numLayers*mesh.numFaces,3);
            volumes = zeros(numLayers*mesh.numFaces,1);
            
            % outer prism layers
            for i = 1:numLayers-1
                
                % split tetrahedra, define vertices of inner facets
                p4 = (numLayers-i)/numLayers*p1;
                p5 = (numLayers-i)/numLayers*p2;
                p6 = (numLayers-i)/numLayers*p3;
                
                % calculated centroids of 3 tets defining prism
                c1 = (p1+p2+p3+p4)/4;
                c2 = (p2+p3+p4+p6)/4;
                c3 = (p2+p4+p5+p6)/4;
                
                % volume of tets defining prism
                vol1 = abs(dot((p1-p4),cross((p2-p4),(p3-p4),2),2)/6);
                vol2 = abs(dot((p4-p6),cross((p2-p6),(p3-p6),2),2)/6);
                vol3 = abs(dot((p2-p6),cross((p4-p6),(p5-p6),2),2)/6);
                
                % indices for prims layer i
                i1 = mesh.numFaces*(i-1)+1;
                i2 = mesh.numFaces*i;
                
                % sum to get prism properties from tetrahedra
                volumes(i1:i2,1) = vol1+vol2+vol3;
                obj.positions(i1:i2,1:3) = (c1.*vol1+c2.*vol2+c3.*vol3)...
                    ./volumes(i1:i2,1);
                
                p1=p4; p2=p5; p3=p6; % reset iterations
                % inner facet --> outer facet
                
            end
            
            % add in the inner core of tetarahera
            i = numLayers;
            i1 = mesh.numFaces*(i-1)+1;
            i2 = mesh.numFaces*i;
            
            obj.positions(i1:i2,1:3)  = (p1+p2+p3)/4;
            obj.positions = obj.positions + mesh.centroid;
            
            volumes(i1:i2,1) = abs(dot((p1),cross(p2,p3,2),2)/6);
            obj.mu = Mu*volumes/sum(volumes);
            obj.numElements = size(obj.mu,1);
        end
        
        function [potential] = potential(obj,p)
        % Gravitational potential using the negative convention
        %------------------------------------------------------------------
        % Inputs:
        %   p --------- Mx3 array of M sample coordinates
        %------------------------------------------------------------------
        % Outputs:
        %   potential - Mx1 gravitational potential at sample sites
        %------------------------------------------------------------------
               
            potential = zeros(size(p,1),1);
            
            for i = 1:size(p,1)
                
                r = obj.positions-p(i,:);
                rinv = 1./sqrt(r(:,1).^2+r(:,2).^2+r(:,3).^2);
                
                potential(i,1) =-obj.mu'*rinv;
                
            end
            
        end
        function [acceleration] = acceleration(obj,p)
        % Gravitational acceleration 
        %------------------------------------------------------------------
        % Inputs:
        %   p ------------ [Mx3] array of M sample coordinates
        %------------------------------------------------------------------
        % Outputs:
        %   acceleration - [Mx3] gravitational acceleration at sample sites
        %------------------------------------------------------------------
              
            acceleration = zeros(size(p,1),3);
            
            for i = 1:size(p,1)
                
                r = obj.positions-p(i,:);
                rinv3 = 1./(r(:,1).^2+r(:,2).^2+r(:,3).^2).^(3/2);
                
                acceleration(i,1:3) = (obj.mu.*rinv3)'*r;
                
            end
            
        end
        function [laplacian] = Laplacian(obj,p)
        % Laplacian of gravitational field. Always zero for point mass 
        %------------------------------------------------------------------
        % Inputs:
        %   p --------- [Mx3] array of M sample coordinates
        %------------------------------------------------------------------
        % Outputs:
        %   laplacian - [Mx1] laplacian at sample sites
        %------------------------------------------------------------------
          
            laplacian = zeros(size(p,1),1);
        end
        function [gravGradient] = gravityGradient(obj,p)
        % Symmetric gravitational gradient tensor stored as a 1x6
        %------------------------------------------------------------------
        % Inputs:
        %   p ------------ [Mx3] array of M sample coordinates
        %------------------------------------------------------------------
        % Outputs:
        %   gravGradient - [Mx6] acceleration gradient at sample sites
        %------------------------------------------------------------------
         
            gravGradient = zeros(size(p,1),6);
            
            for i = 1:size(p,1)
                
                r = obj.positions-p(i,:);
                rinv = 1./sqrt(r(:,1).^2+r(:,2).^2+r(:,3).^2);
                rhat = r.*rinv; 
                muOverR3 = (obj.mu.*rinv.^3);

                gravGradient(i,1) = muOverR3'*(3*rhat(:,1).*rhat(:,1)-1.0);
                gravGradient(i,2) = muOverR3'*(3*rhat(:,1).*rhat(:,2));
                gravGradient(i,3) = muOverR3'*(3*rhat(:,1).*rhat(:,3));
                gravGradient(i,4) = muOverR3'*(3*rhat(:,2).*rhat(:,2)-1.0);
                gravGradient(i,5) = muOverR3'*(3*rhat(:,2).*rhat(:,3));
                gravGradient(i,6) = muOverR3'*(3*rhat(:,3).*rhat(:,3)-1.0);
                
                % VVU = 3*((mu_rinv3.*r)'*(r.*rinv2))-sum(mu_rinv3)*eye(3);
            end
        end
    end
end
