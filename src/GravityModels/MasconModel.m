classdef MasconModel < handle % test why handle is slow
% Mass-concentration model 
%--------------------------------------------------------------------------
%
    properties(GetAccess=public,SetAccess=private)
        frame = 'BFF'  % integrator needs to know if bff or inertial x,y,z
        frameId = 1;   % id for iertial
    end
    properties(GetAccess=public)
        mu;          % G x mass of mascon
        coordinates; % coordinates of mascons
        res;         % lengthscale used to prevent singularities
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
                    
                    obj.initializeSimplePacking( mesh, numMascons, Mu);
                    
                elseif strcmp(distributionType,'extendedTetrahedra')
                    
                    numLayers = max(1, round(numMascons/mesh.numFaces));
                    obj.initializeExtendedTetrahedra(mesh,numLayers,Mu);
                    
                else
                    
                    error('invalid distributionType - available options: simplePacking, extendedTetrahedra')
                    
                end
                
                obj.numElements = size(obj.mu,1);
            end
            
        end
        function initializeFromVolumeMesh(obj, volumeMesh,Mu,quadratureMethod)
        % initializes a mascon distribution from VolumeMesh
        %------------------------------------------------------------------
        % Inputs:
        %   volumeMesh ------- initialize volume mesh object
        %   Mu --------------- standard gravitational parameter
        %   quadratureMethod - 'nodes', 'cells','vertices',
        %                      'excludeSurfaceNodes'
        %   centerOfMass ----- true/false quads at COMs?
        %------------------------------------------------------------------
        
            % error checking and defaults
            %--------------------------------------------------------------
            assert(nargin >= 3, 'requires VolumeMesh and stand grav param as inputs');
            assert(isa(volumeMesh,"VolumeMesh"),"first input must be VolumeMesh")
            assert(isnumeric(Mu),"Second inputqm stand grav param must be number")
            
            quadratureMethod = lower(quadratureMethod);
            if nargin == 4
                assert(any(strcmpi(quadratureMethod,{'vertex','node','cell','excludesurface'})), ...
                       "valid quadrature methods 'vertex' 'node' 'cell' 'excludesurface")
            else
                quadratureMethod = 'vertex';
            end

            switch quadratureMethod
                case 'vertex'
                    [cm,vol] = volumeMesh.vertexCentroids();
                case 'cell'
                    [cm,vol] = volumeMesh.cellCentroids();
                case 'node'
                    [cm,vol] = volumeMesh.nodeCentroids();
                case 'excludesurface'
                    [cm,vol] = volumeMesh.vertexCentroidsExcludingBoundaries();
                otherwise
                    error('that aint right')

            end
           
            obj.calculateResolution(vol);
            obj.coordinates = cm;
            obj.mu = vol/sum(vol)*Mu;
            obj.numElements = size(cm,1);

        end
        function initializeSimplePacking(obj, mesh, numMascons, Mu, latticeType)
        % mascons packed in a primitive cubic lattice.
        %------------------------------------------------------------------
        % Inputs:
        %   mesh -------- SurfaceMesh object or string name of obj file
        %   numMascons -- number of mascon (approx)
        %   Mu ---------- standard gravitational parameter for object
        %------------------------------------------------------------------
            
            if isstring(mesh) || ischar(mesh)
                mesh = SurfaceMesh(mesh); % we need this to be consistent
            elseif ~isa(mesh,'SurfaceMesh')
                error('Mesh input must be SurfaceMesh object or string specifing .obj mesh file')
            end
            if nargin == 3
                Mu = mesh.volume;
            elseif nargin < 3 || nargin > 5
                error('incorrect number of inputs')
            end
            
            if nargin < 5
                latticeType = "pc";
            end
        
            insetSurfaceMesh = mesh.offsetSurfaceMesh(-mesh.resolution/2, ...
                                                       mesh.numVertices);
            % spacing so numInternal works out
            if strcmp(latticeType,"bcc")
                ds  = (2*insetSurfaceMesh.volume/numMascons)^(1/3);
            elseif strcmp(latticeType,"fcc")
                ds  = (4*insetSurfaceMesh.volume/numMascons)^(1/3);
            else
                ds  = (insetSurfaceMesh.volume/numMascons)^(1/3);
            end

            maxExtent=max(mesh.coordinates,[],1) + 0.5*ds;
            minExtent=min(mesh.coordinates,[],1) - 0.5*ds; 

            candidates = createLattice(ds,minExtent,maxExtent,latticeType);

            internalNodes = insetSurfaceMesh.isInside(candidates)==1;

            obj.coordinates = candidates(internalNodes,:);
            obj.mu = Mu*ones(size(obj.coordinates,1),1)/size(obj.coordinates,1);
            obj.numElements = size(obj.mu,1);
            volumes = ones(obj.numElements,1)/obj.numElements*mesh.volume;
            obj.calculateResolution(volumes);

        end
        function initializeExtendedTetrahedra(obj, mesh, Mu, numLayers, method)
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
        %   method ---- 'lumpcore' will make the inner most layer sphere
        %                and treat it as a mascon. Requires numLayer >=2
        %------------------------------------------------------------------
            
            % process inputs
            assert(isa(mesh,"SurfaceMesh"),"input 1 must be a surface mesh")
            assert(nargin>=3,"incorrect number of arguments")
            if nargin < 4
                numLayers = 1;
            end
            if nargin < 5 
                method = "standard";
            end

            % call our general function
            [ centroids,volumes] = createExtendedTetrahedralDistribution( mesh, numLayers, method);
            
            % set our class parameters
            obj.coordinates = centroids;
            obj.mu = Mu*volumes/sum(volumes);
            obj.numElements = size(obj.mu,1);

            % get our resolution
            obj.calculateResolution(volumes);
        end
        
        function c = centroid(obj)
        % moves distributions so centroid is at origin
        %------------------------------------------------------------------
           firstMoment = sum(obj.mu,1);
           secondMoment = sum(obj.mu.*obj.coordinates,1);
           c = secondMoment/firstMoment;
        end
        function center(obj)
        % moves distributions so centroid is at origin
        %------------------------------------------------------------------
           obj.coordinates=obj.coordinates-obj.centroid();
        end

        function calculateResolution(obj,volumes)
            obj.res = 0.5*volumes.^(1/3);
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
                
                r = obj.coordinates-p(i,:);
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
                
                r = obj.coordinates-p(i,:);
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
                
                r = obj.coordinates-p(i,:);
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
        function plot(obj,mesh,markerScaleFactor,sliceAxis)
        % plot slice
        %------------------------------------------------------------------
        
            if nargin < 3
                markerScaleFactor = 10;
            end
            if nargin < 4
                sliceAxis = 2;
            end
            FS=16;
            LW=1;
            MS=8;

            thresh = min(max(obj.coordinates,[],1)-min(obj.coordinates,[],1))*0.125;
            disp(thresh)
            figure
            hold on
        
            % if given a mesh plot it 
            if nargin >= 2
                % plot faces in certain range x
                c = mesh.faceCentroids();
                for i = 1:mesh.numFaces
                    if (c(i,sliceAxis) < thresh && c(i,sliceAxis) > -thresh)
                        pts = [mesh.coordinates(mesh.faces(i,:),:);...
                               mesh.coordinates(mesh.faces(i,1),:)];
                        plot3(pts(:,1),pts(:,2),pts(:,3),'k-')
                    end
                end

            end

            c = obj.coordinates;

            resolution = obj.mu.^(1/3);
            MS = resolution./max(resolution)*markerScaleFactor;
            color = 1.0-(MS-min(MS))/max(max(MS)-min(MS),mean(MS)/10000);
            colorVecMax = [1,1,0];
            colorVecMin = [0,0,1];
            for i = 1:obj.numElements
                if (c(i,sliceAxis) < thresh && c(i,sliceAxis) > -thresh)
                    colorVeci = (1-color(i))*colorVecMax + color(i)*colorVecMin;
                    MSi = MS(i);
                    plot3(obj.coordinates(i,1),...
                          obj.coordinates(i,2),...
                          obj.coordinates(i,3),'ko','MarkerFaceColor',colorVeci,'MarkerSize',MSi)
                end
            end
            set(gcf,'Color',[1,1,1]);
            set(gca,'FontSize',FS)
            set(gca,'TickLabelInterpreter','latex')
            daspect([1,1,1])
            
            if sliceAxis==1
                view([1,0,0])
            elseif sliceAxis==2   
                view([0,1,0])
            elseif sliceAxis==3
                view([0,0,1])
            end

            box off
            axis off
        end
    end
end
