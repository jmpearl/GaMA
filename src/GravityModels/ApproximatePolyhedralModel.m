classdef ApproximatePolyhedralModel < handle
% Quadrature-based gravity model for the homogeneous polyhedron
%==========================================================================
% Approximate polyhedral model using numerical quadrature to approximate
% the gravitational fields of a homogenous polyhedron consisting of
% triangular facets. Algorithm and quadrature formulas are described in
% Pearl 2020.
%--------------------------------------------------------------------------
% Pearl, J.M., Hitt, D.L., "A fast quadrature-based gravity model for the
% homogeneous polyhedron," MNRAS, Vol. 492, pp. 420-430, 2020.
%=========================================================================
% Matlab version requirements:
% auto expand -- matlab 2016b+
% vecnorm ------ matlab 2017b+
%==========================================================================
    properties (GetAccess=public,SetAccess=private)
        frame = 'BFF'; % integrator needs to know if bff or inertial x,y,z
        frameId = 1;   % id for iertial
    end

    properties (GetAccess=public)
        GrhoAn;         % G x density x facet area  vector An
        positions;      % quadrature points
        res;
        quadratureRule  % name of quadrature rule 
        numElements     % number of quadrature points
    end
    
    methods
        function obj = ApproximatePolyhedralModel(mesh,Mu,quadratureRule)
            if nargin==3
                obj.initializeFromMesh(mesh,Mu,quadratureRule);
            elseif nargin==2
                obj.initializeFromMesh(mesh,Mu);
            elseif nargin~=0
                error('incorrect number of arguments in constructor')
            end
            
        end        
        
        function initializeFromMesh(obj,mesh,Mu,quadratureRule)

            if nargin~=4 && nargin~=3
                error(' incorrect number of arguments: (SurfaceMesh,Mu,quadratureRule)')
            end
            
            if isstring(mesh) || ischar(mesh)
                mesh = SurfaceMesh(mesh); % we need this to be consistent
            elseif ~isa(mesh,'SurfaceMesh')
                error('first input must be SurfaceMesh or string with surface mesh name')
            end
            
            if nargin ==3
                [obj.positions,Aq] = mesh.createQuadrature();
                obj.quadratureRule = 'mesh-degree';
            else
                [obj.positions,Aq] = mesh.createQuadrature(quadratureRule);
                obj.quadratureRule = quadratureRule;
            end
            obj.GrhoAn = Mu/mesh.volume*Aq;
            obj.res = sqrt(vecnorm(obj.GrhoAn/Mu*mesh.volume,2,2))/4;
            obj.res = max(obj.res,mesh.resolution/100);
            obj.numElements = size(Aq,1);
            
            
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
            
            potential  = zeros(size(p,1),1);
            for i = 1:size(p,1)
                
                rhat  = obj.positions-p(i,:); 
                rhat = rhat./sqrt(rhat(:,1).^2+rhat(:,2).^2+rhat(:,3).^2);
                %rhat(isnan(rhat)|isinf(rhat))=0;
                potential(i,1) = -1/2 * rhat(:)' * obj.GrhoAn(:); 
                
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
            acceleration  = zeros(size(p,1),3);
            
            for i = 1:size(p,1)
                r = obj.positions-p(i,:);
                rinv = 1./max(sqrt(r(:,1).^2+r(:,2).^2+r(:,3).^2),obj.res);
                acceleration(i,1:3) = -rinv'*obj.GrhoAn;
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

            res2 = obj.res.^2;
            laplacian  = zeros(size(p,1),1);
            for i = 1:size(p,1)
                r = obj.positions-p(i,:); % position of q-points relative to sample point P
                
                rOverRinv3 = r./max(r(:,1).^2+r(:,2).^2+r(:,3).^2,res2).^(3/2);
                laplacian(i,1) =  -obj.GrhoAn(:)'*rOverRinv3(:); % Laplacian
            end
            
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
            
            res2 = obj.res.^2;
            gravGradient  = zeros(size(p,1),6);
            
            for i = 1:size(p,1)
                
                r = obj.positions-p(i,:); % position of q-points relative to sample point P
                rinv3 = 1./max(r(:,1).^2+r(:,2).^2+r(:,3).^2,res2).^(3/2);
                rOverR3 = r.*rinv3; % r/r3
                
                gravGradient(i,1:3) =  -obj.GrhoAn(:,1)'*rOverR3;        % 11, 12, 13
                gravGradient(i,4:5) =  -obj.GrhoAn(:,2)'*rOverR3(:,2:3); %     22, 23
                gravGradient(i,6) =  -obj.GrhoAn(:,3)'*rOverR3(:,3);     %         33
            end
        
        end
 
    end
end



