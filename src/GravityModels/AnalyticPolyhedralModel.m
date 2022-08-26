classdef AnalyticPolyhedralModel < handle
% Werner's  analytic model for the gravitational fields of polyhedra
%==========================================================================
%  Gravitational Model derived from the surface definition of an
% irregularly-shaped homogenous body. This implementation is based on
% Werner 1994, 1996, 2017 and assumes constant density and requires a
% surface mesh composed of triangular facets with outward facing normals
% and consistent winding. Variable names inside methods are selected for
% consistency with nomenclature of Werner & Scheeres 1996
%-------------------------------------------------------------------------%
% Werner, R.A., "The Gravitational Potential of a Homogeneous Polyhedron
% or Don't Cut Corners," Celestial Mechanics and Dynamical Astronomy, Vol.
% 59, pp. 253-278, 1994.
%
% Werner, R.A., Scheeres, D.J., Exterior "Gravitation of a Polyhedron
% Dervied and Compared with Harmonic and Mascon Gravitation
% Representations of Aeroid 4769 Castalia," Celestial Mechanics and
% Dynamical Astronomy, Vol. 65, pp. 313–344, 1996.
%
% Werner, R.A., "The Solid Angle Hidden in Polyhedron Gravitation
% Formulations," Journal of Geodesy, Vol. 97, pp. 307-328, 2017.
%==========================================================================
% Matlab version requirements:
% auto expand -- matlab 2016b+
% vecnorm ------ matlab 2017b+
%==========================================================================
    
    properties (GetAccess=public,SetAccess=private)
        frame = 'BFF'; % integrator needs to know if bff or inertial x,y,z
        frameId = 1;   % id for iertial
    end

    properties(GetAccess=public)
        mu;             % gravitational constant per unit volume
        faceCentroids;  % facet centroids
        faceVertices1;  % coordinates of first vertex
        faceVertices2;  % coordinates of second vertex
        faceVertices3;  % coordinates of third vertex
        edgeVertices1;  % coordinates of unique edge start point
        edgeVertices2;  % coordinates of unique edge end point
        edgeLengths;    % edge lengths
        
        % columns of tensors
        E1; 
        E2;
        E3;
        F1;
        F2;
        F3;
        
        numFaces;    % number of faces
        numEdges;    % number of unqiue edges
        numElements; % number of computational elements (numEdges+numFaces)
    end
    
    methods
        function obj = AnalyticPolyhedralModel(mesh, Mu)
            if nargin == 2
                obj.initializeFromMesh(mesh, Mu);
            elseif nargin ~= 0
                error('incorrect number of arguments in constructor')
            end

        end
        function initializeFromMesh(obj, mesh, Mu)
            
            % change to be consistent w/ approx model nomenclature
            % handle variable input 
            if isstring(mesh) || ischar(mesh)
                mesh = SurfaceMesh(mesh);
            elseif ~ isa(mesh,'SurfaceMesh')
                error('  incorrect mesh input: currently support .obj files or SurfaceMesh objects')
            end
            
            % make sure its degree is 1
            mesh.flatten();
            
            [edges,edgeHalfEdges,~] = mesh.edges();
            edgeFaces = ceil(edgeHalfEdges/3);
            n = mesh.faceNormals();
          
            obj.mu = Mu/mesh.volume;
            obj.faceCentroids = mesh.faceCentroids();
            
            obj.faceVertices1 = mesh.coordinates(mesh.faces(:,1),:);
            obj.faceVertices2 = mesh.coordinates(mesh.faces(:,2),:);
            obj.faceVertices3 = mesh.coordinates(mesh.faces(:,3),:);
            
            obj.edgeVertices1 = mesh.coordinates(edges(:,1),:);
            obj.edgeVertices2 = mesh.coordinates(edges(:,2),:);
            obj.edgeLengths   = vecnorm(obj.edgeVertices1 - obj.edgeVertices2,2,2);
            
            % Symmetric Face Dyads stored as 1x6 for each face
            F = [n(:,1).*n(:,1), n(:,1).*n(:,2), n(:,1).*n(:,3), ...
                                 n(:,2).*n(:,2), n(:,2).*n(:,3), ...
                                                 n(:,3).*n(:,3)];
  
            ne1 = cross(obj.edgeVertices1 - obj.edgeVertices2,n(edgeFaces(:,1),:));
            ne1 = ne1./vecnorm(ne1,2,2);
            direction1 = (obj.edgeVertices1 + obj.edgeVertices2)/2 ...
                         - obj.faceCentroids(edgeFaces(:,1),:);
            direction1 = direction1./vecnorm(direction1,2,2);
            dot1 = dot(ne1,direction1,2);
            ne1(dot1<0.0,:) = -1*ne1(dot1<0.0,:);
            

            ne2 = cross(obj.edgeVertices1 - obj.edgeVertices2,n(edgeFaces(:,2),:));
            ne2 = ne2./vecnorm(ne2,2,2);
            direction2 = (obj.edgeVertices1 + obj.edgeVertices2)/2 ...
                         - obj.faceCentroids(edgeFaces(:,2),:);
            direction2 = direction2 ./vecnorm(direction2,2,2);
            dot2 = dot(ne2,direction2,2);
            ne2(dot2<0.0,:) = -1*ne2(dot2<0.0,:); % edges normal 2

            na = n(edgeFaces(:,1),:); % normal of face A
            nb = n(edgeFaces(:,2),:); % normal of face B
            
            % Edge Dyads stored as 1x9 for each face
            E = [na(:,1).*ne1(:,1), na(:,1).*ne1(:,2), na(:,1).*ne1(:,3),...
                 na(:,2).*ne1(:,1), na(:,2).*ne1(:,2), na(:,2).*ne1(:,3),...
                 na(:,3).*ne1(:,1), na(:,3).*ne1(:,2), na(:,3).*ne1(:,3)] +...
                [nb(:,1).*ne2(:,1), nb(:,1).*ne2(:,2), nb(:,1).*ne2(:,3),...
                 nb(:,2).*ne2(:,1), nb(:,2).*ne2(:,2), nb(:,2).*ne2(:,3),...
                 nb(:,3).*ne2(:,1), nb(:,3).*ne2(:,2), nb(:,3).*ne2(:,3)];
                 
           % saving some time not 
           % futzing with indexing later
           % we'll orient this in the transpose to make things a little
           % faster, but is the confusion of doin it this way worth it??
           obj.E1 = E(:,[1,4,7]);
           obj.E2 = E(:,[2,5,8]);
           obj.E3 = E(:,[3,6,9]);
           
           obj.F1 = F(:,[1,2,3]);
           obj.F2 = F(:,[2,4,5]);
           obj.F3 = F(:,[3,5,6]);
           
           obj.numFaces = size(obj.faceCentroids,1);
           obj.numEdges = size(obj.edgeVertices1,1);
           obj.numElements = obj.numFaces + obj.numEdges;
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
            
            for i=1:size(p,1)
                
                % Summation over Faces
                %---------------------------------------------------------%
                rf = obj.faceCentroids-p(i,:); 
                r1 = obj.faceVertices1-p(i,:); 
                r2 = obj.faceVertices2-p(i,:); 
                r3 = obj.faceVertices3-p(i,:);
                
                r1 = r1./sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2);
                r2 = r2./sqrt(r2(:,1).^2+r2(:,2).^2+r2(:,3).^2);
                r3 = r3./sqrt(r3(:,1).^2+r3(:,2).^2+r3(:,3).^2);
                
                % writing out the box product is a bit faster (matlab 2018a)
                omega = atan2(r1(:,1).*(r2(:,2).*r3(:,3)-r2(:,3).*r3(:,2))+...
                    r1(:,2).*(r2(:,3).*r3(:,1)-r2(:,1).*r3(:,3))+...
                    r1(:,3).*(r2(:,1).*r3(:,2)-r2(:,2).*r3(:,1)),...
                    1+sum(r1.*r2+r3.*(r2+r1),2));
                
                %omega = atan2(sum(r1.*cross(r2,r3),2),...
                %             (1+sum(r1.*r2+r2.*r3+r3.*r1,2)));
                
                omega(isinf(omega)|isnan(omega)) = 0; % for singularities
                
                U_f = dot(rf,(obj.F1.*rf(:,1)+obj.F2.*rf(:,2)+obj.F3.*rf(:,3)),2)'*omega;
                
                % Summation over Edges
                %---------------------------------------------------------%
                ra = obj.edgeVertices1-p(i,:); 
                rb = obj.edgeVertices2-p(i,:);
                
                a = sqrt(ra(:,1).^2+ra(:,2).^2+ra(:,3).^2);
                b = sqrt(rb(:,1).^2+rb(:,2).^2+rb(:,3).^2);
                
                %Le = 2*atanh(obj.edgeLengths./(a+b));
                Le = log((a+b+obj.edgeLengths)./(a+b-obj.edgeLengths)); % log version faster
                                                                        % than tanh version
                Le(isinf(Le)|isnan(Le)) = 0;
                Le = real(Le);
                
                % see if we can get away w/ sym tensors
                U_e = dot(ra,(obj.E1.*ra(:,1)+obj.E2.*ra(:,2)+obj.E3.*ra(:,3)),2)'*Le;
                
                % Total Acceleration
                %---------------------------------------------------------%
                potential(i,1) = obj.mu*(U_f-U_e/2)';
                
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
            
            for i=1:size(p,1)
                
                % Summation over Faces
                %---------------------------------------------------------%
                rf = obj.faceCentroids-p(i,:); 
                r1 = obj.faceVertices1-p(i,:); 
                r2 = obj.faceVertices2-p(i,:); 
                r3 = obj.faceVertices3-p(i,:);
                
                r1 = r1./sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2);
                r2 = r2./sqrt(r2(:,1).^2+r2(:,2).^2+r2(:,3).^2);
                r3 = r3./sqrt(r3(:,1).^2+r3(:,2).^2+r3(:,3).^2);
                
                omega = atan2(r1(:,1).*(r2(:,2).*r3(:,3)-r2(:,3).*r3(:,2))+...
                              r1(:,2).*(r2(:,3).*r3(:,1)-r2(:,1).*r3(:,3))+...
                              r1(:,3).*(r2(:,1).*r3(:,2)-r2(:,2).*r3(:,1)),...
                              1+sum(r1.*r2+r3.*(r2+r1),2));
                          
                %omega = atan2(sum(r1.*cross(r2,r3),2),...
                %    (1+sum((r1+r3).*r2+r3.*r1,2)));
                
                omega(isinf(omega)|isnan(omega)) = 0;
                
                % F11 r1 + F12 r2 + F13 r3 ,
                a_f = (obj.F1.*rf(:,1)+obj.F2.*rf(:,2)+obj.F3.*rf(:,3))'*omega;
                
                % Summation over Edges
                %---------------------------------------------------------%
                ra = obj.edgeVertices1-p(i,:); 
                rb = obj.edgeVertices2-p(i,:);
                
                a = sqrt(ra(:,1).^2+ra(:,2).^2+ra(:,3).^2);
                b = sqrt(rb(:,1).^2+rb(:,2).^2+rb(:,3).^2);
                
                %Le = 2*atanh(obj.edgeLengths./(a+b));
                Le = log((a+b+obj.edgeLengths)./(a+b-obj.edgeLengths)); % log version faster
                                                                        % than tanh version
                Le(isinf(Le)|isnan(Le)) = 0;
                Le = real(Le);
                
                a_e = (obj.E1.*ra(:,1)+obj.E2.*ra(:,2)+obj.E3.*ra(:,3))'*Le;
                
                % Total Acceleration
                %---------------------------------------------------------%
                acceleration(i,1:3) = obj.mu*(2*a_f-a_e)'; % 2 comes from solid angle
            end

        end
        function [laplacian] = Laplacian(obj,p)
        % Laplacian of gravitational field. -4piGrho inside, 0 outside
        %------------------------------------------------------------------
        % Inputs:
        %   p --------- [Mx3] array of M sample coordinates
        %------------------------------------------------------------------
        % Outputs:
        %   laplacian - [Mx1] laplacian at sample sites
        %------------------------------------------------------------------
                
            laplacian = zeros(size(p,1),1);
            
            for i=1:size(p,1)
                
                % Summation over Faces
                %---------------------------------------------------------%
                r1 = obj.faceVertices1-p(i,:); 
                r2 = obj.faceVertices2-p(i,:); 
                r3 = obj.faceVertices3-p(i,:);
                
                r1 = r1./sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2);
                r2 = r2./sqrt(r2(:,1).^2+r2(:,2).^2+r2(:,3).^2);
                r3 = r3./sqrt(r3(:,1).^2+r3(:,2).^2+r3(:,3).^2);
                
                omega = atan2(r1(:,1).*(r2(:,2).*r3(:,3)-r2(:,3).*r3(:,2))+...
                              r1(:,2).*(r2(:,3).*r3(:,1)-r2(:,1).*r3(:,3))+...
                              r1(:,3).*(r2(:,1).*r3(:,2)-r2(:,2).*r3(:,1)),...
                              1+sum(r1.*r2+r3.*(r2+r1),2));
                omega(isinf(omega)|isnan(omega)) = 0;
                
                laplacian(i,1) = -2*obj.mu*sum(omega);
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
             
            gravGradient = zeros(size(p,1),6);
            
            for i=1:size(p,1)
                
                % Summation over Faces
                %---------------------------------------------------------%
                r1 = obj.faceVertices1-p(i,:); 
                r2 = obj.faceVertices2-p(i,:); 
                r3 = obj.faceVertices3-p(i,:);
                
                r1 = r1./sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2);
                r2 = r2./sqrt(r2(:,1).^2+r2(:,2).^2+r2(:,3).^2);
                r3 = r3./sqrt(r3(:,1).^2+r3(:,2).^2+r3(:,3).^2);
                
                omega = atan2(r1(:,1).*(r2(:,2).*r3(:,3)-r2(:,3).*r3(:,2))+...
                              r1(:,2).*(r2(:,3).*r3(:,1)-r2(:,1).*r3(:,3))+...
                              r1(:,3).*(r2(:,1).*r3(:,2)-r2(:,2).*r3(:,1)),...
                              1+sum(r1.*r2+r3.*(r2+r1),2));
                omega(isinf(omega)|isnan(omega)) = 0;

                
                % Summation over Edges
                %---------------------------------------------------------%
                ra = obj.edgeVertices1-p(i,:); 
                rb = obj.edgeVertices2-p(i,:);
                
                a = sqrt(ra(:,1).^2+ra(:,2).^2+ra(:,3).^2); 
                b = sqrt(rb(:,1).^2+rb(:,2).^2+rb(:,3).^2);
                %Le = 2*atanh(obj.edgeLengths./(a+b));
                Le = log((a+b+obj.edgeLengths)./(a+b-obj.edgeLengths)); % log version faster
                                                                        % than tanh version
                Le(isinf(Le)|isnan(Le)) = 0;
                Le = real(Le);
               
                % 2x from solid angle done here
                %---------------------------------------------------------%
                gravGradient(i,1:3) = obj.E1'*Le        - 2*(obj.F1'*omega);
                gravGradient(i,4:5) = obj.E2(:,2:3)'*Le - 2*(obj.F2(:,2:3)'*omega);
                gravGradient(i,6)   = obj.E3(:,3)'*Le   - 2*(obj.F3(:,3)'*omega);
            end 
            gravGradient = obj.mu*gravGradient;
        end
        
    end
end   
