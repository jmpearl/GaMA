classdef SphericalHarmonicModel_OLD
%==========================================================================
% Spheical Harmonic Model graviational model. Potential, acceleration,
% and gravitational gradient calculations are based on Gottlieb 1993. 
% comparison w/ other implementation can be found in Eckman et al 2014. 
% Recursion algothrim for derived Legendre polynomials selected for 
% high degree and order stability based on the work of Lundberg and 
% Schultz 1988. Calculation of harmonic coefficients from polyhedral
% surface meshes and mascon distributions based on Werner 1997.
% -------------------------------------------------------------------------
% Vallado, D.A., "Fundamentals of Astrodynamics and Applications" 
% 4th ed. pp. 543-550 and 596-598, 2013.
%
% Werner, R.A., Spherical Harmonic Coefficients  for the Potential
% of a Constant-Density Polyhedron, Computers & Geoscience,
% Vol. 23, No. 10, pp. 1071-1077, 1997.
%
% Lundberg, J.B., Schultz, B.E., "Recursion Formulas of Legendre Functions
% for Use with Nonsingular Geopotential Models," Journal of Guidance, 
% Vol. 11, No. 1, pp. 31-38, 1988.
%
% Gottlieb, R.G., "Fast Gravity, Gravity Partials, Normalized Gravity,
% Gravity Gradient Torque and Magnetic Field: Derivation, Code and Data,"
% NASA CR-188243, 1993.
%
% Eckman, R.A., Brown, A.J., Adamo, D.R., "Normalization and Implemenation
% of Three Gravitational Acceleration Models," NASA TP-2016-218604, 2014.
%
% BETTER ALGOS FOR CALC OF COEFFS FROM POLYHEDRON
% Jamet O, Thomas E (2004) A linear algorithm for computing the spherical harmonic coefficients of the
% gravitational potential from a constant density polyhedron. In: Proceedings of the second international
% GOCE user workshop, GOCE, The Geoid and Oceanography, ESA-ESRIN, Frascati, Italy, Citeseer,
% pp 8–10
%
% Tsoulis D, Jamet O, Verdun J, Gonindard N (2009) Recursive algorithms for the computation of the poten-
% tial harmonic coefficients of a constant density polyhedron. J Geod 83:925–942
        
%==========================================================================
        
    properties (GetAccess=public)
        % switch C/S to vector format
        C; % cosine coefficients
        S; % sine coefficients
        %Pnn; % derived Legendre Polynomials diagonal
        Ro; % brillouin sphere radius
        Mu; % gravitational parameter
    end
    
    methods (Access=public)
            
%         function obj = SphericalHarmonicModel(Mesh,N,Mu)
%         
%             obj = ConstructFromPolyhedron(obj,Mesh,N,Mu);
%             % obj = ConstructFromMascon(obj,Mesh,N,Mu);
%         end
            
        function [U] = Potential(obj,P)

            % Constants and Preallocation
            %--------------------------------------------------------------
            Cnm = obj.C; % harmonic coeff
            Snm = obj.S; % harmonic coeff
            Ro = obj.Ro; % sphere radius
            
            N = size(Cnm,1)-1; % order of expansion
            r=norm(P); % radial distance
            theta=atan2(P(2),P(1)); % longitude
            phi=asin(P(3)/r); % latitude
            cosphi = cos(phi); % precalc trig functions used multiple times
            sinphi = sin(phi);
            costheta = cos(theta);
            sintheta = sin(theta);
            
            % preallocate to store spherical derivatives for each order (i.e. one
            % allocation for each m-sum in the standard SH equation)
            U = 1;
            gamma=Ro/r;% interative multiplication factor
            
            % recursively calculated terms
            Pnm_v = [1,0,0]; % assoc. Legendre Poly vertical recursion
            cosmtheta = [costheta,1,0];
            sinmtheta = [sintheta,0,0];
            
            A = obj.Mu/r; % Point mass attraction
            
            
            % Zonal Harmonics
            %--------------------------------------------------------------
            for n = 1:N
                
                % preallocate commonly used terms
                n1=n+1;
                
                % recursion of dPnm
                Pnm_v(3) = Pnm_v(1);
                Pnm_v(1) = gamma*(sqrt(((2*n-1)*(2*n+1))/((n)*(n)))*sinphi*Pnm_v(1)...
                    -gamma*sqrt(((2*n+1)/(2*n-3))*((n-1)*(n-1)/((n)*(n))))*Pnm_v(2));
                Pnm_v(2) = Pnm_v(3);
                
                % n < m calc derivatives
                U = U + Pnm_v(1)*Cnm(n1,1);
                
            end
            
            % diagonal recursion of associated Legendre polynomial and derivative
            Pnm_d = gamma*sqrt(3)*cosphi;
            
            
            % Sectoral and Tessoral harmonics
            %--------------------------------------------------------------
            for m = 1:N
                
                % precalc commonly used terms
                m1 = m+1;
                
                % n = m
                U = U + Pnm_d*(Cnm(m1,m1)*cosmtheta(1)+Snm(m1,m1)*sinmtheta(1));
                
                % initialize vertical recursion
                Pnm_v = [Pnm_d,0,0];
                
                for n = m+1:N
                    
                    % recursion of Pnm - vertically
                    n1=n+1;
                    Pnm_v(3) = Pnm_v(1);
                    Pnm_v(1) = gamma*((sqrt(((2*n-1)*(2*n+1))/((n+m)*(n-m))))*sinphi*Pnm_v(1)...
                        -gamma*(sqrt(((2*n+1)/(2*n-3))*((n-m-1)*(n+m-1)/((n+m)*(n-m)))))*Pnm_v(2));
                    Pnm_v(2) = Pnm_v(3);
                    
                    U = U + Pnm_v(1)*(Cnm(n1,m1)*cosmtheta(1)+Snm(n1,m1)*sinmtheta(1));
                    
                end
                
                % diagonal recursion of associated Legendre polynomial and derivative
                Pnm_d = sqrt(1+0.5/(m1))*gamma*cosphi*Pnm_d;
                cosmtheta(3)=cosmtheta(1); sinmtheta(3)=sinmtheta(1);
                cosmtheta(1) = 2*costheta*cosmtheta(1)-cosmtheta(2); cosmtheta(2)=cosmtheta(3);
                sinmtheta(1) = 2*costheta*sinmtheta(1)-sinmtheta(2); sinmtheta(2)=sinmtheta(3);
            end
            
            U=U*A;
            
        end
        
        function [a] = Acceleration(obj,P)
        %------------------------------------------------------------------
        % returns the gravitational acceleration in cartesian coordinates
        % based on an input cartesian coordinate (x,y,z). Recursive 
        % formulation designed to be fast and memory light.
        %------------------------------------------------------------------
        % Inputs:
        %   P - position in 3d cartesian (x,y,z)
        %------------------------------------------------------------------
        % Outputs:
        %   a - acceleration in 3d cartesian (ax,ay,az)
        %------------------------------------------------------------------
        
            % Constants and Preallocation
            %--------------------------------------------------------------
            Cnm = obj.C; % harmonic coeff
            Snm = obj.S; % harmonic coeff
            
            N = size(Cnm,1)-1; % order of expansion
            r = norm(P); % radial distance
            theta= atan2(P(2),P(1)); % longitude
            phi= asin(P(3)/r); % latitude
            cosphi = cos(phi); % precalc trig functions used multiple times
            sinphi = sin(phi);
            costheta = cos(theta);
            sintheta = sin(theta);
            gamma = obj.Ro/r;
            dUdr = 1;
            dUdtheta = 0;
            dUdphi = 0;
            A = obj.Mu/r/r; % Point mass attraction
            
            % recursively calculated terms
            Pnm_v = [1,0,0]; % assoc. Legendre Poly vertical recursion
            dPnm_v = [0,0,0]; % d/dphi assoc. Legendre Poly vertical recursion
            cosmtheta = [costheta,1,0]; % cos(m*theta)
            sinmtheta = [sintheta,0,0]; % sin(m*theta)
            
            % Zonal Harmonics
            %--------------------------------------------------------------
            for n = 1:N
                
                % preallocate commonly used terms
                n1=n+1;
                c1 = sqrt(((2*n-1)*(2*n+1))/((n)*(n)));
                c2 = sqrt(((2*n+1)/(2*n-3))*((n-1)*(n-1)/((n)*(n))));
                
                % recursion of dPnm/dphi
                % (this formulation is faster than using 2x1 iteration... go figure)
                dPnm_v(3) = dPnm_v(1);
                dPnm_v(1) = gamma*(c1*(sinphi*dPnm_v(1)+cosphi*Pnm_v(1))-gamma*c2*dPnm_v(2));
                dPnm_v(2) = dPnm_v(3);
                
                % recursion of dPnm
                Pnm_v(3) = Pnm_v(1);
                Pnm_v(1) = gamma*(c1*sinphi*Pnm_v(1)-gamma*c2*Pnm_v(2));
                Pnm_v(2) = Pnm_v(3);
                
                % n < m calc derivatives
                dUdr = dUdr + n1*Pnm_v(1)*Cnm(n1,1);
                dUdphi = dUdphi + dPnm_v(1)*Cnm(n1,1);
                
            end
            
            % diagonal recursion of associated Legendre polynomial and derivative
            c3 = sqrt(3);
            dPnm_d = -gamma*c3*sinphi;  % this was the issue.........
            Pnm_d = gamma*c3*cosphi;
            
            
            % Sectoral and Tessoral harmonics
            %--------------------------------------------------------------
            for m = 1:N
                
                % precalc commonly used terms
                m1 = m+1; % commonly used index
                A1 = Cnm(m1,m1)*cosmtheta(1)+Snm(m1,m1)*sinmtheta(1); % used twice
                
                % n=m calc derivatives
                dUdr = dUdr + m1*Pnm_d*A1;
                dUdtheta = dUdtheta + m*Pnm_d*(Snm(m1,m1)*cosmtheta(1)-Cnm(m1,m1)*sinmtheta(1));
                dUdphi = dUdphi + dPnm_d*A1;
                
                % initialize vertical recursion
                %disp(Pnm_d)
                %disp(dPnm_d)
                %disp(' ')
                Pnm_v = [Pnm_d,0,0];
                dPnm_v = [dPnm_d,0,0];
                
                for n = m+1:N
                    
                    % preallocate commonly used terms
                    n1=n+1;
                    c1 = sqrt(((2*n-1)*(2*n+1))/((n+m)*(n-m)));
                    c2 = sqrt(((2*n+1)/(2*n-3))*((n-m-1)*(n+m-1)/((n+m)*(n-m))));
                    A1 = Cnm(n1,m1)*cosmtheta(1)+Snm(n1,m1)*sinmtheta(1);
                    
                    % recursion of dPnm/dphi
                    dPnm_v(3) = dPnm_v(1);
                    dPnm_v(1) = gamma*(c1*(sinphi*dPnm_v(1)+cosphi*Pnm_v(1))-gamma*c2*dPnm_v(2));
                    dPnm_v(2) = dPnm_v(3);
                    
                    % recursion of dPnm
                    Pnm_v(3) = Pnm_v(1);
                    Pnm_v(1) = gamma*(c1*sinphi*Pnm_v(1)-gamma*c2*Pnm_v(2));
                    Pnm_v(2) = Pnm_v(3);

                    % n < m calc derivatives
                    dUdr = dUdr + n1*Pnm_v(1)*A1;
                    dUdphi = dUdphi + dPnm_v(1)*A1;
                    dUdtheta = dUdtheta + m*Pnm_v(1)*(Snm(n+1,m+1)*cosmtheta(1)-Cnm(n+1,m+1)*sinmtheta(1));
                    
                end
                %disp([cosmtheta(1)*cosphi^m,sinmtheta(1)*cosphi^m,Pnm_d/cosphi^m])
                % diagonal recursion of associated Legendre polynomial and derivative
                c3 = sqrt(1+0.5/m1); % used twice
                dPnm_d = c3*gamma*(cosphi*dPnm_d-sinphi*Pnm_d);
                Pnm_d = c3*gamma*cosphi*Pnm_d;
                cosmtheta(3)=cosmtheta(1); sinmtheta(3)=sinmtheta(1);
                cosmtheta(1) = 2*costheta*cosmtheta(1)-cosmtheta(2); cosmtheta(2)=cosmtheta(3);
                sinmtheta(1) = 2*costheta*sinmtheta(1)-sinmtheta(2); sinmtheta(2)=sinmtheta(3);
            end
    
            % Convert to Cartesian
            %--------------------------------------------------------------
            a=[0,0,0];
            a(1) = A*(-dUdr/r*P(1)-sintheta/cosphi*dUdtheta-costheta*sinphi*dUdphi);
            a(2) = A*(-dUdr/r*P(2)+costheta/cosphi*dUdtheta-sintheta*sinphi*dUdphi);
            a(3) = A*(-dUdr/r*P(3)+cosphi*dUdphi);
            
        end
        
        function [U] = Potential_gottlieb(obj,P)
        %------------------------------------------------------------------
        % returns the gravitational acceleration in cartesian coordinates
        % based on an input cartesian coordinate (x,y,z). Recursive 
        % formulation designed to be fast and memory light.
        %------------------------------------------------------------------
        % Inputs:
        %   P - position in 3d cartesian (x,y,z)
        %------------------------------------------------------------------
        % Outputs:
        %   a - acceleration in 3d cartesian (ax,ay,az)
        %------------------------------------------------------------------
        
            % Constants and Preallocation
            %--------------------------------------------------------------
            Cnm = obj.C; % harmonic coeff
            Snm = obj.S; % harmonic coeff
            
           
            N = size(Cnm,1)-1; % order of expansion
            r = norm(P); % radial distance
            rinv = 1/r;
            
            U = 1;
            
            drdx=P*rinv;
            gamma = obj.Ro*rinv; 
            A = obj.Mu*rinv; % Point mass attraction
            Pnm_v = [1,0,0]; % assoc. Legendre Poly vertical recursion

            % Zonal Harmonics
            %--------------------------------------------------------------
            for n = 1:N
                
                % preallocate commonly used terms
                n1=n+1;
                c1 = sqrt(((2*n-1)*(2*n+1))/((n)*(n)));
                c2 = sqrt(((2*n+1)/(2*n-3))*((n-1)*(n-1)/((n)*(n))));
                
                % recursion of dPnm
                Pnm_v(3) = Pnm_v(1);
                Pnm_v(1) = gamma*(c1*drdx(3)*Pnm_v(1)-c2*gamma*Pnm_v(2));
                Pnm_v(2) = Pnm_v(3);

                U = U + Pnm_v(1)*Cnm(n1,1);

            end
            
            
            Cm_old = 1.0;
            Cm = drdx(1);
            Sm = drdx(2);
            Pnm_d = sqrt(3)*gamma;
            
            % Sectoral and Tessoral harmonics
            %--------------------------------------------------------------
            for m = 1:N
                
                % precalc commonly used terms
                m1 = m+1; % commonly used index
            
                % n=m calc derivatives
                U = U + Pnm_d*(Cnm(m1,m1)*Cm+Snm(m1,m1)*Sm);
                
                % initialize vertical recursion
                Pnm_v(1) = Pnm_d;
                Pnm_v(2) = 0.0;
                Pnm_v(3) = 0.0;

                for n = m+1:N
                    
                    % preallocate commonly used terms
                    n1=n+1;

                    % recursion of dPnm
                    Pnm_v(3) = Pnm_v(1);
                    Pnm_v(1) = gamma*(sqrt(((2*n-1)*(2*n+1))/((n+m)*(n-m)))*drdx(3)*Pnm_v(1)-gamma*sqrt(((2*n+1)/(2*n-3))*((n-m-1)*(n+m-1)/((n+m)*(n-m))))*Pnm_v(2));
                    Pnm_v(2) = Pnm_v(3);
                    
                    % n < m calc derivatives
                    U = U + Pnm_v(1)*(Cnm(n1,m1)*Cm+Snm(n1,m1)*Sm);
                   
                end

                
                Cm_old = Cm;
                Cm = drdx(1)*Cm-drdx(2)*Sm;
                Sm = drdx(1)*Sm+drdx(2)*Cm_old;
                Pnm_d = gamma*sqrt(1+0.5/m1)*Pnm_d;
            end
            U = U*A;
        end 
        
        function [a] = Acceleration_gottlieb(obj,P)
        %------------------------------------------------------------------
        % returns the gravitational acceleration in cartesian coordinates
        % based on an input cartesian coordinate (x,y,z). Recursive 
        % formulation designed to be fast and memory light.
        %------------------------------------------------------------------
        % Inputs:
        %   P - position in 3d cartesian (x,y,z)
        %------------------------------------------------------------------
        % Outputs:
        %   a - acceleration in 3d cartesian (ax,ay,az)
        %------------------------------------------------------------------
        
            % Constants and Preallocation
            %--------------------------------------------------------------
            Cnm = obj.C; % harmonic coeff
            Snm = obj.S; % harmonic coeff
            
            N = size(Cnm,1)-1; % order of expansion
            r = norm(P); % radial distance
            rinv = 1/r;
            drdx=P*rinv;
            
            dUdr = 1;
            dUdeps = 0;
             a=[0,0,0];
            gamma = obj.Ro*rinv; 
            A = obj.Mu*rinv^2; % Point mass attraction
            Pnm_v = [1,0,0]; % assoc. Legendre Poly vertical recursion

            % Zonal Harmonics
            %--------------------------------------------------------------
            for n = 1:N
                
                % preallocate commonly used terms
                n1=n+1;
                c1 = sqrt(((2*n-1)*(2*n+1))/((n)*(n)));
                c2 = sqrt(((2*n+1)/(2*n-3))*((n-1)*(n-1)/((n)*(n))));
                
                % recursion of dPnm
                Pnm_v(3) = Pnm_v(1);
                Pnm_v(1) = gamma*(c1*drdx(3)*Pnm_v(1)-c2*gamma*Pnm_v(2));
                Pnm_v(2) = Pnm_v(3);

                dUdr = dUdr + n1*Pnm_v(1)*Cnm(n1,1);

            end
            
            D=2.0;
            Cm_old = 1.0;
            Sm_old = 0.0;
            Cm = drdx(1);
            Sm = drdx(2);
            Pnm_d = sqrt(3)*gamma;
            
            % Sectoral and Tessoral harmonics
            %--------------------------------------------------------------
            for m = 1:N
                
                % precalc commonly used terms
                m1 = m+1; % commonly used index
            
                % n=m calc derivatives
                dUdr = dUdr + (m+m1)*Pnm_d*(Cnm(m1,m1)*Cm+Snm(m1,m1)*Sm);
                a(1) = a(1) + Pnm_d*m*(Cnm(m1,m1)*Cm_old+Snm(m1,m1)*Sm_old);
                a(2) = a(2) + Pnm_d*m*(Snm(m1,m1)*Cm_old-Cnm(m1,m1)*Sm_old);
                dUdeps = dUdeps + sqrt((2*m)/D)*Pnm_d*(Cnm(m1,m)*Cm_old+Snm(m1,m)*Sm_old);
                
                % initialize vertical recursion
                Pnm_v(1) = Pnm_d;
                Pnm_v(2) = 0.0;
                Pnm_v(3) = 0.0;

                for n = m+1:N
                    
                    % preallocate commonly used terms
                    n1=n+1;

                    % recursion of dPnm
                    Pnm_v(3) = Pnm_v(1);
                    Pnm_v(1) = gamma*(sqrt(((2*n-1)*(2*n+1))/((n+m)*(n-m)))*drdx(3)*Pnm_v(1)-gamma*sqrt(((2*n+1)/(2*n-3))*((n-m-1)*(n+m-1)/((n+m)*(n-m))))*Pnm_v(2));
                    Pnm_v(2) = Pnm_v(3);
                    
                    % n < m calc derivatives
                    dUdr = dUdr + (m+n+1)*Pnm_v(1)*(Cnm(n1,m1)*Cm+Snm(n1,m1)*Sm);
                    dUdeps = dUdeps +  sqrt((n-m+1)*(n+m)/D)*Pnm_v(1)*(Cnm(n1,m)*Cm_old+Snm(n1,m)*Sm_old);
                    a(1) = a(1) + Pnm_v(1)*m*(Cnm(n1,m1)*Cm_old+Snm(n1,m1)*Sm_old);
                    a(2) = a(2) + Pnm_v(1)*m*(Snm(n1,m1)*Cm_old-Cnm(n1,m1)*Sm_old);
                end

                D=1.0; % to handle horizontal ratio of Nnm w/ kronecker
                Cm_old = Cm;
                Sm_old = Sm;
                Cm = drdx(1)*Cm_old-drdx(2)*Sm_old;
                Sm = drdx(1)*Sm_old+drdx(2)*Cm_old;
                Pnm_d = gamma*sqrt(1+0.5/m1)*Pnm_d;
            end
  
            dUdr=dUdr+drdx(3)*dUdeps;
            a(1) = A*(a(1)-dUdr*drdx(1));
            a(2) = A*(a(2)-dUdr*drdx(2));
            a(3) = A*(a(3)-dUdr*drdx(3)+dUdeps);   
        end       
        
        function [V2U] = Laplacian(obj,P)
            V2U=0;
        end
        %function [VVU] = GravityGradient(obj,P) 
    
        function obj = ConstructFromPolyhedron(obj,mesh,N,Mu)
            
            % Initializes harmonic coefficients
            %--------------------------------------------------------------
            if isa(mesh,'SurfaceMesh')
                M = mesh;
                assert(mesh.degree==1);
            elseif isstring(mesh) || ischar(mesh)
                M = SurfaceMesh(mesh);
            else
                error('  incorrect mesh input: currently support .obj files or SurfaceMesh objects')
            end
            
            pts = M.coordinates; f = M.faces;
            
            clear M
            
            a = max(calc_Mag(pts)); % radius of expansion sphere
            pts = pts/a; % save some compute time later
            
            C = zeros(N+1,N+1);%zeros((N+1)*(N+2)/2,2);
            S = zeros(N+1,N+1);%zeros((N+1)*(N+2)/2,2);
            %obj.Pnn = zeros(N+1,1); obj.Pnn(1)=1; obj.Pnn(2)=-sqrt(3);
            N_store = N^2/2+5/2*N+3; % storage required 
            indices = linspace(1,N_store,N_store)'; % for convienence later
            alpha_d = zeros(N_store,2); % storage for diagonal recursion
            beta_d = zeros(N_store,2); % storage for diagonal recursion
            alpha_v = zeros(N_store,3); % storage for vertical recursion
            beta_v = zeros(N_store,3); % storage for vertical recursion
            ijk = zeros(N_store,3); % stores exponents of X,Y,Z (vertical)
            ijkd = zeros(N_store,3); % stores exponents of X,Y,Z (diagonal)
            fXYZ_v = ones(N_store,3); % factorial of X,Y,Z calculated recursively
            fXYZ_d = ones(N_store,3); % factorial of X,Y,Z calculated recursively
            DETJ=0; % intialize volume
            
            
            % Calculate contribution from each face.
            for i=1:length(f(:,1))
                
                p = pts(f(i,:),:); % vertices of triangular face (0,0,0) other vertex
                detJ = det(p); % determinant = 6*V_tet
                DETJ = DETJ+detJ; % sum determinant to get total volume
                
                % x/a, y/a, z/a, and r/a in Werner 1997
                xa = p(:,1); % x coord of vertices
                ya = p(:,2); % y coord of vertices
                za = p(:,3); % z coord of vertices
                ra = [p(1,:)*p(1,:)';...
                    2*p(1,:)*p(2,:)';...
                    2*p(1,:)*p(3,:)';...
                    p(2,:)*p(2,:)';...
                    2*p(2,:)*p(3,:)';...
                    p(3,:)*p(3,:)']; % [xx;xy;xz;yy;yz;zz]
                
                cd = 1/sqrt(3); % first coeff in diagonal recursion
                
                alpha_d(:,:) = 0.0; % storage for diagonal recursion
                beta_d(:,:) = 0.0; % storage for diagonal recursion
                alpha_v(:,:) = 0.0; % storage for vertical recursion
                beta_v(:,:) = 0.0; % storage for vertical recursion
                
                ijk(:,:) = 0.0; % exponents of x,y,z in monomials (vertical)
                ijkd(:,:) = 0.0; % exponents of x,y,z in monomials (diagonal)
                ijk(1:3,1:3) = eye(3); % start with Ax + By + Cz
                ijkd(1:3,1:3) = eye(3); % start with Ax + By + Cz
                
                fXYZ_v(:,:) = 1;
                fXYZ_d(:,:) = 1;
                fn3 = 24;
                fn3_d = 24;
                % initialize recursion for facet f
                alpha_d(1,1) = 1;
                alpha_v(1,2) = 1;
                beta_d(1,1) = 0;
                beta_v(1,2) = 1;
                N_ijk(1:3) = 1;
                M_ijk(1:2) = 1;
                C(1,1) = C(1,1) + detJ/6.0; % first coeff, really a sanity check
                
                
                
                for n=1:N %---Zonal-Harmonics---
%                     if n>1
%                         n1=n+1;
%                         obj.Pnn(n+1) = -sqrt(1+0.5/(n1))*obj.Pnn(n); % correct of Nnm
%                     end
                    % coefficients for recursion, calculate once
                    c1 = (2*n-1)*sqrt((2*n-1)/((2*n+1)*n^2));
                    c2 = sqrt(((2*n-3)*(n-1)^2)/((2*n+1)*n^2));
                    
                    idx = indices(1:N_ijk(2));
                    idy = idx+n-ijk(idx,1)+1;
                    idz = idy + 1;
                    
                    % monomials degree
                    alpha_v(idx,3) =  c1*za(1).*alpha_v(idx,2);
                    alpha_v(idy,3) =  alpha_v(idy,3)+c1*za(2).*alpha_v(idx,2);
                    alpha_v(idz,3) =  alpha_v(idz,3)+c1*za(3).*alpha_v(idx,2);
                    
                    idxx = indices(1:N_ijk(1));
                    idxy = idxx + n-ijk(idxx,1) + 1;
                    idxz = idxy + 1;
                    idyy = idxx + 2*(n-ijk(idxx,1)) + 3;
                    idyz = idyy + 1;
                    idzz = idyz + 1;
                    
                    N_ijk(1) = N_ijk(2); % update here avoids extra storage
                    N_ijk(2) = (n^2/2+3/2*n+1);
                    
                    alpha_v(idxx,3) =  alpha_v(idxx,3)-c2*ra(1).*alpha_v(idxx,1);
                    alpha_v(idxy,3) =  alpha_v(idxy,3)-c2*ra(2).*alpha_v(idxx,1);
                    alpha_v(idxz,3) =  alpha_v(idxz,3)-c2*ra(3).*alpha_v(idxx,1);
                    alpha_v(idyy,3) =  alpha_v(idyy,3)-c2*ra(4).*alpha_v(idxx,1);
                    alpha_v(idyz,3) =  alpha_v(idyz,3)-c2*ra(5).*alpha_v(idxx,1);
                    alpha_v(idzz,3) =  alpha_v(idzz,3)-c2*ra(6).*alpha_v(idxx,1);
                    
                    % Integrate monimals
                    C(n+1,1) = C(n+1,1) + detJ/fn3...
                        *sum(fXYZ_v(1:N_ijk(2),1).*...
                        fXYZ_v(1:N_ijk(2),2).*...
                        fXYZ_v(1:N_ijk(2),3).*...
                        alpha_v(1:N_ijk(2),3));
                    
                    
                    % increment var for next iter.
                    N_ijk(1) = N_ijk(2);
                    
                    ijk(1:N_ijk(2),1) = ijk(1:N_ijk(2),1)+1;
                    ijk(N_ijk(2)+1:N_ijk(2)+n+2,2) =  ijk(N_ijk(2)+1:N_ijk(2)+n+2,2) + fliplr(indices(1:n+2)')'-1;
                    ijk(N_ijk(2)+1:N_ijk(2)+n+2,3) =  ijk(N_ijk(2)+1:N_ijk(2)+n+2,3) + indices(1:n+2)-1;
                    
                    % handle integration factorials recursively
                    fn3 = fn3*(n+4);
                    fXYZ_v(1:N_ijk(1),1) = fXYZ_v(1:N_ijk(1),1).*ijk(1:N_ijk(1),1);
                    fXYZ_v(N_ijk(2)+1:N_ijk(2)+n,2) =  ijk(N_ijk(2)+1:N_ijk(2)+n,2).*fXYZ_v(N_ijk(2)-n:N_ijk(2)-1,2);
                    fXYZ_v(N_ijk(2)+3:N_ijk(2)+n+2,3) =  ijk(N_ijk(2)+3:N_ijk(2)+n+2,3).*fXYZ_v(N_ijk(2)-n+1:N_ijk(2),3);
                    
                    alpha_v(:,1:2) = alpha_v(:,2:3);
                    alpha_v(:,3) = 0;
                end
                
                
                for m=1:N %---Sectoral-Harmonics---
                    
                    
                    idx = indices(1:M_ijk);
                    idy = idx+m-ijkd(idx,1)+1;
                    idz = idy + 1;
                    
                    % (-1)^m and diagonal recursion of Pnm
                    alpha_d(idx,2) =  alpha_d(idx,2)+cd*(xa(1).*alpha_d(idx,1)-ya(1).*beta_d(idx,1));
                    alpha_d(idy,2) =  alpha_d(idy,2)+cd*(xa(2).*alpha_d(idx,1)-ya(2).*beta_d(idx,1));
                    alpha_d(idz,2) =  alpha_d(idz,2)+cd*(xa(3).*alpha_d(idx,1)-ya(3).*beta_d(idx,1));
                    
                    beta_d(idx,2) =  beta_d(idx,2)+cd*(ya(1).*alpha_d(idx,1)+xa(1).*beta_d(idx,1));
                    beta_d(idy,2) =  beta_d(idy,2)+cd*(ya(2).*alpha_d(idx,1)+xa(2).*beta_d(idx,1));
                    beta_d(idz,2) =  beta_d(idz,2)+cd*(ya(3).*alpha_d(idx,1)+xa(3).*beta_d(idx,1));
                    
                    
                    M_ijk = (m^2/2+3/2*m+1); % update here avoids extra storage
                    
                    C(m+1,m+1) = C(m+1,m+1) + detJ/fn3_d...
                        *sum(fXYZ_d(1:M_ijk,1).*...
                        fXYZ_d(1:M_ijk,2).*...
                        fXYZ_d(1:M_ijk,3).*alpha_d(1:M_ijk,2));
                    
                    S(m+1,m+1) = S(m+1,m+1) + detJ/fn3_d...
                        *sum(fXYZ_d(1:M_ijk,1).*...
                        fXYZ_d(1:M_ijk,2).*...
                        fXYZ_d(1:M_ijk,3).*beta_d(1:M_ijk,2));
                    
                    % initialize vertical recursion
                    
                    ijkd(1:M_ijk,1) = ijkd(1:M_ijk,1)+1;
                    ijkd(M_ijk+1:M_ijk+m+2,2) =  ijkd(M_ijk+1:M_ijk+m+2,2) + fliplr(indices(1:m+2)')'-1;
                    ijkd(M_ijk+1:M_ijk+m+2,3) =  ijkd(M_ijk+1:M_ijk+m+2,3) + indices(1:m+2)-1;
                    
                    fn3_d = fn3_d*(m+4);
                    fXYZ_d(1:M_ijk,1) = fXYZ_d(1:M_ijk,1).*ijkd(1:M_ijk,1);
                    fXYZ_d(M_ijk+1:M_ijk+m,2) =  ijkd(M_ijk+1:M_ijk+m,2).*fXYZ_d(M_ijk-m:M_ijk-1,2);
                    fXYZ_d(M_ijk+3:M_ijk+m+2,3) =  ijkd(M_ijk+3:M_ijk+m+2,3).*fXYZ_d(M_ijk-m+1:M_ijk,3);
                    
                    
                    N_ijk(:) = 1;
                    N_ijk(2) = M_ijk;
                    ijk = ijkd;
                    fXYZ_v(:,:)=fXYZ_d(:,:);
                    fn3 = fn3_d;
                    
                    % initialize vertical recursion
                    alpha_v(:,:) = 0;
                    beta_v(:,:) = 0;
                    alpha_v(:,2)= alpha_d(:,2);
                    beta_v(:,2) = beta_d(:,2);
                    
                    % reset for next diagonal recursion
                    alpha_d(:,1) = alpha_d(:,2);
                    alpha_d(:,2) = 0;
                    
                    beta_d(:,1) = beta_d(:,2);
                    beta_d(:,2) = 0;
                    
                    cd = (2*m+1)/sqrt(2*(m+1)*(2*m+3));
                    
                    for n=m+1:N
                        
                        % coefficients for recursion, calculate once
                        c1 = (2*n-1)*sqrt((2*n-1)/((2*n+1)*(n+m)*(n-m)));
                        c2 = sqrt(((2*n-3)*(n+m-1)*(n-m-1))/((2*n+1)*(n+m)*(n-m)));
                        
                        
                        idx = indices(1:N_ijk(2));
                        idy = idx+n-ijk(idx,1)+1;
                        idz = idy + 1;
                        
                        % monomials degree
                        alpha_v(idx,3) =  c1*za(1).*alpha_v(idx,2);
                        alpha_v(idy,3) =  alpha_v(idy,3)+c1*za(2).*alpha_v(idx,2);
                        alpha_v(idz,3) =  alpha_v(idz,3)+c1*za(3).*alpha_v(idx,2);
                        
                        beta_v(idx,3) =  c1*za(1).*beta_v(idx,2);
                        beta_v(idy,3) =  beta_v(idy,3)+c1*za(2).*beta_v(idx,2);
                        beta_v(idz,3) =  beta_v(idz,3)+c1*za(3).*beta_v(idx,2);
                        
                        idxx = indices(1:N_ijk(1));
                        idxy = idxx + n-ijk(idxx,1) + 1;
                        idxz = idxy + 1;
                        idyy = idxx + 2*(n-ijk(idxx,1)) + 3;
                        idyz = idyy + 1;
                        idzz = idyz + 1;
                        
                        N_ijk(1) = N_ijk(2);
                        N_ijk(2) = (n^2/2+3/2*n+1);
                        
                        alpha_v(idxx,3) =  alpha_v(idxx,3)-c2*ra(1).*alpha_v(idxx,1);
                        alpha_v(idxy,3) =  alpha_v(idxy,3)-c2*ra(2).*alpha_v(idxx,1);
                        alpha_v(idxz,3) =  alpha_v(idxz,3)-c2*ra(3).*alpha_v(idxx,1);
                        alpha_v(idyy,3) =  alpha_v(idyy,3)-c2*ra(4).*alpha_v(idxx,1);
                        alpha_v(idyz,3) =  alpha_v(idyz,3)-c2*ra(5).*alpha_v(idxx,1);
                        alpha_v(idzz,3) =  alpha_v(idzz,3)-c2*ra(6).*alpha_v(idxx,1);
                        
                        beta_v(idxx,3) =  beta_v(idxx,3)-c2*ra(1).*beta_v(idxx,1);
                        beta_v(idxy,3) =  beta_v(idxy,3)-c2*ra(2).*beta_v(idxx,1);
                        beta_v(idxz,3) =  beta_v(idxz,3)-c2*ra(3).*beta_v(idxx,1);
                        beta_v(idyy,3) =  beta_v(idyy,3)-c2*ra(4).*beta_v(idxx,1);
                        beta_v(idyz,3) =  beta_v(idyz,3)-c2*ra(5).*beta_v(idxx,1);
                        beta_v(idzz,3) =  beta_v(idzz,3)-c2*ra(6).*beta_v(idxx,1);
                        
                        % Integrate monimals
                        C(n+1,m+1) = C(n+1,m+1) + detJ/fn3...
                            *sum(fXYZ_v(1:N_ijk(2),1).*...
                            fXYZ_v(1:N_ijk(2),2).*...
                            fXYZ_v(1:N_ijk(2),3).*alpha_v(1:N_ijk(2),3));
                        
                        S(n+1,m+1) = S(n+1,m+1) + detJ/fn3...
                            *sum(fXYZ_v(1:N_ijk(2),1).*...
                            fXYZ_v(1:N_ijk(2),2).*...
                            fXYZ_v(1:N_ijk(2),3).*beta_v(1:N_ijk(2),3));
                        
                        % increment var for next iter.
                        N_ijk(1) = N_ijk(2); N_ijk(2) = N_ijk(2);
                        
                        ijk(1:N_ijk(2),1) = ijk(1:N_ijk(2),1)+1;
                        ijk(N_ijk(2)+1:N_ijk(2)+n+2,2) =  ijk(N_ijk(2)+1:N_ijk(2)+n+2,2) + fliplr(indices(1:n+2)')'-1;
                        ijk(N_ijk(2)+1:N_ijk(2)+n+2,3) =  ijk(N_ijk(2)+1:N_ijk(2)+n+2,3) + indices(1:n+2)-1;
                        
                        fn3 = fn3*(n+4);
                        fXYZ_v(1:N_ijk(1),1) = fXYZ_v(1:N_ijk(1),1).*ijk(1:N_ijk(1),1);
                        fXYZ_v(N_ijk(2)+1:N_ijk(2)+n,2) =  ijk(N_ijk(2)+1:N_ijk(2)+n,2).*fXYZ_v(N_ijk(2)-n:N_ijk(2)-1,2);
                        fXYZ_v(N_ijk(2)+3:N_ijk(2)+n+2,3) =  ijk(N_ijk(2)+3:N_ijk(2)+n+2,3).*fXYZ_v(N_ijk(2)-n+1:N_ijk(2),3);
                        
                        alpha_v(:,1:2) = alpha_v(:,2:3);
                        alpha_v(:,3) = 0;
                        beta_v(:,1:2) = beta_v(:,2:3);
                        beta_v(:,3) = 0;
                    end
                end
            end
            
            C = C/DETJ*6; % divide out that volume
            S = S/DETJ*6; % divide out that volume
            
            % load result into object
            obj.C=C; % cosine ocefficients
            obj.S=S; % sine coefficients
            obj.Ro=a; % radius of sphere
            obj.Mu=Mu; % gravitational parameter
        end
        
        function obj = ConstructFromMascon(obj,pts,mu,Ro,N)
            
            obj.Ro = Ro;
            obj.Mu = sum(mu); 
            obj.C = zeros(N+1,N+1);
            obj.S = zeros(N+1,N+1);
            obj.C(1,1) = 1.0; 
            
            xa = pts(:,1)/Ro;
            ya = pts(:,2)/Ro;
            za = pts(:,3)/Ro;
            ra = calc_Mag(pts)/Ro;
            ra2 = ra.^2;
            mu = mu/obj.Mu;
            
            C_d=zeros(length(pts(:,1)),2);
            S_d=zeros(length(pts(:,1)),1);
            C_v=zeros(length(pts(:,1)),3);
            S_v=zeros(length(pts(:,1)),3);
            C_v(:,1) = 1.0; C_d(:,1) = 1.0;

            
            for n=1:N
                c1 = (2*n-1)*sqrt((2*n-1)/((2*n+1)*n^2));
                c2 = sqrt(((2*n-3)*(n-1)^2)/((2*n+1)*n^2));
                C_v(:,3) = c1*za.*C_v(:,1)-c2*ra2.*C_v(:,2);
                C_v(:,2) = C_v(:,1);
                C_v(:,1) = C_v(:,3);
                obj.C(n+1,1) = (mu'*C_v(:,1));
            end
            
            cd = 1/sqrt(3); % first coeff in diagonal recursion
                
            for m=1:N
                
                C_d(:,2) = cd*(xa.*C_d(:,1)-ya.*S_d(:));
                S_d(:) = cd*(ya.*C_d(:,1)+xa.*S_d(:));
                C_d(:,1)=C_d(:,2);
                
                obj.C(m+1,m+1) = (mu'*C_d(:,1));
                obj.S(m+1,m+1) = (mu'*S_d(:));
                    
                C_v(:,2:3)=0;
                C_v(:,1) = C_d(:,1);
                S_v(:,2:3)=0;
                S_v(:,1) = S_d(:);
                
                cd  = (2*m+1)/sqrt((2*m+2)*(2*m+3));
                for n = m+1:N
                    c1 = (2*n-1)*sqrt((2*n-1)/((2*n+1)*(n+m)*(n-m)));
                    c2 = sqrt(((2*n-3)*(n+m-1)*(n-m-1))/((2*n+1)*(n+m)*(n-m)));
                     
                    C_v(:,3) = c1*za.*C_v(:,1)-c2*ra2.*C_v(:,2);
                    C_v(:,2) = C_v(:,1);
                    C_v(:,1) = C_v(:,3);
                    
                    S_v(:,3) = c1*za.*S_v(:,1)-c2*ra2.*S_v(:,2);
                    S_v(:,2) = S_v(:,1);
                    S_v(:,1) = S_v(:,3);
                    
                    obj.C(n+1,m+1) = (mu'*C_v(:,1));
                    obj.S(n+1,m+1) = (mu'*S_v(:,1));
                end
            end
            
%             for m = 1:N
%                 for n=m+1:N
%                     c1 = (2*n-1)*sqrt((2*n-1)/((2*n+1)*(n+m)*(n-m)));
%                     c2 = sqrt(((2*n-3)*(n+m-1)*(n-m-1))/((2*n+1)*(n+m)*(n-m)));
%                         
%                         
%                 end
%             end
        end
    end
end

