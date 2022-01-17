classdef SphericalHarmonicModel
%==========================================================================
% Spheical Harmonic graviational model. Potential, acceleration,
% and gravitational gradient calculations are based on Gottlieb 1993,
% with the standard recursion algorithm of Lundberg and Schultz 1988.
% This implementation has been validated against the normalized Gottlieb, 
% Lear, and Pines matlab implementations of Eckman and Schultz 2014 for 
% expansions up to degree and order 180. 
%
% Features methods to calculate harmonic coefficients from polyhedral
% surface meshes and mascon distributions according to Werner 1997. Here, 
% Werner's method is prefered to that of Jamet 2004, Tsoulis 2009 due to 
% the singularities present in the recursion of the later.
% -------------------------------------------------------------------------
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
% Jamet O, Thomas E, "A linear algorithm for computing the spherical 
% coefficients of the gravitational potential from a constant density 
% polyhedron," In: Proceedings of the second international GOCE user 
% workshop, GOCE, The Geoid and Oceanography, ESA-ESRIN, Frascati, Italy, 
% Citeseer, pp 8–10, 2004.
%
% Tsoulis D, Jamet O, Verdun J, Gonindard N (2009) "Recursive algorithms 
% for the computation of the potential harmonic coefficients of a constant 
% density polyhedron," J Geod 83:925–942, 2009.       
%==========================================================================
    
    properties (GetAccess=public)
        frame; % integrator needs to know if bff or inertial x,y,z
       
        L; % degree of model
        C; % cosine coefficients
        S; % sine coefficients
        Ro; % brillouin sphere radius
        Mu; % gravitational parameter
    end
    
    methods (Access=public)
        function obj = SphericalHarmonicModel(mesh,Mu,N)
        % provides several options for initializing SH models
        %------------------------------------------------------------------
        % Implementations:
        %
        %   SphericalHarmonicModel(SurfaceMesh,Mu,N) or...
        %   SphericalHarmonicModel('filename.obj',Mu,N)
        %
        %       nominal construction using Werner 1997 algorithm to
        %       construct degree N harmonic model from polyhedron
        %
        %   SphericalHarmonicModel(SurfaceMesh,Mu)
        %
        %       same as above defaults to N=4
        %
        %   SphericalHarmonicModel(MasconModel,Mu,N)
        %
        %       nominal construction using Werner 1997 algorithm to
        %       construct degree N harmonic model from mascons
        %
        %   SphericalHarmonicModel('filename.gfc')
        %
        %       reads coefficients, Mu, Ro from gfc file
        %
        %   SphericalHarmonicModel
        %
        %       blank construction
        %
        %------------------------------------------------------------------
        % Inputs:
        %   mesh - SurfaceMesh, MasconModel, 'filename.obj','filename.gfc'
        %   Mu --- 
        %   N ----
        %------------------------------------------------------------------
            obj.frame = 'BFF';
            
            if nargin==2
                N = 4;
            elseif nargin==1 % need to allow construction from GCF file
                assert(isa(mesh,'MasconModel') || contains(mesh,'.gfc'))
                N=4;
            elseif nargin~=3 && nargin~=0
                error('incorrect number of inputs')
            end
            
            if nargin~=0
                
                if (isstring(mesh) || ischar(mesh))
                    
                    if contains(mesh,'obj')
                        mesh = SurfaceMesh(mesh); 
                    elseif constains(mesh,'.gfc')
                        obj = readGFC(mesh);
                    end
                    
                elseif ~isa(mesh,'SurfaceMesh') && ~isa(mesh,'MasconModel')
                    error('mesh needs to be SurfaceMesh, MasconModel, or name of obj file')
                end
                
                if isa(mesh,'SurfaceMesh')
                    obj = obj.initializeFromMesh(mesh,Mu,N);
                elseif isa(mesh,'MasconModel')
                    obj = obj.initializeFromMascons(mesh,N);
                end
            end
            
        end
        function obj = initializeFromMesh(obj,mesh,Mu,N)
        % calculate harmonic coefficients from polyhedral surface def
        %------------------------------------------------------------------
        % Werner, R.A., Spherical Harmonic Coefficients  for the Potential
        % of a Constant-Density Polyhedron, Computers & Geoscience,
        % Vol. 23, No. 10, pp. 1071-1077, 1997.  
        %------------------------------------------------------------------
        
            % Initializes harmonic coefficients
            %--------------------------------------------------------------
            if isstring(mesh) || ischar(mesh)
                mesh = SurfaceMesh(mesh);
            elseif ~isa(mesh,'SurfaceMesh')
                error('  incorrect mesh input: currently support .obj files or SurfaceMesh objects')
            end
            
            % ensure we're working with polyhedron
            mesh.flatten();
            
            pts = mesh.coordinates;
            f = mesh.faces;
            
            a = max(vecnorm(pts,2,2)); % radius of expansion sphere
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
        function obj = initializeFromMascons(obj,masconModel,N)
        % calculate harmonic coefficients from mascon distribution
        %------------------------------------------------------------------
        % Werner, R.A., Spherical Harmonic Coefficients  for the Potential
        % of a Constant-Density Polyhedron, Computers & Geoscience,
        % Vol. 23, No. 10, pp. 1071-1077, 1997.  
        %------------------------------------------------------------------
          
            obj.Ro = max(vecnorm(masconModel.positions,2,2));
            obj.Mu = sum(masconModel.mu); 
            obj.C = zeros(N+1,N+1);
            obj.S = zeros(N+1,N+1);
            obj.C(1,1) = 1.0; 
            
            xa = masconModel.positions(:,1)/Ro;
            ya = masconModel.positions(:,2)/Ro;
            za = masconModel.positions(:,3)/Ro;
            ra = vecnorm(masconModel.positions)/Ro;
            ra2 = ra.^2;
            mu = masconModel.mu/obj.Mu;
            
            C_d=zeros(masconModel.numElements,2);
            S_d=zeros(masconModel.numElements,1);
            C_v=zeros(masconModel.numElements,3);
            S_v=zeros(masconModel.numElements,3);
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
                
                C_d(:,2) = cd*(xa.*C_d(:,1) - ya.*S_d(:));
                S_d(:)   = cd*(ya.*C_d(:,1) + xa.*S_d(:));
                C_d(:,1) = C_d(:,2);
                
                obj.C(m+1,m+1) = (mu'*C_d(:,1));
                obj.S(m+1,m+1) = (mu'*S_d(:));
                    
                C_v(:,2:3) = 0;
                C_v(:,1)   = C_d(:,1);
                S_v(:,2:3) = 0;
                S_v(:,1)   = S_d(:);
                
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
            

        end
        function obj = readGFC(obj,fileName)
        % reads a GCF file to create spherical harmonic model
            
            obj.C = [];
            obj.S = [];
            
            fid = fopen(fileName,'r');
            tline = fgetl(fid);
            
            while ischar(tline)
                data = strsplit(tline);
                
                if contains(tline,"gfc")
                    ni = round(str2double(data(2)));
                    mi = round(str2double(data(3)));
                    obj.C(ni+1,mi+1) = str2double(data(4));
                    obj.S(ni+1,mi+1) = str2double(data(5));
                    
                elseif contains(tline,"radius")
                    obj.Ro = str2double(data(2));
                    
                elseif contains(tline,"max_degree")
                    ni = str2double(data(2));
                    ob.C = zeros(ni+1,ni+1);
                    ob.S = zeros(ni+1,ni+1);
                    
                elseif contains(tline,"earth_gravity_constant")
                    obj.Mu = str2double(data(2));
                    
                end
                tline = fgetl(fid);
            end
            fclose(fid);
        end
        
        function [potential] = potential(obj,p)
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
            
            r = sqrt(p(:,1).^2+p(:,2).^2+p(:,3).^2); % radial distance
            rinv = 1./r;
            
            potential = ones(size(p,1),1);
            
            drdx=p.*rinv;         % unit vector
            gamma = obj.Ro*rinv;  % ratio of brillioun radius to distance
            A = - obj.Mu*rinv;    % Point mass attraction
            
            % assoc. Legendre Poly vertical recursion
            Pnm_v = zeros(size(p,1),3); 
            Pnm_v(:,1) = 1;
            
            % Zonal Harmonics
            %--------------------------------------------------------------
            for n = 1:N
                
                % preallocate commonly used terms
                n1=n+1;
                c1 = sqrt(((2*n-1)*(2*n+1))/((n)*(n)));
                c2 = sqrt(((2*n+1)/(2*n-3))*((n-1)*(n-1)/((n)*(n))));
                
                % recursion of dPnm
                Pnm_v(:,3) = Pnm_v(:,1);
                Pnm_v(:,1) = gamma.*(c1*drdx(:,3).*Pnm_v(:,1)-c2*gamma.*Pnm_v(:,2));
                Pnm_v(:,2) = Pnm_v(:,3);

                potential = potential + Pnm_v(:,1)*Cnm(n1,1);

            end
            
            
            Cm = drdx(:,1);
            Sm = drdx(:,2);
            Pnm_d = sqrt(3)*gamma;
            
            % Sectoral and Tessoral harmonics
            %--------------------------------------------------------------
            for m = 1:N
                
                % precalc commonly used terms
                m1 = m+1; % commonly used index
            
                % n=m calc derivatives
                potential = potential + Pnm_d.*(Cnm(m1,m1)*Cm+Snm(m1,m1)*Sm);
                
                % initialize vertical recursion
                Pnm_v(:,1) = Pnm_d;
                Pnm_v(:,2) = 0.0;
                Pnm_v(:,3) = 0.0;

                for n = m+1:N
                    
                    % preallocate commonly used terms
                    n1=n+1;

                    % recursion of dPnm
                    Pnm_v(:,3) = Pnm_v(:,1);
                    Pnm_v(:,1) = gamma.*(sqrt(((2*n-1)*(2*n+1))/((n+m)*(n-m)))*drdx(:,3).*Pnm_v(:,1)-gamma.*sqrt(((2*n+1)/(2*n-3))*((n-m-1)*(n+m-1)/((n+m)*(n-m)))).*Pnm_v(:,2));
                    Pnm_v(:,2) = Pnm_v(:,3);
                    
                    % n < m calc derivatives
                    potential = potential + Pnm_v(:,1).*(Cnm(n1,m1)*Cm+Snm(n1,m1)*Sm);
                   
                end

                
                Cm_old = Cm;
                Cm = drdx(:,1).*Cm-drdx(:,2).*Sm;
                Sm = drdx(:,1).*Sm+drdx(:,2).*Cm_old;
                Pnm_d = gamma*sqrt(1+0.5/m1).*Pnm_d;
            end
            potential = potential.*A;
        end 
        function [acceleration] = acceleration(obj,p)
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
                        
            r = sqrt(p(:,1).^2+p(:,2).^2+p(:,3).^2); % radial distance
            rinv = 1./r;
            drdx=p.*rinv;
            
            dUdr = 1;
            dUdeps = 0;
            acceleration=zeros(size(p,1),3);
            gamma = obj.Ro*rinv; 
            A = obj.Mu*rinv.^2; % Point mass attraction
            
            % assoc. Legendre Poly vertical recursion
            Pnm_v = zeros(size(p,1),3); 
            Pnm_v(:,1) = 1;
            
            % Zonal Harmonics
            %--------------------------------------------------------------
            for n = 1:N
                
                % preallocate commonly used terms
                n1=n+1;
                c1 = sqrt(((2*n-1)*(2*n+1))/((n)*(n)));
                c2 = sqrt(((2*n+1)/(2*n-3))*((n-1)*(n-1)/((n)*(n))));
                
                % recursion of dPnm
                Pnm_v(:,3) = Pnm_v(:,1);
                Pnm_v(:,1) = gamma.*(c1*drdx(:,3).*Pnm_v(:,1)-c2*gamma.*Pnm_v(:,2));
                Pnm_v(:,2) = Pnm_v(:,3);

                dUdr = dUdr + n1*Pnm_v(:,1)*Cnm(n1,1);

            end
            
            D=2.0;
            Cm_old = ones(size(p,1),1);
            Sm_old = zeros(size(p,1),1);
            Cm = drdx(:,1);
            Sm = drdx(:,2);
            Pnm_d = sqrt(3)*gamma;
            
            % Sectoral and Tessoral harmonics
            %--------------------------------------------------------------
            for m = 1:N
                
                % precalc commonly used terms
                m1 = m+1; % commonly used index
            
                % n=m calc derivatives
                dUdr = dUdr + (m+m1)*Pnm_d.*(Cnm(m1,m1)*Cm+Snm(m1,m1)*Sm);
                acceleration(:,1) = acceleration(:,1) + Pnm_d.*m.*(Cnm(m1,m1)*Cm_old+Snm(m1,m1)*Sm_old);
                acceleration(:,2) = acceleration(:,2) + Pnm_d.*m.*(Snm(m1,m1)*Cm_old-Cnm(m1,m1)*Sm_old);
                dUdeps = dUdeps + sqrt((2*m)/D)*Pnm_d.*(Cnm(m1,m)*Cm_old+Snm(m1,m)*Sm_old);
                
                % initialize vertical recursion
                Pnm_v(:,1) = Pnm_d;
                Pnm_v(:,2) = 0.0;
                Pnm_v(:,3) = 0.0;

                for n = m+1:N
                    
                    % preallocate commonly used terms
                    n1=n+1;

                    % recursion of dPnm
                    Pnm_v(:,3) = Pnm_v(:,1);
                    Pnm_v(:,1) = gamma.*(sqrt(((2*n-1)*(2*n+1))/((n+m)*(n-m)))*drdx(:,3).*Pnm_v(:,1)-gamma.*sqrt(((2*n+1)/(2*n-3))*((n-m-1)*(n+m-1)/((n+m)*(n-m)))).*Pnm_v(:,2));
                    Pnm_v(:,2) = Pnm_v(:,3);
                    
                    % n < m calc derivatives
                    dUdr = dUdr + (m+n+1)*Pnm_v(:,1).*(Cnm(n1,m1)*Cm+Snm(n1,m1)*Sm);
                    dUdeps = dUdeps +  sqrt((n-m+1)*(n+m)/D)*Pnm_v(:,1).*(Cnm(n1,m)*Cm_old+Snm(n1,m)*Sm_old);
                    acceleration(:,1) = acceleration(:,1) + Pnm_v(:,1).*m.*(Cnm(n1,m1)*Cm_old+Snm(n1,m1)*Sm_old);
                    acceleration(:,2) = acceleration(:,2) + Pnm_v(:,1).*m.*(Snm(n1,m1)*Cm_old-Cnm(n1,m1)*Sm_old);
                end

                D=1.0; % to handle horizontal ratio of Nnm w/ kronecker
                Cm_old = Cm;
                Sm_old = Sm;
                Cm = drdx(:,1).*Cm_old-drdx(:,2).*Sm_old;
                Sm = drdx(:,1).*Sm_old+drdx(:,2).*Cm_old;
                Pnm_d = gamma.*sqrt(1+0.5/m1).*Pnm_d;
            end
  
            dUdr=dUdr+drdx(:,3).*dUdeps;
            acceleration(:,1) = A.*(acceleration(:,1)-dUdr.*drdx(:,1));
            acceleration(:,2) = A.*(acceleration(:,2)-dUdr.*drdx(:,2));
            acceleration(:,3) = A.*(acceleration(:,3)-dUdr.*drdx(:,3)+dUdeps);   
        end       
        function [laplacian] = Laplacian(obj,p)
            laplacian=zeros(size(p,1));
            r = sqrt(p(:,1).^2+p(:,2).^2+p(:,3).^2);
            laplacian(r<obj.Ro) = -4*pi*obj.Mu;
            
        end
        
        %function [VVU] = gravityGradient(obj,P) 
    
    end
end

