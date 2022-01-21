function [ cm,v ] = octreePacking( pts,f,N_refines,ds, C2 )
%=========================================================================%
% Function to create a octree packing mascon distribution within a
% polygonal mesh. Mascons vary in size depending upon their "octree level".
% Process begins with a coarse cubic mesh. Cubes completely inside the
% polygonal mesh are kept, cubes completely outside are discarded, cubes
% that lie on the boundary are sub-divided into 8 smaller cubes. The
% process is then repeated for the new boundary cubes until the maximum
% octree level is obtained.
%
%
% Inputs:
%   pts - mesh vertices (Mx3)
%   f   - mesh connectivity (must be closed triangular surface mesh)
%         orientation of mesh faces must be outward (Nx3)
%   Nr  - Maximum number of splits (i.e. boundary resolution) (integer)
%   ds  - resolution of coarse initial mesh (integer) (1/2 min dimension of
%         asteroid is a good place to start)
%   C2  - 1.0-2.0 institutes more gradual refinement towards boundary. 1.0
%         results in no grading 2.0 more grading.
% Outputs:
%   cm  - coordinates of mascons (Lx3)
%   v   - normalized volume (element volume divided by total volume) (Lx1)
%=========================================================================%

disp('Generating Mascon Distribution...')
disp('    type : Octree Packing')

% Begin: Create Initial Coarse Mesh
%-------------------------------------------------------------------------%
% Admittedly this portion of the code is kind-of clunky and inefficient; 
% but, in comparison to the meat of the function that performs the octree
% refinement, the initial mesh generation takes an insignificant amount of
% time and I really didn't feel like wasting my time trying to make it
% quicker.
%-------------------------------------------------------------------------%

Max=max(pts,[],1); % max in each dimension
Min=min(pts,[],1); % min in each dimension

% x-y-z coordinates
xp = 0; xm = -ds; %pos and neg directions
yp = 0; ym = -ds;
zp = 0; zm = -ds;

for i = 2:100
    if xp(end)<Max(1)
        xp = [xp,xp(i-1)+ds];
    end
    if yp(end)<Max(2)
        yp = [yp,yp(i-1)+ds];
    end
    if zp(end)<Max(3)
        zp = [zp,zp(i-1)+ds];
    end
    if xm(end)>Min(1)
        xm = [xm,xm(i-1)-ds];
    end
    if ym(end)>Min(2)
        ym = [ym,ym(i-1)-ds];
    end
    if zm(end)>Min(3)
        zm = [zm,zm(i-1)-ds];
    end
end

x=[xm,xp];
Nx = length(x);
y=[ym,yp];
Ny = length(y);
z=[zm,zp];
Nz = length(z);

% Create coordinate vectors for mascon positions
X = zeros(Nx*Ny*Nz,1);
Y = zeros(Nx*Ny*Nz,1);
Z = zeros(Nx*Ny*Nz,1);
iter = 1;

for i =1:length(x)
    for j = 1:length(y)
        for k = 1:length(z)
            X(iter) = x(i);
            Y(iter) = y(j);
            Z(iter) = z(k);
            iter = iter +1;
        end
    end
end

% End: Create Initial Coarse Mesh
%-------------------------------------------------------------------------%



% Begin: Octree Refinement and Mascon Distribution Generation
%-------------------------------------------------------------------------%

%Initialize
Octree_Level =[];
cm=[];

for j = 1:N_refines
    % Only keep mascons located within the body
    % location of vertices relative to centroid dependent upon otree level
    Offset = 1/2^j*ds*C2;
    
    if j<N_refines
        
        InOrOutVar = zeros(length(X),8);
        
        for i = 1:length(X)
            
            % Determine if the vertices of each cube are external or
            % internal to the domain. InOrOutVar returns 0 if external and
            % 4 pi (approx 12.5664) if internal. 
            [ InOrOutVar(i,1) ] = InOrOut( pts, f, [X(i)+Offset,Y(i)+Offset,Z(i)+Offset] );
            [ InOrOutVar(i,2) ] = InOrOut( pts, f, [X(i)-Offset,Y(i)+Offset,Z(i)+Offset] );
            [ InOrOutVar(i,3) ] = InOrOut( pts, f, [X(i)+Offset,Y(i)-Offset,Z(i)+Offset] );
            [ InOrOutVar(i,4) ] = InOrOut( pts, f, [X(i)-Offset,Y(i)-Offset,Z(i)+Offset] );
            [ InOrOutVar(i,5) ] = InOrOut( pts, f, [X(i)+Offset,Y(i)+Offset,Z(i)-Offset] );
            [ InOrOutVar(i,6) ] = InOrOut( pts, f, [X(i)-Offset,Y(i)+Offset,Z(i)-Offset] );
            [ InOrOutVar(i,7) ] = InOrOut( pts, f, [X(i)+Offset,Y(i)-Offset,Z(i)-Offset] );
            [ InOrOutVar(i,8) ] = InOrOut( pts, f, [X(i)-Offset,Y(i)-Offset,Z(i)-Offset] );
            
        end
        
    else
        
        InOrOutVar = zeros(length(X),1);
        
        for i = 1:length(X)
            
            % For the last refinement step determine if centroid of cube is
            % internal or external
            [ InOrOutVar(i,1) ] = InOrOut( pts, f, [X(i),Y(i),Z(i)] );
        
        end
        
    end
    
    %Sum InOrOutVar result for all vertices
    IOVar = sum(InOrOutVar,2);
    
    if j<N_refines
        %Add cubes to the octree distribution if they lie completely 
        %inside the body of interest (12.5664*8 = 100.5)
        cm = [cm;X(IOVar>100),Y(IOVar>100),Z(IOVar>100)];
        
        %Record number of cubes at current octree level added to the
        %distribution
        N_add = length(X(IOVar>100));
        
        %Cubes that the domain boundary passes through
        X = X(IOVar>.1 & IOVar<100);
        Y = Y(IOVar>.1 & IOVar<100);
        Z = Z(IOVar>.1 & IOVar<100);
        
        %Subdivide boundary cube into 8 smaller cubes
        C=1/2^(j+1);
        cm_next = [X+C*ds,Y+C*ds,Z+C*ds];
        cm_next = [cm_next;X-C*ds,Y+C*ds,Z+C*ds];
        cm_next = [cm_next;X+C*ds,Y-C*ds,Z+C*ds];
        cm_next = [cm_next;X-C*ds,Y-C*ds,Z+C*ds];
        cm_next = [cm_next;X+C*ds,Y+C*ds,Z-C*ds];
        cm_next = [cm_next;X-C*ds,Y+C*ds,Z-C*ds];
        cm_next = [cm_next;X+C*ds,Y-C*ds,Z-C*ds];
        cm_next = [cm_next;X-C*ds,Y-C*ds,Z-C*ds];
        
        %Redefine x,y,z for next iteration of loop
        X = cm_next(:,1);Y = cm_next(:,2);Z = cm_next(:,3);
        
    else
        
        %In the final refinement step cubes are added to the distribution 
        %if their centriod lies within the domain
        cm = [cm;X(IOVar>12),Y(IOVar>12),Z(IOVar>12)];
        N_add = length(X(IOVar>12));
        
    end
    
    %Octree level of cubes in distribution used to calculate mascon volume
    Octree_Level = [Octree_Level;j*ones(N_add,1)];
    
end

% Calculate volume per element with elements of layer i+1 weighted 1/8 as
% much as the i layer. Then volume is normalized so that it can be
% corrected outside the function.
v = ((1/2).^(Octree_Level-1)).^3;
v = v/sum(v);

% End: Octree Refinement and Mascon Distribution Generation
%-------------------------------------------------------------------------%
end

