function [ cm ] = simplePacking( pts,f,ds )
%=========================================================================%
% Function to create a simple packing mascon distribution within a
% polygonal mesh. All mascon are the same size and equally spaced using a
% uniform grid.
%
% Inputs:
%   pts - mesh vertices (Mx3)
%   f   - mesh connectivity (must be closed triangular surface mesh)
%         orientation of mesh faces must be outward (Nx3)
%   ds  - resolution of mascon distribution (space between elements) (1)
%
% Outputs:
%   cm  - coordinates of mascons (Lx3)
%=========================================================================%

disp('Generating Mascon Distribution...')
disp('    type : Simple Packing')

Max=max(pts,[],1); % max in each dimension
Min=min(pts,[],1); % min in each dimension

Max=max(pts,[],1); % max in each dimension
Min=min(pts,[],1); % min in each dimension

Range=(Max-Min); % ranges in each direction

%Extend slightly past ranges
Min = Min - 0.05*Range;
Max = Max + 0.05*Range;

% Number of elements per dimension
Nx = floor((Max(1)-Min(1))/ds);
Ny = floor((Max(2)-Min(2))/ds);
Nz = floor((Max(3)-Min(3))/ds);

% Correct so that ds is constant
Max = Min+[Nx+1,Ny+1,Nz+1]*ds;
    
% x-y-z coordinates
x = linspace(Min(1),Max(1),Nx+2);
y = linspace(Min(2),Max(2),Ny+2);
z = linspace(Min(3),Max(3),Nz+2);

% Create coordinate vectors for mascon positions
X = zeros((Nx+2)*(Ny+2)*(Nz+2),1);
Y = zeros((Nx+2)*(Ny+2)*(Nz+2),1);
Z = zeros((Nx+2)*(Ny+2)*(Nz+2),1);
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

% Only keep mascons located within the body
InOrOutVar = zeros(length(X),1);
for i = 1:length(X)
    [ InOrOutVar(i) ] = InOrOut( pts, f, [X(i),Y(i),Z(i)] );
end

% package as Nx3 matrix
cm = [X(InOrOutVar>10),Y(InOrOutVar>10),Z(InOrOutVar>10)];

end

