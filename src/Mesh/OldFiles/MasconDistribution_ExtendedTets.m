function [ cm,V ] = MasconDistribution_ExtendedTets( pts, f )
%=========================================================================%
% Mascon distribution method proposed by T.G.G Chanut, S. Aljbaae, and V.
% Carruba in "Mascon Gravitation Model Using A Shaped Polyhedral Source"
% Monthly Notices of the Royal Astronomical Society 450,3742-3749 (2015)
%
% triangular facets of a closed surface mesh are used to create tetrahedra
% that extend down to the center of the body. Each tetrahedron is then
% approximated as a point-mass. In the paper this is refered to as the
% Mascon-1 distribution, which doesn't involve splitting of the tetrahedra
%
% Inputs:
%   pts - mesh vertices (Mx3)
%   f   - connectivity of the mesh (Nx3)
%
% Outputs:
%   cm  - mascon locations (Nx3)
%
%=========================================================================%

disp('Generating Mascon Distribution...')
disp('    type : Extended-Tetrahedra')
% centroid is the average of the 4 vertices with the 4th vertices being at
% [x,y,z] = [0,0,0]
Body_Center = sum((pts(f(:,1),:)+pts(f(:,2),:)+pts(f(:,3),:))/3,1)/length(f(:,1));
pts = [pts(:,1)-Body_Center(1),pts(:,2)-Body_Center(2),pts(:,3)-Body_Center(3)];

cm = (pts(f(:,1),:)+pts(f(:,2),:)+pts(f(:,3),:))/4;
cm=[cm(:,1)+Body_Center(1),cm(:,2)+Body_Center(2),cm(:,3)+Body_Center(3)];

%Calculate Tet Volume
V = abs(dot(pts(f(:,1),:),cross(pts(f(:,3),:),pts(f(:,2),:),2),2)/6);
end

