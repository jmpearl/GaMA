function [ cm,V ] = MasconDistribution_ExtendedTets3( pts, f, N_splits )
%=========================================================================%
% Mascon distribution method proposed by T.G.G Chanut, S. Aljbaae, and V.
% Carruba in "Mascon Gravitation Model Using A Shaped Polyhedral Source"
% Monthly Notices of the Royal Astronomical Society 450,3742-3749 (2015)
%
% triangular facets of a closed surface mesh are used to create tetrahedra
% that extend down to the center of the body. Each tetrahedron is then
% approximated as a point-mass. In the paper this is refered to as the
% Mascon-N distribution, tets are splt into N_splits sections
%
% Inputs:
%   pts - mesh vertices (Mx3)
%   f   - connectivity of the mesh (Nx3)
%   N_splits - number of times tets are split
%
% Outputs:
%   cm  - mascon locations (Nx3)
%
%=========================================================================%

disp('Generating Mascon Distribution...')
disp('    type : Extended-Tetrahedra w/ Splits')
Body_Center = sum((pts(f(:,1),:)+pts(f(:,2),:)+pts(f(:,3),:))/3,1)/length(f(:,1));
pts = [pts(:,1)-Body_Center(1),pts(:,2)-Body_Center(2),pts(:,3)-Body_Center(3)];

p1 = pts(f(:,1),:);
p2 = pts(f(:,2),:);
p3 = pts(f(:,3),:);

Nf = length(f(:,1));

cm = zeros(N_splits*Nf,3);
V = zeros(N_splits*Nf,1);

% storage vars
p11=p1;
p22=p2;
p33=p3;

for i = 1:N_splits-1
    p4 = (N_splits-i)/N_splits*p11;
    p5 = (N_splits-i)/N_splits*p22;
    p6 = (N_splits-i)/N_splits*p33;
    
    c1 = (p1+p2+p3+p4)/4;
    c2 = (p2+p3+p4+p6)/4;
    c3 = (p2+p4+p5+p6)/4;
    
    V1 = abs(dot((p1-p4),cross((p2-p4),(p3-p4),2),2)/6);
    V2 = abs(dot((p4-p6),cross((p2-p6),(p3-p6),2),2)/6);
    V3 = abs(dot((p2-p6),cross((p4-p6),(p5-p6),2),2)/6);
    
    i1 = Nf*(i-1)+1;
    i2 = Nf*i;
    
    V(i1:i2,1) = V1+V2+V3;
    
    cm(i1:i2,1:3) = (c1.*[V1,V1,V1]+c2.*[V2,V2,V2]+c3.*[V3,V3,V3])...
        ./([V(i1:i2,1),V(i1:i2,1),V(i1:i2,1)]);
    
    p1=p4; p2=p5; p3=p6;
end

i=N_splits;
i1 = Nf*(i-1)+1;
i2 = Nf*i;

V(i1:i2,1) = abs(dot((p1),cross((p2),(p3),2),2)/6);
cm(i1:i2,1:3)  = (p1+p2+p3)/4;
cm=[cm(:,1)+Body_Center(1),cm(:,2)+Body_Center(2),cm(:,3)+Body_Center(3)];


