function [pts_q,V_q] = preprocess_Mascon(pts,tets,ModelType)
%==========================================================================
%   Jason M. Pearl
%   University of Vermont
%   College of Engineering and Mathematical Sciences (CEMS)
%
%   jmpearl@uvm.edu | jasonmpearl91@gmail.com
%==========================================================================
% Generates mascon distributions from a tetrahedral mesh discretizing the
% interior of a body. Numerical quadrature is used to distribute mascons
% and three techniques are available to the user: 1st order Gauss-Legendre,
% 1st order Newton-Cotes, and 2nd order Newton-Cotes; with flags: Gauss-1,
% Lagrange-1, and Lagrange-2 respectively.
%-------------------------------------------------------------------------%
% Inputs:
%   pts -------- vertices of trahedral volume mesh (Nvx3)
%   tets ------- connectivity of tetrahedral volume mesh (Ntx4)
%   ModelType -- type of quadrature (string)
%--------------------------------------------------------------------------
% Outputs:
%   pts_q ------ coordinates of quadrature points  (Nqx3)
%   V_q -------- Volume associated w/ each quadrature point (Nqx1)
%--------------------------------------------------------------------------
% Nomenclature : Nv is the number of vertices in the mesh, Nt is the number
%                of tetrahedra, Nq is the number of quadrature points. The
%                number of quadrature points depends upon the number of
%                cells and the degree of exactness of the quadrature.
%==========================================================================

p1 = pts(tets(:,1),:); % vertex 1
p2 = pts(tets(:,2),:); % vertex 2
p3 = pts(tets(:,3),:); % vertex 3
p4 = pts(tets(:,4),:); % vertex 4

V_t = 1/6 * dot(p2-p1,cross(p3-p1,p4-p1),2); % volume of tetrahedra


% 1st order Gauss-Legendre quadrature with a single quadrature point 
% located at the centroid of the element. (Stroud 1971) Tn:1-1
if strcmp('Gauss-1',ModelType)  
    
    pts_q = (p1+p2+p3+p4)/4; % assign quad-point coordinate output
    V_q = V_t; % assign quad-point volume output

% 1st order Newton-Cotes quadrature derived from the integration of 1st
% degree Lagrange interpolating polynomials. Quadrature points are located
% at the vertices of tetrahedra. (Stroud 1971) Tn:1-2
elseif strcmp('Lagrange-1',ModelType)
    
    % determine all tets associated with vertices
    [ tetsOfVertices ] = calc_VertexRelations_v2_Vol(tets);
    
    % volume associated w/ vertex is 1/4 the volume of all associated tets
    V_pts = zeros(size(pts,1),1);
    for i =1:size(pts,1)
        V_pts(i) = (sum(V_t(tetsOfVertices{i})))/4;
    end
    
    pts_q = pts; % assign quad-point coordinate output
    V_q = V_pts; % assign quad-point volume output
 
% 2nd order Newton-Cotes quadrature derived from the integration of 2nd
% degree Lagrange interpolating polynomials. Quadrature points are located
% at the vertices and edge midpoints of tetrahedra. (Stroud 1971) Tn:2-1
elseif strcmp('Lagrange-2',ModelType)
    
    % define relations between edges and tets and vertices and tets
    [ tetsOfVertices ] = calc_VertexRelations_v2_Vol(tets);
    [ edges,tetsOfEdges,~ ] = calc_EdgeRelations_v2_Vol(tets);
    
    % weighted volume to each vertex
    V_pts = zeros(size(pts,1),1);
    for i =1:size(pts,1)
        V_pts(i) = -(sum(V_t(tetsOfVertices{i})))/20;
    end
    
    % weigthed volume to each edge
    V_mpts = zeros(size(edges,1),1);
    for i =1:size(edges,1)
        V_mpts(i) = 4*(sum(V_t(tetsOfEdges{i})))/20;
    end
    
    mpts = (pts(edges(:,1),:)+pts(edges(:,2),:))/2; % midpoints
    
    pts_q = [pts;mpts]; % assign quad-point coordinate output
    V_q = [V_pts;V_mpts]; % assign quad-point volume output
    
end

