function [pts_q,V_q] = preprocess_Mascon_P2_H_L1G1(pts,tet,ptsT,fT)


[ ~,TetsOfFace,~ ] = calc_FaceRelations_v2_Vol(tet);

BCtets=[];
for i =1:length(TetsOfFace)
    if length(TetsOfFace{i})==1
        BCtets = [BCtets,TetsOfFace{i}];
    end
end


% Relations
[ edges,~,EdgesOfTet ] = calc_EdgeRelations_v2_Vol(tet);


% project boundary points
mpts = (pts(edges(:,1),:)+pts(edges(:,2),:))/2;
[edgesBC,ne] = calc_BCedges(tet,pts);
Res = 2*mean(calc_Mag(pts(edges(edgesBC(:),1),:)-pts(edges(edgesBC(:),2),:)));
[pts_new] = project_Points2TriMesh(ptsT,fT,mpts(edgesBC,:),ne,Res);
mpts(edgesBC,:)=pts_new;



[coord,weight] = NewtonCotes_Quadrature_3DTet(8);
[phi,phiu,phiv,phiw] = phi2_3D(coord(:,1),coord(:,2),coord(:,3));
f= zeros(size(phi,1),4);
fBC= zeros(size(phi,1),1);


MassInterior = zeros(size(pts,1),1);

MassBC  = zeros(length(BCtets),1);
pts_q_BC = zeros(length(BCtets),1);

 for k = 1:length(BCtets)
     % midpoint of edges
     
     i = BCtets(k);
     mpts_local = mpts(abs(EdgesOfTet(i,:)),:);
     X = [pts(tet(i,:),:);mpts_local];
     
     % high degree quadrature is used to calculated the weights for
     % quadrature points in a curved mesh. For each high degree quadrature
     % point the determinant of the jacobian is calculated and multiplied
     % by phi.
     
     for j = 1:size(phi,1)
         
         F = [phiu(j,:)*X(:,1),phiv(j,:)*X(:,1),phiw(j,:)*X(:,1);...
             phiu(j,:)*X(:,2),phiv(j,:)*X(:,2),phiw(j,:)*X(:,2);...
             phiu(j,:)*X(:,3),phiv(j,:)*X(:,3),phiw(j,:)*X(:,3)];
         
         fBC(j,1) = abs(det(F));
         
         com(j,1:3) = phi(j,:)*X*(det(F));
         
     end
     
     MassBC(k) = fBC'*weight';
     pts_q_BC(k,1:3) = fBC'*com/MassBC(k);
     
 end
 
 
 for i = 1:size(tet,1)
     
     if sum(BCtets==i)==0
     % midpoint of edges
     mpts_local = mpts(abs(EdgesOfTet(i,:)),:);
     X = [pts(tet(i,:),:);mpts_local];
     
     % high degree quadrature is used to calculated the weights for
     % quadrature points in a curved mesh. For each high degree quadrature
     % point the determinant of the jacobian is calculated and multiplied
     % by phi.
     
     for j = 1:size(phi,1)
         
         F = [phiu(j,:)*X(:,1),phiv(j,:)*X(:,1),phiw(j,:)*X(:,1);...
             phiu(j,:)*X(:,2),phiv(j,:)*X(:,2),phiw(j,:)*X(:,2);...
             phiu(j,:)*X(:,3),phiv(j,:)*X(:,3),phiw(j,:)*X(:,3)];
         
         f(j,1:4) = abs(det(F))*[coord(j,1),coord(j,2),coord(j,3),...
             1-coord(j,1)-coord(j,2)-coord(j,3)];
         
     end

     MassInterior(tet(i,:)) = MassInterior(tet(i,:))+f'*weight';
     end
     
 end
 
 Indices2Remove = find(MassInterior==0);
 
 pts(Indices2Remove,:) =[];
 MassInterior(Indices2Remove,:) =[];
 

 pts_q = [pts;pts_q_BC];
 V_q  = [MassInterior;MassBC];

 
end

