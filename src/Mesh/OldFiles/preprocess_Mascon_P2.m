function [pts_q,V_q] = preprocess_Mascon_P2(pts,tet,ptsT,fT)



% Relations
[ tetsOfVertices ] = calc_VertexRelations_v2_Vol(tet);
[ edges,~,EdgesOfTet ] = calc_EdgeRelations_v2_Vol(tet);


% project boundary points
mpts = (pts(edges(:,1),:)+pts(edges(:,2),:))/2;
[edgesBC,ne] = calc_BCedges(tet,pts);
Res = 2*mean(calc_Mag(pts(edges(edgesBC(:),1),:)-pts(edges(edgesBC(:),2),:)));
[pts_new] = project_Points2TriMesh(ptsT,fT,mpts(edgesBC,:),ne,Res);
mpts(edgesBC,:)=pts_new;



[coord,weight] = NewtonCotes_Quadrature_3DTet(8);
[phi,phiu,phiv,phiw] = phi2_3D(coord(:,1),coord(:,2),coord(:,3));
f= zeros(size(phi,1),size(phi,2));
Massp = zeros(size(pts,1),1);
Massm = zeros(size(edges,1),1);



 for i = 1:size(tet,1)
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
         
         f(j,1:10) = abs(det(F))*phi(j,:);
         
     end
     
     integral = f'*weight';
     Massp(tet(i,:)) =  Massp(tet(i,:))+integral(1:4);
     Massm(abs(EdgesOfTet(i,:))) =  Massm(abs(EdgesOfTet(i,:)))+integral(5:10);
     
 end
 
 pts_q = [pts;mpts];
 V_q  = [Massp;Massm];

end

