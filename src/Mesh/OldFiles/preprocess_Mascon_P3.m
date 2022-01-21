function [pts_q,V_q] = preprocess_Mascon_P3(pts,tet,ptsT,fT)


% Relations
[ edges,~,EdgesOfTet ] = calc_EdgeRelations_v2_Vol(tet);
[ faces,~,FacesOfTet ] = calc_FaceRelations_v2_Vol(tet);


% project boundary points
mpts1 = (2*pts(edges(:,1),:)+pts(edges(:,2),:))/3;
mpts2 = (pts(edges(:,1),:)+2*pts(edges(:,2),:))/3;
fpts  = (pts(faces(:,1),:)+pts(faces(:,2),:)+pts(faces(:,3),:))/3;

[edgesBC,ne] = calc_BCedges(tet,pts);
[facesBC,nf] = calc_BCfaces(tet,pts);

Res = 2*mean(calc_Mag(pts(edges(edgesBC(:),1),:)-pts(edges(edgesBC(:),2),:)));
[mpts_new1] = project_Points2TriMesh(ptsT,fT,mpts1(edgesBC,:),ne,Res);
mpts1(edgesBC,:)=mpts_new1;

[mpts_new2] = project_Points2TriMesh(ptsT,fT,mpts2(edgesBC,:),ne,Res);
mpts2(edgesBC,:)=mpts_new2;

[fpts_new] = project_Points2TriMesh(ptsT,fT,fpts(facesBC,:),nf,Res);
fpts(facesBC,:)=fpts_new;

[coord,weight] = NewtonCotes_Quadrature_3DTet(8);

[phi,phiu,phiv,phiw] = phi3_3D(coord(:,1),coord(:,2),coord(:,3));


f= zeros(size(phi,1),size(phi,2));
Massp = zeros(size(pts,1),1);
Massm1 = zeros(size(edges,1),1);
Massm2 = zeros(size(edges,1),1);
Massf = zeros(size(faces,1),1);
ftp = zeros(4,3);

for i = 1:size(tet,1)
    % midpoint of edges
    for k = 1:6
        if EdgesOfTet(i,k)>0
            epts_local1(k,1:3) = mpts1(EdgesOfTet(i,k),:);
            epts_local2(k,1:3) = mpts2(EdgesOfTet(i,k),:);
        else
            epts_local1(k,1:3) = mpts2(-EdgesOfTet(i,k),:);
            epts_local2(k,1:3) = mpts1(-EdgesOfTet(i,k),:);
        end
    end
    
    % reorder facets to be consistent w/ basis function ordering
    ft = faces(FacesOfTet(i,:),:);
    ftp(ft==tet(i,1))=1;
    ftp(ft==tet(i,2))=2;
    ftp(ft==tet(i,3))=3;
    ftp(ft==tet(i,4))=4;
    fts=sum(ftp,2);
    fti(1)=find(fts==6);
    fti(2)=find(fts==7);
    fti(3)=find(fts==9);
    fti(4)=find(fts==8);
    
    X = [pts(tet(i,:),:);epts_local1;epts_local2;fpts(FacesOfTet(i,fti),:)];
    
    for j = 1:size(phi,1)
        
        F = [phiu(j,:)*X(:,1),phiv(j,:)*X(:,1),phiw(j,:)*X(:,1);...
            phiu(j,:)*X(:,2),phiv(j,:)*X(:,2),phiw(j,:)*X(:,2);...
            phiu(j,:)*X(:,3),phiv(j,:)*X(:,3),phiw(j,:)*X(:,3)];
        
        f(j,1:20) = abs(det(F))*phi(j,:);
        
    end
    %VOLUME(i) = volume*weight';
    integral = f'*weight';
    Massp(tet(i,:)) =  Massp(tet(i,:))+integral(1:4);
    
    for k = 1:6
        if EdgesOfTet(i,k)>0
            Massm1(EdgesOfTet(i,k)) = Massm1(EdgesOfTet(i,k))+integral(k+4);
            Massm2(EdgesOfTet(i,k)) = Massm2(EdgesOfTet(i,k))+integral(k+10);
        else
            Massm1(-EdgesOfTet(i,k)) = Massm1(-EdgesOfTet(i,k))+integral(k+10);
            Massm2(-EdgesOfTet(i,k)) = Massm2(-EdgesOfTet(i,k))+integral(k+4);
        end
    end
    Massf(FacesOfTet(i,fti)) = Massf(FacesOfTet(i,fti)) + integral(17:20);
    
end

pts_q = [pts;mpts1;mpts2;fpts];
V_q  = [Massp;Massm1;Massm2;Massf];
end

