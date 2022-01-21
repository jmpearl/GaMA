Order = 6;
for i = 0:Order
   [u,v,weights] = NewtonCotesTetrahedron(i);
   points = [u,v];
   
   quadratureScheme.points = [u,v];
   quadratureScheme.weights = weights;
   quadratureScheme.degree = i;
   quadratureScheme.name = ['tetrahedronNewtonCotes',num2str(i,'%02.f')];
   quadratureScheme.credit = 'GAMA';
   quadratureScheme.doi = 'None';
   quadratureScheme.geometry = 'tetrahedron';
   save([quadratureScheme.name,'.mat'],'quadratureScheme')
end

