function failed = Test_SurfaceMesh_Bennu(failed)
% Tests several features of the SurfaceMesh class


if isfile('Bennu_TestMesh.obj')
    delete('Bennu_TestMesh.obj')
end


% Test SurfaceMesh class on 1442 facet mesh of Bennu
%--------------------------------------------------------------------------
tol = 1e-14;


meshfile = 'Bennu_1442.obj'; %'Eros_1208.obj';

mesh = SurfaceMesh(meshfile);
mesh.writeOBJ('Bennu_TestMesh.obj',8);

meshTest = SurfaceMesh('Bennu_TestMesh.obj');
totalArea = sum(mesh.faceAreas());
totalVolume = mesh.volume();


disp('================================================================')
disp([' Testing SurfaceMesh -- ','Geometry: ',meshfile])
disp('================================================================')

if sum(vecnorm(mesh.coordinates-meshTest.coordinates,2,2))/...
        sum(vecnorm(meshTest.coordinates,2,2))<tol...
        && isequal(mesh.faces,meshTest.faces)
    disp('    PASSED: read/write obj files')
else
    disp(' ')
    disp('    FAILED: read/write obj files')
    disp(' ')
    failed=true;
end

% Area Vectors
%--------------------------------------------------------------------------
% area vectors should sum to zero within tolerance
%--------------------------------------------------------------------------
a_sum_f = sum(mesh.faceAreaVectors(),1);
a_sum_v = sum(mesh.vertexAreaVectors(),1);
a_sum_e = sum(mesh.edgeAreaVectors(),1);

if any(abs(a_sum_f)/totalArea>tol) || any(abs(a_sum_v)/totalArea>tol) || any(abs(a_sum_e)/totalArea>tol)
    disp(' ')
    disp('    FAILED: area vectors not summing to zero')
    disp(['            sum(Af)= ',num2str(a_sum_f)])
    disp(['            sum(Ae)= ',num2str(a_sum_e/3)])
    disp(['            sum(Av)= ',num2str(a_sum_v/3)])
    disp(' ')
    failed=true;
else
    disp('    PASSED: area vectors sum to zero')
    
end

% Area
%-------------------------------------------------------------------------
% test for consistent total surface area calculations the vertex and edge
% area calculations are only roughly consistent w/ the total surface area
%-------------------------------------------------------------------------
a_sum_f = sum(vecnorm(abs(mesh.faceAreaVectors()),2,2))-totalArea;
a_sum_v = sum(vecnorm(abs(mesh.vertexAreaVectors()),2,2))-totalArea;
a_sum_e = sum(vecnorm(abs(mesh.edgeAreaVectors()),2,2))-totalArea;
if abs(a_sum_f)/totalArea>tol || abs(a_sum_v)/totalArea>0.05 || abs(a_sum_e)/totalArea>0.05
    disp(' ')
    disp('    FAILED: total surface area calculations')
    disp(['            surface area error = ',num2str(100*a_sum_f/totalArea),'%'])
    disp(['            surface area error = ',num2str(100*a_sum_e/totalArea),'%'])
    disp(['            surface area error = ',num2str(100*a_sum_v/totalArea),'%'])
    disp(' ')
    failed=true;
else
    disp('    PASSED: area vectors and areas consistent')
    
end
disp(' ')


% Quadrature Rules
%--------------------------------------------------------------------------
% test quadrature rules area vectors should again sum to zero
% and volumes should be consistent with M.volume()
%--------------------------------------------------------------------------
quadRule = {'G1','L1','L2','G2','B2','L3','B3','O4'};
for i=1:length(quadRule)
    disp(['  Quadrature Rule: ',quadRule{i}])
    [ptsq,Aq] = mesh.createQuadrature(quadRule{i});
    V_i = sum(dot(ptsq,Aq,2))/3;
    A_i = sum(Aq,1);
    if any(abs(A_i)/totalArea>tol)
        disp(' ')
        disp('    FAILED: area vectors not summing to zero')
        disp(['            sum(A)= ',num2str(A_i)])
        failed=true;
    else
        disp('    PASSED: area vectors sum to zero')
    end
    
    if abs(V_i-totalVolume)/totalVolume > tol
        disp(' ')
        disp('    FAILED: inconsistent volume')
        disp(['            V = ',num2str(V_i)])
        failed=true;
    else
        disp('    PASSED: volume consistent')
    end
    disp(' ')
end
end
