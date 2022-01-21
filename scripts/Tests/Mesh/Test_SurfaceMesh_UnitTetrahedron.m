function failed = Test_SurfaceMesh_UnitTetrahedron(failed)
% Tests several features of surface mesh class using unit tetrahedron

disp('=======================================================================')
disp(' Testing SurfaceMesh -- Geometry: unit tetrahedron')
disp('=======================================================================')
disp('      z              we are going to manually enter the vertices and ')
disp('      o 3            connectivity for the simple case of a unit right ')
disp('      |              tetrahedron, then use it to test the various')
disp('      |              meshing utilities')
disp('      | 4     2')
disp('      o-------o y')
disp('     /')
disp('    /')
disp('   o 1')
disp(' x')
disp('------------------------------------------------------------------------')
disp(' pts = [1,0,0; % vertices')
disp('        0,1,0;')
disp('        0,0,1;')
disp('        0,0,0];')
disp('    ')
disp(' f   = [1,2,3; % connectivity')
disp('        1,4,2;')
disp('        1,3,4;')
disp('        2,4,3];')
disp('------------------------------------------------------------------------')

load('unitTetValidationCase.mat')


tol = 1e-14;

mesh = SurfaceMesh;
mesh.initializeFromFaceData(unitTet.vertices,unitTet.faces);


if sum(vecnorm(mesh.coordinates-unitTet.vertices,2,2))/...
        sum(vecnorm(unitTet.vertices,2,2))<tol
    disp('    PASSED: initialization')
else
    disp(' ')
    disp('    FAILED: initialization')
    disp(' ')
    failed=true;
end

if isequal(mesh.faces,unitTet.faces)
    disp('    PASSED: set.faces')
else
    disp(' ')
    disp('    FAILED: set.faces')
    disp(' ')
    failed=true;
end

%face_normals(1,1) =face_normals(1,1) +1e-16;
if sum(vecnorm(unitTet.faceNormals-mesh.faceNormals(),2,2))/ ...
        sum(vecnorm(unitTet.faceNormals,2,2))<tol
    disp('    PASSED: faceNormals()')
else
    disp(' ')
    disp('    FAILED: faceNormals()')
    disp(' ')
    failed=true;
end

if sum(vecnorm(unitTet.faceAreaVectors-mesh.faceAreaVectors(),2,2))/ ...
        sum(vecnorm(unitTet.faceAreaVectors,2,2))<tol
    disp('    PASSED: faceAreaVectors()')
else
    disp(' ')
    disp('    FAILED: faceAreaVectors()')
    disp(' ')
    failed=true;
end

if sum(vecnorm(unitTet.faceAreas-mesh.faceAreas(),2,2))/ ...
        sum(vecnorm(unitTet.faceAreas,2,2))<tol
    disp('    PASSED: faceAreas()')
else
    disp(' ')
    disp('    FAILED: faceAreas()')
    disp(' ')
    failed=true;
end

if sum(vecnorm(unitTet.faceCentroids-mesh.faceCentroids(),2,2))/ ...
        sum(vecnorm(unitTet.faceCentroids,2,2))<tol
    disp('    PASSED: faceCentroids()')
else
    disp(' ')
    disp('    FAILED: faceCentroids()')
    disp(' ')
    failed=true;
end

if (unitTet.volume-mesh.volume)/unitTet.volume<tol
    disp('    PASSED: volume()')
else
    disp(' ')
    disp('    FAILED: volume()')
    disp(' ')
    failed=true;
end

if norm(unitTet.centroid-mesh.centroid)/...
        norm(unitTet.centroid)<tol
    disp('    PASSED: centroid()')
else
    disp(' ')
    disp('    FAILED: centroid()')
    disp(' ')
    failed=true;
end

if isequal(unitTet.uniqueEdges,mesh.edges())
    disp('    PASSED: uniqueEdges()')
else
    disp(' ')
    disp('    FAILED: uniqueEdges()')
    disp(' ')
    failed=true;
end

if isequal(unitTet.faceNeighbors,mesh.faceNeighbors())
    disp('    PASSED: faceNeighbors()')
else
    disp(' ')
    disp('    FAILED: faceNeighbors()')
    disp(' ')
    failed=true;
end

if sum(unitTet.faceAverageMeshAngle-mesh.faceAverageMeshAngles)/ ...
        sum(unitTet.faceAverageMeshAngle)<tol
    disp('    PASSED: faceAverageMeshAngle()')
else
    disp(' ')
    disp('    FAILED: faceAverageMeshAngle()')
    disp(' ')
    failed=true;
end

if sum(vecnorm(unitTet.vertexAreaVectors-mesh.vertexAreaVectors(),2,2))/ ...
        sum(vecnorm(unitTet.vertexAreaVectors,2,2))<tol
    disp('    PASSED: vertexAreaVectors()')
else
    disp(' ')
    disp('    FAILED: vertexAreaVectors()')
    disp(' ')
    failed=true;
end

if sum(vecnorm(unitTet.edgeAreaVectors-mesh.edgeAreaVectors(),2,2))/ ...
        sum(vecnorm(unitTet.edgeAreaVectors,2,2))<tol
    disp('    PASSED: edgeAreaVectors()')
else
    disp(' ')
    disp('    FAILED: edgeAreaVectors()')
    disp(' ')
    failed=true;
end

a_sum_f = sum(mesh.faceAreaVectors(),1);
a_sum_v = sum(mesh.vertexAreaVectors(),1);
a_sum_e = sum(mesh.edgeAreaVectors(),1);

if any(abs(a_sum_f)>tol) || any(abs(a_sum_v)>tol) || any(abs(a_sum_e)>tol)
    disp(' ')
    disp('    FAILED: area vectors not summing to zero')
    disp(' ')
    failed=true;
else
    disp('    PASSED: area vectors sum to zero')
end
disp(' ')
end