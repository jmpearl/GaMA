close all; clear all; clc;

meshFile = 'Eros_46906.obj';      % define surface mesh file
mesh = SurfaceMesh(meshFile);     % create tri surface mesh
meshCoarse = mesh.coarsen(2000);  % coursen to 2000 faces
 
% smooth 10 iteration w/ cotan weights
meshCoarse = meshCoarse.smooth(10,'cotan');

% use min included angle to test if any face-pairs need to be flipped
meshCoarse = meshCoarse.edgeFlipAll();

% refine features based on mesh smoothness
meshCoarse = meshCoarse.refineFeatures(mesh,... % hi-res mesh to project onto
                                       2,...    % number of iterations
                                       9,...    % average mesh angle (degrees)
                                       2);      % max refinement level

     
% test validity of half-edge data structure
% this wont get all invalid meshes but if you get a failure out of this
% there is definitely something very wrong.
meshCoarse.checkValid()

figure(1)
meshCoarse.plot('plotType','alpha')  

% add a fieldData and write a vtk file
alpha = meshCoarse.faceAverageMeshAngles();             % example face field
alphaVertex = meshCoarse.vertexValuesFromFaces(alpha);  % convert to vertex field
alpha = meshCoarse.faceValuesFromVertices(alphaVertex); % example converting back

meshCoarse = meshCoarse.addFaceField(alpha,'faceAngle');
meshCoarse = meshCoarse.addNodeField(alphaVertex,'vertexAngle');
meshCoarse.writeVTK('ExampleRectilinearMesh.vtk')      

% convert mesh to a P2 curvilinear mesh and write vtk file
meshCoarse = meshCoarse.clearFields();     % stored fields can be cleared out
meshCoarse = meshCoarse.setDegree(2);      % 2 3 or 4
meshCoarse = meshCoarse.curve(mesh);       % curve through projection
meshCoarse.writeVTK('ExampleCurvilinearDegree2Mesh.vtk')


meshCoarse.flatten() % mesh can be converted back to rectilinear