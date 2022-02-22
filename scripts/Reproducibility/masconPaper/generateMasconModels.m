function masconModel = generateMasconModels(mesh,meshCoarse,Mu,Ni)


    insetSurfaceMesh = meshCoarse.offsetSurfaceMesh(-meshCoarse.resolution/2, ...
                                                     meshCoarse.numVertices);
    meshCoarseP2 = SurfaceMesh(meshCoarse);
    meshCoarseP2.setDegree(2);
    meshCoarseP2.curve(mesh);

    % coarse mesh
    vmCoarse = VolumeMesh(meshCoarse);
    vmCoarse.initializeFromSimpleLattice(Ni-meshCoarse.numVertices);
    vmCoarse.smooth(5)

    % coarse mesh P2 
    vmCoarseP2 = VolumeMesh(vmCoarse);
    vmCoarseP2.setDegree(2);
    vmCoarseP2.curve(mesh);

    % based on full surface resolution
    vmFine = VolumeMesh(mesh);
    vmFine.initializeFromSimpleLattice(Ni);
    vmFine.smooth(5);

    % octree infill with coarse surface - degree 2
    vmOctreeDegree2 = VolumeMesh(mesh);
    vmOctreeDegree2.initializeFromOctree(round(Ni/4.0),2,2);
    vmOctreeDegree2.smooth(5);
    vmOctreeDegree2.setDegree(2);
    vmOctreeDegree2.curve(mesh);

    % surface interation 1/2 coarsening degree 2
    vmIter = VolumeMesh(meshCoarse);
    vmIter.initializeFromSurfaceIteration(1/2);
    vmIter.smooth(2)
    vmIter.setDegree(2);
    vmIter.curve(mesh);

    % surface interation 1/3 coarsening degree 2
    vmIter2 = VolumeMesh(meshCoarse);
    vmIter2.initializeFromSurfaceIteration(1/3);
    vmIter2.smooth(2)
    vmIter2.setDegree(2);
    vmIter2.curve(mesh);
    
    % surface interation 1/2 full refine surface
    vmIter4 = VolumeMesh(mesh);
    vmIter4.initializeFromSurfaceIteration(1/2,2.0*meshCoarse.numVertices);
    vmIter4.smooth(2)

    % surface interation 1/3 full refine surface
    vmIter5 = VolumeMesh(mesh);
    vmIter5.initializeFromSurfaceIteration(1/3,round(meshCoarse.numVertices*3.0));
    vmIter5.smooth(2)

     % surface interation 1/4 full refine surface
    vmIter6 = VolumeMesh(mesh);
    vmIter6.initializeFromSurfaceIteration(1/4,round(meshCoarse.numVertices*4.0));
    vmIter6.smooth(2)

    % Gravity Models
    %-------------------------------------------------------------------------

    i=1;

    % mascon - lattice packing distributions
    masconModel{i} = MasconModel();
    masconModel{i}.initializeSimplePacking(insetSurfaceMesh,Ni,Mu);
    i=i+1;
    masconModel{i} = MasconModel();
    masconModel{i}.initializeSimplePacking(insetSurfaceMesh,Ni,Mu,'bcc');
    i=i+1;
    masconModel{i} = MasconModel();
    masconModel{i}.initializeSimplePacking(insetSurfaceMesh,Ni,Mu,'fcc');
    i=i+1;

    % mascon - volume mesh P1 vertex quads
    masconModel{i} = MasconModel();
    masconModel{i}.initializeFromVolumeMesh(vmCoarse,Mu,'vertex');
    i=i+1;
    masconModel{i} = MasconModel();
    masconModel{i}.initializeFromVolumeMesh(vmCoarse,Mu,'excludesurface');
    i=i+1;

    % mascon - volume mesh P2 mesh vertex quad
    masconModel{i} = MasconModel();
    masconModel{i}.initializeFromVolumeMesh(vmCoarseP2,Mu,'vertex');
    i=i+1;
    masconModel{i} = MasconModel();
    masconModel{i}.initializeFromVolumeMesh(vmCoarseP2,Mu,'cell');
    i=i+1;
    masconModel{i} = MasconModel();
    masconModel{i}.initializeFromVolumeMesh(vmCoarseP2,Mu,'node');
    i=i+1;
    masconModel{i} = MasconModel();
    masconModel{i}.initializeFromVolumeMesh(vmCoarseP2,Mu,'excludesurface');
    i=i+1;

    % mascon - degree 2 octree vertex quad
    masconModel{i} = MasconModel();
    masconModel{i}.initializeFromVolumeMesh(vmOctreeDegree2,Mu,'excludesurface');
    i=i+1;

    % based on true mesh excluding surface points simple lattice
    masconModel{i} = MasconModel();
    masconModel{i}.initializeFromVolumeMesh(vmFine,Mu,'excludesurface');
    i=i+1;

    % based on true mesh excluding surface points w/ surface iteration
    masconModel{i} = MasconModel();
    masconModel{i}.initializeFromVolumeMesh(vmIter,Mu,'excludesurface');
    i=i+1;
    masconModel{i} = MasconModel();
    masconModel{i}.initializeFromVolumeMesh(vmIter2,Mu,'excludesurface');
    i=i+1;
    masconModel{i} = MasconModel();
    masconModel{i}.initializeFromVolumeMesh(vmIter4,Mu,'excludesurface');
    i=i+1;
    masconModel{i} = MasconModel();
    masconModel{i}.initializeFromVolumeMesh(vmIter5,Mu,'excludesurface');
    i=i+1;
    masconModel{i} = MasconModel();
    masconModel{i}.initializeFromVolumeMesh(vmIter6,Mu,'excludesurface');
    i=i+1;

    % based on extended tet method
    masconModel{i} = MasconModel();
    masconModel{i}.initializeExtendedTetrahedra(meshCoarse,Mu,1);
    i=i+1;
    masconModel{i} = MasconModel();
    masconModel{i}.initializeExtendedTetrahedra(meshCoarse,Mu,3);
    i=i+1;
    masconModel{i} = MasconModel();
    masconModel{i}.initializeExtendedTetrahedra(meshCoarse,Mu,2,'lumpcore');
    i=i+1;
    masconModel{i} = MasconModel();
    masconModel{i}.initializeExtendedTetrahedra(meshCoarseP2,Mu,1);
    i=i+1;
    masconModel{i} = MasconModel();
    masconModel{i}.initializeExtendedTetrahedra(meshCoarseP2,Mu,2,'lumpcore');
    i=i+1;
end