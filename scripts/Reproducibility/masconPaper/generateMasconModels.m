function masconModel = generateMasconModels(mesh,meshCoarse,Mu,Ni)


    insetSurfaceMesh = meshCoarse.offsetSurfaceMesh(-meshCoarse.resolution/2, ...
                                                     meshCoarse.numVertices);

    vmCoarse = VolumeMesh(meshCoarse);
    vmCoarse.initializeFromSimpleLattice(Ni-meshCoarse.numVertices);
    vmCoarse.smooth(5)

    vmFine = VolumeMesh(mesh);
    vmFine.initializeFromSimpleLattice(Ni);
    vmFine.smooth(5);

    vmOctreeDegree2 = VolumeMesh(mesh);
    vmOctreeDegree2.initializeFromOctree(round(Ni/4.0),2,2);
    vmOctreeDegree2.smooth(5);
    vmOctreeDegree2.setDegree(2);
    vmOctreeDegree2.curve(mesh);

    vmIter = VolumeMesh(meshCoarse);
    vmIter.initializeFromSurfaceIteration(1/2);
    vmIter.smooth(2)
    vmIter.setDegree(2);
    vmIter.curve(mesh);

    vmIter2 = VolumeMesh(meshCoarse);
    vmIter2.initializeFromSurfaceIteration(1/3);
    vmIter2.smooth(2)
    vmIter2.setDegree(2);
    vmIter2.curve(mesh);

    vmIter4 = VolumeMesh(mesh);
    vmIter4.initializeFromSurfaceIteration(1/2,meshCoarse.numVertices);
    vmIter4.smooth(2)

    vmIter5 = VolumeMesh(mesh);
    vmIter5.initializeFromSurfaceIteration(1/3,round(meshCoarse.numVertices*1.5));
    vmIter5.smooth(2)

    vmIter6 = VolumeMesh(mesh);
    vmIter6.initializeFromSurfaceIteration(1/4,round(meshCoarse.numVertices*2.0));
    vmIter6.smooth(2)

    % Gravity Models
    %-------------------------------------------------------------------------
    % polyhedral test mesh

    % mascon - packing
    masconModel{1} = MasconModel(insetSurfaceMesh,Mu,Ni);

    % mascon - volume mesh P1 vertex quads
    masconModel{2} = MasconModel();
    masconModel{2}.initializeFromVolumeMesh(vmCoarse,Mu,'vertex');

    % mascon - volume mesh P2 mesh vertex quad
    vmCoarse.setDegree(2);
    vmCoarse.curve(mesh);
    masconModel{3} = MasconModel();
    masconModel{3}.initializeFromVolumeMesh(vmCoarse,Mu,'vertex');
    masconModel{4} = MasconModel();
    masconModel{4}.initializeFromVolumeMesh(vmCoarse,Mu,'cell');
    masconModel{5} = MasconModel();
    masconModel{5}.initializeFromVolumeMesh(vmCoarse,Mu,'node');
    masconModel{6} = MasconModel();
    masconModel{6}.initializeFromVolumeMesh(vmCoarse,Mu,'excludesurface');

    % mascon - degree 2 octree vertex quad
    masconModel{7} = MasconModel();
    masconModel{7}.initializeFromVolumeMesh(vmOctreeDegree2,Mu,'excludesurface');

    % based on true mesh excluding surface points simple lattice
    masconModel{8} = MasconModel();
    masconModel{8}.initializeFromVolumeMesh(vmFine,Mu,'excludesurface');

    % based on true mesh excluding surface points w/ surface iteration
    masconModel{9} = MasconModel();
    masconModel{9}.initializeFromVolumeMesh(vmIter,Mu,'excludesurface');

    masconModel{10} = MasconModel();
    masconModel{10}.initializeFromVolumeMesh(vmIter2,Mu,'excludesurface');

    masconModel{11} = MasconModel();
    masconModel{11}.initializeFromVolumeMesh(vmIter4,Mu,'excludesurface');

    masconModel{12} = MasconModel();
    masconModel{12}.initializeFromVolumeMesh(vmIter5,Mu,'excludesurface');

    masconModel{13} = MasconModel();
    masconModel{13}.initializeFromVolumeMesh(vmIter6,Mu,'excludesurface');

    % based on true mesh excluding surface points
    masconModel{14} = MasconModel();
    masconModel{14}.initializeExtendedTetrahedra(vmCoarse,3,Mu);
    masconModel{15} = MasconModel();
    masconModel{15}.initializeExtendedTetrahedra(vmCoarse,3,Mu);
end