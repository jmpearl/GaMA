clear all; close all; clc;

addpath(genpath('/home/jmpearl/GravityModels/CLEO'))

meshFine = SurfaceMesh('Eros_14024.obj');
mesh = SurfaceMesh('Eros_5176.obj');
mesh.writeVTK('surfaeMesh.vtk')
volMesh = VolumeMesh(mesh);
volMesh.initializeFromOctree(300,2);
volMesh.smooth(4)
volMesh.setDegree(2)
volMesh.curve(meshFine);



volMesh.writeVTK('MeshNoClip.vtk');

