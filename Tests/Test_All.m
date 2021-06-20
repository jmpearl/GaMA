clear all; close all; clc


if isfile('Test_log')
    delete('Test_log')
end

diary Test_log
diary on

failureTracker = false;

failureTracker = Test_SurfaceMesh_UnitTetrahedron(failureTracker);
failureTracker = Test_SurfaceMesh_Bennu(failureTracker);

diary off