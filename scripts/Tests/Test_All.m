clear all; close all; clc
addpath(genpath('../..'))

if isfile('Test_log')
    delete('Test_log')
end

diary Test_log
diary on

failureTracker = false;


% Surface Mesh Class
failureTracker = Test_SurfaceMesh_UnitTetrahedron(failureTracker);
failureTracker = Test_SurfaceMesh_Bennu(failureTracker);

% Volume Mesh Class
failureTracker = Test_VolumeMesh_Bennu(failureTracker);

% Gravity Models -- Spherical Harmonic
failureTracker = Test_SphericalHarmonicModel_Bennu(failureTracker);
failureTracker = Test_SphericalHarmonicModel_Earth(failureTracker);
disp(" ")

if ~failureTracker
    disp("ALL TESTS PASSED")
else
    disp("FAILURE DETECTED")
end

diary off