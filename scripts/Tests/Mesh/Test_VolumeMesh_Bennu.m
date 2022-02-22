function failed = Test_VolumeMesh_Bennu(failed)
% Tests consistency of VolumeMesh classes volume calculations

% Test SurfaceMesh class on 1442 facet mesh of Bennu
%--------------------------------------------------------------------------
tol = 1e-13;


meshfile = 'Bennu_2692.obj';

sm = SurfaceMesh(meshfile);          % construct surface mesh
vm = VolumeMesh(sm);                 % construct volume mesh
vm.initializeFromSimpleLattice(500); % initialize with a lattice 

disp('================================================================')
disp([' Testing VolumeMesh -- ','Geometry: ',meshfile])
disp('================================================================')

% degree 1 mesh
volError = (sm.volume-vm.volume)/vm.volume;
volError = max( volError, (sum(vm.cellVolumes())-vm.volume)/vm.volume);
volError = max( volError, (sum(vm.vertexVolumes())-vm.volume)/vm.volume);
volError = max( volError, (sum(vm.nodeVolumes())-vm.volume)/vm.volume);

if volError < tol
    disp('    PASSED: volume consistency -- degree 1')
else
    disp(' ')
    disp('    FAILED: volume consistency -- degree 1')
    disp(' ')
    failed=true;
end


% degree 2 mesh
vm.setDegree(2);
volError =  (sm.volume-vm.volume)/vm.volume;
volError = max( volError, (sum(vm.cellVolumes())-vm.volume)/vm.volume);
volError = max( volError, (sum(vm.vertexVolumes())-vm.volume)/vm.volume);
volError = max( volError, (sum(vm.nodeVolumes())-vm.volume)/vm.volume);

if volError < tol
    disp('    PASSED: volume consistency -- degree 2')
else
    disp(' ')
    disp('    FAILED: volume consistency -- degree 2')
    disp(' ')
    failed=true;
end

disp(' ')

end
