function failed = Test_SphericalHarmonicModel_Bennu(failed)
disp('====================================================================')
disp('Testing Spherical Harmonic Gravity Model')
disp('====================================================================')
disp('Comparison to algorithms of Eckman R.A., Brown, A.J., Adamo, D.R.,')
disp('"Normalization and Implemenation of Three Gravitational Acceleration')
disp('Models," NASA TP-2016-218604, 2014.')
disp('--------------------------------------------------------------------')
disp(' ')

format long


% Test Cnm Snm calculation from polyhedron 
%--------------------------------------------------------------------------
tol = 1e-14;                        % tolerance for error calculation
Nsample = 1000;                     % number of times to calculate to timing
P = [300,-100,100];                 % location of calculation
N = 8;                              % number of harmonics
Mu = 1.0;                           % standard grav param
meshfile = 'Bennu_1442.obj';        % mesh file

M = SurfaceMesh(meshfile);          % construct surface mesh

SH = SphericalHarmonicModel;        % construct empty SH object
SH = SH.initializeFromMesh(M,Mu,N); % calc our harmonic coeffs
SH.C(1,1)=1;                        % assumed by test functions
SH.C(2,1)=0; SH.C(2,2)=0;           % assumed by test functions
SH.S(2,1)=0; SH.S(2,2)=0;           % assumed by test functions


disp('Body            : Bennu')
disp(['Mesh File       : ',meshfile])
disp(['Degree and Order: ',num2str(N)])
disp(['Sample Location : ',num2str(P),' meters'])
disp(' ')
disp('--------------------------------------------------------------------')
disp('Our Implementation -- normalized Gottlieb') 

tic
for i = 1:Nsample
a = SH.acceleration(P);
end
toc
disp(['Acceleration = ',num2str(a,'%.8e')])


disp('--------------------------------------------------------------------')
disp('Eckman et al. -- normalized Gottlieb') 
tic
for i = 1:Nsample
accel_gottlieb = gottliebnorm(SH.Mu, SH.Ro, P', SH.C, SH.S, N, N, eye(3));
end
toc
disp(['Acceleration = ',num2str(accel_gottlieb','%.8e')])

disp('--------------------------------------------------------------------')
disp('Eckman et al. -- normalized Pines') 
tic
for i = 1:Nsample
accel_pines = pinesnorm(SH.Mu, SH.Ro, P', SH.C, SH.S, N, N, eye(3));
end
toc
disp(['Acceleration = ',num2str(accel_pines','%.8e')])

disp('--------------------------------------------------------------------')

disp('Eckman et al. -- normalized Lear') 
tic
for i =1:Nsample
accel_lear = learnorm(SH.Mu, SH.Ro, P', SH.C, SH.S, N, N, eye(3));
end
toc
disp(['Acceleration = ',num2str(accel_lear','%.8e')])
disp(' ')



% check our error 
maxDiff = norm(accel_gottlieb'-a)/norm(a);
maxDiff = max(maxDiff, norm(accel_pines'-a)/norm(a));
maxDiff = max(maxDiff, norm(accel_lear'-a)/norm(a));

if maxDiff < tol
    disp('    PASSED: SH acceleration comparison Eckman et. al.')
else
    disp(' ')
    disp('    FAILED: SH acceleration comparison Eckman et. al.')
    disp(' ')
    failed=true;
end

disp(' ')

end
