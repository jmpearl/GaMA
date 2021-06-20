clear all; close all; 
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
Nsample = 1000;
P = [300,-100,100];
meshfile = 'Bennu_1442.obj';
M = SurfaceMesh(meshfile);
N = 10;
Mu = 1.0;
SH = SphericalHarmonicModel;
SH = SH.initializeFromMesh(M,Mu,N);
SH.C(1,1)=1;
SH.C(2,1)=0; SH.C(2,2)=0;
SH.S(2,1)=0; SH.S(2,2)=0;


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
disp('====================================================================')

