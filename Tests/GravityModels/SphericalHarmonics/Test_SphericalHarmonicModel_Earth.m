clear all; close all;

disp('====================================================================')
disp('Testing Spherical Harmonic Gravity Model Degree and Order 120')
disp('====================================================================')
disp('Comparison to algorithms of Eckman R.A., Brown, A.J., Adamo, D.R.,')
disp('"Normalization and Implemenation of Three Gravitational Acceleration')
disp('Models," NASA TP-2016-218604, 2014.')
disp('--------------------------------------------------------------------')
disp(' ')
format long


Nsample=1000; % for averaging algo run time
P = 1.0001*[0,0,1]/norm([0,0,1]); % sample point normalized
GFCFile = 'ULux_CHAMP2013s.gfc'; %120x120 harmonic model
SH = SphericalHarmonicModel; % instantiate harmonic model
SH = SH.readGFC(GFCFile); % load in the harmonic data
P=P*SH.Ro; % set the sample point

% Eckman's implementations assume these terms are zero
% (i.e. the expansion is centered). Here we force this
% to  allow for a fair comparison
SH.C(1,1)=1;
SH.C(2,1)=0; SH.C(2,2)=0;
SH.S(2,1)=0; SH.S(2,2)=0;
N = size(SH.C,1)-1; % degree

disp('Body            : Earth')
disp(['Degree and Order: ',num2str(N)])
disp(['Sample Location : ',num2str(P),' meters'])

disp(['GFC File        : ',GFCFile])
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


