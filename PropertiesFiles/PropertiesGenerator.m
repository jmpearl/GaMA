addpath(genpath('/home/jmpearl/GravityModels/cleanedCode/GravModel_OOP'))

AU2m = 1.496e11

bodyProperties.name = 'A433_Eros';
bodyProperties.units = 'kg_m_s_rad';
bodyProperties.mass = 6.687e15;
bodyProperties.density = 2670;
bodyProperties.aphelion = AU2m*1.7825;
bodyProperties.perihelion = AU2m*1.1334;
bodyProperties.rotationRate = 2*pi/(5.270*3600); %1/sec
save('Eros.mat','bodyProperties')

clear bodyProperties
bodyProperties.name = 'A25143_Itokawa';
bodyProperties.units = 'kg_m_s_rad';
bodyProperties.mass = 3.5e10;
bodyProperties.density = 1900;
bodyProperties.aphelion = AU2m*1.6951;
bodyProperties.perihelion = AU2m*0.9532;
bodyProperties.rotationRate = 2*pi/(12.132*3600); %1/sec
save('Itokawa.mat','bodyProperties')

clear bodyProperties
bodyProperties.name = 'A101955_Bennu';
bodyProperties.units = 'kg_m_s_rad';
bodyProperties.mass = 7.329e10;
bodyProperties.density = 1190;
bodyProperties.aphelion = AU2m*1.3559;
bodyProperties.perihelion = AU2m*0.89689;
bodyProperties.rotationRate = 2*pi/(4.296057*3600); %1/sec
save('Bennu.mat','bodyProperties')

clear bodyProperties
bodyProperties.name = 'M1_Phobos';
bodyProperties.units = 'kg_m_s_rad';
bodyProperties.mass = 1.0659e16;
bodyProperties.density = 1876;
bodyProperties.aphelion = AU2m*1.666;
bodyProperties.perihelion = AU2m*1.382;
bodyProperties.rotationRate = 2*pi/(7.6533*3600); %1/sec
save('Phobos.mat','bodyProperties')

clear bodyProperties
bodyProperties.name = 'C67P/Churyumov-Gersimenko';
bodyProperties.units = 'kg_m_s_rad';
bodyProperties.mass = 9.982e12;
bodyProperties.density = 533;
bodyProperties.aphelion = AU2m*5.5829;
bodyProperties.perihelion = AU2m*1.2432;
bodyProperties.rotationRate = 2*pi/(12.4043*3600); %1/sec
save('67P.mat','bodyProperties')

clear bodyProperties
bodyProperties.name = 'A162173_Ryugu';
bodyProperties.units = 'kg_m_s_rad';
bodyProperties.mass = 4.5e11;
bodyProperties.density = 1190;
bodyProperties.aphelion = AU2m*1.4159;
bodyProperties.perihelion = AU2m*0.9633;
bodyProperties.rotationRate = 2*pi/(7.627*3600); %1/sec
save('Ryugu.mat','bodyProperties')