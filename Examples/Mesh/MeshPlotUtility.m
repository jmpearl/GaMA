close all; clear all; clc;

meshFile = 'Eros_46906.obj';           % define surface mesh file
mesh = SurfaceMesh(meshFile);          % create tri surface mesh
meshCoarse = mesh.coarsen(2000);       % coursen to 2000 faces

% Example 1 : default
%--------------------------------------------------------------------------

figure(1)
meshCoarse.plot()  
  

% Example 2 : wireframe
%--------------------------------------------------------------------------

% figure(1)
% meshCoarse.plot('plotType','wireFrame',... 
%                 'edgeColor','r',...
%                  'lineStyle','--')  


% Example 3 : mesh smoothness
%--------------------------------------------------------------------------

% figure(1)
% meshCoarse.plot('plotType','alpha',... 
%                 'edgeColor','k')       % optional edge color
          
% Example 4 : mesh quality
%--------------------------------------------------------------------------

% figure(1)
% meshCoarse.plot('plotType','minIncludedAngle',... 
%                 'edgeColor','k',...               % optional edge color
%                 'faceAlpha',1.0)                  % optional transparency
            
            
% Example 5 : constant face color
% optional face color can be specified w/ RGB vector or matlab color char
%--------------------------------------------------------------------------

% figure(1)
% meshCoarse.plot('plotType','userDefined',... 
%                 'edgeColor','k',...             
%                 'faceColor',[1,1,1],...         
%                 'faceAlpha',0.9)  
         

% Example 6 : variable face color
% face color can be set proportional to a Nfx1 variable here we use the x
% position of the face centroid
%--------------------------------------------------------------------------

% centroids = meshCoarse.faceCentroids();
% colorData = centroids(:,1); 
% figure(1)
% meshCoarse.plot('plotType','userDefined',... 
%                 'edgeColor','k',...             
%                 'faceColor',colorData,...         
%                 'faceAlpha',0.9,...
%                 'lineStyle','-') 
            
 

% Example 7 : constant face color w/ light source 
%--------------------------------------------------------------------------

% figure(1)  % with ray tracing
% meshCoarse.plot('plotType','userDefined',... 
%                 'edgeColor','k',...             
%                 'faceColor','w',...         
%                 'faceAlpha',0.9,...
%                 'lineStyle','-',...
%                 'lighting',[1,0,0],... % specify vec of incoming light ray
%                 'shadowing','on')        % ray tracing on (default)
%             
% 
% figure(2) % without ray tracing
% meshCoarse.plot('plotType','userDefined',... 
%                 'edgeColor','k',...             
%                 'faceColor','w',...         
%                 'faceAlpha',0.9,...
%                 'lineStyle','-',...
%                 'lighting',[1,0,0],... % specify vec of incoming light ray
%                 'shadowing','off')       % override default ray-tracing
%           
%                 
