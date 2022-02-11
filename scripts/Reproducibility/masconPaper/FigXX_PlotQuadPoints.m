%==========================================================================
% Pearl J.M., Hitt, D.L.
%
% ** Note GaMA is under development small variations from published 
% ** results may occur.
%
% plots our distributions of mascons. Uncomment the distribution you want
% to get it to plot. For the extended tetrahedral distribution, ln. 24 
% setNumFaces should be 500 instead of 2000.
%==========================================================================

clear all; clc; close all;

vertices = [0,0,0;...
            1,0,0;...
            0,1,0;...
            0,0,1];

degree2Qpoints = [0,0,0;...
            1,0,0;...
            0,1,0;...
            0,0,1;...
            0,0.5,0;...
            0.5,0,0;...
            0,0,0.5;...
            0.5,0.5,0;...
            0,0.5,0.5;...
            0.5,0.5,0];

edges = [0,0,0;...
         1,0,0;...
         0,1,0;...
         0,0,0;...
         0,0,1;...
         1,0,0;...
         0,0,1;...
         0,1,0];

MS=14;
LW=1.5;
figure(1)
hold on
plot3(edges(:,1),edges(:,2),edges(:,3),'k-','LineWidth',LW)
plot3(vertices(:,1),vertices(:,2),vertices(:,3),'ks','MarkerFaceColor','c','MarkerSize',MS)
set(gcf,'Color',[1,1,1]);
set(gca,'TickLabelInterpreter','latex')
daspect([1,1,1])
view([1,0.2,0.3])
box off
axis off


figure(2)
hold on
plot3(edges(:,1),edges(:,2),edges(:,3),'k-','LineWidth',LW)
plot3(1/4,1/4,1/4,'ks','MarkerFaceColor','c','MarkerSize',MS)
set(gcf,'Color',[1,1,1]);
set(gca,'TickLabelInterpreter','latex')
daspect([1,1,1])
view([1,0.2,0.3])
box off
axis off



figure(3)
hold on
plot3(edges(:,1),edges(:,2),edges(:,3),'k-','LineWidth',LW)
plot3(degree2Qpoints(:,1),degree2Qpoints(:,2),degree2Qpoints(:,3),'ks','MarkerFaceColor','c','MarkerSize',MS)
set(gcf,'Color',[1,1,1]);
set(gca,'TickLabelInterpreter','latex')
daspect([1,1,1])
view([1,0.2,0.3])
box off
axis off


