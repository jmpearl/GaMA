%==========================================================================
% Pearl J.M., Hitt, D.L.
%
% ** Note GaMA is under development small variations from published 
% ** results may occur if using the latest version
%==========================================================================

clear all; clc; close all;

body = 'Eros';
if strcmp(body,'Eros')
    load('Eros_46906_Mascon.mat');
    numElementsMax=22000;
    plotLocation = 'NorthEast';
elseif strcmp(body,'Bennu')
    load('Bennu_200000_Mascon.mat');
    numElementsMax=105000;
    plotLocation = 'SouthWest';
else
    error('invalid body')
end


% Packing Distributions highest altitude
%--------------------------------------------------------------------------
FS = 12;
LW = 0.5;
MS = 10;
MSsmall=3;
% 
% % % Packing Distributions lowest altitude
% % %--------------------------------------------------------------------------
% h=1
% plotName = ['masconPacking_h=',num2str(altitudes(h)),'_',body,'.png']
% figure(2)
% hold on
% plot(numElements{h}(:,1),100*L2{h}(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,2),100*L2{h}(:,2),'ks-','MarkerFaceColor','y','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,3),100*L2{h}(:,3),'ks-','MarkerFaceColor','c','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,4),100*L2{h}(:,4),'ks-','MarkerFaceColor','b','LineWidth',LW,'MarkerSize',MS)
% 
% plot(numElements{h}(:,1),100*L1{h}(:,1)/numSamplesActual(h)*100,'k--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,1),100*Linf{h}(:,1)/numSamplesActual(h)*100,'k--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% 
% plot(numElements{h}(:,2),100*L1{h}(:,2)/numSamplesActual(h)*100,'y--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,3),100*L1{h}(:,3)/numSamplesActual(h)*100,'c--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,4),100*L1{h}(:,4)/numSamplesActual(h)*100,'b--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,2),100*Linf{h}(:,2)/numSamplesActual(h)*100,'y--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,3),100*Linf{h}(:,3)/numSamplesActual(h)*100,'c--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,4),100*Linf{h}(:,4)/numSamplesActual(h)*100,'b--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% 
% plot(numElements{h}(:,1),100*L2{h}(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,2),100*L2{h}(:,2),'ks-','MarkerFaceColor','y','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,3),100*L2{h}(:,3),'ks-','MarkerFaceColor','c','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,4),100*L2{h}(:,4),'ks-','MarkerFaceColor','b','LineWidth',LW,'MarkerSize',MS)
% 
% 
% plot(numElements{h}(:,1),100*L2{h}(:,1),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% plot(numElements{h}(:,2),100*L2{h}(:,2),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% plot(numElements{h}(:,3),100*L2{h}(:,3),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% plot(numElements{h}(:,4),100*L2{h}(:,4),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% 
% legend({['Analytic Poly.'],...
%     ['PC Packing'],...
%     ['BCC Packing'],...
%     ['FCC Packing']},...
%     'interpreter','latex',...
%     'FontSize',FS,...
%     'Location','SouthWest')
% ylabel('error $\%$','interpreter','latex','FontSize',FS)
% xlabel('number of elements','interpreter','latex','FontSize',FS)
% set(gcf,'Color',[1,1,1]);
% set(gca,'FontSize',FS)
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'xscale','log');
% set(gca,'yscale','log');
% axis tight
% grid on
% %set(gca,'xlim',[10,numElementsMax])
% %axis([0.1,10000,0.1,1000])
% %legend boxoff
% box on
% 
% [f1, cdata] = myaa([4, 2]); imwrite(cdata, plotName, 'png');
% 
% close all
% 
% % % Packing Distributions highest altitude
% % %--------------------------------------------------------------------------
% h=10
% plotName = ['masconPacking_h=',num2str(altitudes(h)),'_',body,'.png']
% figure(2)
% hold on
% plot(numElements{h}(:,1),100*L2{h}(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,2),100*L2{h}(:,2),'ks-','MarkerFaceColor','y','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,3),100*L2{h}(:,3),'ks-','MarkerFaceColor','c','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,4),100*L2{h}(:,4),'ks-','MarkerFaceColor','b','LineWidth',LW,'MarkerSize',MS)
% 
% plot(numElements{h}(:,1),100*L1{h}(:,1)/numSamplesActual(h)*100,'k--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,1),100*Linf{h}(:,1)/numSamplesActual(h)*100,'k--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% 
% plot(numElements{h}(:,2),100*L1{h}(:,2)/numSamplesActual(h)*100,'y--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,3),100*L1{h}(:,3)/numSamplesActual(h)*100,'c--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,4),100*L1{h}(:,4)/numSamplesActual(h)*100,'b--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,2),100*Linf{h}(:,2)/numSamplesActual(h)*100,'y--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,3),100*Linf{h}(:,3)/numSamplesActual(h)*100,'c--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,4),100*Linf{h}(:,4)/numSamplesActual(h)*100,'b--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% 
% plot(numElements{h}(:,1),100*L2{h}(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,2),100*L2{h}(:,2),'ks-','MarkerFaceColor','y','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,3),100*L2{h}(:,3),'ks-','MarkerFaceColor','c','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,4),100*L2{h}(:,4),'ks-','MarkerFaceColor','b','LineWidth',LW,'MarkerSize',MS)
% 
% 
% plot(numElements{h}(:,1),100*L2{h}(:,1),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% plot(numElements{h}(:,2),100*L2{h}(:,2),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% plot(numElements{h}(:,3),100*L2{h}(:,3),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% plot(numElements{h}(:,4),100*L2{h}(:,4),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% 
% legend({['Analytic Poly'],...
%     ['PC Packing'],...
%     ['BCC Packing'],...
%     ['FCC Packing']},...
%     'interpreter','latex',...
%     'FontSize',FS,...
%     'Location','SouthWest')
% ylabel('error $\%$','interpreter','latex','FontSize',FS)
% xlabel('number of elements','interpreter','latex','FontSize',FS)
% set(gcf,'Color',[1,1,1]);
% set(gca,'FontSize',FS)
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'xscale','log');
% set(gca,'yscale','log');
% axis tight
% grid on
% %set(gca,'xlim',[10,numElementsMax])
% %axis([0.1,10000,0.1,1000])
% %legend boxoff
% box on
% 
% [f1, cdata] = myaa([4, 2]); imwrite(cdata, plotName, 'png');
% 
% close all
% 
% 
% % Plot extend tet distros suface
% %--------------------------------------------------------------------------
% h=1
% plotName = ['masconExtendedTet_h=',num2str(altitudes(h)),'_',body,'.png']
% figure(3)
% hold on
% plot(numElements{h}(:,1),100*L2{h}(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,18),100*L2{h}(:,18),'k>-','MarkerFaceColor','r','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,19),100*L2{h}(:,19),'k<-','MarkerFaceColor',[1,0.5,0],'LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,20),100*L2{h}(:,20),'kv-','MarkerFaceColor','y','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,21),100*L2{h}(:,21),'k^-','MarkerFaceColor','c','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,22),100*L2{h}(:,22),'ko-','MarkerFaceColor','b','LineWidth',LW,'MarkerSize',MS)
% %
% plot(numElements{h}(:,1),100*L1{h}(:,1)/numSamplesActual(h)*100,'k--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,1),100*Linf{h}(:,1)/numSamplesActual(h)*100,'k--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% %
% plot(numElements{h}(:,18),100*L1{h}(:,18)/numSamplesActual(h)*100,'r--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,19),100*L1{h}(:,19)/numSamplesActual(h)*100,'y--','Color',[1,0.5,0],'LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,20),100*L1{h}(:,20)/numSamplesActual(h)*100,'y--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,21),100*L1{h}(:,21)/numSamplesActual(h)*100,'c--','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,22),100*L1{h}(:,22)/numSamplesActual(h)*100,'b--','LineWidth',LW,'MarkerSize',MS)
% 
% plot(numElements{h}(:,18),100*Linf{h}(:,18)/numSamplesActual(h)*100,'r--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,19),100*Linf{h}(:,19)/numSamplesActual(h)*100,'y--','Color',[1,0.5,0],'LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,20),100*Linf{h}(:,20)/numSamplesActual(h)*100,'y--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,21),100*Linf{h}(:,21)/numSamplesActual(h)*100,'c--','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,22),100*Linf{h}(:,22)/numSamplesActual(h)*100,'b--','LineWidth',LW,'MarkerSize',MS)
% 
% plot(numElements{h}(:,1),100*L2{h}(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,18),100*L2{h}(:,18),'k>-','MarkerFaceColor','r','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,19),100*L2{h}(:,19),'k<-','MarkerFaceColor',[1,0.5,0],'LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,20),100*L2{h}(:,20),'kv-','MarkerFaceColor','y','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,21),100*L2{h}(:,21),'k^-','MarkerFaceColor','c','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,22),100*L2{h}(:,22),'ko-','MarkerFaceColor','b','LineWidth',LW,'MarkerSize',MS)
% 
% plot(numElements{h}(:,1),100*L2{h}(:,1),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% plot(numElements{h}(:,18),100*L2{h}(:,18),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% plot(numElements{h}(:,19),100*L2{h}(:,19),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% plot(numElements{h}(:,20),100*L2{h}(:,20),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% plot(numElements{h}(:,21),100*L2{h}(:,21),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% plot(numElements{h}(:,22),100*L2{h}(:,22),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% 
% legend({['Analytic Poly.'],...
%     ['P1-1'],...
%     ['P1-3'],...
%     ['P1-CC'], ...
%     ['P2-1'], ...
%     ['P2-CC']},...
%     'interpreter','latex',...
%     'FontSize',FS,...
%     'Location',plotLocation)
% ylabel('error $\%$','interpreter','latex','FontSize',FS)
% xlabel('number of elements','interpreter','latex','FontSize',FS)
% set(gcf,'Color',[1,1,1]);
% set(gca,'FontSize',FS)
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'xscale','log');
% set(gca,'yscale','log');
% axis tight
% grid on
% 
% box on
% 
% [f1, cdata] = myaa([4, 2]); imwrite(cdata, plotName, 'png');
% 
% close all
% 
% 
% % Plot extend tet distros r= Rmax
% %--------------------------------------------------------------------------
% h=10
% plotName = ['masconExtendedTet_h=',num2str(altitudes(h)),'_',body,'.png']
% figure(4)
% hold on
% plot(numElements{h}(:,1),100*L2{h}(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,18),100*L2{h}(:,18),'k>-','MarkerFaceColor','r','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,19),100*L2{h}(:,19),'k<-','MarkerFaceColor',[1,0.5,0],'LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,20),100*L2{h}(:,20),'kv-','MarkerFaceColor','y','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,21),100*L2{h}(:,21),'k^-','MarkerFaceColor','c','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,22),100*L2{h}(:,22),'ko-','MarkerFaceColor','b','LineWidth',LW,'MarkerSize',MS)
% %
% plot(numElements{h}(:,1),100*L1{h}(:,1)/numSamplesActual(h)*100,'k--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,1),100*Linf{h}(:,1)/numSamplesActual(h)*100,'k--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% %
% plot(numElements{h}(:,18),100*L1{h}(:,18)/numSamplesActual(h)*100,'r--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,19),100*L1{h}(:,19)/numSamplesActual(h)*100,'y--','Color',[1,0.5,0],'LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,20),100*L1{h}(:,20)/numSamplesActual(h)*100,'y--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,21),100*L1{h}(:,21)/numSamplesActual(h)*100,'c--','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,22),100*L1{h}(:,22)/numSamplesActual(h)*100,'b--','LineWidth',LW,'MarkerSize',MS)
% 
% plot(numElements{h}(:,18),100*Linf{h}(:,18)/numSamplesActual(h)*100,'r--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,19),100*Linf{h}(:,19)/numSamplesActual(h)*100,'y--','Color',[1,0.5,0],'LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,20),100*Linf{h}(:,20)/numSamplesActual(h)*100,'y--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,21),100*Linf{h}(:,21)/numSamplesActual(h)*100,'c--','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,22),100*Linf{h}(:,22)/numSamplesActual(h)*100,'b--','LineWidth',LW,'MarkerSize',MS)
% 
% plot(numElements{h}(:,1),100*L2{h}(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,18),100*L2{h}(:,18),'k>-','MarkerFaceColor','r','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,19),100*L2{h}(:,19),'k<-','MarkerFaceColor',[1,0.5,0],'LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,20),100*L2{h}(:,20),'kv-','MarkerFaceColor','y','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,21),100*L2{h}(:,21),'k^-','MarkerFaceColor','c','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,22),100*L2{h}(:,22),'ko-','MarkerFaceColor','b','LineWidth',LW,'MarkerSize',MS)
% 
% plot(numElements{h}(:,1),100*L2{h}(:,1),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% plot(numElements{h}(:,18),100*L2{h}(:,18),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% plot(numElements{h}(:,19),100*L2{h}(:,19),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% plot(numElements{h}(:,20),100*L2{h}(:,20),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% plot(numElements{h}(:,21),100*L2{h}(:,21),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% plot(numElements{h}(:,22),100*L2{h}(:,22),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% 
% legend({['Analytic Poly'],...
%     ['P1-1'],...
%     ['P1-3'],...
%     ['P1-CC'], ...
%     ['P2-1'], ...
%     ['P2-CC']},...
%     'interpreter','latex',...
%     'FontSize',FS,...
%     'Location','SouthWest')
% ylabel('error $\%$','interpreter','latex','FontSize',FS)
% xlabel('number of elements','interpreter','latex','FontSize',FS)
% set(gcf,'Color',[1,1,1]);
% set(gca,'FontSize',FS)
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'xscale','log');
% set(gca,'yscale','log');
% axis tight
% grid on
% 
% box on
% 
% [f1, cdata] = myaa([4, 2]); imwrite(cdata, plotName, 'png');
% 
% close all
% 
% stop


% Plot quad methods
%--------------------------------------------------------------------------
% h=1
% markerType = {'k<-','k<-','kv-','k^-'};
% markerType2 = {'k--','k--','k--','k--'};
% color = {'r',[1,0.5,0],'y','c'};
% plotIndices = [7,8,9,10];
% 
% plotName = ['masconQuadRule_h=',num2str(altitudes(h)),'_',body,'.png']
% 
% figure(3)
% hold on
% plot(numElements{h}(:,1),100*L2{h}(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% 
% 
% for i = 1:length(markerType)
%     i
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*L2{h}(:,j),markerType{i},'MarkerFaceColor',color{i},'LineWidth',LW,'MarkerSize',MS)
% end
% plot(numElements{h}(:,1),100*L1{h}(:,1)/numSamplesActual(h)*100,'k--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,1),100*Linf{h}(:,1)/numSamplesActual(h)*100,'k--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,1),100*L2{h}(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,1),100*L2{h}(:,1),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% 
% for i = 1:length(markerType)
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*L1{h}(:,j)/numSamplesActual(h)*100,markerType2{i},'Color',color{i},'LineWidth',LW,'MarkerSize',MS)
% end
% for i = 1:length(markerType)
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*Linf{h}(:,j)/numSamplesActual(h)*100,markerType2{i},'Color',color{i},'LineWidth',LW,'MarkerSize',MS)
% end
% for i = 1:length(markerType)
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*L2{h}(:,j),markerType{i},'MarkerFaceColor',color{i},'LineWidth',LW,'MarkerSize',MS)
% end
% for i = 1:length(markerType)
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*L2{h}(:,j),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% end
% 
% legend({['Analytic Poly.'],...
%     ['P2-vertex'],...
%     ['P2-cell'],...
%     ['P2-node (T10)'], ...
%     ['P2-exclude surface']},...
%     'interpreter','latex',...
%     'FontSize',FS,...
%     'Location',plotLocation)
% ylabel('error $\%$','interpreter','latex','FontSize',FS)
% xlabel('number of elements','interpreter','latex','FontSize',FS)
% set(gcf,'Color',[1,1,1]);
% set(gca,'FontSize',FS)
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'xscale','log');
% set(gca,'yscale','log');
% axis tight
% grid on
% 
% box on
% 
% [f1, cdata] = myaa([4, 2]); imwrite(cdata, plotName, 'png');
% 
% close all
% 
% %quad methods max altitude
% %--------------------------------------------------------------------------
% h=10
% markerType = {'k<-','k<-','kv-','k^-'};
% markerType2 = {'k--','k--','k--','k--'};
% color = {'r',[1,0.5,0],'y','c'};
% plotIndices = [7,8,9,10];
% 
% plotName = ['masconQuadRule_h=',num2str(altitudes(h)),'_',body,'.png']
% 
% figure(3)
% hold on
% plot(numElements{h}(:,1),100*L2{h}(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% 
% 
% for i = 1:length(markerType)
%     i
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*L2{h}(:,j),markerType{i},'MarkerFaceColor',color{i},'LineWidth',LW,'MarkerSize',MS)
% end
% 
% plot(numElements{h}(:,1),100*L1{h}(:,1)/numSamplesActual(h)*100,'k--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,1),100*Linf{h}(:,1)/numSamplesActual(h)*100,'k--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,1),100*L2{h}(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,1),100*L2{h}(:,1),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% 
% for i = 1:length(markerType)
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*L1{h}(:,j)/numSamplesActual(h)*100,markerType2{i},'Color',color{i},'LineWidth',LW,'MarkerSize',MS)
% end
% for i = 1:length(markerType)
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*Linf{h}(:,j)/numSamplesActual(h)*100,markerType2{i},'Color',color{i},'LineWidth',LW,'MarkerSize',MS)
% end
% for i = 1:length(markerType)
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*L2{h}(:,j),markerType{i},'MarkerFaceColor',color{i},'LineWidth',LW,'MarkerSize',MS)
% end
% for i = 1:length(markerType)
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*L2{h}(:,j),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% end
% 
% legend({['Analytic Poly.'],...
%     ['P2-vertex'],...
%     ['P2-cell'],...
%     ['P2-node (T10)'], ...
%     ['P2-exclude surface']},...
%     'interpreter','latex',...
%     'FontSize',FS,...
%     'Location',plotLocation)
% ylabel('error $\%$','interpreter','latex','FontSize',FS)
% xlabel('number of elements','interpreter','latex','FontSize',FS)
% set(gcf,'Color',[1,1,1]);
% set(gca,'FontSize',FS)
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'xscale','log');
% set(gca,'yscale','log');
% axis tight
% grid on
% 
% box on
% 
% [f1, cdata] = myaa([4, 2]); imwrite(cdata, plotName, 'png');
% 
% close all
% 
% 
% 
% %surface method min altitude
% %--------------------------------------------------------------------------
% h=1
% markerType = {'k<-','k<-','kv-'};
% markerType2 = {'k--','k--','k--'};
% color = {'r','c','m',};
% plotIndices = [6,10,12];
% 
% plotName = ['masconSurfaceMethod_h=',num2str(altitudes(h)),'_',body,'.png']
% 
% figure(3)
% hold on
% plot(numElements{h}(:,1),100*L2{h}(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% 
% 
% for i = 1:length(markerType)
%     i
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*L2{h}(:,j),markerType{i},'MarkerFaceColor',color{i},'LineWidth',LW,'MarkerSize',MS)
% end
% 
% plot(numElements{h}(:,1),100*L1{h}(:,1)/numSamplesActual(h)*100,'k--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,1),100*Linf{h}(:,1)/numSamplesActual(h)*100,'k--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,1),100*L2{h}(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,1),100*L2{h}(:,1),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% 
% for i = 1:length(markerType)
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*L1{h}(:,j)/numSamplesActual(h)*100,markerType2{i},'Color',color{i},'LineWidth',LW,'MarkerSize',MS)
% end
% for i = 1:length(markerType)
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*Linf{h}(:,j)/numSamplesActual(h)*100,markerType2{i},'Color',color{i},'LineWidth',LW,'MarkerSize',MS)
% end
% for i = 1:length(markerType)
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*L2{h}(:,j),markerType{i},'MarkerFaceColor',color{i},'LineWidth',LW,'MarkerSize',MS)
% end
% for i = 1:length(markerType)
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*L2{h}(:,j),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% end
% 
% legend({['Analytic Poly.'],...
%     ['P1'],...
%     ['P2'],...
%     ['Piecewise']},...
%     'interpreter','latex',...
%     'FontSize',FS,...
%     'Location',plotLocation)
% ylabel('error $\%$','interpreter','latex','FontSize',FS)
% xlabel('number of elements','interpreter','latex','FontSize',FS)
% set(gcf,'Color',[1,1,1]);
% set(gca,'FontSize',FS)
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'xscale','log');
% set(gca,'yscale','log');
% axis tight
% grid on
% 
% box on
% 
% [f1, cdata] = myaa([4, 2]); imwrite(cdata, plotName, 'png');
% 
% close all
% 
% %surface method max altitude
% %--------------------------------------------------------------------------
% h=10
% markerType = {'k<-','k<-','kv-'};
% markerType2 = {'k--','k--','k--'};
% color = {'r','c','m',};
% plotIndices = [6,10,12];
% 
% plotName = ['masconSurfaceMethod_h=',num2str(altitudes(h)),'_',body,'.png']
% 
% figure(3)
% hold on
% plot(numElements{h}(:,1),100*L2{h}(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% 
% 
% for i = 1:length(markerType)
%     i
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*L2{h}(:,j),markerType{i},'MarkerFaceColor',color{i},'LineWidth',LW,'MarkerSize',MS)
% end
% 
% plot(numElements{h}(:,1),100*L1{h}(:,1)/numSamplesActual(h)*100,'k--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,1),100*Linf{h}(:,1)/numSamplesActual(h)*100,'k--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,1),100*L2{h}(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,1),100*L2{h}(:,1),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% 
% for i = 1:length(markerType)
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*L1{h}(:,j)/numSamplesActual(h)*100,markerType2{i},'Color',color{i},'LineWidth',LW,'MarkerSize',MS)
% end
% for i = 1:length(markerType)
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*Linf{h}(:,j)/numSamplesActual(h)*100,markerType2{i},'Color',color{i},'LineWidth',LW,'MarkerSize',MS)
% end
% for i = 1:length(markerType)
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*L2{h}(:,j),markerType{i},'MarkerFaceColor',color{i},'LineWidth',LW,'MarkerSize',MS)
% end
% for i = 1:length(markerType)
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*L2{h}(:,j),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% end
% 
% legend({['Analytic Poly.'],...
%     ['P1'],...
%     ['P2'],...
%     ['Piecewise']},...
%     'interpreter','latex',...
%     'FontSize',FS,...
%     'Location',plotLocation)
% ylabel('error $\%$','interpreter','latex','FontSize',FS)
% xlabel('number of elements','interpreter','latex','FontSize',FS)
% set(gcf,'Color',[1,1,1]);
% set(gca,'FontSize',FS)
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'xscale','log');
% set(gca,'yscale','log');
% axis tight
% grid on
% 
% box on
% 
% [f1, cdata] = myaa([4, 2]); imwrite(cdata, plotName, 'png');
% 
% close all
% 
% 
% 
% 
% 
% 
% %
% %--------------------------------------------------------------------------
% h=1
% markerType = {'ks-','ko-','kv-','kp-','kh-'};
% markerType2 = {'k--','k--','k--','k--','k--','k--','k--','k--'};
% color = {'r','y','c','b','m',};
% plotIndices = [13,14,15,16,17];
% 
% plotName = ['masconIterationMethod_h=',num2str(altitudes(h)),'_',body,'.png']
% 
% figure(3)
% hold on
% plot(numElements{h}(:,1),100*L2{h}(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% 
% 
% for i = 1:length(markerType)
%     i
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*L2{h}(:,j),markerType{i},'MarkerFaceColor',color{i},'LineWidth',LW,'MarkerSize',MS)
% end
% 
% plot(numElements{h}(:,1),100*L1{h}(:,1)/numSamplesActual(h)*100,'k--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,1),100*Linf{h}(:,1)/numSamplesActual(h)*100,'k--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,1),100*L2{h}(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,1),100*L2{h}(:,1),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% 
% for i = 1:length(markerType)
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*L1{h}(:,j)/numSamplesActual(h)*100,markerType2{i},'Color',color{i},'LineWidth',LW,'MarkerSize',MS)
% end
% for i = 1:length(markerType)
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*Linf{h}(:,j)/numSamplesActual(h)*100,markerType2{i},'Color',color{i},'LineWidth',LW,'MarkerSize',MS)
% end
% for i = 1:length(markerType)
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*L2{h}(:,j),markerType{i},'MarkerFaceColor',color{i},'LineWidth',LW,'MarkerSize',MS)
% end
% for i = 1:length(markerType)
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*L2{h}(:,j),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% end
% 
% legend({['Analytic Poly.'],...
%     ['P2-1/2'],...
%     ['P2-1/3'],...
%     ['Fine-1/2'], ...
%     ['Fine-1/3'], ...
%     ['Fine-1/4']},...
%     'interpreter','latex',...
%     'FontSize',FS,...
%     'Location',plotLocation)
% ylabel('error $\%$','interpreter','latex','FontSize',FS)
% xlabel('number of elements','interpreter','latex','FontSize',FS)
% set(gcf,'Color',[1,1,1]);
% set(gca,'FontSize',FS)
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'xscale','log');
% set(gca,'yscale','log');
% axis tight
% newXLim = xlim
% newXLim(1) = max(80,newXLim(1));
% set(gca,'xlim',newXLim)
% grid on
% 
% box on
% 
% [f1, cdata] = myaa([4, 2]); imwrite(cdata, plotName, 'png');
% 
% close all
% 
% % and at r = Rmax
% h=10
% markerType = {'ks-','ko-','kv-','kp-','kh-'};
% markerType2 = {'k--','k--','k--','k--','k--','k--','k--','k--'};
% color = {'r','y','c','b','m',};
% plotIndices = [13,14,15,16,17];
% 
% plotName = ['masconIterationMethod_h=',num2str(altitudes(h)),'_',body,'.png']
% 
% figure(3)
% hold on
% plot(numElements{h}(:,1),100*L2{h}(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% 
% 
% for i = 1:length(markerType)
%     i
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*L2{h}(:,j),markerType{i},'MarkerFaceColor',color{i},'LineWidth',LW,'MarkerSize',MS)
% end
% 
% plot(numElements{h}(:,1),100*L1{h}(:,1)/numSamplesActual(h)*100,'k--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,1),100*Linf{h}(:,1)/numSamplesActual(h)*100,'k--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,1),100*L2{h}(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
% plot(numElements{h}(:,1),100*L2{h}(:,1),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% 
% for i = 1:length(markerType)
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*L1{h}(:,j)/numSamplesActual(h)*100,markerType2{i},'Color',color{i},'LineWidth',LW,'MarkerSize',MS)
% end
% for i = 1:length(markerType)
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*Linf{h}(:,j)/numSamplesActual(h)*100,markerType2{i},'Color',color{i},'LineWidth',LW,'MarkerSize',MS)
% end
% for i = 1:length(markerType)
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*L2{h}(:,j),markerType{i},'MarkerFaceColor',color{i},'LineWidth',LW,'MarkerSize',MS)
% end
% for i = 1:length(markerType)
%     j = plotIndices(i);
%     plot(numElements{h}(:,j),100*L2{h}(:,j),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
% end
% 
% legend({['Analytic Poly.'],...
%     ['P2-1/2'],...
%     ['P2-1/3'],...
%     ['Fine-1/2'], ...
%     ['Fine-1/3'], ...
%     ['Fine-1/4']},...
%     'interpreter','latex',...
%     'FontSize',FS,...
%     'Location',plotLocation)
% ylabel('error $\%$','interpreter','latex','FontSize',FS)
% xlabel('number of elements','interpreter','latex','FontSize',FS)
% set(gcf,'Color',[1,1,1]);
% set(gca,'FontSize',FS)
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'xscale','log');
% set(gca,'yscale','log');
% axis tight
% newXLim = xlim
% newXLim(1) = max(80,newXLim(1));
% set(gca,'xlim',newXLim)
% 
% grid on
% 
% box on
% 
% [f1, cdata] = myaa([4, 2]); imwrite(cdata, plotName, 'png');
% 
% close all


% Plot quad methods
%--------------------------------------------------------------------------
h=1
markerType = {'ko-','kv-','ks-'};
markerType2 = {'k--','k--','k--'};
color = {'r',[1,0.5,0],'c'};
plotIndices = [2,18,14];   
% 2 -- simple packing
% 18-- extend tet 1
% 14 -- P2 layers

plotName = ['crossComparison_h=',num2str(altitudes(h)),'_',body,'.png']

figure(4)
hold on
plot(numElements{h}(:,1),100*L2{h}(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)


for i = 1:length(markerType)
    i
    j = plotIndices(i);
    plot(numElements{h}(:,j),100*L2{h}(:,j),markerType{i},'MarkerFaceColor',color{i},'LineWidth',LW,'MarkerSize',MS)
end
plot(numElements{h}(:,1),100*L1{h}(:,1)/numSamplesActual(h)*100,'k--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
plot(numElements{h}(:,1),100*Linf{h}(:,1)/numSamplesActual(h)*100,'k--','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
plot(numElements{h}(:,1),100*L2{h}(:,1),'kd-','MarkerFaceColor','k','LineWidth',LW,'MarkerSize',MS)
plot(numElements{h}(:,1),100*L2{h}(:,1),'k.','Color',[1,1,1],'MarkerSize',MSsmall)

for i = 1:length(markerType)
    j = plotIndices(i);
    plot(numElements{h}(:,j),100*L1{h}(:,j)/numSamplesActual(h)*100,markerType2{i},'Color',color{i},'LineWidth',LW,'MarkerSize',MS)
end
for i = 1:length(markerType)
    j = plotIndices(i);
    plot(numElements{h}(:,j),100*Linf{h}(:,j)/numSamplesActual(h)*100,markerType2{i},'Color',color{i},'LineWidth',LW,'MarkerSize',MS)
end
for i = 1:length(markerType)
    j = plotIndices(i);
    plot(numElements{h}(:,j),100*L2{h}(:,j),markerType{i},'MarkerFaceColor',color{i},'LineWidth',LW,'MarkerSize',MS)
end
for i = 1:length(markerType)
    j = plotIndices(i);
    plot(numElements{h}(:,j),100*L2{h}(:,j),'k.','Color',[1,1,1],'MarkerSize',MSsmall)
end

legend({['Analytic Poly.'],...
    ['PC Packing'],...
    ['Extend Tet - P1'],...
    ['FVM - offset surfaces - P2-1/3']},...
    'interpreter','latex',...
    'FontSize',FS,...
    'Location','SouthWest')
ylabel('error $\%$','interpreter','latex','FontSize',FS)
xlabel('number of elements','interpreter','latex','FontSize',FS)
set(gcf,'Color',[1,1,1]);
set(gca,'FontSize',FS)
set(gca,'TickLabelInterpreter','latex')
set(gca,'xscale','log');
set(gca,'yscale','log');
axis tight
grid on

box on

[f1, cdata] = myaa([4, 2]); imwrite(cdata, plotName, 'png');

close all
