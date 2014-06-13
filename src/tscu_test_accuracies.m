clear all
close all
clc
addpath('lib/export_fig')
x = load('accuracies.txt');

figure
% NONE vs DTW
d = x(:,2)-x(:,1); % DTW - NONE
snot = find(x(:,5)==0);
s    = find(x(:,5)==1);

plot(sort(d(s)),s,'ks','MarkerSize',12); hold on; plot(sort(d(snot)),snot,'k+','MarkerSize',12);
plot([0 0],[0 41],'k--'); legend('Significant','Other','Location','NorthWest'); ylim([0 41]); xlim([-0.20 0.20]);
ylabel('Dataset index','FontSize',12); xlabel('Difference between accuracies','FontSize',12);
set(gca,'YTick',[]);
title('NONE .vs. DTW');
set(gca,'XTick',[-0.1 0 0.1],'XTickLabel',{'DTW < NONE',' DTW = NONE','DTW > NONE'},'FontSize',10);
export_fig('-pdf','-transparent','tscu_test_accuracies_NONE_vs_DTW.pdf');

figure
% NONE vs CDTW
d = x(:,3)-x(:,1); % CDTW - NONE
snot = find(x(:,6)==0);
s    = find(x(:,6)==1);

plot(sort(d(s)),s,'ks','MarkerSize',12); hold on; plot(sort(d(snot)),snot,'k+','MarkerSize',12);
plot([0 0],[0 41],'k--'); legend('Significant','Other','Location','NorthWest'); ylim([0 41]); xlim([-0.20 0.20]);
ylabel('Dataset index','FontSize',12); xlabel('Difference between accuracies','FontSize',12);
set(gca,'YTick',[]);
title('NONE .vs. CDTW');
set(gca,'XTick',[-0.1 0 0.1],'XTickLabel',{'CDTW < NONE','CDTW = NONE','CDTW > NONE'},'FontSize',10);
export_fig('-pdf','-transparent','tscu_test_accuracies_NONE_vs_CDTW.pdf');

figure
% NONE vs SAGA
snot = find(x(:,7)==0);
s    = find(x(:,7)==1);
d = x(:,4)-x(:,1); % SAGA - NONE
plot(sort(d(s)),s,'ks','MarkerSize',12); hold on; plot(sort(d(snot)),snot,'k+','MarkerSize',12);
plot([0 0],[0 41],'k--'); legend('Significant','Other','Location','NorthWest'); ylim([0 41]); xlim([-0.20 0.20]);
ylabel('Dataset index','FontSize',12); xlabel('Difference between accuracies','FontSize',12);
set(gca,'YTick',[]);
title('NONE .vs. SAGA');
set(gca,'XTick',[-0.1 0 0.1],'XTickLabel',{'SAGA < NONE','SAGA = NONE','SAGA > NONE'},'FontSize',10);
export_fig('-pdf','-transparent','tscu_test_accuracies_NONE_vs_SAGA.pdf');

figure
% DTW vs CDTW
snot = find(x(:,8)==0);
s    = find(x(:,8)==1);
d = x(:,3)-x(:,2); % CDTW-DTW
plot(sort(d(s)),s,'ks','MarkerSize',12); hold on; plot(sort(d(snot)),snot,'k+','MarkerSize',12);
plot([0 0],[0 41],'k--'); legend('Significant','Other','Location','NorthWest'); ylim([0 41]); xlim([-0.20 0.20]);
ylabel('Dataset index','FontSize',12); xlabel('Difference between accuracies','FontSize',12);
set(gca,'YTick',[]);
title('DTW .vs. CDTW');
set(gca,'XTick',[-0.1 0 0.1],'XTickLabel',{'CDTW < DTW','CDTW = DTW','CDTW > DTW'},'FontSize',10);
export_fig('-pdf','-transparent','tscu_test_accuracies_CDTW_vs_DTW.pdf');

figure
% DTW vs SAGA
snot = find(x(:,9)==0);
s    = find(x(:,9)==1);
d = x(:,4)-x(:,2); % SAGA-DTW
plot(sort(d(s)),s,'ks','MarkerSize',12); hold on; plot(sort(d(snot)),snot,'k+','MarkerSize',12);
plot([0 0],[0 41],'k--'); legend('Significant','Other','Location','NorthWest'); ylim([0 41]); xlim([-0.20 0.20]);
ylabel('Dataset index','FontSize',12); xlabel('Difference between accuracies','FontSize',12);
set(gca,'YTick',[]);
title('DTW .vs. SAGA');
set(gca,'XTick',[-0.1 0 0.1],'XTickLabel',{'SAGA < DTW','SAGA = DTW','SAGA > DTW'},'FontSize',10);
export_fig('-pdf','-transparent','tscu_test_accuracies_DTW_vs_SAGA.pdf');

figure
% CDTW vs SAGA
snot = find(x(:,10)==0);
s    = find(x(:,10)==1);
d = x(:,4)-x(:,3); % SAGA-CDTW
plot(sort(d(s)),s,'ks','MarkerSize',12); hold on; plot(sort(d(snot)),snot,'k+','MarkerSize',12);
plot([0 0],[0 41],'k--'); legend('Significant','Other','Location','NorthWest'); ylim([0 41]); xlim([-0.20 0.20]);
ylabel('Dataset index','FontSize',12); xlabel('Difference between accuracies','FontSize',12);
set(gca,'YTick',[]);
title('CDTW .vs. SAGA');
set(gca,'XTick',[-0.1 0 0.1],'XTickLabel',{'SAGA < CDTW','SAGA = CDTW','SAGA > CDTW'},'FontSize',10);
export_fig('-pdf','-transparent','tscu_test_accuracies_CDTW_vs_SAGA.pdf');