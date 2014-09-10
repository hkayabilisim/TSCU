clear all
close all
clc
addpath('lib/export_fig')
tmp=load('results_svm.mat','a'); % load NONE, DTW, CDTW,SAGA
x=tmp.a;

cname = {'NONE','DTW','CDTW','SAGA'};
c = [1 2;1 3;1 4;2 3;2 4;3 4];
for j=1:size(c,1) 
    m1=cname{c(j,1)};
    m2=cname{c(j,2)};
    figure
    % NONE vs DTW means DTW-NONE
    d = x(:,c(j,2))-x(:,c(j,1)); % DTW - NONE
    snot = find(x(:,4+j)==0);
    s    = find(x(:,4+j)==1);
    limx = max(abs(d));
    plot(sort(d(s)),s,'ks','MarkerSize',12);
    hold on;
    plot(sort(d(snot)),snot,'k+','MarkerSize',12);
    plot([0 0],[0 41],'k--');
    %legend('Significant','Other','Location','NorthWest');
    ylim([0 41]); xlim([-limx limx]);
    %ylabel('Dataset index','FontSize',12);
    %xlabel('Difference between accuracies','FontSize',12);
    set(gca,'YTick',[]);
    title(sprintf('%4s .vs. %4s',m1,m2));
    
    set(gca,'XTick',[-0.5 0 0.5],'XTickLabel',...
        {sprintf('%4s < %4s',m2,m1),...
         sprintf('%4s = %4s',m2,m1),...
         sprintf('%4s > %4s',m2,m1)},'FontSize',10);
    export_fig('-pdf','-transparent',...
        sprintf('tscu_test_accuracies_svm_%s_vs_%s.pdf',m1,m2));
    
end