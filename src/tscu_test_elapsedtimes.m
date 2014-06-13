clear all
close all
clc

addpath('lib/export_fig');
%1  2       3   4   5       6    7   8    9
%id	length	trn	tst	trnxtst	NONE DTW CTDW SAGA
x=load('tscu_test_elapsedtimes.txt');
alignments = (x(:,3).*x(:,4));
n          = x(:,2);
[nsorted,nidx]=sort(n);
none = alignments./x(:,6);
dtw  = alignments./x(:,7);
cdtw = alignments./x(:,8);
saga = alignments./x(:,9);
lin  = 100000./nsorted;
quad = (10000./nsorted).^2;
loglog((nsorted),(none(nidx)),'ks','MarkerSize',12); hold on
loglog((nsorted),(dtw(nidx)),'k.','MarkerSize',12);
loglog((nsorted),(cdtw(nidx)),'kd','MarkerSize',12);
loglog((nsorted),(saga(nidx)),'k*','MarkerSize',12); 
%semilogx((nsorted),(lin),'k--');
%semilogx((nsorted),(quad),'k--');

%xlim([4 12])
%ylim([1.5 5])
xlabel('The length of time series in log scale','FontSize',12);
ylabel('# of alignments per second','FontSize',12);
legend('NONE','DTW','CDTW','SAGA','Location','NorthEast')
%set(gca,'XTick',5:11); %'XTickLabel',{10,100,1000});
%set(gca,'YTick',[2 3 4]);

export_fig('-pdf','-transparent','tscu_test_elapsedtimes.pdf');
