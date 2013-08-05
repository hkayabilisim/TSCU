clear all
close all
clc

addpath('../../lib/export_fig');
%trn=load('../../UCR/DiatomSizeReduction/DiatomSizeReduction_TRAIN');
%tst=load('../../UCR/DiatomSizeReduction/DiatomSizeReduction_TEST');
% 
% tscu(trn,tst,'Alignment','SAGA','MATLABPool','karadeniz_64');
% tscu(trn,tst,'Alignment','SAGA','MATLABPool','karadeniz_32');
% tscu(trn,tst,'Alignment','SAGA','MATLABPool','karadeniz_16');
% tscu(trn,tst,'Alignment','SAGA','MATLABPool','karadeniz_8');
% tscu(trn,tst,'Alignment','SAGA','MATLABPool','karadeniz_4');
% tscu(trn,tst,'Alignment','SAGA','MATLABPool','karadeniz_2');


nps   = [1 2 4 8 16 32 64];
times = [5146.74 2573.37 1305.21 690.90 344.25 172.19 87.03];
speedups = times(1)./times;
exptimes = times(1)*ones(1,length(nps))./nps;
obttimes = times;
efficieny = speedups./nps;

figure
plot(log2(nps),log2(speedups),'--bs')
hold on
plot(log2(nps),log2(nps),'--rs')
legend('Obtained','Expected','Location','NorthWest');
xlabel('# of processors');
ylabel('speedup')
title('Speedup (log2/log2 scaled)');
set(gca,'XTickLabel',nps);
set(gca,'YTickLabel',nps);
export_fig('-pdf','-transparent','tscu_test08_speedup.pdf');

figure
plot(log2(nps),log2(obttimes),'--bs')
hold on
plot(log2(nps),log2(exptimes),'--rs')
legend('Obtained','Expected','Location','NorthEast');
xlabel('# of processors (log2 scaled)');
ylabel('log2 of elapsed time (sec)')
title('Elapsed times');
set(gca,'XTickLabel',nps);
%set(gca,'YTickLabel',nps);
%figuresize(15,15,'centimeters');
export_fig('-pdf','-transparent','tscu_test08_elapsed.pdf');

figure
plot(log2(nps),100*efficieny,'--bs')
xlabel('# of processors (log2 scaled)');
ylabel('efficieny')
title('Efficieny');
set(gca,'XTickLabel',nps);
%set(gca,'YTickLabel',nps);
%figuresize(15,15,'centimeters');
%print -dpdf 'Experiment42-efficieny.pdf'
export_fig('-pdf','-transparent','tscu_test08_efficiency.pdf');



