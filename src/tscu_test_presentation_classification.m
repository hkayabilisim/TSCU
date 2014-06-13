clear all
close all
clc

addpath('lib/export_fig/');

trn=load('../../UCR/synthetic_control/synthetic_control_TRAIN');
tst=load('../../UCR/synthetic_control/synthetic_control_TEST');

for i=1:6
   idx = find(trn(:,1)==i); 
   plot(trn(idx(1:10),2:end)','k');
   set(gca,'XTick',[],'YTick',[]);
   export_fig('-pdf','-transparent',sprintf('tscu_test_presentation_classification%d.pdf',i));

end

x = tst(14,2:end);
plot(x,'k');
set(gca,'XTick',[],'YTick',[]);
text(30,0,'?','FontSize',200,'HorizontalAlignment','center','VerticalAlignment','middle');
export_fig('-pdf','-transparent','tscu_test_presentation_classification_unknown.pdf');
