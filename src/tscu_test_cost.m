clear all
close all
clc
%% cost surface
trn=load('../../UCR/FaceFour/FaceFour_TRAIN');
tst=load('../../UCR/FaceFour/FaceFour_TEST');

x = trn(1,2:end);% class 2
y = tst(3,2:end);% class 2
n = length(x);
m = 50;
w1 = linspace(-10,10,m);
w2 = linspace(-10,10,m);
[W1,W2]=meshgrid(w1,w2);
C = zeros(m,m);
for i=1:m
    for j=1:m
        C(i,j) = tscu_saga_cost(x,y,[w1(i) w2(j)],zeros(1,n),zeros(1,n));
    end
end

% [p1,p2,d]=tscu_saga_register(x,y,2,zeros(1,n),zeros(1,n),zeros(1,2),zeros(1,2));
% 
% x_saga = interp1(1:n,x,p1);
% y_saga = interp1(1:n,y,p2);
% 
% [y_aligned,w]=tscu_saga_warp(y,[w1(24),w2(27)]);
%%
surfc(W1,W2,C)
colormap gray
xlabel('c_1'); ylabel('c_2'); zlabel('J');
%title('Surface plot of a cost function');
export_fig('-pdf','-transparent',...
    sprintf('tscu_test_cost.pdf'));

% figure
% subplot(211)
% plot(x,'k'); hold on
% plot(y,'k--')
% legend('x','y');
% title(sprintf('Original time series [d=%8.5f]',norm(x-y)));
% 
% subplot(212)
% plot(x_saga,'k'); hold on
% plot(y_saga,'k--')
% %plot(y_aligned,'k.');
% legend('x','y'); %,'y_manuel');
% title(sprintf('Aligned time series [d=%8.5f]',norm(x_saga-y_saga)));

