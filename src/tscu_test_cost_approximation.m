clear all
close all
clc

integration    = @(s) cumsum(s)/(size(s,1)-1);
monotonization = @(s) integration(exp(s));
normalize      = @(s) (s-s(1))/(s(end)-s(1));

solve_ode_real = @(c,u)     normalize(monotonization(integration(u*c')));
solve_ode_appr = @(c,u,t) normalize(t+integration(integration(u))*c');

%u = [ones(1,50)  , -1*ones(1,50),    zeros(1,100); ...
%     zeros(1,100),  1*ones(1,50), -1*ones(1,50)]';
o = zeros(10,1);
l = ones(10,1);
u = [ l o o o o o o o;
      o l o o o o o o;
      o o l o o o o o;
      o o o l o o o o;
      o o o o l o o o;
      o o o o o l o o;
      o o o o o o l o;
      o o o o o o o l];

n = size(u,1);
t = linspace(0,1,n)'; 
ntest=10000;
e=zeros(ntest,1);
tic
for i=1:ntest
c = 2*rand(1,8)-1;
r=solve_ode_real(c,u);
end
fprintf('Real: %8.2f                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                second\n',toc);

tic
for i=1:ntest
    c = 2*rand(1,8)-1;

a=solve_ode_appr(c,u,t);
end
fprintf('Appr: %8.2f second\n',toc);

% e(i)=norm(a-r)/norm(r);
% 
% subplot(211);
% plot(t,r); hold on;
% plot(t,a,'r'); hold off;
% xlim([0 1]); ylim([0 1]);
% subplot(212);
% plot(e(1:i));
