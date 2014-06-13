clear all
close all
clc

addpath('lib/export_fig');


integration    = @(s) cumsum(s)/(size(s,1)-1);
monotonization = @(s) integration(exp(s));
normalize      = @(s) (s-s(1))/(s(end)-s(1));

nk=20;
o = zeros(nk,1);
l = ones(nk,1);
phi1 =[ l  o  o  o;
       -l  l  o  o;
        o -l  l  o;
        o  o -l  l;
        o  o  o -l];
phi2 =[ l  o  o  o;
        o  l  o  o;
        o  o  l  o;
        o  o  o  l];
phi8old = [l o o o o o o o;
           o l o o o o o o;
           o o l o o o o o;
           o o o l o o o o;
           o o o o l o o o;
           o o o o o l o o;
           o o o o o o l o;
           o o o o o o o l];
phi8new = [l  o  o  o  o  o  o  o;
          -l  l  o  o  o  o  o  o;
           o -l  l  o  o  o  o  o;
           o  o -l  l  o  o  o  o;
           o  o  o -l  l  o  o  o;
           o  o  o  o -l  l  o  o;
           o  o  o  o  o -l  l  o;
           o  o  o  o  o  o -l  l];
       
phi = phi1;    
K = size(phi,2);
n = size(phi,1);

t = linspace(0,1,n)';

varphi = integration(phi);
Phi = integration(varphi);
Beta = zeros(n*K,K);
for i=1:K
    for j=1:K
        Betaij =  integration(integration(phi(:,i)).*integration(phi(:,j)));
        start = (i-1)*n+1;
        stop  = start+n-1;
        Beta(start:stop,j)= Betaij;
    end
end

solve_ode_direct    = @(c) normalize(monotonization(integration(phi*c')));
solve_ode_taylor1st = @(c) normalize(t+Phi*c');
solve_ode_taylor2nd = @(c,n,k) normalize(t+Phi*c' + ...
    0.5*(sum(reshape(reshape(repmat(c,n,1),n*k,1).*(Beta*c'),n,k),2)));

%solve_ode_taylor2nd2 = @(c,n,k) normalize(t+Phi*c' + ...
%    0.5*(reshape(repmat(c,n*n,1),n,n*k)*Beta*c'));
% Direct Taylor1st Taylor2nd
mtest=1000;
errors = zeros(mtest,3);
etimes = zeros(mtest,3);
u_direct = zeros(n,mtest);
u_taylor1st = zeros(n,mtest);
u_taylor2nd = zeros(n,mtest);
for i=1:mtest
    c = 2*rand(1,K)-1;
    
    tic
    u_direct(:,i)=solve_ode_direct(c); 
    etimes(i,1)=toc;
 
    tic
    u_taylor1st(:,i)=solve_ode_taylor1st(c);
    etimes(i,2)=toc;
 
    tic
    u_taylor2nd(:,i)=solve_ode_taylor2nd(c,n,K);
    etimes(i,3)=toc;
    
    

    plot(u_direct(:,i),'k'); 
    hold on; 
    plot(u_taylor1st(:,i),'r');
    plot(u_taylor2nd(:,i),'g');
    
    legend('Direct','1st','2nd','Location','SouthEast');
    hold off
    pause(0.1)
end

%%
close all

figure
plot(t,u_direct(:,2),'k');
hold on;
plot(t,u_taylor1st(:,2),'k*');
plot(t,u_taylor2nd(:,2),'k.');
legend('Direct','1st Taylor','2nd Taylor','Location','SouthEast');
xlabel('model time','FontSize',12);
ylabel('observed time','FontSize',12);
export_fig('-pdf','-transparent','tscu_test_taylor_warpings.pdf');

figure
plot(t,t-u_direct(:,2),'k');
hold on;
plot(t,t-u_taylor1st(:,2),'k*');
plot(t,t-u_taylor2nd(:,2),'k.');
legend('Direct','1st Taylor','2nd Taylor','Location','SouthEast');
xlabel('model time','FontSize',12);
ylabel('observed time - t','FontSize',12);
export_fig('-pdf','-transparent','tscu_test_taylor_warpings_minus_t.pdf');

figure
idx = round(linspace(1,n,20));
boxplot(abs(u_direct(idx,:)-u_taylor1st(idx,:))');
export_fig('-pdf','-transparent','tscu_test_taylor_boxplot1st.pdf');

figure
boxplot(abs(u_direct(idx,:)-u_taylor2nd(idx,:))');
export_fig('-pdf','-transparent','tscu_test_taylor_boxplot2nd.pdf');

figure
plot(phi); title('\phi');
export_fig('-pdf','-transparent','tscu_test_taylor_phi.pdf');

figure
plot(varphi); title('varphi');
export_fig('-pdf','-transparent','tscu_test_taylor_varphi.pdf');

figure
plot(Phi);title('\Phi');
export_fig('-pdf','-transparent','tscu_test_taylor_Phi.pdf');

figure
plot(reshape(Beta,n,K*K)); title('\beta');
export_fig('-pdf','-transparent','tscu_test_taylor_beta.pdf');



