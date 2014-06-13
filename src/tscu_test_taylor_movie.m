clear all
close all
clc

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

mtest=1000;
errors = zeros(mtest,2);
etimes = zeros(mtest,3);
u_direct = zeros(n,mtest);
u_taylor1st = zeros(n,mtest);
u_taylor2nd = zeros(n,mtest);
for i=1:mtest
    c = 4*rand(1,K)-2;
    
    tic
    u_direct(:,i)=solve_ode_direct(c); 
    etimes(i,1)=toc;
 
    tic
    u_taylor1st(:,i)=solve_ode_taylor1st(c);
    etimes(i,2)=toc;
 
    tic
    u_taylor2nd(:,i)=solve_ode_taylor2nd(c,n,K);
    etimes(i,3)=toc;

end
%%
errors(:,1)=(sqrt(sum((u_direct-u_taylor1st).^2))./(sqrt(sum(u_direct.^2))));
errors(:,2)=(sqrt(sum((u_direct-u_taylor2nd).^2))./(sqrt(sum(u_direct.^2))));
%%
close all    
figure; %('Position',[500,500,1200,500]);
for i=1:mtest    
subplot(2,2,[1 3]);
plot(t,u_taylor1st(:,i),'r','LineWidth',1);
hold on;
plot(t,u_taylor2nd(:,i),'g','LineWidth',1);
plot(t,u_direct(:,i),'k','LineWidth',1); 
set(gca,'XTick',[],'YTick',[]);
h_legend=legend('Taylor 1st','Taylor 2nd','Direct','Location','SouthEast');
set(h_legend,'FontSize',14);
hold off
axis equal
box on

subplot(2,2,2);
hist(errors(1:i,1),20);
h_legend=legend('Taylors 1st Relative Error');
set(h_legend,'FontSize',10);
set(gca,'XTick',[0 0.01 0.02],'XTickLabel',{'0','%1','%2'});
xlim([0 0.02]);

subplot(2,2,4);
hist(errors(1:i,2),20);
h_legend=legend('Taylors 2nd Relative Error');
set(h_legend,'FontSize',10);
set(gca,'XTick',[0 0.001 0.002],'XTickLabel',{'0','%0.1','%0.2'});
xlim([0 0.002]);
print('-dpng',sprintf('taylor_movie/taylor_movie%06d.png',i));

end


