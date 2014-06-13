clear all
close all
clc

addpath('lib/export_fig');


ntest = [10 20 40];
u_real=cell(length(ntest),1);
u_direct=cell(length(ntest),1);
t_all=cell(length(ntest),1);
errors = zeros(length(ntest),1);

for i=1:length(ntest)
nhalf=ntest(i);
n=2*nhalf;
K = 2;
t = linspace(0,1,n)';
t1 = linspace(0,0.5,nhalf+1)';
t2 = linspace(0.5,1,nhalf)';

integration    = @(s) cumsum(s)/(size(s,1)-1);
monotonization = @(s) integration(exp(s));
normalize      = @(s) (s-s(1))/(s(end)-s(1));

o = zeros(nhalf,1);
l = ones(nhalf,1);
phi = [ l o ;
        o l ];
Phi = integration(integration(phi));
Beta11 =  integration(integration(phi(:,1)).*integration(phi(:,1)));
Beta12 =  integration(integration(phi(:,1)).*integration(phi(:,2)));
Beta21 =  integration(integration(phi(:,2)).*integration(phi(:,1)));
Beta22 =  integration(integration(phi(:,2)).*integration(phi(:,2)));

solve_ode_direct    = @(c) normalize(monotonization(integration(phi*c')));
solve_ode_taylor1st = @(c) normalize(t+Phi*c');
solve_ode_taylor2nd = @(c) normalize(t+Phi*c' + ...
    0.5*(c(1)*c(1)*Beta11+c(1)*c(2)*Beta12+c(2)*c(1)*Beta21+c(2)*c(2)*Beta22) );

    c1 =  10;
    c2 = -10;
    c = [c1 c2];

    e1 = exp(c1*0.5);
    e2 = exp(c2*0.5);
    a  = -e2-c2*e2*(1-e1)/(c1*e1);
    b2 = 1/(e2^2+a);
    a2 = b2*a;
    b1 = b2*c2*e2/(c1*e1);
    a1 = -b1;
    
    u1 = a1+b1*exp(c1*t1);
    u2 = a2+b2*exp(c2*t2);

    u_real{i} = [u1(1:end-1);u2];

    u_direct{i} = solve_ode_direct(c);
    t_all{i} = t;
    errors(i)=norm(u_direct{i}-u_real{i})/norm(u_real{i});
end

plot(t_all{3},u_real{3},'k');
hold on;
plot(t_all{1},u_direct{1},'k--');
plot(t_all{2},u_direct{2},'k-.');
plot(t_all{3},u_direct{3},'k-*');

xlim([-0.1 1.1])
ylim([-0.1 1.1])
set(gca,'XTick',[0 0.5 1],'YTick',[0 0.5 1])
xlabel('model time','FontSize',12);
ylabel('observed time','FontSize',12);
legend('Analytical',...
    sprintf('n=%-3d [Error: %3.0f%%]',ntest(1)*2,100*errors(1)),...
    sprintf('n=%-3d [Error: %3.0f%%]',ntest(2)*2,100*errors(2)),...
    sprintf('n=%-3d [Error: %3.0f%%]',ntest(3)*2,100*errors(3)),...
    'FontSize',12,'Location','SouthEast');
export_fig('-pdf','-transparent',...
    sprintf('tscu_test_ode.pdf'));