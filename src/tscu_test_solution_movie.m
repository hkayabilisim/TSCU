function tscu_test_solution_movie
close all
clc
tic
k=8;
s=zeros(k,1);
nk=10;
k=length(s);
n=nk*k;
t=linspace(0,1,n)';
maxwalk=1000;
walk=zeros(maxwalk,k);
for i=1:maxwalk
    walk(i,:)=s';
    [w0,w1,w2,w3]=solve_ode(k,n,nk,s); 
    direction = randi([1 k],1,1);
    s(direction)=s(direction)+0.1*(2*rand(1,1)-1);
    plot_ode(i,t,walk,w0,w1,w2,w3);
    pause(0.01)
end
toc
end
function [w0,w1,w2,w3]=solve_ode(k,n,nk,s)

w0=zeros(n,1);
for i=0:k-1
    w0(i*nk+1:i*nk+nk)=s(i+1);    
end
%w0=[10*ones(50,1); -10*ones(50,1)]; %20*sin(10*pi*t);

w1=cumsum(w0)/(n-1);
w2=exp(w1);
w3=cumsum(w2)/(n-1); w3=w3/w3(end);

end

function plot_ode(i,t,walk,w0,w1,w2,w3)
subplot(2,3,1); plot(w0,'k.'); title('w0: sin(t)'); ylim([-5 5]);
subplot(2,3,2); plot(t,w1); title('w1: int(w0)');
subplot(2,3,4); plot(t,w2); title('w2: exp(w1)');
subplot(2,3,5); plot(walk(1:i,1),walk(1:i,2),'k'); title('random walk');
%xlim([-2 2]); ylim([-2 2]); axis equal
subplot(2,3,[3 6]); plot(t,w3,'k'); title('w3: int(w2)'); 
%xlim([0 1]); ylim([0 1]); axis equal
end

