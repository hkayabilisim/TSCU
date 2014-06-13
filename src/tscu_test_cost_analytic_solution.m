clear all
close all
clc

integration    = @(s) cumsum(s)/(size(s,1)-1);
monotonization = @(s) integration(exp(s));
normalization  = @(s) (s-s(1))/(s(end)-s(1));

m = 5;
l = ones(m,1);
o = ones(m,1);

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
       
tt = linspace(0,1,m)';       
phifourier = zeros(m,11);
for i=1:5
    phifourier(:,i)  =sqrt(2)*cos(2*pi*i*tt);
    phifourier(:,2*i)=sqrt(2)*sin(2*pi*i*tt);
end
phifourier(:,end) = ones(m,1);

phi = phi8new;
[n k]=size(phi);
t = linspace(0,1,n)';

x = t.^2;
y = x - t;

Phi = integration(phi);   
varphi = integration(Phi);
Q = varphi'*varphi;
s = 1*((t.^2)-t);
c = (Q\(varphi'*s))';

z = normalization(t+varphi*c');
[path1,path2,d]=tscu_saga_register(x,t,k,zeros(1,n),zeros(1,n),zeros(1,k),zeros(1,k),varphi');
z2 = interp1(1:n,t,path2);

subplot(131);
plot(t,x,'k'); 
hold on
plot(t,t,'r');
legend('x','y');

subplot(132);
plot(t,x,'k'); 
hold on
plot(t,z,'r');
legend('x','y analytic');

subplot(133);
plot(t,x,'k'); 
hold on
plot(t,z2,'r');
legend('x','SAGA');



   
   
