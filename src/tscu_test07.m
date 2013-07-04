n = 100;
t = linspace(0,1,n);
x = sin(t);
y = cos(t);

m=50000;
d = zeros(1,m);

tic
for i=1:m
s = -1+rand(1,8)*2;
d(i) = Jcost1(y,x,s);
end
fprintf('%8f %8.2f\n',sum(d),toc);
