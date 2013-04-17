function [u c d] = ramsay(t,w)
% Solution of ODE for a piecewise defined weight
% function
% t 1 x n : evaluation points
% w 1 x m : curvatures
% out:
% u 1 x n : values on t
% c 1 x m : Coefficients: c
% d 1 x m : Coefficients: d

m = length(w);
n = length(t);

A = zeros(2*m,2*m);
y = zeros(2*m,1);
u = zeros(1,n);
b = linspace(0,1,m+1);

% If the variables are off limits then
% we correct the variables.
w(w<-30) = -30;
w(w>10)  = 10;

eps = 0.001;
for i=1:m-1
    if abs(w(i)) < eps
        A(i,i)     = b(i+1); % ok
        A(i+m-1,i) = 1;  %ok
    else
        A(i,i)     = exp(w(i)*b(i+1)); % ok
        A(i+m-1,i) = w(i)*exp(w(i)*b(i+1)); %ok
    end
    
    if abs(w(i+1)) < eps
        A(i,i+1)     = -b(i+1); %ok
        A(i+m-1,i+1) = -1; %ok
    else
        A(i,i+1)     = -exp(w(i+1)*b(i+1)); %ok
        A(i+m-1,i+1) = -w(i+1)*exp(w(i+1)*b(i+1)); %ok
    end
    
    A(i,m+i) = 1; %ok
    A(i,m+i+1) = -1; %ok
end

if abs(w(1)) > eps
    A(2*m-1,1) = 1; % ok
end
A(2*m-1,m+1)=1; % ok
if abs(w(m)) < eps
    A(2*m,m) = 1; %ok
else
    A(2*m,m) = exp(w(m)); %ok
end
A(2*m,2*m)=1;

%fprintf('rcond: %e min: %8.2f max %8.2f\n',rcond(A),min(w),max(w));
% Right hand side of Ax=y
y(2*m) = 1;

x = (A\y); % solution of Ax=y
c = x(1:m)';
d = x(m+1:2*m)';

for i=1:m
    idx = find(t >= b(i) & t <= b(i+1));
    if abs(w(i)) < eps
        u(idx) = c(i)*t(idx)+d(i);
    else
        u(idx) = c(i)*exp(w(i)*t(idx))+d(i);
    end
end
u(1) = 0;
u(end) = 1;


end