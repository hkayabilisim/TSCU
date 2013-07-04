function d = Jcost0(y,x,s,t)
% %5     linspace
% %75    ramsay
% %75    interp1
%t = linspace(0,1,length(x));
d = norm(interp1(t,y,ramsay(t,s))-x);

%d = norm(ramsay(t,s)-x);
%d = norm(interp1(t,y,y)-x);
%d = norm(y-x);
end
