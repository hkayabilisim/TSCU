addpath('../fdaM')

T = 4;
M = 4;
B = 4;

[p,e,t] = ClipTriPET(T, M, B);

np = size(p,2);

pdemesh(p,e,t)
hold on
plot([0,0],[0,T],'g-')
hold off

ngrid = 101;
x = linspace(0,4,ngrid)';
v = linspace(-4,4,ngrid)';
y = linspace(0,4,ngrid)';

%  plot each basis function

for i=1:np
    ui = zeros(np, 1);
    ui(i) = 1;
    [uxyi,tn,al2,al3] = tri2grid(p, t, ui, x, y);
    uxyi(isnan(uxyi)) = 0;
    surf(x,y,uxyi)
    axis([0,4,0,4,0,1])
    title(['\fontsize{16} Basis function ',num2str(i)])
    pause
end
 
%  plot each diagonal crossproduct function

for i=1:np
    ui = zeros(np, 1);
    ui(i) = 1;
    [uxyi,tn,al2,al3] = tri2grid(p, t, ui, v, y);
    uxyi(isnan(uxyi)) = 0;
    cxyii = ((uxyi*uxyi') - ...
        0.5.*(uxyi(:,1    )*uxyi(:,1    )' + ...
        uxyi(:,ngrid)*uxyi(:,ngrid)'))./(ngrid-1);
    surf(x,y,cxyii)
    axis([0,4,0,4,0,0.12])
    title(['\fontsize{16} Diagonal product ',num2str(i)])
    pause
end

%  set up surface

u = ones(np,1);
% u(6:9) = 0.5;
% u(10:12) = 0.25;

% tn  = x;
% al2 = y;

[uxy,tn,al2,al3] = tri2grid(p, t, u, v, y);

uxy(isnan(uxy))=0;

surf(x,y,uxy*uxy'./(ngrid-1))

contour(uxy)

size(isnan(uxy))

uxy(isnan(uxy)) = 0; 

cxy = uxy*uxy';

surf(cxy)

m = 0;
for i=1:np
    ui = zeros(np, 1);
    ui(i) = 1;
    [uxyi,tn,al2,al3] = tri2grid(p, t, ui, v, y);
    uxyi(isnan(uxyi)) = 0;
    for j=1:i
        if abs(p(1,i)-p(1,j)) < 2
        uj = zeros(np, 1);
        uj(j) = 1;
        [uxyj,tn,al2,al3] = tri2grid(p, t, uj, v, y); 
        uxyj(isnan(uxyj)) = 0;
        cxyij = ((uxyi*uxyj') - ...
            0.5.*(uxyi(:,1    )*uxyj(:,1    )' + ...
                  uxyi(:,ngrid)*uxyj(:,ngrid)'))./(ngrid-1);
        surf(x,y,cxyij)
        axis([0,4,0,4])
        title(['\fontsize{16}',num2str(i),' - ',num2str(j)])
        m = m + 1;
        disp([i,j,m])
        pause
        end
    end
end
 
   
  
 
