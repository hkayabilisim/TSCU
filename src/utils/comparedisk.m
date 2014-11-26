clear all
close all
clc


mnames = {'NONE' ,'DTW'  ,'CDTW','SAGA','PTW'};
mcolor = {'r'    ,'g'    ,'b'   ,'m'   ,'c'  };

for i=1:40
    dnames{i}=num2str(i);
end


%input = [88.0  99.3  98.7  90.3  94.3;...
%         76.0 100.0  99.0  99.0  99.0];
input = load('input.mat');
input = input.input;
bez = @(t,P) ...
  bsxfun(@times,(1-t).^3,P(1,:)) + ...
  bsxfun(@times,3*(1-t).^2.*t,P(2,:)) + ...
  bsxfun(@times,3*(1-t).^1.*t.^2,P(3,:)) + ...
  bsxfun(@times,t.^3,P(4,:));

a=linspace(0,2*pi,length(mnames)*length(dnames)*2);

wh = @(mi,di) (mi-1)*2*length(dnames)+di ;

plot(cos(a),sin(a),'k')
axis equal
hold on;
alpha(0.5);

for i=1:length(mnames)
    text(1.1*cos(a(wh(i,20))),1.1*sin(a(wh(i,20))),mnames{i});
end
    
for i=1:length(mnames)-1
    for j=i+1:length(mnames)
        for k=1:length(dnames)
            fprintf('%20s %5.1f (%5s) %5.1f (%5s)\n',...
                dnames{k},input(k,i),mnames{i},input(k,j),mnames{j});
            nn=400;
            xx=((abs((input(k,i)-input(k,j)))/nn)^(1/3)+1)/2;
            if input(k,j) > input(k,i)
                xx=1-xx;
            end
            t=linspace(0,xx,nn)';
            from = [cos(a(wh(i,k))) sin(a(wh(i,k)))];
            to   = [cos(a(wh(j,k))) sin(a(wh(j,k)))];
            if mod(abs(i-j),3) == 1
                pp = (from+to)/norm(from+to)*(1/2);
            else
                pp = [0 0];
            end
            % From i --> j
            P=[from(1) from(2);...
               pp(1)   pp(2)  ;...
               pp(1)   pp(2)  ;
               to(1)   to(2) ];
            X = bez(t,P);
            patchline(X(:,1),X(:,2),'edgecolor',mcolor{i},'linewidth',2,...
                'edgealpha',0.3);
            % From i --> j
            t=linspace(xx,1,nn)';
            X = bez(t,P);
            patchline(X(:,1),X(:,2),'edgecolor',mcolor{j},'linewidth',2,...
                'edgealpha',0.3);
            %pause(0.01);
        end
    end
end
     





