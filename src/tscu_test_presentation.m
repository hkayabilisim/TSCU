clear all
close all
clc

movie_counter=1;
addpath('lib/export_fig');
% trn = load('../../UCR/ECGFiveDays/ECGFiveDays_TRAIN');
% tst = load('../../UCR/ECGFiveDays/ECGFiveDays_TEST');
% 
% plot(trn(:,2:end)','k')
% set(gca,'XTick',[],'YTick',[]);
% box on
% 
% export_fig('-pdf','-transparent','tscu_test_presentation_ECGFiveDays.pdf');

trn = load('../../UCR/ItalyPowerDemand/ItalyPowerDemand_TRAIN');
tst = load('../../UCR/ItalyPowerDemand/ItalyPowerDemand_TEST');

print_fig = @(s)  export_fig('-pdf','-transparent',sprintf('dtw_movie/dtw_movie_%06d.pdf',s));

x = trn(1,2:end)';
y = trn(2,2:end)';
n = length(x);

d = (repmat(y,1,n)-repmat(x',n,1)).^2;
d = d/max(max(d));

[~,p1,p2]=tscu_dtw(x,y,length(x));
p1 = fliplr(p1);
p2 = fliplr(p2);

c = zeros(n+1,n+1); 
c(:,1)=inf;
c(1,:)=inf;
c(1,1)=0;

tr = zeros(n+1,n+1);

for i=2:n+1
    for j=2:n+1
        [mval,midx]=min([c(i-1,j-1),c(i-1,j),c(i,j-1)]);
        c(i,j)=d(i-1,j-1)+mval;
        tr(i,j)=midx;
    end
end
c = c/max(max(c(2:end,2:end)));
u = tr(2:end,2:end);
v = tr(2:end,2:end);
u(u==1)=-1; v(v==1)= -1; % from south-west
u(u==2)= 0; v(v==2)= -1; % from west
u(u==3)=-1; v(v==3)=  0; % from south
[xx,yy]=meshgrid(1:n,1:n);

% figure
% imagesc(d)
% set(gca,'YDir','normal');
% colormap gray
% 
% figure
% imagesc(c(2:end,2:end));
% set(gca,'YDir','normal');
% colormap gray
% 
% figure
% plot(x,'k');
% hold on;
% plot(y,'k-.');
% legend('x','y');
% 
% figure
% plot(interp1(1:n,x,p1),'k');
% hold on;
% plot(interp1(1:n,y,p2),'k-.');
% legend('x dtw','y dtw');



rechandles = zeros(n,n);
figure('Position', [0, 0, 1024, 768]);
axis equal
box on
xlim([-4 n+1]);ylim([-4 n+1]);
hold on
for i=1:n
    h=plot(1:i,x(1:i)-2,'k'); 
    g=plot(y(1:i)-2,1:i,'k');
    print_fig(movie_counter); movie_counter = movie_counter+1;
    delete(h);
    delete(g);
end
plot(1:n,x-2,'k');
plot(y-2,1:n,'k');
print_fig(movie_counter); movie_counter = movie_counter+1;


for i=1:n
    for j=1:n
        rechandles(i,j)=rectangle('Position',[i-0.5,j-0.5,1,1],...
          'FaceColor',1-d(i,j)*[1 1 1]);
        print_fig(movie_counter); movie_counter = movie_counter+1;
    end
end

rectangle('Position',[0-0.5,0-0.5,1,1],'FaceColor',[1 1 1]);
print_fig(movie_counter); movie_counter = movie_counter+1;


for i=1:n
    rectangle('Position',[0-0.5,i-0.5,1,1],'FaceColor',[0 0 0]);
    rectangle('Position',[i-0.5,0-0.5,1,1],'FaceColor',[0 0 0]);
    print_fig(movie_counter); movie_counter = movie_counter+1;
end

arrows=zeros(n,n);
for i=1:n
    for j=1:n
        southwest = rectangle('Position',[i-1-0.5,j-1-0.5,1,1],'EdgeColor','green');
        west = rectangle('Position',[i-1-0.5,j-0.5,1,1],'EdgeColor','green');
        south = rectangle('Position',[i-0.5,j-1-0.5,1,1],'EdgeColor','green');
        
        arrow_southwest = arrow([i-1,j-1],[i,j],5,'BaseAngle',90,'TipAngle',10,'FaceColor','g','EdgeColor','g');
        arrow_west      = arrow([i-1,  j],[i,j],5,'BaseAngle',90,'TipAngle',10,'FaceColor','g','EdgeColor','g');
        arrow_south     = arrow([i  ,j-1],[i,j],5,'BaseAngle',90,'TipAngle',10,'FaceColor','g','EdgeColor','g');        

        print_fig(movie_counter); movie_counter = movie_counter+1;

        set(rechandles(i,j),'FaceColor',1-c(i+1,j+1)*[1 1 1]);  
        
        switch tr(i+1,j+1)
            case 1
                delete(arrow_west);
                delete(arrow_south);
                arrows(i,j)=arrow_southwest ;
                set(southwest,'FaceColor','r');
            case 2
                delete(arrow_southwest);
                delete(arrow_south);
                arrows(i,j)=arrow_west ;
                set(west,'FaceColor','r');
            case 3
                delete(arrow_southwest);
                delete(arrow_west);
                set(south,'FaceColor','r');
                arrows(i,j)=arrow_south;
        end
        set(arrows(i,j),'FaceColor','r','EdgeColor','r');
        print_fig(movie_counter); movie_counter = movie_counter+1;

        delete(southwest);
        delete(south);
        delete(west);
        
        
    end
end



%%
for k=1:30
    alpha=k*pi/30;
    rot = [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];
    % Delete arrows
    for i=1:n
        for j=1:n
            delete(arrows(i,j));
            switch tr(i+1,j+1)
                case 1
                    start=[i-1 j-1]';
                    stop=[i   j]';
                    
                case 2
                    start=[i-1 j]';
                    stop=[i   j]';
                case 3
                    start=[i j-1]';
                    stop=[i   j]';
            end
            coord = [start stop];
            middle=mean(coord,2);
            ncoord = rot*(coord-[middle middle])+[middle middle];
            arrows(i,j) = arrow(ncoord(:,1)',ncoord(:,2)',5,'BaseAngle',90,'TipAngle',10,'FaceColor','r','EdgeColor','r');
        end
    end
    
print_fig(movie_counter); movie_counter = movie_counter+1;

end

%%
plen=length(p2);
for i=plen:-1:1
    h=plot(p2(i:end),p1(i:end),'LineWidth',4);
    print_fig(movie_counter); movie_counter = movie_counter+1;
    delete(h)
end
h=plot(p2,p1,'LineWidth',4');
print_fig(movie_counter); movie_counter = movie_counter+1;

% figure
% quiver(xx,yy,u,v);
% hold on;
% plot(p1,p2,'r');


% subplot(121);
% plot(trn(trn(:,1)==1,2:end)','k');
% hold on;
% plot(tst(tst(:,1)==1,2:end)','r');
% 
% subplot(122);
% plot(trn(trn(:,1)==2,2:end)','k');
% hold on;
% plot(tst(tst(:,1)==2,2:end)','r');
