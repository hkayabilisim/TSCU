%% TSCU test suite: 01
% This test runs Time Series Classification Utility (TSCU) in default 
% settings on a small toy example.
%
% * Author : Huseyin Kaya
% * Website: <http://timewarping.org>
% * Sources: <https://github.com/hkayabilisim/TSCU>

%% Initialization
% I prefer to clear and close everything to stay out of any nonsense.
clear all
close all
clc

%% Creating a toy example
H = @(s) sin(4*pi*s);

m = 100;
t = linspace(0,1,m);
z=zeros(1,m);
w=zeros(1,m);
s=zeros(1,4);
sbest=zeros(1,4);

%%
[f1,u1]=tscu_saga_warp(H(t),2*[ 0 1 0 0]);
u1 = (u1-1)/(m-1);

%%
[f2,u2]=tscu_saga_warp(H(t),2*[ 0 -1 0 0]);
u2 = (u2-1)/(m-1);

%%
[u1_saga,u2_saga, d_saga]=tscu_saga_register(f1,f2,4,z,w,s,sbest);
%%
[d_dtw,u1_dtw, u2_dtw]=tscu_dtw(f1,f2,100);
%%
%%
u1_dtw=fliplr(u1_dtw);
u2_dtw=fliplr(u2_dtw);
%% Normalize all warping function to [0 1]

u1_dtw = (u1_dtw-1)/(m-1);
u2_dtw = (u2_dtw-1)/(m-1);
u1_saga = (u1_saga-1)/(m-1);
u2_saga = (u2_saga-1)/(m-1);


dtw_len=length(u1_dtw);
t_dtw=linspace(0,1,dtw_len);

u1_inv=interp1(u1,t,t);
u2_inv=interp1(u2,t,t);

%%
close all
clc

u1ou2 = interp1(t,u1,u2);
u1ou2_inv = interp1(u1ou2,t,t);

u2_invou1 = interp1(t,u2_inv,u1);
d2od1 = interp1(t_dtw,u2_dtw,u1_dtw);

u1od1 = interp1(t,u1,u1_dtw);
u2od2 = interp1(t,u2,u2_dtw);
u1os1 = interp1(t,u1,u1_saga);
u2os2 = interp1(t,u2,u2_saga);

plot(t,H(t)); title('Model function');

figure
plot(t,u1,'k'); hold on
plot(t,u2,'b'); legend('u1','u2');
title('GT Warping functions');

figure
plot(t,f1,'k'); hold on
plot(t,f2,'b'); legend('f1','f2');
title('Observed functions');

figure
plot(t    ,u1_inv ,'r'); hold on;
plot(t_dtw,u1_dtw ,'k');
plot(t    ,u1_saga,'b'); 
legend('GT (u1-)','DTW (v1)','SAGA (v1)');
title('First warping');

figure
plot(t    ,u2_inv,  'r'); hold on;
plot(t_dtw,u2_dtw,'k');
plot(t    ,u2_saga,     'b'); hold on;
legend('GT (u2-)','DTW (v2)','SAGA (v2)');
title('Second warping');

% figure
% plot(t    ,u1ou2_inv,     'r'); hold on;
% plot(t_dtw,d2od1,'k')
% legend('GT (u1 o u2)-','DTW (v1- o v2-)-');
% title('Combination: Sacma bir grafik');

figure
plot(t,u2_invou1,'k'); hold on;
plot(t,u2_saga, 'b'); 
legend('GT (u2- o u1)', 'SAGA (v2)');

figure;
subplot(121);
plot(interp1(t,f1,u1_saga),'k');hold on
plot(interp1(t,f2,u2_saga),'r');
plot(H(t),'b');
legend('g1','g2','H');title('SAGA');
subplot(122);
plot(interp1(t,f1,u1_dtw),'k'); hold on
plot(interp1(t,f2,u2_dtw),'r');
plot(H(t_dtw),'b');
legend('g1','g2','H');title('DTW');

figure
subplot(121);
plot(t_dtw,u1od1,'k'); hold on;
plot(t_dtw,u2od2,'r'); 
plot(t_dtw,t_dtw,'b');legend('u1 o v1','u2 o v2','t');
title('DTW')
subplot(122);
plot(t,u1os1,'k'); hold on;
plot(t,u2os2,'r'); 
plot(t,u1,'b');legend('u1 o v1','u2 o v2','u1');
title('SAGA')

%figure
%for i=1:100
%	[f1 u1]=tscu_saga_warp(H(t),[0 0 10*i/100 0]);
%	plot(1:m,u1);
%	pause(0.1);
%end
