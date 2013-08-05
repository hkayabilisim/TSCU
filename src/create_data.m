%% Create synthetic data

%% Start from scratch
clear all
close all
clc
RandStream.setGlobalStream(RandStream('mt19937ar','seed',4));

%% Model function
f1 = @(s) sin(4*pi*s);
f2 = @(s) sin(6*pi*s);

%% The time domain
nt = 60;
t = linspace(0,1,nt);

%% Warping functions
w1 = @(s) s;%0.2+0.1*rand + s.^(1.8+0.4*rand);
w2 = @(s) s;%0.1*rand + s;

%% Time series
n = 6;
tsClass1 = zeros(n,nt);
tsClass2 = zeros(n,nt);

for i=1:n
    tsClass1(i,:)=f1(w1(t));
    tsClass2(i,:)=f2(w2(t));
end

%% Visualization
m = n/2;
trn = [ones(m,1) tsClass1(1:m,:)    ;2*ones(m,1) tsClass2(1:m,:)];
tst = [ones(m,1) tsClass1(m+1:end,:);2*ones(m,1) tsClass2(m+1:end,:)];

figure;
subplot(221); plot(trn(trn(:,1)==1,2:end)','k'); title('1 trn');
subplot(222); plot(trn(trn(:,1)==2,2:end)','k'); title('2 trn');
subplot(223); plot(tst(tst(:,1)==1,2:end)','k'); title('1 tst');
subplot(224); plot(tst(tst(:,1)==2,2:end)','k'); title('2 tst');

%%
%save('BASIC_TRAIN','trn','-ascii');
%save('BASIC_TEST','tst','-ascii');
%%
tscu(trn,tst,'Alignment','None');
tscu(trn,tst,'Alignment','DTW');
tscu(trn,tst,'Alignment','CDTW');
tscu(trn,tst,'Alignment','SAGA');
