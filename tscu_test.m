% Time Series Classification Utility (TSCU)
% This script is written to test the basic functionality of TSCU.
%
clear all
close all
clc

%% Synthetic Control dataset from UCR repository.
fprintf('--\n');
fprintf('-- Loading synthetic_control dataset\n');
fprintf('--\n');
trn=load('synthetic_control_TRAIN');
tst=load('synthetic_control_TEST');


%% Run with default options.
fprintf('--\n');
fprintf('-- Running with default options\n');
fprintf('--\n');
tscu(trn,tst)

%% Run with Dynamic Time Warping
fprintf('--\n');
fprintf('-- Same options except use DTW in alignment\n');
fprintf('--\n');
tscu(trn,tst,'Alignment','DTW');

%% Run with Constrained Dynamic Time Warping
fprintf('--\n');
fprintf('-- Same options except use DTW with bandwidth of 6\n');
fprintf('--\n');
tscu(trn,tst,'Alignment','CDTW','DTWbandwidth',6);


%% Synthetic Control dataset from UCR repository.
clear all
clc
fprintf('--\n');
fprintf('-- Loading Coffee dataset\n');
fprintf('--\n');
trn=load('../UCR/Coffee/Coffee_TRAIN');
tst=load('../UCR/Coffee/Coffee_TEST');

%% Run with SAGA
fprintf('--\n');
fprintf('-- Using SAGA in alignment with debuggin on\n');
fprintf('--\n');
tscu(trn,tst,...
    'Alignment'  ,'SAGA',...
    'LogLevel'   ,'Info',...
    'SAGAOptimizationMethod','Simplex',...
    'SAGACostFunction','Jcost0');

%% MUFIT'in datasi
clear all
close all
clc
fprintf('--\n');
fprintf('-- Loading MUFIT_SMALL t\n');
fprintf('--\n');
trn=load('../UCR/MUFIT/MUFIT_ALL_TRAIN');
tst=load('../UCR/MUFIT/MUFIT_ALL_TEST');

%% Run with default options.
fprintf('--\n');
fprintf('-- Running with default options\n');
fprintf('--\n');
tscu(trn,tst);
%%
%tscu(trn,tst,'Alignment','CDTW','DTWbandwidth',3,'LogLevel'   ,'Debug');
% Size of training set.....................: 207
% Size of testing set......................: 491
% Time series length.......................: 2151
% No alignment 
% Overall Accuracy.........................: 0.668   
% Overall Error............................: 0.332   
% Producer Accuracy........................: 0.013   0.901   0.949   0.242   0.728   
% User Accuracy............................: 0.125   0.912   0.636   0.889   0.615   
% Kappa....................................: 0.537   
% Z-value..................................: 1.943   
% Confusion matrix.........................: 
% Confusion matrix
%           1     2     3     4     5    UA    TO 
%     1     1     4     0     3     0 0.125     8 
%     2     0    73     0     7     0 0.912    80 
%     3     0     0   112    14    50 0.636   176 
%     4     0     1     0     8     0 0.889     9 
%     5    74     3     6     1   134 0.615   218 
%    PA 0.013 0.901 0.949 0.242 0.728 
%    TO    75    81   118    33   184         491 
% 
% Time elapsed (sec).......................: 567.49  
% CDTW
% Overall Accuracy.........................: 0.670   
% Overall Error............................: 0.330   
% Producer Accuracy........................: 0.027   0.728   0.941   0.576   0.750   
% User Accuracy............................: 0.154   0.843   0.712   0.559   0.633   
% Kappa....................................: 0.545   
% Z-value..................................: 1.985   
% Confusion matrix.........................: 
% Confusion matrix
%           1     2     3     4     5    UA    TO 
%     1     2     9     0     2     0 0.154    13 
%     2     0    59     0    11     0 0.843    70 
%     3     0     0   111     1    44 0.712   156 
%     4     0    13     0    19     2 0.559    34 
%     5    73     0     7     0   138 0.633   218 
%    PA 0.027 0.728 0.941 0.576 0.750 
%    TO    75    81   118    33   184         491 
% 
% Time elapsed (sec).......................: 1143.27 
%%
x=rand(1,1000); y=rand(1,1000); s=[-2 2 -2 2 -2 2];
tic
t=linspace(0,1,length(x));
for i=1:100
    norm(interp1(t,y,ramsay(t,s))-x);
end
fprintf('Jcost0 : %8.2f\n',toc);
tic
for i=1:100
    Jcost1(x,y,s);
end
fprintf('Jcost1 : %8.2f\n',toc);

