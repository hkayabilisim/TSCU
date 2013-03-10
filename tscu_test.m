% Time Series Classification Utility (TSCU)
% This script is written to test the basic functionality of TSCU.
%
clear all
close all
clc

%% Synthetic Control dataset from UCR repository.
trn=load('synthetic_control_TRAIN');
tst=load('synthetic_control_TEST');

%% Run with default options.
tscu(trn,tst)

%% Run with Dynamic Time Warping
tscu(trn,tst,'Alignment','DTW')

%% Run with Constrained Dynamic Time Warping
tscu(trn,tst,'Alignment','CDTW','DTWbandwidth',6)


%% Synthetic Control dataset from UCR repository.
trn=load('../UCR/Coffee/Coffee_TRAIN');
tst=load('../UCR/Coffee/Coffee_TEST');

%% Run with SAGA
tscu(trn,tst,...
    'Alignment'  ,'SAGA',...
    'LogLevel'   ,'Debug',...
    'MATLABPool' ,'local',...
    'SAGACostFcn','Jcost0')