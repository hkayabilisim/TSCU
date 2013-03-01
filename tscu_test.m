clear all
close all
clc

% Synthetic Control dataset from UCR repository.
trn=load('synthetic_control_TRAIN');
tst=load('synthetic_control_TEST');

% Run with default options.
tscu(trn,tst)

% Run with Dynamic Time Warping
tscu(trn,tst,'alignment','DTW')

% Run with Constrained Dynamic Time Warping
tscu(trn,tst,'alignment','CDTW','DTWbandwidth',6)