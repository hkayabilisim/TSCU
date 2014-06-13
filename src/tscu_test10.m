%% TSCU test suite: 10
% In this test, we use Synthetic Control dataset provided in UCI Machine
% Learning Repository but I use the one on UCR time series repository 
% (http://www.cs.ucr.edu/~eamonn/time_series_data) since the format is 
% suitable for TSCU. 
%
% In this test, I don't use any alignment algorithm but |NONE| which means
% standard Euclidean distance is used. The classification accuracy of
% this approach is 88\%. It matches the accuracy given on UCR time series
% repository web page.
%
% * Author : Huseyin Kaya
% * Website: <http://timewarping.org>
% * Sources: <https://github.com/hkayabilisim/TSCU>

%% Initialization
% As always, I begin with clearing everything to stay out of any nonsense.
clear all
close all
clc

%% Loading data
% I load synthetic control dataset. In may setup they are in the folder
% '../../UCR/'. You should change this folder fitting your needs.
trn=load('../../UCR/synthetic_control/synthetic_control_TRAIN');
tst=load('../../UCR/synthetic_control/synthetic_control_TEST');

%% Running without alignment
% You don't have to use an alignment method to classify time series as in
% this example. You simply set |Alignment| to |NONE| which in turn uses
% standard Euclidean distance. I also use |DisplayInputData| to see the
% synthetic control dataset.
tscu(trn,tst,'Alignment','NONE','DisplayInputData','yes');
