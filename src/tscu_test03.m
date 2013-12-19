%% Time Series Classification Utility (TSCU) test suite.
% The test runs TSCU with Dynamic Time Warping on a small artifical
% dataset.
%
% * Author : Huseyin Kaya
% * Website: <http://tscu.blogspot.com>
% * Sources: <https://github.com/hkayabilisim/TSCU>
clear all
close all
clc

%% Creating a toy dataset
% Let's create 4 time series with two different classes: 1 and 2. First
% class represents a sine wave, whereas the later represents a cosine wave.
% We also deviced an artifical change within the same class time series by
% warping the time axis with w(t)=t^2.
% 
%
%   Name  Function       Class
%   a     sin(2*pi*t)    1  
%   b     sin(2*pi*t*t)  1
%   c     cos(2*pi*t)    2
%   d     cos(2*pi*t*t)  2
%
% If you have UCR data available, then load it as following:
%
%   trn=load('synthetic_control_TRAIN');
%   tst=load('synthetic_control_TEST');
%
t = linspace(0,1,29);
a=sin(2*pi*t); b=sin(2*pi*t.^2);
c=cos(2*pi*t); d=cos(2*pi*t.^2);
tst = [ 1 a ; 2 c];
trn = [ 1 b ; 2 d];

%% Running TSCU with Dynamic Time Warping
% This time, let's display input data and use DTW as the alignment. Now
% thanks to the alignment you will that there is no misclassification 
% error.
tscu(trn,tst,'Alignment','DTW');
