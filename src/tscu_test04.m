%% Time Series Classification Utility (TSCU) test suite.
% The test runs TSCU in default settings. 
%
% * Author : Huseyin Kaya
% * Website: <http://web.itu.edu.tr/huseyinkaya/tscu>
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



%% Running TSCU
% This time, we want to display a specific alignment between the first
% sample in the training and the first sample in testing.
%
% You will see three figures
%
% * original signals
% * aligned signals
% * warping between the time series.
% * mapping between the time series.
tscu(trn,tst,'Alignment','DTW','DisplayAlignment',{1,1});
