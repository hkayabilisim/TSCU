%% TSCU test suite: 02
% This test runs Time Series Classification Utility (TSCU) with 
% |DisplayInputData| option set to |yes| in order to display the input 
% data. 
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
% Let's create 4 time series with two different classes: sine and cosine.
% We also deviced an artifical change within the same class time series by
% warping the time axis with w(t)=t^2.
%
%   Name  Function       Class  Set
%   ----  --------       -----  --------
%   a     sin(2*pi*t)    1      Training
%   b     sin(2*pi*t*t)  1      Testing
%   c     cos(2*pi*t)    2      Training
%   d     cos(2*pi*t*t)  2      testing
%
% tst and trn vectors contain both the time series and their class labels.
t = linspace(0,1,29);
a = sin(2*pi*t); 
b = sin(2*pi*t.^2);
c = cos(2*pi*t); 
d = cos(2*pi*t.^2);
trn = [ 1 a ; 2 c];
tst = [ 1 b ; 2 d];

%% Running TSCU
% We run TSCU with default settings except |DisplayInputData| option set to 
% yes. In this way, we have the opportunity to observe the input data. You
% can see all the four time series. The first two are the same class
% (since wave) but the second one is sligtly different. Namely the time
% domain of the second time series has warped by using the warping function 
% of t^2. Similary story applied to the other class (cosine wave).
%
% If we had more than one time series for each set (training or testing),
% all of them will be plotted on top of each other.
tscu(trn,tst,'DisplayInputData','yes');

