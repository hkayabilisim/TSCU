%% TSCU test suite: 03
% This test runs Time Series Classification Utility (TSCU) with |Alignment|
% option set to |DTW|.
%
% * Author : Huseyin Kaya
% * Website: <http://timewarping.org>
% * Sources: <https://github.com/hkayabilisim/TSCU>

%% Initialization
% I usually prefer to clear and close everything not to deal with
% unexpected behaviours.
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
% |tst| and |trn| vectors contain both the time series and their class 
% labels. 
t = linspace(0,1,29);
a = sin(2*pi*t); 
b = sin(2*pi*t.^2);
c = cos(2*pi*t); 
d = cos(2*pi*t.^2);
trn = [ 1 a ; 2 c];
tst = [ 1 b ; 2 d];

%% Building the MEX file
% DTW is implemented in C, so we should compile the corresponding function
% by using mex 
mex tscu_dtw.c

%% Running TSCU with DTW
% We run TSCU with default settings except |Alignment| option set to 
% |DTW|. Dynamic Time Warping (DTW) is a powerfull alignment method dating
% back to 1960s. By using DTW we are expecting an improvement in overall
% classification accuracy.
%
% If you look at the confusion matrix, you will see the improvement: no
% misclassification.
%
tscu(trn,tst,'Alignment','DTW');