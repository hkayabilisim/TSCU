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
% We run with default settings which means, no alignment is applied. 
% You will see many info messages. The confusion matrix is one of the important 
% output of TSCU. In this example you will see a misclassification resulting an 
% overall classification accuracy of 0.5.
tscu(trn,tst);

