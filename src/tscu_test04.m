%% TSCU test suite: 04
% The test runs TSCU with Dynamic Time Warping (DTW) as the alignment method.
% This time we set |DisplayAlignment| option to display the 
% aligned time series as well as the warping function obtained by |DTW|
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
% By setting |Alignment| option to |DTW| we choose DTW as the alignment
% method. In order to use |DTW|, you have to compile the mex file |dtw.c|.
% Please read Installation page on the web site.
%
% In addition to |Alignment|, I also use |DisplayAlignment| option to see
% what is going on during alignment. Since there will be too many pairwise
% alignments during the classification, you should determine exactly which
% alignment you want to see. In this example, the alignment between first
% time series in the training set and the first time series in the testing
% set will be shown.
%
% |DisplayAlignment| option creates four different figures:
%
% * original signals: The original signals are displayed.
% * aligned signals: The aligned signals are displayed. Please note that
% the length of alignment signals will be longer than original signals.
% * warping between the time series: This is a path traversing along the
% main diagonal of similarity matrix. The path depends on the alignment
% method. For instance in |DTW|, the path is found by Dynamic
% Programming (please refer to the User Manual). 
% * mapping between the time series: Actually this is just another way of
% visualization of alignment. In this figure the signals are in their
% original forms but the mapping is displayed by drawing lines. In this way
% you can feel the amount of warping. 
tscu(trn,tst,'Alignment','DTW','DisplayAlignment',{1,1});
