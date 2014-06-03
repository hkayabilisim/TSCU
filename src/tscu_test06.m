%% TSCU test suite: 06
% In this example, I demonstrator the |LogLevel| option. Setting it to 
% |Debug|,  we can collect more information during classification and 
% alignment. 
%
% * Author : Huseyin Kaya
% * Website: <http://timewarping.org>
% * Sources: <https://github.com/hkayabilisim/TSCU>

%% Initialization
% As always I clear and close everything.
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
% SAGA is implemented in C, so we should compile the corresponding function
% by using mex 
mex tscu_saga_register.c tscu_saga_util.c

%% Running TSCU with increased log level
% We can icrease the log level by setting |LogLevel| property to |Debug|.
% This option turns on verbose messages. Most of the diagnostic messages 
% are self explantory and don't need explanation. I will cover some of
% them:
% 
% * |Class information|: the number of classes are displayed. You will also
% see all the distances between the training testing objects.
% * |index of testing objects|: these are the indexed of testing objects.
% Since we have two testing objects in this example, you are seeing 1 and
% 2.
% * |index of testing objects (True)|: the class labels of testing objects.
% * |index of testing objects (Estimated)|: estimated class labels of
% testing objects. These informations can be useful to create contingency
% matrix and hence statistical significance analysis such as McNemar test.
% * |closest training objects|: For eact test object, these are the index 
% of closest training objects. It can be valueble for deep inspection of
% the performance of classifier.
% 
% It is not recommended to increase the log level especially for relatively 
% large datasets since you will have to deal with a lot of log messages.
tscu(trn,tst,'Alignment','SAGA','LogLevel','Debug');