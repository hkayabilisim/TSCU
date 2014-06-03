%% TSCU test suite: 07
% In this example, I use Constrained Dynamic Time Warping (CDTW) as the
% alignment method. CDTW restricts the search space of classical DTW to the
% diagonal of the similarity matrix. The width of this band can be adjusted
% by using |DTWbandwidth| option defined in percentage. Increasing it
% toward 100 percent is equivalent to use classical DTW. Whereas decreasing
% it will result a faster and sometimes more accurate alignment. 
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
% DTW is implemented in C, so we should compile the corresponding function
% by using mex 
mex tscu_dtw.c

%% Running Constrained DTW
% Constrained Dynamic Time Warping (CDTW) is nearly identical to |DTW|
% except that the path in the dynamic programing is confined into a narrow 
% band  along the diagonal of the similarity matrix. This is achieved by 
% using |DTWbandwidth| option. It is defined as percentage. So if you set
% it to 50 then the width of the band will be half width of the time
% series. In order to see the effect, you can use |DisplayAlignment|
% options. If |CDTW| is used as the alignment, the borders of this band
% will be shown in the figures.

%% 30% percent 
% In this example we took %30 of the matrix which is enough to get a good
% alignment. The result will be very similary to the one obtained by using
% |DTW| but faster.
tscu(trn,tst,'Alignment','CDTW','DTWbandwidth',30,'DisplayAlignment',{1,1});

%% 10% percent
% If we further narrow the band, then the alignment is getting worse. Why?
% Because the true warping function (sqrt(t)) is out of the bounds. So 
% there is no way to obtain the real warping function.
tscu(trn,tst,'Alignment','CDTW','DTWbandwidth',10,'DisplayAlignment',{1,1});
