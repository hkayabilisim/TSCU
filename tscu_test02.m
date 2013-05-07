%% Time Series Classification Utility (TSCU) test suite.
%
% This test is to demonstrate the effect of alignment
% in classification accuracy.
%
% * Author : Huseyin Kaya
% * Website: <http://web.itu.edu.tr/huseyinkaya/tscu>
% * Sources: <https://github.com/hkayabilisim/TSCU>

%% Loading data
% I'm using Synthetic Control dataset downloaded from 
% UCR Time Series web site 
% (<http://www.cs.ucr.edu/~eamonn/time_series_data/>). 
% You should find files in the TSCU distribution. 
% If not, then go ahead and download them from UCR web site.
trn=load('synthetic_control_TRAIN');
tst=load('synthetic_control_TEST');

%% Without alignment
% If one uses default options, then the overall accuracy for
% synthetic control dataset is 0.88.
tscu(trn,tst);

%% With alignment
% Now alignment is carried out by using Dynamic Time Warping (DTW) 
% which in turn dramatically increases the classification accuracy.
% The new overall accuracy (0.993) is significantly better than the old 
% accuracy. However the classification time increased.
tscu(trn,tst,'Alignment','DTW');
