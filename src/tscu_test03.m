%% Time Series Classification Utility (TSCU) test suite.
%
% This test is to demonstrate the paralel computing option.
% If 'MATLABPool' option is set to 'local', then classification
% is carried out on available cores in the local computer.
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

%% Running TSCU in serial
% By default, it runs in serial.
tscu(trn,tst,'Alignment','DTW');

%% Running TSCU in parallel
% Running with default options but 'MATLABPool' is set to
% 'local' to speed up the calculation by using the available
% cores in the PC.
tscu(trn,tst,'Alignment','DTW','MATLABPool','local');
