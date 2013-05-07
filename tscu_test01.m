%% Time Series Classification Utility (TSCU) test suite.
% The test runs TSCU whose settings are set to default values.
% The data set is Synthetic Control which is provided
% by <http://www.cs.ucr.edu/~eamonn/time_series_data Eamonn Keogh>. 
%
% * Author : Huseyin Kaya
% * Website: <http://web.itu.edu.tr/huseyinkaya/tscu>
% * Sources: <https://github.com/hkayabilisim/TSCU>

%% Loading data
% I'm using the Synthetic Control dataset downloaded from 
% UCR Time Series web site 
% (<http://www.cs.ucr.edu/~eamonn/time_series_data/>).
% The files should be in the TSCU distribution. If not, then go ahead
% and download them from UCR web site.
trn=load('synthetic_control_TRAIN');
tst=load('synthetic_control_TEST');

%% Running TSCU
% TSCU is run with its defaults options. It should run in couple
% of seconds since the alignment is set to 'None' by default. 
% It means standard good old friend Euclidean distance is used to 
% measure how two time series are close to each other. 
tscu(trn,tst);
