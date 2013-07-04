%% Time Series Classification Utility (TSCU) test suite.
%
% The test runs TSCU with using the alignment method SAGA.
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

%% Running SAGA with Jcost0
%tscu(trn,tst,'Alignment','SAGA','SAGACostFunction','Jcost0',...
%	     'MATLABPool','anadolu_dual_64');
%% Running SAGA with Jcost1
tscu(trn,tst,'Alignment','SAGA','SAGACostFunction','Jcost1')
