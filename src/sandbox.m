%% Using SVM in TSCU
% Time Series Classification Utility (TSCU) is a collection of MATLAB(R) 
% and C functions written to create an easy to use framework for 
% classification of time series. If you have a collection of time series 
% that needs to be classified by using Support Vector Machined (SVM)
% then continue reading this tutorial.
% 
% * Author : Huseyin Kaya
% * Website: <http://timewarping.org>
% * Sources: <https://github.com/hkayabilisim/TSCU>
%% Installation
% TSCU is written in MATLAB(R), so there is no setup; just download the
% package and run from MATLAB(R) Command Window. The complete package is
% available for free from <https://github.com/hkayabilisim/TSCU GitHub>. 
% You have these two options:
%
% *Option 1* Use <https://github.com/hkayabilisim/TSCU/archive/master.zip Download 
% ZIP> option to download the package in a zip file. If you choose this way
% you have to download the whole package to obtain the most current version
% of the utility.
%
% *Option 2* Another option is to use a command line to fetch the git repository. In
% this way, it is easy to update the package by using suitable |git|
% options. To check out the repository you can use the following command.
% If you don't have git on your command line, then you should install to 
% the operating system. For further information, please take a look at 
% <http://git-scm.com>.
%
%   # git clone https://github.com/hkayabilisim/TSCU.git
%
% In both methods, you will end up with a directory named TSCU. Open your
% MATLAB(R) Command Window, and go to the TSCU/src directory. Now you are 
% ready to run TSCU. But please be patient. Just read this tutorial and
% follow the step by step instructions.
%% Loading a time series dataset
% For this tutorial we will use 
% <http://www.cs.ucr.edu/~eamonn/time_series_data UCR time series
% repository> which contains more than 40 different datasets.
% You can send an e-mail to
% Dr. Keogh to download all of them. For the time being, we will use 
% Gun_Point, ECG200, Lighting2 and yoga. If you 
% haven't already downloaded it, then go ahead and send and e-mail to
% Dr. Keogh and save them under ../../UCR.

%                     RBF-SVM     RBF-DTW    DTW
%         Data        sigma  C    sigma C    bandwith
data = { 'Gun_Point' ,0.03  ,50,  3    ,20,  0,...
    'ECG200'         ,0.1   ,50,  0.2  ,30,  0,...
    'Lighting2'      ,0.001 ,30,  0.01 ,30,  6,...
    'yoga'           ,0.013 ,90,  2    ,20,  2 ...
    };

fprintf('%-15s %-15s %-15s %-4s %-4s\n',...
        'Dataset','RBF_SVM','DTW_SVM','DTW-1NN','1-NN');
for i=1:3
    dataname  = data{6*(i-1)+1};
    s_rbf_svm = data{6*(i-1)+2};
    c_rbf_svm = data{6*(i-1)+3};
    s_dtw_svm = data{6*(i-1)+4};
    c_dtw_svm = data{6*(i-1)+5};
    bandwith  = data{6*(i-1)+6};

    trn=load(sprintf('../../UCR/%s/%s_TRAIN',dataname,dataname));
    tst=load(sprintf('../../UCR/%s/%s_TEST',dataname,dataname));
    
    RBF_SVM=tscu(trn,tst             ,...
        'Classifier'    ,'SVM'       ,...
        'Alignment'     ,'NONE'      ,...
        'SVMKernel'     ,'gaussian'  ,...
        'SVMSoftMargin' ,c_rbf_svm   ,...
        'SVMGamma'      ,s_rbf_svm   ,...
        'LogLevel'      ,'Emergency');
    
    DTW_SVM =tscu(trn,tst            ,...
        'Classifier'    ,'SVM'       ,...
        'Alignment'     ,'CDTW'      ,...
        'DTWbandwidth'  ,bandwith    ,...
        'SVMKernel'     ,'gaussian'  ,...
        'SVMSoftMargin' ,c_dtw_svm   ,...
        'SVMGamma'      ,s_dtw_svm   ,...
        'LogLevel'      ,'Emergency');   
    
    DTW_NN =tscu(trn,tst             ,...
        'Classifier'    ,'KNN'       ,...
        'Alignment'     ,'CDTW'      ,...
        'DTWbandwidth'  ,bandwith    ,...
        'LogLevel'      ,'Emergency');
    
    NN =tscu(trn,tst                 ,...
        'Classifier'    ,'KNN'       ,...
        'Alignment'     ,'NONE'      ,...
        'LogLevel'      ,'Emergency');

    fprintf('%-15s %4.3f (%5.3f %2d) %4.3f (%5.3f %2d) %4.3f %4.3f\n',...
        dataname,RBF_SVM.perf.OA, s_rbf_svm, c_rbf_svm, ...
        DTW_SVM.perf.OA, s_dtw_svm, c_dtw_svm, ...
        DTW_NN.perf.OA,NN.perf.OA);
end