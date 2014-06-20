%% A Tutorial for Time Series Classification Utility (TSCU)
% Time Series Classification Utility (TSCU) is a collection of MATLAB(R) 
% and C functions written to create an easy to use framework for 
% classification of time series. If you have a collection of time series 
% that needs to be classified, then continue reading this tutorial.
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
%% Compiling MEX functions
% In the *TSCU/src* directory, you will see some MEX functions and their
% precompiled versions for Windows XP SP3 and Mac OS X. If someshow you
% need to compile it again, than you should issue the following commands.
% If everything goes smoothly, then you will see the compiled mex files on
% the same directory. If not, you need to configure MEX. A good starting 
% point is 
% <http://www.mathworks.com/help/matlab/matlab_external/what-you-need-to-build-mex-files.html
% here>.
mex tscu_saga_register.c tscu_saga_util.c
mex tscu_saga_warp.c     tscu_saga_util.c
mex tscu_dtw.c
%% Loading a time series dataset
% For this tutorial we will use 
% <http://www.cs.ucr.edu/~eamonn/time_series_data UCR time series
% repository> which contains more than 40 different datasets.
% You can send an e-mail to
% Dr. Keogh to download all of them. For the time being, we will use the 
% only one available for public access: Synthetic Control Dataset. If you 
% haven't already downloaded it, then go ahead and run the following 
% commands to fetch the dataset into the MATLAB workspace.
% ucr_address='http://www.cs.ucr.edu/~eamonn/time_series_data';
% trnfile='synthetic_control_TRAIN';
% tstfile='synthetic_control_TEST';
% urlwrite([ucr_address '/' trnfile],trnfile);
% urlwrite([ucr_address '/' tstfile],tstfile);
% trn=load(trnfile);
% tst=load(tstfile);
% n = size(trn,2)-1;
%% Signal Alignment via Genetic Algorithm
% Now let's go back to the alignment methods. I demonstrated DTW and
% Constrained DTW algorithms in the beginning of this tutorial. TCSU
% support another alignment algorithm called Signal Alignment via Genetic
% Algorithm (SAGA). It combines the idea of using smooth monotone 
% increasing functions introduced by Ramsay with Genetic Algorithm. It
% finds the most suitable warping function in the alignment by solving the
% alignment problem with Genetic Algorithm. For the details, you can
% consult "SAGA: A novel signal alignment method based on genetic 
% algorithm", H Kaya, S Gunduz-Oguducu, Information Sciences 228, 113-130
% <http://www.sciencedirect.com/science/article/pii/S0020025512007955 
% DOI:10.1016/j.ins.2012.12.012>.
%
% In order to use SAGA, you should set the |Alignment| option to |SAGA|.
% Please note that, SAGA is implemented as a MEX file, so you should first
% compile the required MEX files.
mex tscu_saga_register.c tscu_saga_util.c
mex tscu_saga_warp.c     tscu_saga_util.c
%%
% If there is no problem in compiling MEX file then you can run TSCU. Here 
% I use also |DisplayAlignment' option to see the shape of the warping 
% function TSCU created. I choosed the same training and testing samples so
% that you can compare the warping function with the one obtained by 
% Constrained DTW in the previous sections.
tscu(trn,tst,'Alignment','SAGA','DisplayAlignment',{[42],[142]});
%%
% The first thing to take note of is the smoothness of the warping
% function. Since SAGA uses an Ordinary Differential Equation to create 
% smooth functions, you don't see wildly rapid changes or staircase
% patterns in the warping function.  For this reason, SAGA can not make
% local adjustments in the time axis. This feature becomes a liability in
% this synthetic control dataset because the there are random fluctuations
% in the data especially in the first class. The random noise in the time
% series make it harder to discrimate the first class from others. 
%
% Another aspect of SAGA is its speed. As you see from classification
% time message in the output it is way slower than DTW or CDTW. This is 
% an excepted behaviour becase Genetic Algorithm is known to be a slow 
% optimization solver.
