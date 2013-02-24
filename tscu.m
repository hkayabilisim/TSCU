function e = tscu(x,y,options)
%TSCU Time Series Classification Utility
%   TSCU(X,Y) classifies the time series in testing set Y by using 
%   the time series in training set X with default classification 
%   method (1-NN) and distance metric (Euclidean). Rows of both X and Y
%   corresponds to the time series in training and testing sets, 
%   respectively.
%
%   TSCU(X,Y,OPTIONS) classifies the time series in Y by using
%   traininf set X by using the options parameter.

if nargin==2,
  if size(x,2) ~= size(y,2)
    error('tscu:invalidlength','Length of time series in training and testing sets should be equal');
  end
  options = getDefaultOptions;
elseif nargin==3
    
  
end
end

function options = getDefaultOptions(varargin)
%GETDEFAULTOPTIONS  Default options for the Time Series Classification Utility (TSCU)
%   

%% Debugging
% Display warping functions in each iteration
options.displayBestWarping = 0;
% Display best solution in each iteration
options.displayBestSolution = 0;
% Display best variables in each iteration
options.displayBestVariables = 0;
% Display best alignment in each iteration
options.displayBestAlignment =0;
% Display alignments
options.displayGroupAlignment = 0;
% Display data
options.displayData = 0;
% Number of spline bases
options.m = 8   ;
% Minimum value of a design variable
options.mindesign = -30;
% Maximum value of a design variable
options.maxdesign = 10;
% Default population size in genetic algorithm
options.popsize = 20;
% Maximumum number of iterations:
% Valid for fminsearch,and genetic algorithm (max # of generations) 
options.maxiter = 200;
% Initial range for SAGA. Initial population will be
% created with in the range of -val +val
options.initialrange = 0.1;
% Maximum number of cost function evaluations;
options.maxFunEvals = 500;
% Radius used in constrained DTW. If you want to
% repeat the experimets in eamonn time series, then you should
% multiply it with 2. Set it to infinity for unconstrained DTW.
options.radius = 22;
% Run classification algorithm maxRetry times. In each try
% we obtain different results. We take the one with min error.
options.maxRetry = 1;
% Do k-fold cross-validation wherever possible
% Set it to zero to disable it.
options.crossvalidation = 2;
end

