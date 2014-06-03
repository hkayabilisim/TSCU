%% TSCU test suite: 09
% In this test we user the curve registration technique
% developed by J. O. Ramsay & X. C. Li published on the 
% following paper:
%
% * 1998. Curve registration, Journal of the Royal Statistical Society 
%   Series  B-statistical Methodology, 60(Part 2), 351?363.
%
% Ramsay published the related codes on
% http://www.psych.mcgill.ca/misc/fda. Since the codes are published under 
% GPL, I include them in TSCU's 'lib' folder so that you don't need 
% a separate download. 
%
% * Author : Huseyin Kaya
% * Website: <http://timewarping.org>
% * Sources: <https://github.com/hkayabilisim/TSCU>

%% Initialization
% As always, I begin with clearing everything to stay out of any nonsense.
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

%% Running Curve Registration
% In this example, we use the curve registraion technique of Ramsay & Li.
% The related keywoard is |CREG| for this purpose. The related software is
% already bundled in TSCU so that you don't need to download them. In
% addition to that, I use |DisplayAlignment| option to see an example
% alignment in action. 
%
% The resulting alignment (registration in Ramsay's 
% terminology) is nearly identical to |SAGA| as expected since they both 
% use the same ODE model. Since the curve registration technique uses a
% deterministic optimization routines, we don't expect difference from run
% to run. This is not the case for |SAGA| since it uses Genetic Algorithm
% which result in slight different results for each run.
tscu(trn,tst,'Alignment','CREG','DisplayAlignment',{1,1});
