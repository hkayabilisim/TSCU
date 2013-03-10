function out = tscu(x,y,varargin)
%TSCU Time Series Classification Utility
%   TSCU(X) first divides the set X into training and testing set
%   randomly, then classifies the time series in testing set by using
%   the time series in training set  with default classification
%   method (1-NN) and distance metric (Euclidean). Rows of X
%   corresponds to the time series. First column is the labels of
%   integer.
%
%   TSCU(X,Y) classifies the time series in testing set Y by using
%   the time series in training set X with default classification
%   method (1-NN) and distance metric (Euclidean). Rows of both X and Y
%   corresponds to the time series in training and testing sets,
%   respectively. First column of both X and Y should be labels of
%   integer.
%
%   TSCU(X,Y,'option1',value1,'option2',value2,...) classifies 
%   the time series in Y by using training set X by using the options. 
%
%   Z = TSCU(...) returns output values in the structure Z.
%
tic
options = getOptions;
if nargin == 0
    error('tscu:noinput','Not enough input arguments.');
elseif nargin == 1
    [x y]=divideset(x,options);
    displine('Info','Size of input set',sprintf('%d',size(x,1)),options);
    displine('Info','Dividing input into trn/tst','done',options);
elseif nargin == 2,
    if size(x,2) ~= size(y,2)
        error('tscu:invalidlength','Length of time series in training and testing sets should be equal');
    end
end

optarglength = size(varargin,2);
if mod(optarglength,2) ~= 0
    error('tscu:argerror','The number of optional parameters should be even');
else
    for i=1:2:optarglength
        switch varargin{i}
            case 'Classifier'
                options.classifier = varargin{i+1};
            case 'Alignment'
                options.alignment = varargin{i+1};
            case 'DTWbandwidth'
                options.DTWbandwidth = varargin{i+1};
            case 'LogLevel'
                options.loglevel = varargin{i+1};
            case 'MATLABPool'
                options.MATLABPool = varargin{i+1};
            case 'SAGACostFcn'
                options.SAGACostFcn = varargin{i+1};
        end
    end
end
    
% Opening MATLAB pool
if ~isempty(options.MATLABPool)
    displine('Info','Setting parallel processing','',options);
    matlabpool close
    matlabpool('open',options.MATLABPool);
end

% Display some debug
displine('Info','Size of training set',sprintf('%d',size(x,1)),options);
displine('Info','Size of testing set',sprintf('%d',size(y,1)),options);
displine('Info','Time series length',sprintf('%d',size(x,2)-1),options);
displine('Info','Classification method',options.classifier,options);
displine('Info','Alignment method',options.alignment,options);
if strcmp(options.alignment,'CDTW')
    displine('Info','DTW band width (%)',sprintf('%5.2f',options.DTWbandwidth),options);
end
displine('Info','MATLAB Pool',options.MATLABPool,options);
displine('Info','SAGA cost function',options.SAGACostFcn,options);

% Classification
switch options.classifier
    case '1-NN'
        labels = nnclassifier(x,y,options);
    otherwise
        labels = nnclassifier(x,y,options);
end


% Performance
perf = performance(y(:,1),labels);
displine('Info','Overall Accuracy',sprintf('%-8.3f',perf.OA),options);
displine('Info','Overall Error',sprintf('%-8.3f',perf.error),options);
displine('Info','Producer Accuracy',sprintf('%-8.3f',perf.PA),options);
displine('Info','User Accuracy',sprintf('%-8.3f',perf.UA),options);
displine('Info','Kappa',sprintf('%-8.3f',perf.kappa),options);
displine('Info','Z-value',sprintf('%-8.3f',perf.Z),options);
displine('Debug','Confusion matrix',sprintf('\n%s',perf.confmatdisplay),options);
displine('Info','Time elapsed (sec)',sprintf('%-8.2f',toc),options);

% Returning output
out.labels = labels;
end


function options = getOptions(varargin)
% Available options are:
%
% Classifier: The preferrred classification technique
%   '1-NN'   : 1 Nearest Neighbor
%   default  : '1-NN'
%
% Alignment: Alignment method
%   'None'   : no alignment
%   'DTW'    : Dynamic Time Warping
%   'CDTW'   : Constained Time Warping
%   'SAGA'   : Signal Alignment via Genetic Algorithm
%   default  : 'None'
%
% LogLevel: Log level
%   'Emergency'  : (level 0)
%   'Alert'      : (level 1)
%   'Critical'   : (level 2)
%   'Error'      : (level 3)
%   'Warning'    : (level 4)
%   'Notice'     : (level 5)
%   'Info'       : (level 6) 
%   'Debug'      : (level 7)
%   default      : 'Info'
%
% SAGACostFcn: Cost function used in SAGA
%   'Jcost0'     : MATLAB, interp1
%   'Jcost1'     : MEX, LAPACK,
%   default      : 'Jcost0'
%
% MATLABPool: MATLAB pool used for parallel computing
%   default : '' no pool for parallel computing.
%   
% reportLineWidth: Line width of report lines. Actually it defines
%   the width of the first part of the lines.
%   default = 60
%
% trainingRatio: Ratios of training set to the whole set of training
%   and testing. Defined between (0,1).
%   default = 0.30

options.classifier          = '1-NN';
options.alignment           = 'None';
options.loglevel            = 'Info';
options.reportLineWidth     = 40;
options.trainingRatio       = 0.3;
options.DTWbandwidth        = 6;
options.MATLABPool          = '';
options.SAGACostFcn         = 'Jcost0';

if nargin > 0 && mod(nargin,2) ~= 0
    error('tscu:invalidoption','The number of input variables must be even');
end


end

function displine(l,k,v,o)
%DISPLINE Display a report line
%   DISPLINE(L,K,V,O) display the string K with value V by using
%   the parameters in O if current loglevel is less than or equal to
%   logvel K. For example, DISPLINE('Info','Length',12,options)
%   will display
%   Length.......................: 12 
%   if the current log level option is 'Info' or 'Debug'. 
%   See options.loglevel setting.

if getloglevelindex(l) <= getloglevelindex(o.loglevel)
    fprintf('%s',k);
    for i=length(k):o.reportLineWidth
        fprintf('.');
    end
    fprintf(': %s\n',v);
end
end

function i = getloglevelindex(l)
%GETLOGLEVELINDEX  Conversion of log level string.
%   GETLOGLEVELINDEX(L) gets the corresponding integer for a given
%   log level L.
switch l
    case 'Emergency' 
        i = 0;
    case 'Alert'
        i = 1;
    case 'Critical' 
        i = 2;
    case 'Error'
        i = 3; 
    case 'Warning' 
        i = 4;
    case 'Notice'
        i = 5;
    case 'Info' 
        i = 6;
    case 'Debug'
        i = 7; 
    otherwise
        i = 6;
end
end
function [trn tst]=divideset(x,options)
%DIVIDESET Divides a set of time series into two parts.
%   DIVIDESET(X,OPTIONS) divides the dataset X into training (X) and testing
%   sets (Y) by using the ratios OPTIONS.trainingRatio.
%
%
[m n]=size(x);
labels = x(:,1);

% index of which objects will be reserved for training.
trnidx = zeros(m,1);
% If a label is not an integer, exit with error.
floatlabels=find(abs(round(labels)-labels)>eps);
if numel(floatlabels) > 0
    error('tscu:floatinglabels',...
        'Some labels are not integer. Example: label of object [%d] is %f',floatlabels(1),x(floatlabels(1),1));
end
uniquelabels = unique(labels);
for i=1:length(uniquelabels)
    label=uniquelabels(i);
    labelidx = find(label==labels);
    
    trnlength = round(length(labelidx)*options.trainingRatio);
    
    if trnlength < 1
        warning('tscu:divideset','There is no object left to training set for label %d',label);
        trnlength = 0;
    end
    trnidx(labelidx(1:trnlength))=1;
end
trn = x(trnidx==1,:);
tst = x(trnidx==0,:);
end

function labels = nnclassifier(x,y,options)
%NNCLASSIFIER Nearest Neighbor Classification
%   LABELS = NNCLASSIFIER(X,Y,OPTIONS) classifies the time series in 
%   testing set Y by using the time series in training set X with the
%   nearest neighbor algorithm resulting estimated labels LABELS.
xlabels = x(:,1);
n = size(x,1);
m = size(y,1);
[yIdx xIdx] = meshgrid(1:m,1:n);
alldistances = zeros(1,n*m);

% It may seem awkward not to use two inner loops for a simple
% 1-NN classifier. The reason is that I should use just one loop 
% to easily share the jobs. The need is evident
% if the number of available processors are more than max(n,m).
parfor i = 1 : n*m
        yObject = y(yIdx(i),2:end);
        xObject = x(xIdx(i),2:end);
        switch options.alignment
            case 'None'
                distance = sqrt(sum((xObject - yObject).^2)); 
            case 'DTW'
                distance = dtwalignment(xObject,yObject,options);
            case 'CDTW'
                distance = cdtwalignment(xObject,yObject,options);
            case 'SAGA'
                distance = sagaalignment(xObject,yObject,options);
            otherwise
                distance = sqrt(sum((xObject - yObject).^2));
        end
        alldistances(i) = distance;
        displine('Debug',sprintf('[%5d of %5d] dist(%4d,%4d)',i,n*m,yIdx(i),xIdx(i)),...
            sprintf('%f',distance),options);
end
    
[mindistances mindistanceIdx]=min(reshape(alldistances,n,m));
labels = xlabels(mindistanceIdx);
end


function [distance path1 path2] = sagaalignment(x,y,options)
%SAGAALIGNMENT Signal Alignment via Genetic Algorithm
%   [DISTANCE PATH1 PATH2]=SAGAALIGNMENT(X,Y,OPTIONS) aligns
%   the objects X and Y using available options OPTIONS and returns
%   the distance in DISTANCE and warping paths in PATH1 and PATH2.
        gaoptions = gaoptimset('Generations',200,...
            'TolFun', eps, ...
            'StallGenLimit',40,...
            'Display','off',...
            'PopulationSize',20,...
            'PopInitRange',1*[-1;1]);
        switch options.SAGACostFcn
            case 'Jcost0'
                t=linspace(0,1,length(x));
                J = @(s) norm(interp1(t,y,ramsay3(t,s))-x);
            case 'Jcost1'
                J = @(s) Jcost1(y,x,s);
            otherwise
                t=linspace(0,1,length(x));
                J = @(s) norm(interp1(t,y,ramsay3(t,s))-x);
        end
        [sbest distance]=ga(J,8,[],[],[],[],-30,10,[],gaoptions);
end

function [distance path1 path2] = dtwalignment(x,y,options)
%DTWALIGNMENT Dynamic Time Warping alignment
%   [DISTANCE PATH1 PATH2]=DTWALIGNMENT(X,Y,OPTIONS) aligns
%   the objects X and Y using available options OPTIONS and returns
%   the distance in DISTANCE and warping paths in PATH1 and PATH2.
[distance path1 path2]=dtw(x,y,length(x));
end

function [distance path1 path2] = cdtwalignment(x,y,options)
%CDTWALIGNMENT Constained Dynamic Time Warping alignment
%   [DISTANCE PATH1 PATH2]=CDTWALIGNMENT(X,Y,OPTIONS) aligns
%   the objects X and Y using available options OPTIONS and returns
%   the distance in DISTANCE and warping paths in PATH1 and PATH2.
[distance path1 path2]=dtw(x,y,floor(options.DTWbandwidth*length(x)/100));
end

function perf = performance(truelabels,estimatedlabels)
%PERFORMANCE Evaluates the performance of a classifier.
%   P = PERFORMANCE(TRUELABELS,ESTIMATEDLABELS) calculated
%   performance metrics and store them in the ouput P.
uniqtruelabels = unique(sort(truelabels));
uniqestimatedlabels = unique(sort(estimatedlabels));
n = length(uniqtruelabels);
confmat = zeros(n,n);
for i = 1:length(truelabels)
    truelabel = truelabels(i);
    estimatedlabel = estimatedlabels(i);
    
    ii = find(uniqestimatedlabels == estimatedlabel);
    jj = find(uniqtruelabels == truelabel);
    confmat(ii,jj) = confmat(ii,jj)+1;
end

% sum of row, niplus
niplus  = sum(confmat,2);
njplus  = sum(confmat,1);
nn       = sum(njplus);

% overall, producer's and user's accuracy
OA      = trace(confmat)/nn;
PA      = diag(confmat)./njplus';
UA      = diag(confmat)./niplus ;

% proportion of samples
pij     = confmat./nn;
piplus  = sum(pij,2);
pjplus  = sum(pij,1);

po      = trace(pij);
pc      = sum(piplus.*pjplus');
kappa   = (po - pc )/(1-pc);
theta1  = sum(diag(confmat))/nn;
theta2  = sum(niplus.*njplus')/(nn*nn);
theta3  = sum(diag(confmat).*(niplus+njplus'))/(nn*nn);

theta4  = sum(sum(confmat .* (repmat(niplus',size(confmat,1),1) +...
    repmat(njplus,size(confmat,1),1)).^2,2))/(nn^3);
varK    = ((theta1*(1-theta1)/(1-theta2)^2) + ...
    (2*(1-theta1)*(2*theta1*theta2-theta3))/(1-theta2)^3 + ...
    (((1-theta1)^2)*(theta4-(4*(theta2^2)))) /((1-theta2)^2)) / n;
Z       = kappa/sqrt(varK);

% Display stuff
confmatdisplay = '';
confmatdisplay = sprintf('Confusion matrix\n%s      ',confmatdisplay);
for i=1:n
    confmatdisplay = sprintf('%s%5d ',confmatdisplay,uniqtruelabels(i));
end
confmatdisplay = sprintf('%s%5s %5s ',confmatdisplay,'UA','TO');
confmatdisplay = sprintf('%s\n',confmatdisplay);
for i=1:n
    confmatdisplay=sprintf('%s%5d ',confmatdisplay,uniqtruelabels(i));
    for j=1:n
        confmatdisplay=sprintf('%s%5d ',confmatdisplay,confmat(i,j));
    end
    confmatdisplay=sprintf('%s%5.3f %5d ',confmatdisplay,UA(i),sum(confmat(i,:)));
    confmatdisplay=sprintf('%s\n',confmatdisplay);
end
confmatdisplay=sprintf('%s%5s ',confmatdisplay,'PA');
for j=1:n
        confmatdisplay=sprintf('%s%5.3f ',confmatdisplay,PA(j));
end
confmatdisplay=sprintf('%s\n',confmatdisplay);
confmatdisplay=sprintf('%s%5s ',confmatdisplay,'TO');
for j=1:n
        confmatdisplay=sprintf('%s%5d ',confmatdisplay,sum(confmat(:,j)));
end
confmatdisplay=sprintf('%s%5s %5d ',confmatdisplay,'',sum(confmat(:)));
confmatdisplay=sprintf('%s\n',confmatdisplay);

perf.OA = OA;
perf.PA = PA;
perf.UA = UA;
perf.error = 1-OA;
perf.kappa = kappa;
perf.Z = Z;
perf.confmat = confmat;
perf.confmatdisplay = confmatdisplay;

end

function [u c d] = ramsay3(t,w)
% Solution of ODE for a piecewise defined weight
% function
% t 1 x n : evaluation points
% w 1 x m : curvatures
% out:
% u 1 x n : values on t
% c 1 x m : Coefficients: c
% d 1 x m : Coefficients: d

m = length(w);
n = length(t);

A = zeros(2*m,2*m);
y = zeros(2*m,1);
u = zeros(1,n);
b = linspace(0,1,m+1);

% If the variables are off limits then
% we correct the variables.
w(w<-30) = -30;
w(w>10)  = 10;

eps = 0.001;
for i=1:m-1
    if abs(w(i)) < eps
        A(i,i)     = b(i+1); % ok
        A(i+m-1,i) = 1;  %ok
    else
        A(i,i)     = exp(w(i)*b(i+1)); % ok
        A(i+m-1,i) = w(i)*exp(w(i)*b(i+1)); %ok
    end
    
    if abs(w(i+1)) < eps
        A(i,i+1)     = -b(i+1); %ok
        A(i+m-1,i+1) = -1; %ok
    else
        A(i,i+1)     = -exp(w(i+1)*b(i+1)); %ok
        A(i+m-1,i+1) = -w(i+1)*exp(w(i+1)*b(i+1)); %ok
    end
    
    A(i,m+i) = 1; %ok
    A(i,m+i+1) = -1; %ok
end

if abs(w(1)) > eps
    A(2*m-1,1) = 1; % ok
end
A(2*m-1,m+1)=1; % ok
if abs(w(m)) < eps
    A(2*m,m) = 1; %ok
else
    A(2*m,m) = exp(w(m)); %ok
end
A(2*m,2*m)=1;

%fprintf('rcond: %e min: %8.2f max %8.2f\n',rcond(A),min(w),max(w));
% Right hand side of Ax=y
y(2*m) = 1;

x = (A\y); % solution of Ax=y
c = x(1:m)';
d = x(m+1:2*m)';

for i=1:m
    idx = find(t >= b(i) & t <= b(i+1));
    if abs(w(i)) < eps
        u(idx) = c(i)*t(idx)+d(i);
    else
        u(idx) = c(i)*exp(w(i)*t(idx))+d(i);
    end
end
u(1) = 0;
u(end) = 1;


end