function e = tscu(x,y,varargin)
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
%   the time series in Y by using training set X by using the options 
%
%
options = getOptions;
if nargin == 0
    error('tscu:noinput','Not enough input arguments.');
elseif nargin == 1
    [x y]=divideset(x,options);
    if options.displayReportLines
        displine('Size of input set',sprintf('%d',size(x,1)),options);
        displine('Dividing input into trn/tst','done',options);
    end
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
            case 'classifier'
                options.classifier = varargin{i+1};
            case 'alignment'
                options.alignment = varargin{i+1};
        end
    end
end
    
    

% Display some debug
if options.displayReportLines
    displine('Size of training set',sprintf('%d',size(x,1)),options);
    displine('Size of testing set',sprintf('%d',size(y,1)),options);
    displine('Time series length',sprintf('%d',size(x,2)-1),options);
end

% Classification
switch options.classifier
    case '1-NN'
        labels = nnclassifier(x,y,options);
    otherwise
        labels = nnclassifier(x,y,options);
end


% Performance
perf = performance(y(:,1),labels);
if options.displayReportLines
    displine('Overall Accuracy',sprintf('%-8.3f',perf.OA),options);
    displine('Overall Error',sprintf('%-8.3f',perf.error),options);
    displine('Producer Accuracy',sprintf('%-8.3f',perf.PA),options);
    displine('User Accuracy',sprintf('%-8.3f',perf.UA),options);
    displine('Kappa',sprintf('%-8.3f',perf.kappa),options);
    displine('Z-value',sprintf('%-8.3f',perf.Z),options);
    fprintf('%s',perf.confmatdisplay);
end
end


function options = getOptions(varargin)
% Available options are:
%
% classifier: The preferrred classification technique
%   '1-NN'   : 1 Nearest Neighbor
%   default  : '1-NN'
%
% alignment: Alignment method
%   'None'   : no alignment
%   'DTW'    : Dynamic Time Warping
%   default  : 'None'
%
% displayReportLines: Will debugging info lines be displayed?
%   0        : Do not display report lines
%   1        : Display report lines
%   default  : 1
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
options.displayReportLines  = 1;
options.reportLineWidth     = 40;
options.trainingRatio       = 0.3;

if nargin > 0 && mod(nargin,2) ~= 0
    error('tscu:invalidoption','The number of input variables must be even');
end


end

function displine(k,v,o)
%DISPLINE Display a report line
%   DISPLINE(K,V,O) display the string K with value V by using
%   the parameters in O. For example, DISPLINE('Length',12,options)
%   will display
%   Length.......................: 12

if o.displayReportLines
    fprintf('%s',k);
    for i=length(k):o.reportLineWidth
        fprintf('.');
    end
    fprintf(': %s\n',v);
    
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

for i = 1 : n*m
        yObject = y(yIdx(i),2:end);
        xObject = x(xIdx(i),2:end);
        switch options.alignment
            case 'None'
                distance = sqrt(sum((xObject - yObject).^2)); 
            case 'DTW'
                distance = dtwalignment(xObject,yObject,options);
            otherwise
                distance = sqrt(sum((xObject - yObject).^2));
        end
        alldistances(i) = distance;
end
    
[mindistances mindistanceIdx]=min(reshape(alldistances,n,m));
labels = xlabels(mindistanceIdx);
end

function [distance path1 path2] = dtwalignment(x,y,options)
%DTWALIGNMENT Dynamic Time Warping alignment
%   [DISTANCE PATH1 PATH2]=DTWALIGNMENT(X,Y,OPTIONS) aligns
%   the objects X and Y using available options OPTION and returns
%   the distance in DISTANCE and warping paths in PATH1 and PATH2.
[distance path1 path2]=dtw(x,y,length(x));
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