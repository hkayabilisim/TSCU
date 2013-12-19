function out = tscu(x,y,varargin)
%TSCU Time Series Classification Utility
%   TSCU(X) first divides the set X into training and testing sets
%   randomly, then classifies the time series in the testing set by using
%   the time series in the training set  with default classification
%   method (K-NN) and distance metric (Euclidean). The rows of X
%   corresponds to the time series. The first column defines the labels of
%   the time series.
%
%   TSCU(X,Y) classifies the time series in the testing set Y by using
%   the time series in the training set X with default classification
%   method (K-NN) and distance metric (Euclidean). The rows of both X and Y
%   corresponds to the time series in the training and testing sets,
%   respectively. The first column of both X and Y should be the labels of
%   the time series.
%
%   TSCU(X,Y,'option1',value1,'option2',value2,...) classifies 
%   the time series in Y by using training set X by using the options. 
%   Available options are:
%
%   'Classifier': The preferrred classification technique
%    'K-NN'   : K Nearest Neighbor
%    default  : 'K-NN'
%
%   'Alignment': Alignment method
%    'None'   : no alignment
%    'DTW'    : Dynamic Time Warping
%    'CDTW'   : Constained Time Warping
%    'SAGA'   : Signal Alignment via Genetic Algorithm
%    default  : 'None'
%
%   'LogLevel': Log level
%    'Emergency'  : (level 0)
%    'Alert'      : (level 1)
%    'Critical'   : (level 2)
%    'Error'      : (level 3)
%    'Warning'    : (level 4)
%    'Notice'     : (level 5)
%    'Info'       : (level 6) 
%    'Debug'      : (level 7)
%    default      : 'Info'
%
%   'SAGACostFunction': Cost function used in SAGA
%    'Jcost0'     : Euclidean distance of x to y(w(t))
%    'Jcost1'     : MEX counterpart of Jcost0
%    default      : 'Jcost0'
%
%   'SAGAOptimizationMethod': Optimization technique used in SAGA
%    'GA'         : Genetic Algorithm
%    'Simplex'    : Nelder-Mead Simplex method (fminsearch of MATLAB)
%    default      : 'GA'
%
%   'SAGABaseLength': The number of B-spline bases in ODE
%    default      : 8
%
%   'SAGAInitialSolution': Initial solution in SAGA
%    default      : zero vector with length SAGABaseLength.
%
%   'CrossValidation': An integer specifying how many time the 
%    cross validation takes place.
%    default      : 0 means don't do cross validation
%
%   'MATLABPool': MATLAB pool used for parallel computing
%    'local' : Set it it 'local' if you want to use the processors
%              available in local PC.
%    default : '' no pool for parallel computing.
%   
%   'reportLineWidth': Line width of report lines. Actually it defines
%    the width of the first part of the lines.
%    default : 60
%
%   'trainingRatio': Ratios of training set to the whole set of training
%    and testing. Defined between 0 and 1.
%    default : 0.30
%
%   'DTWbandwidth': It is the width of the Sakoe-Chiba band defined in
%    percentage. Setting it to 100 is the same effect as running DTW.
%    default: 6
%
%   'DisplayInputData': Display both training and testing data
%    'yes'          : plot the input data
%    'no'           : don't display!
%    default : 'no'
%
%   'DisplayAlignment': Display alignment for the specified instances
%    defined as {trnidx, tstidx}
%    default : {[],[]} (means no alignment is displayed)
%    example : tscu(trn,tst,'DisplayAlignment',{[1 3],[1]}) will display the 
%    alignment of test sample 1 to the training samples 1 and 3, so 
%    two alignments will be displayed.
%    
%    The figures will be saved in PDF on the current directory.
%
%   'DumpDistanceMatrix': Dump the distance matrix ta a txt file.
%    'yes' : dump it
%    'no'  : don't!
%    default : 'no'
%
%   Z = TSCU(...) returns output values in the structure Z.
%

% I'm adding this small library of Oliver Woodford to produce PDFs 
% ready for publication It's license allows me to include it. 
% See lib/export_fig/license.txt
addpath('lib/export_fig');

options = getDefaultOptions;
if nargin == 0
    error('tscu:noinput','Not enough input arguments.');
elseif nargin == 1
    [x,y]=divideset(x,options);
    displine('Info','Size of input set',sprintf('%d',size(x,1)),options);
    displine('Info','Dividing input into trn/tst','done',options);
elseif nargin == 2,
    if size(x,2) ~= size(y,2)
        error('tscu:invalidlength',...
	 'Length of time series in training and testing sets should be equal');
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
            case 'SAGACostFunction'
                options.SAGACostFunction = varargin{i+1};
            case 'SAGAOptimizationMethod'
                options.SAGAOptimizationMethod = varargin{i+1};
            case 'SAGABaseLength'
                options.SAGABaseLength = varargin{i+1};
            case 'SAGAInitialSolution'
                options.SAGAInitialSolution = varargin{i+1};
            case 'CrossValidation'
                options.CrossValidation = varargin{i+1};
            case 'DisplayInputData'
		        options.DisplayInputData = varargin{i+1};    
            case 'DisplayAlignment'
                options.DisplayAlignment = varargin{i+1};
            case 'DumpDistanceMatrix'
                options.DumpDistanceMatrix = varargin{i+1};
        end
    end
end

% Opening MATLAB pool
if ~isempty(options.MATLABPool)
    displine('Info','Setting parallel processing','',options);
    matlabpool('open',options.MATLABPool);
end

if options.SAGABaseLength ~= length(options.SAGAInitialSolution)
    displine('Warning','Size of initial solution should be',...
     sprintf('%d',length(options.SAGABaseLength)),options);
end

% Check the optimization method
switch options.SAGAOptimizationMethod
    case 'GA'
    case 'Simplex'
    otherwise
        displine('Warning',sprintf('The method "%s" is not recognized',...
         options.SAGAOptimizationMethod),'GA will be used',options);
        options.SAGAOptimizationMethod = 'GA';
end

% Display some debug
displine('Info','Size of training set',sprintf('%d',size(x,1)),options);
displine('Info','Size of testing set',sprintf('%d',size(y,1)),options);
displine('Info','Time series length',sprintf('%d',size(x,2)-1),options);

displine('Info','Classification method',options.classifier,options);
displine('Info','Alignment method',options.alignment,options);
displine('Info','Displaying input data',options.DisplayInputData,options);

if options.CrossValidation < 1
    displine('Info','No cross validation is chosen',...
        sprintf('%d',options.CrossValidation),options);
else
    displine('Info','Cross validation is not implemented',...
        sprintf('%d',options.CrossValidation),options);
end
if strcmp(options.alignment,'SAGA')
    displine('Info','SAGA number of spline bases',...
     sprintf('%d',options.SAGABaseLength),options);
    displine('Info','SAGA optimization method',...
     options.SAGAOptimizationMethod,options);
    displine('Info','SAGA initial solution',...
     sprintf('%5.2f ',options.SAGAInitialSolution),options);
    displine('Info','SAGA cost function',options.SAGACostFunction,options);
end
if strcmp(options.alignment,'CDTW')
    displine('Info','DTW band width (%)',...
     sprintf('%5.2f',options.DTWbandwidth),options);
end
if ~isempty(options.MATLABPool)
    displine('Info','MATLAB Pool',options.MATLABPool,options);
end

if numel(options.DisplayAlignment{1}) > 0 && numel(options.DisplayAlignment{2}) > 0
    displine('Info','Displaying alignments (trn)',...
        sprintf('%d',options.DisplayAlignment{1}),options);
    displine('Info','Displaying alignments (tst)',...
        sprintf('%d',options.DisplayAlignment{2}),options);
else
    displine('Info','Displaying alignments','none',options);
end
displine('Info','Dumping distance matrix',options.DumpDistanceMatrix,options);


% Displaying Input Data
if strcmp(options.DisplayInputData,'yes')
    displayInputData(x,y,options);
end
displayClassInfo(x,y,options);


% Classification
tic
switch options.classifier
    case 'K-NN'
        labels = nnclassifier(x,y,options);
    otherwise
        labels = nnclassifier(x,y,options);
end
classification_time = toc;

% Performance
perf = performance(y(:,1),labels);
displine('Info','Overall Accuracy',sprintf('%-8.3f',perf.OA),options);
displine('Info','Overall Error',sprintf('%-8.3f',perf.error),options);
displine('Info','Producer Accuracy',sprintf('%-8.3f',perf.PA),options);
displine('Info','User Accuracy',sprintf('%-8.3f',perf.UA),options);
displine('Info','Kappa',sprintf('%-8.3f',perf.kappa),options);
displine('Info','Z-value',sprintf('%-8.3f',perf.Z),options);
displine('Info','Confusion matrix',sprintf('\n%s',perf.confmatdisplay),options);
displine('Info','Classification time (sec)',...
 sprintf('%-8.2f',classification_time),options);

% Closing the MATLAB pool if paralel process
if ~isempty(options.MATLABPool)
    displine('Info','Setting parallel processing','',options);
    matlabpool close
end

% Returning output
out.labels = labels;
displine('Info','The end of TSCU','FINISHED',options);
end

function displayClassInfo(x,y,options)
% Display class information

uniquelabels=unique([x(:,1);y(:,1)]);
% For each class
for i=1:length(uniquelabels)
    trnidx = uniquelabels(i)==x(:,1); 
    tstidx = uniquelabels(i)==y(:,1); 
    displine('Info','Class information',...
        sprintf('%d [TRN:%3d TST:%3d]',uniquelabels(i),sum(trnidx),sum(tstidx)),options);

end
end

function displayInputData(x,y,options)
% Display training and testing sets.

uniquelabels=unique([x(:,1);y(:,1)]);
% For each class
for i=1:length(uniquelabels)
    trnidx = uniquelabels(i)==x(:,1); 
    tstidx = uniquelabels(i)==y(:,1); 
    
    
    figure
    subplot(121);
    plot(x(trnidx,2:end)','k');
    title(sprintf('Class index: %d [TRN:%d]',...
     uniquelabels(i),sum(trnidx)));
    subplot(122);
    plot(y(tstidx,2:end)','k');
    title(sprintf('Class index: %d [TST:%d]',...
     uniquelabels(i),sum(tstidx)));
end
end

function options = getDefaultOptions(varargin)
% Please see help tscu for available options.

options.classifier               = 'K-NN';
options.alignment                = 'None';
options.loglevel                 = 'Info';
options.reportLineWidth          = 40;
options.trainingRatio            = 0.3;
options.DTWbandwidth             = 6;
options.MATLABPool               = '';
options.SAGACostFunction        = 'Jcost0';
options.SAGAOptimizationMethod  = 'GA';
options.SAGABaseLength          = 8;
options.SAGAInitialSolution     = zeros(1,options.SAGABaseLength);
options.DisplayInputData        = 'no';
options.CrossValidation         = 0;
options.DisplayAlignment        = {[],[]};
options.DumpDistanceMatrix      = 'no';
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
    out=sprintf('%s',k);
    for i=length(k):o.reportLineWidth
        out=sprintf('%s.',out);
    end
    out=sprintf('%s: %s',out,v);
    fprintf('%s\n',out)
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
function [trn,tst]=divideset(x,options)
%DIVIDESET Divides a set of time series into two parts.
%   DIVIDESET(X,OPTIONS) divides the dataset X into training (X) and testing
%   sets (Y) by using the ratios OPTIONS.trainingRatio.
%
%
m=size(x,1);
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
ylabels = y(:,1);
n = size(x,1); % training
m = size(y,1); % testing
[yIdx,xIdx] = meshgrid(1:m,1:n);
alldistances = zeros(1,n*m);

Alignment = options.alignment;
DisplayAlignmentMat = zeros(n,m);
DisplayAlignmentMat(options.DisplayAlignment{1},options.DisplayAlignment{2})=1;

% It may seem awkward not to use two inner loops for a simple
% K-NN classifier. The reason is that I should use just one loop 
% to easily distribute the job. 
parfor i = 1 : n*m
        yObject = y(yIdx(i),2:end); 
        xObject = x(xIdx(i),2:end);
        path1 = 1:numel(xObject); 
        path2 = 1:numel(yObject);
        switch Alignment
            case 'None'
                alldistances(i) = sqrt(sum((xObject - yObject).^2)); 
            case 'DTW'
                [alldistances(i), path1, path2] = dtwalignment(xObject,yObject);
            case 'CDTW'
                [alldistances(i), path1, path2] = cdtwalignment(xObject,yObject,options);
            case 'SAGA'
                [alldistances(i), path1, path2] = sagaalignment(xObject,yObject,options);
            otherwise
                alldistances(i) = sqrt(sum((xObject - yObject).^2));
        end
        displine('Debug',sprintf('[%5d of %5d] dist(%4d,%4d)',i,n*m,yIdx(i),xIdx(i)),...
            sprintf('%f',alldistances(i)),options);
        if DisplayAlignmentMat(i)
            displayAlignment(x,y,xIdx(i),yIdx(i),path1,path2,Alignment,options)
        end
end
distancematrix = reshape(alldistances,n,m);
[dummy, mindistanceIdx]=min(distancematrix);
labels = xlabels(mindistanceIdx);
displine('Info','index of testing objects',sprintf('%3d ',1:m),options);
displine('Info','labels of testing objects (True)',sprintf('%3d ',ylabels),options);
displine('Info','labels of testing objects (Estimated)',sprintf('%3d ',labels),options);
displine('Info','closest training objects',sprintf('%3d ',mindistanceIdx),options);

if strcmp(options.DumpDistanceMatrix,'yes')
    save(sprintf('tscu_distancematrix_%s.txt',Alignment),'distancematrix','-ascii');
end

end

function  displayAlignment(x,y,xIdx,yIdx,path1,path2,Alignment,options)
%DISPLAYALIGNMENT Displays the aligned time series
%   DISPLAYALIGNMENT(x,y,xIdx,yIdx,path1,path2,Alignment,options) will
%   display the alignment between the time series x(xIdx,2:end))
%   and y(yIdx,2:end)) by using the determined warpings path1 and path2
%   found by the alignment method Alignment.
    xObject = x(xIdx,2:end);
    yObject = y(yIdx,2:end);
    xAligned = xObject(path1);
    yAligned = yObject(path2);
    
    nx = numel(xObject);
    xmin = min(xObject);
    xmax = max(xObject);

    ny = numel(yObject);
    ymin = min(yObject);
    ymax = max(yObject);
    
    dx = 10;
    if nx < 30
        dx = nx;
    end
    xgrid=round(linspace(1,nx,dx));
    ygrid=round(linspace(1,ny,dx));
    
    relErrBefore = norm(xObject - yObject)/norm(xObject);
    relErrAfter  = norm(xAligned - yAligned)/norm(xAligned);
    
    figure('Visible','on');
    plot(xObject,'b');
    hold on;
    plot(yObject,'r');
    %legend(sprintf('TRN %d',xIdx),sprintf('TST %d',yIdx));
    legend('Training','Testing');
    
    title(sprintf('Originals Distance: %8.5f',relErrBefore));
    export_fig('-pdf','-transparent',...
        sprintf('tscu_alignment_%03d_%03d_%s_before.pdf',xIdx,yIdx,Alignment));

    figure('Visible','on');
    plot(xAligned,'b');
    hold on;
    plot(yAligned,'r');
    %legend(sprintf('TRN %d',xIdx),sprintf('TST %d',yIdx));
    legend('Training','Testing');
    title(sprintf('Alignment (%s) Distance: %8.5f',Alignment,relErrAfter)); 
    export_fig('-pdf','-transparent',...
        sprintf('tscu_alignment_%03d_%03d_%s_after.pdf',xIdx,yIdx,Alignment));
    
    figure('Visible','on');
    plot(path1,path2,'b');
    hold on
    plot(nx*0.20*((xObject-xmin)/(xmax-xmin)-1),'r');
    plot(ny*0.20*((yObject-ymin)/(ymax-ymin)-1),1:ny,'r');
    set(gca,'XTick',[],'YTick',[])
    box on
    plot(xgrid,meshgrid(xgrid,ygrid),'k:');
    plot(meshgrid(xgrid,ygrid),ygrid,'k:');
    
    
    if strcmp(Alignment,'CDTW')
        band=floor(options.DTWbandwidth*nx/100);
        plot([1 1    nx-band nx],[1 band ny      ny],'k');
        plot([1 band nx      nx],[1 1    ny-band ny],'k');
    end
    
    axis equal
    xlim([-nx*0.20 nx+1]);
    ylim([-ny*0.20 ny+1]);
    export_fig('-pdf','-transparent',...
        sprintf('tscu_alignment_%03d_%03d_%s_warping.pdf',xIdx,yIdx,Alignment));
    
    figure('Visible','on');
    for i=round(linspace(1,length(path1),50))
        plot([path1(i) path2(i)], [xObject(path1(i)) yObject(path2(i))+2*xmax],'k'); hold on 
    end
    plot(xObject,'b');
    plot(2*xmax+yObject','b');
    xlim([0 nx]);
    set(gca,'XTick',[],'YTick',[])
    box on
    export_fig('-pdf','-transparent',...
        sprintf('tscu_alignment_%03d_%03d_%s_warpingLines.pdf',xIdx,yIdx,Alignment));    
end

function [distance, path1, path2] = sagaalignment(x,y,options)
%SAGAALIGNMENT Signal Alignment via Genetic Algorithm
%   [DISTANCE PATH1 PATH2]=SAGAALIGNMENT(X,Y,OPTIONS) aligns
%   the objects X and Y using available options OPTIONS and returns
%   the distance in DISTANCE and warping paths in PATH1 and PATH2.
        t=linspace(0,1,length(x));

        switch options.SAGACostFunction
            case 'Jcost0'
                J = @(s) Jcost0(y,x,s,t);
            case 'Jcost1'
                J = @(s) Jcost1(y,x,s);
            otherwise
                displine('Warning',sprintf('Cost function "%s" is not defined. Using ',...
                    options.SAGACostFunction),'Jcost0',options);
                J = @(s) Jcost0(y,x,s,t);
        end
        
        switch options.SAGAOptimizationMethod
            case 'GA'
                gaoptions = gaoptimset('Generations',100,...
                    'TolFun', eps, ...
                    'StallGenLimit',40,...
                    'Display','off',...
                    'PopulationSize',20,...
                    'PopInitRange',1*[-1;1]);
                [sbest, distance]=ga(J,options.SAGABaseLength,[],[],[],[],-10,10,[],gaoptions);
            case 'Simplex'
                [sbest, distance] = fminsearch(J,options.SAGAInitialSolution);
            otherwise
                displine('Warning',sprintf('Optimization function "%s" is not defined. Using',...
                    options.SAGAOptimizationMethod),'GA',options);
                [sbest, distance]=ga(J,options.SAGABaseLength);
        end
        path1 = 1:length(x);
        path2 = round(interp1(t,1:length(y),ramsay(t,sbest)));
        path2(path2<1)=1;
        path2(path2>length(y))=length(y);
end

function [distance,path1,path2] = dtwalignment(x,y)
%DTWALIGNMENT Dynamic Time Warping alignment
%   [DISTANCE PATH1 PATH2]=DTWALIGNMENT(X,Y,OPTIONS) aligns
%   the objects X and Y using available options OPTIONS and returns
%   the distance in DISTANCE and warping paths in PATH1 and PATH2.
[distance,path1, path2]=dtw(x,y,length(x));
path1 = fliplr(path1);
path2 = fliplr(path2);
end

function [distance, path1, path2] = cdtwalignment(x,y,options)
%CDTWALIGNMENT Constained Dynamic Time Warping alignment
%   [DISTANCE PATH1 PATH2]=CDTWALIGNMENT(X,Y,OPTIONS) aligns
%   the objects X and Y using available options OPTIONS and returns
%   the distance in DISTANCE and warping paths in PATH1 and PATH2.
[distance, path1, path2]=dtw(x,y,floor(options.DTWbandwidth*length(x)/100));
path1 = fliplr(path1);
path2 = fliplr(path2);
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

