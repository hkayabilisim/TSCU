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
%    'SVM'    : Linear Support Vector Machines
%    default  : 'K-NN'
%
%   'Alignment': Alignment method
%    'NONE'   : no alignment
%    'DTW'    : Dynamic Time Warping
%    'CDTW'   : Constained Time Warping
%    'SAGA'   : Signal Alignment via Genetic Algorithm
%    'CREG'   : Curve Registration of Ramsay & Silverman
%    default  : 'NONE'
%
%   'SVMKernel': Kernel type of SVM classifier
%    'linear'     : Linear
%    'gaussian'   : Gaussian
%    default      : 'linear'
%
%   'SVMSoftMargin': Soft margin parameter (C) of SVM
%    default      : 8
%
%   'SVMSigma': sigma value in gaussian SVM kernel
%    default      : 1
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
%   'SAGAOptimizationMethod': Optimization technique used in SAGA
%    'GA'         : Genetic Algorithm
%    'GA_MEX'     : A simplified MEX version of Genetic Algorithm
%    'Simplex'    : Nelder-Mead Simplex method (fminsearch of MATLAB)
%    default      : 'GA_MEX'
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
addpath('lib/creg');

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
            case 'SVMKernel'
                options.svmkernel = varargin{i+1};
            case 'SVMSoftMargin'
                options.svmsoftmargin = varargin{i+1};            
            case 'SVMSigma'
                options.svmsigma = varargin{i+1};
            case 'DTWbandwidth'
                options.DTWbandwidth = varargin{i+1};
            case 'LogLevel'
                options.loglevel = varargin{i+1};
            case 'MATLABPool'
                options.MATLABPool = varargin{i+1};
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
    case 'GA_MEX'
    otherwise
        displine('Warning',sprintf('The method "%s" is not recognized',...
            options.SAGAOptimizationMethod),'GA_MEX will be used',options);
        options.SAGAOptimizationMethod = 'GA_MEX';
end

switch options.alignment
    case 'NONE'
        options.alignmentfunction = @nonealignment;
    case 'DTW'
        options.alignmentfunction = @dtwalignment;
    case 'CDTW'
        options.alignmentfunction = @cdtwalignment;
    case 'SAGA'
        options.alignmentfunction = @sagaalignment;
        n = size(x,2)-1;
        k = options.SAGABaseLength;
        options.SAGAz = zeros(1,n);
        options.SAGAw = zeros(1,n);
        options.SAGAs = zeros(1,k);
        options.SAGAsbest = zeros(1,k);
        u = zeros(n,k);
        nk = round(n/k);
        lastpiece = n-nk*k;
        for i=0:k-2
            u(nk*i+1:nk*i+nk,i+1)=1;
        end
        u(n-lastpiece-1:n,k)=1;
        integration    = @(s) cumsum(s)/(size(s,1)-1);
        options.SAGAbmat = integration(integration(u));        
    case 'CREG'
        options.alignmentfunction = @cregalignment;
    otherwise
        options.alignmentfunction = @nonealignment;
end

% Display some debug
displine('Info','Size of training set',sprintf('%d',size(x,1)),options);
displine('Info','Size of testing set',sprintf('%d',size(y,1)),options);
displine('Info','Time series length',sprintf('%d',size(x,2)-1),options);

displine('Info','Classification method',options.classifier,options);
if strcmp(options.classifier,'SVM')
    displine('Info','SVM kernel type',options.svmkernel,options);
    displine('Info','SVM Soft margin',...
        sprintf('%8.5f',options.svmsoftmargin),options);
    if strcmp(options.svmkernel,'gaussian')
        displine('Info','SVM sigma parameter',...
            sprintf('%8.5f',options.svmsigma),options);
    end
end
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
    case 'SVM'
        labels = svmclassifier(x,y,options);
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
out.labels              = labels;
out.truelabels          = y(:,1);
out.classification_time = classification_time;
out.perf                = perf;
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
    subplot(211);
    plot(x(trnidx,2:end)','k');
    title(sprintf('%d time series in training set of class label %d',...
        sum(trnidx),uniquelabels(i)));
    subplot(212);
    plot(y(tstidx,2:end)','k');
    title(sprintf('%d time series in testing set of class label %d',...
        sum(tstidx),uniquelabels(i)));
end
end

function options = getDefaultOptions(varargin)
% Please see help tscu for available options.

options.classifier               = 'K-NN';
options.alignment                = 'NONE';
options.loglevel                 = 'Info';
options.svmkernel                = 'linear';
options.svmsoftmargin            = 8;
options.svmsigma                 = 1;
options.reportLineWidth          = 40;
options.trainingRatio            = 0.3;
options.DTWbandwidth             = 6;
options.MATLABPool               = '';
options.SAGAOptimizationMethod  = 'GA_MEX';
options.SAGABaseLength          = 8;
options.SAGAInitialSolution     = zeros(1,options.SAGABaseLength);
options.DisplayInputData        = 'no';
options.CrossValidation         = 0;
options.DisplayAlignment        = {[],[]};
options.DumpDistanceMatrix      = 'no';
options.alignmentfunction       = @nonealignment;
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
for i = 1 : n*m
    yObject = y(yIdx(i),2:end);
    xObject = x(xIdx(i),2:end);
    path1 = 1:numel(xObject);
    path2 = 1:numel(yObject);
    switch Alignment
        case 'NONE'
            [alldistances(i), path1, path2] = nonealignment(xObject,yObject,options);
        case 'DTW'
            [alldistances(i), path1, path2] = dtwalignment(xObject,yObject,options);
        case 'CDTW'
            [alldistances(i), path1, path2] = cdtwalignment(xObject,yObject,options);
        case 'SAGA'
            [alldistances(i), path1, path2] = sagaalignment(xObject,yObject,options);
        case 'CREG'
            [alldistances(i), path1, path2] = cregalignment(xObject,yObject,options);
        otherwise
            [alldistances(i), path1, path2] = nonealignment(xObject,yObject,options);
    end
    displine('Debug',sprintf('[%5d of %5d] dist(%4d,%4d)',i,n*m,yIdx(i),xIdx(i)),...
        sprintf('%f',alldistances(i)),options);
    if DisplayAlignmentMat(i)
        displayAlignment(x,y,xIdx(i),yIdx(i),path1,path2,Alignment,options)
    end
end
distancematrix = reshape(alldistances,n,m);
[~,mindistanceIdx]=min(distancematrix);
labels = xlabels(mindistanceIdx);
displine('Debug','index of testing objects',sprintf('%3d ',1:m),options);
displine('Debug','labels of testing objects (True)',sprintf('%3d ',ylabels),options);
displine('Debug','labels of testing objects (Estimated)',sprintf('%3d ',labels),options);
displine('Debug','closest training objects',sprintf('%3d ',mindistanceIdx),options);

if strcmp(options.DumpDistanceMatrix,'yes')
    save(sprintf('tscu_distancematrix_%s.txt',Alignment),'distancematrix','-ascii');
end

end

function labels = svmclassifier(x,y,options)
%SVMCLASSIFIER Support Vector Machine Classification
%   LABELS = SVMCLASSIFIER(X,Y,OPTIONS)

classlabels=unique(sort(x(:,1)));
m=length(classlabels);
tstlabels=zeros(size(y,1),(m-1)*m/2);

opt=optimset('maxiter',500,'Algorithm','active-set','display','off');
% One against one approach
k=1;
for i=1:m-1
    for j=i+1:m
        % Training
        class_one=x(x(:,1)==classlabels(i),2:end);
        class_two=x(x(:,1)==classlabels(j),2:end);
        objects=[class_one;class_two];
        truelabels=[-ones(size(class_one,1),1);ones(size(class_two,1),1)];
        displine('Debug','SVM one-against-one',...
            sprintf('class %2d .vs. class %2d',classlabels(i),classlabels(j)),options);
        nx=size(objects,1);
        xlen=size(objects,2);
        R=zeros(nx,nx);
        e=1e-6;
        C=options.svmsoftmargin;    

        switch options.svmkernel
            case 'linear'
                if strcmp(options.alignment,'NONE')
                    R=objects*objects';
                else
                    for ii=1:nx-1
                        for jj=ii+1:nx
                            [distance, path1, path2]=options.alignmentfunction(objects(ii,:),objects(jj,:),options);
                            xAligned = interp1(1:xlen,objects(ii,:),path1);
                            yAligned = interp1(1:xlen,objects(jj,:),path2);
                            displine('Debug','Aligning',sprintf('%5d .vs. %5d',ii,jj),options);
                            R(ii,jj)=xAligned*yAligned';
                            R(jj,ii)=R(ii,jj);
                        end
                    end
                    for ii=1:nx
                        R(ii,ii)=objects(ii,:)*objects(ii,:)';
                    end
                end
            case 'gaussian'
                sigma=options.svmsigma;       % Parameter of the kernel
                D=buffer(sum([kron(objects,ones(nx,1))...
                  - kron(ones(1,nx),objects')'].^2,2),nx,0);
                R=exp(-D/(2*sigma)); % Kernel Matrix               
        end      
        
        Y=diag(truelabels);
        H=Y*R*Y+1e-6*eye(length(truelabels));         %Matrix H regularized
        f=-ones(size(truelabels)); a=truelabels'; K=0;
        Kl=zeros(size(truelabels));
        Ku=C*ones(size(truelabels));
        alpha=quadprog(H,f,[],[],a,K,Kl,Ku,[],opt); %Solver
        
        switch options.svmkernel
            case 'linear'
                w=objects'*(alpha.*truelabels);               
                %ind=find(alpha>e & alpha<=C-e);
                ind=find(alpha<C);
                b=mean(truelabels(ind) - objects(ind,:)*w);
                % Testing
                tstvalues=y(:,2:end)*w+b;
                
                %figure
                %subplot(121)
                %plot(alpha,'.');
                %title(sprintf('%d %d',i,j));
                %subplot(122);
                %plot(eig(H),'.');
                

            case 'gaussian'
                ind=find(alpha>e);
                x_sv=objects(ind,:);    % Extraction of the
                                         % support vectors
                N_SV=length(ind); 
                ind=find(alpha>e & alpha<C-e);
                N_margin=length(ind);
                D=buffer(sum([kron(x_sv,ones(N_margin,1))...
                 - kron(ones(1,N_SV),objects(ind,:)')'].^2,2),N_margin,0);
                R_margin=exp(-D/(2*sigma));
                y_margin=R_margin*(truelabels(ind).*alpha(ind));
                b=mean(truelabels(ind) - y_margin); 
                
                N_test=size(y,1);                
                D=buffer(sum([kron(x_sv,ones(N_test,1))...
                 - kron(ones(1,N_SV),y(:,2:end)')'].^2,2),N_test,0);
                R_test=exp(-D/(2*sigma));
                tstvalues=R_test*(truelabels(ind).*alpha(ind))+b;
        end
        tstlabels(tstvalues< 0,k)=classlabels(i);
        tstlabels(tstvalues>=0,k)=classlabels(j);
        k = k + 1;
    end
end
labels=mode(tstlabels,2);

displine('Debug','index of testing objects',sprintf('%3d ',1:size(y,1)),options);
displine('Debug','labels of testing objects (True)',sprintf('%3d ',y(:,1)),options);
displine('Debug','labels of testing objects (Estimated)',sprintf('%3d ',labels),options);

end

function  displayAlignment(x,y,xIdx,yIdx,path1,path2,Alignment,options)
%DISPLAYALIGNMENT Displays the aligned time series
%   DISPLAYALIGNMENT(x,y,xIdx,yIdx,path1,path2,Alignment,options) will
%   display the alignment between the time series x(xIdx,2:end))
%   and y(yIdx,2:end)) by using the determined warpings path1 and path2
%   found by the alignment method Alignment.
xObject = x(xIdx,2:end);
yObject = y(yIdx,2:end);
xAligned = interp1(1:length(xObject),xObject,path1);
yAligned = interp1(1:length(yObject),yObject,path2);

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
legend('Training','Testing');
title(sprintf('Originals Distance: %8.5f',relErrBefore));
xlim([1 max([nx ny])]);
%export_fig('-pdf','-transparent',...
%    sprintf('tscu_alignment_%03d_%03d_%s_before.pdf',xIdx,yIdx,Alignment));

figure('Visible','on');
plot(xAligned,'b');
hold on;
plot(yAligned,'r');
legend('Training','Testing');
title(sprintf('Alignment (%s) Distance: %8.5f',Alignment,relErrAfter));
xlim([1 length(xAligned)]);
%export_fig('-pdf','-transparent',...
%    sprintf('tscu_alignment_%03d_%03d_%s_after.pdf',xIdx,yIdx,Alignment));

% warping is scaled to [0 1]
figure('Visible','on');
%plot((path1-1)/(length(path1)-1),(path2-1)/(length(path2)-1),'b');
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
%export_fig('-pdf','-transparent',...
%    sprintf('tscu_alignment_%03d_%03d_%s_warping.pdf',xIdx,yIdx,Alignment));

figure('Visible','on');
for i=round(linspace(1,length(path1),50))
    plot([path1(i) path2(i)], ...
        [interp1(1:nx,xObject,path1(i)) interp1(1:ny,yObject,path2(i))+2*xmax],'k');
    hold on
end
plot(xObject,'b');
plot(2*xmax+yObject','b');
xlim([0 nx]);
set(gca,'XTick',[],'YTick',[])
box on
%export_fig('-pdf','-transparent',...
%    sprintf('tscu_alignment_%03d_%03d_%s_warpingLines.pdf',xIdx,yIdx,Alignment));
end

function [distance, path1, path2] = nonealignment(x,y,options)
%NONEALIGNMENT Does nothing but a trivial alignment
%   [DISTANCE PATH1 PATH2]=NONEALIGNMENT(X,Y,OPTIONS) produces
%   a trivial alignment in which warping is a simple line.
path1 = 1:length(x);
path2 = 1:length(y);
distance = sqrt(sum((x - y).^2));
end

function [distance, path1, path2] = sagaalignment(x,y,options)
%SAGAALIGNMENT Signal Alignment via Genetic Algorithm
%   [DISTANCE PATH1 PATH2]=SAGAALIGNMENT(X,Y,OPTIONS) aligns
%   the objects X and Y using available options OPTIONS and returns
%   the distance in DISTANCE and warping paths in PATH1 and PATH2.

J = @(s) tscu_saga_cost(x,y,s,options.SAGAw,options.SAGAz);
switch options.SAGAOptimizationMethod
    case 'GA'
        gaoptions = gaoptimset('Generations',100,...
            'TolFun', eps, ...
            'StallGenLimit',40,...
            'Display','off',...
            'PopulationSize',20,...
            'PopInitRange',1*[-1;1]);
        [sbest, distance]=ga(J,options.SAGABaseLength,[],[],[],[],-2,2,[],gaoptions);
	[~,path2] = tscu_saga_warp(y,sbest);
	path1=1:length(x);
    case 'Simplex'
        [sbest, distance] = fminsearch(J,options.SAGAInitialSolution);
        [~, path2] = tscu_saga_warp(y,sbest);
        path1=1:length(x);
    case 'GA_MEX'
	[path1, path2, distance]=tscu_saga_register(x,y,options.SAGABaseLength,...
    		options.SAGAz,...
    		options.SAGAw,...
    		options.SAGAs,...
    		options.SAGAsbest,...
            options.SAGAbmat');
    otherwise
        displine('Warning',sprintf('Optimization function "%s" is not defined. Using',...
            options.SAGAOptimizationMethod),'GA',options);
	[path1, path2, distance]=tscu_saga_register(x,y,options.SAGABaseLength,...
    		options.SAGAz,...
    		options.SAGAw,...
    		options.SAGAs,...
    		options.SAGAsbest,...
            options.SAGAbmat');
end
end

function [distance, path1, path2] = cregalignment(x,y,options)
%CREGALIGNMENT Curve Registration of Ramsay & Silverman
%   [DISTANCE PATH1 PATH2]=CREGALIGNMENT(X,Y,OPTIONS) aligns
%   the objects X and Y using available options OPTIONS and returns
%   the distance in DISTANCE and warping paths in PATH1 and PATH2.
nx = length(x);
ny = length(y);
nbasis = 10;
norder = 6;
basis = create_bspline_basis([1,nx], nbasis, norder);
xfd = data2fd(x, 1:nx, basis);
yfd = data2fd(y, 1:ny, basis);

nwbasis = 4;
nworder = 3;
wbasis = create_bspline_basis([1,nx], nwbasis, nworder);
Wfd0 = fd(zeros(nwbasis,1), wbasis);

[yfdnew,dummy,warpfd] = registerfd(xfd, yfd, Wfd0);

ynew = eval_fd(yfdnew, 1:ny);
warpvec = eval_mon(1:nx, warpfd);
warpvec = warpvec/max(warpvec);

path1 = 1:nx;
path2 = 1 + warpvec'/max(warpvec)*(nx-1);

distance = sqrt(sum((ynew'-x).^2));
end


function [distance,path1,path2] = dtwalignment(x,y,options)
%DTWALIGNMENT Dynamic Time Warping alignment
%   [DISTANCE PATH1 PATH2]=DTWALIGNMENT(X,Y,OPTIONS) aligns
%   the objects X and Y using available options OPTIONS and returns
%   the distance in DISTANCE and warping paths in PATH1 and PATH2.
[distance,path1, path2]=tscu_dtw(x,y,length(x));
path1 = fliplr(path1);
path2 = fliplr(path2);
end

function [distance, path1, path2] = cdtwalignment(x,y,options)
%CDTWALIGNMENT Constained Dynamic Time Warping alignment
%   [DISTANCE PATH1 PATH2]=CDTWALIGNMENT(X,Y,OPTIONS) aligns
%   the objects X and Y using available options OPTIONS and returns
%   the distance in DISTANCE and warping paths in PATH1 and PATH2.
[distance, path1, path2]=tscu_dtw(x,y,floor(options.DTWbandwidth*length(x)/100));
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
