%%
% Publish tscu_tutorial.m in PDF and HTML format for displaying on the web
% site.
%publish('tscu_tutorial.m','format','pdf' ,'outputDir','../doc/pdf');
%publish('tscu_tutorial.m','format','html','outputDir','../doc/html');
%publish('tscu.m','format','html','outputDir','../doc/html');
close all
%% Using Support Vector Machines (SVM) classification
% From the beginning, we have been using Nearest Neighbor classification
% scheme which is the default classifier in TSCU. However you have an
% alternative: Support Vector Machines (SVM). If you want to give it a try,
% you should use the option |Classifier|. The default kernel type of SVM is
% linear with C=10. 
% 
% The implementation of SVM is based on the MATLAB(R) scripts in the book 
% "Support Vector Machines for Antenna Array Processing and 
% Electromagnetics", Manel Martinez-Ramon, Christos G. Christodoulou, 
% Morgan & Claypool Publishers, 2006. 
% ISSN <http://www.amazon.com/dp/159829024X 159829024X>
tscu(trn,tst,'Classifier','SVM');
%%
% Let's analyze the results. First of all it took around 1 seconds to do 
% all the classification whereas KNN required approximately 10 seconds.
% So it is clearly faster than KNN. But remember the exact running times
% will change depending on the number and speed of the processors you have.
% Don't be surprized if you have very different running times. If you have
% a fast CPU, both SVM and KNN may run under a second. If that is the case
% you may not notice a difference between the two. On the other
% hand, if you have a slow PC, then the difference between the two will be
% dramatic.
%
% The next thing to look at is the accuracy. The overall accuracy of SVM
% 92.7% which is higher than KNN (remember that, it was 88%). If you 
% compare the confusion matrices, you will see that SVM is better from KNN 
% at  discriminating the first class from other classes. However it is 
% still having a hard time to separete the first class from the second as 
% it  makes 10 misclassifications (see the 2nd row and first column entry 
% of the confusion matrix). Other than that SVM is make only 4 
% misclassifications.
%% Comparing SVM with KNN under different alignments
% We compared SVM with KNN without using alignment or technically speaking
% alignment with NONE. In this case, SVM seems better than KNN in terms of
% both speed and accuracy. But what about if I use other alignments such as
% DTW and CDTW? I will run TSCU with different options and prepare a table
% similar to the below.
%
%                           Alignment
%                      -------------------
%  Classification      NONE    DTW    CDTW
%  KNN                 ...     ...    ...
%  SVM                 ...     ...    ...
svm_none=tscu(trn,tst,'Classifier','SVM','LogLevel','Alert');
svm_dtw =tscu(trn,tst,'Classifier','SVM','Alignment','DTW',...
          'LogLevel','Alert');
svm_cdtw=tscu(trn,tst,'Classifier','SVM','Alignment','CDTW',...
          'LogLevel','Alert');
knn_none=tscu(trn,tst,'LogLevel','Alert');
knn_dtw =tscu(trn,tst,'LogLevel','Alert','Alignment','DTW');
knn_cdtw=tscu(trn,tst,'LogLevel','Alert','Alignment','CDTW');
%%
% I have the results. Now, in order to create the table I mentioned above, 
% I have to use these ugly looking commands. Sorry! You don't have to use
% these commands but I prefer to create these kind of text tables. 
fprintf('%12s %-17s\n','','Alignment');
fprintf('%12s %s\n','','-----------------');
fprintf('%12s %-5s %-5s %-5s\n','Classifier','NONE','DTW','CDTW');
fprintf('%12s %3.1f%% %3.1f%% %3.1f%%\n','KNN',100*knn_none.perf.OA,...
    100*knn_dtw.perf.OA,100*knn_cdtw.perf.OA);
fprintf('%12s %3.1f%% %3.1f%% %3.1f%%\n','SVM',100*svm_none.perf.OA,...
    100*svm_dtw.perf.OA,100*svm_cdtw.perf.OA);
%%
% As you see, if I don't use alignment, then SVM is definitely better than
% KNN. But using alignment does not improve the accuracy of SVM.  Contrary 
% it slightly decreases the accuracy. So, I can claim that the
% classification accuracy depends not only on the classification method, 
% but also on the alignment method. The two are tighly coupled. A third 
% dependancy is data. The performance of classification algorithms heavily
% depends on the type of data. In the above example, KNN with NONE 
% alignment is better than SVM with NONE alignment. But you find a dataset
% for which the opposite is true. In summary, the classification accuracy
% depends on three main factors:
%
% * classification method (KNN,SVM)
% * alignment method (NONE,DTW,CDTW)
% * data (Synthetic Control,...)
%% Tuning soft margin parameter of SVM
% Like many other machine learning algorithms, SVM has also some parameters
% that controls the details of the algorithm. The most well known parameter
% is the soft margin parameter; usually denoted by C. It is very useful for
% the cases where it is impossible to obtain a separating hyperlane without
% any misclassification. Setting C to a non-zero value, you basically give
% some flexibility to SVM. In doing so, SVM is allowed to make
% some errors. Increasing the soft margin parameter give higher
% flexibility. For further information, you can take a look at 
% Cortes, C.;Vapnik, V. (1995), "Support-vector networks". Machine Learning
% 20 (3): 273 <http://dx.doi.org/10.1007%2FBF00994018
% doi:10.1007/BF00994018>
%%
% The default soft margin parameter is 8 in TSCU. You can change it by
% using |SVMSoftMargin| option. Here I will try 11 different parameters and
% plot the resulting classification accuracies. It is common to try powers
% of 2, so I will follow the general trend.
% The parameters that I used are C={pow(2,i) where i=-5,...,5}
margins=2.^(-5:5);
accuracies=zeros(1,numel(margins));
for i=1:length(margins)
    out=tscu(trn,tst,'Classifier','SVM','SVMSoftMargin',margins(i),...
        'LogLevel','Alert');
    accuracies(i)=out.perf.OA;
end
%%
% Then I will run tscu in silent mode (remember the |LogLevel| option),
% collect the accuracies and plot them. 
figure
plot(log2(margins),100*accuracies,'o-');
xlabel('Soft margin parameter (log2)');
ylabel('Overall accuracy (%)');
%%%
% Here are the classificaton accuracies over different soft margins. As
% you see, the results didn't change except C=1/8 and C=1/4. But the change
% in the accuracy is around 0.006 which is very small and negligible.
% Therefore, changing the soft margin parameter of SVM didn't
% improve the accuracy at least in this dataset and alignment scheme. 
%% Using nonlinear kernel for SVM
% TSCU allows you to change the default linear kernel used in SVM to a
% nonlinear gaussian kernel. For this, you have assign |SVMKernel| option
% to |gaussian|. 
tscu(trn,tst,'Classifier','SVM','SVMKernel','gaussian');
%%
% You can see the kernel type of SVM in the info messages. Other than that,
% the rest of the messages are same as before except the classification
% accuracy: 50% which is very low compared to the linear SVM.
% The confusion matrix tells us that the first class (normal) is confused 
% with the other classes. 
%% Tuning gaussian kernel
% You can change the gaussian kernel parameter sigma by using the option
% |SVMSigma|. Here I will try the following sigma values:
% s={pow(2,i) where i=-5,...,5}. Then I will run tscu in silent mode 
% (remember the |LogLevel| option) collect the accuracies and plot them. 
sigmas=2.^(-5:5);
accuracies=zeros(1,numel(sigmas));
for i=1:length(sigmas)
    out=tscu(trn,tst,'Classifier','SVM','SVMKernel','gaussian',...
        'SVMSigma',sigmas(i),...
        'LogLevel','Alert');
    accuracies(i)=out.perf.OA;
end
%%
figure
plot(log2(sigmas),100*accuracies,'o-');
xlabel('sigma parameter of Gaussian kernel (log2)');
ylabel('Overall accuracy (%)');
%% 
% As you see, the classification accuracy heavily depends of the sigma 
% parameter. As it is being increased, the overall accuracy increases. The
% largest accuracy (the y-axis of the biggest peak) is 97.67\%.  This is
% very high compared to the linear SVM.
fprintf('Maximum accuracy: %8.2f\n',100*max(accuracies));
