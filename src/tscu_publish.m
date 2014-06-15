%%
% Publish tscu_tutorial.m in PDF and HTML format for displaying on the web
% site.
publish('tscu_tutorial.m','format','pdf' ,'outputDir','../doc/pdf');
publish('tscu_tutorial.m','format','html','outputDir','../doc/html');
%publish('tscu.m','format','html','outputDir','../doc/html');
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