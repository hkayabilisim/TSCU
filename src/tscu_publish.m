%%
% Publish tscu_tutorial.m in PDF and HTML format for displaying on the web
% site.
%publish('tscu_tutorial.m','format','pdf' ,'outputDir','../doc/pdf');
%publish('tscu_tutorial.m','format','html','outputDir','../doc/html');
%publish('tscu.m','format','html','outputDir','../doc/html');
%%
% In order to create the table I mentioned above, I have to use these ugly
% commands. Sorry! 
fprintf('%12s %-17s\n','','Alignment');
fprintf('%12s %s\n','','-----------------');
fprintf('%12s %-5s %-5s %-5s\n','Classifier','NONE','DTW','CDTW');
fprintf('%12s %3.1f%% %3.1f%% %3.1f%%\n','KNN',100*knn_none.perf.OA,...
    100*knn_dtw.perf.OA,100*knn_cdtw.perf.OA);
fprintf('%12s %3.1f%% %3.1f%% %3.1f%%\n','SVM',100*svm_none.perf.OA,...
    100*svm_dtw.perf.OA,100*svm_cdtw.perf.OA);